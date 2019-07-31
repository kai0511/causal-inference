library(grf)
library(MASS)
library(dplyr)
# library(doMC)
# library(conformalInference)
# library(causalLearning)
# library(devtools)
# library(TCGA2STAT)
# library(SummarizedExperiment)
# library(TCGAbiolinks)
# install_github("saberpowers/causalLearning")
library(data.table)
library(survival)
# library(survminer)
library(mppa)
library(doParallel)
library(methods)
library(randomForest)


# source('/exeh_3/kai/R.code/HTE/causal_inference_models.R')

registerDoParallel(5)

impute.survival <- function(surv.time, censor){
    # @surv.time: a vector of times, if the subject is alive, then it's NA; else, it's a number, indicated survival time
    # @censor: The status indicator, normally 0 = alive, 1 = dead, required by "Surv" function
    # need to change censoring status with the longest censoring to non-censor and assume death

    surv.obj <- Surv(surv.time, censor)
    fit.obj <- survfit(surv.obj ~ 1)
    step.surv <- -diff(c(1, fit.obj$surv))
    log.step <- log(fit.obj$time) * step.surv

    # impute the survival time of censored obs.(censor == 0), unknown to us.
    log.times <- log(surv.time)
    # for(i in 1: length(surv.time)){
    #     ind <- which(fit.obj$time == surv.time[i])
    #     if(cens[i] == 0){
    #         log.times[i] <- sum(log.step[fit.obj$time > fit.obj$time[ind]]) / fit.obj$surv[ind]
    #     }
    # }

    imputed.log.surv.times <- sapply(1:length(surv.time), function(i){
        ind <- which(fit.obj$time == surv.time[i])
        if(censor[i] == 0){
            return(sum(log.step[fit.obj$time > fit.obj$time[ind]]) / fit.obj$surv[ind])
        }else{
            return(log.times[i])
        }
    })
    return(imputed.log.surv.times)
}

extract.freq.mutation <- function(dataset, pos, threshold, is.binary = TRUE ){
    # @dataset mutation dataset
    # @pos the position where mutation variables start
    # @threshold threshold for frequent mutation
    if(is.binary == TRUE){
        mutation.count <- apply(dataset[, pos:ncol(dataset)], 2, function(x) sum(!x == 0))
    }else{
        mutation.count <- apply(dataset[, pos:ncol(dataset)], 2, function(x) sum(x >= 0.75))
    }

    # find mutations with mutation count greater than the defined threshold
    freq.mutations <- colnames(dataset)[pos:ncol(dataset)][mutation.count >= threshold]
    print(paste0('number of mutation will be studied:', sum(mutation.count >= threshold)))
    return(freq.mutations)
}

assess.tau.risk <- function(fitted.obj, Y, treatment){
    # calculate tau-risk_R
    tau.pred <- predict(fitted.obj, newdata = NULL, estimate.variance = TRUE, num.threads = 8)
    Y.hat <- fitted.obj[["Y.hat"]]
    W.hat <- fitted.obj[["W.hat"]]
    tau.1 <- Y - Y.hat

    tau.2 <- (treatment - W.hat) * tau.pred$predictions
    risk <- mean((tau.1 - tau.2)^2)
    return(risk)
}

assess.constant.risk <- function(Y, treatment){
    # calculate constant-risk_R
    trainset <- cbind(Y, treatment)
    rf_fitted <- randomForest(Y ~ treatment, data = trainset, ntree = 2000)
    pred <- predict(rf_fitted)

    syn.trainset <- cbind(Y = Y, treatment = 1 - treatment)
    syn.rf_fitted <- randomForest(Y ~ treatment, data = syn.trainset, ntree = 2000)
    syn.pred <- predict(syn.rf_fitted)

    res <- pred - syn.pred
    tau <- sapply(seq(length(treatment)), function(i) if(treatment[i] == 1) res[i] else -res[i])
    Y.hat <- mean(Y)
    propensity.score <- mean(treatment)
    constant.risk <- mean((Y - Y.hat - (treatment - propensity.score) * tau)^2)
    return(constant.risk)
}

bootstrap.testing <- function(covariates, Y, treatment, tau.risk, constant.risk, tau.var, num.strap = 1000){
    # calculate average difference in mean survival time in treatment v.s non-treatment groups
    # including two subsampling testings, so two pvalues vill be returned in a vector.
    subsampled.result <- foreach(i = seq(num.strap), .combine = 'rbind') %dopar%  {
        samp <- sample(dim(covariates)[1], replace = T)
        sampled.covariates <- covariates[samp, ]
        sampled.Y <- Y[samp]
        sampled.treatment <- treatment[samp]
        
        params = tune_causal_forest(sampled.covariates, sampled.Y, sampled.treatment)$params
        tau.forest <- cf(sampled.covariates, sampled.Y, sampled.treatment, params)  
        bootstrap.risk <- assess.tau.risk(tau.forest, sampled.Y, sampled.treatment)
        
        bootstrap.risk.constant <- assess.constant.risk(sampled.Y, sampled.treatment)

        tau.pred <- predict(tau.forest, newdata = NULL, estimate.variance = TRUE, num.threads = 8)
        bootstrap.tau.var <- var(tau.pred$predictions)

        c(bootstrap.risk, bootstrap.risk.constant, bootstrap.tau.var)
    }

    subsampled.risk <- subsampled.result[, 1]
    subsampled.risk.constant <- subsampled.result[, 2]
    subsampled.var <- subsampled.result[, 3]

    theta.hat <- constant.risk - tau.risk
    theta.star <- subsampled.risk.constant - subsampled.risk
    subsampling.pval <- mean(theta.star - theta.hat >= theta.hat)

    var.pval <- mean(subsampled.var - tau.var >= tau.var)

    return(c(subsampling.pval, var.pval))
}

permutate.covariates.testing <- function(covariates, Y, treatment, tau.risk, num.strap = 1000){
    # permute covariates to investigate whether covariates contribute to HTE
    # return tau.risk for each permutation
    perm.risk <- foreach(i = seq(num.strap), .combine = 'c') %dopar%  {
        samp <- sample(dim(covariates)[1], replace = F)
        sampled.covariates <- covariates[samp, ]

        params = tune_causal_forest(sampled.covariates, Y, treatment)$params
        perm.tau.forest <- cf(sampled.covariates, Y, treatment, params)
        perm.risk <- assess.tau.risk(perm.tau.forest, Y, treatment)
        perm.risk
    }
    return(mean(perm.risk - tau.risk < 0))
}

permutated.pval <- function(Y1, Y2, no.simulations = 10000){
    # test fitting with permutating one outcome variable, and a pval will be returned
    set.seed(0)
    risk <- mean((Y1 - Y2)^2)
    
    permutated.risk <- sapply(seq(no.simulations), function(x){
        mean((Y1[sample(length(Y1))] - Y2)^2)
    })
    return(mean(permutated.risk < risk))
}

cf <- function(covarites, Y, W, params){
    # set and run casual forest
    tau.forest <- causal_forest(X = covarites,
                                Y = Y,
                                W = W, 
                                sample.fraction = as.numeric(params["sample.fraction"]), 
                                num.trees = 1000,
                                min.node.size = as.numeric(params["min.node.size"]),  
                                mtry =  as.numeric(params["mtry"]), # default is ceiling(sqrt(num.col.covariates) + 20)
                                alpha = as.numeric(params["alpha"]),
                                imbalance.penalty = as.numeric(params["imbalance.penalty"]),
                                seed = 0,
                                num.threads = 10)
    return(tau.forest)
}

run.hte <- function(mat, mutations, dataset, cancerType, is.binary = TRUE, thres = 0.75){
    # @mat: covariates and treatment assigments, including gender, sex, and mutations etc. For treatment assignment,
    #   should be 0 or 1, can be derived by treating sigle or multiple gene mutation as treatment assigment
    # @mutations: a vector of mutations will be investigated.
    # @dataset: contains all information, including patients and mutations
    
    subsampPValue <- data.frame(mutName = character(),
                                permutated.pval = double(),
                                # permutated.covar.pval = double(),
                                subsampling.pval = double(),
                                subsampling.var.pval = double(),
                                stringsAsFactors = FALSE)
    
    Y <- dataset$imputed.log.times

    for(mut in mutations){
        # since treatment variable is 0 or 1, if the feature considered is binary, then no transformation is neeeded;
        # otherwise, we set values of the feature greater than specific quantile, say 0.75, to 1
        if(is.binary == TRUE){
            treatment <- as.numeric(!mat[, mut] == 0)  # for binary number treatment
            # treatment <- mat[, mut] # for real number treatment assignment
        }else{
            treatment <- as.numeric(mat[, mut] >= quantile(mat[, mut], thres))
        }

        X.covariates <- select(mat, -mut)

	# method <- match.arg(method, c("rfscr", "cf", "surv_bart"))
        params = tune_causal_forest(X.covariates, Y, treatment)$params
	
        tau.forest <- cf(X.covariates, Y, treatment, params)   # run causal forests by default
        tau.prediction <- predict(tau.forest, newdata = NULL, estimate.variance = TRUE, num.threads = 8)

        # compute zval, pval and ajusted.p
        tau.zval <- tau.prediction$predictions/sqrt(tau.prediction$variance.estimates)
        tau.pval <- pnorm(abs(tau.zval), lower.tail = FALSE) * 2
        tau.p.adjust <- p.adjust(tau.pval, method = 'BH')

        overall.pval <- simes.test(tau.pval)
        
        if(overall.pval <= 0.05){ 
            # save the result
            pred.ret <- data.frame(tau.prediction$predictions, tau.zval, tau.pval, tau.p.adjust)
            colnames(pred.ret) <- c('tau.val', 'tau.zval', 'tau.pval', 'tau.p.adjust')
            write.table(pred.ret, paste0(cancerType, '.tau.', mut, '.csv'), sep = ',', quote = F, row.names = F)

            # *****************************************************************************************
            # permutate estimated tau values to validate HET esimation
            # *****************************************************************************************
            # cacluate tau-risk_R
            no.obs <- length(Y)
            Y.hat <- tau.forest[["Y.hat"]]
            W.hat <- tau.forest[["W.hat"]]
            tau.1 <- Y - Y.hat
            tau.2 <- (treatment - W.hat) * tau.prediction$predictions

            # compute permutated p.value for the mutation
            permutated.p.val <- permutated.pval(tau.1, tau.2)
            # permutated.p.val2 <- permutated.pval(Y, Y.hat + tau.2) # discarded
            print(paste0('permutated.p.val of tau.values for ', mut, ':', permutated.p.val))

            # *****************************************************************************************
            # permutate covariates to validate HET esimation
            # *****************************************************************************************
            tau.risk <- assess.tau.risk(tau.forest, Y, treatment)
            # perm.pval <- permutate.covariates.testing(X.covariates, Y, treatment, tau.risk, num.strap = 1000)
            # print(paste0('permutated.p.val of covariates for ', mut, ':', perm.pval))

            # *****************************************************************************************
            # use subsampling to validate HET esimation
            # *****************************************************************************************
            avg.treatment.hete <- mean(Y[treatment == 1])
            avg.control.hete <- mean(Y[treatment == 0])
            # avg.hete <- (avg.control.hete + avg.treatment.hete)/length(treatment)
            tau.hat <- avg.treatment.hete - avg.control.hete
            tau.var <- var(tau.prediction$predictions)

            constant.risk <- assess.constant.risk(Y, treatment) 
            
            pvals <- bootstrap.testing(X.covariates, Y, treatment, tau.risk, constant.risk, tau.var, num.strap = 1000)
            print(paste0('subsampling.pval for ', mut, ':', pvals[1]))
            print(paste0('var.pval for ', mut, ':', pvals[2]))
            
            subsampPValue[nrow(subsampPValue) + 1, ] <- list(mut, 
                                                             permutated.p.val, 
                                                             # perm.pval,
                                                             pvals[1],
                                                             pvals[2])
            
            # extract feature importance and save
            varImp <- variable_importance(tau.forest,  max.depth = 4)
            varImp.ret <- data.frame(variable = colnames(X.covariates),  varImp)
            write.table(varImp.ret, paste0(cancerType, '.varimp.', mut, '.csv'), sep = ',', quote = F, row.names = F)
        }
    }
    return(subsampPValue)
}

cancerType <- 'OV' # cancer type to be studied
mutation.threshold <- 30  # minimum number of mutation to be studied
setwd("/exeh_3/kai/data/heterogenous_treatment_effects")

# load TCGA dataset
dataset <- read.csv('TCGA_survival_CDR.csv', header = TRUE, stringsAsFactors = F)
dataset <- dataset[c(2, 5, 4, 32, 33, 3)]   # select relevant information
dataset$gender <- as.numeric(dataset$gender == 'FEMALE')

# some clearning
dataset <- dataset[dataset$PFI.time!=0,]
dataset[dataset == "#N/A"] <- NA
dataset <- subset(dataset, !is.na(PFI.time) & !is.na(PFI) & !is.na(age_at_initial_pathologic_diagnosis))

# load ICGC dataset
icgcMutation <- read.csv('ICGC/ICGC.covatiates.mutation.csv', header = T, stringsAsFactors = F)
donorInfo <- icgcMutation[c(1:4, 8, 5)]
donorInfo$sex <- as.numeric(donorInfo$sex == 'female')

# change names of two dataset
index <- c('donorId', 'gender', 'age', 'censor', 'PFI.time', 'type')
colnames(dataset) <- index
colnames(donorInfo) <- index

# union two datasets
patientInfo <-rbind(dataset, donorInfo)

# specify cancer type here
patientInfo <- subset(patientInfo, type == cancerType)
surv.times <- as.numeric(as.character(patientInfo$PFI.time))
cens <- as.numeric(patientInfo$censor)

# get imputed log survival times
max.censored <- max(surv.times[cens == 0])
cens[surv.times == max.censored] <- 1
imputed.log.times = impute.survival(surv.times, cens)

# attach imputed.log.times to original dataset
dataset <- cbind(patientInfo, imputed.log.times)

#**********************************************************************************
#                                mutation data
#**********************************************************************************

# load and process TCGA dataset
load("TCGA_PanCanAtlas_No_Of_Mutations_NoSilentMut_intersected.rda")
donorId <- substr(intersected.mut_adj.rev$Tumor_Sample_Barcode, start=1, stop=12)
mut_adj.rev <- cbind(donorId, intersected.mut_adj.rev[2: dim(intersected.mut_adj.rev)[2]])

# obtain ICGC mutation data
icgcMut <- icgcMutation[c(1, 9: dim(icgcMutation)[2])]
colnames(icgcMut)[1] <- 'donorId'

# rbind two dataset
mutDataset <- rbind(mut_adj.rev, icgcMut)

# merge dataset with mutitation dataset by "bcr_patient_barcode"
merged.dataset <- inner_join(dataset, mutDataset, by="donorId")

# find mutations with mutation count greater than the defined mutation.threhold
pos <- which(colnames(merged.dataset) == "TSPAN6")
freq.mutations <- extract.freq.mutation(merged.dataset, pos, mutation.threshold)

mutation.mat <- merged.dataset[, pos:ncol(merged.dataset)]

# adding confonding variables
age <- merged.dataset$age
gender <- merged.dataset$gender

# remove stage variable first because of lots of missing value
# stage <- as.factor(merged.dataset$clinical_stage)
mutation.mat <- cbind(age, gender, mutation.mat)
mutation.mat$age <- as.numeric(as.character(mutation.mat$age))

# generate train id for causal forests
set.seed(1)
obsNumber <- dim(mutation.mat)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)

# run causal forests model
result <- run.hte(mutation.mat, freq.mutations, merged.dataset, cancerType)
write.csv(result, paste0(cancerType, '.mutation.subsampling.result.csv'), sep = ',', quote = F, row.names = F)


#*******************************************************************************
#                            CNV data
#*******************************************************************************
# CNV.dataset <- fread("ISAR_GISTIC.all_thresholded.by_genes_formatted.csv")
# merged.dataset <- inner_join(dataset, CNV.dataset, by="bcr_patient_barcode")
#
# # find mutations with mutation count greater than the defined mutation.threhold
# pos <- which(colnames(merged.dataset) == "X1.2.SBSRNA4")
# CNV.freq.mutations <- extract.freq.mutation(merged.dataset, pos, mutation.threshold, is.binary = FALSE)
#
# CNV.mat <- merged.dataset[, pos:ncol(merged.dataset)]
#
# # adding confonding variables
# age <- as.numeric(as.character(merged.dataset$age_at_initial_pathologic_diagnosis))
# gender <- as.numeric(merged.dataset$gender == "FEMALE")
# because of many missing values in stage, so ignore stage currently
# stage <- as.factor(merged.dataset$clinical_stage)
# CNV.mat <- cbind(age, gender, CNV.mat)
# run.hte(CNV.mat, CNV.freq.mutations, merged.dataset)

# #*******************************************************************************
#                             expression (RNAseq) data
# #*******************************************************************************
# expr.dataset = fread("PANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp_formatted.csv")
# merged.dataset = inner_join(dataset, expr.dataset, by="bcr_patient_barcode")
#
# pos = which(colnames(merged.dataset) == "?|100130426")
# expr.freq.mutations <- extract.freq.mutation(merged.dataset, pos, mutation.threshold, is.binary = FALSE)
#
# expr.mat = merged.dataset[, pos:ncol(merged.dataset)]
#
# # adding confonding variables
# age <- as.numeric(as.character(merged.dataset$age_at_initial_pathologic_diagnosis))
# gender <- as.numeric(merged.dataset$gender == "FEMALE")
# # because of lots of missing values, ignore stage at first
# stage <- as.factor(merged.dataset$clinical_stage)
# exp.mat <- data.frame(age, gender, expr.mat)
# run.hte(expr.mat, expr.freq.mutations, merged.dataset, is.binary = FALSE)
