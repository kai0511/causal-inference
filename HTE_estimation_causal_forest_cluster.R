library(grf)
library(MASS)
library(dplyr)
library(foreach)
# library(doMC)
# library(causalLearning)
# library(devtools)
# library(TCGA2STAT)
library(SummarizedExperiment)
# library(TCGAbiolinks)
# install_github("saberpowers/causalLearning")
library(data.table)
library(survminer)
library(mppa)
library(doParallel)
library(methods)

source('/exeh_3/kai/R.code/HTE/survival.imputation.R')
source('/exeh_3/kai/R.code/HTE/causal_inference_models.R')
source('/exeh_3/kai/R.code/HTE/hte.validation.R')

registerDoParallel(20)

extract.freq.mutation <- function(dataset, pos, is.binary = TRUE, threshold = 10){
    # @dataset mutation dataset
    # @pos the position where mutation variables start
    # @threshold threshold for frequent mutation
    mutation.names <- colnames(dataset)[pos:ncol(dataset)]

    mutation.selected <- sapply(mutation.names, function(x) {
        temp <- dataset[c('type', x)]
        colnames(temp) <- c('type', 'mut')
        
        if(is.binary == TRUE){
            temp$is_mut <- temp$mut != 0
        } else {
            temp$mut >= 0.75
        }
        res <- temp %>% group_by(type) %>% summarise(cnt = sum(is_mut))
        return(all(res$cnt >= threshold))
    })

    # find mutations with mutation count greater than the defined threshold
    freq.mutations <- mutation.names[mutation.selected]
    print(paste0('number of mutation will be studied:', sum(mutation.selected)))
    return(freq.mutations)
}

# if the generic function existis, i.e., append, we should use the same argument names as exist.
# otherwise, we need to create a new one.
# setGeneric("append") # this declear may be omited since append function already exists.
setMethod("append", signature(x = "data.frame", values = "vector"),
    function(x, values) {
        if(ncol(x) != length(values)){
            stop('dimension of data.frame does not match with the length of values.')
        }else if(typeof(x) == 'list'){
            x[nrow(x)+1, ] <- values
            x
        }else{
            x[nrow(x)+1, ] <- as.list(values)
            x
        }
    }
)

compute_stats <- function(tau.pred){
    # compute zval, pval and ajusted.p
    tau.zval <- tau.pred$predictions/sqrt(tau.pred$variance.estimates)
    tau.pval <- pnorm(abs(tau.zval), lower.tail = FALSE) * 2
    tau.p.adjust <- p.adjust(tau.pval, method = 'BH')
    stats <- data.frame(tau.pred$predictions, tau.zval, tau.pval, tau.p.adjust)
    return(stats)
}

correlation_test <- function(x, y, methods, alt = 'greater') {
    print(head(data.frame(x, y)))
    fitted_res <- sapply(methods, function(m) {
        fitted <- cor.test(x, y, method = m, alternative = alt, use = "na.or.complete")
        print(fitted)
        c(fitted$estimate, fitted$p.value)
    })
    print(fitted_res)
    fitted_res
}

stratified_sample <- function(cluster_id){
    obs <- seq(length(cluster_id))
    splits <- rep(FALSE, length(cluster_id))
    clusters <- levels(cluster_id)

    for(cls in clusters){
        idx <- obs[clusters == cls]
        train_id <- sample(idx, floor(length(idx)/2), replace = FALSE)
        splits[train_id] = TRUE
    }
    return(data.frame(obs, splits))
}

split_half_testing <- function(covariates, Y, treatment, cluster_id, no_repeats = 10, n_core = 8){

    obsNumber <- dim(covariates)[1]
    correlation_matrix <- foreach(i = seq(no_repeats), .combine = 'rbind') %dopar%  {
        set.seed(i)

        splitset <- stratified_sample(cluster_id)
        trainId <- splitset$obs[splitset$splits == TRUE]
        train_cluster_id <- cluster_id[splitset$splits == TRUE]
        # trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)

        no.obs <- dim(covariates)[1]
        no.obs.train <- length(trainId)
        no.obs.test <- no.obs - no.obs.train

        tau.forest.train <- cf.tuned(covariates[trainId, ], Y[trainId], treatment[trainId], cluster_id = cluster_id[splitset$splits == TRUE])   
        tau.pred.train.a <- predict(tau.forest.train, newdata = NULL, estimate.variance = TRUE, num.threads = n_core)
        tau.pred.test.a <- predict(tau.forest.train, newdata = covariates[-trainId, ], estimate.variance = TRUE, num.threads = n_core)

        # params = tune_causal_forest(X.covariates[-trainId, ], Y[-trainId], treatment[-trainId])
        tau.forest.test <- cf.tuned(covariates[-trainId, ], Y[-trainId], treatment[-trainId], cluster_id = cluster_id[splitset$splits == FALSE])   
        tau.pred.train.b <- predict(tau.forest.test, newdata = NULL, estimate.variance = TRUE, num.threads = n_core)
        tau.pred.test.b <- predict(tau.forest.test, newdata = covariates[trainId, ], estimate.variance = TRUE, num.threads = n_core)

        # compute z-score, pvalues, and ajusted.pvalues
        tau.a.train.stats <- compute_stats(tau.pred.train.a)
        tau.a.stats <- compute_stats(tau.pred.test.a)
        tau.b.train.stats <- compute_stats(tau.pred.train.b)
        tau.b.stats <- compute_stats(tau.pred.test.b)

        simes_pval.a <- simes.test(tau.a.stats[[1]])
        simes_pval.b <- simes.test(tau.b.stats[[1]])

        partial_simes_pval.a <- simes.partial(floor(no.obs.test * 0.05), tau.a.stats[[1]])
        partial_simes_pval.b <- simes.partial(floor(no.obs.train * 0.05), tau.b.stats[[1]])

        # check the correlation between two predictions from two datasets
        test_res.a <- correlation_test(tau.a.train.stats[[1]], tau.b.stats[[1]], methods = c('pearson', 'kendall', 'spearman'))
        test_res.b <- correlation_test(tau.a.stats[[1]], tau.b.train.stats[[1]], methods = c('pearson', 'kendall', 'spearman'))

        fisher.pval.a <- fisher.exact.test(tau.a.train.stats[[3]], tau.b.stats[[3]])
        fisher.pval.b <- fisher.exact.test(tau.a.stats[[3]], tau.b.train.stats[[3]])

        t.test.pval.a <- quantile.t.test(tau.a.train.stats[[1]], tau.b.stats[[1]])
        t.test.pval.b <- quantile.t.test(tau.a.stats[[1]], tau.b.train.stats[[1]]) 
        c(simes_pval.a, simes_pval.b, partial_simes_pval.a, partial_simes_pval.b, test_res.a, test_res.b, fisher.pval.a, fisher.pval.b, t.test.pval.a, t.test.pval.b) 
    }

    print(correlation_matrix)

    aggregated_result <- sapply(seq(dim(correlation_matrix)[2]), function(i){
        if(i %in% c(5, 7, 9, 11, 13, 15)){
            return(median(correlation_matrix[, i], na.rm = TRUE))
        }else{
            removed_na <- na.omit(correlation_matrix[, i])
            simes_pval <- ifelse(length(removed_na) > 0, simes.test(removed_na), NA)
            return(simes_pval)
        }
    })
    return(aggregated_result)
}

run.hte <- function(mat, mutations, dataset, cancer_type, trainId, is.binary = TRUE, thres = 0.75, n_core = 8){
    # @mat: covariates and treatment assigments, including gender, sex, and mutations etc. For treatment assignment,
    #   should be 0 or 1, can be derived by treating sigle or multiple gene mutation as treatment assigment
    # @mutations: a vector of mutations will be investigated.
    # @dataset: contains all information, including patients and mutations
    
    correlation.test.ret <- data.frame(mutName = character(),
                                       simes.pval.a = double(),
                                       simes.pval.b = double(),
                                       partial.simes.pval.a = double(), 
                                       partial.simes.pval.b = double(),
                                       train.pearson.estimate = double(),
                                       train.pearson.pvalue = double(),
                                       test.pearson.estimate = double(),
                                       test.pearson.pvalue = double(),
                                       train.kendall.estimate = double(),
                                       train.kendall.pvalue = double(), 
                                       test.kendall.estimate = double(),
                                       test.kendall.pvalue = double(), 
                                       train.spearman.estimate = double(),
                                       train.spearman.pvalue = double(),
                                       test.spearman.estimate = double(),
                                       test.spearman.pvalue = double(),
                                       stringsAsFactors = FALSE)

    double.dataset.test.ret <- data.frame(mutName = character(),
                                          train.fisher.pval = double(),
                                          test.fisher.pval = double(),
                                          train.t.test.pval.a = double(),
                                          train.t.test.pval.b = double(),
                                          test.t.test.pval.a = double(),
                                          test.t.test.pval.b = double(),
                                          stringsAsFactors = FALSE)

    calibration.ret <- data.frame(mutName = character(),
                                  simes.pval = double(),
                                  partial.conjunction.pval = double(),
                                  cor.mut.Y.estimate = double(),
                                  cor.mut.Y.pvalue = double(),
                                  permutation.pval = double(),
                                  mean.pred.estimate = double(),
                                  mean.pred.pval = double(),
                                  differential.pred.estimate = double(),
                                  differential.pred.pval = double(),
                                  stringsAsFactors = FALSE)

    permutate.testing.ret <- data.frame(mutName = character(),                                       
                                        var.pval = double(),
                                        risk.pval = double(),
                                        fixed.W.risk.pval = double(),
                                        fixed.YW.risk.pval = double(),
                                        stringsAsFactors = FALSE)
    no.obs <- dim(mat)[1]
    no.obs.train <- length(trainId)
    no.obs.test <- no.obs - no.obs.train

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

        X.covariates <- as.matrix(select(mat, -mut))
        
        print(paste0(c('#', rep('-', 40), ' begin a new mutation ', rep('-', 40)), collapse = ''))
        
        # split whole dataset into two parts, and the idea of validation is similar to prediction strength.
        aggregated_ret <- split_half_testing(X.covariates, Y, treatment, cancer_type)
        print(aggregated_ret)

        # print some result from above 
        print(paste0('mutation name:', mut))
        print(paste0('Fisher extact test pval in trainset:', aggregated_ret[17], '; corresponding value in testset:', aggregated_ret[18]))
        print(paste0('pearson correlation pval of tau predictions in trainset:', aggregated_ret[6], '; pval in testset:', aggregated_ret[12]))
        print(paste0('spearman correlation pval of tau predictions in trainset:', aggregated_ret[10], '; pval value in testset:', aggregated_ret[16]))
        
        # append results to correspoding dataframes
        current_ret <- do.call(c, list(list(mut), as.list(aggregated_ret[1: 16])))
        correlation.test.ret <- append(correlation.test.ret, current_ret)
        two_sample_test_ret <- do.call(c, list(list(mut), aggregated_ret[17: 22]))
        double.dataset.test.ret <- append(double.dataset.test.ret, two_sample_test_ret)
                                                                                              
        # *****************************************************************************************
        # fit causal forest on the whole dataset
        # validate fitting with permutating covariates
        # *****************************************************************************************
        tau.forest <- cf.tuned(X.covariates, Y, treatment, cluster_id = cancer_type)   # run causal forests by default
        tau.prediction <- predict(tau.forest, newdata = NULL, estimate.variance = TRUE, num.threads = n_core)

        # compute zval, pval and ajusted.p
        tau_stats <- compute_stats(tau.prediction)
        simes.pval <- simes.test(tau_stats$tau.pval)
        partial.simes.pval <- simes.partial(floor(no.obs * 0.05), tau_stats[[3]])

        if(simes.pval <= 0.05){ 
            cor.overall <- cor.test(mat[, mut], Y, method = 'pearson', alternative = 'greater', use="na.or.complete")
            # save the result
            pred.ret <- cbind(dataset$donorId, tau_stats) 
            colnames(pred.ret) <- c('donorId', 'tau.val', 'tau.zval', 'tau.pval', 'tau.p.adjust')
            write.csv(pred.ret, paste0('mutation_result/cluster/cluster.tau.', mut, '.csv'), quote = F, row.names = F)

            # permutate estimated tau values to validate HET esimation
            Y.hat <- tau.forest[["Y.hat"]]
            W.hat <- tau.forest[["W.hat"]]
            tau.1 <- Y - Y.hat
            tau.2 <- (treatment - W.hat) * tau.prediction$predictions

            # compute permutated p.value for the mutation
            permutated.p.val <- permutated.pval(tau.1, tau.2)
            # permutated.p.val2 <- permutated.pval(Y, Y.hat + tau.2) # discarded
            print(paste0('p.val by permutating (Y - Y.hat) or [(W - W.hat) * tau] for ', mut, ':', permutated.p.val))

            # test of the calibration with test_calibration from grf 
            # as reported, the pvalue from the method is not acuuray. 
            calibration.fit <- test_calibration(tau.forest)
            print(paste0('mean.pred.estimate of test_calibration:', calibration.fit[1, 1], '; its pval:', calibration.fit[1, 4]))
            print(paste0('differential.pred.estimate of test_calibration:', calibration.fit[2, 1], '; its pval:', calibration.fit[2, 4]))

            # save results
            test_result <- c(calibration.fit[1, c(1, 4)], calibration.fit[2, c(1, 4)])
            test_result <- c(simes.pval, partial.simes.pval, cor.overall$estimate, cor.overall$p.value, permutated.p.val, test_result)
            record <- do.call(c, list(list(mut), as.list(test_result)))
            calibration.ret <- append(calibration.ret, record)

            # valiate hte using by shuffling covariates, and we validate whether the variance of tau after shuffling is likey to 
            # to be less than that observed. 
            # And also, we also use the tau_risk defined by other authors to validate, but in this case, we expect that the tau risk
            # observed is more likey to be smaller tau risk from permutation.
            tau.risk <- assess.explained.tau.risk(tau.forest, Y, treatment)
            fixed.W.tau.risk <- assess.explained.tau.fixed.W.risk(tau.forest, Y, treatment, W.hat)
            fixed.YW.tau.risk <- assess.explained.tau.fixed.YW.risk(tau.forest, Y, Y.hat, treatment, W.hat)


            perm.pvals <- permutate.covariates.testing(X.covariates, 
                                                        Y, Y.hat, 
                                                        treatment, W.hat, 
                                                        tau.risk, 
                                                        fixed.W.tau.risk, 
                                                        fixed.YW.tau.risk, 
                                                        var(tau.prediction$predictions), 
                                                        mut = mut, 
                                                        num.strap = 500,
                                                        num_trees = 1000)

            perm_pval_record <- do.call(c, list(list(mut), as.list(perm.pvals)))
            permutate.testing.ret <- append(permutate.testing.ret, perm_pval_record)    

            print(paste0('pval of variance by permutation covariates:', perm.pvals[1]))
            print(paste0('pval of tau.risk (without fixed YW) by permutation covariates:', perm.pvals[2]))
            print(paste0('pval of tau.risk (fixed W) by permutation covariates:', perm.pvals[3]))
            print(paste0('pval of tau.risk (fixed YW) by permutation covariates:', perm.pvals[4]))

            # extract feature importance and save
            varImp <- variable_importance(tau.forest,  max.depth = 4)
            varImp.ret <- data.frame(variable = colnames(X.covariates),  varImp)
            write.csv(varImp.ret, paste0('mutation_result/cluster/cluster.varimp.', mut, '.csv'), quote = F, row.names = F) 
        }
    }
    return(list(correlation.test.ret, calibration.ret, double.dataset.test.ret, permutate.testing.ret))
}

cancerType <- c('BRCA', 'GBM', 'OV', 'UCEC', 'KIRC', 'HNSC', 'LGG', 'THCA', 'PRAD', 'LUAD', 'LUSC') # cancer types to be included
mutation.threshold <- 50  # minimum number of mutation to be included
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
patientInfo <- subset(patientInfo, !is.na(age))

# specify cancer type here
patientInfo <- subset(patientInfo, type %in% cancerType)
surv.times <- as.numeric(as.character(patientInfo$PFI.time))
cens <- as.numeric(patientInfo$censor)
cluster_id <- as.factor(patientInfo$type)

# get imputed log survival times
# max.censored <- max(surv.times[cens == 0])
# cens[surv.times == max.censored] <- 1 
imputed.log.times = impute.survival(surv.times, cens, cluters = cluster_id)

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
# freq.mutations <- extract.freq.mutation(merged.dataset, pos, threshold = mutation.threshold)
freq.mutations <- c('GABRA3', 'CNTLN', 'PREX2', 'ARHGAP6', 'CTNNA2', 'MYO3B', 'MYO3B', 'TRPS1', 'ARHGAP15', 'GLIS3', 'MTUS2', 'CRB1')

mutation.mat <- merged.dataset[, pos:ncol(merged.dataset)]

# adding confonding variables
age <- merged.dataset$age
gender <- merged.dataset$gender

# remove stage variable first because of lots of missing value
# stage <- as.factor(merged.dataset$clinical_stage)
mutation.mat <- cbind(age, gender, mutation.mat)
mutation.mat$age <- as.numeric(as.character(mutation.mat$age))
cluster_id <- as.factor(merged.dataset$type)

# generate train id for causal forests
set.seed(1)
obsNumber <- dim(mutation.mat)[1]
trainId <- sample(1: obsNumber, floor(obsNumber/2), replace = FALSE)

# split TCGA and TCGC separately
# trainId <- grep('TCGA', dataset$donorId)

# run causal forests model
result <- run.hte(mutation.mat, freq.mutations, merged.dataset, cluster_id, trainId)
write.csv(result[[1]], paste0('mutation_result/cluster/cluster.mutation.correlation.test.result.csv'), quote = F, row.names = F)
write.csv(result[[2]], paste0('mutation_result/cluster/cluster.mutation.calibration.result.csv'), quote = F, row.names = F)
write.csv(result[[3]], paste0('mutation_result/cluster/cluster.mutation.median.t.test.result.csv'), quote = F, row.names = F)
write.csv(result[[4]], paste0('mutation_result/cluster/cluster.mutation.permutate.testing.result.csv'), quote = F, row.names = F)

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



