library(grf)
library(MASS)
library(dplyr)
library(data.table)
library(survival)
library(survminer)
library(doParallel)

registerDoParallel(1)
args = commandArgs(trailingOnly=TRUE)

source('survival.imputation.R')  # load impute.survival function
source('/exeh_3/kai/R.code/HTE/hte.validation.R')     # load validation methods for HTE
source('simulation.scenarios.R') # load different simulation functions
source('causal_inference_models.R')

generate.covariates <- function(n, p = 10){
    X = matrix(rnorm(n*p), nrow = n, ncol = p)
    X[, seq(2, p, by = 2)] = (X[, seq(2, p, by = 2)] > 0)
    colnames(X) <- paste0('v', seq(p))
    return(X)
}

randomized.tx <- function(n){
    rbinom(n, size = 1, p = 0.5) # random treatment assignment
}

biased.tx <- function(mu.x, tau.x){
    diff <- mu.x - tau.x/2
    probs <- exp(diff) / (1 + exp(diff))
    # tx <- sapply(probs, function(i) rbinom(1, size = 1, p = i))
    tx <- rbinom(length(mu.x), 1, probs)
    tx
}

generate.censored.outcome <- function(Y, event.rate){
    cutoff <- quantile(Y, 1 - event.rate)
    cen <- rexp(length(Y), rate = 1/cutoff)
    y.censored <- pmin(Y, cen)
    censor <- as.numeric(Y <= y.censored)
    res <- cbind(y.censored, censor)
    colnames(res) <- c('y.censored', 'censor')
    return(res)
}

generate.scenario <- function(n, p, mu.fun, tau.fun, noise.sd = 1, randomized = TRUE, event.rate = 0.1){
    x <- generate.covariates(n, p)

    mu.x <- do.call(mu.fun, list(x))
    tau.x <- do.call(tau.fun, list(x))
    # mu.x <- f8(x)
    # tau.x <- f1(x)

    if(randomized == TRUE){
        tx <-randomized.tx(n)
    }else{
        tx <- biased.tx(mu.x, tau.x)
    }

    log.y <- mu.x + (tx - 1/2) * tau.x + rnorm(n, sd = noise.sd)
    # log.y <- mu.x + tx  * tau.x + rnorm(n, sd = noise.sd)

    y <- exp(log.y)
    # y <- rexp(n, rate = 1/y.bar)
    res <- generate.censored.outcome(y, event.rate)

    max.censored <- max(res[, 1][res[, 2] == 0])
    res[, 2][res[, 1] == max.censored] <- 1

    imputed.y <- impute.survival(res[, 1], res[, 2])

    d <- cbind(y, tau.x, tx, imputed.y, res[, 2], x)
    # d <- cbind(y, tau.x, tx, log(y), rep(0, length(y)), x)
    colnames(d)[1:5] <- c('y',  'tau', 'tx', 'y.imputed', 'censor')
    return(d)
}

# original sample size: 600
d1 <- generate.scenario(450, 400, 'f8', 'f1')
d9 <- generate.scenario(450, 400, 'f8', 'f1', randomized = FALSE)

d2 <- generate.scenario(450, 400, 'f5', 'f2', noise.sd = 0.25)
d10 <- generate.scenario(450, 400, 'f5', 'f2', noise.sd = 0.25, randomized = FALSE)

# original sample size: 700
d3 <- generate.scenario(550, 300, 'f4', 'f3')
d11 <- generate.scenario(550, 300, 'f4', 'f3', randomized = FALSE)

d4 <- generate.scenario(550, 300, 'f7', 'f4', noise.sd = 0.25)
d12 <- generate.scenario(550, 300, 'f7', 'f4', noise.sd = 0.25, randomized = FALSE)

# original sample size: 800
d5 <- generate.scenario(600, 200, 'f3', 'f5')
d13 <- generate.scenario(600, 200, 'f3', 'f5', randomized = FALSE)

d6 <- generate.scenario(600, 200, 'f1', 'f6')
d14 <- generate.scenario(600, 200, 'f1', 'f6', randomized = FALSE)

# original sample size: 900
d7 <- generate.scenario(700, 100, 'f2', 'f7', noise.sd = 4)
d15 <- generate.scenario(700, 100, 'f2', 'f7', noise.sd = 4, randomized = FALSE)

d8 <- generate.scenario(700, 100, 'f6', 'f8', noise.sd = 4)
d16 <- generate.scenario(700, 100, 'f6', 'f8', noise.sd = 4, randomized = FALSE)

# num <- c(seq(9, 1, -1), seq(16, 10, -1))
# no_repeats <- 200
# num <- c(rep(1, no_repeats), rep(9, no_repeats))

print_false_rate <- function(test_names, false_rates){
    num <- length(test_names)

    for(i in seq(num)){
        print(paste0(test_names[i], ': ', false_rates[i]))
    }
}

# --------------------------------------------------------------------------------------
# run simulation for each scenario for a number of times, defined by no_repeats below
# collect simulation result wiht matrix and save it in local directory
# --------------------------------------------------------------------------------------

if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
    cat('Scenario: ', args[1], '\n')
}

num <- as.numeric(args[1])
no_repeats <- 500 
verbose <- TRUE; is_tuned <- FALSE; is_save <- TRUE; seed <- NULL

col_names <- c('simes.pval', 'partial.simes.pval', 'pearson.estimate','pearson.pvalue',
               'kendall.estimate','kendall.pvalue', 'spearman.estimate','spearman.pvalue',
               'fisher.pval', 't.test.a.pval', 't.test.b.pval') 
calibration_test_names <- c('differential.estimate', 'differential.pvalue')
perm_test_names <- c('variance.pval', 'fixed.YW.pval')

cf.estimator <- ifelse(is_tuned, cf.tuned, cf)

for(i in seq(length(num))){
    
    data <- get(paste0('d', num[i]))
    
    write.csv(data, file = paste0('scenario_', num[i], '_dataset.csv'), quote = FALSE, row.names = F)
   
    # remove rows containing NAs 
    ind <- apply(data, 1, function(x) sum(is.na(x)) == 0)
    data <- data[ind, ]

    covariates <- data[, -c(1:5)]
    real.Y <- data[, 1]
    tau.real <- data[, 2]
    W <- data[, 3]
    imputed.Y <- data[, 4]
    censor <- data[, 5]
    sample_size <- dim(covariates)

    scenario_result <- foreach(j = seq(no_repeats), .combine = 'rbind', .options.multicore = list(set.seed = T)) %do%  {
        # set.seed(j)
        file_prefix <- paste0('simulation_result/scenario_', num[i], '_repeat_', j)
        
        fitted.obj <- cf.estimator(covariates, imputed.Y, W, num_trees = 1000)
        tau.pred <- predict(fitted.obj, newdata = NULL, estimate.variance = TRUE, num.threads = 8)
        
        # print(fitted.obj[["tuning.output"]])
        
        Y.hat <- fitted.obj[["Y.hat"]]
        W.hat <- fitted.obj[["W.hat"]]

        cf.risk = mean((tau.real - tau.pred$predictions)^2)
        tau.var <- var(tau.pred$predictions) 
        permutated.p.val <- permutated.pval(tau.real, (W - W.hat) * tau.pred$predictions) 
        
        Y.hat <- fitted.obj[["Y.hat"]]
        tau.estimated <- imputed.Y - Y.hat
        
        cor.obj <- cor.test(tau.real, tau.pred$predictions, alternative = 'greater')
        obs_calibration_pvals <- unname(test_calibration(fitted.obj)[2, c(1, 4)])
        risk <- mean((tau.real - tau.pred$predictions)^2)

        if(verbose){
            print(paste0('# ------------------------------ begin scenario ', num[i],' ---------------------------------------'))
            print(paste0('censor rate:', 1 - mean(censor)))
            print(paste0('risk fitted by CF in scenario ', num[i], ': ', cf.risk))

            # print(paste0('variance in scenario ', i, ': ', var(tau.real)))    
            # print(paste0('R square in scenario ', i, ': ', 1 - var(tau.real - tau.pred$predictions) / var(tau.real)))
            print(paste0('p.val of permutating either real tau or (W - W.hat) * tau.pred$predictions in scenario ', 
                num[i], ':', permutated.p.val))
            print(paste0('The estimated measure of correlation between real tau and predicted tau in scenario ', num[i], ': ', 
                cor.obj$estimate, '; its pval:', cor.obj$p.value))
            print(paste0('risk of mean squared difference between real tau and predicted tau value in scenario ', num[i], ':', risk))
        }

        # *****************************************************************************************
        # permutate covariates to validate HET esimation
        # *****************************************************************************************
        no.obs <- length(imputed.Y)
        tau.1 <- imputed.Y - Y.hat
        tau.2 <- (W - W.hat) * tau.pred$predictions
        permutated.p.val <- permutated.pval(tau.1, tau.2)

        if(verbose){
            print(paste0('p.val of permutating either (Y - Y.hat) or [(W - W.hat) * tau.pred$predictions] in scenario ', num[i], ':', permutated.p.val))
            print(test_calibration(fitted.obj))
        }
       
        calibration_pvalues <- unname(test_calibration(fitted.obj)[2, c(1, 4)]) 

        # print('# ------------------------ result of calibration by maunual regression -----------------------------------')
        # tau.bar <- mean(tau.pred$predictions)
        # pred.ret <- data.frame(outcome = tau.1, tau.diff = tau.pred$predictions - tau.bar, tau.avg = rep(tau.bar, length(tau.1)))
        # print(summary(lm(outcome ~ tau.diff + tau.avg + 0, data = pred.ret)))
         
         
        print('# ------------------------------ result of permutation covariates ---------------------------------------')
        fixed.YW.tau.risk <- assess.explained.tau.fixed.YW.risk(fitted.obj, imputed.Y, Y.hat, W, W.hat)    
         
        # perm.pvals <- permutate.covariates.testing(covariates, 
        #                                            imputed.Y, 
        #                                            Y.hat, 
        #                                            W, W.hat, 
        #                                            fixed.YW.tau.risk, 
        #                                            tau.var, 
        #                                            is_tuned = is_tuned,
        #                                            is_save = is_save,
        #                                            file_prefix = file_prefix,
        #                                            num_trees = 1000, 
        #                                            num.strap = 500) 
        
        perm.pvals <- adaptive.permutate.covariates.testing(covariates, 
                                                            imputed.Y, 
                                                            Y.hat, 
                                                            W, W.hat, 
                                                            fixed.YW.tau.risk, 
                                                            tau.var,
                                                            is_tuned = is_tuned,
                                                            is_save = is_save,
                                                            file_prefix = file_prefix,
                                                            num_trees = 1000, 
                                                            num.strap = 500) 
        if(verbose){
            print(paste0('p-values of variance testing by permutating covariates: ', perm.pvals[1]))
            print(paste0('p-values of tau risk by permutating covariates with fixing YW: ', perm.pvals[2]))
        }


        # *****************************************************************************************
        # replication testing by spliting dataset into two pieces.
        # *****************************************************************************************
        no.obs <- dim(covariates)[1]
        # seed <- (j + floor(runif(1) * 888))

        pvalues <- split_half_testing(covariates, 
                                      imputed.Y, W, 
                                      binary = T, 
                                      is_save = T, 
                                      is_tuned = is_tuned,
                                      file_prefix = file_prefix, 
                                      col_names = col_names,
                                      seed = seed)

        print(paste0('Fisher extact test pval in trainset:', pvalues[9]))
        print(paste0('pearson correlation pval of tau predictionst:', pvalues[4])) 
        print(paste0('spearman correlation pval of tau predictions in trainset:', pvalues[8]))
        print(paste0('kendall correlation pval of tau predictions in trainset:', pvalues[6]))

        # c(tau.var, fixed.YW.tau.risk, pvalues, NA, NA)
        c(tau.var, fixed.YW.tau.risk, pvalues, calibration_pvalues, perm.pvals)
    }
     
    write.csv(scenario_result[, c(1, 2)], file = paste0('simulation_result/scenario_', num[i], '.observed.var.risk.csv'), quote = F, row.names = F)

    colnames(scenario_result[, -c(1, 2)]) <- c(col_names, calibration_test_names, perm_test_names)
    false_rates <- apply(scenario_result[, -c(1, 2)], 2, function(x) sum(x <= 0.05)/length(x))
    
    write.csv(scenario_result[, -c(1, 2)], paste0('scenario_', num[i], '_result.csv'), quote = F, row.names = F)

    if(verbose){
        print(paste0('# ------------------------------ summary of false discovery rate ', num[i],' ---------------------------------------'))
        print(paste0("dimension of sample size:", sample_size))
        print_false_rate(col_names, false_rates)
    } 
}
