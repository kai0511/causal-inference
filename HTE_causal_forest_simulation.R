library(grf)
library(MASS)
library(dplyr)
library(data.table)
library(survival)
library(survminer)
# library(mppa)
library(doParallel)
library(randomForest)

registerDoParallel(3)

source('survival.imputation.R')  # load impute.survival function
source('hte.validation.R')     # load validation methods for HTE
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
    diff <- (mu.x - tau.x)/2
    probs <- exp(diff) / (1 + exp(diff))
    tx <- sapply(probs, function(i) rbinom(1, size = 1, p = i))
    tx
}

generate.censored.outcome <- function(Y, event.rate){
    cutoff <- quantile(Y, event.rate)
    cen <- rexp(length(Y), rate = 1/cutoff)
    y.censored <- pmin(Y, cen)
    censor <- as.numeric(Y <= y.censored)
    res <- cbind(y.censored, censor)
    colnames(res) <- c('y.censored', 'censor')
    return(res)
}

generate.scenario <- function(n, p, mu.fun, tau.fun, noise.sd = 1, randomized = TRUE, event.rate = 0.2){
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

calculate_mean_surv <- function(rfsrv_obj){
    step_surv <- cbind(rep(1, dim(rfsrv_obj$survival)[1]), rfsrv_obj$survival)
    step_prob <- t(apply(step_surv, 1, function(x) -diff(x)))
    mean_surv_pred <- step_prob %*% rfsrv_obj$time.interest
    # mean_surv_pred <- step_prob %*% diff(c(0, rfsrv_obj$time.interest))
    return(mean_surv_pred)
}

calculate_surv_auc <- function(survival, time.interest){
    step_surv <- cbind(rep(1, dim(survival)[1]), survival)
    step_prob <- t(apply(step_surv, 1, function(x) (x[-1] + x[-length(x)])/2))
    mean_surv_pred <- step_prob %*% diff(c(0, time.interest))
    return(mean_surv_pred)
}

integrate_surv_fun <- function(survival, time.interest){
    step_surv <- cbind(rep(1, dim(survival)[1]), survival)
    timeInterest <- c(0, time.interest)
    meanSurv <- apply(step_surv, 1, function(x, timeInterest){ 
        integrate.xy(timeInterest, x, use.spline= T)}, timeInterest = timeInterest)
    return(meanSurv)
}

fit.cf.splited.datasets <- function(covariates, Y, W, trainId){
    tau.forest.train <- cf.tuned(covariates[trainId, ], Y[trainId], W[trainId], num_trees = 1000)   # run causal forests by default
    tau.pred.train.a <- predict(tau.forest.train, newdata = NULL, estimate.variance = TRUE, num.threads = 8)
    tau.pred.test.a <- predict(tau.forest.train, newdata = covariates[-trainId, ], estimate.variance = TRUE, num.threads = 8)
   
    tau.forest.test <- cf.tuned(covariates[-trainId, ], Y[-trainId], W[-trainId], num_trees = 1000)   # run causal forests by default
    tau.pred.train.b <- predict(tau.forest.test, newdata = NULL, estimate.variance = TRUE, num.threads = 8)
    tau.pred.test.b <- predict(tau.forest.test, newdata = covariates[trainId, ], estimate.variance = TRUE, num.threads = 8)
    
    fitted.models <- list(tau.pred.train.a, tau.pred.test.a, tau.pred.train.b, tau.pred.test.b)
    return(fitted.models)
}


d1 <- generate.scenario(800, 400, 'f8', 'f1')
d9 <- generate.scenario(800, 400, 'f8', 'f1', randomized = FALSE)

d2 <- generate.scenario(800, 400, 'f5', 'f2', noise.sd = 0.25)
d10 <- generate.scenario(800, 400, 'f5', 'f2', noise.sd = 0.25, randomized = FALSE)

d3 <- generate.scenario(800, 300, 'f4', 'f3')
d11 <- generate.scenario(800, 300, 'f4', 'f3', randomized = FALSE)

d4 <- generate.scenario(900, 300, 'f7', 'f4', noise.sd = 0.25)
d12 <- generate.scenario(900, 300, 'f7', 'f4', noise.sd = 0.25, randomized = FALSE)

d5 <- generate.scenario(900, 200, 'f3', 'f5')
d13 <- generate.scenario(800, 200, 'f3', 'f5', randomized = FALSE)

d6 <- generate.scenario(900, 200, 'f1', 'f6')
d14 <- generate.scenario(900, 200, 'f1', 'f6', randomized = FALSE)

d7 <- generate.scenario(1000, 100, 'f2', 'f7', noise.sd = 4)
d15 <- generate.scenario(1000, 100, 'f2', 'f7', noise.sd = 4, randomized = FALSE)

d8 <- generate.scenario(1000, 100, 'f6', 'f8', noise.sd = 4)
d16 <- generate.scenario(1000, 100, 'f6', 'f8', noise.sd = 4, randomized = FALSE)

num <- c(seq(1:8), seq(9, 16))
# no_repeats <- 200
# num <- c(rep(1, no_repeats), rep(9, no_repeats))

print.correlation <- function(fitted.obj, type){
    print(paste0('correlation of ', type,' in :', fitted.obj$estimate, ', its pvalues:', fitted.obj$p.value))
}

print_false_rate <- function(test_names, false_rates){
    num <- length(test_names)

    for(i in seq(num)){
        print(paste0(test_names[i], ': ', false_rates[i]))
    }
}



# --------------------------------------------------------------------------------------
#
# run simulation for each scenario for a number of times, defined by no_repeats below
# collect simulation result wiht matrix and save it in local directory
#
# --------------------------------------------------------------------------------------
no_repeats <- 200
is_print <- FALSE
col_names <- c('variance.p', 'risk.p', 'fixed.W.risk.p', 'fixed.YW.risk.p', 
    'pearson.a.p', 'pearson.b.p', 'spearman.a.p', 'spearman.b.p', 'median.a.p', 'median.b.p')


for(i in seq(length(num))){

    scenario_result <- foreach(i = seq(no_repeats), .combine = 'rbind') %dopar%  {
        
        data <- get(paste0('d', num[i]))
        
        # remove rows containing NAs 
        ind <- apply(data, 1, function(x) sum(is.na(x)) == 0)
        data <- data[ind, ]

        covariates <- data[, -c(1:5)]
        real.Y <- data[, 1]
        tau.real <- data[, 2]
        W <- data[, 3]
        imputed.Y <- data[, 4]
        censor <- data[, 5]   
       

        fitted.obj <- cf.tuned(covariates, imputed.Y, W, num_trees = 1000)
        tau.pred <- predict(fitted.obj, newdata = NULL, estimate.variance = TRUE, num.threads = 8)

        Y.hat <- fitted.obj[["Y.hat"]]
        W.hat <- fitted.obj[["W.hat"]]

        cf.risk = mean((tau.real - tau.pred$predictions)^2)
        
        permutated.p.val <- permutated.pval(tau.real, (W - W.hat) * tau.pred$predictions)
     
        
        Y.hat <- fitted.obj[["Y.hat"]]
        tau.estimated <- imputed.Y - Y.hat
        
        cor.obj <- cor.test(tau.real, tau.pred$predictions, alternative = 'greater')
        risk <- mean((tau.real - tau.pred$predictions)^2)

        if(is_print){
            print(paste0('# ------------------------------ begin scenario ', num[i],' ---------------------------------------'))
            print(paste0('censor rate:', 1 - mean(censor)))
            print(paste0('risk fitted by CF in scenario ', num[i], ': ', cf.risk))

            # print(paste0('variance in scenario ', i, ': ', var(tau.real)))    
            # print(paste0('R square in scenario ', i, ': ', 1 - var(tau.real - tau.pred$predictions) / var(tau.real)))
            print(paste0('p.val of permutating either real tau from simulation or (W - W.hat) * tau.pred$predictions in scenario ', num[i], ':', permutated.p.val))
            print(paste0('The estimated measure of correlation between real tau and predicted tau in scenario ', num[i], ': ', cor.obj$estimate, '; its pval:', cor.obj$p.value ))
            print(paste0('risk of mean squared difference between real tau and predicted tau value in scenario ', num[i], ': ', risk))
        }

        # *****************************************************************************************
        # permutate covariates to validate HET esimation
        # *****************************************************************************************
        no.obs <- length(imputed.Y)
        tau.1 <- imputed.Y - Y.hat
        tau.2 <- (W - W.hat) * tau.pred$predictions
        permutated.p.val <- permutated.pval(tau.1, tau.2)

        if(is_print){
            print(paste0('p.val of permutating either (Y - Y.hat) or [(W - W.hat) * tau.pred$predictions] in scenario ', num[i], ':', permutated.p.val))
            print(test_calibration(fitted.obj))
        }
        
        
        # print('# ------------------------ result of calibration by maunual regression -----------------------------------')
        # tau.bar <- mean(tau.pred$predictions)
        # pred.ret <- data.frame(outcome = tau.1, tau.diff = tau.pred$predictions - tau.bar, tau.avg = rep(tau.bar, length(tau.1)))
        # print(summary(lm(outcome ~ tau.diff + tau.avg + 0, data = pred.ret)))
       
        tau.risk <- assess.explained.tau.risk(fitted.obj, imputed.Y, W)
        fixed.W.tau.risk <- assess.explained.tau.fixed.W.risk(fitted.obj, imputed.Y, W, W.hat)    
        fixed.YW.tau.risk <- assess.explained.tau.fixed.YW.risk(fitted.obj, imputed.Y, Y.hat, W, W.hat)    

        perm.p.vals <- permutate.covariates.testing(covariates, imputed.Y, Y.hat, W, W.hat, tau.risk, fixed.W.tau.risk, fixed.YW.tau.risk, var(tau.pred$predictions), num_trees = 1000, num.strap = 500) 
        
        if(is_print){
            print('# ------------------------------ result of permutation covariates ---------------------------------------')
            print(paste0('p-values of variance testing by permutating covariates: ', perm.p.vals[1]))
            print(paste0('p-values of tau risk by permutating covariates without fixing YW: ', perm.p.vals[2]))
            print(paste0('p-values of tau risk by permutating covariates with fixing W: ', perm.p.vals[3]))
            print(paste0('p-values of tau risk by permutating covariates with fixing YW: ', perm.p.vals[4]))
        }
        

        # *****************************************************************************************
        # permutate treatment to validate HET esimation
        # *****************************************************************************************
        # Y.residual <- imputed.Y - Y.hat - tau.pred$predictions

        # treatment.var <- var(Y.residual[W == 1])
        # control.var <- var(Y.residual[W == 0])
        # obs.t.var <- treatment.var/control.var

        # perm.treament.pval <- permutate.treatment.testing(covariates, imputed.Y, W, obs.t.var, num_trees = 1000, num.strap = 500)
        # print(paste0('p-values of t.var statistics by permutating W:', perm.treament.pval))

        # *****************************************************************************************
        # replication testing by spliting dataset into two pieces.
        # *****************************************************************************************
        no.obs <- dim(covariates)[1]
        trainId <- sample(1:no.obs, floor(no.obs/2), replace = F)
        
        fitted.models <- fit.cf.splited.datasets(covariates, imputed.Y, W, trainId)
        
        # save(fitted.models, file = paste0('fitted.model', i, '.RData')) 
        cor.obj.a <- cor.test(fitted.models[[1]]$predictions, fitted.models[[4]]$predictions, alternative = 'greater', method = 'pearson', use="na.or.complete")
        cor.obj.b <- cor.test(fitted.models[[2]]$predictions, fitted.models[[3]]$predictions, alternative = 'greater', method = 'pearson', use="na.or.complete")
        
        rank.cor.obj.a <- cor.test(fitted.models[[1]]$predictions, fitted.models[[4]]$predictions, alternative = 'greater', method = 'spearman', use="na.or.complete")
        rank.cor.obj.b <- cor.test(fitted.models[[2]]$predictions, fitted.models[[3]]$predictions, alternative = 'greater', method = 'spearman', use="na.or.complete")
        
        

        t.test.pval.a <- median.t.test(fitted.models[[1]]$predictions, fitted.models[[4]]$predictions)
        t.test.pval.b <- median.t.test(fitted.models[[2]]$predictions, fitted.models[[3]]$predictions)

        if(is_print){
            print.correlation(cor.obj.a, 'pearson')
            print.correlation(cor.obj.b, 'pearson')
            print.correlation(rank.cor.obj.a, 'spearman')
            print.correlation(rank.cor.obj.b, 'spearman')
            print(paste0('median test pval in trainset:', t.test.pval.a, '; corresponding pval in testset:', t.test.pval.b))
        }
        
        # connect results for saving
        correlation.res <- c(cor.obj.a$p.value, cor.obj.b$p.value, rank.cor.obj.a$p.value, rank.cor.obj.b$p.value, t.test.pval.a, t.test.pval.b)
        c(perm.p.vals, correlation.res)
    }

    
    colnames(scenario_result) <- col_names
    false_rates <- apply(scenario_result, 2, function(x) sum(x <= 0.05)/length(x))

    if(is_print){
        print(paste0('# ------------------------------ summary of false discovery rate ', num[i],' ---------------------------------------'))
        print_false_rate(test_names, false_rates)
    } 

    write.csv(scenario_result, paste0('scenario_', num[i], '_result.csv'), quote = F, row.names = F)

}
write.csv(res, file = paste0('permutation.pvalues.csv'), quote = F, row.names = F)
