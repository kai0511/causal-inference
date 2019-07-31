require(spatstat)

# some methods for validate causal models
# for bootstrap testing refer to "a comparison of methods for model selection when estimating individual treatment effects"
# permutate.covariates.testing is to test whether covariates contribute to variation in heterogeneous treatment effect.
# for permutated.pval ref. equation 9 to the paper above, we permute y_i - m_hat(x_i) to check our fitting


# Shifted kolmogorov-smirnov statistic
# Calculate KS distance between Y0 and Y1 shifted by estimated tau.

weighted.SKS.stat <- function(Y, W, menbership){
    Y1 = Y[W == 1]
    Y0 = Y[W == 0]
    
    Y1.star = Y1 - mean(Y1)
    Y0.star = Y0 - mean(Y0)
    
    unique.points = c(Y1.star, Y0.star)
    
    # calculate emprical weighted cum distribution
    Fn1 = ewcdf(Y1.star, menbership[W == 1])
    Fn0 = ewcdf(Y0.star, 1 - menbership[W == 0])
    
    difference = Fn1(unique.points) - Fn0(unique.points)
    
    return(max(abs(difference)))
}

weighted.var.stat <- function(Y, W, menbership){
    Y1 <- Y[W == 1]
    Y0 <- Y[W == 0]
    
    Y1.star <- Y1 - mean(Y1)
    Y0.star <- Y0 - mean(Y0)
    
    # calculate emprical weighted cum distribution
    var1 <- mean(Y1.star^2 %*% menbership[W == 1])
    var0 <- mean(Y0.star^2 %*% (1 - menbership[W == 0]))
    
    return(var1 / var0)
}

cf_SKS_test <- function(fitted.obj, repeats = 5000){
    Y <- fitted.obj[["Y.orig"]]
    Y.hat <- fitted.obj[["Y.hat"]]
    residuals <- Y - Y.hat
    
    W <- forest[["W.orig"]]
    # here we regard w.hat as weights, which measure the menbership of each obs. belonging to treatment group
    W.hat <- forest[["W.hat"]]

    # compute weighted SKS.stat 
    t_obs <- weighted.SKS.stat(residuals, W, W.hat)

    perm_obs <- replicate(repeats, {
        perm.residuals = sample(residuals, length(Y))
        t = SKS.stat(perm.residuals, W, W.hat)
        t
    })
    mean(t_obs >= perm_obs)
}

cf_variance_test <- function(fitted.obj, repeats = 5000){
    Y <- fitted.obj[["Y.orig"]]
    Y.hat <- fitted.obj[["Y.hat"]]
    residuals <- Y - Y.hat
    
    W <- forest[["W.orig"]]
    # here we regard w.hat as weights, which measure the menbership of each obs. belonging to treatment group
    W.hat <- forest[["W.hat"]]

    # compute weighted SKS.stat 
    t_obs <- weighted.var.stat(residuals, W, W.hat)

    perm_obs <- replicate(repeats, {
        perm.residuals = sample(residuals, length(Y))
        t = weighted.var.stat(perm.residuals, W, W.hat)
        t
    })
    mean(t_obs >= perm_obs)
}

assess.tau.risk <- function(fitted.obj, Y, treatment, constant.tau = NULL){
    # calculate tau-risk_R
    tau.pred <- predict(fitted.obj, newdata = NULL, estimate.variance = TRUE, num.threads = 8)
    Y.hat <- fitted.obj[["Y.hat"]]
    W.hat <- fitted.obj[["W.hat"]]
    tau.1 <- Y - Y.hat

    if (is.null(constant.tau)){
        tau.2 <- (treatment - W.hat) * tau.pred$predictions
    }else {
        tau.2 <- (treatment - W.hat) * constant.tau
    }
    
    risk <- mean((tau.1 - tau.2)^2)
    return(risk)
}

assess.explained.tau.risk <- function(fitted.obj, Y, treatment, constant.tau = NULL){
    # calculate tau-risk_R
    tau.pred <- predict(fitted.obj, newdata = NULL, estimate.variance = TRUE, num.threads = 8)
    Y.hat <- fitted.obj[["Y.hat"]]
    W.hat <- fitted.obj[["W.hat"]]
    tau.1 <- Y - Y.hat

    if (is.null(constant.tau)){
        mean.tau <- (treatment - W.hat) * mean(tau.pred$predictions)
        tau.2 <- (treatment - W.hat) * tau.pred$predictions
    } else {
        mean.tau <- (treatment - W.hat) * mean(constant.tau)
        tau.2 <- (treatment - W.hat) * constant.tau
    }
    total.risk <- (tau.1 - mean.tau)^2
    unexplained.risk <- (tau.1 - tau.2)^2

    risk.explained <- mean(total.risk - unexplained.risk)
    return(risk.explained)
}

assess.explained.tau.fixed.W.risk <- function(fitted.obj, Y, treatment, W.hat, constant.tau = NULL){
    # calculate tau-risk_R
    tau.pred <- predict(fitted.obj, newdata = NULL, estimate.variance = TRUE, num.threads = 8)
    Y.hat <- fitted.obj[["Y.hat"]]
    #W.hat <- fitted.obj[["W.hat"]]
    tau.1 <- Y - Y.hat

    if (is.null(constant.tau)){
        mean.tau <- (treatment - W.hat) * mean(tau.pred$predictions)
        tau.2 <- (treatment - W.hat) * tau.pred$predictions
    }else {
        mean.tau <- (treatment - W.hat) * mean(constant.tau)
        tau.2 <- (treatment - W.hat) * constant.tau
    }

    risk.total <- (tau.1 - mean.tau)^2
    risk.unexplained <- (tau.1 - tau.2)^2
    risk.explained <- mean(risk.total - risk.unexplained)
    return(risk.explained)
}

assess.explained.tau.fixed.YW.risk <- function(fitted.obj, Y, Y.hat,  treatment, W.hat, constant.tau = NULL){
    # calculate tau-risk_R
    tau.pred <- predict(fitted.obj, newdata = NULL, estimate.variance = TRUE, num.threads = 8)
    # Y.hat <- fitted.obj[["Y.hat"]]
    # W.hat <- fitted.obj[["W.hat"]]
    tau.1 <- Y - Y.hat

    if (is.null(constant.tau)){
        mean.tau <- (treatment - W.hat) * mean(tau.pred$predictions)
        tau.2 <- (treatment - W.hat) * tau.pred$predictions
    }else {
        mean.tau <- (treatment - W.hat) * mean(constant.tau)
        tau.2 <- (treatment - W.hat) * constant.tau
    }
    
    risk.total <- (tau.1 - mean.tau)^2
    risk.unexplained <- (tau.1 - tau.2)^2
    risk.explained <- mean(risk.total - risk.unexplained)
    return(risk.explained)
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

bootstrap.testing <- function(covariates, Y, treatment, tau.risk, constant.risk, tau.var, num_trees = 2000, num.strap = 1000, frac = 1){
    # calculate average difference in mean survival time in treatment v.s non-treatment groups
    # including two subsampling testings, so two pvalues vill be returned in a vector.
    subsampled.result <- foreach(i = seq(num.strap), .combine = 'rbind') %dopar%  {
        samp <- sample(dim(covariates)[1], frac * dim(covariates)[1], replace = F)
        sampled.covariates <- covariates[samp, ]
        sampled.Y <- Y[samp]
        sampled.treatment <- treatment[samp]


        # params <- tune_causal_forest(sampled.covariates, sampled.Y, sampled.treatment)$params
        tau.forest <- cf.tuned(sampled.covariates, sampled.Y, sampled.treatment, num_trees = num_trees)
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

    var.pval <- mean((subsampled.var - tau.var) >= tau.var)

    return(c(subsampling.pval, var.pval))
}

stratified.permutation <- function(cluster_id){
    cluters <- levels(cluster_id)
    obs <- seq(length(cluster_id))

    samp <- rep(0, length(cluster_id))
    for(cls in cluters){
        cls_idx <- obs[cluster_id == cls]
        samp[cls_idx] <- sample(cls_idx, length(cls_idx), replace = F)
    }
    return(samp)
}

permutate.covariates.testing <- function(covariates, Y, Y.hat, treatment, W.hat, tau.risk, fixed.W.tau.risk, fixed.YW.tau.risk, tau.var, cluster_id = NULL, num_trees = 2000, num.strap = 1000, frac = 1, mut = NULL, seed =333){
    # permute covariates to investigate whether covariates contribute to HTE
    # return tau.risk for each permutation
    set.seed(seed)

    perm.risk.mat <- foreach(i = seq(num.strap), .combine = 'rbind') %dopar%  {
        set.seed((i-1) * num.strap + seed)

        # if cluster_id is provided, we do a stratified permutation to keep cluter heterogeneity present
        if(is.null(cluster_id)){
            samp <- sample(dim(covariates)[1], frac * dim(covariates)[1], replace = F)
        }else{
            samp <- stratified.permutation(cluster_id)
        }

        sampled.covariates <- covariates[samp, ]
        perm.tau.forest <- cf.tuned(sampled.covariates, Y, treatment, num_trees = num_trees)
        perm.var = var(perm.tau.forest$predictions)
        perm.risk <- assess.explained.tau.risk(perm.tau.forest, Y, treatment) ##the larger the better
        fixed.W.perm.risk <- assess.explained.tau.fixed.W.risk(perm.tau.forest, Y, treatment, W.hat)
        fixed.YW.perm.risk <- assess.explained.tau.fixed.YW.risk(perm.tau.forest, Y, Y.hat, treatment, W.hat)
        c(perm.var, perm.risk, fixed.W.perm.risk, fixed.YW.perm.risk)
    }

    if(!is.null(mut)){    
        observed_risk <- c(tau.var, tau.risk, fixed.YW.tau.risk)
        write.csv(perm.risk.mat, paste0('mutation_result/', mut, '.permutation.risk.result.fixed.W.csv'), quote = F, row.names = F)
        write.csv(observed_risk, paste0('mutation_result/', mut, '.observed.risk.result_fixed.W.csv'), quote = F, row.names = F)
    }
	
    permute.var.p = mean(perm.risk.mat[,1] >  tau.var)
    permute.risk.p = mean(perm.risk.mat[,2] > tau.risk)
    fixed.W.permute.risk.p = mean(perm.risk.mat[,3] > fixed.W.tau.risk)
    fixed.YW.permute.risk.p = mean(perm.risk.mat[,4] > fixed.YW.tau.risk)
    
    return(c(permute.var.p, permute.risk.p, fixed.W.permute.risk.p, fixed.YW.permute.risk.p))
}


permutate.treatment.testing <- function(covariates, Y, treatment, t.var, cluster_id = NULL, num_trees = 2000, num.strap = 1000, frac = 1, mut = NULL, seed =333){
    # permute treatment to investigate the ratio of variance under treatment and control 
    # return the ratio of variance for each permutation

    perm.t.vars <- foreach(i = seq(num.strap), .combine = 'c') %dopar%  {
        set.seed((i-1) * num.strap + seed)

        # if cluster_id is provided, we do a stratified permutation to keep cluter heterogeneity present
        if(is.null(cluster_id)){
            samp <- sample(dim(covariates)[1], frac * dim(covariates)[1], replace = F)
        }else{
            samp <- stratified.permutation(cluster_id)
        }

        sampled.treatment <- treatment[samp]

        perm.tau.forest <- cf.tuned(covariates, Y, sampled.treatment, num_trees = num_trees)
        tau.pred <- predict(perm.tau.forest, newdata = NULL, estimate.variance = TRUE, num.threads = 8)

        Y.hat <- fitted.obj[["Y.hat"]]
        Y.residual <- Y - Y.hat - perm.tau.forest$predictions

        treatment.var <- var(Y.residual[sampled.treatment == 1])
        control.var <- var(Y.residual[sampled.treatment == 0])
        perm.t.var <- treatment.var/control.var
        perm.t.var
    }

    # Here assume the ratio after permutation is greater than that without permutation because greater unexplananed heterogeneity.
    # so we expect that t.var observed is less than most of t.var after permutation.
    return(mean(perm.t.vars < t.var))
}

permutate.covariates.excessRisk.testing <- function(covariates, Y, treatment, tau.excess.risk, tau.var, num.strap = 1000, frac = 1, seed =333){
    # permute covariates to investigate whether covariates contribute to HTE
    # return tau.risk for each permutation
    set.seed(seed)

    perm.risk.mat <- foreach(i = seq(num.strap), .combine = 'rbind') %dopar%  {
        set.seed((i-1)*num.strap +seed)
        samp <- sample(dim(covariates)[1], frac * dim(covariates)[1], replace = F)
        sampled.covariates <- covariates[samp, ]
        perm.tau.forest <- cf(sampled.covariates, Y, treatment)
        perm.risk <- assess.explained.tau.risk(perm.tau.forest, Y, treatment) ##the larger the better
        perm.var = var(perm.tau.forest$predictions)
        c(perm.risk, perm.var)
    }

    permuteExcessRisk.p = mean(perm.risk.mat[,1] > tau.excess.risk )
    permute.var.p = mean(perm.risk.mat[,2] >  tau.var)
    return(c(permuteExcessRisk.p, permute.var.p))
}

permutate.covariates.fixed.YW.excessRisk.testing <- function(covariates, Y, Y.hat, treatment, W.hat, tau.excess.risk, num.strap = 1000, frac = 1, seed =333){
    # permute covariates to investigate whether covariates contribute to HTE
    # return tau.risk for each permutation
    set.seed(seed)

    perm.risk <- foreach(i = seq(num.strap), .combine = 'c') %dopar% {
        set.seed((i-1)*num.strap +seed)
        samp <- sample(dim(covariates)[1], frac * dim(covariates)[1], replace = F)
        sampled.covariates <- covariates[samp, ]
        perm.tau.forest <- cf(sampled.covariates, Y, treatment)
        perm.risk <- assess.explained.tau.fixed.YW.risk(perm.tau.forest, Y, Y.hat, treatment, W.hat) ##the larger the better 
        perm.risk
    }
    return(mean(perm.risk > tau.excess.risk))
}

fisher.exact.test <- function(x, y, sig = 0.1){
    crosstabular  <- table(as.factor(x > sig), as.factor(y > sig))
    unname(crosstabular)
    if(dim(crosstabular)[1] == 2 && dim(crosstabular)[2] == 2){
        fitted <- fisher.test(crosstabular, alternative = 'greater')
        return(fitted$p.value)
    }else{
        return(NA)
    }
}

median.t.test <- function(x, y){
    x.median <- median(x)
    y.median <- median(y)

    if(sum(x <= y.median) < 10 | sum(x > y.median) < 10){
        x.p.value = NA
    }else{
        x.fit <- t.test(x[x <= y.median], x[x > y.median], alternative = 'less')
        x.p.value <- x.fit$p.value
    }

    if(sum(y <= x.median) < 10 | sum(y > x.median) < 10){
        y.p.value = NA
    }else{
        y.fit <- t.test(y[y <= x.median], y[y > x.median], alternative = 'less')
        y.p.value <- y.fit$p.value
    }
    return(c(x.p.value, y.p.value))
}

quantile.t.test <- function(x, y, quantile = 0.5){
    x.quantile <- quantile(x, quantile)
    y.quantile <- quantile(y, quantile)

    if(sum(x <= y.quantile) < 5 | sum(x > y.quantile) < 5){
        x.p.value = NA
    }else{
        x.fit <- t.test(x[x <= y.quantile], x[x > y.quantile], alternative = 'less')
        x.p.value <- x.fit$p.value
    }

    if(sum(y <= x.quantile) < 5 | sum(y > x.quantile) < 5){
        y.p.value = NA
    }else{
        y.fit <- t.test(y[y <= x.quantile], y[y > x.quantile], alternative = 'less')
        y.p.value <- y.fit$p.value
    }
    return(c(x.p.value, y.p.value))
}

simes.partial <- function(u, pvec)  {  ##equation 2 in the above reference
    n = length(pvec)
    pvec = pvec[order(pvec)]  ##sort pval in ascending order

    corr.p = numeric(n-u+1)

    for (i in 1: (n - u + 1) ){
        corr.p[i] = (n-u+1)/i  * pvec[u-1+i]
    }

    partial.conj.p = min(corr.p)
    return(partial.conj.p)
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
