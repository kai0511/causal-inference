library(grf)
library(BART)
library(ranger)
library(randomForestSRC)
library(randomForest)


cf <- function(covarites, Y, W, cluster_id = NULL, frac_selected = 0.333){
    # set and run casual forest
    tau.forest <- causal_forest(X = covarites,
                                Y = Y,
                                W = W,  # merged$AgeOnset
                                sample.fraction = 0.5,  # When confidence intervals are enabled, the sampling fraction must be less than 0.5.
                                num.trees = 2000,
                                min.node.size = 10,   # by default
                                mtry = round(ncol(covarites)* frac_selected) , # default is ceiling(sqrt(num.col.covariates) + 20)
                                clusters = cluster_id,
                                seed = 0,
                                num.threads = 20)
    return(tau.forest)
}

cf.tuned <- function(covarites, Y, W, cluster_id = NULL, num_trees = 2000, is_tune = TRUE, seed = 0){
    # set and run casual forest
    tau.forest <- causal_forest(X = covarites,
                                Y = Y,
                                W = W,
                                # sample.fraction = 0.5,  # When confidence intervals are enabled, the sampling fraction must be less than 0.5.
                                num.trees = num_trees,
                                clusters = cluster_id, 
                                # min.node.size = 10,   # by default
                                # mtry = round(ncol(covarites)* frac_selected) , # default is ceiling(sqrt(num.col.covariates) + 20)
                                tune.parameters = is_tune,
                                seed = seed,
                                num.threads = 20)
    return(tau.forest)
}


# cf <- function(covarites, Y, W, params){
#     # This function is for grf of old version, and params can be obtain by using tune_causal_forest, which is changed in usage in new version. 
#     tau.forest <- causal_forest(X = covarites,
#                                 Y = Y,
#                                 W = W, 
#                                 sample.fraction = as.numeric(params["sample.fraction"]), 
#                                 num.trees = 1000,
#                                 min.node.size = as.numeric(params["min.node.size"]),  
#                                 mtry =  as.numeric(params["mtry"]), # default is ceiling(sqrt(num.col.covariates) + 20)
#                                 alpha = as.numeric(params["alpha"]),
#                                 imbalance.penalty = as.numeric(params["imbalance.penalty"]),
#                                 seed = 0,
#                                 num.threads = 10)
#     return(tau.forest)
# }


# rfscr_hte <- function(covarites, Y, W, censor){
#     #as required by the method, the censoring variable must be coded with 0 reserved for censoring, 
#     # and (usually) 1 = death (event).
# 
#     treatment = (W == 1)
#     tau_pred <- rep(NA, dim(covarites)[1])
#     colnames(covarites) <- paste0('V', seq(dim(covarites)[2]))
# 
#     surv_data = as.data.frame(cbind(Y, censor, covarites))
# 
#     colnames(surv_data)[1:2] <- c('surv_time', 'censor')
# 
#     rfsrc_treatment <- rfsrc(Surv(surv_time, censor) ~ .,
#                              data = surv_data[treatment, ],
#                              splitrules = 'logrankscore',
#                              importance = FALSE)
# 
#     rfsrc_control <- rfsrc(Surv(surv_time, censor) ~ .,
#                            data = surv_data[!treatment, ],
#                            splitrules = 'logrankscore',
#                            importance = FALSE)
# 
#     rfsrc_treatment_pred = predict(rfsrc_treatment, surv_data[!treatment,])
#     rfsrc_control_pred = predict(rfsrc_control, surv_data[treatment,])
# 
#     rfsrc_treatment_surv <- calculate_surv_auc(rfsrc_treatment_pred$survival, rfsrc_treatment_pred$time.interest)
#     rfsrc_control_surv <- calculate_surv_auc(rfsrc_control_pred$survival, rfsrc_control_pred$time.interest)
# 
#     tau_control <- log(rfsrc_treatment_surv) - log(Y[!treatment])
#     tau_treatment <- log(Y[treatment]) - log(rfsrc_control_surv)
# 
#     tau_pred[treatment] <- tau_treatment
#     tau_pred[!treatment] <- tau_control
# 
#     return(tau_pred)
# }
# 
# ranger_hte <- function(covariates, Y, W, censor){
#     # as required by the method, the censoring variable must be coded with 0 reserved for censoring, 
#     # and (usually) 1 = death (event).
# 
#     treatment = (W == 1)
#     tau_pred <- rep(NA, dim(covariates)[1])
# 
#     colnames(covariates) <- paste0('V', seq(dim(covariates)[2]))
# 
#     surv_data = as.data.frame(cbind(Y, censor, covariates))
# 
#     colnames(surv_data)[1:2] <- c('surv_time', 'censor')
# 
#     rg_treatment <- ranger(Surv(surv_time, censor) ~ .,
#                                  data = surv_data[treatment,])
# 
#     rg_control <- ranger(Surv(surv_time, censor) ~ .,
#                                data = surv_data[!treatment,])
# 
#     rg_treatment_pred = predict(rg_treatment, surv_data[!treatment,])
#     rg_control_pred = predict(rg_control, surv_data[treatment,])
# 
#     rg_treatment_surv <- calculate_surv_auc(rg_treatment_pred$survival, rg_treatment_pred$unique.death.times)
#     rg_control_surv <- calculate_surv_auc(rg_control_pred$survival, rg_control_pred$unique.death.times)
# 
#     tau_control <- log(rg_treatment_surv) - log(Y[!treatment])
#     tau_treatment <- log(Y[treatment]) - log(rg_control_surv)
# 
#     tau_pred[treatment] <- tau_treatment
#     tau_pred[!treatment] <- tau_control
# 
#     return(tau_pred)
# }
# 
# bart.hte <- function(covarites, Y, W, censor, ndpost = 2000, nskip = 2000){
#     # set and run bart with time-to-event outcomes
#     treatment.assgin <- (W == 1)
#     tau.pred <- rep(NA, dim(covarites)[1])
# 
#     bart.treatment.obj <- abart(covarites[treatment.assgin, ],
#                                 Y[treatment.assgin],
#                                 censor[treatment.assgin],
#                                 covarites[!treatment.assgin, ],
#                                 ndpost = ndpost,  # The number of posterior draws returned.
#                                 nskip = nskip)    # Number of MCMC iterations to be treated as burn in
# 
#     bart.control.obj <- abart(covarites[!treatment.assgin, ],
#                               Y[!treatment.assgin],
#                               censor[!treatment.assgin],
#                               covarites[treatment.assgin, ],
#                               ndpost = ndpost,     # same as above
#                               nskip = nskip)       # save as above
# 
#     tau.control <- log(bart.treatment.obj$yhat.test.mean) - log(Y[!treatment.assgin])
#     tau.treatment <- log(Y[treatment.assgin]) - log(bart.control.obj$yhat.test.mean)
#     print(head(bart.treatment.obj$yhat.test.mean))
#     print(head(log(bart.treatment.obj$yhat.test.mean)))
#     tau.pred[treatment.assgin] <- tau.treatment
#     tau.pred[!treatment.assgin] <- tau.control
# 
#     return(tau.pred)
# }
# 
# surv_bart_hte <- function(covariates, Y, W, censor){
#     tau_pred <- rep(NA, dim(covariates)[1])
# 
#     trainset <- cbind(W, covariates)
# 
#     pre.fit <- surv.pre.bart(Y, censor, x.train = trainset, x.test = trainset)
# 
#     unique.times <- pre.fit$K
#     test.covar <- pre.fit$tx.test[, -2]
#     long.W <- pre.fit$tx.test[, 2]
#     unique.death.times <- pre.fit$tx.test[1: unique.times, 1]
#     # long.W <- as.vector(sapply(W, function(x) rep(x, unique.times)))
# 
#     treatment.post <- mc.surv.bart(x.train = covariates[W == 1, ],
#                                    times = Y[W == 1],
#                                    delta = censor[W == 1],
#                                    mc.cores = 8,
#                                    seed = 99)
# 
#     control.post <- mc.surv.bart(x.train = covariates[W != 1, ],
#                                  times = Y[W != 1],
#                                  delta = censor[W != 1],
#                                  mc.cores = 8,
#                                  seed = 99)
# 
#     treatment.pred <- predict(treatment.post, newdata=test.covar[long.W != 1,], mc.cores=8)
#     control.pred <- predict(control.post, newdata=test.covar[long.W == 1,], mc.cores=8)
# 
#     treatment.survival <- t(matrix(treatment.pred$surv.test.mean, nrow = unique.times))
#     control.survival <- t(matrix(control.pred$surv.test.mean, nrow = unique.times))
# 
#     bart_treatment_surv <- calculate_surv_auc(treatment.survival, unique.death.times)
#     bart_control_surv <- calculate_surv_auc(control.survival, unique.death.times)
# 
#     tau.pred[W != 1] <- log(bart_treatment_surv) - log(Y[W != 1])
#     tau.pred[W == 1] <- log(Y[W == 1]) - log(bart_control_surv)
#     return(tau.pred)
# }
# 
