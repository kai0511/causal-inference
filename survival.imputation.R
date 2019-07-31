
library(survival)

impute.survival <- function(surv.time, censor, cluters = NULL){
    # note that cluters should be a factor of length of surv.time
    if(is.null(cluters)){
        max.censored <- max(surv.time[censor == 0])
        censor[surv.time == max.censored] <- 1
        imputed.surv.times <- compute.surv.time(surv.time, censor)
    }else{
        imputed.surv.times <- rep(0, length(cluters))
        cluter_id <- levels(cluters)

        for(type in cluter_id){
            max.censored <- max(surv.time[censor == 0 & cluters == type])
            censor[surv.time == max.censored & cluters == type] <- 1
            
            idx <- (cluters == type)
            idx.surv.time <- compute.surv.time(surv.time[idx], censor[idx])
            imputed.surv.times[idx] <- idx.surv.time
        }
        return(imputed.surv.times)
    }
}

compute.surv.time <- function(surv.time, censor){
    # @surv.time: a vector of times, if the subject is alive, then it's NA; else, it's a number, indicated survival time
    # @censor: The status indicator, normally 0 = alive, 1 = dead, required by "Surv" function
    # need to change censoring status with the longest censoring to non-censor and assume death

    surv.obj <- Surv(surv.time, censor)
    fit.obj <- survfit(surv.obj ~ 1)
    step.surv <- -diff(c(1, fit.obj$surv))
    log.step <- log(fit.obj$time) * step.surv

    # impute the survival time of censored obs.(censor == 0), unknown to us.
    log.times <- log(surv.time)

    surv.times <- sapply(1:length(surv.time), function(i){
        ind <- which(fit.obj$time == surv.time[i])
        
        # It may be caused by a bug of survival package, sometimes we cannot find specific survival time 
        # from "time" list of fitted object by Surv function.
        if (length(ind) == 0){
            return(NA)
        }
        
        if(censor[i] == 0){
            return(sum(log.step[fit.obj$time > fit.obj$time[ind]]) / fit.obj$surv[ind])
        }else{
            return(log.times[i])
        }
    })
    return(surv.times)
}