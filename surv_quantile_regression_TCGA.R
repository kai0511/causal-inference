require(quantreg)
require(dplyr)

for (tau in (1:4)/5){
    fitted.obj <- crq(Curv(surv.time,censor) ~ mat, tau = tau, method = "Pow")
}
