# set eight scenarios for estimating HTE, details refers to Powers et al. 2018

f1 = function(x) {
    rep(0, nrow(x))
}

f2 = function(x) {
    5 * (x[, 1] > 1) - 1 - (5*pnorm(-1) -1)
} 
    

f3 = function(x) {
    5 * x[, 1]
}

f4 = function(x) {
    1 * x[,2] * x[,4] * x[,6] + 
    2 * x[,2] * x[,4] * (1-x[,6]) +
    3 * x[,2] * (1-x[,4]) * x[,6] + 
    4 * x[,2] * (1-x[,4]) * (1-x[,6]) +
    5 * (1-x[,2]) * x[,4] * x[,6] + 
    6 * (1-x[,2]) * x[,4] * (1-x[,6]) +
    7 * (1-x[,2]) * (1-x[,4]) * x[,6] + 
    8 * (1-x[,2]) * (1-x[,4]) * (1-x[,6]) - 
    5 - (-0.5)
}

f5 = function(x) {
    x[, 1] + x[, 3] + x[, 5] + x[, 7] + x[, 8] + x[, 9] -
    0.5 
}

f6 = function(x) {
    4 * (x[,1]>1) * (x[,3]>0) + 4 * (x[,5]>1) * (x[,7]>0) + 2 * x[,8] * x[,9] - 1 - (4*pnorm(-1)-1)
}

f7 = function(x){
  ( x[, 1]^2 + 
    x[, 2] + 
    x[, 3]^2 + 
    x[, 4] + 
    x[, 5]^2 + 
    x[, 6] + 
    x[, 7]^2 +
    x[, 8] + 
    x[, 9]^2) / sqrt(2) - 5 -
    (7/sqrt(2) - 5)
}

f8 = function(x){
    (f4(x) + f5(x)) / sqrt(2)
}    

f9 = function(x) {
    (x[,1])^2 - 1
}

