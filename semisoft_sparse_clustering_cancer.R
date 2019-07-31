
setwd("/exeh_3/kai/data/heterogenous_treatment_effects")




load("TCGA_PanCanAtlas_No_Of_Mutations_NoSilentMut.rda")
mutation <- mut_adj.rev[,-1]
rownames(mutation) <- mut_adj.rev[, 1]



get_similarity_matrix <- function(expr, type="log") {

    ## for raw counts, normalize and log-transform
    if (type == "count") {
        depths = rowSums(expr)
        A = expr %*% t(expr) - diag(depths)
    } else if (type == "log") {
        ## center by genes is helpful
        expr = scale(expr, center=TRUE, scale=FALSE)
        A = expr %*% t(expr)
    } else {
        stop("data type must be either 'log' or 'count'.\n")
    }
    return (A) 
}