require(arules)
require(foreach)
require(doParallel)

registerDoParallel(40)

# Write function to do arules random sampling

RSarules <- function(data, M, ig, rhs, lhs_offset = NULL) {
    # Check the inputs are of appropriate types
    if (!is.matrix(data))
        stop("The dataset must be of matrix type.")

    if (max(data) > 1)
        stop("The dataset must be contain only 0 and 1.")

    if (!is.numeric(rhs)) {
        stop("rhs must be of numeric type.")
    }

    if (is.numeric(rhs)) {
        if (length(rhs) > 1)
            stop("The rhs can only contain one item.")
    }

    if (!is.null(lhs_offset)) {
        if (!is.numeric(lhs_offset)) {
            stop("The offset must be of numeric type.")
        }
    }

    if (is.null(M)) {
        "The sample size must be specified."
    }

    if (!is.null(M)) {
        if (!is.numeric(M)) {
            stop("The number of association rules to be sampled must be of numeric type.")
        }
    }

    if (is.null(ig)) {
        "The value of the tuning parameter must be specified."
    }

    if (!is.null(ig)) {
        if (!is.numeric(ig)) {
            stop("The value for the tuning parameter must be of numeric type.")
        }
    }

    # Rearrage the transaction data set 'dd'
    l_1 <- dim(data)[2]
    colnames(data) <- c(paste0("I", 1:l_1))

    if (is.null(lhs_offset)) {
        dd_col <- append(c(1:l_1)[-rhs], rhs)
    }

    if (!is.null(lhs_offset)) {
        dd_col <- append(c(1:l_1)[-append(lhs_offset, rhs)], rhs)
    }

    dd <- data[, dd_col]

    ## dimenstions of dd
    n <- dim(dd)[1]
    l <- dim(dd)[2]

    # Use functions in the R package 'arules' library(arules) a function supp to
    # compute support of a given itemset
    supp <- function(itemset, dataset) {
        # itemset: a given itemset, which is a vector with 0 and 1.
        # ddt: the transaction dataset
        cnt <- 0
        # if (sum(itemset) != 0) {
        #     for(i in 1: dim(dataset)[1]) {
        #         # check whether it is contained by the row
        #         if (length(which(dataset[i, ] - itemset < 0)) == 0){
        #             cnt <- cnt + 1
        #         }
        #     }
        # }

        if (sum(itemset) != 0) {
            cnt <- sum(apply(dataset, 1, function(x, y) sum((x - y) < 0) == 0, itemset))
        }
        return(cnt / dim(dataset)[1])
    }

    # ***************************************************
    #             Gibbs Sampler
    # ****************************************************

    # 2.1 Initial values
    set.seed(100)
    J <- matrix(c(rbinom(l - 1, 1, 0.5), 0), nrow = 1)
    S <- matrix(c(0), nrow = M, ncol = l)

    ## 2.2 Given J_{-1}, time1 <- proc.time()
    S <- foreach(m = 1: M, .combine = 'rbind' ) %dopar% {
        i <- 1
        J2 <- J
        J2[i] <- 1
        p21 <- supp(J2, dd)

        J2[l] <- 1
        p22 <- supp(J2, dd)

        if (p21 == 0) {
            p1 <- 0  # 1
        }else {
            p1 <- exp(ig * p22 * p22/p21)  # support * confidence
        }

        J3 <- J
        J3[i] <- 0
        p31 <- supp(J3, dd)

        J3[l] <- 1
        p32 <- supp(J3, dd)

        if (p31 == 0) {
            p2 <- 0
        } else {
            p2 <- exp(ig * p32 * p32/p31)
        }

        if (p1 == 0 & p2 == 0)
            p <- 0 else p <- p1/(p1 + p2)

        J[i] <- rbinom(1, 1, p)

        for (i in 2:(l - 1)) {
            if (J[i] == 1 & J[i - 1] == 1) {
                J3 <- J
                J3[i] <- 0
                p31 <- supp(J3, dd)

                J3[l] <- 1
                p32 <- supp(J3, dd)

                if (p31 == 0) {
                  p2 <- 0
                } else {
                  p2 <- exp(ig * p32 * p32/p31)
                }
            } else if (J[i] == 0 & J[i - 1] == 0) {
                J2 <- J
                J2[i] <- 1
                p21 <- supp(J2, dd)

                J2[l] <- 1
                p22 <- supp(J2, dd)

                if (p21 == 0) {
                  p1 <- 0
                } else {
                  p1 <- exp(ig * p22 * p22/p21)
                }
            } else {
                J3 <- J
                J3[i] <- 0
                p31 <- supp(J3, dd)

                J3[l] <- 1
                p32 <- supp(J3, dd)

                if (p31 == 0) {
                  p2 <- 0
                } else {
                  p2 <- exp(ig * p32 * p32/p31)
                }

                J2 <- J
                J2[i] <- 1
                p21 <- supp(J2, dd)

                J2[l] <- 1
                p22 <- supp(J2, dd)

                if (p21 == 0) {
                  p1 <- 0  # 1
                } else {
                  p1 <- exp(ig * p22 * p22/p21)
                }
            }

            if (p1 == 0 & p2 == 0)
                p <- 0 else p <- p1/(p1 + p2)
            J[i] <- rbinom(1, 1, p)
        }
        # S[m, ] <- J
        J
    }
    
    # proc.time()-time1
    colnames(S) <- colnames(dd)
    S1 <- as(S, "itemMatrix")

    # 1) compute the marginal frequency for each item appeared in the random sample
    sf1 <- itemFrequency(S1)
    marginal_frequency <- sort(sf1[which(sf1 > 0)], decreasing = F)  # 67

    # 2) Compute the frequency of the sampled itemsets for lhs
    W1 <- unique(S1)
    W <- as(W1, "matrix")

    coun <- rep(0, dim(W)[1])

    for (i in 1:dim(S)[1]) {
        for (j in 1:dim(W)[1]) {
            if (sum(abs(S[i, ] - W[j, ])) == 0)
                coun[j] <- coun[j] + 1
        }
    }

    frequency <- coun/M

    # 3) compute the C*S for unique sampled association rules
    importance <- NULL
    support <- NULL
    confidence <- NULL

    for (i in 1:dim(W)[1]) {
        rw <- W[i, ]
        p21 <- supp(rw, dd)

        rw[l] <- 1
        p22 <- supp(rw, dd)

        if (p21 == 0) {
            p1 <- 0  # 1
            p2 <- 0
        } else {
            p1 <- p22 * p22/p21
            p2 <- p22/p21
        }

        support <- append(support, p22)
        confidence <- append(confidence, p2)
        importance <- append(importance, p1)
    }

    Ms <- matrix(0, nrow = dim(W)[1], ncol = 4)
    Ms[, 1] <- support
    Ms[, 2] <- confidence
    Ms[, 3] <- importance
    Ms[, 4] <- frequency

    WW <- cbind(Ms, W)
    WW1 <- WW[order(WW[, 3], decreasing = TRUE), ]
    WW2 <- as(WW1[, -c(1:4)], "transactions")

    Me <- WW1[, c(1:4)]
    colnames(Me) <- c("support", "confidence", "importance", "frequencies")

    # W2 <- inspect(WW1)
    list(sampled_items = marginal_frequency, sampled_rules = WW2, measures = Me)
}

