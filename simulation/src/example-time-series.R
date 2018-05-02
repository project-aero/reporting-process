#!/usr/bin/Rscript

library(dplyr)
load("checkpoint-02.rda")

rho_levs <- c(0.04736614, 1)
k_levs <- c(.01, 100)

which(with(observation_des_mat, (near(neg_bin_k, k_levs[1]) | near(neg_bin_k, k_levs[2])) &
                                (near(reporting_prob, rho_levs[1]) | near(reporting_prob, rho_levs[2])) &
                                (near(aggregation_days, 7)))) -> inds

testts <- lapply(inds, function(x) simulated_observations[[2]][[x]]$test[[1]][[1]])
nullts <- lapply(inds, function(x) simulated_observations[[2]][[x]]$null[[1]][[1]])

day <- seq(1, 7300, by = 7)

plotf <- function(m, n, i) {
    main <- paste(c("phi", "reporting_prob"),
                  round(observation_des_mat[i,c("neg_bin_k", "reporting_prob")], digits=3),
                  sep = "=")
    sub <- paste("aggregation_days == 7, infectious_days == 7")
    matplot(x = day, y = cbind(m, n), type = 'h',
            ylab = "Case reports", xlab = "Day", lty = 1, main = main, sub = sub)
    legend("topleft", legend=c("test", "null"), col=1:2, lty = 1)
}

pdf("example-case-reports-ts-high-low.pdf")
mapply(plotf, testts, nullts, inds)
dev.off()

