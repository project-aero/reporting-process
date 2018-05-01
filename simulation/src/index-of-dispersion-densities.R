load("checkpoint-03.rda")


obs_des_levs <- c(85, 87, 89)
observation_des_mat[obs_des_levs, ]

sapply(obs_des_levs, function(x) {
    analyzed_observations[[1]][[x]][[1]]$auc["index_of_dispersion"]
})


tmpf <- function(y) {
    tmpff <- function(x) {
        analyzed_observations[[1]][[y]][[1]]$sct$test$sc[[1]][[x]]$stats$index_of_dispersion
    }
    sapply(1:100, tmpff)
}
idts <- lapply(obs_des_levs, tmpf)

get_quant <- function(x) apply(x, 1, quantile, na.rm = TRUE)

test_quants <- lapply(idts, get_quant)

tmpf <- function(y) {
    tmpff <- function(x) {
        analyzed_observations[[1]][[y]][[1]]$sct$null$sc[[1]][[x]]$stats$index_of_dispersion
    }
    sapply(1:100, tmpff)
}
idts <- lapply(obs_des_levs, tmpf)

null_quants <- lapply(idts, get_quant)

tmpf <- function(x, y, z) {
    matplot(y=t(x), type = 'l', col=1, xlab = "Index of window",
            ylab = "Quantiles of Index of Dispersion",
            main = paste("k =", z), lty = 1,
            sub = "aggregation_days = 7, reporting_prob = 0.19, infectious_days = 30, bandwidth=35" )
    legend("topleft", legend=c("test", "null"), col=1:2, lty = 1)
    matplot(y=t(y), type = 'l', col=2, add = TRUE, lty = 1)
   
}

pdf("quantiles-iod.pdf")
mapply(tmpf, test_quants, null_quants, observation_des_mat$neg_bin_k[obs_des_levs])
dev.off()




