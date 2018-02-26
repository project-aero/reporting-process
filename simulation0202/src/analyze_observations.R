#!/usr/bin/Rscript

library(foreach)

doParallel::registerDoParallel(cores=5)
RNGkind("L'Ecuyer-CMRG")
set.seed(1234)
parallel::mc.reset.stream()

load("checkpoint-02.rda")

levs <- list()
levs$bandwidth <- c(36, 156)
levs$lag <- 1
analysis_des_mat <- do.call(expand.grid, levs)

calc_auc <- function(predictions, is_null){
  test <- !is.na(predictions)
  predictions <- predictions[test]
  is_null <- is_null[test]
  r <- rank(predictions)
  r1 <- sum(r[!is_null])
  n1 <- sum(!is_null)
  n2 <- sum(is_null)
  (r1 - n1 * (n1 + 1) / 2) / (n1 * n2)
}

analyze_observations <- function(bandwidth, lag, obs){
  get_tau <- function(stat_ts, time){
    stats::cor(stat_ts, time, use="pairwise.complete.obs", method="kendall")
  }
  get_ews_cor <- function(reports, bw=bandwidth, l=lag){
      stats <- spaero::get_stats(reports,
                                 center_bandwidth=bw, stat_bandwidth=bw,
                                 center_trend="local_constant",
                                 stat_trend="local_constant",
                                 center_kernel="uniform",
                                 stat_kernel="uniform", lag=l)$stats
      taus <- sapply(stats, get_tau, time=seq_along(reports))
      list(stats=stats, taus=taus)
  }
  loop_manip <- function(x){
    loop_omult <- function(x) lapply(x, get_ews_cor)
    sc <- lapply(x, loop_omult)
    taus <- lapply(sc, function(x) sapply(x, "[[", "taus"))
    taus <- do.call(cbind, taus)
    list(sc=sc, taus=taus)
  }
  sct <- lapply(obs, loop_manip)
  is_null <- c(rep(FALSE, ncol(sct$test$taus)), rep(TRUE, ncol(sct$null$taus)))
  preds <- cbind(sct$test$taus, sct$null$taus)

  df <- data.frame(is_null, t(preds), row.names=1:ncol(preds))
  samp_stat <- function(x, w) {
    apply(x[w, -1], 2, calc_auc, is_null=df$is_null[w])
  }
  bs <- boot::boot(data=df, statistic=samp_stat, R=300)
  bssd <- apply(bs$t, 2, sd, na.rm=TRUE)
  list(sct=sct, auc=bs$t0, auc_stderr=bssd)
}

analyzed_observations <-
  foreach (i=seq(1, nrow(process_des_mat)),
           .options.multicore=list(set.seed=TRUE, preschedule=FALSE)) %:%
    foreach (j=seq(1, nrow(observation_des_mat)),
             .options.multicore=list(set.seed=TRUE, preschedule=FALSE)) %:%
      foreach (m=seq(1, nrow(analysis_des_mat)),
               .options.multicore=list(set.seed=TRUE, preschedule=FALSE)) %dopar% {
                 do.call(analyze_observations,
                         c(analysis_des_mat[m, ],
                           list(obs=simulated_observations[[i]][[j]])))
    }
warnings()

res <- list()
n <- 1
for (i in seq(1, nrow(process_des_mat))){
  pvars <- process_des_mat[i, ]
  for (j in seq(1, nrow(observation_des_mat))){
    ovars <- observation_des_mat[j, ]
    for (m in seq(1, nrow(analysis_des_mat))){
      auc <- analyzed_observations[[i]][[j]][[m]]$auc
      names(auc) <- paste("AUC", names(auc), sep="_")
      aucse <- analyzed_observations[[i]][[j]][[m]]$auc_stderr
      names(aucse) <- paste(names(auc), "stderr", sep="_")
      avars <- analysis_des_mat[m, ]
      res[[n]] <- data.frame(c(pvars, ovars, avars, auc, aucse))
      n <- n + 1
    }
  }
}
res <- do.call(rbind, res)

save.image(file="checkpoint-03.rda")
