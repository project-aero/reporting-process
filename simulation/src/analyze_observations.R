#!/usr/bin/Rscript

library(foreach)

doParallel::registerDoParallel(cores=5)
RNGkind("L'Ecuyer-CMRG")
set.seed(1234)
parallel::mc.reset.stream()



levs <- list()
levs$bandwidth <- c(35, 100)
levs$lag <- c(1)
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

parallel_setup <- function(is_small, seed = 1){
### Due to the large size results of parallelized jobs, the SNOW-style
### cluster backend seems better than multicore approach to
### parallelization. With the cluster approach, we can fork before any
### results are computed to create workers rather than over time as
### results are generated. The forked processes share memory with the
### parent. So as the results accumulate and the parent's memory use
### grows, the memory allocated to each forked process also
### grows. Even though the forked processes don't need to modify the
### results, and thus the memory actually used by them is not great,
### the allocation seems to trigger some mechanism that kills jobs on
### the Olympus cluster. It could be the Linux OOM killer.
  require(foreach)
  if (is_small) {
    clust <- parallel::makeForkCluster(nnodes = 2, outfile = "cluster-outfile")
    doParallel::registerDoParallel(cl = clust, outfile = "cluster-outfile")
  } else {
    clust <- parallel::makeForkCluster(nnodes = 23)
    doParallel::registerDoParallel(cl = clust)
  }
  parallel::clusterSetRNGStream(cl = clust, iseed = seed)
  clust
}

is_small_scale <- Sys.getenv("is_small_scale") == "TRUE"
cluster <- parallel_setup(is_small = is_small_scale)

load("checkpoint-02.rda")
options <- list(preschedule = FALSE) # To keep results <2GB limit

analyzed_observations <-
  foreach (i=seq(1, nrow(process_des_mat)),
           .options.snow = options)) %:%
    foreach (j=seq(1, nrow(observation_des_mat)),
             .options.snow = options) %:%
      foreach (m=seq(1, nrow(analysis_des_mat)),
               .options.snow = options) %dopar% {
                 do.call(analyze_observations,
                         c(analysis_des_mat[m, ],
                           list(obs=simulated_observations[[i]][[j]])))
    }
warnings()

parallel::stopCluster(cluster)

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

save(res, file="res.RData")
