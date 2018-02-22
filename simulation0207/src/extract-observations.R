#!/usr/bin/Rscript

load("checkpoint-02.rda")
set.seed(1)

res <- list()
for (i in seq(1, nrow(process_des_mat))){
  pvars <- process_des_mat[i, ]
  jseq <- sample.int(nrow(observation_des_mat), 5) ## to reduce the size of the data
##  for (j in seq(1, nrow(observation_des_mat))){
  for (j in jseq){
    ovars <- observation_des_mat[j, ]
    sim_obs <- simulated_observations[[i]][[j]]
    stopifnot(length(sim_obs$null) == 1L) ## To simplify code, we
                                          ## assume that there is only
                                          ## one realization of the
                                          ## observation model per
                                          ## realization of the
                                          ## process model.
    stopifnot(length(sim_obs$test) == 1L)
    tdata <- sim_obs$test[[1L]]
    ndata <- sim_obs$null[[1L]]
    make_df <- function(dat, repid, is_emergence){
      t <- seq_along(dat)
      data.frame(c(pvars, ovars), process_replicate_id = repid,
                 is_emergence = is_emergence,
                 cbind(time_index = t, reported_cases = dat))
    }
    tdf <- Map(make_df, dat = tdata, repid = seq_along(tdata),
               is_emergence = TRUE)
    tdn <- Map(make_df, dat = ndata, repid = seq_along(ndata),
               is_emergence = FALSE)
    res <- c(res, tdf, tdn)
  }
}

reports_ts <- do.call(rbind, res)
## Delete potentially confusing column
reports_ts$scenario <- NULL
rownames(reports_ts) <- NULL

write.csv(reports_ts, file="reports_ts.csv")
save(reports_ts, file = "reports_ts.RData")
