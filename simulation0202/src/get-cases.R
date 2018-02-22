#!/usr/bin/Rscript

load("checkpoint-01.rda")

## We'll delete observations that do not occur at multiples of
## ``snapshots`` to save space because they are not used in any
## further analysis

test_snap <- function(x, period){
  mod <- x %% period
  abs(mod) < .Machine$double.eps^.5
}

get_snaps <- function(df) {
  tms <- df$time
  snapshots <- c(7, 30)
  tests <- lapply(snapshots, function(per) test_snap(tms, period = per))
  is_snap <- colSums(do.call(rbind, tests)) > 0
  df[is_snap, ]
}

for(i in seq_len(nrow(process_des_mat))){
  simulated_procs[[i]]$test <- get_snaps(simulated_procs[[i]]$test)
  simulated_procs[[i]]$null <- get_snaps(simulated_procs[[i]]$null)
}

save(simulated_procs, process_des_mat, file="cases.RData")
