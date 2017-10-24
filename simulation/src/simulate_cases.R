#!/usr/bin/Rscript

library(foreach)

doParallel::registerDoParallel(cores=30)
RNGkind("L'Ecuyer-CMRG")
set.seed(3234)
parallel::mc.reset.stream()

levs <- list()
levs$external_forcing <- c(1 / 7, 1 / 30, 1 / 365)
levs$host_lifetime <- c(70 * 365)


levs$infectious_days <- c(30, 7)
levs$observation_days <- c(20) * 365
levs$population_size <- 10^c(6)
levs$process_reps <- 1000

levs$scenario <- c("emergence")
process_des_mat <- do.call(expand.grid, levs)

sample_process <- function(external_forcing, host_lifetime, infectious_days,
                           observation_days, population_size, process_reps,
                           scenario){
    stopifnot(scenario == "emergence")

    times <- seq(0, observation_days)
    params <- c(gamma=1 / infectious_days, mu=1 / host_lifetime,
                d=1 / host_lifetime, eta=external_forcing / population_size,
                beta=0, rho=0.1, S_0=1, I_0=0, R_0=0, N_0=population_size)
    covnul <- data.frame(gamma_t=c(0, 0), mu_t=c(0, 0), d_t=c(0, 0),
                         eta_t=c(0, 0), beta_t=c(0, 0),
                         time=c(0, observation_days))
    beta_final <- (params["gamma"] + params["d"]) / population_size
    covalt <- data.frame(gamma_t=c(0, 0), mu_t=c(0, 0), d_t=c(0, 0),
                         eta_t=c(0, 0), beta_t=c(0, beta_final),
                         time=c(0, observation_days))

    simtest <- spaero::create_simulator(times=times, params=params, covar=covalt)
    simnull <- spaero::create_simulator(times=times, params=params, covar=covnul)

    ret <- list()
    do_sim <- function(obj, nsim=process_reps){
        cols_to_delete <- c("reports", "gamma_t", "mu_t", "d_t", "eta_t")
        ret <- pomp::simulate(obj, nsim=nsim, as.data.frame=TRUE)
        ret[, !colnames(ret) %in% cols_to_delete]
    }
    ret$test <- do_sim(simtest)
    ret$null <- do_sim(simnull)
    ret
}

simulated_procs <- foreach(i=seq(1, nrow(process_des_mat)),
                           .options.multicore=list(set.seed=TRUE)) %dopar%
  do.call(sample_process, process_des_mat[i, ])

save.image(file="checkpoint-01.rda")

## Next we'll output a times series of cases for one of the figures in
## the paper. We'll delete observations that do not occur at multiples
## of ``snapshots`` to save space because they are not used in any
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
