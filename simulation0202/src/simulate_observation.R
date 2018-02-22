#!/usr/bin/Rscript

library(foreach)

doParallel::registerDoParallel(cores=30)
RNGkind("L'Ecuyer-CMRG")
set.seed(3222)
parallel::mc.reset.stream()

load("checkpoint-01.rda")

levs <- list()
levs$aggregation_days <- c(7, 30)
levs$neg_bin_k <- c(0.01, 1, 100)
levs$observation_multiplier <- 1
levs$reporting_prob <- 2^(seq(-8, 0, len=21))
observation_des_mat <- do.call(expand.grid, levs)

sample_observation <- function(aggregation_days, neg_bin_k,
                               observation_multiplier, reporting_prob,
                               procs){
    splt <- function(x) {
        split(x, x$sim)
    }
    aggr <- function(x, p=reporting_prob, m=observation_multiplier,
                     k=neg_bin_k) {
        tots <- aggregate.ts(x$cases, nfrequency=1 / aggregation_days)
        mu <- tots * p
        n <- length(mu)
        rnbinom(n=n, mu=mu, size=k)
    }
    transform <- function(x, n=observation_multiplier){
        sims <- splt(x)
        replicate(n=n, lapply(sims, aggr), simplify=FALSE)
    }
    lapply(procs, transform)
}

simulated_observations <-
  foreach (i=seq(1, nrow(process_des_mat)),
           .options.multicore=list(set.seed=TRUE)) %:%
    foreach (j=seq(1, nrow(observation_des_mat))) %dopar% {
      do.call(sample_observation,
              c(observation_des_mat[j, ], list(simulated_procs[[i]])))
    }

save.image(file="checkpoint-02.rda")
