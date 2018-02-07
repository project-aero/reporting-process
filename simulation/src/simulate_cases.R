#!/usr/bin/Rscript

library(foreach)

doParallel::registerDoParallel(cores=30)
RNGkind("L'Ecuyer-CMRG")
set.seed(3234)
parallel::mc.reset.stream()

levs <- list()
levs$external_forcing <- c(1 / 7)
levs$host_lifetime <- c(70 * 365)


levs$infectious_days <- c(30, 7)
levs$observation_days <- c(10) * 365
levs$population_size <- 10^c(6)
levs$process_reps <- 1000

levs$scenario <- c("emergence")
process_des_mat <- do.call(expand.grid, levs)

EndemicEquilSIR <- function(beta=(R0 * (mu + gamma)), eta=17/5e4,
                            gamma=365/22, mu=1/50, p=0,  R0=17, verbose=FALSE) {
      ## Computes the endemic equilibrium of an SIR model with immigration
      ##
      ## Args:
      ##   beta: numeric. The transmission rate.
      ##   eta: numeric. The rate of infection from outside.
      ##   gamma: numeric. The recovery rate.
      ##   mu: numeric. The birth rate.
      ##   p: numeric. The vaccination uptake.
      ##   R0: numeric. The basic reproduction number.
      ##
      ## Returns:
      ##   A list with numeric elements S, I, and R, coresponding to the
      ##   equilibrium fractions of the population in the
      ##   susceptible, infected, and removed states.
      stopifnot(c(beta, eta, gamma, p, R0) >= 0)
      stopifnot(p <= 1)
      a <- - beta * (gamma + mu)
      b <- beta * mu * (1 - p) - (gamma + mu) * (eta + mu)
      c <- mu * (1 - p) * eta
      eq <- (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
      i.star <- ifelse(p == 1, 0, eq)
      s.star <- ifelse(p == 1, 0, mu * (1 - p)/ (beta * i.star + eta + mu))
      if (verbose) {
        ds.star <- mu *(1 - p) - beta * s.star * i.star - eta * s.star - mu * s.star
        di.star <- beta * s.star * i.star + eta * s.star - (gamma + mu) * i.star
        cat('dS = ', ds.star, '\n')
        cat('dI = ', di.star, '\n')
      }
      return(list(S=s.star, I=i.star, R=1 - i.star - s.star))
}

sample_process <- function(external_forcing, host_lifetime, infectious_days,
                           observation_days, population_size, process_reps,
                           scenario){
    stopifnot(scenario == "emergence")

    times <- seq(0, observation_days)

    params <- c(gamma=1 / infectious_days, mu=1 / host_lifetime,
                d=1 / host_lifetime, eta=external_forcing / population_size,
                beta=0, rho=0.1, S_0=1, I_0=0, R_0=0, N_0=population_size)

    beta_critical <- (params["gamma"] + params["d"]) / population_size
    initial_fraction_critical <- 0.5
    equil <- EndemicEquilSIR(beta = beta_critical * initial_fraction_critical,
                             eta = params["eta"], gamma = params["gamma"],
                             mu = params["mu"], p = 0)
    params["S_0"] <- equil[["S"]]
    params["I_0"] <- equil[["I"]]
    params["R_0"] <- equil[["R"]]

    covnul <- data.frame(gamma_t=c(0, 0), mu_t=c(0, 0), d_t=c(0, 0),
                         eta_t=c(0, 0), beta_t=rep(beta_critical, 2) * 0.5,
                         time=c(0, observation_days))
    covalt <- data.frame(gamma_t=c(0, 0), mu_t=c(0, 0), d_t=c(0, 0),
                         eta_t=c(0, 0),
                         beta_t=c(beta_critical * 0.5, beta_critical),
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
