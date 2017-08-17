usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  library(p, character.only = TRUE)
}

sample_process <- function(external_forcing = 1 / 7, host_lifetime = 70 * 365,
                           infectious_days = 7, observation_days = 20 * 365,
                           sampling_interval = 1, population_size = 1e6,
                           process_reps = 1){
  ## Runs model and returns sample of state variables and flow from I to R, which is the cases column
  ## Default parameters are those used in the simulation study.
  
  times <- seq(0, observation_days, by = sampling_interval)
  params <- c(gamma=1 / infectious_days, mu=1 / host_lifetime,
              d=1 / host_lifetime, eta=external_forcing / population_size,
              beta=0, rho=0.1, S_0=1, I_0=0, R_0=0, N_0=population_size)
  beta_final <- (params["gamma"] + params["d"]) / population_size
  covalt <- data.frame(gamma_t=c(0, 0), mu_t=c(0, 0), d_t=c(0, 0),
                       eta_t=c(0, 0), beta_t=c(0, beta_final),
                       time=c(0, observation_days))
  
  simtest <- spaero::create_simulator(times=times, params=params, covar=covalt)
  
  ret <- list()
  do_sim <- function(obj, nsim=process_reps){
    cols_to_delete <- c("reports", "gamma_t", "mu_t", "d_t", "eta_t")
    ret <- pomp::simulate(obj, nsim=nsim, as.data.frame=TRUE)
    ret[, !colnames(ret) %in% cols_to_delete]
  }
  do_sim(simtest)
}

# new sample_observation() function taking a time series as input
sample_observation <- function (cases, sampling_interval = 1, tau = 1, reporting_prob = 1, dispersion_parameter = 100) {
  tots <- aggregate.ts(cases, nfrequency = 1/(tau/sampling_interval))
  mu <- tots * reporting_prob
  n <- length(mu)
  sampled <- rnbinom(n = n, mu = mu, size = dispersion_parameter)
  ts(sampled,start=start(tots),end=end(tots),frequency=frequency(tots))
}

# sample_observation <- function(cases, reporting_prob, dispersion_parameter, aggregation_days){
#   ## Aggregates time series of cases and adds random reporting error.
#   ## Note, a high values of the dispersion parameters _reduces_ the dispersion
#   ts(cases, start = toe)
#   tots <- aggregate.ts(cases, nfrequency=1 / aggregation_days)
#   mu <- tots * reporting_prob
#   n <- length(mu)
#   rnbinom(n=n, mu=mu, size=dispersion_parameter)
# }