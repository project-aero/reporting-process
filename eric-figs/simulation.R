# This script assumes the working directory is set to "eric-figs/"

# Dependencies and helper functions in header.R

source("./R/header.R")

# Random Number Generator state recall in a robust, portable way, for reproducibility

# Restore random number generator state from file
load("./data/RNGdata20170602_163629.Rdata")
RNGversion(RNGversion) # enforces RNG version
do.call("RNGkind",as.list(RNGkind)) # enforces RNG kind
.Random.seed <- RNGseed # restore seed to global environment

# Run Simulations

# simulation parameters
observation_days = 20 * 365
sampling_interval <- 1 # simulation sampling interval in days
population_size = 1e6
host_lifetime <- 70 * 365
infectious_days <- 7  # infectious period in days
external_forcing = 1 / 7
gamma <- 1/infectious_days
d <- 1/host_lifetime

sim_params <- list(
  observation_days = observation_days,
  sampling_interval = sampling_interval,
  population_size = population_size,
  host_lifetime = host_lifetime,
  infectious_days = infectious_days,
  external_forcing = external_forcing,
  gamma = gamma,
  d = d
)
save(sim_params,file="./output/data/sim_params.Rda")

# analysis parameters

analysis_params <- list(
  center_trend = "local_constant",
  stat_trend = "local_constant",
  center_kernel = "uniform",
  stat_kernel = "uniform" ,
  center_bandwidth = 100,
  stat_bandwidth = 100,
  lag = 1
)
save(analysis_params,file="./output/data/analysis_params.Rda")

# do single simulation

sim.7di <- sample_process(observation_days = observation_days,
                          sampling_interval = sampling_interval,
                          population_size = population_size,
                          host_lifetime = host_lifetime,
                          infectious_days  = infectious_days,
                          external_forcing = external_forcing)

save(sim.7di,file="./output/data/sim.7di.Rda")


# Stats for number of infected

## init nested dataframe (using tidyr definition of nested dataframe)
sim.7di.I <- list(interval=NA,lag=NA,bandwidth=NA,data=list(NULL),stats=list(NULL),taus=list(NULL))
attr(sim.7di.I, "row.names") <- 1L
class(sim.7di.I) <- c("tbl_df", "data.frame")

## resample and get stats, save results in nested dataframe
tmpfnc <- function(params) {
  i <- params['id']
  analysis_params$center_bandwidth <- analysis_params$stat_bandwidth <- params['bandwidth']
  analysis_params$lag <- params['lag']
  # Resample data, calc stats
  tmp.ts <- data[data$time %in% seq(0,max(data$time),params['resampling_interval']),]
  tmp.ts.zoo <- zoo(tmp.ts$I,tmp.ts$time)
  tmp.stats <- analysis(tmp.ts.zoo,analysis_params)
  tmp.stats.zoo <- zooreg(as.data.frame(tmp.stats$stats), start=0, deltat = params['resampling_interval'])
  tmp.taus <- as.data.frame(tmp.stats$taus)
  sim.7di.I[i,] <<- list(
    interval=params['resampling_interval'],
    lag=params['lag'],
    bandwidth=params['bandwidth'],
    data=list(tmp.ts.zoo),
    stats=list(tmp.stats.zoo),
    taus=list(tmp.taus)
  )
}

tmpgrid <- rbind(
  expand.grid(
    resampling_interval=1,
    lag=c(1,7,30),
    bandwidth=c(35,100)
  ),
  expand.grid(
    resampling_interval=c(7,30),
    lag=1,
    bandwidth=c(35,100)
  )
)
tmpgrid <- cbind(id=as.integer(rownames(tmpgrid)), tmpgrid)
data <- sim.7di # data used by tmpfunction

# do simulation
start.time <- proc.time()
cat("About to write to sim.7di.I (this may take a few minutes)...\n")
invisible(apply(tmpgrid, MARGIN=1, tmpfnc)) # get stats
save(sim.7di.I,file="./output/data/sim.7di.I.Rda")
cat("sim.7di.I saved to ./output/data/sim.7di.I.Rda\n","elapsed time =",as.numeric((proc.time()-start)[3]),"seconds")

# Reporting:

## init nested dataframe (using tidyr definition of nested dataframe)
reports.7di <- list(
  tau=NA,lag=NA,bandwidth=NA,rho=NA,disp=NA,
  data=list(NULL),stats=list(NULL),taus=list(NULL)
)
attr(reports.7di, "row.names") <- 1L
class(reports.7di) <- c("tbl_df", "data.frame")

## aggregate cases and get stats, save results in nested dataframe
tmpfnc <- function(params) {
  i <- params['id']
  analysis_params$center_bandwidth <- analysis_params$stat_bandwidth <- params['bandwidth']
  analysis_params$lag <- params['lag']
  # aggregate
  if(is.na(params['disp'])) {
    tmp.ts.zoo <- zoo(aggregate.ts(data, nfrequency=1 / (params['tau']/sampling_interval)))
  }else{
    tmp.ts.zoo <- zoo(
      sample_observation(data,
                         sampling_interval = sampling_interval,
                         tau = params['tau'],
                         reporting_prob = params['rho'],
                         dispersion_parameter = params['disp']
      )
    )
  }
  # calc stats
  tmp.stats <- analysis(tmp.ts.zoo,analysis_params)
  tmp.stats.zoo <- zooreg(as.data.frame(tmp.stats$stats), 
                          start=0, 
                          frequency = 1/(params['tau']/sampling_interval)
  )
  tmp.taus <- as.data.frame(tmp.stats$taus)
  
  # add stats to nested df
  i <- params['id']
  reports.7di[i,] <<- list(
    tau = params['tau'],
    lag = params['lag'],
    bandwidth = params['bandwidth'],
    rho = params['rho'],
    disp = params['disp'],
    data = list(tmp.ts.zoo),
    stats = list(tmp.stats.zoo),
    taus = list(tmp.taus)
  )
}

tmpgrid <- rbind(
  expand.grid(
    tau = c(7,30),
    lag=1,
    bandwidth=c(35,100),
    rho = 1,
    disp = NA
  ),
  expand.grid(
    tau = c(7,30),
    lag=1,
    bandwidth=c(35,100),
    rho = 2^seq(from = 0, to = -8, by = -0.4),
    disp = 10^seq(2, -2, by = -1)
  )
)
tmpgrid <- cbind(id=as.integer(rownames(tmpgrid)), tmpgrid)
data <- sim.7di.cases # data used by tmpfnc

# do imperfect reporting
start.time <- proc.time()
cat("About to write reports.7di (this may take a few minutes)... \n")
invisible(apply(tmpgrid, MARGIN=1, tmpfnc))
save(reports.7di,file="./output/data/reports.7di.Rda")
cat("reports.7di saved to ./output/data/reports.7di.Rda\n","elapsed time =",as.numeric((proc.time()-start)[3]),"seconds")

