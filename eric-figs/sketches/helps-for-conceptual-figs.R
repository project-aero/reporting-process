## sample_process() and sample_observation() function definitions moved to ../R/helpers.R

## Example usage


## Get data frame of states of the model as a function of time and
## flows from I to R for the preceding time interval. Those flows are
## in the column named cases.

set.seed(1)
head(sample_process())
set.seed(1) ## The results should match if the seed is the same
head(sample_process())

## Use a finer sampling grid. This can slow down things considerably.
head(sample_process(sampling_interval = 0.1))


## Get a time series of reports based on passing a time series of
## cases through the observation model.

sim <- sample_process()
reports <- sample_observation(
  sim$cases, 
  reporting_prob = 0.1,
  dispersion_parameter = 100, 
  aggregation_days = 1  # multiple of sampling intervals
  )

## Get distribution of errors for a given number of cases

reports <- sample_observation(cases = rep(10, 1000), reporting_prob = 0.1,
                              dispersion_parameter = 100, aggregation_days = 7)
hist(reports)

## plot time series of I

plot(sim$I, type='l')

## plot time series of reports

barplot(reports)

## scenario #1: aggregation period = infection period

sampling_interval <- 1 # simulation sampling interval in days
tau <- 7  # reporting interval in days
infectious_days <- 7  # infectious period in days

sim.7dayinfectious <- sample_process(sampling_interval=sampling_interval, infectious_days = infectious_days)
reports.7day <- sample_observation(
  sim.7dayinfectious$cases, 
  aggregation_days = tau/sampling_interval,
  reporting_prob = 0.1,
  dispersion_parameter = 100 
)

plot(sim.7dayinfectious$I, type='l')
barplot(reports.7day)

## Summary of 4 senarios to plot:
  
# reporting_prob: .005 to 1.0
# dispersion_parameter: .001 to 100

  