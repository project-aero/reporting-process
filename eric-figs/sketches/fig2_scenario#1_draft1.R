# Random Number Generator state recall in a robust, portable way, for reproducibility

# Restore random number generator state from file
load("./data/RNGstate20170602_151509.Rdata")
RNGversion(RNGversion) # enforces RNG version
do.call("RNGkind",as.list(RNGkind)) # enforces RNG kind
.Random.seed <- RNGseed # restore seed

# parameters
observation_days = 20 * 365
sampling_interval <- 1 # simulation sampling interval in days
tau <- 7  # reporting interval in days
infectious_days <- 7  # infectious period in days
reporting_prob <- 1.0
dispersion_parameter <- 100

sim.7di <- sample_process(observation_days = observation_days,
                          sampling_interval=sampling_interval, 
                          infectious_days = infectious_days)

# reports.7di.7dr.hrp.ldisp <- sample_observation(sim.7di$cases, 
#                                                 aggregation_days = tau/sampling_interval,
#                                                 reporting_prob = reporting_prob,
#                                                 dispersion_parameter = dispersion_parameter
#                                                 )

reports.7di.7dr.perfect <- aggregate.ts(sim.7di$cases, 
                                        nfrequency=1 / (7/sampling_interval)
)
reports.7di.30dr.perfect <- aggregate.ts(sim.7di$cases, 
                                         nfrequency=1 / (30/sampling_interval)
)

# plot(sim.7dayinfectious$I, type='l')
# barplot(reports.7day)
# max(sim.7dayinfectious$I)
# max(reports.7day)

# ## plot cases (flow from I to R) vs time.
# plot(cases~time, data=sim.7dayinfectious, xlim=c(0,365), ylim=(c(1,max(sim.7dayinfectious$cases))))


# caselist <- sim.7di[sim.7di$cases > 0,][, c("time", "cases")]
# caselist.e <- caselist[rep(seq(nrow(caselist)), caselist$cases),]

caselist <- sim.7di[rep(seq(nrow(sim.7di)), sim.7di$cases), c("time", "cases")]
caselist$cases <- caselist$cases != 0

transmissionlist <- sim.7di[rep(seq(nrow(sim.7di)), sim.7di$transmissions), c("time", "transmissions")]
transmissionlist$transmissions <- transmissionlist$transmissions != 0
rownames(transmissionlist) <- NULL

deathlist <- sim.7di[rep(seq(nrow(sim.7di)), sim.7di$deathsI), c("time", "deathsI")]
deathlist$deathsI <- deathlist$deathsI != 0

caseanddeathlist <- merge(caselist,deathlist,all=TRUE) # renumbered



## random heights for bars

# heights <- rnorm(n=sum(sim.7dayinfectious$cases),mean=0,sd=1.5)
# heights <- runif(n=sum(sim.7dayinfectious$cases),min=-4.4,max=4.5)

height.levels <- 35

heights <- as.vector(
  replicate(
    ceiling(sum(transmissionlist$transmissions)/height.levels), # number of calls
    c(
      sample(1:height.levels, height.levels, replace=F)+runif(1,-.25,.25),
      sample(1:height.levels, height.levels, replace=F)+runif(1,-.25,.25)+.5
    )
  )
)[1:length(transmissionlist$time)]

# Plotting window (days)
plotxmin <- 7*52*15 # days
plotxlength <- 180 # days
plotxmax <- plotxmin+plotxlength # days

# Plot Cases as segments

par(fig=c(0,1,0,.6))

plot(0,0, type='n', axes=FALSE, ann=FALSE, yaxt="n",
     xlim=c(plotxmin,plotxmax), ylim=c(min(heights)-3,max(heights)+3))
# segments(
#   caselist$time,heights, # x0,y0
#   caselist$time-infectious_days,heights #x1,y1
#   )

# plot cases with fixed infection period
# for(i in 1:length(caselist$cases)){
#     heightvector <- heights[i]+runif(caselist$cases[i],-5,5) ## runif = KLUDGE
#     segments(
#         caselist$time[i]-infectious_days,heightvector, # x0,y0
#         caselist$time[i],heightvector #x1,y1
#     )
#     points(x=rep(caselist$time[i],caselist$cases[i]),y=heightvector, pch=20)
# }

# # Plot cases with actual transmission times (onsets), FIFO pairing
## Each case time is paired with the earliest available transmission time.
# for(i in 1:length(transmissionlist$time)){
#     # heightvector <- heights[i]+runif(caselist$cases[i],-5,5) ## runif = KLUDGE
#     segments(
#       transmissionlist$time[i],heights[i], # x0,y0
#       # x1 = recovery, death or end of obs period:
#       x1 = ifelse(!is.na(caseanddeathlist$time[i]),caseanddeathlist$time[i],observation_days), 
#       y1 = heights[i]
#       )
#     points(
#       x = ifelse(!is.na(caseanddeathlist$time[i]),caseanddeathlist$time[i],observation_days),
#       y = heights[i], 
#       pch = ifelse(caseanddeathlist$cases[i]==TRUE,20,NA)
#       )
# }

# plot cases with actual transmission times (onsets), random pairing
## Each case time is paired with a random transmission time <= case time.
## For some cases, case time == transmission time
tmp_transmissionlist <- transmissionlist
for(i in 1:length(transmissionlist$time)){  
  # heightvector <- heights[i]+runif(caselist$cases[i],-5,5) ## runif = KLUDGE
  
  offsettime <- caseanddeathlist$time[i]
  onsetindex <- sample(1:nrow(tmp_transmissionlist[tmp_transmissionlist$time<=offsettime,]),1)
  # onsettime <- sample(transmissionlist$time[1:"firstindexwhosetimeis less than time of R[i]"],1)
  
  segments(
    x0 = tmp_transmissionlist$time[onsetindex],
    y0 = heights[i],
    # x1 = recovery, death or end of obs period:
    x1 = ifelse(!is.na(caseanddeathlist$time[i]),caseanddeathlist$time[i],observation_days), 
    y1 = heights[i]
  )
  points(
    x = ifelse(!is.na(caseanddeathlist$time[i]),caseanddeathlist$time[i],observation_days),
    y = heights[i], 
    pch = ifelse(caseanddeathlist$cases[i]==TRUE,20,NA)
  )
  
  tmp_transmissionlist <- tmp_transmissionlist[-onsetindex, , drop = FALSE]    
  
}
axis(side = 1, # 1 specifies bottom axis
     at = seq(plotxmin,plotxmax,7), # vector of tick locations (in data units).
     labels = seq(plotxmin,plotxmax,7), # vector of tick labels
     las = 1, # label orientation (0:parall., 1:horiz., 2:perp., 3:vert.)
     tck = -.05, # tick length (fraction of plot width, neg to draw outside plot)
     pos = -1, # axis position, in units of y data
     lty = "solid",
     lwd = 1, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
     # col.ticks = NULL, # tick line color
)# abline(v=caselist$time) # case times
title(xlab="day")

# plot observation windows
abline(v=(7*0:(observation_days/7))-.1, col=rgb(1,0,0,.5))
abline(v=(30*0:(observation_days/30))-.1, col=rgb(0,0,1,.5), lty=2)


# plot time series

par(fig=c(0,1,.35,1), new=TRUE)

# plot time series of reports (7 day)

plot(reports.7di.7dr.perfect, type='s', lwd=3, xaxt="n", xlab="", ylab="cases", col=rgb(1,0,0,.5),
     xlim=c(plotxmin+1,plotxmax+1),  # KLUDGE: xlim=c(plotxmin+1,plotxmax+1) compensates for type='s'
     ylim=c(0,max(window(reports.7di.30dr.perfect,plotxmin-30,plotxmax+30)))  # max 30day reports in window
)  

# plot time series of reports (30 day)
lines(reports.7di.30dr.perfect, type='s', lwd=3, col=rgb(0,0,1,.5))

# plot time series of cases
lines(sim.7di$I, type='s')

