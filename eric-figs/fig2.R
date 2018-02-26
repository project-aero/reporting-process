# This script assumes the working directory is set to "eric-figs/"

# Dependencies and helper functions in header.R

source("./R/header.R")

# Restore random number generator state from file
load("./data/RNGdata20170602_163629.Rdata")
RNGversion(RNGversion) # enforces RNG version
do.call("RNGkind",as.list(RNGkind)) # enforces RNG kind
.Random.seed <- RNGseed # restore seed to global environment

# Data

load("./output/data/sim.7di.Rda")
load("./output/data/sim.7di.I.Rda")
load("./output/data/reports.7di.Rda")
load("./output/data/sim_params.Rda")
load("./output/data/analysis_params.Rda")

## Subsetted Data

reports.7di.30dr.perfect <- reports.7di[which(
  reports.7di$tau==30 & 
    reports.7di$bandwidth==35 & 
    reports.7di$rho==1 &
    is.na(reports.7di$disp)
), ]$data[[1]]

reports.7di.7dr.perfect <- reports.7di[which(
  reports.7di$tau==7 & 
    reports.7di$bandwidth==35 & 
    reports.7di$rho==1 &
    is.na(reports.7di$disp)
), ]$data[[1]]

## Data must be coerced to a time series first, 
## setting start time to the minimum time in dataframe (0 in this case), 
## otherwise aggregate.ts() will coerce the data to a time series starting at time 1.  
sim.7di.cases <- ts(sim.7di$cases, start = sim.7di$time[1])  

# extract caselist (under perfect reporting) (I --> R)
caselist <- sim.7di[rep(seq(nrow(sim.7di)), sim.7di$cases), c("time", "cases")]
caselist$cases <- caselist$cases != 0

# extract list of transmissions (S --> I)
transmissionlist <- sim.7di[rep(seq(nrow(sim.7di)), sim.7di$transmissions), c("time", "transmissions")]
transmissionlist$transmissions <- transmissionlist$transmissions != 0
rownames(transmissionlist) <- NULL
transmissionlist$ID <- as.integer(rownames(transmissionlist))

# extract list of deaths of I
deathlist <- sim.7di[rep(seq(nrow(sim.7di)), sim.7di$deathsI), c("time", "deathsI")]
deathlist$deathsI <- deathlist$deathsI != 0

# merged list of deaths of I (I --> death) and "cases" (I --> R)
caseanddeathlist <- merge(caselist,deathlist,all=TRUE) # renumbered

# Scales & Ranges

# Set plotting window (days)
plotxmin <- 7*870 # days (Week 870)
plotxlength <- 7*12 # days (12 weeks)
plotxmax <- plotxmin+plotxlength # days
windows.7day <- seq(from = 0, to = sim_params$observation_days, by = 7) # aggregation period boundaries
windows.30day <- seq(from = 0, to = sim_params$observation_days, by = 30) # aggregation period boundaries

# Top panel y limits
plotymin <- 0
plotymax <- max(
  window(reports.7di.30dr.perfect,start=plotxmin-30,end=plotxmax+30),
  window(reports.7di.7dr.perfect,start=plotxmin-30,end=plotxmax+30),
  window(sim.7di$I,start=plotxmin-30,end=plotxmax+30)
)

# Bottom panel Y limits
min.transmissionID <- as.integer(rownames(head(transmissionlist[transmissionlist$time >= plotxmin,],1)))
min.transmissionID <- floor(min.transmissionID/10)*10

max.transmissionID <- as.integer(rownames(tail(transmissionlist[transmissionlist$time <= plotxmax,],1)))
max.transmissionID <- ceiling(max.transmissionID/10)*10

# Top panel Y ticks
top.y.tick.interval <- 10
top.y.ticks <- top.y.tick.interval*seq(from = 0, to = ceiling(plotymax/top.y.tick.interval), by = 1)

## Visual Variables

color.7day <- unipalette['lightblue']
color.30day <- unipalette['red']
color.I <- unipalette['black']
color.grid <- rgb(0,0,0,.25)
lty.7day <- "dotted" # for vertical dividers
lty.30day <- "dashed" # for vertical dividers
lwd.7day <- 1.5
lwd.30day <- 3
lwd.I <- 1

## divider parameters as list
dividers.7day <- list(v=windows.7day, col=color.7day, lty=lty.7day, lwd=lwd.7day, xpd=FALSE)
dividers.30day <- list(v=windows.30day, col=color.30day, lty=lty.30day, lwd=lwd.30day, xpd=FALSE)

# Typography

font.family <- "Times" # PLOS publications require either Times, Arial or Symbol.  only Times is included in R
font.sizes <- seq(from = 8, # publisher's minimum point size (points)
                  to = 12, # publisher's maximum point size (points) 
                  length.out = 5)
font.size.normal <- mean(font.sizes)
font.scales <- font.sizes/mean(font.sizes)
names(font.scales) <- names(font.sizes) <- c("XS", "S", "M", "L", "XL")

# Figure dimensions

figure.widths <- c(min=2.63, page=7.5, column=5.2) # in inches, as defined by publisher
figure.heights <- c(min=1, page=8.75) # in inches, as defined by publisher

# Margins and Figure Bounds

margins = c(4,5,4,8)+0.1
bottom.pannel.height <- .75
bottom.pannel <- c(0,1,0,bottom.pannel.height) # panel bounds: x0,x1,y0,y1 as fraction of figure region
top.pannel <- c(0,1,bottom.pannel.height-.2,1) # panel bounds: x0,x1,y0,y1 as fraction of figure region

# PDF output

pdf(
  file = "./output/plots/fig2.pdf",
  title = "Figure 2", # displayed in title bar of PDF Reader
  width = figure.widths['column'], # inches.  Must fit publisher's min and max figure dimensions
  height = figure.heights['page']*.7, # inches.  Must fit publisher's min and max figure dimensions
  family = font.family, 
  pointsize = font.size.normal # default size of text (points).
)

# init figure

par(mar=margins)
par(lend="butt")

# Bottom Panel (plot of individual infections as multiple line segments)

## init plot
par(fig=bottom.pannel)
plot(0,0, type='n', axes=FALSE, ann=FALSE, yaxt="n",
     xlim=c(plotxmin,plotxmax), ylim=c(min.transmissionID,max.transmissionID))

## plot cases with actual transmission times (onsets), random pairing
## Each case time is paired with a random transmission time <= case time.
## For some cases, case time == transmission time
tmp_transmissionlist <- transmissionlist

## Plot segments
for(i in 1:length(transmissionlist$time)){  
  
  offsettime <- caseanddeathlist$time[i]
  onsetindex <- sample(1:nrow(tmp_transmissionlist[tmp_transmissionlist$time<=offsettime,]),1)
  id <- tmp_transmissionlist$ID[onsetindex]
  xstart <- tmp_transmissionlist$time[onsetindex]
  # recovery, death or end of obs period:
  xend <- ifelse(!is.na(caseanddeathlist$time[i]),caseanddeathlist$time[i],sim_params$observation_days) 
  
  # draw segment
  segments(
    x0 = xstart, y0 = id,
    x1 = xend, y1 = id
  )
  symbols(x = xend, y = id, add=TRUE, inches = FALSE, circles = .3, bg='black')
  
  tmp_transmissionlist <- tmp_transmissionlist[-onsetindex, , drop = FALSE]
}

## plot observation window dividers
do.call(abline, dividers.7day)
do.call(abline, dividers.30day)

# Axes
axis(side = 1, # 1 specifies bottom axis - minor ticks
     at = seq(plotxmin-7,plotxmax+7,7), # tick locations (in data units).
     labels = FALSE, # vector of tick labels
     las = 1, # label orientation (0:parall., 1:horiz., 2:perp., 3:vert.)
     tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
     line = .5, # number of lines into margin at which axis is drawn 
     lty = "solid",
     lwd = 0, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
)
tmp.at <- seq(plotxmin,plotxmax+7,28)
tmp.lab <- (tmp.at-plotxmin)/7
axis(side = 1, # 1 specifies bottom axis - major ticks
     at = tmp.at, # tick locations (in data units).
     labels = tmp.lab,
     tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
     line = .5, # number of lines into margin at which axis is drawn 
     lty = "solid",
     lwd = 0, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
)
tmp.at <- seq(min.transmissionID,max.transmissionID,20)
tmp.lab <- tmp.at - min.transmissionID
axis(side = 2, # 1 specifies left axis
     at = tmp.at, # tick locations (in data units).
     labels = tmp.lab, # vector of tick labels
     las = 1, # label orientation (0:parall., 1:horiz., 2:perp., 3:vert.)
     tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
     line = .5, # number of lines into margin at which axis is drawn 
     lty = "solid",
     lwd = 0, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
)
title(xlab = "Week", line = 2.5)
title(ylab = "Case ID", line = 3)

# Top Panel (initialize plot of time series)

## init plot
par(fig=top.pannel, new=TRUE)

## plot time series of reports (7 day)
plot(reports.7di.7dr.perfect, type='s', axes = FALSE, xlab="", bty='n', ylab="", 
     lwd=lwd.7day, 
     col=color.7day,
     xlim=c(plotxmin,plotxmax),  
     ylim=c(plotymin,plotymax) # max 30day reports in window 
)  

## plot time series of reports (30 day)
lines(reports.7di.30dr.perfect, type='s', lwd=lwd.30day, col=color.30day)

## plot time series of true cases
lines(sim.7di$I, type='S', lwd=lwd.I, col=color.I) # note: type="S" not "s"

## grid lines
abline(
  h = top.y.ticks,
  col = color.grid, 
  lty = 'solid',
  lwd = 1
)

## plot observation window dividers
do.call(abline, dividers.7day)
do.call(abline, dividers.30day)

## Axes
axis(side = 1, # 1 specifies bottom axis - minor ticks
     at = seq(plotxmin-7,plotxmax+7,7), # tick locations (in data units).
     labels = FALSE, # vector of tick labels
     las = 1, # label orientation (0:parall., 1:horiz., 2:perp., 3:vert.)
     tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
     line = .5, # number of lines into margin at which axis is drawn 
     lty = "solid",
     lwd = 0, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
)
axis(side = 3, # 1 specifies bottom axis - minor ticks
     at = seq(plotxmin-7,plotxmax+7,7), # tick locations (in data units).
     labels = FALSE, # vector of tick labels
     las = 1, # label orientation (0:parall., 1:horiz., 2:perp., 3:vert.)
     tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
     line = .5, # number of lines into margin at which axis is drawn 
     lty = "solid",
     lwd = 0, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
)
mtext(text="Week", side=3, line = 2)

tmp.at <- seq(plotxmin,plotxmax+7,28)
tmp.lab <- (tmp.at-plotxmin)/7
axis(side = 3, # 1 specifies bottom axis - major ticks
     at = tmp.at, # tick locations (in data units).
     labels = tmp.lab,
     tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
     line = .5, # number of lines into margin at which axis is drawn 
     lty = "solid",
     lwd = 0, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
)
axis(side = 2, # 1 specifies left axis
     at = top.y.ticks, # tick locations (in data units).
     labels = top.y.ticks, # vector of tick labels
     las = 1, # label orientation (0:parall., 1:horiz., 2:perp., 3:vert.)
     tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
     line = .5, # number of lines into margin at which axis is drawn 
     lty = "solid",
     lwd = 0, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
)
title(ylab = "Number", line = 3)

## annotations

## Legend
par(lend="butt")
legend("topleft", xpd=NA, inset=c(1.01,0), xjust=0, yjust=0, cex=font.scales["XS"], y.intersp = 1.5,
       legend=c("Daily Snapshots\nof Number Infected", "Weekly Reports","Monthly Reports"),
       col=c(color.I,color.7day,color.30day),
       lwd=c(lwd.I,lwd.7day,lwd.30day),
       lty=1,
       bty='n')

# Close PDF
invisible(dev.off())