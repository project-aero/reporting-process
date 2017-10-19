# This script assumes the working directory is set to "eric-figs/"

# Dependencies and helper functions in header.R

source("./R/header.R")

# Data

load("./output/data/sim.7di.Rda")
load("./output/data/sim.7di.I.Rda")
load("./output/data/reports.7di.Rda")
load("./output/data/sim_params.Rda")
load("./output/data/analysis_params.Rda")

# Extract Case Reports

reports.7di.7dr.perfect <- reports.7di[which(
  reports.7di$tau==7 & 
    reports.7di$bandwidth==35 & 
    reports.7di$rho==1 &
    is.na(reports.7di$disp)
), ]$data[[1]]

reports.7di.7dr.imperfect <- reports.7di[which(
  reports.7di$tau==7 & 
    reports.7di$bandwidth==35 & 
    reports.7di$rho==.25 &
    reports.7di$disp==0.1
), ]$data[[1]]

# Set plotting window

plotxmin <- 7*960 # days
plotxlength <-  7*30 # days

plotxmax <- plotxmin+plotxlength # days

plotymin <- 0
plotymax <- max(
  window(reports.7di.7dr.perfect,start=plotxmin,end=plotxmax),
  window(reports.7di.7dr.imperfect,start=plotxmin,end=plotxmax),
  window(sim.7di$I,start=plotxmin,end=plotxmax)
)

# Y ticks

y.tick.interval <- 10
y.ticks <- seq(from = plotymin, to = plotymax, by = y.tick.interval)

## Visual Variables

color.7day <- unipalette['lightblue']
color.7day.imperfect <- unipalette['red']
color.I <- unipalette['black']
color.grid <- rgb(0,0,0,.25)
lty.7day <- "dotted" # for vertical dividers
lty.30day <- "dashed" # for vertical dividers
lwd.7day <- 1.5
lwd.30day <- 3
lwd.I <- 1

## About Devices: Plot is made to the PDF device.  After the PDF device is closed, the figure is read in as an image and converted  to PNG.

# PDF output

path <- "./output/plots/fig3.pdf"
pdf(
  file = path,
  title = "Figure 3",
  width = 7.25, height=3
)

# plot

## init
par(mar=c(4,5,1,7.5)+0.1)
plot(0,0, type='n', axes=FALSE, ann=FALSE, yaxt="n",
     xlim=c(plotxmin,plotxmax), ylim=c(plotymin,plotymax))

## plot time series of number of cases aggregated weekly, with imperfect reporting
lines(reports.7di.7dr.imperfect, type='s', pch=20, lwd=lwd.7day, col=color.7day.imperfect)

## plot time series of true number of cases aggregated weekly
lines(reports.7di.7dr.perfect, type='s', pch=20, lwd=lwd.7day, col=color.7day)

## plot time series of number of infected
lines(sim.7di$I, type='S', pch=20, lwd=lwd.I, col=color.I) # note: type="S" not "s"

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
tmp.at <- seq(plotxmin,plotxmax+7,28)
tmp.lab <- (tmp.at-plotxmin)/7
axis(side = 1, # 1 specifies bottom axis - major ticks
     at = tmp.at, # tick locations (in data units).
     labels = tmp.lab,
     tcl = -1, # tick length (fraction of plot width, neg to draw outside plot)
     line = .5, # number of lines into margin at which axis is drawn 
     lty = "solid",
     lwd = 0, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
)

axis(side = 2, # 1 specifies left axis
     at = y.ticks, # tick locations (in data units).
     labels = y.ticks, # vector of tick labels
     las = 1, # label orientation (0:parall., 1:horiz., 2:perp., 3:vert.)
     tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
     line = .5, # number of lines into margin at which axis is drawn 
     lty = "solid",
     lwd = 0, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
)

## grid lines
abline(
  h = y.ticks,
  col = color.grid, 
  lty = 'solid',
  lwd = 1
)

## annotations
title(ylab = "Number", line = 3)
title(xlab="Week", line = 2.5)

## Legend
par(lend="butt")
legend("topleft", xpd=NA, inset=c(1.01,0), xjust=0, yjust=0, cex=.65, y.intersp = 3.0,
       legend=c("Daily Snapshots\nof Number Infected", 
                "Case Reports\n(Perfect Reporting)",
                "Case Reports\n(Imperfect Reporting)"),
       col=c(color.I,color.7day,color.7day.imperfect),
       lwd=c(lwd.I,lwd.7day,lwd.7day),
       lty=1,
       bty='n')

## Close PDF
invisible(dev.off())