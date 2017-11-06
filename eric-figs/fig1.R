# This script assumes the working directory is set to "eric-figs/"

# Dependencies and helper functions in header.R

source("./R/header.R")

# Restore random number generator state from file
load("./data/RNGdata20170602_163629.Rdata")
RNGversion(RNGversion) # enforces RNG version
do.call("RNGkind",as.list(RNGkind)) # enforces RNG kind
.Random.seed <- RNGseed # restore seed to global environment

# Data

load("./output/data/sim.7di.I.Rda")
load("./output/data/reports.7di.Rda")

## Subsetted data
data.7day.snap <- sim.7di.I[which(
  sim.7di.I$interval==7 & 
    sim.7di.I$bandwidth==35 & 
    sim.7di.I$lag==1
), ]
data.7day.imperfect <- reports.7di[which(
  reports.7di$tau==7 & 
    reports.7di$bandwidth==35 & 
    reports.7di$rho==2^-3.2 & 
    reports.7di$disp==100
), ]
data.30day.imperfect <- reports.7di[which(
  reports.7di$tau==30 & 
    reports.7di$bandwidth==35 & 
    reports.7di$rho==2^-3.2 & 
    reports.7di$disp==1
), ]

# Scales & Ranges

## Set plotting window
plotxmin <- 10*365  # days
plotxlength <-  10*365  # days
plotxmax <- plotxmin+plotxlength # days

plotymin <- 0
plotymax <- max(
  window(data.7day.snap$data[[1]],start=plotxmin,end=plotxmax),
  window(data.7day.imperfect$data[[1]],start=plotxmin,end=plotxmax),
  window(data.30day.imperfect$data[[1]],start=plotxmin,end=plotxmax)
)

# Y ticks
y.tick.interval <- 10 
y.ticks <- seq(from = plotymin, to = plotymax, by = y.tick.interval)

# Y ticks
top.y.ticks <- seq(-.5,1,.5)

# Visual Variables

## Colors - REVIEW
color.7day.snap <- unipalette['black']
color.7day.snap.fill <- rgb(0,0,0,.1)
color.7day.imperfect <- unipalette['pink']
color.30day.imperfect <- unipalette['red']
color.crit <- unipalette['green']

## Line types
lty.7day <- "dotted" # for vertical dividers
lty.30day <- "dashed" # for vertical dividers
lty.crit <- "dashed" # for verticals
lwd.7day <- .5
lwd.30day <- 1.5
lwd.I <- .5
lwd.crit <- 4

# Typography

font.family <- "Times"
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
top.pannel <- c(0,1,.45,1) # panel bounds: x0,x1,y0,y1 as fraction of figure region
bottom.pannel <- c(0,1,0,.65) # panel bounds: x0,x1,y0,y1 as fraction of figure region

# PDF output
pdf(
  file = "./output/plots/fig1.pdf",
  title = "Figure 1", # displayed in title bar of PDF Reader
  width = figure.widths['page'], # inches.  Must fit publisher's min and max figure dimensions
  height = figure.heights['page']*.7, # inches.  Must fit publisher's min and max figure dimensions
  family = font.family, 
  pointsize = font.size.normal # default size of text (points).
)

# init figure
par(mar=margins)
par(lend="butt")

# Top Panel

## Initialize plot of stats
par(fig=top.pannel)
plot(0,0, type='n', axes=FALSE, ann=FALSE, yaxt="n",
     xlim=c(plotxmin,plotxmax), ylim=range(top.y.ticks))

## grid
abline(h=top.y.ticks,col='grey')

## weekly snapshots of I, lag 1
lines(data.7day.snap$stats[[1]]$autocorrelation,
      col=color.7day.snap, lwd=lwd.I)

## weekly reports
lines(data.7day.imperfect$stats[[1]]$autocorrelation,
      col=color.7day.imperfect, lwd=lwd.7day)

## monthly reports
lines(data.30day.imperfect$stats[[1]]$autocorrelation,
      col=color.30day.imperfect, lwd=lwd.30day)  

## Annotation
title(ylab="Autocorrelation", line = 3)
text(x=19.75*365,
     y=c(.925,.3,-.1),
     adj=1,
     cex=font.scales['M'],
     col=c(color.7day.snap,
           color.7day.imperfect,
           color.30day.imperfect),
     labels=c(
       TeX(sprintf("$\\tau = $ %.2f", data.7day.snap$taus[[1]]$autocorrelation)),
       TeX(sprintf("$\\tau = $ %.2f", data.7day.imperfect$taus[[1]]$autocorrelation)),
       TeX(sprintf("$\\tau = $ %.2f", data.30day.imperfect$taus[[1]]$autocorrelation))
     )
)

## plot critical threshold
abline(
  v=20*365, col=color.crit, lty=lty.crit, lwd = lwd.crit
)

## Axes

axis(side = 1, # 1 specifies bottom axis - major ticks
     at = seq(plotxmin,plotxmax,365), # tick locations (in data units).
     labels = FALSE,
     tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
     line = 1, # axis position (in lines)
     lty = "solid",
     lwd = 0, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
)
axis(side = 3, # 1 specifies bottom axis - major ticks
     at = seq(plotxmin,plotxmax,365), # tick locations (in data units).
     labels = seq(plotxmin,plotxmax,365)/365,
     tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
     line = .5, # axis position (in lines)
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
     line = .5, # axis position, in lines
     lty = "solid",
     lwd = 0, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
)
mtext(text="Year", side=3, line = 2.5)

## Legend
legend("topleft", xpd=NA, inset=c(1.01,0), xjust=0, yjust=0, cex=font.scales['XS'], y.intersp = 3,
       legend=c(
         "Weekly Snapshots\nof Number Infected", 
         "Weekly Reports", 
         "Monthly Reports"
       ),
       col=c(
         # color.ac.I.lag1,
         color.7day.snap,
         color.7day.imperfect,
         color.30day.imperfect
       ),
       lwd=c(
         # lwd.I,
         lwd.I,
         lwd.7day,
         lwd.30day),
       lty=1,
       bty='n')

# Bottom Panel

## init plot
par(fig=bottom.pannel, new=TRUE)
plot(0,0, type='n', axes=FALSE, ann=FALSE, yaxt="n",
     xlim=c(plotxmin,plotxmax), ylim=c(plotymin,plotymax))

## plot time series of weekly snapshots of number infected
poly <- data.7day.snap$data[[1]]
x <- time(poly)
y <- as.numeric(poly)
poly.x <- c(head(x,1),x,tail(x,1))
# poly.y <- c(0,y,0)
poly.x <- c(rbind(x,x)) ## double each time
poly.y <- c(0,
            c(rbind(y[2:length(y)],y[2:length(y)])),
            0)
polygon(x=poly.x,y=poly.y, lwd=.1, border=color.7day.snap, col=color.7day.snap.fill)

## plot time series of cases, aggregated Weekly (imperfect reporting), filled staircase
poly <- data.7day.imperfect$data[[1]]
x <- time(poly) 
y <- as.numeric(poly)
poly.x <- c(rbind(x,x)) ## double each time
poly.y <- c(0,
            c(rbind(y[2:length(y)],y[2:length(y)])),
            0)
polygon(x=poly.x,y=poly.y, lwd=lwd.7day, border=color.7day.imperfect, col=color.7day.imperfect)

## plot time series of cases, aggregated Monthly (imperfect reporting), staircase
lines(data.30day.imperfect$data[[1]], 
      type='s', pch=20, lwd=lwd.30day, col=color.30day.imperfect)

## plot critical threshold
abline(
  v=20*365, col=color.crit, lty=lty.crit, lwd = lwd.crit
)

## Axes

axis(side = 1, # 1 specifies bottom axis - major ticks
     at = seq(plotxmin,plotxmax,365), # tick locations (in data units).
     labels = seq(plotxmin,plotxmax,365)/365,
     tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
     line = .5, # axis position, in lines
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
     line = .5, # axis position, in lines
     lty = "solid",
     lwd = 0, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
)

## grid lines
abline(
  h = y.ticks,
  col = rgb(0,0,0,.25), 
  lty = 'solid',
  lwd = 1
)

## annotations
title(ylab = "Number", line = 3)
title(xlab = "Year", line = 2.5)

text(365*20, plotymax, xpd=NA, cex=font.scales['M'], labels="epidemic threshold", adj=c(1,1.5), srt=90, col=color.crit)

## Legend
par(lend="butt")
legend("topleft", xpd=NA, inset=c(1.01,0), xjust=0, yjust=0, cex=font.scales['XS'], y.intersp = 1.5,
       legend=c(
         "Weekly Snapshots\nof Number Infected", 
         "Weekly Reports", 
         "Monthly Reports"
       ),
       col=c(color.7day.snap.fill, color.7day.imperfect, color.30day.imperfect),
       lty=1,
       lwd=c(9,9,lwd.30day),
       bty='n')

# close PDF
invisible(dev.off())