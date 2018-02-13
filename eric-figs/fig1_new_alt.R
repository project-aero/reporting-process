# NOTES -------------------------------------------------------------------


# Level 2 data is missing values

# Use data.gov data with download script


# HEADER ------------------------------------------------------------------

# This script assumes the working directory is set to "eric-figs/"

# Dependencies and helper functions in header.R

source("./R/header.R")

# Add these to Header.R
library(lubridate)
library(plotly)
library(padr)
library(imputeTS)

# Set UTC time
localTZ <- Sys.timezone()
Sys.setenv(TZ='UTC')

## To Closeout script

# Sys.setenv(TZ=localTZ)

# DATA --------------------------------------------------------------------

# Project Tycho level 1 data (from data.gov)

# local copy of Project Tycho level 1 data
tychoL1 <- read.csv('data/ProjectTycho_Level1_v1.0.0.csv', stringsAsFactors = F, na.strings = "\\N")

# # repository copy of Project Tycho level 1 data
# tychoL1 <- read.csv(url("https://www.healthdata.gov/sites/default/files/ProjectTycho_Level1_v1.0.0.csv", stringsAsFactors = F, na.strings = "\\N"))
# # Metadata: 
# browseURL('https://catalog.data.gov/harvest/object/cb73ca20-e127-4b96-8a48-646d0d4a606f')

# add column for epiweek enddate; calculate epiweek enddate from epiweek
# cdcweekToDate() calculates the start date from the 6 digit epiweek code (yyyyww)
# this code is inefficeint and shoudl be vectorized.
tychoL1 <- cbind(enddate=NA,tychoL1)
tychoL1$enddate <- do.call("c",
                           lapply(tychoL1$epi_week, cdcweekToDate, weekday=6)
                           )

# extract CA Measles cases and pad missing weeks
tychoL1.measles.CA <- pad(
  tychoL1[(tychoL1$disease == 'MEASLES') & (tychoL1$state=='CA'), c("enddate","epi_week","cases")],
  interval = "week"
  )

# impute cases with Kalman filter
tychoL1.measles.CA$cases.imp <- round(na.kalman(as.numeric(tychoL1.measles.CA.pad$cases)))

# convert to zoo time series
tychoL1.measles.CA.cases.imp.zoo <- zoo(tychoL1.measles.CA$cases.imp,tychoL1.measles.CA$enddate)


# Set windows -----------------------------------------------------

# date of peak of epidemic
peakdate <- cdcweekToDate(year=1990, week=34, weekday=6)


bandwidth_weeks <- 104
bandwidth <- bandwidth_weeks * 7  # 104 weeks = 104*7 days
# statsend <- peakdate-bandwidth
statsend <- cdcweekToDate(year=1989, week=34, weekday=6)
plotxmin.year <- 1984
plotxmax.year <- 1992

# cdcweekToDate() calculates the date from the epiweek (year and week)
plotxmin <- cdcweekToDate(year=plotxmin.year,week=plotxmin.week, weekday=6)
plotxmax <- cdcweekToDate(year=plotxmax.year,week=plotxmin.week, weekday=6)

plotymin <- 0
plotymax <- max(
  window(tychoL1.measles.CA.cases.imp.zoo,start=plotxmin,end=plotxmax),
  window(tychoL1.measles.CA.cases.imp.zoo.2wk,start=plotxmin,end=peakdate),
  window(tychoL1.measles.CA.cases.imp.zoo.4wk,start=plotxmin,end=peakdate)
)


# Window Data -------------------------------------------------------------


tychoL1.measles.CA.cases.imp.zoo.statswindow <- window(tychoL1.measles.CA.cases.imp.zoo,
                                                  start=plotxmin-bandwidth,
                                                  end=statsend)


# aggregate biweekly and 4-weekly -----------------------------------------


tmpfcn <- function(x,agg.interval=1) {
  index <- index(x)
  interval <- mean(diff(index))
  agg.groups <- rep(seq(index[length(index)], index[1], -interval*agg.interval), each=agg.interval)
  agg.groups <- rev(agg.groups[1:length(x)])
  out <- aggregate(x, by = agg.groups, sum)
  return(out)
}
tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk <- tmpfcn(tychoL1.measles.CA.cases.imp.zoo.statswindow,2)
tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk <- tmpfcn(tychoL1.measles.CA.cases.imp.zoo.statswindow,4)

tychoL1.measles.CA.cases.imp.zoo.2wk <- tmpfcn(tychoL1.measles.CA.cases.imp.zoo,2)
tychoL1.measles.CA.cases.imp.zoo.4wk <- tmpfcn(tychoL1.measles.CA.cases.imp.zoo,4)


# Get Stats ---------------------------------------------------------------


analysis_params <- list(
  center_trend = "local_constant",
  stat_trend = "local_constant",
  center_kernel = "uniform",
  stat_kernel = "uniform" ,
  center_bandwidth = bandwidth_weeks,
  stat_bandwidth = bandwidth_weeks,
  lag = 1
)

# weekly reports, bandwidth = 104 weeks
analysis_params$center_bandwidth <- analysis_params$stat_bandwidth <- bandwidth_weeks
# No edge
tychoL1.measles.CA.cases.imp.zoo.statswindow.stats <- analysis(tychoL1.measles.CA.cases.imp.zoo.statswindow, analysis_params)
# complete
tychoL1.measles.CA.cases.imp.zoo.stats <- analysis(tychoL1.measles.CA.cases.imp.zoo, analysis_params)

# bi-weekly reports, bandwidth = 52 bi-weeks = 104 weeks
analysis_params$center_bandwidth <- analysis_params$stat_bandwidth <- bandwidth_weeks/2
# No edge
tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk.stats <- analysis(tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk, analysis_params)
# complete
tychoL1.measles.CA.cases.imp.zoo.2wk.stats <- analysis(tychoL1.measles.CA.cases.imp.zoo.2wk, analysis_params)

# four-weekly reports, bandwidth = 26 four-weeks = 104 weeks
analysis_params$center_bandwidth <- analysis_params$stat_bandwidth <- bandwidth_weeks/4 
# No edge
tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk.stats <- analysis(tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk, analysis_params)
# complete
tychoL1.measles.CA.cases.imp.zoo.4wk.stats <- analysis(tychoL1.measles.CA.cases.imp.zoo.4wk, analysis_params)


# Visual Variables --------------------------------------------------------

# Y ticks
y.tick.interval <- 250 
y.ticks <- seq(from = plotymin, to = plotymax, by = y.tick.interval)

# Y ticks
top.y.ticks <- seq(-1,1,.5)

## Colors - REVIEW
color.7day.snap <- unipalette['black']
color.7day.snap.fill <- rgb(0,0,0,.1)
color.7day.imperfect <- unipalette['pink']
color.30day.imperfect <- unipalette['red']
color.crit <- unipalette['green']
color.shaded <- rgb(0,0,0,.1)

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




# Init Figure -------------------------------------------------------------

# PDF output
pdf(
  file = "./output/plots/fig1_revised_alt.pdf",
  title = "Figure 1", # displayed in title bar of PDF Reader
  width = figure.widths['page'], # inches.  Must fit publisher's min and max figure dimensions
  height = figure.heights['page']*.7, # inches.  Must fit publisher's min and max figure dimensions
  family = font.family, 
  pointsize = font.size.normal # default size of text (points).
)

# init figure region
par(mar=margins)
par(lend="butt")


# Top Panel ---------------------------------------------------------------

## Initialize plot of stats
par(fig=top.pannel)
plot(0,0, type='n', axes=FALSE, ann=FALSE, yaxt="n",
     xlim=c(plotxmin,plotxmax), ylim=range(top.y.ticks))

## grid
abline(h=top.y.ticks,col='grey')

## weekly reports
tmp.zoo <- zoo(tychoL1.measles.CA.cases.imp.zoo.statswindow.stats$stats$autocorrelation, 
               time(tychoL1.measles.CA.cases.imp.zoo.statswindow))
lines(window(tmp.zoo, end=statsend), 
      col=color.7day.snap, lwd=lwd.I)

## bi-weekly reports
tmp.zoo <- zoo(tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk.stats$stats$autocorrelation, 
               time(tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk))
lines(window(tmp.zoo, end=statsend), 
      col=color.7day.imperfect, lwd=lwd.7day)

## four-weekly reports
tmp.zoo <- zoo(tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk.stats$stats$autocorrelation, 
               time(tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk))
lines(window(tmp.zoo, end=statsend),
      col=color.30day.imperfect, lwd=lwd.30day)

## Annotation
title(ylab="Autocorrelation", line = 3)

# ### tuas
# text(x=19.75*365,
#      y=c(.925,.3,-.1),
#      adj=1,
#      cex=font.scales['M'],
#      col=c(color.7day.snap,
#            color.7day.imperfect,
#            color.30day.imperfect),
#      labels=c(
#        TeX(sprintf("$\\tau = $ %.2f", data.7day.snap$taus[[1]]$autocorrelation)),
#        TeX(sprintf("$\\tau = $ %.2f", data.7day.imperfect$taus[[1]]$autocorrelation)),
#        TeX(sprintf("$\\tau = $ %.2f", data.30day.imperfect$taus[[1]]$autocorrelation))
#      )
# )

## shade unanalyzed times 

rect(statsend, 0, plotxmax+bandwidth, 1, border = F, col = color.shaded)

## Axes

axis(side = 1, # 1 specifies bottom axis - major ticks
     at = seq(plotxmin,plotxmax,"years"), # tick locations (in data units).
     labels = FALSE,
     tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
     line = 1, # axis position (in lines)
     lty = "solid",
     lwd = 0, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
)
axis(side = 3, # 1 specifies bottom axis - major ticks
     at = seq(plotxmin,plotxmax,"years"), # tick locations (in data units).
     labels = year(seq(plotxmin,plotxmax,"years")),
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
         "Four-weekly\nReports",
         "Bi-weekly\nReports", 
         "Weekly\nReports"
       ),
       col=c(
         # color.ac.I.lag1,
         color.30day.imperfect,
         color.7day.imperfect,
         color.7day.snap
       ),
       lwd=c(
         # lwd.I,
         lwd.30day,
         lwd.7day,
         lwd.I),
       lty=1,
       bty='n')


# Bottom Panel ------------------------------------------------------------

## init plot
par(fig=bottom.pannel, new=TRUE)
plot(0,0, type='n', axes=FALSE, ann=FALSE, yaxt="n",
     xlim=c(plotxmin,plotxmax), ylim=c(plotymin,plotymax))

# 
# tychoL1.measles.CA.cases.imp.zoo
# tychoL1.measles.CA.cases.imp.zoo.statswindow
# tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk
# tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk

## plot time series of weekly snapshots of number infected
filledLine <- function (data, lwd=.1, border, col) {
  data <- tychoL1.measles.CA.cases.imp.zoo
  x <- time(data)
  y <- as.numeric(data)
  data.x <- c(head(x,1),x,tail(x,1))
  # poly.y <- c(0,y,0)
  data.x <- c(rbind(x,x)) ## double each time
  data.y <- c(0,
              c(rbind(y[2:length(y)],y[2:length(y)])),
              0)
  polygon(x=data.x,y=data.y, lwd=lwd, border=border, col=col)
}

lines(tychoL1.measles.CA.cases.imp.zoo.4wk, col=color.30day.imperfect, type='S')
lines(tychoL1.measles.CA.cases.imp.zoo.2wk, col=color.7day.imperfect, type='S')
filledLine(data=tychoL1.measles.CA.cases.imp.zoo,border=color.7day.snap, col=color.7day.snap.fill)


# filledLine(data=tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk,border=color.30day.imperfect, col=color.30day.imperfect)
# filledLine(data=tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk,border=color.7day.imperfect, col=color.7day.imperfect)

# filledLine(data=tmp.zoo,border=color.7day.snap, col="black")

# lines(tychoL1.measles.CA.cases.imp.zoo, col="black", type='S', lwd=.25)

# 
# 
# 
# poly <- tychoL1.measles.CA.cases.imp.zoo
# x <- time(poly)
# y <- as.numeric(poly)
# poly.x <- c(head(x,1),x,tail(x,1))
# # poly.y <- c(0,y,0)
# poly.x <- c(rbind(x,x)) ## double each time
# poly.y <- c(0,
#             c(rbind(y[2:length(y)],y[2:length(y)])),
#             0)
# polygon(x=poly.x,y=poly.y, lwd=.1, border=color.7day.snap, col=color.7day.snap.fill)
# 
# ## plot time series of cases, aggregated Weekly (imperfect reporting), filled staircase
# poly <- data.7day.imperfect$data[[1]]
# x <- time(poly) 
# y <- as.numeric(poly)
# poly.x <- c(rbind(x,x)) ## double each time
# poly.y <- c(0,
#             c(rbind(y[2:length(y)],y[2:length(y)])),
#             0)
# polygon(x=poly.x,y=poly.y, lwd=lwd.7day, border=color.7day.imperfect, col=color.7day.imperfect)
# 
# ## plot time series of cases, aggregated Monthly (imperfect reporting), staircase
# lines(data.30day.imperfect$data[[1]], 
#       type='s', pch=20, lwd=lwd.30day, col=color.30day.imperfect)

# shade unplotted data region

rect(statsend, 0, plotxmax+bandwidth, plotymax, border = F, col = color.shaded)

## Axes

axis(side = 1, # 1 specifies bottom axis - major ticks
     at = seq(plotxmin,plotxmax,"year"), # tick locations (in data units).
     labels = year(seq(plotxmin,plotxmax,"year")),
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
# 
# ## Legend
# par(lend="butt")
# legend("topleft", xpd=NA, inset=c(1.01,0), xjust=0, yjust=0, cex=font.scales['XS'], y.intersp = 1.5,
#        legend=c(
#          "Four-weekly Reports", 
#          "Bi-weekly Reports", 
#          "Weekly Reports"
#        ),
#        col=c(color.30day.imperfect,color.7day.imperfect,color.7day.snap.fill),
#        lty=1,
#        lwd=c(lwd.30day,lwd.7day,lwd.7day.snap),
#        bty='n')

# Close PDF ---------------------------------------------------------------

invisible(dev.off())
