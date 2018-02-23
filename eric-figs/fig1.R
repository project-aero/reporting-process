# HEADER ------------------------------------------------------------------

# This script assumes the working directory is set to "eric-figs/"

## Dependencies and helper functions in header.R
source("./R/header.R")

## Additional packages for managing real data
library(lubridate)
library(padr)
library(imputeTS)

## Set UTC time
localTZ <- Sys.timezone()
Sys.setenv(TZ='UTC')

# DATA --------------------------------------------------------------------

# Project Tycho level 1 data

## local copy of Project Tycho level 1 data from:
  # browseURL:("https://www.healthdata.gov/sites/default/files/ProjectTycho_Level1_v1.0.0.csv")
  # Metadata: 
  # browseURL('https://catalog.data.gov/harvest/object/cb73ca20-e127-4b96-8a48-646d0d4a606f')

tychoL1 <- read.csv('data/ProjectTycho_Level1_v1.0.0.csv', stringsAsFactors = F, na.strings = "\\N")

# The original copy of Project Tycho level 1 data may be loaded directly using the line of code below:
# tychoL1 <- read.csv(url("https://www.healthdata.gov/sites/default/files/ProjectTycho_Level1_v1.0.0.csv", stringsAsFactors = F, na.strings = "\\N"))

## Add column for epiweek enddate; calculate epiweek enddate from epiweek
  # Note: cdcweekToDate() calculates the start date from the 6 digit epiweek code (yyyyww)
  # Note: cdcweekToDate() is not a vectorized function and may take several minutes to run on a large dataset.
tychoL1 <- cbind(enddate=NA,tychoL1)
tychoL1$enddate <- do.call("c",
                           lapply(tychoL1$epi_week, cdcweekToDate, weekday=6)
)

## Extract CA Measles cases and pad missing weeks with NA's
tychoL1.measles.CA <- pad(
  tychoL1[(tychoL1$disease == 'MEASLES') & (tychoL1$state=='CA'), c("enddate","epi_week","cases")],
  interval = "week"
)

## Impute cases with Kalman filter
tychoL1.measles.CA$cases.imp <- round(na.kalman(as.numeric(tychoL1.measles.CA$cases)))

## Convert to zoo time series
tychoL1.measles.CA.cases.imp.zoo <- zoo(tychoL1.measles.CA$cases.imp,tychoL1.measles.CA$enddate)

# Set windows -----------------------------------------------------

## Statistics bandwidth
bandwidth_weeks <- 26
bandwidth <- bandwidth_weeks * 7  # 104 weeks = 104*7 days

## date of peak of epidemic
peakdate <- cdcweekToDate(year=1990, week=34, weekday=6)

## date prior to peak of outbreak from which past data is used to calculate statistics.
obsdate <- cdcweekToDate(year=1990, week=34-31, weekday=6)

## end date for plot of stats
statsend <- obsdate

## set range of years for overall plat
plotxmin.year <- 1984
plotxmax.year <- 1992

## cdcweekToDate() calculates the date from the epiweek (year and week)
plotxmin <- cdcweekToDate(year=plotxmin.year, week=1, weekday=6)
plotxmax <- cdcweekToDate(year=plotxmax.year, week=1, weekday=6)

# Window data -------------------------------------------------------------

## window data to plot region
tychoL1.measles.CA.cases.imp.zoo.plotwindow <- window(tychoL1.measles.CA.cases.imp.zoo,
                                                       start=plotxmin,
                                                       end=plotxmax)

# Aggregate biweekly and 4-weekly -----------------------------------------

tmpfcn <- function(x,agg.interval=1) {
  index <- index(x)
  interval <- mean(diff(index))
  agg.groups <- rep(seq(index[1], index[length(index)], interval*agg.interval), each=agg.interval)
  agg.groups <- agg.groups[1:length(x)]
  out <- aggregate(x, by = agg.groups, sum)
  return(out)
}

# aggregate windowed dataset
tychoL1.measles.CA.cases.imp.zoo.plotwindow.2wk <- tmpfcn(tychoL1.measles.CA.cases.imp.zoo.plotwindow,2)
tychoL1.measles.CA.cases.imp.zoo.plotwindow.4wk <- tmpfcn(tychoL1.measles.CA.cases.imp.zoo.plotwindow,4)


# Window aggregated data for analysis -------------------------------------

## window to end at a chosen "observation date."
  # This windowed data is used below to calculate statistics

tychoL1.measles.CA.cases.imp.zoo.statswindow <- window(tychoL1.measles.CA.cases.imp.zoo.plotwindow,
                                                       end=obsdate)
tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk <- window(tychoL1.measles.CA.cases.imp.zoo.plotwindow.2wk,
                                                       end=obsdate)
tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk <- window(tychoL1.measles.CA.cases.imp.zoo.plotwindow.4wk,
                                                       end=obsdate)

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

# weekly reports
analysis_params$center_bandwidth <- analysis_params$stat_bandwidth <- bandwidth_weeks
tychoL1.measles.CA.cases.imp.zoo.statswindow.stats <- analysis(tychoL1.measles.CA.cases.imp.zoo.statswindow, analysis_params)

# bi-weekly reports
analysis_params$center_bandwidth <- analysis_params$stat_bandwidth <- bandwidth_weeks/2
tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk.stats <- analysis(tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk, analysis_params)

# four-weekly reports
analysis_params$center_bandwidth <- analysis_params$stat_bandwidth <- bandwidth_weeks/4 
tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk.stats <- analysis(tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk, analysis_params)


# Set Y ranges ------------------------------------------------------------

cases.min <- 0
cases.max <- max(
  window(tychoL1.measles.CA.cases.imp.zoo.plotwindow,start=plotxmin,end=plotxmax),
  window(tychoL1.measles.CA.cases.imp.zoo.plotwindow.2wk,start=plotxmin,end=plotxmax),
  window(tychoL1.measles.CA.cases.imp.zoo.plotwindow.4wk,start=plotxmin,end=plotxmax)
)
cases.tick.interval <- 250 
cases.ticks <- seq(from = cases.min, by = cases.tick.interval, length.out = ceiling(cases.max/cases.tick.interval)+1)

var.min <- 0
var.max.actual <- max(
  tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk.stats$stats$variance,
  tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk.stats$stats$variance,
  tychoL1.measles.CA.cases.imp.zoo.statswindow.stats$stats$variance
)
var.max <- 10000 # manually set plot limit

var.tick.interval <- 2000
var.ticks <- seq(from = var.min, to = var.max, by = var.tick.interval)

# Visual Variables --------------------------------------------------------

## Colors
color.1wk <- unipalette['black']
color.2wk <- unipalette['lightblue']
color.4wk <- unipalette['red']
color.omitted <- rgb(0,0,0,.1)

## Line widths
lwd.1wk <- 0.5
lwd.2wk <- 1.0
lwd.4wk <- 2.0

## Filled line border widths and fill alpha
border.1wk <- 0.5
border.2wk <- 1.0
border.4wk <- 2.0
fill.alpha.1wk <- 0.7
fill.alpha.2wk <- 0.2
fill.alpha.4wk <- 0.0

## Typography
font.family <- "Times"
font.sizes <- seq(from = 8, # publisher's minimum point size (points)
                  to = 12, # publisher's maximum point size (points) 
                  length.out = 5)
font.size.normal <- mean(font.sizes)
font.scales <- font.sizes/mean(font.sizes)
names(font.scales) <- names(font.sizes) <- c("XS", "S", "M", "L", "XL")

## Figure dimensions
figure.widths <- c(min=2.63, page=7.5, column=5.2) # in inches, as defined by publisher
figure.heights <- c(min=1, page=8.75) # in inches, as defined by publisher

## Margins and Figure Bounds
margins = c(4,5,4,8)+0.1
top.pannel <- c(0,1,.45,1) # panel bounds: x0,x1,y0,y1 as fraction of figure region
bottom.pannel <- c(0,1,0,.6) # panel bounds: x0,x1,y0,y1 as fraction of figure region

# Init Figure -------------------------------------------------------------

## PDF output
pdf(
  file = "./output/plots/fig1.pdf",
  title = "Figure 1", # displayed in title bar of PDF Reader
  width = figure.widths['page'], # inches.  Must fit publisher's min and max figure dimensions
  height = figure.heights['page']*.7, # inches.  Must fit publisher's min and max figure dimensions # change to .7
  family = font.family, 
  pointsize = font.size.normal # default size of text (points).
)

## init figure region
par(mar=margins)
par(lend="butt", ljoin="mitre")

# Top Panel (Variance) ----------------------------------------------------

## Initialize plot of stats
par(fig=top.pannel)
plot(0,0, type='n', axes=FALSE, ann=FALSE, yaxt="n",xaxs='i',
     xlim=c(plotxmin,plotxmax), ylim=range(var.ticks))

## Shade unanalyzed times 
rect(obsdate, var.min, plotxmax+bandwidth, var.max, border = F, col = color.omitted)

## Grid
abline(h = var.ticks, col = 'grey')

## weekly reports
tmp.zoo <- zoo(tychoL1.measles.CA.cases.imp.zoo.statswindow.stats$stats$variance, 
               time(tychoL1.measles.CA.cases.imp.zoo.statswindow))
lines(tmp.zoo, col=color.1wk, lwd=lwd.1wk)

## bi-weekly reports
tmp.zoo <- zoo(tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk.stats$stats$variance, 
               time(tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk))
lines(tmp.zoo, col=color.2wk, lwd=lwd.2wk)

## four-weekly reports
tmp.zoo <- zoo(tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk.stats$stats$variance, 
               time(tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk))
lines(tmp.zoo, col=color.4wk, lwd=lwd.4wk)

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
     at = var.ticks, # tick locations (in data units).
     labels = var.ticks, # vector of tick labels
     las = 1, # label orientation (0:parall., 1:horiz., 2:perp., 3:vert.)
     tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
     line = .5, # axis position, in lines
     lty = "solid",
     lwd = 0, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
)

## Axis Labels
mtext(text="Year", side=3, line = 2.5)
title(ylab="Variance", line = 4)

## Taus (Tau depends on all plotted values of the variance)
text(x=c(statsend-104*7,statsend-18*7,statsend+6*7),
     y=c(9500,8500,5500),
     adj=0,
     cex=font.scales['M'],
     col=c(color.4wk,
           color.2wk,
           color.1wk),
     labels=c(
       TeX(sprintf("$\\tau = $ %.2f", tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk.stats$taus$variance)),
       TeX(sprintf("$\\tau = $ %.2f", tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk.stats$taus$variance)),
       TeX(sprintf("$\\tau = $ %.2f", tychoL1.measles.CA.cases.imp.zoo.statswindow.stats$taus$variance))
     )
)

## Legend
legend("topleft", xpd=NA, inset=c(1.01,0), xjust=0, yjust=0, cex=font.scales['XS'], y.intersp = 3,
       legend=c("Four-weekly\nReports", "Bi-weekly\nReports", "Weekly\nReports"),
       col=c(color.4wk,color.2wk,color.1wk),
       lwd=c(lwd.4wk,lwd.2wk,lwd.1wk),
       lty=1,
       bty='n'
       )

# Bottom Panel ------------------------------------------------------------

## Initialize plot of cases
par(fig=bottom.pannel, new=TRUE)
plot(0,0, type='n', axes=FALSE, ann=FALSE, yaxt="n", xaxs='i', yaxs = 'i',
     xlim=c(plotxmin,plotxmax), ylim=range(cases.ticks))

## shade unanalyzed times 
rect(obsdate, range(cases.ticks)[1], plotxmax+bandwidth, range(cases.ticks)[2], border = F, col = color.omitted)

## Grid
abline(h = cases.ticks, col = 'grey')

## Plot time series of cases, aggreaged weekly, bi-weekly, four-weekly

filledLine <- function (data, lwd, color, fill.alpha) {
  x <- time(data)
  y <- as.numeric(data)
  data.x <- c(head(x,1),x,tail(x,1))
  data.x <- c(rbind(x,x)) ## double each time
  data.y <- c(0,
              c(rbind(y[2:length(y)],y[2:length(y)])),
              0)
  polygon(x=data.x,y=data.y, lwd=lwd, border=color, col=alpha(color,fill.alpha))
}

# four-weekly
filledLine(data=tychoL1.measles.CA.cases.imp.zoo.plotwindow.4wk, lwd = border.4wk, color=color.4wk, fill.alpha=fill.alpha.4wk)
# bi-weekly
filledLine(data=tychoL1.measles.CA.cases.imp.zoo.plotwindow.2wk, lwd = border.2wk, color=color.2wk, fill.alpha=fill.alpha.2wk)
# weekly
filledLine(data=tychoL1.measles.CA.cases.imp.zoo.plotwindow, lwd = border.1wk, color=color.1wk, fill.alpha=fill.alpha.1wk)

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
     at = cases.ticks, # tick locations (in data units).
     labels = cases.ticks, # vector of tick labels
     las = 1, # label orientation (0:parall., 1:horiz., 2:perp., 3:vert.)
     tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
     line = .5, # axis position, in lines
     lty = "solid",
     lwd = 0, # axis line weight
     lwd.ticks = 1, # tick line weight
     col = "grey" # axis line color
)

## Axis Labels
title(ylab = "Number", line = 4)
title(xlab = "Year", line = 2.5)

# Close PDF ---------------------------------------------------------------

invisible(dev.off())


# FOOTER ------------------------------------------------------------------

## reset Timezone to previous stored local timezone
Sys.setenv(TZ=localTZ)

