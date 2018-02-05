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
tychoL1 <- read.csv(url("https://www.healthdata.gov/sites/default/files/ProjectTycho_Level1_v1.0.0.csv", stringsAsFactors = F))

# Metadata:
# browseURL('https://catalog.data.gov/harvest/object/cb73ca20-e127-4b96-8a48-646d0d4a606f')
# browseURL('https://jsoneditoronline.org/?url=https%3A%2F%2Fcatalog.data.gov%2Fharvest%2Fobject%2Fcb73ca20-e127-4b96-8a48-646d0d4a606f')


# local copy of Project Tycho level 1 data
tychoL1 <- read.csv('data/ProjectTycho_Level1_v1.0.0.csv', stringsAsFactors = F)

# add column for epiweek enddate to use as index
tychoL1 <- cbind(enddate=NA,tychoL1)
tychoL1$enddate <- as.Date(sapply(tychoL1$epi_week, cdcweekToDate))+6

# extract CA Measles time series 
tychoL1.measles.CA <- tychoL1[(tychoL1$disease == 'MEASLES') & (tychoL1$state=='CA'), ]

# pad missing weeks
tychoL1.measles.CA.pad <- pad(tychoL1.measles.CA, interval = "week")

# impute cases with Kalman filter
# cases <- as.numeric(tychoL1.measles.CA.pad$cases)
tychoL1.measles.CA.pad$cases.imp <- round(na.kalman(as.numeric(tychoL1.measles.CA.pad$cases)))
# convert to zoo time series
tychoL1.measles.CA.imp.zoo <- zoo(cases.imp,as.numeric(tychoL1.measles.CA$enddate))

# Test Plots -------------------------------------------------------------------


# test plot, padded with NA's
plot(y = tychoL1.measles.CA.pad$cases, x = tychoL1.measles.CA.pad$enddate, type='l',
     ylim = c(0,1200)
)
plot(y = tychoL1.measles.CA.pad$cases.imp, x = tychoL1.measles.CA.pad$enddate, type='l',
     ylim = c(0,1200)
)

plot(y = tychoL1.measles.CA.pad$cases.imp, col=2, x = tychoL1.measles.CA.pad$enddate, type='l',
     xlim = c(cdcweekToDate(198801), cdcweekToDate(199101)),
     ylim = c(0,1200)
)
lines(y = tychoL1.measles.CA.pad$cases, col=1, lwd=1.5, x = tychoL1.measles.CA.pad$enddate, type='l')


barplot(tychoL1.measles.CA.pad$cases.imp, col=2, x = tychoL1.measles.CA.pad$enddate, type='l',
     xlim = c(cdcweekToDate(198801), cdcweekToDate(199101)),
     ylim = c(0,1200)
)
lines(y = tychoL1.measles.CA.pad$cases, col=1, lwd=1.5, x = tychoL1.measles.CA.pad$enddate, type='l')

# plot imputed case time series as zoo 
tmpzoo <- window(tychoL1.measles.CA.imp.zoo, start=cdcweekToDate(199001),end=cdcweekToDate(199401))
plot.zoo(tmpzoo, type='l', ylim = c(0,1200))

# # initial experiment to set week = 1/52(year)
# date <- measles$YEAR+(measles$WEEK-1)/52
# 
# CA <- measles$CALIFORNIA
# plot(date, CA, type='l', 
#      xlim = c(cdcweekToDate(198101), cdcweekToDate(199101)),
#      ylim = c(0,1200)
#      )