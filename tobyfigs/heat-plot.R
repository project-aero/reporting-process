#Script for heat plot 
# 3-#2017-07-05 -- Toby

library(colorspace)
library(lattice)
library(fields)
library(dplyr)
library(Cairo)
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





## Eric's AUC gradient ##
# Given error of ~0.5, set number of levels to 10.  For a higher precision, use nlevels = 20 
nlevels = 100  # number of bins = number of colors
# Color scale with {colorspace}
AUC.colors <- diverge_hcl(
  n=nlevels,  # number of colors
  h = c(45, 225),  # Hues (low, hi)
  c = 100, # fixed Chroma or Chroma Range (edges, center)
  l = c(90, 10),  # Lightness range (edges, center) 
  power = 1  # exponent
)
#Gray colour used in paper:
gray_colour <- rgb(0,0,0,.25)
dl <- seq(0,1,0.01)

## AUC data ##
load("res.RData")
res <- res[which(res$reporting_prob > 0.01),]


heat_map_plot <- function(bw){
    
  ##Filter for bandwidth Choose either 35 or 100
  res <- filter(res, bandwidth == bw)
  
  ## Results for 7 day recovery
  res.7.wkly <- res[ which(res$infectious_days==7 
                           & res$aggregation_days == 7), ]
  
  res.7.mon <- res[ which(res$infectious_days==7 
                          & res$aggregation_days == 30), ]
  
  ## Results for 30 day recovery
  res.30.wkly <- res[ which(res$infectious_days==30 
                            & res$aggregation_days == 7), ]
  
  res.30.mon <- res[ which(res$infectious_days==30 
                           & res$aggregation_days == 30), ]
  
  ## function to organize data for heat plots ##
  org_data <- function(dat){
    binLevs <- 3
    ewsLevs <- ncol(dat[,14:23])
    repLevs <- length(unique(dat$reporting_prob))
    
    tmpList <- list() # For each EWS
    # Make each level of dispersion a new column
    for(i in 1:ewsLevs){
      cl <- 13 + i
      tmp <- dat[, c(1:13, cl)]
      tmp=filter(tmp, neg_bin_k!=.1)
      tmp=filter(tmp, neg_bin_k!=10)
      tmp <- arrange(tmp, neg_bin_k)
      
      tmp <- matrix(tmp[,14], ncol=repLevs, byrow = T)
      row.names(tmp)=c("100",  "1" , ".01")
      
      colnames(tmp) <- round(unique(dat$reporting_prob), 5)
      tmpList[[i]] <- tmp
    }
    
    # format for heat plot
    newDat <- do.call(rbind, tmpList)
    newDat <- t(newDat)
    
    kurt=newDat[,c(28:30)]
    skew=newDat[,c(25:27)]
    cov =newDat[,c(22:24)]
    ind =newDat[,c(19:21)]
    men =newDat[,c(16:18)]
    dec =newDat[,c(13:15)]
    acf =newDat[,c(10:12)]
    acv =newDat[,c(7:9)]
    vdif =newDat[,c(4:6)]
    var =newDat[,c(1:3)]
    ord = cbind(men, var, vdif,
                ind, acv, acf,
                dec, cov, skew,
                kurt)
    return(ord)
  }
  
  yLabs <- c("Mean", "Variance", "1st Diff. Var.", 
             "Ind. of Dis.", "Autocovar.", "Autocorr.",
             "Decay time", "Coeff. Var", "Skewness",
             "Kurtosis")
  
  # For plotting
  nRepProb <- length(unique(res$reporting_prob))
  xVals <- seq(0,nRepProb-1)/(nRepProb-1)  # seq(0,20)/20
  xLabVals <- xVals[seq(1, length(xVals), 6)]
  xLabs <- round(unique(res$reporting_prob), 2)
  xLabsNew <- xLabs[seq(1, length(xLabs), 6)]
  yVals <- seq(2.5,50, by=5)/50
  
  addLines <- function(clr="white", lwd=1.5){
    abline(h=18/200, col=clr, lwd=lwd)
    abline(h=39/200, col=clr, lwd=lwd)
    abline(h=59/200, col=clr, lwd=lwd)
    abline(h=79.5/200, col=clr, lwd=lwd)
    abline(h=100/200, col=clr, lwd=lwd)
    abline(h=120/200, col=clr, lwd=lwd)
    abline(h=141/200, col=clr, lwd=lwd)
    abline(h=161/200, col=clr, lwd=lwd)
    abline(h=181.5/200, col=clr, lwd=lwd)
  }
  
  dispLabs3 <- function(){
    #dispersion labels
    axis(4, at = seq(46.5,50, length.out = 3)/50, 
         labels = rep("", 3), tck=-.1, col.ticks = gray_colour,
         col = "white", pos = 1.05,family = font.family, cex=font.scales['S'])
    mtext(side=4, "less dispersed", line=.65, las=1, family = font.family, cex=font.scales['XS'], at=1.07, xpd=T)
    mtext(side=4, "100", line=3, las=1, family = font.family, cex=font.scales['XS'], at=1.01, xpd=T)
    mtext(side=4, " 1 ", line=3.1, las=1, family = font.family, cex=font.scales['XS'], at=.965, xpd=T)
    mtext(side=4, "0.01", line=3, las=1,family = font.family, cex=font.scales['XS'], at=.92, xpd=T)
    mtext(side=4, "more dispersed", line=.51, las=1,family = font.family, cex=font.scales['XS'], at=.86, xpd=T)
  }
  
  addLines3 <- function(clr="white", lwd=1.5){
    abline(h=17/200, col=clr, lwd=lwd)
    abline(h=38/200, col=clr, lwd=lwd)
    abline(h=59/200, col=clr, lwd=lwd)
    abline(h=79.5/200, col=clr, lwd=lwd)
    abline(h=100/200, col=clr, lwd=lwd)
    abline(h=120/200, col=clr, lwd=lwd)
    abline(h=141/200, col=clr, lwd=lwd)
    abline(h=162/200, col=clr, lwd=lwd)
    abline(h=182/200, col=clr, lwd=lwd)
  }
  
  # Data for heat plots
  H5 <- org_data(res.7.wkly)
  H6 <- org_data(res.7.mon)
  
  H7 <- org_data(res.30.wkly)
  H8 <- org_data(res.30.mon)
  
  # PDF output
  ## specify 
  tiff(
    file = paste("heat-plot",bw,".tiff", sep=""),
    #type="tiff",
    title = "Fig. 5", # displayed in title bar of PDF Reader
    width = figure.widths['column'], # full width, in inches
    height = figure.heights['page']*.7, # 70% of full height, in inches
    units = "in",
    compression="lzw", 
    family = font.family, 
    pointsize = font.size.normal, # default size of text (points).
    res = 300
  )

  
  
  # separate heat plots
  layout(matrix(c(1, 2,
                  3, 4,
                  5, 5), nrow=3, byrow=TRUE),
         heights = c(1,1,.25))
  
  par(oma=c(3,2,2,2), mar=c(1,5,3,3))
  
  image(H5, col=AUC.colors, xlab="", 
        ylab="", axes=F, main="",
        add.expr= abline(h=5/50, col="white"))
  axis(3, at=xLabVals, labels=xLabsNew, col.ticks = gray_colour,
       col = "white", pos = 1.03, family = font.family, cex.axis=font.scales['XL'])
  axis(3, at=xVals, labels=FALSE, tck=-0.05, col.ticks = gray_colour, 
       col = "white", pos = 1.03, family = font.family, cex.axis=font.scales['L'])
  axis(2, at=yVals, labels=yLabs, las=1, col.ticks = gray_colour, 
       col = "white", pos = -0.05, family = font.family, cex.axis=font.scales['L'])
  mtext(side=3, text="Aggregated weekly", line=3, cex=font.scales['L'])
  mtext(side=2, expression(paste(gamma, "=1/7")), line=5.5, cex=font.scales['L'])
  addLines3()
  dispLabs3()
  
  par(mar=c(1,4,3,4))
  image(H6, col=AUC.colors, xlab="", 
        ylab="", axes=F)
  axis(3, at=xLabVals, labels=xLabsNew, col.ticks = gray_colour,
       col = "white", pos = 1.03,family = font.family, cex.axis=font.scales['M'])
  axis(3, at=xVals, labels=FALSE, tck=-0.05, col.ticks = gray_colour,
       col = "white", pos = 1.03, family = font.family, cex.axis=font.scales['M'])
  mtext(side=3, text="Aggregated monthly", line=3, cex=font.scales['L'])
  axis(4, at=yVals, labels=yLabs, las=1, col.ticks = gray_colour,
       col = "white", pos =  1.05, family = font.family, cex.axis=font.scales['M'])
  addLines3()
  axis(2, at = seq(46.5,50,length.out = 3)/50, 
       labels = rep("", 3), tck=-.1,col.ticks = gray_colour,
       col = "white", pos = -0.05, family = font.family, cex.axis=font.scales['M'])
  
  
  par(mar=c(4,5,0,3))
  image(H7, col=AUC.colors, xlab="", 
        ylab="", axes=F)
  axis(1, at=xLabVals, labels=xLabsNew, col.ticks = gray_colour,
       col = "white", pos = -0.03, family = font.family, cex.axis=font.scales['M'])
  axis(1, at=xVals, labels=FALSE, tck=-0.05, col.ticks = gray_colour,
       col = "white", pos = -0.03, family = font.family, cex.axis=font.scales['M'])
  axis(2, at=yVals, labels=yLabs, las=1, col.ticks = gray_colour,
       col = "white",pos = -0.05, family = font.family, cex.axis=font.scales['M'])
  mtext(side=2, expression(paste(gamma, "=1/30")), line=5.5)
  mtext(side=1, "Reporting probability", line=3, cex=font.scales['M'])
  addLines3()
  dispLabs3()
  
  
  par(mar=c(4,4,0,4))
  image(H8, col=AUC.colors, xlab="", 
        ylab="", axes=F)
  axis(1, at=xLabVals, labels=xLabsNew, col.ticks = gray_colour,
       col = "white", pos = -0.03,family = font.family, cex.axis=font.scales['M'])
  axis(1, at=xVals, labels=FALSE, tck=-0.05, col.ticks = gray_colour,
       col = "white", pos = -0.03, family = font.family, cex.axis=font.scales['M'])
  axis(4, at=yVals, labels=yLabs, las=1, col.ticks = gray_colour,
       col = "white", pos = 1.05, family = font.family, cex.axis=font.scales['M'])
  mtext(side=1, "Reporting probability", line=3, cex=font.scales['M'])
  addLines3()
  axis(2, at = seq(46.5,50, length.out = 3)/50, 
       labels = rep("", 3), tck=-.1, col.ticks = gray_colour, 
       col = "white", pos = -0.05, family = font.family, cex.axis=font.scales['M'])
  
  
  par(mar=c(3,7,1.2,7)) 
  legend_data <- matrix(dl, nrow = length(dl), ncol = 1)
  image(legend_data, col=AUC.colors, xlab="", 
              ylab="", axes=F, useRaster = T)
  axis(1, at=seq(0.0,1,0.2), labels=seq(0.0,1,0.2), col.ticks = gray_colour, 
       col = "white", pos = -1.5, family = font.family, cex.axis=font.scales['M'])
  mtext(side=1, "AUC", line=3, , cex=font.scales['M'])
  dev.off()
}

heat_map_plot(35)
heat_map_plot(100)
