#Script for heat plot 
# 3-#2017-07-05 -- Toby

library(colorspace)
library(lattice)
library(fields)
library(dplyr)

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

##Filter for windowsize
res <- filter(res, bandwidth == 100)

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
xLabs <- round(unique(res$reporting_prob), 2)
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
       labels = rep("", 3), tck=-.1, col.ticks = gray_colour,col = "white", pos = 1.05)
  mtext(side=4, "less dispersed", line=.95, las=1, cex=.6, at=1.05, xpd=T)
  
  mtext(side=4, "100", line=3, las=1, cex=.4, at=1, xpd=T)
  mtext(side=4, " 1 ", line=3.1, las=1, cex=.4, at=.965, xpd=T)
  mtext(side=4, "0.01", line=3, las=1, cex=.4, at=.93, xpd=T)
  
  mtext(side=4, "more dispersed", line=.91, las=1, cex=.6, at=.88, xpd=T)
  
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


# separate heat plots
layout(matrix(c(1, 2,
                3, 4,
                5, 5), nrow=3, byrow=TRUE),
       heights = c(1,1,.25))
layout.show(n=5)

par(oma=c(3,2,2,2), mar=c(1,5,3,3))
#par(mfrow=c(2,2), oma=c(4,6.5,4,4.5), mar=c(1,1,1,1))
image(H5, col=AUC.colors, xlab="", 
      ylab="", axes=F, main="",
      add.expr= abline(h=5/50, col="white"))
axis(3, at=xVals, labels=xLabs, col.ticks = gray_colour, col = "white", pos = 1.03)
axis(2, at=yVals, labels=yLabs, las=1, col.ticks = gray_colour, col = "white", pos = -0.05)
mtext(side=3, text="Aggregated weekly", line=3)
mtext(side=2, expression(paste(gamma, "=1/7")), line=5.5)
addLines3()
dispLabs3()
#box()

#arrows(x0=1.2, x1=1.2, y0=.87, y1=1.06, xpd=T, code=3, length = .05)

par(mar=c(1,4,3,4))
image(H6, col=AUC.colors, xlab="", 
      ylab="", axes=F)
axis(3, at=xVals, labels=xLabs, col.ticks = gray_colour,col = "white", pos = 1.03)
mtext(side=3, text="Aggregated monthly", line=3)
axis(4, at=yVals, labels=yLabs, las=1, col.ticks = gray_colour, col = "white", pos =  1.05)
addLines3()
axis(2, at = seq(46.5,50,length.out = 3)/50, 
     labels = rep("", 3), tck=-.1,col.ticks = gray_colour, col = "white", pos = -0.05)
#box()

par(mar=c(4,5,0,3))
image(H7, col=AUC.colors, xlab="", 
      ylab="", axes=F)
axis(1, at=xVals, labels=xLabs, col.ticks = gray_colour, col = "white",pos = -0.03)
axis(2, at=yVals, labels=yLabs, las=1, col.ticks = gray_colour, col = "white",pos = -0.05)
mtext(side=2, expression(paste(gamma, "=1/30")), line=5.5)
mtext(side=1, "Reporting probability", line=3)
addLines3()
#axis(4, at = seq(46.5,50, length.out = 3)/50, labels = rep("", 3), col.ticks = gray_colour)
dispLabs3()
#box()

par(mar=c(4,4,0,4))
image(H8, col=AUC.colors, xlab="", 
      ylab="", axes=F)
axis(1, at=xVals, labels=xLabs, col.ticks = gray_colour, col = "white", pos = -0.03)
axis(4, at=yVals, labels=yLabs, las=1, col.ticks = gray_colour, col = "white", pos = 1.05)
mtext(side=1, "Reporting probability", line=3)
addLines3()
axis(2, at = seq(46.5,50, length.out = 3)/50, 
     labels = rep("", 3), tck=-.1, col.ticks = gray_colour, col = "white", pos = -0.05)
#box()

par(mar=c(1.2,7,1.2,7)) 
#par(mar= c(0,0,0,0), oma=c(1.,0,0,0))

legend_data <- matrix(dl, nrow = length(dl), ncol = 1)
legend_data <- as.matrix(blur(as.im(legend_data), sigma=6))
image(legend_data, col=AUC.colors, xlab="", 
            ylab="", axes=F, useRaster = T)
axis(1, at=seq(0.0,1,0.2), labels=seq(0.0,1,0.2), col.ticks = gray_colour, col = "white", pos = -1.5)
mtext(side=1, "AUC", line=3)
#image.plot( zlim=c(0,1.0), legend.only=TRUE, horizontal=TRUE,
   #         col=AUC.colors, smallplot=c(.25,.8,0.7,.85), legend.lab = "AUC",
    #        axis.args = c(col = "white", col.ticks = gray_colour, pos =0.),border = NA )

dev.copy(pdf,"heat-plot.pdf",width=4.52,height=5.19)
dev.off()
