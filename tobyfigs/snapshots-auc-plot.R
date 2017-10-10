library(ggplot2)
library(ggthemes) ##
library(tidyverse)
#library(plyr)
library(reshape2)
library(colorspace)

auc_data <- read.csv("snapshots-auc.csv")
auc_data$`Infectious period` <- as.factor(auc_data$Infectious.period)
levels(auc_data$`Infectious period`) <- c(expression(paste(gamma, "=1/30")),expression(paste(gamma, "=1/7")))

auc_data$`Aggregation period` <- as.factor(auc_data$Aggregation.period)
levels(auc_data$`Aggregation period`) <- c("Monthly~snapshots","Weekly~snapshots")

auc_data$variable <- factor(auc_data$variable, levels = levels(auc_data$variable)[c(7,10,9,5,2,1,4,3,8,6)])


auc_plot <- function(df){
  
  ## Eric's AUC gradient ##
  # Given error of ~0.5, set number of levels to 10.  For a higher precision, use nlevels = 20 
  nlevels = 100  # number of bins = number of colors
  # Color scale with {colorspace}
  AUC_colors <- diverge_hcl(
    n=nlevels,  # number of colors
    h = c(45, 225),  # Hues (low, hi)
    c = 100, # fixed Chroma or Chroma Range (edges, center)
    l = c(90, 10),  # Lightness range (edges, center) 
    power = 1  # exponent
  )
  names(AUC_colors) <- seq(0,1,1/(nlevels-1))
  
  
  ggplot(df) + 
    # geom_pointrange(aes(x=variable,y=AUC, ymin=AUC-AUC_err, ymax= AUC+AUC_err, color=variable)) +
    geom_bar(aes(x=variable,y=AUC-0.5, fill = AUC, color = AUC),stat="identity"  ) + #, fill = "#088E7C") +
    facet_grid(`Infectious period`~`Aggregation period`,labeller = label_parsed)+
    geom_rangeframe(colour ="black") +
    scale_fill_gradientn(limits = c(0,1),colours=AUC_colors) +
    scale_color_gradientn(limits = c(0,1),colours=AUC_colors) +
    scale_y_continuous(name = "AUC",limits = c(-0.5,0.5),labels=c("0.0","0.25","0.5","0.75","1.0"))
  #scale_x_continuous(name = "")+#,labels=c("Variance","Variance convexity","Autocovariance",
  #                                      "Autocorrelation","Decay time", "Mean",
  #                                      "Index of Dispersion", "Coefficient of variation",
  #                                      "Skewness", "Kurtosis"))+
  # theme_minimal() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

par(font.lab = 3)
text_color <- "black"
background_color <- "white"
font_chosen <- "Helvetica"
auc_plot(auc_data) +
  theme(text = element_text(color=text_color, size = 18, family=font_chosen),
        title = element_text(size = 14),
        line = element_line(color=text_color),
        rect = element_rect(color="black"),
        axis.title.x=element_blank(),
        axis.title.y = element_text(color=text_color,face = "plain", 
                                    family=font_chosen, size = 14),
        axis.text.x = element_text(color = text_color, angle = 45,
                                   vjust = 1, hjust=1),#family=font_chosen, size = 14),
        axis.text.y = element_text(color = text_color, family=font_chosen),
        axis.ticks = element_line(color=rgb(0,0,0,.25),line),
        axis.line = element_blank(),
        legend.text = element_text(color = text_color, family=font_chosen),
        legend.title = element_text(color = text_color, family = font_chosen),
        panel.grid.major = element_blank(),
        panel.grid.major.y = element_line(color = rgb(0,0,0,.25)),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(color=NA,fill=NA),
        panel.spacing = unit(2, "lines"),
        strip.text.y= element_text(color = text_color,angle = 270, family=font_chosen),
        strip.text.x= element_text(color = text_color, family=font_chosen),
        strip.background = element_rect(fill=NA),
        plot.background = element_rect(color=NA, fill=background_color),
        plot.margin = unit(c(1,1,1,1), "cm"))


ggsave("snapshots-auc.pdf",width =1.2*8.64,height=1.2*4.20)


