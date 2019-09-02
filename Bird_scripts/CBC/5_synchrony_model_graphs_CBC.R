###########################################################
## Title: Synchrony Model Graphs BTO CBC data  
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: January 2018
##########################################################

rm(list=ls()) # clear R

library(ggplot2)
library(gridExtra)

####################################
########## FCI GRAPHS ##############
####################################

### read in results tables ###
# results_final_all_spp <- read.csv("../Results/Bird_results/results_final_all_spp_CBC.csv", header=TRUE)
# results_final_sp <- read.csv("../Results/Bird_results/results_final_sp_CBC.csv", header=TRUE)
# results_final_spec <- read.csv("../Results/Bird_results/results_final_spec_CBC.csv", header=TRUE) 
results_final_all_spp <- read.csv("../Results/Bird_results/results_final_all_spp_climate_CBC.csv", header=TRUE)

## FCI plot with error bars and sacled to 100
FCI_plot_scaled_error <- ggplot(results_final_all_spp, aes(x = parameter, y = rescaled_FCI)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_sd, ymax = rescaled_FCI + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Population synchrony") +
  #scale_y_continuous(breaks=seq(0,180,20)) +
  scale_x_continuous(breaks=seq(1985,1996,3)) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
FCI_plot_scaled_error

ggsave("../Graphs/Connectivity_plots/FCI_all_spp_CBC_scaled_error_climate.png", plot = FCI_plot_scaled_error, width=5, height=5)

##### method to obtain smoothed values with 95% confidence intervals for indicator
loess_model <- loess(rescaled_FCI~parameter, data=results_final_all_spp)
loess_predict <- predict(loess_model, se=T)
FCI_indicator_plot <- data.frame("fitted"=loess_predict$fit, "lwl"=(loess_predict$fit-1.96*loess_predict$se.fit),
                                 "upl"=(loess_predict$fit+1.96*loess_predict$se.fit))
### merge with results_final_all_spp file
FCI_indicator_plot <- cbind(FCI_indicator_plot, results_final_all_spp)
## remove columns not needed
FCI_indicator_plot <- FCI_indicator_plot[-c(4:6,9:10)]
## reorder columns
FCI_indicator_plot <- FCI_indicator_plot[c(4,5,1,3,2)]
names(FCI_indicator_plot) <-  c("year", "unsmoothed", "smoothed", "upperCI", "lowerCI")
## save file
write.csv(FCI_indicator_plot, file="../Connectivity fiche/Data files/CBC/CBC trend data.csv", row.names=FALSE)

