###########################################################
## Title: Synchrony Model Graphs BTO BBS data  
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: May 2018
##########################################################

rm(list=ls()) # clear R

library(ggplot2)
library(dplyr)

####################################
########## FCI GRAPHS ##############
####################################

### read in results tables ###
# results_final_all_spp_BBS_1 <- read.csv("../Results/Bird_results/results_final_all_spp_BBS_1.csv", header=TRUE)
# results_final_all_spp_BBS_2 <- read.csv("../Results/Bird_results/results_final_all_spp_BBS_2.csv", header=TRUE)
# results_final_all_spp_BBS_3 <- read.csv("../Results/Bird_results/results_final_all_spp_BBS_3.csv", header=TRUE)
# results_final_all_spp_BBS_4 <- read.csv("../Results/Bird_results/results_final_all_spp_BBS_4.csv", header=TRUE)
# results_final_all_spp_BBS_5 <- read.csv("../Results/Bird_results/results_final_all_spp_BBS_5.csv", header=TRUE)
results_final_all_spp_BBS <- read.csv("../Results/Bird_results/results_final_all_spp_BBS.csv", header=TRUE)

## FCI plot with error bars and sacled to 100
FCI_plot_scaled_error <- ggplot(results_final_all_spp_BBS, aes(x = parameter, y = rescaled_FCI)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_sd, ymax = rescaled_FCI + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Population synchrony") +
  scale_y_continuous(breaks=seq(-20,180,20)) +
  scale_x_continuous(breaks=seq(1999,2012,3)) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
FCI_plot_scaled_error

ggsave("../Graphs/Connectivity_plots/FCI_all_spp_BBS_scaled_error.png", plot = FCI_plot_scaled_error, width=5, height=5)

##### method to obtain smoothed values with 95% confidence intervals for indicator
loess_model <- loess(rescaled_FCI~parameter, data=results_final_all_spp_BBS)
loess_predict <- predict(loess_model, se=T)
FCI_indicator_plot <- data.frame("fitted"=loess_predict$fit, "lwl"=(loess_predict$fit-1.96*loess_predict$se.fit),
                                 "upl"=(loess_predict$fit+1.96*loess_predict$se.fit))
### merge with results_final_all_spp file
FCI_indicator_plot <- cbind(FCI_indicator_plot, results_final_all_spp_BBS)
## remove columns not needed
FCI_indicator_plot <- FCI_indicator_plot[-c(4:6,9:10)]
## reorder columns
FCI_indicator_plot <- FCI_indicator_plot[c(4,5,1,3,2)]
names(FCI_indicator_plot) <-  c("year", "unsmoothed", "smoothed", "upperCI", "lowerCI")
## save file
write.csv(FCI_indicator_plot, file="../Connectivity fiche/Data files/BBS/BBS trend data.csv", row.names=FALSE)

#### plot with 5 runs of BBS 
results_final_all_spp_BBS_1$Run <- 1
results_final_all_spp_BBS_2$Run <- 2
results_final_all_spp_BBS_3$Run <- 3
results_final_all_spp_BBS_4$Run <- 4
results_final_all_spp_BBS_5$Run <- 5
results_final_all_spp_BBS <- rbind(results_final_all_spp_BBS_1, results_final_all_spp_BBS_2, results_final_all_spp_BBS_3, results_final_all_spp_BBS_4, results_final_all_spp_BBS_5)
## stagger error bars
results_final_all_spp_BBS <- results_final_all_spp_BBS %>% mutate(parameter = ifelse(Run=="2", parameter+0.125, parameter))
results_final_all_spp_BBS <- results_final_all_spp_BBS %>% mutate(parameter = ifelse(Run=="3", parameter+0.25, parameter))
results_final_all_spp_BBS <- results_final_all_spp_BBS %>% mutate(parameter = ifelse(Run=="4", parameter+0.375, parameter))
results_final_all_spp_BBS <- results_final_all_spp_BBS %>% mutate(parameter = ifelse(Run=="5", parameter+0.5, parameter))
results_final_all_spp_BBS$Run <- as.factor(results_final_all_spp_BBS$Run)

FCI_BBS_runs <- ggplot(results_final_all_spp_BBS, aes(x = parameter, y = rescaled_FCI)) +
  stat_smooth(aes(group=Run, color=Run), method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_ci, ymax = rescaled_FCI + rescaled_ci), width=0.2, size = 0.5) +
  geom_point(size=2, aes(color=Run)) + 
  labs(x = "Mid-year of moving window", y = "Population synchrony") +
  scale_y_continuous(breaks=seq(-20,180,20)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_color_brewer(palette="Set1")
FCI_BBS_runs

ggsave("../Graphs/Connectivity_plots/FCI_BBS_runs.png", plot = FCI_BBS_runs, width=7, height=5)
