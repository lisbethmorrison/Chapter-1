###########################################################
## Title: Moving window functional connectivity graphs 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: September 2017
##########################################################


#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script") ####

rm(list=ls()) # clear R

library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggpubr)
library(gridExtra)

## input data ##
results_final_all_spp <- read.csv("../Results/Butterfly_results/results_final_all_spp.csv", header=TRUE)
results_final_sp <- read.csv("../Results/Butterfly_results/results_final_sp.csv", header=TRUE)
results_final_spec <- read.csv("../Results/Butterfly_results/results_final_spec.csv", header=TRUE)
results_final_climate <- read.csv("../Results/Butterfly_results/results_final_all_spp_climate.csv", header=TRUE)

## read in specialism data
spp_data <- read.csv("../Data/UKBMS_data/UKBMS_UKspecieslist.csv", header=TRUE)
spp_data <- spp_data[c(1,5)]
results_final_sp <- merge(results_final_sp, spp_data, by.x="sp", by.y="BMSCODE") ## merge to get strategy info 

#################################
## PLOT GRAPH FOR ALL SPECIES ##
#################################


## FCI plot scaled to 100 with smoothed line and SE error bars ##
FCI_plot_scaled <- ggplot(results_final_climate, aes(x = parameter, y = rescaled_FCI)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_sd, ymax = rescaled_FCI + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Population synchrony") +
  #scale_y_continuous(breaks=seq(40,160,10)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 11)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
FCI_plot_scaled
ggsave("../Graphs/Connectivity_plots/FCI_plot_all_spp_scaled_climate.png", plot = FCI_plot_scaled, width=7, height=5)

############### plot with all 3 temporal synchrony graphs for ms #########################
results_final_all_spp_BBS <- read.csv("../Results/Bird_results/results_final_all_spp_BBS.csv", header=TRUE)
results_final_all_spp_CBC <- read.csv("../Results/Bird_results/results_final_all_spp_CBC.csv", header=TRUE)

FCI_BBS <- ggplot(results_final_all_spp_BBS, aes(x = parameter, y = rescaled_FCI)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_sd, ymax = rescaled_FCI + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Population synchrony") +
  scale_y_continuous(breaks=seq(-20,180,20)) +
  scale_x_continuous(breaks=seq(1999,2012,3)) +
  theme_bw() +
  theme(text = element_text(size = 11)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_BBS

FCI_CBC <- ggplot(results_final_all_spp_CBC, aes(x = parameter, y = rescaled_FCI)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_sd, ymax = rescaled_FCI + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Population synchrony") +
  scale_y_continuous(breaks=seq(0,200,20)) +
  scale_x_continuous(breaks=seq(1985,1996,3)) +
  theme_bw() +
  theme(text = element_text(size = 11)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_CBC

#######  add in percentage change bar graphs for ms ########
UKBMS_percentages <- read.csv("../Results/Butterfly_results/ukbms_percentage_trends.csv", header=TRUE)
CBC_percentages <- read.csv("../Results/Bird_results/cbc_percentage_trends.csv", header=TRUE)
BBS_percentages <- read.csv("../Results/Bird_results/bbs_percentage_trends.csv", header=TRUE)

CBC_percentages$comparison <- "1985-1996"
BBS_percentages$comparison <- "1999-2012"

## bind the two bird datasets together
model_comp_birds <- rbind(BBS_percentages, CBC_percentages)
## change levels 
model_comp_birds$comparison <- factor(model_comp_birds$comparison, levels=c("1985-1996", "1999-2012"))
model_comp_birds$change <- factor(model_comp_birds$change, levels=c("Increase", "No change", "Decrease"))
UKBMS_percentages$comparison <- factor(UKBMS_percentages$comparison, levels=c("1985-2000", "2000-2012", "1985-2012"))
UKBMS_percentages$change <- factor(UKBMS_percentages$change, levels=c("Increase", "No change", "Decrease"))

butterfly <- ggplot(data=UKBMS_percentages, aes(x=comparison, y=percentage, fill=change)) +
  geom_bar(stat="identity", width=0.7) +
  labs(y="Percentage of species", x="", fill="") +
  scale_fill_manual(values = c("#339900", "#999999", "#990000")) +
  theme_bw() +
  theme(axis.text.x=element_text(colour = "black")) +
  theme(axis.text.y=element_text(colour = "black")) +
  scale_y_continuous(breaks = seq(0,100,10), expand = c(0, 0)) +
  theme(legend.position="bottom",legend.direction="horizontal") +
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  theme(text = element_text(size = 11), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
butterfly
## without overall change
UKBMS_percentages <- UKBMS_percentages[!UKBMS_percentages$comparison == "1985-2012",]
butterfly2 <- ggplot(data=UKBMS_percentages, aes(x=comparison, y=percentage, fill=change)) +
  geom_bar(stat="identity", width=0.7) +
  labs(y="Percentage of species", x="", fill="") +
  scale_fill_manual(values = c("#339900", "#999999", "#990000")) +
  theme_bw() +
  theme(axis.text.x=element_text(colour = "black")) +
  theme(axis.text.y=element_text(colour = "black")) +
  scale_y_continuous(breaks = seq(0,100,10), expand = c(0, 0)) +
  theme(legend.position="bottom",legend.direction="horizontal") +
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  theme(text = element_text(size = 11), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
butterfly2

birds <- ggplot(data=model_comp_birds, aes(x=comparison, y=percentage, fill=change)) +
  geom_bar(stat="identity", width=0.7) +
  labs(y="Percentage of species", x="", fill="") +
  scale_fill_manual(values = c("#339900", "#999999", "#990000")) +
  theme_bw() +
  theme(axis.text.x=element_text(colour = "black")) +
  theme(axis.text.y=element_text(colour = "black")) +
  guides(fill = guide_legend(override.aes = list(size = 6))) +
  scale_y_continuous(breaks = seq(0,100,10), expand = c(0, 0)) +
  theme(legend.position="top") +
  theme(text = element_text(size = 11), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
birds

#### put all 5 graphs together
## create percentage plots with shared legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(butterfly)
butterfly2 <- butterfly2 + theme(legend.position="none")
birds <- birds + theme(legend.position="none")

percentages <- grid.arrange(butterfly2, legend, birds, ncol=1, nrow = 3, widths = c(2.7), heights = c(2.5, 0.4, 2.5))

png("../Graphs/FINAL/Figure2.png", height = 150, width = 220, units = "mm", res = 300)
grid.arrange(FCI_plot_scaled, percentages,
             FCI_CBC, FCI_BBS, 
             layout_matrix=cbind(c(1,3), c(1,4), c(2,2)), nrow=2)
dev.off()

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
write.csv(FCI_indicator_plot, file="../Connectivity fiche/Data files/UKBMS/C2i data.csv", row.names=FALSE)
 

# #######################################
# ## PLOT SOME GRAPHS FOR PRESENTATION ##
# #######################################
# 
# #### individual FCI plot for presentation
# FCI_plot_93 <- ggplot(results_final_sp, aes(x = parameter, y = rescaled_FCI)) +
#   stat_smooth(data=subset(results_final_sp, sp=="93"), colour="black", method=loess, se=FALSE) +
#   geom_errorbar(data=subset(results_final_sp, sp=="93"), aes(ymin = rescaled_FCI - rescaled_ci, ymax = rescaled_FCI + rescaled_ci), width=0.2, size = 0.5) +
#   geom_point(data=subset(results_final_sp, sp=="93"), size=2) + 
#   labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
#   scale_y_continuous(breaks=seq(40,250,20)) +
#   scale_x_continuous(breaks=seq(1985,2012,3)) +
#   geom_hline(yintercept = 100, linetype = "dashed") +
#   theme_bw() +
#   theme(text = element_text(size = 12)) +
#   labs(size=3) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# FCI_plot_93
# ## species 23 = speckled wood (wider countryside woodland species)
# 
# ggsave("../Graphs/Presentation/Spp93_FCI_plot.png", plot = FCI_plot_93, width=7, height=5)
# 
# FCI_plot_98 <- ggplot(results_final_sp, aes(x = parameter, y = rescaled_FCI)) +
#   stat_smooth(data=subset(results_final_sp, sp=="98"), colour="black", method=loess, se=FALSE) +
#   geom_errorbar(data=subset(results_final_sp, sp=="98"), aes(ymin = rescaled_FCI - rescaled_ci, ymax = rescaled_FCI + rescaled_ci), width=0.2, size = 0.5) +
#   geom_point(data=subset(results_final_sp, sp=="98"), size=2) + 
#   labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
#   scale_y_continuous(breaks=seq(40,200,20)) +
#   scale_x_continuous(breaks=seq(1985,2012,3)) +
#   geom_hline(yintercept = 100, linetype = "dashed") +
#   theme_bw() +
#   theme(text = element_text(size = 12)) +
#   labs(size=3) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# FCI_plot_98
# ## species 98 = large white (wider countryside garden and hedgerow species)
# 
# ggsave("../Graphs/Presentation/Spp98_FCI_plot.png", plot = FCI_plot_98, width=7, height=5)

