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

## input data ##
results_final_all_spp <- read.csv("../Results/Butterfly_results/results_final_all_spp_no_zeros2.csv", header=TRUE)
results_final_sp <- read.csv("../Results/Butterfly_results/results_final_sp_no_zeros2.csv", header=TRUE)
results_final_spec <- read.csv("../Results/Butterfly_results/results_final_spec_no_zeros2.csv", header=TRUE)
## climate data
results_table_winter_rain <- read.csv("../Results/Climate_results/winter_rainfall_synchrony.csv", header=TRUE)
results_table_spring_rain <- read.csv("../Results/Climate_results/spring_rainfall_synchrony.csv", header=TRUE)
results_table_summer_rain <- read.csv("../Results/Climate_results/summer_rainfall_synchrony.csv", header=TRUE)
results_table_autumn_rain <- read.csv("../Results/Climate_results/autumn_rainfall_synchrony.csv", header=TRUE)
results_table_winter_temp <- read.csv("../Results/Climate_results/winter_temp_synchrony.csv", header=TRUE)
results_table_spring_temp <- read.csv("../Results/Climate_results/spring_temp_synchrony.csv", header=TRUE)
results_table_summer_temp <- read.csv("../Results/Climate_results/summer_temp_synchrony.csv", header=TRUE)
results_table_autumn_temp <- read.csv("../Results/Climate_results/autumn_temp_synchrony.csv", header=TRUE)

## read in strategy data
spp_data <- read.csv("../Data/UKBMS_data/UKBMS_UKspecieslist.csv", header=TRUE)
spp_data <- spp_data[c(1,5)]
results_final_sp <- merge(results_final_sp, spp_data, by.x="sp", by.y="BMSCODE") ## merge to get strategy info 

#################################
## PLOT GRAPH FOR ALL SPECIES ##
#################################

## FCI plot WITHOUT scaling to 100 with smoothed line and SD error bars ##
FCI_plot_errorbars <- ggplot(results_final_all_spp, aes(x = parameter, y = FCI)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = FCI - SD, ymax = FCI + SD), width = .4) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(-0.2,0.4,0.05)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_plot_errorbars

ggsave("../Graphs/Connectivity_plots/FCI_plot_all_spp_errorbars_nomigrants.png", plot = FCI_plot_errorbars, width=7, height=5)

## FCI plot without scaled to 100 with smoothed line and 95% CI shaded error ## 
FCI_plot_unscaled <- ggplot(results_final_all_spp, aes(x = parameter, y = FCI)) +
  stat_smooth(colour="black", method=loess) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(0.1,0.36,0.02)) +
  scale_x_continuous(breaks=seq(1985,2012,2)) +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_plot_unscaled

ggsave("../Graphs/Connectivity_plots/FCI_plot_all_spp.png", plot = FCI_plot_unscaled, width=7, height=5)

## FCI plot scaled to 100 with smoothed line and CI error bars ##
FCI_plot_scaled <- ggplot(results_final_all_spp, aes(x = parameter, y = rescaled_FCI)) +
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
ggsave("../Graphs/Connectivity_plots/FCI_plot_all_spp_scaled3.png", plot = FCI_plot_scaled, width=7, height=5)
#ggsave("../Graphs/Presentation/FCI_UKBMS.png", plot = FCI_plot_scaled, width=7, height=5)

### add season to datasets
results_table_autumn_rain$category <- "autumn"
results_table_spring_rain$category <- "spring"
results_table_summer_rain$category <- "summer"
results_table_winter_rain$category <- "winter"

results_table_autumn_temp$category <- "autumn"
results_table_spring_temp$category <- "spring"
results_table_summer_temp$category <- "summer"
results_table_winter_temp$category <- "winter"

results_final_all_spp$category <- "population"
names(results_final_all_spp)[1] <- "synchrony"
names(results_final_all_spp)[5] <- "rescaled_sync"

## merge files
rainfall_synchrony <- rbind(results_table_autumn_rain, results_table_spring_rain, results_table_summer_rain, results_table_winter_rain)
temp_synchrony <- rbind(results_table_autumn_temp, results_table_spring_temp, results_table_summer_temp, results_table_winter_temp)
rainfall_synchrony$category <- factor(rainfall_synchrony$category, levels=c("spring", "summer", "autumn", "winter"))
temp_synchrony$category <- factor(temp_synchrony$category, levels=c("spring", "summer", "autumn", "winter"))

## plot rainfall synchrony
rainfall_sync_plot <- ggplot(rainfall_synchrony, aes(x = parameter, y = rescaled_sync, group=category)) +
  stat_smooth(aes(colour=category), method=loess, se=FALSE) +
  geom_errorbar(aes(colour=category, ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5, position=myjit) +
  geom_point(size=2, aes(colour=category), position=myjit) + 
  labs(x = "Mid-year of moving window", y = "Rainfall synchrony") +
  scale_y_continuous(breaks=seq(90,110,0.01)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  scale_color_manual(values=c("tomato1", "yellowgreen", "cyan3", "medium purple")) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
rainfall_sync_plot
ggsave("../Graphs/Climate/rainfall_synchrony.png", plot = rainfall_sync_plot, width=10, height=8)

## plot temperature synchrony
temp_sync_plot <- ggplot(temp_synchrony, aes(x = parameter, y = rescaled_sync, group=category)) +
  stat_smooth(aes(colour=category), method=loess, se=FALSE) +
  geom_errorbar(aes(colour=category, ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5, position=myjit) +
  geom_point(size=2, aes(colour=category), position=myjit) + 
  labs(x = "Mid-year of moving window", y = "Temperture synchrony") +
  scale_y_continuous(breaks=seq(90,110,0.01)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  scale_color_manual(values=c("tomato1", "yellowgreen", "cyan3", "medium purple")) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
temp_sync_plot
ggsave("../Graphs/Climate/temperature_synchrony.png", plot = temp_sync_plot, width=10, height=8)

## combine population and temperature synchrony data
pop_temp_sync <- rbind(results_final_all_spp, temp_synchrony)
levels(as.factor(pop_temp_sync$category))
pop_temp_sync$category <- factor(pop_temp_sync$category, levels=c("population", "spring", "summer", "autumn", "winter"))
levels(pop_temp_sync$category)

## plot graph
pop_temp_plot <- ggplot(pop_temp_sync, aes(x = parameter, y = rescaled_sync, group=category)) +
  stat_smooth(aes(colour=category), method=loess, se=FALSE) +
  geom_errorbar(aes(colour=category, ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5, position=myjit) +
  geom_point(size=2, aes(colour=category), position=myjit) + 
  labs(x = "Mid-year of moving window", y = "Synchrony") +
  scale_y_continuous(breaks=seq(20,120,10)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  scale_color_manual(values=c("black", "tomato1", "yellowgreen", "cyan3", "medium purple")) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
pop_temp_plot
ggsave("../Graphs/Climate/population_temperature_synchrony2.png", plot = pop_temp_plot, width=10, height=8)

pop_rain_sync <- rbind(results_final_all_spp, rainfall_synchrony)
levels(as.factor(pop_rain_sync$category))
pop_rain_sync$category <- factor(pop_rain_sync$category, levels=c("population", "spring", "summer", "autumn", "winter"))
levels(pop_rain_sync$category)

## plot graph
pop_rain_plot <- ggplot(pop_rain_sync, aes(x = parameter, y = rescaled_sync, group=category)) +
  stat_smooth(aes(colour=category), method=loess, se=FALSE) +
  geom_errorbar(aes(colour=category, ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5, position=myjit) +
  geom_point(size=2, aes(colour=category), position=myjit) + 
  labs(x = "Mid-year of moving window", y = "Synchrony") +
  scale_y_continuous(breaks=seq(20,120,10)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  scale_color_manual(values=c("black", "tomato1", "yellowgreen", "cyan3", "medium purple")) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
pop_rain_plot
ggsave("../Graphs/Climate/population_rainfall_synchrony2.png", plot = pop_rain_plot, width=10, height=8)

############ jitter code #################
myjit <- ggproto("fixJitter", PositionDodge,
                 width = 0.5,
                 dodge.width = 0.15,
                 jit = NULL,
                 compute_panel =  function (self, data, params, scales) 
                 {
                   
                   #Generate Jitter if not yet
                   if(is.null(self$jit) ) {
                     self$jit <-jitter(rep(0, nrow(data)), amount=self$dodge.width)
                   }
                   
                   data <- ggproto_parent(PositionDodge, self)$compute_panel(data, params, scales)
                   
                   data$x <- data$x + self$jit
                   #For proper error extensions
                   if("xmin" %in% colnames(data)) data$xmin <- data$xmin + self$jit
                   if("xmax" %in% colnames(data)) data$xmax <- data$xmax + self$jit
                   data
                 } )




############### plot with all 3 temporal synchrony graphs for ms #########################
results_final_all_spp_BBS <- read.csv("../Results/Bird_results/results_final_all_spp_BBS_final.csv", header=TRUE)
results_final_all_spp_CBC <- read.csv("../Results/Bird_results/results_final_all_spp_CBC_no_zeros2_correct.csv", header=TRUE)

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

library(ggpubr)

png("../Graphs/FINAL/FigureS1_2.png", height = 250, width = 220, units = "mm", res = 300)
ggarrange(FCI_plot_scaled,                                                 
      ggarrange(FCI_BBS, FCI_CBC, ncol = 2, labels = c("(b)", "(c)"), font.label = list(size = 12, color ="black")), 
      nrow = 2, labels = "(a)", font.label = list(size = 10, color ="black")) 
dev.off()

####### also add in percentage change bar graphs for ms ########
UKBMS_percentages <- read.csv("../Results/Butterfly_results/model_comp_percentages_no_zeros2.csv", header=TRUE)
CBC_percentages <- read.csv("../Results/Bird_results/cbc_overall_trends_correct.csv", header=TRUE)
BBS_percentages <- read.csv("../Results/Bird_results/bbs_overall_trends.csv", header=TRUE)

CBC_percentages$comparison <- "1985-1996"
BBS_percentages$comparison <- "1999-2012"

model_comp_birds <- rbind(BBS_percentages, CBC_percentages)
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

#### put all 5 graphs together somehow....
library(ggpubr)
library(gridExtra)

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

png("../Graphs/FINAL/Figure2_final.png", height = 150, width = 220, units = "mm", res = 300)
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
 
############# plot of woodland species (one line) FCI only #################
woodland_sp <- ggplot(subset(results_final_hab, habitat %in% ("Woodland")), aes(x = parameter, y = FCI)) +
  stat_smooth(colour="black", method=loess) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(0.05,0.4,0.02)) +
  scale_x_continuous(breaks=seq(1985,2012,2)) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size = 10))
woodland_sp

ggsave("../Graphs/Connectivity_plots/Woodland_FCI_spp.png", plot = woodland_sp, width=7, height=5)

############# Woodland species with error bars and unscaled #################
FCI_woodland_errorbars <- ggplot(subset(results_final_hab, habitat %in% ("Woodland")), aes(x = parameter, y = FCI)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = FCI - SD, ymax = FCI + SD), width = .4) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(0.05,0.4,0.02)) +
  scale_x_continuous(breaks=seq(1985,2012,2)) +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_woodland_errorbars

ggsave("../Graphs/Connectivity_plots/FCI_woodland_errorbars.png", plot = FCI_woodland_errorbars, width=7, height=5)


########### woodland species with error bars and scaled to 100 ################
FCI_woodland_scaled <- ggplot(results_final_all_spp_wood, aes(x = parameter, y = rescaled_FCI)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_ci, ymax = rescaled_FCI + rescaled_ci), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(20,220,20)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_woodland_scaled

ggsave("../Graphs/Connectivity_plots/FCI_woodland_scaled.png", plot = FCI_woodland_scaled, width=7, height=5)


############################
## PLOT GRAPHS BY SPECIES ##
############################

#### individual FCI plot for presentation
FCI_plot_93 <- ggplot(results_final_sp, aes(x = parameter, y = rescaled_FCI)) +
  stat_smooth(data=subset(results_final_sp, sp=="93"), colour="black", method=loess, se=FALSE) +
  geom_errorbar(data=subset(results_final_sp, sp=="93"), aes(ymin = rescaled_FCI - rescaled_ci, ymax = rescaled_FCI + rescaled_ci), width=0.2, size = 0.5) +
  geom_point(data=subset(results_final_sp, sp=="93"), size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(40,250,20)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_plot_93
## species 23 = speckled wood (wider countryside woodland species)

ggsave("../Graphs/Presentation/Spp93_FCI_plot.png", plot = FCI_plot_93, width=7, height=5)

FCI_plot_98 <- ggplot(results_final_sp, aes(x = parameter, y = rescaled_FCI)) +
  stat_smooth(data=subset(results_final_sp, sp=="98"), colour="black", method=loess, se=FALSE) +
  geom_errorbar(data=subset(results_final_sp, sp=="98"), aes(ymin = rescaled_FCI - rescaled_ci, ymax = rescaled_FCI + rescaled_ci), width=0.2, size = 0.5) +
  geom_point(data=subset(results_final_sp, sp=="98"), size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(40,200,20)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_plot_98
## species 98 = large white (wider countryside garden and hedgerow species)

ggsave("../Graphs/Presentation/Spp98_FCI_plot.png", plot = FCI_plot_98, width=7, height=5)

#########################################################################################################
#########################################################################################################

### FCI for each individual species unscaled
FCI_35_spp_unscaled <- ggplot(results_final_sp, aes(x = parameter, y = FCI, group = COMMON_NAME)) +
  geom_line(aes(colour=COMMON_NAME), lwd=0.5) +
  geom_point(size=1, aes(colour=COMMON_NAME)) + 
  labs(x = "Year", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(-1,1,0.1)) +
  scale_x_continuous(breaks=seq(1985,2012,2)) +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(colour = "Species", size=2) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_35_spp_unscaled

ggsave("../Graphs/Connectivity_plots/FCI_plot_35spp_unscaled.png", plot = FCI_35_spp_unscaled, width=7, height=5)

## CI for each individual species scaled to 100
FCI_35_spp_scaled <- ggplot(results_final_sp, aes(x = parameter, y = rescaled_FCI, group = COMMON_NAME)) +
  geom_line(aes(colour=COMMON_NAME), lwd=0.5) +
  geom_point(size=1, aes(colour=COMMON_NAME)) + 
  labs(x = "Year", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(-200,10000,100)) +
  scale_x_continuous(breaks=seq(1985,2012,2)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(colour = "Species", size=2)
FCI_35_spp_scaled

ggsave("../Graphs/Connectivity_plots/FCI_35_spp_scaled.png", plot = FCI_35_spp_scaled, width=10, height=8)

## FCI plot coloured by habitat type and scaled to 100 ##
SD_hab_plot <- ggplot(results_final_sp, aes(x = parameter, y = rescaled_FCI, group = sp)) +
  geom_line(aes(colour=HABITAT), lwd=0.5) +
  geom_point(size=1, aes(colour=HABITAT)) + 
  labs(x = "Year", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(-200,10000,100)) +
  scale_x_continuous(breaks=seq(1985,2012,2)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 10))
SD_hab_plot

ggsave("../Graphs/Habitat_plots/FCI_habitat_plot_35spp.png", plot = SD_hab_plot, width=7, height=5)

## FCI plot coloured by habitat type and unscaled ##
SD_hab_plot_unscaled <- ggplot(results_final_sp, aes(x = parameter, y = FCI, group = sp)) +
  geom_line(aes(colour=HABITAT), lwd=0.5) +
  geom_point(size=1, aes(colour=HABITAT)) + 
  labs(x = "Year", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(-0.5,1,0.2)) +
  scale_x_continuous(breaks=seq(1985,2012,2)) +
  theme_bw() +
  theme(text = element_text(size = 10))
SD_hab_plot_unscaled

ggsave("../Graphs/Habitat_plots/FCI_habitat_plot_unscaled_35spp.png", plot = SD_hab_plot_unscaled, width=7, height=5)

## PLOT WOODLAND SPECIES GRAPH ##
FCI_woodland <- ggplot(subset(results_final_sp, HABITAT %in% ("Woodland")), aes(x = parameter, y = FCI, group = sp)) +
  stat_smooth(aes(colour=COMMON_NAME), method=loess, level=0.95, se=FALSE) +
  geom_point(aes(colour=COMMON_NAME),size=1) + 
  labs(x = "Year", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(-0.2,0.4,0.1)) +
  scale_x_continuous(breaks=seq(1985,2012,2)) +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(colour = "Woodland Species", size=2)
FCI_woodland

ggsave("../Graphs/Habitat_plots/FCI_woodland_plot.png", plot = FCI_woodland, width=7, height=5)

## PLOT WOODLAND SPECIES GRAPH SCALED TO 100 ##
FCI_woodland_scaled <- ggplot(subset(results_final_sp, HABITAT %in% ("Woodland")), aes(x = parameter, y = rescaled_FCI, group = sp)) +
  stat_smooth(aes(colour=COMMON_NAME), method=loess, level=0.95, se=FALSE) +
  geom_point(aes(colour=COMMON_NAME),size=1) + 
  labs(x = "Year", y = "Functional connectivity index") +
  #scale_y_continuous(breaks=seq(-1500,1000,200)) +
  scale_x_continuous(breaks=seq(1985,2012,2)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(colour = "Woodland Species", size=2)
FCI_woodland_scaled

ggsave("../Graphs/Habitat_plots/FCI_woodland_plot_sacled.png", plot = FCI_woodland_scaled, width=7, height=5)

## PLOT GRASSLAND SPECIES GRAPH ##
FCI_grassland <- ggplot(subset(results_final_sp, HABITAT %in% ("Grassland")), aes(x = parameter, y = rescaled_FCI, group = sp)) +
  geom_line(aes(color=COMMON_NAME)) +
  geom_point(size=1, aes(color=COMMON_NAME)) + 
  labs(x = "Year", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(-100,1100,100)) +
  scale_x_continuous(breaks=seq(1985,2012,2)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(colour = "Grassland Species", size=2)
FCI_grassland

ggsave("../Graphs/Habitat_plots/FCI_grassland_plot.png", plot = FCI_grassland, width=7, height=5)

## PLOT GARDEN AND HEDGEROW SPECIES GRAPH ##
FCI_garden <- ggplot(subset(results_final_sp, HABITAT %in% ("Garden and hedgerow")), aes(x = parameter, y = rescaled_FCI, group = sp)) +
  geom_line(aes(color=COMMON_NAME)) +
  geom_point(size=1, aes(color=COMMON_NAME)) + 
  labs(x = "Year", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(-100,400,50)) +
  scale_x_continuous(breaks=seq(1985,2012,2)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(colour = "Garden and hedgerow Species", size=2)
FCI_garden

ggsave("../Graphs/Habitat_plots/FCI_garden_plot.png", plot = FCI_garden, width=7, height=5)

## PLOT HEATHLAND SPECIES GRAPH ##
FCI_heath <- ggplot(subset(results_final_sp, HABITAT %in% ("Heathland")), aes(x = parameter, y = rescaled_FCI, group = sp)) +
  geom_line(aes(color=COMMON_NAME)) +
  geom_point(size=1, aes(color=COMMON_NAME)) + 
  labs(x = "Year", y = "Functional connectivity index") +
  scale_y_continuous(breaks=c(20,100,10)) +
  scale_x_continuous(breaks=seq(1985,2012,2)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(colour = "Heathland Species", size=2)
FCI_heath

ggsave("../Graphs/Habitat_plots/FCI_heath_plot.png", plot = FCI_heath, width=7, height=5)

##### plot woodland species by strategy ######
## subset to get just woodland species

## plot wider countryside species
FCI_wider_country <- ggplot(subset(results_final_sp_woodland, STRATEGY %in% ("Wider countryside sp")), aes(x = parameter, y = rescaled_FCI, group = sp)) +
  stat_smooth(aes(color=COMMON_NAME), method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_sd, ymax = rescaled_FCI + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=1, color="black") + 
  labs(x = "Year", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(0,210,20)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(colour = "Wider countryside species", size=2)
FCI_wider_country

ggsave("../Graphs/Habitat_plots/FCI_wider_country.png", plot = FCI_wider_country, width=7, height=5)


FCI_specialist <- ggplot(subset(results_final_sp_woodland, STRATEGY %in% ("Habitat specialist")), aes(x = parameter, y = rescaled_FCI, group = sp)) +
  stat_smooth(aes(color=COMMON_NAME), method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_sd, ymax = rescaled_FCI + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=1, color="black") + 
  labs(x = "Year", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(-400,300,50)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(colour = "Habitat specialist species", size=2)
FCI_specialist

ggsave("../Graphs/Habitat_plots/FCI_specialist.png", plot = FCI_specialist, width=7, height=5)

############################
## PLOT GRAPHS BY HABITAT ##
############################

### CI plot with coloured shaded CI error 
FCI_CI_habitat_plot <- ggplot(results_final_hab, aes(x = parameter, y = rescaled_FCI, group = habitat)) +
  geom_ribbon(aes(ymin = rescaled_FCI - rescaled_ci, ymax = rescaled_FCI + rescaled_ci), alpha=0.2, lwd=0.5) +
  geom_line(aes(linetype=habitat)) +
  geom_point(size=1) + 
  labs(x = "Year", y = "Functional Connectivity index") +
  scale_y_continuous(breaks=seq(40,200,20)) +
  scale_x_continuous(breaks=seq(1985,2012,2)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(linetype = "Habitat", size=2)
FCI_CI_habitat_plot

ggsave("../Graphs/Habitat_plots/FCI_CI_habitat_plot.png", plot = FCI_CI_habitat_plot, width=7, height=5)

## main habitat plot
hab_plot <- ggplot(results_final_hab, aes(x = parameter, y = FCI, group = habitat)) +
  stat_smooth(aes(linetype=habitat), colour="black", method=loess) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(0.1,0.6,0.05)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  labs(aes(linetype = "Habitat Type"), size=3) +
  theme(legend.position = c(0.68,0.8)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
hab_plot

ggsave("../Graphs/Habitat_plots/FCI_habitat_plot.png", plot = hab_plot, width=7, height=5)

#### habitat plot with error bars and smoothed line sacled ####
FCI_habitat_error_scaled <- ggplot(results_final_hab, aes(x = parameter, y = rescaled_FCI, group = habitat)) +
  stat_smooth(aes(linetype=habitat), colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_ci, ymax = rescaled_FCI + rescaled_ci), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(20,220,20)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  theme(legend.position = c(0.7,0.8)) +
  labs(aes(linetype = "Habitat"), size=3) +
  scale_linetype_manual(values=c("solid", "twodash", "dotted")) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_habitat_error_scaled

ggsave("../Graphs/Habitat_plots/FCI_habitat_error_scaled.png", plot = FCI_habitat_error_scaled, width=7, height=5)

#############################
## PLOT GRAPHS BY STRATEGY ##
#############################

## take out regular migrants
results_final_strat$strategy <- as.character(results_final_strat$strategy)
results_final_strat <- results_final_strat[!results_final_strat$strategy=="Regular migrants    ",]

library(dplyr)
results_final_strat <- results_final_strat %>% mutate(parameter = ifelse(strategy=="Wider countryside", parameter+0.25, parameter))

strat_plot <- ggplot(results_final_strat, aes(x = parameter, y = FCI, group = strategy)) +
  stat_smooth(aes(linetype=strategy), colour="black", method=loess) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(0.1,0.6,0.05)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  labs(aes(linetype = "Strategy"), size=3) +
  theme(legend.position = c(0.8,0.6)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
strat_plot

ggsave("../Graphs/Connectivity_plots/FCI_strategy_plot.png", plot = strat_plot, width=7, height=5)

FCI_strat_error_scaled <- ggplot(results_final_strat, aes(x = parameter, y = rescaled_FCI, group = strategy)) +
  stat_smooth(aes(linetype=strategy), colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_ci, ymax = rescaled_FCI + rescaled_ci), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(20,200,20)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  labs(aes(linetype = "Strategy"), size=3) +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_strat_error_scaled

ggsave("../Graphs/Connectivity_plots/FCI_strat_error_scaled.png", plot = FCI_strat_error_scaled, width=7, height=5)

FCI_strat_error_unscaled <- ggplot(results_final_strat, aes(x = parameter, y = FCI, group = strategy)) +
  stat_smooth(aes(linetype=strategy), colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = FCI - SD, ymax = FCI + SD), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(0.1,0.8,0.1)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  labs(aes(linetype = "Strategy"), size=3) +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_strat_error_unscaled

ggsave("../Graphs/Connectivity_plots/FCI_strat_error_unscaled.png", plot = FCI_strat_error_unscaled, width=7, height=5)

##### same plots on woodland species only

## add 0.25 to year for wider countryside species to offset error bars
library(dplyr)
results_final_strat_wood <- results_final_strat_wood %>% mutate(parameter = ifelse(strategy=="Wider countryside", parameter+0.25, parameter))

FCI_strat_wood_error_scaled <- ggplot(results_final_strat_wood, aes(x = parameter, y = rescaled_FCI, group = strategy)) +
  stat_smooth(aes(linetype=strategy), colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_ci, ymax = rescaled_FCI + rescaled_ci), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(0,260,20)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  labs(aes(linetype = "Strategy"), size=3) +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  theme(legend.title = element_blank()) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_strat_wood_error_scaled

ggsave("../Graphs/Connectivity_plots/FCI_strat_error_scaled_wood.png", plot = FCI_strat_wood_error_scaled, width=7, height=5)

FCI_strat_wood_error_unscaled <- ggplot(results_final_strat_wood, aes(x = parameter, y = FCI, group = strategy)) +
  stat_smooth(aes(linetype=strategy), colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = FCI - SD, ymax = FCI + SD), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(0,0.8,0.05)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  labs(aes(linetype = "Strategy"), size=3) +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  theme(legend.position = c(0.6,0.87)) +
  theme(legend.title = element_blank()) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_strat_wood_error_unscaled

ggsave("../Graphs/Connectivity_plots/FCI_strat_error_unscaled_wood.png", plot = FCI_strat_wood_error_unscaled, width=7, height=5)


#############################
## PLOT GRAPHS BY MOBILITY ##
#############################

### CI plot with coloured shaded CI error 
CI_plot2 <- ggplot(results_final_mobility, aes(x = parameter, y = rescaled_FCI, group = mobility.score)) +
  geom_ribbon(aes(ymin = rescaled_FCI - rescaled_ci, ymax = rescaled_FCI + rescaled_ci), alpha=0.2, lwd=0.5) +
  geom_line(aes(linetype=mobility.score)) +
  geom_point(size=1) + 
  labs(x = "Year", y = "Connectivity index") +
  scale_y_continuous(breaks=seq(-50,150,25)) +
  scale_x_continuous(breaks=seq(1985,2012,2)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(linetype = "Mobility", size=2)
CI_plot2

ggsave("../Graphs/Mobility_plots/FCI_CI_mobility_plot.png", plot = CI_plot2, width=7, height=5)

## SD plot with smooth lines and NOT scaled to 100 ## 
SD_plot2 <- ggplot(results_final_mobility, aes(x = parameter, y = rescaled_FCI, group = mobility.score)) +
  geom_ribbon(aes(ymin = rescaled_FCI - rescaled_sd, ymax = rescaled_FCI + rescaled_sd), alpha=0.2, lwd=0.5) +
  geom_line(aes(linetype=mobility.score)) +
  geom_point(size=1) + 
  labs(x = "Year", y = "Connectivity index") +
  scale_y_continuous(breaks=seq(-50,150,25)) +
  scale_x_continuous(breaks=seq(1985,2012,2)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  labs(linetype = "Mobility", size=2)
SD_plot2

ggsave("../Graphs/Habitat_plots/FCI_SD_mobility_plot.png", plot = SD_plot2, width=7, height=5)
