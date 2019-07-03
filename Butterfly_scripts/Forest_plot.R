#############################################################################
## Title: Create coefficient and CI plots for synchrony manuscript
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: November 2018
#############################################################################

rm(list=ls()) # clear R
library(plyr)
library(dplyr)
library(ggplot2)
library(plotrix)
library(lme4)
library(arm)
options(scipen=999)

############## Graph for average synchrony ###################
## read data
## butterflies
ukbms_strat <- read.csv( "../Results/Model_outputs/UKBMS/average_spec_ukbms.csv", header=TRUE)
ukbms_mob <- read.csv("../Results/Model_outputs/UKBMS/average_mob_ukbms.csv", header=TRUE)
ukbms_common <- read.csv("../Results/Model_outputs/UKBMS/average_abund_ukbms.csv", header=TRUE)
ukbms_sti <- read.csv( "../Results/Model_outputs/UKBMS/average_STIgroup_ukbms.csv", header=TRUE)
ukbms_fixed <- read.csv("../Results/Model_outputs/UKBMS/fixed_effect_results_ukbms.csv", header=TRUE)
## cbc birds
cbc_strat <- read.csv( "../Results/Model_outputs/CBC/average_spec_cbc.csv", header=TRUE)
cbc_mob <- read.csv("../Results/Model_outputs/CBC/average_mob_cbc.csv", header=TRUE)
cbc_common <- read.csv("../Results/Model_outputs/CBC/average_abund_cbc.csv", header=TRUE)
cbc_sti <- read.csv("../Results/Model_outputs/CBC/average_STIgroup_cbc.csv", header=TRUE)
cbc_fixed <- read.csv("../Results/Model_outputs/CBC/fixed_effect_results_cbc.csv", header=TRUE)
## bbs birds
bbs_strat <- read.csv( "../Results/Model_outputs/BBS/average_spec_bbs.csv", header=TRUE)
bbs_mob <- read.csv("../Results/Model_outputs/BBS/average_mob_bbs.csv", header=TRUE)
bbs_common <- read.csv("../Results/Model_outputs/BBS/average_abund_bbs.csv", header=TRUE)
bbs_sti <- read.csv("../Results/Model_outputs/BBS/average_STIgroup_bbs.csv", header=TRUE)
bbs_fixed <- read.csv("../Results/Model_outputs/BBS/fixed_effect_results_bbs.csv", header=TRUE)

### remove fixed & random effects from species attribute models
ukbms_strat <- ukbms_strat[-c(1:31),] 
ukbms_mob <- ukbms_mob[-c(1:31),] 
ukbms_common <- ukbms_common[-c(1:31),]
ukbms_sti <- ukbms_sti[-c(1:31),]
cbc_strat <- cbc_strat[-c(1:15),] 
cbc_mob <- cbc_mob[-c(1:5),] 
cbc_common <- cbc_common[-c(1:5),] 
cbc_sti <- cbc_sti[-c(1:15),]
bbs_strat <- bbs_strat[-c(1:17),] 
bbs_mob <- bbs_mob[-c(1:17),] 
bbs_common <- bbs_common[-c(1:17),] 
bbs_sti <- bbs_sti[-c(1:17),]
bbs_fixed <- bbs_fixed[-c(1,5:17),]

## create column for scheme name
ukbms_strat$Scheme <- "UKBMS"
ukbms_mob$Scheme <- "UKBMS"
ukbms_common$Scheme <- "UKBMS"
ukbms_sti$Scheme <- "UKBMS"
ukbms_fixed$Scheme <- "UKBMS"
cbc_strat$Scheme <- "CBC"
cbc_mob$Scheme <- "CBC"
cbc_common$Scheme <- "CBC"
cbc_sti$Scheme <- "CBC"
cbc_fixed$Scheme <- "CBC"
bbs_strat$Scheme <- "BBS"
bbs_mob$Scheme <- "BBS"
bbs_common$Scheme <- "BBS"
bbs_sti$Scheme <- "BBS"
bbs_fixed$Scheme <- "BBS"

################### correct way to plot graph #####################
## merge datasets
names(ukbms_fixed)[6] <- "X"
ukbms_fixed <- ukbms_fixed[,c(6,1,2,3,4,5,7)]
names(cbc_fixed)[6] <- "X"
cbc_fixed <- cbc_fixed[,c(6,1,2,3,4,5,7)]

average_synchrony <- rbind(ukbms_strat, ukbms_common, ukbms_mob, ukbms_sti, ukbms_fixed, cbc_strat, cbc_common, 
                           cbc_mob, cbc_sti, cbc_fixed, bbs_strat, bbs_common, bbs_mob, bbs_sti, bbs_fixed)
colnames(average_synchrony)[1] <- "vars"
rownames(average_synchrony) <- 1:nrow(average_synchrony)

### Change viarable names
lookup <- c("specialismwider.countryside"="Specialism", "specialismgeneralist"="Specialism", "specialismspecialist"="Specialism",
            "pop_est_stand"="Abundance", "pop_estimate_log"="Abundance", "pop_estimate_log"="Abundance",
            "mobility_score2high"="Mobility", "Breeding_AM_score2high"="Mobility", "Breeding_AM_score22"="Mobility",
            "mean_northing"="Mean northing", "distance"="Distance", "renk_hab_sim"="Habitat similarity", "hab_sim"="Habitat similarity",
            "specialismspecialist"="Specialism", "hab_sim1"="Habitat similarity", "STI_score2"="STI")
average_synchrony <- average_synchrony %>% mutate(vars=lookup[as.character(vars)])

## order alphabetically
average_synchrony <- average_synchrony[order(average_synchrony$vars),]
## calculate CI's from SE values
average_synchrony$CI <- average_synchrony$Std..Error*1.96

average_synchrony$Scheme <- as.factor(average_synchrony$Scheme)
levels(average_synchrony$Scheme)
average_synchrony$Scheme <- factor(average_synchrony$Scheme, levels=c("UKBMS", "CBC", "BBS"))

######## FINAL PLOT!!! ########
png("../Graphs/FINAL/Figure3_2.png", height = 100, width = 173, units = "mm", res = 300)
ggplot(data=average_synchrony,aes(x=vars,y=Estimate, group=Scheme),position=position_dodge(width=0.5)) +
  geom_point(aes(shape=Scheme), size=3,position=position_dodge(width=0.5)) +
  scale_x_discrete(limits=c("Specialism", "Mobility", "Abundance", "STI", "Mean northing",  "Distance","Habitat similarity")) +
  geom_errorbar(aes(ymin=Estimate-CI,ymax=Estimate+CI), position=position_dodge(width=0.5),width=0.2) +
  geom_hline(yintercept=0,linetype="dashed") +
  scale_y_continuous(breaks=seq(-0.1,0.2,0.05)) +
  labs(x = "Fixed effects", y = "Standardised coefficient") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

## also save as pdf
figure3 <- ggplot(data=average_synchrony,aes(x=vars,y=Estimate, group=Scheme),position=position_dodge(width=0.5)) +
  geom_point(aes(shape=Scheme), size=3,position=position_dodge(width=0.5)) +
  scale_x_discrete(limits=c("Specialism", "Mobility", "Abundance", "STI", "Mean northing",  "Distance","Habitat similarity")) +
  geom_errorbar(aes(ymin=Estimate-CI,ymax=Estimate+CI), position=position_dodge(width=0.5),width=0.2) +
  geom_hline(yintercept=0,linetype="dashed") +
  scale_y_continuous(breaks=seq(-0.1,0.2,0.05)) +
  labs(x = "Fixed effects", y = "Standardised coefficient") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave(filename="../Graphs/FINAL/Figure3_2.pdf", plot=figure3, width=173, height=100, units="mm")


############# Plot for BES presentation ############# 

## remove northing, distance and hab sim
average_synchrony2 <- average_synchrony[!grepl("Habitat similarity", average_synchrony$vars),]
average_synchrony2 <- average_synchrony2[!grepl("Mean northing", average_synchrony2$vars),]
average_synchrony2 <- average_synchrony2[!grepl("Distance", average_synchrony2$vars),]

### plot graph again with transparent background for powerpoint
png("../Graphs/FINAL/Presentation1.png", width = 4 * 800,
    height = 4 * 800, res = 600)

ggplot(data=average_synchrony2,aes(x=vars,y=Estimate, group=Scheme),position=position_dodge(width=0.5)) +
  geom_point(aes(shape=Scheme), size=6,position=position_dodge(width=0.5), colour="white") +
  scale_x_discrete(limits=c("Specialism", "Mobility", "Abundance")) +
  geom_errorbar(aes(ymin=Estimate-CI,ymax=Estimate+CI), colour="white", position=position_dodge(width=0.5),width=0.4) +
  geom_hline(yintercept=0,linetype="dashed", colour="white") +
  scale_y_continuous(breaks=seq(-0.1,0.2,0.05)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        legend.position="none",
        legend.key = element_blank(),
        axis.line = element_line(colour = "white"),
        axis.ticks = element_blank(),
        axis.text.y = element_text(colour="white", size=16),
        axis.text.x = element_text(colour="white", size=16),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA))
ggsave("../Graphs/FINAL/Presentation1.png", bg = "transparent")
dev.off()