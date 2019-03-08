###########################################################
## Title: Synchrony Model Graphs BTO CBC data  
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: January 2018
##########################################################

rm(list=ls()) # clear R

library(ggplot2)

####################################
########## FCI GRAPHS ##############
####################################

### read in results tables ###
results_final_all_spp <- read.csv("../Results/Bird_results/results_final_all_spp_CBC.csv", header=TRUE)
results_final_sp <- read.csv("../Results/Bird_results/results_final_spp_CBC.csv", header=TRUE)
results_final_spec <- read.csv("../Results/Bird_results/results_final_spec_CBC.csv", header=TRUE) 

#### one line for all species with shaded smoothed line and unscaled ####
FCI_plot_unscaled <- ggplot(results_final_all_spp, aes(x = parameter, y = FCI)) +
  stat_smooth(colour="black", method=loess, level=0.95) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(0.02,0.06,0.01)) +
  scale_x_continuous(breaks=seq(1985,1996,1)) +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_plot_unscaled

ggsave("../Graphs/Connectivity_plots/FCI_all_spp_CBC_unscaled.png", plot = FCI_plot_unscaled, width=4, height=5)

## FCI plot with error bars and unscaled
FCI_plot_errorbars <- ggplot(results_final_all_spp, aes(x = parameter, y = FCI)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = FCI - SD, ymax = FCI + SD), width = .4) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(0.02,0.08,0.01)) +
  scale_x_continuous(breaks=seq(1985,1996,1)) +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_plot_errorbars

ggsave("../Graphs/Connectivity_plots/FCI_all_spp_CBC_error_unscaled.png", plot = FCI_plot_errorbars, width=7, height=5)

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
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_plot_scaled_error

ggsave("../Graphs/Connectivity_plots/FCI_all_spp_CBC_scaled_error2.png", plot = FCI_plot_scaled_error, width=5, height=5)

#### plot for specialist/generalist ####
results_final_spec$Strategy <- as.factor(results_final_spec$Strategy)
## change 0 to specialist and 1 to generalist
results_final_spec$Strategy <- gsub("0", "Specialist", results_final_spec$Strategy)
results_final_spec$Strategy <- gsub("1", "Generalist", results_final_spec$Strategy)

## add 0.25 to generalist year to offset error bars
library(dplyr)
results_final_spec <- results_final_spec %>% mutate(parameter = ifelse(Strategy=="Generalist", parameter+0.25, parameter))

FCI_plot_spec <- ggplot(results_final_spec, aes(x = parameter, y = rescaled_FCI, group = Strategy)) +
  stat_smooth(aes(linetype=Strategy), colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_ci, ymax = rescaled_FCI + rescaled_ci), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(-100,500,50)) +
  scale_x_continuous(breaks=seq(1985,1996,1)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  labs(aes(linetype=Generalist),size=3) +
  theme(legend.title = element_blank()) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_plot_spec

ggsave("../Graphs/Connectivity_plots/FCI_generalist_specialist_CBC_scaled.png", plot = FCI_plot_spec, width=7, height=5)

FCI_plot_spec_unscaled <- ggplot(results_final_spec, aes(x = parameter, y = FCI, group = Strategy)) +
  stat_smooth(aes(linetype=Strategy), colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = FCI - SD, ymax = FCI + SD), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Functional connectivity index") +
  scale_y_continuous(breaks=seq(0.02,0.1,0.01)) +
  scale_x_continuous(breaks=seq(1985,1996,1)) +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  labs(aes(linetype=Generalist),size=3) +
  theme(legend.position = c(0.8,0.8)) +
  theme(legend.title = element_blank()) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
FCI_plot_spec_unscaled

ggsave("../Graphs/Connectivity_plots/FCI_generalist_specialist_BTO_unscaled.png", plot = FCI_plot_spec_unscaled, width=7, height=5)

#################### maps of CBC sites ######################

#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script") ####

rm(list=ls()) # clear R

## add data
pair_attr <- read.csv("../Data/Bird_sync_data/pair_attr_BTO.csv", header=TRUE)
site_data <- read.csv("../Data/BTO_data/pair_attr_mean_north_dist_hab_sim.csv", header = TRUE)

### save data on included sites for plot
# create the following: site, spp for which it is included for. Do this for 3 start.years (1980, 1985 & 1991)
## TP1 = 1980 DATA 
## TP2 = 1985 DATA 
## TP3 = 1991 DATA 
inc_sites_T1 <- NULL
inc_sites_T2 <- NULL
inc_sites_T3 <- NULL
SPECIES_T1 <- NULL
SPECIES_T2 <- NULL
SPECIES_T3 <- NULL

for (i in unique(pair_attr$spp)){
  temp_attr <- pair_attr[pair_attr$spp==i,]
  inc_sites_T1 <- c(inc_sites_T1, unique(c(temp_attr[temp_attr$start.year=="1980", "site1"], temp_attr[temp_attr$start.year=="1980", "site2"])))
  inc_sites_T2 <- c(inc_sites_T2, unique(c(temp_attr[temp_attr$start.year=="1985", "site1"], temp_attr[temp_attr$start.year=="1985", "site2"])))
  inc_sites_T3 <- c(inc_sites_T3, unique(c(temp_attr[temp_attr$start.year=="1991", "site1"], temp_attr[temp_attr$start.year=="1991", "site2"])))
  SPECIES_T1 <- c(SPECIES_T1, rep(i, length(unique(c(temp_attr[temp_attr$start.year=="1980", "site1"], temp_attr[temp_attr$start.year=="1980", "site2"])))))
  SPECIES_T2 <- c(SPECIES_T2, rep(i, length(unique(c(temp_attr[temp_attr$start.year=="1985", "site1"], temp_attr[temp_attr$start.year=="1985", "site2"])))))
  SPECIES_T3 <- c(SPECIES_T3, rep(i, length(unique(c(temp_attr[temp_attr$start.year=="1991", "site1"], temp_attr[temp_attr$start.year=="1991", "site2"])))))
}

T1_map_data <- data.frame(site = inc_sites_T1, spp = SPECIES_T1)
T2_map_data <- data.frame(site = inc_sites_T2, spp = SPECIES_T2)
T3_map_data <- data.frame(site = inc_sites_T3, spp = SPECIES_T3)

T1_map <- data.frame(site = names(table(T1_map_data$site)), spp_count = as.numeric(table(T1_map_data$site)))
T2_map <- data.frame(site = names(table(T2_map_data$site)), spp_count = as.numeric(table(T2_map_data$site)))
T3_map <- data.frame(site = names(table(T3_map_data$site)), spp_count = as.numeric(table(T3_map_data$site)))

# merge in the site xy data #
temp_site_data <- site_data[,c("site_a", "site_a_NORTH", "site_a_EAST")]
names(temp_site_data) <- c("site", "north", "east")
temp_site_data2 <- site_data[,c("site_b", "site_b_NORTH", "site_b_EAST")]
names(temp_site_data2) <- c("site", "north", "east")

T1_map <- merge(T1_map, unique(rbind(temp_site_data, temp_site_data2)))
T2_map <- merge(T2_map, unique(rbind(temp_site_data, temp_site_data2)))
T3_map <- merge(T3_map, unique(rbind(temp_site_data, temp_site_data2)))

install.packages("blighty")
library(blighty)

par(mfrow=c(1,3)) 

x1 <- (T1_map$east/1000)
y1 <- (T1_map$north/1000)
x2 <- (T2_map$east/1000)
y2 <- (T2_map$north/1000)
x3 <- (T3_map$east/1000)
y3 <- (T3_map$north/1000)

## T1 map 
T1.map <- blighty(place="set.British.Isles") +
  points(x1,y1, col="black", pch=19) +
  title(main = "1980 sites", cex.main=2)

## T2 map
T2.map <- blighty(place="set.British.Isles") +
  points(x2,y2, col="black", pch=19)

## T3 map
T3.map <- blighty(place="set.British.Isles") +
  points(x3,y3, col="black", pch=19)

png("../Graphs/Maps_of_sites_CBC.png", width=1000, height=600)
par(mfrow=c(1,3))
T1.map <- blighty(place="set.British.Isles") +
  points(x1,y1, col="black", pch=19) +
  title(main = "1980 sites", cex.main=2)
T2.map <- blighty(place="set.British.Isles") +
  points(x2,y2, col="black", pch=19) +
  title(main = "1985 sites", cex.main=2)
T3.map <- blighty(place="set.British.Isles") +
  points(x3,y3, col="black", pch=19) +
  title(main = "1990 sites", cex.main=2)
dev.off()
