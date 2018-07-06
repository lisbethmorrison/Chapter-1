###########################################################
## Title: Moving window functional connectivity graphs BBS
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: March 2018
##########################################################

rm(list=ls()) # clear R

library(ggplot2)

## add final data 
final_data <- read.csv("../Data/Bird_sync_data/final_pair_data_all_spp_BBS.csv", header=TRUE)

head(final_data)
summary(final_data)
length(unique(final_data$spp)) ## 10 species

## basic plot of raw synchrony against mid year of moving window 
## to get an idea of the trend

## take the average synchrony value across all sites for each moving window
bbs_species <- unique(final_data$spp)
bbs_synchrony_results <- NULL

for (g in bbs_species){ # loop for each species #
  
  species.tab<-final_data[final_data$spp==g,] 
  head(species.tab)
  print(paste("species",g))
  
  year.list <- unique(final_data$mid.year)   
  
  for (i in year.list){
    
    species.10.yr.data<-species.tab[species.tab$mid.year==i,]
    
    species<-g
    year<-i
    mean_synchrony <- mean(species.10.yr.data$lag0)
    results.temp<-data.frame(year,mean_synchrony,species)
    bbs_synchrony_results<-rbind(bbs_synchrony_results,results.temp)
    
  } # end i in year.list
} # end g in ukbms_species

bbs_synchrony_results$species <- as.factor(bbs_synchrony_results$species)
bbs_synchrony_results$year <- as.factor(bbs_synchrony_results$year)

## merge in species name data
common_name <- read.csv("../Data/BTO_data/woodland_generalist_specialist.csv", header=TRUE)
bbs_synchrony_results <- merge(bbs_synchrony_results, common_name, by.x="species", by.y="species_code")
## change species to common_name
names(bbs_synchrony_results)[4] <- "common_name"

synchrony_bbs <- ggplot(bbs_synchrony_results, aes(x=year, y=mean_synchrony, group=common_name)) +
  stat_smooth(aes(colour=common_name), method=loess, se=FALSE) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Raw population synchrony") +
  scale_y_continuous(breaks=seq(0,0.2,0.05)) +
  scale_x_continuous(breaks=seq(1999,2012,3)) +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
synchrony_bbs
