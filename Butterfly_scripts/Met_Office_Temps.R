###########################################################
## Title: Trend of Met Office temperatures over time  
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: October 2017
##########################################################


#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script") ####

require(ggplot2)

rm(list=ls())

## load met office temperature data ##
met_office_temp <- read.csv("../Data/Temp_data/Met_Office_Seasonal_Temps.csv", header=T)
head(met_office_temp)

##### PLOT SEASONAL TEMPERATURES OVER TIME USING MIN AND MAX AS CONFIDENCE INTERVALS ######

#### plot spring temp over time ####
spring_temp <- ggplot(met_office_temp, aes(x=year, y=mean_MAM)) +
  geom_point(size=1) +
  geom_line(lwd=0.5) +
  geom_ribbon(aes(ymin=min_MAM, ymax=max_MAM), linetype=2, alpha=0.1, lwd=0.5) +
  labs(x = "Year", y = "Spring temperature") +
  scale_y_continuous(breaks=seq(3,16,2)) +
  scale_x_continuous(breaks=seq(1980,2016,2)) +
  theme_bw()
spring_temp

#### plot summer temp over time ####
summer_temp <- ggplot(met_office_temp, aes(x=year, y=mean_JJA)) +
  geom_point(size=1) +
  geom_line(lwd=0.5) +
  geom_ribbon(aes(ymin=min_JJA, ymax=max_JJA), linetype=2, alpha=0.1, lwd=0.5) +
  labs(x = "Year", y = "Summer temperature") +
  scale_y_continuous(breaks=seq(10,25,2)) +
  scale_x_continuous(breaks=seq(1980,2016,2)) +
  theme_bw()
summer_temp

#### plot autumn temp over time ####
autumn_temp <- ggplot(met_office_temp, aes(x=year, y=mean_SON)) +
  geom_point(size=1) +
  geom_line(lwd=0.5) +
  geom_ribbon(aes(ymin=min_SON, ymax=max_SON), linetype=2, alpha=0.1, lwd=0.5) +
  labs(x = "Year", y = "Autumn temperature") +
  scale_y_continuous(breaks=seq(5,17,2)) +
  scale_x_continuous(breaks=seq(1980,2016,2)) +
  theme_bw()
autumn_temp

#### plot winter temp over time ####
winter_temp <- ggplot(met_office_temp, aes(x=year, y=mean_DJF)) +
  geom_point(size=1) +
  geom_line(lwd=0.5) +
  geom_ribbon(aes(ymin=min_DJF, ymax=max_DJF), linetype=2, alpha=0.1, lwd=0.5) +
  labs(x = "Year", y = "Winter temperature") +
  scale_y_continuous(breaks=seq(0,10,2)) +
  scale_x_continuous(breaks=seq(1980,2016,2)) +
  theme_bw()
winter_temp

require(gridExtra)
temp_graphs <- grid.arrange(spring_temp, summer_temp, autumn_temp, winter_temp, ncol=2)

ggsave("../Graphs/Met_office_temp_plots/met_office_temp_graphs.png", plot=temp_graphs, width=15, height=8)

spring_model <- lm(mean_MAM ~ year, data=met_office_temp)
summary(spring_model)
anova(spring_model)

summer_model <- lm(mean_JJA ~ year, data=met_office_temp)
summary(summer_model)
anova(summer_model)

autumn_model <- lm(mean_SON ~ year, data=met_office_temp)
summary(autumn_model)
anova(autumn_model)

winter_model <- lm(mean_DJF ~ year, data=met_office_temp)
summary(winter_model)
anova(winter_model)
