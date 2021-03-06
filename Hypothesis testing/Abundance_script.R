#############################################################################
## Title: Abundance testing for woodland birds and butterflies
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: April 2018
#############################################################################

rm(list=ls()) # clear R
library(plyr) # load packages
library(dplyr) 
library(lme4)
library(ggplot2)
library(lmerTest)
library(plotrix)
options(scipen=999)


################################# AVERAGE SYNCHRONY #####################################
 
###################################
#### Abundance model for UKBMS ####
###################################

## read data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pop_climate_synchrony.csv", header=TRUE) # butterfly pair attribute data
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pop_climate_synchrony_CBC.csv", header=TRUE) # CBC pair attribute data
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pop_climate_synchrony_BBS.csv", header=TRUE) # BBS pair attribute data
bird_common <- read.csv("../Data/BTO_data/pop_estimates_birds.csv", header=TRUE)
WCBS_data <- read.csv("../Data/UKBMS_data/WCBS_data.csv", header=TRUE)

######### UKBMS ###########
## remove unneccessary columns from WCBS data
WCBS_data <- WCBS_data[-c(3,5:7,9:11,13:15,17:19,21:23,25:27,29:31,33:35,37:38)] ## left with abundance from 09-17
## take the average abundance over time
WCBS_data$average_abundance <- rowMeans(subset(WCBS_data, select=c(2:10), na.rm=TRUE)) ## this is our measure of 'commonness'
## save dataframe
write.csv(WCBS_data, file="../Data/Butterfly_sync_data/Average_abundance_UKBMS.csv", row.names=FALSE)

## merge the two datasets
pair_attr <- merge(pair_attr, WCBS_data, by.x="spp", by.y="Species_code", all=FALSE)
length(unique(pair_attr$spp)) # 31 species (Grizzled Skipper doesn't have abundance data)
summary(pair_attr)

## rescale pop estimate
pair_attr$pop_est_stand <- (pair_attr$average_abundance - mean(na.omit(pair_attr$average_abundance)))/sd(na.omit(pair_attr$average_abundance))

pair_attr$mid.year <- as.factor(pair_attr$mid.year)
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)

### full model with average_abundance as a measure of commonness
common_model_ukbms2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain +
                              winter_temp + autumn_temp + spring_temp + summer_temp + pop_est_stand + (1|pair.id) + (1|spp), data = pair_attr)

qqnorm(residuals(common_model_ukbms2))
qqline(residuals(common_model_ukbms2))
plot(common_model_ukbms2)

summary(common_model_ukbms2)
anova(common_model_ukbms2)
## non-significant
results_table_abund_ukbms <- data.frame(summary(common_model_ukbms2)$coefficients[,1:5]) ## 31 species
write.csv(results_table_abund_ukbms, file = "../Results/Model_outputs/UKBMS/average_abund_ukbms.csv", row.names=TRUE)

#################################
#### Abundance model for CBC ####
#################################

## merge the two datasets
pair_attr_CBC <- merge(pair_attr_CBC, bird_common, by.x="spp", by.y="species_code", all=FALSE)

## rescale pop estimate
pair_attr_CBC$pop_estimate_stand <- (pair_attr_CBC$pop_estimate - mean(na.omit(pair_attr_CBC$pop_estimate)))/sd(na.omit(pair_attr_CBC$pop_estimate))

### main model with pop_estimate as a measure of commonness
common_model_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + summer_temp + pop_estimate_stand + (1|family) + (1|spp) + (1|pair.id), data = pair_attr_CBC)
## this model fails to converge
## run model without family random effect and check if coefficients are still similar
common_model_cbc2 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + summer_temp + pop_estimate_stand + (1|spp) + (1|pair.id), data = pair_attr_CBC)

summary(common_model_cbc) ## pop estimate is significant (p=0.0128)
summary(common_model_cbc2) ## pop estimate is significant (p=0.0105)
## both models are very similar
## save results from one with family

results_table_abund_cbc <- data.frame(summary(common_model_cbc)$coefficients[,1:5]) ## 26 species
write.csv(results_table_abund_cbc, file = "../Results/Model_outputs/CBC/average_abund_cbc.csv", row.names=TRUE)

#################################
#### Abundance model for BBS ####
#################################

## merge the two datasets
pair_attr_BBS <- merge(pair_attr_BBS, bird_common, by.x="spp", by.y="species_code", all=FALSE)

## rescale pop estimate
pair_attr_BBS$pop_estimate_stand <- (pair_attr_BBS$pop_estimate - mean(na.omit(pair_attr_BBS$pop_estimate)))/sd(na.omit(pair_attr_BBS$pop_estimate))

### full model with pop_estimate as a measure of commonness
common_model_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + 
                           pop_estimate_stand + (1|pair.id) + (1|spp), data = pair_attr_BBS)
anova(common_model_bbs)
summary(common_model_bbs)
## non-significant
## save result
results_table_abund_bbs <- data.frame(summary(common_model_bbs)$coefficients[,1:5]) ## 24 species
write.csv(results_table_abund_bbs, file = "../Results/Model_outputs/BBS/average_abund_bbs.csv", row.names=TRUE)


#######################################################################################################
####################################################################################################### 
#######################################################################################################

################################# CHANGE IN SYNCHRONY #####################################

###################################
#### Abundance model for UKBMS ####
###################################

### add abundance data
abundance_data <- read.csv("../Data/UKBMS_data/Collated_Indices_2016.csv", header=TRUE)
## add synchrony data
results_final_sp <- read.csv("../Results/Butterfly_results/results_final_sp.csv", header=TRUE)

## removed columns not needed (leaving species and habitat info)
results_final_sp2 <- results_final_sp[-c(2:8)]
results_final_sp2 <- unique(results_final_sp2)

## merge files together
abundance_data <- merge(abundance_data, results_final_sp2, by.x="Species.code", by.y="sp")
length(unique(abundance_data$Common.name)) ## 32 species
## 3 species missing from the abundance data - clouded yellow, red admiral and painted lady (all migrant species)
## remove some columns not needed
abundance_data <- abundance_data[-c(2,5,7,9:10)]
abundance_data <- droplevels(abundance_data)
length(unique(abundance_data$Species.code))
####### do species significantly change in abundance between early and late years?

## subset two abundance data frames
## abundance t1 with years 1980-1989 (mid-year=1985)
## abundance t2 with years 1995-2004 (mid-year=2000)
abundance_t1 <- subset(abundance_data[(abundance_data$Year>=1980) & (abundance_data$Year<=1989),])
abundance_t1$Year <- "TP1"
abundance_t2 <- subset(abundance_data[(abundance_data$Year>=1995) & (abundance_data$Year<=2004),])
abundance_t2$Year <- "TP2"
## rbind together
abundance_final <- rbind(abundance_t1, abundance_t2)
str(abundance_final)
################## T TEST ON EACH SPECIES BETWEEN EARLY AND LATE YEARS ############################
## calculates t test between early and late, produces p value, and produces mean of early and mean of late years
abundance_results = grouped_df(abundance_final, vars="Species.code") %>% summarise(p=t.test(Collated.Index[which(Year=="TP1")], 
                   Collated.Index[which(Year=="TP2")])$p.value, mean_x=mean(Collated.Index[which(Year=="TP1")]), 
                   mean_y=mean(Collated.Index[which(Year=="TP2")]))

## add in column to say whether species t test is significant or not
abundance_results$significance <- ifelse(abundance_results$p<0.05, "yes", "no")

## mean difference column of mean_y[LATE] - mean_x[EARLY]
## therefore negative difference == species has decrease in abundance between early and late years
abundance_results$mean_change_abund <- abundance_results$mean_y - abundance_results$mean_x
abundance_results$ab_change_85_00 <- ifelse(abundance_results$mean_change_abund>0, "increase", "decrease")                                                               
## save abundance results
write.csv(abundance_results, file="../Results/Butterfly_results/abundance_results_85_00.csv", row.names=FALSE)
abundance_results <- read.csv("../Results/Butterfly_results/abundance_results_85_00.csv", header=TRUE)

################ Interaction between mid-year and abundance change ####################
## UKBMS butterflies
pair_attr <- read.csv("../Data/Butterfly_sync_data/pop_climate_synchrony.csv", header=TRUE) ## MUST READ IN DATA AGAIN
## subset to only look at mid years 1985 and 2000
pair_attr_1985 <- pair_attr[pair_attr$mid.year==1984.5,]
pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_ukbms <- rbind(pair_attr_1985, pair_attr_2000) ## 2 years

pair_attr_ukbms$mid.year <- as.factor(pair_attr_ukbms$mid.year)
pair_attr_ukbms$pair.id <- as.character(pair_attr_ukbms$pair.id)
pair_attr_ukbms$spp <- as.factor(pair_attr_ukbms$spp)

## merge with abundance data 
pair_attr_ukbms <- merge(pair_attr_ukbms, abundance_results, by.x="spp", by.y="Species.code", all=FALSE)

## full model with interaction between mid year and mobility
abund_model1 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + summer_rain + 
                       winter_temp + autumn_temp + spring_temp + summer_temp + mid.year*ab_change_85_00 + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
summary(abund_model1)
anova(abund_model1)
## not significant (p=0.22)
results_table_abund_ukbms <- data.frame(summary(abund_model1)$coefficients[,1:5]) ## 32 species
write.csv(results_table_abund_ukbms, file = "../Results/Model_outputs/UKBMS/change_abund_85_00_ukbms.csv", row.names=TRUE)

################ same again but with abundance from 00-12 ################ 

## subset two abundance data frames
## abundance t3 with years 1995-2004 (mid-year=2000)
## abundance t4 with years 2007-2016 (mid-year=2012)
abundance_t3 <- subset(abundance_data[(abundance_data$Year>=1995) & (abundance_data$Year<=2004),])
abundance_t3$Year <- "TP3"
abundance_t4 <- subset(abundance_data[(abundance_data$Year>=2007) & (abundance_data$Year<=2016),])
abundance_t4$Year <- "TP4"
## rbind together
abundance_final <- rbind(abundance_t3, abundance_t4)
str(abundance_final)
################## T TEST ON EACH SPECIES BETWEEN EARLY AND LATE YEARS ############################
## calculates t test between early and late, produces p value, and produces mean of early and mean of late years
abundance_results = grouped_df(abundance_final, vars="Species.code") %>% summarise(p=t.test(Collated.Index[which(Year=="TP3")], 
                    Collated.Index[which(Year=="TP4")])$p.value, mean_x=mean(Collated.Index[which(Year=="TP3")]), 
                    mean_y=mean(Collated.Index[which(Year=="TP4")]))

## add in column to say whether species t test is significant or not
abundance_results$significance <- ifelse(abundance_results$p<0.05, "yes", "no")

## mean difference column of mean_y[LATE] - mean_x[EARLY]
## therefore negative difference == species has decrease in abundance between early and late years
abundance_results$mean_change_abund <- abundance_results$mean_y - abundance_results$mean_x
abundance_results$ab_change_00_12 <- ifelse(abundance_results$mean_change_abund>0, "increase", "decrease")                                                               
## save abundance results
write.csv(abundance_results, file="../Results/Butterfly_results/abundance_results_00_12.csv", row.names=FALSE)
abundance_results <- read.csv("../Results/Butterfly_results/abundance_results_00_12.csv", header=TRUE)

################ Interaction between mid-year and abundance change ####################
## UKBMS butterflies
## subset to only look at mid years 1985 and 2000
pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]
pair_attr_ukbms <- rbind(pair_attr_2000, pair_attr_2012) ## 2 years

pair_attr_ukbms$mid.year <- as.numeric(pair_attr_ukbms$mid.year)
pair_attr_ukbms$pair.id <- as.character(pair_attr_ukbms$pair.id)
pair_attr_ukbms$spp <- as.factor(pair_attr_ukbms$spp)

## merge with abundance data 
pair_attr_ukbms <- merge(pair_attr_ukbms, abundance_results, by.x="spp", by.y="Species.code", all=FALSE)

## full model with interaction between mid year and mobility
abund_model2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + summer_rain + 
                       winter_temp + autumn_temp + spring_temp + summer_temp + mid.year*ab_change_00_12 + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
summary(abund_model2)
anova(abund_model2)
## very significant interaction (p<0.00001)
results_table_abund_ukbms <- data.frame(summary(abund_model2)$coefficients[,1:5]) ## 32 species
write.csv(results_table_abund_ukbms, file = "../Results/Model_outputs/UKBMS/change_abund_00_12_ukbms.csv", row.names=TRUE)
############################################# PLOT REUSLTS ############################################# 

############### LATE MODEL ############### 
## run model without abundance for each species and extract the slope (mid.year) coefficient for each
results_table<-NULL
for (i in unique(pair_attr_ukbms$spp)){
  print(i)
  ## if statement: if if the mixed effects model doesn't throw an error (levels of pairID), then run an lmer, else run a linear model 
  loadError=F
  a=try({spec_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain + 
                              winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id), data=pair_attr_ukbms[pair_attr_ukbms$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain + 
                         winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id), data=pair_attr_ukbms[pair_attr_ukbms$spp==i,])
  }else{
    spec_ukbms <- lm(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain + 
                       winter_temp + autumn_temp + spring_temp + summer_temp, data=pair_attr_ukbms[pair_attr_ukbms$spp==i,])
  }
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(spec_ukbms)$coefficients[,1:2],i)
  results_table_temp$parameter <- paste(row.names(results_table_temp))
  rownames(results_table_temp) <- 1:nrow(results_table_temp)
  results_table_temp <- results_table_temp[grep("mid.year", results_table_temp$parameter),]
  results_table<-rbind(results_table,results_table_temp)
}

## merge in this dataframe with abundance info 
abundance <- pair_attr_ukbms[,c(1,43)]
abundance <- unique(abundance)
results_table <- merge(results_table, abundance, by.x="i", by.y="spp")
## remove parameter column
results_table <- results_table[,-c(4)]
## re-name some columns
names(results_table) <- c("species", "slope", "SE", "Abundance")
results_table$Abundance <- revalue(results_table$Abundance, c("increase"="Increase"))
results_table$Abundance <- revalue(results_table$Abundance, c("decrease"="Decrease"))

########## calculate slope and SE from main model (with specialism)
### plot points and SE ontop of raw data
results_table_abund_ukbms <- results_table_abund_ukbms[,-c(3:5)]
## leave rows with year and interaction
results_table_abund_ukbms$Abundance <- paste(row.names(results_table_abund_ukbms))
rownames(results_table_abund_ukbms) <- 1:nrow(results_table_abund_ukbms)
results_table_abund_ukbms <- results_table_abund_ukbms[c(13,15),]
results_table_abund_ukbms$Abundance <- revalue(results_table_abund_ukbms$Abundance, c("mid.year"="Decrease"))
results_table_abund_ukbms$Abundance <- revalue(results_table_abund_ukbms$Abundance, c("mid.year:ab_change_00_12increase"="Increase"))

## create new dataframe
decrease_slope <- results_table_abund_ukbms[1,1]
increase_slope <- sum(results_table_abund_ukbms$Estimate)
decrease_SE <- results_table_abund_ukbms[1,2]
interaction_SE <- results_table_abund_ukbms[2,2]
increase_SE <- sqrt((decrease_SE)^2) + ((interaction_SE)^2)
model_summary <- data.frame(Abundance=c("Decrease", "Increase"), slope=c(decrease_slope, increase_slope), SE=c(decrease_SE, increase_SE))

model_summary$Abundance <- factor(model_summary$Abundance, levels=c("Increase", "Decrease"))
results_table$Abundance <- factor(results_table$Abundance, levels=c("Increase", "Decrease"))

png("../Graphs/Abundance/Abundance_change_ukbms_00_12.png", height = 150, width = 180, units = "mm", res = 300)
abund <- ggplot(mapping=aes(x=Abundance, y=slope, ymin=slope-SE, ymax=slope+SE)) +
  geom_point(data=results_table, size=2, color="lightgrey", position=myjit) +
  geom_errorbar(data=results_table, color="lightgrey", position=myjit, width=0.1) +
  geom_point(data=model_summary, size=4) +
  geom_errorbar(data=model_summary, width=0.1) +
  labs(x="Change in abundance", y="Change in population synchrony 2000-2012") +
  scale_y_continuous(breaks = seq(-0.02,0.03,0.01)) +
  theme_bw() +
  theme(text = element_text(size = 10), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black"))
abund
dev.off()

readPNG("../Graphs/Specialism/Specialism_change_predicted_ukbms_85_00.png", native = FALSE, info = FALSE)

library(ggpubr)
png("../Graphs/FINAL/Figure4.png", height = 200, width = 150, units = "mm", res = 300)
ggarrange(spec, abund, 
          labels = c("(a)", "(b)"), font.label = list(size = 10, color ="black"),
          ncol = 1, nrow = 2)
dev.off()

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


###################################################################################
###################################################################################
################################ WOODLAND BIRDS ###################################
###################################################################################
###################################################################################

## Woodland bird change in abundance
## add abundance files
file_names <- list.files("../Data/BTO_data/CBC_BBS_abundance", full.names=TRUE) # where all the files are stored
bird_abundance <- do.call(rbind,lapply(file_names,read.csv)) # create new dataframe with all species rbind together

## remove columns not needed (unsmoothed and CI columns)
bird_abundance <- bird_abundance[-c(2,4,5)]
## remove years before 1980
bird_abundance <- bird_abundance[bird_abundance$year>=1980,]

#################################
#### Abundance model for CBC ####
#################################

## for cbc only interested in years 1980-1989 and 1991-2000
bird_abundance_cbc_1 <- subset(bird_abundance[(bird_abundance$year>=1980) & (bird_abundance$year<=1989),])
bird_abundance_cbc_1$time <- "Early"
bird_abundance_cbc_2 <- subset(bird_abundance[(bird_abundance$year>=1991) & (bird_abundance$year<=2000),])
bird_abundance_cbc_2$time <- "Late"
bird_abundance_cbc <- rbind(bird_abundance_cbc_1, bird_abundance_cbc_2)

## now do t test on difference in abundance between early and late time periods
## calculates t test between early and late, produces p value, and produces mean of early and mean of late years
abundance_results_cbc = bird_abundance_cbc %>% group_by(species_code) %>% summarise(p=t.test(sm[which(time=="Early")], 
                        sm[which(time=="Late")])$p.value, mean_x=mean(sm[which(time=="Early")]), mean_y=mean(sm[which(time=="Late")]))

## add in column to say whether species t test is significant or not
abundance_results_cbc$significance <- ifelse(abundance_results_cbc$p<0.05, "yes", "no")

## mean difference column of mean_y[LATE] - mean_x[EARLY]
## therefore negative difference == species has decrease in abundance between early and late years
abundance_results_cbc$mean_abund_change <- abundance_results_cbc$mean_y - abundance_results_cbc$mean_x
abundance_results_cbc$ab_change_85_96 <- ifelse(abundance_results_cbc$mean_abund_change>0, "increase", "decrease")                                                               

## save abundance results
write.csv(abundance_results_cbc, file="../Results/Bird_results/abundance_results_85_96.csv", row.names=FALSE)
abundance_results_cbc <- read.csv("../Results/Bird_results/abundance_results_85_96.csv", header=TRUE)

####### Interaction with mid-year and abundance change
pair_attr_cbc_1985 <- pair_attr_CBC[pair_attr_CBC$mid.year==1984.5,]
pair_attr_cbc_1996 <- pair_attr_CBC[pair_attr_CBC$mid.year==1995.5,]
pair_attr_cbc <- rbind(pair_attr_cbc_1985, pair_attr_cbc_1996)

pair_attr_cbc$mid.year <- as.factor(pair_attr_cbc$mid.year)
pair_attr_cbc$pair.id <- as.character(pair_attr_cbc$pair.id)
pair_attr_cbc$spp <- as.factor(pair_attr_cbc$spp)

## merge abundance data
pair_attr_cbc <- merge(pair_attr_cbc, abundance_results_cbc, by.x="spp", by.y="species_code", all=FALSE)
length(unique(pair_attr_cbc$spp)) # 22 species 
#pair_attr_cbc <- pair_attr_cbc[-c(21:25)]

###### run model with strategy and year interaction
abund_model_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + summer_temp + mid.year*ab_change_85_96 + (1|family) + (1|pair.id) + (1|spp), data = pair_attr_cbc)
summary(abund_model_cbc) ## interaction is non-significant (p=0.0601)
## singular fit warning (species random effect has zero variance)
## check if similar to model without species random effect
# abund_model_cbc2 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*ab_change_85_96 + (1|family) + (1|pair.id), data = pair_attr_cbc)
# summary(abund_model_cbc2)
## results are very similar even when species random effect removed
## save results from model with family

anova(abund_model_cbc)
results_table_abund_cbc <- data.frame(summary(abund_model_cbc)$coefficients[,1:5]) ## 22 species
write.csv(results_table_abund_cbc, file = "../Results/Model_outputs/CBC/change_abund_cbc.csv", row.names=TRUE)

################################################################################
################################# BBS BIRDS ####################################
################################################################################

## for bbs only interested in years 1994-2003 and 2007-2016
## but only have data to 2011 for now
bird_abundance_bbs_1 <- subset(bird_abundance[(bird_abundance$year>=1994) & (bird_abundance$year<=2003),])
bird_abundance_bbs_1$time <- "Early"
bird_abundance_bbs_2 <- subset(bird_abundance[(bird_abundance$year>=2007) & (bird_abundance$year<=2016),]) 
bird_abundance_bbs_2$time <- "Late"
bird_abundance_bbs <- rbind(bird_abundance_bbs_1, bird_abundance_bbs_2)
## take out year 2001 => no data (foot&mouth)
bird_abundance_bbs<-bird_abundance_bbs[!(bird_abundance_bbs$year=="2001"),]

## now do t test on difference in abundance between early and late time periods
## calculates t test between early and late, produces p value, and produces mean of early and mean of late years
bird_abundance_bbs$species_code <- as.factor(bird_abundance_bbs$species_code)
abundance_results_bbs <- grouped_df(bird_abundance_bbs, vars="species_code") %>% summarise(p=t.test(sm[which(time=="Early")], sm[which(time=="Late")])$p.value, mean_x=mean(sm[which(time=="Early")]), mean_y=mean(sm[which(time=="Late")]))

## add in column to say whether species t test is significant or not
abundance_results_bbs$significance <- ifelse(abundance_results_bbs$p<0.05, "yes", "no")

## mean difference column of mean_y[LATE] - mean_x[EARLY]
## therefore negative difference == species has decrease in abundance between early and late years
abundance_results_bbs$mean_abund_change <- abundance_results_bbs$mean_y - abundance_results_bbs$mean_x
abundance_results_bbs$ab_change_99_12 <- ifelse(abundance_results_bbs$mean_abund_change>0, "increase", "decrease")                                                               
## save abundance results
write.csv(abundance_results_bbs, file="../Results/Bird_results/abundance_results_99_12.csv", row.names=FALSE)
abundance_results_bbs <- read.csv("../Results/Bird_results/abundance_results_99_12.csv", header=TRUE)

#### BBS birds
pair_attr_bbs_1999 <- pair_attr_BBS[pair_attr_BBS$mid.year==1998.5,]
pair_attr_bbs_2012 <- pair_attr_BBS[pair_attr_BBS$mid.year==2011.5,]
pair_attr_bbs <- rbind(pair_attr_bbs_1999, pair_attr_bbs_2012)

pair_attr_bbs$mid.year <- as.factor(pair_attr_bbs$mid.year)
pair_attr_bbs$pair.id <- as.character(pair_attr_bbs$pair.id)
pair_attr_bbs$spp <- as.factor(pair_attr_bbs$spp)

## merge with abundance data
pair_attr_bbs <- merge(pair_attr_bbs, abundance_results_bbs, by.x="spp", by.y="species_code", all=FALSE)

###### run model with abund and year interaction
abund_model_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + 
                          mid.year*ab_change_99_12 + (1|pair.id) + (1|spp), data = pair_attr_bbs)
summary(abund_model_bbs)
## intercept is mid.year 1999
anova(abund_model_bbs)
## mid.year*mean_abund_change interaction is non-significant (0.26)
results_table_abund_bbs <- data.frame(summary(abund_model_bbs)$coefficients[,1:5]) ## 18 species
write.csv(results_table_abund_bbs, file = "../Results/Model_outputs/BBS/change_abund_bbs.csv", row.names=TRUE)
