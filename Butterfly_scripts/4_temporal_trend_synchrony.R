###########################################################
## Title: Moving window temporal trend in functional connectivity
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: August 2017
##########################################################

#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script") ####

rm(list=ls()) # clear R

library(lme4)
library(lmerTest)

################################
## CREATE PAIR ATTRIBUTE DATA ##
################################

## add data ##
final_pair_data <- read.csv("../Data/Butterfly_sync_data/final_pair_data_all_spp.csv", header = TRUE) # pair attr synchrony data 37 species
site_data <- read.csv("../Data/UKBMS_data/pair_attr_mean_northing_dist_sim.csv", header = TRUE) # pair attr site data
trait_data <- read.csv("../Data/UKBMS_data/species.traits.full.table.csv", header = TRUE) # mobility trait data for each species

## merge to site data for each species ##
pair_attr_1 <- merge(final_pair_data, site_data, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))  # merge the site comparisons in one direction - site a to a, b to b
site_data_reverse <- site_data
names(site_data_reverse)[1:6] <- c("site_b", "site_a", "site_b_EAST", "site_b_NORTH", "site_a_EAST", "site_a_NORTH")
pair_attr_2 <- merge(final_pair_data, site_data_reverse, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))	# merge the site comparisons in the same direction using site_data_reverse
pair_attr <- rbind(pair_attr_1, pair_attr_2) # combine the two datasets
length(unique(pair_attr$spp))## 37 species

pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr <- pair_attr[!pair_attr$spp==121,] # DROPPED AS ITS A DUPLICATE OF SMALL SKIPPER
pair_attr <- pair_attr[!pair_attr$spp==34,] # MIGRANT (CLOUDED YELLOW)
pair_attr <- pair_attr[!pair_attr$spp==122,] # MIGRANT (RED ADMIRAL)
pair_attr <- pair_attr[!pair_attr$spp==123,] # MIGRANT (PAINTED LADY)

pair_attr$spp <- droplevels(pair_attr$spp)
length(unique(pair_attr$spp)) ## 33 species

# centre and standardize the distance and northing variables
### standardise by mean and sd
pair_attr$distance <- (pair_attr$distance - mean(na.omit(pair_attr$distance)))/sd(na.omit(pair_attr$distance))
pair_attr$mean_northing <- (pair_attr$mean_northing - mean(na.omit(pair_attr$mean_northing)))/sd(na.omit(pair_attr$mean_northing))
pair_attr$renk_hab_sim <- (pair_attr$renk_hab_sim - mean(na.omit(pair_attr$renk_hab_sim)))/sd(na.omit(pair_attr$renk_hab_sim))

var(pair_attr$distance) # =1
sd(pair_attr$distance) # =1
var(pair_attr$mean_northing) # =1
sd(pair_attr$mean_northing) # =1
var(pair_attr$renk_hab_sim) # =1
sd(pair_attr$renk_hab_sim) # =1

# check colinearity
cor.test(pair_attr$mean_northing, pair_attr$distance) # -0.045
cor.test(pair_attr$mean_northing, pair_attr$renk_hab_sim) # 0.097
summary(lm(mean_northing ~ end.year, data = pair_attr)) # 0.001
cor.test(pair_attr$distance, pair_attr$renk_hab_sim) # -0.08
summary(lm(distance ~ end.year, data = pair_attr)) # 0.002
summary(lm(renk_hab_sim ~ end.year, data = pair_attr)) # 4.796e-05

## merge with trait data 
pair_attr <- merge(pair_attr, trait_data, by.x="spp", by.y="species")
summary(pair_attr)
length(unique(pair_attr$spp)) # 33 species
pair_attr <- droplevels(pair_attr)
## remove columns not needed
pair_attr <- pair_attr[-c(16, 18:46, 48:50, 52:54)]
summary(pair_attr)
## 33 species - NA's are there because spp 100 has no mobility data

## create pair.id variable and make species and year variables factors ## 
str(pair_attr)
pair_attr$pair.id <- paste("ID", pair_attr$site1, pair_attr$site2, sep = "_")
pair_attr$end.year <- as.factor(pair_attr$end.year)
pair_attr$mid.year <- as.factor(pair_attr$mid.year)
pair_attr$start.year <- as.factor(pair_attr$start.year)
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)

write.csv(pair_attr, file = "../Data/Butterfly_sync_data/pair_attr.csv", row.names = FALSE) # save pair_attr file 

####################################################################################################
################################## RUN SYNCHRONY MODELS ############################################
####################################################################################################

pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE) # load up pair attribute data to save time

## make sure correct variables are factors ## 
str(pair_attr)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$start.year <- as.factor(pair_attr$start.year)
pair_attr$end.year <- as.factor(pair_attr$end.year)
pair_attr$mid.year <- as.factor(pair_attr$mid.year)
pair_attr$pair.id <- as.character(pair_attr$pair.id)

###############################################
### run the synchrony model for all species ###
###############################################
length(unique(pair_attr$spp)) ## 33 species

##  model to produce aggregate FCI model for all 33 species - no intercept and no migrants ##
all_spp_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp)-1, data = pair_attr)
summary(all_spp_model) 
anova(all_spp_model)

# results_table_full_model_ukbms <- data.frame(summary(all_spp_model)$coefficients[,1:5])
# write.csv(results_table_full_model_ukbms, file = "../Results/Model_outputs/full_model_ukbms.csv", row.names=TRUE)

## save and plot results ## 
results_table_all_spp <- data.frame(summary(all_spp_model)$coefficients[,1:3])

## change names and add in parameter column ##
names(results_table_all_spp) <- c("FCI", "SD", "t")
results_table_all_spp$parameter <- paste(row.names(results_table_all_spp))
rownames(results_table_all_spp) <- 1:nrow(results_table_all_spp)

## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_tab4 <- NULL
results_tab1 <- results_table_all_spp[grep("mean_northing", results_table_all_spp$parameter),]
results_tab2 <- results_table_all_spp[grep("distance", results_table_all_spp$parameter),]
results_tab3 <- results_table_all_spp[grep("renk_hab_sim", results_table_all_spp$parameter),]
results_tab4 <- rbind(results_tab4, results_tab1, results_tab2, results_tab3)
results_table_all_spp <- results_table_all_spp[!results_table_all_spp$parameter%in%results_tab4$parameter,]

## change parameter names to year
results_table_all_spp$parameter <- rep(1985:2012)

### rescale estimate, SD and CI ### 
results_table_all_spp$rescaled_FCI <- results_table_all_spp$FCI*(100/results_table_all_spp$FCI[1])
results_table_all_spp$rescaled_sd <- results_table_all_spp$SD*(100/results_table_all_spp$FCI[1])
results_table_all_spp$rescaled_ci <- results_table_all_spp$rescaled_sd*1.96

## save final results table ##
write.csv(results_table_all_spp, file = "../Results/Butterfly_results/results_final_all_spp.csv", row.names=FALSE)


##############################################
## model for all species without covariates ##
##############################################
all_spp_model2 <- lmer(lag0 ~ mid.year + (1|pair.id) + (1|spp), data = pair_attr)
summary(all_spp_model2) 
anova(all_spp_model2)

## AIC of full model against minimal model
AIC(all_spp_model, all_spp_model2)
## all_spp_model = 2037132
## all_spp_model2 = 2037505
## more complicated model has lower AIC (by 373)

################################################
### run the synchrony model for each species ###
################################################

results_table_sp<-NULL
for (i in unique(pair_attr$spp)){
  print(i)
  
  species_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id)-1, data = pair_attr[pair_attr$spp==i,])
  summary(species_model)
  anova(species_model)

### save and plot the results ###
results_table_temp <- data.frame(summary(species_model)$coefficients[,1:3],i)
results_table_sp<-rbind(results_table_sp,results_table_temp)
}

## change names and add in parameter column ##
names(results_table_sp) <- c("FCI", "SD", "t","sp")
results_table_sp$parameter <- paste(row.names(results_table_sp))
rownames(results_table_sp) <- 1:nrow(results_table_sp)

## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_table4 <- NULL
results_table1 <- results_table_sp[grep("mean_northing", results_table_sp$parameter),]
results_table2 <- results_table_sp[grep("distance", results_table_sp$parameter),]
results_table3 <- results_table_sp[grep("renk_hab_sim", results_table_sp$parameter),]
results_table4 <- rbind(results_table4, results_table1, results_table2, results_table3)
results_table_sp <- results_table_sp[!results_table_sp$parameter%in%results_table4$parameter,]

## change parameter names to year
results_table_sp$parameter <- rep(1985:2012)

### rescale estimate, SD and CI for each species 
results_final_sp <- NULL
for (i in unique(results_table_sp$sp)){

results_temp_sp <- results_table_sp[results_table_sp$sp==i,]  
  
results_temp_sp$rescaled_FCI <- results_temp_sp$FCI*(100/results_temp_sp$FCI[1])
results_temp_sp$rescaled_sd <- results_temp_sp$SD*(100/results_temp_sp$FCI[1])
results_temp_sp$rescaled_ci <- results_temp_sp$rescaled_sd*1.96

results_final_sp <- rbind(results_final_sp, results_temp_sp)

}

## Add in common names of species to final results table ##
species_info <- pair_attr[c(1,16:18)]
species_info <- unique(species_info)
results_final_sp <- merge(results_final_sp, species_info, by.x="sp", by.y="spp")

## save final results table ##
write.csv(results_final_sp, file="../Results/Butterfly_results/results_final_sp.csv", row.names=FALSE)


#########################################################################
### run the synchrony model for each strategy (specialist/generalist) ###
#########################################################################

results_table_spec<-NULL
for (g in unique(pair_attr$specialism)){
  
  spec_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp)-1, data = pair_attr[pair_attr$specialism==g,])
  summary(spec_model)
  anova(spec_model)
  
  ### save and plot the results ###
  results_table_spec_temp <- data.frame(summary(spec_model)$coefficients[,1:3],g)
  results_table_spec<-rbind(results_table_spec,results_table_spec_temp)
}

## change names and add in parameter column ##
names(results_table_spec) <- c("FCI", "SD", "t","specialism")
results_table_spec$parameter <- paste(row.names(results_table_spec))
rownames(results_table_spec) <- 1:nrow(results_table_spec)

## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_table_strat4 <- NULL
results_table_strat1 <- results_table_spec[grep("mean_northing", results_table_spec$parameter),]
results_table_strat2 <- results_table_spec[grep("distance", results_table_spec$parameter),]
results_table_strat3 <- results_table_spec[grep("renk_hab_sim", results_table_spec$parameter),]
results_table_strat4 <- rbind(results_table_strat4, results_table_strat1, results_table_strat2, results_table_strat3)
results_table_spec <- results_table_spec[!results_table_spec$parameter%in%results_table_strat4$parameter,]

## change parameter names to year
results_table_spec$parameter <- rep(1985:2012)

### rescale estimate, SD and CI for each species 
results_final_spec <- NULL
for (g in unique(results_table_spec$specialism)){
  
  results_temp_spec <- results_table_spec[results_table_spec$specialism==g,]  
  
  results_temp_spec$rescaled_FCI <- results_temp_spec$FCI*(100/results_temp_spec$FCI[1])
  results_temp_spec$rescaled_sd <- results_temp_spec$SD*(100/results_temp_spec$FCI[1])
  results_temp_spec$rescaled_ci <- results_temp_spec$rescaled_sd*1.96
  
  results_final_spec <- rbind(results_final_spec, results_temp_spec)
  
}

## save final results table	##
write.csv(results_final_spec, file = "../Results/Butterfly_results/results_final_spec.csv", row.names = FALSE) 

# ##############################################################################
# ##################### WOODLAND SPECIES ANALYSIS ##############################
# ##############################################################################
# 
# ### subset pair_attr to woodland species only
# pair_attr_woodland <- subset(pair_attr[pair_attr$HABITAT=="Woodland",])
# 
# ##  model to produce aggregate FCI model for all 9 woodland species - no intercept ##
# all_spp_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp)-1, data = pair_attr_woodland)
# summary(all_spp_model) 
# anova(all_spp_model)
# 
# ## save and plot results ## 
# results_table_all_spp <- data.frame(summary(all_spp_model)$coefficients[,1:3])
# 
# ## change names and add in parameter column ##
# names(results_table_all_spp) <- c("FCI", "SD", "t")
# results_table_all_spp$parameter <- paste(row.names(results_table_all_spp))
# rownames(results_table_all_spp) <- 1:nrow(results_table_all_spp)
# 
# ## take out 3 columns for each species: mean northing, distance and renk_hab_sim
# results_tab4 <- NULL
# results_tab1 <- results_table_all_spp[grep("mean_northing", results_table_all_spp$parameter),]
# results_tab2 <- results_table_all_spp[grep("distance", results_table_all_spp$parameter),]
# results_tab3 <- results_table_all_spp[grep("renk_hab_sim", results_table_all_spp$parameter),]
# results_tab4 <- rbind(results_tab4, results_tab1, results_tab2, results_tab3)
# results_table_all_spp <- results_table_all_spp[!results_table_all_spp$parameter%in%results_tab4$parameter,]
# 
# ## change parameter names to year
# results_table_all_spp$parameter <- rep(1985:2012)
# 
# ### rescale estimate, SD and CI ### 
# results_table_all_spp$rescaled_FCI <- results_table_all_spp$FCI*(100/results_table_all_spp$FCI[1])
# results_table_all_spp$rescaled_sd <- results_table_all_spp$SD*(100/results_table_all_spp$FCI[1])
# results_table_all_spp$rescaled_ci <- results_table_all_spp$rescaled_sd*1.96
# 
# ## save final results table ##
# write.csv(results_table_all_spp, file = "../Results/Butterfly_results/results_final_all_spp_woodland.csv", row.names=FALSE)
# 
# ### STRATEGY MODEL ##
# 
# results_table_strat<-NULL
# for (g in unique(pair_attr_woodland$STRATEGY)){
#   
#   strategy_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp)-1, data = pair_attr_woodland[pair_attr_woodland$STRATEGY==g,])
#   summary(strategy_model)
#   anova(strategy_model)
#   
#   ### save and plot the results ###
#   results_table_strat_temp <- data.frame(summary(strategy_model)$coefficients[,1:3],g)
#   results_table_strat<-rbind(results_table_strat,results_table_strat_temp)
# }
# 
# ## change names and add in parameter column ##
# names(results_table_strat) <- c("FCI", "SD", "t","strategy")
# results_table_strat$parameter <- paste(row.names(results_table_strat))
# rownames(results_table_strat) <- 1:nrow(results_table_strat)
# 
# ## take out 3 columns for each species: mean northing, distance and renk_hab_sim
# results_table_strat4 <- NULL
# results_table_strat1 <- results_table_strat[grep("mean_northing", results_table_strat$parameter),]
# results_table_strat2 <- results_table_strat[grep("distance", results_table_strat$parameter),]
# results_table_strat3 <- results_table_strat[grep("renk_hab_sim", results_table_strat$parameter),]
# results_table_strat4 <- rbind(results_table_strat4, results_table_strat1, results_table_strat2, results_table_strat3)
# results_table_strat <- results_table_strat[!results_table_strat$parameter%in%results_table_strat4$parameter,]
# 
# ## change parameter names to year
# results_table_strat$parameter <- rep(1985:2012)
# 
# ### rescale estimate, SD and CI for each species 
# results_final_strat <- NULL
# for (g in unique(results_table_strat$strategy)){
#   
#   results_temp_strat <- results_table_strat[results_table_strat$strategy==g,]  
#   
#   results_temp_strat$rescaled_FCI <- results_temp_strat$FCI*(100/results_temp_strat$FCI[1])
#   results_temp_strat$rescaled_sd <- results_temp_strat$SD*(100/results_temp_strat$FCI[1])
#   results_temp_strat$rescaled_ci <- results_temp_strat$rescaled_sd*1.96
#   
#   results_final_strat <- rbind(results_final_strat, results_temp_strat)
#   
# }
# 
# results_final_strat$strategy <- as.character(results_final_strat$strategy)
# results_final_strat$strategy <- replace(results_final_strat$strategy, results_final_strat$strategy=="Wider countryside sp", "Wider countryside")
# 
# ## save final results table	##
# write.csv(results_final_strat, file = "../Results/Butterfly_results/results_final_strat_woodland.csv", row.names = FALSE) 
# 
# #######################################################
# ### run the synchrony model for each mobility score ###
# #######################################################
# 
# ## split species into 3 groups based on mobility aility ## 
# pair_attr$mobility.wil <- as.numeric(pair_attr$mobility.wil)
# pair_attr$mobility.score <- cut(pair_attr$mobility.wil, 3, labels=FALSE)
# pair_attr$mobility.score <- as.factor(pair_attr$mobility.score)
# 
# ## run synchrony model for each mobility score ##
# results_table_mobility<-NULL
# for (k in unique(pair_attr$mobility.score)){
#   
#   mobility_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), data = pair_attr[pair_attr$mobility.score==k,])
#   summary(mobility_model)
#   anova(mobility_model)
#   
#   ### save and plot the results ###
#   results_table_mobility_temp <- data.frame(summary(mobility_model)$coefficients[,1:3],k)
#   results_table_mobility<-rbind(results_table_mobility,results_table_mobility_temp)
# }
# 
# names(results_table_mobility) <- c("Estimate", "SD", "t","mobility.score")
# results_table_mobility$parameter <- paste(row.names(results_table_mobility))
# rownames(results_table_mobility) <- 1:nrow(results_table_mobility)
# 
# ## take out 3 columns for each species: mean northing, distance and renk_hab_sim
# results_table_mob4 <- NULL
# results_table_mob1 <- results_table_mobility[grep("mean_northing", results_table_mobility$parameter),]
# results_table_mob2 <- results_table_mobility[grep("distance", results_table_mobility$parameter),]
# results_table_mob3 <- results_table_mobility[grep("renk_hab_sim", results_table_mobility$parameter),]
# results_table_mob4 <- rbind(results_table_mob4, results_table_mob1, results_table_mob2, results_table_mob3)
# results_table_mobility <- results_table_mobility[!results_table_mobility$parameter%in%results_table_mob4$parameter,]
# 
# ## change parameter names to year
# results_table_mobility$parameter <- rep(1985:2012)
# 
# #### calculate FCI for each year relative to intercept (=1985)
# library(dplyr)
# results_table_mobility <- results_table_mobility %>% group_by(mobility.score) %>% mutate(FCI = Estimate + Estimate[parameter==1985]) # this doubles the intercept (1985) value
# results_table_mobility$FCI[results_table_mobility$parameter==1985] <- results_table_mobility$Estimate[results_table_mobility$parameter==1985] # change the first value back to the intercept (1985) estimate value
# 
# ### rescale estimate, SD and CI for each species 
# results_final_mobility <- NULL
# for (k in unique(results_table_mobility$mobility.score)){
#   
#   results_temp_mob <- results_table_mobility[results_table_mobility$mobility.score==k,]  
#   
#   results_temp_mob$rescaled_FCI <- results_temp_mob$FCI*(100/results_temp_mob$FCI[1])
#   results_temp_mob$rescaled_sd <- results_temp_mob$SD*(100/results_temp_mob$FCI[1])
#   results_temp_mob$rescaled_ci <- results_temp_mob$rescaled_sd*1.96
#   
#   results_final_mobility <- rbind(results_final_mobility, results_temp_mob)
#   
# }
# 
# ## change mobility score from numbers to low, medium and high ## 
# results_final_mobility$mobility.score <- as.character(results_final_mobility$mobility.score)
# results_final_mobility$mobility.score[results_final_mobility$mobility.score==1]<-"low"
# results_final_mobility$mobility.score[results_final_mobility$mobility.score==2]<-"medium"
# results_final_mobility$mobility.score[results_final_mobility$mobility.score==3]<-"high"
# 
# write.csv(results_final_mobility, file = "../Results/Butterfly_results/results_final_mobility.csv", row.names = FALSE) # save table	
# 
# ## save pair_attr data with mobility score ##
# write.csv(pair_attr, file = "../Data/Butterfly_sync_data/pair_attr.csv", row.names = FALSE) # save pair_attr file 
# 
