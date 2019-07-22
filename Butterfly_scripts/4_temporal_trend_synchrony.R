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
library(broom)
library(MuMIn)

################################
## CREATE PAIR ATTRIBUTE DATA ##
################################

## add data ##
final_pair_data <- read.csv("../Data/Butterfly_sync_data/final_pair_data_all_spp.csv", header = TRUE) # pair attr synchrony data 35 species
site_data <- read.csv("../Data/UKBMS_data/pair_attr_mean_northing_dist_sim.csv", header = TRUE) # pair attr site data
trait_data <- read.csv("../Data/UKBMS_data/species.traits.full.table.csv", header = TRUE) # mobility trait data for each species

## merge to site data for each species ##
pair_attr_1 <- merge(final_pair_data, site_data, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))  # merge the site comparisons in one direction - site a to a, b to b
site_data_reverse <- site_data
names(site_data_reverse)[1:6] <- c("site_b", "site_a", "site_b_EAST", "site_b_NORTH", "site_a_EAST", "site_a_NORTH")
pair_attr_2 <- merge(final_pair_data, site_data_reverse, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))	# merge the site comparisons in the same direction using site_data_reverse
pair_attr <- rbind(pair_attr_1, pair_attr_2) # combine the two datasets
length(unique(pair_attr$spp))## 35 species

## remove migrant amd duplicated species
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr <- pair_attr[!pair_attr$spp==121,] # DROPPED AS ITS A DUPLICATE OF SMALL SKIPPER
pair_attr <- pair_attr[!pair_attr$spp==122,] # MIGRANT (RED ADMIRAL)
pair_attr <- pair_attr[!pair_attr$spp==123,] # MIGRANT (PAINTED LADY)

pair_attr$spp <- droplevels(pair_attr$spp)
length(unique(pair_attr$spp)) ## 32 species

# centre and standardize the distance and northing variables
### standardise by mean and sd
pair_attr$distance <- (pair_attr$distance - mean(na.omit(pair_attr$distance)))/sd(na.omit(pair_attr$distance))
pair_attr$mean_northing <- (pair_attr$mean_northing - mean(na.omit(pair_attr$mean_northing)))/sd(na.omit(pair_attr$mean_northing))
pair_attr$renk_hab_sim <- (pair_attr$renk_hab_sim - mean(na.omit(pair_attr$renk_hab_sim)))/sd(na.omit(pair_attr$renk_hab_sim))

## double check SD is 1
sd(pair_attr$distance) # =1
sd(pair_attr$mean_northing) # =1
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
length(unique(pair_attr$spp)) # 32 species
pair_attr <- droplevels(pair_attr)
## remove columns not needed
pair_attr <- pair_attr[-c(9:10, 18, 22:48, 50:52, 54:56)]
summary(pair_attr)
length(unique(pair_attr$spp)) # 32 species
## 32 species - NA's are there because spp 100 has no mobility data

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

## read in data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE) # pair attribute synchrony data for 35 species

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

length(unique(pair_attr$spp)) ## 32 species

## first check whether family and genus is significant in main model (phylogenetic checks)
### run model with family to test for a relationship with synchrony
all_spp_model_family <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + family + (1|pair.id) + (1|spp), data = pair_attr)
summary(all_spp_model_family) 
anova(all_spp_model_family) ## family non-significant

### run model with genus to test for a relationship with synchrony
all_spp_model_genus <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + genus + (1|pair.id) + (1|spp), data = pair_attr)
summary(all_spp_model_genus) 
anova(all_spp_model_genus) ## genus non-significant

## run model with intercept to get fixed effect results
## model with intercept
all_spp_model_int <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), data = pair_attr)

### save results for northing, distance and hab sim
fixed_results <- data.frame(summary(all_spp_model_int)$coefficients[,1:5])
fixed_results$parameter <- paste(row.names(fixed_results))
rownames(fixed_results) <- 1:nrow(fixed_results)
r.squaredGLMM(all_spp_model_int)
## R^2 marginal (fixed) = 0.020, R^2 conditional (random) = 0.222
## remove mid.year rows
fixed_results <- fixed_results[-c(1,5:31),]
## save results
write.csv(fixed_results, file = "../Results/Model_outputs/UKBMS/fixed_effect_results_ukbms.csv", row.names=FALSE)

##  model to produce aggregate synchrony values for all 32 species - no intercept
all_spp_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp)-1, data = pair_attr)
summary(all_spp_model) 
anova(all_spp_model)

# save results ## 
results_table_all_spp <- data.frame(summary(all_spp_model)$coefficients[,1:3])
## change names and add in parameter column ##
names(results_table_all_spp) <- c("FCI", "SD", "t")
results_table_all_spp$parameter <- paste(row.names(results_table_all_spp))
rownames(results_table_all_spp) <- 1:nrow(results_table_all_spp)
## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_table_all_spp <- results_table_all_spp[grep("mid.year", results_table_all_spp$parameter),]
## change parameter names to year
results_table_all_spp$parameter <- rep(1985:2012)
### rescale estimate, SD and CI ### 
results_table_all_spp$rescaled_FCI <- results_table_all_spp$FCI*(100/results_table_all_spp$FCI[1])
results_table_all_spp$rescaled_sd <- results_table_all_spp$SD*(100/results_table_all_spp$FCI[1])
results_table_all_spp$rescaled_ci <- results_table_all_spp$rescaled_sd*1.96

## save final results table ##
write.csv(results_table_all_spp, file = "../Results/Butterfly_results/results_final_all_spp.csv", row.names=FALSE)

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
results_table_sp <- results_table_sp[grep("mid.year", results_table_sp$parameter),]
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
species_info <- pair_attr[c(1,16)]
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


#######################################################################################################################
########################################### CLIMATE AND SYNCHRONY ANALYSIS ############################################ 
#######################################################################################################################

## read in climate data
final_pair_data_rain <- read.csv("../Data/MetOffice_data/final_pair_data_mean_rainfall.csv", header=TRUE)
final_pair_data_temp <- read.csv("../Data/MetOffice_data/final_pair_data_mean_temp.csv", header=TRUE)

## summer, autumn and winter temperature and rainfall (these all decline then increase over time)
final_pair_data_sum_temp <- final_pair_data_temp[final_pair_data_temp$season=="c",] ## summer
final_pair_data_aut_temp <- final_pair_data_temp[final_pair_data_temp$season=="d",] ## autumn
final_pair_data_win_temp <- final_pair_data_temp[final_pair_data_temp$season=="a",] ## winter

names(final_pair_data_sum_temp)[3] <- "lag0_sum_temp"
names(final_pair_data_aut_temp)[3] <- "lag0_aut_temp"
names(final_pair_data_win_temp)[3] <- "lag0_win_temp"

## remove some columns (season, start and end year)
final_pair_data_sum_temp <- subset(final_pair_data_sum_temp, select = -c(5:7))
final_pair_data_aut_temp <- subset(final_pair_data_aut_temp, select = -c(5:7))
final_pair_data_win_temp <- subset(final_pair_data_win_temp, select = -c(5:7))

### merge in summer temp
pair_attr_temp <- pair_attr
pair_attr_1 <- merge(pair_attr_temp, final_pair_data_sum_temp, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"), all=FALSE)  # merge the site comparisons in one direction - site a to a, b to b
summer_reverse <- final_pair_data_sum_temp
names(summer_reverse)[1:2] <- c("site2", "site1")
pair_attr_2 <- merge(pair_attr_temp, summer_reverse, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"), all=FALSE)	# merge the site comparisons in the same direction using site_data_reverse
pair_attr_temp <- rbind(pair_attr_1, pair_attr_2) # combine the two datasets
pair_attr_temp <- unique(pair_attr_temp)
length(unique(pair_attr_temp$spp))## 32 species
length(unique(pair_attr_temp$site1)) # 454
length(unique(pair_attr_temp$site2)) # 473 
### merge in autumn temp
pair_attr_1 <- merge(pair_attr_temp, final_pair_data_aut_temp, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"))  # merge the site comparisons in one direction - site a to a, b to b
autumn_reverse <- final_pair_data_aut_temp
names(autumn_reverse)[1:2] <- c("site2", "site1")
pair_attr_2 <- merge(pair_attr_temp, autumn_reverse, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"))	# merge the site comparisons in the same direction using site_data_reverse
pair_attr_temp <- rbind(pair_attr_1, pair_attr_2) # combine the two datasets
length(unique(pair_attr_temp$spp))## 32 species
length(unique(pair_attr_temp$site1)) # 2482
length(unique(pair_attr_temp$site2)) # 2471 
### merge in winter temp
pair_attr_1 <- merge(pair_attr_temp, final_pair_data_win_temp, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"))  # merge the site comparisons in one direction - site a to a, b to b
winter_reverse <- final_pair_data_win_temp
names(winter_reverse)[1:2] <- c("site2", "site1")
pair_attr_2 <- merge(pair_attr_temp, winter_reverse, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"))	# merge the site comparisons in the same direction using site_data_reverse
pair_attr_temp <- rbind(pair_attr_1, pair_attr_2) # combine the two datasets
length(unique(pair_attr_temp$spp))## 32 species
length(unique(pair_attr_temp$site1)) # 2482
length(unique(pair_attr_temp$site2)) # 2471 


### run model with summer, autumn and winter temperature
all_spp_model_temp <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + lag0_sum_temp + lag0_aut_temp + lag0_win_temp + (1|pair.id) + (1|spp), data = pair_attr_temp)
summary(all_spp_model_temp) 
anova(all_spp_model_temp) ## summer and autumn temperature are significant, winter is non-significant
## summer is positive and autumn is negative

## same model without intercept
all_spp_model_temp2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + lag0_sum_temp + lag0_aut_temp + lag0_win_temp + (1|pair.id) + (1|spp)-1, data = pair_attr_temp)

pop_temp_model <- data.frame(summary(all_spp_model_temp2)$coefficients[,1:3])
names(pop_temp_model) <- c("FCI", "SD", "t")
pop_temp_model$parameter <- paste(row.names(pop_temp_model))
rownames(pop_temp_model) <- 1:nrow(pop_temp_model)

pop_temp_model <- pop_temp_model[grep("mid.year", pop_temp_model$parameter),]

## change parameter names to year
pop_temp_model$parameter <- rep(1985:2012)

### rescale estimate, SD and CI ### 
pop_temp_model$rescaled_FCI <- pop_temp_model$FCI*(100/pop_temp_model$FCI[1])
pop_temp_model$rescaled_sd <- pop_temp_model$SD*(100/pop_temp_model$FCI[1])
pop_temp_model$rescaled_ci <- pop_temp_model$rescaled_sd*1.96

## save final results table ##
# write.csv(pop_temp_model, file = "../Results/Butterfly_results/results_final_all_spp_temperature.csv", row.names=FALSE)




