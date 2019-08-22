###########################################################
## Title: Moving window temporal trend in functional connectivity BBS BTO data  
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: May 2018
##########################################################


#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script/Bird_scripts") ####

rm(list=ls()) # clear R

library(lme4)
library(lmerTest)

################################
## CREATE PAIR ATTRIBUTE DATA ##
################################
## add data ##
final_pair_data <- read.csv("../Data/Bird_sync_data/final_pair_data_all_spp_BBS_1.csv", header = TRUE) # final pair data 1 (can choose any file as all v similar)
site_data <- read.csv("../Data/BTO_data/pair_attr_mean_north_dist_hab_sim_BBS.csv", header = TRUE) # pair attr site data
## family/genus data
species_traits <- read.csv("../Data/BTO_data/woodland_generalist_specialist.csv", header=TRUE)

### remove species from final_pair_data which skipped for SOME years (i.e. do not have synchrony data for whole time series)
## species which were skipped for SOME years:
## 456, 431, 468
final_pair_data <- subset(final_pair_data, spp!="456")
final_pair_data <- subset(final_pair_data, spp!="431")
final_pair_data <- subset(final_pair_data, spp!="468")
length(unique(final_pair_data$spp)) # 24

## merge to site data for each species ##
pair_attr_1 <- merge(final_pair_data, site_data, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))  # merge the site comparisons in one direction - site a to a, b to b
site_data_reverse <- site_data
names(site_data_reverse)[1:6] <- c("site_b", "site_a", "site_b_EAST", "site_b_NORTH", "site_a_EAST", "site_a_NORTH")
pair_attr_2 <- merge(final_pair_data, site_data_reverse, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))	# merge the site comparisons in the same direction using site_data_reverse
pair_attr <- rbind(pair_attr_1, pair_attr_2) # combine the two datasets

# centre and standardize the distance and northing variables
pair_attr[c(13:15)] <- lapply(pair_attr[c(13:15)], function(pair_attr) c(scale(pair_attr, center = TRUE, scale = TRUE))) 
# 
# pair_attr$distance <- pair_attr$distance - mean(na.omit(pair_attr$distance))
# pair_attr$distance <- pair_attr$distance/(max(abs(na.omit(pair_attr$distance))))
# 
# pair_attr$mean_northing <- pair_attr$mean_northing - mean(na.omit(pair_attr$mean_northing))
# pair_attr$mean_northing <- pair_attr$mean_northing/(max(abs(na.omit(pair_attr$mean_northing))))
# 
# pair_attr$renk_hab_sim <- pair_attr$renk_hab_sim - mean(na.omit(pair_attr$renk_hab_sim))
# pair_attr$renk_hab_sim <- pair_attr$renk_hab_sim/(max(abs(na.omit(pair_attr$renk_hab_sim))))

# check colinearity
cor.test(pair_attr$mean_northing, pair_attr$distance) # -0.006
summary(lm(mean_northing ~ renk_hab_sim, data = pair_attr)) # 0.02
summary(lm(mean_northing ~ end.year, data = pair_attr)) # 0.001
summary(lm(distance ~ renk_hab_sim, data = pair_attr)) # 6.128e-06
summary(lm(distance ~ end.year, data = pair_attr)) # 0.000009
summary(lm(renk_hab_sim ~ end.year, data = pair_attr)) # 0.0038

## create pair.id variable and make species and year variables factors ## 
pair_attr$pair.id <- paste("ID", pair_attr$site1, pair_attr$site2, sep = "_")
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$end.year <- as.factor(pair_attr$end.year)
pair_attr$mid.year <- as.factor(pair_attr$mid.year)

## merge in specialism, family and genus info
specialism <- read.csv("../Data/BTO_data/woodland_generalist_specialist.csv", header=TRUE)
pair_attr <- merge(pair_attr, specialism, by.x="spp", by.y="species_code")
length(unique(pair_attr$spp)) ## 24 species

### save final pair_attr file with standardised variables
write.csv(pair_attr, file = "../Data/Bird_sync_data/pair_attr_BBS.csv", row.names = FALSE) # save pair_attr file 

##################################################
## CREATE PAIR ATTRIBUTE DATA WITH CLIMATE DATA ##
##################################################

rm(list=ls()) # clear R
library(lme4)
library(lmerTest)
library(data.table)

## read in data
final_pair_data_rain <- read.csv("../Data/MetOffice_data/final_pair_data_mean_rainfall_BBS.csv", header=TRUE)
final_pair_data_temp <- read.csv("../Data/MetOffice_data/final_pair_data_mean_temp_BBS.csv", header=TRUE)  

## change dataframe to wide format
rain_wide <- reshape(final_pair_data_rain, idvar=c("site1", "site2", "mid.year"), timevar="season", v.names="lag0", direction="wide", sep="_")
rain_wide <- rain_wide[,-c(4:5)] ## remove start and end year columns
names(rain_wide)[4:7] <- c("winter_rain","autumn_rain", "spring_rain", "summer_rain") ## rename sync columns (a, d, b, c)

temp_wide <- reshape(final_pair_data_temp, idvar=c("site1", "site2", "mid.year"), timevar="season", v.names="lag0", direction="wide", sep="_")
temp_wide <- temp_wide[,-c(4:5)] ## remove start and end year columns
names(temp_wide)[4:7] <- c("spring_temp","autumn_temp", "winter_temp", "summer_temp") ## rename sync columns (b,d,a,c)

## put into one dataframe
climate_wide <- merge(rain_wide, temp_wide, by=c("site1", "site2", "mid.year"))
## check correlation between climate synchrony values
climate_cor <- cor(climate_wide[c(4:11)])
climate_cor ## highest correlation is 0.59 - keep all variables


## save climate_wide 
write.csv(climate_wide, file = "../Data/MetOffice_data/climate_synchrony_final_BBS.csv", row.names=FALSE)

###################################################################################
### run model with pop sync and all 8 climate varaibles to pick important ones ###
###################################################################################

climate_wide <- read.csv("../Data/MetOffice_data/climate_synchrony_final_BBS.csv", header=TRUE) ## climate synchrony data
pair_attr <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE) ## read in pop sync data for butterflies
site_data <- read.csv("../Data/MetOffice_data/site_list_5km_BBS.csv", header=TRUE) ## read in site list with 1km and 5km easting and northing

pair_attr$site1 <- as.factor(pair_attr$site1)
pair_attr$site2 <- as.factor(pair_attr$site2)
climate_wide$site1 <- as.factor(climate_wide$site1)
climate_wide$site2 <- as.factor(climate_wide$site2)

## merge climate wide and pair_attr then add in 5k info later
pair_attr1 <- merge(pair_attr, climate_wide, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"), all=TRUE)
climate_wide_rev <- climate_wide
names(climate_wide_rev)[1:2] <- c("site_b", "site_a")
pair_attr2 <- merge(pair_attr, climate_wide_rev, by.x=c("site1", "site2", "mid.year"), by.y=c("site_a", "site_b", "mid.year"), all=TRUE)

pair_attr1$site1 <- as.factor(pair_attr1$site1)
pair_attr1$site2 <- as.factor(pair_attr1$site2)
pair_attr2$site1 <- as.factor(pair_attr2$site1)
pair_attr2$site2 <- as.factor(pair_attr2$site2)

## remove NAs from climate data (sites which don't have climate info - 15 sites)
pair_attr1 <- pair_attr1 %>% 
  filter_at(vars(winter_rain, autumn_rain, spring_rain, summer_rain, winter_temp, autumn_temp, spring_temp, summer_temp), any_vars(!is.na(.)))
pair_attr2 <- pair_attr2 %>% 
  filter_at(vars(winter_rain, autumn_rain, spring_rain, summer_rain, winter_temp, autumn_temp, spring_temp, summer_temp), any_vars(!is.na(.)))

## bind the two dataframes together and remove NAs from lag0
pair_attr_climate <- rbind(pair_attr1, pair_attr2)
pair_attr_climate <- pair_attr_climate[!is.na(pair_attr_climate$lag0),]
## less rows than pair_attr because there are 2490 sites, not 2499

## check that unique sites == 2490
site1 <- unique(pair_attr_climate[,1, drop=FALSE])
site2 <- unique(pair_attr_climate[,2, drop=FALSE])
colnames(site1)[1] <- "site"
colnames(site2)[1] <- "site"
site_list2 <- rbind(site1, site2)
site_list2 <- unique(site_list2) ## 2490 sites

## create unique site_5k number using 5k easting and northing
DT = data.table(site_data)
site_data <- DT[, site_5k := as.numeric(factor(paste(east.5k, north.5k), levels = unique(paste(east.5k, north.5k))))]
site_data <- as.data.frame(site_data)
## remove original site number
site_data2 <- site_data[,-1]

## add in 5km easting and northing data
pair_attr_test1 <- merge(site_data2, pair_attr_climate, by.x=c("east", "north"), by.y=c("site_a_EAST", "site_a_NORTH"))
## rename columns
names(pair_attr_test1)[1:4] <- c("site_a_EAST", "site_a_NORTH", "site_5k_a_EAST", "site_5k_a_NORTH")
names(pair_attr_test1)[5] <- "site_5k_a"
pair_attr_test1 <- merge(site_data2, pair_attr_test1, by.x=c("east", "north"), by.y=c("site_b_EAST", "site_b_NORTH"))
names(pair_attr_test1)[1:5] <- c("site_b_EAST", "site_b_NORTH", "site_5k_b_EAST", "site_5k_b_NORTH", "site_5k_b")
## remove duplicated rows (now the same length as pair_attr_climate)
pair_attr_test2 <- unique(pair_attr_test1)
## reorder columns
pair_attr_climate <- pair_attr_test2[,c(11,12,6,7,1,2,10,5,8,9,3,4,13:34)]
## create new 5km pairID column
pair_attr_climate$pair.id_5k <- paste("ID", pair_attr_climate$site_5k_a, pair_attr_climate$site_5k_b, sep = "_")
length(unique(pair_attr_climate$pair.id_5k)) ## 197,293
length(unique(pair_attr_climate$pair.id)) ## 321,137

## save pop + climate sync file
write.csv(pair_attr_climate, file = "../Data/Bird_sync_data/pop_climate_synchrony_BBS.csv", row.names=FALSE)

###### run model with pop sync as response and climate variables as explanatory
pair_attr_climate$spp <- as.factor(pair_attr_climate$spp)
pair_attr_climate$mid.year <- as.factor(pair_attr_climate$mid.year)
pair_attr_climate$pair.id <- as.character(pair_attr_climate$pair.id)
pair_attr_climate$pair.id_5k <- as.character(pair_attr_climate$pair.id_5k)

pop_climate <- lmer(lag0 ~ mid.year + winter_rain + autumn_rain + spring_rain + summer_rain + winter_temp + autumn_temp + 
                      spring_temp + summer_temp + (1|spp) + (1|pair.id) + (1|pair.id_5k), data=pair_attr_climate)
summary(pop_climate)
anova(pop_climate)
## winter, autumn and spring rain are significant (all positive)
## save results
pop_climate_results <- data.frame(summary(pop_climate)$coefficients[,1:5])
pop_climate_results$parameter <- paste(row.names(pop_climate_results))
rownames(pop_climate_results) <- 1:nrow(pop_climate_results)
## remove mid.year rows
pop_climate_results <- pop_climate_results[-c(1:14),]
## save results
write.csv(pop_climate_results, file = "../Results/Model_outputs/BBS/pop_climate_results_BBS.csv", row.names=FALSE)
## winter, autumn and spring rainfall need to be included in future analyses

####################################################################################################
################################## RUN SYNCHRONY MODELS ############################################
####################################################################################################

pair_attr <- read.csv(file = "../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE) # read in pair_attr file 
pair_attr_climate <- read.csv(file = "../Data/Bird_sync_data/pop_climate_synchrony_BBS.csv", header=TRUE) # read in pair_attr file 

## make sure correct variables are factors ## 
str(pair_attr)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$mid.year <- as.factor(pair_attr$mid.year)
pair_attr$pair.id <- as.character(pair_attr$pair.id)

pair_attr_climate$spp <- as.factor(pair_attr_climate$spp)
pair_attr_climate$mid.year <- as.factor(pair_attr_climate$mid.year)
pair_attr_climate$pair.id <- as.character(pair_attr_climate$pair.id)

################################################
##  model to produce one line for all species ##
################################################

### run model with family to test for a relationship with synchrony
all_spp_model_family <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + 
                               mid.year + family + (1|pair.id) + (1|pair.id_5k) + (1|spp), data = pair_attr_climate)
summary(all_spp_model_family) 
anova(all_spp_model_family) ## non-significant (p=0.78)

### run model with genus to test for a relationship with synchrony
all_spp_model_genus <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + 
                              mid.year + genus + (1|pair.id) + (1|pair.id_5k) + (1|spp), data = pair_attr_climate)
summary(all_spp_model_genus) 
anova(all_spp_model_genus) ## non-significant (p=0.79)

## run model with intercept to get fixed effect results
all_spp_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + 
                        mid.year + (1|pair.id) + (1|pair.id_5k) + (1|spp), data = pair_attr_climate)
### save results for northing, distance and hab sim
fixed_results <- data.frame(summary(all_spp_model)$coefficients[,1:5])
fixed_results$parameter <- paste(row.names(fixed_results))
rownames(fixed_results) <- 1:nrow(fixed_results)
## remove mid.year rows and intercept
fixed_results <- fixed_results[-c(1,5:20),]
## save results
write.csv(fixed_results, file = "../Results/Model_outputs/BBS/fixed_effect_results_climate_bbs.csv", row.names=FALSE)

##  model to produce aggregate synchrony values for all 32 species - no intercept
all_spp_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + 
                        mid.year + (1|pair.id) + (1|pair.id_5k) + (1|spp)-1, data = pair_attr_climate)
## convergence warnings
relgrad <- with(all_spp_model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) ## 0.00011 (ignore warnings if <0.001)

summary(all_spp_model) 
anova(all_spp_model)

## save results ## 
results_table_all_spp <- data.frame(summary(all_spp_model)$coefficients[,1:3])

## change names and add in parameter column ##
names(results_table_all_spp) <- c("FCI", "SD", "t")
results_table_all_spp$parameter <- paste(row.names(results_table_all_spp))
rownames(results_table_all_spp) <- 1:nrow(results_table_all_spp)
## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_table_all_spp <- results_table_all_spp[grep("mid.year", results_table_all_spp$parameter),]
## change parameter names to year
results_table_all_spp$parameter <- rep(1999:2012)
### rescale estimate, SD and CI ### 
results_table_all_spp$rescaled_FCI <- results_table_all_spp$FCI*(100/results_table_all_spp$FCI[1])
results_table_all_spp$rescaled_sd <- results_table_all_spp$SD*(100/results_table_all_spp$FCI[1])
results_table_all_spp$rescaled_ci <- results_table_all_spp$rescaled_sd*1.96

## save final results table ##
write.csv(results_table_all_spp, file = "../Results/Bird_results/results_final_all_spp_climate_BBS.csv", row.names=FALSE)

##################################################
##### model to produce one line per species #####
#################################################

results_table_sp<-NULL
for (i in unique(pair_attr_climate$spp)){
  print(i)
  
  species_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + 
                          mid.year + (1|pair.id_5k) + (1|pair.id)-1, data = pair_attr_climate[pair_attr_climate$spp==i,])
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
results_table3 <- results_table_sp[grep("hab_sim", results_table_sp$parameter),]
results_table4 <- rbind(results_table4, results_table1, results_table2, results_table3)
results_table_sp <- results_table_sp[!results_table_sp$parameter%in%results_table4$parameter,]

## change parameter names to year
results_table_sp$parameter <- rep(1999:2012)

results_final_sp <- NULL
for (i in unique(results_table_sp$sp)){
  
  results_temp_sp <- results_table_sp[results_table_sp$sp==i,]  
  
  results_temp_sp$rescaled_FCI <- results_temp_sp$FCI*(100/results_temp_sp$FCI[1])
  results_temp_sp$rescaled_sd <- results_temp_sp$SD*(100/results_temp_sp$FCI[1])
  results_temp_sp$rescaled_ci <- results_temp_sp$rescaled_sd*1.96
  
  results_final_sp <- rbind(results_final_sp, results_temp_sp)
  
}

### merge with species common name info ###
species_names <- read.csv("../Data/BTO_data/woodland_generalist_specialist.csv", header=TRUE)
results_final_sp <- merge(results_final_sp, species_names, by.x="sp", by.y="species_code")

## save species final results table ##
write.csv(results_final_sp, file = "../Results/Bird_results/results_final_sp_climate_BBS.csv", row.names=FALSE)




#######################################################################################################################
########################################### CLIMATE AND SYNCHRONY ANALYSIS ############################################ 
#######################################################################################################################

## climate data
final_pair_data_temp <- read.csv("../Data/MetOffice_data/final_pair_data_mean_temp_BBS.csv", header=TRUE)  
final_pair_data_rain <- read.csv("../Data/MetOffice_data/final_pair_data_mean_rainfall_BBS.csv", header=TRUE)

## merge in climate data
## spring and winter temperature
final_pair_data_spr_temp <- final_pair_data_temp[final_pair_data_temp$season=="b",] ## spring
final_pair_data_win_temp <- final_pair_data_temp[final_pair_data_temp$season=="a",] ## winter
names(final_pair_data_spr_temp)[3] <- "lag0_spr_temp"
names(final_pair_data_win_temp)[3] <- "lag0_win_temp"
## autumn and winter rainfall
final_pair_data_aut_rain <- final_pair_data_rain[final_pair_data_rain$season=="d",] ## autumn
final_pair_data_win_rain <- final_pair_data_rain[final_pair_data_rain$season=="a",] ## winter
names(final_pair_data_aut_rain)[3] <- "lag0_aut_rain"
names(final_pair_data_win_rain)[3] <- "lag0_win_rain"
## remove some columns (season, start and end year)
final_pair_data_spr_temp <- subset(final_pair_data_spr_temp, select = -c(5:7))
final_pair_data_win_temp <- subset(final_pair_data_win_temp, select = -c(5:7))
final_pair_data_aut_rain <- subset(final_pair_data_aut_rain, select = -c(5:7))
final_pair_data_win_rain <- subset(final_pair_data_win_rain, select = -c(5:7))

### merge in spring temp
pair_attr_temp <- pair_attr
pair_attr_1 <- merge(pair_attr_temp, final_pair_data_spr_temp, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"), all=FALSE)  # merge the site comparisons in one direction - site a to a, b to b
spring_reverse <- final_pair_data_spr_temp
names(spring_reverse)[1:2] <- c("site2", "site1")
pair_attr_2 <- merge(pair_attr_temp, spring_reverse, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"), all=FALSE)	# merge the site comparisons in the same direction using site_data_reverse
pair_attr_temp <- rbind(pair_attr_1, pair_attr_2) # combine the two datasets
pair_attr_temp <- unique(pair_attr_temp)
length(unique(pair_attr_temp$spp))## 24 species
length(unique(pair_attr_temp$site1)) # 2482
length(unique(pair_attr_temp$site2)) # 2471 
### merge in winter temp
pair_attr_1 <- merge(pair_attr_temp, final_pair_data_win_temp, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"))  # merge the site comparisons in one direction - site a to a, b to b
winter_reverse <- final_pair_data_win_temp
names(winter_reverse)[1:2] <- c("site2", "site1")
pair_attr_2 <- merge(pair_attr_temp, winter_reverse, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"))	# merge the site comparisons in the same direction using site_data_reverse
pair_attr_temp <- rbind(pair_attr_1, pair_attr_2) # combine the two datasets
pair_attr_temp <- unique(pair_attr_temp)
length(unique(pair_attr_temp$spp))## 24 species
length(unique(pair_attr_temp$site1)) # 2482
length(unique(pair_attr_temp$site2)) # 2471 

### merge in autumn rain
pair_attr_rain <- pair_attr
pair_attr_1 <- merge(pair_attr_rain, final_pair_data_aut_rain, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"))  # merge the site comparisons in one direction - site a to a, b to b
autumn_reverse <- final_pair_data_aut_rain
names(autumn_reverse)[1:2] <- c("site2", "site1")
pair_attr_2 <- merge(pair_attr_rain, autumn_reverse, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"))	# merge the site comparisons in the same direction using site_data_reverse
pair_attr_rain <- rbind(pair_attr_1, pair_attr_2) # combine the two datasets
pair_attr_rain <- unique(pair_attr_rain)
length(unique(pair_attr_rain$spp))## 24 species
length(unique(pair_attr_rain$site1)) # 2482
length(unique(pair_attr_rain$site2)) # 2471
### merge in winter rain
pair_attr_1 <- merge(pair_attr_rain, final_pair_data_win_rain, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"))  # merge the site comparisons in one direction - site a to a, b to b
winter_reverse <- final_pair_data_win_rain
names(winter_reverse)[1:2] <- c("site2", "site1")
pair_attr_2 <- merge(pair_attr_rain, winter_reverse, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"))	# merge the site comparisons in the same direction using site_data_reverse
pair_attr_rain <- rbind(pair_attr_1, pair_attr_2) # combine the two datasets
pair_attr_rain <- unique(pair_attr_rain)
length(unique(pair_attr_rain$spp))## 24 species
length(unique(pair_attr_rain$site1)) # 2482
length(unique(pair_attr_rain$site2)) # 2471

summary(pair_attr_rain)
summary(pair_attr_temp)

pair_attr_rain$pair.id <- as.character(pair_attr_rain$pair.id)
pair_attr_rain$spp <- as.factor(pair_attr_rain$spp)
pair_attr_rain$end.year <- as.factor(pair_attr_rain$end.year)
pair_attr_rain$mid.year <- as.factor(pair_attr_rain$mid.year)

pair_attr_temp$pair.id <- as.character(pair_attr_temp$pair.id)
pair_attr_temp$spp <- as.factor(pair_attr_temp$spp)
pair_attr_temp$end.year <- as.factor(pair_attr_temp$end.year)
pair_attr_temp$mid.year <- as.factor(pair_attr_temp$mid.year)

######## temperature model
all_spp_model_temp <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + lag0_spr_temp + lag0_win_temp + (1|pair.id) + (1|spp)-1, data = pair_attr_temp)
summary(all_spp_model_temp) 
anova(all_spp_model_temp) ## both temp variables are non-significant

pop_temp_model <- data.frame(summary(all_spp_model_temp)$coefficients[,1:3])
names(pop_temp_model) <- c("FCI", "SD", "t")
pop_temp_model$parameter <- paste(row.names(pop_temp_model))
rownames(pop_temp_model) <- 1:nrow(pop_temp_model)

results_tab5 <- NULL
results_tab1 <- pop_temp_model[grep("mean_northing", pop_temp_model$parameter),]
results_tab2 <- pop_temp_model[grep("distance", pop_temp_model$parameter),]
results_tab3 <- pop_temp_model[grep("renk_hab_sim", pop_temp_model$parameter),]
results_tab4 <- pop_temp_model[grep("lag0_spr_temp", pop_temp_model$parameter),]
results_tab5 <- pop_temp_model[grep("lag0_win_temp", pop_temp_model$parameter),]
results_tab6 <- rbind(results_tab4, results_tab1, results_tab2, results_tab3, results_tab4, results_tab5)
pop_temp_model <- pop_temp_model[!pop_temp_model$parameter%in%results_tab6$parameter,]

## change parameter names to year
pop_temp_model$parameter <- rep(1999:2012)

### rescale estimate, SD and CI ### 
pop_temp_model$rescaled_FCI <- pop_temp_model$FCI*(100/pop_temp_model$FCI[1])
pop_temp_model$rescaled_sd <- pop_temp_model$SD*(100/pop_temp_model$FCI[1])
pop_temp_model$rescaled_ci <- pop_temp_model$rescaled_sd*1.96

## save final results table ##
write.csv(pop_temp_model, file = "../Results/Bird_results/results_final_all_spp_temp_BBS.csv", row.names=FALSE)

## graph
FCI_plot_BBS_temp <- ggplot(pop_temp_model, aes(x = parameter, y = rescaled_FCI)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_sd, ymax = rescaled_FCI + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Population synchrony") +
  #scale_y_continuous(breaks=seq(40,160,10)) +
  scale_x_continuous(breaks=seq(1999,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
FCI_plot_BBS_temp
ggsave("../Graphs/Connectivity_plots/FCI_plot_temp_BBS.png", plot = FCI_plot_BBS_temp, width=7, height=5)

######## rainfall model
all_spp_model_rain <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + lag0_aut_rain + lag0_win_rain + (1|pair.id) + (1|spp)-1, data = pair_attr_rain)
summary(all_spp_model_rain) 
anova(all_spp_model_rain) ## both climate varaibles are non-significant

pop_rain_model <- data.frame(summary(all_spp_model_rain)$coefficients[,1:3])
names(pop_rain_model) <- c("FCI", "SD", "t")
pop_rain_model$parameter <- paste(row.names(pop_rain_model))
rownames(pop_rain_model) <- 1:nrow(pop_rain_model)

results_tab5 <- NULL
results_tab1 <- pop_rain_model[grep("mean_northing", pop_rain_model$parameter),]
results_tab2 <- pop_rain_model[grep("distance", pop_rain_model$parameter),]
results_tab3 <- pop_rain_model[grep("renk_hab_sim", pop_rain_model$parameter),]
results_tab4 <- pop_rain_model[grep("lag0_aut_rain", pop_rain_model$parameter),]
results_tab5 <- pop_rain_model[grep("lag0_win_rain", pop_rain_model$parameter),]
results_tab6 <- rbind(results_tab4, results_tab1, results_tab2, results_tab3, results_tab4, results_tab5)
pop_rain_model <- pop_rain_model[!pop_rain_model$parameter%in%results_tab6$parameter,]

## change parameter names to year
pop_rain_model$parameter <- rep(1999:2012)

### rescale estimate, SD and CI ### 
pop_rain_model$rescaled_FCI <- pop_rain_model$FCI*(100/pop_rain_model$FCI[1])
pop_rain_model$rescaled_sd <- pop_rain_model$SD*(100/pop_rain_model$FCI[1])
pop_rain_model$rescaled_ci <- pop_rain_model$rescaled_sd*1.96

## save final results table ##
write.csv(pop_rain_model, file = "../Results/Bird_results/results_final_all_spp_rain_BBS.csv", row.names=FALSE)

## graph
FCI_plot_BBS_rain <- ggplot(pop_rain_model, aes(x = parameter, y = rescaled_FCI)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_sd, ymax = rescaled_FCI + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Population synchrony") +
  #scale_y_continuous(breaks=seq(40,160,10)) +
  scale_x_continuous(breaks=seq(1999,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
FCI_plot_BBS_rain
ggsave("../Graphs/Connectivity_plots/FCI_plot_rain_BBS.png", plot = FCI_plot_BBS_rain, width=7, height=5)
