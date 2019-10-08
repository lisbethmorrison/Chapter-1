###########################################################
## Title: Moving window temporal trend in functional connectivity CBC BTO data  
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: November 2017
##########################################################


#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script/Bird_scripts") ####

rm(list=ls()) # clear R

library(lme4)
library(lmerTest)
options(scipen=999)

################################
## CREATE PAIR ATTRIBUTE DATA ##
################################

## add data ##
final_pair_data <- read.csv("../Data/Bird_sync_data/final_pair_data_all_spp_CBC.csv", header = TRUE) # pair attr synchrony data
site_data <- read.csv("../Data/BTO_data/pair_attr_mean_north_dist_hab_sim_CBC.csv", header = TRUE) # pair attr site data

## merge to site data for each species ##
pair_attr_1 <- merge(final_pair_data, site_data, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))  # merge the site comparisons in one direction - site a to a, b to b
site_data_reverse <- site_data
names(site_data_reverse)[1:6] <- c("site_b", "site_a", "site_b_EAST", "site_b_NORTH", "site_a_EAST", "site_a_NORTH")
pair_attr_2 <- merge(final_pair_data, site_data_reverse, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))	# merge the site comparisons in the same direction using site_data_reverse
pair_attr <- rbind(pair_attr_1, pair_attr_2) # combine the two datasets

# centre and standardize the distance and northing variables
pair_attr[c(17:18)] <- lapply(pair_attr[c(17:18)], function(pair_attr) c(scale(pair_attr, center = TRUE, scale = TRUE))) 

# pair_attr$distance <- pair_attr$distance - mean(na.omit(pair_attr$distance))
# pair_attr$distance <- pair_attr$distance/(max(abs(na.omit(pair_attr$distance))))
# 
# pair_attr$mean_northing <- pair_attr$mean_northing - mean(na.omit(pair_attr$mean_northing))
# pair_attr$mean_northing <- pair_attr$mean_northing/(max(abs(na.omit(pair_attr$mean_northing))))

# check colinearity
cor.test(pair_attr$mean_northing, pair_attr$distance) # -0.002
summary(lm(mean_northing ~ hab_sim, data = pair_attr)) # 0.02
summary(lm(mean_northing ~ end.year, data = pair_attr)) # 0.001
summary(lm(distance ~ hab_sim, data = pair_attr)) # 6.128e-06
summary(lm(distance ~ end.year, data = pair_attr)) # 0.000009
summary(lm(hab_sim ~ end.year, data = pair_attr)) # 0.0038

## create pair.id variable and make species and year variables factors ## 
pair_attr$pair.id <- paste("ID", pair_attr$site1, pair_attr$site2, sep = "_")
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$end.year <- as.factor(pair_attr$end.year)
pair_attr$mid.year <- as.factor(pair_attr$mid.year)
pair_attr$hab_sim <- as.factor(pair_attr$hab_sim)

pair_attr <- subset(pair_attr, spp!= "370") # all sites have hab_sim = 1, spp = nightingale 

length(unique(pair_attr$spp)) # 26 species

## merge in trait data (specialism, family and genus)
specialism <- read.csv("../Data/BTO_data/woodland_generalist_specialist.csv", header=TRUE)
pair_attr <- merge(pair_attr, specialism, by.x="spp", by.y="species_code")
length(unique(pair_attr$spp)) ## 26 species

write.csv(pair_attr, file = "../Data/Bird_sync_data/pair_attr_CBC.csv", row.names = FALSE) # save pair_attr file (31 spp)

##################################################
## CREATE PAIR ATTRIBUTE DATA WITH CLIMATE DATA ##
##################################################

rm(list=ls()) # clear R
library(lme4)
library(lmerTest)
library(data.table)

## read in data
final_pair_data_rain <- read.csv("../Data/MetOffice_data/final_pair_data_mean_rainfall_CBC.csv", header=TRUE)
final_pair_data_temp <- read.csv("../Data/MetOffice_data/final_pair_data_mean_temp_CBC.csv", header=TRUE)  

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
climate_cor ## highest correlation is 0.6 - keep all variables


## save climate_wide 
write.csv(climate_wide, file = "../Data/MetOffice_data/climate_synchrony_final_CBC.csv", row.names=FALSE)

###################################################################################
### run model with pop sync and all 8 climate varaibles to pick important ones ###
###################################################################################

climate_wide <- read.csv("../Data/MetOffice_data/climate_synchrony_final_CBC.csv", header=TRUE) ## climate synchrony data
pair_attr <- read.csv("../Data/Bird_sync_data/pair_attr_CBC.csv", header=TRUE) ## read in pop sync data for butterflies
site_data <- read.csv("../Data/MetOffice_data/site_list_5km_CBC.csv", header=TRUE) ## read in site list with 1km and 5km easting and northing
## 106 CBC sites
## 96 5km squares with climate data

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

## check that unique sites == 106
site1 <- unique(pair_attr_climate[,1, drop=FALSE])
site2 <- unique(pair_attr_climate[,2, drop=FALSE])
colnames(site1)[1] <- "site"
colnames(site2)[1] <- "site"
site_list2 <- rbind(site1, site2)
site_list2 <- unique(site_list2) ## 106 sites

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
pair_attr_climate <- pair_attr_test2[,c(11,12,6,7,1,2,10,5,8,9,3,4,13:37)]
## create new 5km pairID column
pair_attr_climate$pair.id_5k <- paste("ID", pair_attr_climate$site_5k_a, pair_attr_climate$site_5k_b, sep = "_")
length(unique(pair_attr_climate$pair.id_5k)) ## 1599
length(unique(pair_attr_climate$pair.id)) ## 1923

## save pop + climate sync file
write.csv(pair_attr_climate, file = "../Data/Bird_sync_data/pop_climate_synchrony_CBC.csv", row.names=FALSE)

###### run model with pop sync as response and climate variables as explanatory
pair_attr_climate$spp <- as.factor(pair_attr_climate$spp)
pair_attr_climate$mid.year <- as.factor(pair_attr_climate$mid.year)
pair_attr_climate$pair.id <- as.character(pair_attr_climate$pair.id)
pair_attr_climate$pair.id_5k <- as.character(pair_attr_climate$pair.id_5k)

pop_climate <- lmer(lag0 ~ mid.year + winter_rain + autumn_rain + spring_rain + summer_rain + winter_temp + autumn_temp + 
                      spring_temp + summer_temp + (1|spp) + (1|pair.id) + (1|pair.id_5k), data=pair_attr_climate)
summary(pop_climate)
anova(pop_climate)
## summer temp is significant (positive) (p=0.023) 
## save results
pop_climate_results <- data.frame(summary(pop_climate)$coefficients[,1:5])
pop_climate_results$parameter <- paste(row.names(pop_climate_results))
rownames(pop_climate_results) <- 1:nrow(pop_climate_results)
## remove mid.year rows
pop_climate_results <- pop_climate_results[-c(1:12),]
## save results
write.csv(pop_climate_results, file = "../Results/Model_outputs/CBC/pop_climate_results_CBC.csv", row.names=FALSE)
## only summer temperature needs to be added to future analyses

####################################################################################################
################################## RUN SYNCHRONY MODELS ############################################
####################################################################################################

# pair_attr <- read.csv("../Data/Bird_sync_data/pair_attr_CBC.csv", header=TRUE) ## read in data
## and read in pop and climate synchrony
pair_attr <- read.csv("../Data/Bird_sync_data/pop_climate_synchrony_CBC.csv", header=TRUE) 

## check all variables are correct format
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$end.year <- as.factor(pair_attr$end.year)
pair_attr$mid.year <- as.factor(pair_attr$mid.year)
str(pair_attr)

###############################################
### run the synchrony model for all species ###
###############################################

## first check whether family and genus is significant in main model (phylogenetic checks)

### run model with family to test for a relationship with synchrony
all_spp_model_family <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + family + summer_temp + (1|pair.id) + (1|spp), data = pair_attr)
summary(all_spp_model_family) 
anova(all_spp_model_family) ## family significant (p=0.042)

all_spp_model_genus <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + genus + summer_temp + (1|pair.id) + (1|spp), data = pair_attr)
summary(all_spp_model_genus) 
anova(all_spp_model_genus) ## genus non-significant (p=0.21)

## run model with intercept to get fixed effect results
## model with intercept
all_spp_model_int <- lmer(lag0 ~ mean_northing + distance + hab_sim + summer_temp + mid.year + (1|pair.id) + (1|spp), data = pair_attr)
summary(all_spp_model_int)

### save results for northing, distance and hab sim
fixed_results <- data.frame(summary(all_spp_model_int)$coefficients[,1:5])
fixed_results$parameter <- paste(row.names(fixed_results))
rownames(fixed_results) <- 1:nrow(fixed_results)
## remove mid.year rows
fixed_results <- fixed_results[-c(1,5:16),]
## save results
write.csv(fixed_results, file = "../Results/Model_outputs/CBC/fixed_effect_results_cbc.csv", row.names=FALSE)

##  model to produce aggregate synchrony values for all 26 species - no intercept (with climate data - summer temp is significant)
all_spp_model <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + summer_temp + (1|pair.id) + (1|spp)-1, data = pair_attr)
summary(all_spp_model) 
anova(all_spp_model)

### check model fit ###
hist(residuals(all_spp_model))
qqnorm(resid(all_spp_model))
qqline(resid(all_spp_model)) # not too bad

## save results ## 
results_table_all_spp <- data.frame(summary(all_spp_model)$coefficients[,1:3])

## change names and add in parameter column ##
names(results_table_all_spp) <- c("FCI", "SD", "t")
results_table_all_spp$parameter <- paste(row.names(results_table_all_spp))
rownames(results_table_all_spp) <- 1:nrow(results_table_all_spp)

## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_table_all_spp <- results_table_all_spp[grep("mid.year", results_table_all_spp$parameter),]

## change parameter names to year
results_table_all_spp$parameter <- rep(1985:1996)

### rescale estimate, SD and CI ### 
results_table_all_spp$rescaled_FCI <- results_table_all_spp$FCI*(100/results_table_all_spp$FCI[1])
results_table_all_spp$rescaled_sd <- results_table_all_spp$SD*(100/results_table_all_spp$FCI[1])
results_table_all_spp$rescaled_ci <- results_table_all_spp$rescaled_sd*1.96

## save final results table ##
write.csv(results_table_all_spp, file = "../Results/Bird_results/results_final_all_spp_climate_CBC.csv", row.names=FALSE)

## run model without with all fixed effects and then minus climate varaible to compare R2
climate_model1 <- lmer(lag0 ~ mean_northing + distance + hab_sim + summer_temp + mid.year + (1|pair.id) + (1|spp), data = pair_attr)
climate_model2 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + (1|pair.id) + (1|spp), data = pair_attr)
## calc R2 and take difference (marginal == fixed effects)
r.squaredGLMM(climate_model1)
r.squaredGLMM(climate_model2)
x = 0.0008511256 - 0.000607122
x ## 0.00024
x2 = x*100 ## 0.024%


##################################################
##### model to produce one line per species #####
#################################################

results_table_sp<-NULL
for (i in unique(pair_attr_climate$spp)){
  
  species_model <- lmer(lag0 ~ mean_northing + distance + hab_sim + summer_temp + mid.year + (1|pair.id)-1, data = pair_attr_climate[pair_attr_climate$spp==i,])
  print(i)
  summary(species_model)
  anova(species_model)
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(species_model)$coefficients[,1:3],i)
  results_table_sp<-rbind(results_table_sp,results_table_temp)
}

## lots of errors when use pair.id_5k random effect - so this is taken out (permutation tests run on this anyway so shouldn't be a problem)
## change names and add in parameter column ##
names(results_table_sp) <- c("FCI", "SD", "t","sp")
results_table_sp$parameter <- paste(row.names(results_table_sp))
rownames(results_table_sp) <- 1:nrow(results_table_sp)
## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_table_sp <- results_table_sp[grep("mid.year", results_table_sp$parameter),]
## change parameter names to year
results_table_sp$parameter <- rep(1985:1996)

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
length(unique(results_final_sp$sp)) # 26 species
## save species final results table ##
write.csv(results_final_sp, file = "../Results/Bird_results/results_final_sp_climate_CBC.csv", row.names=FALSE)

# #########################################################
# ### run the synchrony model for generalist/specialist ###
# #########################################################
# 
# ## merge with generalist/specialist data ##
# str(pair_attr)
# pair_attr$specialism <- as.factor(pair_attr$specialism)
# 
# results_table_spec<-NULL
# for (j in unique(pair_attr$specialism)){
#   
#   specialist_model <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + (1|pair.id) + (1|spp)-1, data = pair_attr[pair_attr$specialism==j,])
#   summary(specialist_model)
#   anova(specialist_model)
#   
#   ### save and plot the results ###
#   results_table_spec_temp <- data.frame(summary(specialist_model)$coefficients[,1:3],j)
#   results_table_spec<-rbind(results_table_spec,results_table_spec_temp)
# }
# 
# ## change names and add in parameter column ##
# names(results_table_spec) <- c("FCI", "SD", "t","Strategy")
# results_table_spec$parameter <- paste(row.names(results_table_spec))
# rownames(results_table_spec) <- 1:nrow(results_table_spec)
# 
# ## take out 3 columns for each species: mean northing, distance and renk_hab_sim
# results_table_spec <- results_table_spec[grep("mid.year", results_table_spec$parameter),]
# 
# ## change parameter names to year
# results_table_spec$parameter <- rep(1985:1996)
# 
# ## rescale FCI, SD and CI
# results_final_spec <- NULL
# for (k in unique(results_table_spec$Strategy)){
#   
#   results_temp_spec <- results_table_spec[results_table_spec$Strategy==k,]  
#   
#   results_temp_spec$rescaled_FCI <- results_temp_spec$FCI*(100/results_temp_spec$FCI[1])
#   results_temp_spec$rescaled_sd <- results_temp_spec$SD*(100/results_temp_spec$FCI[1])
#   results_temp_spec$rescaled_ci <- results_temp_spec$rescaled_sd*1.96
#   
#   results_final_spec <- rbind(results_final_spec, results_temp_spec)
#   
# }
# 
# ## save species specialism final results
# write.csv(results_final_spec, file = "../Results/Bird_results/results_final_spec_CBC.csv", row.names=FALSE)


