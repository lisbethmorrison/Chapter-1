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

## remove frequency columns
final_pair_data <- final_pair_data[,-c(9,10)]

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
pair_attr[c(13:15)] <- lapply(pair_attr[c(13:15)], function(pair_attr) c(scale(pair_attr, center = TRUE, scale = TRUE))) 

## double check SD is 1 and mean = 0
sd(pair_attr$distance) # =1
sd(pair_attr$mean_northing) # =1
sd(pair_attr$renk_hab_sim) # =1
mean(pair_attr$distance) # =0
mean(pair_attr$mean_northing) # =0
mean(pair_attr$renk_hab_sim) # =0

## merge with trait data 
pair_attr <- merge(pair_attr, trait_data, by.x="spp", by.y="species")
summary(pair_attr)
length(unique(pair_attr$spp)) # 32 species
pair_attr <- droplevels(pair_attr)
## remove columns not needed
pair_attr <- pair_attr[-c(20:46,48:50,52:54)]
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

##################################################
## CREATE PAIR ATTRIBUTE DATA WITH CLIMATE DATA ##
##################################################

rm(list=ls()) # clear R
library(lme4)
library(lmerTest)
library(data.table)

## read in data
final_pair_data_rain <- read.csv("../Data/MetOffice_data/final_pair_data_mean_rainfall.csv", header=TRUE)
final_pair_data_temp <- read.csv("../Data/MetOffice_data/final_pair_data_mean_temp.csv", header=TRUE)  

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
climate_cor ## highest correlation is 0.56 - keep all variables


## save climate_wide 
write.csv(climate_wide, file = "../Data/MetOffice_data/climate_synchrony_final.csv", row.names=FALSE)

###################################################################################
### run model with pop sync and all 8 climate varaibles to pick important ones ###
###################################################################################

climate_wide <- read.csv("../Data/MetOffice_data/climate_synchrony_final.csv", header=TRUE) ## climate synchrony data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE) ## read in pop sync data for butterflies
site_data <- read.csv("../Data/MetOffice_data/site_list_5km.csv", header=TRUE) ## read in site list with 1km and 5km easting and northing

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

## check that unique sites == 686
site1 <- unique(pair_attr_climate[,1, drop=FALSE])
site2 <- unique(pair_attr_climate[,2, drop=FALSE])
colnames(site1)[1] <- "site"
colnames(site2)[1] <- "site"
site_list2 <- rbind(site1, site2)
site_list2 <- unique(site_list2) ## 686 sites

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
pair_attr_climate <- pair_attr_test2[,c(11,12,6,7,1,2,10,5,8,9,3,4,13:36)]
## create new 5km pairID column
pair_attr_climate$pair.id_5k <- paste("ID", pair_attr_climate$site_5k_a, pair_attr_climate$site_5k_b, sep = "_")
length(unique(pair_attr_climate$pair.id_5k)) ## 23,265
length(unique(pair_attr_climate$pair.id)) ## 37,616

## save pop + climate sync file
write.csv(pair_attr_climate, file = "../Data/Butterfly_sync_data/pop_climate_synchrony.csv", row.names=FALSE)
pair_attr_climate <- read.csv("../Data/Butterfly_sync_data/pop_climate_synchrony.csv", header=TRUE)

###### run model with pop sync as response and climate variables as explanatory
pair_attr_climate$spp <- as.factor(pair_attr_climate$spp)
pair_attr_climate$mid.year <- as.factor(pair_attr_climate$mid.year)
pair_attr_climate$pair.id <- as.character(pair_attr_climate$pair.id)
pair_attr_climate$pair.id_5k <- as.character(pair_attr_climate$pair.id_5k)

pop_climate <- lmer(lag0 ~ mid.year + winter_rain + autumn_rain + spring_rain + summer_rain + winter_temp + autumn_temp + 
                      spring_temp + summer_temp + (1|spp) + (1|pair.id) + (1|pair.id_5k), data=pair_attr_climate)
summary(pop_climate)

r.squaredGLMM(pop_climate) ## R^2 marginal = 0.017 (1.7%)

## run model again with only

anova(pop_climate)
## all 8 are significant
## most positive
## autumn temperature negative
## save results
pop_climate_results <- data.frame(summary(pop_climate)$coefficients[,1:5])
pop_climate_results$parameter <- paste(row.names(pop_climate_results))
rownames(pop_climate_results) <- 1:nrow(pop_climate_results)
## remove mid.year rows
pop_climate_results <- pop_climate_results[-c(1:28),]
## save results
write.csv(pop_climate_results, file = "../Results/Model_outputs/UKBMS/pop_climate_results_ukbms.csv", row.names=FALSE)
## all 8 variables need to be added to future analyses


####################################################################################################
################################## RUN SYNCHRONY MODELS ############################################
####################################################################################################

## use pop & climate sync data

## read in data
# pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE) # pair attribute synchrony data for 35 species
## and read in pop and climate synchrony
pair_attr <- read.csv("../Data/Butterfly_sync_data/pop_climate_synchrony.csv", header=TRUE)  # pop synchrony for 32 species and climate synchrony for 8 variables

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
all_spp_model_family <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + 
                               summer_rain + winter_temp + autumn_temp + spring_temp + summer_temp + family + (1|pair.id) + (1|spp), data = pair_attr)
summary(all_spp_model_family) 
anova(all_spp_model_family) ## family non-significant

### run model with genus to test for a relationship with synchrony
all_spp_model_genus <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + 
                              summer_rain + winter_temp + autumn_temp + spring_temp + summer_temp + genus + (1|pair.id) + (1|spp), data = pair_attr)
summary(all_spp_model_genus) 
anova(all_spp_model_genus) ## genus non-significant

## run model with intercept to get fixed effect results
## model with intercept
all_spp_model_int <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + 
                            summer_rain + winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr)

### save results for northing, distance and hab sim
fixed_results <- data.frame(summary(all_spp_model_int)$coefficients[,1:5])
fixed_results$parameter <- paste(row.names(fixed_results))
rownames(fixed_results) <- 1:nrow(fixed_results)
## remove mid.year rows
fixed_results <- fixed_results[-c(1,5:39),]
## save results
write.csv(fixed_results, file = "../Results/Model_outputs/UKBMS/fixed_effect_results_ukbms.csv", row.names=FALSE)

## run model without with all fixed effects and then minus climate varaible to compare R2
climate_model1 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + 
                         summer_rain + winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr)
climate_model2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), data = pair_attr)
## calc R2 and take difference (marginal == fixed effects)
r.squaredGLMM(climate_model1)
r.squaredGLMM(climate_model2)
x = 0.02077165 - 0.02042361
x2 = x*100 ## 0.035% or 0.00035 proportion

###########################################################################################
##  model to produce aggregate synchrony values for all 32 species - no intercept
all_spp_model2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + 
                        summer_rain + winter_temp + autumn_temp + spring_temp + summer_temp + mid.year + (1|pair.id)
                        + (1|pair.id_5k) + (1|spp)-1, data = pair_attr_climate)
## convergence warnings when 5k pairID random effect is in the model
relgrad <- with(all_spp_model2@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) ## 0.000093 (ignore warnings if <0.001)

summary(all_spp_model2) 
anova(all_spp_model2)

hist(residuals(all_spp_model2))
## normally distributed
qqnorm(resid(all_spp_model))
qqline(resid(all_spp_model)) ## looks good

# save results ## 
results_table_all_spp <- data.frame(summary(all_spp_model2)$coefficients[,1:3])
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
write.csv(results_table_all_spp, file = "../Results/Butterfly_results/results_final_all_spp_climate.csv", row.names=FALSE)

################################################
### run the synchrony model for each species ###
################################################

results_table_sp<-NULL
for (i in unique(pair_attr_climate$spp)){
  print(i)
  
  species_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + 
                          summer_rain + winter_temp + autumn_temp + spring_temp + summer_temp + mid.year + 
                          (1|pair.id_5k) + (1|pair.id)-1, data = pair_attr_climate[pair_attr_climate$spp==i,])
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
species_info <- pair_attr_climate[c(14,22)]
species_info <- unique(species_info) ## 32 species
results_final_sp <- merge(results_final_sp, species_info, by.x="sp", by.y="spp")

## save final results table ##
write.csv(results_final_sp, file="../Results/Butterfly_results/results_final_sp_climate.csv", row.names=FALSE)


# #########################################################################
# ### run the synchrony model for each strategy (specialist/generalist) ###
# #########################################################################
# 
# results_table_spec<-NULL
# for (g in unique(pair_attr$specialism)){
#   
#   spec_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp)-1, data = pair_attr[pair_attr$specialism==g,])
#   summary(spec_model)
#   anova(spec_model)
#   
#   ### save and plot the results ###
#   results_table_spec_temp <- data.frame(summary(spec_model)$coefficients[,1:3],g)
#   results_table_spec<-rbind(results_table_spec,results_table_spec_temp)
# }
# 
# ## change names and add in parameter column ##
# names(results_table_spec) <- c("FCI", "SD", "t","specialism")
# results_table_spec$parameter <- paste(row.names(results_table_spec))
# rownames(results_table_spec) <- 1:nrow(results_table_spec)
# 
# ## take out 3 columns for each species: mean northing, distance and renk_hab_sim
# results_table_strat4 <- NULL
# results_table_strat1 <- results_table_spec[grep("mean_northing", results_table_spec$parameter),]
# results_table_strat2 <- results_table_spec[grep("distance", results_table_spec$parameter),]
# results_table_strat3 <- results_table_spec[grep("renk_hab_sim", results_table_spec$parameter),]
# results_table_strat4 <- rbind(results_table_strat4, results_table_strat1, results_table_strat2, results_table_strat3)
# results_table_spec <- results_table_spec[!results_table_spec$parameter%in%results_table_strat4$parameter,]
# 
# ## change parameter names to year
# results_table_spec$parameter <- rep(1985:2012)
# 
# ### rescale estimate, SD and CI for each species 
# results_final_spec <- NULL
# for (g in unique(results_table_spec$specialism)){
#   
#   results_temp_spec <- results_table_spec[results_table_spec$specialism==g,]  
#   
#   results_temp_spec$rescaled_FCI <- results_temp_spec$FCI*(100/results_temp_spec$FCI[1])
#   results_temp_spec$rescaled_sd <- results_temp_spec$SD*(100/results_temp_spec$FCI[1])
#   results_temp_spec$rescaled_ci <- results_temp_spec$rescaled_sd*1.96
#   
#   results_final_spec <- rbind(results_final_spec, results_temp_spec)
#   
# }
# 
# ## save final results table	##
# write.csv(results_final_spec, file = "../Results/Butterfly_results/results_final_spec.csv", row.names = FALSE) 
# 
