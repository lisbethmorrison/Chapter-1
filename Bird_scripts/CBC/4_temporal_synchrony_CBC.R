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
#final_pair_data <- read.csv("../Data/Bird_sync_data/final_pair_data_all_spp_CBC_zeros.csv", header = TRUE) # pair attr synchrony data
final_pair_data <- read.csv("../Data/Bird_sync_data/final_pair_data_all_spp_CBC_no_zeros2_correct.csv", header = TRUE) # pair attr synchrony data
site_data <- read.csv("../Data/BTO_data/pair_attr_mean_north_dist_hab_sim_CBC.csv", header = TRUE) # pair attr site data

## merge to site data for each species ##
pair_attr_1 <- merge(final_pair_data, site_data, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))  # merge the site comparisons in one direction - site a to a, b to b
site_data_reverse <- site_data
names(site_data_reverse)[1:6] <- c("site_b", "site_a", "site_b_EAST", "site_b_NORTH", "site_a_EAST", "site_a_NORTH")
pair_attr_2 <- merge(final_pair_data, site_data_reverse, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))	# merge the site comparisons in the same direction using site_data_reverse
pair_attr <- rbind(pair_attr_1, pair_attr_2) # combine the two datasets

# centre and standardize the distance and northing variables
### should we standardise and centre here?! - standardise by mean and sd?
pair_attr$distance <- pair_attr$distance - mean(na.omit(pair_attr$distance))
pair_attr$distance <- pair_attr$distance/(max(abs(na.omit(pair_attr$distance))))

pair_attr$mean_northing <- pair_attr$mean_northing - mean(na.omit(pair_attr$mean_northing))
pair_attr$mean_northing <- pair_attr$mean_northing/(max(abs(na.omit(pair_attr$mean_northing))))

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

pair_attr <- subset(pair_attr, spp!= "323") # all sites have hab_sim = 1, spp = lesser spotted woodpecker
pair_attr <- subset(pair_attr, spp!= "370") # all sites have hab_sim = 1, spp = nightingale 

length(unique(pair_attr$spp)) # 26 species

## merge in trait data (specialism and dispersal distance)
specialism <- read.csv("../Data/BTO_data/woodland_generalist_specialist.csv", header=TRUE)
pair_attr <- merge(pair_attr, specialism, by.x="spp", by.y="species_code")
length(unique(pair_attr$spp)) ## 26 species

write.csv(pair_attr, file = "../Data/Bird_sync_data/pair_attr_CBC_no_zeros2_correct.csv", row.names = FALSE) # save pair_attr file (31 spp)


########## READ IN PAIR ATTRIBUTE FILE ################
pair_attr <- read.csv("../Data/Bird_sync_data/pair_attr_CBC_no_zeros2_correct.csv", header=TRUE)
## climate data
final_pair_data_temp <- read.csv("../Data/MetOffice_data/final_pair_data_mean_temp_CBC.csv", header=TRUE)  

## merge in climate data
## just summer and autumn temperature (this is the only significant change in climate data)
final_pair_data_sum_temp <- final_pair_data_temp[final_pair_data_temp$season=="c",] ## summer
final_pair_data_aut_temp <- final_pair_data_temp[final_pair_data_temp$season=="d",] ## autumn
names(final_pair_data_sum_temp)[3] <- "lag0_summer"
names(final_pair_data_aut_temp)[3] <- "lag0_autumn"

## remove some columns (season, start and end year)
final_pair_data_sum_temp <- subset(final_pair_data_sum_temp, select = -c(5:7))
final_pair_data_aut_temp <- subset(final_pair_data_aut_temp, select = -c(5:7))

### merge in summer temp
pair_attr_1 <- merge(pair_attr, final_pair_data_sum_temp, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"))  # merge the site comparisons in one direction - site a to a, b to b
summer_reverse <- final_pair_data_sum_temp
names(summer_reverse)[1:2] <- c("site2", "site1")
pair_attr_2 <- merge(pair_attr, summer_reverse, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"))	# merge the site comparisons in the same direction using site_data_reverse
pair_attr <- rbind(pair_attr_1, pair_attr_2) # combine the two datasets
length(unique(pair_attr$spp))## 29 species
length(unique(pair_attr$site1)) # 109
length(unique(pair_attr$site2)) # 109 
## merge in autumn temp
pair_attr_1 <- merge(pair_attr, final_pair_data_aut_temp, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"))  # merge the site comparisons in one direction - site a to a, b to b
autumn_reverse <- final_pair_data_aut_temp
names(autumn_reverse)[1:2] <- c("site2", "site1")
pair_attr_2 <- merge(pair_attr, autumn_reverse, by.x=c("site1", "site2", "mid.year"), by.y=c("site1", "site2", "mid.year"))	# merge the site comparisons in the same direction using site_data_reverse
pair_attr <- rbind(pair_attr_1, pair_attr_2) # combine the two datasets
length(unique(pair_attr$spp))## 29 species
length(unique(pair_attr$site1)) # 109
length(unique(pair_attr$site2)) # 109 

##
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$end.year <- as.factor(pair_attr$end.year)
pair_attr$mid.year <- as.factor(pair_attr$mid.year)
str(pair_attr)

################################################
##  model to produce one line for all species ##
################################################

### run model with spring rainfall 
all_spp_model_temp <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + lag0_summer + lag0_autumn + (1|pair.id) + (1|spp)-1, data = pair_attr)
summary(all_spp_model_temp) 
anova(all_spp_model_temp) ## summer and autumn temp are non-significant

pop_temp_model <- data.frame(summary(all_spp_model_temp)$coefficients[,1:3])
names(pop_temp_model) <- c("FCI", "SD", "t")
pop_temp_model$parameter <- paste(row.names(pop_temp_model))
rownames(pop_temp_model) <- 1:nrow(pop_temp_model)

results_tab5 <- NULL
results_tab1 <- pop_temp_model[grep("mean_northing", pop_temp_model$parameter),]
results_tab2 <- pop_temp_model[grep("distance", pop_temp_model$parameter),]
results_tab3 <- pop_temp_model[grep("hab_sim", pop_temp_model$parameter),]
results_tab4 <- pop_temp_model[grep("lag0_summer", pop_temp_model$parameter),]
results_tab5 <- pop_temp_model[grep("lag0_autumn", pop_temp_model$parameter),]
results_tab6 <- rbind(results_tab4, results_tab1, results_tab2, results_tab3, results_tab4, results_tab5)
pop_temp_model <- pop_temp_model[!pop_temp_model$parameter%in%results_tab6$parameter,]

## change parameter names to year
pop_temp_model$parameter <- rep(1985:1996)

### rescale estimate, SD and CI ### 
pop_temp_model$rescaled_FCI <- pop_temp_model$FCI*(100/pop_temp_model$FCI[1])
pop_temp_model$rescaled_sd <- pop_temp_model$SD*(100/pop_temp_model$FCI[1])
pop_temp_model$rescaled_ci <- pop_temp_model$rescaled_sd*1.96

## save final results table ##
write.csv(pop_temp_model, file = "../Results/Bird_results/results_final_all_spp_temp_CBC.csv", row.names=FALSE)

## graph
FCI_plot_CBC_temp <- ggplot(pop_temp_model, aes(x = parameter, y = rescaled_FCI)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_FCI - rescaled_sd, ymax = rescaled_FCI + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Population synchrony") +
  #scale_y_continuous(breaks=seq(40,160,10)) +
  scale_x_continuous(breaks=seq(1985,1996,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
FCI_plot_CBC_temp
ggsave("../Graphs/Connectivity_plots/FCI_plot_temp_CBC.png", plot = FCI_plot_CBC_temp, width=7, height=5)

############################## merge in species data (family) ###############################
## subset family and genus info
## family/genus data
species_traits <- read.csv("../Data/BTO_data/woodland_generalist_specialist.csv", header=TRUE)
species_traits <- species_traits[,c(2,3,4)]
pair_attr <- merge(pair_attr, species_traits, by.x="spp", by.y="species_code", all=FALSE)
pair_attr <- droplevels(pair_attr)
length(unique(pair_attr$spp)) # 26 species
summary(pair_attr)

### run model with family to test for a relationship with synchrony
all_spp_model_family <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + family + (1|pair.id) + (1|spp), data = pair_attr)
summary(all_spp_model_family) 
anova(all_spp_model_family) ## family significant (p=0.042)

all_spp_model_genus <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + genus + (1|pair.id) + (1|spp), data = pair_attr)
summary(all_spp_model_genus) 
anova(all_spp_model_genus) ## genus non-significant (p=0.21)

## model with intercept
start_time <- Sys.time()
all_spp_model_int <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + (1|pair.id) + (1|spp), data = pair_attr)
summary(all_spp_model_int)
end_time <- Sys.time()
end_time - start_time # 20.5 seconds

### save results for northing, distance and hab sim
fixed_results <- data.frame(summary(all_spp_model_int)$coefficients[,1:5])
fixed_results$parameter <- paste(row.names(fixed_results))
rownames(fixed_results) <- 1:nrow(fixed_results)
## remove mid.year rows
fixed_results <- fixed_results[-c(1,5:15),]
## save results
write.csv(fixed_results, file = "../Results/Model_outputs/fixed_effect_results_cbc_no_zeros2_correct.csv", row.names=FALSE)

## model without intercept
all_spp_model <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + (1|pair.id) + (1|spp)-1, data = pair_attr)
summary(all_spp_model) 
anova(all_spp_model)

### check model fit ###
plot(all_spp_model)
qqnorm(resid(all_spp_model))
qqline(resid(all_spp_model)) # not too bad

## save results ## 
results_table_all_spp <- data.frame(summary(all_spp_model)$coefficients[,1:3])

## change names and add in parameter column ##
names(results_table_all_spp) <- c("FCI", "SD", "t")
results_table_all_spp$parameter <- paste(row.names(results_table_all_spp))
rownames(results_table_all_spp) <- 1:nrow(results_table_all_spp)

## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_tab4 <- NULL
results_tab1 <- results_table_all_spp[grep("mean_northing", results_table_all_spp$parameter),]
results_tab2 <- results_table_all_spp[grep("distance", results_table_all_spp$parameter),]
results_tab3 <- results_table_all_spp[grep("hab_sim", results_table_all_spp$parameter),]
results_tab4 <- rbind(results_tab4, results_tab1, results_tab2, results_tab3)
results_table_all_spp <- results_table_all_spp[!results_table_all_spp$parameter%in%results_tab4$parameter,]

## change parameter names to year
results_table_all_spp$parameter <- rep(1985:1996)

### rescale estimate, SD and CI ### 
results_table_all_spp$rescaled_FCI <- results_table_all_spp$FCI*(100/results_table_all_spp$FCI[1])
results_table_all_spp$rescaled_sd <- results_table_all_spp$SD*(100/results_table_all_spp$FCI[1])
results_table_all_spp$rescaled_ci <- results_table_all_spp$rescaled_sd*1.96

## save final results table ##
write.csv(results_table_all_spp, file = "../Results/Bird_results/results_final_all_spp_CBC_no_zeros2_correct.csv", row.names=FALSE)

## same model but without covariates
all_spp_model2 <- lmer(lag0 ~ mid.year + (1|pair.id) + (1|spp), data = pair_attr)
summary(all_spp_model2) 
anova(all_spp_model2)

## AIc between the full and minimal model
AIC(all_spp_model, all_spp_model2)
## all_spp_model = 60373.99
## all_spp_model2 = 60347.78
## minimal model is better (by 26.21)

##################################################
##### model to produce one line per species #####
#################################################

results_table_sp<-NULL
for (i in unique(pair_attr$spp)){
  
  species_model <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + (1|pair.id)-1, data = pair_attr[pair_attr$spp==i,])
  print(i)
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
length(unique(results_final_sp$sp)) # 29 species
## save species final results table ##
write.csv(results_final_sp, file = "../Results/Bird_results/results_final_spp_CBC_no_zeros2.csv", row.names=FALSE)

#########################################################
### run the synchrony model for generalist/specialist ###
#########################################################

## merge with generalist/specialist data ##
str(pair_attr)
pair_attr$specialism <- as.factor(pair_attr$specialism)

results_table_spec<-NULL
for (j in unique(pair_attr$specialism)){
  
  specialist_model <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + (1|pair.id) + (1|spp)-1, data = pair_attr[pair_attr$specialism==j,])
  summary(specialist_model)
  anova(specialist_model)
  
  ### save and plot the results ###
  results_table_spec_temp <- data.frame(summary(specialist_model)$coefficients[,1:3],j)
  results_table_spec<-rbind(results_table_spec,results_table_spec_temp)
}

## change names and add in parameter column ##
names(results_table_spec) <- c("FCI", "SD", "t","Strategy")
results_table_spec$parameter <- paste(row.names(results_table_spec))
rownames(results_table_spec) <- 1:nrow(results_table_spec)

## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_table_spec4 <- NULL
results_table_spec1 <- results_table_spec[grep("mean_northing", results_table_spec$parameter),]
results_table_spec2 <- results_table_spec[grep("distance", results_table_spec$parameter),]
results_table_spec3 <- results_table_spec[grep("hab_sim", results_table_spec$parameter),]
results_table_spec4 <- rbind(results_table_spec4, results_table_spec1, results_table_spec2, results_table_spec3)
results_table_spec <- results_table_spec[!results_table_spec$parameter%in%results_table_spec4$parameter,]

## change parameter names to year
results_table_spec$parameter <- rep(1985:1996)

## rescale FCI, SD and CI
results_final_spec <- NULL
for (k in unique(results_table_spec$Strategy)){
  
  results_temp_spec <- results_table_spec[results_table_spec$Strategy==k,]  
  
  results_temp_spec$rescaled_FCI <- results_temp_spec$FCI*(100/results_temp_spec$FCI[1])
  results_temp_spec$rescaled_sd <- results_temp_spec$SD*(100/results_temp_spec$FCI[1])
  results_temp_spec$rescaled_ci <- results_temp_spec$rescaled_sd*1.96
  
  results_final_spec <- rbind(results_final_spec, results_temp_spec)
  
}

## save species specialism final results
write.csv(results_final_spec, file = "../Results/Bird_results/results_final_spec_CBC_no_zeros2.csv", row.names=FALSE)

