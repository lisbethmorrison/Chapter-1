###########################################################
## Title: Functional Connectivity Index and Abundance analysis 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: November 2017
##########################################################

#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script") ####
options(scipen=999)
require(ggplot2)

rm(list=ls()) # clear R

## add data ##
abundance_data <- read.csv("../Data/Bird_sync_data/abundance_data_CBC.csv", header = TRUE) # mean abundance per species
fci_data <- read.csv("../Results/Bird_results/results_final_spp.csv", header = TRUE) # fci data per species
species_names <- read.csv("../Data/BTO_data/woodland_generalist_specialist.csv", header=TRUE)

abundance_data <- subset(abundance_data, species!= "323") # all sites have hab_sim = 1, spp = lesser spotted woodpecker
abundance_data <- subset(abundance_data, species!= "370") # all sites have hab_sim = 1, spp = nightingale 

## add 0.5 to mid.year of abundance data so they can be merged together ##
abundance_data$mid.year <- abundance_data$mid.year+0.5

## merge fci and abundance data by the mid.year ##
fci_abundance <- merge(abundance_data, fci_data, by.x=c("species", "mid.year"), by.y=c("sp", "parameter"), all.x=TRUE)
fci_abundance = subset(fci_abundance, select = -c(3:4,6:8, 10:12))
colnames(fci_abundance)[5] <- "common_name"
colnames(fci_abundance)[2] <- "year"

## log abundance data ##
fci_abundance$log.abundance <- log10(fci_abundance$abundance.index)

## save FCI and abundance data
write.csv(fci_abundance, file = "../Data/Bird_sync_data/FCI_abundance_data.csv", row.names=FALSE)

fci_abundance <- read.csv("../Data/Bird_sync_data/FCI_abundance_data.csv", header=TRUE)
str(fci_abundance)

## linear model of all species together ##
fci_abund_model <- lm(FCI ~ abundance.index, data=fci_abundance)
anova(fci_abund_model)
summary(fci_abund_model)

## log linear model of all species together ## 
fci_abund_model_log <- lm(FCI ~ log.abundance, data=fci_abundance)
anova(fci_abund_model_log)
summary(fci_abund_model_log)

par(mfrow=c(2,2))
plot(fci_abund_model)
plot(fci_abund_model_log) ## looks better 

### linear model of FCI against abundance with loop for each species ###
fci_abund_results <- NULL
for (i in unique(fci_abundance$species)){
  
  fci_abund_model_sp <- lm(FCI ~ log.abundance, data=fci_abundance[fci_abundance$species==i,])
  summary(fci_abund_model_sp)
  
  fci_abund_temp <- data.frame(summary(fci_abund_model_sp)$coefficients[,1:4], summary(fci_abund_model_sp)$r.squared,i)
  fci_abund_results <- rbind(fci_abund_results, fci_abund_temp)
  
} # end i in species

## change column names and row names ##
names(fci_abund_results) <- c("Estimate", "SD", "t", "p value", "R^2", "species_code")
fci_abund_results$parameter <- paste(row.names(fci_abund_results))

## remove log.abundance values, just keep intercept ones ## 
results_tab1 <- fci_abund_results[grep("(Intercept)", fci_abund_results$parameter),]
fci_abund_results <- fci_abund_results[!fci_abund_results$parameter%in%results_tab1$parameter,]
rownames(fci_abund_results) <- 1:nrow(fci_abund_results)
fci_abund_results = subset(fci_abund_results, select = -c(parameter))

## merge with species common name ##
fci_abund_final <- merge(fci_abund_results, species_names, by.x="species_code", by.y="species_code")
fci_abund_final = subset(fci_abund_final, select = -c(generalist, specialist))

###### SPECIES 462, 467, 361, 347, 428, 365, 123, 503, 453, 320, 450, 298, 468, 368, 485 and 470 ARE SIGNIFICANT ######
###### 16 species in total: 7 w/positive relationship, 9 w/negative relationship ######

write.csv(fci_abund_final, file="../Results/Bird_results/fci_abund_results.csv", row.names=FALSE)

