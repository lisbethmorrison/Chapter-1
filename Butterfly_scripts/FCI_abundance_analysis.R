###########################################################
## Title: Functional Connectivity Index and Abundance analysis 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: September 2017
##########################################################

#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script") ####
options(scipen=999)
require(ggplot2)

rm(list=ls()) # clear R

## add data ##
abundance_data <- read.csv("../Data/Butterfly_sync_data/abundance_data.csv", header = TRUE) # mean abundance per species
fci_data <- read.csv("../Results/Butterfly_results/results_final_sp.csv", header = TRUE) # fci data per species

## remove 2 species ##
abundance_data <- abundance_data[!abundance_data$species==100,] # DROP SPECIES 100 (Small White) as no mobility data
abundance_data <- abundance_data[!abundance_data$species==121,] # DROP SPECIES 121 as its a duplicate of small skipper

## add 0.5 to mid.year of abundance data so they can be merged together ##
abundance_data$mid.year <- abundance_data$mid.year+0.5

## merge fci and abundance data by the mid.year ##
fci_abundance <- merge(abundance_data, fci_data, by.x=c("species", "mid.year"), by.y=c("sp", "parameter"), all.x=TRUE)
fci_abundance = subset(fci_abundance, select = -c(3:4,6:8,10:12))
fci_abundance$species <- as.factor(fci_abundance$species)

## log abundance data ##
fci_abundance$log.abundance <- log10(fci_abundance$abundance.index)

colnames(fci_abundance)[2] <- "year"

## save FCI and abundance data
write.csv(fci_abundance, file = "../Data/Butterfly_sync_data/FCI_abundance_data.csv", row.names=FALSE)

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
plot(fci_abund_model_log) ## much better

### linear model of FCI against abundance with loop for each species ###
fci_abund_results <- NULL
for (i in unique(fci_abundance$species)){
  
  fci_abund_model_sp <- lm(FCI ~ log.abundance, data=fci_abundance[fci_abundance$species==i,])
  summary(fci_abund_model_sp)

  fci_abund_temp <- data.frame(summary(fci_abund_model_sp)$coefficients[,1:4], summary(fci_abund_model_sp)$r.squared,i)
  fci_abund_results <- rbind(fci_abund_results, fci_abund_temp)
  
} # end i in species

## change column names and row names ##
names(fci_abund_results) <- c("Estimate", "SD", "t", "p value", "R^2", "species")
fci_abund_results$parameter <- paste(row.names(fci_abund_results))

## remove log.abundance values, just keep intercept ones ## 
results_tab1 <- fci_abund_results[grep("(Intercept)", fci_abund_results$parameter),]
fci_abund_results <- fci_abund_results[!fci_abund_results$parameter%in%results_tab1$parameter,]
rownames(fci_abund_results) <- 1:nrow(fci_abund_results)
fci_abund_results = subset(fci_abund_results, select = -c(parameter))

## merge results with fci_abundance common name data ##
fci_abund_final <- merge(fci_abundance, fci_abund_results, by.x="species", by.y="species")
fci_abund_final = fci_abund_final[!duplicated(fci_abund_final$species),]
fci_abund_final = subset(fci_abund_final, select = -c(abundance.index))

## reorder columns ## 
fci_abund_final <- fci_abund_final[,c(1,4,2,5,3,6,7,8,9,10)]

###### SPECIES 68, 120, 94, 123, 2, 29, 17, 98, 27, 18, 75, 64, 4, 78, 20, 116, 99 and 23 ARE SIGNIFICANT ######
###### 18 species in total: 12 w/positive relationship, 6 w/negative relationship ######

write.csv(fci_abund_final, file="../Results/Butterfly_results/fci_abund_results.csv", row.names=FALSE)

###### proportion test #######
35*0.05
## = 1.75, rounds up to 2 ##
prop.test(x=c(18,2), n=c(35,35))
# p = 0.00007229

par(mfrow=c(2,2))
plot(fci_abund_model_sp) ## doesn't look great


############## PLOTS 2 SPECIES ###############

######## PLOT THE TWO SIGNIFICANT SPECIES: WHITE ADMIRAL AND PEARL-BORDERED FRITILLARY #########
## White Admiral ##
spp18 <- ggplot(fci_abundance, aes(x=log.abundance, y=FCI)) +
  geom_point(data=subset(fci_abundance,species=="18"), size=6) +
  geom_smooth(data=subset(fci_abundance,species=="18"),se=FALSE,method=lm,level=0.95,col="red") +
  xlab(expression(paste(log[10], sep=" ", "abundance index"))) +
  ylab(expression("Functional connectivity index")) +
  scale_y_continuous(breaks=seq(-0.2,0.3,0.1)) +
  scale_x_continuous(breaks=seq(0.5,1.5,0.05)) +
  theme(axis.title.x = element_text(size=20)) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text = element_text(size=25))
spp18

spp98 <- ggplot(fci_abundance, aes(x=log.abundance, y=FCI)) +
  geom_point(data=subset(fci_abundance,species=="98"), size=6) +
  geom_smooth(data=subset(fci_abundance,species=="98"),se=FALSE,method=lm,level=0.95,col="red") +
  xlab(expression(paste(log[10], sep=" ", "abundance index"))) +
  ylab(expression("Functional connectivity index")) +
  scale_y_continuous(breaks=seq(-0.2,0.8,0.1)) +
  scale_x_continuous(breaks=seq(1.5,1.85,0.05)) +
  theme(axis.title.x = element_text(size=20)) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text = element_text(size=25))
spp98


######################################
######################################
##### FCI ABUNDANCE WOODLAND SPP #####
######################################
######################################

## create new woodland fci_abundance data frame 
woodland_fci_abund <- fci_abundance[fci_abundance$HABITAT=="Woodland",]
## remove habitat column as not needed anymore
woodland_fci_abund <- subset(woodland_fci_abund, select=-c(HABITAT))
summary(woodland_fci_abund)

## linear model of all species together ##
wood_fci_abund_model <- lm(FCI ~ abundance.index, data=woodland_fci_abund)
anova(wood_fci_abund_model)
summary(wood_fci_abund_model)

## log linear model of all species together ## 
wood_fci_abund_model_log <- lm(FCI ~ log.abundance, data=woodland_fci_abund)
anova(wood_fci_abund_model_log)
summary(wood_fci_abund_model_log)

## check model ##
par(mfrow=c(2,2))
plot(wood_fci_abund_model)
plot(wood_fci_abund_model_log) # looks better than the above one 

### linear model of FCI against abundance with loop for each species ###
wood_fci_abund_results <- NULL
for (i in unique(woodland_fci_abund$species)){
  
  wood_fci_abund_model_sp <- lm(FCI ~ log.abundance, data=woodland_fci_abund[woodland_fci_abund$species==i,])
  summary(wood_fci_abund_model_sp)
  
  fci_abund_temp <- data.frame(summary(wood_fci_abund_model_sp)$coefficients[,1:4], summary(wood_fci_abund_model_sp)$r.squared,i)
  wood_fci_abund_results <- rbind(wood_fci_abund_results, fci_abund_temp)
  
} # end i in species

## change column names and row names ##
names(wood_fci_abund_results) <- c("Estimate", "SD", "t", "p value", "R^2", "species")
wood_fci_abund_results$parameter <- paste(row.names(wood_fci_abund_results))

## remove log.abundance values, just keep intercept ones ## 
results_tab1 <- wood_fci_abund_results[grep("(Intercept)", wood_fci_abund_results$parameter),]
wood_fci_abund_results <- wood_fci_abund_results[!wood_fci_abund_results$parameter%in%results_tab1$parameter,]
rownames(wood_fci_abund_results) <- 1:nrow(wood_fci_abund_results)
wood_fci_abund_results = subset(wood_fci_abund_results, select = -c(parameter))

## merge results with fci_abundance common name data ##
wood_fci_abund_final <- merge(fci_abundance, wood_fci_abund_results, by.x="species", by.y="species")
wood_fci_abund_final = wood_fci_abund_final[!duplicated(wood_fci_abund_final$species),]
wood_fci_abund_final = subset(wood_fci_abund_final, select = -c(2:4,6:7))

## save final results
write.csv(wood_fci_abund_final, file = "../Results/Butterfly_results/wood_fci_abund_results.csv", row.names=FALSE)

