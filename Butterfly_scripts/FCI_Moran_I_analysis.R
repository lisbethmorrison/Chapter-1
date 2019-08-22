###########################################################
## Title: Functional Connectivity and mean Moran's I analysis
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: September 2017
##########################################################


#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script") ####

rm(list=ls()) # clear R

library(ggplot2)

## read in data ##
results_table_sp <- read.csv("../Results/Butterfly_results/results_final_sp.csv", header=TRUE)
moran.results <- read.csv("../Results/Butterfly_results/moran.i.mean.results.csv", header=TRUE)

## remove mid.years 2009-2012  - moran's i data only goes up to 2008 ## 
results_table_sp <- results_table_sp[results_table_sp$parameter<2009,]

## merge moran's i data with fci data by year to create fci_moran ##
fci_moran <- merge(results_table_sp,moran.results,by.x="parameter", by.y="mid.year", all=TRUE)

## change column names ##
colnames(fci_moran)[1] <- "mid_year"
colnames(fci_moran)[3] <- "FCI"
colnames(fci_moran)[11] <- "mean_moran_i"

## remove unused columns 
fci_moran <- fci_moran[-c(4:8,10)]

## plot of FCI against mean moran value ##
fci_moran$sp <- as.factor(fci_moran$sp)

## mode of fci against mean moran's i for each species ##
## loop through each species ## 
fci_moran_results <- NULL
for(i in unique(fci_moran$sp)){
  
  fci_moran_model <- lm(FCI ~ mean_moran_i, data=fci_moran[fci_moran$sp==i,])
  summary(fci_moran_model)
  
  fci_moran_temp <- data.frame(summary(fci_moran_model)$coefficients[,1:4],i)
  fci_moran_results <- rbind(fci_moran_results, fci_moran_temp)
  
} # end i in species.list

## change column names and row names ##
names(fci_moran_results) <- c("Estimate", "SD", "t", "p value", "species")
fci_moran_results$parameter <- paste(row.names(fci_moran_results))

## remove moran i values, just keep intercept ones ## 
results_tab1 <- fci_moran_results[grep("(Intercept)", fci_moran_results$parameter),]
fci_moran_results <- fci_moran_results[!fci_moran_results$parameter%in%results_tab1$parameter,]
rownames(fci_moran_results) <- 1:nrow(fci_moran_results)
fci_moran_results = subset(fci_moran_results, select = -c(parameter))

write.csv(fci_moran_results, file = "../Results/Butterfly_results/fci_moran_results.csv", row.names=FALSE)

##### 4 SIGNIFICANT SPECIES: 48, 23, 46 AND 106 ######

## check model ## 
par(mfrow = c(2, 2))
plot(fci_moran_model)
