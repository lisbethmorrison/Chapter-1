###########################################################
## Title: Moving window climate synchrony calculation  
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: March 2019
##########################################################

rm(list=ls()) # clear R
library(dplyr)

### add all files 
file_names <- list.files("../Data/MetOffice_data/5km_gridded_mean_temperature", full.names=TRUE) # where all the files are stored
mean_temp <- do.call(cbind,lapply(file_names,read.csv)) # create new dataframe with all species rbind together

head(mean_temp)
summary(mean_temp)
## remove years <1980
