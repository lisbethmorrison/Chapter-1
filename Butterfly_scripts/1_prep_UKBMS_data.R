###########################################################
## Title: Prep UKBMS Data
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: August 2017
##########################################################

#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script") ####

rm(list=ls()) # clear R

### add data
b_data <- read.csv("../../Data/UKBMS_data/ukbms_sindex.csv", sep = ",", header = TRUE) # UKBMS species data
species_codes <- read.csv("../../Data/UKBMS_data/ukbms_species.csv", header = TRUE) # UKBMS species code data

## merge b_data with species code data to exclude moths ##
b_data <- b_data[b_data$SPECIES%in%species_codes$SPECIES,]

### data cleaning ###

# remove negative SINDEX values
b_data <- b_data[b_data$SINDEX>=0,]

# drop records based on BROOD = 0, take total annual counts
b_data <- b_data[b_data$BROOD == 0,]

# remove pilot years (>1975)
b_data <- b_data[b_data$YEAR>1975,]

# drop egg counts
b_data <- b_data[b_data$SITE<10000,]

# take the max value where multiple sindex for the same year site species combo
b_data <- aggregate(SINDEX ~ SPECIES + BROOD + SITE + YEAR, FUN = max, data = b_data)

### loop through species - estimate number of sites for each year ###
spp_data <- NULL
spp_list <- unique(b_data$SPECIES) # species to loop through (ALL SPECIES = 84)
year.list <- min(b_data$YEAR):max(b_data$YEAR)

for (i in spp_list){ # loop through all species
  temp <- b_data[b_data$SPECIES==i,]

num.data <- NULL
  
  for (j in year.list) { # loop through year - estimate the number of survey data for each year #
  num.data <- c(num.data, nrow(temp[temp$YEAR==j,]))
  } # end j in year.list

prop_data <- data.frame(year.list, num.data)
prop_data$SPECIES <- i   
spp_data<- rbind(spp_data, prop_data)
} # end i in spp_list
    

### loop through species - estimate years of data for each site ###
spp_data2 <- NULL
spp_list <- unique(b_data$SPECIES) # species to loop through (ALL SPECIES = 84)
site.list <- unique(b_data$SITE)

for (i in spp_list){ # loop through all species
  temp2 <- b_data[b_data$SPECIES==i,]
  
  num.data2 <- NULL
  
  for (j in site.list) { # loop through site - estimate the number of years of data for each site #
    num.data2 <- c(num.data2, nrow(temp2[temp2$SITE==j,]))
  } # end j in year.list
  
  prop_data2 <- data.frame(site.list, num.data2)
  prop_data2$SPECIES <- i   
  spp_data2<- rbind(spp_data2, prop_data2)
} # end i in spp_list

write.csv(spp_data, file = "../../Data/Butterfly_sync_data/years_data_spp.csv", row.names = FALSE)
write.csv(spp_data2, file = "../../Data/Butterfly_sync_data/site_data_spp.csv", row.names = FALSE)
write.csv(b_data, file = "../../Data/Butterfly_sync_data/b_data.csv", row.names = FALSE)

  
  
  
