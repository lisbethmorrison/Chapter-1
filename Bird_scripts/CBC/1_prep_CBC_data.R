###########################################################
## Title: Prep CBC BTO data 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: November 2017
##########################################################


#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script/Bird_scripts") ####

rm(list=ls()) # clear R

## add BTO data ##
cbc_data <- read.csv("../Data/BTO_data/wCBC61-08.csv", header=TRUE)

## add generalist/specialist woodland bird data
woodland_spp <- read.csv("../Data/BTO_data/woodland_generalist_specialist.csv", header=TRUE)

## add habitat code data 
hab_codes <- read.csv("../Data/BTO_data/CBC_woodland_codes.csv", header=TRUE)

#### merge data together by species code, so only left with woodland species ####
woodland_cbc <- merge(woodland_spp, cbc_data, by.x="species_code", by.y="bou")

### remove unnecessary columns for the synchrony analysis ###
woodland_cbc <- subset(woodland_cbc, select = -c(Spname, Note))

## remove NAs in data from the Count column ##
woodland_cbc <- na.omit(woodland_cbc)
summary(woodland_cbc) # check NAs have been removed 

# take the max value where multiple sindex for the same year site species combo
woodland_cbc <- aggregate(Count ~ species_code + Pcode + Year + Gridref + hab, FUN = max, data = woodland_cbc)

####################################################################
##### remove incorrect Gridrefs (keep those with 1km grid ref) #####
####################################################################

woodland_cbc$Gridref <- as.character(woodland_cbc$Gridref) # first change Gridref into a character
woodland_cbc$Grid <- substring(woodland_cbc$Gridref, 1, 2) # create a Grid column with the first 2 characters
woodland_cbc <- woodland_cbc[-grep("\\.", woodland_cbc$Grid),] # remove rows with full stops in Grid column (so only have gridrefs with 2 letters at the beginning, not one) 
woodland_cbc$Gridref <- gsub("\\.", "", woodland_cbc$Gridref) # remove full stops in GridRef data (i.e. those with full stops in middle of gridref) 
woodland_cbc <- woodland_cbc[nchar(woodland_cbc$Gridref) >= 6, ] # remove gridrefs with less than 6 characters
woodland_cbc$Gridref <- substring(woodland_cbc$Gridref, 1, 6) # remove the last 2 characters from each gridref, so left with 1km gridrefs
woodland_cbc <- subset(woodland_cbc, select = -(Grid)) # remove Grid column 

## merge with habitat code info ##
woodland_cbc <- merge(woodland_cbc, hab_codes, by.x="hab", by.y="Habitat.Code")

## remove columns not needed, and just leave H2 habitat code (0,1,4,7,9)
woodland_cbc <- subset(woodland_cbc, select= -c(hab, H1, H1.Description, H2.Notes, H3, H3.Description))

## add in generalist/specialist data
woodland_gen_spec <- read.csv("../Data/BTO_data/woodland_generalist_specialist.csv", header=TRUE)
woodland_cbc <- merge(woodland_cbc, woodland_gen_spec, by="species_code")

#### save file ####
write.csv(woodland_cbc, file="../Data/Bird_sync_data/cbc_woodland_birds.csv", row.names=FALSE)

#### loop through species to estimate number of years of data for each site ####
spp_data <- NULL
spp_list <- unique(woodland_cbc$species) # create unique species list
site_list <- unique(woodland_cbc$Pcode) # create unique site list

for (i in spp_list){ # loop through each species
  temp <- woodland_cbc[woodland_cbc$species==i,]
  
  num.data <- NULL
  
  for (j in site_list){ # loop through each site
    num.data <- c(num.data, nrow(temp[temp$Pcode==j,]))
  } # end j in site_list
  
  prop_data <- data.frame(site_list, num.data)
  prop_data$species <- i
  spp_data <- rbind(spp_data, prop_data)
} # end i in spp_list

#### save spp_data = the number of years of data for each site and each species ####
write.csv(spp_data, file="../Data/Bird_sync_data/site_year_species_CBC.csv", row.names=FALSE)

#### loop through species to estimate number of sites per year ####
spp_data2 <- NULL
spp_list <- unique(woodland_cbc$species) # create unique species list
year_list <- min(woodland_cbc$Year):max(woodland_cbc$Year)

for (i in spp_list){ # loop through each species
  temp2 <- woodland_cbc[woodland_cbc$species==i,]
  
  num.data2 <- NULL
  
  for (j in year_list){ # loop through each site
    num.data2 <- c(num.data2, nrow(temp2[temp2$Year==j,]))
  } # end j in site_list
  
  prop_data2 <- data.frame(year_list, num.data2)
  prop_data2$species <- i
  spp_data2 <- rbind(spp_data2, prop_data2)
} # end i in spp_list

#### save spp_data2 = the number of sites for each year and each species ####
write.csv(spp_data2, file="../Data/Bird_sync_data/year_site_species_CBC.csv", row.names=FALSE)
