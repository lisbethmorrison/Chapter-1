###########################################################
## Title: Prep BBS BTO data 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: January 2018
##########################################################

rm(list=ls()) # clear R

## add BBS data
bbs_data1 <- read.csv("../Data/BTO_data/BBS_woodland_data1.csv", header=TRUE)
bbs_data2 <- read.csv("../Data/BTO_data/BBS_woodland_data2.csv", header=TRUE)
bbs_data3 <- read.csv("../Data/BTO_data/BBS_woodland_data3.csv", header=TRUE)

## add the two datasets together
bbs_data <- rbind(bbs_data1, bbs_data2, bbs_data3)

summary(bbs_data)
head(bbs_data)
str(bbs_data)

## remove negative values
bbs_data <- bbs_data[bbs_data$TOT>=0,] # removes 3 values

## remove NAs in data from the Count column ##
bbs_data <- na.omit(bbs_data)
summary(bbs_data) # check NAs have been removed 

## add in species codes
cbc_codes <- read.csv("../Data/BTO_data/woodland_generalist_specialist.csv", header=TRUE)
bbs_data <- merge(bbs_data, cbc_codes, by.x="ENGLISH_NAME", by.y="species", all.x=TRUE)
summary(bbs_data)
str(bbs_data)

## create site_code column
bbs_data <- transform(bbs_data, site_code=as.numeric(factor(GRIDREF)))

# take the max value where multiple totals for the same year site species combo (i.e. the max value across the 2 (early and late) site visits per year)
bbs_data <- aggregate(TOT ~ species_code + site_code + YEAR + ENGLISH_NAME + GRIDREF, FUN = max, data = bbs_data)

## save data
write.csv(bbs_data, file="../Data/Bird_sync_data/bbs_woodland_birds.csv", row.names=FALSE)


#### loop through species to estimate number of years of data for each site ####
spp_data <- NULL
spp_list <- unique(bbs_data$species_code) # create unique species list
site_list <- unique(bbs_data$site_code) # create unique site list

for (i in spp_list){ # loop through each species
  temp <- bbs_data[bbs_data$species_code==i,]
  
  num.data <- NULL
  
  for (j in site_list){ # loop through each site
    num.data <- c(num.data, nrow(temp[temp$site_code==j,]))
  } # end j in site_list
  
  prop_data <- data.frame(site_list, num.data)
  prop_data$species <- i
  spp_data <- rbind(spp_data, prop_data)
} # end i in spp_list

#### save spp_data = the number of years of data for each site and each species ####
write.csv(spp_data, file="../Data/Bird_sync_data/site_year_species_BBS.csv", row.names=FALSE)

#### loop through species to estimate number of sites per year ####
spp_data2 <- NULL
spp_list <- unique(bbs_data$species) # create unique species list
year_list <- min(bbs_data$YEAR):max(bbs_data$YEAR)

for (i in spp_list){ # loop through each species
  temp2 <- bbs_data[bbs_data$species_code==i,]
  
  num.data2 <- NULL
  
  for (j in year_list){ # loop through each site
    num.data2 <- c(num.data2, nrow(temp2[temp2$YEAR==j,]))
  } # end j in site_list
  
  prop_data2 <- data.frame(year_list, num.data2)
  prop_data2$species <- i
  spp_data2 <- rbind(spp_data2, prop_data2)
} # end i in spp_list

#### save spp_data2 = the number of sites for each year and each species ####
write.csv(spp_data2, file="../Data/Bird_sync_data/year_site_species_BBS.csv", row.names=FALSE)


### remove sites which have less than 10 years of data
## make site_years with site_code and the minimum year of each site
min_years <- aggregate(bbs_data$YEAR, list(bbs_data$site_code), min) 
colnames(min_years) <- c("site_code", "min_year")
## then create a column with max year
max_years <- aggregate(bbs_data$YEAR, list(bbs_data$site_code), max) 
colnames(max_years) <- c("site_code", "max_year")
## merge the two together
site_years <- merge(min_years, max_years, by="site_code")
## create a column of max-min years
site_years$year_diff <- site_years$max_year - site_years$min_year

## remove sites which have a year_diff of less than 10
site_years <- site_years[!site_years$year_diff<10,]
## this removes 2551 sites

## save file
write.csv(site_years, file="../Data/Bird_sync_data/sites_10_yrs_data.csv", row.names=FALSE)

