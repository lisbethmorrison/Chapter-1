###########################################################
## Title: Calculate BTO site data info  
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: November 2017
##########################################################

rm(list=ls()) # clear R

load("Bird_scripts/gr_let2num.Rdata") # load function to convert gridrefs to easting and northing

####################################
#### add woodland bird cbc data ####
####################################

woodland_cbc <- read.csv("../Data/Bird_sync_data/cbc_woodland_birds.csv", header=TRUE)

head(woodland_cbc)
summary(woodland_cbc)

##################################################
#### change grid refs to easting and northing ####
##################################################

#### use function "gr_let2num" from R script "CONVERT_gridrefs_to_east_norths_Colin_H" ####
easting_northing_values <- gr_let2num(woodland_cbc$Gridref)
woodland_cbc <- cbind(woodland_cbc, easting_northing_values)

#### create unique site list #####
site_data <- woodland_cbc[,c("Pcode", "Gridref", "EASTING", "NORTHING")] # create temp dataframe with only site number, easting and northing
site_data <- unique(site_data) # take out duplicated sites, left with 411 sites
names(site_data)[1:4] <- c("Site_code", "Gridref", "Easting", "Northing")

#### save site list info ####
write.csv(site_data, file = "../Data/BTO_data/site_list_east_north.csv", row.names=FALSE)

## save woodland_cbc data which has 411 sites
write.csv(woodland_cbc, file = "../Data/Bird_sync_data/cbc_woodland_birds.csv", row.names=FALSE)


###############################################################
### merge woodland_cbc gridref with hab_data gridref ###
## therefore hab_data only includes sites we are looking at ##
# how many sites do we have attribute data for
# hab_data <- read.table("../Data/Land_cover_data/BBS.CBC.UKBMS.complete.landcover.data.all.scales.soil.dem.configuration.txt", header=TRUE)
# 
# ## subset BBS data and 500m buffer only
# hab_data <- hab_data[hab_data$Surv=="CBC",]
# hab_data <- hab_data[hab_data$buffer=="500",]
# 
# length(unique(woodland_cbc$Gridref)[unique(site_data$Easting)%in%unique(hab_data$POINT_X)]) # 39 sites that have site attribute data
# 
# correct_sites <- merge(site_data, hab_data, by.x="Easting", by.y="POINT_X")
#### not enough site attribute data to calculate renk hab sim for CBC data
## use binary code instead
############################################################

######################################################
#### create list of all pairwise site comparisons ####
######################################################
site.list <- unique(as.character(site_data$Site_code))
num.ts <- length(site.list)

pair.list <- NULL

pair.list <- t(combn(site.list,2))
pair_attr <- data.frame(pair.list)
nrow(pair_attr) ## 84255 rows 
names(pair_attr) <- c("site_a", "site_b")

## add xy data for site_a and site_b
pair_attr <- merge(pair_attr, site_data, by.x="site_a", by.y="Site_code")
names(pair_attr)[4:5] <- c("site_a_EAST", "site_a_NORTH") # rename them
pair_attr <- merge(pair_attr, site_data, by.x="site_b", by.y="Site_code")
names(pair_attr)[7:8] <- c("site_b_EAST", "site_b_NORTH") # rename them
## remove grid ref columns
pair_attr <- pair_attr[-c(3,6)]
# sort out column order #
pair_attr <- pair_attr[,c("site_a", "site_b", "site_a_EAST", "site_a_NORTH", "site_b_EAST", "site_b_NORTH")]

# ## Add habitat code info ##
# ## create site code vector which has site number and habitat code ##
# site_codes <- woodland_cbc[,-c(1,3:5,7:12)]
# site_codes <- unique(site_codes) ## 411 rows
# 
# pair_attr <- merge(pair_attr, site_codes, by.x="site_a", by.y="Pcode") # merge by site_a
# names(pair_attr)[7] <- "site_a_hab_code" # rename the column
# pair_attr <- merge(pair_attr, site_codes, by.x="site_b", by.y="Pcode") # merge by site_b
# names(pair_attr)[8] <- "site_b_hab_code" # rename the column
# 
# ## sort column order
# pair_attr <- pair_attr[,c("site_a", "site_b", "site_a_EAST", "site_a_NORTH", "site_b_EAST", "site_b_NORTH", "site_a_hab_code", "site_b_hab_code")]
# 
# #### create binary habitat similarity column ####
# pair_attr$hab_sim <- ifelse(pair_attr$site_a_hab_code==pair_attr$site_b_hab_code, 1, 0) # matching habitat type = 1, not matching = 0

# take average northing between the pairwise comparisons #
pair_attr$mean_northing <-rowMeans(subset(pair_attr, select = c(site_a_NORTH, site_b_NORTH)), na.rm = TRUE)

### Estimate distance between sites ###
pair_attr$distance <- sqrt((pair_attr$site_a_NORTH-pair_attr$site_b_NORTH)^2 + (pair_attr$site_a_EAST-pair_attr$site_b_EAST)^2) ## pythagoras equation
pair_attr$distance <- pair_attr$distance/1000
## this gives distance in km

summary(pair_attr)

### save pair_attr ###
#write.csv(pair_attr, file="../Data/BTO_data/pair_attr_mean_north_dist_CBC.csv", row.names=FALSE)

###### Calculate habitat similarity index

## read in hab_data
hab_data <- read.table("../Data/Land_cover_data/BBS.CBC.UKBMS.complete.landcover.data.all.scales.soil.dem.configuration.txt", header=TRUE)
final_data <- read.csv("../Data/Bird_sync_data/final_pair_data_all_spp_CBC.csv", header = TRUE)

## subset BBS data and 500m buffer only
hab_data <- hab_data[hab_data$Surv=="CBC",]
hab_data <- hab_data[hab_data$buffer=="500",]

## remove columns not needed
hab_data <- hab_data[-c(2:3,6:16,30:85)]

## create total_count column
hab_data$total_count <- rowSums(hab_data[4:16])

# divide all columns by total land area
habs <- c("A","BgRo","Br","BW","C","CW","F","G","H","M","S","R","UG")

for (i in habs){
  hab_data[,i] <- hab_data[,i]/hab_data[,"total_count"]
}

## merge hab_data with bbs_site_data
## to have site numbers in hab_data
## this removes some sites which have habitat data
## sites removed are those which have less than 10 years of data
hab_data <- merge(hab_data, site_data, by.x="siteno.gref", by.y="Gridref", all=FALSE) ## 23 sites with habitat data
## remove duplicated easting and northing columns
hab_data <- hab_data[-c(2:3)] ## only 23 sites which have habitat data...
## re-order columns
# hab_data <- hab_data[,c(1,15,16,17,18,2:14)]
# 
# pair_attr$renk_hab_sim <- -9999
# row.names(hab_data) = hab_data$site_code
# 
# all_mins = pmin(hab_data[as.character(pair_attr$site_a),6:18], hab_data[as.character(pair_attr$site_b),6:18])
# pair_attr[,"renk_hab_sim"] = rowSums(all_mins)
# 
# summary(pair_attr)
# ### 513684 NA's => result of 152 sites without habitat data
# 
# ## remove NA's
# pair_attr <- na.omit(pair_attr) ## 5456556 rows :) 
# 
# ## save file
# ## this is site data for sites with at least 10 years of data and those with habitat data to calc renk_hab_sim
# write.csv(pair_attr, file="../Data/BTO_data/pair_attr_mean_north_dist_sim_BBS.csv", row.names=FALSE)
