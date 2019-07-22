###########################################################
## Title: Calculate BBS site data info => change loops into functions
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: February 2018
##########################################################

rm(list=ls()) # clear R

load("Bird_scripts/gr_let2num.Rdata") # load function to convert gridrefs to easting and northing

#### add woodland bird bbs data ####
woodland_bbs <- read.csv("../Data/Bird_sync_data/bbs_woodland_birds.csv", header=TRUE)

head(woodland_bbs)
summary(woodland_bbs)

##################################################
#### change grid refs to easting and northing ####
##################################################

## create bbs_site_data with only site code and gridref columns
bbs_site_data <- woodland_bbs[,c("site_code", "GRIDREF")]
bbs_site_data <- unique(bbs_site_data)
## there are irish Gridrefs (start with I*), these cannot be coverted so need taken out
bbs_site_data <- bbs_site_data[!grepl("I", bbs_site_data$GRIDREF),]

#### use function "gr_let2num" from R script "CONVERT_gridrefs_to_east_norths_Colin_H" ####
bbs_east_north <- gr_let2num(bbs_site_data$GRIDREF)
bbs_site_data <- cbind(bbs_site_data, bbs_east_north) ## 5919 sites

#### save site list info ####
write.csv(bbs_site_data, file = "../Data/BTO_data/site_list_east_north_bbs.csv", row.names=FALSE)

## read in the site_years data which has sites with at least 10 years of data
site_years <- read.csv("../Data/Bird_sync_data/sites_10_yrs_data.csv", header=TRUE)

## merge woodland_bbs with site_years (i.e. take out sites with less than 10 years of data)
bbs_site_data <- merge(bbs_site_data, site_years, by="site_code", all.y=TRUE)
## left with 3588 sites
summary(bbs_site_data)
## 1932 NAs (these are the Irish sites which need taken out)
bbs_site_data <- na.omit(bbs_site_data)
summary(bbs_site_data)
## 3456 sites

## remove min_year, max_yaer and year_diff columns
bbs_site_data <- subset(bbs_site_data, select = -c(5:7))

######################################################
#### create list of all pairwise site comparisons ####
######################################################

site.list <- unique(as.integer(bbs_site_data$site_code))

pair.list <- t(combn(site.list,2))
pair_attr <- data.frame(pair.list)
nrow(pair_attr)
## ALL SITES = 17514321 ROWS 
## exactly the same result (but much quicker) as the for loop version
names(pair_attr) <- c("site_a", "site_b")

## add xy data for site_a and site_b
pair_attr <- merge(pair_attr, bbs_site_data, by.x="site_a", by.y="site_code")
names(pair_attr)[4:5] <- c("site_a_EAST", "site_a_NORTH") # rename them
pair_attr <- merge(pair_attr, bbs_site_data, by.x="site_b", by.y="site_code")
names(pair_attr)[7:8] <- c("site_b_EAST", "site_b_NORTH") # rename them
pair_attr <- pair_attr[-c(3,6)] # remove gridref columns

# sort out column order #
pair_attr <- pair_attr[,c("site_a", "site_b", "site_a_EAST", "site_a_NORTH", "site_b_EAST", "site_b_NORTH")]
summary(pair_attr)

# take average northing between the pairwise comparisons #
pair_attr$mean_northing <-rowMeans(subset(pair_attr, select = c(site_a_NORTH, site_b_NORTH)), na.rm = TRUE)

### Estimate distance between sites ###
pair_attr$distance <- sqrt((pair_attr$site_a_NORTH-pair_attr$site_b_NORTH)^2 + (pair_attr$site_a_EAST-pair_attr$site_b_EAST)^2) ## pythagoras equation
pair_attr$distance <- pair_attr$distance/1000
## this gives distance in km

summary(pair_attr)

write.csv(pair_attr, file="../Data/BTO_data/pair_attr_mean_north_dist_BBS.csv", row.names=FALSE)

###### Calculate habitat similarity index

## read in hab_data
hab_data <- read.table("../Data/Land_cover_data/BBS.CBC.UKBMS.complete.landcover.data.all.scales.soil.dem.configuration.txt", header=TRUE)

## subset BBS data and 500m buffer only
hab_data <- hab_data[hab_data$Surv=="BBS",]
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
hab_data <- merge(hab_data, bbs_site_data, by.x="siteno.gref", by.y="GRIDREF", all=FALSE) ## 3304 sites with habitat data
## remove duplicated easting and northing columns
hab_data <- hab_data[-c(2:3)]
## re-order columns
hab_data <- hab_data[,c(1,15,16,17,18,2:14)]

pair_attr$renk_hab_sim <- -9999
row.names(hab_data) = hab_data$site_code

### calculate renk hab sim
all_mins = pmin(hab_data[as.character(pair_attr$site_a),6:18], hab_data[as.character(pair_attr$site_b),6:18])
pair_attr[,"renk_hab_sim"] = rowSums(all_mins)

summary(pair_attr)
### 513684 NA's => result of 152 sites without habitat data

## remove NA's
pair_attr <- na.omit(pair_attr) ## 5456556 rows :) 

## save file
## this is site data for sites with at least 10 years of data and those with habitat data to calc renk_hab_sim
write.csv(pair_attr, file="../Data/BTO_data/pair_attr_mean_north_dist_hab_sim_BBS.csv", row.names=FALSE)


