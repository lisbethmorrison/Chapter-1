###########################################################
## Title: Calculate BBS site data info
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: January 2018
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

## read in the site_years data which has sites at least 5 years old
site_years <- read.csv("../Data/Bird_sync_data/sites_10_yrs_data.csv", header=TRUE)

## merge woodland_bbs with site_years (i.e. take out sites younger than 5 years)
bbs_site_data <- merge(bbs_site_data, site_years, by="site_code", all.y=TRUE)
## left with 3588 sites
summary(bbs_site_data)
## 1932 NAs (these are the Irish sites which need taken out)
bbs_site_data <- na.omit(bbs_site_data)
summary(bbs_site_data)
## 3456 sites 

## remove min_year, max_yaer and year_diff columns
bbs_site_data <- subset(bbs_site_data, select = -c(5:7)) 

## randomly select 25% of the sites (repeat this 4 times, so end up with 4 pair_attr files)
library(dplyr)
bbs_site_data_4 <- bbs_site_data %>% sample_frac(.25)

######################################################
#### create list of all pairwise site comparisons ####
######################################################
site.list <- unique(as.character(bbs_site_data_4$site_code))
num.ts <- length(site.list)

pair.list2 <- NULL
for(i in 1:(num.ts-1)){
  print(i)
  pair.list2 <- rbind(pair.list2, cbind(rep(site.list[i], num.ts-i), (site.list[(i+1):num.ts])))
}
nrow(pair.list2)   ## 839,160 for each 25% random sample of data (plus excluding sites younger than 5 years old)
## 17,514,321 rows (if all sites included) or 14,453,376 if only sites older than 5 years included

pair_attr2 <- data.frame(pair.list2)
names(pair_attr2) <- c("site_a", "site_b")

## add xy data for site_a and site_b
pair_attr2 <- merge(pair_attr2, bbs_site_data, by.x="site_a", by.y="site_code")
names(pair_attr2)[4:5] <- c("site_a_EAST", "site_a_NORTH") # rename them
pair_attr2 <- merge(pair_attr2, bbs_site_data, by.x="site_b", by.y="site_code")
names(pair_attr2)[7:8] <- c("site_b_EAST", "site_b_NORTH") # rename them
pair_attr2 <- pair_attr2[-c(3,6)] # remove gridref columns

# sort out column order #
pair_attr2 <- pair_attr2[,c("site_a", "site_b", "site_a_EAST", "site_a_NORTH", "site_b_EAST", "site_b_NORTH")]

# take average northing between the pairwise comparisons #
pair_attr2$mean_northing <- -9999

for (i in 1:nrow(pair_attr2)){
  print(i)
  pair_attr2$mean_northing[i] <- mean(c(pair_attr2[i,"site_a_NORTH"], pair_attr2[i,"site_b_NORTH"]))
}

### Estimate distance between sites ###
# use distance matrix for this?!
temp <- bbs_site_data[,c("EASTING", "NORTHING")]
temp <- as.matrix(temp) # check this

#Run spDists on a matrix of site_x, site_y.  It estimates the euclidean distance between all pairwise combinations.
install.packages("sp")
library(sp)
dist_mat <- spDists(temp, temp, longlat = FALSE)
rownames(dist_mat) <- bbs_site_data_4$site_code
colnames(dist_mat) <- bbs_site_data_4$site_code
dist_mat <- as.data.frame(dist_mat)

# add in distance column to pair attr table
pair_attr_4$distance <- -9999

for (i in 1:nrow(pair_attr_4)){
  print(i)
  pair_attr_4$distance[i] <- dist_mat[rownames(dist_mat)==as.character(pair_attr_4[i,"site_a"]),colnames(dist_mat)==as.character(pair_attr_4[i, "site_b"])]
}

### save pair_attr ###
write.csv(pair_attr_4, file="../Data/BTO_data/pair_attr_mean_north_dist_BBS_4.csv", row.names=FALSE)

