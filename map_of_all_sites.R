###########################################################
## Title: MMap of all sites for UKBMS, CBC and BBS
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: November 2018
##########################################################

##### NOTE: Blighty can only be used on old versions of R (3.0.0)

library(blighty)
library(gridExtra)

## add data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE)
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC_no_zeros2.csv", header=TRUE)
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE)

###########################
## number of UKBMS sites ##
###########################

site1 <- unique(subset(pair_attr[c(2,9,10)]))
site2 <- unique(subset(pair_attr[c(3,11,12)]))
colnames(site1) <- c("site", "east", "north")
colnames(site2) <- c("site", "east", "north")
site_list <- rbind(site1, site2)
site_list <- unique(site_list) ## 709 sites

### divide by 10000 to work with Blighty
site_list$east <- site_list$east/1000
site_list$north <- site_list$north/1000

### ukbms map of all sites (709)
UKBMS_map <- blighty(place="set.British.Isles") +
  points(site_list$east,site_list$north, col="black", pch=19)

#########################
## number of CBC sites ##
#########################

site1 <- unique(subset(pair_attr_CBC[c(2,10,11)]))
site2 <- unique(subset(pair_attr_CBC[c(3,12,13)]))
colnames(site1) <- c("site", "east", "north")
colnames(site2) <- c("site", "east", "north")
site_list_CBC <- rbind(site1, site2)
site_list_CBC <- unique(site_list_CBC) ## 109 sites

### divide by 10000 to work with Blighty
site_list_CBC$east <- site_list_CBC$east/1000
site_list_CBC$north <- site_list_CBC$north/1000

### ukbms map of all sites (109)
CBC_map <- blighty(place="set.British.Isles") +
  points(site_list_CBC$east,site_list_CBC$north, col="black", pch=19)

#########################
## number of BBS sites ##
#########################

site1 <- unique(subset(pair_attr_BBS[c(2,9,10)]))
site2 <- unique(subset(pair_attr_BBS[c(3,11,12)]))
colnames(site1) <- c("site", "east", "north")
colnames(site2) <- c("site", "east", "north")
site_list_BBS <- rbind(site1, site2)
site_list_BBS <- unique(site_list_BBS) ## 109 sites

### divide by 10000 to work with Blighty
site_list_BBS$east <- site_list_BBS$east/1000
site_list_BBS$north <- site_list_BBS$north/1000

### ukbms map of all sites (109)
BBS_map <- blighty(place="set.British.Isles") +
  points(site_list_BBS$east,site_list_BBS$north, col="black", pch=19)

##### plot all 3 maps together
png("../Graphs/FINAL/FigureS1.png", width=1000, height=400)
par(mfrow=c(1,3))
UKBMS_map <- blighty(place="set.British.Isles") +
  points(site_list$east,site_list$north, col="black", pch=19)
CBC_map <- blighty(place="set.British.Isles") +
  points(site_list_CBC$east,site_list_CBC$north, col="black", pch=19)
BBS_map <- blighty(place="set.British.Isles") +
  points(site_list_BBS$east,site_list_BBS$north, col="black", pch=19)
dev.off()
