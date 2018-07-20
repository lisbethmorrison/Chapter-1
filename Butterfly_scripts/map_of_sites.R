###########################################################
## Title: Map of UKBMS sites use in FCI calculation
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: November 2017
##########################################################

#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script") ####

rm(list=ls()) # clear R

## add data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE)
site_data <- read.csv("../Data/UKBMS_data/pair_attr_mean_northing_dist_sim.csv", header = TRUE)
  
### save data on included sites for plot
# create the following: site, spp for which it is included for. Do this for each TP.
## TP1 = 1980 DATA i.e. early
## TP2 = 1995 DATA i.e. mid
## TP3 = 2007 DATA i.e. late
inc_sites_T1 <- NULL
inc_sites_T2 <- NULL
inc_sites_T3 <- NULL
SPECIES_T1 <- NULL
SPECIES_T2 <- NULL
SPECIES_T3 <- NULL

for (i in unique(pair_attr$spp)){
  temp_attr <- pair_attr[pair_attr$spp==i,]
  inc_sites_T1 <- c(inc_sites_T1, unique(c(temp_attr[temp_attr$start.year=="1980", "site1"], temp_attr[temp_attr$start.year=="1980", "site2"])))
  inc_sites_T2 <- c(inc_sites_T2, unique(c(temp_attr[temp_attr$start.year=="1995", "site1"], temp_attr[temp_attr$start.year=="1995", "site2"])))
  inc_sites_T3 <- c(inc_sites_T3, unique(c(temp_attr[temp_attr$start.year=="2007", "site1"], temp_attr[temp_attr$start.year=="2007", "site2"])))
  SPECIES_T1 <- c(SPECIES_T1, rep(i, length(unique(c(temp_attr[temp_attr$start.year=="1980", "site1"], temp_attr[temp_attr$start.year=="1980", "site2"])))))
  SPECIES_T2 <- c(SPECIES_T2, rep(i, length(unique(c(temp_attr[temp_attr$start.year=="1995", "site1"], temp_attr[temp_attr$start.year=="1995", "site2"])))))
  SPECIES_T3 <- c(SPECIES_T3, rep(i, length(unique(c(temp_attr[temp_attr$start.year=="2007", "site1"], temp_attr[temp_attr$start.year=="2007", "site2"])))))
}

T1_map_data <- data.frame(site = inc_sites_T1, spp = SPECIES_T1)
T2_map_data <- data.frame(site = inc_sites_T2, spp = SPECIES_T2)
T3_map_data <- data.frame(site = inc_sites_T3, spp = SPECIES_T3)

T1_map <- data.frame(site = names(table(T1_map_data$site)), spp_count = as.numeric(table(T1_map_data$site)))
T2_map <- data.frame(site = names(table(T2_map_data$site)), spp_count = as.numeric(table(T2_map_data$site)))
T3_map <- data.frame(site = names(table(T3_map_data$site)), spp_count = as.numeric(table(T3_map_data$site)))

# merge in the site xy data #
temp_site_data <- site_data[,c("site_a", "site_a_NORTH", "site_a_EAST")]
names(temp_site_data) <- c("site", "north", "east")
temp_site_data2 <- site_data[,c("site_b", "site_b_NORTH", "site_b_EAST")]
names(temp_site_data2) <- c("site", "north", "east")

T1_map <- merge(T1_map, unique(rbind(temp_site_data, temp_site_data2)))
T2_map <- merge(T2_map, unique(rbind(temp_site_data, temp_site_data2)))
T3_map <- merge(T3_map, unique(rbind(temp_site_data, temp_site_data2)))

install.packages("blighty")
library(blighty)

par(mfrow=c(1,3)) 

x1 <- (T1_map$east/1000)
y1 <- (T1_map$north/1000)
x2 <- (T2_map$east/1000)
y2 <- (T2_map$north/1000)
x3 <- (T3_map$east/1000)
y3 <- (T3_map$north/1000)

## T1 map 
T1.map <- blighty(place="set.British.Isles") +
  points(x1,y1, col="black", pch=19)

## T2 map
T2.map <- blighty(place="set.British.Isles") +
  points(x2,y2, col="black", pch=19)

## T3 map
T3.map <- blighty(place="set.British.Isles") +
  points(x3,y3, col="black", pch=19)

png("../Graphs/Maps_of_sites_UKBMS.png", width=1000, height=500)
par(mfrow=c(1,3))
T1.map <- blighty(place="set.British.Isles") +
  points(x1,y1, col="black", pch=19) +
  title(main = "1985", cex.main=2)
T2.map <- blighty(place="set.British.Isles") +
  points(x2,y2, col="black", pch=19) +
  title(main = "2000", cex.main=2)
T3.map <- blighty(place="set.British.Isles") +
  points(x3,y3, col="black", pch=19) +
  title(main = "2012", cex.main=2)
dev.off()

##### map of woodland sites only #####

pair_attr <- subset(pair_attr[pair_attr$HABITAT=="Woodland",])

### save data on included sites for plot
# create the following: site, spp for which it is included for. Do this for each TP.
## TP1 = 1980 DATA i.e. early
## TP2 = 1995 DATA i.e. mid
## TP3 = 2007 DATA i.e. late
inc_sites_T1 <- NULL
inc_sites_T2 <- NULL
inc_sites_T3 <- NULL
SPECIES_T1 <- NULL
SPECIES_T2 <- NULL
SPECIES_T3 <- NULL

for (i in unique(pair_attr$spp)){
  temp_attr <- pair_attr[pair_attr$spp==i,]
  inc_sites_T1 <- c(inc_sites_T1, unique(c(temp_attr[temp_attr$start.year=="1980", "site1"], temp_attr[temp_attr$start.year=="1980", "site2"])))
  inc_sites_T2 <- c(inc_sites_T2, unique(c(temp_attr[temp_attr$start.year=="1995", "site1"], temp_attr[temp_attr$start.year=="1995", "site2"])))
  inc_sites_T3 <- c(inc_sites_T3, unique(c(temp_attr[temp_attr$start.year=="2007", "site1"], temp_attr[temp_attr$start.year=="2007", "site2"])))
  SPECIES_T1 <- c(SPECIES_T1, rep(i, length(unique(c(temp_attr[temp_attr$start.year=="1980", "site1"], temp_attr[temp_attr$start.year=="1980", "site2"])))))
  SPECIES_T2 <- c(SPECIES_T2, rep(i, length(unique(c(temp_attr[temp_attr$start.year=="1995", "site1"], temp_attr[temp_attr$start.year=="1995", "site2"])))))
  SPECIES_T3 <- c(SPECIES_T3, rep(i, length(unique(c(temp_attr[temp_attr$start.year=="2007", "site1"], temp_attr[temp_attr$start.year=="2007", "site2"])))))
}

T1_map_data <- data.frame(site = inc_sites_T1, spp = SPECIES_T1)
T2_map_data <- data.frame(site = inc_sites_T2, spp = SPECIES_T2)
T3_map_data <- data.frame(site = inc_sites_T3, spp = SPECIES_T3)

T1_map <- data.frame(site = names(table(T1_map_data$site)), spp_count = as.numeric(table(T1_map_data$site)))
T2_map <- data.frame(site = names(table(T2_map_data$site)), spp_count = as.numeric(table(T2_map_data$site)))
T3_map <- data.frame(site = names(table(T3_map_data$site)), spp_count = as.numeric(table(T3_map_data$site)))

# merge in the site xy data #
temp_site_data <- site_data[,c("site_a", "site_a_NORTH", "site_a_EAST")]
names(temp_site_data) <- c("site", "north", "east")
temp_site_data2 <- site_data[,c("site_b", "site_b_NORTH", "site_b_EAST")]
names(temp_site_data2) <- c("site", "north", "east")

T1_map <- merge(T1_map, unique(rbind(temp_site_data, temp_site_data2)))
T2_map <- merge(T2_map, unique(rbind(temp_site_data, temp_site_data2)))
T3_map <- merge(T3_map, unique(rbind(temp_site_data, temp_site_data2)))

library(blighty)

par(mfrow=c(1,3)) 

x1 <- (T1_map$east/1000)
y1 <- (T1_map$north/1000)
x2 <- (T2_map$east/1000)
y2 <- (T2_map$north/1000)
x3 <- (T3_map$east/1000)
y3 <- (T3_map$north/1000)

## T1 map 
T1.map <- blighty(place="set.British.Isles") +
  points(x1,y1, col="black", pch=19)

## T2 map
T2.map <- blighty(place="set.British.Isles") +
  points(x2,y2, col="black", pch=19)

## T3 map
T3.map <- blighty(place="set.British.Isles") +
  points(x3,y3, col="black", pch=19)

png("../Graphs/Maps_of_woodland_sites_UKBMS.png", width=1000, height=500)
par(mfrow=c(1,3))
T1.map <- blighty(place="set.British.Isles") +
  points(x1,y1, col="black", pch=19) +
  title(main = "(a) 1985", cex.main=2)
T2.map <- blighty(place="set.British.Isles") +
  points(x2,y2, col="black", pch=19) +
  title(main = "(b) 2000", cex.main=2)
T3.map <- blighty(place="set.British.Isles") +
  points(x3,y3, col="black", pch=19) +
  title(main = "(c) 2012", cex.main=2)
dev.off()


#### number of woodland sites
pair_attr_wood <- subset(pair_attr[pair_attr$HABITAT=="Woodland",])

site1_wood <- as.data.frame(unique(pair_attr_wood$site1))
site2_wood <- as.data.frame(unique(pair_attr_wood$site2))
colnames(site1_wood)[1] <- "site"
colnames(site2_wood)[1] <- "site"
site_list <- rbind(site1_wood, site2_wood)
site_list <- unique(site_list)

## number of all sites
site1 <- as.data.frame(unique(pair_attr$site1))
site2 <- as.data.frame(unique(pair_attr$site2))
colnames(site1)[1] <- "site"
colnames(site2)[1] <- "site"
site_list <- rbind(site1, site2)
site_list <- unique(site_list) ## 709 sites

pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC.csv", header=TRUE) 
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS_final.csv", header=TRUE) 

## number of CBC sites
site3 <- as.data.frame(unique(pair_attr_CBC$site1))
site4 <- as.data.frame(unique(pair_attr_CBC$site2))
colnames(site3)[1] <- "site"
colnames(site4)[1] <- "site"
site_list2 <- rbind(site3, site4)
site_list2 <- unique(site_list2) ## 106 sites

## number of BBS sites
site5 <- as.data.frame(unique(pair_attr_BBS$site1))
site6 <- as.data.frame(unique(pair_attr_BBS$site2))
colnames(site5)[1] <- "site"
colnames(site6)[1] <- "site"
site_list3 <- rbind(site5, site6)
site_list3 <- unique(site_list3) ## 2499 sites



library(ggmap)
library(ggplot2)
UK <- map_data(map = "world", region = "UK") # changed map to "world"
ggplot(data = UK, aes(x = long, y = lat, group = group)) + 
  geom_polygon() +
  coord_map()

d <- data.frame(lat=c(50.659631, 50.607213, 50.608129),
                lon=c(3.09319, 3.011473, 3.031529))
UK + geom_point(data=d, aes(x=lon, y=lat), color="red", size=30, alpha=0.5)

UK_map <- get_map(location = c(-2.65, 53.7), zoom = 5, maptype = "hybrid")
UK_map <- ggmap(ggmap=UK_map, extent = "device", legend = "right")

UK_map_1 <- UK_map +
        ggmap(ggmap=UK_map, extent = "device", legend = "right") +
        geom_point(data = T1_map_coord, aes(x=X, y=Y), size  = 2)

UK_map_2 <- UK_map +
  ggmap(ggmap=UK_map, extent = "device", legend = "right") +
  geom_point(data = T2_map_coord, aes(x=X, y=Y), size  = 2)

UK_map_3 <- UK_map +
  ggmap(ggmap=UK_map, extent = "device", legend = "right") +
  geom_point(data = T3_map_coord, aes(x=X, y=Y), size  = 2)



multiplot(UK_map + geom_point(data = T1_map_coord, 
                              aes(x=X, y=Y), size = 0.8), 
          UK_map + geom_point(data = T2_map_coord, 
                              aes(x=X, y=Y), size = 0.8), 
          UK_map + geom_point(data = T3_map_coord,
                              aes(x=X, y=Y), size = 0.8), cols = 3)







multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}















library(tidyverse)
library(sf)
#> Linking to GEOS 3.6.1, GDAL 2.2.3, proj.4 4.9.3
T1_map_coord <- T1_map %>%
  st_as_sf(coords = c("east", "north"), crs = 27700) %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as_tibble()

T2_map_coord <- T2_map %>%
  st_as_sf(coords = c("east", "north"), crs = 27700) %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as_tibble()

T3_map_coord <- T3_map %>%
  st_as_sf(coords = c("east", "north"), crs = 27700) %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as_tibble()
