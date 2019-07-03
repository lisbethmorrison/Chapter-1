###########################################################
## Title: Map of UKBMS sites use in FCI calculation
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: November 2017
##########################################################

#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script") ####

rm(list=ls()) # clear R


library(ggmap)
library(maps)
library(ggplot2)
library(gganimate)
library(tidyverse)
library(sf)
library(gifski)
library(mapproj)

## add data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE)
site_data <- read.csv("../Data/UKBMS_data/pair_attr_mean_northing_dist_sim.csv", header = TRUE)
  
#### gif for pint talk
## site number and lat long only
site1 <- unique(subset(site_data[c(1,3,4)]))
site2 <- unique(subset(site_data[c(2,5,6)]))
colnames(site1) <- c("site", "easting", "northing")
colnames(site2) <- c("site", "easting", "northing")
site_list <- rbind(site1, site2)
site_list <- unique(site_list)

lat_long <- site_list %>%
  st_as_sf(coords = c("easting", "northing"), crs = 27700) %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as_tibble()
site_list <- cbind(site_list, lat_long)
site_list <- transform( site_list, site = sample(site))

UK <- map_data("world") %>% filter(region=="UK")

p <- 
ggplot() +
  geom_polygon(data = UK, aes(x=long, y = lat, group = group, cumulative=TRUE), fill="grey44", alpha=0.3) +
  geom_point(data=site_list, aes(x=X, y=Y, cumulative=TRUE), size = 3, shape = 21, fill="cadetblue4", alpha = 0.7) +
  theme_void() + coord_map() +
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0, size = 16, vjust=0)) +
  transition_manual(site, cumulative = TRUE)
anim1 <- animate(p, width = 400, height = 500)
anim_save("../pintofscience.gif", animation = anim1)

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












######################### MAPS USING GGMAP ###############################
## load data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE)
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC.csv", header=TRUE) 
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE) 

## create initial UK map
UK_map <- get_map(location = c(-2.65, 53.7), zoom = 5, maptype = "satellite")
UK_map <- ggmap(ggmap=UK_map, extent = "device", legend = "right")
UK_map

## number of UKBMS sites
site1 <- unique(subset(pair_attr[c(2,9,10)]))
site2 <- unique(subset(pair_attr[c(3,11,12)]))
colnames(site1) <- c("site", "east", "north")
colnames(site2) <- c("site", "east", "north")
site_list <- rbind(site1, site2)
site_list <- unique(site_list) ## 709 sites
## change easting/northing to coordinates
UKBMS_map_coord <- site_list %>%
  st_as_sf(coords = c("east", "north"), crs = 27700) %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as_tibble()
## UKBMS sites as add points
UKBMS_map <- UK_map +
  geom_point(data = UKBMS_map_coord, aes(x=X, y=Y), size  = 1.5)
UKBMS_map

## number of CBC sites
site3 <- unique(subset(pair_attr_CBC[c(2,9,10)]))
site4 <- unique(subset(pair_attr_CBC[c(3,11,12)]))
colnames(site3) <- c("site", "east", "north")
colnames(site4) <- c("site", "east", "north")
site_list2 <- rbind(site3, site4)
site_list2 <- unique(site_list2) ## 106 sites
## change easting/northing to coordinates
CBC_map_coord <- site_list2 %>%
  st_as_sf(coords = c("east", "north"), crs = 27700) %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as_tibble()
## CBC sites as add points
CBC_map <- UK_map +
  geom_point(data = CBC_map_coord, aes(x=X, y=Y), size  = 1.5)
CBC_map


## number of BBS sites
site5 <- unique(subset(pair_attr_BBS[c(2,9,10)]))
site6 <- unique(subset(pair_attr_BBS[c(3,11,12)]))
colnames(site5) <- c("site", "east", "north")
colnames(site6) <- c("site", "east", "north")
site_list3 <- rbind(site5, site6)
site_list3 <- unique(site_list3) ## 2499 sites
## change easting/northing to coordinates
BBS_map_coord <- site_list3 %>%
  st_as_sf(coords = c("east", "north"), crs = 27700) %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as_tibble()
## CBC sites as add points
BBS_map <- UK_map +
  geom_point(data = BBS_map_coord, aes(x=X, y=Y), size  = 1)
BBS_map

## multiplot function to plot maps on same file
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


## save all 3 maps as one
png("../Graphs/all_maps.png", height = 200, width = 250, units = "mm", res = 300)
multiplot(UK_map + geom_point(data = UKBMS_map_coord, 
                              aes(x=X, y=Y), size = 0.3), 
          UK_map + geom_point(data = CBC_map_coord, 
                              aes(x=X, y=Y), size = 0.3), 
          UK_map + geom_point(data = BBS_map_coord,
                              aes(x=X, y=Y), size = 0.3), cols = 3)
dev.off()






