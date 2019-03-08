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

