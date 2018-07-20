###########################################################
## Title: Moving window calcualte population synchrony BBS BTO data
## User: Lisbeth Morrison
## email: l.morrison@pgr.reading.ac.uk
## Date: January 2018
##########################################################

#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script") ####

############## SCRIPT FOR THE FIRST 25% OF RANDOM SITE DATA ####################

rm(list=ls()) # clear R

library(dplyr)

### add data
final_data <- read.csv("../Data/Bird_sync_data/final_data_all_spp_BBS.csv", header=TRUE) ## 31 spp
woodland_bbs <- read.csv("../Data/Bird_sync_data/bbs_woodland_birds.csv", header=TRUE)
site_data <- read.csv("../Data/BTO_data/pair_attr_mean_north_dist_BBS_10yrs.csv", header=TRUE) 

### Create "good species list" which meet a minimum criteria of 2 filters
woodland_bbs<-woodland_bbs[woodland_bbs$YEAR>=2000&woodland_bbs$YEAR<=2016,]
woodland_bbs$species_code <- as.factor(woodland_bbs$species_code)
bbs_summary_tab<-with(woodland_bbs,table(species_code,YEAR)) # 31 species
bbs_summary_good_years<-NULL

for (i in 1:nrow(bbs_summary_tab)){
  no.good.years<-length(bbs_summary_tab[i,][bbs_summary_tab[i,]>=50])  # filter 1: 50 sites/year
  sp<-row.names(bbs_summary_tab)[i]
  results.temp<-data.frame(sp,no.good.years)
  bbs_summary_good_years<-rbind(bbs_summary_good_years,results.temp)
}
bbs_summary_good_years
# filter 2: only select species with >75% 'good years' (i.e. with more than 50 sites/year)
length(2000:2016)*0.75 # 12.75
good.species.list<-bbs_summary_good_years$sp[bbs_summary_good_years$no.good.years>12.75]
good.species.list
good.species.list <- droplevels(good.species.list)

##### removes 2 species ==> willow tit and wood warbler (these get removed during synchrony analysis anyway)

#############################
#### calculate synchrony ####
#############################

spp.list <- unique(final_data$name)

final_pair_data <- NULL
final_summ_stats <- NULL
final_summary <- NULL
### split based on species ###
for(g in spp.list[1]){

# loop through spp.list 
  spp_data <- final_data[final_data$name==g,]
  total_comp <- NULL
  print(paste("species",g))

  year.list<-1994:2007

  for (i in year.list[1]){ # loop through years
    start.year<-i
    mid.year<-i+4.5
    print(paste("mid.year=",mid.year))
    end.year<-i+9
    species.10.yr.data<-spp_data[spp_data$year>=start.year&spp_data$year<=end.year,] ## create dataframe with data for 10 year moving window

    ################################
    # create a matrix to be filled #
    site.list <- unique(species.10.yr.data$site)
    year.list.temp <- min(species.10.yr.data$year):max(species.10.yr.data$year)
    Grow.matrix<-matrix(c(species.10.yr.data$gr), nrow=length(year.list.temp)) ## create growth matrix with growth rate data
    rownames(Grow.matrix)<-year.list.temp ## moving window years as rows
    colnames(Grow.matrix)<-site.list ## site names as column names

    num.ts <- length(site.list)   # number of time series
    TS <- matrix(Grow.matrix, ncol=num.ts)    # stores Grow data in a matrix - just removes row names
    # create site match table used to correct sites names after calculating synchrony
    site_match <- data.frame(TS_name = 1:ncol(TS), site_name = colnames(Grow.matrix))

    ###############################
    # cross-correlation functions #

    # calculates combination between all the pairs of sites
    pair.list <- t(combn(1:num.ts,2))
    nrow(pair.list) 

    colnames(pair.list) <- c("site1","site2")

    #### Add in the real site names here, merge with site_data, then cut the table to only include pairs where distance<100km
    
    ###  correct site numbers ###
    pair.list <- merge(pair.list, site_match, by.x="site1", by.y="TS_name")
    pair.list <- pair.list[,c("site_name", "site2")] ## get rid of site1 (which are not real site names)
    names(pair.list) <- gsub("site_name", "site1", names(pair.list))
    pair.list <- merge(pair.list, site_match, by.x="site2", by.y="TS_name")
    pair.list <- pair.list[,c("site1", "site_name")] ## get rid of site2 (which are not real site names)
    names(pair.list) <- gsub("site_name", "site2", names(pair.list)) ## now site1 and site2 are REAL site names

    ## merge pair.list with site data
    site_data[,c(1,2)] <- sapply(site_data[,c(1,2)], as.factor)
    site_data_pair1 <- merge(pair.list, site_data, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))
    site_data_reverse <- site_data
    names(site_data_reverse)[1:6] <- c("site_b", "site_a", "site_b_EAST", "site_b_NORTH", "site_a_EAST", "site_a_NORTH")
    site_data_pair2 <- merge(pair.list, site_data_reverse, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))
    site_data_pair <- rbind(site_data_pair1, site_data_pair2)
    
    # length(unique(pair.list$site1)[!unique(pair.list$site1)%in%unique(site_data$site_a)])   # how many sites are missing attribute data

    ## remove site pair comparisons more than 100km apart
    site_data_pair <- site_data_pair[site_data_pair$distance <= 100,] 

    ## merge pair.list with site_data pair so pair.list only contains sites <100km apart
    pair.list<-merge(pair.list,site_data_pair,by.x=c("site1", "site2"),by.y=c("site1", "site2"))[,1:2] ## pair list becomes same length as site_data_pair

    # change site names back to numbers
    pair.list <- merge(pair.list, site_match, by.x="site1", by.y="site_name")
    pair.list <- pair.list[,c("TS_name", "site2")]
    names(pair.list) <- gsub("TS_name", "site1", names(pair.list))
    pair.list <- merge(pair.list, site_match, by.x="site2", by.y="site_name")
    pair.list <- pair.list[,c("site1", "TS_name")]
    names(pair.list) <- gsub("TS_name", "site2", names(pair.list))
    
    # ### take random 25% of sites
    # pair.list <- pair.list[sample(0.25*(nrow(pair.list))),]

    # calculate cross-correlation functions at different time lags
    # (NB - we are only be interested in the current year, i.e. lag = 0)
    max.lag <- 0 # maximum time lag over which cross-correlations are assessed
    CCF <- data.frame(matrix(NA, ncol = 2*max.lag+2, nrow = nrow(pair.list)))	## matrix to store cross-correlation coefficients
    ## bind CCF with pair.list
    CCF <- cbind(CCF, pair.list)
    names(CCF) <- c(paste("lag", -max.lag:max.lag, sep=""), "numYears", "site1", "site2") ## adds column names

    for(k in 1:nrow(pair.list)) { ## loop through each row in pair.list
      ## calculate number of years to base correlation on (how much data is in common between pairs of sites)
      CCF$numYears[k] <- length(na.omit((TS[,pair.list[k,1]]+TS[,pair.list[k,2]]))) ## add up number of times when sum can be done (i.e. neither is NA) 
    }
    
    summary <- data.frame(spp=g, year=i+4.5, no.sites=nrow(CCF), sampled.sites=10000)
    summary$skipped=ifelse(summary$no.sites<2, "yes", "no")
    summary$sampled=ifelse(summary$no.sites>=10000, "yes", "no")
    if(summary$sampled=="yes"){
      summary$percentage=(summary$sampled.sites/summary$no.sites)*100
    }else{
      summary$percentage=NA
    }
    final_summary <- rbind(final_summary,summary)
  }
}

write.csv(final_summary, file="../Data/Bird_sync_data/BBS_summary_8_31.csv", row.names=FALSE) 

    CCF1 <- CCF
    ## rank by numYears
    CCF1 <- arrange(CCF1, desc(numYears))
    ## chop out rows where numYears<7
    CCF1 <- CCF1[CCF1$numYears>6,]
    
    # if statement to skip species which won't run
    if (nrow(CCF1)<2){
      print(paste("skip species", g, "year", i+4.5))
      next
    }
    ## nrow of CCF
    nrow(CCF1) 
    CCF1$site1 <- as.numeric(CCF1$site1)
    CCF1$site2 <- as.numeric(CCF1$site2)

    if (nrow(CCF1)>=10000){
      CCF1 <- CCF1[sample(nrow(CCF1), 10000),]
      print(paste("species", g, "year", i+4.5, "sampled"))
    }else{
      print(paste("species", g, "year", i+4.5, "not sampled"))
    }

    pair.list <- CCF1[,c(3,4)] ## subset pair.list to have only sites that have been sampled
    
    for(k in 1:nrow(pair.list)) {
      ## calculation of CCF
      try(CCF1$lag0[k] <- ccf(TS[,pair.list[k,1]], TS[,pair.list[k,2]], lag.max=max.lag, na.action=na.exclude, plot=F, type="correlation")$acf, silent=T) ## use try to prevent crashes if no data
    } ## this will take a long time!!! ## end k in nrow(pair.list)
  
    pair.attr <- CCF1 ## matrix to hold attributes of pairs
    # colnames(pair.attr) <- c("site1","site2")
    # pair.attr <- cbind(pair.attr, CCF1)     ### add in correlation scores and number of comparisons at each site
    # head(pair.attr)

    ###  correct site numbers ### 
    pair.attr <- merge(pair.attr, site_match, by.x="site1", by.y="TS_name")  
    pair.attr <- pair.attr[,c("site_name", "site2", "lag0", "numYears")]
    names(pair.attr) <- gsub("site_name", "site1", names(pair.attr))
    pair.attr <- merge(pair.attr, site_match, by.x="site2", by.y="TS_name")  
    pair.attr <- pair.attr[,c("site1", "site_name", "lag0", "numYears")]
    names(pair.attr) <- gsub("site_name", "site2", names(pair.attr))

    # new_pair_attr <- pair.attr[pair.attr$numYears>6,] #### THIS IS THE LINE OF CODE WHICH STOPS SPECIES WORKING ####
      pair.attr <- pair.attr[!is.na(pair.attr$lag0),]
      pair.attr$mid.year<-mid.year
      pair.attr$start.year<-start.year
      pair.attr$end.year<-end.year

      # create a table of all pair-wise comparisons, this will be built up as we move through years
      all_pair_attr <- NULL
      all_pair_attr <- rbind(all_pair_attr, pair.attr)

      ### record the number of good pair-wise comparisons and number of unique sites
      temp_summ_stats <- data.frame(total_comps = nrow(pair.attr), uni_sites = length(unique(c(pair.attr$site1, pair.attr$site2))))
      tp_summ_stats <- NULL
      tp_summ_stats <- rbind(tp_summ_stats, temp_summ_stats)
      tp_summ_stats$mid.year<-mid.year
      tp_summ_stats$start.year<-start.year
      tp_summ_stats$end.year<-end.year

      ### add in species info
      tp_summ_stats$spp <- g
      final_summ_stats <- rbind(final_summ_stats, tp_summ_stats)

      all_pair_attr$spp <- g
      final_pair_data <- rbind(final_pair_data, all_pair_attr)

  } # end i in year
} # end g in species

### merge in predictor variables
head(final_pair_data)
head(final_summ_stats)

summary(final_pair_data)

write.csv(final_pair_data, file="../Data/Bird_sync_data/final_pair_data_all_spp_BBS_4th_run.csv", row.names=FALSE) ## save final pair data for all 33 species

write.csv(final_summ_stats, file="../Data/Bird_sync_data/final_summ_stats_all_spp_BBS_4th_run.csv", row.names=FALSE) ## save final summ stats for all 33 species

### species which were skipped for ALL years: 
## 467, 298, 450, 123 ==> 27 species were NOT skipped

## species which were skipped for SOME years:
## 456, 431, 468

############################ left with 24 SPECIES with data for all years ############################ 

### species which were sampled: (for some OR all years)
## 320, 322, 464, 463, 403, 368, 365, 451, 503, 428, 453, 361, 398

