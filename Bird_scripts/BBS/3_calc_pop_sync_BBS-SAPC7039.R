###########################################################
## Title: Moving window calcualte population synchrony BBS BTO data  
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: January 2018
##########################################################

#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script") ####

############## SCRIPT FOR THE FIRST 25% OF RANDOM SITE DATA ####################

rm(list=ls()) # clear R

### add data
final_data <- read.csv("../Data/Bird_sync_data/final_data_all_spp_BBS.csv", header=TRUE) ## 31 spp
woodland_bbs <- read.csv("../Data/Bird_sync_data/bbs_woodland_birds.csv", header=TRUE)
site_data_1 <- read.csv("../Data/BTO_data/pair_attr_mean_north_dist_BBS_1.csv", header=TRUE)
site_data_2 <- read.csv("../Data/BTO_data/pair_attr_mean_north_dist_BBS_2.csv", header=TRUE)
site_data_3 <- read.csv("../Data/BTO_data/pair_attr_mean_north_dist_BBS_3.csv", header=TRUE)
site_data_4 <- read.csv("../Data/BTO_data/pair_attr_mean_north_dist_BBS_4.csv", header=TRUE)

### try and merge the two site_data files 
site_data <- rbind(site_data_1, site_data_2, site_data_3, site_data_4)
site_data <- unique(site_data)

# ### Create "good species list" which meet a minimum criteria of 2 filters
# woodland_bbs<-woodland_bbs[woodland_bbs$YEAR>=1980&woodland_bbs$YEAR<=2016,]
# bbs_summary_tab<-with(woodland_bbs,table(species_code,YEAR))
# bbs_summary_good_years<-NULL
# 
# for (i in 1:nrow(bbs_summary_tab)){
#   no.good.years<-length(bbs_summary_tab[i,][bbs_summary_tab[i,]>=50])  # filter 1: 50 sites/year
#   sp<-row.names(bbs_summary_tab)[i]
#   results.temp<-data.frame(sp,no.good.years)
#   bbs_summary_good_years<-rbind(bbs_summary_good_years,results.temp)
# }
# bbs_summary_good_years
# # filter 2: only select species with >75% 'good years' (i.e. with more than 50 sites/year)
# length(1994:2016)*0.75 # 17.25
# good.species.list<-bbs_summary_good_years$sp[bbs_summary_good_years$no.good.years>17.25]
# good.species.list
# 
# ############## CONCLUSION => ALL 31 SPECIES MEET THE ABOVE CRITERIA, SO ISN'T NEEDED ##################

#############################
#### calculate synchrony ####
#############################

spp.list <- unique(final_data$name)

final_pair_data <- NULL
final_summ_stats <- NULL

### split based on species ###
for (g in spp.list){ # loop through spp.list ### spp 123, 511, 456, 450, 431 are causing issues
  spp_data <- final_data[final_data$name==g,]
  total_comp <- NULL
  print(paste("species",g))    
  
  year.list<-1994:2007
  
  for (i in year.list){ #loop through years
    start.year<-i
    mid.year<-i+4.5
    print(paste("mid.year=",mid.year)) 
    end.year<-i+9
    species.10.yr.data<-spp_data[spp_data$year>=start.year&spp_data$year<=end.year,]
    
    ################################
    # create a matrix to be filled #
    site.list <- unique(species.10.yr.data$site)
    site.list
    year.list.temp <- min(species.10.yr.data$year):max(species.10.yr.data$year)
    year.list.temp
    Grow.matrix<-matrix(c(species.10.yr.data$gr), nrow=length(year.list.temp))
    Grow.matrix
    ncol(Grow.matrix)
    rownames(Grow.matrix)<-year.list.temp
    colnames(Grow.matrix)<-site.list
    Grow.matrix
    
    length(site.list)
    
    num.ts <- length(site.list)   # number of time series
    TS <- matrix(Grow.matrix, ncol=num.ts)    # stores Grow data in a matrix     (simply removes names!)
    TS
    # create site match table used to correct sites names after calculating synchrony
    site_match <- data.frame(TS_name = 1:ncol(TS), site_name = colnames(Grow.matrix))
    
    ###############################
    # cross-correlation functions #
    
    pair.list <- NULL             # calculates combination between all the pairs of sites
    for(k in 1:(num.ts-1)){
      pair.list <- rbind(pair.list, cbind(rep(k, num.ts-k), (k+1):num.ts)) ## species 508 has only 2 sites, so this line of code doesn't work...
    } # end k in num.ts-1
    nrow(pair.list) 
    
    colnames(pair.list) <- c("site1","site2")
    
    #### NOTE- ADD IN THE REAL SITE NAMES HERE, MERGE WITH SITE_DATA TABLE, THEN  CUT DOWN THE TABLE TO ONLY INCLUDE PAIRS WHERE distance < 100000 (100km)) ####
    ## e.g.   merge(data1,data2,by.x=c(col1,col2),by.y=c(site_a,site_b)) ##
    
    ###  correct site numbers ### 
    pair.list <- merge(pair.list, site_match, by.x="site1", by.y="TS_name")  
    pair.list <- pair.list[,c("site_name", "site2")]
    names(pair.list) <- gsub("site_name", "site1", names(pair.list))
    pair.list <- merge(pair.list, site_match, by.x="site2", by.y="TS_name")  
    pair.list <- pair.list[,c("site1", "site_name")]
    names(pair.list) <- gsub("site_name", "site2", names(pair.list))
    
    ## merge pair.list with site data
    site_data_1[,c(1,2)] <- sapply(site_data_1[,c(1,2)], as.factor)
    site_data_pair1 <- merge(pair.list, site_data_1, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))
    site_data_pair2 <- merge(pair.list, site_data_1, by.x=c("site1", "site2"), by.y=c("site_b", "site_a"))
    site_data_pair <- rbind(site_data_pair1, site_data_pair2)

    # unique(pair.list$site1)[!unique(pair.list$site1)%in%unique(site_data$site_a)]   # how many sites are missing attribute data

    ## remove site pair comparisons more than 100km apart
    site_data_pair <- site_data_pair[site_data_pair$distance <= 100000,]

    pair.list<-merge(pair.list,site_data_pair,by.x=c("site1", "site2"),by.y=c("site1", "site2"))[,1:2]

    
    ## calculate cross-correlation functions at different time lags
    ## (NB - we are only be interested in the current year, i.e. lag = 0)
    max.lag <- 0 # maximum time lag over which cross-correlations are assessed
    CCF <- data.frame(matrix(NA, ncol = 2*max.lag+2, nrow = nrow(pair.list)))	## matrix to store cross-correlation coefficients
    names(CCF) <- c(paste("lag", -max.lag:max.lag, sep=""), "numYears") ## adds column names
    
    for(k in 1:nrow(pair.list)){
      ## calculate number of years to base correlation on...
      CCF$numYears[k] <- length(na.omit((TS[,pair.list[k,1]]+TS[,pair.list[k,2]]))) ## add up number of times when sum can be done (i.e. neither is NA) then only include site with > 6 common years of data
      
      ## calculation of CCF...
      try(CCF[k,-ncol(CCF)] <- ccf(TS[,pair.list[k,1]], TS[,pair.list[k,2]], lag.max=max.lag, na.action=na.exclude, plot=F, type="correlation")$acf, silent=T) ## use try to prevent crashes if no data
    } ## this will take a long time!!! ## end k in nrow(pair.list)
    
    #hist(CCF$lag0)
    #hist(CCF$numYears)
    
    pair.attr <- pair.list ## matrix to hold attributes of pairs...
    colnames(pair.attr) <- c("site1","site2")
    pair.attr <- cbind(pair.attr, CCF)     ### add in correlation scores and number of comparisons at each site.
    head(pair.attr)
    
    ### drop sites comparisons with <7 years of survey data in common 
    if(nrow(pair.attr[pair.attr$numYears>6,])==0){
      print(paste("skip year"))
      next
    }else{
    new_pair_attr <- pair.attr[pair.attr$numYears>6,] 
    new_pair_attr <- new_pair_attr[!is.na(new_pair_attr$lag0),]
    new_pair_attr$mid.year<-mid.year
    new_pair_attr$start.year<-start.year
    new_pair_attr$end.year<-end.year
    
    # create a table of all pair-wise comparisons, this will be built up as we move through years
    all_pair_attr <- NULL
    all_pair_attr <- rbind(all_pair_attr, new_pair_attr) 
    
    ### record the number of good pair-wise comparisons and number of unique sites
    temp_summ_stats <- data.frame(total_comps = nrow(new_pair_attr), uni_sites = length(unique(c(new_pair_attr$site1, new_pair_attr$site2))))
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
}
### merge in predictor variables
head(final_pair_data)
head(final_summ_stats)

summary(final_pair_data)

write.csv(final_pair_data, file="../Data/Bird_sync_data/final_pair_data_all_spp_BBS_1.csv", row.names=FALSE) ## save final pair data for all 33 species

write.csv(final_summ_stats, file="../Data/Bird_sync_data/final_summ_stats_all_spp_BBS_1.csv", row.names=FALSE) ## save final summ stats for all 33 species


############## SCRIPT FOR THE SECOND 25% OF RANDOM SITE DATA ####################

rm(list=ls()) # clear R

### add data
final_data <- read.csv("../Data/Bird_sync_data/final_data_all_spp_BBS.csv", header=TRUE) ## 31 spp
woodland_bbs <- read.csv("../Data/Bird_sync_data/bbs_woodland_birds.csv", header=TRUE)
site_data_2 <- read.csv("../Data/BTO_data/pair_attr_mean_north_dist_BBS_2.csv", header=TRUE)

#############################
#### calculate synchrony ####
#############################

spp.list <- unique(final_data$name)

final_pair_data <- NULL
final_summ_stats <- NULL

### split based on species ###
for (g in spp.list){ # loop through spp.list 
  spp_data <- final_data[final_data$name==g,]
  total_comp <- NULL
  print(paste("species",g))    
  
  year.list<-1994:2007
  
  for (i in year.list){ # loop through years
    start.year<-i
    mid.year<-i+4.5
    print(paste("mid.year=",mid.year)) 
    end.year<-i+9
    species.10.yr.data<-spp_data[spp_data$year>=start.year&spp_data$year<=end.year,]
    
    ################################
    # create a matrix to be filled #
    site.list <- unique(species.10.yr.data$site)
    site.list
    year.list.temp <- min(species.10.yr.data$year):max(species.10.yr.data$year)
    year.list.temp
    Grow.matrix<-matrix(c(species.10.yr.data$gr), nrow=length(year.list.temp))
    Grow.matrix
    ncol(Grow.matrix)
    rownames(Grow.matrix)<-year.list.temp
    colnames(Grow.matrix)<-site.list
    Grow.matrix
    
    length(site.list)
    
    num.ts <- length(site.list)   # number of time series
    TS <- matrix(Grow.matrix, ncol=num.ts)    # stores Grow data in a matrix     (simply removes names!)
    TS
    # create site match table used to correct sites names after calculating synchrony
    site_match <- data.frame(TS_name = 1:ncol(TS), site_name = colnames(Grow.matrix))
    
    ###############################
    # cross-correlation functions #
    
    pair.list <- NULL             # calculates combination between all the pairs of sites
    for(k in 1:(num.ts-1)){
      pair.list <- rbind(pair.list, cbind(rep(k, num.ts-k), (k+1):num.ts)) ## species 508 has only 2 sites, so this line of code doesn't work...
    } # end k in num.ts-1
    nrow(pair.list) 
    
    colnames(pair.list) <- c("site1","site2")
    
    #### NOTE- ADD IN THE REAL SITE NAMES HERE, MERGE WITH SITE_DATA TABLE, THEN  CUT DOWN THE TABLE TO ONLY INCLUDE PAIRS WHERE distance < 100000 (100km)) ####
    ## e.g.   merge(data1,data2,by.x=c(col1,col2),by.y=c(site_a,site_b)) ##
    
    ###  correct site numbers ### 
    pair.list <- merge(pair.list, site_match, by.x="site1", by.y="TS_name")  
    pair.list <- pair.list[,c("site_name", "site2")]
    names(pair.list) <- gsub("site_name", "site1", names(pair.list))
    pair.list <- merge(pair.list, site_match, by.x="site2", by.y="TS_name")  
    pair.list <- pair.list[,c("site1", "site_name")]
    names(pair.list) <- gsub("site_name", "site2", names(pair.list))
    
    ## merge pair.list with site data
    site_data_2[,c(1,2)] <- sapply(site_data_2[,c(1,2)], as.factor)
    site_data_pair1 <- merge(pair.list, site_data_2, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))
    site_data_pair2 <- merge(pair.list, site_data_2, by.x=c("site1", "site2"), by.y=c("site_b", "site_a"))
    site_data_pair <- rbind(site_data_pair1, site_data_pair2)
    
    # unique(pair.list$site1)[!unique(pair.list$site1)%in%unique(site_data$site_a)]   # how many sites are missing attribute data
    
    ## remove site pair comparisons more than 100km apart
    site_data_pair <- site_data_pair[site_data_pair$distance <= 100000,]
    
    pair.list<-merge(pair.list,site_data_pair,by.x=c("site1", "site2"),by.y=c("site1", "site2"))[,1:2]
    
    
    ## calculate cross-correlation functions at different time lags
    ## (NB - we are only be interested in the current year, i.e. lag = 0)
    max.lag <- 0 # maximum time lag over which cross-correlations are assessed
    CCF <- data.frame(matrix(NA, ncol = 2*max.lag+2, nrow = nrow(pair.list)))	## matrix to store cross-correlation coefficients
    names(CCF) <- c(paste("lag", -max.lag:max.lag, sep=""), "numYears") ## adds column names
    
    for(k in 1:nrow(pair.list)){
      ## calculate number of years to base correlation on...
      CCF$numYears[k] <- length(na.omit((TS[,pair.list[k,1]]+TS[,pair.list[k,2]]))) ## add up number of times when sum can be done (i.e. neither is NA) then only include site with > 6 common years of data
      
      ## calculation of CCF...
      try(CCF[k,-ncol(CCF)] <- ccf(TS[,pair.list[k,1]], TS[,pair.list[k,2]], lag.max=max.lag, na.action=na.exclude, plot=F, type="correlation")$acf, silent=T) ## use try to prevent crashes if no data
    } ## this will take a long time!!! ## end k in nrow(pair.list)
    
    #hist(CCF$lag0)
    #hist(CCF$numYears)
    
    pair.attr <- pair.list ## matrix to hold attributes of pairs...
    colnames(pair.attr) <- c("site1","site2")
    pair.attr <- cbind(pair.attr, CCF)     ### add in correlation scores and number of comparisons at each site.
    head(pair.attr)
    
    ### drop sites comparisons with <7 years of survey data in common 
    new_pair_attr <- pair.attr[pair.attr$numYears>6,] 
    new_pair_attr <- new_pair_attr[!is.na(new_pair_attr$lag0),]
    new_pair_attr$mid.year<-mid.year
    new_pair_attr$start.year<-start.year
    new_pair_attr$end.year<-end.year
    
    # create a table of all pair-wise comparisons, this will be built up as we move through years
    all_pair_attr <- NULL
    all_pair_attr <- rbind(all_pair_attr, new_pair_attr) 
    
    ### record the number of good pair-wise comparisons and number of unique sites
    temp_summ_stats <- data.frame(total_comps = nrow(new_pair_attr), uni_sites = length(unique(c(new_pair_attr$site1, new_pair_attr$site2))))
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

write.csv(final_pair_data, file="../Data/Bird_sync_data/final_pair_data_all_spp_BBS_2.csv", row.names=FALSE) ## save final pair data for all 33 species

write.csv(final_summ_stats, file="../Data/Bird_sync_data/final_summ_stats_all_spp_BBS_2.csv", row.names=FALSE) ## save final summ stats for all 33 species


# ###########################
# ## abundance calculation ##
# ###########################
# woodland_bbs <- read.csv("../Data/Bird_sync_data/bbs_woodland_birds.csv", header=TRUE)
# 
# spp.list <- unique(woodland_bbs$species_code)
# 
# 
# abundance.results <- NULL
# 
# for (g in spp.list){ # loop for each species #
#   
#   species.tab<-woodland_bbs[woodland_bbs$species_code==g,] 
#   head(species.tab)
#   print(paste("species",g))
#   
#   # create a table to assess how much data in each year for that species
#   # then select only years that fulfill a minumum data criteria
#   # allocate those years to 'year.list'     
#   
#   year.list <-1994:2007   # temp until above steps are complete
#   
#   for (i in year.list){
#     
#     start.year<-i
#     mid.year<-i+4.5
#     print(paste("mid.year=",mid.year))
#     end.year<-i+9
#     species.10.yr.data<-species.tab[species.tab$Year>=start.year&species.tab$Year<=end.year,]
#     
#     species<-g
#     abundance.index <- mean(species.10.yr.data$TOT)
#     results.temp<-data.frame(start.year,mid.year,end.year,abundance.index,species)
#     abundance.results<-rbind(abundance.results,results.temp)
#     
#   }
# }
# 
# write.csv(abundance.results, file = "../Data/Bird_sync_data/abundance_data_BBS.csv", row.names = FALSE)
