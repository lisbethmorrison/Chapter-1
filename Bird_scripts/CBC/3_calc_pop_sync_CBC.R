###########################################################
## Title: Moving window calcualte population synchrony CBC BTO data  
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: November 2017
##########################################################

#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script") ####

rm(list=ls()) # clear R
library(dplyr)

### add data
#final_data <- read.csv("../Data/Bird_sync_data/final_data_all_spp_CBC_zeros.csv", header=TRUE) 
#final_data <- read.csv("../Data/Bird_sync_data/final_data_all_spp_CBC_no_zeros.csv", header=TRUE) 
final_data <- read.csv("../Data/Bird_sync_data/final_data_all_spp_CBC_no_zeros2.csv", header=TRUE) 
woodland_cbc <- read.csv("../Data/Bird_sync_data/cbc_woodland_birds.csv", header=TRUE)
site_data <- read.csv("../Data/BTO_data/pair_attr_mean_north_dist_hab_sim_CBC.csv", header=TRUE)

### Create "good species list" which meet a minimum criteria of 2 filters
woodland_cbc<-woodland_cbc[woodland_cbc$Year>=1980&woodland_cbc$Year<=2000,]
woodland_cbc$species_code <- as.factor(woodland_cbc$species_code)
cbc_summary_tab<-with(woodland_cbc,table(species_code,Year))
cbc_summary_good_years<-NULL

for (i in 1:nrow(cbc_summary_tab)){
  no.good.years<-length(cbc_summary_tab[i,][cbc_summary_tab[i,]>=50])  # filter 1: 50 sites/year
  sp<-row.names(cbc_summary_tab)[i]
  results.temp<-data.frame(sp,no.good.years)
  cbc_summary_good_years<-rbind(cbc_summary_good_years,results.temp)
}
cbc_summary_good_years
# filter 2: only select species with >75% 'good years' (i.e. with more than 50 sites/year)
length(1980:2000)*0.75 # 15.75
good.species.list<-cbc_summary_good_years$sp[cbc_summary_good_years$no.good.years>15.75]
good.species.list
good.species.list <- droplevels(good.species.list)

############## CONCLUSION => ALL 38 SPECIES MEET THE ABOVE CRITERIA, SO ISN'T NEEDED ##################

#############################
#### calculate synchrony ####
#############################

final_data <- final_data[!final_data$name==508,] # removed siskin as not enough data (stops at first year, 1984.5)
final_data <- final_data[!final_data$name==460,] # remove pied flycatcher as very little data (stopped working at 1991.5)
final_data <- final_data[!final_data$name==522,] # removed hawfinch as not enough data (stops at first year, 1984.5)

spp.list <- unique(final_data$name)

final_pair_data <- NULL
final_summ_stats <- NULL

### split based on species ###
for (g in spp.list){ # loop through spp.list 
  spp_data <- final_data[final_data$name==g,]
  total_comp <- NULL
  print(paste("species",g))    
  
  year.list<-1980:1991
  
  for (i in year.list){ # loop through years
    start.year<-i
    mid.year<-i+4.5
    print(paste("mid.year=",mid.year)) 
    end.year<-i+9
    species.10.yr.data<-spp_data[spp_data$year>=start.year&spp_data$year<=end.year,]
    
    ################################
    # create a matrix to be filled #
    site.list <- unique(species.10.yr.data$site)
    year.list.temp <- min(species.10.yr.data$year):max(species.10.yr.data$year)
    Grow.matrix<-matrix(c(species.10.yr.data$gr), nrow=length(year.list.temp))
    ncol(Grow.matrix)
    rownames(Grow.matrix)<-year.list.temp
    colnames(Grow.matrix)<-site.list

    length(site.list)
    
    num.ts <- length(site.list)   # number of time series
    TS <- matrix(Grow.matrix, ncol=num.ts)    # stores Grow data in a matrix     (simply removes names!)
    # create site match table used to correct sites names after calculating synchrony
    site_match <- data.frame(TS_name = 1:ncol(TS), site_name = colnames(Grow.matrix))
    
    ###############################
    # cross-correlation functions #
    
    # calculates combination between all the pairs of sites
    pair.list <- t(combn(1:num.ts,2))
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
    site_data[,c(1,2)] <- sapply(site_data[,c(1,2)], as.factor)
    site_data_pair1 <- merge(pair.list, site_data, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))
    site_data_reverse <- site_data
    names(site_data_reverse)[1:6] <- c("site_b", "site_a", "site_b_EAST", "site_b_NORTH", "site_a_EAST", "site_a_NORTH")
    site_data_pair2 <- merge(pair.list, site_data_reverse, by.x=c("site1", "site2"), by.y=c("site_a", "site_b"))
    site_data_pair <- rbind(site_data_pair1, site_data_pair2)
    
    # unique(pair.list$site1)[!unique(pair.list$site1)%in%unique(site_data$site_a)]   # how many sites are missing attribute data
    
    ## remove site pair comparisons more than 100km apart
    site_data_pair <- site_data_pair[site_data_pair$distance <= 100,]
    
    pair.list<-merge(pair.list,site_data_pair,by.x=c("site1", "site2"),by.y=c("site1", "site2"))[,1:2] 
    
    # change site names back to numbers
    pair.list <- merge(pair.list, site_match, by.x="site1", by.y="site_name")
    pair.list <- pair.list[,c("TS_name", "site2")]
    names(pair.list) <- gsub("TS_name", "site1", names(pair.list))
    pair.list <- merge(pair.list, site_match, by.x="site2", by.y="site_name")
    pair.list <- pair.list[,c("site1", "TS_name")]
    names(pair.list) <- gsub("TS_name", "site2", names(pair.list))
    
    ## calculate cross-correlation functions at different time lags
    ## (NB - we are only be interested in the current year, i.e. lag = 0)
    max.lag <- 0 # maximum time lag over which cross-correlations are assessed
    CCF <- data.frame(matrix(NA, ncol = 2*max.lag+2, nrow = nrow(pair.list)))	## matrix to store cross-correlation coefficients
    names(CCF) <- c(paste("lag", -max.lag:max.lag, sep=""), "numYears") ## adds column names
    
    for(k in 1:nrow(pair.list)){
      # print(paste("row number=",k,"out of",nrow(pair.list))) 
      ## calculate number of years to base correlation on...
      CCF$numYears[k] <- length(na.omit((TS[,pair.list[k,1]]+TS[,pair.list[k,2]]))) ## add up number of times when sum can be done (i.e. neither is NA) then only include site with > 6 common years of data
      
      ## attempt calculation of CCF...
      try(CCF[k,-ncol(CCF)] <- ccf(TS[,pair.list[k,1]], TS[,pair.list[k,2]], lag.max=max.lag, na.action=na.exclude, plot=F, type="correlation")$acf, silent=T) ## use try to prevent crashes if no data
    } ## this will take a long time!!!
    
    ## rank by numYears
    CCF1 <- CCF
    CCF1 <- arrange(CCF1, desc(numYears))
    ## chop out rows where numYears<7
    CCF1 <- CCF1[CCF1$numYears>6,]
    CCF1 <- na.omit(CCF1)

    # if statement to skip species which won't run
    if (nrow(CCF1)<1){
      print(paste("skip species", g, "year", i+4.5))
      next
    }

    pair.attr <- pair.list ## matrix to hold attributes of pairs...
    colnames(pair.attr) <- c("site1","site2")
    pair.attr <- cbind(pair.attr, CCF)     ### add in correlation scores and number of comparisons at each site.
    head(pair.attr)
    
    ###  correct site numbers ### 
    pair.attr <- merge(pair.attr, site_match, by.x="site1", by.y="TS_name")  
    pair.attr <- pair.attr[,c("site_name", "site2", "lag0", "numYears")]
    names(pair.attr) <- gsub("site_name", "site1", names(pair.attr))
    pair.attr <- merge(pair.attr, site_match, by.x="site2", by.y="TS_name")  
    pair.attr <- pair.attr[,c("site1", "site_name", "lag0", "numYears")]
    names(pair.attr) <- gsub("site_name", "site2", names(pair.attr))
    
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

head(final_pair_data)
head(final_summ_stats)

length(unique(final_pair_data$spp)) ## 30 species

final_summ_stats$spp <- as.factor(final_summ_stats$spp) ## 32 species
final_pair_data_summ <- final_summ_stats %>% 
  group_by(spp) %>% 
  summarise(n = n()) ## nrow of each species
## only include species with complete time series (i.e. nrow=12)
final_pair_data_summ <- final_pair_data_summ[final_pair_data_summ$n>=12,] ## 27 species (5 species removed)
## merge back into final_pair_data
final_pair_data <- merge(final_pair_data, final_pair_data_summ, by="spp", all=FALSE)
length(unique(final_pair_data$spp)) ## 27 species

write.csv(final_pair_data, file="../Data/Bird_sync_data/final_pair_data_all_spp_CBC_no_zeros2_correct.csv", row.names=FALSE) ## save final pair data for all 33 species

write.csv(final_summ_stats, file="../Data/Bird_sync_data/final_summ_stats_all_spp_CBC_no_zeros2_correct.csv", row.names=FALSE) ## save final summ stats for all 33 species

########################## NOT USED ANYMORE ##############################

# ###########################
# ## abundance calculation ##
# ###########################
# woodland_cbc <- read.csv("../Data/Bird_sync_data/cbc_woodland_birds.csv", header=TRUE)
# 
# woodland_cbc <- woodland_cbc[!woodland_cbc$species_code==508,] # removed siskin as not enough data (stops at first year, 1984.5)
# woodland_cbc <- woodland_cbc[!woodland_cbc$species_code==460,] # remove pied flycatcher as very little data (stopped working at 1991.5)
# woodland_cbc <- woodland_cbc[!woodland_cbc$species_code==522,] # removed hawfinch as not enough data (stops at first year, 1984.5)
# woodland_cbc <- woodland_cbc[!woodland_cbc$species_code==57,] # removed capercaillie => gets removed at growth rate stage
# woodland_cbc <- woodland_cbc[!woodland_cbc$species_code==515,] # removed common crossbill =>gets removed at growth rate stage
# 
# spp.list <- unique(woodland_cbc$species_code)
# 
# 
# abundance.results <- NULL
# 
# for (g in spp.list){ # loop for each species #
#   
#   species.tab<-woodland_cbc[woodland_cbc$species_code==g,] 
#   head(species.tab)
#   print(paste("species",g))
#   
#   # create a table to assess how much data in each year for that species
#   # then select only years that fulfill a minumum data criteria
#   # allocate those years to 'year.list'     
#   
#   year.list <-1980:1991   # temp until above steps are complete
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
#     abundance.index <- mean(species.10.yr.data$Count)
#     results.temp<-data.frame(start.year,mid.year,end.year,abundance.index,species)
#     abundance.results<-rbind(abundance.results,results.temp)
#     
#   }
# }
# 
# write.csv(abundance.results, file = "../Data/Bird_sync_data/abundance_data_BTO.csv", row.names = FALSE)
