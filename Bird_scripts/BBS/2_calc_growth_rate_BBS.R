###########################################################
## Title: Moving window calcualte growth rate BBS BTO data  
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: January 2018
##########################################################

rm(list=ls()) # clear R

## add data
woodland_bbs <- read.csv("../Data/Bird_sync_data/bbs_woodland_birds.csv", header=TRUE)

woodland_bbs <- unique(woodland_bbs) # drop duplicate rows

woodland_bbs$TOT <- woodland_bbs$TOT+1 # add 1 to each annual count to avoid problem of logging zeros 

####################################
#### calculate growth rate data ####
####################################

final_data <- NULL
zero_count_data <- NULL
spp_list <- unique(woodland_bbs$species_code)

for (g in spp_list){ # loop for each species #
  
  species.tab<-woodland_bbs[woodland_bbs$species_code==g,]
  head(species.tab)
  print(paste("species",g))
  
  new_bird_final<-NULL
  
  # id site list for 10yr moving window
  site.list <- unique(species.tab$site_code)
  gr <- NULL
  new_bird <- NULL
  
  for (j in site.list){
    temp.table <- species.tab[species.tab$site_code==j,] # make a mini table of the site in question
    gr <- NULL # clear gr column
    TOT <- NULL
    TOT <- as.numeric(TOT)
    
    # create a mini table per species, how many years are NOT!! zero counts (1 counts here as +1 above)#
    zero_count_data <- rbind(zero_count_data, data.frame(species_code = g, site_code = j, good_years = nrow(na.omit(temp.table[temp.table$TOT>1,])))) 
    
    for (i in (min(species.tab$YEAR):max(species.tab$YEAR))){ # loop through each year
      # convert to growth rate
      if(length(temp.table[temp.table$YEAR==i,"TOT"]) + length(temp.table[temp.table$YEAR==(i-1),"TOT"]) < 2){		# if there is less than two years of data
        gr <- c(gr,NA) # gr gets NA
      } else {
        gr <- c(gr, (log(temp.table[temp.table$YEAR==i,"TOT"]) - log(temp.table[temp.table$YEAR==(i-1),"TOT"]))) # else gr gets the log growth rate calculation
      }
      
      # take sindex and add NA if no count made
      if(length(temp.table[temp.table$YEAR==i,"TOT"]) == 1){
        TOT <- c(TOT, temp.table[temp.table$YEAR==i,"TOT"])
      } else {
        TOT <- c(TOT, NA)
      }
      
    }  # end i in  min(species.tab$Year):max(species.tab$Year)
    
    # sort out any NaNs, or Infs - make them NAs
    gr <- as.numeric(gsub("NaN", "NA", gr))
    gr <- as.numeric(gsub("Inf", "NA", gr))
    gr <- as.numeric(gsub("-Inf", "NA", gr))
    
    # create a mini table per site
    new_data <- data.frame(site = rep(j,length(min(species.tab$YEAR):max(species.tab$YEAR))), year = min(species.tab$YEAR):max(species.tab$YEAR), gr = gr, TOT = TOT)
    new_bird <- rbind(new_bird, new_data) # rbind the full growth rates per site per year with NAs added for the missing years
    
  }  # end j in site list
  
  new_bird_final <- rbind(new_bird_final, new_bird)
  new_bird_final$name <- g   # add in species code
  final_data <- rbind(final_data, new_bird_final)   #  all species data
  
} # end g in species

## 31 species and 6139 sites

### drop sites with <50% zero counts ###
good_year_data <- zero_count_data[zero_count_data$good_years>5,]
final_data$rec_id <- paste(final_data$name, final_data$site, sep="_")
good_year_data$rec_id <- paste(good_year_data$species_code, good_year_data$site_code, sep="_")

final_data <- final_data[final_data$rec_id%in%good_year_data$rec_id,]
## left with 31 species (none removed) and 4114 sites (2025 sites removed)

write.csv(final_data, file = "../Data/Bird_sync_data/final_data_all_spp_BBS.csv", row.names = FALSE)

