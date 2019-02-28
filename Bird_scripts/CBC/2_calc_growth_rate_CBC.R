###########################################################
## Title: Moving window calcualte growth rate CBC BTO data  
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: November 2017
##########################################################


#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script/Bird_scripts") ####

rm(list=ls()) # clear R

## add data
woodland_cbc <- read.csv("../Data/Bird_sync_data/cbc_woodland_birds.csv", header=TRUE) ## contains ALL 38 woodland species 

### delete unnessecary columns ###
woodland_cbc <- subset(woodland_cbc, select = -c(generalist, specialist, Area, qual, Gridref, hab, Hab1, malt))

woodland_cbc$Count <- woodland_cbc$Count+1 # add 1 to each annual count to avoid problem of logging zeros 

####################################
#### calculate growth rate data ####
####################################

final_data <- NULL
zero_count_data <- NULL
spp_list <- unique(woodland_cbc$species_code)

for (g in spp_list){ # loop for each species #
  
  species.tab<-woodland_cbc[woodland_cbc$species_code==g,]
  head(species.tab)
  print(paste("species",g))
  
  new_bird_final<-NULL

  # id site list for 10yr moving window
  site.list <- unique(species.tab$Pcode)
  gr <- NULL
  new_bird <- NULL
  
  for (j in site.list){
    temp.table <- species.tab[species.tab$Pcode==j,] # make a mini table of the site in question
    gr <- NULL # clear gr column
    Count <- NULL
    
    # create a mini table per species, how many years are NOT!! zero counts (1 counts here as +1 above)#
    zero_count_data <- rbind(zero_count_data, data.frame(species_code = g, Pcode = j, good_years = nrow(na.omit(temp.table[temp.table$Count>1,])))) 
    
    for (i in (min(species.tab$Year):max(species.tab$Year))){ # loop through each year
      # convert to growth rate
      if(length(temp.table[temp.table$Year==i,"Count"]) + length(temp.table[temp.table$Year==(i-1),"Count"]) < 2){		# if there is less than two years of data
        gr <- c(gr,NA) # gr gets NA
      } else if (temp.table[temp.table$Year==(i-1), "Count"] == 1){ ## if previous year == 1 (i.e. no count made), gr = NA to avoid synchrony being calculated on string of zeros
        gr <- c(gr, NA) # gr gets NA
      } else {
        gr <- c(gr, (log(temp.table[temp.table$Year==i,"Count"]) - log(temp.table[temp.table$Year==(i-1),"Count"]))) # else gr gets the log growth rate calculation
      }
      
      # take sindex and add NA if no count made
      if(length(temp.table[temp.table$Year==i,"Count"]) == 1){
        Count <- c(Count, temp.table[temp.table$Year==i,"Count"])
      } else {
        Count <- c(Count, NA)
      }
      
    }  # end i in  min(species.tab$Year):max(species.tab$Year)
    
    # sort out any NaNs, or Infs - make them NAs
    gr <- as.numeric(gsub("NaN", "NA", gr))
    gr <- as.numeric(gsub("Inf", "NA", gr))
    gr <- as.numeric(gsub("-Inf", "NA", gr))
    
    # create a mini table per site
    new_data <- data.frame(site = rep(j,length(min(species.tab$Year):max(species.tab$Year))), year = min(species.tab$Year):max(species.tab$Year), gr = gr, Count = Count)
    new_bird <- rbind(new_bird, new_data) # rbind the full growth rates per site per year with NAs added for the missing years
    
  }  # end j in site list
  
  new_bird_final <- rbind(new_bird_final, new_bird)
  new_bird_final$name <- g   # add in species code
  final_data <- rbind(final_data, new_bird_final)   #  all species data
  
} # end g in species

## 38 species and 411 sites

### drop sites with >50% zero counts ###
good_year_data <- zero_count_data[zero_count_data$good_years>5,] ## dataframe with species & site combo with more than 5 years of non-zero counts
final_data$rec_id <- paste(final_data$name, final_data$site, sep="_")
good_year_data$rec_id <- paste(good_year_data$species_code, good_year_data$Pcode, sep="_")

final_data <- final_data[final_data$rec_id%in%good_year_data$rec_id,]
## NOTE: this process removes 2 species from the analysis: spp 57 (Capercaillie) & 515 (Common crossbill)
## Also removes 231 sites (left with 180 sites)

write.csv(final_data, file = "../Data/Bird_sync_data/final_data_all_spp_CBC.csv", row.names = FALSE)

