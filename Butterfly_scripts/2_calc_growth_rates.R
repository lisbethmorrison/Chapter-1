###########################################################
## Title: Moving window growth rate calculation
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: August 2017
##########################################################


#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script/Butterfly_scripts") ####

rm(list=ls()) # clear R

### add data
b_data <- read.csv("../Data/Butterfly_sync_data/b_data.csv", header = TRUE) # add butterfly count data

b_data <- unique(b_data) # drop duplicate rows

b_data$SINDEX<-b_data$SINDEX+1  # add 1 to each annual count to avoid problem of logging zeros

#### calculate growth rate data ####

final_data <- NULL
zero_count_data <- NULL

for (g in unique(b_data$SPECIES)){ # loop for each species #

species.tab<-b_data[b_data$SPECIES==g,] 
head(species.tab)
print(paste("species",g))

new_butterfly_final<-NULL

# id site list for 10yr moving window
site.list <- unique(species.tab$SITE)
gr <- NULL
new_butterfly <- NULL

for (j in site.list){
  temp.table <- species.tab[species.tab$SITE==j,] # make a mini table of the site in question
	gr <- NULL # clear gr column
	SINDEX <- NULL
	
	# create a mini table per species, how many years are NOT!! zero counts (1 counts here as +1 above)#
	zero_count_data <- rbind(zero_count_data, data.frame(SPECIES = g, SITE = j, good_years = nrow(na.omit(temp.table[temp.table$SINDEX>1,])))) 
	
	for (i in (min(species.tab$YEAR):max(species.tab$YEAR))){ # loop through each year
	  print(paste("Year=", i))
		# convert to growth rate
		if(length(temp.table[temp.table$YEAR==i,"SINDEX"]) + length(temp.table[temp.table$YEAR==(i-1),"SINDEX"]) < 2){		# if there is less than two years of data
			gr <- c(gr,NA) # gr gets NA
		} else if (temp.table[temp.table$YEAR==(i-1), "SINDEX"] == 1){ ## if previous year == 1 (i.e. no count made), gr = NA to avoid synchrony being calculated on string of zeros
		  gr <- c(gr, NA) # gr gets NA
		} else {
			gr <- c(gr, (log(temp.table[temp.table$YEAR==i,"SINDEX"]) - log(temp.table[temp.table$YEAR==(i-1),"SINDEX"]))) # else gr gets the log growth rate calculation
		}
		
		# take sindex and add NA if no count made
		if(length(temp.table[temp.table$YEAR==i,"SINDEX"]) == 1){
			SINDEX <- c(SINDEX, temp.table[temp.table$YEAR==i,"SINDEX"])
		} else {
			SINDEX <- c(SINDEX, NA)
		}
	  
	}  # end i in  min(species.tab$YEAR):max(species.tab$YEAR)
	
	# sort out any NaNs, or Infs - make them NAs
	gr <- as.numeric(gsub("NaN", "NA", gr))
	gr <- as.numeric(gsub("Inf", "NA", gr))
	gr <- as.numeric(gsub("-Inf", "NA", gr))
	
	# create a mini table per site
	new_data <- data.frame(site = rep(j,length(min(species.tab$YEAR):max(species.tab$YEAR))), year = min(species.tab$YEAR):max(species.tab$YEAR), gr = gr, SINDEX = SINDEX)
	new_butterfly <- rbind(new_butterfly, new_data) # rbind the full growth rates per site per year with NAs added for the missing years

	}  # end j in site list

new_butterfly_final <- rbind(new_butterfly_final, new_butterfly)
new_butterfly_final$name <- g   # add in species code
final_data <- rbind(final_data, new_butterfly_final)   #  all species data

} # end g in species

write.csv(final_data, file = "../Data/Butterfly_sync_data/final_data_all_spp.csv", row.names = FALSE)
length(unique(final_data$name)) # 60 species
length(unique(final_data$site)) # 2415 sites
