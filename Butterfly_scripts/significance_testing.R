###########################################################
## Title: Testing model significance & creating barcharts
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: November 2017
##########################################################

#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script") ####

rm(list=ls()) # clear R

library(lme4)
library(lmerTest)
library(dplyr)
library(foreach)
library(doParallel)
library(lmerTest)
library(lme4)
options(scipen=999)

## load data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE)

length(unique(pair_attr$spp)) # 32 species

###### create 3 new pair_attr files whicih compares early, late and overall
pair_attr_1985 <- pair_attr[pair_attr$mid.year==1984.5,]
pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]

pair_attr_early <- rbind(pair_attr_1985, pair_attr_2000) # comparison of early years 1985 & 1985
pair_attr_late <- rbind(pair_attr_2000, pair_attr_2012) # comparison of late years 1985 & 2012
pair_attr_overall <- rbind(pair_attr_1985, pair_attr_2012) # comparison of overall years 1985 & 2012

## centre and standardise continuous variables by mean and SD
## early
pair_attr_early$distance <- (pair_attr_early$distance - mean(na.omit(pair_attr_early$distance)))/sd(na.omit(pair_attr_early$distance))
pair_attr_early$mean_northing <- (pair_attr_early$mean_northing - mean(na.omit(pair_attr_early$mean_northing)))/sd(na.omit(pair_attr_early$mean_northing))
pair_attr_early$renk_hab_sim <- (pair_attr_early$renk_hab_sim - mean(na.omit(pair_attr_early$renk_hab_sim)))/sd(na.omit(pair_attr_early$renk_hab_sim))
## late
pair_attr_late$distance <- (pair_attr_late$distance - mean(na.omit(pair_attr_late$distance)))/sd(na.omit(pair_attr_late$distance))
pair_attr_late$mean_northing <- (pair_attr_late$mean_northing - mean(na.omit(pair_attr_late$mean_northing)))/sd(na.omit(pair_attr_late$mean_northing))
pair_attr_late$renk_hab_sim <- (pair_attr_late$renk_hab_sim - mean(na.omit(pair_attr_late$renk_hab_sim)))/sd(na.omit(pair_attr_late$renk_hab_sim))
## overall
pair_attr_overall$distance <- (pair_attr_overall$distance - mean(na.omit(pair_attr_overall$distance)))/sd(na.omit(pair_attr_overall$distance))
pair_attr_overall$mean_northing <- (pair_attr_overall$mean_northing - mean(na.omit(pair_attr_overall$mean_northing)))/sd(na.omit(pair_attr_overall$mean_northing))
pair_attr_overall$renk_hab_sim <- (pair_attr_overall$renk_hab_sim - mean(na.omit(pair_attr_overall$renk_hab_sim)))/sd(na.omit(pair_attr_overall$renk_hab_sim))

#### write files
write.csv(pair_attr_early, file = "../Data/Butterfly_sync_data/pair_attr_early.csv", row.names = FALSE) # save pair_attr_early file 
write.csv(pair_attr_late, file = "../Data/Butterfly_sync_data/pair_attr_late.csv", row.names = FALSE) # save pair_attr_late file 
write.csv(pair_attr_overall, file = "../Data/Butterfly_sync_data/pair_attr_overall.csv", row.names = FALSE) # save pair_attr_overall file 

## read in files to save time 
pair_attr_early <- read.csv("../Data/Butterfly_sync_data/pair_attr_early.csv", header=TRUE)
pair_attr_late <- read.csv("../Data/Butterfly_sync_data/pair_attr_late.csv", header=TRUE)
pair_attr_overall <- read.csv("../Data/Butterfly_sync_data/pair_attr_overall.csv", header=TRUE)

#### change variables to factors/characters
pair_attr_early$mid.year <- as.factor(pair_attr_early$mid.year)
pair_attr_early$pair.id <- as.character(pair_attr_early$pair.id)
pair_attr_early$spp <- as.factor(pair_attr_early$spp)

pair_attr_late$mid.year <- as.factor(pair_attr_late$mid.year)
pair_attr_late$pair.id <- as.character(pair_attr_late$pair.id)
pair_attr_late$spp <- as.factor(pair_attr_late$spp)

pair_attr_overall$mid.year <- as.factor(pair_attr_overall$mid.year)
pair_attr_overall$pair.id <- as.character(pair_attr_overall$pair.id)
pair_attr_overall$spp <- as.factor(pair_attr_overall$spp)


##################################
## EARLY MODEL FOR EACH SPECIES ##
##################################

## run lmer for each species and save F values to be used in permutation tests

spp.list.early <- unique(pair_attr_early$spp)

results_table_early<-NULL
F_result_final <- NULL
for (i in spp.list.early){
  print(i)
  
  ## create unique pair_attr for each species
  p_attr <- pair_attr_early[pair_attr_early$spp==i,]
  
# if loop (if nrow when start.year=1980 is <20 the skip that species)
  if(nrow(p_attr[p_attr$start.year=="1980",]) < 20){
    print(paste("skip species",i))
    next 
    
  } else {
  
  early_model <- (lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data = pair_attr_early[pair_attr_early$spp==i,]))
  summary(early_model)
  anova(early_model)
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(early_model)$coefficients[,1:5],i) ## save model results
  F_result_table <- data.frame(anova(early_model)[5], i) ## also save model F values to use in permutation tests
  results_table_early <-rbind(results_table_early,results_table_temp)
  F_result_final <-rbind(F_result_final,F_result_table)
  
  }
    
}

## if else function removes 1 species when nrow of 1980 is set at <20:
## species 110, 116, 119, 17, 18, 23, 27 and 48
## 24 species 

## save results
names(results_table_early) <- c("Estimate", "SD", "df", "t","p_value", "species")
results_table_early$parameter <- paste(row.names(results_table_early))
rownames(results_table_early) <- 1:nrow(results_table_early)
## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_table_early <- results_table_early[grep("mid.year", results_table_early$parameter),]
## change parameter names to year of interest
results_table_early$year <- "1985/2000" 
## make final table
results_final_early <- results_table_early[,c(1,6,8)]

##### modify F result table
colnames(F_result_final)[2] <- "j"
F_result_final$i <- 0 ## make i column with zeros 
F_result_final$parameter <- paste(row.names(F_result_final)) ## move row.names to parameter column
rownames(F_result_final) <- 1:nrow(F_result_final) ## change row names to numbers
## take out 3 columns for each species: mean northing, distance and renk_hab_sim
F_result_final <- F_result_final[grep("mid.year", F_result_final$parameter),]

## save true model results
write.csv(F_result_final, file = "../Results/Model_outputs/UKBMS/true_F_values_spp_85_00.csv", row.names=FALSE)
F_result_final <- read.csv("../Results/Model_outputs/UKBMS/true_F_values_spp_85_00.csv", header=TRUE)

## permutation tests for each species

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
F_result_final$j <- as.factor(F_result_final$j)
spp_list <- unique(F_result_final$j) ## only run this on the species that work (24)
## merge spp with pair_attr_early
spp_list <- as.data.frame(spp_list)
pair_attr_early <- merge(pair_attr_early, spp_list, by.x="spp", by.y="spp_list", all=FALSE)
pair_attr_early <- droplevels(pair_attr_early)
length(unique(pair_attr_early$spp))
spp_list <- unique(F_result_final$j) ## only run this on the species that work (24)
str(spp_list)

para_start_time = Sys.time()
ukbms_spp_para <- foreach (i=1:n_sims,  .combine=rbind, .packages=c('lme4', 'foreach')) %dopar% { ## loop through each permutation
  print(i)
  foreach (j=spp_list, .combine=rbind, .packages=c('lme4', 'foreach')) %dopar% { ## loop through each species
    print(j)
    pair_attr_early$mid.year_shuffle <- sample(pair_attr_early$mid.year) ## randomly shuffle mid year variable
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year_shuffle + (1|pair.id), data = pair_attr_early[pair_attr_early$spp==j,])
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i, j=j)
  }
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 1.7hrs for 999 runs

### save results
ukbms_spp_para$parameter <- paste(row.names(ukbms_spp_para)) ## move row.names to parameter column
rownames(ukbms_spp_para) <- 1:nrow(ukbms_spp_para) ## change row names to numbers
ukbms_spp_para <- ukbms_spp_para[,-c(1:3)] ## remove unnecessary columns
## only keep rows with interaction
ukbms_spp_para <- ukbms_spp_para[grep("mid.year_shuffle", ukbms_spp_para$parameter),]
final_results_table <- rbind(F_result_final, ukbms_spp_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_F_values_spp_85_00.csv", row.names=FALSE)
## read in file
perm_ukbms_sp <- read.csv("../Results/Model_outputs/UKBMS/perm_F_values_spp_85_00.csv", header=TRUE) 

## calculate p value for each species
## Calculate p value
number_of_permutations <- 1000
spp_list <- unique(perm_ukbms_sp$j)
p_values_final <- NULL
## loop through each species
for(j in spp_list){
  
  perm_ukbms_sp2 <- perm_ukbms_sp[perm_ukbms_sp$j==j,] ## select species of interest
  true_ukbms_sp <- perm_ukbms_sp2[perm_ukbms_sp2$i==0,] ## obtain true F value (iteration = 0)
  perm_ukbms_sp2 <- perm_ukbms_sp2[!perm_ukbms_sp2$i==0,] ## remove true value to calc. p value
  diff.observed <- true_ukbms_sp$F.value ## true F value
  
  # P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
  pvalue = sum(abs(perm_ukbms_sp2$F.value) >= abs(diff.observed)) / number_of_permutations
  p_values_temp <- data.frame(pvalue=pvalue, species=j)
  p_values_final <- rbind(p_values_final, p_values_temp)
}
## save file
write.csv(p_values_final, file = "../Results/Model_outputs/UKBMS/perm_p_values_spp_85_00.csv", row.names=FALSE)
p_values_early <- read.csv("../Results/Model_outputs/UKBMS/perm_p_values_spp_85_00.csv", header=TRUE)

## merge permutation p values with main model output (need results_final_early to know direction of significant change)
final_results_early <- merge(results_final_early, p_values_early, by="species")

## classify each species as either no change (insignificant p value), decreasing (negative estimte), or increasing (positive estimate)
final_results_early <- final_results_early %>% group_by(species) %>% mutate(change = ifelse(pvalue>0.05, "No change", ifelse(pvalue<0.05 & Estimate<0, "Decrease", "Increase")))

nrow(final_results_early[final_results_early$change=="No change",]) ## 6 species unchanged
nrow(final_results_early[final_results_early$change=="Decrease",]) ## 17 species decreasing
nrow(final_results_early[final_results_early$change=="Increase",]) ## 1 species increasing

write.csv(final_results_early, file="../Results/Butterfly_results/ukbms_spp_trends_85_00.csv", row.names=FALSE)

final_results_early2 <- data.frame(change=unique(final_results_early$change), no_species=c(17,6,1), total_species=24)
final_results_early2$percentage <- (final_results_early2$no_species / final_results_early2$total_species)*100

## save file
write.csv(final_results_early2, file="../Results/Butterfly_results/ukbms_percentage_trends_85_00.csv", row.names=FALSE)


##################################
## LATE MODEL FOR EACH SPECIES ##
#################################

spp.list.late <- unique(pair_attr_late$spp)
results_table_late<-NULL
F_result_final <- NULL
for (i in spp.list.late){
  print(i)
  
  ## create unique pair_attr for each species
  p_attr <- pair_attr_late[pair_attr_late$spp==i,]
  
  # if loop (if nrow when start.year=2007 is <20 the skip that species)
  if(nrow(p_attr[p_attr$start.year=="2007",]) < 20){
    print(paste("skip species",i))
    next 
    
  } else if(length(unique(p_attr$pair.id))==nrow(p_attr)) {
    print(paste("skip species",i))
    next
  }
  
  late_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data = pair_attr_late[pair_attr_late$spp==i,])
  summary(late_model)
  anova(late_model)
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(late_model)$coefficients[,1:5],i)
  results_table_late <-rbind(results_table_late,results_table_temp)
  F_result_table <- data.frame(anova(late_model)[5], i) ## also save model F values to use in permutation tests
  F_result_final <-rbind(F_result_final,F_result_table)
  
}

## if else function removes 2 species when nrow of 2007 is set at <20:
## species 18 (pearl-bordered fritillary) skipped
## 31 species

names(results_table_late) <- c("Estimate", "SD", "df", "t","p_value", "species")
results_table_late$parameter <- paste(row.names(results_table_late))
rownames(results_table_late) <- 1:nrow(results_table_late)
## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_table_late <- results_table_late[grep("mid.year", results_table_late$parameter),]
## change parameter names to year of interest (mid.year = 2012/2011.5)
results_table_late$year <- "2000/2012"  
results_final_late <- results_table_late[,c(1,6,8)]


##### modify F result table
colnames(F_result_final)[2] <- "j"
F_result_final$i <- 0 ## make i column with zeros 
F_result_final$parameter <- paste(row.names(F_result_final)) ## move row.names to parameter column
rownames(F_result_final) <- 1:nrow(F_result_final) ## change row names to numbers
## take out 3 columns for each species: mean northing, distance and renk_hab_sim
F_result_final <- F_result_final[grep("mid.year", F_result_final$parameter),]

## save true model results
write.csv(F_result_final, file = "../Results/Model_outputs/UKBMS/true_F_values_spp_00_12.csv", row.names=FALSE)
F_result_final <- read.csv("../Results/Model_outputs/UKBMS/true_F_values_spp_00_12.csv", header=TRUE)

## permutation tests for each species

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
F_result_final$j <- as.factor(F_result_final$j)
spp_list <- unique(F_result_final$j) ## only run this on the species that work (31)
## merge spp with pair_attr_late
spp_list <- as.data.frame(spp_list)
pair_attr_late <- merge(pair_attr_late, spp_list, by.x="spp", by.y="spp_list", all=FALSE)
pair_attr_late <- droplevels(pair_attr_late)
length(unique(pair_attr_late$spp))
spp_list <- unique(F_result_final$j) ## only run this on the species that work (31)
str(spp_list)

para_start_time = Sys.time()
ukbms_spp_para <- foreach (i=1:n_sims,  .combine=rbind, .packages=c('lme4', 'foreach')) %dopar% { ## loop through each permutation
  print(i)
  foreach (j=spp_list, .combine=rbind, .packages=c('lme4', 'foreach')) %dopar% { ## loop through each species
    print(j)
    pair_attr_late$mid.year_shuffle <- sample(pair_attr_late$mid.year) ## randomly shuffle mid year variable
    model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year_shuffle + (1|pair.id), data = pair_attr_late[pair_attr_late$spp==j,])
    ## run model with shuffled variable
    ## save results
    anoresult<-anova(model)
    data.frame(anoresult, i=i, j=j)
  }
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 9.618 minutes for 999 runs!

### save results
ukbms_spp_para$parameter <- paste(row.names(ukbms_spp_para)) ## move row.names to parameter column
rownames(ukbms_spp_para) <- 1:nrow(ukbms_spp_para) ## change row names to numbers
ukbms_spp_para <- ukbms_spp_para[,-c(1:3)] ## remove unnecessary columns
## only keep rows with interaction
ukbms_spp_para <- ukbms_spp_para[grep("mid.year_shuffle", ukbms_spp_para$parameter),]
final_results_table <- rbind(F_result_final, ukbms_spp_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_F_values_spp_00_12.csv", row.names=FALSE)
## read in file
perm_ukbms_sp <- read.csv("../Results/Model_outputs/UKBMS/perm_F_values_spp_00_12.csv", header=TRUE) 

## calculate p value for each species
## Calculate p value
number_of_permutations <- 1000
spp_list <- unique(perm_ukbms_sp$j)
p_values_final <- NULL
## loop through each species
for(j in spp_list){
  
  perm_ukbms_sp2 <- perm_ukbms_sp[perm_ukbms_sp$j==j,] ## select species of interest
  true_ukbms_sp <- perm_ukbms_sp2[perm_ukbms_sp2$i==0,] ## obtain true F value (iteration = 0)
  perm_ukbms_sp2 <- perm_ukbms_sp2[!perm_ukbms_sp2$i==0,] ## remove true value to calc. p value
  diff.observed <- true_ukbms_sp$F.value ## true F value
  
  # P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
  pvalue = sum(abs(perm_ukbms_sp2$F.value) >= abs(diff.observed)) / number_of_permutations
  p_values_temp <- data.frame(pvalue=pvalue, species=j)
  p_values_final <- rbind(p_values_final, p_values_temp)
}
## save file
write.csv(p_values_final, file = "../Results/Model_outputs/UKBMS/perm_p_values_spp_00_12.csv", row.names=FALSE)
p_values_late <- read.csv("../Results/Model_outputs/UKBMS/perm_p_values_spp_00_12.csv", header=TRUE)

## merge permutation p values with main model output
final_results_late <- merge(results_final_late, p_values_late, by="species")

## classify each species as either no change (insignificant p value), decreasing (negative estimte), or increasing (positive estimate)
final_results_late <- final_results_late %>% group_by(species) %>% mutate(change = ifelse(pvalue>0.05, "No change", ifelse(pvalue<0.05 & Estimate<0, "Decrease", "Increase")))

nrow(final_results_late[final_results_late$change=="No change",]) ## 5 species unchanged
nrow(final_results_late[final_results_late$change=="Decrease",]) ## 3 species decreasing
nrow(final_results_late[final_results_late$change=="Increase",]) ## 23 species increasing

write.csv(final_results_late, file="../Results/Butterfly_results/ukbms_spp_trends_00_12.csv", row.names=FALSE)

final_results_late2 <- data.frame(change=unique(final_results_late$change), no_species=c(23,5,3), total_species=31)
final_results_late2$percentage <- (final_results_late2$no_species / final_results_late2$total_species)*100

## save file
write.csv(final_results_late2, file="../Results/Butterfly_results/ukbms_percentage_trends_00_12.csv", row.names=FALSE)



####################################
## OVERALL MODEL FOR EACH SPECIES ##
####################################

spp.list.overall <- unique(pair_attr_overall$spp)
results_table_overall<-NULL
F_result_final <- NULL
for (i in spp.list.overall){
  print(i)
  
  ## create unique pair_attr for each species
  p_attr <- pair_attr_overall[pair_attr_overall$spp==i,]
  
  # if loop (if a is < 30 and b is < 30)
  if((nrow(p_attr[p_attr$start.year=="1980",]) < 20 || nrow(p_attr[p_attr$start.year=="2007",]) < 20)){
  print(paste("skip species",i))
    
    next 
    
  } else {
    
  if((length(unique(p_attr$pair.id)))==(nrow(p_attr))){
    print(paste("skip species",i))
    
    next
    
  } else {
  
  overall_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data = pair_attr_overall[pair_attr_overall$spp==i,])
  summary(overall_model)
  anova(overall_model)
  

  ### save and plot the results ###
  results_table_temp <- data.frame(summary(overall_model)$coefficients[,1:5],i)
  results_table_overall <-rbind(results_table_overall,results_table_temp)
  F_result_table <- data.frame(anova(overall_model)[5], i) ## also save model F values to use in permutation tests
  F_result_final <-rbind(F_result_final,F_result_table)
  
  }

  }
}
## if else function removes 6 species when nrow of 1980 and 2007 is set at <20:
## spp. 110, 116, 119, 17, 18, 23, 27, 4, 48, 64 and 94 skipped
## 21 species left

names(results_table_overall) <- c("Estimate", "SD", "df", "t","p_value", "species")
results_table_overall$parameter <- paste(row.names(results_table_overall))
rownames(results_table_overall) <- 1:nrow(results_table_overall)
## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_table_overall <- results_table_overall[grep("mid.year", results_table_overall$parameter),]
## change parameter names to year of interest (mid.year = 2012/2011.5)
results_table_overall$year <- "1985/2012"   
results_final_overall <- results_table_overall[,c(1,6,8)]

##### modify F result table
colnames(F_result_final)[2] <- "j"
F_result_final$i <- 0 ## make i column with zeros 
F_result_final$parameter <- paste(row.names(F_result_final)) ## move row.names to parameter column
rownames(F_result_final) <- 1:nrow(F_result_final) ## change row names to numbers
## take out 3 columns for each species: mean northing, distance and renk_hab_sim
F_result_final <- F_result_final[grep("mid.year", F_result_final$parameter),]

## save true model results
write.csv(F_result_final, file = "../Results/Model_outputs/UKBMS/true_F_values_spp_85_12.csv", row.names=FALSE)
F_result_final <- read.csv("../Results/Model_outputs/UKBMS/true_F_values_spp_85_12.csv", header=TRUE)

## permutation tests for each species

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
F_result_final$j <- as.factor(F_result_final$j)
spp_list <- unique(F_result_final$j) ## only run this on the species that work (21)
## merge spp with pair_attr_overall
spp_list <- as.data.frame(spp_list)
pair_attr_overall <- merge(pair_attr_overall, spp_list, by.x="spp", by.y="spp_list", all=FALSE)
pair_attr_overall <- droplevels(pair_attr_overall)
length(unique(pair_attr_overall$spp))
spp_list <- unique(F_result_final$j) ## only run this on the species that work (21)
str(spp_list)

para_start_time = Sys.time()
ukbms_spp_para <- foreach (i=1:n_sims,  .combine=rbind, .packages=c('lme4', 'foreach')) %dopar% { ## loop through each permutation
  print(i)
  foreach (j=spp_list, .combine=rbind, .packages=c('lme4', 'foreach')) %dopar% { ## loop through each species
    print(j)
    pair_attr_overall$mid.year_shuffle <- sample(pair_attr_overall$mid.year) ## randomly shuffle mid year variable
    model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year_shuffle + (1|pair.id), data = pair_attr_overall[pair_attr_overall$spp==j,])
    ## run model with shuffled variable
    ## save results
    anoresult<-anova(model)
    data.frame(anoresult, i=i, j=j)
  }
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 1hr for 999 runs

### save results
ukbms_spp_para$parameter <- paste(row.names(ukbms_spp_para)) ## move row.names to parameter column
rownames(ukbms_spp_para) <- 1:nrow(ukbms_spp_para) ## change row names to numbers
ukbms_spp_para <- ukbms_spp_para[,-c(1:3)] ## remove unnecessary columns
## only keep rows with interaction
ukbms_spp_para <- ukbms_spp_para[grep("mid.year_shuffle", ukbms_spp_para$parameter),]
final_results_table <- rbind(F_result_final, ukbms_spp_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_F_values_spp_85_12.csv", row.names=FALSE)
## read in file
perm_ukbms_sp <- read.csv("../Results/Model_outputs/UKBMS/perm_F_values_spp_85_12.csv", header=TRUE) 

## calculate p value for each species
## Calculate p value
number_of_permutations <- 1000
spp_list <- unique(perm_ukbms_sp$j)
p_values_final <- NULL
## loop through each species
for(j in spp_list){
  
  perm_ukbms_sp2 <- perm_ukbms_sp[perm_ukbms_sp$j==j,] ## select species of interest
  true_ukbms_sp <- perm_ukbms_sp2[perm_ukbms_sp2$i==0,] ## obtain true F value (iteration = 0)
  perm_ukbms_sp2 <- perm_ukbms_sp2[!perm_ukbms_sp2$i==0,] ## remove true value to calc. p value
  diff.observed <- true_ukbms_sp$F.value ## true F value
  
  # P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
  pvalue = sum(abs(perm_ukbms_sp2$F.value) >= abs(diff.observed)) / number_of_permutations
  p_values_temp <- data.frame(pvalue=pvalue, species=j)
  p_values_final <- rbind(p_values_final, p_values_temp)
}
## save file
write.csv(p_values_final, file = "../Results/Model_outputs/UKBMS/perm_p_values_spp_85_12.csv", row.names=FALSE)
p_values_overall <- read.csv("../Results/Model_outputs/UKBMS/perm_p_values_spp_85_12.csv", header=TRUE)

## merge permutation p values with main model output
final_results_overall <- merge(results_final_overall, p_values_overall, by="species")

## classify each species as either no change (insignificant p value), decreasing (negative estimte), or increasing (positive estimate)
final_results_overall <- final_results_overall %>% group_by(species) %>% mutate(change = ifelse(pvalue>0.05, "No change", ifelse(pvalue<0.05 & Estimate<0, "Decrease", "Increase")))

nrow(final_results_overall[final_results_overall$change=="No change",]) ## 6 species unchanged
nrow(final_results_overall[final_results_overall$change=="Decrease",]) ## 6 species decreasing
nrow(final_results_overall[final_results_overall$change=="Increase",]) ## 9 species increasing

write.csv(final_results_overall, file="../Results/Butterfly_results/ukbms_spp_trends_85_12.csv", row.names=FALSE)

final_results_overall2 <- data.frame(change=unique(final_results_overall$change), no_species=c(9,6,6), total_species=21)
final_results_overall2$percentage <- (final_results_overall2$no_species / final_results_overall2$total_species)*100

## save file
write.csv(final_results_overall2, file="../Results/Butterfly_results/ukbms_percentage_trends_85_12.csv", row.names=FALSE)




######################### PLOT GRAPHS #############################

#### save comparison model results for each comparison for each species #####
results_final_early <- read.csv("../Results/Butterfly_results/ukbms_percentage_trends_85_00.csv", header=TRUE)
results_final_late <- read.csv("../Results/Butterfly_results/ukbms_percentage_trends_00_12.csv", header=TRUE)
results_final_overall <- read.csv("../Results/Butterfly_results/ukbms_percentage_trends_85_12.csv", header=TRUE)

## remove columns 2 and 3 of each data frame (so just left with change and percentages)
results_final_early <- results_final_early[,-c(2:3)]
results_final_late <- results_final_late[,-c(2:3)]
results_final_overall <- results_final_overall[,-c(2:3)]

results_final_early$comparison <- "1985-2000"
results_final_late$comparison <- "2000-2012"
results_final_overall$comparison <- "1985-2012"

## rbind all results together 
model_comp_results <- rbind(results_final_early, results_final_late, results_final_overall)
## save file
write.csv(model_comp_results, file="../Results/Butterfly_results/ukbms_percentage_trends.csv", row.names=FALSE)

## now produce graph

library(ggplot2)

## change levels of change so they go in correct order
levels(model_comp_results$change)
model_comp_results$change <- factor(model_comp_results$change, levels=c("Increase", "No change", "Decrease"))
levels(model_comp_results$change)
model_comp_results$comparison <- as.factor(model_comp_results$comparison)
levels(model_comp_results$comparison)
model_comp_results$comparison <- factor(model_comp_results$comparison, levels=c("1985-2000", "2000-2012", "1985-2012"))
levels(model_comp_results$comparison)

## pretty graph
png("../Graphs/Percentage_trends/ukbms_percentage_trends.png", height = 100, width = 120, units = "mm", res = 300)
ggplot(data=model_comp_results, aes(x=comparison, y=percentage, fill=change)) +
  geom_bar(stat="identity", width=0.4) +
  labs(y="Percentage of species", x="", fill="") +
  scale_fill_manual(values = c("#339900", "#999999", "#990000")) +
  theme_bw() +
  theme(axis.text.x=element_text(colour = "black")) +
  theme(axis.text.y=element_text(colour = "black")) +
  scale_y_continuous(breaks = seq(0,100,10), expand = c(0, 0)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
dev.off()

## same plot but without overall change
## remove overall change
model_comp_results <- model_comp_results[!model_comp_results$comparison == "1985-2012",]
png("../Graphs/Model_comps/ukbms_percentage_trends2.png", height = 100, width = 120, units = "mm", res = 300)
ggplot(data=model_comp_results, aes(x=comparison, y=percentage, fill=change)) +
  geom_bar(stat="identity", width=0.4) +
  labs(y="Percentage of species", x="", fill="") +
  scale_fill_manual(values = c("#339900", "#999999", "#990000")) +
  theme_bw() +
  theme(axis.text.x=element_text(colour = "black")) +
  theme(axis.text.y=element_text(colour = "black")) +
  scale_y_continuous(breaks = seq(0,100,20), expand = c(0, 0)) +
  theme(text = element_text(size = 16)) +
  theme(legend.margin=margin(c(-15,40,-5,-20)), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
dev.off()
