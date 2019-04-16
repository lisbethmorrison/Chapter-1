##########################################
##### LIST OF MODELS WHICH NEED BOOTSTRAPPED:

#### UKBMS
## Change + specialism
## Average + mobility
## Change + mobility
## Change + abundance 00-12
## Change + STI

#### BBS
## Change + specialism
## Change + mobility
## Change + STI

rm(list=ls()) # clear R

## install packages
library(foreach)
library(doParallel)
library(lmerTest)
library(lme4)

###########################################################
####################### SPECIALISM ######################## 
###########################################################

## read in data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr_no_zeros2.csv", header=TRUE) 
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE) 

##### UKBMS (change over time)
#### subset only years 1985, 2000 and 2012
pair_attr_1985 <- pair_attr[pair_attr$mid.year==1984.5,]
pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]
pair_attr_ukbms <- rbind(pair_attr_1985, pair_attr_2000, pair_attr_2012) ## 3 years and 33 species
pair_attr_ukbms <- droplevels(pair_attr_ukbms)

pair_attr_ukbms$mid.year <- as.factor(pair_attr_ukbms$mid.year)
pair_attr_ukbms$pair.id <- as.character(pair_attr_ukbms$pair.id)
pair_attr_ukbms$spp <- as.factor(pair_attr_ukbms$spp)
pair_attr_ukbms$specialism <- as.factor(pair_attr_ukbms$specialism)

## run true model
spec_model_ukbms2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*specialism + (1|pair.id) + (1|spp), data = pair_attr_ukbms)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(spec_model_ukbms2)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## remove rows with mean northing, distance, renk hab sim, specialism and mid year F values (only interested in specialism*midyear interaction)
main_result_table <- main_result_table[ !(main_result_table$parameter %in% c("mean_northing", "distance", "renk_hab_sim", "mid.year", "specialism")), ]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_spec_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_ukbms$spec_shuffle <- sample(pair_attr_ukbms$specialism) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*spec_shuffle + (1|pair.id) + (1|spp), data = pair_attr_ukbms)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 9.618 minutes for 999 runs!

### save results
ukbms_spec_para$parameter <- paste(row.names(ukbms_spec_para)) ## move row.names to parameter column
rownames(ukbms_spec_para) <- 1:nrow(ukbms_spec_para) ## change row names to numbers
ukbms_spec_para <- ukbms_spec_para[,-c(1:3)] ## remove unnecessary columns
## only keep rows with interaction
ukbms_spec_para <- ukbms_spec_para[grep("mid.year:spec_shuffle", ukbms_spec_para$parameter),]
final_results_table <- rbind(main_result_table, ukbms_spec_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_change_spec_ukbms.csv", row.names=TRUE)
## read in file
perm_change_spec_ukbms <- read.csv("../Results/Model_outputs/UKBMS/perm_change_spec_ukbms.csv", header=TRUE) 
  
## Calculate p value
number_of_permutations <- 1000
true_change_spec_ukbms <- perm_change_spec_ukbms[perm_change_spec_ukbms$i==0,]
perm_change_spec_ukbms <- perm_change_spec_ukbms[!perm_change_spec_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_spec_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_spec_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant

####################################################
####################################################
#### BBS (change over time)
pair_attr_bbs_1999 <- pair_attr_BBS[pair_attr_BBS$mid.year==1998.5,]
pair_attr_bbs_2012 <- pair_attr_BBS[pair_attr_BBS$mid.year==2011.5,]
pair_attr_bbs <- rbind(pair_attr_bbs_1999, pair_attr_bbs_2012)

pair_attr_bbs$mid.year <- as.factor(pair_attr_bbs$mid.year)
pair_attr_bbs$pair.id <- as.character(pair_attr_bbs$pair.id)
pair_attr_bbs$spp <- as.factor(pair_attr_bbs$spp)
pair_attr_bbs$specialism <- as.factor(pair_attr_bbs$specialism)

## run true model
spec_model_bbs2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*specialism + (1|pair.id) + (1|spp), data = pair_attr_bbs)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(spec_model_bbs2)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## remove rows with mean northing, distance, renk hab sim, specialism and mid year F values (only interested in specialism*midyear interaction)
main_result_table <- main_result_table[ !(main_result_table$parameter %in% c("mean_northing", "distance", "renk_hab_sim", "mid.year", "specialism")), ]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
bbs_spec_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_bbs$spec_shuffle <- sample(pair_attr_bbs$specialism) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*spec_shuffle + (1|pair.id) + (1|spp), data = pair_attr_bbs)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 9.618 minutes for 999 runs!

### save results
bbs_spec_para$parameter <- paste(row.names(bbs_spec_para)) ## move row.names to parameter column
rownames(bbs_spec_para) <- 1:nrow(bbs_spec_para) ## change row names to numbers
bbs_spec_para <- bbs_spec_para[,-c(1:3)] ## remove unnecessary columns
## only interested in specialism interaction
bbs_spec_para <- bbs_spec_para[grep("mid.year:spec_shuffle", bbs_spec_para$parameter),]

final_results_table <- rbind(main_result_table, bbs_spec_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/BBS/perm_change_spec_bbs.csv", row.names=TRUE)
## read in file
perm_change_spec_bbs <- read.csv("../Results/Model_outputs/BBS/perm_change_spec_bbs.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_spec_bbs <- perm_change_spec_bbs[perm_change_spec_bbs$i==0,]
perm_change_spec_bbs <- perm_change_spec_bbs[!perm_change_spec_bbs$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_spec_bbs$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_spec_bbs$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0.056 non-significant

###########################################################
######################## MOBILITY ######################### 
###########################################################

## read in data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr_no_zeros2.csv", header=TRUE) 
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE) 
bird_dispersal <- read.csv("../Data/Woodland_bird_dispersal_Paradis1998.csv", header=TRUE)

## UKBMS (Average synchrony)
pair_attr <- na.omit(pair_attr)
length(unique(pair_attr$spp)) # 31 species

pair_attr$mid.year <- as.factor(pair_attr$mid.year)
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)

## run model with 2 groups (high and low) of mobility 
pair_attr$mobility_wil <- as.numeric(pair_attr$mobility_wil)
pair_attr$mobility_score2 <- cut(pair_attr$mobility_wil, 2, labels=c("low", "high"))
pair_attr$mobility_score2 <- as.factor(pair_attr$mobility_score2)

## run true model
mobility_model4 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(mobility_model4)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## remove rows with mean northing, distance, renk hab sim and mid year F values (only interested in mobility main effect)
main_result_table <- main_result_table[ !(main_result_table$parameter %in% c("mean_northing", "distance", "renk_hab_sim", "mid.year")), ]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_mob_para1 <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr$mob_shuffle <- sample(pair_attr$mobility_score2) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + mob_shuffle + (1|pair.id) + (1|spp), data = pair_attr)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 9.618 minutes for 999 runs!

### save results
ukbms_mob_para1$parameter <- paste(row.names(ukbms_mob_para1)) ## move row.names to parameter column
rownames(ukbms_mob_para1) <- 1:nrow(ukbms_mob_para1) ## change row names to numbers
ukbms_mob_para1 <- ukbms_mob_para1[,-c(1:3)] ## remove unnecessary columns
## only interested in mobility main effect
ukbms_mob_para1 <- ukbms_mob_para1[grep("mob_shuffle", ukbms_mob_para1$parameter),]
final_results_table <- rbind(main_result_table, ukbms_mob_para1) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_average_mob_ukbms.csv", row.names=TRUE)
## read in file
perm_average_mob_ukbms <- read.csv("../Results/Model_outputs/UKBMS/perm_average_mob_ukbms.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_average_mob_ukbms <- perm_average_mob_ukbms[perm_average_mob_ukbms$i==0,]
perm_average_mob_ukbms <- perm_average_mob_ukbms[!perm_average_mob_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_average_mob_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_average_mob_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0.019 significant

####################################################
####################################################
## UKBMS (Change in synchrony)
pair_attr_1985 <- pair_attr[pair_attr$mid.year==1984.5,]
pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]
pair_attr_ukbms <- rbind(pair_attr_1985, pair_attr_2000, pair_attr_2012) ## 3 years and 32 species
## remove NA's to make sure small white is taken out
pair_attr_ukbms <- na.omit(pair_attr_ukbms)
length(unique(pair_attr_ukbms$spp)) # 31 species

pair_attr_ukbms$mid.year <- as.factor(pair_attr_ukbms$mid.year)
pair_attr_ukbms$pair.id <- as.character(pair_attr_ukbms$pair.id)
pair_attr_ukbms$spp <- as.factor(pair_attr_ukbms$spp)
## run model with 2 mobility groups (low and high)
pair_attr_ukbms$mobility_wil <- as.numeric(pair_attr_ukbms$mobility_wil)
pair_attr_ukbms$mobility_score2 <- cut(pair_attr_ukbms$mobility_wil, 2, labels=FALSE)
pair_attr_ukbms$mobility_score2 <- as.factor(pair_attr_ukbms$mobility_score2)

## run true model
mobility_model7 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(mobility_model7)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## remove rows with mean northing, distance, renk hab sim, mobility and mid year F values (only interested in mobility*midyear interaction)
main_result_table <- main_result_table[ !(main_result_table$parameter %in% c("mean_northing", "distance", "renk_hab_sim", "mid.year", "mobility_score2")), ]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_mob_para2 <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_ukbms$mob_shuffle <- sample(pair_attr_ukbms$mobility_score2) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mob_shuffle + (1|pair.id) + (1|spp), data = pair_attr_ukbms)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 9.618 minutes for 999 runs!

### save results
ukbms_mob_para2$parameter <- paste(row.names(ukbms_mob_para2)) ## move row.names to parameter column
rownames(ukbms_mob_para2) <- 1:nrow(ukbms_mob_para2) ## change row names to numbers
ukbms_mob_para2 <- ukbms_mob_para2[,-c(1:3)] ## remove unnecessary columns
## only interested in mobility*midyear interaction
ukbms_mob_para2 <- ukbms_mob_para2[grep("mid.year:mob_shuffle", ukbms_mob_para2$parameter),]

final_results_table <- rbind(main_result_table, ukbms_mob_para2) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_change_mob_ukbms.csv", row.names=TRUE)
## read in file
perm_change_mob_ukbms <- read.csv("../Results/Model_outputs/UKBMS/perm_change_mob_ukbms.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_mob_ukbms <- perm_change_mob_ukbms[perm_change_mob_ukbms$i==0,]
perm_change_mob_ukbms <- perm_change_mob_ukbms[!perm_change_mob_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_mob_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_mob_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant

####################################################
####################################################
#### BBS (Change in synchrony)
pair_attr_BBS <- merge(pair_attr_BBS, bird_dispersal, by.x="spp", by.y="Species_code", all=FALSE)
length(unique(pair_attr_BBS$spp)) # 18 species until remove NAs for Breeding dispersal (==17)

pair_attr_bbs_1999 <- pair_attr_BBS[pair_attr_BBS$mid.year==1998.5,]
pair_attr_bbs_2012 <- pair_attr_BBS[pair_attr_BBS$mid.year==2011.5,]
pair_attr_bbs <- rbind(pair_attr_bbs_1999, pair_attr_bbs_2012)
pair_attr_bbs <- na.omit(pair_attr_bbs) # remove NA's (species without data)

pair_attr_bbs$mid.year <- as.factor(pair_attr_bbs$mid.year)
pair_attr_bbs$pair.id <- as.character(pair_attr_bbs$pair.id)
pair_attr_bbs$spp <- as.factor(pair_attr_bbs$spp)

## run model with 2 mobility groups (low and high)
pair_attr_bbs$Breeding_AM <- as.numeric(pair_attr_bbs$Breeding_AM)
pair_attr_bbs$Breeding_AM_score2 <- cut(pair_attr_bbs$Breeding_AM, 2, labels=FALSE)
pair_attr_bbs$Breeding_AM_score2 <- as.factor(pair_attr_bbs$Breeding_AM_score2)

## run true model
dispersal_model_bbs3 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*Breeding_AM_score2 + (1|spp) + (1|pair.id), data=pair_attr_bbs)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(dispersal_model_bbs3)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## remove rows with mean northing, distance, renk hab sim, mobility and mid year F values (only interested in mobility*midyear interaction)
main_result_table <- main_result_table[ !(main_result_table$parameter %in% c("mean_northing", "distance", "renk_hab_sim", "mid.year", "Breeding_AM_score2")), ]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
bbs_mob_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_bbs$mob_shuffle <- sample(pair_attr_bbs$Breeding_AM_score2) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mob_shuffle + (1|pair.id) + (1|spp), data = pair_attr_bbs)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 9.618 minutes for 999 runs!

### save results
bbs_mob_para$parameter <- paste(row.names(bbs_mob_para)) ## move row.names to parameter column
rownames(bbs_mob_para) <- 1:nrow(bbs_mob_para) ## change row names to numbers
bbs_mob_para <- bbs_mob_para[,-c(1:3)] ## remove unnecessary columns
## only interested in midyear*mobility interaction
bbs_mob_para <- bbs_mob_para[grep("mid.year:mob_shuffle", bbs_mob_para$parameter),]

final_results_table <- rbind(main_result_table, bbs_mob_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/BBS/perm_change_mob_bbs.csv", row.names=TRUE)
## read in file
perm_change_mob_bbs <- read.csv("../Results/Model_outputs/BBS/perm_change_mob_bbs.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_mob_bbs <- perm_change_mob_bbs[perm_change_mob_bbs$i==0,]
perm_change_mob_bbs <- perm_change_mob_bbs[!perm_change_mob_bbs$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_mob_bbs$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_mob_bbs$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant

##############################################################
######################### ABUNDANCE ########################## 
##############################################################

#### UKBMS (Change in synchrony 00-12)
## read in data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr_no_zeros2.csv", header=TRUE) 
abundance_results <- read.csv("../Results/Butterfly_results/abundance_results_00_12.csv", header=TRUE)

pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]
pair_attr_ukbms <- rbind(pair_attr_2000, pair_attr_2012) ## 2 years

pair_attr_ukbms$mid.year <- as.factor(pair_attr_ukbms$mid.year)
pair_attr_ukbms$pair.id <- as.character(pair_attr_ukbms$pair.id)
pair_attr_ukbms$spp <- as.factor(pair_attr_ukbms$spp)

## rescale variables 
pair_attr_ukbms$distance <- (pair_attr_ukbms$distance - mean(na.omit(pair_attr_ukbms$distance)))/sd(na.omit(pair_attr_ukbms$distance))
pair_attr_ukbms$mean_northing <- (pair_attr_ukbms$mean_northing - mean(na.omit(pair_attr_ukbms$mean_northing)))/sd(na.omit(pair_attr_ukbms$mean_northing))
pair_attr_ukbms$renk_hab_sim <- (pair_attr_ukbms$renk_hab_sim - mean(na.omit(pair_attr_ukbms$renk_hab_sim)))/sd(na.omit(pair_attr_ukbms$renk_hab_sim))

## merge with abundance data 
pair_attr_ukbms <- merge(pair_attr_ukbms, abundance_results, by.x="spp", by.y="Species.code", all=FALSE)

## true model
abund_model2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*ab_change_00_12 + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(abund_model2)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## remove rows with mean northing, distance, renk hab sim, mobility and mid year F values (only interested in abundance*midyear interaction)
main_result_table <- main_result_table[ !(main_result_table$parameter %in% c("mean_northing", "distance", "renk_hab_sim", "mid.year", "ab_change_00_12")), ]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_abund_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_ukbms$abund_shuffle <- sample(pair_attr_ukbms$ab_change_00_12) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*abund_shuffle + (1|pair.id) + (1|spp), data = pair_attr_ukbms)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 9.618 minutes for 999 runs!

### save results
ukbms_abund_para$parameter <- paste(row.names(ukbms_abund_para)) ## move row.names to parameter column
rownames(ukbms_abund_para) <- 1:nrow(ukbms_abund_para) ## change row names to numbers
ukbms_abund_para <- ukbms_abund_para[,-c(1:3)] ## remove unnecessary columns
## only interested in midyear*mobility interaction
ukbms_abund_para <- ukbms_abund_para[grep("mid.year:abund_shuffle", ukbms_abund_para$parameter),]

final_results_table <- rbind(main_result_table, ukbms_abund_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_change_abund_ukbms.csv", row.names=TRUE)
## read in file
perm_change_abund_ukbms <- read.csv("../Results/Model_outputs/UKBMS/perm_change_abund_ukbms.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_abund_ukbms <- perm_change_abund_ukbms[perm_change_abund_ukbms$i==0,]
perm_change_abund_ukbms <- perm_change_abund_ukbms[!perm_change_abund_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_abund_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_abund_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant

##############################################################
############################ STI ############################# 
##############################################################

## read in data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr_no_zeros2.csv", header=TRUE) 
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE) 
bird_STI <- read.csv("../Data/birds_STI.csv", header=TRUE)
butterfly_STI <- read.csv("../Data/butterflies_STI.csv", header=TRUE)

#### UKBMS
pair_attr <- merge(pair_attr, butterfly_STI, by.x="spp", by.y="species_code") ## 119 has no STI data
length(unique(pair_attr$spp)) # 31 species

pair_attr_1985 <- pair_attr[pair_attr$mid.year==1984.5,]
pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]
pair_attr_ukbms <- rbind(pair_attr_1985, pair_attr_2000, pair_attr_2012)
length(unique(pair_attr_ukbms$spp)) ## 3 years and 31 species

pair_attr_ukbms$mid.year <- as.factor(pair_attr_ukbms$mid.year)
pair_attr_ukbms$pair.id <- as.character(pair_attr_ukbms$pair.id)
pair_attr_ukbms$spp <- as.factor(pair_attr_ukbms$spp)

## true model
climate_ukbms_change <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(climate_ukbms_change)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## remove rows with mean northing, distance, renk hab sim, STI and mid year F values (only interested in STI*midyear interaction)
main_result_table <- main_result_table[ !(main_result_table$parameter %in% c("mean_northing", "distance", "renk_hab_sim", "mid.year", "STI")), ]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_sti_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_ukbms$STI_shuffle <- sample(pair_attr_ukbms$STI) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI_shuffle + (1|pair.id) + (1|spp), data = pair_attr_ukbms)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 9.618 minutes for 999 runs!

### save results
ukbms_sti_para$parameter <- paste(row.names(ukbms_sti_para)) ## move row.names to parameter column
rownames(ukbms_sti_para) <- 1:nrow(ukbms_sti_para) ## change row names to numbers
ukbms_sti_para <- ukbms_sti_para[,-c(1:3)] ## remove unnecessary columns
## only interested in midyear*mobility interaction
ukbms_sti_para <- ukbms_sti_para[grep("mid.year:STI_shuffle", ukbms_sti_para$parameter),]

final_results_table <- rbind(main_result_table, ukbms_sti_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_change_STI_ukbms.csv", row.names=TRUE)
## read in file
perm_change_STI_ukbms <- read.csv("../Results/Model_outputs/UKBMS/perm_change_STI_ukbms.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_STI_ukbms <- perm_change_STI_ukbms[perm_change_STI_ukbms$i==0,]
perm_change_STI_ukbms <- perm_change_STI_ukbms[!perm_change_STI_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_STI_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_STI_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant

#### BBS
pair_attr_BBS <- merge(pair_attr_BBS, bird_STI, by.x="spp", by.y="species_code") ## 119 has no STI data
length(unique(pair_attr_BBS$spp)) # 24 species

pair_attr_bbs_1999 <- pair_attr_BBS[pair_attr_BBS$mid.year==1998.5,]
pair_attr_bbs_2012 <- pair_attr_BBS[pair_attr_BBS$mid.year==2011.5,]
pair_attr_bbs <- rbind(pair_attr_bbs_1999, pair_attr_bbs_2012)

pair_attr_bbs$mid.year <- as.factor(pair_attr_bbs$mid.year)
pair_attr_bbs$pair.id <- as.character(pair_attr_bbs$pair.id)
pair_attr_bbs$spp <- as.factor(pair_attr_bbs$spp)
length(unique(pair_attr_bbs$spp)) # 24 species

## true model
climate_bbs_change <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI + (1|pair.id) + (1|spp), data = pair_attr_bbs)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(climate_bbs_change)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## remove rows with mean northing, distance, renk hab sim, STI and mid year F values (only interested in STI*midyear interaction)
main_result_table <- main_result_table[ !(main_result_table$parameter %in% c("mean_northing", "distance", "renk_hab_sim", "mid.year", "STI")), ]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
bbs_sti_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_bbs$STI_shuffle <- sample(pair_attr_bbs$STI) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI_shuffle + (1|pair.id) + (1|spp), data = pair_attr_bbs)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 9.618 minutes for 999 runs!

### save results
bbs_sti_para$parameter <- paste(row.names(bbs_sti_para)) ## move row.names to parameter column
rownames(bbs_sti_para) <- 1:nrow(bbs_sti_para) ## change row names to numbers
bbs_sti_para <- bbs_sti_para[,-c(1:3)] ## remove unnecessary columns
## only interested in midyear*mobility interaction
bbs_sti_para <- bbs_sti_para[grep("mid.year:STI_shuffle", bbs_sti_para$parameter),]

final_results_table <- rbind(main_result_table, bbs_sti_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/BBS/perm_change_STI_bbs.csv", row.names=TRUE)
## read in file
perm_change_STI_bbs <- read.csv("../Results/Model_outputs/BBS/perm_change_STI_bbs.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_STI_bbs <- perm_change_STI_bbs[perm_change_STI_bbs$i==0,]
perm_change_STI_bbs <- perm_change_STI_bbs[!perm_change_STI_bbs$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_STI_bbs$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_STI_bbs$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant

