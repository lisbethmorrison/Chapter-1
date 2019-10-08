#############################################################################
## Title: Permutation tests for significant results
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: January 2019
#############################################################################

##########################################
##### LIST OF MODELS WHICH NEED BOOTSTRAPPED:

#### UKBMS
## Change + specialism DONE
## Average + mobility DONE
## Change + mobility 00-12 DONE
## Change + abundance 00-12 DONE
## mean northing DONE
## distance DONE
## habitat similarity DONE

#### BBS
## Change + specialism DONE
## Change + mobility DONE
## habitat similarity DONE

#### CBC
## Change + specialism DONE
## average + abundance DONE

rm(list=ls()) # clear R

## install packages
library(foreach)
library(doParallel)
library(lmerTest)
library(lme4)
#library(doSNOW)

###########################################################
###################### FIXED EFFECTS ###################### 
###########################################################

## read in data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pop_climate_synchrony.csv", header=TRUE) 
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pop_climate_synchrony_BBS.csv", header=TRUE) 

##### UKBMS (northing, distance and habitat similarity)
pair_attr$mid.year <- as.factor(pair_attr$mid.year)
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)

############## NORTHING #################
## run true model
fixed_model1 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain +
                       winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr)

true_result_table <- data.frame(anova(fixed_model1)[5]) ## save anova table from main model
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## keep rows with mean northing
true_result_table <- true_result_table[grep("mean_northing", true_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_perm_north <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr$northing_shuffle <- sample(pair_attr$mean_northing) ## randomly shuffle northing varaible
  model <- lmer(lag0 ~ northing_shuffle + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain +
                  winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr)
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
ukbms_perm_north$parameter <- paste(row.names(ukbms_perm_north)) ## move row.names to parameter column
rownames(ukbms_perm_north) <- 1:nrow(ukbms_perm_north) ## change row names to numbers
ukbms_perm_north <- ukbms_perm_north[,-c(1:3)] ## remove unnecessary columns
## only interested in northing shuffle main effect
ukbms_perm_north <- ukbms_perm_north[grep("northing_shuffle", ukbms_perm_north$parameter),]
final_results_table <- rbind(true_result_table, ukbms_perm_north) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_mean_north_ukbms.csv", row.names=TRUE)
## read in file
perm_mean_north_ukbms <- read.csv("../Results/Model_outputs/UKBMS/perm_mean_north_ukbms.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_mean_north_ukbms <- perm_mean_north_ukbms[perm_mean_north_ukbms$i==0,]
perm_mean_north_ukbms <- perm_mean_north_ukbms[!perm_mean_north_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_mean_north_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_mean_north_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant

############## DISTANCE #################
## run true model
fixed_model1 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain +
                       winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr)

true_result_table <- data.frame(anova(fixed_model1)[5]) ## save anova table from main model
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## keep rows with distance
true_result_table <- true_result_table[grep("distance", true_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores <- parallel::detectCores()
cl <- makeCluster(cores[1]-1) # not to overload your computer
registerDoParallel(cl)

n_sims <- 999

## this code creates a progress percentage bar
pb <- txtProgressBar(max = n_sims, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress) 

para_start_time = Sys.time()
ukbms_perm_distance <- foreach (i=1:n_sims, .options.snow = opts, .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr$distance_shuffle <- sample(pair_attr$distance) ## randomly shuffle northing varaible
  model <- lmer(lag0 ~ distance_shuffle + mean_northing + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain +
                  winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
close(pb)
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 21 hours

### save results
ukbms_perm_distance$parameter <- paste(row.names(ukbms_perm_distance)) ## move row.names to parameter column
rownames(ukbms_perm_distance) <- 1:nrow(ukbms_perm_distance) ## change row names to numbers
ukbms_perm_distance <- ukbms_perm_distance[,-c(1:3)] ## remove unnecessary columns
## only interested in distance shuffle main effect
ukbms_perm_distance <- ukbms_perm_distance[grep("distance_shuffle", ukbms_perm_distance$parameter),]
final_results_table <- rbind(true_result_table, ukbms_perm_distance) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_distance_ukbms.csv", row.names=FALSE)
## read in file
perm_distance_ukbms <- read.csv("../Results/Model_outputs/UKBMS/perm_distance_ukbms.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_distance_ukbms <- perm_distance_ukbms[perm_distance_ukbms$i==0,]
perm_distance_ukbms <- perm_distance_ukbms[!perm_distance_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_distance_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_distance_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant

############## HABITAT SIMILARITY #################
## run true model
## don't need to run the model again (use the same one as northing, but extract distance)
true_result_table <- data.frame(anova(fixed_model1)[5]) ## save anova table from main model
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## keep rows with renk_hab_sim
true_result_table <- true_result_table[grep("renk_hab_sim", true_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores <- parallel::detectCores()
cl <- makeCluster(cores[1]-1) # not to overload your computer
registerDoParallel(cl)

n_sims <- 999

## this code creates a progress percentage bar
pb <- txtProgressBar(max = n_sims, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress) 

para_start_time = Sys.time()
ukbms_perm_hab <- foreach (i=1:n_sims, .options.snow = opts, .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr$hab_shuffle <- sample(pair_attr$renk_hab_sim) ## randomly shuffle northing varaible
  model <- lmer(lag0 ~ hab_shuffle + mean_northing + distance + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain +
                  winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
close(pb)
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 14 hours

### save results
ukbms_perm_hab$parameter <- paste(row.names(ukbms_perm_hab)) ## move row.names to parameter column
rownames(ukbms_perm_hab) <- 1:nrow(ukbms_perm_hab) ## change row names to numbers
ukbms_perm_hab <- ukbms_perm_hab[,-c(1:3)] ## remove unnecessary columns
## only interested in mobility main effect
ukbms_perm_hab <- ukbms_perm_hab[grep("hab_shuffle", ukbms_perm_hab$parameter),]
final_results_table <- rbind(true_result_table, ukbms_perm_hab) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_hab_ukbms.csv", row.names=FALSE)
## read in file
perm_hab_ukbms <- read.csv("../Results/Model_outputs/UKBMS/perm_hab_ukbms.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_hab_ukbms <- perm_hab_ukbms[perm_hab_ukbms$i==0,]
perm_hab_ukbms <- perm_hab_ukbms[!perm_hab_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_hab_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_hab_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant

#############################################################
#############################################################
##### BBS (habitat similarity)
pair_attr_BBS$mid.year <- as.factor(pair_attr_BBS$mid.year)
pair_attr_BBS$pair.id <- as.character(pair_attr_BBS$pair.id)
pair_attr_BBS$spp <- as.factor(pair_attr_BBS$spp)

############## HABITAT SIMILARITY #################
## run true model
fixed_model_BBS <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + 
                       (1|pair.id) + (1|spp), data = pair_attr_BBS)

true_result_table <- data.frame(anova(fixed_model_BBS)[5]) ## save anova table from main model
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## keep rows with renk_hab_sim
true_result_table <- true_result_table[grep("renk_hab_sim", true_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores <- parallel::detectCores()
cl <- makeCluster(cores[1]-1) # not to overload your computer
registerDoParallel(cl)

n_sims <- 999

## this code creates a progress percentage bar
pb <- txtProgressBar(max = n_sims, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

para_start_time = Sys.time()
bbs_perm_hab <- foreach (i=1:n_sims, .options.snow = opts, .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_BBS$hab_shuffle <- sample(pair_attr_BBS$renk_hab_sim) ## randomly shuffle habitat varaible
  model <- lmer(lag0 ~ hab_shuffle + mean_northing + distance + mid.year + winter_rain + autumn_rain + spring_rain + 
                  (1|pair.id) + (1|spp), data = pair_attr_BBS)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
close(pb)
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 13 hours

### save results
bbs_perm_hab$parameter <- paste(row.names(bbs_perm_hab)) ## move row.names to parameter column
rownames(bbs_perm_hab) <- 1:nrow(bbs_perm_hab) ## change row names to numbers
bbs_perm_hab <- bbs_perm_hab[,-c(1:3)] ## remove unnecessary columns
## only interested in mobility main effect
bbs_perm_hab <- bbs_perm_hab[grep("hab_shuffle", bbs_perm_hab$parameter),]
final_results_table <- rbind(true_result_table, bbs_perm_hab) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/BBS/perm_hab_bbs.csv", row.names=FALSE)
## read in file
perm_hab_bbs <- read.csv("../Results/Model_outputs/BBS/perm_hab_bbs.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_hab_bbs <- perm_hab_bbs[perm_hab_bbs$i==0,]
perm_hab_bbs <- perm_hab_bbs[!perm_hab_bbs$i==0,] ## remove true value to calc. p value
diff.observed <- true_hab_bbs$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_hab_bbs$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant 

###########################################################
####################### SPECIALISM ######################## 
###########################################################

## read in data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pop_climate_synchrony.csv", header=TRUE) # butterfly pair attribute data
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pop_climate_synchrony_CBC.csv", header=TRUE) # CBC pair attribute data
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pop_climate_synchrony_BBS.csv", header=TRUE) # BBS pair attribute data

##### UKBMS (change over time)
#### subset only years 1985, 2000 and 2012
pair_attr_1985 <- pair_attr[pair_attr$mid.year==1984.5,]
pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]
pair_attr_early <- rbind(pair_attr_1985, pair_attr_2000) ## 2 years and 33 species
pair_attr_late <- rbind(pair_attr_2000, pair_attr_2012) ## 2 years and 33 species
pair_attr_early <- droplevels(pair_attr_early)
pair_attr_late <- droplevels(pair_attr_late)

pair_attr_early$mid.year <- as.factor(pair_attr_early$mid.year)
pair_attr_early$pair.id <- as.character(pair_attr_early$pair.id)
pair_attr_early$spp <- as.factor(pair_attr_early$spp)
pair_attr_early$specialism <- as.factor(pair_attr_early$specialism)

pair_attr_late$mid.year <- as.factor(pair_attr_late$mid.year)
pair_attr_late$pair.id <- as.character(pair_attr_late$pair.id)
pair_attr_late$spp <- as.factor(pair_attr_late$spp)
pair_attr_late$specialism <- as.factor(pair_attr_late$specialism)

##### EARLY ######
## run true model
spec_model_ukbms2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*specialism + winter_rain + autumn_rain + spring_rain + summer_rain +
                            winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr_early)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(spec_model_ukbms2)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with specialism*midyear interaction
main_result_table <- main_result_table[grep("mid.year:specialism", main_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_spec_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_early$spec_shuffle <- sample(pair_attr_early$specialism) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*spec_shuffle + winter_rain + autumn_rain + spring_rain + summer_rain +
                  winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr_early)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 31 minutes for 999 runs

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
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_change_spec_ukbms_85_00.csv", row.names=FALSE)
## read in file
perm_change_spec_ukbms <- read.csv("../Results/Model_outputs/UKBMS/perm_change_spec_ukbms_85_00.csv", header=TRUE) 
  
## Calculate p value
number_of_permutations <- 1000
true_change_spec_ukbms <- perm_change_spec_ukbms[perm_change_spec_ukbms$i==0,]
perm_change_spec_ukbms <- perm_change_spec_ukbms[!perm_change_spec_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_spec_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_spec_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant

##### LATE ######
## run true model
spec_model_ukbms2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*specialism + winter_rain + autumn_rain + spring_rain + summer_rain +
                            winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr_late)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(spec_model_ukbms2)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with specialism*midyear interaction
main_result_table <- main_result_table[grep("mid.year:specialism", main_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_spec_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_late$spec_shuffle <- sample(pair_attr_late$specialism) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*spec_shuffle + winter_rain + autumn_rain + spring_rain + summer_rain +
                  winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr_late)
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
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_change_spec_ukbms_00_12.csv", row.names=FALSE)
## read in file
perm_change_spec_ukbms <- read.csv("../Results/Model_outputs/UKBMS/perm_change_spec_ukbms_00_12.csv", header=TRUE) 

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
spec_model_bbs2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*specialism + winter_rain + autumn_rain + spring_rain + 
                          (1|pair.id) + (1|spp), data = pair_attr_bbs)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(spec_model_bbs2)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with specialism*midyear interaction
main_result_table <- main_result_table[grep("mid.year:specialism", main_result_table$parameter),]

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
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*spec_shuffle + winter_rain + autumn_rain + spring_rain + 
                  (1|pair.id) + (1|spp), data = pair_attr_bbs)
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
pvalue ## 0.068 non-significant

####################################################
####################################################
#### CBC (change over time)
pair_attr_cbc_1985 <- pair_attr_CBC[pair_attr_CBC$mid.year==1984.5,]
pair_attr_cbc_1996 <- pair_attr_CBC[pair_attr_CBC$mid.year==1995.5,]
pair_attr_cbc <- rbind(pair_attr_cbc_1985, pair_attr_cbc_1996)

pair_attr_cbc$mid.year <- as.factor(pair_attr_cbc$mid.year)
pair_attr_cbc$pair.id <- as.character(pair_attr_cbc$pair.id)
pair_attr_cbc$spp <- as.factor(pair_attr_cbc$spp)
pair_attr_cbc$specialism <- as.factor(pair_attr_cbc$specialism)

## run true model
spec_model_cbc2 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*specialism + summer_temp + (1|family) + (1|pair.id) + (1|spp), data = pair_attr_cbc)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(spec_model_cbc2)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with specialism*midyear interaction
main_result_table <- main_result_table[grep("mid.year:specialism", main_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
cbc_spec_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_cbc$spec_shuffle <- sample(pair_attr_cbc$specialism) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*spec_shuffle + summer_temp + (1|family) + (1|pair.id) + (1|spp), data = pair_attr_cbc)
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
cbc_spec_para$parameter <- paste(row.names(cbc_spec_para)) ## move row.names to parameter column
rownames(cbc_spec_para) <- 1:nrow(cbc_spec_para) ## change row names to numbers
cbc_spec_para <- cbc_spec_para[,-c(1:3)] ## remove unnecessary columns
## only interested in specialism interaction
cbc_spec_para <- cbc_spec_para[grep("mid.year:spec_shuffle", cbc_spec_para$parameter),]

final_results_table <- rbind(main_result_table, cbc_spec_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/CBC/perm_change_spec_cbc.csv", row.names=TRUE)
## read in file
perm_change_spec_cbc <- read.csv("../Results/Model_outputs/CBC/perm_change_spec_cbc.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_spec_cbc <- perm_change_spec_cbc[perm_change_spec_cbc$i==0,]
perm_change_spec_cbc <- perm_change_spec_cbc[!perm_change_spec_cbc$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_spec_cbc$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_spec_cbc$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0


###########################################################
######################## MOBILITY ######################### 
###########################################################

## read in data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pop_climate_synchrony.csv", header=TRUE) # butterfly pair attribute data
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pop_climate_synchrony_BBS.csv", header=TRUE) # BBS pair attribute data
bird_dispersal <- read.csv("../Data/Woodland_bird_dispersal_Paradis1998.csv", header=TRUE)

## UKBMS (Average synchrony)
pair_attr <- na.omit(pair_attr)
length(unique(pair_attr$spp)) # 31 species

pair_attr$mid.year <- as.factor(pair_attr$mid.year)
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)

## run true model
mobility_model4 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + mobility_wil + winter_rain + autumn_rain + spring_rain + summer_rain + 
                          winter_temp + autumn_temp + spring_temp + summer_temp + (1|spp) + (1|pair.id), data=pair_attr)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(mobility_model4)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with mobility main effect
main_result_table <- main_result_table[grep("mobility_wil", main_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_mob_para1 <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr$mob_shuffle <- sample(pair_attr$mobility_wil) ## randomly shuffle mobility varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + mob_shuffle + winter_rain + autumn_rain + spring_rain + summer_rain + 
                  winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr)
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
pvalue ## 0.004 significant

####################################################
####################################################
## UKBMS (Change in synchrony)
pair_attr_1985 <- pair_attr[pair_attr$mid.year==1984.5,]
pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]
pair_attr_early <- rbind(pair_attr_1985, pair_attr_2000) ## 2 years and 32 species
pair_attr_late <- rbind(pair_attr_2000, pair_attr_2012) ## 2 years and 32 species
## remove NA's to make sure small white is taken out
pair_attr_early <- na.omit(pair_attr_early)
pair_attr_late <- na.omit(pair_attr_late)

pair_attr_early$mid.year <- as.factor(pair_attr_early$mid.year)
pair_attr_early$pair.id <- as.character(pair_attr_early$pair.id)
pair_attr_early$spp <- as.factor(pair_attr_early$spp)

pair_attr_late$mid.year <- as.factor(pair_attr_late$mid.year)
pair_attr_late$pair.id <- as.character(pair_attr_late$pair.id)
pair_attr_late$spp <- as.factor(pair_attr_late$spp)

###############################
# ###### EARLY #######
# ## run true model
# mobility_model7 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mobility_wil + winter_rain + autumn_rain + spring_rain + summer_rain + 
#                           winter_temp + autumn_temp + spring_temp + summer_temp + (1|spp) + (1|pair.id), data=pair_attr_early)
# ## save true model results (to merge in with bootstrapped models later)
# main_result_table <- data.frame(anova(mobility_model7)[5]) ## save anova table from main model
# main_result_table$i <- 0 ## make i column with zeros 
# main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
# rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
# ## keep rows with mobility*midyear interaction
# main_result_table <- main_result_table[grep("mid.year:mobility_wil", main_result_table$parameter),]
# 
# ## run 999 permutation tests
# ## Set up number of cores to run on (7)
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)
# 
# n_sims <- 999
# para_start_time = Sys.time()
# ukbms_mob_para2 <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
#   print(i)
#   pair_attr_early$mob_shuffle <- sample(pair_attr_early$mobility_wil) ## randomly shuffle mobility varaible
#   model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mob_shuffle + winter_rain + autumn_rain + spring_rain + summer_rain + 
#                   winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr_early)
#   ## run model with shuffled variable
#   ## save results
#   anoresult<-anova(model)
#   data.frame(anoresult, i=i)
# }
# stopCluster(cl)
# para_end_time = Sys.time()
# para_run_time = para_end_time - para_start_time
# print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 9.618 minutes for 999 runs!
# 
# ### save results
# ukbms_mob_para2$parameter <- paste(row.names(ukbms_mob_para2)) ## move row.names to parameter column
# rownames(ukbms_mob_para2) <- 1:nrow(ukbms_mob_para2) ## change row names to numbers
# ukbms_mob_para2 <- ukbms_mob_para2[,-c(1:3)] ## remove unnecessary columns
# ## only interested in mobility*midyear interaction
# ukbms_mob_para2 <- ukbms_mob_para2[grep("mid.year:mob_shuffle", ukbms_mob_para2$parameter),]
# 
# final_results_table <- rbind(main_result_table, ukbms_mob_para2) ## bind the two data frames together
# 
# F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
# hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)
# 
# ## save file
# write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_change_mob_ukbms_85_00.csv", row.names=TRUE)
# ## read in file
# perm_change_mob_ukbms <- read.csv("../Results/Model_outputs/UKBMS/perm_change_mob_ukbms_85_00.csv", header=TRUE) 
# 
# ## Calculate p value
# number_of_permutations <- 1000
# true_change_mob_ukbms <- perm_change_mob_ukbms[perm_change_mob_ukbms$i==0,]
# perm_change_mob_ukbms <- perm_change_mob_ukbms[!perm_change_mob_ukbms$i==0,] ## remove true value to calc. p value
# diff.observed <- true_change_mob_ukbms$F.value ## true F value
# 
# # P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
# pvalue = sum(abs(perm_change_mob_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
# pvalue ## 0 significant
###############################

###### LATE #######
## run true model
mobility_model7 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mobility_wil + winter_rain + autumn_rain + spring_rain + summer_rain + 
                          winter_temp + autumn_temp + spring_temp + summer_temp + (1|spp) + (1|pair.id), data=pair_attr_late)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(mobility_model7)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with mobility*midyear interaction
main_result_table <- main_result_table[grep("mid.year:mobility_wil", main_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
ukbms_mob_para2 <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_late$mob_shuffle <- sample(pair_attr_late$mobility_wil) ## randomly shuffle mobility varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mob_shuffle + winter_rain + autumn_rain + spring_rain + summer_rain + 
                  winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr_late)
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
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_change_mob_ukbms_00_12.csv", row.names=TRUE)
## read in file
perm_change_mob_ukbms <- read.csv("../Results/Model_outputs/UKBMS/perm_change_mob_ukbms_00_12.csv", header=TRUE) 

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

## run true model
dispersal_model_bbs3 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*Breeding_AM + winter_rain + autumn_rain + spring_rain + 
                               (1|spp) + (1|pair.id), data=pair_attr_bbs)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(dispersal_model_bbs3)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with mobility*midyear interaction
main_result_table <- main_result_table[grep("mid.year:Breeding_AM", main_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

n_sims <- 999
para_start_time = Sys.time()
bbs_mob_para <- foreach (i=1:n_sims,  .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_bbs$mob_shuffle <- sample(pair_attr_bbs$Breeding_AM) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mob_shuffle + winter_rain + autumn_rain + spring_rain + 
                  (1|pair.id) + (1|spp), data = pair_attr_bbs)
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
pair_attr <- read.csv("../Data/Butterfly_sync_data/pop_climate_synchrony.csv", header=TRUE) # butterfly pair attribute data
abundance_results <- read.csv("../Results/Butterfly_results/abundance_results_00_12.csv", header=TRUE)

pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]
pair_attr_ukbms <- rbind(pair_attr_2000, pair_attr_2012) ## 2 years

pair_attr_ukbms$mid.year <- as.factor(pair_attr_ukbms$mid.year)
pair_attr_ukbms$pair.id <- as.character(pair_attr_ukbms$pair.id)
pair_attr_ukbms$spp <- as.factor(pair_attr_ukbms$spp)

## merge with abundance data 
pair_attr_ukbms <- merge(pair_attr_ukbms, abundance_results, by.x="spp", by.y="Species.code", all=FALSE)

## true model
abund_model2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*ab_change_00_12 + winter_rain + autumn_rain + spring_rain + summer_rain + 
                       winter_temp + autumn_temp + spring_temp + summer_temp + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(abund_model2)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with abundance*midyear interaction
main_result_table <- main_result_table[grep("mid.year:ab_change_00_12", main_result_table$parameter),]

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
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*abund_shuffle + winter_rain + autumn_rain + spring_rain + summer_rain + 
                  winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr_ukbms)
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

##### Average abundance (CBC)
## read in data
bird_common <- read.csv("../Data/BTO_data/pop_estimates_birds.csv", header=TRUE)
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pop_climate_synchrony_CBC.csv", header=TRUE) # CBC pair attribute data

## merge the two datasets
pair_attr_CBC <- merge(pair_attr_CBC, bird_common, by.x="spp", by.y="species_code", all=FALSE)
## rescale variables
pair_attr_CBC$pop_estimate_stand <- (pair_attr_CBC$pop_estimate - mean(na.omit(pair_attr_CBC$pop_estimate)))/sd(na.omit(pair_attr_CBC$pop_estimate))

## true model
common_model_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + pop_estimate_stand + summer_temp + (1|pair.id) + (1|spp), data = pair_attr_CBC)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(common_model_cbc)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## keep rows with abundance
main_result_table <- main_result_table[grep("pop_estimate_stand", main_result_table$parameter),]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores <- parallel::detectCores()
cl <- makeCluster(cores[1]-1) # not to overload your computer
registerDoParallel(cl)

n_sims <- 999

pb <- txtProgressBar(max = n_sims, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress) 

para_start_time = Sys.time()
cbc_abund_para <- foreach (i=1:n_sims, .options.snow = opts, .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_CBC$abund_shuffle <- sample(pair_attr_CBC$pop_estimate_stand) ## randomly shuffle abundance varaible
  model <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + abund_shuffle + summer_temp + (1|pair.id) + (1|spp), data = pair_attr_CBC)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 35 minutes for 999 runs

### save results
cbc_abund_para$parameter <- paste(row.names(cbc_abund_para)) ## move row.names to parameter column
rownames(cbc_abund_para) <- 1:nrow(cbc_abund_para) ## change row names to numbers
cbc_abund_para <- cbc_abund_para[,-c(1:3)] ## remove unnecessary columns
## only interested in abund shuffle
cbc_abund_para <- cbc_abund_para[grep("abund_shuffle", cbc_abund_para$parameter),]

final_results_table <- rbind(main_result_table, cbc_abund_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/CBC/perm_average_abund_cbc.csv", row.names=TRUE)
## read in file
perm_average_abund_cbc <- read.csv("../Results/Model_outputs/CBC/perm_average_abund_cbc.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_average_abund_cbc <- perm_average_abund_cbc[perm_average_abund_cbc$i==0,]
perm_average_abund_cbc <- perm_average_abund_cbc[!perm_average_abund_cbc$i==0,] ## remove true value to calc. p value
diff.observed <- true_average_abund_cbc$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_average_abund_cbc$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0.008

