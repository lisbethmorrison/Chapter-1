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
library(doSNOW)

###########################################################
###################### FIXED EFFECTS ###################### 
###########################################################

## read in data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr_no_zeros2.csv", header=TRUE) 
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE) 
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC_no_zeros2.csv", header=TRUE)

##### UKBMS (northing, distance and habitat similarity)
pair_attr$mid.year <- as.factor(pair_attr$mid.year)
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)

############## NORTHING #################
## run true model
fixed_model1 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), data = pair_attr)

true_result_table <- data.frame(anova(fixed_model1)[5]) ## save anova table from main model
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## remove rows with distance, renk hab sim and mid year F values (only interested in mean northing)
true_result_table <- true_result_table[ !(true_result_table$parameter %in% c("distance", "renk_hab_sim", "mid.year")), ]

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
  model <- lmer(lag0 ~ northing_shuffle + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), data = pair_attr)
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
## only interested in mobility main effect
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
pvalue ## 0.019 significant

############## DISTANCE #################
## run true model
fixed_model1 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), data = pair_attr)

true_result_table <- data.frame(anova(fixed_model1)[5]) ## save anova table from main model
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## remove rows with mean_northing, renk hab sim and mid year F values (only interested in distance)
true_result_table <- true_result_table[ !(true_result_table$parameter %in% c("mean_northing", "renk_hab_sim", "mid.year")), ]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores <- parallel::detectCores()
cl <- makeSOCKcluster(cores[1]-1) # not to overload your computer
registerDoSNOW(cl)

n_sims <- 999

## this code creates a progress percentage bar
pb <- txtProgressBar(max = n_sims, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress) 

para_start_time = Sys.time()
ukbms_perm_distance <- foreach (i=1:n_sims, .options.snow = opts, .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr$distance_shuffle <- sample(pair_attr$distance) ## randomly shuffle northing varaible
  model <- lmer(lag0 ~ distance_shuffle + mean_northing + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), data = pair_attr)
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
ukbms_perm_distance$parameter <- paste(row.names(ukbms_perm_distance)) ## move row.names to parameter column
rownames(ukbms_perm_distance) <- 1:nrow(ukbms_perm_distance) ## change row names to numbers
ukbms_perm_distance <- ukbms_perm_distance[,-c(1:3)] ## remove unnecessary columns
## only interested in mobility main effect
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
pvalue ## 0.019 significant

############## HABITAT SIMILARITY #################
## run true model
## don't need to run the model again (use the same one as northing, but extract distance)
true_result_table <- data.frame(anova(fixed_model1)[5]) ## save anova table from main model
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## remove rows with mean_northing, distance and mid year F values (only interested in renk_hab_sim)
true_result_table <- true_result_table[ !(true_result_table$parameter %in% c("mean_northing", "distance", "mid.year")), ]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores <- parallel::detectCores()
cl <- makeSOCKcluster(cores[1]-1) # not to overload your computer
registerDoSNOW(cl)

n_sims <- 999

## this code creates a progress percentage bar
pb <- txtProgressBar(max = n_sims, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress) 

para_start_time = Sys.time()
ukbms_perm_hab <- foreach (i=1:n_sims, .options.snow = opts, .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr$hab_shuffle <- sample(pair_attr$renk_hab_sim) ## randomly shuffle northing varaible
  model <- lmer(lag0 ~ hab_shuffle + mean_northing + distance + mid.year + (1|pair.id) + (1|spp), data = pair_attr)
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
pvalue ## 0.019 significant

#############################################################
#############################################################
##### BBS (distance and habitat similarity)
pair_attr_BBS$mid.year <- as.factor(pair_attr_BBS$mid.year)
pair_attr_BBS$pair.id <- as.character(pair_attr_BBS$pair.id)
pair_attr_BBS$spp <- as.factor(pair_attr_BBS$spp)

############## DISTANCE #################
## run true model
fixed_model_BBS <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), data = pair_attr_BBS)

## don't need to run the model again (use the same one as northing, but extract distance)
true_result_table <- data.frame(anova(fixed_model1)[5]) ## save anova table from main model
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## remove rows with mean_northing, renk hab sim and mid year F values (only interested in distance)
true_result_table <- true_result_table[ !(true_result_table$parameter %in% c("mean_northing", "renk_hab_sim", "mid.year")), ]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores <- parallel::detectCores()
cl <- makeSOCKcluster(cores[1]-1) # not to overload your computer
registerDoSNOW(cl)

n_sims <- 999

## this code creates a progress percentage bar
pb <- txtProgressBar(max = n_sims, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress) 

para_start_time = Sys.time()
bbs_perm_distance <- foreach (i=1:n_sims, .options.snow = opts, .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_BBS$distance_shuffle <- sample(pair_attr_BBS$distance) ## randomly shuffle northing varaible
  model <- lmer(lag0 ~ distance_shuffle + mean_northing + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), data = pair_attr_BBS)
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
bbs_perm_distance$parameter <- paste(row.names(bbs_perm_distance)) ## move row.names to parameter column
rownames(bbs_perm_distance) <- 1:nrow(bbs_perm_distance) ## change row names to numbers
bbs_perm_distance <- bbs_perm_distance[,-c(1:3)] ## remove unnecessary columns
## only interested in mobility main effect
bbs_perm_distance <- bbs_perm_distance[grep("distance_shuffle", bbs_perm_distance$parameter),]
final_results_table <- rbind(true_result_table, bbs_perm_distance) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/BBS/perm_distance_bbs.csv", row.names=FALSE)
## read in file
perm_distance_bbs <- read.csv("../Results/Model_outputs/BBS/perm_distance_bbs.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_distance_bbs <- perm_distance_bbs[perm_distance_bbs$i==0,]
perm_distance_bbs <- perm_distance_bbs[!perm_distance_bbs$i==0,] ## remove true value to calc. p value
diff.observed <- true_distance_bbs$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_distance_bbs$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0.019 significant

############## HABITAT SIMILARITY #################
## run true model
## don't need to run the model again (use the same one as northing, but extract distance)
true_result_table <- data.frame(anova(fixed_model_BBS)[5]) ## save anova table from main model
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## remove rows with mean_northing, distance and mid year F values (only interested in renk_hab_sim)
true_result_table <- true_result_table[ !(true_result_table$parameter %in% c("mean_northing", "distance", "mid.year")), ]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores <- parallel::detectCores()
cl <- makeSOCKcluster(cores[1]-1) # not to overload your computer
registerDoSNOW(cl)

n_sims <- 999

## this code creates a progress percentage bar
pb <- txtProgressBar(max = n_sims, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

para_start_time = Sys.time()
bbs_perm_hab <- foreach (i=1:n_sims, .options.snow = opts, .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_BBS$hab_shuffle <- sample(pair_attr_BBS$renk_hab_sim) ## randomly shuffle northing varaible
  model <- lmer(lag0 ~ hab_shuffle + mean_northing + distance + mid.year + (1|pair.id) + (1|spp), data = pair_attr_BBS)
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
pvalue ## 0.019 significant

#############################################################
#############################################################
##### CBC (distance)
pair_attr_CBC$mid.year <- as.factor(pair_attr_CBC$mid.year)
pair_attr_CBC$pair.id <- as.character(pair_attr_CBC$pair.id)
pair_attr_CBC$spp <- as.factor(pair_attr_CBC$spp)

############## DISTANCE #################
## run true model
fixed_model_CBC <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + (1|pair.id) + (1|spp), data = pair_attr_CBC)

## don't need to run the model again (use the same one as northing, but extract distance)
true_result_table <- data.frame(anova(fixed_model_CBC)[5]) ## save anova table from main model
true_result_table$i <- 0 ## make i column with zeros 
true_result_table$parameter <- paste(row.names(true_result_table)) ## move row.names to parameter column
rownames(true_result_table) <- 1:nrow(true_result_table) ## change row names to numbers
## remove rows with mean_northing, renk hab sim and mid year F values (only interested in distance)
true_result_table <- true_result_table[!(true_result_table$parameter %in% c("mean_northing", "hab_sim", "mid.year")), ]

## run 999 permutation tests

n_sims <- 999
cores <- parallel::detectCores()
cl <- makeSOCKcluster(cores[1]-1) # not to overload your computer
registerDoSNOW(cl)

## this code creates a progress percentage bar
pb <- txtProgressBar(max = n_sims, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress) 

para_start_time = Sys.time()
cbc_perm_distance <- foreach (i=1:n_sims, .options.snow = opts, .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_CBC$distance_shuffle <- sample(pair_attr_CBC$distance) ## randomly shuffle northing varaible
  model <- lmer(lag0 ~ distance_shuffle + mean_northing + hab_sim + mid.year + (1|pair.id) + (1|spp), data = pair_attr_CBC)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
close(pb)
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 46 minutes for 999 runs!

### save results
cbc_perm_distance$parameter <- paste(row.names(cbc_perm_distance)) ## move row.names to parameter column
rownames(cbc_perm_distance) <- 1:nrow(cbc_perm_distance) ## change row names to numbers
cbc_perm_distance <- cbc_perm_distance[,-c(1:3)] ## remove unnecessary columns
## only interested in mobility main effect
cbc_perm_distance <- cbc_perm_distance[grep("distance_shuffle", cbc_perm_distance$parameter),]
final_results_table <- rbind(true_result_table, cbc_perm_distance) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/CBC/perm_distance_cbc.csv", row.names=TRUE)
## read in file
perm_distance_cbc <- read.csv("../Results/Model_outputs/CBC/perm_distance_cbc.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_distance_cbc <- perm_distance_cbc[perm_distance_cbc$i==0,]
perm_distance_cbc <- perm_distance_cbc[!perm_distance_cbc$i==0,] ## remove true value to calc. p value
diff.observed <- true_distance_cbc$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_distance_cbc$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant

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
spec_model_ukbms2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*specialism + (1|pair.id) + (1|spp), data = pair_attr_early)
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
  pair_attr_early$spec_shuffle <- sample(pair_attr_early$specialism) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*spec_shuffle + (1|pair.id) + (1|spp), data = pair_attr_early)
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
spec_model_ukbms2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*specialism + (1|pair.id) + (1|spp), data = pair_attr_late)
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
  pair_attr_late$spec_shuffle <- sample(pair_attr_late$specialism) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*spec_shuffle + (1|pair.id) + (1|spp), data = pair_attr_late)
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
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC_no_zeros2.csv", header=TRUE)
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
pair_attr_early <- rbind(pair_attr_1985, pair_attr_2000) ## 2 years and 32 species
pair_attr_late <- rbind(pair_attr_2000, pair_attr_2012) ## 2 years and 32 species
## remove NA's to make sure small white is taken out
pair_attr_early <- na.omit(pair_attr_early)
pair_attr_late <- na.omit(pair_attr_late)

pair_attr_early$mid.year <- as.factor(pair_attr_early$mid.year)
pair_attr_early$pair.id <- as.character(pair_attr_early$pair.id)
pair_attr_early$spp <- as.factor(pair_attr_early$spp)
pair_attr_early$mobility_wil <- as.numeric(pair_attr_early$mobility_wil)

pair_attr_late$mid.year <- as.factor(pair_attr_late$mid.year)
pair_attr_late$pair.id <- as.character(pair_attr_late$pair.id)
pair_attr_late$spp <- as.factor(pair_attr_late$spp)
pair_attr_late$mobility_wil <- as.numeric(pair_attr_late$mobility_wil)

## run model with 2 mobility groups (low and high)
pair_attr_early$mobility_wil <- as.numeric(pair_attr_early$mobility_wil)
pair_attr_early$mobility_score2 <- cut(pair_attr_early$mobility_wil, 2, labels=FALSE)
pair_attr_early$mobility_score2 <- as.factor(pair_attr_early$mobility_score2)

pair_attr_late$mobility_wil <- as.numeric(pair_attr_late$mobility_wil)
pair_attr_late$mobility_score2 <- cut(pair_attr_late$mobility_wil, 2, labels=FALSE)
pair_attr_late$mobility_score2 <- as.factor(pair_attr_late$mobility_score2)

###### EARLY #######
## run true model
mobility_model7 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr_early)
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
  pair_attr_early$mob_shuffle <- sample(pair_attr_early$mobility_score2) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mob_shuffle + (1|pair.id) + (1|spp), data = pair_attr_early)
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
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_change_mob_ukbms_85_00.csv", row.names=TRUE)
## read in file
perm_change_mob_ukbms <- read.csv("../Results/Model_outputs/UKBMS/perm_change_mob_ukbms_85_00.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_mob_ukbms <- perm_change_mob_ukbms[perm_change_mob_ukbms$i==0,]
perm_change_mob_ukbms <- perm_change_mob_ukbms[!perm_change_mob_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_mob_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_mob_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant

###### LATE #######
## run true model
mobility_model7 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr_late)
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
  pair_attr_late$mob_shuffle <- sample(pair_attr_late$mobility_score2) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mob_shuffle + (1|pair.id) + (1|spp), data = pair_attr_late)
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


####################################################
####################################################
#### CBC (Change in synchrony)
pair_attr_CBC <- merge(pair_attr_CBC, bird_dispersal, by.x="spp", by.y="Species_code", all=FALSE)
length(unique(pair_attr_CBC$spp)) # 22 species until remove NAs for Breeding dispersal 

pair_attr_cbc_1985 <- pair_attr_CBC[pair_attr_CBC$mid.year==1984.5,]
pair_attr_cbc_1996 <- pair_attr_CBC[pair_attr_CBC$mid.year==1995.5,]
pair_attr_cbc <- rbind(pair_attr_cbc_1985, pair_attr_cbc_1996)
pair_attr_cbc <- na.omit(pair_attr_cbc) # remove NA's (species without data)
length(unique(pair_attr_cbc$spp)) # 21 species 

pair_attr_cbc$mid.year <- as.factor(pair_attr_cbc$mid.year)
pair_attr_cbc$pair.id <- as.character(pair_attr_cbc$pair.id)
pair_attr_cbc$spp <- as.factor(pair_attr_cbc$spp)

## run model with 2 mobility groups (low and high)
pair_attr_cbc$Breeding_AM <- as.numeric(pair_attr_cbc$Breeding_AM)
pair_attr_cbc$Breeding_AM_score2 <- cut(pair_attr_cbc$Breeding_AM, 2, labels=FALSE)
pair_attr_cbc$Breeding_AM_score2 <- as.factor(pair_attr_cbc$Breeding_AM_score2)

## run true model
dispersal_model_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*Breeding_AM_score2 + (1|spp) + (1|pair.id), data=pair_attr_cbc)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(dispersal_model_cbc)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## remove rows with mean northing, distance, renk hab sim, mobility and mid year F values (only interested in mobility*midyear interaction)
main_result_table <- main_result_table[ !(main_result_table$parameter %in% c("mean_northing", "distance", "hab_sim", "mid.year", "Breeding_AM_score2")), ]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores <- parallel::detectCores()
cl <- makeSOCKcluster(cores[1]-1) # not to overload your computer
registerDoSNOW(cl)

n_sims <- 999

pb <- txtProgressBar(max = n_sims, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress) 

para_start_time = Sys.time()
cbc_mob_para <- foreach (i=1:n_sims, .options.snow = opts, .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_cbc$mob_shuffle <- sample(pair_attr_cbc$Breeding_AM_score2) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*mob_shuffle + (1|pair.id) + (1|spp), data = pair_attr_cbc)
  ## run model with shuffled variable
  ## save results
  anoresult<-anova(model)
  data.frame(anoresult, i=i)
}
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 6.94 minutes for 999 runs

### save results
cbc_mob_para$parameter <- paste(row.names(cbc_mob_para)) ## move row.names to parameter column
rownames(cbc_mob_para) <- 1:nrow(cbc_mob_para) ## change row names to numbers
cbc_mob_para <- cbc_mob_para[,-c(1:3)] ## remove unnecessary columns
## only interested in midyear*mobility interaction
cbc_mob_para <- cbc_mob_para[grep("mid.year:mob_shuffle", cbc_mob_para$parameter),]

final_results_table <- rbind(main_result_table, cbc_mob_para) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/CBC/perm_change_mob_cbc.csv", row.names=TRUE)
## read in file
perm_change_mob_cbc <- read.csv("../Results/Model_outputs/CBC/perm_change_mob_cbc.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_mob_cbc <- perm_change_mob_cbc[perm_change_mob_cbc$i==0,]
perm_change_mob_cbc <- perm_change_mob_cbc[!perm_change_mob_cbc$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_mob_cbc$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_mob_cbc$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0.058 non-significant


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

##### Average abundance (CBC)
## read in data
bird_common <- read.csv("../Data/BTO_data/pop_estimates_birds.csv", header=TRUE)
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC_no_zeros2.csv", header=TRUE)

## merge the two datasets
pair_attr_CBC <- merge(pair_attr_CBC, bird_common, by.x="spp", by.y="species_code", all=FALSE)
## rescale variables
pair_attr_CBC$pop_estimate_log <- log(pair_attr_CBC$pop_estimate)

## true model
common_model_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + pop_estimate_log + (1|pair.id) + (1|spp), data = pair_attr_CBC)
## save true model results (to merge in with bootstrapped models later)
main_result_table <- data.frame(anova(common_model_cbc)[5]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## remove rows with mean northing, distance, renk hab sim and mid year F values (only interested in abundance)
main_result_table <- main_result_table[ !(main_result_table$parameter %in% c("mean_northing", "distance", "hab_sim", "mid.year")), ]

## run 999 permutation tests
## Set up number of cores to run on (7)
cores <- parallel::detectCores()
cl <- makeSOCKcluster(cores[1]-1) # not to overload your computer
registerDoSNOW(cl)

n_sims <- 999

pb <- txtProgressBar(max = n_sims, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress) 

para_start_time = Sys.time()
cbc_abund_para <- foreach (i=1:n_sims, .options.snow = opts, .combine=rbind, .packages='lme4') %dopar% {
  print(i)
  pair_attr_CBC$abund_shuffle <- sample(pair_attr_CBC$pop_estimate_log) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + abund_shuffle + (1|pair.id) + (1|spp), data = pair_attr_CBC)
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
pvalue ## 0.092 non-significant

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
pair_attr_early <- rbind(pair_attr_1985, pair_attr_2000)
pair_attr_late <- rbind(pair_attr_2000, pair_attr_2012)
length(unique(pair_attr_early$spp)) ## 2 years and 31 species
length(unique(pair_attr_late$spp)) ## 2 years and 31 species

pair_attr_early$mid.year <- as.factor(pair_attr_early$mid.year)
pair_attr_early$pair.id <- as.character(pair_attr_early$pair.id)
pair_attr_early$spp <- as.factor(pair_attr_early$spp)

pair_attr_late$mid.year <- as.factor(pair_attr_late$mid.year)
pair_attr_late$pair.id <- as.character(pair_attr_late$pair.id)
pair_attr_late$spp <- as.factor(pair_attr_late$spp)

##### EARLY #####
## true model
climate_ukbms_change <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI + (1|spp) + (1|pair.id), data=pair_attr_early)
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
  pair_attr_early$STI_shuffle <- sample(pair_attr_early$STI) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI_shuffle + (1|pair.id) + (1|spp), data = pair_attr_early)
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
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_change_STI_ukbms_85_00.csv", row.names=TRUE)
## read in file
perm_change_STI_ukbms <- read.csv("../Results/Model_outputs/UKBMS/perm_change_STI_ukbms_85_00.csv", header=TRUE) 

## Calculate p value
number_of_permutations <- 1000
true_change_STI_ukbms <- perm_change_STI_ukbms[perm_change_STI_ukbms$i==0,]
perm_change_STI_ukbms <- perm_change_STI_ukbms[!perm_change_STI_ukbms$i==0,] ## remove true value to calc. p value
diff.observed <- true_change_STI_ukbms$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(perm_change_STI_ukbms$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 significant

##### LATE #####
## true model
climate_ukbms_change <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI + (1|spp) + (1|pair.id), data=pair_attr_late)
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
  pair_attr_late$STI_shuffle <- sample(pair_attr_late$STI) ## randomly shuffle specialism varaible
  model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI_shuffle + (1|pair.id) + (1|spp), data = pair_attr_late)
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
write.csv(final_results_table, file = "../Results/Model_outputs/UKBMS/perm_change_STI_ukbms_00_12.csv", row.names=TRUE)
## read in file
perm_change_STI_ukbms <- read.csv("../Results/Model_outputs/UKBMS/perm_change_STI_ukbms_00_12.csv", header=TRUE) 

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
pvalue ## 0.04 significant

