###################################################################
## Title: Testing model significance & creating barcharts BBS data
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: May 2018
###################################################################

rm(list=ls()) # clear R

library(lme4)
library(lmerTest)
library(dplyr)
library(foreach)
library(doParallel)
library(doSNOW)
options(scipen=999)


## load data
pair_attr <- read.csv("../Data/Bird_sync_data/pop_climate_synchrony_BBS.csv", header=TRUE) # BBS pair attribute data

### just want to compare 1985 with 1996 
## to see if there has been an increase/decrease/no change in connectivity

## create new pair_attr file with just those years
pair_attr_1999 <- pair_attr[pair_attr$mid.year==1998.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]
pair_attr_comp <- rbind(pair_attr_1999, pair_attr_2012) ## rbind them together
summary(pair_attr_comp)
unique(pair_attr_comp$mid.year) ## 1998.5 and 2011.5

## write file to save time
write.csv(pair_attr_comp, file = "../Data/Bird_sync_data/pair_attr_comp_BBS.csv", row.names = FALSE) 

## read in file
pair_attr_comp <- read.csv("../Data/Bird_sync_data/pair_attr_comp_BBS.csv", header=TRUE)

#### change variables to factors/characters
str(pair_attr_comp)
pair_attr_comp$mid.year <- as.factor(pair_attr_comp$mid.year)
pair_attr_comp$pair.id <- as.character(pair_attr_comp$pair.id)
pair_attr_comp$spp <- as.factor(pair_attr_comp$spp)

## centre and standardise continuous variables by mean and SD
pair_attr_comp$distance <- (pair_attr_comp$distance - mean(na.omit(pair_attr_comp$distance)))/sd(na.omit(pair_attr_comp$distance))
pair_attr_comp$mean_northing <- (pair_attr_comp$mean_northing - mean(na.omit(pair_attr_comp$mean_northing)))/sd(na.omit(pair_attr_comp$mean_northing))
pair_attr_comp$renk_hab_sim <- (pair_attr_comp$renk_hab_sim - mean(na.omit(pair_attr_comp$renk_hab_sim)))/sd(na.omit(pair_attr_comp$renk_hab_sim))

############################################
#### model comparison for each species ####
############################################

spp.list<- unique(pair_attr_comp$spp)

results_table<-NULL
F_result_final <- NULL
for (i in spp.list){
  print(i)
  
  ## create unique pair_attr for each species
  p_attr <- pair_attr_comp[pair_attr_comp$spp==i,]
  
 # if loop (if nrow when start.year=1980 is <20 then skip that species)
  if(nrow(p_attr[p_attr$mid.year=="1998.5",]) < 3){
    print(paste("skip species",i))
    next

  } else {

  model_comp <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + mid.year + (1|pair.id), data = pair_attr_comp[pair_attr_comp$spp==i,])
  summary(model_comp)
  anova(model_comp)
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(model_comp)$coefficients[,1:5],i)
  results_table <-rbind(results_table,results_table_temp)
  F_result_table <- data.frame(anova(model_comp)[5], i) ## also save model F values to use in permutation tests
  F_result_final <-rbind(F_result_final,F_result_table)
  
  } ## species 472 does work ==> 23 species in total now 
}

## save results
names(results_table) <- c("Estimate", "SD", "df", "t","p_value", "species")
results_table$parameter <- paste(row.names(results_table))
rownames(results_table) <- 1:nrow(results_table)
## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_table <- results_table[grep("mid.year", results_table$parameter),]
## change parameter names to year of interest 
results_table$parameter <- 1999/2012
## make final table
results_final <- results_table[,c(1,6)]

##### modify F result table
colnames(F_result_final)[2] <- "j"
F_result_final$i <- 0 ## make i column with zeros 
F_result_final$parameter <- paste(row.names(F_result_final)) ## move row.names to parameter column
rownames(F_result_final) <- 1:nrow(F_result_final) ## change row names to numbers
## take out 3 columns for each species: mean northing, distance and renk_hab_sim
F_result_final <- F_result_final[grep("mid.year", F_result_final$parameter),]

## save true model results
write.csv(F_result_final, file = "../Results/Model_outputs/BBS/true_F_values_spp.csv", row.names=FALSE)
F_result_final <- read.csv("../Results/Model_outputs/BBS/true_F_values_spp.csv", header=TRUE)

## permutation tests for each species
## run 999 permutation tests
n_sims <- 999

## only run on species which have enough data (same ones as above)
F_result_final$j <- as.factor(F_result_final$j)
spp_list <- unique(F_result_final$j) ## only run this on the species that work (23)
## merge spp with pair_attr_comp
spp_list <- as.data.frame(spp_list)
pair_attr_comp <- merge(pair_attr_comp, spp_list, by.x="spp", by.y="spp_list", all=FALSE)
pair_attr_comp <- droplevels(pair_attr_comp)
length(unique(pair_attr_comp$spp)) ## 23
spp_list <- unique(F_result_final$j) ## only run this on the species that work (23)

## Set up number of cores to run on (7)
cores <- parallel::detectCores()
cl <- makeCluster(cores[1]-1) # not to overload your computer
registerDoSNOW(cl)

## this code creates a progress percentage bar
pb <- txtProgressBar(max = n_sims, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress) 
## run parallel
para_start_time = Sys.time()
bbs_spp_para <- foreach (i=1:n_sims,  .combine=rbind, .options.snow = opts, .packages=c('lme4', 'foreach')) %dopar% { ## loop through each permutation
  print(i)
  foreach (j=spp_list, .combine=rbind, .options.snow = opts, .packages=c('lme4', 'foreach')) %dopar% { ## loop through each species
    print(j)
    pair_attr_comp$mid.year_shuffle <- sample(pair_attr_comp$mid.year) ## randomly shuffle mid year variable
    model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + mid.year_shuffle + (1|pair.id), data = pair_attr_comp[pair_attr_comp$spp==j,])
    ## run model with shuffled variable
    ## save results
    anoresult<-anova(model)
    data.frame(anoresult, i=i, j=j)
  }
}
close(pb)
stopCluster(cl)
para_end_time = Sys.time()
para_run_time = para_end_time - para_start_time
print(paste0("TOTAL RUN TIME: ", para_run_time)) ## 5.9 hours for 999 runs
#### end

### save results
bbs_spp_para$parameter <- paste(row.names(bbs_spp_para)) ## move row.names to parameter column
rownames(bbs_spp_para) <- 1:nrow(bbs_spp_para) ## change row names to numbers
bbs_spp_para <- bbs_spp_para[,-c(1:3)] ## remove unnecessary columns
## only keep rows with interaction
bbs_spp_para <- bbs_spp_para[grep("mid.year_shuffle", bbs_spp_para$parameter),]
final_results_table <- rbind(F_result_final, bbs_spp_para) ## bind the two data frames together

final_results_table_429 <- final_results_table[final_results_table$j==429,]
F_value <- with(final_results_table_429, final_results_table_429$F.value[final_results_table_429$i==0]) ## true F value from main model
hist(final_results_table_429$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/BBS/perm_F_values_spp.csv", row.names=FALSE)
## read in file
perm_F_values_spp <- read.csv("../Results/Model_outputs/BBS/perm_F_values_spp.csv", header=TRUE) 

## calculate p value for each species
## Calculate p value
number_of_permutations <- 1000
spp_list <- unique(perm_F_values_spp$j)
p_values_final <- NULL
## loop through each species
for(j in spp_list){
  
  perm_F_values_spp2 <- perm_F_values_spp[perm_F_values_spp$j==j,] ## select species of interest
  true_bbs_spp <- perm_F_values_spp2[perm_F_values_spp2$i==0,] ## obtain true F value (iteration = 0)
  perm_F_values_spp2 <- perm_F_values_spp2[!perm_F_values_spp2$i==0,] ## remove true value to calc. p value
  diff.observed <- true_bbs_spp$F.value ## true F value
  
  # P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
  pvalue = sum(abs(perm_F_values_spp2$F.value) >= abs(diff.observed)) / number_of_permutations
  p_values_temp <- data.frame(pvalue=pvalue, species=j)
  p_values_final <- rbind(p_values_final, p_values_temp)
}
## save file
write.csv(p_values_final, file = "../Results/Model_outputs/BBS/perm_p_values_spp.csv", row.names=FALSE)
p_values_final <- read.csv("../Results/Model_outputs/BBS/perm_p_values_spp.csv", header=TRUE)

## merge permutation p values with main model output
final_results_bbs <- merge(results_final, p_values_final, by="species")

## classify each species as either no change (insignificant p value), decreasing (negative estimte), or increasing (positive estimate)
final_results_bbs <- final_results_bbs %>% group_by(species) %>% mutate(change = ifelse(pvalue>0.05, "No change", ifelse(pvalue<0.05 & Estimate<0, "Decrease", "Increase")))

nrow(final_results_bbs[final_results_bbs$change=="No change",]) ## 17 species unchanged
nrow(final_results_bbs[final_results_bbs$change=="Decrease",]) ## 4 species decreasing
nrow(final_results_bbs[final_results_bbs$change=="Increase",]) ## 2 species increasing

write.csv(final_results_bbs, file="../Results/Bird_results/bbs_spp_trends.csv", row.names=FALSE)

final_results_bbs2 <- data.frame(change=unique(final_results_bbs$change), no_species=c(17,4,2), total_species=23)
final_results_bbs2$percentage <- (final_results_bbs2$no_species / final_results_bbs2$total_species)*100

## save file
write.csv(final_results_bbs2, file="../Results/Bird_results/bbs_percentage_trends.csv", row.names=FALSE)
