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
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr_no_zeros2.csv", header=TRUE)

##### Comparing between early, mid and late years:
## EARLY: start.year = 1980
## MID: start.year = 1995
## LATE: start.year = 2007

length(unique(pair_attr$spp)) # 32 species

###### create 3 new pair_attr files whicih compares early, late and overall
pair_attr_1980 <- pair_attr[pair_attr$start.year==1980,]
pair_attr_1995 <- pair_attr[pair_attr$start.year==1995,]
pair_attr_2007 <- pair_attr[pair_attr$start.year==2007,]

pair_attr_early <- rbind(pair_attr_1980, pair_attr_1995) # comparison of early years 1980 & 1995
pair_attr_late <- rbind(pair_attr_1995, pair_attr_2007) # comparison of late years 1995 & 2007
pair_attr_overall <- rbind(pair_attr_1980, pair_attr_2007) # overall comparison from 1980 to 2007

## check new pair_attr files only contain 2 sets of start.years
unique(pair_attr_early$start.year) ## 1980 1995
unique(pair_attr_late$start.year) ## 1995 2007
unique(pair_attr_overall$start.year) ## 1980 2007

## centre and standardise continuous variables by mean and SD
##early
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
write.csv(pair_attr_early, file = "../Data/Butterfly_sync_data/pair_attr_early_no_zeros2.csv", row.names = FALSE) # save pair_attr_early file 
write.csv(pair_attr_late, file = "../Data/Butterfly_sync_data/pair_attr_late_no_zeros2.csv", row.names = FALSE) # save pair_attr_late file 
write.csv(pair_attr_overall, file = "../Data/Butterfly_sync_data/pair_attr_overall_no_zeros2.csv", row.names = FALSE) # save pair_attr_overall file 

## read in files to save time 
pair_attr_early <- read.csv("../Data/Butterfly_sync_data/pair_attr_early_no_zeros2.csv", header=TRUE)
pair_attr_late <- read.csv("../Data/Butterfly_sync_data/pair_attr_late_no_zeros2.csv", header=TRUE)
pair_attr_overall <- read.csv("../Data/Butterfly_sync_data/pair_attr_overall_no_zeros2.csv", header=TRUE)

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


# ##############################
# ######## run models ##########
# ############################## 
# 
# ## likelihood ratio test to test if year has a significant OVERALL effect ##
# 
# ## EARLY COMPARISON ##
# early_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), REML=FALSE, data = pair_attr_early)
# early_model_null <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + (1|pair.id) + (1|spp), REML=FALSE, data = pair_attr_early)
# ## model comparison
# model_comp_early <- anova(early_model_null, early_model, test="Chisq")
# ## add comparison column
# model_comp_early$comparison <- rep("Early comparison")
# ## create dataframe with results
# results_early_all_spp <- data.frame(summary(early_model)$coefficients[,1:5])
# names(results_early_all_spp) <- c("Estimate", "SD", "df", "t","p_value")
# results_early_all_spp$parameter <- paste(row.names(results_early_all_spp))
# rownames(results_early_all_spp) <- 1:nrow(results_early_all_spp)
# 
# ## take out 3 columns for each species: mean northing, distance and renk_hab_sim
# results_table5 <- NULL
# results_table1 <- results_early_all_spp[grep("mean_northing", results_early_all_spp$parameter),]
# results_table2 <- results_early_all_spp[grep("distance", results_early_all_spp$parameter),]
# results_table3 <- results_early_all_spp[grep("hab_sim", results_early_all_spp$parameter),]
# results_table4 <- results_early_all_spp[grep("(Intercept)", results_early_all_spp$parameter),]
# results_table5 <- rbind(results_table5, results_table4, results_table1, results_table2, results_table3)
# results_early_all_spp <- results_early_all_spp[!results_early_all_spp$parameter%in%results_table5$parameter,]
# 
# ## create year to the comparison (years 1985 and 2000)
# results_early_all_spp$year <- "1985/2000"
# ## remove parameter column
# results_early_all_spp <- results_early_all_spp[-6]
# results_early_all_spp$comparison <- "Early short term"
# results_early_all_spp$trend <- "significant decrease"
# 
# ## LATE COMPARISON ##
# late_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), REML=FALSE, data = pair_attr_late)
# late_model_null <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + (1|pair.id) + (1|spp), REML=FALSE, data = pair_attr_late)
# ## model comparison
# model_comp_late <- anova(late_model_null, late_model, test="Chisq")
# ## add comparison column
# model_comp_late$comparison <- rep("Late comparison")
# ## create dataframe with results
# results_late_all_spp <- data.frame(summary(late_model)$coefficients[,1:5])
# names(results_late_all_spp) <- c("Estimate", "SD", "df", "t","p_value")
# results_late_all_spp$parameter <- paste(row.names(results_late_all_spp))
# rownames(results_late_all_spp) <- 1:nrow(results_late_all_spp)
# 
# ## take out 3 columns for each species: mean northing, distance and renk_hab_sim
# results_table5 <- NULL
# results_table1 <- results_late_all_spp[grep("mean_northing", results_late_all_spp$parameter),]
# results_table2 <- results_late_all_spp[grep("distance", results_late_all_spp$parameter),]
# results_table3 <- results_late_all_spp[grep("hab_sim", results_late_all_spp$parameter),]
# results_table4 <- results_late_all_spp[grep("(Intercept)", results_late_all_spp$parameter),]
# results_table5 <- rbind(results_table5, results_table4, results_table1, results_table2, results_table3)
# results_late_all_spp <- results_late_all_spp[!results_late_all_spp$parameter%in%results_table5$parameter,]
# 
# ## create year to the comparison (years 1985 and 2000)
# results_late_all_spp$year <- "2000/2012"
# ## remove parameter column
# results_late_all_spp <- results_late_all_spp[-6]
# results_late_all_spp$comparison <- "Late short term"
# results_late_all_spp$trend <- "significant increase"
# 
# ## OVERALL COMPARISON ##
# overall_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), REML=FALSE, data = pair_attr_overall)
# overall_model_null <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + (1|pair.id) + (1|spp), REML=FALSE, data = pair_attr_overall)
# ## model comparison
# model_comp_overall <- anova(overall_model_null, overall_model, test="Chisq")
# ## add comparison column
# model_comp_overall$comparison <- rep("Overall comparison")
# ## create dataframe with results
# results_overall_all_spp <- data.frame(summary(overall_model)$coefficients[,1:5])
# names(results_overall_all_spp) <- c("Estimate", "SD", "df", "t","p_value")
# results_overall_all_spp$parameter <- paste(row.names(results_overall_all_spp))
# rownames(results_overall_all_spp) <- 1:nrow(results_overall_all_spp)
# 
# ## take out 3 columns for each species: mean northing, distance and renk_hab_sim
# results_table5 <- NULL
# results_table1 <- results_overall_all_spp[grep("mean_northing", results_overall_all_spp$parameter),]
# results_table2 <- results_overall_all_spp[grep("distance", results_overall_all_spp$parameter),]
# results_table3 <- results_overall_all_spp[grep("hab_sim", results_overall_all_spp$parameter),]
# results_table4 <- results_overall_all_spp[grep("(Intercept)", results_overall_all_spp$parameter),]
# results_table5 <- rbind(results_table5, results_table4, results_table1, results_table2, results_table3)
# results_overall_all_spp <- results_overall_all_spp[!results_overall_all_spp$parameter%in%results_table5$parameter,]
# 
# ## create year to the comparison (years 1985 and 2000)
# results_overall_all_spp$year <- "1985/2012"
# ## remove parameter column
# results_overall_all_spp <- results_overall_all_spp[-6]
# results_overall_all_spp$comparison <- "Long term"
# results_overall_all_spp$trend <- "significant increase"
# 
# ## save model comparison results for all species
# model_comp_all_spp <- data.frame(rbind(results_early_all_spp, results_late_all_spp, results_overall_all_spp))
# write.csv(model_comp_all_spp, file="../Results/Butterfly_results/model_comp_all_spp.csv", row.names=FALSE)
# 
# #### save log likelihood test results for each comparison for all species #####
# ratio_test_all_spp <- data.frame(rbind(model_comp_early, model_comp_late, model_comp_overall))
# ratio_test_all_spp <- ratio_test_all_spp[-c(1,3,5),]
# rownames(ratio_test_all_spp) <- 1:nrow(ratio_test_all_spp)
# 
# write.csv(ratio_test_all_spp, file="../Results/Butterfly_results/ratio_test_all_spp.csv", row.names=FALSE)
# #### results = very significant for all 3 comparisons #####
# 
# ##################################################################################################
# ####################### SAME AS ABOVE BUT FOR WOODLAND SPECIES ONLY ##############################
# ##################################################################################################
# 
# ### create new pair_attr files with woodland species only
# pair_attr_early_wood <- pair_attr_early[pair_attr_early$HABITAT=="Woodland",]
# pair_attr_early_wood <- droplevels(pair_attr_early_wood)
# 
# pair_attr_late_wood <- pair_attr_late[pair_attr_late$HABITAT=="Woodland",]
# pair_attr_late_wood <- droplevels(pair_attr_late_wood)
# 
# pair_attr_overall_wood <- pair_attr_overall[pair_attr_overall$HABITAT=="Woodland",]
# pair_attr_overall_wood <- droplevels(pair_attr_overall_wood)
# 
# ## EARLY COMPARISON ##
# early_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), REML=FALSE, data = pair_attr_early_wood)
# early_model_null <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + (1|pair.id) + (1|spp), REML=FALSE, data = pair_attr_early_wood)
# ## model comparison
# model_comp_early <- anova(early_model_null, early_model, test="Chisq")
# ## add comparison column
# model_comp_early$comparison <- rep("Early comparison")
# ## create dataframe with results
# results_early_all_spp <- data.frame(summary(early_model)$coefficients[,1:5])
# names(results_early_all_spp) <- c("Estimate", "SD", "df", "t","p_value")
# results_early_all_spp$parameter <- paste(row.names(results_early_all_spp))
# rownames(results_early_all_spp) <- 1:nrow(results_early_all_spp)
# 
# ## take out 3 columns for each species: mean northing, distance and renk_hab_sim
# results_table5 <- NULL
# results_table1 <- results_early_all_spp[grep("mean_northing", results_early_all_spp$parameter),]
# results_table2 <- results_early_all_spp[grep("distance", results_early_all_spp$parameter),]
# results_table3 <- results_early_all_spp[grep("hab_sim", results_early_all_spp$parameter),]
# results_table4 <- results_early_all_spp[grep("(Intercept)", results_early_all_spp$parameter),]
# results_table5 <- rbind(results_table5, results_table4, results_table1, results_table2, results_table3)
# results_early_all_spp <- results_early_all_spp[!results_early_all_spp$parameter%in%results_table5$parameter,]
# 
# ## create year to the comparison (years 1985 and 2000)
# results_early_all_spp$year <- "1985/2000"
# ## remove parameter column
# results_early_all_spp <- results_early_all_spp[-6]
# results_early_all_spp$comparison <- "Early short term"
# results_early_all_spp$trend <- "significant decrease"
# 
# ## LATE COMPARISON ##
# late_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), REML=FALSE, data = pair_attr_late_wood)
# late_model_null <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + (1|pair.id) + (1|spp), REML=FALSE, data = pair_attr_late_wood)
# ## model comparison
# model_comp_late <- anova(late_model_null, late_model, test="Chisq")
# ## add comparison column
# model_comp_late$comparison <- rep("Late comparison")
# ## create dataframe with results
# results_late_all_spp <- data.frame(summary(late_model)$coefficients[,1:5])
# names(results_late_all_spp) <- c("Estimate", "SD", "df", "t","p_value")
# results_late_all_spp$parameter <- paste(row.names(results_late_all_spp))
# rownames(results_late_all_spp) <- 1:nrow(results_late_all_spp)
# 
# ## take out 3 columns for each species: mean northing, distance and renk_hab_sim
# results_table5 <- NULL
# results_table1 <- results_late_all_spp[grep("mean_northing", results_late_all_spp$parameter),]
# results_table2 <- results_late_all_spp[grep("distance", results_late_all_spp$parameter),]
# results_table3 <- results_late_all_spp[grep("hab_sim", results_late_all_spp$parameter),]
# results_table4 <- results_late_all_spp[grep("(Intercept)", results_late_all_spp$parameter),]
# results_table5 <- rbind(results_table5, results_table4, results_table1, results_table2, results_table3)
# results_late_all_spp <- results_late_all_spp[!results_late_all_spp$parameter%in%results_table5$parameter,]
# 
# ## create year to the comparison (years 1985 and 2000)
# results_late_all_spp$year <- "2000/2012"
# ## remove parameter column
# results_late_all_spp <- results_late_all_spp[-6]
# results_late_all_spp$comparison <- "Late short term"
# results_late_all_spp$trend <- "significant increase"
# 
# ## OVERALL COMPARISON ##
# overall_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), REML=FALSE, data = pair_attr_overall_wood)
# overall_model_null <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), REML=FALSE, data = pair_attr_overall_wood)
# ## model comparison
# model_comp_overall <- anova(overall_model_null, overall_model, test="Chisq")
# ## add comparison column
# model_comp_overall$comparison <- rep("Overall comparison")
# ## create dataframe with results
# results_overall_all_spp <- data.frame(summary(overall_model)$coefficients[,1:5])
# names(results_overall_all_spp) <- c("Estimate", "SD", "df", "t","p_value")
# results_overall_all_spp$parameter <- paste(row.names(results_overall_all_spp))
# rownames(results_overall_all_spp) <- 1:nrow(results_overall_all_spp)
# 
# ## take out 3 columns for each species: mean northing, distance and renk_hab_sim
# results_table5 <- NULL
# results_table1 <- results_overall_all_spp[grep("mean_northing", results_overall_all_spp$parameter),]
# results_table2 <- results_overall_all_spp[grep("distance", results_overall_all_spp$parameter),]
# results_table3 <- results_overall_all_spp[grep("hab_sim", results_overall_all_spp$parameter),]
# results_table4 <- results_overall_all_spp[grep("(Intercept)", results_overall_all_spp$parameter),]
# results_table5 <- rbind(results_table5, results_table4, results_table1, results_table2, results_table3)
# results_overall_all_spp <- results_overall_all_spp[!results_overall_all_spp$parameter%in%results_table5$parameter,]
# 
# ## create year to the comparison (years 1985 and 2000)
# results_overall_all_spp$year <- "1985/2012"
# ## remove parameter column
# results_overall_all_spp <- results_overall_all_spp[-6]
# results_overall_all_spp$comparison <- "Long term"
# results_overall_all_spp$trend <- "significant increase"
# 
# ## save model comparison results for all species
# model_comp_all_spp_wood <- data.frame(rbind(results_early_all_spp, results_late_all_spp, results_overall_all_spp))
# write.csv(model_comp_all_spp_wood, file="../Results/Butterfly_results/model_comp_all_spp_woodland.csv", row.names=FALSE)


##################################
## EARLY MODEL FOR EACH SPECIES ##
##################################

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

names(results_table_early) <- c("Estimate", "SD", "df", "t","p_value", "species")
results_table_early$parameter <- paste(row.names(results_table_early))
rownames(results_table_early) <- 1:nrow(results_table_early)

## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_table_early <- results_table_early[grep("mid.year", results_table_early$parameter),]

## change parameter names to year of interest (mid.year = 2000/1999.5)
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

## merge permutation p values with main model output
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
write.csv(final_results_early2, file="../Results/Butterfly_results/ukbms_overall_trends_85_00.csv", row.names=FALSE)


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
write.csv(final_results_late2, file="../Results/Butterfly_results/ukbms_overall_trends_00_12.csv", row.names=FALSE)



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
write.csv(final_results_overall2, file="../Results/Butterfly_results/ukbms_spp_trends_85_12.csv", row.names=FALSE)






#### save comparison model results for each comparison for each species #####
results_final_early <- read.csv("../Results/Butterfly_results/ukbms_overall_trends_85_00.csv", header=TRUE)
results_final_late <- read.csv("../Results/Butterfly_results/ukbms_overall_trends_00_12.csv", header=TRUE)
results_final_overall <- read.csv("../Results/Butterfly_results/ukbms_overall_trends_85_12.csv", header=TRUE)

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
write.csv(model_comp_results, file="../Results/Butterfly_results/model_comp_percentages_no_zeros2.csv", row.names=FALSE)

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
png("../Graphs/Model_comps/Model_comp_results_all_spp_no_zeros2.png", height = 100, width = 120, units = "mm", res = 300)
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
png("../Graphs/Model_comps/Model_comp_results_all_spp3_no_zeros2.png", height = 100, width = 120, units = "mm", res = 300)
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
# 
# #### produce graph for percentage of specialists and generalists 
# ## add specialist/generalist data
# spec_gen <- read.csv("../Data/UKBMS_data/UKBMS_UKspecieslist.csv", header=TRUE)
# 
# ## merge datasets
# results_table_early <- merge(results_table_early, spec_gen, by.x="species", by.y="BMSCODE")
# results_table_early <- results_table_early[-c(6,7,9,10)] # remove unecessary columns
# results_table_early_wc <- results_table_early[results_table_early$STRATEGY=="Wider countryside sp",]
# results_table_early_hs <- results_table_early[results_table_early$STRATEGY=="Habitat specialist",]
# 
# results_table_late <- merge(results_table_late, spec_gen, by.x="species", by.y="BMSCODE")
# results_table_late <- results_table_late[-c(6,7,9,10)] # remove unecessary columns
# results_table_late_wc <- results_table_late[results_table_late$STRATEGY=="Wider countryside sp",]
# results_table_late_hs <- results_table_late[results_table_late$STRATEGY=="Habitat specialist",]
# 
# results_table_overall <- merge(results_table_overall, spec_gen, by.x="species", by.y="BMSCODE")
# results_table_overall <- results_table_overall[-c(6,7,9,10)] # remove unecessary columns
# results_table_overall_wc <- results_table_overall[results_table_overall$STRATEGY=="Wider countryside sp",]
# results_table_overall_hs <- results_table_overall[results_table_overall$STRATEGY=="Habitat specialist",]
# 
# ## produce percentages of species (EARLY WIDER COUNTRYSIDE)
# results_final_early_wc <- data.frame(comparison="Early short term", change=unique(results_table_early_wc$change), strategy="Wider countryside sp", no_species=c(13,6,2), total_species=21)
# nrow(results_table_early_wc[results_table_early_wc$change=="No change",]) ## 6 wider countryside species no change
# nrow(results_table_early_wc[results_table_early_wc$change=="Decrease",]) ## 13 wider countryside species decreasing
# nrow(results_table_early_wc[results_table_early_wc$change=="Increase",]) ## 2 wider counytrside species increasing
# results_final_early_wc$percentage <- (results_final_early_wc$no_species / results_final_early_wc$total_species)*100
# 
# ## produce percentages of species (EARLY SPECIALISTS)
# results_final_early_hs <- data.frame(comparison="Early short term", change=unique(results_table_early_wc$change), strategy="Habitat specialist", no_species=c(0,6,2), total_species=8)
# nrow(results_table_early_hs[results_table_early_hs$change=="No change",]) ## 6 wider countryside species no change
# nrow(results_table_early_hs[results_table_early_hs$change=="Decrease",]) ## 0 wider countryside species decreasing
# nrow(results_table_early_hs[results_table_early_hs$change=="Increase",]) ## 2 wider counytrside species increasing
# results_final_early_hs$percentage <- (results_final_early_hs$no_species / results_final_early_hs$total_species)*100
# 
# ## produce percentages of species (LATE WIDER COUNTRYSIDE)
# results_final_late_wc <- data.frame(comparison="Late short term", change=unique(results_table_late_wc$change), strategy="Wider countryside sp", no_species=c(3,17,2), total_species=22)
# nrow(results_table_late_wc[results_table_late_wc$change=="No change",]) ## 3 wider countryside species no change
# nrow(results_table_late_wc[results_table_late_wc$change=="Decrease",]) ## 2 wider countryside species decreasing
# nrow(results_table_late_wc[results_table_late_wc$change=="Increase",]) ## 17 wider counytrside species increasing
# results_final_late_wc$percentage <- (results_final_late_wc$no_species / results_final_late_wc$total_species)*100
# 
# ## produce percentages of species (LATE SPECIALISTS)
# results_final_late_hs <- data.frame(comparison="Late short term", change=unique(results_table_late_hs$change), strategy="Habitat specialist", no_species=c(2,5,1), total_species=8)
# nrow(results_table_late_hs[results_table_late_hs$change=="No change",]) ## 2 wider countryside species no change
# nrow(results_table_late_hs[results_table_late_hs$change=="Decrease",]) ## 1 wider countryside species decreasing
# nrow(results_table_late_hs[results_table_late_hs$change=="Increase",]) ## 5 wider counytrside species increasing
# results_final_late_hs$percentage <- (results_final_late_hs$no_species / results_final_late_hs$total_species)*100
# 
# ## produce percentages of species (OVERALL WIDER COUNTRYSIDE)
# results_final_overall_wc <- data.frame(comparison="Long term", change=unique(results_table_overall_wc$change), strategy="Wider countryside sp", no_species=c(5,9,6), total_species=20)
# nrow(results_table_overall_wc[results_table_overall_wc$change=="No change",]) ## 6 wider countryside species no change
# nrow(results_table_overall_wc[results_table_overall_wc$change=="Decrease",]) ## 5 wider countryside species decreasing
# nrow(results_table_overall_wc[results_table_overall_wc$change=="Increase",]) ## 9 wider counytrside species increasing
# results_final_overall_wc$percentage <- (results_final_overall_wc$no_species / results_final_overall_wc$total_species)*100
# 
# ## produce percentages of species (OVERALL SPECIALISTS)
# results_final_overall_hs <- data.frame(comparison="Long term", change=unique(results_table_overall_hs$change), strategy="Habitat specialist", no_species=c(1,1,5), total_species=7)
# nrow(results_table_overall_hs[results_table_overall_hs$change=="No change",]) ## 5 wider countryside species no change
# nrow(results_table_overall_hs[results_table_overall_hs$change=="Decrease",]) ## 1 wider countryside species decreasing
# nrow(results_table_overall_hs[results_table_overall_hs$change=="Increase",]) ## 1 wider counytrside species increasing
# results_final_overall_hs$percentage <- (results_final_overall_hs$no_species / results_final_overall_hs$total_species)*100
# 
# ## rbind all results together 
# model_comp_results_wc <- rbind(results_final_early_wc, results_final_late_wc, results_final_overall_wc)
# model_comp_results_hs <- rbind(results_final_early_hs, results_final_late_hs, results_final_overall_hs)
# ## save file
# write.csv(model_comp_results_woodland, file="../Results/Butterfly_results/model_comp_percentages_woodland.csv", row.names=FALSE)
# 
# ####### plot 2 graphs - one for wider countryside and one for specialists ##########
# library(ggplot2)
# 
# ## change levels of change so they go in correct order
# levels(model_comp_results_wc$change)
# model_comp_results_wc$change <- factor(model_comp_results_wc$change, levels=c("Increase", "No change", "Decrease"))
# levels(model_comp_results_wc$change)
# 
# levels(model_comp_results_hs$change)
# model_comp_results_hs$change <- factor(model_comp_results_hs$change, levels=c("Increase", "No change", "Decrease"))
# levels(model_comp_results_hs$change)
# 
# ## wider countryside species
# png("../Graphs/Model_comps/Model_comp_results_wc.png", height = 100, width = 120, units = "mm", res = 300)
# ggplot(data=model_comp_results_wc, aes(x=comparison, y=percentage, fill=change)) +
#   geom_bar(stat="identity", width=0.4) +
#   labs(y="Percentage of species", x="", fill="") +
#   scale_fill_manual(values = c("#339900", "#999999", "#990000")) +
#   theme_bw() +
#   theme(axis.text.x=element_text(colour = "black")) +
#   theme(axis.text.y=element_text(colour = "black")) +
#   scale_y_continuous(breaks = seq(0,100,10), expand = c(0, 0)) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
# dev.off()
# 
# ## habitat specialist species
# png("../Graphs/Model_comps/Model_comp_results_hs.png", height = 100, width = 120, units = "mm", res = 300)
# ggplot(data=model_comp_results_hs, aes(x=comparison, y=percentage, fill=change)) +
#   geom_bar(stat="identity", width=0.4) +
#   labs(y="Percentage of species", x="", fill="") +
#   scale_fill_manual(values = c("#339900", "#999999", "#990000")) +
#   theme_bw() +
#   theme(axis.text.x=element_text(colour = "black")) +
#   theme(axis.text.y=element_text(colour = "black")) +
#   scale_y_continuous(breaks = seq(0,100,10), expand = c(0, 0)) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
# dev.off()

# ##########################################################
# ##### Make the same graph for woodland species only ######
# ##########################################################
# 
# ## merge the results with habitat data
# habitat_data <- read.csv("../Data/UKBMS_data/UKBMS_UKspecieslist.csv", header = TRUE) # habitat data for each species
# results_table_early <- merge(results_table_early, habitat_data, by.x="species", by.y="BMSCODE")
# results_table_late <- merge(results_table_late, habitat_data, by.x="species", by.y="BMSCODE")
# results_table_overall <- merge(results_table_overall, habitat_data, by.x="species", by.y="BMSCODE")
# 
# ## subset woodland species only
# results_table_early <- results_table_early[results_table_early$HABITAT=="Woodland",] ## 7 out of 9
# results_table_late <- results_table_late[results_table_late$HABITAT=="Woodland",] ## 7 out of 9
# results_table_overall <- results_table_overall[results_table_overall$HABITAT=="Woodland",] ## 5 out of 9
# 
# ## leave 4 columns: species code, change, common name and strategy
# results_table_early <- results_table_early[,c(1,4:5,8)]
# results_table_late <- results_table_late[,c(1,4:5,8)]
# results_table_overall <- results_table_overall[,c(1,4:5,8)]
# 
# ## create early comparison data frame
# results_final_early <- data.frame(comparison="Early short term", change=unique(results_table_early$change), no_species=c(1,4,2), total_species=7)
# nrow(results_table_early[results_table_early$change=="No change",]) ## 4 species unchanged
# nrow(results_table_early[results_table_early$change=="Decrease",]) ## 2 species decreasing
# nrow(results_table_early[results_table_early$change=="Increase",]) ## 1 species increasing
# results_final_early$percentage <- (results_final_early$no_species / results_final_early$total_species)*100
# 
# ## create late comparison data frame
# results_final_late <- data.frame(comparison="Late short term", change=unique(results_table_late$change), no_species=c(5,1,1), total_species=7)
# nrow(results_table_late[results_table_late$change=="No change",]) ## 1 species unchanged
# nrow(results_table_late[results_table_late$change=="Decrease",]) ## 1 species decreasing
# nrow(results_table_late[results_table_late$change=="Increase",]) ## 5 species increasing
# results_final_late$percentage <- (results_final_late$no_species / results_final_late$total_species)*100
# 
# ## create overall comparison data frame
# results_final_overall <- data.frame(comparison="Long term", change=unique(results_table_overall$change), no_species=c(2,1,2), total_species=5)
# nrow(results_table_overall[results_table_overall$change=="No change",]) ## 2 species unchanged
# nrow(results_table_overall[results_table_overall$change=="Decrease",]) ## 1 species decreasing
# nrow(results_table_overall[results_table_overall$change=="Increase",]) ## 2 species increasing
# results_final_overall$percentage <- (results_final_overall$no_species / results_final_overall$total_species)*100
# 
# ## rbind all results together 
# model_comp_results_woodland <- rbind(results_final_early, results_final_late, results_final_overall)
# ## save file
# write.csv(model_comp_results_woodland, file="../Results/Butterfly_results/model_comp_percentages_woodland.csv", row.names=FALSE)
# 
# ## now produce graph
# 
# library(ggplot2)
# 
# ## change levels of change so they go in correct order
# levels(model_comp_results_woodland$change)
# model_comp_results_woodland$change <- factor(model_comp_results_woodland$change, levels=c("Increase", "No change", "Decrease"))
# levels(model_comp_results_woodland$change)
# 
# ## pretty graph
# png("../Graphs/Model_comps/Model_comp_results_woodland.png", height = 100, width = 120, units = "mm", res = 300)
# ggplot(data=model_comp_results_woodland, aes(x=comparison, y=percentage, fill=change)) +
#   geom_bar(stat="identity", width=0.4) +
#   labs(y="Percentage of species", x="", fill="") +
#   scale_fill_manual(values = c("#339900", "#999999", "#990000")) +
#   theme_bw() +
#   theme(axis.text.x=element_text(colour = "black")) +
#   theme(axis.text.y=element_text(colour = "black")) +
#   scale_y_continuous(breaks = seq(0,100,10), expand = c(0, 0)) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
# dev.off()
