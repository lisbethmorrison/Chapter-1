#############################################################################
## Title: Abundance testing for woodland birds and butterflies
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: April 2018
#############################################################################

rm(list=ls()) # clear R
library(plyr) # load packages
library(dplyr) 
library(lme4)
library(ggplot2)
library(lmerTest)
library(plotrix)
options(scipen=999)


################################# AVERAGE SYNCHRONY #####################################
 
###################################
#### Abundance model for UKBMS ####
###################################

## read data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr_no_zeros2.csv", header=TRUE)
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC_no_zeros2.csv", header=TRUE)
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE)
bird_common <- read.csv("../Data/BTO_data/pop_estimates_birds.csv", header=TRUE)
WCBS_data <- read.csv("../Data/UKBMS_data/WCBS_data.csv", header=TRUE)

######### UKBMS ###########
## remove unneccessary columns from WCBS data
WCBS_data <- WCBS_data[-c(3,5:7,9:11,13:15,17:19,21:23,25:27,29:31,33:35,37:38)] ## left with abundance from 09-17
## take the average abundance over time
WCBS_data$average_abundance <- rowMeans(subset(WCBS_data, select=c(2:10), na.rm=TRUE)) ## this is our measure of 'commonness'
## save dataframe
write.csv(WCBS_data, file="../Data/Butterfly_sync_data/Average_abundance_UKBMS.csv", row.names=FALSE)

## merge the two datasets
pair_attr <- merge(pair_attr, WCBS_data, by.x="spp", by.y="Species_code", all=FALSE)
length(unique(pair_attr$spp)) # 31 species (Grizzled Skipper doesn't have abundance data)
summary(pair_attr)

## run model
## scale variables
pair_attr$distance <- (pair_attr$distance - mean(na.omit(pair_attr$distance)))/sd(na.omit(pair_attr$distance))
pair_attr$mean_northing <- (pair_attr$mean_northing - mean(na.omit(pair_attr$mean_northing)))/sd(na.omit(pair_attr$mean_northing))
pair_attr$renk_hab_sim <- (pair_attr$renk_hab_sim - mean(na.omit(pair_attr$renk_hab_sim)))/sd(na.omit(pair_attr$renk_hab_sim))
pair_attr$pop_est_stand <- (pair_attr$average_abundance - mean(na.omit(pair_attr$average_abundance)))/sd(na.omit(pair_attr$average_abundance))

pair_attr$mid.year <- as.factor(pair_attr$mid.year)
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)

### full model with average_abundance as a measure of commonness
common_model_ukbms2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + pop_est_stand + (1|pair.id) + (1|spp), data = pair_attr)
summary(common_model_ukbms2)
anova(common_model_ukbms2)
## non-significant
results_table_abund_ukbms <- data.frame(summary(common_model_ukbms2)$coefficients[,1:5]) ## 31 species
write.csv(results_table_abund_ukbms, file = "../Results/Model_outputs/UKBMS/average_abund_ukbms.csv", row.names=TRUE)


# ## predict data to plot graphs
# newdata_ukbms <- expand.grid(mean_northing=mean(pair_attr$mean_northing), distance=mean(pair_attr$distance), 
#                              renk_hab_sim=mean(pair_attr$renk_hab_sim), pair.id=sample(pair_attr$pair.id,100),
#                              mid.year=mean(pair_attr$mid.year), spp=unique(pair_attr$spp),
#                              pop_est_stand=unique(pair_attr$pop_est_stand))
# newdata_ukbms$lag0 <- predict(common_model_ukbms2, newdata=newdata_ukbms, re.form=NA)
# 
# mm <- model.matrix(terms(common_model_ukbms2), newdata_ukbms)
# pvar <- diag(mm %*% tcrossprod(vcov(common_model_ukbms2),mm))
# tvar <- pvar+VarCorr(common_model_ukbms2)$spp[1]+VarCorr(common_model_ukbms2)$pair.id[1]
# cmult <- 2
# 
# newdata_ukbms <- data.frame(
#   newdata_ukbms
#   , plo = newdata_ukbms$lag0-1.96*sqrt(pvar)
#   , phi = newdata_ukbms$lag0+1.96*sqrt(pvar)
#   , tlo = newdata_ukbms$lag0-1.96*sqrt(tvar)
#   , thi = newdata_ukbms$lag0+1.96*sqrt(tvar)
# )
# 
# ## run model without abundance and species random effect to get residuals from graph
# model_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data = pair_attr)
# pair_attr$residuals <- resid(model_ukbms) ## save residuals
# 
# ## create new dataframe to calculate mean, SD and SE of residuals for each species
# summary_ukbms <- pair_attr %>% group_by(spp, average_abundance) %>% 
#   summarise_at(vars(residuals), funs(mean,std.error,sd))
# 
# ## plot graph with raw data residuals (+SE error bars) and fitted line
# png("../Graphs/Abundance/Common_average_predicted_ukbms.png", height = 80, width = 120, units = "mm", res = 300)
# ggplot(summary_ukbms, aes(x = average_abundance, y = mean)) +
#   geom_point(size = 1, colour="grey") +
#   geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey") +
#   geom_line(data=newdata_ukbms, aes(x=average_abundance, y=lag0), colour="black", lwd=0.5) +
#   labs(x="Average population abundance", y="Population synchrony") +
#   theme_bw() +
#   theme(text = element_text(size = 10), panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# dev.off()
# 

#################################
#### Abundance model for CBC ####
#################################

## merge the two datasets
pair_attr_CBC <- merge(pair_attr_CBC, bird_common, by.x="spp", by.y="species_code", all=FALSE)

## rescale variables
pair_attr_CBC$pop_estimate_log <- log(pair_attr_CBC$pop_estimate)
# pair_attr_CBC$pop_estimate_sqrt <- sqrt(pair_attr_CBC$pop_estimate)
# pair_attr_CBC$pop_estimate_standardise <- (pair_attr_CBC$pop_estimate - mean(na.omit(pair_attr_CBC$pop_estimate)))/sd(na.omit(pair_attr_CBC$pop_estimate))
# pair_attr_CBC$distance <- (pair_attr_CBC$distance - mean(na.omit(pair_attr_CBC$distance)))/sd(na.omit(pair_attr_CBC$distance))
# pair_attr_CBC$mean_northing <- (pair_attr_CBC$mean_northing - mean(na.omit(pair_attr_CBC$mean_northing)))/sd(na.omit(pair_attr_CBC$mean_northing))
# ## remove columns not needed (region, season, unit, year)
# pair_attr_CBC <- pair_attr_CBC[-c(21:22,24:26)]
# summary(pair_attr_CBC)

hist(pair_attr_CBC$pop_estimate)
hist(pair_attr_CBC$pop_estimate_standardise)
hist(pair_attr_CBC$pop_estimate_log)
hist(pair_attr_CBC$pop_estimate_sqrt)

### main model with pop_estimate as a measure of commonness
start_time <- Sys.time()
common_model_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + pop_estimate_log + (1|pair.id) + (1|spp), data = pair_attr_CBC)
end_time <- Sys.time()
end_time - start_time ## 26.14 seconds
summary(common_model_cbc) ## non-significant
anova(common_model_cbc) 
results_table_abund_cbc <- data.frame(summary(common_model_cbc)$coefficients[,1:5]) ## 29 species
write.csv(results_table_abund_cbc, file = "../Results/Model_outputs/CBC/average_abund_cbc.csv", row.names=TRUE)

## save main (true) model results
main_result_table <- data.frame(anova(common_model_cbc)[,5:6]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## remove rows with mean northing, distance, hab sim and mid year F values (only interested in abundance F values)
main_result_table <- main_result_table[ !(main_result_table$parameter %in% c("mean_northing", "distance", "hab_sim", "mid.year")), ]

## run model 999 times
perm_cbc_abund <- NULL
start_time <- Sys.time()
for (i in 1:999){
  print(i)
  pair_attr_CBC$abund_shuffle <- sample(pair_attr_CBC$pop_estimate_log) ## randomly shuffle pop_estimate variable
  common_model_cbc_shuffle <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + abund_shuffle + (1|pair.id) + (1|spp), data = pair_attr_CBC)
  ## run model with shuffled variable
  ## save results
  results_table_temp <- data.frame(anova(common_model_cbc_shuffle)[,5:6],i)
  perm_cbc_abund<-rbind(perm_cbc_abund,results_table_temp)
  }
end_time <- Sys.time()
end_time - start_time ## 6.3 hours (999 runs)

final_results_table <- rbind(main_result_table, perm_cbc_abund) ## bind the two data frames together

final_results_table$parameter <- paste(row.names(final_results_table)) ## move row.names to parameter column
rownames(final_results_table) <- 1:nrow(final_results_table) ## change row names to numbers

## remove rows with mean northing, distance, hab sim and mid year F values (only interested in abundance F values)
final_results_table <- final_results_table[!grep("mean_northing", final_results_table$parameter),]
final_results_table <- final_results_table[!grep("distance", final_results_table$parameter),]
final_results_table <- final_results_table[!grep("hab_sim", final_results_table$parameter),]
final_results_table <- final_results_table[!grep("mid.year", final_results_table$parameter),]

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/perm_average_abund_cbc.csv", row.names=TRUE)

## Calculate p value
number_of_permutations <- 1000
final_results_table2 <- final_results_table[!final_results_table$i==0,] ## remove true value to calc. p value
diff.observed <- main_result_table$F.value ## true F value
  
# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(final_results_table2$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0.005


# ## predict new data to plot graph
# newdata_cbc <- expand.grid(mean_northing=mean(pair_attr_CBC$mean_northing), distance=mean(pair_attr_CBC$distance), 
#                        hab_sim=mean(pair_attr_CBC$hab_sim), pair.id=sample(pair_attr_CBC$pair.id,10),
#                        mid.year=mean(pair_attr_CBC$mid.year), spp=unique(pair_attr_CBC$spp),
#                        pop_estimate_log=unique(pair_attr_CBC$pop_estimate_log))
# newdata_cbc$lag0 <- predict(common_model_cbc, newdata=newdata_cbc, re.form=NA)
# 
# mm2 <- model.matrix(terms(common_model_cbc), newdata_cbc)
# pvar2 <- diag(mm2 %*% tcrossprod(vcov(common_model_cbc),mm2))
# tvar2 <- pvar2+VarCorr(common_model_cbc)$spp[1]+VarCorr(common_model_cbc)$pair.id[1]
# cmult <- 2
# 
# newdata_cbc <- data.frame(
#   newdata_cbc
#   , plo = newdata_cbc$lag0-1.96*sqrt(pvar2)
#   , phi = newdata_cbc$lag0+1.96*sqrt(pvar2)
#   , tlo = newdata_cbc$lag0-1.96*sqrt(tvar2)
#   , thi = newdata_cbc$lag0+1.96*sqrt(tvar2)
# )
# 
# ## model without abundance and species random effect to get residuals for graph
# model_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + (1|pair.id), data = pair_attr_CBC)
# pair_attr_CBC$residuals <- resid(model_cbc)
# 
# group_by(pair_attr_CBC, pop_estimate_log) %>% summarize(m = mean(lag0)) ## low mean = 0.0256, high mean = 0.0307
# ## put mean of each group into pair_attr dataframe
# pair_attr_CBC <- ddply(pair_attr_CBC, "pop_estimate_log", transform, abund_mean = mean(lag0))
# ## add mean to each residual
# pair_attr_CBC$residuals2 <- pair_attr_CBC$residuals + pair_attr_CBC$abund_mean
# 
# ## create new dataframe to calculate mean, SD and SE of residuals for each species
# summary_cbc <- pair_attr_CBC %>% group_by(spp, pop_estimate_log) %>% 
#   summarise_at(vars(residuals2), funs(mean,std.error,sd))
# 
# ## plot graph with raw data residuals (+SE error bars) and fitted line
# png("../Graphs/Abundance/Common_average_predicted_cbc.png", height = 100, width = 110, units = "mm", res = 300)
# ggplot(summary_cbc, aes(x = pop_estimate_log, y = mean)) +
#   geom_point(size = 2, colour="grey66") +
#   geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), width=0.2, colour="grey66") +
#   geom_line(data=newdata_cbc, aes(x=pop_estimate_log, y=lag0), colour="black", lwd=1) +
#   labs(x="(log) Average population abundance", y="Population synchrony") +
#   theme_bw() +
#   theme(text = element_text(size = 8), panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# dev.off()
# 
# ## plot graph with fitted line and confidence interval shaded error
# png("../Graphs/Abundance/Common_average_predicted_cbc2.png", height = 100, width = 120, units = "mm", res = 300)
# ggplot(newdata_cbc, aes(x=pop_estimate, y=lag0)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = plo, ymax = phi), alpha=0.2, lwd=0.5) +
#   labs(x="Average population abundance", y="Population synchrony") +
#   theme_bw() +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# dev.off()

#################################
#### Abundance model for BBS ####
#################################

## merge the two datasets
pair_attr_BBS <- merge(pair_attr_BBS, bird_common, by.x="spp", by.y="species_code", all=FALSE)

## remove columns not needed (region, season, unit, year)
pair_attr_BBS <- pair_attr_BBS[-c(19:20,22:24)]
summary(pair_attr_BBS)

## rescale variables
pair_attr_BBS$pop_estimate_log <- log(pair_attr_BBS$pop_estimate)
pair_attr_BBS$distance <- (pair_attr_BBS$distance - mean(na.omit(pair_attr_BBS$distance)))/sd(na.omit(pair_attr_BBS$distance))
pair_attr_BBS$mean_northing <- (pair_attr_BBS$mean_northing - mean(na.omit(pair_attr_BBS$mean_northing)))/sd(na.omit(pair_attr_BBS$mean_northing))
pair_attr_BBS$renk_hab_sim <- (pair_attr_BBS$renk_hab_sim - mean(na.omit(pair_attr_BBS$renk_hab_sim)))/sd(na.omit(pair_attr_BBS$renk_hab_sim))

### full model with pop_estimate as a measure of commonness
common_model_bbs2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + pop_estimate_log + (1|pair.id) + (1|spp), data = pair_attr_BBS)
anova(common_model_bbs2)
summary(common_model_bbs2)
## non-significant

### split average abundance into two groups
pair_attr_BBS$pop_estimate <- as.numeric(pair_attr_BBS$pop_estimate)
pair_attr_BBS$rare_common <- cut(pair_attr_BBS$pop_estimate, 2, labels=c("rare", "common"))
pair_attr_BBS$rare_common <- as.factor(pair_attr_BBS$rare_common)

library(Hmisc) # cut2
pair_attr_BBS$rare_common2 <- cut2(pair_attr_BBS$pop_estimate, g=2)

## rare common by species
rare_common_bbs  <- pair_attr_BBS[c(17,19,20,21)]
rare_common_bbs <- unique(rare_common_bbs)
rare_common_bbs <- arrange(rare_common_bbs, as.numeric(pop_estimate))
rare_common_bbs$difference <- c(0, diff(rare_common_bbs$pop_estimate))

mean(rare_common_bbs$pop_estimate[rare_common_bbs$rare_common2=="[2500000,7700000]"]) # 5050000
mean(rare_common_bbs$pop_estimate[rare_common_bbs$rare_common2=="[   4300,2500000)"]) # 592461.1
mean(rare_common_bbs$pop_estimate[rare_common_bbs$rare_common=="common"]) # 6100000
mean(rare_common_bbs$pop_estimate[rare_common_bbs$rare_common=="rare"]) # 828215
## 5271785 difference for cut
## 4457538.9 difference for cut2
## cut has largest difference between means == better option

common_model_bbs2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + rare_common + (1|pair.id) + (1|spp), data = pair_attr_BBS)
summary(common_model_bbs2)
anova(common_model_bbs2)
### still not significant (p=0.425)

#######################################################################################################
####################################################################################################### 
#######################################################################################################

################################# CHANGE IN SYNCHRONY #####################################

###################################
#### Abundance model for UKBMS ####
###################################

### add abundance data
abundance_data <- read.csv("../Data/UKBMS_data/Collated_Indices_2016.csv", header=TRUE)
## add synchrony data
results_final_sp <- read.csv("../Results/Butterfly_results/results_final_sp_no_zeros2.csv", header=TRUE)

## removed columns not needed (leaving species and habitat info)
results_final_sp2 <- results_final_sp[-c(2:8)]
results_final_sp2 <- unique(results_final_sp2)

## merge files together
abundance_data <- merge(abundance_data, results_final_sp2, by.x="Species.code", by.y="sp")
length(unique(abundance_data$Common.name)) ## 32 species
## 3 species missing from the abundance data - clouded yellow, red admiral and painted lady (all migrant species)
## remove some columns not needed
abundance_data <- abundance_data[-c(2,5,7,9:10)]
abundance_data <- droplevels(abundance_data)
length(unique(abundance_data$Species.code))
####### do species significantly change in abundance between early and late years?

## subset two abundance data frames
## abundance t1 with years 1980-1989 (mid-year=1985)
## abundance t2 with years 1995-2004 (mid-year=2000)
abundance_t1 <- subset(abundance_data[(abundance_data$Year>=1980) & (abundance_data$Year<=1989),])
abundance_t1$Year <- "TP1"
abundance_t2 <- subset(abundance_data[(abundance_data$Year>=1995) & (abundance_data$Year<=2004),])
abundance_t2$Year <- "TP2"
## rbind together
abundance_final <- rbind(abundance_t1, abundance_t2)
str(abundance_final)
################## T TEST ON EACH SPECIES BETWEEN EARLY AND LATE YEARS ############################
## calculates t test between early and late, produces p value, and produces mean of early and mean of late years
abundance_results = grouped_df(abundance_final, vars="Species.code") %>% summarise(p=t.test(Collated.Index[which(Year=="TP1")], 
                   Collated.Index[which(Year=="TP2")])$p.value, mean_x=mean(Collated.Index[which(Year=="TP1")]), 
                   mean_y=mean(Collated.Index[which(Year=="TP2")]))

## add in column to say whether species t test is significant or not
abundance_results$significance <- ifelse(abundance_results$p<0.05, "yes", "no")

## mean difference column of mean_y[LATE] - mean_x[EARLY]
## therefore negative difference == species has decrease in abundance between early and late years
abundance_results$mean_change_abund <- abundance_results$mean_y - abundance_results$mean_x
abundance_results$ab_change_85_00 <- ifelse(abundance_results$mean_change_abund>0, "increase", "decrease")                                                               
## save abundance results
write.csv(abundance_results, file="../Results/Butterfly_results/abundance_results_85_00.csv", row.names=FALSE)
abundance_results <- read.csv("../Results/Butterfly_results/abundance_results_85_00.csv", header=TRUE)

################ Interaction between mid-year and abundance change ####################
## UKBMS butterflies
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr_no_zeros2.csv", header=TRUE) ## MUST READ IN DATA AGAIN
## subset to only look at mid years 1985 and 2000
pair_attr_1985 <- pair_attr[pair_attr$mid.year==1984.5,]
pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_ukbms <- rbind(pair_attr_1985, pair_attr_2000) ## 2 years

pair_attr_ukbms$mid.year <- as.factor(pair_attr_ukbms$mid.year)
pair_attr_ukbms$pair.id <- as.character(pair_attr_ukbms$pair.id)
pair_attr_ukbms$spp <- as.factor(pair_attr_ukbms$spp)

## rescale variables 
pair_attr_ukbms$distance <- (pair_attr_ukbms$distance - mean(na.omit(pair_attr_ukbms$distance)))/sd(na.omit(pair_attr_ukbms$distance))
pair_attr_ukbms$mean_northing <- (pair_attr_ukbms$mean_northing - mean(na.omit(pair_attr_ukbms$mean_northing)))/sd(na.omit(pair_attr_ukbms$mean_northing))
pair_attr_ukbms$renk_hab_sim <- (pair_attr_ukbms$renk_hab_sim - mean(na.omit(pair_attr_ukbms$renk_hab_sim)))/sd(na.omit(pair_attr_ukbms$renk_hab_sim))

## merge with abundance data 
pair_attr_ukbms <- merge(pair_attr_ukbms, abundance_results, by.x="spp", by.y="Species.code", all=FALSE)

## full model with interaction between mid year and mobility
abund_model1 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*ab_change_85_00 + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
summary(abund_model1)
anova(abund_model1)
## not significant (p=0.13)
results_table_abund_ukbms <- data.frame(summary(abund_model1)$coefficients[,1:5]) ## 32 species
write.csv(results_table_abund_ukbms, file = "../Results/Model_outputs/UKBMS/change_abund_85_00_ukbms.csv", row.names=TRUE)

################ same again but with abundance from 00-12 ################ 

## subset two abundance data frames
## abundance t3 with years 1995-2004 (mid-year=2000)
## abundance t4 with years 2007-2016 (mid-year=2012)
abundance_t3 <- subset(abundance_data[(abundance_data$Year>=1995) & (abundance_data$Year<=2004),])
abundance_t3$Year <- "TP3"
abundance_t4 <- subset(abundance_data[(abundance_data$Year>=2007) & (abundance_data$Year<=2016),])
abundance_t4$Year <- "TP4"
## rbind together
abundance_final <- rbind(abundance_t3, abundance_t4)
str(abundance_final)
################## T TEST ON EACH SPECIES BETWEEN EARLY AND LATE YEARS ############################
## calculates t test between early and late, produces p value, and produces mean of early and mean of late years
abundance_results = grouped_df(abundance_final, vars="Species.code") %>% summarise(p=t.test(Collated.Index[which(Year=="TP3")], 
                    Collated.Index[which(Year=="TP4")])$p.value, mean_x=mean(Collated.Index[which(Year=="TP3")]), 
                    mean_y=mean(Collated.Index[which(Year=="TP4")]))

## add in column to say whether species t test is significant or not
abundance_results$significance <- ifelse(abundance_results$p<0.05, "yes", "no")

## mean difference column of mean_y[LATE] - mean_x[EARLY]
## therefore negative difference == species has decrease in abundance between early and late years
abundance_results$mean_change_abund <- abundance_results$mean_y - abundance_results$mean_x
abundance_results$ab_change_00_12 <- ifelse(abundance_results$mean_change_abund>0, "increase", "decrease")                                                               
## save abundance results
write.csv(abundance_results, file="../Results/Butterfly_results/abundance_results_00_12.csv", row.names=FALSE)
abundance_results <- read.csv("../Results/Butterfly_results/abundance_results_00_12.csv", header=TRUE)
################ Interaction between mid-year and abundance change ####################
## UKBMS butterflies
## subset to only look at mid years 1985 and 2000
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

## full model with interaction between mid year and mobility
start_time <- Sys.time()
abund_model2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*ab_change_00_12 + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
end_time <- Sys.time()
end_time - start_time ## 22.23 seconds

summary(abund_model2)
anova(abund_model2)
## very significant interaction (p<0.00001)
results_table_abund_ukbms <- data.frame(summary(abund_model2)$coefficients[,1:5]) ## 32 species
write.csv(results_table_abund_ukbms, file = "../Results/Model_outputs/UKBMS/change_abund_00_12_ukbms.csv", row.names=TRUE)

### predict new data to plot graph
pair_attr_ukbms <- droplevels(pair_attr_ukbms)
newdata_ukbms <- expand.grid(mean_northing=mean(pair_attr_ukbms$mean_northing), distance=mean(pair_attr_ukbms$distance), 
                             renk_hab_sim=mean(pair_attr_ukbms$renk_hab_sim), pair.id=sample(pair_attr_ukbms$pair.id,100),
                             spp=unique(pair_attr_ukbms$spp), mid.year=unique(pair_attr_ukbms$mid.year), 
                             ab_change_00_12=unique(pair_attr_ukbms$ab_change_00_12))
newdata_ukbms$lag0 <- predict(abund_model2, newdata=newdata_ukbms, re.form=NA)


mm <- model.matrix(terms(abund_model2), newdata_ukbms)
pvar <- diag(mm %*% tcrossprod(vcov(abund_model2),mm))
tvar <- pvar+VarCorr(abund_model2)$spp[1]+VarCorr(abund_model2)$pair.id[1]
cmult <- 2

newdata_ukbms <- data.frame(
  newdata_ukbms
  , plo = newdata_ukbms$lag0-1.96*sqrt(pvar)
  , phi = newdata_ukbms$lag0+1.96*sqrt(pvar)
  , tlo = newdata_ukbms$lag0-1.96*sqrt(tvar)
  , thi = newdata_ukbms$lag0+1.96*sqrt(tvar)
)

## run model without abundance and species random effect to obtain residuals
abund_model_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data=pair_attr_ukbms)
pair_attr_ukbms$residuals <- resid(abund_model_ukbms)
group_by(pair_attr_ukbms, ab_change_00_12) %>% summarize(m = mean(lag0)) ## decrease=0.25, increase=0.286
## put mean of each group into pair_attr dataframe
pair_attr_ukbms <- ddply(pair_attr_ukbms, "ab_change_00_12", transform, abund_mean = mean(lag0))
## add mean to each residual
pair_attr_ukbms$residuals2 <- pair_attr_ukbms$residuals + pair_attr_ukbms$abund_mean

## create new dataframe to calculate mean, SD and SE of residuals for each species
summary_ukbms2 <- pair_attr_ukbms %>% group_by(spp, ab_change_00_12, mid.year) %>% 
  summarise_at(vars(residuals2), funs(mean,std.error,sd))

## change year values and abundance variable names
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1999.5"="2000"))
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("2011.5"="2012"))
colnames(newdata_ukbms)[7] <- "Abundance_change"
newdata_ukbms$Abundance_change <- revalue(newdata_ukbms$Abundance_change, c("decrease"="Decrease"))
newdata_ukbms$Abundance_change <- revalue(newdata_ukbms$Abundance_change, c("increase"="Increase"))
newdata_ukbms$Abundance_change <- factor(newdata_ukbms$Abundance_change, levels=c("Increase", "Decrease"))
## same for summary dataframe
summary_ukbms2$mid.year <- revalue(summary_ukbms2$mid.year, c("1999.5"="2000"))
summary_ukbms2$mid.year <- revalue(summary_ukbms2$mid.year, c("2011.5"="2012"))
colnames(summary_ukbms2)[2] <- "Abundance_change"
summary_ukbms2$Abundance_change <- revalue(summary_ukbms2$Abundance_change, c("decrease"="Decrease"))
summary_ukbms2$Abundance_change <- revalue(summary_ukbms2$Abundance_change, c("increase"="Increase"))
summary_ukbms2$Abundance_change <- factor(summary_ukbms2$Abundance_change, levels=c("Increase", "Decrease"))

## plot graph with raw data residuals (+SE error bars) and fitted line
png("../Graphs/Abundance/Abundance_change_predicted_ukbms_00_12.png", height = 150, width = 180, units = "mm", res = 300)
butterfly_abund <- ggplot(summary_ukbms2, aes(x = mid.year, y = mean, group=Abundance_change)) +
  geom_point(aes(shape=Abundance_change), colour="grey66", size = 3, position=myjit) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey66", position=myjit, width=0.2) +
  geom_line(data=newdata_ukbms, aes(x=mid.year, y=lag0, linetype=Abundance_change), colour="black", lwd=1) +
  labs(x="Mid year of moving window", y="Population synchrony") +
  theme_bw() +
  theme(legend.key.width = unit(0.8,"cm"), legend.key = element_rect(size = 2), text = element_text(size = 12), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.margin=margin(c(-10,20,-10,0)),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  scale_linetype_manual(name="Change in abundance",
                        labels=c("Increase", "Decrease"), values=c(1,2)) +
  scale_shape_manual(name="", 
                     labels=c("Increase", "Decrease"), values=c(16,4)) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(linetype = guide_legend(override.aes = list(size = 0.5)))

butterfly_abund

dev.off()

myjit <- ggproto("fixJitter", PositionDodge,
                 width = 0.5,
                 dodge.width = 0.15,
                 jit = NULL,
                 compute_panel =  function (self, data, params, scales) 
                 {
                   
                   #Generate Jitter if not yet
                   if(is.null(self$jit) ) {
                     self$jit <-jitter(rep(0, nrow(data)), amount=self$dodge.width)
                   }
                   
                   data <- ggproto_parent(PositionDodge, self)$compute_panel(data, params, scales)
                   
                   data$x <- data$x + self$jit
                   #For proper error extensions
                   if("xmin" %in% colnames(data)) data$xmin <- data$xmin + self$jit
                   if("xmax" %in% colnames(data)) data$xmax <- data$xmax + self$jit
                   data
                 } )


###################################################################################
###################################################################################
################################ WOODLAND BIRDS ###################################
###################################################################################
###################################################################################

## Woodland bird change in abundance
## add abundance files
file_names <- list.files("../Data/BTO_data/CBC_BBS_abundance", full.names=TRUE) # where all the files are stored
bird_abundance <- do.call(rbind,lapply(file_names,read.csv)) # create new dataframe with all species rbind together

## remove columns not needed (unsmoothed and CI columns)
bird_abundance <- bird_abundance[-c(2,4,5)]
## remove years before 1980
bird_abundance <- bird_abundance[bird_abundance$year>=1980,]

#################################
#### Abundance model for CBC ####
#################################

## for cbc only interested in years 1980-1989 and 1991-2000
bird_abundance_cbc_1 <- subset(bird_abundance[(bird_abundance$year>=1980) & (bird_abundance$year<=1989),])
bird_abundance_cbc_1$time <- "Early"
bird_abundance_cbc_2 <- subset(bird_abundance[(bird_abundance$year>=1991) & (bird_abundance$year<=2000),])
bird_abundance_cbc_2$time <- "Late"
bird_abundance_cbc <- rbind(bird_abundance_cbc_1, bird_abundance_cbc_2)

## now do t test on difference in abundance between early and late time periods
## calculates t test between early and late, produces p value, and produces mean of early and mean of late years
abundance_results_cbc = bird_abundance_cbc %>% group_by(species_code) %>% summarise(p=t.test(sm[which(time=="Early")], 
                        sm[which(time=="Late")])$p.value, mean_x=mean(sm[which(time=="Early")]), mean_y=mean(sm[which(time=="Late")]))

## add in column to say whether species t test is significant or not
abundance_results_cbc$significance <- ifelse(abundance_results_cbc$p<0.05, "yes", "no")

## mean difference column of mean_y[LATE] - mean_x[EARLY]
## therefore negative difference == species has decrease in abundance between early and late years
abundance_results_cbc$mean_abund_change <- abundance_results_cbc$mean_y - abundance_results_cbc$mean_x
abundance_results_cbc$ab_change_85_96 <- ifelse(abundance_results_cbc$mean_abund_change>0, "increase", "decrease")                                                               

## save abundance results
write.csv(abundance_results_cbc, file="../Results/Bird_results/abundance_results_85_96.csv", row.names=FALSE)
abundance_results_cbc <- read.csv("../Results/Bird_results/abundance_results_85_96.csv", header=TRUE)

####### Interaction with mid-year and abundance change
pair_attr_cbc_1985 <- pair_attr_CBC[pair_attr_CBC$mid.year==1984.5,]
pair_attr_cbc_1996 <- pair_attr_CBC[pair_attr_CBC$mid.year==1995.5,]
pair_attr_cbc <- rbind(pair_attr_cbc_1985, pair_attr_cbc_1996)

pair_attr_cbc$mid.year <- as.factor(pair_attr_cbc$mid.year)
pair_attr_cbc$pair.id <- as.character(pair_attr_cbc$pair.id)
pair_attr_cbc$spp <- as.factor(pair_attr_cbc$spp)

## merge abundance data
pair_attr_cbc <- merge(pair_attr_cbc, abundance_results_cbc, by.x="spp", by.y="species_code", all=FALSE)
length(unique(pair_attr_cbc$spp)) # 23 species 
pair_attr_cbc <- pair_attr_cbc[-c(21:25)]

###### run model with strategy and year interaction
start_time <- Sys.time()
abund_model_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*ab_change_85_96 + (1|pair.id) + (1|spp), data = pair_attr_cbc)
end_time <- Sys.time()
end_time - start_time ## 3.994 seconds

summary(abund_model_cbc)
anova(abund_model_cbc)
results_table_abund_cbc <- data.frame(summary(abund_model_cbc)$coefficients[,1:5]) ## 23 species
write.csv(results_table_abund_cbc, file = "../Results/Model_outputs/CBC/change_abund_cbc.csv", row.names=TRUE)

## save main (true) model results
main_result_table <- data.frame(anova(abund_model_cbc)[,5:6]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## remove rows with mean northing, distance, hab sim and mid year F values (only interested in abundance F values)
main_result_table <- main_result_table[ !(main_result_table$parameter %in% c("mean_northing", "distance", "hab_sim", "mid.year", "ab_change_85_96")), ]

## run model 999 times
perm_cbc_abund <- NULL
start_time <- Sys.time()
for (i in 1:999){
  print(i)
  pair_attr_cbc$abund_shuffle <- sample(pair_attr_cbc$ab_change_85_96) ## randomly shuffle abundance change varaible
  abund_model_cbc_shuffle <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*abund_shuffle + (1|pair.id) + (1|spp), data = pair_attr_cbc)
  ## run model with shuffled variable
  ## save results
  results_table_temp <- data.frame(anova(abund_model_cbc_shuffle)[,5:6],i)
  perm_cbc_abund<-rbind(perm_cbc_abund,results_table_temp)
}
end_time <- Sys.time()
end_time - start_time ## 59mins to do 999 runs

perm_cbc_abund$parameter <- paste(row.names(perm_cbc_abund)) ## move row.names to parameter column
rownames(perm_cbc_abund) <- 1:nrow(perm_cbc_abund) ## change row names to numbers
## remove rows with mean northing, distance, hab sim and mid year F values (only interested in abundance F values)
perm_cbc_abund <- perm_cbc_abund[-grep("mean_northing", perm_cbc_abund$parameter),]
perm_cbc_abund <- perm_cbc_abund[-grep("distance", perm_cbc_abund$parameter),]
perm_cbc_abund <- perm_cbc_abund[-grep("hab_sim", perm_cbc_abund$parameter),]
perm_cbc_abund <- perm_cbc_abund[-grep("mid.year", perm_cbc_abund$parameter),]

final_results_table <- rbind(main_result_table, perm_cbc_abund) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/perm_change_abund_cbc.csv", row.names=TRUE)

## Calculate p value
number_of_permutations <- 1000
final_results_table2 <- final_results_table[!final_results_table$i==0,] ## remove true value to calc. p value
diff.observed <- main_result_table$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(final_results_table2$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0.024

################################################################################
################################# BBS BIRDS ####################################
################################################################################

## for bbs only interested in years 1994-2003 and 2007-2016
## but only have data to 2011 for now
bird_abundance_bbs_1 <- subset(bird_abundance[(bird_abundance$year>=1994) & (bird_abundance$year<=2003),])
bird_abundance_bbs_1$time <- "Early"
bird_abundance_bbs_2 <- subset(bird_abundance[(bird_abundance$year>=2007) & (bird_abundance$year<=2016),]) 
bird_abundance_bbs_2$time <- "Late"
bird_abundance_bbs <- rbind(bird_abundance_bbs_1, bird_abundance_bbs_2)
## take out year 2001 => no data (foot&mouth)
bird_abundance_bbs<-bird_abundance_bbs[!(bird_abundance_bbs$year=="2001"),]

## now do t test on difference in abundance between early and late time periods
## calculates t test between early and late, produces p value, and produces mean of early and mean of late years
abundance_results_bbs = grouped_df(bird_abundance_bbs, vars="species_code") %>% summarise(p=t.test(sm[which(time=="Early")], sm[which(time=="Late")])$p.value, mean_x=mean(sm[which(time=="Early")]), mean_y=mean(sm[which(time=="Late")]))

## add in column to say whether species t test is significant or not
abundance_results_bbs$significance <- ifelse(abundance_results_bbs$p<0.05, "yes", "no")

## mean difference column of mean_y[LATE] - mean_x[EARLY]
## therefore negative difference == species has decrease in abundance between early and late years
abundance_results_bbs$mean_abund_change <- abundance_results_bbs$mean_y - abundance_results_bbs$mean_x
abundance_results_bbs$ab_change_99_12 <- ifelse(abundance_results_bbs$mean_abund_change>0, "increase", "decrease")                                                               

#### BBS birds
pair_attr_bbs_1999 <- pair_attr_BBS[pair_attr_BBS$mid.year==1998.5,]
pair_attr_bbs_2012 <- pair_attr_BBS[pair_attr_BBS$mid.year==2011.5,]
pair_attr_bbs <- rbind(pair_attr_bbs_1999, pair_attr_bbs_2012)

pair_attr_bbs$mid.year <- as.factor(pair_attr_bbs$mid.year)
pair_attr_bbs$pair.id <- as.character(pair_attr_bbs$pair.id)
pair_attr_bbs$spp <- as.factor(pair_attr_bbs$spp)

## merge with abundance data
pair_attr_bbs <- merge(pair_attr_bbs, abundance_results_bbs, by.x="spp", by.y="species_code", all=FALSE)
pair_attr_bbs <- pair_attr_bbs[-c(19:23)]

###### run model with strategy and year interaction
abund_model_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*ab_change_99_12 + (1|pair.id) + (1|spp), data = pair_attr_bbs)
summary(abund_model_bbs)
## intercept is mid.year 1999
anova(abund_model_bbs)
## mid.year*mean_abund_change interaction is non-significant 
