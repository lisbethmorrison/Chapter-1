#############################################################################
## Title: Abundance testing for woodland birds and butterflies
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: April 2018
#############################################################################

rm(list=ls()) # clear R
library(plyr)
library(dplyr) # load packages
library(lme4)
library(ggplot2)
library(lmerTest)
options(scipen=999)

# remove.packages(c("lme4", "data.table"))
# install.packages('Rcpp', dependencies = TRUE)
# install.packages('lme4', dependencies = TRUE)
# install.packages('data.table', dependencies = TRUE)


################################# AVERAGE ABUNDANCE IN FULL MODEL #####################################
 
## read data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE)
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC.csv", header=TRUE)
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE)
bird_common <- read.csv("../Data/BTO_data/pop_estimates_birds.csv", header=TRUE)
WCBS_data <- read.csv("../Data/UKBMS_data/WCBS_data.csv", header=TRUE)

######### UKBMS ###########
## remove unneccessary columns from WCBS data
WCBS_data <- WCBS_data[-c(3,5:7,9:11,13:15,17:19,21:23,25:27,29:31,33:35,37:38)] ## left with abundance from 09-17
## take the average abundance over time
WCBS_data$average_abundance <- rowMeans(subset(WCBS_data, select=c(2:10), na.rm=TRUE)) ## this is our measure of 'commonness'

## merge the two datasets
pair_attr <- merge(pair_attr, WCBS_data, by.x="spp", by.y="Species_code", all=FALSE)
length(unique(pair_attr$spp)) # 32 species (Grizzled Skipper doesn't have abundance data)
summary(pair_attr)

## run model
## scale average_abundance first (and other variables)
pair_attr$average_abundance <- (pair_attr$average_abundance - mean(na.omit(pair_attr$average_abundance)))/sd(na.omit(pair_attr$average_abundance))
pair_attr$distance <- (pair_attr$distance - mean(na.omit(pair_attr$distance)))/sd(na.omit(pair_attr$distance))
pair_attr$mean_northing <- (pair_attr$mean_northing - mean(na.omit(pair_attr$mean_northing)))/sd(na.omit(pair_attr$mean_northing))
pair_attr$renk_hab_sim <- (pair_attr$renk_hab_sim - mean(na.omit(pair_attr$renk_hab_sim)))/sd(na.omit(pair_attr$renk_hab_sim))

model_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), data = pair_attr)
pair_attr$residuals <- resid(model_ukbms)
pair_attr$fit <- fitted(model_ukbms)

qqnorm(resid(model_ukbms)) # A quantile normal plot - good for checking normality
qqline(resid(model_ukbms))

library(plotrix)
summary <- pair_attr %>% group_by(spp, average_abundance) %>% 
  summarise_at(vars(residuals), funs(mean,std.error))
# summary <- summary[!summary$spp == "18", ]
ggplot(summary, aes(x = average_abundance, y = mean)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error))

### full model with average_abundance as a measure of commonness
common_model_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + average_abundance + (1|pair.id) + (1|spp), data = pair_attr)
summary(common_model_ukbms)
anova(common_model_ukbms)
library(MuMIn)
r.squaredGLMM(common_model_ukbms)
## significant, p=0.039 (positive relationship)

## plot graph
newdata_ukbms <- expand.grid(mean_northing=mean(pair_attr$mean_northing), distance=mean(pair_attr$distance), 
                           renk_hab_sim=mean(pair_attr$renk_hab_sim), pair.id=sample(pair_attr$pair.id,10),
                           mid.year=mean(pair_attr$mid.year), spp=sample(pair_attr$spp,10),
                           average_abundance=unique(pair_attr$average_abundance))
newdata_ukbms$lag0 <- predict(common_model_ukbms, newdata=newdata_ukbms, re.form=NA)

mm <- model.matrix(terms(common_model_ukbms), newdata_ukbms)
pvar <- diag(mm %*% tcrossprod(vcov(common_model_ukbms),mm))
tvar <- pvar+VarCorr(common_model_ukbms)$spp[1]+VarCorr(common_model_ukbms)$pair.id[1]
cmult <- 2

newdata_ukbms <- data.frame(
  newdata_ukbms
  , plo = newdata_ukbms$lag0-1.96*sqrt(pvar)
  , phi = newdata_ukbms$lag0+1.96*sqrt(pvar)
  , tlo = newdata_ukbms$lag0-1.96*sqrt(tvar)
  , thi = newdata_ukbms$lag0+1.96*sqrt(tvar)
)

png("../Graphs/Abundance/Common_average_predicted_ukbms.png", height = 100, width = 120, units = "mm", res = 300)
ggplot(newdata_ukbms, aes(x=average_abundance, y=lag0)) +
  geom_line() +
  geom_ribbon(aes(ymin = plo, ymax = phi), alpha=0.2, lwd=0.5) +
  labs(x="Population abundance", y="Population synchrony") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()


######### CBC BIRDS #############
## merge the two datasets
pair_attr_CBC <- merge(pair_attr_CBC, bird_common, by.x="spp", by.y="species_code", all=FALSE)

## rescale variables
#pair_attr_CBC$pop_estimate <- (pair_attr_CBC$pop_estimate - mean(na.omit(pair_attr_CBC$pop_estimate)))/sd(na.omit(pair_attr_CBC$pop_estimate))
pair_attr_CBC$pop_estimate_log <- log10(pair_attr_CBC$pop_estimate)
pair_attr_CBC$distance <- (pair_attr_CBC$distance - mean(na.omit(pair_attr_CBC$distance)))/sd(na.omit(pair_attr_CBC$distance))
pair_attr_CBC$mean_northing <- (pair_attr_CBC$mean_northing - mean(na.omit(pair_attr_CBC$mean_northing)))/sd(na.omit(pair_attr_CBC$mean_northing))
## remove columns not needed (region, season, unit, year)
pair_attr_CBC <- pair_attr_CBC[-c(21:22,24:26)]
summary(pair_attr_CBC)

model_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + (1|pair.id) + (1|spp), data = pair_attr_CBC)
pair_attr_CBC$residuals <- resid(model_cbc)

library(plotrix)
summary <- pair_attr_CBC %>% group_by(spp, pop_estimate_log) %>% 
  summarise_at(vars(residuals), funs(mean,std.error))
# summary <- summary[!summary$spp == "18", ]
ggplot(summary, aes(x = pop_estimate_log, y = mean)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error))


### full model with pop_estimate as a measure of commonness
common_model_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + pop_estimate_log + (1|pair.id) + (1|spp), data = pair_attr_CBC)
summary(common_model_cbc)
anova(common_model_cbc)
## significant (p=0.0018)

## plot graph
newdata_cbc <- expand.grid(mean_northing=mean(pair_attr_CBC$mean_northing), distance=mean(pair_attr_CBC$distance), 
                       hab_sim=mean(pair_attr_CBC$hab_sim), pair.id=sample(pair_attr_CBC$pair.id,10),
                       mid.year=mean(pair_attr_CBC$mid.year), spp=sample(pair_attr_CBC$spp,10),
                       pop_estimate_log=unique(pair_attr_CBC$pop_estimate_log))
newdata_cbc$lag0 <- predict(common_model_cbc, newdata=newdata_cbc, re.form=NA)

mm2 <- model.matrix(terms(common_model_cbc), newdata_cbc)
pvar2 <- diag(mm2 %*% tcrossprod(vcov(common_model_cbc),mm2))
tvar2 <- pvar2+VarCorr(common_model_cbc)$spp[1]+VarCorr(common_model_cbc)$pair.id[1]
cmult <- 2

newdata_cbc <- data.frame(
  newdata_cbc
  , plo = newdata_cbc$lag0-1.96*sqrt(pvar2)
  , phi = newdata_cbc$lag0+1.96*sqrt(pvar2)
  , tlo = newdata_cbc$lag0-1.96*sqrt(tvar2)
  , thi = newdata_cbc$lag0+1.96*sqrt(tvar2)
)

png("../Graphs/Abundance/Common_average_predicted_cbc.png", height = 100, width = 120, units = "mm", res = 300)
ggplot(newdata_cbc, aes(x=pop_estimate_log, y=lag0)) +
  geom_line() +
  geom_ribbon(aes(ymin = plo, ymax = phi), alpha=0.2, lwd=0.5) +
  labs(x="Population Estimate", y="Population synchrony") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

# 
# ## split commonness into 3 groups and run model again 
# pair_attr_CBC$pop_estimate <- as.numeric(pair_attr_CBC$pop_estimate)
# pair_attr_CBC$common_score1 <- cut(pair_attr_CBC$pop_estimate, 3, labels=FALSE)
# pair_attr_CBC$common_score1 <- as.factor(pair_attr_CBC$common_score1)
# 
# ### full model with common_score1 of 3 levels
# common_model_cbc2 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + common_score1 + (1|pair.id) + (1|spp), data = pair_attr_CBC)
# summary(common_model_cbc2)
# anova(common_model_cbc2)
# ## significant p=0.017
# 
# ## run again with 2 groups 
# pair_attr_CBC$pop_estimate <- as.numeric(pair_attr_CBC$pop_estimate)
# pair_attr_CBC$common_score2 <- cut(pair_attr_CBC$pop_estimate, 2, labels=FALSE)
# pair_attr_CBC$common_score2 <- as.factor(pair_attr_CBC$common_score2)
# 
# ### full model with common_score1 of 2 levels
# common_model_cbc3 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + common_score2 + (1|pair.id) + (1|spp), data = pair_attr_CBC)
# summary(common_model_cbc3)
# anova(common_model_cbc3)
# ## still significant p=0.014
# 
# AIC(common_model_cbc,common_model_cbc2, common_model_cbc3) ## model 3 has lower AIC than model 2 (2 groups better than 3)
# ## model with raw common values still performs the best

######### BBS BIRDS #############
## merge the two datasets
pair_attr_BBS <- merge(pair_attr_BBS, bird_common, by.x="spp", by.y="species_code", all=FALSE)

## remove columns not needed (region, season, unit, year)
pair_attr_BBS <- pair_attr_BBS[-c(19:20,22:24)]
summary(pair_attr_BBS)

## rescale variables
pair_attr_BBS$pop_estimate <- (pair_attr_BBS$pop_estimate - mean(na.omit(pair_attr_BBS$pop_estimate)))/sd(na.omit(pair_attr_BBS$pop_estimate))
pair_attr_BBS$distance <- (pair_attr_BBS$distance - mean(na.omit(pair_attr_BBS$distance)))/sd(na.omit(pair_attr_BBS$distance))
pair_attr_BBS$mean_northing <- (pair_attr_BBS$mean_northing - mean(na.omit(pair_attr_BBS$mean_northing)))/sd(na.omit(pair_attr_BBS$mean_northing))
pair_attr_BBS$renk_hab_sim <- (pair_attr_BBS$renk_hab_sim - mean(na.omit(pair_attr_BBS$renk_hab_sim)))/sd(na.omit(pair_attr_BBS$renk_hab_sim))

### full model with pop_estimate as a measure of commonness
common_model_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + pop_estimate + (1|pair.id) + (1|spp), data = pair_attr_BBS)
anova(common_model_bbs)
summary(common_model_bbs)
## non-significant

### split commonness into 3 groups and run model again 
pair_attr_BBS$pop_estimate <- as.numeric(pair_attr_BBS$pop_estimate)
pair_attr_BBS$common_score1 <- cut(pair_attr_BBS$pop_estimate, 3, labels=FALSE)
pair_attr_BBS$common_score1 <- as.factor(pair_attr_BBS$common_score1)

### full model with common_score1 of 3 levels
common_model_bbs2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + common_score1 + (1|pair.id) + (1|spp), data = pair_attr_BBS)
summary(common_model_bbs2)
anova(common_model_bbs2)
## non-significant

## run again with 2 groups 
pair_attr_BBS$pop_estimate <- as.numeric(pair_attr_BBS$pop_estimate)
pair_attr_BBS$common_score2 <- cut(pair_attr_BBS$pop_estimate, 2, labels=FALSE)
pair_attr_BBS$common_score2 <- as.factor(pair_attr_BBS$common_score2)

### full model with common_score1 of 2 levels
common_model_bbs3 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + common_score2 + (1|pair.id) + (1|spp), data = pair_attr_BBS)
summary(common_model_bbs3)
anova(common_model_bbs3)
## non-significant

AIC(common_model_bbs2, common_model_bbs3) ## model 3 has lower AIC (2 groups better than 3)

#######################################################################################################
#################################### CALCULATE CHANGE IN ABUNDANCE #################################### 
#######################################################################################################

####################### ALL BUTTERFLIES #########################

### add abundance data
abundance_data <- read.csv("../Data/UKBMS_data/Collated_Indices_2016.csv", header=TRUE)
## add synchrony data
results_final_sp <- read.csv("../Results/Butterfly_results/results_final_sp.csv", header=TRUE)

## removed columns not needed (leaving species and habitat info)
results_final_sp2 <- results_final_sp[-c(2:8)]
results_final_sp2 <- unique(results_final_sp2)

## merge files together
abundance_data <- merge(abundance_data, results_final_sp2, by.x="Species.code", by.y="sp")
length(unique(abundance_data$Common.name)) ## 33 species
## 3 species missing from the abundance data - clouded yellow, red admiral and painted lady (all migrant species)
## remove some columns not needed
abundance_data <- abundance_data[-c(2,5,7,9:10)]
abundance_data <- droplevels(abundance_data)
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

################ Interaction between mid-year and abundance change ####################
## UKBMS butterflies
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
## not significant (p=0.64)

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
abund_model2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*ab_change_00_12 + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
summary(abund_model2)
anova(abund_model2)
## very significant interaction (p<0.00001)
r.squaredGLMM(abund_model2)

### plot this result using predicted values
newdata_ukbms <- expand.grid(mean_northing=mean(pair_attr_ukbms$mean_northing), distance=mean(pair_attr_ukbms$distance), 
                       renk_hab_sim=mean(pair_attr_ukbms$renk_hab_sim), pair.id=sample(pair_attr_ukbms$pair.id,10),
                       spp=sample(pair_attr_ukbms$spp,10), mid.year=unique(pair_attr_ukbms$mid.year), 
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

## change year values
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1999.5"="2000"))
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("2011.5"="2012"))
colnames(newdata_ukbms)[7] <- "Abundance_change"
newdata_ukbms$Abundance_change <- revalue(newdata_ukbms$Abundance_change, c("decrease"="Decrease"))
newdata_ukbms$Abundance_change <- revalue(newdata_ukbms$Abundance_change, c("increase"="Increase"))
levels(newdata_ukbms$Abundance_change)
newdata_ukbms$Abundance_change <- factor(newdata_ukbms$Abundance_change, levels=c("Increase", "Decrease"))
levels(newdata_ukbms$Abundance_change)

pd <- position_dodge(0.1)

png("../Graphs/Abundance/Abundance_change_predicted_ukbms_00_12.png", height = 100, width = 120, units = "mm", res = 300)
ggplot(newdata_ukbms, aes(x=mid.year, y=lag0, group=Abundance_change)) +
  geom_point(position=pd) +
  geom_line(aes(linetype=Abundance_change), size=1, position=pd) +
  labs(x="Year", y="Population synchrony") +
  geom_linerange(aes(ymin=plo, ymax=phi), position=pd) +
  labs(linetype="Change in abundance") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()


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

################################################################################
################################# CBC BIRDS ####################################
################################################################################

## for cbc only interested in years 1980-1989 and 1991-2000
bird_abundance_cbc_1 <- subset(bird_abundance[(bird_abundance$year>=1980) & (bird_abundance$year<=1989),])
bird_abundance_cbc_1$time <- "Early"
bird_abundance_cbc_2 <- subset(bird_abundance[(bird_abundance$year>=1991) & (bird_abundance$year<=2000),])
bird_abundance_cbc_2$time <- "Late"
bird_abundance_cbc <- rbind(bird_abundance_cbc_1, bird_abundance_cbc_2)

## now do t test on difference in abundance between early and late time periods
## calculates t test between early and late, produces p value, and produces mean of early and mean of late years
abundance_results_cbc = grouped_df(bird_abundance_cbc, vars="species_code") %>% summarise(p=t.test(sm[which(time=="Early")], sm[which(time=="Late")])$p.value, mean_x=mean(sm[which(time=="Early")]), mean_y=mean(sm[which(time=="Late")]))

## add in column to say whether species t test is significant or not
abundance_results_cbc$significance <- ifelse(abundance_results_cbc$p<0.05, "yes", "no")

## mean difference column of mean_y[LATE] - mean_x[EARLY]
## therefore negative difference == species has decrease in abundance between early and late years
abundance_results_cbc$mean_abund_change <- abundance_results_cbc$mean_y - abundance_results_cbc$mean_x
abundance_results_cbc$ab_change_85_96 <- ifelse(abundance_results_cbc$mean_abund_change>0, "increase", "decrease")                                                               

####### Interaction with mid-year and abundance change
pair_attr_cbc_1985 <- pair_attr_CBC[pair_attr_CBC$mid.year==1984.5,]
pair_attr_cbc_1996 <- pair_attr_CBC[pair_attr_CBC$mid.year==1995.5,]
pair_attr_cbc <- rbind(pair_attr_cbc_1985, pair_attr_cbc_1996)

pair_attr_cbc$mid.year <- as.factor(pair_attr_cbc$mid.year)
pair_attr_cbc$pair.id <- as.character(pair_attr_cbc$pair.id)
pair_attr_cbc$spp <- as.factor(pair_attr_cbc$spp)

## merge abundance data
pair_attr_cbc <- merge(pair_attr_cbc, abundance_results_cbc, by.x="spp", by.y="species_code", all=FALSE)
pair_attr_cbc <- pair_attr_cbc[-c(21:25)]

###### run model with strategy and year interaction
abund_model_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*ab_change_85_96 + (1|pair.id) + (1|spp), data = pair_attr_cbc)
summary(abund_model_cbc)
anova(abund_model_cbc)
## mid.year*ab_change_85_96 interaction is NON-significant overall

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
## mid.year*mean_abund_change interaction is very significant overall

## save model output
results_table_abund_bbs <- data.frame(summary(abund_model_bbs)$coefficients[,1:5])
write.csv(results_table_abund_bbs, file = "../Results/Model_outputs/change_abund_bbs.csv", row.names=TRUE)

### plot result using predicted values
newdata_bbs <- expand.grid(mean_northing=mean(pair_attr_bbs$mean_northing), distance=mean(pair_attr_bbs$distance), renk_hab_sim=mean(pair_attr_bbs$renk_hab_sim),
                       mid.year=unique(pair_attr_bbs$mid.year), ab_change_99_12=unique(pair_attr_bbs$ab_change_99_12),
                       pair.id=sample(pair_attr_bbs$pair.id,10), spp=sample(pair_attr_bbs$spp,10))
newdata_bbs$lag0 <- predict(abund_model_bbs, newdata=newdata_bbs, re.form=NA)

mm2 <- model.matrix(terms(abund_model_bbs), newdata_bbs)
pvar2 <- diag(mm2 %*% tcrossprod(vcov(abund_model_bbs),mm2))
tvar2 <- pvar2+VarCorr(abund_model_bbs)$spp[1]+VarCorr(abund_model_bbs)$pair.id[1]
cmult <- 2

newdata_bbs <- data.frame(
  newdata_bbs
  , plo = newdata_bbs$lag0-1.96*sqrt(pvar2)
  , phi = newdata_bbs$lag0+1.96*sqrt(pvar2)
  , tlo = newdata_bbs$lag0-1.96*sqrt(tvar2)
  , thi = newdata_bbs$lag0+1.96*sqrt(tvar2)
)

## change year values
newdata_bbs$mid.year <- revalue(newdata_bbs$mid.year, c("1998.5"="1999"))
newdata_bbs$mid.year <- revalue(newdata_bbs$mid.year, c("2011.5"="2012"))
colnames(newdata_bbs)[5] <- "Abundance_change"
newdata_bbs$Abundance_change <- revalue(newdata_bbs$Abundance_change, c("decrease"="Decrease"))
newdata_bbs$Abundance_change <- revalue(newdata_bbs$Abundance_change, c("increase"="Increase"))

pd <- position_dodge(0.1)

png("../Graphs/Abundance/Abundance_change_predicted_bbs.png", height = 100, width = 120, units = "mm", res = 300)
ggplot(newdata_bbs, aes(x=mid.year, y=lag0, group=Abundance_change)) +
  geom_point(position=pd) +
  geom_line(aes(linetype=Abundance_change), size=1, position=pd) +
  labs(x="Year", y="Population synchrony") +
  geom_linerange(aes(ymin=plo, ymax=phi), position=pd) +
  labs(linetype="Change in abundance") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

