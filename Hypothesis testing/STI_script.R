#############################################################################
## Title: STI and synchrony analysis
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: Feb 2019
#############################################################################

rm(list=ls()) # clear R
library(lme4)
library(lmerTest)
library(ggsignif)
library(ggplot2)
library(plyr)
library(dplyr)
library(plotrix)
library(parallel)
library(MuMIn)
options(scipen=999)

## read in synchrony data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr_no_zeros2.csv", header=TRUE) 
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC_no_zeros2_correct.csv", header=TRUE) 
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE) 
bird_STI <- read.csv("../Data/birds_STI.csv", header=TRUE)
butterfly_STI <- read.csv("../Data/butterflies_STI.csv", header=TRUE)
bird_dispersal <- read.csv("../Data/Woodland_bird_dispersal_Paradis1998.csv", header=TRUE)
## bird phylogeny info 
species_traits <- read.csv("../Data/BTO_data/woodland_generalist_specialist.csv", header=TRUE)

############### AVERAGE SYNCHORNY #################

#### UKBMS
## merge datasets
pair_attr <- merge(pair_attr, butterfly_STI, by.x="spp", by.y="species_code") ## 119 has no STI data
length(unique(pair_attr$spp)) # 31 species
summary(pair_attr)
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$mid.year <- as.factor(pair_attr$mid.year)

## run model with STI as fixed effect
climate_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + STI + (1|spp) + (1|pair.id), data=pair_attr)
summary(climate_ukbms) ## STI non-significant (p=0.38)
results_table_STI_ukbms <- data.frame(summary(climate_ukbms)$coefficients[,1:5]) ## 31 species
write.csv(results_table_STI_ukbms, file = "../Results/Model_outputs/UKBMS/average_STI_ukbms.csv", row.names=TRUE)

## run model with STI as two groups 
pair_attr$STI <- as.numeric(pair_attr$STI)
pair_attr$STI_score <- cut(pair_attr$STI, 2, labels=FALSE)
pair_attr$STI_score <- as.factor(pair_attr$STI_score)

climate_ukbms2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + STI_score + (1|spp) + (1|pair.id), data=pair_attr)
summary(climate_ukbms2) ## STI non-significant (p=0.155)
results_table_STI_ukbms <- data.frame(summary(climate_ukbms2)$coefficients[,1:5]) ## 31 species
write.csv(results_table_STI_ukbms, file = "../Results/Model_outputs/UKBMS/average_STIgroup_ukbms.csv", row.names=TRUE)

### plot average synchrony result anyway
## plot synchrony as a quadratic variable 
climate_ukbms_poly <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + poly(STI,2) + (1|spp) + (1|pair.id), data=pair_attr)
r.squaredGLMM(climate_ukbms_poly)
## predict new data to plot graph
newdata_ukbms <- expand.grid(mean_northing=mean(pair_attr$mean_northing), distance=mean(pair_attr$distance), 
                           renk_hab_sim=mean(pair_attr$renk_hab_sim), pair.id=sample(pair_attr$pair.id,1),
                           mid.year=mean(pair_attr$mid.year), spp=unique(pair_attr$spp),
                           STI=unique(pair_attr$STI))
newdata_ukbms$lag0 <- predict(climate_ukbms_poly, newdata=newdata_ukbms, re.form=NA)

## model without STI and species random effect to get residuals for graph
model_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data = pair_attr)
pair_attr$residuals <- resid(model_ukbms)

group_by(pair_attr, STI) %>% summarize(m = mean(lag0))
## put mean of each group into pair_attr dataframe
pair_attr <- ddply(pair_attr, "STI", transform, STI_mean = mean(lag0))
## add mean to each residual
pair_attr$residuals2 <- pair_attr$residuals + pair_attr$STI_mean

## create new dataframe to calculate mean, SD and SE of residuals for each species
summary_ukbms <- pair_attr %>% group_by(spp, STI) %>% 
  summarise_at(vars(residuals2), funs(mean,std.error,sd))

## plot graph with raw data residuals (+SE error bars) and polynomial line line
png("../Graphs/STI/STI_average_predicted_ukbms.png", height = 100, width = 110, units = "mm", res = 300)
ggplot(summary_ukbms, aes(x = STI, y = mean)) +
  geom_point(size = 2, colour="grey66") +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), width=0.2, colour="grey66") +
  geom_line(data=newdata_ukbms, aes(x=STI, y=lag0), colour="black", lwd=1) +
  labs(x="Species Temperature Index", y="Population synchrony") +
  theme_bw() +
  theme(text = element_text(size = 6), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

## test with mobility in the model too
## remove NA's (small white has no mobility data)
pair_attr <- na.omit(pair_attr)
## split mobility into two groups
pair_attr$mobility_wil <- as.numeric(pair_attr$mobility_wil)
pair_attr$mobility_score2 <- cut(pair_attr$mobility_wil, 2, labels=c("low", "high"))
pair_attr$mobility_score2 <- as.factor(pair_attr$mobility_score2)
length(unique(pair_attr$spp)) # 30 species

climate_mob_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + STI + mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr)
summary(climate_mob_ukbms) ## STI non-significant (p=0.24)
anova(climate_mob_ukbms) ## mobility is still significant (p=0.02)
results_table_STI_mob_ukbms <- data.frame(summary(climate_mob_ukbms)$coefficients[,1:5]) ## 30 species
write.csv(results_table_STI_mob_ukbms, file = "../Results/Model_outputs/UKBMS/average_STI_mob_ukbms.csv", row.names=TRUE)

##### CBC
# 
# ## merge in phylogeny info 
# species_traits <- species_traits[,c(2,3,4)]
# pair_attr_CBC <- merge(pair_attr_CBC, species_traits, by.x="spp", by.y="species_code", all=FALSE)
# pair_attr_CBC <- droplevels(pair_attr_CBC)
# length(unique(pair_attr_CBC$spp)) # 29 species

## merge datasets
pair_attr_CBC <- merge(pair_attr_CBC, bird_STI, by.x="spp", by.y="species_code") ## 31 species
length(unique(pair_attr_CBC$spp)) # 26 species
str(pair_attr_CBC)

pair_attr_CBC$pair.id <- as.character(pair_attr_CBC$pair.id)
pair_attr_CBC$spp <- as.factor(pair_attr_CBC$spp)
pair_attr_CBC$mid.year <- as.factor(pair_attr_CBC$mid.year)

climate_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + STI + (1|spp) + (1|pair.id), data=pair_attr_CBC)
summary(climate_cbc) ## STI non-significant (p=0.96) (still non-significant with family, p=0.124)
anova(climate_cbc)
results_table_STI_cbc <- data.frame(summary(climate_cbc)$coefficients[,1:5]) ## 26 species
write.csv(results_table_STI_cbc, file = "../Results/Model_outputs/CBC/average_STI_cbc_correct.csv", row.names=TRUE)

## same model with STI as group
pair_attr_CBC$STI <- as.numeric(pair_attr_CBC$STI)
pair_attr_CBC$STI_score <- cut(pair_attr_CBC$STI, 2, labels=FALSE)
pair_attr_CBC$STI_score <- as.factor(pair_attr_CBC$STI_score)

climate_cbc2 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + STI_score + (1|family) + (1|spp) + (1|pair.id), data=pair_attr_CBC)
## no errors
summary(climate_cbc2) ## STI non-significant (p=0.0534) 
anova(climate_cbc2)
results_table_STI_cbc <- data.frame(summary(climate_cbc2)$coefficients[,1:5]) ## 26 species
write.csv(results_table_STI_cbc, file = "../Results/Model_outputs/CBC/average_STIgroup_cbc_correct.csv", row.names=TRUE)

##### BBS
## merge datasets
pair_attr_BBS <- merge(pair_attr_BBS, bird_STI, by.x="spp", by.y="species_code") ## 24 species
length(unique(pair_attr_BBS$spp)) # 24 species

pair_attr_BBS$mid.year <- as.factor(pair_attr_BBS$mid.year)
pair_attr_BBS$pair.id <- as.character(pair_attr_BBS$pair.id)
pair_attr_BBS$spp <- as.factor(pair_attr_BBS$spp)

climate_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + STI + (1|spp) + (1|pair.id), data=pair_attr_BBS)
summary(climate_bbs) ## STI non-significant (p=0.099)
results_table_STI_bbs <- data.frame(summary(climate_bbs)$coefficients[,1:5])
write.csv(results_table_STI_bbs, file = "../Results/Model_outputs/BBS/average_STI_bbs.csv", row.names=TRUE)


## same model with STI as group
pair_attr_BBS$STI <- as.numeric(pair_attr_BBS$STI)
pair_attr_BBS$STI_score <- cut(pair_attr_BBS$STI, 2, labels=FALSE)
pair_attr_BBS$STI_score <- as.factor(pair_attr_BBS$STI_score)

climate_bbs2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + STI_score + (1|spp) + (1|pair.id), data=pair_attr_BBS)
summary(climate_bbs2) ## STI score is significant (p=0.034)
results_table_STI_bbs <- data.frame(summary(climate_bbs2)$coefficients[,1:5])
write.csv(results_table_STI_bbs, file = "../Results/Model_outputs/BBS/average_STIgroup_bbs.csv", row.names=TRUE)

## plot result
## predict new data to plot graph
pair_attr_BBS$mid.year <- as.factor(pair_attr_BBS$mid.year)
newdata_bbs <- expand.grid(mean_northing=mean(pair_attr_BBS$mean_northing), distance=mean(pair_attr_BBS$distance), 
                             renk_hab_sim=mean(pair_attr_BBS$renk_hab_sim), pair.id=sample(pair_attr_BBS$pair.id,1),
                             mid.year=unique(pair_attr_BBS$mid.year), spp=unique(pair_attr_BBS$spp),
                             STI_score=unique(pair_attr_BBS$STI_score))
newdata_bbs$lag0 <- predict(climate_bbs2, newdata=newdata_bbs, re.form=NA)
mm <- model.matrix(terms(climate_bbs2), newdata_bbs)
pvar <- diag(mm %*% tcrossprod(vcov(climate_bbs2),mm))
tvar <- pvar+VarCorr(climate_bbs2)$spp[1]+VarCorr(climate_bbs2)$pair.id[1]
cmult <- 2

newdata_bbs <- data.frame(
  newdata_bbs
  , plo = newdata_bbs$lag0-1.96*sqrt(pvar)
  , phi = newdata_bbs$lag0+1.96*sqrt(pvar)
  , tlo = newdata_bbs$lag0-1.96*sqrt(tvar)
  , thi = newdata_bbs$lag0+1.96*sqrt(tvar)
)

## model without STI and species random effect to get residuals for graph
model_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data = pair_attr_BBS)
pair_attr_BBS$residuals <- resid(model_bbs)

group_by(pair_attr_BBS, STI_score) %>% summarize(m = mean(lag0))
## put mean of each group into pair_attr dataframe
pair_attr_BBS <- ddply(pair_attr_BBS, "STI_score", transform, STI_mean = mean(lag0))
## add mean to each residual
pair_attr_BBS$residuals2 <- pair_attr_BBS$residuals + pair_attr_BBS$STI_mean

## create new dataframe to calculate mean, SD and SE of residuals for each species
summary_bbs <- pair_attr_BBS %>% group_by(spp, STI_score) %>% 
  summarise_at(vars(residuals2), funs(mean,std.error,sd))

## plot graph with raw data residuals (+SE error bars) and polynomial line line
png("../Graphs/STI/STI_average_predicted_ukbms.png", height = 100, width = 110, units = "mm", res = 300)
ggplot(summary_bbs, aes(x = STI_score, y = mean)) +
  geom_point(size = 2, colour="grey66") +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), width=0.2, colour="grey66") +
  geom_line(data=newdata_bbs, aes(x=STI_score, y=lag0), colour="black", lwd=1) +
  labs(x="Species Temperature Index", y="Population synchrony") +
  theme_bw() +
  theme(text = element_text(size = 6), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()


## does this change when mobility is added to the model? 
## This model doesn't run with STI score - the species with low STI is removed as it doesn't have mobility data

# pair_attr_BBS <- merge(pair_attr_BBS, bird_dispersal, by.x="spp", by.y="Species_code", all=FALSE)
# pair_attr_BBS <- na.omit(pair_attr_BBS)
# length(unique(pair_attr_BBS$spp)) # 17 species 
# ## split breeding_AM into 2 groups 
# pair_attr_BBS$Breeding_AM <- as.numeric(pair_attr_BBS$Breeding_AM)
# pair_attr_BBS$Breeding_AM_score2 <- cut(pair_attr_BBS$Breeding_AM, 2, labels=FALSE)
# pair_attr_BBS$Breeding_AM_score2 <- as.factor(pair_attr_BBS$Breeding_AM_score2)
# 
# climate_mob_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + STI_score + Breeding_AM_score2 + (1|spp) + (1|pair.id), data=pair_attr_BBS)
# summary(climate_mob_bbs) ## STI and mobility non-significant

############### CHANGE IN SYNCHORNY #################

## read in pair_attr file again to make sure all species are there
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr_no_zeros2.csv", header=TRUE) 
pair_attr <- merge(pair_attr, butterfly_STI, by.x="spp", by.y="species_code") ## 119 has no STI data
length(unique(pair_attr$spp)) # 31 species
summary(pair_attr)
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$mid.year <- as.factor(pair_attr$mid.year)

#### UKBMS
## subset to only look at mid years 1985, 2000 and 2012
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

## EARLY MODEL
climate_ukbms_change <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI + (1|spp) + (1|pair.id), data=pair_attr_early)
summary(climate_ukbms_change)
anova(climate_ukbms_change) ## interaction is very significant (p<0.0001)
results_table_STI_change_ukbms <- data.frame(summary(climate_ukbms_change)$coefficients[,1:5]) ## 31 species
write.csv(results_table_STI_change_ukbms, file = "../Results/Model_outputs/UKBMS/change_STI_ukbms_85_00.csv", row.names=TRUE)

## LATE MODEL
climate_ukbms_change2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI + (1|spp) + (1|pair.id), data=pair_attr_late)
summary(climate_ukbms_change2)
anova(climate_ukbms_change2) ## interaction is very significant (p<0.0001)
results_table_STI_change_ukbms <- data.frame(summary(climate_ukbms_change2)$coefficients[,1:5]) ## 31 species
write.csv(results_table_STI_change_ukbms, file = "../Results/Model_outputs/UKBMS/change_STI_ukbms_00_12.csv", row.names=TRUE)

## plot mid year and STI interaction result for EARLY MODEL
## split STI into three groups (low, med, high) to plot interaction result
pair_attr_early$STI <- as.numeric(pair_attr_early$STI)
pair_attr_early$STI_score <- cut(pair_attr_early$STI, 2, labels=FALSE)
pair_attr_early$STI_score <- as.factor(pair_attr_early$STI_score)

## run model
climate_ukbms_change3 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI_score + (1|pair.id) + (1|spp), data = pair_attr_early)
summary(climate_ukbms_change3)
anova(climate_ukbms_change3) ## significant (p<0.0001)
results_table_STI_change_ukbms <- data.frame(summary(climate_ukbms_change3)$coefficients[,1:5]) ## 31 species
write.csv(results_table_STI_change_ukbms, file = "../Results/Model_outputs/UKBMS/change_STIgroup_ukbms_85_00.csv", row.names=TRUE)

### predict new data
pair_attr_early$mid.year <- factor(pair_attr_early$mid.year)
newdata_ukbms <- expand.grid(mean_northing=mean(pair_attr_early$mean_northing), distance=mean(pair_attr_early$distance), 
                           renk_hab_sim=mean(pair_attr_early$renk_hab_sim), mid.year=unique(pair_attr_early$mid.year), 
                           pair.id=sample(pair_attr_early$pair.id,10), spp=(unique(pair_attr_early$spp)),
                           STI_score=unique(pair_attr_early$STI_score))
newdata_ukbms$lag0 <- predict(climate_ukbms_change3, newdata=newdata_ukbms, re.form=NA)

mm2 <- model.matrix(terms(climate_ukbms_change3), newdata_ukbms)
pvar2 <- diag(mm2 %*% tcrossprod(vcov(climate_ukbms_change3),mm2))
tvar2 <- pvar2+VarCorr(climate_ukbms_change3)$spp[1]+VarCorr(climate_ukbms_change3)$pair.id[1]
cmult <- 2

newdata_ukbms <- data.frame(
  newdata_ukbms
  , plo = newdata_ukbms$lag0-1.96*sqrt(pvar2)
  , phi = newdata_ukbms$lag0+1.96*sqrt(pvar2)
  , tlo = newdata_ukbms$lag0-1.96*sqrt(tvar2)
  , thi = newdata_ukbms$lag0+1.96*sqrt(tvar2)
)

## run model without STI and species random effect to obtain residuals
STI_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data=pair_attr_early)
pair_attr_early$residuals <- resid(STI_ukbms)
group_by(pair_attr_early, STI_score) %>% summarize(m = mean(lag0)) ## low mean = 0.0256, high mean = 0.0307
## put mean of each group into pair_attr dataframe
pair_attr_early <- ddply(pair_attr_early, "STI_score", transform, STI_mean = mean(lag0))
## add mean to each residual
pair_attr_early$residuals2 <- pair_attr_early$residuals + pair_attr_early$STI_mean

## create new dataframe to calculate mean, SD and SE of residuals for each species
summary_ukbms <- pair_attr_early %>% group_by(spp, STI_score, mid.year) %>% 
  summarise_at(vars(residuals2), funs(mean,std.error))
### 3 species have NA values for standard error (18, 27 and 116)
## these species only have one data point - cannot calculate standard error

## change values
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1984.5"="1985"))
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1999.5"="2000"))
colnames(newdata_ukbms)[7] <- "STI"
newdata_ukbms$STI <- as.factor(newdata_ukbms$STI)
newdata_ukbms$STI <- revalue(newdata_ukbms$STI, c("1"="Low"))
newdata_ukbms$STI <- revalue(newdata_ukbms$STI, c("2"="High"))
## same for summary dataframe
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("1984.5"="1985"))
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("1999.5"="2000"))
colnames(summary_ukbms)[2] <- "STI"
summary_ukbms$STI <- as.factor(summary_ukbms$STI)
summary_ukbms$STI <- revalue(summary_ukbms$STI, c("1"="Low"))
summary_ukbms$STI <- revalue(summary_ukbms$STI, c("2"="High"))
## reorder levels
levels(summary_ukbms$STI)
summary_ukbms$STI <- factor(summary_ukbms$STI, levels=c("High", "Low"))
levels(summary_ukbms$STI)
levels(newdata_ukbms$STI)
newdata_ukbms$STI <- factor(newdata_ukbms$STI, levels=c("High", "Low"))
levels(newdata_ukbms$STI)

png("../Graphs/STI/STI_change_predicted_ukbms_85_00.png", height = 120, width = 180, units = "mm", res = 300)
ggplot(summary_ukbms, aes(x = mid.year, y = mean, group=STI)) +
  geom_point(aes(shape=STI), colour="grey66", size = 3, position=myjit) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey66", width=0.2, position=myjit) +
  geom_line(data=newdata_ukbms, aes(x=mid.year, y=lag0, linetype=STI), lwd=1) +
  labs(x="Mid year of moving window", y="Population synchrony") +
  scale_y_continuous(breaks=seq(-2,0.8,0.1)) +
  theme_bw() +
  theme(legend.key.width = unit(0.8,"cm"), legend.key = element_rect(size = 2), text = element_text(size = 12), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.margin=margin(c(-10,20,-90,0)),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  scale_linetype_manual(name=" ",
                        labels=c("High", "Low"), values=c(1,2)) +
  scale_shape_manual(name="STI", 
                     labels=c("High", "Low"), values=c(16,4)) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(linetype = guide_legend(override.aes = list(size = 0.5)))
dev.off()

## plot mid year and STI interaction result for LATE MODEL
## split STI into three groups (low, med, high) to plot interaction result
pair_attr_late$STI <- as.numeric(pair_attr_late$STI)
pair_attr_late$STI_score <- cut(pair_attr_late$STI, 2, labels=FALSE)
pair_attr_late$STI_score <- as.factor(pair_attr_late$STI_score)

## run model
climate_ukbms_change4 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI_score + (1|pair.id) + (1|spp), data = pair_attr_late)
summary(climate_ukbms_change4)
anova(climate_ukbms_change4) ## significant (p<0.0001)
results_table_STI_change_ukbms <- data.frame(summary(climate_ukbms_change4)$coefficients[,1:5]) ## 31 species
write.csv(results_table_STI_change_ukbms, file = "../Results/Model_outputs/UKBMS/change_STIgroup_ukbms_00_12.csv", row.names=TRUE)

### predict new data
pair_attr_late$mid.year <- factor(pair_attr_late$mid.year)
newdata_ukbms <- expand.grid(mean_northing=mean(pair_attr_late$mean_northing), distance=mean(pair_attr_late$distance), 
                             renk_hab_sim=mean(pair_attr_late$renk_hab_sim), mid.year=unique(pair_attr_late$mid.year), 
                             pair.id=sample(pair_attr_late$pair.id,10), spp=(unique(pair_attr_late$spp)),
                             STI_score=unique(pair_attr_late$STI_score))
newdata_ukbms$lag0 <- predict(climate_ukbms_change4, newdata=newdata_ukbms, re.form=NA)

mm2 <- model.matrix(terms(climate_ukbms_change4), newdata_ukbms)
pvar2 <- diag(mm2 %*% tcrossprod(vcov(climate_ukbms_change4),mm2))
tvar2 <- pvar2+VarCorr(climate_ukbms_change4)$spp[1]+VarCorr(climate_ukbms_change4)$pair.id[1]
cmult <- 2

newdata_ukbms <- data.frame(
  newdata_ukbms
  , plo = newdata_ukbms$lag0-1.96*sqrt(pvar2)
  , phi = newdata_ukbms$lag0+1.96*sqrt(pvar2)
  , tlo = newdata_ukbms$lag0-1.96*sqrt(tvar2)
  , thi = newdata_ukbms$lag0+1.96*sqrt(tvar2)
)

## run model without STI and species random effect to obtain residuals
STI_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data=pair_attr_late)
pair_attr_late$residuals <- resid(STI_ukbms)
group_by(pair_attr_late, STI_score) %>% summarize(m = mean(lag0)) ## low mean = 0.0256, high mean = 0.0307
## put mean of each group into pair_attr dataframe
pair_attr_late <- ddply(pair_attr_late, "STI_score", transform, STI_mean = mean(lag0))
## add mean to each residual
pair_attr_late$residuals2 <- pair_attr_late$residuals + pair_attr_late$STI_mean

## create new dataframe to calculate mean, SD and SE of residuals for each species
summary_ukbms <- pair_attr_late %>% group_by(spp, STI_score, mid.year) %>% 
  summarise_at(vars(residuals2), funs(mean,std.error))
### 3 species have NA values for standard error (18, 27 and 116)
## these species only have one data point - cannot calculate standard error

## change values
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("2011.5"="2012"))
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1999.5"="2000"))
colnames(newdata_ukbms)[7] <- "STI"
newdata_ukbms$STI <- as.factor(newdata_ukbms$STI)
newdata_ukbms$STI <- revalue(newdata_ukbms$STI, c("1"="Low"))
newdata_ukbms$STI <- revalue(newdata_ukbms$STI, c("2"="High"))
## same for summary dataframe
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("2011.5"="2012"))
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("1999.5"="2000"))
colnames(summary_ukbms)[2] <- "STI"
summary_ukbms$STI <- as.factor(summary_ukbms$STI)
summary_ukbms$STI <- revalue(summary_ukbms$STI, c("1"="Low"))
summary_ukbms$STI <- revalue(summary_ukbms$STI, c("2"="High"))
## reorder levels
levels(summary_ukbms$STI)
summary_ukbms$STI <- factor(summary_ukbms$STI, levels=c("High", "Low"))
levels(summary_ukbms$STI)
levels(newdata_ukbms$STI)
newdata_ukbms$STI <- factor(newdata_ukbms$STI, levels=c("High", "Low"))
levels(newdata_ukbms$STI)

png("../Graphs/STI/STI_change_predicted_ukbms_00_12.png", height = 120, width = 180, units = "mm", res = 300)
ggplot(summary_ukbms, aes(x = mid.year, y = mean, group=STI)) +
  geom_point(aes(shape=STI), colour="grey66", size = 3, position=myjit) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey66", width=0.2, position=myjit) +
  geom_line(data=newdata_ukbms, aes(x=mid.year, y=lag0, linetype=STI), lwd=1) +
  labs(x="Mid year of moving window", y="Population synchrony") +
  scale_y_continuous(breaks=seq(-2,0.8,0.1)) +
  theme_bw() +
  theme(legend.key.width = unit(0.8,"cm"), legend.key = element_rect(size = 2), text = element_text(size = 12), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.margin=margin(c(-10,20,-90,0)),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  scale_linetype_manual(name=" ",
                        labels=c("High", "Low"), values=c(1,2)) +
  scale_shape_manual(name="STI", 
                     labels=c("High", "Low"), values=c(16,4)) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(linetype = guide_legend(override.aes = list(size = 0.5)))
dev.off()


myjit <- ggproto("fixJitter", PositionDodge,
                 width = 0.2,
                 dodge.width = 0.35,
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


############ Run models again with midyear*mobility interaction

## remove NAs from early and late datasets (as one species has no mobility data)
pair_attr_early <- na.omit(pair_attr_early)
length(unique(pair_attr_early$spp)) ## 30 species
pair_attr_late <- na.omit(pair_attr_late)
length(unique(pair_attr_late$spp)) ## 30 species

## split mobility into two groups
pair_attr_early$mobility_wil <- as.numeric(pair_attr_early$mobility_wil)
pair_attr_early$mobility_score2 <- cut(pair_attr_early$mobility_wil, 2, labels=c("low", "high"))
pair_attr_early$mobility_score2 <- as.factor(pair_attr_early$mobility_score2)

pair_attr_late$mobility_wil <- as.numeric(pair_attr_late$mobility_wil)
pair_attr_late$mobility_score2 <- cut(pair_attr_late$mobility_wil, 2, labels=c("low", "high"))
pair_attr_late$mobility_score2 <- as.factor(pair_attr_late$mobility_score2)

## EARLY MODEL
climate_ukbms_change_mob <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI_score + mid.year*mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr_early)
summary(climate_ukbms_change_mob) ## year*STI significant, year*mobility non-significant
anova(climate_ukbms_change_mob) 
results_table_STI_change_mob_ukbms <- data.frame(summary(climate_ukbms_change_mob)$coefficients[,1:5])
write.csv(results_table_STI_change_mob_ukbms, file = "../Results/Model_outputs/UKBMS/change_mob_STI_ukbms_85_00.csv", row.names=TRUE)
## early plot with mid.year*STI only is the MAM => so don't need to re-plot this model

## LATE MODEL
climate_ukbms_change_mob2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI_score + mid.year*mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr_late)
summary(climate_ukbms_change_mob2) ## both interactions are very significant
anova(climate_ukbms_change_mob2) 
results_table_STI_change_mob_ukbms <- data.frame(summary(climate_ukbms_change_mob2)$coefficients[,1:5])
write.csv(results_table_STI_change_mob_ukbms, file = "../Results/Model_outputs/UKBMS/change_mob_STI_ukbms_00_12.csv", row.names=TRUE)
## this it the MAM - need to re-plot this model

### predict new data
pair_attr_late$mid.year <- factor(pair_attr_late$mid.year)
newdata_ukbms <- expand.grid(mean_northing=mean(pair_attr_late$mean_northing), distance=mean(pair_attr_late$distance), 
                             renk_hab_sim=mean(pair_attr_late$renk_hab_sim), mid.year=unique(pair_attr_late$mid.year), 
                             pair.id=sample(pair_attr_late$pair.id,10), spp=(unique(pair_attr_late$spp)), 
                             mobility_score2=unique(pair_attr_late$mobility_score2), STI_score=unique(pair_attr_late$STI_score))
newdata_ukbms$lag0 <- predict(climate_ukbms_change_mob2, newdata=newdata_ukbms, re.form=NA)

mm2 <- model.matrix(terms(climate_ukbms_change_mob2), newdata_ukbms)
pvar2 <- diag(mm2 %*% tcrossprod(vcov(climate_ukbms_change_mob2),mm2))
tvar2 <- pvar2+VarCorr(climate_ukbms_change_mob2)$spp[1]+VarCorr(climate_ukbms_change_mob2)$pair.id[1]
cmult <- 2

newdata_ukbms <- data.frame(
  newdata_ukbms
  , plo = newdata_ukbms$lag0-1.96*sqrt(pvar2)
  , phi = newdata_ukbms$lag0+1.96*sqrt(pvar2)
  , tlo = newdata_ukbms$lag0-1.96*sqrt(tvar2)
  , thi = newdata_ukbms$lag0+1.96*sqrt(tvar2)
)

## run model without STI and species random effect to obtain residuals
STI_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data=pair_attr_late)
pair_attr_late$residuals <- resid(STI_ukbms)
group_by(pair_attr_late, STI_score) %>% summarize(m = mean(lag0)) ## low mean = 0.0256, high mean = 0.0307
## put mean of each group into pair_attr dataframe
pair_attr_late <- ddply(pair_attr_late, "STI_score", transform, STI_mean = mean(lag0))
## add mean to each residual
pair_attr_late$residuals2 <- pair_attr_late$residuals + pair_attr_late$STI_mean

## create new dataframe to calculate mean, SD and SE of residuals for each species
summary_ukbms <- pair_attr_late %>% group_by(spp, STI_score, mid.year) %>% 
  summarise_at(vars(residuals2), funs(mean,std.error))
### 3 species have NA values for standard error (18, 27 and 116)
## these species only have one data point - cannot calculate standard error

## change values
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1999.5"="2000"))
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("2011.5"="2012"))
colnames(newdata_ukbms)[8] <- "STI"
newdata_ukbms$STI <- as.factor(newdata_ukbms$STI)
newdata_ukbms$STI <- revalue(newdata_ukbms$STI, c("1"="Low"))
newdata_ukbms$STI <- revalue(newdata_ukbms$STI, c("2"="High"))
## same for summary dataframe
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("1999.5"="2000"))
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("2011.5"="2012"))
colnames(summary_ukbms)[2] <- "STI"
summary_ukbms$STI <- as.factor(summary_ukbms$STI)
summary_ukbms$STI <- revalue(summary_ukbms$STI, c("1"="Low"))
summary_ukbms$STI <- revalue(summary_ukbms$STI, c("2"="High"))
## reorder levels
levels(summary_ukbms$STI)
summary_ukbms$STI <- factor(summary_ukbms$STI, levels=c("High", "Low"))
levels(summary_ukbms$STI)
levels(newdata_ukbms$STI)
newdata_ukbms$STI <- factor(newdata_ukbms$STI, levels=c("High", "Low"))
levels(newdata_ukbms$STI)

png("../Graphs/STI/STI_mob_change_predicted_ukbms.png", height = 120, width = 180, units = "mm", res = 300)
ggplot(summary_ukbms, aes(x = mid.year, y = mean, group=STI)) +
  geom_point(aes(shape=STI), colour="grey66", size = 3, position=myjit) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey66", width=0.2, position=myjit) +
  geom_line(data=newdata_ukbms, aes(x=mid.year, y=lag0, linetype=STI), lwd=1) +
  labs(x="Mid year of moving window", y="Population synchrony") +
  scale_y_continuous(breaks=seq(-2,0.8,0.1)) +
  theme_bw() +
  theme(legend.key.width = unit(0.8,"cm"), legend.key = element_rect(size = 2), text = element_text(size = 12), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.margin=margin(c(-10,20,-90,0)),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  scale_linetype_manual(name=" ",
                        labels=c("High", "Low"), values=c(1,2)) +
  scale_shape_manual(name="STI", 
                     labels=c("High", "Low"), values=c(16,4)) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(linetype = guide_legend(override.aes = list(size = 0.5)))
dev.off()


myjit <- ggproto("fixJitter", PositionDodge,
                 width = 0.2,
                 dodge.width = 0.35,
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


#### CBC
## read in pair_attr_CBC again
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC_no_zeros2_correct.csv", header=TRUE) 
pair_attr_CBC <- merge(pair_attr_CBC, bird_STI, by.x="spp", by.y="species_code") ## 31 species
length(unique(pair_attr_CBC$spp)) # 26 species
pair_attr_CBC$pair.id <- as.character(pair_attr_CBC$pair.id)
pair_attr_CBC$spp <- as.factor(pair_attr_CBC$spp)
pair_attr_CBC$mid.year <- as.factor(pair_attr_CBC$mid.year)

pair_attr_cbc_1985 <- pair_attr_CBC[pair_attr_CBC$mid.year==1984.5,]
pair_attr_cbc_1996 <- pair_attr_CBC[pair_attr_CBC$mid.year==1995.5,]
pair_attr_cbc <- rbind(pair_attr_cbc_1985, pair_attr_cbc_1996)
length(unique(pair_attr_cbc$spp)) ## 26 species

pair_attr_cbc$mid.year <- as.factor(pair_attr_cbc$mid.year)
pair_attr_cbc$pair.id <- as.character(pair_attr_cbc$pair.id)
pair_attr_cbc$spp <- as.factor(pair_attr_cbc$spp)

###### run model with STI and year interaction
climate_cbc_change <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*STI + (1|pair.id) + (1|spp), data = pair_attr_cbc)
summary(climate_cbc_change)
anova(climate_cbc_change) ## interaction is non-significant (p=0.84) (still non-significant with family and genus added)
results_table_STI_change_cbc <- data.frame(summary(climate_cbc_change)$coefficients[,1:5]) ## 29 species
write.csv(results_table_STI_change_cbc, file = "../Results/Model_outputs/CBC/change_STI_cbc.csv", row.names=TRUE)

## same model with STI as group
pair_attr_cbc$STI <- as.numeric(pair_attr_cbc$STI)
pair_attr_cbc$STI_score <- cut(pair_attr_cbc$STI, 2, labels=FALSE)
pair_attr_cbc$STI_score <- as.factor(pair_attr_cbc$STI_score)

climate_cbc_change2 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*STI_score + (1|family) + (1|pair.id) + (1|spp), data = pair_attr_cbc)
## no errors
summary(climate_cbc_change2)
anova(climate_cbc_change2) ## interaction is non-significant (p=0.342) 
results_table_STI_change_cbc <- data.frame(summary(climate_cbc_change2)$coefficients[,1:5]) ## 26 species
write.csv(results_table_STI_change_cbc, file = "../Results/Model_outputs/CBC/change_STIgroup_cbc_correct.csv", row.names=TRUE)

## is mobility still significant with STI in the model?
pair_attr_cbc <- merge(pair_attr_cbc, bird_dispersal, by.x="spp", by.y="Species_code", all=FALSE)
pair_attr_cbc <- na.omit(pair_attr_cbc)
length(unique(pair_attr_cbc$spp)) # 21 species 
## split breeding_AM into 2 groups 
pair_attr_cbc$Breeding_AM <- as.numeric(pair_attr_cbc$Breeding_AM)
pair_attr_cbc$Breeding_AM_score2 <- cut(pair_attr_cbc$Breeding_AM, 2, labels=FALSE)
pair_attr_cbc$Breeding_AM_score2 <- as.factor(pair_attr_cbc$Breeding_AM_score2)

climate_mob_cbc_change <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*Breeding_AM_score2 +STI + (1|pair.id) + (1|spp), data = pair_attr_cbc)
summary(climate_mob_cbc_change)
anova(climate_mob_cbc_change) ## interaction and STI non-significant 
results_table_STI_mob_change_cbc <- data.frame(summary(climate_mob_cbc_change)$coefficients[,1:5]) ## 29 species
write.csv(results_table_STI_change_cbc, file = "../Results/Model_outputs/CBC/change_mob_STI_cbc.csv", row.names=TRUE)

## is STI still significant with mobility in the model?
climate_mob_cbc_change2 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*STI + Breeding_AM_score2 + (1|pair.id) + (1|spp), data = pair_attr_cbc)
summary(climate_mob_cbc_change2)
anova(climate_mob_cbc_change2) ## interaction and mobility non-significant 
results_table_STI_mob_change_cbc <- data.frame(summary(climate_mob_cbc_change2)$coefficients[,1:5]) ## 29 species
write.csv(results_table_STI_mob_change_cbc, file = "../Results/Model_outputs/CBC/change_STI_mob_cbc.csv", row.names=TRUE)

#### BBS
## read in pair_attr file again to make sure all species are there
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE) 
pair_attr_BBS <- merge(pair_attr_BBS, bird_STI, by.x="spp", by.y="species_code") ## 119 has no STI data
length(unique(pair_attr_BBS$spp)) # 24 species
pair_attr_BBS$mid.year <- as.factor(pair_attr_BBS$mid.year)
pair_attr_BBS$pair.id <- as.character(pair_attr_BBS$pair.id)
pair_attr_BBS$spp <- as.factor(pair_attr_BBS$spp)

pair_attr_bbs_1999 <- pair_attr_BBS[pair_attr_BBS$mid.year==1998.5,]
pair_attr_bbs_2012 <- pair_attr_BBS[pair_attr_BBS$mid.year==2011.5,]
pair_attr_bbs <- rbind(pair_attr_bbs_1999, pair_attr_bbs_2012)

pair_attr_bbs$mid.year <- as.factor(pair_attr_bbs$mid.year)
pair_attr_bbs$pair.id <- as.character(pair_attr_bbs$pair.id)
pair_attr_bbs$spp <- as.factor(pair_attr_bbs$spp)
length(unique(pair_attr_bbs$spp)) # 24 species

###### run model with STI and year interaction
climate_bbs_change <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI + (1|pair.id) + (1|spp), data = pair_attr_bbs)
summary(climate_bbs_change)
anova(climate_bbs_change) ## interaction is significant (p=0.0332)
results_table_STI_change_bbs <- data.frame(summary(climate_bbs_change)$coefficients[,1:5]) ## 24 species
write.csv(results_table_STI_change_bbs, file = "../Results/Model_outputs/BBS/change_STI_bbs.csv", row.names=TRUE)
#results_table_STI_change_bbs <- read.csv("../Results/Model_outputs/BBS/change_STI_bbs.csv", header=TRUE)

## split STI into two groups (low and high) to plot midyear*STI interaction result
pair_attr_bbs$STI <- as.numeric(pair_attr_bbs$STI)
pair_attr_bbs$STI_score <- cut(pair_attr_bbs$STI, 2, labels=FALSE)
pair_attr_bbs$STI_score <- as.factor(pair_attr_bbs$STI_score)

## run model
climate_bbs_change2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI_score + (1|pair.id) + (1|spp), data = pair_attr_bbs)
summary(climate_bbs_change2)
anova(climate_bbs_change2) ## interaction is not significant 
results_table_STI_change_bbs <- data.frame(summary(climate_bbs_change2)$coefficients[,1:5]) ## 24 species
write.csv(results_table_STI_change_bbs, file = "../Results/Model_outputs/BBS/change_STIgroup_bbs.csv", row.names=TRUE)

### predict new data
pair_attr_bbs$mid.year <- factor(pair_attr_bbs$mid.year)
newdata_bbs <- expand.grid(mean_northing=mean(pair_attr_bbs$mean_northing), distance=mean(pair_attr_bbs$distance), 
                           renk_hab_sim=mean(pair_attr_bbs$renk_hab_sim), mid.year=unique(pair_attr_bbs$mid.year), 
                           pair.id=sample(pair_attr_bbs$pair.id,10), spp=unique(pair_attr_bbs$spp),
                           STI_score=unique(pair_attr_bbs$STI_score))
newdata_bbs$lag0 <- predict(climate_bbs_change2, newdata=newdata_bbs, re.form=NA)

mm2 <- model.matrix(terms(climate_bbs_change2), newdata_bbs)
pvar2 <- diag(mm2 %*% tcrossprod(vcov(climate_bbs_change2),mm2))
tvar2 <- pvar2+VarCorr(climate_bbs_change2)$spp[1]+VarCorr(climate_bbs_change2)$pair.id[1]
cmult <- 2

newdata_bbs <- data.frame(
  newdata_bbs
  , plo = newdata_bbs$lag0-1.96*sqrt(pvar2)
  , phi = newdata_bbs$lag0+1.96*sqrt(pvar2)
  , tlo = newdata_bbs$lag0-1.96*sqrt(tvar2)
  , thi = newdata_bbs$lag0+1.96*sqrt(tvar2)
)

## run model without STI and species random effect to obtain residuals
STI_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data=pair_attr_bbs)
pair_attr_bbs$residuals <- resid(STI_bbs)
group_by(pair_attr_bbs, STI_score) %>% summarize(m = mean(lag0)) ## low mean = 0.0256, high mean = 0.0307
## put mean of each group into pair_attr dataframe
pair_attr_bbs <- ddply(pair_attr_bbs, "STI_score", transform, STI_mean = mean(lag0))
## add mean to each residual
pair_attr_bbs$residuals2 <- pair_attr_bbs$residuals + pair_attr_bbs$STI_mean

## create new dataframe to calculate mean, SD and SE of residuals for each species
summary_bbs <- pair_attr_bbs %>% group_by(spp, STI_score, mid.year) %>% 
  summarise_at(vars(residuals2), funs(mean,std.error))

## change values
newdata_bbs$mid.year <- revalue(newdata_bbs$mid.year, c("1998.5"="1999"))
newdata_bbs$mid.year <- revalue(newdata_bbs$mid.year, c("2011.5"="2012"))
colnames(newdata_bbs)[7] <- "STI"
newdata_bbs$STI <- as.factor(newdata_bbs$STI)
newdata_bbs$STI <- revalue(newdata_bbs$STI, c("1"="Low"))
newdata_bbs$STI <- revalue(newdata_bbs$STI, c("2"="High"))
## same for summary dataframe
summary_bbs$mid.year <- revalue(summary_bbs$mid.year, c("1998.5"="1999"))
summary_bbs$mid.year <- revalue(summary_bbs$mid.year, c("2011.5"="2012"))
colnames(summary_bbs)[2] <- "STI"
summary_bbs$STI <- as.factor(summary_bbs$STI)
summary_bbs$STI <- revalue(summary_bbs$STI, c("1"="Low"))
summary_bbs$STI <- revalue(summary_bbs$STI, c("2"="High"))
## reorder levels
levels(summary_bbs$STI)
summary_bbs$STI <- factor(summary_bbs$STI, levels=c("High", "Low"))
levels(summary_bbs$STI)
levels(newdata_bbs$STI)
newdata_bbs$STI <- factor(newdata_bbs$STI, levels=c("High", "Low"))
levels(newdata_bbs$STI)

pd <- position_dodge(0.1)
png("../Graphs/STI/STI_change_predicted_bbs.png", height = 150, width = 180, units = "mm", res = 300)
bbs_sti <- ggplot(summary_bbs, aes(x = mid.year, y = mean, group=STI)) +
  geom_point(aes(shape=STI), colour="grey66", size = 3, position=myjit) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey66", width=0.2, position=myjit) +
  geom_line(data=newdata_bbs, aes(x=mid.year, y=lag0, linetype=STI), lwd=1) +
  labs(x="Mid year of moving window", y="Population synchrony") +
  scale_y_continuous(breaks=seq(-2,0.6,0.1)) +
  theme_bw() +
  theme(legend.key.width = unit(0.8,"cm"), legend.key = element_rect(size = 2), text = element_text(size = 12), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.margin=margin(c(-10,20,-90,0)),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  scale_linetype_manual(name=" ",
                        labels=c("High", "Low"), values=c(1,2)) +
  scale_shape_manual(name="STI", 
                     labels=c("High", "Low"), values=c(16,4)) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(linetype = guide_legend(override.aes = list(size = 0.5)))
bbs_sti
dev.off()

## include STI graph for manuscript
library(ggpubr)
png("../Graphs/FINAL/Figure5_2.png", height = 200, width = 150, units = "mm", res = 300)
ggarrange(bbs_mob, bbs_sti, 
          labels = c("(a)", "(b)"), font.label = list(size = 10, color ="black"),
          ncol = 1, nrow = 2)
dev.off()

myjit <- ggproto("fixJitter", PositionDodge,
                 width = 0.2,
                 dodge.width = 0.35,
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


## is mobility still significant with STI in the model?
pair_attr_bbs <- merge(pair_attr_bbs, bird_dispersal, by.x="spp", by.y="Species_code", all=FALSE)
pair_attr_bbs <- na.omit(pair_attr_bbs)
length(unique(pair_attr_bbs$spp)) # 17 species

## split breeding_AM into 3 groups and run model again
pair_attr_bbs$Breeding_AM <- as.numeric(pair_attr_bbs$Breeding_AM)
pair_attr_bbs$Breeding_AM_score2 <- cut(pair_attr_bbs$Breeding_AM, 2, labels=FALSE)
pair_attr_bbs$Breeding_AM_score2 <- as.factor(pair_attr_bbs$Breeding_AM_score2)

climate_mob_change_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*Breeding_AM_score2 + STI + (1|pair.id) + (1|spp), data = pair_attr_bbs)
anova(climate_mob_change_bbs) ## midyear*mobility interaction is still significant (p=0.00001)
results_table_change_mob_STI <- data.frame(summary(climate_mob_change_bbs)$coefficients[,1:5])
write.csv(results_table_change_mob_STI, file = "../Results/Model_outputs/BBS/change_mob_STI_bbs.csv", row.names=TRUE) # 24 species

## is STI still significant with mobility in the model?
climate_mob_change_bbs2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI + Breeding_AM_score2 + (1|pair.id) + (1|spp), data = pair_attr_bbs)
anova(climate_mob_change_bbs2) ## midyear*STI interaction now not significant
## mobility is having an influence on STI changes over time
## maybe more mobile species have low STI?
results_table_change_mob_STI <- data.frame(summary(climate_mob_change_bbs2)$coefficients[,1:5])
write.csv(results_table_change_mob_STI, file = "../Results/Model_outputs/BBS/change_STI_mob_bbs.csv", row.names=TRUE) # 24 species
