#############################################################################
## Title: Mobility and average synchrony for butterflies and woodland birds
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: April 2018
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
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE) 
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC.csv", header=TRUE) 
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE) 
wingspan <- read.csv("../Data/UKBMS_data/butterflywingspans.csv", header=TRUE)
bird_dispersal <- read.csv("../Data/Woodland_bird_dispersal_Paradis1998.csv", header=TRUE)

################################### AVERAGE SYNCRHONY ################################### 

##############################testing wingspand and synchrony first##############################
# ## merge pair attr and wingspan data
# pair_attr_1 <- merge(pair_attr, wingspan, by.x="spp", by.y="BMScode")
# 
# pair_attr_1$mid.year <- as.factor(pair_attr_1$mid.year)
# pair_attr_1$pair.id <- as.character(pair_attr_1$pair.id)
# pair_attr_1$spp <- as.factor(pair_attr_1$spp)
# 
# ## change wingspan into 3 groups (small, med, large)
# pair_attr_1$Wingspan <- as.numeric(pair_attr_1$Wingspan)
# pair_attr_1$Wingspan_group <- cut(pair_attr_1$Wingspan, 3, labels=FALSE)
# pair_attr_1$Wingspan_group <- as.factor(pair_attr_1$Wingspan_group)
# 
# wing_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Wingspan + (1|spp) + (1|pair.id), data=pair_attr_1)
# summary(wing_model)
# anova(wing_model)
# ## wingspan is non-significant (p=0.33)
# 
# wing_model2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Wingspan + Family + (1|spp) + (1|pair.id), data=pair_attr_1)
# summary(wing_model2)
# anova(wing_model2)
# ## neither family nor wingspan are significant 
# 
# wing_model3 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Wingspan + Sub.family + (1|spp) + (1|pair.id), data=pair_attr_1)
# summary(wing_model3)
# anova(wing_model3)
# ## neither sub-family nor wingspan are significant
# 
# wing_model4 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (Wingspan|Sub.family) + (1|spp) + (1|pair.id), data=pair_attr_1)
# anova(wing_model4)
# ## doesn't work ==> model failed to converge
# 
# wing_model5 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Wingspan + (1|Sub.family) + (1|spp) + (1|pair.id), data=pair_attr_1)
# anova(wing_model5)
# ## wingspan still not significant 
# 
# wing_model6 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Wingspan*Sub.family + (1|spp) + (1|pair.id), data=pair_attr_1)
# anova(wing_model6)
# ## interaction is significant
# ## but warning: fixed-effect model matrix is rank deficient so dropping 3 columns / coefficients
# ## this is because 3 sub-families only have one wingspan, so can't model this
# ## Coliadinae, Limenitidinae, and Lycaeninae
# ## remove these sub-families and run model again
# pair_attr_1 <- pair_attr_1[!pair_attr_1$Sub.family=="Coliadinae",]
# pair_attr_1 <- pair_attr_1[!pair_attr_1$Sub.family=="Limenitidinae",]
# pair_attr_1 <- pair_attr_1[!pair_attr_1$Sub.family=="Lycaeninae",]
# ## run model again
# wing_model7 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Wingspan*Sub.family + (1|spp) + (1|pair.id), data=pair_attr_1)
# anova(wing_model7)
# summary(wing_model7)
# 
# ## plot result
# newdata_wing <- expand.grid(mean_northing=mean(pair_attr_1$mean_northing), distance=mean(pair_attr_1$distance), 
#                              renk_hab_sim=mean(pair_attr_1$renk_hab_sim), pair.id=sample(pair_attr_1$pair.id,100),
#                              spp=unique(pair_attr_1$spp), mid.year=unique(pair_attr_1$mid.year), 
#                              Wingspan=unique(pair_attr_1$Wingspan), Sub.family=unique(pair_attr_1$Sub.family))
# 
# newdata_wing$lag0 <- predict(wing_model7, newdata=newdata_wing, re.form=NA)
# 
# mm <- model.matrix(terms(wing_model7), newdata_wing)
# pvar <- diag(mm %*% tcrossprod(vcov(wing_model7),mm))
# tvar <- pvar+VarCorr(wing_model7)$spp[1]+VarCorr(wing_model7)$pair.id[1]
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
# ## run model without mobility or species random effect to obtain residuals
# wing_mod <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data=pair_attr_1)
# pair_attr_1$residuals <- resid(wing_mod, type="pearson")
# ## check mean synchrony of each group
# group_by(pair_attr_1, Wingspan) %>% summarize(m = mean(lag0)) ## WC mean = 0.265, HS mean = 0.195
# ## put mean of each group into pair_attr dataframe
# pair_attr_1 <- ddply(pair_attr_1, "Wingspan", transform, wing_mean = mean(lag0))
# ## add mean to each residual
# pair_attr_1$residuals2 <- pair_attr_1$residuals + pair_attr_1$wing_mean
# 
# ## create dataframe which calculates mean, SD and SE of residuals for each species
# summary_wing2 <- pair_attr_1 %>% group_by(spp, Wingspan, Sub.family) %>% 
#   summarise_at(vars(residuals2), funs(mean,std.error))
# 
# summary_wing2 <- summary_wing2[!summary_wing2$Sub.family=="Coliadinae",]
# summary_wing2 <- summary_wing2[!summary_wing2$Sub.family=="Limenitidinae",]
# summary_wing2 <- summary_wing2[!summary_wing2$Sub.family=="Lycaeninae",]
# 
# ggplot(summary_wing2, aes(x = Wingspan, y = mean, colour=Sub.family)) +
#   geom_point(aes(colour=Sub.family), size = 2) +
#   #geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error, colour=Sub.family), width=0.2) +
#   geom_line(data=newdata_wing, aes(x=Wingspan, y=lag0, colour=Sub.family), lwd=1) +
#   theme_bw() +
#   theme(text = element_text(size = 12), panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 
# ## AIC of all models - which wingspan/sub-family model performs best?
# AIC(wing_model3, wing_model4, wing_model5, wing_model6)
# ## for now wing_model5 has the lowest AIC

##################################
#### Mobility model for UKBMS ####
##################################

## remove NA's ==> removes Small White butterfly which doesn't have mobility data
pair_attr <- na.omit(pair_attr)
length(unique(pair_attr$spp)) # 31 species

pair_attr$mid.year <- as.factor(pair_attr$mid.year)
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)

### model with mobility score of 3 levels (1=low, 2=medium, 3=high)
pair_attr$mobility_wil <- as.numeric(pair_attr$mobility_wil)
pair_attr$mobility_score <- cut(pair_attr$mobility_wil, 3, labels=FALSE)
pair_attr$mobility_score <- as.factor(pair_attr$mobility_score)

mobility_model2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + mobility_score + (1|spp) + (1|pair.id), data=pair_attr)
summary(mobility_model2)
anova(mobility_model2)
### significant - higher mobility score == higher average synchrony
## 3 is significantly different from 1 (p=0.000684)
## 2 is not significantly different from 1 (p=0.133)

## run again with score=2 as intercept
levels(pair_attr$mobility_score)
pair_attr$mobility_score <- factor(pair_attr$mobility_score, levels=c("2", "3", "1"))
levels(pair_attr$mobility_score)

mobility_model3 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + mobility_score + (1|spp) + (1|pair.id), data=pair_attr)
summary(mobility_model3)
anova(mobility_model3)
## mobility still significant 
## 3 is significantly different from 2 (p=0.0137)
## 1 still not significantly different from 2 (p=0.133)

## run model with 2 groups (high and low) of mobility 
pair_attr$mobility_wil <- as.numeric(pair_attr$mobility_wil)
pair_attr$mobility_score2 <- cut(pair_attr$mobility_wil, 2, labels=c("low", "high"))
pair_attr$mobility_score2 <- as.factor(pair_attr$mobility_score2)

mobility_model4 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr)
summary(mobility_model4)
anova(mobility_model4)
## significant (p=0.011)
## mobility score 2 has higher average synchrony than score 1
r.squaredGLMM(mobility_model4)
## R^2 marginal (fixed effects) = 0.029, R^2 conditional (random effects) = 0.224
results_table_mob_ukbms <- data.frame(summary(mobility_model4)$coefficients[,1:5]) ## 31 species
write.csv(results_table_mob_ukbms, file = "../Results/Model_outputs/UKBMS/average_mob_ukbms.csv", row.names=TRUE)

##################################
#### Mobility model for CBC ####
##################################

## merge with bird dispersal data
length(unique(bird_dispersal$Species_code)) ## only 23 species have dispersal data

## merge pair_attr with bird dispersal data
pair_attr_CBC <- merge(pair_attr_CBC, bird_dispersal, by.x="spp", by.y="Species_code", all=FALSE)
pair_attr_CBC <- na.omit(pair_attr_CBC) ## remove NAs (one species == redstart)
length(unique(pair_attr_CBC$spp)) ## 20 species

## breeding arithmetic mean dispersal
dispersal_model_cbc2 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + Breeding_AM + (1|spp) + (1|pair.id), data=pair_attr_CBC)
summary(dispersal_model_cbc2)
anova(dispersal_model_cbc2)
## not significant 

## split breeding_AM into 3 groups and run model again
pair_attr_CBC$Breeding_AM <- as.numeric(pair_attr_CBC$Breeding_AM)
pair_attr_CBC$Breeding_AM_score <- cut(pair_attr_CBC$Breeding_AM, 3, labels=FALSE)
pair_attr_CBC$Breeding_AM_score <- as.factor(pair_attr_CBC$Breeding_AM_score)

dispersal_model_cbc3 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + Breeding_AM_score + (1|spp) + (1|pair.id), data=pair_attr_CBC)
summary(dispersal_model_cbc3)
anova(dispersal_model_cbc3)
## still not significant

## split breeding_AM into 2 groups and compare
pair_attr_CBC$Breeding_AM <- as.numeric(pair_attr_CBC$Breeding_AM)
pair_attr_CBC$Breeding_AM_score2 <- cut(pair_attr_CBC$Breeding_AM, 2, labels=c("low", "high"))
pair_attr_CBC$Breeding_AM_score2 <- as.factor(pair_attr_CBC$Breeding_AM_score2)

pair_attr_CBC <- droplevels(pair_attr_CBC)

dispersal_model_cbc4 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + Breeding_AM_score2 + (1|family) + (1|spp) + (1|pair.id), data=pair_attr_CBC)
## no errors
summary(dispersal_model_cbc4)
anova(dispersal_model_cbc4) ## dispersal score is non-significant (p=0.63)
results_table_mob_cbc <- data.frame(summary(dispersal_model_cbc4)$coefficients[,1:5]) ## 20 species
write.csv(results_table_mob_cbc, file = "../Results/Model_outputs/CBC/average_mob_cbc.csv", row.names=TRUE)

################################
#### Mobility model for BBS ####
################################

## merge pair_attr with bird dispersal data
pair_attr_BBS <- merge(pair_attr_BBS, bird_dispersal, by.x="spp", by.y="Species_code", all=FALSE)
length(unique(pair_attr_BBS$spp)) # 18 species until remove NAs for Breeding dispersal (==17)

## natal arithmetic mean dispersal
dispersal_model_bbs1 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Natal_AM + (1|spp) + (1|pair.id), data=pair_attr_BBS)
summary(dispersal_model_bbs1)
anova(dispersal_model_bbs1)
## not significant

## breeding arithmetic mean dispersal
## remove NAs (one species == redstart)
pair_attr_BBS <- na.omit(pair_attr_BBS)
summary(pair_attr_BBS)  

dispersal_model_bbs2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Breeding_AM + (1|spp) + (1|pair.id), data=pair_attr_BBS)
summary(dispersal_model_bbs2)
anova(dispersal_model_bbs2)
## not significant 

## split breeding_AM into 3 groups and run model again
pair_attr_BBS$Breeding_AM <- as.numeric(pair_attr_BBS$Breeding_AM)
pair_attr_BBS$Breeding_AM_score <- cut(pair_attr_BBS$Breeding_AM, 3, labels=FALSE)
pair_attr_BBS$Breeding_AM_score <- as.factor(pair_attr_BBS$Breeding_AM_score)

dispersal_model_bbs3 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Breeding_AM_score + (1|spp) + (1|pair.id), data=pair_attr_BBS)
summary(dispersal_model_bbs3)
anova(dispersal_model_bbs3)
## still not significant

## split breeding_AM into 2 groups and compare
pair_attr_BBS$Breeding_AM <- as.numeric(pair_attr_BBS$Breeding_AM)
pair_attr_BBS$Breeding_AM_score2 <- cut(pair_attr_BBS$Breeding_AM, 2, labels=FALSE)
pair_attr_BBS$Breeding_AM_score2 <- as.factor(pair_attr_BBS$Breeding_AM_score2)

dispersal_model_bbs4 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Breeding_AM_score2 + (1|spp) + (1|pair.id), data=pair_attr_BBS)
summary(dispersal_model_bbs4)
anova(dispersal_model_bbs4)
## still not significant
results_table_mob_bbs <- data.frame(summary(dispersal_model_bbs4)$coefficients[,1:5]) ## 17 species
write.csv(results_table_mob_bbs, file = "../Results/Model_outputs/BBS/average_mob_bbs.csv", row.names=TRUE)


############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

############### CHANGE IN SYNCHORNY #################

##################################
#### Mobility model for UKBMS ####
##################################

## subset to only look at mid years 1985, 2000 and 2012
pair_attr_1985 <- pair_attr[pair_attr$mid.year==1984.5,]
pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]
pair_attr_early <- rbind(pair_attr_1985, pair_attr_2000) ## 2 years and 32 species
pair_attr_late <- rbind(pair_attr_2000, pair_attr_2012) ## 2 years and 32 species
## remove NA's to make sure small white is taken out
pair_attr_early <- na.omit(pair_attr_early)
pair_attr_late <- na.omit(pair_attr_late)

length(unique(pair_attr_early$spp)) # 31 species
length(unique(pair_attr_late$spp)) # 31 species

# ##### test with wingspan data ######
# pair_attr_ukbms1 <- merge(pair_attr_ukbms, wingspan, by.x="spp", by.y="BMScode")
# 
# pair_attr_ukbms1$mid.year <- as.factor(pair_attr_ukbms1$mid.year)
# pair_attr_ukbms1$pair.id <- as.character(pair_attr_ukbms1$pair.id)
# pair_attr_ukbms1$spp <- as.factor(pair_attr_ukbms1$spp)
# 
# wing_model4 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*Wingspan + (1|spp) + (1|pair.id), data=pair_attr_ukbms1)
# summary(wing_model4)
# anova(wing_model4)
# ## same result as with mobility data ==> butterflies with larger wingspan
# ## increase in syncrhony between 00-12, but not 85-00
#####
pair_attr_early$mid.year <- as.factor(pair_attr_early$mid.year)
pair_attr_early$pair.id <- as.character(pair_attr_early$pair.id)
pair_attr_early$spp <- as.factor(pair_attr_early$spp)
pair_attr_early$mobility_wil <- as.numeric(pair_attr_early$mobility_wil)

pair_attr_late$mid.year <- as.factor(pair_attr_late$mid.year)
pair_attr_late$pair.id <- as.character(pair_attr_late$pair.id)
pair_attr_late$spp <- as.factor(pair_attr_late$spp)
pair_attr_late$mobility_wil <- as.numeric(pair_attr_late$mobility_wil)

## rescale variables 
pair_attr_early$distance <- (pair_attr_early$distance - mean(na.omit(pair_attr_early$distance)))/sd(na.omit(pair_attr_early$distance))
pair_attr_early$mean_northing <- (pair_attr_early$mean_northing - mean(na.omit(pair_attr_early$mean_northing)))/sd(na.omit(pair_attr_early$mean_northing))
pair_attr_early$renk_hab_sim <- (pair_attr_early$renk_hab_sim - mean(na.omit(pair_attr_early$renk_hab_sim)))/sd(na.omit(pair_attr_early$renk_hab_sim))

pair_attr_late$distance <- (pair_attr_late$distance - mean(na.omit(pair_attr_late$distance)))/sd(na.omit(pair_attr_late$distance))
pair_attr_late$mean_northing <- (pair_attr_late$mean_northing - mean(na.omit(pair_attr_late$mean_northing)))/sd(na.omit(pair_attr_late$mean_northing))
pair_attr_late$renk_hab_sim <- (pair_attr_late$renk_hab_sim - mean(na.omit(pair_attr_late$renk_hab_sim)))/sd(na.omit(pair_attr_late$renk_hab_sim))

## full model with interaction between mid year and mobility
mobility_model5 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mobility_wil + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
summary(mobility_model5)
anova(mobility_model5)
## overall interaction very significant
## so the difference in mobility betwen the years differs depending on mobility ability

## run model with three mobility groups (low, med, high)
pair_attr_ukbms$mobility_wil <- as.numeric(pair_attr_ukbms$mobility_wil)
pair_attr_ukbms$mobility_score <- cut(pair_attr_ukbms$mobility_wil, 3, labels=FALSE)
pair_attr_ukbms$mobility_score <- as.factor(pair_attr_ukbms$mobility_score)

mobility_model6 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mobility_score + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
summary(mobility_model6)
anova(mobility_model6)
## interaction is very significant 

## run model with 2 mobility groups (low and high)
pair_attr_early$mobility_wil <- as.numeric(pair_attr_early$mobility_wil)
pair_attr_early$mobility_score2 <- cut(pair_attr_early$mobility_wil, 2, labels=FALSE)
pair_attr_early$mobility_score2 <- as.factor(pair_attr_early$mobility_score2)

pair_attr_late$mobility_wil <- as.numeric(pair_attr_late$mobility_wil)
pair_attr_late$mobility_score2 <- cut(pair_attr_late$mobility_wil, 2, labels=FALSE)
pair_attr_late$mobility_score2 <- as.factor(pair_attr_late$mobility_score2)

## EARLY MODEL
mobility_model7 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr_early)
summary(mobility_model7) ## both years are significant
anova(mobility_model7)
## interaction is  significant 
## save model output
results_table_mobility_ukbms <- data.frame(summary(mobility_model7)$coefficients[,1:5]) ## 31 species
write.csv(results_table_mobility_ukbms, file = "../Results/Model_outputs/UKBMS/change_mobility_ukbms_85_00.csv", row.names=TRUE)

## LATE MODEL
mobility_model8 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr_late)
summary(mobility_model8) ## both years are significant
anova(mobility_model8)
## interaction is  significant 
## save model output
results_table_mobility_ukbms <- data.frame(summary(mobility_model8)$coefficients[,1:5]) ## 31 species
write.csv(results_table_mobility_ukbms, file = "../Results/Model_outputs/UKBMS/change_mobility_ukbms_00_12.csv", row.names=TRUE)

### predict new data for EARLY model
pair_attr_early <- droplevels(pair_attr_early)
newdata_ukbms <- expand.grid(mean_northing=mean(pair_attr_early$mean_northing), distance=mean(pair_attr_early$distance), 
                       renk_hab_sim=mean(pair_attr_early$renk_hab_sim), pair.id=sample(pair_attr_early$pair.id,1),
                       spp=unique(pair_attr_early$spp), mid.year=unique(pair_attr_early$mid.year), 
                       mobility_score2=unique(pair_attr_early$mobility_score2))

newdata_ukbms$lag0 <- predict(mobility_model7, newdata=newdata_ukbms, re.form=NA)

mm <- model.matrix(terms(mobility_model7), newdata_ukbms)
pvar <- diag(mm %*% tcrossprod(vcov(mobility_model7),mm))
tvar <- pvar+VarCorr(mobility_model7)$spp[1]+VarCorr(mobility_model7)$pair.id[1]
cmult <- 2

newdata_ukbms <- data.frame(
  newdata_ukbms
  , plo = newdata_ukbms$lag0-1.96*sqrt(pvar)
  , phi = newdata_ukbms$lag0+1.96*sqrt(pvar)
  , tlo = newdata_ukbms$lag0-1.96*sqrt(tvar)
  , thi = newdata_ukbms$lag0+1.96*sqrt(tvar)
)

## run model without mobility or species random effect to obtain residuals
mob_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data=pair_attr_early)
pair_attr_early$residuals <- resid(mob_model, type="pearson")
group_by(pair_attr_early, mobility_score2) %>% summarize(m = mean(lag0)) ## decrease=0.25, increase=0.286
## put mean of each group into pair_attr dataframe
pair_attr_early <- ddply(pair_attr_early, "mobility_score2", transform, mob_mean = mean(lag0))
## add mean to each residual
pair_attr_early$residuals2 <- pair_attr_early$residuals + pair_attr_early$mob_mean

## create dataframe which calculates mean, SD and SE of residuals for each species
summary_ukbms <- pair_attr_early %>% group_by(spp, mobility_score2, mid.year) %>% 
  summarise_at(vars(residuals2), funs(mean,std.error))

## change values
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1984.5"="1985"))
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1999.5"="2000"))
colnames(newdata_ukbms)[7] <- "Mobility"
newdata_ukbms$Mobility <- revalue(newdata_ukbms$Mobility, c("1"="Low"))
newdata_ukbms$Mobility <- revalue(newdata_ukbms$Mobility, c("2"="High"))
## same for summary file
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("1984.5"="1985"))
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("1999.5"="2000"))
colnames(summary_ukbms)[2] <- "Mobility"
summary_ukbms$Mobility <- revalue(summary_ukbms$Mobility, c("1"="Low"))
summary_ukbms$Mobility <- revalue(summary_ukbms$Mobility, c("2"="High"))
## reorder levels
levels(summary_ukbms$Mobility)
summary_ukbms$Mobility <- factor(summary_ukbms$Mobility, levels=c("High", "Low"))
levels(summary_ukbms$Mobility)
levels(newdata_ukbms$Mobility)
newdata_ukbms$Mobility <- factor(newdata_ukbms$Mobility, levels=c("High", "Low"))
levels(newdata_ukbms$Mobility)

## plot graph with raw data residuals (+SE error bars) and fitted lines
png("../Graphs/Mobility/Mobility_change_predicted_ukbms_85_00.png", height = 150, width = 180, units = "mm", res = 300)
plot1<-ggplot(summary_ukbms, aes(x = mid.year, y = mean, group=Mobility)) +
  geom_point(aes(shape=Mobility), colour="grey", size = 3, position=myjit) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey", width=0.1, position=myjit) +
  geom_line(data=newdata_ukbms, aes(x=mid.year, y=lag0, linetype=Mobility), lwd=1) +
  labs(x="Mid year of moving window", y="Population synchrony") +
  theme_bw() +
  theme(legend.key.width = unit(0.8,"cm"), legend.key = element_rect(size = 2), text = element_text(size = 12), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.margin=margin(c(-10,20,-90,0)),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  scale_linetype_manual(name="Mobility",
                        labels=c("High", "Low"), values=c(1,2)) +
  scale_shape_manual(name=" ", 
                     labels=c("High", "Low"), values=c(16,4)) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(linetype = guide_legend(override.aes = list(size = 0.5)))
plot1
dev.off()

### predict new data for LATE model
pair_attr_late <- droplevels(pair_attr_late)
newdata_ukbms <- expand.grid(mean_northing=mean(pair_attr_late$mean_northing), distance=mean(pair_attr_late$distance), 
                             renk_hab_sim=mean(pair_attr_late$renk_hab_sim), pair.id=sample(pair_attr_late$pair.id,1),
                             spp=unique(pair_attr_late$spp), mid.year=unique(pair_attr_late$mid.year), 
                             mobility_score2=unique(pair_attr_late$mobility_score2))

newdata_ukbms$lag0 <- predict(mobility_model8, newdata=newdata_ukbms, re.form=NA)

mm <- model.matrix(terms(mobility_model8), newdata_ukbms)
pvar <- diag(mm %*% tcrossprod(vcov(mobility_model8),mm))
tvar <- pvar+VarCorr(mobility_model8)$spp[1]+VarCorr(mobility_model8)$pair.id[1]
cmult <- 2

newdata_ukbms <- data.frame(
  newdata_ukbms
  , plo = newdata_ukbms$lag0-1.96*sqrt(pvar)
  , phi = newdata_ukbms$lag0+1.96*sqrt(pvar)
  , tlo = newdata_ukbms$lag0-1.96*sqrt(tvar)
  , thi = newdata_ukbms$lag0+1.96*sqrt(tvar)
)

## run model without mobility or species random effect to obtain residuals
mob_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data=pair_attr_late)
pair_attr_late$residuals <- resid(mob_model, type="pearson")
group_by(pair_attr_late, mobility_score2) %>% summarize(m = mean(lag0)) ## decrease=0.25, increase=0.286
## put mean of each group into pair_attr dataframe
pair_attr_late <- ddply(pair_attr_late, "mobility_score2", transform, mob_mean = mean(lag0))
## add mean to each residual
pair_attr_late$residuals2 <- pair_attr_late$residuals + pair_attr_late$mob_mean

## create dataframe which calculates mean, SD and SE of residuals for each species
summary_ukbms <- pair_attr_late %>% group_by(spp, mobility_score2, mid.year) %>% 
  summarise_at(vars(residuals2), funs(mean,std.error))

## change values
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("2011.5"="2012"))
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1999.5"="2000"))
colnames(newdata_ukbms)[7] <- "Mobility"
newdata_ukbms$Mobility <- revalue(newdata_ukbms$Mobility, c("1"="Low"))
newdata_ukbms$Mobility <- revalue(newdata_ukbms$Mobility, c("2"="High"))
## same for summary file
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("2011.5"="2012"))
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("1999.5"="2000"))
colnames(summary_ukbms)[2] <- "Mobility"
summary_ukbms$Mobility <- revalue(summary_ukbms$Mobility, c("1"="Low"))
summary_ukbms$Mobility <- revalue(summary_ukbms$Mobility, c("2"="High"))
## reorder levels
levels(summary_ukbms$Mobility)
summary_ukbms$Mobility <- factor(summary_ukbms$Mobility, levels=c("High", "Low"))
levels(summary_ukbms$Mobility)
levels(newdata_ukbms$Mobility)
newdata_ukbms$Mobility <- factor(newdata_ukbms$Mobility, levels=c("High", "Low"))
levels(newdata_ukbms$Mobility)

## plot graph with raw data residuals (+SE error bars) and fitted lines
png("../Graphs/Mobility/Mobility_change_predicted_ukbms_00_12.png", height = 150, width = 180, units = "mm", res = 300)
plot2<-ggplot(summary_ukbms, aes(x = mid.year, y = mean, group=Mobility)) +
  geom_point(aes(shape=Mobility), colour="grey", size = 3, position=myjit) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey", width=0.1, position=myjit) +
  geom_line(data=newdata_ukbms, aes(x=mid.year, y=lag0, linetype=Mobility), lwd=1) +
  labs(x="Mid year of moving window", y="Population synchrony") +
  theme_bw() +
  theme(legend.key.width = unit(0.8,"cm"), legend.key = element_rect(size = 2), text = element_text(size = 12), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.margin=margin(c(-10,20,-90,0)),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  scale_linetype_manual(name="Mobility",
                        labels=c("High", "Low"), values=c(1,2)) +
  scale_shape_manual(name=" ", 
                     labels=c("High", "Low"), values=c(16,4)) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(linetype = guide_legend(override.aes = list(size = 0.5)))
plot2
dev.off()


############ jitter code #################
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

################################
#### Mobility model for CBC ####
################################

pair_attr_cbc_1985 <- pair_attr_CBC[pair_attr_CBC$mid.year==1984.5,]
pair_attr_cbc_1996 <- pair_attr_CBC[pair_attr_CBC$mid.year==1995.5,]
pair_attr_cbc <- rbind(pair_attr_cbc_1985, pair_attr_cbc_1996)

pair_attr_cbc$mid.year <- as.factor(pair_attr_cbc$mid.year)
pair_attr_cbc$pair.id <- as.character(pair_attr_cbc$pair.id)
pair_attr_cbc$spp <- as.factor(pair_attr_cbc$spp)

###### run model with strategy and year interaction
dispersal_model_cbc1 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*Breeding_AM + (1|pair.id) + (1|spp), data = pair_attr_cbc)
summary(dispersal_model_cbc1)
anova(dispersal_model_cbc1)
## mid.year*dispersal interaction is NON-significant overall (p=0.15)

## run model with three mobility groups (low, med, high)
pair_attr_cbc$Breeding_AM <- as.numeric(pair_attr_cbc$Breeding_AM)
pair_attr_cbc$Breeding_AM_score <- cut(pair_attr_cbc$Breeding_AM, 3, labels=FALSE)
pair_attr_cbc$Breeding_AM_score <- as.factor(pair_attr_cbc$Breeding_AM_score)

dispersal_model_cbc2 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*Breeding_AM_score + (1|spp) + (1|pair.id), data=pair_attr_cbc)
summary(dispersal_model_cbc2)
anova(dispersal_model_cbc2)

## run model with 2 mobility groups (low and high)
pair_attr_cbc$Breeding_AM <- as.numeric(pair_attr_cbc$Breeding_AM)
pair_attr_cbc$Breeding_AM_score2 <- cut(pair_attr_cbc$Breeding_AM, 2, labels=FALSE)
pair_attr_cbc$Breeding_AM_score2 <- as.factor(pair_attr_cbc$Breeding_AM_score2)
pair_attr_cbc <- na.omit(pair_attr_cbc)

dispersal_model_cbc3 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*Breeding_AM_score2 + (1|family) + (1|spp) + (1|pair.id), data=pair_attr_cbc)
## singular fit error (species RE has variance very close to 0)
summary(dispersal_model_cbc3)
## run model without species random effect and check if results are similar
dispersal_model_cbc4 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*Breeding_AM_score2 + (1|family) + (1|pair.id), data=pair_attr_cbc)
summary(dispersal_model_cbc4)
## results very similar
## save results from model with species and family RE

anova(dispersal_model_cbc3) ## interaction is non-significant (p=0.116)
## save model output
results_table_dispersal_cbc <- data.frame(summary(dispersal_model_cbc3)$coefficients[,1:5]) ## 20 species
write.csv(results_table_dispersal_cbc, file = "../Results/Model_outputs/CBC/change_mob_cbc.csv", row.names=TRUE)

###################################################################
####################### BBS birds #################################
###################################################################

pair_attr_bbs_1999 <- pair_attr_BBS[pair_attr_BBS$mid.year==1998.5,]
pair_attr_bbs_2012 <- pair_attr_BBS[pair_attr_BBS$mid.year==2011.5,]
pair_attr_bbs <- rbind(pair_attr_bbs_1999, pair_attr_bbs_2012)
pair_attr_bbs <- na.omit(pair_attr_bbs) # remove NA's (species without data)

pair_attr_bbs$mid.year <- as.factor(pair_attr_bbs$mid.year)
pair_attr_bbs$pair.id <- as.character(pair_attr_bbs$pair.id)
pair_attr_bbs$spp <- as.factor(pair_attr_bbs$spp)

###### run model with mobility and year interaction
dispersal_model_bbs1 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*Natal_AM + (1|pair.id) + (1|spp), data = pair_attr_bbs)
summary(dispersal_model_bbs1)
## intercept is mid.year 1999
anova(dispersal_model_bbs1)
## mid.year*dispersal interaction is significant overall

## run model with three mobility groups (low, med, high)
pair_attr_bbs$Breeding_AM <- as.numeric(pair_attr_bbs$Breeding_AM)
pair_attr_bbs$Breeding_AM_score <- cut(pair_attr_bbs$Breeding_AM, 3, labels=FALSE)
pair_attr_bbs$Breeding_AM_score <- as.factor(pair_attr_bbs$Breeding_AM_score)

## remove NAs (one species == redstart)
pair_attr_bbs <- na.omit(pair_attr_bbs)
summary(pair_attr_bbs)  

dispersal_model_bbs2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*Breeding_AM_score + (1|spp) + (1|pair.id), data=pair_attr_bbs)
summary(dispersal_model_bbs2)
anova(dispersal_model_bbs2)
## interaction is very significant (p<0.000001)

## run model with 2 mobility groups (low and high)
pair_attr_bbs$Breeding_AM <- as.numeric(pair_attr_bbs$Breeding_AM)
pair_attr_bbs$Breeding_AM_score2 <- cut(pair_attr_bbs$Breeding_AM, 2, labels=FALSE)
pair_attr_bbs$Breeding_AM_score2 <- as.factor(pair_attr_bbs$Breeding_AM_score2)

start_time <- Sys.time()
dispersal_model_bbs3 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*Breeding_AM_score2 + (1|spp) + (1|pair.id), data=pair_attr_bbs)
end_time <- Sys.time()
end_time - start_time ## 34 seconds

summary(dispersal_model_bbs3)
anova(dispersal_model_bbs3)
## interaction is very significant (p<0.000001)
## save model output
results_table_dispersal_bbs <- data.frame(summary(dispersal_model_bbs3)$coefficients[,1:5])
write.csv(results_table_dispersal_bbs, file = "../Results/Model_outputs/BBS/change_mob_bbs.csv", row.names=TRUE)

### predict new data
newdata_bbs <- expand.grid(mean_northing=mean(pair_attr_bbs$mean_northing), distance=mean(pair_attr_bbs$distance), 
                      renk_hab_sim=mean(pair_attr_bbs$renk_hab_sim), mid.year=unique(pair_attr_bbs$mid.year), 
                      pair.id=sample(pair_attr_bbs$pair.id,100), spp=unique(pair_attr_bbs$spp),
                      Breeding_AM_score2=unique(pair_attr_bbs$Breeding_AM_score2))
newdata_bbs$lag0 <- predict(dispersal_model_bbs3, newdata=newdata_bbs, re.form=NA)

mm3 <- model.matrix(terms(dispersal_model_bbs3), newdata_bbs)
pvar3 <- diag(mm3 %*% tcrossprod(vcov(dispersal_model_bbs3),mm3))
tvar3 <- pvar3+VarCorr(dispersal_model_bbs3)$spp[1]+VarCorr(dispersal_model_bbs3)$pair.id[1]
cmult <- 2

newdata_bbs <- data.frame(
  newdata_bbs
  , plo = newdata_bbs$lag0-1.96*sqrt(pvar3)
  , phi = newdata_bbs$lag0+1.96*sqrt(pvar3)
  , tlo = newdata_bbs$lag0-1.96*sqrt(tvar3)
  , thi = newdata_bbs$lag0+1.96*sqrt(tvar3)
)

## run model again without dispersal or species random effect to obtain residuals
mob_model_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data=pair_attr_bbs)
pair_attr_bbs$residuals <- resid(mob_model_bbs)

group_by(pair_attr_bbs, Breeding_AM_score2) %>% summarize(m = mean(lag0)) ## low mean = 0.0256, high mean = 0.0307
## put mean of each group into pair_attr dataframe
pair_attr_bbs <- ddply(pair_attr_bbs, "Breeding_AM_score2", transform, dispersal_mean = mean(lag0))
## add mean to each residual
pair_attr_bbs$residuals2 <- pair_attr_bbs$residuals + pair_attr_bbs$dispersal_mean

## create new dataframe to calculate mean, SD and SE of residuals for each species
summary_bbs <- pair_attr_bbs %>% group_by(spp, Breeding_AM_score2, mid.year) %>% 
  summarise_at(vars(residuals2), funs(mean,std.error))

## change year values
newdata_bbs$mid.year <- revalue(newdata_bbs$mid.year, c("1998.5"="1999"))
newdata_bbs$mid.year <- revalue(newdata_bbs$mid.year, c("2011.5"="2012"))
colnames(newdata_bbs)[7] <- "Dispersal"
newdata_bbs$Dispersal <- revalue(newdata_bbs$Dispersal, c("1"="Low"))
newdata_bbs$Dispersal <- revalue(newdata_bbs$Dispersal, c("2"="High"))
## same for summary file
summary_bbs$mid.year <- revalue(summary_bbs$mid.year, c("1998.5"="1999"))
summary_bbs$mid.year <- revalue(summary_bbs$mid.year, c("2011.5"="2012"))
colnames(summary_bbs)[2] <- "Dispersal"
summary_bbs$Dispersal <- revalue(summary_bbs$Dispersal, c("1"="Low"))
summary_bbs$Dispersal <- revalue(summary_bbs$Dispersal, c("2"="High"))
## reorder levels
levels(summary_bbs$Dispersal)
summary_bbs$Dispersal <- factor(summary_bbs$Dispersal, levels=c("High", "Low"))
levels(summary_bbs$Dispersal)
levels(newdata_bbs$Dispersal)
newdata_bbs$Dispersal <- factor(newdata_bbs$Dispersal, levels=c("High", "Low"))
levels(newdata_bbs$Dispersal)

## plot graph with raw data residuals (+SE errorbars) and fitted line
png("../Graphs/Mobility/Mobility_change_predicted_bbs.png", height = 100, width = 110, units = "mm", res = 300)
bbs_mob <- ggplot(summary_bbs, aes(x = mid.year, y = mean, group=Dispersal)) +
  geom_point(aes(shape=Dispersal), colour="grey66", size = 3, position=myjit) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey66", width=0.2, position=myjit) +
  geom_line(data=newdata_bbs, aes(x=mid.year, y=lag0, linetype=Dispersal), lwd=1) +
  labs(x="Mid year of moving window", y="Population synchrony") +
  scale_y_continuous(breaks=seq(-0.2,0.3,0.05)) +
  theme_bw() +
  theme(legend.key.width = unit(0.8,"cm"), legend.key = element_rect(size = 2), text = element_text(size = 12), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.margin=margin(c(-10,20,-80,0)),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  scale_linetype_manual(name=" ",
                        labels=c("High", "Low"), values=c(1,2)) +
  scale_shape_manual(name="Mobility", 
                     labels=c("High", "Low"), values=c(16,4)) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(linetype = guide_legend(override.aes = list(size = 0.5)))
bbs_mob
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

