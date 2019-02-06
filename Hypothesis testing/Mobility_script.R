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
options(scipen=999)

## read in synchrony data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE) 
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC.csv", header=TRUE) 
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE) 
wingspan <- read.csv("../Data/UKBMS_data/butterflywingspans.csv", header=TRUE)
bird_dispersal <- read.csv("../Data/Woodland_bird_dispersal_Paradis1998.csv", header=TRUE)

################################### AVERAGE SYNCRHONY ################################### 

## merge pair attr and wingspan data
pair_attr_1 <- merge(pair_attr, wingspan, by.x="spp", by.y="BMScode")

pair_attr_1$mid.year <- as.factor(pair_attr_1$mid.year)
pair_attr_1$pair.id <- as.character(pair_attr_1$pair.id)
pair_attr_1$spp <- as.factor(pair_attr_1$spp)

## change wingspan into 3 groups (small, med, large)
pair_attr_1$Wingspan <- as.numeric(pair_attr_1$Wingspan)
pair_attr_1$Wingspan_group <- cut(pair_attr_1$Wingspan, 3, labels=FALSE)
pair_attr_1$Wingspan_group <- as.factor(pair_attr_1$Wingspan_group)

wing_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Wingspan + (1|spp) + (1|pair.id), data=pair_attr_1)
summary(wing_model)
anova(wing_model)
## wingspan is non-significant (p=0.33)

wing_model2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Wingspan + Family + (1|spp) + (1|pair.id), data=pair_attr_1)
summary(wing_model2)
anova(wing_model2)
## neither family nor wingspan are significant 

wing_model3 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Wingspan + Sub.family + (1|spp) + (1|pair.id), data=pair_attr_1)
summary(wing_model3)
anova(wing_model3)
## neither sub-family nor wingspan are significant

wing_model4 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (Wingspan|Sub.family) + (1|spp) + (1|pair.id), data=pair_attr_1)
anova(wing_model4)
## doesn't work ==> model failed to converge

wing_model5 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Wingspan + (1|Sub.family) + (1|spp) + (1|pair.id), data=pair_attr_1)
anova(wing_model5)
## wingspan still not significant 

wing_model6 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Wingspan*Sub.family + (1|spp) + (1|pair.id), data=pair_attr_1)
anova(wing_model6)
## interaction is significant
## but warning: fixed-effect model matrix is rank deficient so dropping 3 columns / coefficients
## this is because 3 sub-families only have one wingspan, so can't model this
## Coliadinae, Limenitidinae, and Lycaeninae
## remove these sub-families and run model again
pair_attr_1 <- pair_attr_1[!pair_attr_1$Sub.family=="Coliadinae",]
pair_attr_1 <- pair_attr_1[!pair_attr_1$Sub.family=="Limenitidinae",]
pair_attr_1 <- pair_attr_1[!pair_attr_1$Sub.family=="Lycaeninae",]
## run model again
wing_model7 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Wingspan*Sub.family + (1|spp) + (1|pair.id), data=pair_attr_1)
anova(wing_model7)
summary(wing_model7)

## plot result
newdata_wing <- expand.grid(mean_northing=mean(pair_attr_1$mean_northing), distance=mean(pair_attr_1$distance), 
                             renk_hab_sim=mean(pair_attr_1$renk_hab_sim), pair.id=sample(pair_attr_1$pair.id,100),
                             spp=unique(pair_attr_1$spp), mid.year=unique(pair_attr_1$mid.year), 
                             Wingspan=unique(pair_attr_1$Wingspan), Sub.family=unique(pair_attr_1$Sub.family))

newdata_wing$lag0 <- predict(wing_model7, newdata=newdata_wing, re.form=NA)

mm <- model.matrix(terms(wing_model7), newdata_wing)
pvar <- diag(mm %*% tcrossprod(vcov(wing_model7),mm))
tvar <- pvar+VarCorr(wing_model7)$spp[1]+VarCorr(wing_model7)$pair.id[1]
cmult <- 2

newdata_ukbms <- data.frame(
  newdata_ukbms
  , plo = newdata_ukbms$lag0-1.96*sqrt(pvar)
  , phi = newdata_ukbms$lag0+1.96*sqrt(pvar)
  , tlo = newdata_ukbms$lag0-1.96*sqrt(tvar)
  , thi = newdata_ukbms$lag0+1.96*sqrt(tvar)
)

## run model without mobility or species random effect to obtain residuals
wing_mod <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data=pair_attr_1)
pair_attr_1$residuals <- resid(wing_mod, type="pearson")
## check mean synchrony of each group
group_by(pair_attr_1, Wingspan) %>% summarize(m = mean(lag0)) ## WC mean = 0.265, HS mean = 0.195
## put mean of each group into pair_attr dataframe
pair_attr_1 <- ddply(pair_attr_1, "Wingspan", transform, wing_mean = mean(lag0))
## add mean to each residual
pair_attr_1$residuals2 <- pair_attr_1$residuals + pair_attr_1$wing_mean

## create dataframe which calculates mean, SD and SE of residuals for each species
summary_wing2 <- pair_attr_1 %>% group_by(spp, Wingspan, Sub.family) %>% 
  summarise_at(vars(residuals2), funs(mean,std.error))

summary_wing2 <- summary_wing2[!summary_wing2$Sub.family=="Coliadinae",]
summary_wing2 <- summary_wing2[!summary_wing2$Sub.family=="Limenitidinae",]
summary_wing2 <- summary_wing2[!summary_wing2$Sub.family=="Lycaeninae",]

ggplot(summary_wing2, aes(x = Wingspan, y = mean, colour=Sub.family)) +
  geom_point(aes(colour=Sub.family), size = 2) +
  #geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error, colour=Sub.family), width=0.2) +
  geom_line(data=newdata_wing, aes(x=Wingspan, y=lag0, colour=Sub.family), lwd=1) +
  theme_bw() +
  theme(text = element_text(size = 12), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

## AIC of all models - which wingspan/sub-family model performs best?
AIC(wing_model3, wing_model4, wing_model5, wing_model6)
## for now wing_model5 has the lowest AIC

##################################
#### Mobility model for UKBMS ####
##################################

## remove NA's ==> removes Small White butterfly which doesn't have mobility data
pair_attr <- na.omit(pair_attr)
length(unique(pair_attr$spp)) # 32 species

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

start_time <- Sys.time()
mobility_model4 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr)
end_time <- Sys.time()
end_time - start_time ## 7.64mins 

summary(mobility_model4)
anova(mobility_model4)
## significant (p=0.011)
## mobility score 2 has higher average synchrony than score 1

#### compare models with 3 groups and 2 groups
AIC(mobility_model2, mobility_model4)
## AIC values are:
## mobility_model2 - 1916320
## mobility_model4 - 1916319
############################# very similar so just use 2 groups ############################# 

## don't plot graph because mobility and change in synchrony is also significant


##################################
#### Mobility model for CBC ####
##################################

## read in bird dispersal data
length(unique(bird_dispersal$Species_code)) ## only 23 species have dispersal data

## merge pair_attr with bird dispersal data
pair_attr_CBC <- merge(pair_attr_CBC, bird_dispersal, by.x="spp", by.y="Species_code", all=FALSE)

## natal arithmetic mean dispersal
dispersal_model_cbc1 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + Natal_AM + (1|spp) + (1|pair.id), data=pair_attr_CBC)
summary(dispersal_model1)
anova(dispersal_model1)
## not significant

## breeding arithmetic mean dispersal
## remove NAs (one species == redstart)
pair_attr_CBC <- na.omit(pair_attr_CBC)
summary(pair_attr_CBC)  

dispersal_model_cbc2 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + Breeding_AM + (1|spp) + (1|pair.id), data=pair_attr_CBC)
summary(dispersal_model2)
anova(dispersal_model2)
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

dispersal_model_cbc4 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + Breeding_AM_score2 + (1|spp) + (1|pair.id), data=pair_attr_CBC)
summary(dispersal_model_cbc4)
anova(dispersal_model_cbc4)
## still not significant

## compare models
AIC(dispersal_model_cbc3, dispersal_model_cbc4)
## model 4 slightly lower AIC => therefore performs better as two groups (low and high dispersal)

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

## compare models
AIC(dispersal_model_bbs3, dispersal_model_bbs4)
## model 4 slightly lower AIC => therefore performs better as two groups (low and high dispersal)


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
pair_attr_ukbms <- rbind(pair_attr_1985, pair_attr_2000, pair_attr_2012) ## 3 years and 32 species
## remove NA's to make sure small white is taken out
pair_attr_ukbms <- na.omit(pair_attr_ukbms)
length(unique(pair_attr_ukbms$spp)) # 32 species

## merge wingspan data
pair_attr_ukbms1 <- merge(pair_attr_ukbms, wingspan, by.x="spp", by.y="BMScode")

pair_attr_ukbms1$mid.year <- as.factor(pair_attr_ukbms1$mid.year)
pair_attr_ukbms1$pair.id <- as.character(pair_attr_ukbms1$pair.id)
pair_attr_ukbms1$spp <- as.factor(pair_attr_ukbms1$spp)

wing_model4 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*Wingspan + (1|spp) + (1|pair.id), data=pair_attr_ukbms1)
summary(wing_model4)
anova(wing_model4)
## same result as with mobility data ==> butterflies with larger wingspan
## increase in syncrhony between 00-12, but not 85-00

pair_attr_ukbms$mid.year <- as.factor(pair_attr_ukbms$mid.year)
pair_attr_ukbms$pair.id <- as.character(pair_attr_ukbms$pair.id)
pair_attr_ukbms$spp <- as.factor(pair_attr_ukbms$spp)
pair_attr_ukbms$mobility_wil <- as.numeric(pair_attr_ukbms$mobility_wil)

## rescale variables 
pair_attr_ukbms$distance <- (pair_attr_ukbms$distance - mean(na.omit(pair_attr_ukbms$distance)))/sd(na.omit(pair_attr_ukbms$distance))
pair_attr_ukbms$mean_northing <- (pair_attr_ukbms$mean_northing - mean(na.omit(pair_attr_ukbms$mean_northing)))/sd(na.omit(pair_attr_ukbms$mean_northing))
pair_attr_ukbms$renk_hab_sim <- (pair_attr_ukbms$renk_hab_sim - mean(na.omit(pair_attr_ukbms$renk_hab_sim)))/sd(na.omit(pair_attr_ukbms$renk_hab_sim))

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
## interaction is very significant (p=<0.00001)

## run model with 2 mobility groups (low and high)
pair_attr_ukbms$mobility_wil <- as.numeric(pair_attr_ukbms$mobility_wil)
pair_attr_ukbms$mobility_score2 <- cut(pair_attr_ukbms$mobility_wil, 2, labels=FALSE)
pair_attr_ukbms$mobility_score2 <- as.factor(pair_attr_ukbms$mobility_score2)

start_time <- Sys.time()
mobility_model7 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
end_time <- Sys.time()
end_time - start_time ## 22.9 seconds

summary(mobility_model7)
anova(mobility_model7)
## interaction is  significant (p=0.00399)
## save model output
results_table_mobility_ukbms <- data.frame(summary(mobility_model7)$coefficients[,1:5])
write.csv(results_table_mobility_ukbms, file = "../Results/Model_outputs/change_mobility_ukbms_nomigrants.csv", row.names=TRUE)

#### compare models with 3 groups and 2 groups
AIC(mobility_model6, mobility_model7)
## AIC values are:
## mobility_model6 - 222980
## mobility_model7 - 223378.2
## group of 3 is better...
## but still continue with 2 groups of mobility as easier to interpret


### predict new data
newdata_ukbms <- expand.grid(mean_northing=mean(pair_attr_ukbms$mean_northing), distance=mean(pair_attr_ukbms$distance), 
                       renk_hab_sim=mean(pair_attr_ukbms$renk_hab_sim), pair.id=sample(pair_attr_ukbms$pair.id,10),
                       spp=sample(pair_attr_ukbms$spp,10), mid.year=unique(pair_attr_ukbms$mid.year), 
                       mobility_score2=unique(pair_attr_ukbms$mobility_score2))

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
mob_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data=pair_attr_ukbms)
pair_attr_ukbms$residuals <- resid(mob_model, type="pearson")

## create dataframe which calculates mean, SD and SE of residuals for each species
summary_ukbms <- pair_attr_ukbms %>% group_by(spp, mobility_score2, mid.year) %>% 
  summarise_at(vars(residuals), funs(mean,std.error))

## change values
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1984.5"="1985"))
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1999.5"="2000"))
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("2011.5"="2012"))
colnames(newdata_ukbms)[7] <- "Mobility"
newdata_ukbms$Mobility <- revalue(newdata_ukbms$Mobility, c("1"="Low"))
newdata_ukbms$Mobility <- revalue(newdata_ukbms$Mobility, c("2"="High"))
## same for summary file
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("1984.5"="1985"))
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("1999.5"="2000"))
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("2011.5"="2012"))
colnames(summary_ukbms)[2] <- "Mobility"
summary_ukbms$Mobility <- revalue(summary_ukbms$Mobility, c("1"="Low"))
summary_ukbms$Mobility <- revalue(summary_ukbms$Mobility, c("2"="High"))

## plot graph with raw data residuals (+SE error bars) and fitted lines
png("../Graphs/Mobility/Mobility_change_predicted_ukbms.png", height = 80, width = 120, units = "mm", res = 300)
pd <- position_dodge(0.1)
ggplot(summary_ukbms, aes(x = mid.year, y = mean, group=Mobility)) +
  geom_point(aes(shape=Mobility), colour="grey", size = 2, position=pd) +
  scale_shape_manual(values=c(16,4)) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey", width=0.1, position=pd) +
  geom_line(data=newdata_ukbms, aes(x=mid.year, y=lag0, linetype=Mobility), lwd=1) +
  labs(x="Mid year of moving window", y="Population synchrony") +
  theme_bw() +
  theme(text = element_text(size = 10), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

## plot graph with fitted line and confidence interval error bars
png("../Graphs/Mobility/Mobility_change_predicted_ukbms2.png", height = 100, width = 120, units = "mm", res = 300)
ggplot(newdata_ukbms, aes(x=mid.year, y=lag0, group=Mobility)) +
  geom_point(position=pd) +
  geom_line(aes(linetype=Mobility), size=1, position=pd) +
  labs(x="Year", y="Population synchrony") +
  geom_linerange(aes(ymin=plo, ymax=phi), position=pd) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

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
## interaction is significant (p=0.0354)

## run model with 2 mobility groups (low and high)
pair_attr_cbc$Breeding_AM <- as.numeric(pair_attr_cbc$Breeding_AM)
pair_attr_cbc$Breeding_AM_score2 <- cut(pair_attr_cbc$Breeding_AM, 2, labels=FALSE)
pair_attr_cbc$Breeding_AM_score2 <- as.factor(pair_attr_cbc$Breeding_AM_score2)
pair_attr_cbc <- na.omit(pair_attr_cbc)

start_time <- Sys.time()
dispersal_model_cbc3 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*Breeding_AM_score2 + (1|spp) + (1|pair.id), data=pair_attr_cbc)
end_time <- Sys.time()
end_time - start_time ## 4.87 seconds

summary(dispersal_model_cbc3)
anova(dispersal_model_cbc3)
## interaction is significant (p=0.0418)
## save model output
results_table_dispersal_cbc <- data.frame(summary(dispersal_model_cbc3)$coefficients[,1:5])
write.csv(results_table_dispersal_cbc, file = "../Results/Model_outputs/change_dispersal_cbc.csv", row.names=TRUE)

## save main (true) model results
main_result_table <- data.frame(anova(dispersal_model_cbc3)[,5:6]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## remove rows with mean northing, distance, hab sim and mid year F values (only interested in abundance F values)
main_result_table <- main_result_table[ !(main_result_table$parameter %in% c("mean_northing", "distance", "hab_sim", "mid.year", "Breeding_AM_score2")), ]

library(foreach)
library(doParallel)
library(snow)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

### parallel ###
## run model 999 times
perm_cbc_mob <- NULL
system.time(
foreach (i=1:trials, .packages='lme4', .combine = combine) %dopar% {
  print(i)
  pair_attr_cbc$mob_shuffle <- sample(pair_attr_cbc$Breeding_AM_score2) ## randomly shuffle abundance change varaible
  abund_model_cbc_shuffle <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*mob_shuffle + (1|pair.id) + (1|spp), data = pair_attr_cbc)
  ## run model with shuffled variable
  ## save results
  results_table_temp <- data.frame(anova(abund_model_cbc_shuffle)[,5:6],i)
  perm_cbc_mob<-rbind(perm_cbc_mob,results_table_temp)
  }
)
stopCluster(cl)


### non-parallel ###
## run model 999 times
perm_cbc_mob <- NULL
system.time(
  for (i in 1:9) {
    print(i)
  pair_attr_cbc$mob_shuffle <- sample(pair_attr_cbc$Breeding_AM_score2) ## randomly shuffle abundance change varaible
  abund_model_cbc_shuffle <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*mob_shuffle + (1|pair.id) + (1|spp), data = pair_attr_cbc)
  ## run model with shuffled variable
  ## save results
  results_table_temp <- data.frame(anova(abund_model_cbc_shuffle)[,5:6],i)
  perm_cbc_mob<-rbind(perm_cbc_mob,results_table_temp)
  }
) ## 33.16 seconds (9 runs)

perm_cbc_mob$parameter <- paste(row.names(perm_cbc_mob)) ## move row.names to parameter column
rownames(perm_cbc_mob) <- 1:nrow(perm_cbc_mob) ## change row names to numbers
## remove rows with mean northing, distance, hab sim and mid year F values (only interested in abundance F values)
perm_cbc_mob <- perm_cbc_mob[-grep("mean_northing", perm_cbc_mob$parameter),]
perm_cbc_mob <- perm_cbc_mob[-grep("distance", perm_cbc_mob$parameter),]
perm_cbc_mob <- perm_cbc_mob[-grep("hab_sim", perm_cbc_mob$parameter),]
perm_cbc_mob <- perm_cbc_mob[-grep("mid.year", perm_cbc_mob$parameter),]

final_results_table <- rbind(main_result_table, perm_cbc_mob) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/perm_change_mob_cbc.csv", row.names=TRUE)

## Calculate p value
number_of_permutations <- 1000
final_results_table2 <- final_results_table[!final_results_table$i==0,] ## remove true value to calc. p value
diff.observed <- main_result_table$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(final_results_table2$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0.023 significant


### predict new data
newdata_cbc <- expand.grid(mean_northing=mean(pair_attr_cbc$mean_northing), distance=mean(pair_attr_cbc$distance), 
                      hab_sim=sample(pair_attr_cbc$hab_sim,1), mid.year=unique(pair_attr_cbc$mid.year), 
                      pair.id=sample(pair_attr_cbc$pair.id,10), spp=sample(pair_attr_cbc$spp,10),
                      Breeding_AM_score2=unique(pair_attr_cbc$Breeding_AM_score2))
newdata_cbc$lag0 <- predict(dispersal_model_cbc3, newdata=newdata_cbc, re.form=NA)

mm2 <- model.matrix(terms(dispersal_model_cbc3), newdata_cbc)
pvar2 <- diag(mm2 %*% tcrossprod(vcov(dispersal_model_cbc3),mm2))
tvar2 <- pvar2+VarCorr(dispersal_model_cbc3)$spp[1]+VarCorr(dispersal_model_cbc3)$pair.id[1]
cmult <- 2

newdata_cbc <- data.frame(
  newdata_cbc
  , plo = newdata_cbc$lag0-1.96*sqrt(pvar2)
  , phi = newdata_cbc$lag0+1.96*sqrt(pvar2)
  , tlo = newdata_cbc$lag0-1.96*sqrt(tvar2)
  , thi = newdata_cbc$lag0+1.96*sqrt(tvar2)
)

## run model without dispersal and species random effect to obtain residuals
dispersal_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + (1|pair.id), data=pair_attr_cbc)
pair_attr_cbc$residuals <- resid(dispersal_cbc)

## create new dataframe to calculate mean, SD and SE of residuals for each species
summary_cbc <- pair_attr_cbc %>% group_by(spp, Breeding_AM_score2, mid.year) %>% 
  summarise_at(vars(residuals), funs(mean,std.error))

## change values
newdata_cbc$mid.year <- revalue(newdata_cbc$mid.year, c("1984.5"="1985"))
newdata_cbc$mid.year <- revalue(newdata_cbc$mid.year, c("1995.5"="1996"))
colnames(newdata_cbc)[7] <- "Dispersal"
newdata_cbc$Dispersal <- revalue(newdata_cbc$Dispersal, c("1"="Low"))
newdata_cbc$Dispersal <- revalue(newdata_cbc$Dispersal, c("2"="High"))
## same for summary dataframe
summary_cbc$mid.year <- revalue(summary_cbc$mid.year, c("1984.5"="1985"))
summary_cbc$mid.year <- revalue(summary_cbc$mid.year, c("1995.5"="1996"))
colnames(summary_cbc)[2] <- "Dispersal"
summary_cbc$Dispersal <- revalue(summary_cbc$Dispersal, c("1"="Low"))
summary_cbc$Dispersal <- revalue(summary_cbc$Dispersal, c("2"="High"))

png("../Graphs/Mobility/Mobility_change_predicted_cbc.png", height = 80, width = 120, units = "mm", res = 300)
pd <- position_dodge(0.1)
ggplot(summary_cbc, aes(x = mid.year, y = mean, group=Dispersal)) +
  geom_point(aes(shape=Dispersal), colour="grey", size = 2, position=pd) +
  scale_shape_manual(values=c(16,4)) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey", width=0.1, position=pd) +
  geom_line(data=newdata_cbc, aes(x=mid.year, y=lag0, linetype=Dispersal), lwd=1) +
  labs(x="Mid year of moving window", y="Population synchrony") +
  theme_bw() +
  theme(text = element_text(size = 10), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

png("../Graphs/Mobility/Mobility_change_predicted_cbc2.png", height = 100, width = 120, units = "mm", res = 300)
pd <- position_dodge(0.1)
ggplot(newdata_cbc, aes(x=mid.year, y=lag0, group=Dispersal)) +
  geom_point(position=pd) +
  geom_line(aes(linetype=Dispersal), size=1, position=pd) +
  labs(x="Year", y="Population synchrony") +
  geom_linerange(aes(ymin=plo, ymax=phi), position=pd) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

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
write.csv(results_table_dispersal_bbs, file = "../Results/Model_outputs/change_dispersal_bbs.csv", row.names=TRUE)

## save main (true) model results
main_result_table <- data.frame(anova(dispersal_model_bbs3)[,5:6]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## remove rows with mean northing, distance, hab sim and mid year F values (only interested in abundance F values)
main_result_table <- main_result_table[ !(main_result_table$parameter %in% c("mean_northing", "distance", "renk_hab_sim", "mid.year", "Breeding_AM_score2")), ]

## run model 999 times
perm_bbs_mob <- NULL
start_time <- Sys.time()
for (i in 1:999){
  print(i)
  pair_attr_bbs$mob_shuffle <- sample(pair_attr_bbs$Breeding_AM_score2) ## randomly shuffle abundance change varaible
  abund_model_bbs_shuffle <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mob_shuffle + (1|pair.id) + (1|spp), data = pair_attr_bbs)
  ## run model with shuffled variable
  ## save results
  results_table_temp <- data.frame(anova(abund_model_bbs_shuffle)[,5:6],i)
  perm_bbs_mob<-rbind(perm_bbs_mob,results_table_temp)
}
end_time <- Sys.time()
end_time - start_time ## 8.81 hours to do 999 runs

perm_bbs_mob$parameter <- paste(row.names(perm_bbs_mob)) ## move row.names to parameter column
rownames(perm_bbs_mob) <- 1:nrow(perm_bbs_mob) ## change row names to numbers
## remove rows with mean northing, distance, hab sim and mid year F values (only interested in abundance F values)
perm_bbs_mob <- perm_bbs_mob[-grep("mean_northing", perm_bbs_mob$parameter),]
perm_bbs_mob <- perm_bbs_mob[-grep("distance", perm_bbs_mob$parameter),]
perm_bbs_mob <- perm_bbs_mob[-grep("hab_sim", perm_bbs_mob$parameter),]
perm_bbs_mob <- perm_bbs_mob[-grep("mid.year", perm_bbs_mob$parameter),]

final_results_table <- rbind(main_result_table, perm_bbs_mob) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/perm_change_mob_bbs.csv", row.names=TRUE)

## Calculate p value
number_of_permutations <- 1000
final_results_table2 <- final_results_table[!final_results_table$i==0,] ## remove true value to calc. p value
diff.observed <- main_result_table$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(final_results_table2$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 (exactly) significant


#### compare models with 3 groups and 2 groups
AIC(dispersal_model_bbs2, dispersal_model_bbs3)
## AIC values are:
## mobility_model6 - 770224.6
## mobility_model7 - 770223.5
## group of 2 is better :) 

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
ggplot(summary_bbs, aes(x = mid.year, y = mean, group=Dispersal)) +
  geom_point(aes(shape=Dispersal), colour="grey66", size = 2, position=myjit) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey66", width=0.2, position=myjit) +
  geom_line(data=newdata_bbs, aes(x=mid.year, y=lag0, linetype=Dispersal), lwd=0.5) +
  labs(x="Mid year of moving window", y="Population synchrony") +
  scale_y_continuous(breaks=seq(-0.2,0.3,0.05)) +
  theme_bw() +
  theme(legend.key.width = unit(0.8,"cm"), legend.key = element_rect(size = 2), text = element_text(size = 8), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.margin=margin(c(-10,20,-80,0)),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  scale_linetype_manual(name=" ",
                        labels=c("High", "Low"), values=c(1,2)) +
  scale_shape_manual(name="Mobility", 
                     labels=c("High", "Low"), values=c(16,4)) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(linetype = guide_legend(override.aes = list(size = 0.5)))
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

## plot graph with fitted line and confidence interval error bars
png("../Graphs/Mobility/Mobility_change_predicted_bbs2.png", height = 100, width = 120, units = "mm", res = 300)
ggplot(newdata_bbs, aes(x=mid.year, y=lag0, group=Dispersal)) +
  geom_point(position=pd) +
  geom_line(aes(linetype=Dispersal), size=1, position=pd) +
  labs(x="Year", y="Population synchrony") +
  geom_linerange(aes(ymin=plo, ymax=phi), position=pd) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

