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
options(scipen=999)

## read in synchrony data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE) 
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC.csv", header=TRUE) 
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE) 
traits <- read.csv("../Data/UKBMS_data/species.traits.full.table.csv", header=TRUE)
######## all butterflies with mobitliy data (=32 species) #########

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
pair_attr$mobility_score2 <- cut(pair_attr$mobility_wil, 2, labels=FALSE)
pair_attr$mobility_score2 <- as.factor(pair_attr$mobility_score2)

mobility_model4 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr)
summary(mobility_model4)
anova(mobility_model4)
## significant (p=0.011)
## mobility score 2 has higher average synchrony than score 1

#### compare models with 3 groups and 2 groups
AIC(mobility_model2, mobility_model4)
## AIC values are:
## mobility_model2 - 1916320
## mobility_model4 - 1916319
## very similar so just use 2 groups

## need to plot graph... 2 groups of mobility (low and high)
## won't use this graph if mobility and change in synchrony is significant
results_final_ukbms <- read.csv("../Results/Butterfly_results/results_final_sp.csv", header=TRUE)
## merge mobility data 
results_final_ukbms <- merge(results_final_ukbms, mobility, by.x="sp", by.y="species", all=FALSE)
## create mobility score
str(results_final_ukbms)
results_final_ukbms$mobility_wil <- as.numeric(results_final_ukbms$mobility_wil)
results_final_ukbms$mobility_score <- cut(results_final_ukbms$mobility_wil, 2, labels=FALSE)
results_final_ukbms$mobility_score <- as.factor(results_final_ukbms$mobility_score)

results_final_ukbms$mobility_score <- revalue(results_final_ukbms$mobility_score, c("1"="Low"))
results_final_ukbms$mobility_score <- revalue(results_final_ukbms$mobility_score, c("2"="High"))

## plot graph
png("../Graphs/Mobility/Mobility_score_average_boxplots_ukbms.png", height = 120, width = 120, units = "mm", res = 300)
ggplot(results_final_ukbms, aes(x=mobility_score, y=FCI)) +
  geom_boxplot() +
  labs(x="Mobility", y="Functional connectivity index") +
  scale_y_continuous(breaks=seq(-0.2,1.2,0.1)) +
  #geom_signif(comparisons = list(c("Low", "High")), map_signif_level = TRUE, textsize = 4, annotation=( "**"), step_increase = 0.1) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

########### CBC birds #################

## read in bird dispersal data
bird_dispersal <- read.csv("../Data/Woodland_bird_dispersal_Paradis1998.csv", header=TRUE)
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
pair_attr_CBC$Breeding_AM_score2 <- cut(pair_attr_CBC$Breeding_AM, 2, labels=FALSE)
pair_attr_CBC$Breeding_AM_score2 <- as.factor(pair_attr_CBC$Breeding_AM_score2)

dispersal_model_cbc4 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + Breeding_AM_score2 + (1|spp) + (1|pair.id), data=pair_attr_CBC)
summary(dispersal_model_cbc4)
anova(dispersal_model_cbc4)
## still not significant

## compare models
AIC(dispersal_model_cbc3, dispersal_model_cbc4)
## model 4 slightly lower AIC => therefore performs better as two groups (low and high dispersal)

########### BBS birds #################

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

## UKBMS butterflies
## subset to only look at mid years 1985, 2000 and 2012
pair_attr_1985 <- pair_attr[pair_attr$mid.year==1984.5,]
pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]
pair_attr_ukbms <- rbind(pair_attr_1985, pair_attr_2000, pair_attr_2012) ## 3 years and 32 species

pair_attr_ukbms$mid.year <- as.factor(pair_attr_ukbms$mid.year)
pair_attr_ukbms$pair.id <- as.character(pair_attr_ukbms$pair.id)
pair_attr_ukbms$spp <- as.factor(pair_attr_ukbms$spp)
pair_attr_ukbms$mobility_wil <- as.numeric(pair_attr_ukbms$mobility_wil)

## rescale variables 
pair_attr_ukbms$distance <- (pair_attr_ukbms$distance - mean(na.omit(pair_attr_ukbms$distance)))/sd(na.omit(pair_attr_ukbms$distance))
pair_attr_ukbms$mean_northing <- (pair_attr_ukbms$mean_northing - mean(na.omit(pair_attr_ukbms$mean_northing)))/sd(na.omit(pair_attr_ukbms$mean_northing))
pair_attr_ukbms$renk_hab_sim <- (pair_attr_ukbms$renk_hab_sim - mean(na.omit(pair_attr_ukbms$renk_hab_sim)))/sd(na.omit(pair_attr_ukbms$renk_hab_sim))

mob_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
pair_attr_ukbms$residuals <- resid(mob_model, type="pearson")
pair_attr_ukbms$predicted <- predict(mob_model)

library(plotrix)
summary4 <- pair_attr_ukbms %>% group_by(spp, mobility_score2, mid.year) %>% 
  summarise_at(vars(predicted), funs(mean,std.error))

ggplot(summary2, aes(x = mid.year, y = mean), group_by(mobility_score2)) +
  geom_point(aes(color=mobility_score2), size = 1) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error, colour=mobility_score2, width=0.06))
  #geom_line(aes(colour=mobility_score2), size=1) 
plot(mobility_model7)
qqnorm(sresid, cex=1.8, pch=10)
qqline(sresid, lty=2, lwd=2)

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

mobility_model7 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
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
## but used group of 2 above so continue with that?

#### plot graph
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

colnames(newdata_ukbms)[7] <- "Mobility"

## change year values
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1984.5"="1985"))
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1999.5"="2000"))
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("2011.5"="2012"))
newdata_ukbms$Mobility <- revalue(newdata_ukbms$Mobility, c("1"="Low"))
newdata_ukbms$Mobility <- revalue(newdata_ukbms$Mobility, c("2"="High"))
pd <- position_dodge(0.1)

png("../Graphs/Mobility/Mobility_change_predicted_ukbms.png", height = 100, width = 120, units = "mm", res = 300)
ggplot(newdata_ukbms, aes(x=mid.year, y=lag0, group=Mobility)) +
  geom_point(position=pd) +
  geom_line(aes(linetype=Mobility), size=1, position=pd) +
  labs(x="Year", y="Population synchrony") +
  geom_linerange(aes(ymin=plo, ymax=phi), position=pd) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

## CBC birds
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

mob_model <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + (1|spp) + (1|pair.id), data=pair_attr_cbc)
pair_attr_cbc$residuals <- resid(mob_model)
pair_attr_cbc$predicted <- predict(mob_model)

library(plotrix)
summary2 <- pair_attr_cbc %>% group_by(spp, Breeding_AM_score2, mid.year) %>% 
  summarise_at(vars(predicted), funs(mean,std.error))

ggplot(summary2, aes(x = mid.year, y = mean), group_by(Breeding_AM_score2)) +
  geom_point(aes(color=Breeding_AM_score2), size = 1) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error, colour=Breeding_AM_score2, width=0.06))


#geom_line(aes(colour=mobility_score2), size=1) 
plot(mobility_model7)
qqnorm(sresid, cex=1.8, pch=10)
qqline(sresid, lty=2, lwd=2)


dispersal_model_cbc3 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*Breeding_AM_score2 + (1|spp) + (1|pair.id), data=pair_attr_cbc)
summary(dispersal_model_cbc3)
anova(dispersal_model_cbc3)
## interaction is significant (p=0.0418)
## save model output
results_table_dispersal_cbc <- data.frame(summary(dispersal_model_cbc3)$coefficients[,1:5])
write.csv(results_table_dispersal_cbc, file = "../Results/Model_outputs/change_dispersal_cbc.csv", row.names=TRUE)

#### compare models with 3 groups and 2 groups
AIC(dispersal_model_cbc2, dispersal_model_cbc3)
## AIC values are:
## mobility_model6 - 12779.42
## mobility_model7 - 12770.79
## group of 2 is better

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

## change year values
newdata_cbc$mid.year <- revalue(newdata_cbc$mid.year, c("1984.5"="1985"))
newdata_cbc$mid.year <- revalue(newdata_cbc$mid.year, c("1995.5"="1996"))
colnames(newdata_cbc)[7] <- "Dispersal"
newdata_cbc$Dispersal <- revalue(newdata_cbc$Dispersal, c("1"="Low"))
newdata_cbc$Dispersal <- revalue(newdata_cbc$Dispersal, c("2"="High"))

pd <- position_dodge(0.1)
png("../Graphs/Mobility/Mobility_change_predicted_cbc.png", height = 100, width = 120, units = "mm", res = 300)
ggplot(newdata_cbc, aes(x=mid.year, y=lag0, group=Dispersal)) +
  geom_point(position=pd) +
  geom_line(aes(linetype=Dispersal), size=1, position=pd) +
  labs(x="Year", y="Population synchrony") +
  geom_linerange(aes(ymin=plo, ymax=phi), position=pd) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

## BBS birds
pair_attr_bbs_1999 <- pair_attr_BBS[pair_attr_BBS$mid.year==1998.5,]
pair_attr_bbs_2012 <- pair_attr_BBS[pair_attr_BBS$mid.year==2011.5,]
pair_attr_bbs <- rbind(pair_attr_bbs_1999, pair_attr_bbs_2012)
pair_attr_bbs <- na.omit(pair_attr_bbs) # remove NA's (species without data)

pair_attr_bbs$mid.year <- as.factor(pair_attr_bbs$mid.year)
pair_attr_bbs$pair.id <- as.character(pair_attr_bbs$pair.id)
pair_attr_bbs$spp <- as.factor(pair_attr_bbs$spp)

###### run model with strategy and year interaction
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

dispersal_model_bbs3 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*Breeding_AM_score2 + (1|spp) + (1|pair.id), data=pair_attr_bbs)
summary(dispersal_model_bbs3)
anova(dispersal_model_bbs3)
## interaction is very significant (p<0.000001)
## save model output
results_table_dispersal_bbs <- data.frame(summary(dispersal_model_bbs3)$coefficients[,1:5])
write.csv(results_table_dispersal_bbs, file = "../Results/Model_outputs/change_dispersal_bbs.csv", row.names=TRUE)

#### compare models with 3 groups and 2 groups
AIC(dispersal_model_bbs2, dispersal_model_bbs3)
## AIC values are:
## mobility_model6 - 770224.6
## mobility_model7 - 770223.5
## group of 2 is better :) 

### predict new data
newdata_bbs <- expand.grid(mean_northing=mean(pair_attr_bbs$mean_northing), distance=mean(pair_attr_bbs$distance), 
                      renk_hab_sim=mean(pair_attr_bbs$renk_hab_sim), mid.year=unique(pair_attr_bbs$mid.year), 
                      pair.id=sample(pair_attr_bbs$pair.id,10), spp=sample(pair_attr_bbs$spp,10),
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

## change year values
newdata_bbs$mid.year <- revalue(newdata_bbs$mid.year, c("1998.5"="1999"))
newdata_bbs$mid.year <- revalue(newdata_bbs$mid.year, c("2011.5"="2012"))
colnames(newdata_bbs)[7] <- "Dispersal"
newdata_bbs$Dispersal <- revalue(newdata_bbs$Dispersal, c("1"="Low"))
newdata_bbs$Dispersal <- revalue(newdata_bbs$Dispersal, c("2"="High"))

png("../Graphs/Mobility/Mobility_change_predicted_bbs.png", height = 100, width = 120, units = "mm", res = 300)
ggplot(newdata_bbs, aes(x=mid.year, y=lag0, group=Dispersal)) +
  geom_point(position=pd) +
  geom_line(aes(linetype=Dispersal), size=1, position=pd) +
  labs(x="Year", y="Population synchrony") +
  geom_linerange(aes(ymin=plo, ymax=phi), position=pd) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

