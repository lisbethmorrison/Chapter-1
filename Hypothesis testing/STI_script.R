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
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC_no_zeros2.csv", header=TRUE) 
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE) 
bird_STI <- read.csv("../Data/birds_STI.csv", header=TRUE)
butterfly_STI <- read.csv("../Data/butterflies_STI.csv", header=TRUE)
bird_dispersal <- read.csv("../Data/Woodland_bird_dispersal_Paradis1998.csv", header=TRUE)

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
## merge datasets
pair_attr_CBC <- merge(pair_attr_CBC, bird_STI, by.x="spp", by.y="species_code") ## 31 species
length(unique(pair_attr_CBC$spp)) # 29 species
str(pair_attr_CBC)

pair_attr_CBC$pair.id <- as.character(pair_attr_CBC$pair.id)
pair_attr_CBC$spp <- as.factor(pair_attr_CBC$spp)
pair_attr_CBC$mid.year <- as.factor(pair_attr_CBC$mid.year)

climate_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + STI + (1|spp) + (1|pair.id), data=pair_attr_CBC)
summary(climate_cbc) ## STI non-significant (p=0.58)
results_table_STI_cbc <- data.frame(summary(climate_cbc)$coefficients[,1:5]) ## 29 species
write.csv(results_table_STI_cbc, file = "../Results/Model_outputs/CBC/average_STI_cbc.csv", row.names=TRUE)

## does this change when mobility is added to the model?
pair_attr_CBC <- merge(pair_attr_CBC, bird_dispersal, by.x="spp", by.y="Species_code", all=FALSE)
pair_attr_CBC <- na.omit(pair_attr_CBC)
length(unique(pair_attr_CBC$spp)) # 21 species 
## split breeding_AM into 2 groups 
pair_attr_CBC$Breeding_AM <- as.numeric(pair_attr_CBC$Breeding_AM)
pair_attr_CBC$Breeding_AM_score2 <- cut(pair_attr_CBC$Breeding_AM, 2, labels=FALSE)
pair_attr_CBC$Breeding_AM_score2 <- as.factor(pair_attr_CBC$Breeding_AM_score2)

climate_mob_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + STI + Breeding_AM_score2 + (1|spp) + (1|pair.id), data=pair_attr_CBC)
summary(climate_mob_cbc) ## STI and mobility non-significant

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

## does this change when mobility is added to the model?
pair_attr_BBS <- merge(pair_attr_BBS, bird_dispersal, by.x="spp", by.y="Species_code", all=FALSE)
pair_attr_BBS <- na.omit(pair_attr_BBS)
length(unique(pair_attr_BBS$spp)) # 17 species 
## split breeding_AM into 2 groups 
pair_attr_BBS$Breeding_AM <- as.numeric(pair_attr_BBS$Breeding_AM)
pair_attr_BBS$Breeding_AM_score2 <- cut(pair_attr_BBS$Breeding_AM, 2, labels=FALSE)
pair_attr_BBS$Breeding_AM_score2 <- as.factor(pair_attr_BBS$Breeding_AM_score2)

climate_mob_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + STI + Breeding_AM_score2 + (1|spp) + (1|pair.id), data=pair_attr_BBS)
summary(climate_mob_bbs) ## STI and mobility non-significant

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
pair_attr_ukbms <- rbind(pair_attr_1985, pair_attr_2000, pair_attr_2012)
length(unique(pair_attr_ukbms$spp)) ## 3 years and 31 species

pair_attr_ukbms$mid.year <- as.factor(pair_attr_ukbms$mid.year)
pair_attr_ukbms$pair.id <- as.character(pair_attr_ukbms$pair.id)
pair_attr_ukbms$spp <- as.factor(pair_attr_ukbms$spp)

climate_ukbms_change <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
summary(climate_ukbms_change)
anova(climate_ukbms_change) ## interaction is very significant (p<0.0001)
results_table_STI_change_ukbms <- data.frame(summary(climate_ukbms_change)$coefficients[,1:5]) ## 31 species
write.csv(results_table_STI_change_ukbms, file = "../Results/Model_outputs/UKBMS/change_STI_ukbms.csv", row.names=TRUE)

## plot mid year and STI interaction result
## split STI into three groups (low, med, high) to plot interaction result
pair_attr_ukbms$STI <- as.numeric(pair_attr_ukbms$STI)
pair_attr_ukbms$STI_score <- cut(pair_attr_ukbms$STI, 2, labels=FALSE)
pair_attr_ukbms$STI_score <- as.factor(pair_attr_ukbms$STI_score)

## run model
climate_ukbms_change2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI_score + (1|pair.id) + (1|spp), data = pair_attr_ukbms)
summary(climate_ukbms_change2)
anova(climate_ukbms_change2) ## significant (p<0.0001)

### predict new data
pair_attr_ukbms$mid.year <- factor(pair_attr_ukbms$mid.year)
newdata_ukbms <- expand.grid(mean_northing=mean(pair_attr_ukbms$mean_northing), distance=mean(pair_attr_ukbms$distance), 
                           renk_hab_sim=mean(pair_attr_ukbms$renk_hab_sim), mid.year=unique(pair_attr_ukbms$mid.year), 
                           pair.id=sample(pair_attr_ukbms$pair.id,10), spp=(unique(pair_attr_ukbms$spp)),
                           STI_score=unique(pair_attr_ukbms$STI_score))
newdata_ukbms$lag0 <- predict(climate_ukbms_change2, newdata=newdata_ukbms, re.form=NA)

mm2 <- model.matrix(terms(climate_ukbms_change2), newdata_ukbms)
pvar2 <- diag(mm2 %*% tcrossprod(vcov(climate_ukbms_change2),mm2))
tvar2 <- pvar2+VarCorr(climate_ukbms_change2)$spp[1]+VarCorr(climate_ukbms_change2)$pair.id[1]
cmult <- 2

newdata_ukbms <- data.frame(
  newdata_ukbms
  , plo = newdata_ukbms$lag0-1.96*sqrt(pvar2)
  , phi = newdata_ukbms$lag0+1.96*sqrt(pvar2)
  , tlo = newdata_ukbms$lag0-1.96*sqrt(tvar2)
  , thi = newdata_ukbms$lag0+1.96*sqrt(tvar2)
)

## run model without STI and species random effect to obtain residuals
STI_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data=pair_attr_ukbms)
pair_attr_ukbms$residuals <- resid(STI_ukbms)
group_by(pair_attr_ukbms, STI_score) %>% summarize(m = mean(lag0)) ## low mean = 0.0256, high mean = 0.0307
## put mean of each group into pair_attr dataframe
pair_attr_ukbms <- ddply(pair_attr_ukbms, "STI_score", transform, STI_mean = mean(lag0))
## add mean to each residual
pair_attr_ukbms$residuals2 <- pair_attr_ukbms$residuals + pair_attr_ukbms$STI_mean

## create new dataframe to calculate mean, SD and SE of residuals for each species
summary_ukbms <- pair_attr_ukbms %>% group_by(spp, STI_score, mid.year) %>% 
  summarise_at(vars(residuals2), funs(mean,std.error))
### 3 species have NA values for standard error (18, 27 and 116)
## these species only have one data point - cannot calculate standard error

## change values
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1984.5"="1985"))
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1999.5"="2000"))
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("2011.5"="2012"))
colnames(newdata_ukbms)[7] <- "STI"
newdata_ukbms$STI <- as.factor(newdata_ukbms$STI)
newdata_ukbms$STI <- revalue(newdata_ukbms$STI, c("1"="Low"))
newdata_ukbms$STI <- revalue(newdata_ukbms$STI, c("2"="High"))
## same for summary dataframe
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("1984.5"="1985"))
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

pd <- position_dodge(0.1)
png("../Graphs/STI/STI_change_predicted_ukbms.png", height = 120, width = 180, units = "mm", res = 300)
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

## remove NAs from pair_attr_ukbms (as one species has no mobility data)
pair_attr_ukbms <- na.omit(pair_attr_ukbms)
length(unique(pair_attr_ukbms$spp)) ## 30 species
## split mobility into two groups
pair_attr_ukbms$mobility_wil <- as.numeric(pair_attr_ukbms$mobility_wil)
pair_attr_ukbms$mobility_score2 <- cut(pair_attr_ukbms$mobility_wil, 2, labels=c("low", "high"))
pair_attr_ukbms$mobility_score2 <- as.factor(pair_attr_ukbms$mobility_score2)

## is mobility still significant with STI in model?
climate_ukbms_change_mob <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*mobility_score2 + STI + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
summary(climate_ukbms_change_mob) ## both years still significant
anova(climate_ukbms_change_mob) ## midyear*mobility interaction still very significant, but STI is not significant
results_table_STI_change_mob_ukbms <- data.frame(summary(climate_ukbms_change_mob)$coefficients[,1:5])
write.csv(results_table_STI_change_mob_ukbms, file = "../Results/Model_outputs/UKBMS/change_mob_STI_ukbms.csv", row.names=TRUE)

## is STI still significant with mobility in model?
climate_ukbms_change_mob2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI + mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
summary(climate_ukbms_change_mob2) ## both years still significant 
anova(climate_ukbms_change_mob2) ## midyear*STI interaction and mobility is significant
results_table_STI_change_mob_ukbms <- data.frame(summary(climate_ukbms_change_mob2)$coefficients[,1:5])
write.csv(results_table_STI_change_mob_ukbms, file = "../Results/Model_outputs/UKBMS/change_STI_mob_ukbms.csv", row.names=TRUE)

## also plot this result ==> MAM would include mobility 
## run model with STI score and mobility
climate_ukbms_change_mob3 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI_score + mobility_score2 + (1|pair.id) + (1|spp), data = pair_attr_ukbms)
summary(climate_ukbms_change_mob3)
anova(climate_ukbms_change_mob3) ## STI_score*midyear and mobility are significant

### predict new data
pair_attr_ukbms$mid.year <- factor(pair_attr_ukbms$mid.year)
newdata_ukbms <- expand.grid(mean_northing=mean(pair_attr_ukbms$mean_northing), distance=mean(pair_attr_ukbms$distance), 
                             renk_hab_sim=mean(pair_attr_ukbms$renk_hab_sim), mid.year=unique(pair_attr_ukbms$mid.year), 
                             pair.id=sample(pair_attr_ukbms$pair.id,10), spp=(unique(pair_attr_ukbms$spp)), 
                             mobility_score2=sample(pair_attr_ukbms$mobility_score2,1), STI_score=unique(pair_attr_ukbms$STI_score))
newdata_ukbms$lag0 <- predict(climate_ukbms_change_mob3, newdata=newdata_ukbms, re.form=NA)

mm2 <- model.matrix(terms(climate_ukbms_change_mob3), newdata_ukbms)
pvar2 <- diag(mm2 %*% tcrossprod(vcov(climate_ukbms_change_mob3),mm2))
tvar2 <- pvar2+VarCorr(climate_ukbms_change2)$spp[1]+VarCorr(climate_ukbms_change2)$pair.id[1]
cmult <- 2

newdata_ukbms <- data.frame(
  newdata_ukbms
  , plo = newdata_ukbms$lag0-1.96*sqrt(pvar2)
  , phi = newdata_ukbms$lag0+1.96*sqrt(pvar2)
  , tlo = newdata_ukbms$lag0-1.96*sqrt(tvar2)
  , thi = newdata_ukbms$lag0+1.96*sqrt(tvar2)
)

## run model without STI and species random effect to obtain residuals
STI_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + mobility_score2 + (1|pair.id), data=pair_attr_ukbms)
pair_attr_ukbms$residuals <- resid(STI_ukbms)
group_by(pair_attr_ukbms, STI_score) %>% summarize(m = mean(lag0)) ## low mean = 0.0256, high mean = 0.0307
## put mean of each group into pair_attr dataframe
pair_attr_ukbms <- ddply(pair_attr_ukbms, "STI_score", transform, STI_mean = mean(lag0))
## add mean to each residual
pair_attr_ukbms$residuals2 <- pair_attr_ukbms$residuals + pair_attr_ukbms$STI_mean

## create new dataframe to calculate mean, SD and SE of residuals for each species
summary_ukbms <- pair_attr_ukbms %>% group_by(spp, STI_score, mid.year) %>% 
  summarise_at(vars(residuals2), funs(mean,std.error))
### 3 species have NA values for standard error (18, 27 and 116)
## these species only have one data point - cannot calculate standard error

## change values
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1984.5"="1985"))
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("1999.5"="2000"))
newdata_ukbms$mid.year <- revalue(newdata_ukbms$mid.year, c("2011.5"="2012"))
colnames(newdata_ukbms)[8] <- "STI"
newdata_ukbms$STI <- as.factor(newdata_ukbms$STI)
newdata_ukbms$STI <- revalue(newdata_ukbms$STI, c("1"="Low"))
newdata_ukbms$STI <- revalue(newdata_ukbms$STI, c("2"="High"))
## same for summary dataframe
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("1984.5"="1985"))
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

pd <- position_dodge(0.1)
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
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC_no_zeros2.csv", header=TRUE) 
pair_attr_CBC <- merge(pair_attr_CBC, bird_STI, by.x="spp", by.y="species_code") ## 31 species
length(unique(pair_attr_CBC$spp)) # 29 species
pair_attr_CBC$pair.id <- as.character(pair_attr_CBC$pair.id)
pair_attr_CBC$spp <- as.factor(pair_attr_CBC$spp)
pair_attr_CBC$mid.year <- as.factor(pair_attr_CBC$mid.year)

pair_attr_cbc_1985 <- pair_attr_CBC[pair_attr_CBC$mid.year==1984.5,]
pair_attr_cbc_1996 <- pair_attr_CBC[pair_attr_CBC$mid.year==1995.5,]
pair_attr_cbc <- rbind(pair_attr_cbc_1985, pair_attr_cbc_1996)
length(unique(pair_attr_cbc$spp)) ## 29 species

pair_attr_cbc$mid.year <- as.factor(pair_attr_cbc$mid.year)
pair_attr_cbc$pair.id <- as.character(pair_attr_cbc$pair.id)
pair_attr_cbc$spp <- as.factor(pair_attr_cbc$spp)

###### run model with STI and year interaction
climate_cbc_change <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*STI + (1|pair.id) + (1|spp), data = pair_attr_cbc)
summary(climate_cbc_change)
anova(climate_cbc_change) ## interaction is non-significant (p=0.84)
results_table_STI_change_cbc <- data.frame(summary(climate_cbc_change)$coefficients[,1:5]) ## 29 species
write.csv(results_table_STI_change_cbc, file = "../Results/Model_outputs/CBC/change_STI_cbc.csv", row.names=TRUE)

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

## split STI into two groups (low and high) to plot midyear*STI interaction result
pair_attr_bbs$STI <- as.numeric(pair_attr_bbs$STI)
pair_attr_bbs$STI_score <- cut(pair_attr_bbs$STI, 2, labels=FALSE)
pair_attr_bbs$STI_score <- as.factor(pair_attr_bbs$STI_score)

## run model
climate_bbs_change2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI_score + (1|pair.id) + (1|spp), data = pair_attr_bbs)
summary(climate_bbs_change2)
anova(climate_bbs_change2) ## not significant 

### predict new data
pair_attr_bbs$mid.year <- factor(pair_attr_bbs$mid.year)
newdata_bbs <- expand.grid(mean_northing=mean(pair_attr_bbs$mean_northing), distance=mean(pair_attr_bbs$distance), 
                           renk_hab_sim=mean(pair_attr_bbs$renk_hab_sim), mid.year=unique(pair_attr_bbs$mid.year), 
                           pair.id=sample(pair_attr_bbs$pair.id,10), spp=sample(pair_attr_bbs$spp,10),
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
ggplot(summary_bbs, aes(x = mid.year, y = mean, group=STI)) +
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
