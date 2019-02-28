#############################################################################
## Title: Position in climatic space and synchrony analysis
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
options(scipen=999)

## read in synchrony data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE) 
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC.csv", header=TRUE) 
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE) 
bird_STI <- read.csv("../Data/birds_STI.csv", header=TRUE)
butterfly_STI <- read.csv("../Data/butterflies_STI.csv", header=TRUE)

############### AVERAGE SYNCHORNY #################

#### UKBMS
## merge datasets
pair_attr <- merge(pair_attr, butterfly_STI, by.x="spp", by.y="species_code") ## 32 species
## No STI data for Essex skipper
str(pair_attr)

pair_attr$mid.year <- as.integer(pair_attr$mid.year)
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)

start_time <- Sys.time()
climate_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + STI + (1|spp) + (1|pair.id), data=pair_attr)
end_time <- Sys.time()
end_time - start_time ## 4.5mins
summary(climate_ukbms) ## STI non-significant (p=0.316)

##### CBC
## merge datasets
pair_attr_CBC <- merge(pair_attr_CBC, bird_STI, by.x="spp", by.y="species_code") ## 31 species
str(pair_attr_CBC)

pair_attr_CBC$mid.year <- as.integer(pair_attr_CBC$mid.year)
pair_attr_CBC$pair.id <- as.character(pair_attr_CBC$pair.id)
pair_attr_CBC$spp <- as.factor(pair_attr_CBC$spp)

start_time <- Sys.time()
climate_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + STI + (1|spp) + (1|pair.id), data=pair_attr_CBC)
end_time <- Sys.time()
end_time - start_time ## 24.36 seconds
summary(climate_cbc) ## STI non-significant (p=0.098)

##### BBS
## merge datasets
pair_attr_BBS <- merge(pair_attr_BBS, bird_STI, by.x="spp", by.y="species_code") ## 24 species
str(pair_attr_BBS)

pair_attr_BBS$mid.year <- as.integer(pair_attr_BBS$mid.year)
pair_attr_BBS$pair.id <- as.character(pair_attr_BBS$pair.id)
pair_attr_BBS$spp <- as.factor(pair_attr_BBS$spp)

start_time <- Sys.time()
climate_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + STI + (1|spp) + (1|pair.id), data=pair_attr_BBS)
end_time <- Sys.time()
end_time - start_time ## 24.36 seconds
summary(climate_bbs) ## STI non-significant (p=0.882)
############### CHANGE IN SYNCHORNY #################

#### UKBMS
## subset to only look at mid years 1985, 2000 and 2012
pair_attr_1985 <- pair_attr[pair_attr$mid.year==1984.5,]
pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]
pair_attr_ukbms <- rbind(pair_attr_1985, pair_attr_2000, pair_attr_2012) ## 3 years and 32 species

pair_attr_ukbms$mid.year <- as.factor(pair_attr_ukbms$mid.year)
pair_attr_ukbms$pair.id <- as.character(pair_attr_ukbms$pair.id)
pair_attr_ukbms$spp <- as.factor(pair_attr_ukbms$spp)

start_time <- Sys.time()
climate_ukbms_change <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI + (1|spp) + (1|pair.id), data=pair_attr_ukbms)
end_time <- Sys.time()
end_time - start_time ## 14.18secs
summary(climate_ukbms_change)
anova(climate_ukbms_change) ## interaction is very significant (p<0.0001)

## predict new data to plot graph
newdata_ukbms <- expand.grid(mean_northing=mean(pair_attr_ukbms$mean_northing), distance=mean(pair_attr_ukbms$distance), 
                           renk_hab_sim=mean(pair_attr_ukbms$renk_hab_sim), pair.id=sample(pair_attr_ukbms$pair.id,10),
                           mid.year=unique(pair_attr_ukbms$mid.year), spp=unique(pair_attr_ukbms$spp),
                           STI=unique(pair_attr_ukbms$STI))
newdata_ukbms$lag0 <- predict(climate_ukbms_change, newdata=newdata_ukbms, re.form=NA)

mm2 <- model.matrix(terms(climate_ukbms_change), newdata_ukbms)
pvar2 <- diag(mm2 %*% tcrossprod(vcov(climate_ukbms_change),mm2))
tvar2 <- pvar2+VarCorr(climate_ukbms_change)$spp[1]+VarCorr(climate_ukbms_change)$pair.id[1]
cmult <- 2

newdata_ukbms <- data.frame(
  newdata_ukbms
  , plo = newdata_ukbms$lag0-1.96*sqrt(pvar2)
  , phi = newdata_ukbms$lag0+1.96*sqrt(pvar2)
  , tlo = newdata_ukbms$lag0-1.96*sqrt(tvar2)
  , thi = newdata_ukbms$lag0+1.96*sqrt(tvar2)
)

## model without abundance and species random effect to get residuals for graph
model_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data = pair_attr_ukbms)
pair_attr_ukbms$residuals <- resid(model_ukbms)

group_by(pair_attr_ukbms, STI) %>% summarize(m = mean(lag0)) 
## put mean of each group into pair_attr dataframe
pair_attr_ukbms <- ddply(pair_attr_ukbms, "STI", transform, STI_mean = mean(lag0))
## add mean to each residual
pair_attr_ukbms$residuals2 <- pair_attr_ukbms$residuals + pair_attr_ukbms$STI_mean

## create new dataframe to calculate mean, SD and SE of residuals for each species
summary_ukbms <- pair_attr_ukbms %>% group_by(spp, STI, mid.year) %>% 
  summarise_at(vars(residuals2), funs(mean,std.error,sd))

## plot graph with raw data residuals (+SE error bars) and fitted line
png("../Graphs/Abundance/Common_average_predicted_cbc.png", height = 100, width = 110, units = "mm", res = 300)
ggplot(summary_ukbms, aes(x = STI, y = mean, colour=mid.year)) +
  geom_point(size = 2) +
  #geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), width=0.2, colour="grey66") +
  #geom_line(data=newdata_cbc, aes(x=pop_estimate_log, y=lag0), colour="black", lwd=1) +
  #labs(x="(log) Average population abundance", y="Population synchrony") +
  theme_bw() +
  theme(text = element_text(size = 8), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()




















#### CBC
pair_attr_cbc_1985 <- pair_attr_CBC[pair_attr_CBC$mid.year==1984.5,]
pair_attr_cbc_1996 <- pair_attr_CBC[pair_attr_CBC$mid.year==1995.5,]
pair_attr_cbc <- rbind(pair_attr_cbc_1985, pair_attr_cbc_1996)

pair_attr_cbc$mid.year <- as.factor(pair_attr_cbc$mid.year)
pair_attr_cbc$pair.id <- as.character(pair_attr_cbc$pair.id)
pair_attr_cbc$spp <- as.factor(pair_attr_cbc$spp)

###### run model with STI and year interaction
start_time <- Sys.time()
climate_cbc_change <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*STI + (1|pair.id) + (1|spp), data = pair_attr_cbc)
end_time <- Sys.time()
end_time - start_time ## 5.886secs
summary(climate_cbc_change)
anova(climate_cbc_change) ## interaction is non-significant (p=0.9763)

#### BBS
pair_attr_bbs_1999 <- pair_attr_BBS[pair_attr_BBS$mid.year==1998.5,]
pair_attr_bbs_2012 <- pair_attr_BBS[pair_attr_BBS$mid.year==2011.5,]
pair_attr_bbs <- rbind(pair_attr_bbs_1999, pair_attr_bbs_2012)

pair_attr_bbs$mid.year <- as.factor(pair_attr_bbs$mid.year)
pair_attr_bbs$pair.id <- as.character(pair_attr_bbs$pair.id)
pair_attr_bbs$spp <- as.factor(pair_attr_bbs$spp)

###### run model with STI and year interaction
start_time <- Sys.time()
climate_bbs_change <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*STI + (1|pair.id) + (1|spp), data = pair_attr_bbs)
end_time <- Sys.time()
end_time - start_time ## 36.953
summary(climate_bbs_change)
anova(climate_bbs_change) ## interaction is significant (p=0.0332)




