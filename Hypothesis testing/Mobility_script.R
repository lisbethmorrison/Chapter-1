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
pair_attr <- read.csv("../Data/Butterfly_sync_data/pop_climate_synchrony.csv", header=TRUE) # butterfly pair attribute data
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pop_climate_synchrony_CBC.csv", header=TRUE) # CBC pair attribute data
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pop_climate_synchrony_BBS.csv", header=TRUE) # BBS pair attribute data
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

str(pair_attr)
pair_attr$mid.year <- as.factor(pair_attr$mid.year)
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$mobility_wil <- as.numeric(pair_attr$mobility_wil)

# ## scale mobility_wil values
# pair_attr[26] <- lapply(pair_attr[26], function(pair_attr) c(scale(pair_attr, center = TRUE, scale = TRUE))) 

### run mobility as continuous varaible
mobility_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain + 
                         winter_temp + autumn_temp + spring_temp + summer_temp + mobility_wil + (1|spp) + (1|pair.id), data=pair_attr)
summary(mobility_model)
anova(mobility_model) ## significant (p=0.01)
results_table_mob_ukbms <- data.frame(summary(mobility_model)$coefficients[,1:5]) ## 31 species
write.csv(results_table_mob_ukbms, file = "../Results/Model_outputs/UKBMS/average_mob_ukbms.csv", row.names=TRUE)

# ## run model with 2 groups (high and low) of mobility 
# pair_attr$mobility_wil <- as.numeric(pair_attr$mobility_wil)
# pair_attr$mobility_score2 <- cut(pair_attr$mobility_wil, 2, labels=c("low", "high"))
# pair_attr$mobility_score2 <- as.factor(pair_attr$mobility_score2)
# 
# mobility_model4 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + mobility_score2 + (1|spp) + (1|pair.id), data=pair_attr)
# summary(mobility_model4)
# anova(mobility_model4)
# ## significant (p=0.011)
# ## mobility score 2 has higher average synchrony than score 1
# r.squaredGLMM(mobility_model4)
# ## R^2 marginal (fixed effects) = 0.029, R^2 conditional (random effects) = 0.224
# results_table_mob_ukbms <- data.frame(summary(mobility_model4)$coefficients[,1:5]) ## 31 species
# write.csv(results_table_mob_ukbms, file = "../Results/Model_outputs/UKBMS/average_mob_ukbms.csv", row.names=TRUE)

##################################
#### Mobility model for CBC ####
##################################

## merge with bird dispersal data
length(unique(bird_dispersal$Species_code)) ## only 23 species have dispersal data

## merge pair_attr with bird dispersal data
pair_attr_CBC <- merge(pair_attr_CBC, bird_dispersal, by.x="spp", by.y="Species_code", all=FALSE)
pair_attr_CBC <- na.omit(pair_attr_CBC) ## remove NAs (one species == redstart)
length(unique(pair_attr_CBC$spp)) ## 20 species

## scale dispersal data
pair_attr_CBC[42] <- lapply(pair_attr_CBC[42], function(pair_attr_CBC) c(scale(pair_attr_CBC, center = TRUE, scale = TRUE))) 

## breeding arithmetic mean dispersal
dispersal_model_cbc2 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + summer_temp + Breeding_AM + (1|spp) + (1|pair.id), data=pair_attr_CBC)
summary(dispersal_model_cbc2)
anova(dispersal_model_cbc2)
## not significant 
results_table_mob_cbc <- data.frame(summary(dispersal_model_cbc2)$coefficients[,1:5]) ## 20 species
write.csv(results_table_mob_cbc, file = "../Results/Model_outputs/CBC/average_mob_cbc.csv", row.names=TRUE)

# ## split breeding_AM into 3 groups and run model again
# pair_attr_CBC$Breeding_AM <- as.numeric(pair_attr_CBC$Breeding_AM)
# pair_attr_CBC$Breeding_AM_score <- cut(pair_attr_CBC$Breeding_AM, 3, labels=FALSE)
# pair_attr_CBC$Breeding_AM_score <- as.factor(pair_attr_CBC$Breeding_AM_score)
# 
# dispersal_model_cbc3 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + Breeding_AM_score + (1|spp) + (1|pair.id), data=pair_attr_CBC)
# summary(dispersal_model_cbc3)
# anova(dispersal_model_cbc3)
# ## still not significant
# 
# ## split breeding_AM into 2 groups and compare
# pair_attr_CBC$Breeding_AM <- as.numeric(pair_attr_CBC$Breeding_AM)
# pair_attr_CBC$Breeding_AM_score2 <- cut(pair_attr_CBC$Breeding_AM, 2, labels=c("low", "high"))
# pair_attr_CBC$Breeding_AM_score2 <- as.factor(pair_attr_CBC$Breeding_AM_score2)
# 
# pair_attr_CBC <- droplevels(pair_attr_CBC)
# 
# dispersal_model_cbc4 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + Breeding_AM_score2 + (1|family) + (1|spp) + (1|pair.id), data=pair_attr_CBC)
# ## no errors
# summary(dispersal_model_cbc4)
# anova(dispersal_model_cbc4) ## dispersal score is non-significant (p=0.63)
# results_table_mob_cbc <- data.frame(summary(dispersal_model_cbc4)$coefficients[,1:5]) ## 20 species
# write.csv(results_table_mob_cbc, file = "../Results/Model_outputs/CBC/average_mob_cbc.csv", row.names=TRUE)

################################
#### Mobility model for BBS ####
################################

## merge pair_attr with bird dispersal data
pair_attr_BBS <- merge(pair_attr_BBS, bird_dispersal, by.x="spp", by.y="Species_code", all=FALSE)
## breeding arithmetic mean dispersal
## remove NAs (one species == redstart)
pair_attr_BBS <- na.omit(pair_attr_BBS)
length(unique(pair_attr_BBS$spp)) # 17spp

## scale dispersal data
pair_attr_BBS[39] <- lapply(pair_attr_BBS[39], function(pair_attr_BBS) c(scale(pair_attr_BBS, center = TRUE, scale = TRUE))) 

dispersal_model_bbs2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + spring_rain + autumn_rain +
                               Breeding_AM + (1|spp) + (1|pair.id), data=pair_attr_BBS)
summary(dispersal_model_bbs2)
anova(dispersal_model_bbs2)
## not significant 
results_table_mob_bbs <- data.frame(summary(dispersal_model_bbs2)$coefficients[,1:5]) ## 17 species
write.csv(results_table_mob_bbs, file = "../Results/Model_outputs/BBS/average_mob_bbs.csv", row.names=TRUE)

# ## split breeding_AM into 3 groups and run model again
# pair_attr_BBS$Breeding_AM <- as.numeric(pair_attr_BBS$Breeding_AM)
# pair_attr_BBS$Breeding_AM_score <- cut(pair_attr_BBS$Breeding_AM, 3, labels=FALSE)
# pair_attr_BBS$Breeding_AM_score <- as.factor(pair_attr_BBS$Breeding_AM_score)
# 
# dispersal_model_bbs3 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Breeding_AM_score + (1|spp) + (1|pair.id), data=pair_attr_BBS)
# summary(dispersal_model_bbs3)
# anova(dispersal_model_bbs3)
# ## still not significant
# 
# ## split breeding_AM into 2 groups and compare
# pair_attr_BBS$Breeding_AM <- as.numeric(pair_attr_BBS$Breeding_AM)
# pair_attr_BBS$Breeding_AM_score2 <- cut(pair_attr_BBS$Breeding_AM, 2, labels=FALSE)
# pair_attr_BBS$Breeding_AM_score2 <- as.factor(pair_attr_BBS$Breeding_AM_score2)
# 
# dispersal_model_bbs4 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + Breeding_AM_score2 + (1|spp) + (1|pair.id), data=pair_attr_BBS)
# summary(dispersal_model_bbs4)
# anova(dispersal_model_bbs4)
# ## still not significant
# results_table_mob_bbs <- data.frame(summary(dispersal_model_bbs4)$coefficients[,1:5]) ## 17 species
# write.csv(results_table_mob_bbs, file = "../Results/Model_outputs/BBS/average_mob_bbs.csv", row.names=TRUE)


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

## rescale climate variables
pair_attr_early[29:36] <- lapply(pair_attr_early[29:36], function(pair_attr_early) c(scale(pair_attr_early, center = TRUE, scale = TRUE))) 
pair_attr_late[29:36] <- lapply(pair_attr_late[29:36], function(pair_attr_late) c(scale(pair_attr_late, center = TRUE, scale = TRUE))) 

## full model with interaction between mid year and mobility EARLY
mobility_model5 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + summer_rain + 
                          winter_temp + autumn_temp + spring_temp + summer_temp + mid.year*mobility_wil + (1|spp) + (1|pair.id), data=pair_attr_early)
summary(mobility_model5)
anova(mobility_model5)
## interaction non-significant (p=0.65)
results_table_mobility_ukbms <- data.frame(summary(mobility_model5)$coefficients[,1:5]) ## 31 species
write.csv(results_table_mobility_ukbms, file = "../Results/Model_outputs/UKBMS/change_mobility_ukbms_85_00.csv", row.names=TRUE)

## full model with interaction between mid year and mobility LATE
mobility_model6 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + summer_rain + 
                          winter_temp + autumn_temp + spring_temp + summer_temp + mid.year*mobility_wil + (1|spp) + (1|pair.id), data=pair_attr_late)
summary(mobility_model6)
anova(mobility_model6)
## overall interaction very significant (p<0.0001)
results_table_mobility_ukbms <- data.frame(summary(mobility_model6)$coefficients[,1:5]) ## 31 species
write.csv(results_table_mobility_ukbms, file = "../Results/Model_outputs/UKBMS/change_mobility_ukbms_00_12.csv", row.names=TRUE)


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

## scale dispersal data
pair_attr_cbc[42] <- lapply(pair_attr_cbc[42], function(pair_attr_cbc) c(scale(pair_attr_cbc, center = TRUE, scale = TRUE))) 

###### run model with strategy and year interaction
dispersal_model_cbc1 <- lmer(lag0 ~ mean_northing + distance + hab_sim + summer_temp + mid.year*Breeding_AM + (1|pair.id) + (1|spp), data = pair_attr_cbc)
summary(dispersal_model_cbc1)
anova(dispersal_model_cbc1)
## mid.year*dispersal interaction is NON-significant overall (0.19)
results_table_dispersal_cbc <- data.frame(summary(dispersal_model_cbc1)$coefficients[,1:5]) ## 20 species
write.csv(results_table_dispersal_cbc, file = "../Results/Model_outputs/CBC/change_mob_cbc.csv", row.names=TRUE)

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

pair_attr_bbs$mid.year <- as.numeric(pair_attr_bbs$mid.year)
pair_attr_bbs$pair.id <- as.character(pair_attr_bbs$pair.id)
pair_attr_bbs$spp <- as.factor(pair_attr_bbs$spp)

pair_attr_bbs[39] <- lapply(pair_attr_bbs[39], function(pair_attr_bbs) c(scale(pair_attr_bbs, center = TRUE, scale = TRUE))) 

###### run model with mobility and year interaction
dispersal_model_bbs1 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + spring_rain + autumn_rain + 
                               mid.year*Breeding_AM + (1|pair.id) + (1|spp), data = pair_attr_bbs)
summary(dispersal_model_bbs1)
## intercept is mid.year 1999
anova(dispersal_model_bbs1)
## mid.year*dispersal interaction is significant overall
results_table_dispersal_bbs <- data.frame(summary(dispersal_model_bbs1)$coefficients[,1:5])
write.csv(results_table_dispersal_bbs, file = "../Results/Model_outputs/BBS/change_mob_bbs.csv", row.names=TRUE)

# ## run model with three mobility groups (low, med, high)
# pair_attr_bbs$Breeding_AM <- as.numeric(pair_attr_bbs$Breeding_AM)
# pair_attr_bbs$Breeding_AM_score <- cut(pair_attr_bbs$Breeding_AM, 3, labels=FALSE)
# pair_attr_bbs$Breeding_AM_score <- as.factor(pair_attr_bbs$Breeding_AM_score)
# 
# ## remove NAs (one species == redstart)
# pair_attr_bbs <- na.omit(pair_attr_bbs)
# summary(pair_attr_bbs)  
# 
# dispersal_model_bbs2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*Breeding_AM_score + (1|spp) + (1|pair.id), data=pair_attr_bbs)
# summary(dispersal_model_bbs2)
# anova(dispersal_model_bbs2)
# ## interaction is very significant (p<0.000001)
# 
# ## run model with 2 mobility groups (low and high)
# pair_attr_bbs$Breeding_AM <- as.numeric(pair_attr_bbs$Breeding_AM)
# pair_attr_bbs$Breeding_AM_score2 <- cut(pair_attr_bbs$Breeding_AM, 2, labels=FALSE)
# pair_attr_bbs$Breeding_AM_score2 <- as.factor(pair_attr_bbs$Breeding_AM_score2)
# 
# start_time <- Sys.time()
# dispersal_model_bbs3 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*Breeding_AM_score2 + (1|spp) + (1|pair.id), data=pair_attr_bbs)
# end_time <- Sys.time()
# end_time - start_time ## 34 seconds
# 
# summary(dispersal_model_bbs3)
# anova(dispersal_model_bbs3)
# ## interaction is very significant (p<0.000001)
# ## save model output
# results_table_dispersal_bbs <- data.frame(summary(dispersal_model_bbs3)$coefficients[,1:5])
# write.csv(results_table_dispersal_bbs, file = "../Results/Model_outputs/BBS/change_mob_bbs.csv", row.names=TRUE)

############################################# PLOT REUSLT ############################################# 

## run model without mobility for each species and extract the slope (mid.year) coefficient for each
results_table<-NULL
for (i in unique(pair_attr_bbs$spp)){
  print(i)
  ## if statement: if if the mixed effects model doesn't throw an error (levels of pairID), then run an lmer, else run a linear model 
  loadError=F
  a=try({mob_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + spring_rain + autumn_rain +
                              (1|pair.id), data=pair_attr_bbs[pair_attr_bbs$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    mob_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + spring_rain + autumn_rain
                       + (1|pair.id), data=pair_attr_bbs[pair_attr_bbs$spp==i,])
  }else{
    mob_bbs <- lm(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + spring_rain + autumn_rain,
                     data=pair_attr_bbs[pair_attr_bbs$spp==i,])
  }
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(mob_bbs)$coefficients[,1:2],i)
  results_table_temp$parameter <- paste(row.names(results_table_temp))
  rownames(results_table_temp) <- 1:nrow(results_table_temp)
  results_table_temp <- results_table_temp[grep("mid.year", results_table_temp$parameter),]
  results_table<-rbind(results_table,results_table_temp)
}

## merge in this dataframe with specialism info (spec/gen)
mobility <- pair_attr_bbs[,c(1,39)]
mobility <- unique(mobility)
results_table <- merge(results_table, mobility, by.x="i", by.y="spp")
## remove parameter column
results_table <- results_table[,-c(4)]
## re-name some columns
names(results_table) <- c("species", "slope", "SE", "Mobility")


########## calculate slope and SE from main model (with mobility as x axis)
### plot points and SE ontop of raw data
results_table_dispersal_bbs <- results_table_dispersal_bbs[,-c(3:5)]
## leave rows with year and interaction
results_table_dispersal_bbs$Mobility <- paste(row.names(results_table_dispersal_bbs))
rownames(results_table_dispersal_bbs) <- 1:nrow(results_table_dispersal_bbs)
results_table_dispersal_bbs <- results_table_dispersal_bbs[c(8,10),]

mid.year_slope <- results_table_dispersal_bbs[1,1]
interaction_slope <- results_table_dispersal_bbs[1,2]

at.x2<-unique(results_table$Mobility)
slopes <- mid.year_slope + interaction_slope * at.x2
slopes

estvar<-vcov(dispersal_model_bbs1); dispersal_model_bbs1.vcov<-as.data.frame(as.matrix(estvar))
var.b1<-dispersal_model_bbs1.vcov["mid.year","mid.year"]
var.b3<-dispersal_model_bbs1.vcov["mid.year:Breeding_AM","mid.year:Breeding_AM"]
cov.b1.b3<-dispersal_model_bbs1.vcov["mid.year","mid.year:Breeding_AM"]

SEs <- rep(NA, length(at.x2))
for (i in 1:length(at.x2)){
  j <- at.x2[i]  
  SEs[i] <- sqrt(var.b1 + var.b3 * j^2 + 2*j* cov.b1.b3)
}
model_summary <- data.frame(Mobility=at.x2, slope=slopes, SE=SEs)

png("../Graphs/Mobility/Mobility_change_bbs_99_12.png", height = 150, width = 180, units = "mm", res = 300)
ggplot(mapping=aes(x=Mobility, y=slope, ymin=slope-SE, ymax=slope+SE)) +
  geom_point(data=results_table, size=2, color="lightgrey") +
  geom_errorbar(data=results_table, color="lightgrey", width=0.1) +
  geom_ribbon(data=model_summary, fill = "grey70") +
  geom_line(data=model_summary, size=1) +
  labs(x="Standardised dispersal distance", y="Change in population synchrony 1999-2012") +
  scale_x_continuous(breaks = seq(-2,3,0.5)) +
  theme_bw() +
  theme(text = element_text(size = 12), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black"))
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

