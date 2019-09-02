#############################################################################
## Title: Testing for significant difference between gen/spec of woodland birds and butterflies
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: April 2018
#############################################################################

rm(list=ls()) # clear R
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggsignif)
library(plyr)
library(dplyr)
library(plotrix)
options(scipen=999)

## read data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pop_climate_synchrony.csv", header=TRUE) # butterfly pair attribute data
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pop_climate_synchrony_CBC.csv", header=TRUE) # CBC pair attribute data
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pop_climate_synchrony_BBS.csv", header=TRUE) # BBS pair attribute data

##############################################
#### specialism model for ALL butterflies ####
##############################################

pair_attr$mid.year <- as.factor(pair_attr$mid.year)
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$specialism <- as.factor(pair_attr$specialism)

spec_model_ukbms1 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + summer_rain + winter_temp +
                            autumn_temp + spring_temp + summer_temp + mid.year + specialism + (1|pair.id) + (1|spp), data = pair_attr)
summary(spec_model_ukbms1)
anova(spec_model_ukbms1)
## non-significant
results_table_average_spec <- data.frame(summary(spec_model_ukbms1)$coefficients[,1:5]) ## 32 species
write.csv(results_table_average_spec, file = "../Results/Model_outputs/UKBMS/average_spec_ukbms.csv", row.names=TRUE) # 32 species

###############################################
########### Same for woodland birds ###########
###############################################

## CBC ##
str(pair_attr_CBC)
pair_attr_CBC$mid.year <- as.factor(pair_attr_CBC$mid.year)
pair_attr_CBC$pair.id <- as.character(pair_attr_CBC$pair.id)
pair_attr_CBC$spp <- as.factor(pair_attr_CBC$spp)

## run model with specialism as fixed categorical variable
strategy_model_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + summer_temp + specialism + (1|family) + (1|pair.id) + (1|spp), data = pair_attr_CBC)
## no errors
summary(strategy_model_cbc)
anova(strategy_model_cbc) ## specialism is non-significant (p=0.09)
results_table_average_spec <- data.frame(summary(strategy_model_cbc)$coefficients[,1:5]) ## 26 species
write.csv(results_table_average_spec, file = "../Results/Model_outputs/CBC/average_spec_cbc.csv", row.names=TRUE) # 26 species

## BBS ##
str(pair_attr_BBS)
pair_attr_BBS$mid.year <- as.factor(pair_attr_BBS$mid.year)
pair_attr_BBS$pair.id <- as.character(pair_attr_BBS$pair.id)
pair_attr_BBS$spp <- as.factor(pair_attr_BBS$spp)

## run model with strategy as fixed categorical variable
strategy_model_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + 
                             specialism + (1|pair.id) + (1|spp), data = pair_attr_BBS)
summary(strategy_model_bbs)
anova(strategy_model_bbs)
## non-significant (p=0.83)
results_table_average_spec <- data.frame(summary(strategy_model_bbs)$coefficients[,1:5]) 
write.csv(results_table_average_spec, file = "../Results/Model_outputs/BBS/average_spec_bbs.csv", row.names=TRUE) # 24 species

#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################

### Now test for change in synchony against specialism 

############################### interaction in model #################################
pair_attr <- read.csv("../Data/Butterfly_sync_data/pop_climate_synchrony.csv", header=TRUE)
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pop_climate_synchrony_CBC.csv", header=TRUE)
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pop_climate_synchrony_BBS.csv", header=TRUE)

#### subset only years 1985, 2000 and 2012
pair_attr_1985 <- pair_attr[pair_attr$mid.year==1984.5,]
pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]
pair_attr_early <- rbind(pair_attr_1985, pair_attr_2000) ## 3 years and 33 species
pair_attr_late <- rbind(pair_attr_2000, pair_attr_2012) ## 3 years and 33 species
pair_attr_early <- droplevels(pair_attr_early)
pair_attr_late <- droplevels(pair_attr_late)

pair_attr_early$mid.year <- as.numeric(pair_attr_early$mid.year)
pair_attr_early$pair.id <- as.character(pair_attr_early$pair.id)
pair_attr_early$spp <- as.factor(pair_attr_early$spp)
pair_attr_early$specialism <- as.factor(pair_attr_early$specialism)

pair_attr_late$mid.year <- as.numeric(pair_attr_late$mid.year)
pair_attr_late$pair.id <- as.character(pair_attr_late$pair.id)
pair_attr_late$spp <- as.factor(pair_attr_late$spp)
pair_attr_late$specialism <- as.factor(pair_attr_late$specialism)

###### run model with strategy and year interaction [EARLY]
spec_model_ukbms2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*specialism + winter_rain + autumn_rain + spring_rain + summer_rain + 
                            winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr_early)
summary(spec_model_ukbms2) ## interaction is significant
anova(spec_model_ukbms2)
## mid.year*strategy interaction is significant overall
## save model output
results_table_spec_ukbms <- data.frame(summary(spec_model_ukbms2)$coefficients[,1:5]) ## 32 species
write.csv(results_table_spec_ukbms, file = "../Results/Model_outputs/UKBMS/change_spec_ukbms_85_00.csv", row.names=TRUE)

###### run model with strategy and year interaction [LATE]
spec_model_ukbms3 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*specialism + winter_rain + autumn_rain + spring_rain + summer_rain + 
                            winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id) + (1|spp), data = pair_attr_late)
summary(spec_model_ukbms3) ## interaction is significant
anova(spec_model_ukbms3)
## mid.year*strategy interaction is significant overall
## save model output
results_table_spec_ukbms2 <- data.frame(summary(spec_model_ukbms3)$coefficients[,1:5]) ## 32 species
write.csv(results_table_spec_ukbms2, file = "../Results/Model_outputs/UKBMS/change_spec_ukbms_00_12.csv", row.names=TRUE)

############################################# PLOT REUSLTS ############################################# 

############### EARLY MODEL ############### 
## run model without specialism for each species and extract the slope (mid.year) coefficient for each
results_table<-NULL
for (i in unique(pair_attr_early$spp)){
  print(i)
  ## if statement: if if the mixed effects model doesn't throw an error (levels of pairID), then run an lmer, else run a linear model 
  loadError=F
  a=try({spec_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain + 
                              winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id), data=pair_attr_early[pair_attr_early$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain + 
                         winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id), data=pair_attr_early[pair_attr_early$spp==i,])
  }else{
    spec_ukbms <- lm(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain + 
                         winter_temp + autumn_temp + spring_temp + summer_temp, data=pair_attr_early[pair_attr_early$spp==i,])
  }

    ### save and plot the results ###
    results_table_temp <- data.frame(summary(spec_ukbms)$coefficients[,1:2],i)
    results_table_temp$parameter <- paste(row.names(results_table_temp))
    rownames(results_table_temp) <- 1:nrow(results_table_temp)
    results_table_temp <- results_table_temp[grep("mid.year", results_table_temp$parameter),]
    results_table<-rbind(results_table,results_table_temp)
}

## merge in this dataframe with specialism info (spec/gen)
specialism <- pair_attr_early[,c(14,26)]
specialism <- unique(specialism)
results_table <- merge(results_table, specialism, by.x="i", by.y="spp")
## remove parameter column
results_table <- results_table[,-c(4)]
## re-name some columns
names(results_table) <- c("species", "slope", "SE", "Specialism")
results_table$Specialism <- revalue(results_table$Specialism, c("specialist"="Specialist"))
results_table$Specialism <- revalue(results_table$Specialism, c("wider.countryside"="Generalist"))

########## calculate slope and SE from main model (with specialism)
### plot points and SE ontop of raw data
results_table_spec_ukbms <- results_table_spec_ukbms[,-c(3:5)]
## leave rows with year and interaction
results_table_spec_ukbms$Specialism <- paste(row.names(results_table_spec_ukbms))
rownames(results_table_spec_ukbms) <- 1:nrow(results_table_spec_ukbms)
results_table_spec_ukbms <- results_table_spec_ukbms[-c(1:4,6:14),]
results_table_spec_ukbms$Specialism <- revalue(results_table_spec_ukbms$Specialism, c("mid.year"="Specialist"))
results_table_spec_ukbms$Specialism <- revalue(results_table_spec_ukbms$Specialism, c("mid.year:specialismwider.countryside"="Generalist"))

## create new dataframe
specialist_slope <- results_table_spec_ukbms[1,1]
generalist_slope <- sum(results_table_spec_ukbms$Estimate)
specialist_SE <- results_table_spec_ukbms[1,2]
interaction_SE <- results_table_spec_ukbms[2,2]
generalist_SE <- sqrt((specialist_SE)^2) + ((interaction_SE)^2)
model_summary <- data.frame(Specialism=c("Specialist", "Generalist"), slope=c(specialist_slope, generalist_slope), SE=c(specialist_SE, generalist_SE))

model_summary$Specialism <- factor(model_summary$Specialism, levels=c("Generalist", "Specialist"))
results_table$Specialism <- factor(results_table$Specialism, levels=c("Generalist", "Specialist"))

png("../Graphs/Specialism/Specialism_change_ukbms_85_00.png", height = 150, width = 180, units = "mm", res = 300)
spec <- ggplot(mapping=aes(x=Specialism, y=slope, ymin=slope-SE, ymax=slope+SE)) +
  geom_point(data=results_table, size=2, color="lightgrey", position=myjit) +
  geom_errorbar(data=results_table, color="lightgrey", position=myjit, width=0.1) +
  geom_point(data=model_summary, size=4) +
  geom_errorbar(data=model_summary, width=0.1) +
  labs(x="Specialism", y="Change in population synchrony 1985-2000") +
  #scale_y_continuous(breaks = seq(-1.8,0.8,0.4)) +
  theme_bw() +
  theme(text = element_text(size = 12), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black"))
spec
## generalists have a more negative change in syncrhony than specialists 
dev.off()


############### LATE MODEL ############### 
## run model without specialism for each species and extract the slope (mid.year) coefficient for each
results_table<-NULL
for (i in unique(pair_attr_late$spp)){
  print(i)
  ## if statement: if if the mixed effects model doesn't throw an error (levels of pairID), then run an lmer, else run a linear model 
  loadError=F
  a=try({spec_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain + 
                              winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])})
  loadError <- (is(a, 'try-error')|is(a,'error'))
  if(loadError==F){
    spec_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain + 
                         winter_temp + autumn_temp + spring_temp + summer_temp + (1|pair.id), data=pair_attr_late[pair_attr_late$spp==i,])
  }else{
    spec_ukbms <- lm(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + winter_rain + autumn_rain + spring_rain + summer_rain + 
                       winter_temp + autumn_temp + spring_temp + summer_temp, data=pair_attr_late[pair_attr_late$spp==i,])
  }
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(spec_ukbms)$coefficients[,1:2],i)
  results_table_temp$parameter <- paste(row.names(results_table_temp))
  rownames(results_table_temp) <- 1:nrow(results_table_temp)
  results_table_temp <- results_table_temp[grep("mid.year", results_table_temp$parameter),]
  results_table<-rbind(results_table,results_table_temp)
}

## merge in this dataframe with specialism info (spec/gen)
specialism <- pair_attr_late[,c(14,26)]
specialism <- unique(specialism)
results_table <- merge(results_table, specialism, by.x="i", by.y="spp")
## remove parameter column
results_table <- results_table[,-c(4)]
## re-name some columns
names(results_table) <- c("species", "slope", "SE", "Specialism")
results_table$Specialism <- revalue(results_table$Specialism, c("specialist"="Specialist"))
results_table$Specialism <- revalue(results_table$Specialism, c("wider.countryside"="Generalist"))

########## calculate slope and SE from main model (with specialism)
### plot points and SE ontop of raw data
results_table_spec_ukbms2 <- results_table_spec_ukbms2[,-c(3:5)]
## leave rows with year and interaction
results_table_spec_ukbms2$Specialism <- paste(row.names(results_table_spec_ukbms2))
rownames(results_table_spec_ukbms2) <- 1:nrow(results_table_spec_ukbms2)
results_table_spec_ukbms2 <- results_table_spec_ukbms2[-c(1:4,6:14),]
results_table_spec_ukbms2$Specialism <- revalue(results_table_spec_ukbms2$Specialism, c("mid.year"="Specialist"))
results_table_spec_ukbms2$Specialism <- revalue(results_table_spec_ukbms2$Specialism, c("mid.year:specialismwider.countryside"="Generalist"))

## create new dataframe
specialist_slope <- results_table_spec_ukbms2[1,1]
generalist_slope <- sum(results_table_spec_ukbms2$Estimate)
specialist_SE <- results_table_spec_ukbms2[1,2]
interaction_SE <- results_table_spec_ukbms2[2,2]
generalist_SE <- sqrt((specialist_SE)^2) + ((interaction_SE)^2)
model_summary <- data.frame(Specialism=c("Specialist", "Generalist"), slope=c(specialist_slope, generalist_slope), SE=c(specialist_SE, generalist_SE))

model_summary$Specialism <- factor(model_summary$Specialism, levels=c("Generalist", "Specialist"))
results_table$Specialism <- factor(results_table$Specialism, levels=c("Generalist", "Specialist"))

png("../Graphs/Specialism/Specialism_change_ukbms_00_12.png", height = 150, width = 180, units = "mm", res = 300)
ggplot(mapping=aes(x=Specialism, y=slope, ymin=slope-SE, ymax=slope+SE)) +
  geom_point(data=results_table, size=2, color="lightgrey", position=myjit) +
  geom_errorbar(data=results_table, color="lightgrey", position=myjit, width=0.1) +
  geom_point(data=model_summary, size=4) +
  geom_errorbar(data=model_summary, width=0.1) +
  labs(x="Specialism", y="Change in population synchrony 1985-2000") +
  #scale_y_continuous(breaks = seq(-0.3,0.3,0.05)) +
  theme_bw() +
  theme(text = element_text(size = 12), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black"))
## generalists have a more negative change in syncrhony than specialists 
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

##### cbc dataset #####

############################### interaction in model #################################
#### subset only years 1985 and 2012
pair_attr_cbc_1985 <- pair_attr_CBC[pair_attr_CBC$mid.year==1984.5,]
pair_attr_cbc_1996 <- pair_attr_CBC[pair_attr_CBC$mid.year==1995.5,]
pair_attr_cbc <- rbind(pair_attr_cbc_1985, pair_attr_cbc_1996)
pair_attr_cbc <- droplevels(pair_attr_cbc)

pair_attr_cbc$mid.year <- as.factor(pair_attr_cbc$mid.year)
pair_attr_cbc$pair.id <- as.character(pair_attr_cbc$pair.id)
pair_attr_cbc$spp <- as.factor(pair_attr_cbc$spp)
pair_attr_cbc$specialism <- as.factor(pair_attr_cbc$specialism)
pair_attr_cbc$hab_sim <- as.factor(pair_attr_cbc$hab_sim)

###### run model with strategy and year interaction
spec_model_cbc2 <- lmer(lag0 ~ mean_northing + distance + hab_sim + summer_temp + mid.year*specialism + (1|family) + (1|pair.id) + (1|spp), data = pair_attr_cbc)
## no errors
summary(spec_model_cbc2) ## interaction is significant (p=0.00057)
anova(spec_model_cbc2)
## save model output
results_table_strategy_cbc <- data.frame(summary(spec_model_cbc2)$coefficients[,1:5]) ## 26 species
write.csv(results_table_strategy_cbc, file = "../Results/Model_outputs/CBC/change_spec_cbc.csv", row.names=TRUE)

## predict new data
newdata_cbc <- expand.grid(mean_northing=mean(pair_attr_cbc$mean_northing), distance=mean(pair_attr_cbc$distance), hab_sim=sample(pair_attr_cbc$hab_sim,1),
                       mid.year=unique(pair_attr_cbc$mid.year), specialism=unique(pair_attr_cbc$specialism), family=sample(pair_attr_cbc$family,10),
                       pair.id=sample(pair_attr_cbc$pair.id,10), spp=sample(pair_attr_cbc$spp,10))

newdata_cbc$lag0 <- predict(spec_model_cbc2, newdata=newdata_cbc, re.form=NA)

mm <- model.matrix(terms(spec_model_cbc2), newdata_cbc)
pvar1 <- diag(mm %*% tcrossprod(vcov(spec_model_cbc2),mm))
tvar1 <- pvar1+VarCorr(spec_model_cbc2)$spp[1]+VarCorr(spec_model_cbc2)$pair.id[1]
cmult <- 2

newdata_cbc <- data.frame(
  newdata_cbc
  , plo = newdata_cbc$lag0-1.96*sqrt(pvar1)
  , phi = newdata_cbc$lag0+1.96*sqrt(pvar1)
  , tlo = newdata_cbc$lag0-1.96*sqrt(tvar1)
  , thi = newdata_cbc$lag0+1.96*sqrt(tvar1)
)

## run model without specialism or species random effect to obtain residuals
spec_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + (1|family) + (1|pair.id), data=pair_attr_cbc)
pair_attr_cbc$residuals <- resid(spec_cbc)

## create dataframe which calculates mean, SD and SE of residuals for each species
summary_cbc <- pair_attr_cbc %>% group_by(spp, specialism, mid.year) %>%
  summarise_at(vars(residuals), funs(mean,std.error,sd))

## change values
newdata_cbc$mid.year <- revalue(newdata_cbc$mid.year, c("1984.5"="1985"))
newdata_cbc$mid.year <- revalue(newdata_cbc$mid.year, c("1995.5"="1996"))
colnames(newdata_cbc)[5] <- "Specialism"
newdata_cbc$Specialism <- revalue(newdata_cbc$Specialism, c("generalist"="Generalist"))
newdata_cbc$Specialism <- revalue(newdata_cbc$Specialism, c("specialist"="Specialist"))
## same for summary dataframe
summary_cbc$mid.year <- revalue(summary_cbc$mid.year, c("1984.5"="1985"))
summary_cbc$mid.year <- revalue(summary_cbc$mid.year, c("1995.5"="1996"))
colnames(summary_cbc)[2] <- "Specialism"
summary_cbc$Specialism <- revalue(summary_cbc$Specialism, c("generalist"="Generalist"))
summary_cbc$Specialism <- revalue(summary_cbc$Specialism, c("specialist"="Specialist"))

## plot graph with raw data residuals (+SE error bars) and fitted lines
png("../Graphs/Specialism/Specialism_change_predicted_cbc.png", height = 150, width = 180, units = "mm", res = 300)
ggplot(summary_cbc, aes(x = mid.year, y = mean, group=Specialism)) +
  geom_point(aes(shape=Specialism), colour="grey66", size = 3, position = myjit) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey66", width=0.2, position = myjit) +
  geom_line(data=newdata_cbc, aes(x=mid.year, y=lag0, linetype=Specialism), colour="black", lwd=1) +
  labs(x="Mid year of moving window", y="Population synchrony") +
  theme_bw() +
  theme(legend.key.width = unit(0.8,"cm"), legend.key = element_rect(size = 2), text = element_text(size = 12), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.margin=margin(c(-10,20,-10,0)),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  scale_linetype_manual(name=" ",
                        labels=c("Generalist", "Specialist"), values=c(1,2)) +
  scale_shape_manual(name="Biotype specialism", 
                     labels=c("Generalist", "Specialist"), values=c(16,4)) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(linetype = guide_legend(override.aes = list(size = 0.5)))
dev.off()

## plot shows that specialists are increasing in synchrony, whereas generalists are relatively constant over time

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

##### bbs dataset #####
############################### interaction in model #################################
#### subset only years 1985 and 2012
pair_attr_bbs_1999 <- pair_attr_BBS[pair_attr_BBS$mid.year==1998.5,]
pair_attr_bbs_2012 <- pair_attr_BBS[pair_attr_BBS$mid.year==2011.5,]
pair_attr_bbs <- rbind(pair_attr_bbs_1999, pair_attr_bbs_2012)

pair_attr_bbs$mid.year <- as.factor(pair_attr_bbs$mid.year)
pair_attr_bbs$pair.id <- as.character(pair_attr_bbs$pair.id)
pair_attr_bbs$spp <- as.factor(pair_attr_bbs$spp)
pair_attr_bbs$specialism <- as.factor(pair_attr_bbs$specialism)

###### run model with strategy and year interaction
spec_model_bbs2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + winter_rain + autumn_rain + spring_rain + mid.year*specialism + (1|pair.id) + (1|spp), data = pair_attr_bbs)
summary(spec_model_bbs2)
anova(spec_model_bbs2) ## just significant! (p=0.045)
results_table_strategy_bbs <- data.frame(summary(spec_model_bbs2)$coefficients[,1:5]) ## 24 species
write.csv(results_table_strategy_bbs, file = "../Results/Model_outputs/BBS/change_spec_bbs.csv", row.names=TRUE)
