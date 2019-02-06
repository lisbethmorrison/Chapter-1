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
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE) # butterfly pair attribute data
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC.csv", header=TRUE) # CBC pair attribute data
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE) # BBS pair attribute data

##############################################
#### specialism model for ALL butterflies ####
##############################################

pair_attr$end.year <- as.factor(pair_attr$end.year)
pair_attr$mid.year <- as.factor(pair_attr$mid.year)
pair_attr$start.year <- as.factor(pair_attr$start.year)
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$specialism <- as.factor(pair_attr$specialism)

start_time <- Sys.time()
spec_model_ukbms1 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + specialism + (1|pair.id) + (1|spp), data = pair_attr)
end_time <- Sys.time()
end_time - start_time ## 7.36 mins

## try and reshuffle specialism variable
pair_attr$specialism_shuffle <- sample(pair_attr$specialism)

### run 9 models 
 

summary(spec_model_ukbms1)
## specialists are the intercept
anova(spec_model_ukbms1)
### WC higher average synchrony compared to HS

###############################################
########### Same for woodland birds ###########
###############################################

## CBC ##
str(pair_attr_CBC)

pair_attr_CBC$end.year <- as.factor(pair_attr_CBC$end.year)
pair_attr_CBC$mid.year <- as.factor(pair_attr_CBC$mid.year)
pair_attr_CBC$start.year <- as.factor(pair_attr_CBC$start.year)
pair_attr_CBC$pair.id <- as.character(pair_attr_CBC$pair.id)
pair_attr_CBC$spp <- as.factor(pair_attr_CBC$spp)

## run model with specialism as fixed categorical variable
strategy_model_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + specialism + (1|pair.id) + (1|spp), data = pair_attr_CBC)
summary(strategy_model_cbc)
anova(strategy_model_cbc)
## non-significant (p=-0.63)

## check fit
sresid <- resid(strategy_model_cbc, type="pearson")
hist(sresid)

## BBS ##
str(pair_attr_BBS)

pair_attr_BBS$distance <- pair_attr_BBS$distance - mean(na.omit(pair_attr_BBS$distance))
pair_attr_BBS$distance <- pair_attr_BBS$distance/(max(abs(na.omit(pair_attr_BBS$distance))))

pair_attr_BBS$mean_northing <- pair_attr_BBS$mean_northing - mean(na.omit(pair_attr_BBS$mean_northing))
pair_attr_BBS$mean_northing <- pair_attr_BBS$mean_northing/(max(abs(na.omit(pair_attr_BBS$mean_northing))))

pair_attr_BBS$end.year <- as.factor(pair_attr_BBS$end.year)
pair_attr_BBS$mid.year <- as.factor(pair_attr_BBS$mid.year)
pair_attr_BBS$start.year <- as.factor(pair_attr_BBS$start.year)
pair_attr_BBS$pair.id <- as.character(pair_attr_BBS$pair.id)
pair_attr_BBS$spp <- as.factor(pair_attr_BBS$spp)

## run model with strategy as fixed categorical variable
strategy_model_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + specialism + (1|pair.id) + (1|spp), data = pair_attr_BBS)
summary(strategy_model_bbs)
anova(strategy_model_bbs)
## non-significant (p=0.82)

## check fit
sresid <- resid(strategy_model_bbs, type="pearson")
hist(sresid)


#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################

### Now test for change in synchony against specialism 

############################### interaction in model #################################
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE)
pair_attr_CBC <- read.csv("../Data/Bird_sync_data/pair_attr_CBC.csv", header=TRUE)
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS_1.csv", header=TRUE)

#### subset only years 1985, 2000 and 2012
pair_attr_1985 <- pair_attr[pair_attr$mid.year==1984.5,]
pair_attr_2000 <- pair_attr[pair_attr$mid.year==1999.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]
pair_attr_ukbms <- rbind(pair_attr_1985, pair_attr_2000, pair_attr_2012) ## 3 years and 33 species
pair_attr_ukbms <- droplevels(pair_attr_ukbms)

pair_attr_ukbms$mid.year <- as.factor(pair_attr_ukbms$mid.year)
pair_attr_ukbms$pair.id <- as.character(pair_attr_ukbms$pair.id)
pair_attr_ukbms$spp <- as.factor(pair_attr_ukbms$spp)
pair_attr_ukbms$specialism <- as.factor(pair_attr_ukbms$specialism)

###### run model with strategy and year interaction (without migrants)
start_time <- Sys.time()
spec_model_ukbms2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*specialism + (1|pair.id) + (1|spp), data = pair_attr_ukbms)
end_time <- Sys.time()
end_time - start_time ## 22 seconds

summary(spec_model_ukbms2)
## intercept is mid.year 1985 and wider countryside species
anova(spec_model_ukbms2)
## mid.year*strategy interaction is significant overall
## save model output
results_table_spec_ukbms <- data.frame(summary(spec_model_ukbms2)$coefficients[,1:5])
write.csv(results_table_strategy_ukbms, file = "../Results/Model_outputs/change_spec_ukbms.csv", row.names=TRUE)

## predict new data
newdata <- expand.grid(mean_northing=mean(pair_attr_ukbms$mean_northing), distance=mean(pair_attr_ukbms$distance), 
              renk_hab_sim=mean(pair_attr_ukbms$renk_hab_sim), mid.year=unique(pair_attr_ukbms$mid.year), 
              specialism=unique(pair_attr_ukbms$specialism), pair.id=sample(pair_attr_ukbms$pair.id,100), 
              spp=unique(pair_attr_ukbms$spp))
                         
newdata$lag0 <- predict(spec_model_ukbms2, newdata=newdata, re.form=NA)

mm2 <- model.matrix(terms(spec_model_ukbms2), newdata)
pvar2 <- diag(mm2 %*% tcrossprod(vcov(spec_model_ukbms2),mm2))
tvar2 <- pvar2+VarCorr(spec_model_ukbms2)$spp[1]+VarCorr(spec_model_ukbms2)$pair.id[1]
cmult <- 2

newdata <- data.frame(
  newdata
  , plo = newdata$lag0-1.96*sqrt(pvar2)
  , phi = newdata$lag0+1.96*sqrt(pvar2)
  , tlo = newdata$lag0-1.96*sqrt(tvar2)
  , thi = newdata$lag0+1.96*sqrt(tvar2)
)

## run model without specialism or species random effect to obtain residuals
spec_ukbms <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data=pair_attr_ukbms)
pair_attr_ukbms$residuals <- resid(spec_ukbms)
## check mean synchrony of each group
group_by(pair_attr_ukbms, specialism) %>% summarize(m = mean(lag0)) ## WC mean = 0.265, HS mean = 0.195
## put mean of each group into pair_attr dataframe
pair_attr_ukbms <- ddply(pair_attr_ukbms, "specialism", transform, spec_mean = mean(lag0))
## add mean to each residual
pair_attr_ukbms$residuals2 <- pair_attr_ukbms$residuals + pair_attr_ukbms$spec_mean


## create dataframe which calculates mean, SD and SE of residuals for each species
summary_ukbms <- pair_attr_ukbms %>% group_by(spp, specialism, mid.year) %>% 
  summarise_at(vars(residuals2), funs(mean,std.error,sd))

## change year values
newdata$mid.year <- revalue(newdata$mid.year, c("1984.5"="1985"))
newdata$mid.year <- revalue(newdata$mid.year, c("1999.5"="2000"))
newdata$mid.year <- revalue(newdata$mid.year, c("2011.5"="2012"))
## change strategy heading
colnames(newdata)[5] <- "Specialism"
newdata$Specialism <- revalue(newdata$Specialism, c("specialist"="Specialist"))
newdata$Specialism <- revalue(newdata$Specialism, c("wider.countryside"="Generalist"))
newdata$Specialism <- factor(newdata$Specialism, levels=c("Generalist", "Specialist"))
## same as above for summary dataframe
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("1984.5"="1985"))
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("1999.5"="2000"))
summary_ukbms$mid.year <- revalue(summary_ukbms$mid.year, c("2011.5"="2012"))
## change strategy heading
colnames(summary_ukbms)[2] <- "Specialism"
summary_ukbms$Specialism <- revalue(summary_ukbms$Specialism, c("specialist"="Specialist"))
summary_ukbms$Specialism <- revalue(summary_ukbms$Specialism, c("wider.countryside"="Generalist"))
summary_ukbms$Specialism <- factor(summary_ukbms$Specialism, levels=c("Generalist", "Specialist"))

## plot graph with raw data residuals (+SE error bars) and fitted lines
png("../Graphs/Specialism/Specialism_change_predicted_ukbms.png", height = 150, width = 180, units = "mm", res = 300)
butterfly_spec <- ggplot(summary_ukbms, aes(x = mid.year, y = mean, group=Specialism)) +
  geom_point(aes(shape=Specialism), colour="grey66", size = 2, position = myjit) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey66", width=0.2, position = myjit) +
  geom_line(data=newdata, aes(x=mid.year, y=lag0, linetype=Specialism), colour="black", lwd=0.5) +
  labs(x="Mid year of moving window", y="Population synchrony") +
  theme_bw() +
  theme(legend.key.width = unit(0.8,"cm"), legend.key = element_rect(size = 2), text = element_text(size = 8), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.margin=margin(c(-15,20,-5,0)),
        axis.text.x=element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  scale_shape_manual(name="Biotype specialism", 
                     labels=c("Generalist", "Specialist"), values=c(16,4)) +
  scale_linetype_manual(name=" ",
                        labels=c("Generalist", "Specialist"), values=c(1,2)) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(linetype = guide_legend(override.aes = list(size = 0.5)))
  
butterfly_spec

dev.off()

install.packages("ggpubr")
library(ggpubr)
png("../Graphs/FINAL/Figure3.png", height = 200, width = 110, units = "mm", res = 300)
ggarrange(butterfly_spec, butterfly_abund, 
          labels = c("(a)", "(b)"), font.label = list(size = 10, color ="black"),
          ncol = 1, nrow = 2)
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

## plot graph with fitted line and confidence interval error bars
pd <- position_dodge(0.1)
png("../Graphs/Specialism/Specialism_change_predicted_ukbms2.png", height = 100, width = 120, units = "mm", res = 300)
ggplot(newdata, aes(x=mid.year, y=lag0, group=Specialism)) +
  geom_point(position=pd) +
  geom_line(aes(linetype=Specialism), size=1, position=pd) +
  labs(x="Year", y="Population synchrony") +
  theme_bw() +
  geom_linerange(aes(ymin=plo, ymax=phi), position=pd) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()


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
spec_model_cbc2 <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year*specialism + (1|pair.id) + (1|spp), data = pair_attr_cbc)
summary(spec_model_cbc2)
## intercept is mid.year 1985 and generalists
anova(spec_model_cbc2)
## mid.year*strategy interaction is NON-significant overall
## save model output
results_table_strategy_cbc <- data.frame(summary(spec_model_cbc2)$coefficients[,1:5])
write.csv(results_table_strategy_cbc, file = "../Results/Model_outputs/change_strategy_cbc.csv", row.names=TRUE)

## predict new data
newdata_cbc <- expand.grid(mean_northing=mean(pair_attr_cbc$mean_northing), distance=mean(pair_attr_cbc$distance), hab_sim=sample(pair_attr_cbc$hab_sim,1),
                       mid.year=unique(pair_attr_cbc$mid.year), specialism=unique(pair_attr_cbc$specialism),
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
spec_cbc <- lmer(lag0 ~ mean_northing + distance + hab_sim + mid.year + (1|pair.id), data=pair_attr_cbc)
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
png("../Graphs/Specialism/Specialism_change_predicted_cbc.png", height = 80, width = 120, units = "mm", res = 300)
pd <- position_dodge(0.2)
ggplot(summary_cbc, aes(x = mid.year, y = mean, group=Specialism)) +
  geom_point(aes(shape=Specialism), colour="grey", size = 2, position=pd) +
  scale_shape_manual(values=c(16,4)) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey", width=0.1, position=pd) +
  geom_line(data=newdata_cbc, aes(x=mid.year, y=lag0, linetype=Specialism), lwd=1) +
  labs(x="Mid year of moving window", y="Population synchrony") +
  theme_bw() +
  theme(text = element_text(size = 10), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

## plot graph with fitted line and confidence interval error bars
png("../Graphs/Specialism/Specialism_change_predicted_cbc2.png", height = 100, width = 120, units = "mm", res = 300)
pd <- position_dodge(0.1)
ggplot(newdata_cbc, aes(x=mid.year, y=lag0, group=Specialism)) +
  geom_point(position=pd) +
  geom_line(aes(linetype=Specialism), size=1, position=pd) +
  labs(x="Year", y="Population synchrony") +
  theme_bw() +
  geom_linerange(aes(ymin=plo, ymax=phi), position=pd) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
## plot shows that generalists are declining in connectivity over time, whereas specialists are increasing

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
start_time <- Sys.time()
spec_model_bbs2 <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*specialism + (1|pair.id) + (1|spp), data = pair_attr_bbs)
end_time <- Sys.time()
end_time - start_time ## 39 seconds

## save main (true) model results
main_result_table <- data.frame(anova(spec_model_bbs2)[,5:6]) ## save anova table from main model
main_result_table$i <- 0 ## make i column with zeros 
main_result_table$parameter <- paste(row.names(main_result_table)) ## move row.names to parameter column
rownames(main_result_table) <- 1:nrow(main_result_table) ## change row names to numbers
## remove rows with mean northing, distance, hab sim and mid year F values (only interested in abundance F values)
main_result_table <- main_result_table[ !(main_result_table$parameter %in% c("mean_northing", "distance", "renk_hab_sim", "mid.year", "specialism")), ]

## run model 999 times
perm_bbs_spec <- NULL
start_time <- Sys.time()
for (i in 1:999){
  print(i)
  pair_attr_bbs$spec_shuffle <- sample(pair_attr_bbs$specialism) ## randomly shuffle abundance change varaible
  spec_model_bbs_shuffle <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year*spec_shuffle + (1|pair.id) + (1|spp), data = pair_attr_bbs)
  ## run model with shuffled variable
  ## save results
  results_table_temp <- data.frame(anova(spec_model_bbs_shuffle)[,5:6],i)
  perm_bbs_spec<-rbind(perm_bbs_spec,results_table_temp)
}
end_time <- Sys.time()
end_time - start_time ## 8.81 hours to do 999 runs

perm_bbs_spec$parameter <- paste(row.names(perm_bbs_spec)) ## move row.names to parameter column
rownames(perm_bbs_spec) <- 1:nrow(perm_bbs_spec) ## change row names to numbers
## remove rows with mean northing, distance, hab sim and mid year F values (only interested in abundance F values)
perm_bbs_spec <- perm_bbs_spec[-grep("mean_northing", perm_bbs_spec$parameter),]
perm_bbs_spec <- perm_bbs_spec[-grep("distance", perm_bbs_spec$parameter),]
perm_bbs_spec <- perm_bbs_spec[-grep("hab_sim", perm_bbs_spec$parameter),]
perm_bbs_spec <- perm_bbs_spec[-grep("mid.year", perm_bbs_spec$parameter),]

final_results_table <- rbind(main_result_table, perm_bbs_spec) ## bind the two data frames together

F_value <- with(final_results_table, final_results_table$F.value[final_results_table$i==0]) ## true F value from main model
hist(final_results_table$F.value) + abline(v=F_value, col="red") ## plot distribution of F values with vertical line (true F value)

## save file
write.csv(final_results_table, file = "../Results/Model_outputs/perm_change_spec_bbs.csv", row.names=TRUE)

## Calculate p value
number_of_permutations <- 1000
final_results_table2 <- final_results_table[!final_results_table$i==0,] ## remove true value to calc. p value
diff.observed <- main_result_table$F.value ## true F value

# P-value is the fraction of how many times the permuted difference is equal or more extreme than the observed difference
pvalue = sum(abs(final_results_table2$F.value) >= abs(diff.observed)) / number_of_permutations
pvalue ## 0 (exactly) significant







summary(spec_model_bbs2)
## intercept is mid.year 1999 and generalists
anova(spec_model_bbs2)
## mid.year*strategy interaction is significant overall (just)
## save model output
results_table_strategy_bbs <- data.frame(summary(spec_model_bbs2)$coefficients[,1:5])
write.csv(results_table_strategy_bbs, file = "../Results/Model_outputs/change_strategy_bbs.csv", row.names=TRUE)

## predict new data
newdata_bbs <- expand.grid(mean_northing=mean(pair_attr_bbs$mean_northing), distance=mean(pair_attr_bbs$distance), renk_hab_sim=mean(pair_attr_bbs$renk_hab_sim),
                       mid.year=unique(pair_attr_bbs$mid.year), specialism=unique(pair_attr_bbs$specialism),
                       pair.id=sample(pair_attr_bbs$pair.id,10), spp=sample(pair_attr_bbs$spp,10))
newdata_bbs$lag0 <- predict(spec_model_bbs2, newdata=newdata_bbs, re.form=NA)

mm3 <- model.matrix(terms(spec_model_bbs2), newdata_bbs)
pvar3 <- diag(mm3 %*% tcrossprod(vcov(spec_model_bbs2),mm3))
tvar3 <- pvar3+VarCorr(spec_model_bbs2)$spp[1]+VarCorr(spec_model_bbs2)$pair.id[1]
cmult <- 2

newdata_bbs <- data.frame(
  newdata_bbs
  , plo = newdata_bbs$lag0-1.96*sqrt(pvar3)
  , phi = newdata_bbs$lag0+1.96*sqrt(pvar3)
  , tlo = newdata_bbs$lag0-1.96*sqrt(tvar3)
  , thi = newdata_bbs$lag0+1.96*sqrt(tvar3)
)

## run model without specialism or species random effect to obtain residuals
spec_bbs <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), data=pair_attr_bbs)
pair_attr_bbs$residuals <- resid(spec_bbs)

## create dataframe which calculates mean, SD and SE of residuals for each species
summary_bbs <- pair_attr_bbs %>% group_by(spp, specialism, mid.year) %>% 
  summarise_at(vars(residuals), funs(mean,std.error,sd))

## change values
newdata_bbs$mid.year <- revalue(newdata_bbs$mid.year, c("1998.5"="1999"))
newdata_bbs$mid.year <- revalue(newdata_bbs$mid.year, c("2011.5"="2012"))
colnames(newdata_bbs)[5] <- "Specialism"
newdata_bbs$Specialism <- revalue(newdata_bbs$Specialism, c("generalist"="Generalist"))
newdata_bbs$Specialism <- revalue(newdata_bbs$Specialism, c("specialist"="Specialist"))
## same as above for summary dataframe
summary_bbs$mid.year <- revalue(summary_bbs$mid.year, c("1998.5"="1999"))
summary_bbs$mid.year <- revalue(summary_bbs$mid.year, c("2011.5"="2012"))
colnames(summary_bbs)[2] <- "Specialism"
summary_bbs$Specialism <- revalue(summary_bbs$Specialism, c("generalist"="Generalist"))
summary_bbs$Specialism <- revalue(summary_bbs$Specialism, c("specialist"="Specialist"))

## plot graph with raw data residuals (+SE error bars) and fitted lines
png("../Graphs/Specialism/Specialism_change_predicted_bbs.png", height = 80, width = 120, units = "mm", res = 300)
pd <- position_dodge(0.2)
ggplot(summary_bbs, aes(x = mid.year, y = mean, group=Specialism)) +
  geom_point(aes(shape=Specialism), colour="grey", size = 2, position=pd) +
  scale_shape_manual(values=c(16,4)) +
  geom_errorbar(aes(ymin = mean-std.error, ymax = mean+std.error), colour="grey", width=0.1, position=pd) +
  geom_line(data=newdata_bbs, aes(x=mid.year, y=lag0, linetype=Specialism), lwd=1) +
  labs(x="Mid year of moving window", y="Population synchrony") +
  theme_bw() +
  theme(text = element_text(size = 10), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

## plot graph with fitted line and confidence interval error bars
png("../Graphs/Specialism/Specialism_change_predicted_bbs2.png", height = 100, width = 120, units = "mm", res = 300)
ggplot(newdata_bbs, aes(x=mid.year, y=lag0, group=Strategy)) +
  geom_point(position=pd) +
  geom_line(aes(linetype=Strategy), size=1, position=pd) +
  labs(x="Year", y="Population synchrony") +
  theme_bw() +
  geom_linerange(aes(ymin=plo, ymax=phi), position=pd) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
## generalists have declined in connectivity more steeply than specialists 

