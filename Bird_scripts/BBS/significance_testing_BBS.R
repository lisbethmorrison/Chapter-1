###################################################################
## Title: Testing model significance & creating barcharts BBS data
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: May 2018
###################################################################

rm(list=ls()) # clear R

library(lme4)
library(lmerTest)
library(dplyr)
options(scipen=999)


## load data
pair_attr <- read.csv("../Data/Bird_sync_data/pair_attr_BBS_final.csv", header=TRUE)

### just want to compare 1985 with 1996 
## to see if there has been an increase/decrease/no change in connectivity

## create new pair_attr file with just those years
pair_attr_1999 <- pair_attr[pair_attr$mid.year==1998.5,]
pair_attr_2012 <- pair_attr[pair_attr$mid.year==2011.5,]
pair_attr_comp <- rbind(pair_attr_1999, pair_attr_2012) ## rbind them together
summary(pair_attr_comp)
unique(pair_attr_comp$mid.year) ## 1998.5 and 2011.5

## centre and standardise continuous variables by mean and SD
pair_attr_comp$distance <- (pair_attr_comp$distance - mean(na.omit(pair_attr_comp$distance)))/sd(na.omit(pair_attr_comp$distance))
pair_attr_comp$mean_northing <- (pair_attr_comp$mean_northing - mean(na.omit(pair_attr_comp$mean_northing)))/sd(na.omit(pair_attr_comp$mean_northing))

## write file to save time
write.csv(pair_attr_comp, file = "../Data/Bird_sync_data/pair_attr_comp_BBS.csv", row.names = FALSE) 

## read in file
pair_attr_comp <- read.csv("../Data/Bird_sync_data/pair_attr_comp_BBS.csv", header=TRUE)

#### change variables to factors/characters
str(pair_attr_comp)
pair_attr_comp$pair.id <- paste("ID", pair_attr_comp$site1, pair_attr_comp$site2, sep = "_")
pair_attr_comp$end.year <- as.factor(pair_attr_comp$end.year)
pair_attr_comp$mid.year <- as.factor(pair_attr_comp$mid.year)
pair_attr_comp$start.year <- as.factor(pair_attr_comp$start.year)
pair_attr_comp$pair.id <- as.character(pair_attr_comp$pair.id)
pair_attr_comp$spp <- as.factor(pair_attr_comp$spp)

##############################
######## run models ##########
############################## 

## likelihood ratio test to test if year has a significant OVERALL effect ##

overall_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id) + (1|spp), REML=FALSE, data = pair_attr_comp)
overall_model_null <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + (1|pair.id) + (1|spp), REML=FALSE, data = pair_attr_comp)
## model comparison
model_comp <- anova(overall_model_null, overall_model, test="Chisq")
## year has a significant effect on year (p<0.01)
results_overall_all_spp <- data.frame(summary(overall_model)$coefficients[,1:5])
names(results_overall_all_spp) <- c("Estimate", "SD", "df", "t","p_value")
results_overall_all_spp$parameter <- paste(row.names(results_overall_all_spp))
rownames(results_overall_all_spp) <- 1:nrow(results_overall_all_spp)

## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_table5 <- NULL
results_table1 <- results_overall_all_spp[grep("mean_northing", results_overall_all_spp$parameter),]
results_table2 <- results_overall_all_spp[grep("distance", results_overall_all_spp$parameter),]
results_table3 <- results_overall_all_spp[grep("hab_sim", results_overall_all_spp$parameter),]
results_table4 <- results_overall_all_spp[grep("(Intercept)", results_overall_all_spp$parameter),]
results_table5 <- rbind(results_table5, results_table4, results_table1, results_table2, results_table3)
results_overall_all_spp <- results_overall_all_spp[!results_overall_all_spp$parameter%in%results_table5$parameter,]

## create year to the comparison (years 1985 and 1996)
results_overall_all_spp$year <- "1999/2012"
results_overall_all_spp$trend <- "significant decrease"
#### connectivity is significantly declining between 1999 and 2012 ==> but only a small decline

############################################
#### model comparison for each species ####
############################################

spp.list<- unique(pair_attr_comp$spp)

results_table<-NULL
for (i in spp.list){
  print(i)
  
  ## create unique pair_attr for each species
  p_attr <- pair_attr_comp[pair_attr_comp$spp==i,]
  
 # if loop (if nrow when start.year=1980 is <20 then skip that species)
  if(nrow(p_attr[p_attr$mid.year=="1998.5",]) < 3){
    print(paste("skip species",i))
    next

  } else {

  model_comp <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + mid.year + (1|pair.id), REML=FALSE, data = pair_attr_comp[pair_attr_comp$spp==i,])
  summary(model_comp)
  anova(model_comp)
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(model_comp)$coefficients[,1:5],i)
  results_table <-rbind(results_table,results_table_temp)
  
  } ## species 472 does work ==> 23 species in total now 
}

## Add REML=FALSE because you should use ML (maximum likelihood) when comparing models that differ in their fixed effects

names(results_table) <- c("Estimate", "SD", "df", "t","p_value", "species")
results_table$parameter <- paste(row.names(results_table))
rownames(results_table) <- 1:nrow(results_table)

## take out 3 columns for each species: mean northing, distance and renk_hab_sim
results_table5 <- NULL
results_table1 <- results_table[grep("mean_northing", results_table$parameter),]
results_table2 <- results_table[grep("distance", results_table$parameter),]
results_table3 <- results_table[grep("hab_sim", results_table$parameter),]
results_table4 <- results_table[grep("(Intercept)", results_table$parameter),]
results_table5 <- rbind(results_table5, results_table4, results_table1, results_table2, results_table3)
results_table <- results_table[!results_table$parameter%in%results_table5$parameter,]

## change parameter names to year of interest 
results_table$parameter <- 1999/2012

## classify each species as either no change (insignificant p value), decreasing (negative estimte), or increasing (positive estimate)
results_table <- results_table %>% group_by(species) %>% mutate(change = ifelse(p_value>0.05, "No change", ifelse(p_value<0.05 & Estimate<0, "Decrease", "Increase")))

nrow(results_table[results_table$change=="No change",]) ## 6 species unchanged
nrow(results_table[results_table$change=="Decrease",]) ## 13 species decreasing
nrow(results_table[results_table$change=="Increase",]) ## 4 species increasing

## make final table
results_final <- results_table[,c(1,5,6,8)]

## save results
write.csv(results_final, file="../Results/Bird_results/model_comp_results_BBS.csv", row.names=FALSE)

results_final_overall <- data.frame(change=unique(results_final$change), no_species=c(4,13,6), total_species=23)
results_final_overall$percentage <- (results_final_overall$no_species / results_final_overall$total_species)*100

## save file
write.csv(results_final_overall, file="../Results/Bird_results/model_comp_percentages_BBS.csv", row.names=FALSE)

levels(results_final_overall$change)
results_final_overall$change <- factor(results_final_overall$change, levels=c("Increase", "No change", "Decrease"))
levels(results_final_overall$change)

######## input the CBC model comparison percentages to create barchart

model_comp_CBC <- read.csv("../Results/Bird_results/model_comp_percentages_CBC.csv", header=TRUE)
model_comp_BBS <- read.csv("../Results/Bird_results/model_comp_percentages_BBS.csv", header=TRUE)

## add in scheme column
model_comp_CBC$comparison <- "1985-1996"
model_comp_BBS$comparison <- "1999-2012"

model_comp <- rbind(model_comp_BBS, model_comp_CBC)
model_comp$comparison <- factor(model_comp$comparison, levels=c("1985-1996", "1999-2012"))
model_comp$change <- factor(model_comp$change, levels=c("Increase", "No change", "Decrease"))

## pretty graph
png("../Graphs/Model_comps/Model_comp_results_CBC_BBS.png", height = 100, width = 120, units = "mm", res = 300)
ggplot(data=model_comp, aes(x=comparison, y=percentage, fill=change)) +
  geom_bar(stat="identity", width=0.4) +
  labs(y="Percentage of species", x="", fill="") +
  scale_fill_manual(values = c("#339900", "#999999", "#990000")) +
  theme_bw() +
  theme(axis.text.x=element_text(colour = "black")) +
  theme(axis.text.y=element_text(colour = "black")) +
  scale_y_continuous(breaks = seq(0,100,10), expand = c(0, 0)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
dev.off()



