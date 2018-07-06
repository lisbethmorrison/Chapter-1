###########################################################
## Title: Testing larvae length and average synchrony for butterflies
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: May 2018
##########################################################

rm(list=ls()) # clear R

library(lme4)
library(lmerTest)

### add data
species_traits <- read.csv("../Data/UKBMS_data/species.traits.full.table.csv", header=TRUE) ## add trait data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE) 

## create larvae_size dataframe
larvae_size <- species_traits[-c(1:2,4,6:38,40)] ## remove columns not needed - left with species code, family and larvae length
larvae_size <- na.omit(larvae_size) ## remove NAs
summary(larvae_size) ## check data

## merge larvae_size with pair_attr
pair_attr <- merge(pair_attr, larvae_size, by.x="spp", by.y="species", all=TRUE)
pair_attr <- na.omit(pair_attr) ## remove NAs
length(unique(pair_attr$spp)) ## 34 species (spp 119 doesn't have larvae length data)

pair_attr$mid.year <- as.factor(pair_attr$mid.year)
pair_attr$pair.id <- as.character(pair_attr$pair.id)
pair_attr$spp <- as.factor(pair_attr$spp)
pair_attr$family <- as.character(pair_attr$family)
pair_attr$larvae_length_mm <- as.numeric(pair_attr$larvae_length_mm)

## run model with larvae length and family interaction
larvae_model <- lmer(lag0 ~ mean_northing + distance + renk_hab_sim + family*larvae_length_mm + (1|spp) + (1|pair.id), data=pair_attr)
anova(larvae_model) ## non-significant 
summary(larvae_model)
