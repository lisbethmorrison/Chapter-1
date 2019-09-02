############################# Packages Required ###############################

rm(list=ls()) # clear R

library(tidyverse)
library(ggpubr)
library(gridExtra)


############################### Read in Data ##################################

# List of 14 species of interest
species_list <- read.csv("../IBM_UKBMS/species_list.csv")
# Traits database of 14 species
trait_data <- read.csv("../IBM_UKBMS/UKBMS_16spp_traits.csv")
# Site index scores for all years of UKBMS up to 2016 for UK species
s_index <-read.csv("../IBM_UKBMS/ukbms_sindex.csv")
# CEH landcover map data
gis_data <- read.table("../IBM_UKBMS/BBS.CBC.UKBMS.complete.landcover.data.all.scales.soil.dem.configuration.txt",header=T)

# change all column heading to lower case
names(s_index)<-tolower(names(s_index))

############################# Data Manipulation ###############################


#### ~ Site Selection ####

# Remove all site index data before 2007 - only want last 10 years of data 
# 2007-2016 inc. Only want where brood = 0 and where s_index is 0 or greater
s_index <- s_index %>%
  filter(year >= 2007) %>%
  filter(brood == 0) %>%
  filter(sindex >= 0)

## merge in species_list so sindex only has the 14 species we're interested in
s_index <- merge(s_index, species_list, by="species", all=FALSE)
length(unique(s_index$species)) ## 16 species

s_index$site <- as.factor(s_index$site)
s_index$species <- as.factor(s_index$species)
s_index$year <- as.factor(s_index$year)
s_index$sindex <- as.numeric(s_index$sindex)

# Create list of potential sites by finding out how many sites have 10 years of 
# counts for all species
possible_sites <- s_index %>% count(site,species)

# Remove any sites that have less than 5 years for all species
possible_sites2 <- possible_sites %>% 
  group_by(species) %>%
  filter(n>=5)
# ## now some sites don't have all 14 species
# ## remove sites which don't have all 14 species
# ## number of rows per site (needs to == 14)
# possible_sites3 <- possible_sites2 %>%
#   group_by(site) %>% 
#   summarise(n = n())
# possible_sites3 <- subset(possible_sites3, n==14) ## 102 sites

# Create list of sites
site_list <- unique(possible_sites2$site) # 1895 sites

# Filter gis data to only include UKBMS sites with 500m buffer
# Select only woodland (combine broadleaf and coniferous), arable, chalk and 
# urban
## only 505 (out of 1895) have gis data
gis_data <- gis_data %>%
  filter(Surv == "UKBMS" & buffer == "2000") %>%
  filter(siteno.gref %in% site_list) %>%
  mutate(W = BW + CW) %>%
  mutate(non_chalkG = LM + NG + UDG + G)
  # select(siteno.gref, A, W, LC, G) %>%
  #rename(arable = A, woodland = W, chalk_grassland = LC)
## select columns
gis_data <- gis_data[,c(1,17:29,70:73,86:87)]


## calc percentage
gis_data$total_count <- rowSums(gis_data[2:20])
gis_data <- gis_data %>%
  mutate_at(vars(A:non_chalkG) , funs(P = ./gis_data$total_count * 100))
## remove columns with raw data
gis_data <- gis_data[,-c(2:21)]
## re-name columns
names(gis_data)[2] <- "arable"
names(gis_data)[19] <- "woodland"
names(gis_data)[15] <- "chalk_grassland"
names(gis_data)[20] <- "non_chalk_grassland"
## remove other columns
gis_data <- gis_data[,c(1:2,15,19:20)]

## from each habitat type, choose the top 5 sites (end up with 20 sites chosen)
top_sites <- gis_data %>% gather(key=habitat_type,value=value,-siteno.gref) %>% #convert to long format
  group_by(habitat_type) %>% #group by habitat type
  top_n(5,value) %>% #select the top 5 sites for each habitat type
  arrange(habitat_type,-value) #sort by habitat type and descending value
top_sites <- droplevels(top_sites)
## two arable sites have exact same percentage
## remove one of them (keep site = 173)
top_sites <- top_sites[-c(6),]
top_sites <- droplevels(top_sites)
# Create list of sites
selected_sites <- unique(top_sites$siteno.gref) ## 20 unique sites

#### ~ Species Selection ####

# Filter s_index data to include include the 20 selected sites
species_data <- merge(s_index, top_sites, by.x="site", by.y="siteno.gref", all=FALSE) # 20 sites and 14 species
## remove columns not needed (brood and value)
species_data <- species_data[,-c(3,8)]

# Calculate proportions of each species at each site type
## mean is calculated across 5 sites for each species and over time (at least 5 years of data)
grouped_species <- species_data %>%
  group_by(site, common_name) %>%
  summarise (mean = mean(sindex)) %>%
  mutate(proportion = mean / sum(mean))
## complete dataset so each site has each species (fill these with zeros)
grouped_species <- as.data.frame(complete(grouped_species,common_name,site,fill=list(mean=0,proportion=0)))
site_match <- species_data[,c(1,6)]
site_match <- unique(site_match)
## merge this back in 
grouped_species <- merge(grouped_species, site_match, by="site")
## now take the mean proportion for each habitat type
grouped_species2 <- grouped_species %>%
  group_by(habitat_type, common_name) %>%
  summarise(proportion=mean(proportion))

# To check if the above has worked - all proportions per site type should sum
# to 1
test1 <- filter(grouped_species2, habitat_type == "arable")
sum(test1$proportion)
test2 <- filter(grouped_species2, habitat_type == "woodland")
sum(test2$proportion)
test3 <- filter(grouped_species2, habitat_type == "chalk_grassland")
sum(test3$proportion)
test4 <- filter(grouped_species2, habitat_type == "non_chalk_grassland")
sum(test4$proportion)

# Create a list of final species by selecting columns from grouped_species
final_species <- grouped_species2 %>%
  select(habitat_type, common_name, proportion)

# Add traits data to get whether wider countryside or specialist#
traits <- trait_data[,c(2,32)] ## common name and specialism
final_species <- merge(final_species, traits, by = "common_name")

#####################################
# # Get list of habitat types
# habitat_type <- unique(final_species$habitat_type)
# 
# # Create list of specialism types
# specialism <- c("wider.countryside", "specialist")
# 
# # Create empty dataframes
# specialist_species <- NULL
# wider_species <- NULL
# 
# # Loop finds the specialist species with the highest proportion at each site 
# # type
# for(i in site_type){
#   specialists <- NULL
#   another_loop <- final_species %>%
#     filter(specialism == "specialist") %>%
#     filter(site_type == i)
#   another_loop <- arrange(another_loop, desc(proportion)) 
#   another_loop <- slice(another_loop, 1)
#   specialists <- rbind(another_loop, specialists)
#   specialist_species <- rbind(specialists, specialist_species)
# }
# 
# # As above but for wider countryside and top 2 species
# for(i in site_type){
#   wider_country <- NULL
#   another_loop <- final_species %>%
#     filter(specialism == "wider.countryside") %>%
#     filter(site_type == i)
#   another_loop <- arrange(another_loop, desc(proportion)) 
#   another_loop <- slice(another_loop, 1:2)
#   wider_country <- rbind(another_loop, wider_country)
#   wider_species <- rbind(wider_country, wider_species)
# }
# 
# # Combine specialiss and wider species
# chosen_species <- rbind(specialist_species, wider_species)
# 
# # Create list of the species
# chosen_species <- unique(chosen_species$common_name)
# 
# # Using the new list of species filter the dataset grouped_species
# plot_data <- as.data.frame(grouped_species %>%
#                              filter(common_name %in% chosen_species))
# 
# # Species that score zero not included in data, need to add them in
# a <- as.data.frame(t(c("arable", "Marsh fritillary", 0,0)))
# names(a) <- c("site_type", "common_name", "sum", "proportion")
# 
# b <- as.data.frame(t(c("grassland", "Chalk-hill blue", 0,0)))
# names(b) <- c("site_type", "common_name", "sum", "proportion")
# 
# c <- as.data.frame(t(c("woodland", "Marsh fritillary", 0,0)))
# names(c) <- c("site_type", "common_name", "sum", "proportion")
# 
# plot_data <- rbind(plot_data, a, b, c)
#####################################

# Set proportion to numeric rather than character
final_species$proportion <- as.numeric(final_species$proportion)

# Round proportion to 3dp
final_species$proportion <- round(final_species$proportion, digits = 3)


################################# Plotting ####################################


#### ~ Arable ####


# Filter plot data to only include arable sites
arable <- final_species %>%
  filter(habitat_type == "arable")

# Plot it
arable_plot <- ggplot(arable,
                      aes(common_name, proportion, fill = specialism)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Species") +
  labs(title = "arable") +
  scale_fill_discrete(name = "Specialism",
                      labels = c("Specialist", "Wider Countryside"))


arable_plot

## choose 8 species for Luke's poster
arable <- arable[c(2:4,7,9:11,15),] 

## save plot
png("../IBM_UKBMS/arable_proportions.png", height = 250, width = 300, units = "mm", res = 300)
arable_plot2 <- ggplot(arable,
                      aes(common_name, proportion, fill = specialism)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Species") +
  ylab("Proportion") +
  ylim(0,0.45) +
  labs(title = "Arable") +
  scale_fill_discrete(name = "Specialism",
                      labels = c("Specialist", "Wider Countryside")) +
  theme(axis.text.x=element_text(size=20, angle=-45, hjust=0), axis.text.y=element_text(size=20), axis.title=element_text(size=20),
        legend.title=element_text(size=20), legend.key.size = unit(0.9, "cm"),
        legend.text=element_text(size=20))
arable_plot2
dev.off()

#### ~ Woodland ####


woodland <- final_species %>%
  filter(habitat_type == "woodland")

woodland_plot <-
  ggplot(woodland, aes(common_name, proportion, fill = specialism)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Species") +
  labs(title = "woodland") +
  scale_fill_discrete(name = "Specialism",
                      labels = c("Specialist", "Wider Countryside"))

woodland_plot


#### ~ Chalk_Grassland ####


chalk_grassland <- final_species %>%
  filter(habitat_type == "chalk_grassland")

chalk_grassland_plot <-
  ggplot(chalk_grassland,
         aes(common_name, proportion, fill = specialism)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Species") +
  labs(title = "chalk grassland") +
  scale_fill_discrete(name = "Specialism",
                      labels = c("Specialist", "Wider Countryside"))

chalk_grassland_plot

## choose 8 species for Luke's poster
chalk_grassland <- chalk_grassland[c(2:4,7,9:11,15),] 

## save plot
png("../IBM_UKBMS/chalk_grass_proportions.png", height = 250, width = 300, units = "mm", res = 300)
chalk_grassland_plot2 <- ggplot(chalk_grassland,
                       aes(common_name, proportion, fill = specialism)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Species") +
  ylab("Proportion") +
  ylim(0,0.45) +
  labs(title = "Chalk grassland") +
  scale_fill_discrete(name = "Specialism",
                      labels = c("Specialist", "Wider Countryside")) +
  theme(axis.text.x=element_text(size=20, angle=-45, hjust=0), axis.text.y=element_text(size=20), axis.title=element_text(size=20),
        legend.title=element_text(size=20), legend.key.size = unit(0.9, "cm"),
        legend.text=element_text(size=20))
chalk_grassland_plot2
dev.off()


#### ~ Grassland ####


grassland <- final_species %>%
  filter(habitat_type == "non_chalk_grassland")

grassland_plot <- ggplot(grassland, aes(common_name, proportion, fill = specialism)) +
  geom_bar(stat="identity") +
  theme_classic() +
  xlab("Species") +
  labs(title = "non_chalk_grassland") +
  scale_fill_discrete(name = "Specialism", labels = c("Specialist", "Wider Countryside")) 

grassland_plot


#### ~ All Plots ####


ggarrange(
  arable_plot,
  woodland_plot,
  chalk_grassland_plot,
  grassland_plot,
  ncol = 2,
  nrow = 2,
  common.legend = TRUE,
  legend = "bottom"
)


