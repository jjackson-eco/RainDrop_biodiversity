####################################################
##                                                ##
##             RainDrop biodiversity              ##
##                                                ##
##        RainDrop simple mixed effects model     ##
##                                                ##
##                 Nov 5th 2021                   ##
##                                                ##
####################################################

## Assistance with Andy's teaching to provide a mixed-effect model example of a biomass effect.

rm(list = ls())
options(width = 100)

library(tidyverse)
library(lme4)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

biomass <- read_csv("../../RainDropRobotics/Data/Raindrop_biomass_2016-20.csv")
glimpse(biomass)


## simplify biomass data for teaching - look at years and total biomass across both seasons
biomass %>%
  group_by(year, block, treatment) %>%
  summarise(total_biomass_g = sum(biomass_g)) %>%
  ggplot(aes(x = treatment, y = total_biomass_g)) +
  geom_boxplot() +
  facet_wrap(~year)

# Looks like 2018 is a good example of the trend
biomass2018 <- biomass %>%
  group_by(year, block, treatment) %>%
  summarise(total_biomass_g = sum(biomass_g)) %>%
  ungroup() %>%
  filter(year == 2018)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. simple model ####

biomass_model <- lmer(total_biomass_g ~ treatment + (1|block),
                      data = biomass2018)

summary(biomass_model)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. The data ####

write.csv(biomass2018, file = "../../RainDropRobotics/Data/biomass_total_2018.csv",row.names = F)

