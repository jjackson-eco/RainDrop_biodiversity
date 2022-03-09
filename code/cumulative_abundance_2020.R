####################################################
##                                                ##
##             RainDrop biodiversity              ##
##                                                ##
##   Cumulative species abundances for 2020 data  ##
##                                                ##
##                Aug 12th 2021                   ##
##                                                ##
####################################################

## Looking at the

rm(list = ls())
options(width = 100)

library(tidyverse)
library(patchwork)
library(viridis)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

load("../../RainDropRobotics/Data/raindrop_biodiversity_2016_2020.RData",
     verbose = TRUE)

# colours
raindrop_colours <-
  tibble(treatment = c("Ambient", "Control", "Drought", "Irrigated"),
         num = rep(100,4),
         colour = c("#61D94E", "#BFBFBF", "#EB8344", "#6ECCFF"))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Wrangle data ####

cumulative_abundance2020 <- percent_cover %>%
  # only june 2020 data and only species-level
  filter(year == 2020 &
         species_level == "Yes" &
         month == "June") %>%
  group_by(block, treatment) %>%
  # arrange by the most to least abundant species
  arrange(-percent_cover, .by_group = TRUE) %>%
  # get the species rank, relative and cumulative abundance
  mutate(species_order = 1:n(),
         relative_abundance = percent_cover/sum(percent_cover),
         cumulative_abundance = cumsum(relative_abundance))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Plot ####

p1 <- cumulative_abundance2020 %>%
  ggplot(aes(x = species_order, y = cumulative_abundance, colour = treatment)) +
  geom_hline(yintercept = 0.8) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = raindrop_colours$colour, guide = "none") +
  facet_grid(rows = vars(block), cols = vars(treatment)) +
  labs(x = "Number of species", y = "Cumulative abundance",title = "2020") +
  theme_bw(base_size = 12)

ggsave(plot = p1, filename = "output/cumulative_abundance2020.pdf",
       width = 20, height = 27, units = "cm")




