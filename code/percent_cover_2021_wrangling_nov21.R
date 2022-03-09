####################################################
##                                                ##
##             RainDrop biodiversity              ##
##                                                ##
##      Data cleaning and wrangling for 2021      ##
##                                                ##
##                Nov 12th 2021                   ##
##                                                ##
####################################################

## Loading transposed cover data for 2021 and wrangling it to tidy data

rm(list = ls())
options(width = 100)

library(tidyverse)
library(patchwork)
library(viridis)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

raw_cover_2021 <- read_csv("../../Raindrop Data to Share/Percent_Cover/RainDrop_cover_2021_nontidy.csv")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Wrangling ####

percentcover_2021 <- raw_cover_2021 %>%
  pivot_longer(-c(Date, Block, Treatment), # main wrangling part, this is a big pivot of the table in one line of code
               names_to = "species_name",
               values_to = "percent_cover") %>%
  rename_all(tolower) %>% # all column names to lower case - useful
  filter(is.na(percent_cover) == FALSE) %>%
  # these functions below are more advaced but useful for tidying up data
  mutate(species_level = if_else(species_name %in% c("Bryophytes", "Bareground") == TRUE |
                                   grepl("sp[.]", species_name) == TRUE, 0, 1),
         date = as.Date(date, format = "%d/%m/%Y"),
         treatment = case_when(
           treatment == "Control" ~ "Ambient",
           treatment == "P. Control" ~ "Control",
           TRUE ~ treatment
         ))


##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Save for Phil ####

# RData files are super useful
save(percentcover_2021,file = "../../../Supervision/Phil_Fernandes/percent_cover_2021.RData")







