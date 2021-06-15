####################################################
##                                                ##
##             RainDrop biodiversity              ##
##                                                ##
##           Total biodiversity models            ##
##                                                ##
##                Apr 15th 2021                   ##
##                                                ##
####################################################

## Building Bayesian heirarchical models of total biodiversity indices with respect to
## the rainfall treatments, including block and annual effects. Simple effects first.

rm(list = ls())
options(width = 100)

library(tidyverse)
library(brms)
library(patchwork)
library(viridis)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load and wrangle data ####

# load data
load("../../RainDropRobotics/Data/raindrop_biodiversity_2016_2020.RData",
     verbose = TRUE)

# Converting to total biomass - Here we are considering TRUE 0 observations - need to chat to Andy about this
totbiomass <- biomass %>%
  group_by(year, harvest, block, treatment) %>%
  summarise(tot_biomass = sum(biomass_g),
            log_tot_biomass = log(tot_biomass + 1)) %>%
  ungroup() %>%
  mutate(year_f = as.factor(year),
         year_s = as.numeric(scale(year)))

## group biomass - Focussing on only graminoids, forbs and legumes
biomass_gr <- biomass %>%
  filter(group %in% c("Forbs", "Graminoids", "Legumes")) %>%
  mutate(biomass_log = log(biomass_g + 1),
         year_f = as.factor(year),
         year_s = as.numeric(scale(year)))


# Percent cover - total biodiversity indices
percent_cover <- percent_cover %>%
  drop_na(percent_cover) %>%
  filter(species_level == "Yes")

diversity <- percent_cover %>%
  # first convert percentages to a proportion
  mutate(proportion = percent_cover/100) %>%
  group_by(year, month, block, treatment) %>%
  # diversity indices for each group
  summarise(richness = n(),
            simpsons = sum(proportion^2),
            shannon = -sum(proportion * log(proportion))) %>%
  ungroup() %>%
  mutate(year_f = as.factor(year),
         year_s = as.numeric(scale(year)),
         shannon = as.numeric(scale(shannon)))

#%>%
  # pivot_longer(cols = c(richness, simpsons, shannon),
  #              names_to = "index") %>%
  # mutate(index_lab = case_when(
  #   index == "richness" ~ "Species richness",
  #   index == "simpsons" ~ "Simpson's index",
  #   index == "shannon" ~ "Shannon-Wiener index"
  # ))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Total biomass models ####

## Simpler model - year as random effect
set.seed(666)
totbiomass_base <- brm(
  log_tot_biomass ~ 1 + harvest + (1|block/treatment) + (1|year_f),
  data = totbiomass, family = gaussian(),
  prior = c(
    prior(normal(3.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(3), class = sd, group = "block:treatment"),
    prior(exponential(3), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
  )

## Treatment effect
set.seed(666)
totbiomass_treatment <- brm(
  log_tot_biomass ~ 1 + treatment + harvest + (1|block/treatment) + (1|year_f),
  data = totbiomass, family = gaussian(),
  prior = c(
    prior(normal(3.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(3), class = sd, group = "block:treatment"),
    prior(exponential(3), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year linear effect
set.seed(666)
totbiomass_year_linear <- brm(
  log_tot_biomass ~ 1 + treatment + harvest + year_s + (1|block/treatment) + (1|year_f),
  data = totbiomass, family = gaussian(),
  prior = c(
    prior(normal(3.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(3), class = sd, group = "block:treatment"),
    prior(exponential(3), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year interaction with harvest
set.seed(666)
totbiomass_year_harvest <- brm(
  log_tot_biomass ~ 1 + treatment + year_s*harvest + (1|block/treatment) + (1|year_f),
  data = totbiomass, family = gaussian(),
  prior = c(
    prior(normal(3.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(3), class = sd, group = "block:treatment"),
    prior(exponential(3), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## treatment interaction with year
set.seed(666)
totbiomass_year_treat <- brm(
  log_tot_biomass ~ 1 + treatment*year_s + harvest + (1|block/treatment) + (1|year_f),
  data = totbiomass, family = gaussian(),
  prior = c(
    prior(normal(3.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(3), class = sd, group = "block:treatment"),
    prior(exponential(3), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year as an autocorrelated term
set.seed(666)
totbiomass_year_auto <- brm(
  log_tot_biomass ~ 1 + treatment + ar(gr = year_s, p = 1) + harvest + (1|block/treatment),
  data = totbiomass, family = gaussian(),
  prior = c(
    prior(normal(3.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(3), class = sd, group = "block:treatment")),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Model comparisons
totbiomass_base <- add_criterion(totbiomass_base, criterion = c("loo","waic"))
totbiomass_treatment <- add_criterion(totbiomass_treatment, criterion = c("loo","waic"))
totbiomass_year_linear <- add_criterion(totbiomass_year_linear, criterion = c("loo","waic"))
totbiomass_year_harvest <- add_criterion(totbiomass_year_harvest, criterion = c("loo","waic"))
totbiomass_year_treat <- add_criterion(totbiomass_year_treat, criterion = c("loo","waic"))
totbiomass_year_auto <- add_criterion(totbiomass_year_auto, criterion = c("loo","waic"))

mod_comp_totbiomass <- as.data.frame(loo_compare(totbiomass_base, totbiomass_treatment,
                                      totbiomass_year_linear, totbiomass_year_treat,
                                      totbiomass_year_harvest, totbiomass_year_auto,
                                      criterion = "loo"))
mod_comp_totbiomass

save(totbiomass_treatment, mod_comp_totbiomass, file = "data/totbiomass_models.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Group-level biomass ####

# Here more conservative prior for the beta terms (sd = 0.5) to reflect more elaborate predictors

## Simpler model - year as random effect
set.seed(666)
biomass_base <- brm(
  biomass_log ~ 1 + harvest + (1|block/treatment) + (1|year_f),
  data = biomass_gr, family = gaussian(),
  prior = c(
    prior(normal(1, 1), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(3), class = sd, group = "block:treatment"),
    prior(exponential(3), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## treatment effect
set.seed(666)
biomass_treatment <- brm(
  biomass_log ~ 1 + treatment + harvest + (1|block/treatment) + (1|year_f),
  data = biomass_gr, family = gaussian(),
  prior = c(
    prior(normal(1, 1), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(3), class = sd, group = "block:treatment"),
    prior(exponential(3), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## group effect
set.seed(666)
biomass_group <- brm(
  biomass_log ~ 1 + group + harvest + (1|block/treatment) + (1|year_f),
  data = biomass_gr, family = gaussian(),
  prior = c(
    prior(normal(1, 1), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(3), class = sd, group = "block:treatment"),
    prior(exponential(3), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## treatment and group effects
set.seed(666)
biomass_treatment_group <- brm(
  biomass_log ~ 1 + treatment + group + harvest + (1|block/treatment) + (1|year_f),
  data = biomass_gr, family = gaussian(),
  prior = c(
    prior(normal(1, 1), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(3), class = sd, group = "block:treatment"),
    prior(exponential(3), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## treatment and group and interaction
set.seed(666)
biomass_treatment_group_int <- brm(
  biomass_log ~ 1 + treatment*group + harvest + (1|block/treatment) + (1|year_f),
  data = biomass_gr, family = gaussian(),
  prior = c(
    prior(normal(1, 1), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(3), class = sd, group = "block:treatment"),
    prior(exponential(3), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## treatment, group and year(linear) interactions
set.seed(666)
biomass_treatment_group_year <- brm(
  biomass_log ~ 1 + treatment*group*year_s + harvest + (1|block/treatment) + (1|year_f),
  data = biomass_gr, family = gaussian(),
  prior = c(
    prior(normal(1, 1), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(3), class = sd, group = "block:treatment"),
    prior(exponential(3), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Model comparisons
biomass_base <- add_criterion(biomass_base, criterion = c("loo","waic"))
biomass_treatment <- add_criterion(biomass_treatment, criterion = c("loo","waic"))
biomass_group <- add_criterion(biomass_group, criterion = c("loo","waic"))
biomass_treatment_group <- add_criterion(biomass_treatment_group, criterion = c("loo","waic"))
biomass_treatment_group_int <- add_criterion(biomass_treatment_group_int, criterion = c("loo","waic"))
biomass_treatment_group_year <- add_criterion(biomass_treatment_group_year, criterion = c("loo","waic"))

mod_comp_biomass <- as.data.frame(loo_compare(biomass_base,
                          biomass_treatment,
                          biomass_group,
                          biomass_treatment_group,
                          biomass_treatment_group_int,
                          biomass_treatment_group_year,
                          criterion = "loo"))

mod_comp_biomass

save(biomass_treatment_group_int, mod_comp_biomass, file = "data/biomass_models.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Shannon-Wiener index across treatments ####

## Base model - year as a random effect
set.seed(666)
shannon_base <- brm(
  shannon ~ 1 + month + (1|block/treatment) + (1|year_f),
  data = diversity, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment"),
    prior(exponential(5), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Simple treatment effect
set.seed(666)
shannon_treatment <- brm(
  shannon ~ 1 + month + treatment + (1|block/treatment) + (1|year_f),
  data = diversity, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment"),
    prior(exponential(5), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year linear effect
set.seed(666)
shannon_year_linear <- brm(
  shannon ~ 1 + month + treatment + year_s + (1|block/treatment) + (1|year_f),
  data = diversity, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment"),
    prior(exponential(5), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year treatment interaction
set.seed(666)
shannon_treatment_year_int <- brm(
  shannon ~ 1 + month + treatment*year_s + (1|block/treatment) + (1|year_f),
  data = diversity, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment"),
    prior(exponential(5), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year as an autocorrelated term
set.seed(666)
shannon_treatment_year_auto <- brm(
  shannon ~ 1 + month + treatment + ar(gr = year_s, p = 1) + (1|block/treatment) + (1|year_f),
  data = diversity, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment"),
    prior(exponential(5), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Model comparisons
shannon_base <- add_criterion(shannon_base, criterion = c("loo","waic"))
shannon_treatment <- add_criterion(shannon_treatment, criterion = c("loo","waic"))
shannon_year_linear <- add_criterion(shannon_year_linear, criterion = c("loo","waic"))
shannon_treatment_year_int <- add_criterion(shannon_treatment_year_int, criterion = c("loo","waic"))
shannon_treatment_year_auto <- add_criterion(shannon_treatment_year_auto, criterion = c("loo","waic"))

mod_comp_shannon <- as.data.frame(loo_compare(shannon_base, shannon_treatment,
                                              shannon_year_linear,
                                              shannon_treatment_year_int,
                                              shannon_treatment_year_auto,
                                              criterion = "loo"))
mod_comp_shannon

save(shannon_treatment_year_int, mod_comp_shannon, file = "data/shannon_models.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Species richness across treatments ####

# looking at a better prior for the intercept
hist(log(diversity$richness))
hist(rnorm(1000, 3,0.25))

## Base model - year as a random effect
set.seed(666)
richness_base <- brm(
  richness ~ 1 + month + (1|block/treatment) + (1|year_f),
  data = diversity, family = poisson(link = "log"),
  prior = c(
    prior(normal(3, 0.25), class =  Intercept),
    prior(normal(0, 0.25), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(8), class = sd, group = "block:treatment"),
    prior(exponential(8), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.97),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Simple treatment effect
set.seed(666)
richness_treatment <- brm(
  richness ~ 1 + month + treatment + (1|block/treatment) + (1|year_f),
  data = diversity, family = poisson(link = "log"),
  prior = c(
    prior(normal(3, 0.25), class =  Intercept),
    prior(normal(0, 0.25), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(8), class = sd, group = "block:treatment"),
    prior(exponential(8), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.97),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year linear effect
set.seed(666)
richness_year_linear <- brm(
  richness ~ 1 + month + treatment + year_s + (1|block/treatment) + (1|year_f),
  data = diversity, family = poisson(link = "log"),
  prior = c(
    prior(normal(3, 0.25), class =  Intercept),
    prior(normal(0, 0.25), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(8), class = sd, group = "block:treatment"),
    prior(exponential(8), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.97),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year treatment interaction
set.seed(666)
richness_treatment_year_int <- brm(
  richness ~ 1 + month + treatment*year_s + (1|block/treatment) + (1|year_f),
  data = diversity, family = poisson(link = "log"),
  prior = c(
    prior(normal(3, 0.25), class =  Intercept),
    prior(normal(0, 0.25), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(8), class = sd, group = "block:treatment"),
    prior(exponential(8), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year as an autocorrelated term
set.seed(666)
richness_treatment_year_auto <- brm(
  richness ~ 1 + month + treatment + ar(gr = year_s, p = 1) + (1|block/treatment) + (1|year_f),
  data = diversity, family = poisson(link = "log"),
  prior = c(
    prior(normal(3, 0.25), class =  Intercept),
    prior(normal(0, 0.25), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(8), class = sd, group = "block:treatment"),
    prior(exponential(8), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Model comparisons
richness_base <- add_criterion(richness_base, criterion = c("loo","waic"))
richness_treatment <- add_criterion(richness_treatment, criterion = c("loo","waic"))
richness_year_linear <- add_criterion(richness_year_linear, criterion = c("loo","waic"))
richness_treatment_year_int <- add_criterion(richness_treatment_year_int, criterion = c("loo","waic"))
richness_treatment_year_auto <- add_criterion(richness_treatment_year_auto, criterion = c("loo","waic"))

mod_comp_richness <- as.data.frame(loo_compare(richness_base, richness_treatment,
                                              richness_year_linear,
                                              richness_treatment_year_int,
                                              richness_treatment_year_auto,
                                              criterion = "loo"))
mod_comp_richness

save(richness_year_linear, mod_comp_richness, file = "data/richness_models.RData")


