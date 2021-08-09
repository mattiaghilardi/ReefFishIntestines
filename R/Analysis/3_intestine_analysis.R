#### ------------ This script includes the analysis used to test how reef fish intestinal  ------------ ####
#### ------------         morphology relates to phylogeny, body size, elongation,          ------------ ####
#### ------------         trophic level, the presence of a stomach, and durophagy          ------------ ####


#### ------------ Load the isotope dataset ------------ ####

stomach_durophagy <- read.csv("Data/stomach_durophagy_fish_Moorea.csv")

#### ------------ Analysis ------------ ####

int_moor_traits <- int_moor_traits %>%
  left_join(stomach_durophagy[,3:5]) %>%
  mutate(stomach = factor(stomach, levels = c(0, 1), labels = c('Absent', 'Present')),
         durophagy = factor(durophagy, levels = c(0, 1), labels = c('Not durophagous', 'Durophagous')))

# Select species with at least 3 sampled individuals
int_moor_traits_f <- int_moor_traits %>%
  group_by(species) %>%
  mutate(n = n()) %>%
  filter(n > 2) %>%
  ungroup() %>%
  mutate(id = as.factor(as.character(id)), # 1208 observations
         family = as.factor(as.character(family)), # 31 families
         genus = as.factor(as.character(genus)), # 72 genera
         species = as.factor(as.character(species))) # 142 species

# Summarise data
int_moor_traits_f_summary <- int_moor_traits_f %>%
  group_by(family, genus, species) %>%
  mutate(min_tl = min(tl),
         max_tl = max(tl)) %>%
  summarise_at(vars(weight:troph, elongation, MaxLengthTL, n, min_tl, max_tl), mean, na.rm = TRUE) %>%
  ungroup() %>%
  select(1:5, 15, 16, 6:14)

# Table with number of individual per habitat (Appendix S2: Table 1)
ind_per_habitat <- int_moor_traits_f %>% 
  group_by(family, species, habitat) %>% 
  count() %>%
  pivot_wider(names_from = habitat, values_from = n) %>% 
  mutate(across(where(is.numeric), ~replace_na(.x, 0)),
         Total = sum(across(where(is.numeric))),
         species = recode(species, "Acanthurus nigroris" = "Acanthurus nigros")) %>% # Rename Acanthurus nigroris as A. nigros
  rename(Family = family, Species = species)
# Save table
write.table(ind_per_habitat, file = "Output/Tables/Individuals per habitat.txt", sep = ";", quote = FALSE, row.names = FALSE)

# Download the phylogeny from The Fish Tree of Life 
# and create a covariance matrix for brms models

# Check species names
check_name_fishtree(unique(as.character(int_moor_traits_f$species)), sampled = TRUE)
# All species names are correct
# These species do not have genetic data:
# [1] "Epinephelus hexagonatus" "Hemiramphus depauperatus" "Exallias brevis"  

# Use function "fishtree_complete_phylogeny" (100 trees)
tree_int_complete <- fishtree_complete_phylogeny(species = unique(int_moor_traits_f$species))

# Download also phylogentetic tree of sampled species only, will be used for figure 1
tree_int <- fishtree_phylogeny(species = unique(int_moor_traits_f$species))

# Save trees
ape::write.tree(tree_int, file = "Output/Trees/tree_intestine")
ape::write.tree(tree_int_complete, file = "Output/Trees/tree_intestine_complete")

# Phylogenetic covariance matrix 
corphy_int <- phylo_cov(tree_int_complete)
corphy_int_m <- corphy_int$mean

# Add column "phylo" to the dataset 
int_moor_traits_f$phylo <- gsub(" ", "_", int_moor_traits_f$species)

# Prepare the data
int_moor_traits_f <-
  int_moor_traits_f %>%
  
  # natural log transformation
  mutate(il_log = log(intestine_length),
         id_log = log(intestine_diameter),
         is_log = log(intestine_surface),
         sl_log = log(sl),
         elon_log = log(elongation)) %>%
  
  # to investigate interspecific relationship between intestinal traits and body size
  # calculate mean standard length per species
  group_by(phylo) %>% 
  mutate(spec_mean_sl_log = mean(sl_log),
         
         # to investigate intraspecific relationship between intestinal traits and body size
         # calculate intraspecific variability of standard length
         within_spec_sl_log = sl_log - spec_mean_sl_log)%>%
  ungroup()

## Run phylogenetic hierarchical linear models in brms

### INTESTINAL LENGTH

## Model 1: include interspecific relationship with body size, trophic level and phylogeny (random effect)
## Include a further random effect on the intercept to account for any specific effect unrelated to phylogeny (1 |species)
## Scale and center predictors to:
## 1 - investigate their relative contribution
## 2 - provide a meaningful interpretation of the intercept (value at mean body size and mean trophic level)
## Use family "student" that is less influenced by outliers

# Set priors 
# Uninformative normal priors with mean of 0 for intercept, slopes
# Set a normal prior with mean of 0 also for sd because all predictors are standardised
# These will be used for all models
prior <- c(prior(normal(0, 10), class = Intercept),
           prior(normal(0, 5), class = b),
           prior(normal(0, 5), class = sd))

m1_leng <- brm(il_log ~ scale(spec_mean_sl_log) + scale(troph) + (1 |gr(phylo, cov = corphy_int_m)) + (1 |species),
               data = int_moor_traits_f, family = student(), data2 = list(corphy_int_m = corphy_int_m),
               prior = prior, warmup = 2000, iter = 8000, chains = 4, cores = 4,
               inits = "random", file = "Output/Models/m1_leng")

# Model checking
pp_check(m1_leng, type = "dens_overlay", nsamples = 100) # ok
plot(m1_leng, ask = FALSE) # ok
loo_m1_leng <- loo(m1_leng) # All pareto_k < 0.5

# Model summary
summary(m1_leng, priors = TRUE) # positive relationship with body size and negative with trophic level
posterior_summary(m1_leng)[1:7,] 

# Explained variance (R2)
bayes_R2(m1_leng) # 88.7% [88.2%, 89.2%]

# Phylogenetic signal
hyp <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
(hyp <- hypothesis(m1_leng, hyp, class = NULL)) # 0.89[0.80, 0.93]
plot(hyp)


## Model 2: add intraspecific relationship with standard length
m2_leng <- update(
  m1_leng, formula = ~ . + scale(within_spec_sl_log) + (0 + scale(within_spec_sl_log) |species),
  newdata = int_moor_traits_f, warmup = 2000, iter = 8000, chains = 4, cores = 4,
  inits = "random", file = "Output/Models/m2_leng")

# Model checking
pp_check(m2_leng, type = "dens_overlay", nsamples = 100) # ok
plot(m2_leng, ask = FALSE) # ok
loo_m2_leng <- loo(m2_leng) # All pareto_k < 0.7

# Model comparison
loo_compare(loo_m1_leng, loo_m2_leng) # model 2 is better

# Model summary
summary(m2_leng, priors = TRUE) # variability in scaling with body size (however many species have only few replicates and the results for them are not reliable)
posterior_summary(m2_leng)[1:9,]

# Get unscaled slopes
round(fixef(m2_leng)[2,]/sd(int_moor_traits_f$spec_mean_sl_log), 3) # 0.860[0.694, 1.025]
round(fixef(m2_leng)[3,]/sd(int_moor_traits_f$troph), 3) # -0.456[-0.625, -0.297]
round(fixef(m2_leng)[4,]/sd(int_moor_traits_f$within_spec_sl_log), 3) # 0.825[0.713, 0.935]

# Explained variance (R2)
bayes_R2(m2_leng) # 92.8% [92.5, 93.2] --> increased

# Phylogenetic signal
hyp <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
(hyp <- hypothesis(m2_leng, hyp, class = NULL)) # 0.92[0.83, 0.95]
plot(hyp)


## Model 3: add elongation and an interaction between stomach and durophagy
m3_leng <- update(
  m2_leng, formula = ~ . + scale(elon_log) + stomach*durophagy,
  newdata = int_moor_traits_f, warmup = 2000, iter = 8000, chains = 4, cores = 4,
  inits = "random", file = "Output/Models/m3_leng")

# Model checking
pp_check(m3_leng, type = "dens_overlay", nsamples = 100) # ok
plot(m3_leng, ask = FALSE) # ok
loo_m3_leng <- loo(m3_leng) # All pareto_k < 0.7

# Model comparison
loo_compare(loo_m1_leng, loo_m2_leng, loo_m3_leng) 
# No noteworthy difference between model 2 and 3 because elongation, stomach and durophagy are highly phylogenetically conserved
# Use full model to show all relationships

# Model summary
summary(m3_leng, priors = TRUE) # negative relationship with body elongation
posterior_summary(m3_leng)[1:13,]

# Plot interaction between stomach and durophagy
plot(conditional_effects(m3_leng, effects = "stomach:durophagy"))

# Get unscaled slopes
round(fixef(m3_leng)[2,]/sd(int_moor_traits_f$spec_mean_sl_log), 3) # 0.972[0.818, 1.128]
round(fixef(m3_leng)[3,]/sd(int_moor_traits_f$troph), 3) # -0.379[-0.531, -0.236]
round(fixef(m3_leng)[4,]/sd(int_moor_traits_f$within_spec_sl_log), 3) # 0.826[0.716, 0.936]
round(fixef(m3_leng)[5,]/sd(int_moor_traits_f$elon_log), 3) # -0.784[-1.047, -0.516]

# Explained variance (R2)
bayes_R2(m3_leng) # 92.8% [92.4, 93.1] --> same as model 2

# Phylogenetic signal
hyp <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
(hyp <- hypothesis(m3_leng, hyp, class = NULL)) # 0.90[0.80, 0.94]
plot(hyp)




### INTESTINAL DIAMETER

## Model 1
m1_diam <- brm(id_log ~ scale(spec_mean_sl_log) + scale(troph) + (1 |gr(phylo, cov = corphy_int_m)) + (1 |species),
               data = int_moor_traits_f, family = student(), data2 = list(corphy_int_m = corphy_int_m),
               prior = prior, warmup = 2000, iter = 8000, chains = 4, cores = 4,
               inits = "random", file = "Output/Models/m1_diam")

# Model checking
pp_check(m1_diam, type = "dens_overlay", nsamples = 100) # ok
plot(m1_diam, ask = FALSE) # ok
loo_m1_diam <- loo(m1_diam) # All pareto_k < 0.5

# Model summary
summary(m1_diam, priors = TRUE) # intestinal diameter also decrease with trophic level
posterior_summary(m1_diam)[1:7,]
# Explained variance (R2)
bayes_R2(m1_diam) # 79.0% [77.8%, 80.0%]

# Phylogenetic signal
hyp <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
(hyp <- hypothesis(m1_diam, hyp, class = NULL)) # 0.6[0.39, 0.77]
plot(hyp)


## Model 2
m2_diam <- update(
  m1_diam, formula = ~ . + scale(within_spec_sl_log) + (0 + scale(within_spec_sl_log) |species),
  newdata = int_moor_traits_f, warmup = 2000, iter = 8000, chains = 4, cores = 4,
  inits = "random", file = "Output/Models/m2_diam")

# Model checking
pp_check(m2_diam, type = "dens_overlay", nsamples = 100) # ok
plot(m2_diam, ask = FALSE) # ok
loo_m2_diam <- loo(m2_diam) # All pareto_k < 0.7

# Model comparison
loo_compare(loo_m1_diam, loo_m2_diam) # model 2 is better

# Model summary
summary(m2_diam, priors = TRUE)
posterior_summary(m2_diam)[1:9,]

# Get unscaled slopes
round(fixef(m2_diam)[2,]/sd(int_moor_traits_f$spec_mean_sl_log), 3) # 0.838[0.720, 0.956]
round(fixef(m2_diam)[3,]/sd(int_moor_traits_f$troph), 3) # -0.166[-0.276, -0.052]
round(fixef(m2_diam)[4,]/sd(int_moor_traits_f$within_spec_sl_log), 3) # 0.717[0.616, 0.819]

# Explained variance (R2)
bayes_R2(m2_diam) # 84.6% [83.6%, 85.4%]

# Phylogenetic signal
hyp <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
(hyp <- hypothesis(m2_diam, hyp, class = NULL)) # 0.66[0.44, 0.83]
plot(hyp)


## Model 3
m3_diam <- update(
  m2_diam, formula = ~ . + scale(elon_log) + stomach*durophagy,
  newdata = int_moor_traits_f, warmup = 2000, iter = 8000, chains = 4, cores = 4,
  inits = "random", file = "Output/Models/m3_diam")

# Model checking
pp_check(m3_diam, type = "dens_overlay", nsamples = 100) # ok
plot(m3_diam, ask = FALSE) # ok
loo_m3_diam <- loo(m3_diam) # All pareto_k < 0.7

# Model comparison
loo_compare(loo_m1_diam, loo_m2_diam, loo_m3_diam) 
# Same as intestinal length, no noteworthy difference between model 2 and 3

# Model summary
summary(m3_diam, priors = TRUE) # negative relationship with body elongation
posterior_summary(m3_diam)[1:13,]

# Plot interaction between stomach and durophagy
plot(conditional_effects(m3_diam, effects = "stomach:durophagy"))

# Get unscaled slopes
round(fixef(m3_diam)[2,]/sd(int_moor_traits_f$spec_mean_sl_log), 3) # 0.931[0.828, 1.034]
round(fixef(m3_diam)[3,]/sd(int_moor_traits_f$troph), 3) # -0.165[-0.250, -0.074]
round(fixef(m3_diam)[4,]/sd(int_moor_traits_f$within_spec_sl_log), 3) # 0.718[0.615, 0.819]
round(fixef(m3_diam)[5,]/sd(int_moor_traits_f$elon_log), 3) # -0.415[-0.558, -0.276]

# Explained variance (R2)
bayes_R2(m3_diam) # 84.6% [83.6%, 85.5%]

# Phylogenetic signal
hyp <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
(hyp <- hypothesis(m3_diam, hyp, class = NULL)) # 0.34[0.12, 0.59]
plot(hyp)




### INTESTINAL SURFACE

## Model 1
m1_surf <- brm(is_log ~ scale(spec_mean_sl_log) + scale(troph) + (1 |gr(phylo, cov = corphy_int_m)) + (1 |species),
               data = int_moor_traits_f, family = student(), data2 = list(corphy_int_m = corphy_int_m),
               prior = prior, warmup = 2000, iter = 8000, chains = 4, cores = 4,
               inits = "random", file = "Output/Models/m1_surf")

# Model checking
pp_check(m1_surf, type = "dens_overlay", nsamples = 100) # ok
plot(m1_surf, ask = FALSE) # ok
loo_m1_surf <- loo(m1_surf) # All pareto_k < 0.5

# Model summary
summary(m1_surf, priors = TRUE)
posterior_summary(m1_surf)[1:7,]

# Explained variance (R2)
bayes_R2(m1_surf) # 85.1% [84.4%, 85.7%]

# Phylogenetic signal
hyp <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
(hyp <- hypothesis(m1_surf, hyp, class = NULL)) # 0.82[0.67, 0.89]
plot(hyp)


## Model 2
m2_surf <- update(
  m1_surf, formula = ~ . + scale(within_spec_sl_log) + (0 + scale(within_spec_sl_log) |species),
  newdata = int_moor_traits_f, warmup = 2000, iter = 8000, chains = 4, cores = 4,
  inits = "random", control = list(max_treedepth = 12), file = "Output/Models/m2_surf")

# Model checking
pp_check(m2_surf, type = "dens_overlay", nsamples = 100) # ok
plot(m2_surf, ask = FALSE) # ok
loo_m2_surf <- loo(m2_surf) # All pareto_k < 0.7

# Model comparison
loo_compare(loo_m1_surf, loo_m2_surf) # model 2 is better

# Model summary
summary(m2_surf, priors = TRUE)
posterior_summary(m2_surf)[1:9,]

# Get unscaled slopes
round(fixef(m2_surf)[2,]/sd(int_moor_traits_f$spec_mean_sl_log), 3) # 1.682[1.441, 1.917]
round(fixef(m2_surf)[3,]/sd(int_moor_traits_f$troph), 3) # -0.588[-0.858, -0.344]
round(fixef(m2_surf)[4,]/sd(int_moor_traits_f$within_spec_sl_log), 3) # 1.540[1.360, 1.713]

# Explained variance (R2)
bayes_R2(m2_surf) # 91.9% [91.5%, 92.3%]

# Phylogenetic signal
hyp <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
(hyp <- hypothesis(m2_surf, hyp, class = NULL)) # 0.87[0.72, 0.93]
plot(hyp)


## Model 3
m3_surf <- update(
  m2_surf, formula = ~ . + scale(elon_log) + stomach*durophagy,
  newdata = int_moor_traits_f, warmup = 2000, iter = 8000, chains = 4, cores = 4,
  inits = "random", file = "Output/Models/m3_surf")

# Model checking
pp_check(m3_surf, type = "dens_overlay", nsamples = 100) # ok
plot(m3_surf, ask = FALSE) # ok
loo_m3_surf <- loo(m3_surf) # All pareto_k < 0.7

# Model comparison
loo_compare(loo_m1_surf, loo_m2_surf, loo_m3_surf) 
# Same as above, no noteworthy difference between model 2 and 3

# Model summary
summary(m3_surf, priors = TRUE) # negative relationship with body elongation
posterior_summary(m3_surf)[1:13,]

# Plot interaction between stomach and durophagy
plot(conditional_effects(m3_surf, effects = "stomach:durophagy"))

# Get unscaled slopes
round(fixef(m3_surf)[2,]/sd(int_moor_traits_f$spec_mean_sl_log), 3) # 1.870[1.650, 2.087]
round(fixef(m3_surf)[3,]/sd(int_moor_traits_f$troph), 3) # -0.548[-0.805, -0.305]
round(fixef(m3_surf)[4,]/sd(int_moor_traits_f$within_spec_sl_log), 3) # 1.539[1.362, 1.713]
round(fixef(m3_surf)[5,]/sd(int_moor_traits_f$elon_log), 3) # -1.195[-1.544, -0.852]

# Explained variance (R2)
bayes_R2(m3_surf) # 91.9% [91.5%, 92.3%]

# Phylogenetic signal
hyp <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
(hyp <- hypothesis(m3_surf, hyp, class = NULL)) # 0.76[0.50, 0.90]
plot(hyp)




### Summary table of the 3 models (Appendix S2: Table 2)
int_model_summary <- left_join(as_tibble(posterior_summary(m3_leng)[1:13, -2], rownames = "parameter"), 
          as_tibble(posterior_summary(m3_diam)[1:13, -2], rownames = "parameter"), 
          by = "parameter", suffix = c(".l", ".d")) %>%
  left_join(as_tibble(posterior_summary(m3_surf)[1:13, -2], rownames = "parameter") %>%
              rename_with( ~ paste0(.x, ".s"), .cols = 2:4),
            by = "parameter") %>%
  mutate(across(-1, round, 2))

# Save table
write.table(int_model_summary, file = "Output/Tables/Summary intestine models.txt", sep = ";", row.names = FALSE)




### Calculate % of variation of intestinal traits in the range of trophic level relative to the lowest trophic level

# Get min and max value of trophic level
min_troph_moor <- tibble(troph = min(int_moor_traits_f$troph),
                         spec_mean_sl_log = mean(int_moor_traits_f$spec_mean_sl_log),
                         within_spec_sl_log = mean(int_moor_traits_f$within_spec_sl_log),
                         elon_log = mean(int_moor_traits_f$elon_log),
                         stomach = "Present",
                         durophagy = "Not durophagous")
                        
max_troph_moor <- tibble(troph = max(int_moor_traits_f$troph),
                         spec_mean_sl_log = mean(int_moor_traits_f$spec_mean_sl_log),
                         within_spec_sl_log = mean(int_moor_traits_f$within_spec_sl_log),
                         elon_log = mean(int_moor_traits_f$elon_log),
                         stomach = "Present",
                         durophagy = "Not durophagous")

# Here we need exp = TRUE in "percent_variation()" because all intestinal traits were log-transformed before analysis
# Intestinal length --> -59.46%
percent_variation(m3_leng, min_troph_moor, max_troph_moor, exp = TRUE)

# Intestinal diameter --> -32.49%
percent_variation(m3_diam, min_troph_moor, max_troph_moor, exp = TRUE)

# Intestinal surface --> -72.87%
percent_variation(m3_surf, min_troph_moor, max_troph_moor, exp = TRUE)




### Scaling parameters for each trait and species

# Extract scaling values for each species

# INTESTINAL LENGTH
slope_leng <- m3_leng %>%
  tidybayes::spread_draws(b_scalewithin_spec_sl_log, r_species[species,term]|term) %>% # get draws of random intercept and slope for species
  mutate(spec_slope = (b_scalewithin_spec_sl_log + scalewithin_spec_sl_log)/sd(int_moor_traits_f$within_spec_sl_log)) %>% # unscale slope values
  group_by(species) %>% 
  tidybayes::median_qi(spec_slope, .width = c(.80, .95)) %>% # get median and 95% credible interval per species
  ungroup() %>%
  mutate(species = gsub("\\.", " ", species)) # adjust species names

# Add number of replicates per species
slope_leng <- left_join(slope_leng, int_moor_traits_f_summary[, c("species", "n")])


# INTESTINAL DIAMETER
slope_diam <- m3_diam %>%
  tidybayes::spread_draws(b_scalewithin_spec_sl_log, r_species[species,term]|term) %>%
  mutate(spec_slope = (b_scalewithin_spec_sl_log + scalewithin_spec_sl_log)/sd(int_moor_traits_f$within_spec_sl_log)) %>%
  group_by(species) %>%
  tidybayes::median_qi(spec_slope, .width = c(.80, .95)) %>%
  ungroup() %>%
  mutate(species = gsub("\\.", " ", species))

# Add number of replicates per species
slope_diam <- left_join(slope_diam, int_moor_traits_f_summary[,c("species", "n")])


# INTESTINAL SURFACE
slope_surf <- m3_surf %>%
  tidybayes::spread_draws(b_scalewithin_spec_sl_log, r_species[species,term]|term) %>%
  mutate(spec_slope = (b_scalewithin_spec_sl_log + scalewithin_spec_sl_log)/sd(int_moor_traits_f$within_spec_sl_log)) %>%
  group_by(species) %>%
  tidybayes::median_qi(spec_slope, .width = c(.80, .95)) %>%
  ungroup() %>%
  mutate(species = gsub("\\.", " ", species))

# Add number of replicates per species
slope_surf <- left_join(slope_surf, int_moor_traits_f_summary[,c("species", "n")])


# Create table with summary information (median and 95% CI)
slope_leng <- slope_leng %>% mutate(summary = paste0(formatC(spec_slope, 2, format = "f"), " (", formatC(.lower, 2, format = "f"), ", ", formatC(.upper, 2, format = "f"), ")"))
slope_diam <- slope_diam %>% mutate(summary = paste0(formatC(spec_slope, 2, format = "f"), " (", formatC(.lower, 2, format = "f"), ", ", formatC(.upper, 2, format = "f"), ")"))
slope_surf <- slope_surf %>% mutate(summary = paste0(formatC(spec_slope, 2, format = "f"), " (", formatC(.lower, 2, format = "f"), ", ", formatC(.upper, 2, format = "f"), ")"))

slope_summary <- tibble(Species = slope_leng[slope_leng$.width==0.95, "species", drop = TRUE], 
                        `Length slope` = slope_leng[slope_leng$.width==0.95, "summary", drop = TRUE], 
                        `Diameter slope` = slope_diam[slope_diam$.width==0.95, "summary", drop = TRUE], 
                        `Surface slope` = slope_surf[slope_surf$.width==0.95, "summary", drop = TRUE], 
                        n = slope_leng[slope_leng$.width==0.95, "n", drop = TRUE])

slope_summary <- left_join(slope_summary, int_moor_traits_f_summary[, c("family", "species", "min_tl", "max_tl")], by = c("Species"="species")) %>%
  mutate(`TL range (mm)` = paste0(min_tl, "-", max_tl), # Size range
         Species = recode(Species, "Acanthurus nigroris" = "Acanthurus nigros")) %>% # Rename Acanthurus nigroris as A. nigros
  rename(Family = family) %>%
  arrange(Family, Species) %>%
  select(6, 1:5, 9)

# Save table (Appendix S2: Table 5)
write.table(slope_summary, file = "Output/Tables/Summary species slopes.txt", sep = ";", row.names = FALSE)

# Select species that meet 3 requirements:
# 1 - have the 95% CI above 0;
# 2 - have at least 10 sampled individuals;
# 3 - have a minimum size range of 25% of the maximum length from FishBase.

# Intestinal length
slope_leng_selection <- left_join(slope_leng, int_moor_traits_f_summary[, c("species", "min_tl", "max_tl", "MaxLengthTL")]) %>%
  mutate(min_tl = min_tl/10,
         max_tl = max_tl/10,
         perc_max_tl = ((max_tl-min_tl)/MaxLengthTL)*100) %>%
  filter(.lower > 0 & n > 9 & perc_max_tl >= 25) %>% 
  group_by(species) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  select(-c(min_tl, max_tl, MaxLengthTL))
# considering 80%CI: 2 species show negative allometry and 1 positive allometry out of 19 species (16%)
# considering 95%CI: 1 species show negative allometry and 1 positive allometry out of 19 species (11%)

# Intestinal diameter
slope_diam_selection <- left_join(slope_diam, int_moor_traits_f_summary[, c("species", "min_tl", "max_tl", "MaxLengthTL")]) %>%
  mutate(min_tl = min_tl/10,
         max_tl = max_tl/10,
         perc_max_tl = ((max_tl-min_tl)/MaxLengthTL)*100) %>%
  filter(.lower > 0 & n > 9 & perc_max_tl >= 25) %>% 
  group_by(species) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  select(-c(min_tl, max_tl, MaxLengthTL)) 
# considering 80%CI: 8 species show negative allometry out of 18 species (44%)
# considering 95%CI: 2 species show negative allometry out of 18 species (11%)

# Intestinal surface
slope_surf_selection <- left_join(slope_surf, int_moor_traits_f_summary[, c("species", "min_tl", "max_tl", "MaxLengthTL")]) %>%
  mutate(min_tl = min_tl/10,
         max_tl = max_tl/10,
         perc_max_tl = ((max_tl-min_tl)/MaxLengthTL)*100) %>%
  filter(.lower > 0 & n > 9 & perc_max_tl >= 25) %>% 
  group_by(species) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  select(-c(min_tl, max_tl, MaxLengthTL)) 
# considering 80%CI: 8 species show negative allometry and 1 positive allometry out of 20 species (45%)
# considering 95%CI: 5 species show negative allometry out of 20 speciesn (25%)

# These tibbles will be used to create figure 4

#### ------------ END ------------ ####
