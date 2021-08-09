#### ------------ This script includes the analysis and figure in appendix 1:  ------------ ####
#### ------------        relationship between delta15N and trophic level       ------------ ####


#### ------------ Load the isotope dataset ------------ ####

d15n_moor <- read.csv("Data/SI_fish_Moorea.csv")

#### ------------ Analysis ------------ ####

# Select species with at least 3 replicates
d15n_moor <- d15n_moor %>%
  group_by(species) %>%
  mutate(n = n()) %>%
  filter(n > 2) %>%
  ungroup() %>%
  mutate(id = as.factor(as.character(id)), # 700 observations
         family = as.factor(as.character(family)), # 22 families
         genus = as.factor(as.character(genus)), # 51 genera
         species = as.factor(as.character(species)), # 83 species
         location = as.factor(as.character(location))) # Moorea

# Summarise data
d15n_summary <- d15n_moor %>%
  group_by(family, genus, species) %>%
  summarise_at(vars(sl:n), mean, na.rm = T) %>%
  ungroup()

# Natural log-transform standard length
d15n_moor$sl_log <- log(d15n_moor$sl)

# Download the phylogeny from The Fish Tree of Life 
# and create a covariance matrix for brms model

# Check species names
check_name_fishtree(unique(as.character(d15n_moor$species)), sampled = TRUE)
# All species names are correct
# These species do not have genetic data:
# [1] "Epinephelus hexagonatus" 

# Use function "fishtree_complete_phylogeny" (100 trees)
tree_d15n <- fishtree_complete_phylogeny(species = unique(d15n_moor$species))

# Save trees
ape::write.tree(tree_d15n, file = "Output/Trees/tree_d15n")

# Phylogenetic covariance matrix 
corphy_d15n <- phylo_cov(tree_d15n)
corphy_d15n_m <- corphy_d15n$mean

# Add column "phylo" to the dataset 
d15n_moor$phylo <- gsub(" ", "_", d15n_moor$species)

# Run a phylogenetic hierarchical linear model for delta15N in brms that include:
# - Body size and trophic level as fixed effects
# - 2 random effects on the intercept that account for species-level effects, both related (phylo) and unrelated (species) to phylogeny
# - Random slope for body size to account for intraspecific scaling
# This method allows to investigate the relationship with trophic level after accounting for body size, phylogeny and other species-level effects
m1_d15n <- brm(delta15N ~ scale(sl_log) + scale(troph) + (1 |gr(phylo, cov = corphy_d15n_m)) + (1 |species) + (0 + scale(sl_log) |species), 
               data = d15n_moor, family = student(), data2 = list(corphy_d15n_m = corphy_d15n_m),
               prior = c(prior(normal(0, 15), class = Intercept),
                         prior(normal(0, 5), class = b),
                         prior(normal(0, 5), class = sd)),
               sample_prior = TRUE, warmup = 1000, iter = 4000, chains = 4,  cores = 4, 
               inits = "random", file = "Output/Models/m1_d15n")

# Model checking
pp_check(m1_d15n, type = "dens_overlay", nsamples = 100) # ok
plot(m1_d15n, N = 8) # ok
loo(m1_d15n) # All pareto_k < 0.7

# Model summary
summary(m1_d15n, priors = TRUE) # strong positive relationship with trophic level

# Appendix S1: table 1
d15n_model_summary <- as_tibble(posterior_summary(m1_d15n)[1:8,-2], rownames = "parameter") %>%
  mutate(across(-1, round, 2))
write.table(d15n_model_summary, file = "Output/Tables/Summary isotope model.txt", sep = ";", row.names = FALSE)

# Get unscaled slopes
round(fixef(m1_d15n)[2,]/sd(d15n_moor$sl_log), 3) # 0.793[0.446, 1.138]
round(fixef(m1_d15n)[3,]/sd(d15n_moor$troph), 3) # 0.970[0.419, 1.513] --> for each unit increase in trophic level, delta15N increases by one unit on average

# Explained variance
bayes_R2(m1_d15n) # 82.2% [81.0%, 83.2%]

# Phylogenetic signal
hyp <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
(hyp <- hypothesis(m1_d15n, hyp, class = NULL)) # 0.74[0.40, 0.94]
plot(hyp)

# Calculate % of variation of delta15N in the range of trophic level relative to lowest trophic level

# Get min and max value of trophic level
min_troph_d15n <- tibble(troph = min(d15n_moor$troph), sl_log = mean(d15n_moor$sl_log))
max_troph_d15n <- tibble(troph = max(d15n_moor$troph), sl_log = mean(d15n_moor$sl_log))

# Here we do not need to exponentiate because delta15N was not log-transformed before analysis
percent_variation(m1_d15n, min_troph_d15n, max_troph_d15n, exp = FALSE)
# 27.58%

#### ------------ END ------------ ####
