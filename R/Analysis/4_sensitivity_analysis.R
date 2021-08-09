#### ------------ This script includes the sensitivity analysis to assess the robustness of our results ------------ ####


#### ------------ Analysis ------------ ####

### 1 - Repeat analysis only for species with at least 5 sampled individuals

# Filter dataset
int_moor_traits_f_5 <- filter(int_moor_traits_f, n > 4)
length(unique(int_moor_traits_f_5$species)) # 122 species - removed 20
length(unique(int_moor_traits_f_5$family)) # 26 families - removed 5

# Download the phylogeny from The Fish Tree of Life 
# and create a covariance matrix for brms models

# Check species names
check_name_fishtree(unique(as.character(int_moor_traits_f_5$species)), sampled = TRUE)
# All species names are correct
# These species do not have genetic data:
# [1] "Epinephelus hexagonatus" "Exallias brevis"  

# Download multiphylo
tree_int_5 <- fishtree_complete_phylogeny(species = unique(int_moor_traits_f_5$species))

# Phylogenetic covariance matrix 
corphy_int_5 <- phylo_cov(tree_int_5)
corphy_int_m_5 <- corphy_int_5$mean

# Run models

# INTESTINAL LENGTH
m_leng_5 <- brm(il_log ~ scale(spec_mean_sl_log) + scale(troph) + scale(within_spec_sl_log) + scale(elon_log) + stomach*durophagy +
                  (1 |gr(phylo, cov = corphy_int_m_5)) + (1 |species) + (0 + scale(within_spec_sl_log) |species),
                data = int_moor_traits_f_5, family = student(), data2 = list(corphy_int_m_5 = corphy_int_m_5),
                prior = prior, warmup = 2000, iter = 8000, chains = 4, cores = 4,
                inits = "random", file = "Output/Models/m_leng_5")

# Checks and summary
pp_check(m_leng_5, type = "dens_overlay", nsamples = 100) # ok
plot(m_leng_5, ask = FALSE) # ok
loo(m_leng_5) # All pareto_k < 0.5
summary(m_leng_5, priors = TRUE)
bayes_R2(m_leng_5) # 92.7% [92.4%, 93.1%]

# Get unscaled slopes
round(fixef(m_leng_5)[2,]/sd(int_moor_traits_f_5$spec_mean_sl_log), 3) # 0.951[0.771, 1.127]
round(fixef(m_leng_5)[3,]/sd(int_moor_traits_f_5$troph), 3) # -0.269[-0.434, -0.109]
round(fixef(m_leng_5)[4,]/sd(int_moor_traits_f_5$within_spec_sl_log), 3) # 0.859[0.747, 0.970]
round(fixef(m_leng_5)[5,]/sd(int_moor_traits_f_5$elon_log), 3) # -0.822[-1.114, -0.523]

# INTESTINAL DIAMETER
m_diam_5 <- brm(id_log ~ scale(spec_mean_sl_log) + scale(troph) + scale(within_spec_sl_log) + scale(elon_log) + stomach*durophagy +
                  (1 |gr(phylo, cov = corphy_int_m_5)) + (1 |species) + (0 + scale(within_spec_sl_log) |species),
                data = int_moor_traits_f_5, family = student(), data2 = list(corphy_int_m_5 = corphy_int_m_5),
                prior = prior, warmup = 2000, iter = 8000, chains = 4, cores = 4,
                inits = "random", file = "Output/Models/m_diam_5")

# Checks and summary
pp_check(m_diam_5, type = "dens_overlay", nsamples = 100) # ok
plot(m_diam_5, ask = FALSE) # ok
loo(m_diam_5) # All pareto_k < 0.7
summary(m_diam_5, priors = TRUE)
bayes_R2(m_diam_5) # 84.6% [83.6%, 85.4%]

# Get unscaled slopes
round(fixef(m_diam_5)[2,]/sd(int_moor_traits_f_5$spec_mean_sl_log), 3) # 0.948[0.827, 1.070]
round(fixef(m_diam_5)[3,]/sd(int_moor_traits_f_5$troph), 3) # -0.162[-0.262, -0.049]
round(fixef(m_diam_5)[4,]/sd(int_moor_traits_f_5$within_spec_sl_log), 3) # 0.722[0.616, 0.830]
round(fixef(m_diam_5)[5,]/sd(int_moor_traits_f_5$elon_log), 3) # -0.423[-0.587, -0.264]

# INTESTINAL SURFACE
m_surf_5 <- brm(is_log ~ scale(spec_mean_sl_log) + scale(troph) + scale(within_spec_sl_log) + scale(elon_log) + stomach*durophagy +
                  (1 |gr(phylo, cov = corphy_int_m_5)) + (1 |species) + (0 + scale(within_spec_sl_log) |species),
                data = int_moor_traits_f_5, family = student(), data2 = list(corphy_int_m_5 = corphy_int_m_5),
                prior = prior, warmup = 2000, iter = 8000, chains = 4, cores = 4,
                inits = "random", file = "Output/Models/m_surf_5")

# Checks and summary
pp_check(m_surf_5, type = "dens_overlay", nsamples = 100) # ok
plot(m_surf_5, ask = FALSE) # ok
loo(m_surf_5) # All pareto_k < 0.7
summary(m_surf_5, priors = TRUE)
bayes_R2(m_surf_5) # 91.7% [91.3%, 92.1%]

# Get unscaled slopes
round(fixef(m_surf_5)[2,]/sd(int_moor_traits_f_5$spec_mean_sl_log), 3) # 1.855[1.609, 2.109]
round(fixef(m_surf_5)[3,]/sd(int_moor_traits_f_5$troph), 3) # -0.396[-0.716, -0.124]
round(fixef(m_surf_5)[4,]/sd(int_moor_traits_f_5$within_spec_sl_log), 3) # 1.578[1.397, 1.756]
round(fixef(m_surf_5)[5,]/sd(int_moor_traits_f_5$elon_log), 3) # -1.245[-1.662, -0.839]


# Summary table of the 3 models (Appendix S2: Table 3)
int_model_5_summary <- left_join(as_tibble(posterior_summary(m_leng_5)[1:13, -2], rownames = "parameter"), 
                               as_tibble(posterior_summary(m_diam_5)[1:13, -2], rownames = "parameter"), 
                               by = "parameter", suffix = c(".l", ".d")) %>%
  left_join(as_tibble(posterior_summary(m_surf_5)[1:13, -2], rownames = "parameter") %>%
              rename_with( ~ paste0(.x, ".s"), .cols = 2:4),
            by = "parameter") %>%
  mutate(across(-1, round, 2))

# Save table
write.table(int_model_5_summary, file = "Output/Tables/Summary intestine models sensitivity 5.txt", sep = ";", row.names = FALSE)


# Calculate % of variation of intestinal traits in the range of trophic level relative to the lowest trophic level
# Get min and max value of trophic level
min_troph_moor_5 <- tibble(troph = min(int_moor_traits_f_5$troph), 
                           spec_mean_sl_log = mean(int_moor_traits_f_5$spec_mean_sl_log),
                           within_spec_sl_log = mean(int_moor_traits_f_5$within_spec_sl_log),
                           elon_log = mean(int_moor_traits_f$elon_log),
                           stomach = "Absent",
                           durophagy = "Not durophagous")
max_troph_moor_5 <- tibble(troph = max(int_moor_traits_f_5$troph), 
                           spec_mean_sl_log = mean(int_moor_traits_f_5$spec_mean_sl_log),
                           within_spec_sl_log = mean(int_moor_traits_f_5$within_spec_sl_log),
                           elon_log = mean(int_moor_traits_f$elon_log),
                           stomach = "Absent",
                           durophagy = "Not durophagous")

# Here we need exp = TRUE in "percent_variation()" because all intestinal traits were log-transformed before analysis
# Intestinal length --> -47.22%
percent_variation(m_leng_5, min_troph_moor_5, max_troph_moor_5, exp = TRUE)

# Intestinal diameter --> -32.01%
percent_variation(m_diam_5, min_troph_moor_5, max_troph_moor_5, exp = TRUE)

# Intestinal surface --> -61.02%
percent_variation(m_surf_5, min_troph_moor_5, max_troph_moor_5, exp = TRUE)

# Results for all 3 models are consistent with the main analysis



### 2 - Repeat analysis only for species with at least 10 sampled individuals

# Filter dataset
int_moor_traits_f_8 <- filter(int_moor_traits_f, n > 7)
length(unique(int_moor_traits_f_8$species)) # 69 species
length(unique(int_moor_traits_f_8$family)) # 17 families

# The dataset has been largely reduced

# Download the phylogeny from The Fish Tree of Life 
# and create a covariance matrix for brms models

# Check species names
check_name_fishtree(unique(as.character(int_moor_traits_f_8$species)), sampled = TRUE)
# All species names are correct
# All species have genetic data

# Download phylogeny
tree_int_8 <- fishtree_phylogeny(species = unique(int_moor_traits_f_8$species))

# Phylogenetic covariance matrix 
corphy_int_8 <- phylo_cov(tree_int_8)

# Run models

# INTESTINAL LENGTH
m_leng_8 <- brm(il_log ~ scale(spec_mean_sl_log) + scale(troph) + scale(within_spec_sl_log) + scale(elon_log) + stomach*durophagy +
                  (1 |gr(phylo, cov = corphy_int_8)) + (1 |species) + (0 + scale(within_spec_sl_log) |species),
                data = int_moor_traits_f_8, family = student(), data2 = list(corphy_int_8 = corphy_int_8),
                prior = prior, warmup = 2000, iter = 8000, chains = 4, cores = 4,
                inits = "random", control = list(adapt_delta = 0.95, max_treedepth = 15), file = "Output/Models/m_leng_8")

# Checks and summary
pp_check(m_leng_8, type = "dens_overlay", nsamples = 100) # ok
plot(m_leng_8, ask = FALSE) # ok
loo(m_leng_8) # All pareto_k < 0.5
summary(m_leng_8, priors = TRUE)
bayes_R2(m_leng_8) # 93.1% [92.7%, 93.4%]

# Get unscaled slopes
round(fixef(m_leng_8)[2,]/sd(int_moor_traits_f_8$spec_mean_sl_log), 3) # 0.837[0.501, 1.164]
round(fixef(m_leng_8)[3,]/sd(int_moor_traits_f_8$troph), 3) # -0.318[-0.535, -0.125]
round(fixef(m_leng_8)[4,]/sd(int_moor_traits_f_8$within_spec_sl_log), 3) # 0.835[0.693, 0.974]
round(fixef(m_leng_8)[5,]/sd(int_moor_traits_f_8$elon_log), 3) # -0.802[-1.188, -0.422]

# INTESTINAL DIAMETER
m_diam_8 <- brm(id_log ~ scale(spec_mean_sl_log) + scale(troph) + scale(within_spec_sl_log) + scale(elon_log) + stomach*durophagy +
                  (1 |gr(phylo, cov = corphy_int_8)) + (1 |species) + (0 + scale(within_spec_sl_log) |species),
                data = int_moor_traits_f_8, family = student(), data2 = list(corphy_int_8 = corphy_int_8),
                prior = prior, warmup = 2000, iter = 8000, chains = 4, cores = 4,
                inits = "random", control = list(adapt_delta = 0.95, max_treedepth = 15), file = "Output/Models/m_diam_8")

# Checks and summary
pp_check(m_diam_8, type = "dens_overlay", nsamples = 100) # ok
plot(m_diam_8, ask = FALSE) # ok
loo(m_diam_8) # All pareto_k < 0.5
summary(m_diam_8, priors = TRUE)
bayes_R2(m_diam_8) # 82.2% [80.9%, 83.3%]

# Get unscaled slopes
round(fixef(m_diam_8)[2,]/sd(int_moor_traits_f_8$spec_mean_sl_log), 3) # 0.995[0.747, 1.245]
round(fixef(m_diam_8)[3,]/sd(int_moor_traits_f_8$troph), 3) # -0.197[-0.322, -0.051]
round(fixef(m_diam_8)[4,]/sd(int_moor_traits_f_8$within_spec_sl_log), 3) # 0.699[0.567, 0.828]
round(fixef(m_diam_8)[5,]/sd(int_moor_traits_f_8$elon_log), 3) # -0.411[-0.643, -0.186]

# INTESTINAL SURFACE
m_surf_8 <- brm(is_log ~ scale(spec_mean_sl_log) + scale(troph) + scale(within_spec_sl_log) + scale(elon_log) + stomach*durophagy +
                  (1 |gr(phylo, cov = corphy_int_8)) + (1 |species) + (0 + scale(within_spec_sl_log) |species),
                data = int_moor_traits_f_8, family = student(), data2 = list(corphy_int_8 = corphy_int_8),
                prior = prior, warmup = 2000, iter = 8000, chains = 4, cores = 4,
                inits = "random", control = list(adapt_delta = 0.95, max_treedepth = 15), file = "Output/Models/m_surf_8")

# Checks and summary
pp_check(m_surf_8, type = "dens_overlay", nsamples = 100) # ok
plot(m_surf_8, ask = FALSE) # ok
loo(m_surf_8) # All pareto_k < 0.7
summary(m_surf_8, priors = TRUE)
bayes_R2(m_surf_8) # 90.4% [89.8%, 90.9%]

# Get unscaled slopes
round(fixef(m_surf_8)[2,]/sd(int_moor_traits_f_8$spec_mean_sl_log), 3) # 1.770[1.301, 2.231]
round(fixef(m_surf_8)[3,]/sd(int_moor_traits_f_8$troph), 3) # -0.481[-0.848, -0.165]
round(fixef(m_surf_8)[4,]/sd(int_moor_traits_f_8$within_spec_sl_log), 3) # 1.529[1.305, 1.746]
round(fixef(m_surf_8)[5,]/sd(int_moor_traits_f_8$elon_log), 3) # -1.213[-1.747, -0.700]


# Summary table of the 3 models (Appendix S2: Table 4)
int_model_8_summary <- left_join(as_tibble(posterior_summary(m_leng_8)[1:13, -2], rownames = "parameter"), 
                                 as_tibble(posterior_summary(m_diam_8)[1:13, -2], rownames = "parameter"), 
                                 by = "parameter", suffix = c(".l", ".d")) %>%
  left_join(as_tibble(posterior_summary(m_surf_8)[1:13, -2], rownames = "parameter") %>%
              rename_with( ~ paste0(.x, ".s"), .cols = 2:4),
            by = "parameter") %>%
  mutate(across(-1, round, 2))

# Save table
write.table(int_model_8_summary, file = "Output/Tables/Summary intestine models sensitivity 8.txt", sep = ";", row.names = FALSE)

#### ------------ END ------------ ####