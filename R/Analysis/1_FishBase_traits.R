#### --------- This script downloads the trait data from FishBase --------- ####


#### ------------ Load the dataset ------------ ####

int_moor <- read.csv("Data/intestine_fish_Moorea.csv")

#### ------------ Download data from FishBase ------------ ####

# Select species in the dataset
species <- unique(int_moor$species)

# Validate species names
species_val <- validate_names(species)

# Check invalid species names
setdiff(species, species_val)
# "Myripristis speA"   "Centropyge loricula" "Acanthurus nigros"  

# These 3 species have non-valid names in FishBase

# We can drop the first but we retain the other two
int_moor <- int_moor[!(int_moor$species == "Myripristis speA"),]

# Acanthurus nigros is considered different from A. nigroris following Randall et al. 2011
# Centropyge loricula is accepted as C. loriculus
# Rename these two species to get data from FishBase, need this later for the phylogeny as well
int_moor$species[int_moor$species == "Centropyge loricula"] <- "Centropyge loriculus"
int_moor$species[int_moor$species == "Acanthurus nigros"] <- "Acanthurus nigroris"

# Reselect species in the dataset
species <- levels(factor(int_moor$species))

# Validate species names
species_val <- validate_names(species)

# Check invalid species names
setdiff(species, species_val) # ok

# TROPHIC LEVEL

# Retrieve trophic level data
troph <- ecology(species, fields = c("Species", "FoodTroph", "DietTroph"))

# Check for NAs
nrow(troph[is.na(troph$FoodTroph),]) # 22 species
nrow(troph[is.na(troph$DietTroph),]) # 115 species

# Use FoodTroph values and replace NAs with the genus-level or family-level average of FoodTroph
troph <- trophic_level_FB(species, type = "food items")

# ELONGATION

# Retrieve elongation data (standard length/body depth)
elongation <- elongation_FB(species)

# MAX LENGTH

# Retrieve max length data
max_length <- estimate(species, fields = c("Species", "MaxLengthTL")) %>% 
  rename(species = Species)

# Join traits to the dataset
int_moor_traits <- left_join(int_moor, troph) %>% left_join(elongation) %>% left_join(max_length)

# Save csv
write.csv(int_moor_traits, "Data/intestine_fish_Moorea_FB.traits.csv", row.names = FALSE)

#### ------------ END ------------ ####
