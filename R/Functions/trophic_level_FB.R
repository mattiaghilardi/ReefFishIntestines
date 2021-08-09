# Get trophic level from FishBase
# Borrow some code from fishflux::trophic_level
# Get trophic level at species, genus or family level

# Input:
# - species_list: a vector of taxonomic names in the form "Genus species"
# - type: one of "diet", "food items", "both".
#         - "diet" --> for trophic level based on diet studies only. These are available only for some species. Note that for some families there are no diet studies available.
#         - "food items" --> for trophic level based on food items only. 
#         - "both" --> use the available value, "diet" or "food items", or the mean of the two when both are available.

# Value:
# A tibble with 3 columns: 
# - "species": the taxonomic name
# - "troph": the value retrieved from FishBase
# - "level_troph": the taxonomic level at which the value has been retrieved


trophic_level_FB <- function(species_list, type = c("diet", "food items", "both")) {
  
  # Check type
  if (length(type) > 1) {
    stop("'type' must be one")
  }
  
  if (type != "diet" & type != "food items" & type != "both") {
    stop("'type' is incorrect. Please use one of 'diet', 'food items' or 'both'")
  }
  
  # Check species names
  sp_correct <- suppressWarnings(rfishbase::validate_names(species_list))
  sp_error   <- species_list[!(species_list %in% sp_correct)]
  if (length(sp_error) > 0) {
    stop("Species name is incorrect: ", paste(sp_error, collapse = ", "))
  }
  
  # Get trophic level at the species, genus or family level for each species
  all_troph <- lapply(species_list, function(x){
    
    # Retrieve taxonomic info
    taxo <- rfishbase::load_taxa()
    taxo <- as.data.frame(taxo)
    
    # Get trophic level for all species within the family
    fam <- taxo[taxo$Species==x, "Family"]
    fam_all <- rfishbase::species_list(Family = fam)
    tl  <- if (type == "diet") {
      rfishbase::ecology(fam_all, fields = c("Species", "DietTroph"))
    } else if (type == "food items") {
      rfishbase::ecology(fam_all, fields = c("Species", "FoodTroph"))
    } else {
      rfishbase::ecology(fam_all, fields = c("Species", "DietTroph", "FoodTroph"))
    }
    
    # Join taxonomic info to the trophic levels
    tl <- suppressMessages(dplyr::left_join(tl, taxo[, c("Family", "Genus", "Species")]))
    
    # Compute species, genus and family level trophic lavel
    tl$troph_sp <- if (type == "diet") {
      tl$DietTroph
    } else if (type == "food items") {
      tl$FoodTroph
    } else {
      rowMeans(tl[, c("DietTroph", "FoodTroph")], na.rm = TRUE)
    }
    
    tl <- dplyr::group_by(tl, Genus)
    tl <- dplyr::mutate(tl, troph_gen = mean(troph_sp, na.rm = TRUE))
    tl <- dplyr::ungroup(tl)
    tl$troph_fam <- mean(tl$troph_sp, na.rm = TRUE)
    
    # Retain only the desired species
    tl <- tl[tl$Species==x,]
    
    # Some species have different entries for different stocks of the same species
    # Retain only the first
    # Data is either identical to the first or simply missing in the additional stocks 
    tl <- tl[1,]
    
    # Retain the trophic level at the lowest possible taxonomic level
    # Add taxonomic level of trophic level's value
    tl <- if (!is.na(tl$troph_sp)) {
      dplyr::mutate(tl,
                    troph = troph_sp,
                    level = "species")
      } else if (is.na(tl$troph_sp) & !is.na(tl$troph_gen)) {
        dplyr::mutate(tl,
                      troph = troph_gen,
                      level = "genus")
        } else {
          dplyr::mutate(tl,
                        troph = troph_fam,
                        level = "family")
        }
    
    # Round trophic level
    tl$troph <- round(tl$troph, 2)
    
    # Tibble for each species
    dplyr::tibble(species = x,
                  troph = tl$troph,
                  level_troph = tl$level)
    
  })
  
  # List to tibble
  all_troph <- dplyr::bind_rows(all_troph)
  all_troph
  
}
 
