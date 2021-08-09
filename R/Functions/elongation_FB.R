# Get elongation from FishBase
# Get average elongation at species, genus or family level

# Input:
# - species_list: a vector of taxonomic names in the form "Genus species"

# Value: 
# A tibble with 3 columns: 
# - "species": the taxonomic name
# - "elongation": the value retrieved from FishBase
# - "level_elon": the taxonomic level at which the value has been retrieved

elongation_FB <- function(species_list) {
  
  # Check species names
  sp_correct <- suppressWarnings(rfishbase::validate_names(species_list))
  sp_error   <- species_list[!(species_list %in% sp_correct)]
  if (length(sp_error) > 0) {
    stop("Species name is incorrect: ", paste(sp_error, collapse = ", "))
  }
  
  # Get elongation at the species, genus or family level for each species
  all_elon <- lapply(species_list, function(x){
    
    # Retrieve taxonomic info
    taxo <- rfishbase::load_taxa()
    taxo <- as.data.frame(taxo)

    # Get elongation for all species within the family
    fam  <- taxo[taxo$Species==x, "Family"]
    fam_all  <- rfishbase::species_list(Family = fam)
    elon <- rfishbase::morphometrics(fam_all, fields = c("Species", "SL", "BD"))
    elon$SL <- as.numeric(elon$SL)
    elon <- dplyr::mutate(elon, elon = SL/BD)
    
    # Join taxonomic info to the aspect ratios
    elon <- suppressMessages(dplyr::left_join(elon, taxo[, c("Family", "Genus", "Species")]))
    
    # Compute species, genus and family level aspect ratio
    elon <- dplyr::group_by(elon, Species)
    elon <- dplyr::mutate(elon, elon_sp = mean(elon, na.rm = TRUE))
    elon <- dplyr::ungroup(elon)
    elon <- dplyr::group_by(elon, Genus)
    elon <- dplyr::mutate(elon, elon_gen = mean(elon_sp, na.rm = TRUE))
    elon <- dplyr::ungroup(elon)
    elon$elon_fam <- mean(elon$elon_sp, na.rm = TRUE)
    
    # Retain only the desired species
    elon <- elon[elon$Species==x,]
    
    # Retain only the first row
    elon <- elon[1,]
    
    # Retain the elongation at the lowest possible taxonomic level
    # Add taxonomic level of elongation's value
    elon <- if (!is.na(elon$elon_sp)) {
      dplyr::mutate(elon,
                    elongation = elon_sp,
                    level = "species")
    } else if (is.na(elon$elon_sp) & !is.na(elon$elon_gen)) {
      dplyr::mutate(elon,
                    elongation = elon_gen,
                    level = "genus")
    } else {
      dplyr::mutate(elon,
                    elongation = elon_fam,
                    level = "family")
    }
    
    # Round elongation
    elon$elongation <- round(elon$elongation, 2)
    
    # Tibble for each species
    dplyr::tibble(species = x,
                  elongation = elon$elongation,
                  level_elon = elon$level)
    
  })
  
  # List to tibble
  all_elon <- dplyr::bind_rows(all_elon)
  all_elon
  
}
