# A function to check species names and sampled species in the Fish Tree of Life

# Input: 
# - species_list: a list of species in the form "Genus species"
# - sampled: LOGICAL; if the list has to be checked for unsampled species (i.e. those without genetic data, see Rabosky et al. 2018 for detailed methodology) 

# Value: a message and, if sampled=TRUE, a vector of species without genetic data

check_name_fishtree <- function(species_list, sampled = FALSE) {
  
  # Retrieve all species in the Fish tree of Life
  sp_fishtree <- fishtree::fishtree_taxonomy("Actinopteri")[[1]]
  
  # Check names
  sp_error <- species_list[!(species_list %in% sp_fishtree$species)]
  
  if (length(sp_error) == 0) {
    message("All species names are correct")
  } else {
    stop("Species name is incorrect or not present in the Fish Tree of Life: ", paste(sp_error, collapse = ", "))
  }
  
  # If sampled=TRUE check sampled species
  if (isTRUE(sampled)) {
    sp_not_sampled <- species_list[!(species_list %in% sp_fishtree$sampled_species)]
    if (length(sp_not_sampled) == 0) {
      message("All species have genetic data")
    } else {
      message("These species do not have genetic data:")
      sp_not_sampled
    }
  }
  
}
