# Phylogenetic covariance matrix for class 'phylo' or 'multiPhylo'
# Code written by Nina MD Schiettekatte and modified by Mattia Ghilardi

# Input: 
# - phy: object of class 'phylo' or 'multiPhylo'
# - scale: LOGICAL; if the matrix has to be normalised to a scale 0-1 (defaults to TRUE) 
# - mc.cores: the number of cores to use for parallel::mclapply(), i.e. at most how many child processes will be run simultaneously. Must be exactly 1 on Windows (which uses the master process); default to 1

# Value: 
# - class 'phylo' --> a numeric matrix with the names of the tips as colnames and rownames
# - class 'multiphylo' --> a list with the mean and standard deviation covariance matrices

phylo_cov <- function(phy, scale = TRUE, mc.cores = 1) {
  
  # Check tree 
  if (!inherits(phy, "phylo") & !inherits(phy, "multiPhylo")) {
    stop("Tree must be of class 'phylo' or 'multiPhylo'")
  }
  
  # Covariance matrix for phylo object
  if (inherits(phy, "phylo")) {
    phy2 <- phy %>% ape::chronoMPL()
    inv.phylo <- MCMCglmm::inverseA(phy2, nodes = "TIPS", scale = scale)
    A <- base::solve(inv.phylo$Ainv)
    rownames(A) <- rownames(inv.phylo$Ainv)
    A[sort(rownames(A)), sort(rownames(A))]
    return(A)
  } 
  
  # Covariance matrix for multiphylo object
  if (inherits(phy, "multiPhylo")) {
    
    # Covariance matrix for each tree
    corphy_list <- parallel::mclapply(phy, function(x){
      phy2 <- x %>% ape::chronoMPL()
      inv.phylo <- MCMCglmm::inverseA(phy2, nodes = "TIPS", scale = scale)
      A <- base::solve(inv.phylo$Ainv)
      rownames(A) <- rownames(inv.phylo$Ainv)
      A[sort(rownames(A)), sort(rownames(A))]
      return(A)
    }, mc.cores = mc.cores)
    
    # Transform list of covariance matrices to 3d array
    corphy_array <- simplify2array(corphy_list)
    
    # Summarise these matrices
    corphy_m <- apply(corphy_array, 1:2, mean)
    corphy_sd <- apply(corphy_array, 1:2, sd)
    
    # return list
    list(mean = corphy_m, sd = corphy_sd)
  }
}
