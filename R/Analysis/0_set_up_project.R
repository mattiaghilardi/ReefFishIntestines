#### ------------ This script sets up the project ------------ ####

# 1 - create "Output" folder and subfolders
# 2 - load the required packages
# 3 - load the supporting functions
# 4 - set theme for all plots

#### ------------ Set up output folders ------------ ####

dir.create("Output")
dir.create("Output/Models")
dir.create("Output/Trees")
dir.create("Output/Figures")
dir.create("Output/Tables")

#### ------------ Load all packages and dependencies required for the analysis ------------ ####

# Packages' list
pkg_list <- c("dplyr", "tidyr", "forcats", "rfishbase", "parallel", 
              "MCMCglmm", "ape", "fishtree", "brms", "tidybayes", "bayesplot",
              "extrafont", "ggplot2", "ggpubr", "ggtree", "ggtext", "ggridges",
              "RColorBrewer", "ggnewscale", "cowplot", "fishualize", "ggrepel", "ggimage")

# Install packages if not already installed
new_pkgs <- pkg_list[!(pkg_list %in% installed.packages()[,"Package"])]
if(length(new_pkgs) > 0) install.packages(new_pkgs)

# Load packages
lapply(pkg_list, require, character.only = TRUE)

#### ------------ Source supporting functions ------------ ####

functions <- list.files("R/Functions")
sapply(functions, function(x) source(paste0("R/Functions/", x)))

#### ------------ Set plot's theme ------------ ####

font_import() # only do this one time - it takes few minutes
loadfonts(device = "win")
theme_set(theme_bw(base_size = 12, base_family = "serif") + 
            theme(panel.grid = element_blank(),
                  axis.text = element_text(color = "black")))

#### ------------ END ------------ ####
