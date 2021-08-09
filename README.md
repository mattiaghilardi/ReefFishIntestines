# Phylogeny, body morphology, and trophic level shape intestinal traits in coral reef fishes

This repository holds all the code and data needed to reproduce the manuscript:

Ghilardi, M., Schiettekatte, N. M. D., Casey, J. M., Brandl, S. J., Degregori, S., Mercière, A., Morat, F., Letourneur, Y., Bejarano, S., Parravicini, V. (2021) Phylogeny, body morphology, and trophic level shape intestinal traits in coral reef fishes. *Ecology and Evolution* (in press).

When using code or data from this project, please cite it as:

Ghilardi, M., Schiettekatte, N. M. D., Casey, J. M., Brandl, S. J., Degregori, S., Mercière, A., Morat, F., Letourneur, Y., Bejarano, S., Parravicini, V. (2021) Data and code of accepted version of manuscript: Phylogeny, body morphology, and trophic level shape intestinal traits in coral reef fishes (Ecology and Evolution). *Zenodo*. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5172790.svg)](https://doi.org/10.5281/zenodo.5172790)



# Instructions

All modelling and plots were produced in R using [Stan](https://mc-stan.org/) through the [brms](https://paul-buerkner.github.io/brms/) package.

If you would like to replicate our analysis, run the code in the folder [`R/Analysis`](https://github.com/mattiaghilardi/ReefFishIntestines/blob/main/R/Analysis) in this order:

- [0_set_up_project.R](https://github.com/mattiaghilardi/ReefFishIntestines/blob/main/R/Analysis/0_set_up_project.R)
- [1_FishBase_traits.R](https://github.com/mattiaghilardi/ReefFishIntestines/blob/main/R/Analysis/1_FishBase_traits.R)
- [2_isotope_analysis.R](https://github.com/mattiaghilardi/ReefFishIntestines/blob/main/R/Analysis/2_isotope_analysis.R)
- [3_intestine_analysis.R](https://github.com/mattiaghilardi/ReefFishIntestines/blob/main/R/Analysis/3_intestine_analysis.R)
- [4_sensitivity_analysis.R](https://github.com/mattiaghilardi/ReefFishIntestines/blob/main/R/Analysis/4_sensitivity_analysis.R)
- [5_figures.R](https://github.com/mattiaghilardi/ReefFishIntestines/blob/main/R/Analysis/5_figures.R)

This will reproduce the entire analysis from start to finish.

The output will be saved in four separate folders:

📁 **Models**: for the Bayesian models   
📁 **Trees**: for the phylogenetic trees  
📁 **Figures**: for all the figures in the manuscript and supporting information  
📁 **Tables**: for all supplementary tables

Note that running the script `1_FishBase_traits.R` may generate slightly different values for the traits from those used in the paper due to updates in FishBase. Also, taxonomic names are continuously revised in FishBase and this may lead to errors when running the script. The complete dataset used in the paper, including traits' values, is included in `Data` ([*intestine_fish_Moorea_FB.traits.csv*](https://github.com/mattiaghilardi/ReefFishIntestines/blob/main/Data/intestine_fish_Moorea_FB.traits.csv)).

## The following software and associated packages were used:

```{r}
R version 4.0.2 (2020-06-22)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19043)

Matrix products: default

locale:
[1] LC_COLLATE=English_Germany.1252  LC_CTYPE=English_Germany.1252   
[3] LC_MONETARY=English_Germany.1252 LC_NUMERIC=C                    
[5] LC_TIME=English_Germany.1252    

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] ggimage_0.2.8      ggrepel_0.8.2      fishualize_0.2.0   cowplot_1.0.0     
 [5] ggnewscale_0.4.3   RColorBrewer_1.1-2 ggridges_0.5.2     ggtext_0.1.0      
 [9] ggtree_2.2.4       ggpubr_0.4.0       ggplot2_3.3.2      extrafont_0.17    
[13] bayesplot_1.7.2    tidybayes_2.1.1    brms_2.14.4        Rcpp_1.0.7        
[17] fishtree_0.3.2     MCMCglmm_2.29      ape_5.4-1          coda_0.19-3       
[21] Matrix_1.2-18      rfishbase_3.0.4    forcats_0.5.0      tidyr_1.1.2       
[25] dplyr_1.0.2       

loaded via a namespace (and not attached):
  [1] readxl_1.3.1         backports_1.1.9      plyr_1.8.6          
  [4] igraph_1.2.5         lazyeval_0.2.2       splines_4.0.2       
  [7] svUnit_1.0.3         crosstalk_1.1.0.1    rstantools_2.1.1    
 [10] inline_0.3.15        digest_0.6.25        htmltools_0.5.1.1   
 [13] magick_2.4.0         rsconnect_0.8.16     fansi_0.4.1         
 [16] magrittr_1.5         memoise_1.1.0        openxlsx_4.1.5      
 [19] readr_1.3.1          RcppParallel_5.0.2   matrixStats_0.56.0  
 [22] xts_0.12-0           extrafontdb_1.0      prettyunits_1.1.1   
 [25] colorspace_1.4-1     ggdist_2.2.0         haven_2.3.1         
 [28] xfun_0.16            callr_3.5.1          crayon_1.3.4        
 [31] jsonlite_1.7.1       lme4_1.1-23          zoo_1.8-8           
 [34] glue_1.4.2           gtable_0.3.0         V8_3.2.0            
 [37] distributional_0.2.0 car_3.0-9            pkgbuild_1.1.0      
 [40] Rttf2pt1_1.3.8       rstan_2.21.2         abind_1.4-5         
 [43] scales_1.1.1         mvtnorm_1.1-1        rstatix_0.6.0       
 [46] miniUI_0.1.1.1       gridtext_0.1.1       xtable_1.8-4        
 [49] gridGraphics_0.5-0   tidytree_0.3.3       foreign_0.8-80      
 [52] stats4_4.0.2         StanHeaders_2.21.0-6 DT_0.15             
 [55] htmlwidgets_1.5.1    httr_1.4.2           threejs_0.3.3       
 [58] arrayhelpers_1.1-0   ellipsis_0.3.1       pkgconfig_2.0.3     
 [61] loo_2.3.1            farver_2.0.3         ggplotify_0.0.5     
 [64] tidyselect_1.1.0     rlang_0.4.11         reshape2_1.4.4      
 [67] later_1.1.0.1        munsell_0.5.0        cellranger_1.1.0    
 [70] tools_4.0.2          cli_2.0.2            generics_0.1.0      
 [73] broom_0.7.0          stringr_1.4.0        fastmap_1.0.1       
 [76] processx_3.4.5       zip_2.1.1            purrr_0.3.4         
 [79] gh_1.1.0             nlme_3.1-148         mime_0.9            
 [82] projpred_2.0.2       aplot_0.0.5          xml2_1.3.2          
 [85] compiler_4.0.2       shinythemes_1.1.2    rstudioapi_0.11     
 [88] png_0.1-7            curl_4.3             gamm4_0.2-6         
 [91] ggsignif_0.6.0       treeio_1.12.0        tibble_3.0.3        
 [94] statmod_1.4.34       stringi_1.4.6        ps_1.3.4            
 [97] Brobdingnag_1.2-6    cubature_2.0.4.1     lattice_0.20-41     
[100] nloptr_1.2.2.2       markdown_1.1         shinyjs_1.1         
[103] tensorA_0.36.1       vctrs_0.3.3          pillar_1.4.6        
[106] lifecycle_0.2.0      BiocManager_1.30.10  bridgesampling_1.0-0
[109] data.table_1.13.0    corpcor_1.6.9        patchwork_1.0.1     
[112] httpuv_1.5.4         R6_2.4.1             promises_1.1.1      
[115] gridExtra_2.3        rio_0.5.16           codetools_0.2-16    
[118] boot_1.3-25          colourpicker_1.0     MASS_7.3-51.6       
[121] gtools_3.8.2         assertthat_0.2.1     withr_2.2.0         
[124] shinystan_2.5.0      mgcv_1.8-31          hms_0.5.3           
[127] grid_4.0.2           minqa_1.2.4          rvcheck_0.1.8       
[130] carData_3.0-4        shiny_1.5.0          base64enc_0.1-3     
[133] dygraphs_1.1.1.6     tinytex_0.25
```

## How to download this project for people not familiar with GitHub:

Click on the green button `clone or download` in the [project's main page on GitHub](https://github.com/mattiaghilardi/ReefFishIntestines) and then click on `Download ZIP`.

# Acknowledgements
Thanks to Clause Wilke for the function [align_legend()](https://github.com/clauswilke/dviz.supp/blob/master/R/align_legend.R) that helped adjusting legends' alignment in Fig. 2.

# Bug reporting
Please [report any bug and issue](https://github.com/mattiaghilardi/ReefFishIntestines/issues)
