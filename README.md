# Exhaustive ADDIS procedures

This repository contains all the code used for the simulations and real data applications in the paper "Exhaustive ADDIS procedures". In the aforementioned paper, we introduced exhaustive ADDIS procedures for online multiple testing that uniformly improve the existing ADDIS procedures by Tian & Ramdas (2021). All procedures control the FWER strongly under independence of the null $p$-values. 



Files:

exhaustive_procedures.R: Contains all applied online procedures. 

plot_generating_function.R: Implements the simulation setting of Section 5.

plot_creator.R: Generates Figures 3-5 of the paper and saves the results in the results folder. 

plot_real_data.R: Real data application based on the data in IMPC_data.R. Generates Figure 6 of the paper and saves it to the results folder.

IMPC_data.csv: Contains 5000 observations from a full data set available at a Zenodo repository organized by Robertson et al. (2018) (https://doi.org/10.5281/zenodo.2396572). 

results folder: Contains all the results that were provided in the paper.



Version information:
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22621)

Matrix products: default

locale:
[1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8   
[3] LC_MONETARY=German_Germany.utf8 LC_NUMERIC=C                   
[5] LC_TIME=German_Germany.utf8    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] patchwork_1.1.2 ggplot2_3.4.2  

loaded via a namespace (and not attached):
 [1] fansi_1.0.4      withr_2.5.0      dplyr_1.1.3      utf8_1.2.3       grid_4.2.2      
 [6] R6_2.5.1         lifecycle_1.0.3  gtable_0.3.3     magrittr_2.0.3   scales_1.2.1    
[11] pillar_1.9.0     rlang_1.1.1      cli_3.6.0        rstudioapi_0.14  generics_0.1.3  
[16] vctrs_0.6.4      tools_4.2.2      glue_1.6.2       munsell_0.5.0    compiler_4.2.2  
[21] pkgconfig_2.0.3  colorspace_2.1-0 tidyselect_1.2.0 tibble_3.2.1  

