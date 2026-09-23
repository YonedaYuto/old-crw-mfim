# R session information

The analyses were run in R 4.6.1. The two sessions below were printed by
`utils::sessionInfo()` at the end of the scripts named in each heading. They
differ in platform: the first script was run on Ubuntu 24.04.4 LTS and the
second on Windows 11 x64. The scripts also print session information at the end
of their logs in `output/`, which is not included in this repository.

## 1. `script/07_1_sens_ranktest.R` (rank-based sensitivity analysis)

- Platform: Ubuntu 24.04.4 LTS (x86_64-pc-linux-gnu)
- Attached package: rankFD 0.1.1

```text
R version 4.6.1 (2026-06-24)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 24.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C          
 [3] LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8    
 [5] LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C             
 [9] LC_ADDRESS=C           LC_TELEPHONE=C        
[11] LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

time zone: UTC
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods  
[7] base     

other attached packages:
[1] rankFD_0.1.1

loaded via a namespace (and not attached):
 [1] codetools_0.2-20  multcomp_1.4-32   Matrix_1.7-5     
 [4] lattice_0.22-9    TH.data_1.1-5     splines_4.6.1    
 [7] matrixStats_1.5.0 modeltools_0.2-24 zoo_1.9-0        
[10] libcoin_1.0-13    parallel_4.6.1    stats4_4.6.1     
[13] mvtnorm_1.4-2     grid_4.6.1        sandwich_3.1-3   
[16] compiler_4.6.1    tools_4.6.1       coin_1.4-5       
[19] survival_3.8-6    MASS_7.3-65  
```

## 2. `script/12_figures_tables.R` (figures and tables)

- Platform: Windows 11 x64 (build 26200)
- Attached packages: ordinal 2026.7-26, ggplot2 4.0.3, patchwork 1.3.2 and splines (base R)

```text
R version 4.6.1 (2026-06-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=Japanese_Japan.utf8  LC_CTYPE=Japanese_Japan.utf8    LC_MONETARY=Japanese_Japan.utf8
[4] LC_NUMERIC=C                    LC_TIME=Japanese_Japan.utf8    

time zone: Asia/Tokyo
tzcode source: internal

attached base packages:
[1] splines   stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] patchwork_1.3.2   ggplot2_4.0.3     ordinal_2026.7-26

loaded via a namespace (and not attached):
 [1] generics_0.1.4          tidyr_1.3.2             fontLiberation_0.1.0    xml2_1.6.0             
 [5] lattice_0.22-9          digest_0.6.39           magrittr_2.0.5          evaluate_1.0.5         
 [9] grid_4.6.1              RColorBrewer_1.1-3      flextable_0.10.1        fastmap_1.2.0          
[13] Matrix_1.7-5            zip_3.0.2               purrr_1.2.2             scales_1.4.0           
[17] fontBitstreamVera_0.1.1 numDeriv_2016.8-1.1     textshaping_1.0.5       cli_3.6.6              
[21] rlang_1.2.0             fontquiver_0.2.1        withr_3.0.3             otel_0.2.0             
[25] gdtools_0.5.1           parallel_4.6.1          tools_4.6.1             officer_0.7.6          
[29] ucminf_1.2.3            uuid_1.2-2              dplyr_1.2.1             vctrs_0.7.3            
[33] R6_2.6.1                stats4_4.6.1            lifecycle_1.0.5         magick_2.9.1           
[37] MASS_7.3-65             ragg_1.5.2              pkgconfig_2.0.3         pillar_1.11.1          
[41] gtable_0.3.6            data.table_1.18.6.1     glue_1.8.1              Rcpp_1.1.2             
[45] systemfonts_1.3.2       xfun_0.60               tibble_3.3.1            tidyselect_1.2.1       
[49] rstudioapi_0.19.0       knitr_1.52              farver_2.1.2            nlme_3.1-169           
[53] htmltools_0.5.9         labeling_0.4.3          rmarkdown_2.32          compiler_4.6.1         
[57] S7_0.2.2                askpass_1.2.1           openssl_2.4.2      
```
