Common Variant GWAS MAGMA Preparation
================
Amalya Murrill
2026-04-18

This script takes the downloaded Grove et al. 2019 ASD GWAS summary
statistics and prepares them for the MAGMA (Multi-marker Analysis of
GenoMic Annotation) calculation of gene-level association scores. QC is
run to determine if claims made in the methods of the paper match the
data output. Indels are handled poorly by MAGMA and are thus filtered
from the data. Finally, the required columns are reformatted in the
correct order and written in a text file for MAGMA. Dependencies and
other session info are listed below. Here, the common variant arm of the
pipeline is handled, while the rare variants are dealt with separately.

``` r
knitr::opts_chunk$set(warning = FALSE, message = FALSE)
```

## Libraries

``` r
library(tidyverse)
```

## Data Import

``` r
#read in
common <- read_tsv("../data/raw/iPSYCH-PGC_ASD_Nov2017.gz", col_names = TRUE)
head(common)
```

    ## # A tibble: 6 × 9
    ##     CHR SNP               BP A1    A2     INFO    OR     SE       P
    ##   <dbl> <chr>          <dbl> <chr> <chr> <dbl> <dbl>  <dbl>   <dbl>
    ## 1     8 rs62513865 101592213 T     C     0.949 1.01  0.027  0.809  
    ## 2     8 rs79643588 106973048 A     G     0.997 1.02  0.024  0.461  
    ## 3     8 rs17396518 108690829 T     G     0.987 0.961 0.014  0.00465
    ## 4     8 rs983166   108681675 A     C     0.998 0.980 0.0139 0.145  
    ## 5     8 rs28842593 103044620 T     C     0.857 0.996 0.0203 0.842  
    ## 6     8 rs35107696 109712249 A     AT    0.999 1.01  0.0165 0.430

## Quality Control

``` r
#num SNPs
nrow(common) #9112386 total SNPs
```

    ## [1] 9112386

``` r
#autosomes only?
common %>% distinct(CHR) %>% nrow() #22 chromosomes
```

    ## [1] 22

``` r
common %>% distinct(CHR) %>% arrange(CHR) #each present
```

    ## # A tibble: 22 × 1
    ##      CHR
    ##    <dbl>
    ##  1     1
    ##  2     2
    ##  3     3
    ##  4     4
    ##  5     5
    ##  6     6
    ##  7     7
    ##  8     8
    ##  9     9
    ## 10    10
    ## # ℹ 12 more rows

``` r
#duplicates?
common %>% filter(duplicated(SNP)) %>% nrow() #0, none
```

    ## [1] 0

``` r
#significant hits?
common %>% filter(P < 5e-8) %>% nrow() #93
```

    ## [1] 93

``` r
#INFO >= 0.7?
common %>% filter(INFO < 0.7) %>% nrow() #0, none
```

    ## [1] 0

``` r
#num NA
sum(is.na(common)) #0, no NA
```

    ## [1] 0

## Indel Removal

``` r
common_clean <- common %>%
  filter(nchar(A1) == 1, nchar(A2) == 1)
nrow(common_clean) #8058797 single nucleotide SNPs
```

    ## [1] 8058797

``` r
common_clean
```

    ## # A tibble: 8,058,797 × 9
    ##      CHR SNP               BP A1    A2     INFO    OR     SE       P
    ##    <dbl> <chr>          <dbl> <chr> <chr> <dbl> <dbl>  <dbl>   <dbl>
    ##  1     8 rs62513865 101592213 T     C     0.949 1.01  0.027  0.809  
    ##  2     8 rs79643588 106973048 A     G     0.997 1.02  0.024  0.461  
    ##  3     8 rs17396518 108690829 T     G     0.987 0.961 0.014  0.00465
    ##  4     8 rs983166   108681675 A     C     0.998 0.980 0.0139 0.145  
    ##  5     8 rs28842593 103044620 T     C     0.857 0.996 0.0203 0.842  
    ##  6     8 rs7014597  104152280 C     G     0.993 1.02  0.0182 0.293  
    ##  7     8 rs3134156  100479917 T     C     0.998 0.982 0.019  0.353  
    ##  8     8 rs6980591  103144592 A     C     0.997 1.04  0.0167 0.00839
    ##  9     8 rs72670434 108166508 A     T     0.985 1.01  0.0147 0.390  
    ## 10     8 rs10955343 105201080 T     C     1     0.994 0.0144 0.659  
    ## # ℹ 8,058,787 more rows

## MAGMA Formatting

``` r
common_MAGMA <- common_clean %>%
  mutate(NOBS = INFO*46350) %>%  #Total sample size from Grove et al. 2019 (n = 46,350, from Fig. 2 caption) using R-squared value for per-SNP N-effect
  dplyr::select(SNP, CHR, BP, P, NOBS)

common_MAGMA
```

    ## # A tibble: 8,058,797 × 5
    ##    SNP          CHR        BP       P   NOBS
    ##    <chr>      <dbl>     <dbl>   <dbl>  <dbl>
    ##  1 rs62513865     8 101592213 0.809   43986.
    ##  2 rs79643588     8 106973048 0.461   46211.
    ##  3 rs17396518     8 108690829 0.00465 45747.
    ##  4 rs983166       8 108681675 0.145   46257.
    ##  5 rs28842593     8 103044620 0.842   39722.
    ##  6 rs7014597      8 104152280 0.293   46026.
    ##  7 rs3134156      8 100479917 0.353   46257.
    ##  8 rs6980591      8 103144592 0.00839 46211.
    ##  9 rs72670434     8 108166508 0.390   45655.
    ## 10 rs10955343     8 105201080 0.659   46350 
    ## # ℹ 8,058,787 more rows

## Output

``` r
write_tsv(common_MAGMA, "../data/processed/common_MAGMA.txt")
```

## Session Info

``` r
sessionInfo()
```

    ## R version 4.4.3 (2025-02-28)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: AlmaLinux 8.10 (Cerulean Leopard)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: FlexiBLAS NETLIB;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: America/New_York
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4    
    ##  [5] purrr_1.0.4     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1   
    ##  [9] ggplot2_3.5.1   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bit_4.6.0         gtable_0.3.6      crayon_1.5.3      compiler_4.4.3   
    ##  [5] tidyselect_1.2.1  parallel_4.4.3    scales_1.3.0      yaml_2.3.11      
    ##  [9] fastmap_1.2.0     R6_2.6.1          generics_0.1.4    knitr_1.49       
    ## [13] munsell_0.5.1     pillar_1.10.1     tzdb_0.4.0        rlang_1.1.6      
    ## [17] utf8_1.2.6        stringi_1.8.7     xfun_0.54         bit64_4.6.0-1    
    ## [21] timechange_0.3.0  cli_3.6.5         withr_3.0.2       magrittr_2.0.4   
    ## [25] digest_0.6.39     grid_4.4.3        vroom_1.6.5       rstudioapi_0.17.1
    ## [29] hms_1.1.3         lifecycle_1.0.4   vctrs_0.6.5       evaluate_1.0.5   
    ## [33] glue_1.8.0        colorspace_2.1-2  rmarkdown_2.29    tools_4.4.3      
    ## [37] pkgconfig_2.0.3   htmltools_0.5.8.1
