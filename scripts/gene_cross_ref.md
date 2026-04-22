SFARI Gene List Cross-Reference
================
Amalya Murrill
2026-04-22

The SFARI database collects information about genes empirically
implicated in ASD and categorizes them, specifically by how
well-verified their mutations are in ASD risk. Our annotated genes were
cross-referenced against this SFARI gene list, identifying 33% as
SFARI-annotated.

``` r
knitr::opts_chunk$set(warning = FALSE, message = FALSE)
```

## Libraries

``` r
library(tidyverse)
```

## Import SFARI Genes

``` r
#import SFARI gene database data
sfari_data <- read.csv("../../data/raw/SFARI_genes.csv", header = TRUE)
colnames(sfari_data)
```

    ##  [1] "status"            "gene.symbol"       "gene.name"        
    ##  [4] "ensembl.id"        "chromosome"        "genetic.category" 
    ##  [7] "gene.score"        "syndromic"         "eagle"            
    ## [10] "number.of.reports"

``` r
sfari_data <- sfari_data %>% dplyr::select(gene.symbol, gene.name, genetic.category, gene.score)
head(sfari_data)
```

    ##   gene.symbol                                            gene.name
    ## 1        ABAT                     4-aminobutyrate aminotransferase
    ## 2      ABCA10 ATP-binding cassette, sub-family A (ABC1), member 10
    ## 3      ABCA13           ATP binding cassette subfamily A member 13
    ## 4       ABCA2            ATP binding cassette subfamily A member 2
    ## 5       ABCA7  ATP-binding cassette, sub-family A (ABC1), member 7
    ## 6       ABCE1            ATP binding cassette subfamily E member 1
    ##                                 genetic.category gene.score
    ## 1 Rare Single Gene Mutation, Genetic Association          2
    ## 2                      Rare Single Gene Mutation          2
    ## 3          Rare Single Gene Mutation, Functional          2
    ## 4                      Rare Single Gene Mutation          3
    ## 5                      Rare Single Gene Mutation          2
    ## 6                      Rare Single Gene Mutation          1

## Import Annotated Pathways

``` r
anno_results <- read.csv("../../results/loeuf_anno_results.csv", header = TRUE)
head(anno_results) 
```

    ##   variant                                                   pathway       padj
    ## 1    rare                                         GOBP_NEURON_DEATH 0.08450704
    ## 2    rare                  GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH 0.08450704
    ## 3    rare GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT 0.08450704
    ## 4    rare        GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION 0.08450704
    ## 5    rare                  GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS 0.08450704
    ## 6    rare                  GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS 0.09208973
    ##                                                                                                                                                                                                  leadingEdge
    ## 1                                                                                                             SCN2A;SYNGAP1;ADNP;GRIN2B;CTNNB1;STXBP1;GABRB2;CORO1A;EPHB1;NF1;THRB;MAP2K7;HSP90AB1;GSK3B;UBB
    ## 2                                                                                                                                        GRIN2B;CTNNB1;NF1;MAP2K7;GSK3B;SRPK2;FCGR2B;ATF4;FOS;GRN;GRIK5;PRNP
    ## 3                                                                                                                                                                                   SYNGAP1;PTEN;BCL11A;GFAP
    ## 4                                                        PTEN;SHANK3;SIN3A;TBR1;LDB1;EPHB1;HSP90AB1;UBB;SATB2;SUFU;GABRB1;NFIB;EPHB2;NR4A2;HES5;HSP90AA1;FEZF2;ISL2;RAPGEF2;RAC3;RORA;NDNF;ARHGAP35;DCC;MAP2
    ## 5                                                                                                                                                   ADNP;DSCAM;RUFY3;PLXNA2;PLXNA1;MAP2K1;STK11;LIMK1;PLXNB1
    ## 6 ADNP;SHANK3;CTNNB1;DSCAM;GFAP;RUFY3;RELN;PLXNA2;ITPKA;CUX2;EPHB2;PPP1CC;PTPRD;PLXNA1;CAPRIN1;MAP2K1;STK11;SOX10;MTOR;LIMK1;PLXNB1;DICER1;ACTR2;VEGFC;XRCC5;HIF1A;PTPRZ1;KIT;TIAM1;SEMA4D;NPTN;PLXNB2;NDEL1
    ##   converge_diverge percent_constrained
    ## 1        divergent            86.66667
    ## 2       convergent            58.33333
    ## 3        divergent            50.00000
    ## 4       convergent            68.00000
    ## 5       convergent            55.55556
    ## 6       convergent            78.78788

## Join SFARI Data

``` r
#prepare annotated pathway data
anno_paths <- anno_results %>% 
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
  mutate(gene = sapply(leadingEdge, str_split ,";"))%>% 
  dplyr::select(-leadingEdge) %>% 
  unnest(gene) %>%
  relocate(c(gene, pathway), .before = 1)

nrow(anno_paths) # 1422 individual gene-pathway groupings
```

    ## [1] 1422

``` r
sfari_data_rename <- sfari_data %>%
  dplyr::rename(gene = gene.symbol)

cross_ref_inner <- anno_paths %>%
  dplyr::inner_join(sfari_data_rename, by = "gene")
cross_ref_inner
```

    ## # A tibble: 470 × 9
    ##    gene    pathway variant   padj converge_diverge percent_constrained gene.name
    ##    <chr>   <chr>   <chr>    <dbl> <chr>                          <dbl> <chr>    
    ##  1 SCN2A   GOBP_N… rare    0.0845 divergent                       86.7 sodium c…
    ##  2 SYNGAP1 GOBP_N… rare    0.0845 divergent                       86.7 synaptic…
    ##  3 ADNP    GOBP_N… rare    0.0845 divergent                       86.7 Activity…
    ##  4 GRIN2B  GOBP_N… rare    0.0845 divergent                       86.7 glutamat…
    ##  5 CTNNB1  GOBP_N… rare    0.0845 divergent                       86.7 catenin …
    ##  6 STXBP1  GOBP_N… rare    0.0845 divergent                       86.7 Syntaxin…
    ##  7 GABRB2  GOBP_N… rare    0.0845 divergent                       86.7 gamma-am…
    ##  8 CORO1A  GOBP_N… rare    0.0845 divergent                       86.7 coronin …
    ##  9 EPHB1   GOBP_N… rare    0.0845 divergent                       86.7 EPH rece…
    ## 10 NF1     GOBP_N… rare    0.0845 divergent                       86.7 neurofib…
    ## # ℹ 460 more rows
    ## # ℹ 2 more variables: genetic.category <chr>, gene.score <int>

``` r
nrow(cross_ref_inner) # 470 individual gene-pathway groupings, 33% match with SFARI
```

    ## [1] 470

``` r
cross_ref_outer <- anno_paths %>%
  dplyr::left_join(sfari_data_rename, by = "gene")
cross_ref_outer %>% filter(is.na(gene.name)) # perhaps genes of future interest?
```

    ## # A tibble: 952 × 9
    ##    gene    pathway variant   padj converge_diverge percent_constrained gene.name
    ##    <chr>   <chr>   <chr>    <dbl> <chr>                          <dbl> <chr>    
    ##  1 THRB    GOBP_N… rare    0.0845 divergent                       86.7 <NA>     
    ##  2 MAP2K7  GOBP_N… rare    0.0845 divergent                       86.7 <NA>     
    ##  3 HSP90A… GOBP_N… rare    0.0845 divergent                       86.7 <NA>     
    ##  4 GSK3B   GOBP_N… rare    0.0845 divergent                       86.7 <NA>     
    ##  5 UBB     GOBP_N… rare    0.0845 divergent                       86.7 <NA>     
    ##  6 MAP2K7  GOBP_P… rare    0.0845 convergent                      58.3 <NA>     
    ##  7 GSK3B   GOBP_P… rare    0.0845 convergent                      58.3 <NA>     
    ##  8 SRPK2   GOBP_P… rare    0.0845 convergent                      58.3 <NA>     
    ##  9 FCGR2B  GOBP_P… rare    0.0845 convergent                      58.3 <NA>     
    ## 10 ATF4    GOBP_P… rare    0.0845 convergent                      58.3 <NA>     
    ## # ℹ 942 more rows
    ## # ℹ 2 more variables: genetic.category <chr>, gene.score <int>

## Explore Results

``` r
#top 5 by padj
cross_ref_inner %>% arrange(padj) %>% distinct(pathway, .keep_all = TRUE) %>% head(5)
```

    ## # A tibble: 5 × 9
    ##   gene    pathway  variant   padj converge_diverge percent_constrained gene.name
    ##   <chr>   <chr>    <chr>    <dbl> <chr>                          <dbl> <chr>    
    ## 1 HDAC4   GOBP_CH… common  0.0212 divergent                       53.8 histone …
    ## 2 PLXNA4  GOBP_NE… common  0.0672 convergent                      30   Plexin A4
    ## 3 SCN2A   GOBP_NE… rare    0.0845 divergent                       86.7 sodium c…
    ## 4 GRIN2B  GOBP_PO… rare    0.0845 convergent                      58.3 glutamat…
    ## 5 SYNGAP1 GOBP_NE… rare    0.0845 divergent                       50   synaptic…
    ## # ℹ 2 more variables: genetic.category <chr>, gene.score <int>

``` r
#top 5 by percent constrained
cross_ref_inner %>% arrange(desc(percent_constrained)) %>% distinct(pathway, .keep_all = TRUE) %>% head(5)
```

    ## # A tibble: 5 × 9
    ##   gene   pathway   variant   padj converge_diverge percent_constrained gene.name
    ##   <chr>  <chr>     <chr>    <dbl> <chr>                          <dbl> <chr>    
    ## 1 SCN2A  GOBP_NEU… rare    0.0845 divergent                       86.7 sodium c…
    ## 2 CHD8   GOBP_CHR… rare    0.0967 divergent                       85.7 chromodo…
    ## 3 DSCAM  GOBP_NEU… rare    0.106  convergent                      85.7 Down syn…
    ## 4 PLXNA4 GOBP_CEN… common  0.144  convergent                      80   Plexin A4
    ## 5 ADNP   GOBP_POS… rare    0.0921 convergent                      78.8 Activity…
    ## # ℹ 2 more variables: genetic.category <chr>, gene.score <int>

``` r
#how often genes appear (hub genes)
cross_ref_inner %>% count(gene, sort = TRUE)
```

    ## # A tibble: 183 × 2
    ##    gene       n
    ##    <chr>  <int>
    ##  1 CTNNB1    12
    ##  2 EPHB2     11
    ##  3 EPHB1     10
    ##  4 NF1       10
    ##  5 PTEN      10
    ##  6 SYNJ1     10
    ##  7 PPFIA3     9
    ##  8 UNC13A     9
    ##  9 ADCY1      8
    ## 10 CADPS      8
    ## # ℹ 173 more rows
