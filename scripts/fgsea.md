Rare Variant Preparation and FGSEA
================
Amalya Murrill
2026-04-19

This script performs fgsea pathway enrichment analysis on both the rare
variant (Satterstrom et al. 2020 TADA+ q-values) and common variant
(MAGMA Z-statistics) arms. A harmonization problem exists between the
two summary statistics; as a provisional measure, rare variant genes
below the first quartile of -log10(q-values) were retained to reduce
zero-inflation. Gene IDs are standardized to Entrez, ranked, and fed
into fgsea using ASD-relevant GO-BP and KEGG pathway sets. Convergent
and divergent pathways are identified through qualitative overlap and a
MAD-based quantitative classification on NES differences. An FDR
threshold of 0.25 serves as a placeholder pending formal harmonization.
Classified pathways are exported for downstream gene-level and
constraint analysis.

``` r
knitr::opts_chunk$set(warning = FALSE, message = FALSE)
```

## Libraries

``` r
library(tidyverse)
library(ggplot2)
library(org.Hs.eg.db)
library(fgsea)
library(msigdbr)
library(patchwork)
```

## Data Import: Rare Variants

``` r
rare <- read_tsv("../data/raw/ASC_gene_results.tsv.bgz", col_names= TRUE)
#quick check
head(rare) # two junk gene names
```

    ## # A tibble: 6 × 15
    ##   gene_id         group xcase_dn_ptv xcont_dn_ptv xcase_dn_misb xcont_dn_misb
    ##   <chr>           <chr>        <dbl>        <dbl>         <dbl>         <dbl>
    ## 1 .               All              0            0             0             0
    ## 2 0               All              0            0             0             0
    ## 3 ENSG00000000419 All              0            0             0             0
    ## 4 ENSG00000000457 All              0            0             0             0
    ## 5 ENSG00000000460 All              0            0             0             0
    ## 6 ENSG00000000938 All              0            0             0             0
    ## # ℹ 9 more variables: xcase_dn_misa <dbl>, xcont_dn_misa <dbl>,
    ## #   xcase_dbs_ptv <dbl>, xcont_dbs_ptv <dbl>, xcase_swe_ptv <dbl>,
    ## #   xcont_swe_ptv <dbl>, xcase_tut <dbl>, xcont_tut <dbl>, qval <dbl>

``` r
colnames(rare)
```

    ##  [1] "gene_id"       "group"         "xcase_dn_ptv"  "xcont_dn_ptv" 
    ##  [5] "xcase_dn_misb" "xcont_dn_misb" "xcase_dn_misa" "xcont_dn_misa"
    ##  [9] "xcase_dbs_ptv" "xcont_dbs_ptv" "xcase_swe_ptv" "xcont_swe_ptv"
    ## [13] "xcase_tut"     "xcont_tut"     "qval"

``` r
nrow(rare) # 17247 genes
```

    ## [1] 17347

``` r
#remove junk gene names
rare <- rare %>% filter(startsWith(gene_id, "ENSG"))
head(rare)
```

    ## # A tibble: 6 × 15
    ##   gene_id         group xcase_dn_ptv xcont_dn_ptv xcase_dn_misb xcont_dn_misb
    ##   <chr>           <chr>        <dbl>        <dbl>         <dbl>         <dbl>
    ## 1 ENSG00000000419 All              0            0             0             0
    ## 2 ENSG00000000457 All              0            0             0             0
    ## 3 ENSG00000000460 All              0            0             0             0
    ## 4 ENSG00000000938 All              0            0             0             0
    ## 5 ENSG00000000971 All              0            0             0             0
    ## 6 ENSG00000001036 All              0            0             0             0
    ## # ℹ 9 more variables: xcase_dn_misa <dbl>, xcont_dn_misa <dbl>,
    ## #   xcase_dbs_ptv <dbl>, xcont_dbs_ptv <dbl>, xcase_swe_ptv <dbl>,
    ## #   xcont_swe_ptv <dbl>, xcase_tut <dbl>, xcont_tut <dbl>, qval <dbl>

``` r
nrow(rare) # 17345, 2 bad genes removed
```

    ## [1] 17345

## Data Import: Common Variant MAGMA Output

``` r
magma_common <- read_table("../data/processed/common_MAGMA_genes.genes.out")
#quick check
head(magma_common)
```

    ## # A tibble: 6 × 9
    ##     GENE   CHR  START   STOP NSNPS NPARAM     N    ZSTAT      P
    ##    <dbl> <dbl>  <dbl>  <dbl> <dbl>  <dbl> <dbl>    <dbl>  <dbl>
    ## 1 148398     1 859993 879961     5      1 40195 -0.480   0.684 
    ## 2  26155     1 879583 894679    30      2 39915  0.00248 0.499 
    ## 3 339451     1 895967 901099    18      2 38249 -0.0244  0.510 
    ## 4  84069     1 901872 910488     6      3 34616  0.925   0.178 
    ## 5  84808     1 910579 917473    18      4 38898  2.10    0.0180
    ## 6  57801     1 934342 936608    11      2 39671  0.835   0.202

``` r
colnames(magma_common)
```

    ## [1] "GENE"   "CHR"    "START"  "STOP"   "NSNPS"  "NPARAM" "N"      "ZSTAT" 
    ## [9] "P"

``` r
nrow(magma_common) # 18072 genes
```

    ## [1] 18072

## Visualize Distributions

``` r
#plot the rare
rare_dist <- rare %>% ggplot(aes(x = qval)) +
  geom_histogram(fill = "#F4D03F",alpha = 0.9, bins = 30) +
  xlab("TADA+ q-value") +
  theme_classic()

rare_dist
```

![](fgsea_files/figure-gfm/stat-distributions-1.png)<!-- -->

``` r
#plot the common
common_dist <- magma_common %>% ggplot(aes(x = ZSTAT)) +
  geom_histogram(fill = "lightblue", alpha = 0.9, bins = 50) +
  xlab("MAGMA Z Statistic") +
  theme_classic()

common_dist
```

![](fgsea_files/figure-gfm/stat-distributions-2.png)<!-- -->

``` r
#put them together
distributions <- common_dist + rare_dist

distributions
```

![](fgsea_files/figure-gfm/stat-distributions-3.png)<!-- -->

``` r
#plot -log10rare
rare_log_quants <- rare %>% dplyr::select(qval) %>% pull() %>% quantile()

rare_dist_log <- rare %>% ggplot(aes(x = -log10(qval))) +
  geom_histogram(fill = "#F4D03F",alpha = 0.9, bins = 80) +
  geom_vline(xintercept = rare_log_quants, color = "grey", linetype = "dashed" ) +
  xlim(0.5, 2.5) +
  xlab("TADA+ -log10(q-value)") +
  theme_classic()

rare_dist_log # remove everything before first quartile
```

![](fgsea_files/figure-gfm/stat-distributions-4.png)<!-- -->

``` r
rare_log_quants
```

    ##      0%     25%     50%     75%    100% 
    ## 0.00000 0.87522 0.91261 0.92517 0.93140

``` r
#-log10(q) = 0.87522, so 
qval_thresh <- 10^(-0.87522) # 0.13
```

### Chi-Squared Goodness of Fit

``` r
#use chi-square GOF to see threshold of deviation from uniform distribution in this range 
#qval_start <- rare %>% arrange(qval) %>% dplyr::select(qval) %>% pull()
#test <- list(p.value = 0)
#threshold <- 0.05

#while(test$p.value < 0.05 & threshold < 0.5){
#threshold <- threshold + 0.01
#qval_current <- qval_start
#keep <- qval_current <= threshold
#qval_current <- qval_current[keep]
#binned_qvals<- cut(qval_current, breaks = 50)
#observed <- table(binned_qvals) 
#test <- chisq.test(observed)
#if(test$p.value > 0.05){
  #print(paste0("threshold =",threshold))
#}
#}
```

## Apply Data-Driven Filter to Rare Arm

``` r
rare_filt <- rare %>% filter(qval < qval_thresh)
#quick check
nrow(rare_filt) # from 17247 genes to 118
```

    ## [1] 118

``` r
head(rare_filt)
```

    ## # A tibble: 6 × 15
    ##   gene_id         group xcase_dn_ptv xcont_dn_ptv xcase_dn_misb xcont_dn_misb
    ##   <chr>           <chr>        <dbl>        <dbl>         <dbl>         <dbl>
    ## 1 ENSG00000005339 All              0            0             1             0
    ## 2 ENSG00000005483 All              2            0             0             0
    ## 3 ENSG00000011485 All              0            0             0             0
    ## 4 ENSG00000021574 All              2            0             0             0
    ## 5 ENSG00000038382 All              1            0             1             0
    ## 6 ENSG00000042753 All              0            0             2             0
    ## # ℹ 9 more variables: xcase_dn_misa <dbl>, xcont_dn_misa <dbl>,
    ## #   xcase_dbs_ptv <dbl>, xcont_dbs_ptv <dbl>, xcase_swe_ptv <dbl>,
    ## #   xcont_swe_ptv <dbl>, xcase_tut <dbl>, xcont_tut <dbl>, qval <dbl>

## Make Rare Ranked Vector

``` r
rare_tib <- rare %>% mutate(neglog10q = -log10(qval)) %>% dplyr::select(gene_id, neglog10q) %>% arrange(desc(neglog10q))
rare_tib 
```

    ## # A tibble: 17,345 × 2
    ##    gene_id         neglog10q
    ##    <chr>               <dbl>
    ##  1 ENSG00000100888    Inf   
    ##  2 ENSG00000136531    Inf   
    ##  3 ENSG00000197283    Inf   
    ##  4 ENSG00000101126     14.1 
    ##  5 ENSG00000114861     11.8 
    ##  6 ENSG00000143442      9.96
    ##  7 ENSG00000049618      9.59
    ##  8 ENSG00000110066      9.26
    ##  9 ENSG00000157540      9.09
    ## 10 ENSG00000157103      8.72
    ## # ℹ 17,335 more rows

``` r
# look at data to check for NAs or Inf before making ranked vector
sum(is.na(rare_tib)) # 0 NA values
```

    ## [1] 0

``` r
sum(is.infinite(rare_tib$neglog10q)) # 3 Inf values
```

    ## [1] 3

``` r
#go back and correct for Inf values using qval floor of 1e-300 (winsorize)
rare_tib_adj <- rare %>% mutate(neglog10q = -log10(pmax(qval, 1e-300))) %>% dplyr::select(gene_id, neglog10q) %>% arrange(desc(neglog10q))

head(rare_tib_adj)
```

    ## # A tibble: 6 × 2
    ##   gene_id         neglog10q
    ##   <chr>               <dbl>
    ## 1 ENSG00000100888    300   
    ## 2 ENSG00000136531    300   
    ## 3 ENSG00000197283    300   
    ## 4 ENSG00000101126     14.1 
    ## 5 ENSG00000114861     11.8 
    ## 6 ENSG00000143442      9.96

``` r
#now rank for fgsea
rare_ranked <- deframe(rare_tib_adj)
head(rare_ranked)
```

    ## ENSG00000100888 ENSG00000136531 ENSG00000197283 ENSG00000101126 ENSG00000114861 
    ##      300.000000      300.000000      300.000000       14.069509       11.751487 
    ## ENSG00000143442 
    ##        9.961499

``` r
#convert to Entrez IDs for fgsea
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = names(rare_ranked),
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

#how many are NA?
sum(is.na(entrez_ids)) # 110
```

    ## [1] 108

``` r
#subset by non-NA values
keep_idx <- !is.na(entrez_ids)
entrez_ids <- entrez_ids[keep_idx]
rare_ranked <- rare_ranked[keep_idx]
#rename
names(rare_ranked) <- entrez_ids

#deduplicate
rare_ranked_df <- enframe(rare_ranked, name = "gene_id", value = "neglog10q") %>% group_by(gene_id) %>% slice_max(neglog10q, n = 1, with_ties = FALSE) %>% ungroup()

rare_ranked <-deframe(rare_ranked_df)
  
head(rare_ranked)
```

    ##          1         10        100       1000      10000      10001 
    ## 0.04586366 0.09229312 0.04716398 0.03637742 0.03789122 0.03342184

## Make Common Ranked Vector

``` r
common_tib <- magma_common %>% dplyr::select(GENE, ZSTAT) %>% arrange(desc(ZSTAT)) #use ZSTAT to preserve direction
head(common_tib)
```

    ## # A tibble: 6 × 2
    ##     GENE ZSTAT
    ##    <dbl> <dbl>
    ## 1  22803  5.53
    ## 2  55857  5.27
    ## 3   3781  5.15
    ## 4 140733  4.82
    ## 5   2550  4.57
    ## 6   4137  4.51

``` r
# look at data to check for NAs or Inf before making ranked vector
sum(is.na(common_tib)) # 0 NA values
```

    ## [1] 0

``` r
sum(ZSTAT = 0) # no 0 ZSTAT values
```

    ## [1] 0

``` r
#now rank for fgsea
common_ranked <- deframe(common_tib)
head(common_ranked)
```

    ##  22803  55857   3781 140733   2550   4137 
    ## 5.5322 5.2678 5.1470 4.8219 4.5739 4.5124

## Pathway Selection

``` r
go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")

head(go_bp)
```

    ## # A tibble: 6 × 15
    ##   gs_cat gs_subcat gs_name                  gene_symbol entrez_gene ensembl_gene
    ##   <chr>  <chr>     <chr>                    <chr>             <int> <chr>       
    ## 1 C5     GO:BP     GOBP_10_FORMYLTETRAHYDR… AASDHPPT          60496 ENSG0000014…
    ## 2 C5     GO:BP     GOBP_10_FORMYLTETRAHYDR… ALDH1L1           10840 ENSG0000014…
    ## 3 C5     GO:BP     GOBP_10_FORMYLTETRAHYDR… ALDH1L2          160428 ENSG0000013…
    ## 4 C5     GO:BP     GOBP_10_FORMYLTETRAHYDR… MTHFD1             4522 ENSG0000010…
    ## 5 C5     GO:BP     GOBP_10_FORMYLTETRAHYDR… MTHFD1L           25902 ENSG0000012…
    ## 6 C5     GO:BP     GOBP_10_FORMYLTETRAHYDR… MTHFD2L          441024 ENSG0000016…
    ## # ℹ 9 more variables: human_gene_symbol <chr>, human_entrez_gene <int>,
    ## #   human_ensembl_gene <chr>, gs_id <chr>, gs_pmid <chr>, gs_geoid <chr>,
    ## #   gs_exact_source <chr>, gs_url <chr>, gs_description <chr>

``` r
head(kegg)
```

    ## # A tibble: 6 × 15
    ##   gs_cat gs_subcat gs_name               gene_symbol entrez_gene ensembl_gene   
    ##   <chr>  <chr>     <chr>                 <chr>             <int> <chr>          
    ## 1 C2     CP:KEGG   KEGG_ABC_TRANSPORTERS ABCA1                19 ENSG00000165029
    ## 2 C2     CP:KEGG   KEGG_ABC_TRANSPORTERS ABCA10            10349 ENSG00000154263
    ## 3 C2     CP:KEGG   KEGG_ABC_TRANSPORTERS ABCA12            26154 ENSG00000144452
    ## 4 C2     CP:KEGG   KEGG_ABC_TRANSPORTERS ABCA13           154664 ENSG00000179869
    ## 5 C2     CP:KEGG   KEGG_ABC_TRANSPORTERS ABCA2                20 ENSG00000107331
    ## 6 C2     CP:KEGG   KEGG_ABC_TRANSPORTERS ABCA3                21 ENSG00000167972
    ## # ℹ 9 more variables: human_gene_symbol <chr>, human_entrez_gene <int>,
    ## #   human_ensembl_gene <chr>, gs_id <chr>, gs_pmid <chr>, gs_geoid <chr>,
    ## #   gs_exact_source <chr>, gs_url <chr>, gs_description <chr>

``` r
# filter pathways by canonical terms: SYNAP, GABA, EXCIT, INHIB, CHROMATIN, NEURO, AXON, GLUTAMAT, MTOR (metabolism)

name_filt_pattern <- "SYNAP|GABA|EXCIT|INHIB|CHROMATIN|NEURO|AXON|GLUTAMAT|MTOR"
go_bp_filt <- go_bp %>% filter(grepl(name_filt_pattern, gs_name))
head(go_bp_filt)
```

    ## # A tibble: 6 × 15
    ##   gs_cat gs_subcat gs_name                  gene_symbol entrez_gene ensembl_gene
    ##   <chr>  <chr>     <chr>                    <chr>             <int> <chr>       
    ## 1 C5     GO:BP     GOBP_ADENYLATE_CYCLASE_… ADCY5               111 ENSG0000017…
    ## 2 C5     GO:BP     GOBP_ADENYLATE_CYCLASE_… DRD2               1813 ENSG0000014…
    ## 3 C5     GO:BP     GOBP_ADENYLATE_CYCLASE_… DRD3               1814 ENSG0000015…
    ## 4 C5     GO:BP     GOBP_ADENYLATE_CYCLASE_… DRD4               1815 ENSG0000006…
    ## 5 C5     GO:BP     GOBP_ADENYLATE_CYCLASE_… DRD4               1815 ENSG0000027…
    ## 6 C5     GO:BP     GOBP_ADENYLATE_CYCLASE_… FLNA               2316 ENSG0000019…
    ## # ℹ 9 more variables: human_gene_symbol <chr>, human_entrez_gene <int>,
    ## #   human_ensembl_gene <chr>, gs_id <chr>, gs_pmid <chr>, gs_geoid <chr>,
    ## #   gs_exact_source <chr>, gs_url <chr>, gs_description <chr>

``` r
kegg_filt <- kegg %>% filter(grepl(name_filt_pattern, gs_name))
head (kegg_filt)
```

    ## # A tibble: 6 × 15
    ##   gs_cat gs_subcat gs_name                  gene_symbol entrez_gene ensembl_gene
    ##   <chr>  <chr>     <chr>                    <chr>             <int> <chr>       
    ## 1 C2     CP:KEGG   KEGG_ALANINE_ASPARTATE_… ABAT                 18 ENSG0000018…
    ## 2 C2     CP:KEGG   KEGG_ALANINE_ASPARTATE_… ACY3              91703 ENSG0000013…
    ## 3 C2     CP:KEGG   KEGG_ALANINE_ASPARTATE_… ADSL                158 ENSG0000023…
    ## 4 C2     CP:KEGG   KEGG_ALANINE_ASPARTATE_… ADSS1            122622 ENSG0000018…
    ## 5 C2     CP:KEGG   KEGG_ALANINE_ASPARTATE_… ADSS2               159 ENSG0000003…
    ## 6 C2     CP:KEGG   KEGG_ALANINE_ASPARTATE_… AGXT                189 ENSG0000017…
    ## # ℹ 9 more variables: human_gene_symbol <chr>, human_entrez_gene <int>,
    ## #   human_ensembl_gene <chr>, gs_id <chr>, gs_pmid <chr>, gs_geoid <chr>,
    ## #   gs_exact_source <chr>, gs_url <chr>, gs_description <chr>

## List Format Conversion

``` r
# merge two dataframes and make list
pathways <- rbind(go_bp_filt, kegg_filt) %>% split(x = as.character(.$entrez_gene), f = .$gs_name)
head(pathways)
```

    ## $GOBP_ADENYLATE_CYCLASE_INHIBITING_DOPAMINE_RECEPTOR_SIGNALING_PATHWAY
    ## [1] "111"   "1813"  "1814"  "1815"  "1815"  "2316"  "10419"
    ## 
    ## $GOBP_ADENYLATE_CYCLASE_INHIBITING_G_PROTEIN_COUPLED_ACETYLCHOLINE_RECEPTOR_SIGNALING_PATHWAY
    ## [1] "1128"  "1129"  "1131"  "1132"  "1133"  "11255" "59340" "4988" 
    ## 
    ## $GOBP_ADENYLATE_CYCLASE_INHIBITING_G_PROTEIN_COUPLED_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY
    ## [1] "2899" "2911" "2912" "2913" "2914" "2915" "2916" "2917" "2918"
    ## 
    ## $GOBP_ADENYLATE_CYCLASE_INHIBITING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY
    ##  [1] "111"    "112"    "134"    "150"    "9590"   "9495"   "846"    "1128"  
    ##  [9] "1129"   "1131"   "1132"   "1133"   "1325"   "1813"   "1814"   "1815"  
    ## [17] "1815"   "1906"   "1909"   "2865"   "2316"   "2358"   "2550"   "2550"  
    ## [25] "2550"   "2550"   "2550"   "2550"   "2550"   "9568"   "2770"   "2771"  
    ## [33] "2773"   "2779"   "2781"   "11245"  "2861"   "9283"   "2899"   "2911"  
    ## [41] "2912"   "2913"   "2914"   "2915"   "2916"   "2917"   "2918"   "11255" 
    ## [49] "59340"  "3350"   "3351"   "3352"   "3354"   "3355"   "3360"   "3361"  
    ## [57] "3640"   "1902"   "8685"   "2847"   "406922" "407034" "4543"   "4883"  
    ## [65] "4886"   "4887"   "4985"   "4986"   "4987"   "4987"   "4988"   "165140"
    ## [73] "5028"   "64805"  "5064"   "5138"   "10419"  "5660"   "768239" "11251" 
    ## [81] "5996"   "5997"   "60626"  "1901"   "1903"   "6752"  
    ## 
    ## $GOBP_ADENYLATE_CYCLASE_INHIBITING_SEROTONIN_RECEPTOR_SIGNALING_PATHWAY
    ## [1] "3350" "3351" "3352" "3354" "3355" "3360" "3361"
    ## 
    ## $GOBP_AMINO_ACID_NEUROTRANSMITTER_REUPTAKE
    ##  [1] "477"  "2055" "2055" "3688" "3766" "8864" "6529" "6538" "6539" "6540"

## FGSEA

``` r
set.seed(42) #set seed for reproducibility
rare_fgseaRes <- fgsea(
  pathways = pathways,
  stats = rare_ranked,
  minSize = 15,
  maxSize = 500,
  scoreType = c("pos") # from TADA+ qvals, all positive
)

nrow(rare_fgseaRes[order(pval), ]) # 182 
```

    ## [1] 182

``` r
rare_pathways <- rare_fgseaRes[order(pval), ]
```

### Collapse Rare Variant Pathways

``` r
rare_collapsed <- fgsea::collapsePathways(
  fgseaRes = rare_fgseaRes,
  pathways = pathways,
  stats = rare_ranked
)
rare_collapsed$mainPathways
```

    ##  [1] "GOBP_ENSHEATHMENT_OF_NEURONS"                     
    ##  [2] "GOBP_NEGATIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION"
    ##  [3] "GOBP_NEUROMUSCULAR_PROCESS"                       
    ##  [4] "GOBP_NEURONAL_ACTION_POTENTIAL"                   
    ##  [5] "GOBP_NEURON_APOPTOTIC_PROCESS"                    
    ##  [6] "GOBP_NEURON_PROJECTION_EXTENSION"                 
    ##  [7] "GOBP_NEURON_PROJECTION_ORGANIZATION"              
    ##  [8] "GOBP_POSTSYNAPSE_ORGANIZATION"                    
    ##  [9] "GOBP_REGULATION_OF_AXONOGENESIS"                  
    ## [10] "GOBP_REGULATION_OF_CHROMATIN_ORGANIZATION"        
    ## [11] "GOBP_REGULATION_OF_SYNAPTIC_PLASTICITY"           
    ## [12] "GOBP_SYNAPSE_ASSEMBLY"

``` r
#write.csv(as.data.frame(rare_collapsed$parentPathways), "../results/rare_collapsed_pathway_mapping.csv") #for summary
```

``` r
set.seed(42) #set seed for reproducibility
common_fgseaRes <- fgsea(
  pathways = pathways,
  stats = common_ranked,
  minSize = 15,
  maxSize = 500
)

nrow(common_fgseaRes[order(pval), ]) # 184
```

    ## [1] 184

``` r
common_pathways <- common_fgseaRes[order(pval), ]
```

### Collapse Common Variant Pathways

``` r
common_collapsed <- fgsea::collapsePathways(
  fgseaRes = common_fgseaRes,
  pathways = pathways,
  stats = common_ranked
)
common_collapsed$mainPathways
```

    ## [1] "GOBP_CHROMATIN_REMODELING"                                              
    ## [2] "GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE"
    ## [3] "GOBP_NEUROTRANSMITTER_SECRETION"                                        
    ## [4] "GOBP_SYNAPSE_ORGANIZATION"                                              
    ## [5] "KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM"                        
    ## [6] "KEGG_MTOR_SIGNALING_PATHWAY"

``` r
#write.csv(as.data.frame(common_collapsed$parentPathways), "../results/common_collapsed_pathway_mapping.csv") # for summary
```

### Visualize Pathway Statistics

``` r
#plot the rare paths
rare_path_padjs <- rare_pathways %>% ggplot(aes(x = padj)) +
  geom_histogram(fill = "lightpink",alpha = 0.9, bins = 30) +
  xlab("Rare Pathway Padj") +
  theme_classic()

#plot the common paths
common_path_padjs <- common_pathways %>% ggplot(aes(x = padj)) +
  geom_histogram(fill = "lightgreen", alpha = 0.9, bins = 50) +
  xlab("Common Pathway Padj") +
  theme_classic()


#put them together
pathway_padjs <- common_path_padjs + rare_path_padjs

pathway_padjs
```

![](fgsea_files/figure-gfm/plot-pathway-padjs-1.png)<!-- -->

## Combine Pathway Findings

``` r
rare_pathways <- rare_pathways %>% mutate(variant = rep("rare", nrow(rare_pathways))) %>% relocate(variant)
common_pathways <- common_pathways %>% mutate(variant = rep("common", nrow(common_pathways))) %>% relocate(variant)
all_paths <- rbind(rare_pathways, common_pathways)

#are common and rare both there?
head(all_paths)
```

    ##    variant                                                       pathway
    ##     <char>                                                        <char>
    ## 1:    rare                         GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE
    ## 2:    rare                                 GOBP_NEURON_APOPTOTIC_PROCESS
    ## 3:    rare                                             GOBP_NEURON_DEATH
    ## 4:    rare GOBP_POSITIVE_REGULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL
    ## 5:    rare                              GOBP_EXCITATORY_SYNAPSE_ASSEMBLY
    ## 6:    rare                    GOBP_POSITIVE_REGULATION_OF_AXON_EXTENSION
    ##            pval       padj   log2err        ES      NES  size  leadingEdge
    ##           <num>      <num>     <num>     <num>    <num> <int>       <list>
    ## 1: 0.0002756957 0.05017661 0.4984931 0.9924235 2.534766    19  8831, 22941
    ## 2: 0.0007189713 0.06542639 0.4772708 0.9695995 2.284397   218 6326, 88....
    ## 3: 0.0023493938 0.08450704 0.4317077 0.9595688 2.212238   325 6326, 88....
    ## 4: 0.0025713793 0.08450704 0.4317077 0.9160108 2.344457    23 5728, 85....
    ## 5: 0.0027933652 0.08450704 0.4317077 0.9278930 2.379915    24 5728, 85....
    ## 6: 0.0042014678 0.08450704 0.4070179 0.9051298 2.323169    31 23394, 1....

``` r
tail(all_paths)
```

    ##    variant
    ##     <char>
    ## 1:  common
    ## 2:  common
    ## 3:  common
    ## 4:  common
    ## 5:  common
    ## 6:  common
    ##                                                                           pathway
    ##                                                                            <char>
    ## 1:                                       GOBP_L_GLUTAMATE_TRANSMEMBRANE_TRANSPORT
    ## 2:                                               GOBP_EXCITATORY_SYNAPSE_ASSEMBLY
    ## 3: GOBP_ADENYLATE_CYCLASE_INHIBITING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY
    ## 4:                      GOBP_POSITIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT
    ## 5:                                  GOBP_REGULATION_OF_SYNAPTIC_VESICLE_RECYCLING
    ## 6:                                 GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_BALANCE
    ##         pval      padj    log2err         ES        NES  size  leadingEdge
    ##        <num>     <num>      <num>      <num>      <num> <int>       <list>
    ## 1: 0.9832168 0.9980822 0.02970646  0.1814376  0.5753444    28 10550, 7....
    ## 2: 0.9857347 0.9980822 0.03056081  0.1747264  0.5405930    25 64101, 7....
    ## 3: 0.9900374 0.9980822 0.02322483  0.1533360  0.5893707    75 2550, 29....
    ## 4: 0.9909808 0.9980822 0.01699711  0.1600589  0.6670175   139 51517, 1....
    ## 5: 0.9926579 0.9980822 0.03163697  0.1610146  0.4673024    19 6622, 60....
    ## 6: 1.0000000 1.0000000 0.08266464 -0.1427050 -0.6048106    49 7401, 64....

``` r
#QC
nrow(all_paths) # 366
```

    ## [1] 366

``` r
sum(is.na(all_paths)) # 0
```

    ## [1] 0

## Qualitatively Overlap Pathways

``` r
common_paths_main <- all_paths %>% filter(variant == "common") %>% pull(pathway)
rare_paths_main <- all_paths %>% filter(variant == "rare")  %>% pull(pathway)
overlap <- common_paths_main %in% rare_paths_main
overlap_paths <- common_paths_main[overlap]
overlap_paths # 182
```

    ##   [1] "GOBP_CHROMATIN_REMODELING"                                                     
    ##   [2] "GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE"       
    ##   [3] "GOBP_SYNAPSE_ORGANIZATION"                                                     
    ##   [4] "GOBP_REGULATION_OF_NEUROGENESIS"                                               
    ##   [5] "GOBP_NEUROTRANSMITTER_SECRETION"                                               
    ##   [6] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION"                            
    ##   [7] "KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM"                               
    ##   [8] "GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY"                                        
    ##   [9] "KEGG_MTOR_SIGNALING_PATHWAY"                                                   
    ##  [10] "GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE"                                    
    ##  [11] "GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH"                                      
    ##  [12] "GOBP_NEURON_PROJECTION_GUIDANCE"                                               
    ##  [13] "GOBP_NEURON_DEATH"                                                             
    ##  [14] "GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS"                                      
    ##  [15] "GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE"          
    ##  [16] "GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS"                    
    ##  [17] "GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS"                                      
    ##  [18] "GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS"                                    
    ##  [19] "GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                     
    ##  [20] "GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS"                                              
    ##  [21] "KEGG_AXON_GUIDANCE"                                                            
    ##  [22] "GOBP_SYNAPTIC_VESICLE_PRIMING"                                                 
    ##  [23] "GOBP_LONG_TERM_SYNAPTIC_DEPRESSION"                                            
    ##  [24] "GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE"                                         
    ##  [25] "GOBP_REGULATION_OF_PROTEIN_LOCALIZATION_TO_SYNAPSE"                            
    ##  [26] "GOBP_NEUROTRANSMITTER_TRANSPORT"                                               
    ##  [27] "GOBP_REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY"                              
    ##  [28] "GOBP_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT"                                 
    ##  [29] "GOBP_FOREBRAIN_GENERATION_OF_NEURONS"                                          
    ##  [30] "GOBP_NEGATIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION"                             
    ##  [31] "GOBP_REGULATION_OF_AXONOGENESIS"                                               
    ##  [32] "GOBP_AXONAL_TRANSPORT"                                                         
    ##  [33] "GOBP_ANTEROGRADE_AXONAL_TRANSPORT"                                             
    ##  [34] "GOBP_RETINAL_GANGLION_CELL_AXON_GUIDANCE"                                      
    ##  [35] "GOBP_GLUTAMATE_SECRETION"                                                      
    ##  [36] "GOBP_NEGATIVE_REGULATION_OF_NEURON_DIFFERENTIATION"                            
    ##  [37] "GOBP_NEURON_PROJECTION_ORGANIZATION"                                           
    ##  [38] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_AXONOGENESIS"                               
    ##  [39] "GOBP_NEURON_FATE_COMMITMENT"                                                   
    ##  [40] "GOBP_SYNAPSE_MATURATION"                                                       
    ##  [41] "GOBP_NEURON_PROJECTION_EXTENSION"                                              
    ##  [42] "GOBP_POSTSYNAPSE_ORGANIZATION"                                                 
    ##  [43] "GOBP_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                              
    ##  [44] "GOBP_REGULATION_OF_GLUTAMATE_SECRETION"                                        
    ##  [45] "GOBP_DNA_REPLICATION_DEPENDENT_CHROMATIN_ORGANIZATION"                         
    ##  [46] "GOBP_RETROGRADE_AXONAL_TRANSPORT"                                              
    ##  [47] "GOBP_NEURON_DEATH_IN_RESPONSE_TO_OXIDATIVE_STRESS"                             
    ##  [48] "GOBP_NEURON_APOPTOTIC_PROCESS"                                                 
    ##  [49] "GOBP_PRESYNAPTIC_ENDOCYTOSIS"                                                  
    ##  [50] "GOBP_AXON_DEVELOPMENT"                                                         
    ##  [51] "GOBP_NEUROBLAST_PROLIFERATION"                                                 
    ##  [52] "GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION"                                    
    ##  [53] "GOBP_AXONAL_TRANSPORT_OF_MITOCHONDRION"                                        
    ##  [54] "GOBP_POSITIVE_REGULATION_OF_NEURON_APOPTOTIC_PROCESS"                          
    ##  [55] "GOBP_NEGATIVE_REGULATION_OF_AXONOGENESIS"                                      
    ##  [56] "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION"                                  
    ##  [57] "GOBP_FOREBRAIN_NEURON_DIFFERENTIATION"                                         
    ##  [58] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DEVELOPMENT"                                
    ##  [59] "GOBP_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_POTENTIAL"                            
    ##  [60] "GOBP_MOTOR_NEURON_AXON_GUIDANCE"                                               
    ##  [61] "GOBP_AXON_EXTENSION"                                                           
    ##  [62] "GOBP_PROTEIN_LOCALIZATION_TO_SYNAPSE"                                          
    ##  [63] "GOBP_SYNAPSE_ASSEMBLY"                                                         
    ##  [64] "GOBP_REGULATION_OF_NEURONAL_SYNAPTIC_PLASTICITY"                               
    ##  [65] "GOBP_REGULATION_OF_NEURON_DIFFERENTIATION"                                     
    ##  [66] "GOBP_SYNAPTIC_VESICLE_LOCALIZATION"                                            
    ##  [67] "GOBP_POSITIVE_REGULATION_OF_NEUROBLAST_PROLIFERATION"                          
    ##  [68] "GOBP_POSITIVE_REGULATION_OF_NEURON_DIFFERENTIATION"                            
    ##  [69] "GOBP_POSITIVE_REGULATION_OF_CHROMATIN_ORGANIZATION"                            
    ##  [70] "GOBP_SYNAPTIC_MEMBRANE_ADHESION"                                               
    ##  [71] "GOBP_POSITIVE_REGULATION_OF_NEUROTRANSMITTER_SECRETION"                        
    ##  [72] "GOBP_SYNAPTIC_VESICLE_TRANSPORT"                                               
    ##  [73] "GOBP_POSITIVE_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT"                        
    ##  [74] "GOBP_SYNAPTIC_VESICLE_MEMBRANE_ORGANIZATION"                                   
    ##  [75] "GOBP_POSTSYNAPTIC_MEMBRANE_ORGANIZATION"                                       
    ##  [76] "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY"                                           
    ##  [77] "GOBP_NEUROEPITHELIAL_CELL_DIFFERENTIATION"                                     
    ##  [78] "GOBP_REGULATION_OF_NEUROBLAST_PROLIFERATION"                                   
    ##  [79] "GOBP_GLUTAMATE_METABOLIC_PROCESS"                                              
    ##  [80] "GOBP_REGULATION_OF_SYNAPTIC_PLASTICITY"                                        
    ##  [81] "GOBP_POSITIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION"                             
    ##  [82] "GOBP_SPINAL_CORD_MOTOR_NEURON_DIFFERENTIATION"                                 
    ##  [83] "GOBP_NEURONAL_STEM_CELL_POPULATION_MAINTENANCE"                                
    ##  [84] "GOBP_REGULATION_OF_LONG_TERM_NEURONAL_SYNAPTIC_PLASTICITY"                     
    ##  [85] "GOBP_LONG_TERM_SYNAPTIC_POTENTIATION"                                          
    ##  [86] "GOBP_AXONEMAL_DYNEIN_COMPLEX_ASSEMBLY"                                         
    ##  [87] "GOBP_NEUROPEPTIDE_SIGNALING_PATHWAY"                                           
    ##  [88] "GOBP_NEURON_RECOGNITION"                                                       
    ##  [89] "GOBP_SYNAPTIC_VESICLE_RECYCLING"                                               
    ##  [90] "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING"                                   
    ##  [91] "GOBP_POSITIVE_REGULATION_OF_NEUROINFLAMMATORY_RESPONSE"                        
    ##  [92] "GOBP_NEGATIVE_REGULATION_OF_NEURON_DEATH"                                      
    ##  [93] "GOBP_NEURON_MIGRATION"                                                         
    ##  [94] "GOBP_NEUROTROPHIN_TRK_RECEPTOR_SIGNALING_PATHWAY"                              
    ##  [95] "GOBP_ENSHEATHMENT_OF_NEURONS"                                                  
    ##  [96] "GOBP_NEURON_CELLULAR_HOMEOSTASIS"                                              
    ##  [97] "GOBP_RESPONSE_TO_LEUKEMIA_INHIBITORY_FACTOR"                                   
    ##  [98] "GOBP_HETEROCHROMATIN_ORGANIZATION"                                             
    ##  [99] "GOBP_NEUROMUSCULAR_JUNCTION_DEVELOPMENT"                                       
    ## [100] "GOBP_MODULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL"                          
    ## [101] "GOBP_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY"                                     
    ## [102] "GOBP_REGULATION_OF_LONG_TERM_SYNAPTIC_POTENTIATION"                            
    ## [103] "GOBP_REGULATION_OF_SYNAPTIC_VESICLE_EXOCYTOSIS"                                
    ## [104] "GOBP_POSTSYNAPTIC_SPECIALIZATION_ORGANIZATION"                                 
    ## [105] "GOBP_PROTEIN_LOCALIZATION_TO_POSTSYNAPSE"                                      
    ## [106] "GOBP_REGULATION_OF_NEUROTRANSMITTER_UPTAKE"                                    
    ## [107] "GOBP_AXON_ENSHEATHMENT_IN_CENTRAL_NERVOUS_SYSTEM"                              
    ## [108] "GOBP_CEREBRAL_CORTEX_NEURON_DIFFERENTIATION"                                   
    ## [109] "GOBP_SYNAPTIC_TRANSMISSION_CHOLINERGIC"                                        
    ## [110] "GOBP_REGULATION_OF_POSTSYNAPSE_ORGANIZATION"                                   
    ## [111] "GOBP_MODIFICATION_OF_SYNAPTIC_STRUCTURE"                                       
    ## [112] "GOBP_NEURON_FATE_SPECIFICATION"                                                
    ## [113] "GOBP_POSITIVE_REGULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL"                 
    ## [114] "GOBP_REGULATION_OF_RECEPTOR_LOCALIZATION_TO_SYNAPSE"                           
    ## [115] "GOBP_SYNAPTIC_VESICLE_CYTOSKELETAL_TRANSPORT"                                  
    ## [116] "GOBP_G_PROTEIN_COUPLED_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY"                   
    ## [117] "GOBP_DOPAMINERGIC_NEURON_DIFFERENTIATION"                                      
    ## [118] "GOBP_NEUROMUSCULAR_SYNAPTIC_TRANSMISSION"                                      
    ## [119] "GOBP_NEUROTRANSMITTER_REUPTAKE"                                                
    ## [120] "GOBP_REGULATION_OF_SYNAPSE_ASSEMBLY"                                           
    ## [121] "GOBP_DNA_REPLICATION_INDEPENDENT_CHROMATIN_ORGANIZATION"                       
    ## [122] "GOBP_PROTEIN_LOCALIZATION_TO_CHROMATIN"                                        
    ## [123] "GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_POSTURE"                                
    ## [124] "GOBP_NEGATIVE_REGULATION_OF_NEURON_APOPTOTIC_PROCESS"                          
    ## [125] "GOBP_INHIBITORY_SYNAPSE_ASSEMBLY"                                              
    ## [126] "GOBP_AXONAL_FASCICULATION"                                                     
    ## [127] "GOBP_NEURON_PROJECTION_ARBORIZATION"                                           
    ## [128] "GOBP_NEUROINFLAMMATORY_RESPONSE"                                               
    ## [129] "GOBP_PRESYNAPSE_ORGANIZATION"                                                  
    ## [130] "GOBP_NEURON_MATURATION"                                                        
    ## [131] "GOBP_NEUROMUSCULAR_PROCESS"                                                    
    ## [132] "GOBP_NEUROTROPHIN_SIGNALING_PATHWAY"                                           
    ## [133] "GOBP_REGULATION_OF_PRESYNAPSE_ORGANIZATION"                                    
    ## [134] "GOBP_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_NEUROTRANSMITTER_RECEPTOR_LEVELS"     
    ## [135] "GOBP_POSITIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"               
    ## [136] "GOBP_NEURONAL_ACTION_POTENTIAL"                                                
    ## [137] "GOBP_FOREBRAIN_NEURON_DEVELOPMENT"                                             
    ## [138] "GOBP_FACULTATIVE_HETEROCHROMATIN_ASSEMBLY"                                     
    ## [139] "GOBP_CHROMATIN_DISASSEMBLY"                                                    
    ## [140] "GOBP_REGULATION_OF_NEUROTRANSMITTER_RECEPTOR_ACTIVITY"                         
    ## [141] "GOBP_POSTSYNAPTIC_SIGNAL_TRANSDUCTION"                                         
    ## [142] "GOBP_REGULATION_OF_DNA_METHYLATION_DEPENDENT_HETEROCHROMATIN_ASSEMBLY"         
    ## [143] "GOBP_DNA_METHYLATION_DEPENDENT_HETEROCHROMATIN_ASSEMBLY"                       
    ## [144] "GOBP_POSITIVE_REGULATION_OF_LONG_TERM_SYNAPTIC_POTENTIATION"                   
    ## [145] "GOBP_CHEMICAL_SYNAPTIC_TRANSMISSION_POSTSYNAPTIC"                              
    ## [146] "GOBP_NEGATIVE_REGULATION_OF_OXIDATIVE_STRESS_INDUCED_NEURON_DEATH"             
    ## [147] "GOBP_REGULATION_OF_CHROMATIN_ASSEMBLY"                                         
    ## [148] "GOBP_REGULATION_OF_CHROMATIN_BINDING"                                          
    ## [149] "GOBP_REGULATION_OF_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY"                          
    ## [150] "GOBP_MIDBRAIN_DOPAMINERGIC_NEURON_DIFFERENTIATION"                             
    ## [151] "GOBP_AXONEME_ASSEMBLY"                                                         
    ## [152] "GOBP_NEUROTRANSMITTER_METABOLIC_PROCESS"                                       
    ## [153] "GOBP_SYNAPTONEMAL_COMPLEX_ORGANIZATION"                                        
    ## [154] "GOBP_REGULATION_OF_CHROMATIN_ORGANIZATION"                                     
    ## [155] "GOBP_REGULATION_OF_NEURON_PROJECTION_REGENERATION"                             
    ## [156] "GOBP_SYNAPTIC_TRANSMISSION_GABAERGIC"                                          
    ## [157] "GOBP_NEURON_PROJECTION_REGENERATION"                                           
    ## [158] "GOBP_RECEPTOR_LOCALIZATION_TO_SYNAPSE"                                         
    ## [159] "GOBP_SYNAPTIC_TRANSMISSION_DOPAMINERGIC"                                       
    ## [160] "GOBP_NEUROTRANSMITTER_UPTAKE"                                                  
    ## [161] "GOBP_POSITIVE_REGULATION_OF_NEURON_MIGRATION"                                  
    ## [162] "GOBP_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"                                      
    ## [163] "GOBP_RNA_MEDIATED_GENE_SILENCING_BY_INHIBITION_OF_TRANSLATION"                 
    ## [164] "GOBP_RESPONSE_TO_AXON_INJURY"                                                  
    ## [165] "GOBP_POSTSYNAPSE_ASSEMBLY"                                                     
    ## [166] "GOBP_POSITIVE_REGULATION_OF_AXON_EXTENSION"                                    
    ## [167] "GOBP_POSITIVE_REGULATION_OF_SYNAPSE_ASSEMBLY"                                  
    ## [168] "GOBP_CALCIUM_ION_REGULATED_EXOCYTOSIS_OF_NEUROTRANSMITTER"                     
    ## [169] "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GABAERGIC"                            
    ## [170] "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"                        
    ## [171] "GOBP_L_GLUTAMATE_IMPORT_ACROSS_PLASMA_MEMBRANE"                                
    ## [172] "GOBP_REGULATION_OF_NEURON_MIGRATION"                                           
    ## [173] "GOBP_NEUROTRANSMITTER_RECEPTOR_TRANSPORT"                                      
    ## [174] "GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_REGENERATION"                    
    ## [175] "GOBP_NEUROTRANSMITTER_RECEPTOR_INTERNALIZATION"                                
    ## [176] "GOBP_POSTSYNAPTIC_SPECIALIZATION_ASSEMBLY"                                     
    ## [177] "GOBP_L_GLUTAMATE_TRANSMEMBRANE_TRANSPORT"                                      
    ## [178] "GOBP_EXCITATORY_SYNAPSE_ASSEMBLY"                                              
    ## [179] "GOBP_ADENYLATE_CYCLASE_INHIBITING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY"
    ## [180] "GOBP_POSITIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                     
    ## [181] "GOBP_REGULATION_OF_SYNAPTIC_VESICLE_RECYCLING"                                 
    ## [182] "GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_BALANCE"

``` r
not_in_overlap <- common_paths_main[!overlap]
not_in_overlap # 2 divergent by definition
```

    ## [1] "GOBP_NEUROTRANSMITTER_RECEPTOR_TRANSPORT_TO_PLASMA_MEMBRANE"       
    ## [2] "GOBP_REGULATION_OF_POSTSYNAPTIC_NEUROTRANSMITTER_RECEPTOR_ACTIVITY"

``` r
all_paths_overlap <- all_paths %>% filter(pathway %in% overlap_paths)
all_paths_overlap
```

    ##      variant
    ##       <char>
    ##   1:    rare
    ##   2:    rare
    ##   3:    rare
    ##   4:    rare
    ##   5:    rare
    ##  ---        
    ## 360:  common
    ## 361:  common
    ## 362:  common
    ## 363:  common
    ## 364:  common
    ##                                                                             pathway
    ##                                                                              <char>
    ##   1:                                          GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE
    ##   2:                                                  GOBP_NEURON_APOPTOTIC_PROCESS
    ##   3:                                                              GOBP_NEURON_DEATH
    ##   4:                  GOBP_POSITIVE_REGULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL
    ##   5:                                               GOBP_EXCITATORY_SYNAPSE_ASSEMBLY
    ##  ---                                                                               
    ## 360:                                               GOBP_EXCITATORY_SYNAPSE_ASSEMBLY
    ## 361: GOBP_ADENYLATE_CYCLASE_INHIBITING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY
    ## 362:                      GOBP_POSITIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT
    ## 363:                                  GOBP_REGULATION_OF_SYNAPTIC_VESICLE_RECYCLING
    ## 364:                                 GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_BALANCE
    ##              pval       padj    log2err         ES        NES  size
    ##             <num>      <num>      <num>      <num>      <num> <int>
    ##   1: 0.0002756957 0.05017661 0.49849311  0.9924235  2.5347658    19
    ##   2: 0.0007189713 0.06542639 0.47727082  0.9695995  2.2843970   218
    ##   3: 0.0023493938 0.08450704 0.43170770  0.9595688  2.2122378   325
    ##   4: 0.0025713793 0.08450704 0.43170770  0.9160108  2.3444566    23
    ##   5: 0.0027933652 0.08450704 0.43170770  0.9278930  2.3799155    24
    ##  ---                                                               
    ## 360: 0.9857346648 0.99808222 0.03056081  0.1747264  0.5405930    25
    ## 361: 0.9900373599 0.99808222 0.02322483  0.1533360  0.5893707    75
    ## 362: 0.9909808343 0.99808222 0.01699711  0.1600589  0.6670175   139
    ## 363: 0.9926578561 0.99808222 0.03163697  0.1610146  0.4673024    19
    ## 364: 1.0000000000 1.00000000 0.08266464 -0.1427050 -0.6048106    49
    ##       leadingEdge
    ##            <list>
    ##   1:  8831, 22941
    ##   2: 6326, 88....
    ##   3: 6326, 88....
    ##   4: 5728, 85....
    ##   5: 5728, 85....
    ##  ---             
    ## 360: 64101, 7....
    ## 361: 2550, 29....
    ## 362: 51517, 1....
    ## 363: 6622, 60....
    ## 364: 7401, 64....

``` r
nrow(all_paths_overlap) # 364
```

    ## [1] 364

``` r
all_paths_overlap %>% filter(variant == "rare") %>% nrow() # 182
```

    ## [1] 182

``` r
all_paths_overlap %>% filter(variant == "common") %>% nrow() # 182
```

    ## [1] 182

## Visualize NES Distribution Calculate Difference in NES

``` r
NES_hist <- ggplot(all_paths_overlap, aes(x = NES, fill = variant)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  scale_fill_manual(values = c("lightblue", "#F4D03F"))+
  theme_classic()

NES_hist # shows right skew for rare, will use MAD for defining convergence
```

![](fgsea_files/figure-gfm/visualize-NES-1.png)<!-- -->

``` r
#ggsave("NES_distributions.png", plot = NES_hist, path = "../figures", width = 14, height = 10, dpi = 300)
```

``` r
variant_path_NES_padj <- all_paths_overlap %>%
  dplyr::select(variant, pathway, NES, padj) %>%
  pivot_wider(
    names_from = variant,
    values_from = c(NES, padj),
    names_glue = "{variant}_{.value}"
  )

variant_path_NES_padj
```

    ## # A tibble: 182 × 5
    ##    pathway                             rare_NES common_NES rare_padj common_padj
    ##    <chr>                                  <dbl>      <dbl>     <dbl>       <dbl>
    ##  1 GOBP_MAINTENANCE_OF_SYNAPSE_STRUCT…     2.53      1.52     0.0502       0.312
    ##  2 GOBP_NEURON_APOPTOTIC_PROCESS           2.28      1.19     0.0654       0.441
    ##  3 GOBP_NEURON_DEATH                       2.21      1.33     0.0845       0.140
    ##  4 GOBP_POSITIVE_REGULATION_OF_EXCITA…     2.34     -0.999    0.0845       0.692
    ##  5 GOBP_EXCITATORY_SYNAPSE_ASSEMBLY        2.38      0.541    0.0845       0.998
    ##  6 GOBP_POSITIVE_REGULATION_OF_AXON_E…     2.32      0.738    0.0845       0.961
    ##  7 GOBP_REGULATION_OF_NEURONAL_SYNAPT…     2.53      1.23     0.0845       0.510
    ##  8 GOBP_POSTSYNAPTIC_SPECIALIZATION_O…     2.54      1.06     0.0845       0.666
    ##  9 GOBP_REGULATION_OF_LONG_TERM_NEURO…     2.54      1.18     0.0845       0.510
    ## 10 GOBP_POSTSYNAPTIC_SPECIALIZATION_A…     2.39      0.533    0.0845       0.998
    ## # ℹ 172 more rows

## Calculate Unfiltered MAD

``` r
variant_path_NES_padj <- variant_path_NES_padj %>% mutate(
  NES_diff = rare_NES - common_NES, 
  MAD_unfilt = mad(NES_diff),
  median_NES_diff = median(NES_diff)) %>% relocate(NES_diff, .before = MAD_unfilt)
variant_path_NES_padj
```

    ## # A tibble: 182 × 8
    ##    pathway         rare_NES common_NES rare_padj common_padj NES_diff MAD_unfilt
    ##    <chr>              <dbl>      <dbl>     <dbl>       <dbl>    <dbl>      <dbl>
    ##  1 GOBP_MAINTENAN…     2.53      1.52     0.0502       0.312    1.02       0.581
    ##  2 GOBP_NEURON_AP…     2.28      1.19     0.0654       0.441    1.09       0.581
    ##  3 GOBP_NEURON_DE…     2.21      1.33     0.0845       0.140    0.881      0.581
    ##  4 GOBP_POSITIVE_…     2.34     -0.999    0.0845       0.692    3.34       0.581
    ##  5 GOBP_EXCITATOR…     2.38      0.541    0.0845       0.998    1.84       0.581
    ##  6 GOBP_POSITIVE_…     2.32      0.738    0.0845       0.961    1.59       0.581
    ##  7 GOBP_REGULATIO…     2.53      1.23     0.0845       0.510    1.29       0.581
    ##  8 GOBP_POSTSYNAP…     2.54      1.06     0.0845       0.666    1.48       0.581
    ##  9 GOBP_REGULATIO…     2.54      1.18     0.0845       0.510    1.36       0.581
    ## 10 GOBP_POSTSYNAP…     2.39      0.533    0.0845       0.998    1.86       0.581
    ## # ℹ 172 more rows
    ## # ℹ 1 more variable: median_NES_diff <dbl>

``` r
#unfiltered MA for robustness check: 0.5844998
```

## Unfiltered Convergence (or Divergence) Classification

``` r
#threshold of 1 MAD away from the median for divergence:

#unfiltered median ± 1 MAD = 0.75 ± 0.58
#unfiltered convergence zone: 0.17 to 1.33

most_divergent <- variant_path_NES_padj %>% filter(NES_diff == max(NES_diff)) %>% dplyr::select(pathway, NES_diff)
most_convergent <- variant_path_NES_padj %>% filter(abs(NES_diff) == min(abs(NES_diff))) %>% dplyr::select(pathway, NES_diff)

most_divergent
```

    ## # A tibble: 1 × 2
    ##   pathway                                              NES_diff
    ##   <chr>                                                   <dbl>
    ## 1 GOBP_MODULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL     3.36

``` r
most_convergent
```

    ## # A tibble: 1 × 2
    ##   pathway                                 NES_diff
    ##   <chr>                                      <dbl>
    ## 1 GOBP_NEUROMUSCULAR_JUNCTION_DEVELOPMENT  -0.0150

``` r
classify_paths <- variant_path_NES_padj %>% mutate(
  converge_diverge =
    case_when(
      (NES_diff > 1.33 | NES_diff < 0.17) ~ "divergent",
      (NES_diff <= 1.33 & NES_diff >= 0.17) ~ "convergent"
    )
)
classify_paths
```

    ## # A tibble: 182 × 9
    ##    pathway         rare_NES common_NES rare_padj common_padj NES_diff MAD_unfilt
    ##    <chr>              <dbl>      <dbl>     <dbl>       <dbl>    <dbl>      <dbl>
    ##  1 GOBP_MAINTENAN…     2.53      1.52     0.0502       0.312    1.02       0.581
    ##  2 GOBP_NEURON_AP…     2.28      1.19     0.0654       0.441    1.09       0.581
    ##  3 GOBP_NEURON_DE…     2.21      1.33     0.0845       0.140    0.881      0.581
    ##  4 GOBP_POSITIVE_…     2.34     -0.999    0.0845       0.692    3.34       0.581
    ##  5 GOBP_EXCITATOR…     2.38      0.541    0.0845       0.998    1.84       0.581
    ##  6 GOBP_POSITIVE_…     2.32      0.738    0.0845       0.961    1.59       0.581
    ##  7 GOBP_REGULATIO…     2.53      1.23     0.0845       0.510    1.29       0.581
    ##  8 GOBP_POSTSYNAP…     2.54      1.06     0.0845       0.666    1.48       0.581
    ##  9 GOBP_REGULATIO…     2.54      1.18     0.0845       0.510    1.36       0.581
    ## 10 GOBP_POSTSYNAP…     2.39      0.533    0.0845       0.998    1.86       0.581
    ## # ℹ 172 more rows
    ## # ℹ 2 more variables: median_NES_diff <dbl>, converge_diverge <chr>

``` r
divergent_unfilt <- classify_paths %>% filter(converge_diverge == "divergent")
divergent_unfilt
```

    ## # A tibble: 61 × 9
    ##    pathway         rare_NES common_NES rare_padj common_padj NES_diff MAD_unfilt
    ##    <chr>              <dbl>      <dbl>     <dbl>       <dbl>    <dbl>      <dbl>
    ##  1 GOBP_POSITIVE_…     2.34     -0.999    0.0845       0.692     3.34      0.581
    ##  2 GOBP_EXCITATOR…     2.38      0.541    0.0845       0.998     1.84      0.581
    ##  3 GOBP_POSITIVE_…     2.32      0.738    0.0845       0.961     1.59      0.581
    ##  4 GOBP_POSTSYNAP…     2.54      1.06     0.0845       0.666     1.48      0.581
    ##  5 GOBP_REGULATIO…     2.54      1.18     0.0845       0.510     1.36      0.581
    ##  6 GOBP_POSTSYNAP…     2.39      0.533    0.0845       0.998     1.86      0.581
    ##  7 GOBP_POSITIVE_…     2.26      0.912    0.0845       0.780     1.35      0.581
    ##  8 GOBP_SYNAPTIC_…     2.19      0.789    0.0845       0.909     1.40      0.581
    ##  9 GOBP_NEURONAL_…     2.53      0.917    0.0845       0.780     1.62      0.581
    ## 10 GOBP_POSTSYNAP…     2.37      0.721    0.0845       0.956     1.65      0.581
    ## # ℹ 51 more rows
    ## # ℹ 2 more variables: median_NES_diff <dbl>, converge_diverge <chr>

``` r
nrow(divergent_unfilt)
```

    ## [1] 61

``` r
convergent_unfilt <- classify_paths %>% filter(converge_diverge == "convergent")
convergent_unfilt
```

    ## # A tibble: 121 × 9
    ##    pathway         rare_NES common_NES rare_padj common_padj NES_diff MAD_unfilt
    ##    <chr>              <dbl>      <dbl>     <dbl>       <dbl>    <dbl>      <dbl>
    ##  1 GOBP_MAINTENAN…     2.53      1.52     0.0502       0.312    1.02       0.581
    ##  2 GOBP_NEURON_AP…     2.28      1.19     0.0654       0.441    1.09       0.581
    ##  3 GOBP_NEURON_DE…     2.21      1.33     0.0845       0.140    0.881      0.581
    ##  4 GOBP_REGULATIO…     2.53      1.23     0.0845       0.510    1.29       0.581
    ##  5 GOBP_SYNAPSE_M…     2.27      1.40     0.0845       0.413    0.876      0.581
    ##  6 GOBP_POSTSYNAP…     2.24      1.21     0.0845       0.510    1.03       0.581
    ##  7 GOBP_NEGATIVE_…     2.47      1.28     0.0845       0.446    1.18       0.581
    ##  8 GOBP_GLUTAMATE…     2.12      1.08     0.0845       0.627    1.04       0.581
    ##  9 GOBP_REGULATIO…     2.11      0.797    0.0845       0.904    1.31       0.581
    ## 10 GOBP_NEGATIVE_…     2.28      1.39     0.0845       0.325    0.897      0.581
    ## # ℹ 111 more rows
    ## # ℹ 2 more variables: median_NES_diff <dbl>, converge_diverge <chr>

``` r
nrow(convergent_unfilt)
```

    ## [1] 121

``` r
#determine padj filter by visualizing NES difference against -log10padj

rare_NES_padj <- ggplot(classify_paths, aes(x = NES_diff, y = -log10(rare_padj), color = converge_diverge)) +
  geom_jitter(alpha = 0.8) +
  scale_color_manual(values = c("convergent" = "lightblue", "divergent" = "#F4D03F")) +
  xlab("Difference in NES")+
  ylab("Rare Variant -log10(padj")+
  theme_classic()

common_NES_padj <- ggplot(classify_paths, aes(x = NES_diff, y = -log10(common_padj), color = converge_diverge)) +
  geom_jitter(alpha = 0.8) +
  scale_color_manual(values = c("convergent" = "lightblue", "divergent" = "#F4D03F")) +
  xlab("Difference in NES")+
  ylab("Common Variant -log10(padj")+
  theme_classic()

NES_padj <- rare_NES_padj + common_NES_padj

NES_padj
```

![](fgsea_files/figure-gfm/visualize-NES-padj-1.png)<!-- -->

``` r
#ggsave("NES_distributions.png", plot = NES_hist, path = "../figures", width = 14, height = 10, dpi = 300)
```

## Apply FDR Pathway Filter

``` r
# apply a filter of padj < 0.25 to rare and common arms
rare_sig <- rare_pathways %>% filter(padj < 0.25) %>% pull(pathway)
common_sig <- common_pathways %>% filter(padj < 0.25) %>% pull(pathway)
both_sig <- intersect(rare_sig, common_sig)
all_paths_filt <- all_paths_overlap %>% filter(pathway %in% both_sig)
nrow(all_paths_filt)
```

    ## [1] 38

``` r
setdiff(rare_sig, common_sig)  #rare only
```

    ##   [1] "GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE"                                         
    ##   [2] "GOBP_NEURON_APOPTOTIC_PROCESS"                                                 
    ##   [3] "GOBP_POSITIVE_REGULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL"                 
    ##   [4] "GOBP_EXCITATORY_SYNAPSE_ASSEMBLY"                                              
    ##   [5] "GOBP_POSITIVE_REGULATION_OF_AXON_EXTENSION"                                    
    ##   [6] "GOBP_REGULATION_OF_NEURONAL_SYNAPTIC_PLASTICITY"                               
    ##   [7] "GOBP_POSTSYNAPTIC_SPECIALIZATION_ORGANIZATION"                                 
    ##   [8] "GOBP_REGULATION_OF_LONG_TERM_NEURONAL_SYNAPTIC_PLASTICITY"                     
    ##   [9] "GOBP_POSTSYNAPTIC_SPECIALIZATION_ASSEMBLY"                                     
    ##  [10] "GOBP_SYNAPSE_MATURATION"                                                       
    ##  [11] "GOBP_POSITIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"               
    ##  [12] "GOBP_SYNAPTIC_TRANSMISSION_GABAERGIC"                                          
    ##  [13] "GOBP_NEURONAL_ACTION_POTENTIAL"                                                
    ##  [14] "GOBP_POSTSYNAPSE_ASSEMBLY"                                                     
    ##  [15] "GOBP_POSTSYNAPTIC_MEMBRANE_ORGANIZATION"                                       
    ##  [16] "GOBP_NEGATIVE_REGULATION_OF_AXONOGENESIS"                                      
    ##  [17] "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GABAERGIC"                            
    ##  [18] "GOBP_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY"                                     
    ##  [19] "GOBP_POSITIVE_REGULATION_OF_SYNAPSE_ASSEMBLY"                                  
    ##  [20] "GOBP_MODULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL"                          
    ##  [21] "GOBP_REGULATION_OF_CHROMATIN_ORGANIZATION"                                     
    ##  [22] "GOBP_NEGATIVE_REGULATION_OF_NEURON_APOPTOTIC_PROCESS"                          
    ##  [23] "GOBP_NEGATIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION"                             
    ##  [24] "GOBP_REGULATION_OF_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY"                          
    ##  [25] "GOBP_REGULATION_OF_POSTSYNAPSE_ORGANIZATION"                                   
    ##  [26] "GOBP_POSITIVE_REGULATION_OF_LONG_TERM_SYNAPTIC_POTENTIATION"                   
    ##  [27] "GOBP_LONG_TERM_SYNAPTIC_POTENTIATION"                                          
    ##  [28] "GOBP_MIDBRAIN_DOPAMINERGIC_NEURON_DIFFERENTIATION"                             
    ##  [29] "GOBP_LONG_TERM_SYNAPTIC_DEPRESSION"                                            
    ##  [30] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_AXONOGENESIS"                               
    ##  [31] "GOBP_CHEMICAL_SYNAPTIC_TRANSMISSION_POSTSYNAPTIC"                              
    ##  [32] "GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_REGENERATION"                    
    ##  [33] "GOBP_FOREBRAIN_GENERATION_OF_NEURONS"                                          
    ##  [34] "GOBP_AXON_ENSHEATHMENT_IN_CENTRAL_NERVOUS_SYSTEM"                              
    ##  [35] "GOBP_FOREBRAIN_NEURON_DIFFERENTIATION"                                         
    ##  [36] "GOBP_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_POTENTIAL"                            
    ##  [37] "GOBP_DOPAMINERGIC_NEURON_DIFFERENTIATION"                                      
    ##  [38] "GOBP_INHIBITORY_SYNAPSE_ASSEMBLY"                                              
    ##  [39] "GOBP_NEGATIVE_REGULATION_OF_OXIDATIVE_STRESS_INDUCED_NEURON_DEATH"             
    ##  [40] "GOBP_NEUROMUSCULAR_PROCESS"                                                    
    ##  [41] "GOBP_REGULATION_OF_NEURON_PROJECTION_REGENERATION"                             
    ##  [42] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DEVELOPMENT"                                
    ##  [43] "GOBP_NEUROTRANSMITTER_REUPTAKE"                                                
    ##  [44] "GOBP_NEURON_PROJECTION_ORGANIZATION"                                           
    ##  [45] "GOBP_ENSHEATHMENT_OF_NEURONS"                                                  
    ##  [46] "GOBP_POSITIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                     
    ##  [47] "GOBP_POSITIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION"                             
    ##  [48] "GOBP_REGULATION_OF_AXONOGENESIS"                                               
    ##  [49] "GOBP_RETROGRADE_AXONAL_TRANSPORT"                                              
    ##  [50] "GOBP_POSTSYNAPSE_ORGANIZATION"                                                 
    ##  [51] "GOBP_PRESYNAPSE_ORGANIZATION"                                                  
    ##  [52] "GOBP_REGULATION_OF_LONG_TERM_SYNAPTIC_POTENTIATION"                            
    ##  [53] "GOBP_NEUROTRANSMITTER_UPTAKE"                                                  
    ##  [54] "GOBP_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"                                      
    ##  [55] "GOBP_SYNAPTIC_VESICLE_LOCALIZATION"                                            
    ##  [56] "GOBP_REGULATION_OF_SYNAPSE_ASSEMBLY"                                           
    ##  [57] "GOBP_NEGATIVE_REGULATION_OF_NEURON_DEATH"                                      
    ##  [58] "GOBP_NEURON_PROJECTION_EXTENSION"                                              
    ##  [59] "GOBP_SYNAPSE_ASSEMBLY"                                                         
    ##  [60] "GOBP_NEURON_FATE_COMMITMENT"                                                   
    ##  [61] "GOBP_NEURON_PROJECTION_REGENERATION"                                           
    ##  [62] "GOBP_POSITIVE_REGULATION_OF_NEURON_APOPTOTIC_PROCESS"                          
    ##  [63] "GOBP_AXON_EXTENSION"                                                           
    ##  [64] "GOBP_NEURON_FATE_SPECIFICATION"                                                
    ##  [65] "GOBP_REGULATION_OF_SYNAPTIC_PLASTICITY"                                        
    ##  [66] "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"                        
    ##  [67] "GOBP_POSITIVE_REGULATION_OF_NEUROBLAST_PROLIFERATION"                          
    ##  [68] "GOBP_SYNAPTIC_MEMBRANE_ADHESION"                                               
    ##  [69] "GOBP_REGULATION_OF_NEUROTRANSMITTER_RECEPTOR_ACTIVITY"                         
    ##  [70] "GOBP_REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY"                              
    ##  [71] "GOBP_REGULATION_OF_NEUROBLAST_PROLIFERATION"                                   
    ##  [72] "GOBP_PROTEIN_LOCALIZATION_TO_SYNAPSE"                                          
    ##  [73] "GOBP_HETEROCHROMATIN_ORGANIZATION"                                             
    ##  [74] "GOBP_NEURON_DEATH_IN_RESPONSE_TO_OXIDATIVE_STRESS"                             
    ##  [75] "GOBP_FACULTATIVE_HETEROCHROMATIN_ASSEMBLY"                                     
    ##  [76] "GOBP_NEUROTRANSMITTER_TRANSPORT"                                               
    ##  [77] "KEGG_AXON_GUIDANCE"                                                            
    ##  [78] "GOBP_AXONAL_TRANSPORT"                                                         
    ##  [79] "GOBP_POSITIVE_REGULATION_OF_CHROMATIN_ORGANIZATION"                            
    ##  [80] "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING"                                   
    ##  [81] "GOBP_DNA_METHYLATION_DEPENDENT_HETEROCHROMATIN_ASSEMBLY"                       
    ##  [82] "GOBP_REGULATION_OF_CHROMATIN_ASSEMBLY"                                         
    ##  [83] "GOBP_RESPONSE_TO_AXON_INJURY"                                                  
    ##  [84] "GOBP_SYNAPTIC_VESICLE_MEMBRANE_ORGANIZATION"                                   
    ##  [85] "GOBP_SYNAPTIC_VESICLE_TRANSPORT"                                               
    ##  [86] "GOBP_REGULATION_OF_NEURON_DIFFERENTIATION"                                     
    ##  [87] "GOBP_CALCIUM_ION_REGULATED_EXOCYTOSIS_OF_NEUROTRANSMITTER"                     
    ##  [88] "GOBP_NEURON_CELLULAR_HOMEOSTASIS"                                              
    ##  [89] "GOBP_NEURON_RECOGNITION"                                                       
    ##  [90] "GOBP_GLUTAMATE_SECRETION"                                                      
    ##  [91] "GOBP_POSITIVE_REGULATION_OF_NEUROTRANSMITTER_SECRETION"                        
    ##  [92] "GOBP_NEUROBLAST_PROLIFERATION"                                                 
    ##  [93] "GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_POSTURE"                                
    ##  [94] "GOBP_SYNAPTIC_VESICLE_PRIMING"                                                 
    ##  [95] "GOBP_AXONAL_TRANSPORT_OF_MITOCHONDRION"                                        
    ##  [96] "GOBP_CHROMATIN_DISASSEMBLY"                                                    
    ##  [97] "GOBP_NEURON_MIGRATION"                                                         
    ##  [98] "GOBP_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                              
    ##  [99] "GOBP_ANTEROGRADE_AXONAL_TRANSPORT"                                             
    ## [100] "GOBP_CEREBRAL_CORTEX_NEURON_DIFFERENTIATION"                                   
    ## [101] "GOBP_NEUROMUSCULAR_SYNAPTIC_TRANSMISSION"                                      
    ## [102] "GOBP_NEURON_MATURATION"                                                        
    ## [103] "GOBP_AXON_DEVELOPMENT"                                                         
    ## [104] "GOBP_POSITIVE_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT"                        
    ## [105] "GOBP_RETINAL_GANGLION_CELL_AXON_GUIDANCE"                                      
    ## [106] "GOBP_REGULATION_OF_NEUROTRANSMITTER_UPTAKE"                                    
    ## [107] "GOBP_RNA_MEDIATED_GENE_SILENCING_BY_INHIBITION_OF_TRANSLATION"                 
    ## [108] "GOBP_REGULATION_OF_NEURON_MIGRATION"                                           
    ## [109] "GOBP_REGULATION_OF_GLUTAMATE_SECRETION"                                        
    ## [110] "GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_BALANCE"                                
    ## [111] "GOBP_POSITIVE_REGULATION_OF_NEURON_DIFFERENTIATION"                            
    ## [112] "GOBP_REGULATION_OF_PRESYNAPSE_ORGANIZATION"                                    
    ## [113] "GOBP_NEURON_PROJECTION_ARBORIZATION"                                           
    ## [114] "GOBP_NEUROTROPHIN_SIGNALING_PATHWAY"                                           
    ## [115] "GOBP_POSITIVE_REGULATION_OF_NEURON_MIGRATION"                                  
    ## [116] "GOBP_RESPONSE_TO_LEUKEMIA_INHIBITORY_FACTOR"                                   
    ## [117] "GOBP_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT"                                 
    ## [118] "GOBP_NEUROTRANSMITTER_RECEPTOR_INTERNALIZATION"                                
    ## [119] "GOBP_NEURONAL_STEM_CELL_POPULATION_MAINTENANCE"                                
    ## [120] "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION"                                  
    ## [121] "GOBP_RECEPTOR_LOCALIZATION_TO_SYNAPSE"                                         
    ## [122] "GOBP_NEUROTRANSMITTER_RECEPTOR_TRANSPORT"                                      
    ## [123] "GOBP_REGULATION_OF_PROTEIN_LOCALIZATION_TO_SYNAPSE"                            
    ## [124] "GOBP_REGULATION_OF_RECEPTOR_LOCALIZATION_TO_SYNAPSE"                           
    ## [125] "GOBP_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_NEUROTRANSMITTER_RECEPTOR_LEVELS"     
    ## [126] "GOBP_FOREBRAIN_NEURON_DEVELOPMENT"                                             
    ## [127] "GOBP_AXONEMAL_DYNEIN_COMPLEX_ASSEMBLY"                                         
    ## [128] "GOBP_SYNAPTIC_VESICLE_RECYCLING"                                               
    ## [129] "GOBP_NEGATIVE_REGULATION_OF_NEURON_DIFFERENTIATION"                            
    ## [130] "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY"                                           
    ## [131] "GOBP_PRESYNAPTIC_ENDOCYTOSIS"                                                  
    ## [132] "GOBP_SYNAPTIC_VESICLE_CYTOSKELETAL_TRANSPORT"                                  
    ## [133] "GOBP_PROTEIN_LOCALIZATION_TO_POSTSYNAPSE"                                      
    ## [134] "GOBP_ADENYLATE_CYCLASE_INHIBITING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY"
    ## [135] "GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION"                                    
    ## [136] "GOBP_PROTEIN_LOCALIZATION_TO_CHROMATIN"                                        
    ## [137] "GOBP_REGULATION_OF_CHROMATIN_BINDING"                                          
    ## [138] "GOBP_NEUROTROPHIN_TRK_RECEPTOR_SIGNALING_PATHWAY"                              
    ## [139] "GOBP_AXONAL_FASCICULATION"                                                     
    ## [140] "GOBP_POSTSYNAPTIC_SIGNAL_TRANSDUCTION"                                         
    ## [141] "GOBP_REGULATION_OF_DNA_METHYLATION_DEPENDENT_HETEROCHROMATIN_ASSEMBLY"         
    ## [142] "GOBP_AXONEME_ASSEMBLY"

``` r
setdiff(common_sig, rare_sig)  #common only
```

    ## [1] "KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM"

``` r
all_paths_filt_export <- all_paths_filt %>% mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
#write.csv(all_paths_filt_export, "../results/all_path_data.csv")
```

``` r
variant_path_NES <- all_paths_filt %>%
  dplyr::select(variant, pathway, NES) %>%
  pivot_wider(
    names_from = variant,
    values_from = NES,
    names_glue = "{variant}_NES"
  )

variant_path_NES
```

    ## # A tibble: 19 × 3
    ##    pathway                                                   rare_NES common_NES
    ##    <chr>                                                        <dbl>      <dbl>
    ##  1 GOBP_NEURON_DEATH                                             2.21       1.33
    ##  2 GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH                      1.98       1.54
    ##  3 GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT     2.34       1.39
    ##  4 GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION            1.79       1.50
    ##  5 GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS                      2.03       1.47
    ##  6 GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS                      1.76       1.43
    ##  7 GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE                    1.60       1.45
    ##  8 GOBP_CHROMATIN_REMODELING                                     2.22       1.62
    ##  9 GOBP_NEURON_PROJECTION_GUIDANCE                               1.55       1.39
    ## 10 GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS                    1.56       1.38
    ## 11 GOBP_SYNAPSE_ORGANIZATION                                     2.11       1.40
    ## 12 GOBP_REGULATION_OF_NEUROGENESIS                               2.13       1.36
    ## 13 GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY                        1.49       1.46
    ## 14 GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJ…     1.80       1.82
    ## 15 GOBP_NEUROTRANSMITTER_SECRETION                               1.48       1.49
    ## 16 GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENES…     1.65       1.60
    ## 17 GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_A…     1.45       1.64
    ## 18 GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS                              1.34       1.43
    ## 19 KEGG_MTOR_SIGNALING_PATHWAY                                   1.38       1.60

``` r
#sum(is.na(variant_path_NES))
```

## Calculate Filtered MAD

``` r
variant_path_NES <- variant_path_NES %>% mutate(
  NES_diff = rare_NES - common_NES, 
  MAD_filt = mad(NES_diff),
  median_NES_diff = median(NES_diff)) %>% relocate(NES_diff, .before = MAD_filt)

variant_path_NES
```

    ## # A tibble: 19 × 6
    ##    pathway                 rare_NES common_NES NES_diff MAD_filt median_NES_diff
    ##    <chr>                      <dbl>      <dbl>    <dbl>    <dbl>           <dbl>
    ##  1 GOBP_NEURON_DEATH           2.21       1.33   0.881     0.397           0.175
    ##  2 GOBP_POSITIVE_REGULATI…     1.98       1.54   0.443     0.397           0.175
    ##  3 GOBP_NEGATIVE_REGULATI…     2.34       1.39   0.959     0.397           0.175
    ##  4 GOBP_CENTRAL_NERVOUS_S…     1.79       1.50   0.298     0.397           0.175
    ##  5 GOBP_POSITIVE_REGULATI…     2.03       1.47   0.558     0.397           0.175
    ##  6 GOBP_POSITIVE_REGULATI…     1.76       1.43   0.325     0.397           0.175
    ##  7 GOBP_VESICLE_MEDIATED_…     1.60       1.45   0.154     0.397           0.175
    ##  8 GOBP_CHROMATIN_REMODEL…     2.22       1.62   0.603     0.397           0.175
    ##  9 GOBP_NEURON_PROJECTION…     1.55       1.39   0.161     0.397           0.175
    ## 10 GOBP_REGULATION_OF_NEU…     1.56       1.38   0.175     0.397           0.175
    ## 11 GOBP_SYNAPSE_ORGANIZAT…     2.11       1.40   0.709     0.397           0.175
    ## 12 GOBP_REGULATION_OF_NEU…     2.13       1.36   0.775     0.397           0.175
    ## 13 GOBP_CHROMATIN_ASSEMBL…     1.49       1.46   0.0333    0.397           0.175
    ## 14 GOBP_NEURON_PROJECTION…     1.80       1.82  -0.0247    0.397           0.175
    ## 15 GOBP_NEUROTRANSMITTER_…     1.48       1.49  -0.0159    0.397           0.175
    ## 16 GOBP_CENTRAL_NERVOUS_S…     1.65       1.60   0.0508    0.397           0.175
    ## 17 GOBP_NEGATIVE_REGULATI…     1.45       1.64  -0.190     0.397           0.175
    ## 18 GOBP_SYNAPTIC_VESICLE_…     1.34       1.43  -0.0951    0.397           0.175
    ## 19 KEGG_MTOR_SIGNALING_PA…     1.38       1.60  -0.226     0.397           0.175

``` r
#unfiltered MAD for robustness check: 0.5844998
mad_filt <- 0.3576003
median_filt <- 0.2000136
```

## Re-Classify

``` r
#threshold of 1 MAD away from the median for divergence:
#new median ± 1 MAD = 0.20 ± 0.36
#new convergence zone: -0.16 to 0.56

most_divergent <- variant_path_NES %>% filter(NES_diff == max(NES_diff)) %>% dplyr::select(pathway, NES_diff)
most_convergent <- variant_path_NES %>% filter(abs(NES_diff) == min(abs(NES_diff))) %>% dplyr::select(pathway, NES_diff)

most_divergent
```

    ## # A tibble: 1 × 2
    ##   pathway                                                   NES_diff
    ##   <chr>                                                        <dbl>
    ## 1 GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT    0.959

``` r
most_convergent
```

    ## # A tibble: 1 × 2
    ##   pathway                         NES_diff
    ##   <chr>                              <dbl>
    ## 1 GOBP_NEUROTRANSMITTER_SECRETION  -0.0159

``` r
classify_paths_filt <- variant_path_NES %>% mutate(
  converge_diverge =
    case_when(
      (NES_diff > 0.56 | NES_diff < -0.16) ~ "divergent",
      (NES_diff <= 0.56 & NES_diff >= -0.16) ~ "convergent"
    )
)
classify_paths_filt
```

    ## # A tibble: 19 × 7
    ##    pathway                 rare_NES common_NES NES_diff MAD_filt median_NES_diff
    ##    <chr>                      <dbl>      <dbl>    <dbl>    <dbl>           <dbl>
    ##  1 GOBP_NEURON_DEATH           2.21       1.33   0.881     0.397           0.175
    ##  2 GOBP_POSITIVE_REGULATI…     1.98       1.54   0.443     0.397           0.175
    ##  3 GOBP_NEGATIVE_REGULATI…     2.34       1.39   0.959     0.397           0.175
    ##  4 GOBP_CENTRAL_NERVOUS_S…     1.79       1.50   0.298     0.397           0.175
    ##  5 GOBP_POSITIVE_REGULATI…     2.03       1.47   0.558     0.397           0.175
    ##  6 GOBP_POSITIVE_REGULATI…     1.76       1.43   0.325     0.397           0.175
    ##  7 GOBP_VESICLE_MEDIATED_…     1.60       1.45   0.154     0.397           0.175
    ##  8 GOBP_CHROMATIN_REMODEL…     2.22       1.62   0.603     0.397           0.175
    ##  9 GOBP_NEURON_PROJECTION…     1.55       1.39   0.161     0.397           0.175
    ## 10 GOBP_REGULATION_OF_NEU…     1.56       1.38   0.175     0.397           0.175
    ## 11 GOBP_SYNAPSE_ORGANIZAT…     2.11       1.40   0.709     0.397           0.175
    ## 12 GOBP_REGULATION_OF_NEU…     2.13       1.36   0.775     0.397           0.175
    ## 13 GOBP_CHROMATIN_ASSEMBL…     1.49       1.46   0.0333    0.397           0.175
    ## 14 GOBP_NEURON_PROJECTION…     1.80       1.82  -0.0247    0.397           0.175
    ## 15 GOBP_NEUROTRANSMITTER_…     1.48       1.49  -0.0159    0.397           0.175
    ## 16 GOBP_CENTRAL_NERVOUS_S…     1.65       1.60   0.0508    0.397           0.175
    ## 17 GOBP_NEGATIVE_REGULATI…     1.45       1.64  -0.190     0.397           0.175
    ## 18 GOBP_SYNAPTIC_VESICLE_…     1.34       1.43  -0.0951    0.397           0.175
    ## 19 KEGG_MTOR_SIGNALING_PA…     1.38       1.60  -0.226     0.397           0.175
    ## # ℹ 1 more variable: converge_diverge <chr>

``` r
divergent <- classify_paths_filt %>% filter(converge_diverge == "divergent")
divergent
```

    ## # A tibble: 7 × 7
    ##   pathway rare_NES common_NES NES_diff MAD_filt median_NES_diff converge_diverge
    ##   <chr>      <dbl>      <dbl>    <dbl>    <dbl>           <dbl> <chr>           
    ## 1 GOBP_N…     2.21       1.33    0.881    0.397           0.175 divergent       
    ## 2 GOBP_N…     2.34       1.39    0.959    0.397           0.175 divergent       
    ## 3 GOBP_C…     2.22       1.62    0.603    0.397           0.175 divergent       
    ## 4 GOBP_S…     2.11       1.40    0.709    0.397           0.175 divergent       
    ## 5 GOBP_R…     2.13       1.36    0.775    0.397           0.175 divergent       
    ## 6 GOBP_N…     1.45       1.64   -0.190    0.397           0.175 divergent       
    ## 7 KEGG_M…     1.38       1.60   -0.226    0.397           0.175 divergent

``` r
nrow(divergent) # 6
```

    ## [1] 7

``` r
convergent <- classify_paths_filt %>% filter(converge_diverge == "convergent")
convergent
```

    ## # A tibble: 12 × 7
    ##    pathway                 rare_NES common_NES NES_diff MAD_filt median_NES_diff
    ##    <chr>                      <dbl>      <dbl>    <dbl>    <dbl>           <dbl>
    ##  1 GOBP_POSITIVE_REGULATI…     1.98       1.54   0.443     0.397           0.175
    ##  2 GOBP_CENTRAL_NERVOUS_S…     1.79       1.50   0.298     0.397           0.175
    ##  3 GOBP_POSITIVE_REGULATI…     2.03       1.47   0.558     0.397           0.175
    ##  4 GOBP_POSITIVE_REGULATI…     1.76       1.43   0.325     0.397           0.175
    ##  5 GOBP_VESICLE_MEDIATED_…     1.60       1.45   0.154     0.397           0.175
    ##  6 GOBP_NEURON_PROJECTION…     1.55       1.39   0.161     0.397           0.175
    ##  7 GOBP_REGULATION_OF_NEU…     1.56       1.38   0.175     0.397           0.175
    ##  8 GOBP_CHROMATIN_ASSEMBL…     1.49       1.46   0.0333    0.397           0.175
    ##  9 GOBP_NEURON_PROJECTION…     1.80       1.82  -0.0247    0.397           0.175
    ## 10 GOBP_NEUROTRANSMITTER_…     1.48       1.49  -0.0159    0.397           0.175
    ## 11 GOBP_CENTRAL_NERVOUS_S…     1.65       1.60   0.0508    0.397           0.175
    ## 12 GOBP_SYNAPTIC_VESICLE_…     1.34       1.43  -0.0951    0.397           0.175
    ## # ℹ 1 more variable: converge_diverge <chr>

``` r
nrow(convergent) # 12
```

    ## [1] 12

``` r
#inner join to add converge_diverge onto all_paths_overlap
classify_all_paths_overlap_filt <- inner_join(all_paths_filt, classify_paths_filt)
classify_all_paths_overlap_filt
```

    ##     variant
    ##      <char>
    ##  1:    rare
    ##  2:    rare
    ##  3:    rare
    ##  4:    rare
    ##  5:    rare
    ##  6:    rare
    ##  7:    rare
    ##  8:    rare
    ##  9:    rare
    ## 10:    rare
    ## 11:    rare
    ## 12:    rare
    ## 13:    rare
    ## 14:    rare
    ## 15:    rare
    ## 16:    rare
    ## 17:    rare
    ## 18:    rare
    ## 19:    rare
    ## 20:  common
    ## 21:  common
    ## 22:  common
    ## 23:  common
    ## 24:  common
    ## 25:  common
    ## 26:  common
    ## 27:  common
    ## 28:  common
    ## 29:  common
    ## 30:  common
    ## 31:  common
    ## 32:  common
    ## 33:  common
    ## 34:  common
    ## 35:  common
    ## 36:  common
    ## 37:  common
    ## 38:  common
    ##     variant
    ##                                                                     pathway
    ##                                                                      <char>
    ##  1:                                                       GOBP_NEURON_DEATH
    ##  2:                                GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH
    ##  3:               GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT
    ##  4:                      GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION
    ##  5:                                GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS
    ##  6:                                GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS
    ##  7:                              GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE
    ##  8:                                               GOBP_CHROMATIN_REMODELING
    ##  9:                                         GOBP_NEURON_PROJECTION_GUIDANCE
    ## 10:                              GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS
    ## 11:                                               GOBP_SYNAPSE_ORGANIZATION
    ## 12:                                         GOBP_REGULATION_OF_NEUROGENESIS
    ## 13:                                  GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY
    ## 14: GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE
    ## 15:                                         GOBP_NEUROTRANSMITTER_SECRETION
    ## 16:              GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS
    ## 17:    GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE
    ## 18:                                        GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS
    ## 19:                                             KEGG_MTOR_SIGNALING_PATHWAY
    ## 20:                                               GOBP_CHROMATIN_REMODELING
    ## 21: GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE
    ## 22:                                               GOBP_SYNAPSE_ORGANIZATION
    ## 23:                                         GOBP_REGULATION_OF_NEUROGENESIS
    ## 24:                                         GOBP_NEUROTRANSMITTER_SECRETION
    ## 25:                      GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION
    ## 26:                                  GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY
    ## 27:                                             KEGG_MTOR_SIGNALING_PATHWAY
    ## 28:                              GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE
    ## 29:                                GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH
    ## 30:                                         GOBP_NEURON_PROJECTION_GUIDANCE
    ## 31:                                                       GOBP_NEURON_DEATH
    ## 32:                                GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS
    ## 33:    GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE
    ## 34:              GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS
    ## 35:                                GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS
    ## 36:                              GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS
    ## 37:               GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT
    ## 38:                                        GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS
    ##                                                                     pathway
    ##             pval       padj   log2err        ES      NES  size  leadingEdge
    ##            <num>      <num>     <num>     <num>    <num> <int>       <list>
    ##  1: 0.0023493938 0.08450704 0.4317077 0.9595688 2.212238   325 6326, 88....
    ##  2: 0.0238811681 0.08450704 0.3524879 0.7955304 1.981162    80 2904, 14....
    ##  3: 0.0249750250 0.08450704 0.2878571 0.9651832 2.344443   123 8831, 57....
    ##  4: 0.0299700300 0.08450704 0.2616635 0.7493933 1.793925   157 5728, 85....
    ##  5: 0.0329046665 0.08450704 0.3217759 0.8138417 2.029486    66 23394, 1....
    ##  6: 0.0389610390 0.09208973 0.2279872 0.7457535 1.759605   199 23394, 8....
    ##  7: 0.0409590410 0.09557110 0.2220560 0.6767406 1.603133   183 5728, 14....
    ##  8: 0.0429570430 0.09671180 0.2165428 0.9450217 2.224850   224 57680, 2....
    ##  9: 0.0449550450 0.09671180 0.2114002 0.6573949 1.549044   217 1826, 93....
    ## 10: 0.0479520480 0.09671180 0.2042948 0.6611369 1.559953   199 6529, 93....
    ## 11: 0.0509490509 0.09760766 0.1978220 0.9249977 2.113272   370 8831, 23....
    ## 12: 0.0569430569 0.10363636 0.1864326 0.9269535 2.134081   326 8831, 23....
    ## 13: 0.0589410589 0.10606061 0.1830239 0.6311725 1.489802   196 11011, 1....
    ## 14: 0.0629370629 0.10606061 0.1766943 0.6973979 1.795029    33 1826, 65....
    ## 15: 0.0879120879 0.12800000 0.1473312 0.6095070 1.475186   136 9378, 68....
    ## 16: 0.0899100899 0.12884753 0.1455161 0.6443596 1.648539    25 2047, 47....
    ## 17: 0.1498501499 0.18808777 0.1088201 0.5653401 1.450018    24 6585, 64....
    ## 18: 0.1578421578 0.19524100 0.1055209 0.5450988 1.336301   101 6812, 27....
    ## 19: 0.1598401598 0.19524100 0.1047328 0.5405839 1.376318    48 208, 253....
    ## 20: 0.0001152112 0.02119885 0.5384341 0.3704143 1.621940   227 9759, 21....
    ## 21: 0.0007306983 0.06722425 0.4772708 0.5462061 1.819726    36 7473, 91....
    ## 22: 0.0016485947 0.10111381 0.4550599 0.3106873 1.404026   384 4137, 11....
    ## 23: 0.0036537739 0.12821204 0.4317077 0.3032521 1.359135   338 4137, 74....
    ## 24: 0.0041033213 0.12821204 0.4070179 0.3592682 1.491053   137 783, 990....
    ## 25: 0.0046781638 0.12821204 0.4070179 0.3530124 1.495717   160 4137, 91....
    ## 26: 0.0057849919 0.12978981 0.4070179 0.3368143 1.456460   191 3010, 83....
    ## 27: 0.0063484148 0.12978981 0.4070179 0.4540769 1.602714    49 2475, 60....
    ## 28: 0.0088819892 0.13977450 0.3807304 0.3350995 1.449486   186 783, 990....
    ## 29: 0.0090875811 0.13977450 0.3807304 0.3951873 1.538347    81 4137, 67....
    ## 30: 0.0096026847 0.13977450 0.3807304 0.3180250 1.387951   221 7473, 91....
    ## 31: 0.0098753721 0.13977450 0.3807304 0.2979639 1.331642   327 4137, 67....
    ## 32: 0.0109958509 0.14404240 0.3807304 0.3294199 1.434285   208 4137, 74....
    ## 33: 0.0131143100 0.14404240 0.3807304 0.5234268 1.640478    26 7473, 10....
    ## 34: 0.0137656389 0.14404240 0.3807304 0.5163934 1.597690    25 91584, 1....
    ## 35: 0.0140840022 0.14404240 0.3807304 0.3855854 1.471753    71 4137, 74....
    ## 36: 0.0140911043 0.14404240 0.3807304 0.3190564 1.384651   201 783, 990....
    ## 37: 0.0193602626 0.18748886 0.3524879 0.3374266 1.385743   128 7473, 23....
    ## 38: 0.0215557747 0.19831313 0.3524879 0.3579009 1.431405   101 783, 990....
    ##             pval       padj   log2err        ES      NES  size  leadingEdge
    ##     rare_NES common_NES    NES_diff  MAD_filt median_NES_diff converge_diverge
    ##        <num>      <num>       <num>     <num>           <num>           <char>
    ##  1: 2.212238   1.331642  0.88059543 0.3966149       0.1753014        divergent
    ##  2: 1.981162   1.538347  0.44281445 0.3966149       0.1753014       convergent
    ##  3: 2.344443   1.385743  0.95870057 0.3966149       0.1753014        divergent
    ##  4: 1.793925   1.495717  0.29820775 0.3966149       0.1753014       convergent
    ##  5: 2.029486   1.471753  0.55773251 0.3966149       0.1753014       convergent
    ##  6: 1.759605   1.434285  0.32532045 0.3966149       0.1753014       convergent
    ##  7: 1.603133   1.449486  0.15364703 0.3966149       0.1753014       convergent
    ##  8: 2.224850   1.621940  0.60291002 0.3966149       0.1753014        divergent
    ##  9: 1.549044   1.387951  0.16109275 0.3966149       0.1753014       convergent
    ## 10: 1.559953   1.384651  0.17530136 0.3966149       0.1753014       convergent
    ## 11: 2.113272   1.404026  0.70924642 0.3966149       0.1753014        divergent
    ## 12: 2.134081   1.359135  0.77494577 0.3966149       0.1753014        divergent
    ## 13: 1.489802   1.456460  0.03334152 0.3966149       0.1753014       convergent
    ## 14: 1.795029   1.819726 -0.02469702 0.3966149       0.1753014       convergent
    ## 15: 1.475186   1.491053 -0.01586671 0.3966149       0.1753014       convergent
    ## 16: 1.648539   1.597690  0.05084856 0.3966149       0.1753014       convergent
    ## 17: 1.450018   1.640478 -0.19045985 0.3966149       0.1753014        divergent
    ## 18: 1.336301   1.431405 -0.09510492 0.3966149       0.1753014       convergent
    ## 19: 1.376318   1.602714 -0.22639549 0.3966149       0.1753014        divergent
    ## 20: 2.224850   1.621940  0.60291002 0.3966149       0.1753014        divergent
    ## 21: 1.795029   1.819726 -0.02469702 0.3966149       0.1753014       convergent
    ## 22: 2.113272   1.404026  0.70924642 0.3966149       0.1753014        divergent
    ## 23: 2.134081   1.359135  0.77494577 0.3966149       0.1753014        divergent
    ## 24: 1.475186   1.491053 -0.01586671 0.3966149       0.1753014       convergent
    ## 25: 1.793925   1.495717  0.29820775 0.3966149       0.1753014       convergent
    ## 26: 1.489802   1.456460  0.03334152 0.3966149       0.1753014       convergent
    ## 27: 1.376318   1.602714 -0.22639549 0.3966149       0.1753014        divergent
    ## 28: 1.603133   1.449486  0.15364703 0.3966149       0.1753014       convergent
    ## 29: 1.981162   1.538347  0.44281445 0.3966149       0.1753014       convergent
    ## 30: 1.549044   1.387951  0.16109275 0.3966149       0.1753014       convergent
    ## 31: 2.212238   1.331642  0.88059543 0.3966149       0.1753014        divergent
    ## 32: 1.759605   1.434285  0.32532045 0.3966149       0.1753014       convergent
    ## 33: 1.450018   1.640478 -0.19045985 0.3966149       0.1753014        divergent
    ## 34: 1.648539   1.597690  0.05084856 0.3966149       0.1753014       convergent
    ## 35: 2.029486   1.471753  0.55773251 0.3966149       0.1753014       convergent
    ## 36: 1.559953   1.384651  0.17530136 0.3966149       0.1753014       convergent
    ## 37: 2.344443   1.385743  0.95870057 0.3966149       0.1753014        divergent
    ## 38: 1.336301   1.431405 -0.09510492 0.3966149       0.1753014       convergent
    ##     rare_NES common_NES    NES_diff  MAD_filt median_NES_diff converge_diverge

``` r
p <- ggplot(classify_all_paths_overlap_filt, aes(x = NES, y = reorder(pathway, NES), fill = variant)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("common" = "lightblue", "rare" = "#F4D03F"))+
  facet_wrap(~ converge_diverge)+
  labs(
    title = "Convergence and Divergence in Enrichment Across Variant Type in Filtered Dataset",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    fill = "Variant Type"
  ) +
  theme_minimal()

p
```

![](fgsea_files/figure-gfm/visualize-summary-1.png)<!-- -->

``` r
#ggsave("converge_diverge_enrichment.png", plot = p2, path = "../figures", width = 14, height = 10, dpi = 300)
```

``` r
#save results
path_results <- classify_paths_filt %>% dplyr::select(pathway, converge_diverge)
#write.csv(path_results, "../results/path_results.csv")
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
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] patchwork_1.3.0      msigdbr_7.5.1        fgsea_1.32.4        
    ##  [4] org.Hs.eg.db_3.20.0  AnnotationDbi_1.68.0 IRanges_2.40.1      
    ##  [7] S4Vectors_0.44.0     Biobase_2.66.0       BiocGenerics_0.52.0 
    ## [10] lubridate_1.9.4      forcats_1.0.0        stringr_1.5.1       
    ## [13] dplyr_1.1.4          purrr_1.0.4          readr_2.1.5         
    ## [16] tidyr_1.3.1          tibble_3.2.1         ggplot2_3.5.1       
    ## [19] tidyverse_2.0.0     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] fastmatch_1.1-6         KEGGREST_1.46.0         babelgene_22.9         
    ##  [4] gtable_0.3.6            xfun_0.54               lattice_0.22-7         
    ##  [7] tzdb_0.4.0              vctrs_0.6.5             tools_4.4.3            
    ## [10] generics_0.1.4          parallel_4.4.3          RSQLite_2.3.9          
    ## [13] blob_1.2.4              pkgconfig_2.0.3         Matrix_1.7-4           
    ## [16] data.table_1.17.8       lifecycle_1.0.4         GenomeInfoDbData_1.2.13
    ## [19] farver_2.1.2            compiler_4.4.3          Biostrings_2.74.1      
    ## [22] munsell_0.5.1           codetools_0.2-20        GenomeInfoDb_1.42.3    
    ## [25] htmltools_0.5.8.1       yaml_2.3.11             pillar_1.10.1          
    ## [28] crayon_1.5.3            BiocParallel_1.40.2     cachem_1.1.0           
    ## [31] tidyselect_1.2.1        digest_0.6.39           stringi_1.8.7          
    ## [34] labeling_0.4.3          cowplot_1.1.3           fastmap_1.2.0          
    ## [37] grid_4.4.3              colorspace_2.1-2        cli_3.6.5              
    ## [40] magrittr_2.0.4          utf8_1.2.6              withr_3.0.2            
    ## [43] scales_1.3.0            UCSC.utils_1.2.0        bit64_4.6.0-1          
    ## [46] timechange_0.3.0        rmarkdown_2.29          XVector_0.46.0         
    ## [49] httr_1.4.7              bit_4.6.0               png_0.1-8              
    ## [52] hms_1.1.3               memoise_2.0.1           evaluate_1.0.5         
    ## [55] knitr_1.49              rlang_1.1.6             Rcpp_1.1.0             
    ## [58] glue_1.8.0              DBI_1.2.3               vroom_1.6.5            
    ## [61] rstudioapi_0.17.1       jsonlite_2.0.0          R6_2.6.1               
    ## [64] zlibbioc_1.52.0
