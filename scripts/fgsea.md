Rare Variant Preparation and FGSEA
================
Amalya Murrill
2026-04-18

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
    ##  [6] "GOBP_NEURON_FATE_SPECIFICATION"                           
    ##  [7] "GOBP_NEURON_PROJECTION_EXTENSION"                         
    ##  [8] "GOBP_NEURON_PROJECTION_ORGANIZATION"                      
    ##  [9] "GOBP_REGULATION_OF_CHROMATIN_ORGANIZATION"                
    ## [10] "GOBP_REGULATION_OF_LONG_TERM_NEURONAL_SYNAPTIC_PLASTICITY"
    ## [11] "GOBP_REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY"         
    ## [12] "GOBP_SYNAPSE_ASSEMBLY"

``` r
write.csv(as.data.frame(rare_collapsed$parentPathways), "../results/rare_collapsed_pathway_mapping.csv") #for summary
```

``` r
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

    ## [1] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION"                     
    ## [2] "GOBP_CHROMATIN_REMODELING"                                              
    ## [3] "GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE"
    ## [4] "GOBP_NEUROTRANSMITTER_SECRETION"                                        
    ## [5] "GOBP_SYNAPSE_ORGANIZATION"                                              
    ## [6] "GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE"                             
    ## [7] "KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM"                        
    ## [8] "KEGG_MTOR_SIGNALING_PATHWAY"

``` r
write.csv(as.data.frame(common_collapsed$parentPathways), "../results/common_collapsed_pathway_mapping.csv") # for summary
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
    ## 1:    rare                                 GOBP_NEURON_APOPTOTIC_PROCESS
    ## 2:    rare                                     GOBP_POSTSYNAPSE_ASSEMBLY
    ## 3:    rare                         GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE
    ## 4:    rare GOBP_POSITIVE_REGULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL
    ## 5:    rare                      GOBP_NEGATIVE_REGULATION_OF_AXONOGENESIS
    ## 6:    rare                 GOBP_POSTSYNAPTIC_SPECIALIZATION_ORGANIZATION
    ##            pval       padj   log2err        ES      NES  size  leadingEdge
    ##           <num>      <num>     <num>     <num>    <num> <int>       <list>
    ## 1: 0.0003713898 0.03583278 0.4984931 0.9695995 2.274621   218 6326, 88....
    ## 2: 0.0006731776 0.03583278 0.4772708 0.9237861 2.334302    26 5728, 85....
    ## 3: 0.0010400696 0.03583278 0.4550599 0.9924235 2.491460    19  8831, 22941
    ## 4: 0.0010768852 0.03583278 0.4550599 0.9160108 2.314084    23 5728, 85....
    ## 5: 0.0015370831 0.03583278 0.4550599 0.9848314 2.426876    56   8831, 5728
    ## 6: 0.0020534137 0.03583278 0.4317077 0.9876555 2.498870    32 8831, 57....

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
    ## 1:                      GOBP_POSITIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT
    ## 2:                                               GOBP_EXCITATORY_SYNAPSE_ASSEMBLY
    ## 3:             GOBP_REGULATION_OF_POSTSYNAPTIC_NEUROTRANSMITTER_RECEPTOR_ACTIVITY
    ## 4:                                 GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_BALANCE
    ## 5:                                  GOBP_REGULATION_OF_SYNAPTIC_VESICLE_RECYCLING
    ## 6: GOBP_ADENYLATE_CYCLASE_INHIBITING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY
    ##         pval     padj    log2err         ES        NES  size  leadingEdge
    ##        <num>    <num>      <num>      <num>      <num> <int>       <list>
    ## 1: 0.9867110 0.996337 0.01603044  0.1600589  0.6605510   139 51517, 1....
    ## 2: 0.9873950 0.996337 0.02956455  0.1747264  0.5243340    25 64101, 7....
    ## 3: 0.9895105 0.996337 0.07271411 -0.1793849 -0.5899148    17 254263, ....
    ## 4: 0.9910314 0.996337 0.08578444 -0.1427050 -0.6179620    49 7401, 64....
    ## 5: 0.9958506 0.996337 0.02850386  0.1610146  0.4586781    19 6622, 60....
    ## 6: 0.9963370 0.996337 0.02172400  0.1533360  0.5784461    75 2550, 29....

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
    ##   [4] "GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE"                                    
    ##   [5] "KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM"                               
    ##   [6] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION"                            
    ##   [7] "GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY"                                        
    ##   [8] "GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS"                                      
    ##   [9] "GOBP_NEUROTRANSMITTER_SECRETION"                                               
    ##  [10] "GOBP_REGULATION_OF_NEUROGENESIS"                                               
    ##  [11] "GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH"                                      
    ##  [12] "GOBP_NEURON_PROJECTION_GUIDANCE"                                               
    ##  [13] "GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS"                                    
    ##  [14] "GOBP_NEURON_DEATH"                                                             
    ##  [15] "KEGG_MTOR_SIGNALING_PATHWAY"                                                   
    ##  [16] "GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE"          
    ##  [17] "GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS"                                              
    ##  [18] "GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                     
    ##  [19] "KEGG_AXON_GUIDANCE"                                                            
    ##  [20] "GOBP_FOREBRAIN_GENERATION_OF_NEURONS"                                          
    ##  [21] "GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS"                    
    ##  [22] "GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS"                                      
    ##  [23] "GOBP_NEUROTRANSMITTER_TRANSPORT"                                               
    ##  [24] "GOBP_SYNAPTIC_VESICLE_PRIMING"                                                 
    ##  [25] "GOBP_REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY"                              
    ##  [26] "GOBP_LONG_TERM_SYNAPTIC_DEPRESSION"                                            
    ##  [27] "GOBP_ANTEROGRADE_AXONAL_TRANSPORT"                                             
    ##  [28] "GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE"                                         
    ##  [29] "GOBP_REGULATION_OF_PROTEIN_LOCALIZATION_TO_SYNAPSE"                            
    ##  [30] "GOBP_REGULATION_OF_AXONOGENESIS"                                               
    ##  [31] "GOBP_AXONAL_TRANSPORT"                                                         
    ##  [32] "GOBP_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT"                                 
    ##  [33] "GOBP_NEGATIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION"                             
    ##  [34] "GOBP_RETINAL_GANGLION_CELL_AXON_GUIDANCE"                                      
    ##  [35] "GOBP_NEURON_FATE_COMMITMENT"                                                   
    ##  [36] "GOBP_GLUTAMATE_SECRETION"                                                      
    ##  [37] "GOBP_NEURON_PROJECTION_ORGANIZATION"                                           
    ##  [38] "GOBP_NEGATIVE_REGULATION_OF_NEURON_DIFFERENTIATION"                            
    ##  [39] "GOBP_NEURON_PROJECTION_EXTENSION"                                              
    ##  [40] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_AXONOGENESIS"                               
    ##  [41] "GOBP_SYNAPSE_MATURATION"                                                       
    ##  [42] "GOBP_NEUROBLAST_PROLIFERATION"                                                 
    ##  [43] "GOBP_PRESYNAPTIC_ENDOCYTOSIS"                                                  
    ##  [44] "GOBP_NEGATIVE_REGULATION_OF_AXONOGENESIS"                                      
    ##  [45] "GOBP_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                              
    ##  [46] "GOBP_POSITIVE_REGULATION_OF_NEURON_APOPTOTIC_PROCESS"                          
    ##  [47] "GOBP_REGULATION_OF_GLUTAMATE_SECRETION"                                        
    ##  [48] "GOBP_POSTSYNAPSE_ORGANIZATION"                                                 
    ##  [49] "GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION"                                    
    ##  [50] "GOBP_AXON_DEVELOPMENT"                                                         
    ##  [51] "GOBP_NEURON_DEATH_IN_RESPONSE_TO_OXIDATIVE_STRESS"                             
    ##  [52] "GOBP_DNA_REPLICATION_DEPENDENT_CHROMATIN_ORGANIZATION"                         
    ##  [53] "GOBP_AXONAL_TRANSPORT_OF_MITOCHONDRION"                                        
    ##  [54] "GOBP_NEURON_APOPTOTIC_PROCESS"                                                 
    ##  [55] "GOBP_RETROGRADE_AXONAL_TRANSPORT"                                              
    ##  [56] "GOBP_FOREBRAIN_NEURON_DIFFERENTIATION"                                         
    ##  [57] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DEVELOPMENT"                                
    ##  [58] "GOBP_MOTOR_NEURON_AXON_GUIDANCE"                                               
    ##  [59] "GOBP_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_POTENTIAL"                            
    ##  [60] "GOBP_SYNAPSE_ASSEMBLY"                                                         
    ##  [61] "GOBP_NEUROEPITHELIAL_CELL_DIFFERENTIATION"                                     
    ##  [62] "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION"                                  
    ##  [63] "GOBP_PROTEIN_LOCALIZATION_TO_SYNAPSE"                                          
    ##  [64] "GOBP_REGULATION_OF_NEURONAL_SYNAPTIC_PLASTICITY"                               
    ##  [65] "GOBP_SYNAPTIC_VESICLE_LOCALIZATION"                                            
    ##  [66] "GOBP_SYNAPTIC_VESICLE_TRANSPORT"                                               
    ##  [67] "GOBP_REGULATION_OF_NEURON_DIFFERENTIATION"                                     
    ##  [68] "GOBP_AXON_EXTENSION"                                                           
    ##  [69] "GOBP_POSITIVE_REGULATION_OF_NEUROINFLAMMATORY_RESPONSE"                        
    ##  [70] "GOBP_POSITIVE_REGULATION_OF_CHROMATIN_ORGANIZATION"                            
    ##  [71] "GOBP_POSITIVE_REGULATION_OF_NEUROBLAST_PROLIFERATION"                          
    ##  [72] "GOBP_GLUTAMATE_METABOLIC_PROCESS"                                              
    ##  [73] "GOBP_POSITIVE_REGULATION_OF_NEUROTRANSMITTER_SECRETION"                        
    ##  [74] "GOBP_POSTSYNAPTIC_MEMBRANE_ORGANIZATION"                                       
    ##  [75] "GOBP_SYNAPTIC_MEMBRANE_ADHESION"                                               
    ##  [76] "GOBP_POSITIVE_REGULATION_OF_NEURON_DIFFERENTIATION"                            
    ##  [77] "GOBP_REGULATION_OF_NEUROBLAST_PROLIFERATION"                                   
    ##  [78] "GOBP_REGULATION_OF_SYNAPTIC_PLASTICITY"                                        
    ##  [79] "GOBP_POSITIVE_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT"                        
    ##  [80] "GOBP_POSITIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION"                             
    ##  [81] "GOBP_AXONEMAL_DYNEIN_COMPLEX_ASSEMBLY"                                         
    ##  [82] "GOBP_NEURON_RECOGNITION"                                                       
    ##  [83] "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY"                                           
    ##  [84] "GOBP_LONG_TERM_SYNAPTIC_POTENTIATION"                                          
    ##  [85] "GOBP_SYNAPTIC_VESICLE_MEMBRANE_ORGANIZATION"                                   
    ##  [86] "GOBP_SPINAL_CORD_MOTOR_NEURON_DIFFERENTIATION"                                 
    ##  [87] "GOBP_NEUROPEPTIDE_SIGNALING_PATHWAY"                                           
    ##  [88] "GOBP_SYNAPTIC_VESICLE_RECYCLING"                                               
    ##  [89] "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING"                                   
    ##  [90] "GOBP_NEURONAL_STEM_CELL_POPULATION_MAINTENANCE"                                
    ##  [91] "GOBP_ENSHEATHMENT_OF_NEURONS"                                                  
    ##  [92] "GOBP_NEGATIVE_REGULATION_OF_NEURON_DEATH"                                      
    ##  [93] "GOBP_REGULATION_OF_LONG_TERM_NEURONAL_SYNAPTIC_PLASTICITY"                     
    ##  [94] "GOBP_NEURON_MIGRATION"                                                         
    ##  [95] "GOBP_HETEROCHROMATIN_ORGANIZATION"                                             
    ##  [96] "GOBP_NEUROMUSCULAR_JUNCTION_DEVELOPMENT"                                       
    ##  [97] "GOBP_NEUROTROPHIN_TRK_RECEPTOR_SIGNALING_PATHWAY"                              
    ##  [98] "GOBP_MODULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL"                          
    ##  [99] "GOBP_REGULATION_OF_LONG_TERM_SYNAPTIC_POTENTIATION"                            
    ## [100] "GOBP_NEURON_CELLULAR_HOMEOSTASIS"                                              
    ## [101] "GOBP_RESPONSE_TO_LEUKEMIA_INHIBITORY_FACTOR"                                   
    ## [102] "GOBP_POSITIVE_REGULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL"                 
    ## [103] "GOBP_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY"                                     
    ## [104] "GOBP_REGULATION_OF_SYNAPTIC_VESICLE_EXOCYTOSIS"                                
    ## [105] "GOBP_REGULATION_OF_NEUROTRANSMITTER_UPTAKE"                                    
    ## [106] "GOBP_PROTEIN_LOCALIZATION_TO_POSTSYNAPSE"                                      
    ## [107] "GOBP_POSTSYNAPTIC_SPECIALIZATION_ORGANIZATION"                                 
    ## [108] "GOBP_SYNAPTIC_TRANSMISSION_CHOLINERGIC"                                        
    ## [109] "GOBP_AXON_ENSHEATHMENT_IN_CENTRAL_NERVOUS_SYSTEM"                              
    ## [110] "GOBP_CEREBRAL_CORTEX_NEURON_DIFFERENTIATION"                                   
    ## [111] "GOBP_NEURON_FATE_SPECIFICATION"                                                
    ## [112] "GOBP_REGULATION_OF_RECEPTOR_LOCALIZATION_TO_SYNAPSE"                           
    ## [113] "GOBP_REGULATION_OF_POSTSYNAPSE_ORGANIZATION"                                   
    ## [114] "GOBP_G_PROTEIN_COUPLED_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY"                   
    ## [115] "GOBP_MODIFICATION_OF_SYNAPTIC_STRUCTURE"                                       
    ## [116] "GOBP_DNA_REPLICATION_INDEPENDENT_CHROMATIN_ORGANIZATION"                       
    ## [117] "GOBP_SYNAPTIC_VESICLE_CYTOSKELETAL_TRANSPORT"                                  
    ## [118] "GOBP_NEURON_MATURATION"                                                        
    ## [119] "GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_POSTURE"                                
    ## [120] "GOBP_REGULATION_OF_SYNAPSE_ASSEMBLY"                                           
    ## [121] "GOBP_DOPAMINERGIC_NEURON_DIFFERENTIATION"                                      
    ## [122] "GOBP_PROTEIN_LOCALIZATION_TO_CHROMATIN"                                        
    ## [123] "GOBP_NEURON_PROJECTION_ARBORIZATION"                                           
    ## [124] "GOBP_NEUROTRANSMITTER_REUPTAKE"                                                
    ## [125] "GOBP_NEGATIVE_REGULATION_OF_NEURON_APOPTOTIC_PROCESS"                          
    ## [126] "GOBP_NEUROMUSCULAR_SYNAPTIC_TRANSMISSION"                                      
    ## [127] "GOBP_INHIBITORY_SYNAPSE_ASSEMBLY"                                              
    ## [128] "GOBP_NEUROINFLAMMATORY_RESPONSE"                                               
    ## [129] "GOBP_PRESYNAPSE_ORGANIZATION"                                                  
    ## [130] "GOBP_AXONAL_FASCICULATION"                                                     
    ## [131] "GOBP_REGULATION_OF_PRESYNAPSE_ORGANIZATION"                                    
    ## [132] "GOBP_NEUROTROPHIN_SIGNALING_PATHWAY"                                           
    ## [133] "GOBP_NEUROMUSCULAR_PROCESS"                                                    
    ## [134] "GOBP_NEURONAL_ACTION_POTENTIAL"                                                
    ## [135] "GOBP_FOREBRAIN_NEURON_DEVELOPMENT"                                             
    ## [136] "GOBP_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_NEUROTRANSMITTER_RECEPTOR_LEVELS"     
    ## [137] "GOBP_POSITIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"               
    ## [138] "GOBP_CHROMATIN_DISASSEMBLY"                                                    
    ## [139] "GOBP_FACULTATIVE_HETEROCHROMATIN_ASSEMBLY"                                     
    ## [140] "GOBP_DNA_METHYLATION_DEPENDENT_HETEROCHROMATIN_ASSEMBLY"                       
    ## [141] "GOBP_REGULATION_OF_NEUROTRANSMITTER_RECEPTOR_ACTIVITY"                         
    ## [142] "GOBP_REGULATION_OF_DNA_METHYLATION_DEPENDENT_HETEROCHROMATIN_ASSEMBLY"         
    ## [143] "GOBP_POSTSYNAPTIC_SIGNAL_TRANSDUCTION"                                         
    ## [144] "GOBP_POSITIVE_REGULATION_OF_LONG_TERM_SYNAPTIC_POTENTIATION"                   
    ## [145] "GOBP_CHEMICAL_SYNAPTIC_TRANSMISSION_POSTSYNAPTIC"                              
    ## [146] "GOBP_NEGATIVE_REGULATION_OF_OXIDATIVE_STRESS_INDUCED_NEURON_DEATH"             
    ## [147] "GOBP_REGULATION_OF_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY"                          
    ## [148] "GOBP_NEURON_PROJECTION_REGENERATION"                                           
    ## [149] "GOBP_REGULATION_OF_CHROMATIN_BINDING"                                          
    ## [150] "GOBP_REGULATION_OF_NEURON_PROJECTION_REGENERATION"                             
    ## [151] "GOBP_MIDBRAIN_DOPAMINERGIC_NEURON_DIFFERENTIATION"                             
    ## [152] "GOBP_REGULATION_OF_CHROMATIN_ASSEMBLY"                                         
    ## [153] "GOBP_REGULATION_OF_CHROMATIN_ORGANIZATION"                                     
    ## [154] "GOBP_AXONEME_ASSEMBLY"                                                         
    ## [155] "GOBP_RECEPTOR_LOCALIZATION_TO_SYNAPSE"                                         
    ## [156] "GOBP_NEUROTRANSMITTER_METABOLIC_PROCESS"                                       
    ## [157] "GOBP_SYNAPTONEMAL_COMPLEX_ORGANIZATION"                                        
    ## [158] "GOBP_SYNAPTIC_TRANSMISSION_GABAERGIC"                                          
    ## [159] "GOBP_POSITIVE_REGULATION_OF_NEURON_MIGRATION"                                  
    ## [160] "GOBP_NEUROTRANSMITTER_UPTAKE"                                                  
    ## [161] "GOBP_SYNAPTIC_TRANSMISSION_DOPAMINERGIC"                                       
    ## [162] "GOBP_RESPONSE_TO_AXON_INJURY"                                                  
    ## [163] "GOBP_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"                                      
    ## [164] "GOBP_RNA_MEDIATED_GENE_SILENCING_BY_INHIBITION_OF_TRANSLATION"                 
    ## [165] "GOBP_POSITIVE_REGULATION_OF_AXON_EXTENSION"                                    
    ## [166] "GOBP_POSTSYNAPSE_ASSEMBLY"                                                     
    ## [167] "GOBP_CALCIUM_ION_REGULATED_EXOCYTOSIS_OF_NEUROTRANSMITTER"                     
    ## [168] "GOBP_POSITIVE_REGULATION_OF_SYNAPSE_ASSEMBLY"                                  
    ## [169] "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GABAERGIC"                            
    ## [170] "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"                        
    ## [171] "GOBP_REGULATION_OF_NEURON_MIGRATION"                                           
    ## [172] "GOBP_NEUROTRANSMITTER_RECEPTOR_TRANSPORT"                                      
    ## [173] "GOBP_L_GLUTAMATE_IMPORT_ACROSS_PLASMA_MEMBRANE"                                
    ## [174] "GOBP_NEUROTRANSMITTER_RECEPTOR_INTERNALIZATION"                                
    ## [175] "GOBP_L_GLUTAMATE_TRANSMEMBRANE_TRANSPORT"                                      
    ## [176] "GOBP_POSTSYNAPTIC_SPECIALIZATION_ASSEMBLY"                                     
    ## [177] "GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_REGENERATION"                    
    ## [178] "GOBP_POSITIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                     
    ## [179] "GOBP_EXCITATORY_SYNAPSE_ASSEMBLY"                                              
    ## [180] "GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_BALANCE"                                
    ## [181] "GOBP_REGULATION_OF_SYNAPTIC_VESICLE_RECYCLING"                                 
    ## [182] "GOBP_ADENYLATE_CYCLASE_INHIBITING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY"

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
    ##   1:                                                  GOBP_NEURON_APOPTOTIC_PROCESS
    ##   2:                                                      GOBP_POSTSYNAPSE_ASSEMBLY
    ##   3:                                          GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE
    ##   4:                  GOBP_POSITIVE_REGULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL
    ##   5:                                       GOBP_NEGATIVE_REGULATION_OF_AXONOGENESIS
    ##  ---                                                                               
    ## 360:                      GOBP_POSITIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT
    ## 361:                                               GOBP_EXCITATORY_SYNAPSE_ASSEMBLY
    ## 362:                                 GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_BALANCE
    ## 363:                                  GOBP_REGULATION_OF_SYNAPTIC_VESICLE_RECYCLING
    ## 364: GOBP_ADENYLATE_CYCLASE_INHIBITING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY
    ##              pval       padj    log2err         ES        NES  size
    ##             <num>      <num>      <num>      <num>      <num> <int>
    ##   1: 0.0003713898 0.03583278 0.49849311  0.9695995  2.2746212   218
    ##   2: 0.0006731776 0.03583278 0.47727082  0.9237861  2.3343020    26
    ##   3: 0.0010400696 0.03583278 0.45505987  0.9924235  2.4914598    19
    ##   4: 0.0010768852 0.03583278 0.45505987  0.9160108  2.3140843    23
    ##   5: 0.0015370831 0.03583278 0.45505987  0.9848314  2.4268765    56
    ##  ---                                                               
    ## 360: 0.9867109635 0.99633700 0.01603044  0.1600589  0.6605510   139
    ## 361: 0.9873949580 0.99633700 0.02956455  0.1747264  0.5243340    25
    ## 362: 0.9910313901 0.99633700 0.08578444 -0.1427050 -0.6179620    49
    ## 363: 0.9958506224 0.99633700 0.02850386  0.1610146  0.4586781    19
    ## 364: 0.9963369963 0.99633700 0.02172400  0.1533360  0.5784461    75
    ##       leadingEdge
    ##            <list>
    ##   1: 6326, 88....
    ##   2: 5728, 85....
    ##   3:  8831, 22941
    ##   4: 5728, 85....
    ##   5:   8831, 5728
    ##  ---             
    ## 360: 51517, 1....
    ## 361: 64101, 7....
    ## 362: 7401, 64....
    ## 363: 6622, 60....
    ## 364: 2550, 29....

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
    ##  1 GOBP_NEURON_APOPTOTIC_PROCESS           2.27      1.19     0.0358       0.438
    ##  2 GOBP_POSTSYNAPSE_ASSEMBLY               2.33      0.698    0.0358       0.977
    ##  3 GOBP_MAINTENANCE_OF_SYNAPSE_STRUCT…     2.49      1.49     0.0358       0.345
    ##  4 GOBP_POSITIVE_REGULATION_OF_EXCITA…     2.31     -1.06     0.0358       0.636
    ##  5 GOBP_NEGATIVE_REGULATION_OF_AXONOG…     2.43      1.27     0.0358       0.433
    ##  6 GOBP_POSTSYNAPTIC_SPECIALIZATION_O…     2.50      1.04     0.0358       0.695
    ##  7 GOBP_REGULATION_OF_LONG_TERM_NEURO…     2.50      1.15     0.0358       0.570
    ##  8 GOBP_LONG_TERM_SYNAPTIC_POTENTIATI…     2.11      1.13     0.0358       0.537
    ##  9 GOBP_NEURON_DEATH                       2.20      1.33     0.0358       0.174
    ## 10 GOBP_NEURONAL_ACTION_POTENTIAL          2.49      0.908    0.0358       0.837
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
    ##  1 GOBP_NEURON_AP…     2.27      1.19     0.0358       0.438    1.09       0.563
    ##  2 GOBP_POSTSYNAP…     2.33      0.698    0.0358       0.977    1.64       0.563
    ##  3 GOBP_MAINTENAN…     2.49      1.49     0.0358       0.345    1.00       0.563
    ##  4 GOBP_POSITIVE_…     2.31     -1.06     0.0358       0.636    3.37       0.563
    ##  5 GOBP_NEGATIVE_…     2.43      1.27     0.0358       0.433    1.15       0.563
    ##  6 GOBP_POSTSYNAP…     2.50      1.04     0.0358       0.695    1.45       0.563
    ##  7 GOBP_REGULATIO…     2.50      1.15     0.0358       0.570    1.36       0.563
    ##  8 GOBP_LONG_TERM…     2.11      1.13     0.0358       0.537    0.979      0.563
    ##  9 GOBP_NEURON_DE…     2.20      1.33     0.0358       0.174    0.878      0.563
    ## 10 GOBP_NEURONAL_…     2.49      0.908    0.0358       0.837    1.58       0.563
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
    ##   pathway                                                       NES_diff
    ##   <chr>                                                            <dbl>
    ## 1 GOBP_POSITIVE_REGULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL     3.37

``` r
most_convergent
```

    ## # A tibble: 1 × 2
    ##   pathway                                    NES_diff
    ##   <chr>                                         <dbl>
    ## 1 GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION 0.000771

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
    ##  1 GOBP_NEURON_AP…     2.27      1.19     0.0358       0.438    1.09       0.563
    ##  2 GOBP_POSTSYNAP…     2.33      0.698    0.0358       0.977    1.64       0.563
    ##  3 GOBP_MAINTENAN…     2.49      1.49     0.0358       0.345    1.00       0.563
    ##  4 GOBP_POSITIVE_…     2.31     -1.06     0.0358       0.636    3.37       0.563
    ##  5 GOBP_NEGATIVE_…     2.43      1.27     0.0358       0.433    1.15       0.563
    ##  6 GOBP_POSTSYNAP…     2.50      1.04     0.0358       0.695    1.45       0.563
    ##  7 GOBP_REGULATIO…     2.50      1.15     0.0358       0.570    1.36       0.563
    ##  8 GOBP_LONG_TERM…     2.11      1.13     0.0358       0.537    0.979      0.563
    ##  9 GOBP_NEURON_DE…     2.20      1.33     0.0358       0.174    0.878      0.563
    ## 10 GOBP_NEURONAL_…     2.49      0.908    0.0358       0.837    1.58       0.563
    ## # ℹ 172 more rows
    ## # ℹ 2 more variables: median_NES_diff <dbl>, converge_diverge <chr>

``` r
divergent_unfilt <- classify_paths %>% filter(converge_diverge == "divergent")
divergent_unfilt
```

    ## # A tibble: 61 × 9
    ##    pathway         rare_NES common_NES rare_padj common_padj NES_diff MAD_unfilt
    ##    <chr>              <dbl>      <dbl>     <dbl>       <dbl>    <dbl>      <dbl>
    ##  1 GOBP_POSTSYNAP…     2.33      0.698    0.0358       0.977     1.64      0.563
    ##  2 GOBP_POSITIVE_…     2.31     -1.06     0.0358       0.636     3.37      0.563
    ##  3 GOBP_POSTSYNAP…     2.50      1.04     0.0358       0.695     1.45      0.563
    ##  4 GOBP_REGULATIO…     2.50      1.15     0.0358       0.570     1.36      0.563
    ##  5 GOBP_NEURONAL_…     2.49      0.908    0.0358       0.837     1.58      0.563
    ##  6 GOBP_POSTSYNAP…     2.35      0.523    0.0358       0.996     1.83      0.563
    ##  7 GOBP_POSITIVE_…     2.28      0.731    0.0358       0.965     1.55      0.563
    ##  8 GOBP_EXCITATOR…     2.35      0.524    0.0358       0.996     1.83      0.563
    ##  9 GOBP_REGULATIO…     2.27      0.642    0.0388       0.996     1.63      0.563
    ## 10 GOBP_MODULATIO…     2.26     -1.08     0.0410       0.628     3.33      0.563
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
    ##  1 GOBP_NEURON_AP…     2.27      1.19     0.0358       0.438    1.09       0.563
    ##  2 GOBP_MAINTENAN…     2.49      1.49     0.0358       0.345    1.00       0.563
    ##  3 GOBP_NEGATIVE_…     2.43      1.27     0.0358       0.433    1.15       0.563
    ##  4 GOBP_LONG_TERM…     2.11      1.13     0.0358       0.537    0.979      0.563
    ##  5 GOBP_NEURON_DE…     2.20      1.33     0.0358       0.174    0.878      0.563
    ##  6 GOBP_POSITIVE_…     2.00      1.44     0.0358       0.345    0.558      0.563
    ##  7 GOBP_REGULATIO…     2.48      1.23     0.0410       0.463    1.25       0.563
    ##  8 GOBP_REGULATIO…     1.95      1.19     0.0680       0.463    0.759      0.563
    ##  9 GOBP_POSITIVE_…     1.95      1.52     0.0680       0.174    0.429      0.563
    ## 10 GOBP_CHEMICAL_…     2.01      0.872    0.0774       0.899    1.14       0.563
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

    ## [1] 30

``` r
setdiff(rare_sig, common_sig)  #rare only
```

    ##   [1] "GOBP_NEURON_APOPTOTIC_PROCESS"                                                 
    ##   [2] "GOBP_POSTSYNAPSE_ASSEMBLY"                                                     
    ##   [3] "GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE"                                         
    ##   [4] "GOBP_POSITIVE_REGULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL"                 
    ##   [5] "GOBP_NEGATIVE_REGULATION_OF_AXONOGENESIS"                                      
    ##   [6] "GOBP_POSTSYNAPTIC_SPECIALIZATION_ORGANIZATION"                                 
    ##   [7] "GOBP_REGULATION_OF_LONG_TERM_NEURONAL_SYNAPTIC_PLASTICITY"                     
    ##   [8] "GOBP_LONG_TERM_SYNAPTIC_POTENTIATION"                                          
    ##   [9] "GOBP_NEURONAL_ACTION_POTENTIAL"                                                
    ##  [10] "GOBP_POSTSYNAPTIC_SPECIALIZATION_ASSEMBLY"                                     
    ##  [11] "GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS"                                      
    ##  [12] "GOBP_POSITIVE_REGULATION_OF_AXON_EXTENSION"                                    
    ##  [13] "GOBP_EXCITATORY_SYNAPSE_ASSEMBLY"                                              
    ##  [14] "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GABAERGIC"                            
    ##  [15] "GOBP_MODULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL"                          
    ##  [16] "GOBP_REGULATION_OF_NEURONAL_SYNAPTIC_PLASTICITY"                               
    ##  [17] "GOBP_SYNAPTIC_TRANSMISSION_GABAERGIC"                                          
    ##  [18] "GOBP_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_POTENTIAL"                            
    ##  [19] "GOBP_POSITIVE_REGULATION_OF_SYNAPSE_ASSEMBLY"                                  
    ##  [20] "GOBP_CHEMICAL_SYNAPTIC_TRANSMISSION_POSTSYNAPTIC"                              
    ##  [21] "GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_REGENERATION"                    
    ##  [22] "GOBP_REGULATION_OF_POSTSYNAPSE_ORGANIZATION"                                   
    ##  [23] "GOBP_POSTSYNAPTIC_MEMBRANE_ORGANIZATION"                                       
    ##  [24] "GOBP_SYNAPTIC_VESICLE_LOCALIZATION"                                            
    ##  [25] "GOBP_POSITIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"               
    ##  [26] "GOBP_NEGATIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION"                             
    ##  [27] "GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                     
    ##  [28] "GOBP_NEURON_PROJECTION_ORGANIZATION"                                           
    ##  [29] "GOBP_NEURON_FATE_SPECIFICATION"                                                
    ##  [30] "GOBP_SYNAPSE_MATURATION"                                                       
    ##  [31] "GOBP_FOREBRAIN_GENERATION_OF_NEURONS"                                          
    ##  [32] "GOBP_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY"                                     
    ##  [33] "GOBP_LONG_TERM_SYNAPTIC_DEPRESSION"                                            
    ##  [34] "GOBP_REGULATION_OF_SYNAPSE_ASSEMBLY"                                           
    ##  [35] "GOBP_NEGATIVE_REGULATION_OF_NEURON_APOPTOTIC_PROCESS"                          
    ##  [36] "GOBP_REGULATION_OF_CHROMATIN_ORGANIZATION"                                     
    ##  [37] "GOBP_FOREBRAIN_NEURON_DIFFERENTIATION"                                         
    ##  [38] "GOBP_AXON_ENSHEATHMENT_IN_CENTRAL_NERVOUS_SYSTEM"                              
    ##  [39] "GOBP_INHIBITORY_SYNAPSE_ASSEMBLY"                                              
    ##  [40] "GOBP_REGULATION_OF_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY"                          
    ##  [41] "GOBP_NEUROTRANSMITTER_UPTAKE"                                                  
    ##  [42] "GOBP_AXON_EXTENSION"                                                           
    ##  [43] "GOBP_DOPAMINERGIC_NEURON_DIFFERENTIATION"                                      
    ##  [44] "GOBP_REGULATION_OF_LONG_TERM_SYNAPTIC_POTENTIATION"                            
    ##  [45] "GOBP_REGULATION_OF_NEURON_PROJECTION_REGENERATION"                             
    ##  [46] "GOBP_NEGATIVE_REGULATION_OF_OXIDATIVE_STRESS_INDUCED_NEURON_DEATH"             
    ##  [47] "GOBP_POSITIVE_REGULATION_OF_LONG_TERM_SYNAPTIC_POTENTIATION"                   
    ##  [48] "GOBP_NEUROTRANSMITTER_REUPTAKE"                                                
    ##  [49] "GOBP_PRESYNAPSE_ORGANIZATION"                                                  
    ##  [50] "GOBP_ENSHEATHMENT_OF_NEURONS"                                                  
    ##  [51] "GOBP_NEUROMUSCULAR_PROCESS"                                                    
    ##  [52] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_AXONOGENESIS"                               
    ##  [53] "GOBP_POSITIVE_REGULATION_OF_NEURON_APOPTOTIC_PROCESS"                          
    ##  [54] "GOBP_REGULATION_OF_AXONOGENESIS"                                               
    ##  [55] "GOBP_RETROGRADE_AXONAL_TRANSPORT"                                              
    ##  [56] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DEVELOPMENT"                                
    ##  [57] "GOBP_NEURON_FATE_COMMITMENT"                                                   
    ##  [58] "GOBP_NEURON_PROJECTION_REGENERATION"                                           
    ##  [59] "GOBP_MIDBRAIN_DOPAMINERGIC_NEURON_DIFFERENTIATION"                             
    ##  [60] "GOBP_POSITIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                     
    ##  [61] "GOBP_POSITIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION"                             
    ##  [62] "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"                        
    ##  [63] "GOBP_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"                                      
    ##  [64] "GOBP_NEGATIVE_REGULATION_OF_NEURON_DEATH"                                      
    ##  [65] "GOBP_POSTSYNAPSE_ORGANIZATION"                                                 
    ##  [66] "GOBP_REGULATION_OF_NEUROTRANSMITTER_RECEPTOR_ACTIVITY"                         
    ##  [67] "GOBP_POSITIVE_REGULATION_OF_NEUROBLAST_PROLIFERATION"                          
    ##  [68] "GOBP_SYNAPSE_ASSEMBLY"                                                         
    ##  [69] "GOBP_SYNAPTIC_MEMBRANE_ADHESION"                                               
    ##  [70] "GOBP_NEURON_PROJECTION_EXTENSION"                                              
    ##  [71] "GOBP_REGULATION_OF_SYNAPTIC_PLASTICITY"                                        
    ##  [72] "GOBP_REGULATION_OF_NEUROBLAST_PROLIFERATION"                                   
    ##  [73] "GOBP_PROTEIN_LOCALIZATION_TO_SYNAPSE"                                          
    ##  [74] "GOBP_REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY"                              
    ##  [75] "GOBP_NEURON_DEATH_IN_RESPONSE_TO_OXIDATIVE_STRESS"                             
    ##  [76] "GOBP_AXONAL_TRANSPORT"                                                         
    ##  [77] "GOBP_FACULTATIVE_HETEROCHROMATIN_ASSEMBLY"                                     
    ##  [78] "GOBP_HETEROCHROMATIN_ORGANIZATION"                                             
    ##  [79] "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING"                                   
    ##  [80] "GOBP_NEUROTRANSMITTER_TRANSPORT"                                               
    ##  [81] "GOBP_RESPONSE_TO_AXON_INJURY"                                                  
    ##  [82] "KEGG_AXON_GUIDANCE"                                                            
    ##  [83] "GOBP_POSITIVE_REGULATION_OF_CHROMATIN_ORGANIZATION"                            
    ##  [84] "GOBP_REGULATION_OF_CHROMATIN_ASSEMBLY"                                         
    ##  [85] "GOBP_DNA_METHYLATION_DEPENDENT_HETEROCHROMATIN_ASSEMBLY"                       
    ##  [86] "GOBP_SYNAPTIC_VESICLE_MEMBRANE_ORGANIZATION"                                   
    ##  [87] "GOBP_NEURON_RECOGNITION"                                                       
    ##  [88] "GOBP_SYNAPTIC_VESICLE_TRANSPORT"                                               
    ##  [89] "GOBP_NEUROBLAST_PROLIFERATION"                                                 
    ##  [90] "GOBP_REGULATION_OF_NEURON_DIFFERENTIATION"                                     
    ##  [91] "GOBP_CALCIUM_ION_REGULATED_EXOCYTOSIS_OF_NEUROTRANSMITTER"                     
    ##  [92] "GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_POSTURE"                                
    ##  [93] "GOBP_SYNAPTIC_VESICLE_PRIMING"                                                 
    ##  [94] "GOBP_AXONAL_TRANSPORT_OF_MITOCHONDRION"                                        
    ##  [95] "GOBP_POSITIVE_REGULATION_OF_NEUROTRANSMITTER_SECRETION"                        
    ##  [96] "GOBP_NEURON_CELLULAR_HOMEOSTASIS"                                              
    ##  [97] "GOBP_GLUTAMATE_SECRETION"                                                      
    ##  [98] "GOBP_ANTEROGRADE_AXONAL_TRANSPORT"                                             
    ##  [99] "GOBP_CHROMATIN_DISASSEMBLY"                                                    
    ## [100] "GOBP_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                              
    ## [101] "GOBP_NEURON_MATURATION"                                                        
    ## [102] "GOBP_CEREBRAL_CORTEX_NEURON_DIFFERENTIATION"                                   
    ## [103] "GOBP_NEUROMUSCULAR_SYNAPTIC_TRANSMISSION"                                      
    ## [104] "GOBP_NEURON_MIGRATION"                                                         
    ## [105] "GOBP_AXON_DEVELOPMENT"                                                         
    ## [106] "GOBP_POSITIVE_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT"                        
    ## [107] "GOBP_RETINAL_GANGLION_CELL_AXON_GUIDANCE"                                      
    ## [108] "GOBP_REGULATION_OF_NEUROTRANSMITTER_UPTAKE"                                    
    ## [109] "GOBP_REGULATION_OF_GLUTAMATE_SECRETION"                                        
    ## [110] "GOBP_RNA_MEDIATED_GENE_SILENCING_BY_INHIBITION_OF_TRANSLATION"                 
    ## [111] "GOBP_POSITIVE_REGULATION_OF_NEURON_DIFFERENTIATION"                            
    ## [112] "GOBP_REGULATION_OF_NEURON_MIGRATION"                                           
    ## [113] "GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS"                    
    ## [114] "GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_BALANCE"                                
    ## [115] "GOBP_REGULATION_OF_PRESYNAPSE_ORGANIZATION"                                    
    ## [116] "GOBP_NEURON_PROJECTION_ARBORIZATION"                                           
    ## [117] "GOBP_POSITIVE_REGULATION_OF_NEURON_MIGRATION"                                  
    ## [118] "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION"                                  
    ## [119] "GOBP_RESPONSE_TO_LEUKEMIA_INHIBITORY_FACTOR"                                   
    ## [120] "GOBP_NEUROTROPHIN_SIGNALING_PATHWAY"                                           
    ## [121] "GOBP_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT"                                 
    ## [122] "GOBP_NEUROTRANSMITTER_RECEPTOR_INTERNALIZATION"                                
    ## [123] "GOBP_NEURONAL_STEM_CELL_POPULATION_MAINTENANCE"                                
    ## [124] "GOBP_REGULATION_OF_PROTEIN_LOCALIZATION_TO_SYNAPSE"                            
    ## [125] "GOBP_NEUROTRANSMITTER_RECEPTOR_TRANSPORT"                                      
    ## [126] "GOBP_FOREBRAIN_NEURON_DEVELOPMENT"                                             
    ## [127] "GOBP_RECEPTOR_LOCALIZATION_TO_SYNAPSE"                                         
    ## [128] "GOBP_SYNAPTIC_VESICLE_RECYCLING"                                               
    ## [129] "GOBP_REGULATION_OF_RECEPTOR_LOCALIZATION_TO_SYNAPSE"                           
    ## [130] "GOBP_NEGATIVE_REGULATION_OF_NEURON_DIFFERENTIATION"                            
    ## [131] "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY"                                           
    ## [132] "GOBP_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_NEUROTRANSMITTER_RECEPTOR_LEVELS"     
    ## [133] "GOBP_SYNAPTIC_VESICLE_CYTOSKELETAL_TRANSPORT"                                  
    ## [134] "GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS"                                              
    ## [135] "GOBP_AXONEMAL_DYNEIN_COMPLEX_ASSEMBLY"                                         
    ## [136] "GOBP_PRESYNAPTIC_ENDOCYTOSIS"                                                  
    ## [137] "GOBP_ADENYLATE_CYCLASE_INHIBITING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY"
    ## [138] "GOBP_REGULATION_OF_CHROMATIN_BINDING"                                          
    ## [139] "GOBP_NEUROTROPHIN_TRK_RECEPTOR_SIGNALING_PATHWAY"                              
    ## [140] "GOBP_AXONAL_FASCICULATION"

``` r
setdiff(common_sig, rare_sig)  #common only
```

    ## [1] "KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM"

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

    ## # A tibble: 15 × 3
    ##    pathway                                                   rare_NES common_NES
    ##    <chr>                                                        <dbl>      <dbl>
    ##  1 GOBP_NEURON_DEATH                                             2.20       1.33
    ##  2 GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH                      1.95       1.52
    ##  3 GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION            1.78       1.49
    ##  4 GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS                      1.76       1.42
    ##  5 GOBP_CHROMATIN_REMODELING                                     2.22       1.61
    ##  6 GOBP_SYNAPSE_ORGANIZATION                                     2.12       1.39
    ##  7 GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS                    1.56       1.37
    ##  8 GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE                    1.60       1.43
    ##  9 GOBP_NEURON_PROJECTION_GUIDANCE                               1.54       1.38
    ## 10 GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY                        1.49       1.44
    ## 11 GOBP_REGULATION_OF_NEUROGENESIS                               2.13       1.35
    ## 12 GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJ…     1.76       1.80
    ## 13 GOBP_NEUROTRANSMITTER_SECRETION                               1.46       1.48
    ## 14 GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_A…     1.43       1.59
    ## 15 KEGG_MTOR_SIGNALING_PATHWAY                                   1.34       1.60

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

    ## # A tibble: 15 × 6
    ##    pathway                 rare_NES common_NES NES_diff MAD_filt median_NES_diff
    ##    <chr>                      <dbl>      <dbl>    <dbl>    <dbl>           <dbl>
    ##  1 GOBP_NEURON_DEATH           2.20       1.33   0.878     0.335           0.189
    ##  2 GOBP_POSITIVE_REGULATI…     1.95       1.52   0.429     0.335           0.189
    ##  3 GOBP_CENTRAL_NERVOUS_S…     1.78       1.49   0.297     0.335           0.189
    ##  4 GOBP_POSITIVE_REGULATI…     1.76       1.42   0.337     0.335           0.189
    ##  5 GOBP_CHROMATIN_REMODEL…     2.22       1.61   0.607     0.335           0.189
    ##  6 GOBP_SYNAPSE_ORGANIZAT…     2.12       1.39   0.723     0.335           0.189
    ##  7 GOBP_REGULATION_OF_NEU…     1.56       1.37   0.189     0.335           0.189
    ##  8 GOBP_VESICLE_MEDIATED_…     1.60       1.43   0.166     0.335           0.189
    ##  9 GOBP_NEURON_PROJECTION…     1.54       1.38   0.161     0.335           0.189
    ## 10 GOBP_CHROMATIN_ASSEMBL…     1.49       1.44   0.0466    0.335           0.189
    ## 11 GOBP_REGULATION_OF_NEU…     2.13       1.35   0.777     0.335           0.189
    ## 12 GOBP_NEURON_PROJECTION…     1.76       1.80  -0.0373    0.335           0.189
    ## 13 GOBP_NEUROTRANSMITTER_…     1.46       1.48  -0.0203    0.335           0.189
    ## 14 GOBP_NEGATIVE_REGULATI…     1.43       1.59  -0.156     0.335           0.189
    ## 15 KEGG_MTOR_SIGNALING_PA…     1.34       1.60  -0.259     0.335           0.189

``` r
#unfiltered MA for robustness check: 0.5844998
mad_filt <- 0.2611755
median_filt <- 0.1677036
```

## Re-Classify

``` r
#threshold of 1 MAD away from the median for divergence:

#new median ± 1 MAD = 0.17 ± 0.26
#new convergence zone: -0.09 to 0.33

most_divergent <- variant_path_NES %>% filter(NES_diff == max(NES_diff)) %>% dplyr::select(pathway, NES_diff)
most_convergent <- variant_path_NES %>% filter(abs(NES_diff) == min(abs(NES_diff))) %>% dplyr::select(pathway, NES_diff)

most_divergent
```

    ## # A tibble: 1 × 2
    ##   pathway           NES_diff
    ##   <chr>                <dbl>
    ## 1 GOBP_NEURON_DEATH    0.878

``` r
most_convergent
```

    ## # A tibble: 1 × 2
    ##   pathway                         NES_diff
    ##   <chr>                              <dbl>
    ## 1 GOBP_NEUROTRANSMITTER_SECRETION  -0.0203

``` r
classify_paths_filt <- variant_path_NES %>% mutate(
  converge_diverge =
    case_when(
      (NES_diff > .33 | NES_diff < -0.09) ~ "divergent",
      (NES_diff <= .33 & NES_diff >= -0.09) ~ "convergent"
    )
)
classify_paths_filt
```

    ## # A tibble: 15 × 7
    ##    pathway                 rare_NES common_NES NES_diff MAD_filt median_NES_diff
    ##    <chr>                      <dbl>      <dbl>    <dbl>    <dbl>           <dbl>
    ##  1 GOBP_NEURON_DEATH           2.20       1.33   0.878     0.335           0.189
    ##  2 GOBP_POSITIVE_REGULATI…     1.95       1.52   0.429     0.335           0.189
    ##  3 GOBP_CENTRAL_NERVOUS_S…     1.78       1.49   0.297     0.335           0.189
    ##  4 GOBP_POSITIVE_REGULATI…     1.76       1.42   0.337     0.335           0.189
    ##  5 GOBP_CHROMATIN_REMODEL…     2.22       1.61   0.607     0.335           0.189
    ##  6 GOBP_SYNAPSE_ORGANIZAT…     2.12       1.39   0.723     0.335           0.189
    ##  7 GOBP_REGULATION_OF_NEU…     1.56       1.37   0.189     0.335           0.189
    ##  8 GOBP_VESICLE_MEDIATED_…     1.60       1.43   0.166     0.335           0.189
    ##  9 GOBP_NEURON_PROJECTION…     1.54       1.38   0.161     0.335           0.189
    ## 10 GOBP_CHROMATIN_ASSEMBL…     1.49       1.44   0.0466    0.335           0.189
    ## 11 GOBP_REGULATION_OF_NEU…     2.13       1.35   0.777     0.335           0.189
    ## 12 GOBP_NEURON_PROJECTION…     1.76       1.80  -0.0373    0.335           0.189
    ## 13 GOBP_NEUROTRANSMITTER_…     1.46       1.48  -0.0203    0.335           0.189
    ## 14 GOBP_NEGATIVE_REGULATI…     1.43       1.59  -0.156     0.335           0.189
    ## 15 KEGG_MTOR_SIGNALING_PA…     1.34       1.60  -0.259     0.335           0.189
    ## # ℹ 1 more variable: converge_diverge <chr>

``` r
divergent <- classify_paths_filt %>% filter(converge_diverge == "divergent")
divergent
```

    ## # A tibble: 8 × 7
    ##   pathway rare_NES common_NES NES_diff MAD_filt median_NES_diff converge_diverge
    ##   <chr>      <dbl>      <dbl>    <dbl>    <dbl>           <dbl> <chr>           
    ## 1 GOBP_N…     2.20       1.33    0.878    0.335           0.189 divergent       
    ## 2 GOBP_P…     1.95       1.52    0.429    0.335           0.189 divergent       
    ## 3 GOBP_P…     1.76       1.42    0.337    0.335           0.189 divergent       
    ## 4 GOBP_C…     2.22       1.61    0.607    0.335           0.189 divergent       
    ## 5 GOBP_S…     2.12       1.39    0.723    0.335           0.189 divergent       
    ## 6 GOBP_R…     2.13       1.35    0.777    0.335           0.189 divergent       
    ## 7 GOBP_N…     1.43       1.59   -0.156    0.335           0.189 divergent       
    ## 8 KEGG_M…     1.34       1.60   -0.259    0.335           0.189 divergent

``` r
nrow(divergent) # 8
```

    ## [1] 8

``` r
convergent <- classify_paths_filt %>% filter(converge_diverge == "convergent")
convergent
```

    ## # A tibble: 7 × 7
    ##   pathway rare_NES common_NES NES_diff MAD_filt median_NES_diff converge_diverge
    ##   <chr>      <dbl>      <dbl>    <dbl>    <dbl>           <dbl> <chr>           
    ## 1 GOBP_C…     1.78       1.49   0.297     0.335           0.189 convergent      
    ## 2 GOBP_R…     1.56       1.37   0.189     0.335           0.189 convergent      
    ## 3 GOBP_V…     1.60       1.43   0.166     0.335           0.189 convergent      
    ## 4 GOBP_N…     1.54       1.38   0.161     0.335           0.189 convergent      
    ## 5 GOBP_C…     1.49       1.44   0.0466    0.335           0.189 convergent      
    ## 6 GOBP_N…     1.76       1.80  -0.0373    0.335           0.189 convergent      
    ## 7 GOBP_N…     1.46       1.48  -0.0203    0.335           0.189 convergent

``` r
nrow(convergent) # 7
```

    ## [1] 7

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
    ## 16:  common
    ## 17:  common
    ## 18:  common
    ## 19:  common
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
    ##     variant
    ##                                                                     pathway
    ##                                                                      <char>
    ##  1:                                                       GOBP_NEURON_DEATH
    ##  2:                                GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH
    ##  3:                      GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION
    ##  4:                                GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS
    ##  5:                                               GOBP_CHROMATIN_REMODELING
    ##  6:                                               GOBP_SYNAPSE_ORGANIZATION
    ##  7:                              GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS
    ##  8:                              GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE
    ##  9:                                         GOBP_NEURON_PROJECTION_GUIDANCE
    ## 10:                                  GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY
    ## 11:                                         GOBP_REGULATION_OF_NEUROGENESIS
    ## 12: GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE
    ## 13:                                         GOBP_NEUROTRANSMITTER_SECRETION
    ## 14:    GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE
    ## 15:                                             KEGG_MTOR_SIGNALING_PATHWAY
    ## 16:                                               GOBP_CHROMATIN_REMODELING
    ## 17: GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE
    ## 18:                                               GOBP_SYNAPSE_ORGANIZATION
    ## 19:                              GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE
    ## 20:                      GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION
    ## 21:                                  GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY
    ## 22:                                GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS
    ## 23:                                         GOBP_NEUROTRANSMITTER_SECRETION
    ## 24:                                         GOBP_REGULATION_OF_NEUROGENESIS
    ## 25:                                GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH
    ## 26:                                         GOBP_NEURON_PROJECTION_GUIDANCE
    ## 27:                              GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS
    ## 28:                                                       GOBP_NEURON_DEATH
    ## 29:                                             KEGG_MTOR_SIGNALING_PATHWAY
    ## 30:    GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE
    ##                                                                     pathway
    ##             pval       padj    log2err        ES      NES  size  leadingEdge
    ##            <num>      <num>      <num>     <num>    <num> <int>       <list>
    ##  1: 2.312396e-03 0.03583278 0.43170770 0.9595688 2.203600   325 6326, 88....
    ##  2: 7.697156e-03 0.06800396 0.40701792 0.7955304 1.946403    80 2904, 14....
    ##  3: 3.996004e-02 0.09828010 0.22496609 0.7493933 1.781518   157 5728, 85....
    ##  4: 4.295704e-02 0.10023310 0.21654284 0.7457535 1.759636   199 23394, 8....
    ##  5: 4.695305e-02 0.10053476 0.20658792 0.9450217 2.215015   224 57680, 2....
    ##  6: 4.695305e-02 0.10053476 0.20658792 0.9249977 2.116172   370 8831, 23....
    ##  7: 4.995005e-02 0.10214505 0.19991523 0.6611369 1.559979   199 6529, 93....
    ##  8: 4.995005e-02 0.10214505 0.19991523 0.6767406 1.595984   183 5728, 14....
    ##  9: 5.494505e-02 0.10496720 0.19002331 0.6573949 1.542288   217 1826, 93....
    ## 10: 5.894106e-02 0.10909091 0.18302394 0.6311725 1.487743   196 11011, 1....
    ## 11: 6.393606e-02 0.11297440 0.17520405 0.9269535 2.128731   326 8831, 23....
    ## 12: 6.593407e-02 0.11384877 0.17232434 0.6973979 1.758611    33 1826, 65....
    ## 13: 8.691309e-02 0.13072878 0.14826150 0.6095070 1.460584   136 9378, 68....
    ## 14: 1.718282e-01 0.20988408 0.10027911 0.5653401 1.432189    24 6585, 64....
    ## 15: 1.918082e-01 0.23118603 0.09374654 0.5405839 1.344756    48 208, 253....
    ## 16: 9.687278e-05 0.01782459 0.53843410 0.3704143 1.607976   227 9759, 21....
    ## 17: 1.305539e-03 0.12010963 0.45505987 0.5462061 1.795908    36 7473, 91....
    ## 18: 2.040575e-03 0.12515526 0.43170770 0.3106873 1.392875   384 4137, 11....
    ## 19: 4.906584e-03 0.17026099 0.40701792 0.3350995 1.430100   186 783, 990....
    ## 20: 6.777987e-03 0.17026099 0.40701792 0.3530124 1.485000   160 4137, 91....
    ## 21: 8.089132e-03 0.17026099 0.38073040 0.3368143 1.441185   191 3010, 83....
    ## 22: 8.814512e-03 0.17026099 0.38073040 0.3294199 1.422300   208 4137, 74....
    ## 23: 9.069706e-03 0.17026099 0.38073040 0.3592682 1.480884   137 783, 990....
    ## 24: 9.253315e-03 0.17026099 0.38073040 0.3032521 1.351234   338 4137, 74....
    ## 25: 1.092743e-02 0.17372614 0.38073040 0.3951873 1.517007    81 4137, 67....
    ## 26: 1.280778e-02 0.17372614 0.38073040 0.3180250 1.381370   221 7473, 91....
    ## 27: 1.381113e-02 0.17372614 0.38073040 0.3190564 1.371380   201 783, 990....
    ## 28: 1.408892e-02 0.17372614 0.38073040 0.2979639 1.325979   327 4137, 67....
    ## 29: 1.507598e-02 0.17372614 0.38073040 0.4540769 1.604229    49 2475, 60....
    ## 30: 1.510662e-02 0.17372614 0.38073040 0.5234268 1.588080    26 7473, 10....
    ##             pval       padj    log2err        ES      NES  size  leadingEdge
    ##     rare_NES common_NES    NES_diff  MAD_filt median_NES_diff converge_diverge
    ##        <num>      <num>       <num>     <num>           <num>           <char>
    ##  1: 2.203600   1.325979  0.87762091 0.3349146       0.1885995        divergent
    ##  2: 1.946403   1.517007  0.42939546 0.3349146       0.1885995        divergent
    ##  3: 1.781518   1.485000  0.29651800 0.3349146       0.1885995       convergent
    ##  4: 1.759636   1.422300  0.33733565 0.3349146       0.1885995        divergent
    ##  5: 2.215015   1.607976  0.60703960 0.3349146       0.1885995        divergent
    ##  6: 2.116172   1.392875  0.72329751 0.3349146       0.1885995        divergent
    ##  7: 1.559979   1.371380  0.18859952 0.3349146       0.1885995       convergent
    ##  8: 1.595984   1.430100  0.16588390 0.3349146       0.1885995       convergent
    ##  9: 1.542288   1.381370  0.16091758 0.3349146       0.1885995       convergent
    ## 10: 1.487743   1.441185  0.04655805 0.3349146       0.1885995       convergent
    ## 11: 2.128731   1.351234  0.77749650 0.3349146       0.1885995        divergent
    ## 12: 1.758611   1.795908 -0.03729727 0.3349146       0.1885995       convergent
    ## 13: 1.460584   1.480884 -0.02030026 0.3349146       0.1885995       convergent
    ## 14: 1.432189   1.588080 -0.15589156 0.3349146       0.1885995        divergent
    ## 15: 1.344756   1.604229 -0.25947279 0.3349146       0.1885995        divergent
    ## 16: 2.215015   1.607976  0.60703960 0.3349146       0.1885995        divergent
    ## 17: 1.758611   1.795908 -0.03729727 0.3349146       0.1885995       convergent
    ## 18: 2.116172   1.392875  0.72329751 0.3349146       0.1885995        divergent
    ## 19: 1.595984   1.430100  0.16588390 0.3349146       0.1885995       convergent
    ## 20: 1.781518   1.485000  0.29651800 0.3349146       0.1885995       convergent
    ## 21: 1.487743   1.441185  0.04655805 0.3349146       0.1885995       convergent
    ## 22: 1.759636   1.422300  0.33733565 0.3349146       0.1885995        divergent
    ## 23: 1.460584   1.480884 -0.02030026 0.3349146       0.1885995       convergent
    ## 24: 2.128731   1.351234  0.77749650 0.3349146       0.1885995        divergent
    ## 25: 1.946403   1.517007  0.42939546 0.3349146       0.1885995        divergent
    ## 26: 1.542288   1.381370  0.16091758 0.3349146       0.1885995       convergent
    ## 27: 1.559979   1.371380  0.18859952 0.3349146       0.1885995       convergent
    ## 28: 2.203600   1.325979  0.87762091 0.3349146       0.1885995        divergent
    ## 29: 1.344756   1.604229 -0.25947279 0.3349146       0.1885995        divergent
    ## 30: 1.432189   1.588080 -0.15589156 0.3349146       0.1885995        divergent
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
write.csv(path_results, "../results/path_results.csv")
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
