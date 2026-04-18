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
    ##  [2] "GOBP_INHIBITORY_SYNAPSE_ASSEMBLY"                         
    ##  [3] "GOBP_LONG_TERM_SYNAPTIC_POTENTIATION"                     
    ##  [4] "GOBP_NEGATIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION"        
    ##  [5] "GOBP_NEUROMUSCULAR_PROCESS"                               
    ##  [6] "GOBP_NEURONAL_ACTION_POTENTIAL"                           
    ##  [7] "GOBP_NEURON_DEATH"                                        
    ##  [8] "GOBP_NEURON_FATE_SPECIFICATION"                           
    ##  [9] "GOBP_NEURON_PROJECTION_ORGANIZATION"                      
    ## [10] "GOBP_NEUROTRANSMITTER_UPTAKE"                             
    ## [11] "GOBP_POSTSYNAPTIC_SPECIALIZATION_ORGANIZATION"            
    ## [12] "GOBP_REGULATION_OF_AXONOGENESIS"                          
    ## [13] "GOBP_REGULATION_OF_CHROMATIN_ORGANIZATION"                
    ## [14] "GOBP_REGULATION_OF_LONG_TERM_NEURONAL_SYNAPTIC_PLASTICITY"
    ## [15] "GOBP_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_POTENTIAL"       
    ## [16] "GOBP_REGULATION_OF_SYNAPSE_ASSEMBLY"                      
    ## [17] "GOBP_SYNAPTIC_VESICLE_LOCALIZATION"

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

    ##    variant                                   pathway         pval       padj
    ##     <char>                                    <char>        <num>      <num>
    ## 1:    rare             GOBP_NEURON_APOPTOTIC_PROCESS 0.0002210583 0.04023261
    ## 2:    rare                         GOBP_NEURON_DEATH 0.0010584774 0.08975834
    ## 3:    rare GOBP_POSTSYNAPTIC_SPECIALIZATION_ASSEMBLY 0.0024603865 0.08975834
    ## 4:    rare  GOBP_NEGATIVE_REGULATION_OF_AXONOGENESIS 0.0035333199 0.08975834
    ## 5:    rare     GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE 0.0049450725 0.08975834
    ## 6:    rare          GOBP_EXCITATORY_SYNAPSE_ASSEMBLY 0.0064322875 0.08975834
    ##      log2err        ES      NES  size  leadingEdge
    ##        <num>     <num>    <num> <int>       <list>
    ## 1: 0.5188481 0.9695995 2.319869   218 6326, 88....
    ## 2: 0.4550599 0.9595688 2.242748   325 6326, 88....
    ## 3: 0.4317077 0.9353871 2.427934    20 5728, 85....
    ## 4: 0.4317077 0.9848314 2.521540    56   8831, 5728
    ## 5: 0.4070179 0.9924235 2.577443    19  8831, 22941
    ## 6: 0.4070179 0.9278930 2.407059    24 5728, 85....

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
    ## 1:                                      GOBP_POSTSYNAPTIC_SPECIALIZATION_ASSEMBLY
    ## 2: GOBP_ADENYLATE_CYCLASE_INHIBITING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY
    ## 3:                     GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_REGENERATION
    ## 4:                                               GOBP_EXCITATORY_SYNAPSE_ASSEMBLY
    ## 5:                                  GOBP_REGULATION_OF_SYNAPTIC_VESICLE_RECYCLING
    ## 6:                                 GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_BALANCE
    ##         pval  padj    log2err         ES        NES  size  leadingEdge
    ##        <num> <num>      <num>      <num>      <num> <int>       <list>
    ## 1: 0.9852507     1 0.03221532  0.1815711  0.5370878    20 78999, 2....
    ## 2: 0.9900990     1 0.02286940  0.1533360  0.5919384    75 2550, 29....
    ## 3: 0.9939302     1 0.03316264  0.1884915  0.5118334    15 5728, 71....
    ## 4: 0.9942693     1 0.03034673  0.1747264  0.5406954    25 64101, 7....
    ## 5: 0.9970238     1 0.03207042  0.1610146  0.4701332    19 6622, 60....
    ## 6: 1.0000000     1 0.07977059 -0.1427050 -0.6197831    49 7401, 64....

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
    ##   [4] "KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM"                               
    ##   [5] "GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY"                                        
    ##   [6] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION"                            
    ##   [7] "GOBP_REGULATION_OF_NEUROGENESIS"                                               
    ##   [8] "KEGG_MTOR_SIGNALING_PATHWAY"                                                   
    ##   [9] "GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS"                                      
    ##  [10] "GOBP_NEUROTRANSMITTER_SECRETION"                                               
    ##  [11] "GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH"                                      
    ##  [12] "GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE"          
    ##  [13] "GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE"                                    
    ##  [14] "GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS"                                    
    ##  [15] "GOBP_NEURON_DEATH"                                                             
    ##  [16] "GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS"                                      
    ##  [17] "GOBP_NEURON_PROJECTION_GUIDANCE"                                               
    ##  [18] "KEGG_AXON_GUIDANCE"                                                            
    ##  [19] "GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS"                    
    ##  [20] "GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                     
    ##  [21] "GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS"                                              
    ##  [22] "GOBP_FOREBRAIN_GENERATION_OF_NEURONS"                                          
    ##  [23] "GOBP_NEUROTRANSMITTER_TRANSPORT"                                               
    ##  [24] "GOBP_ANTEROGRADE_AXONAL_TRANSPORT"                                             
    ##  [25] "GOBP_SYNAPTIC_VESICLE_PRIMING"                                                 
    ##  [26] "GOBP_REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY"                              
    ##  [27] "GOBP_LONG_TERM_SYNAPTIC_DEPRESSION"                                            
    ##  [28] "GOBP_AXONAL_TRANSPORT"                                                         
    ##  [29] "GOBP_NEGATIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION"                             
    ##  [30] "GOBP_REGULATION_OF_PROTEIN_LOCALIZATION_TO_SYNAPSE"                            
    ##  [31] "GOBP_REGULATION_OF_AXONOGENESIS"                                               
    ##  [32] "GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE"                                         
    ##  [33] "GOBP_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT"                                 
    ##  [34] "GOBP_RETINAL_GANGLION_CELL_AXON_GUIDANCE"                                      
    ##  [35] "GOBP_NEURON_FATE_COMMITMENT"                                                   
    ##  [36] "GOBP_NEGATIVE_REGULATION_OF_NEURON_DIFFERENTIATION"                            
    ##  [37] "GOBP_NEURON_PROJECTION_ORGANIZATION"                                           
    ##  [38] "GOBP_GLUTAMATE_SECRETION"                                                      
    ##  [39] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_AXONOGENESIS"                               
    ##  [40] "GOBP_SYNAPSE_MATURATION"                                                       
    ##  [41] "GOBP_NEURON_PROJECTION_EXTENSION"                                              
    ##  [42] "GOBP_RETROGRADE_AXONAL_TRANSPORT"                                              
    ##  [43] "GOBP_PRESYNAPTIC_ENDOCYTOSIS"                                                  
    ##  [44] "GOBP_NEURON_DEATH_IN_RESPONSE_TO_OXIDATIVE_STRESS"                             
    ##  [45] "GOBP_NEUROBLAST_PROLIFERATION"                                                 
    ##  [46] "GOBP_REGULATION_OF_GLUTAMATE_SECRETION"                                        
    ##  [47] "GOBP_POSITIVE_REGULATION_OF_NEURON_APOPTOTIC_PROCESS"                          
    ##  [48] "GOBP_DNA_REPLICATION_DEPENDENT_CHROMATIN_ORGANIZATION"                         
    ##  [49] "GOBP_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                              
    ##  [50] "GOBP_POSTSYNAPSE_ORGANIZATION"                                                 
    ##  [51] "GOBP_NEGATIVE_REGULATION_OF_AXONOGENESIS"                                      
    ##  [52] "GOBP_FOREBRAIN_NEURON_DIFFERENTIATION"                                         
    ##  [53] "GOBP_AXON_DEVELOPMENT"                                                         
    ##  [54] "GOBP_NEURON_APOPTOTIC_PROCESS"                                                 
    ##  [55] "GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION"                                    
    ##  [56] "GOBP_AXONAL_TRANSPORT_OF_MITOCHONDRION"                                        
    ##  [57] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DEVELOPMENT"                                
    ##  [58] "GOBP_MOTOR_NEURON_AXON_GUIDANCE"                                               
    ##  [59] "GOBP_PROTEIN_LOCALIZATION_TO_SYNAPSE"                                          
    ##  [60] "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION"                                  
    ##  [61] "GOBP_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_POTENTIAL"                            
    ##  [62] "GOBP_AXON_EXTENSION"                                                           
    ##  [63] "GOBP_SYNAPTIC_VESICLE_LOCALIZATION"                                            
    ##  [64] "GOBP_SYNAPSE_ASSEMBLY"                                                         
    ##  [65] "GOBP_REGULATION_OF_NEURONAL_SYNAPTIC_PLASTICITY"                               
    ##  [66] "GOBP_POSITIVE_REGULATION_OF_NEURON_DIFFERENTIATION"                            
    ##  [67] "GOBP_REGULATION_OF_NEURON_DIFFERENTIATION"                                     
    ##  [68] "GOBP_POSITIVE_REGULATION_OF_NEUROBLAST_PROLIFERATION"                          
    ##  [69] "GOBP_POSITIVE_REGULATION_OF_CHROMATIN_ORGANIZATION"                            
    ##  [70] "GOBP_SYNAPTIC_MEMBRANE_ADHESION"                                               
    ##  [71] "GOBP_NEUROEPITHELIAL_CELL_DIFFERENTIATION"                                     
    ##  [72] "GOBP_SYNAPTIC_VESICLE_TRANSPORT"                                               
    ##  [73] "GOBP_POSITIVE_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT"                        
    ##  [74] "GOBP_POSITIVE_REGULATION_OF_NEUROTRANSMITTER_SECRETION"                        
    ##  [75] "GOBP_POSTSYNAPTIC_MEMBRANE_ORGANIZATION"                                       
    ##  [76] "GOBP_GLUTAMATE_METABOLIC_PROCESS"                                              
    ##  [77] "GOBP_REGULATION_OF_NEUROBLAST_PROLIFERATION"                                   
    ##  [78] "GOBP_LONG_TERM_SYNAPTIC_POTENTIATION"                                          
    ##  [79] "GOBP_REGULATION_OF_SYNAPTIC_PLASTICITY"                                        
    ##  [80] "GOBP_SYNAPTIC_VESICLE_MEMBRANE_ORGANIZATION"                                   
    ##  [81] "GOBP_POSITIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION"                             
    ##  [82] "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY"                                           
    ##  [83] "GOBP_NEUROPEPTIDE_SIGNALING_PATHWAY"                                           
    ##  [84] "GOBP_REGULATION_OF_LONG_TERM_NEURONAL_SYNAPTIC_PLASTICITY"                     
    ##  [85] "GOBP_AXONEMAL_DYNEIN_COMPLEX_ASSEMBLY"                                         
    ##  [86] "GOBP_NEURON_RECOGNITION"                                                       
    ##  [87] "GOBP_SPINAL_CORD_MOTOR_NEURON_DIFFERENTIATION"                                 
    ##  [88] "GOBP_NEURONAL_STEM_CELL_POPULATION_MAINTENANCE"                                
    ##  [89] "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING"                                   
    ##  [90] "GOBP_SYNAPTIC_VESICLE_RECYCLING"                                               
    ##  [91] "GOBP_POSITIVE_REGULATION_OF_NEUROINFLAMMATORY_RESPONSE"                        
    ##  [92] "GOBP_NEGATIVE_REGULATION_OF_NEURON_DEATH"                                      
    ##  [93] "GOBP_ENSHEATHMENT_OF_NEURONS"                                                  
    ##  [94] "GOBP_NEURON_MIGRATION"                                                         
    ##  [95] "GOBP_MODULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL"                          
    ##  [96] "GOBP_NEUROTROPHIN_TRK_RECEPTOR_SIGNALING_PATHWAY"                              
    ##  [97] "GOBP_NEUROMUSCULAR_JUNCTION_DEVELOPMENT"                                       
    ##  [98] "GOBP_HETEROCHROMATIN_ORGANIZATION"                                             
    ##  [99] "GOBP_RESPONSE_TO_LEUKEMIA_INHIBITORY_FACTOR"                                   
    ## [100] "GOBP_NEURON_CELLULAR_HOMEOSTASIS"                                              
    ## [101] "GOBP_REGULATION_OF_LONG_TERM_SYNAPTIC_POTENTIATION"                            
    ## [102] "GOBP_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY"                                     
    ## [103] "GOBP_REGULATION_OF_SYNAPTIC_VESICLE_EXOCYTOSIS"                                
    ## [104] "GOBP_PROTEIN_LOCALIZATION_TO_POSTSYNAPSE"                                      
    ## [105] "GOBP_POSITIVE_REGULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL"                 
    ## [106] "GOBP_POSTSYNAPTIC_SPECIALIZATION_ORGANIZATION"                                 
    ## [107] "GOBP_REGULATION_OF_NEUROTRANSMITTER_UPTAKE"                                    
    ## [108] "GOBP_REGULATION_OF_POSTSYNAPSE_ORGANIZATION"                                   
    ## [109] "GOBP_AXON_ENSHEATHMENT_IN_CENTRAL_NERVOUS_SYSTEM"                              
    ## [110] "GOBP_CEREBRAL_CORTEX_NEURON_DIFFERENTIATION"                                   
    ## [111] "GOBP_SYNAPTIC_TRANSMISSION_CHOLINERGIC"                                        
    ## [112] "GOBP_REGULATION_OF_RECEPTOR_LOCALIZATION_TO_SYNAPSE"                           
    ## [113] "GOBP_NEURON_FATE_SPECIFICATION"                                                
    ## [114] "GOBP_G_PROTEIN_COUPLED_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY"                   
    ## [115] "GOBP_MODIFICATION_OF_SYNAPTIC_STRUCTURE"                                       
    ## [116] "GOBP_SYNAPTIC_VESICLE_CYTOSKELETAL_TRANSPORT"                                  
    ## [117] "GOBP_REGULATION_OF_SYNAPSE_ASSEMBLY"                                           
    ## [118] "GOBP_DNA_REPLICATION_INDEPENDENT_CHROMATIN_ORGANIZATION"                       
    ## [119] "GOBP_NEUROTRANSMITTER_REUPTAKE"                                                
    ## [120] "GOBP_DOPAMINERGIC_NEURON_DIFFERENTIATION"                                      
    ## [121] "GOBP_PROTEIN_LOCALIZATION_TO_CHROMATIN"                                        
    ## [122] "GOBP_NEGATIVE_REGULATION_OF_NEURON_APOPTOTIC_PROCESS"                          
    ## [123] "GOBP_NEUROMUSCULAR_SYNAPTIC_TRANSMISSION"                                      
    ## [124] "GOBP_NEURON_PROJECTION_ARBORIZATION"                                           
    ## [125] "GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_POSTURE"                                
    ## [126] "GOBP_INHIBITORY_SYNAPSE_ASSEMBLY"                                              
    ## [127] "GOBP_AXONAL_FASCICULATION"                                                     
    ## [128] "GOBP_NEUROINFLAMMATORY_RESPONSE"                                               
    ## [129] "GOBP_PRESYNAPSE_ORGANIZATION"                                                  
    ## [130] "GOBP_NEURON_MATURATION"                                                        
    ## [131] "GOBP_NEUROMUSCULAR_PROCESS"                                                    
    ## [132] "GOBP_NEUROTROPHIN_SIGNALING_PATHWAY"                                           
    ## [133] "GOBP_REGULATION_OF_PRESYNAPSE_ORGANIZATION"                                    
    ## [134] "GOBP_NEURONAL_ACTION_POTENTIAL"                                                
    ## [135] "GOBP_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_NEUROTRANSMITTER_RECEPTOR_LEVELS"     
    ## [136] "GOBP_POSITIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"               
    ## [137] "GOBP_FOREBRAIN_NEURON_DEVELOPMENT"                                             
    ## [138] "GOBP_FACULTATIVE_HETEROCHROMATIN_ASSEMBLY"                                     
    ## [139] "GOBP_CHROMATIN_DISASSEMBLY"                                                    
    ## [140] "GOBP_REGULATION_OF_NEUROTRANSMITTER_RECEPTOR_ACTIVITY"                         
    ## [141] "GOBP_POSTSYNAPTIC_SIGNAL_TRANSDUCTION"                                         
    ## [142] "GOBP_REGULATION_OF_DNA_METHYLATION_DEPENDENT_HETEROCHROMATIN_ASSEMBLY"         
    ## [143] "GOBP_CHEMICAL_SYNAPTIC_TRANSMISSION_POSTSYNAPTIC"                              
    ## [144] "GOBP_DNA_METHYLATION_DEPENDENT_HETEROCHROMATIN_ASSEMBLY"                       
    ## [145] "GOBP_AXONEME_ASSEMBLY"                                                         
    ## [146] "GOBP_POSITIVE_REGULATION_OF_LONG_TERM_SYNAPTIC_POTENTIATION"                   
    ## [147] "GOBP_NEGATIVE_REGULATION_OF_OXIDATIVE_STRESS_INDUCED_NEURON_DEATH"             
    ## [148] "GOBP_NEURON_PROJECTION_REGENERATION"                                           
    ## [149] "GOBP_REGULATION_OF_CHROMATIN_ASSEMBLY"                                         
    ## [150] "GOBP_REGULATION_OF_CHROMATIN_BINDING"                                          
    ## [151] "GOBP_REGULATION_OF_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY"                          
    ## [152] "GOBP_REGULATION_OF_NEURON_PROJECTION_REGENERATION"                             
    ## [153] "GOBP_MIDBRAIN_DOPAMINERGIC_NEURON_DIFFERENTIATION"                             
    ## [154] "GOBP_NEUROTRANSMITTER_METABOLIC_PROCESS"                                       
    ## [155] "GOBP_SYNAPTONEMAL_COMPLEX_ORGANIZATION"                                        
    ## [156] "GOBP_RESPONSE_TO_AXON_INJURY"                                                  
    ## [157] "GOBP_RECEPTOR_LOCALIZATION_TO_SYNAPSE"                                         
    ## [158] "GOBP_SYNAPTIC_TRANSMISSION_GABAERGIC"                                          
    ## [159] "GOBP_POSITIVE_REGULATION_OF_NEURON_MIGRATION"                                  
    ## [160] "GOBP_REGULATION_OF_CHROMATIN_ORGANIZATION"                                     
    ## [161] "GOBP_SYNAPTIC_TRANSMISSION_DOPAMINERGIC"                                       
    ## [162] "GOBP_POSITIVE_REGULATION_OF_AXON_EXTENSION"                                    
    ## [163] "GOBP_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"                                      
    ## [164] "GOBP_NEUROTRANSMITTER_UPTAKE"                                                  
    ## [165] "GOBP_POSTSYNAPSE_ASSEMBLY"                                                     
    ## [166] "GOBP_RNA_MEDIATED_GENE_SILENCING_BY_INHIBITION_OF_TRANSLATION"                 
    ## [167] "GOBP_POSITIVE_REGULATION_OF_SYNAPSE_ASSEMBLY"                                  
    ## [168] "GOBP_CALCIUM_ION_REGULATED_EXOCYTOSIS_OF_NEUROTRANSMITTER"                     
    ## [169] "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"                        
    ## [170] "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GABAERGIC"                            
    ## [171] "GOBP_REGULATION_OF_NEURON_MIGRATION"                                           
    ## [172] "GOBP_NEUROTRANSMITTER_RECEPTOR_TRANSPORT"                                      
    ## [173] "GOBP_L_GLUTAMATE_TRANSMEMBRANE_TRANSPORT"                                      
    ## [174] "GOBP_L_GLUTAMATE_IMPORT_ACROSS_PLASMA_MEMBRANE"                                
    ## [175] "GOBP_POSITIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                     
    ## [176] "GOBP_NEUROTRANSMITTER_RECEPTOR_INTERNALIZATION"                                
    ## [177] "GOBP_POSTSYNAPTIC_SPECIALIZATION_ASSEMBLY"                                     
    ## [178] "GOBP_ADENYLATE_CYCLASE_INHIBITING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY"
    ## [179] "GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_REGENERATION"                    
    ## [180] "GOBP_EXCITATORY_SYNAPSE_ASSEMBLY"                                              
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
    ##   1:                                                  GOBP_NEURON_APOPTOTIC_PROCESS
    ##   2:                                                              GOBP_NEURON_DEATH
    ##   3:                                      GOBP_POSTSYNAPTIC_SPECIALIZATION_ASSEMBLY
    ##   4:                                       GOBP_NEGATIVE_REGULATION_OF_AXONOGENESIS
    ##   5:                                          GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE
    ##  ---                                                                               
    ## 360: GOBP_ADENYLATE_CYCLASE_INHIBITING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY
    ## 361:                     GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_REGENERATION
    ## 362:                                               GOBP_EXCITATORY_SYNAPSE_ASSEMBLY
    ## 363:                                  GOBP_REGULATION_OF_SYNAPTIC_VESICLE_RECYCLING
    ## 364:                                 GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_BALANCE
    ##              pval       padj    log2err         ES        NES  size
    ##             <num>      <num>      <num>      <num>      <num> <int>
    ##   1: 0.0002210583 0.04023261 0.51884808  0.9695995  2.3198685   218
    ##   2: 0.0010584774 0.08975834 0.45505987  0.9595688  2.2427484   325
    ##   3: 0.0024603865 0.08975834 0.43170770  0.9353871  2.4279337    20
    ##   4: 0.0035333199 0.08975834 0.43170770  0.9848314  2.5215405    56
    ##   5: 0.0049450725 0.08975834 0.40701792  0.9924235  2.5774428    19
    ##  ---                                                               
    ## 360: 0.9900990099 1.00000000 0.02286940  0.1533360  0.5919384    75
    ## 361: 0.9939301973 1.00000000 0.03316264  0.1884915  0.5118334    15
    ## 362: 0.9942693410 1.00000000 0.03034673  0.1747264  0.5406954    25
    ## 363: 0.9970238095 1.00000000 0.03207042  0.1610146  0.4701332    19
    ## 364: 1.0000000000 1.00000000 0.07977059 -0.1427050 -0.6197831    49
    ##       leadingEdge
    ##            <list>
    ##   1: 6326, 88....
    ##   2: 6326, 88....
    ##   3: 5728, 85....
    ##   4:   8831, 5728
    ##   5:  8831, 22941
    ##  ---             
    ## 360: 2550, 29....
    ## 361: 5728, 71....
    ## 362: 64101, 7....
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
    ##  1 GOBP_NEURON_APOPTOTIC_PROCESS           2.32      1.20     0.0402       0.381
    ##  2 GOBP_NEURON_DEATH                       2.24      1.34     0.0898       0.150
    ##  3 GOBP_POSTSYNAPTIC_SPECIALIZATION_A…     2.43      0.537    0.0898       1    
    ##  4 GOBP_NEGATIVE_REGULATION_OF_AXONOG…     2.52      1.30     0.0898       0.370
    ##  5 GOBP_MAINTENANCE_OF_SYNAPSE_STRUCT…     2.58      1.53     0.0898       0.279
    ##  6 GOBP_EXCITATORY_SYNAPSE_ASSEMBLY        2.41      0.541    0.0898       1    
    ##  7 GOBP_POSITIVE_REGULATION_OF_AXON_E…     2.34      0.744    0.0898       0.942
    ##  8 GOBP_POSTSYNAPTIC_SPECIALIZATION_O…     2.55      1.07     0.0898       0.667
    ##  9 GOBP_POSITIVE_REGULATION_OF_NEURON…     1.86      0.671    0.0898       1    
    ## 10 GOBP_POSTSYNAPSE_ASSEMBLY               2.40      0.721    0.0898       0.942
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
    ##  1 GOBP_NEURON_AP…     2.32      1.20     0.0402       0.381    1.12       0.586
    ##  2 GOBP_NEURON_DE…     2.24      1.34     0.0898       0.150    0.901      0.586
    ##  3 GOBP_POSTSYNAP…     2.43      0.537    0.0898       1        1.89       0.586
    ##  4 GOBP_NEGATIVE_…     2.52      1.30     0.0898       0.370    1.22       0.586
    ##  5 GOBP_MAINTENAN…     2.58      1.53     0.0898       0.279    1.05       0.586
    ##  6 GOBP_EXCITATOR…     2.41      0.541    0.0898       1        1.87       0.586
    ##  7 GOBP_POSITIVE_…     2.34      0.744    0.0898       0.942    1.60       0.586
    ##  8 GOBP_POSTSYNAP…     2.55      1.07     0.0898       0.667    1.49       0.586
    ##  9 GOBP_POSITIVE_…     1.86      0.671    0.0898       1        1.19       0.586
    ## 10 GOBP_POSTSYNAP…     2.40      0.721    0.0898       0.942    1.68       0.586
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
    ## 1 GOBP_POSITIVE_REGULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL     3.41

``` r
most_convergent
```

    ## # A tibble: 1 × 2
    ##   pathway                         NES_diff
    ##   <chr>                              <dbl>
    ## 1 GOBP_NEUROTRANSMITTER_SECRETION  0.00154

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
    ##  1 GOBP_NEURON_AP…     2.32      1.20     0.0402       0.381    1.12       0.586
    ##  2 GOBP_NEURON_DE…     2.24      1.34     0.0898       0.150    0.901      0.586
    ##  3 GOBP_POSTSYNAP…     2.43      0.537    0.0898       1        1.89       0.586
    ##  4 GOBP_NEGATIVE_…     2.52      1.30     0.0898       0.370    1.22       0.586
    ##  5 GOBP_MAINTENAN…     2.58      1.53     0.0898       0.279    1.05       0.586
    ##  6 GOBP_EXCITATOR…     2.41      0.541    0.0898       1        1.87       0.586
    ##  7 GOBP_POSITIVE_…     2.34      0.744    0.0898       0.942    1.60       0.586
    ##  8 GOBP_POSTSYNAP…     2.55      1.07     0.0898       0.667    1.49       0.586
    ##  9 GOBP_POSITIVE_…     1.86      0.671    0.0898       1        1.19       0.586
    ## 10 GOBP_POSTSYNAP…     2.40      0.721    0.0898       0.942    1.68       0.586
    ## # ℹ 172 more rows
    ## # ℹ 2 more variables: median_NES_diff <dbl>, converge_diverge <chr>

``` r
divergent_unfilt <- classify_paths %>% filter(converge_diverge == "divergent")
divergent_unfilt
```

    ## # A tibble: 60 × 9
    ##    pathway         rare_NES common_NES rare_padj common_padj NES_diff MAD_unfilt
    ##    <chr>              <dbl>      <dbl>     <dbl>       <dbl>    <dbl>      <dbl>
    ##  1 GOBP_POSTSYNAP…     2.43      0.537    0.0898       1         1.89      0.586
    ##  2 GOBP_EXCITATOR…     2.41      0.541    0.0898       1         1.87      0.586
    ##  3 GOBP_POSITIVE_…     2.34      0.744    0.0898       0.942     1.60      0.586
    ##  4 GOBP_POSTSYNAP…     2.55      1.07     0.0898       0.667     1.49      0.586
    ##  5 GOBP_POSTSYNAP…     2.40      0.721    0.0898       0.942     1.68      0.586
    ##  6 GOBP_POSITIVE_…     2.22      0.721    0.0898       0.969     1.50      0.586
    ##  7 GOBP_REGULATIO…     2.57      1.18     0.0898       0.509     1.38      0.586
    ##  8 GOBP_NEURONAL_…     2.55      0.923    0.0898       0.799     1.62      0.586
    ##  9 GOBP_SYNAPTIC_…     2.20      0.793    0.0898       0.910     1.41      0.586
    ## 10 GOBP_REGULATIO…     2.32      0.653    0.0898       1         1.67      0.586
    ## # ℹ 50 more rows
    ## # ℹ 2 more variables: median_NES_diff <dbl>, converge_diverge <chr>

``` r
nrow(divergent_unfilt)
```

    ## [1] 60

``` r
convergent_unfilt <- classify_paths %>% filter(converge_diverge == "convergent")
convergent_unfilt
```

    ## # A tibble: 122 × 9
    ##    pathway         rare_NES common_NES rare_padj common_padj NES_diff MAD_unfilt
    ##    <chr>              <dbl>      <dbl>     <dbl>       <dbl>    <dbl>      <dbl>
    ##  1 GOBP_NEURON_AP…     2.32      1.20     0.0402       0.381    1.12       0.586
    ##  2 GOBP_NEURON_DE…     2.24      1.34     0.0898       0.150    0.901      0.586
    ##  3 GOBP_NEGATIVE_…     2.52      1.30     0.0898       0.370    1.22       0.586
    ##  4 GOBP_MAINTENAN…     2.58      1.53     0.0898       0.279    1.05       0.586
    ##  5 GOBP_POSITIVE_…     1.86      0.671    0.0898       1        1.19       0.586
    ##  6 GOBP_POSITIVE_…     2.07      1.47     0.0898       0.150    0.606      0.586
    ##  7 GOBP_REGULATIO…     2.54      1.24     0.0898       0.447    1.30       0.586
    ##  8 GOBP_CHEMICAL_…     2.08      0.891    0.0898       0.862    1.19       0.586
    ##  9 GOBP_NEGATIVE_…     2.33      1.38     0.0898       0.275    0.948      0.586
    ## 10 GOBP_LONG_TERM…     2.18      1.15     0.0898       0.497    1.03       0.586
    ## # ℹ 112 more rows
    ## # ℹ 2 more variables: median_NES_diff <dbl>, converge_diverge <chr>

``` r
nrow(convergent_unfilt)
```

    ## [1] 122

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

    ##   [1] "GOBP_NEURON_APOPTOTIC_PROCESS"                                                 
    ##   [2] "GOBP_POSTSYNAPTIC_SPECIALIZATION_ASSEMBLY"                                     
    ##   [3] "GOBP_NEGATIVE_REGULATION_OF_AXONOGENESIS"                                      
    ##   [4] "GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE"                                         
    ##   [5] "GOBP_EXCITATORY_SYNAPSE_ASSEMBLY"                                              
    ##   [6] "GOBP_POSITIVE_REGULATION_OF_AXON_EXTENSION"                                    
    ##   [7] "GOBP_POSTSYNAPTIC_SPECIALIZATION_ORGANIZATION"                                 
    ##   [8] "GOBP_POSITIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                     
    ##   [9] "GOBP_POSTSYNAPSE_ASSEMBLY"                                                     
    ##  [10] "GOBP_POSITIVE_REGULATION_OF_SYNAPSE_ASSEMBLY"                                  
    ##  [11] "GOBP_REGULATION_OF_LONG_TERM_NEURONAL_SYNAPTIC_PLASTICITY"                     
    ##  [12] "GOBP_NEURONAL_ACTION_POTENTIAL"                                                
    ##  [13] "GOBP_SYNAPTIC_TRANSMISSION_GABAERGIC"                                          
    ##  [14] "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GABAERGIC"                            
    ##  [15] "GOBP_MODULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL"                          
    ##  [16] "GOBP_REGULATION_OF_NEURONAL_SYNAPTIC_PLASTICITY"                               
    ##  [17] "GOBP_POSITIVE_REGULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL"                 
    ##  [18] "GOBP_CHEMICAL_SYNAPTIC_TRANSMISSION_POSTSYNAPTIC"                              
    ##  [19] "GOBP_NEGATIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION"                             
    ##  [20] "GOBP_LONG_TERM_SYNAPTIC_POTENTIATION"                                          
    ##  [21] "GOBP_NEURON_PROJECTION_ORGANIZATION"                                           
    ##  [22] "GOBP_POSITIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"               
    ##  [23] "GOBP_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY"                                     
    ##  [24] "GOBP_REGULATION_OF_POSTSYNAPSE_ORGANIZATION"                                   
    ##  [25] "GOBP_REGULATION_OF_SYNAPSE_ASSEMBLY"                                           
    ##  [26] "GOBP_ENSHEATHMENT_OF_NEURONS"                                                  
    ##  [27] "GOBP_POSITIVE_REGULATION_OF_LONG_TERM_SYNAPTIC_POTENTIATION"                   
    ##  [28] "GOBP_LONG_TERM_SYNAPTIC_DEPRESSION"                                            
    ##  [29] "GOBP_NEUROMUSCULAR_PROCESS"                                                    
    ##  [30] "GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_REGENERATION"                    
    ##  [31] "GOBP_REGULATION_OF_CHROMATIN_ORGANIZATION"                                     
    ##  [32] "GOBP_SYNAPSE_MATURATION"                                                       
    ##  [33] "GOBP_RETROGRADE_AXONAL_TRANSPORT"                                              
    ##  [34] "GOBP_POSTSYNAPSE_ORGANIZATION"                                                 
    ##  [35] "GOBP_POSTSYNAPTIC_MEMBRANE_ORGANIZATION"                                       
    ##  [36] "GOBP_AXON_ENSHEATHMENT_IN_CENTRAL_NERVOUS_SYSTEM"                              
    ##  [37] "GOBP_PRESYNAPSE_ORGANIZATION"                                                  
    ##  [38] "GOBP_REGULATION_OF_SYNAPTIC_PLASTICITY"                                        
    ##  [39] "GOBP_AXON_EXTENSION"                                                           
    ##  [40] "GOBP_NEGATIVE_REGULATION_OF_NEURON_APOPTOTIC_PROCESS"                          
    ##  [41] "GOBP_POSITIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION"                             
    ##  [42] "GOBP_NEURON_FATE_SPECIFICATION"                                                
    ##  [43] "GOBP_POSITIVE_REGULATION_OF_NEURON_APOPTOTIC_PROCESS"                          
    ##  [44] "GOBP_SYNAPTIC_VESICLE_LOCALIZATION"                                            
    ##  [45] "GOBP_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"                                      
    ##  [46] "GOBP_NEURON_FATE_COMMITMENT"                                                   
    ##  [47] "GOBP_NEURON_PROJECTION_REGENERATION"                                           
    ##  [48] "GOBP_MIDBRAIN_DOPAMINERGIC_NEURON_DIFFERENTIATION"                             
    ##  [49] "GOBP_NEGATIVE_REGULATION_OF_NEURON_DEATH"                                      
    ##  [50] "GOBP_POSITIVE_REGULATION_OF_NEUROBLAST_PROLIFERATION"                          
    ##  [51] "GOBP_REGULATION_OF_NEUROBLAST_PROLIFERATION"                                   
    ##  [52] "GOBP_REGULATION_OF_NEUROTRANSMITTER_RECEPTOR_ACTIVITY"                         
    ##  [53] "GOBP_REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY"                              
    ##  [54] "GOBP_SYNAPTIC_MEMBRANE_ADHESION"                                               
    ##  [55] "GOBP_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_POTENTIAL"                            
    ##  [56] "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC"                        
    ##  [57] "GOBP_FOREBRAIN_GENERATION_OF_NEURONS"                                          
    ##  [58] "GOBP_FOREBRAIN_NEURON_DIFFERENTIATION"                                         
    ##  [59] "GOBP_NEURON_DEATH_IN_RESPONSE_TO_OXIDATIVE_STRESS"                             
    ##  [60] "GOBP_PROTEIN_LOCALIZATION_TO_SYNAPSE"                                          
    ##  [61] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DEVELOPMENT"                                
    ##  [62] "GOBP_NEURON_PROJECTION_EXTENSION"                                              
    ##  [63] "GOBP_FACULTATIVE_HETEROCHROMATIN_ASSEMBLY"                                     
    ##  [64] "GOBP_HETEROCHROMATIN_ORGANIZATION"                                             
    ##  [65] "GOBP_NEGATIVE_REGULATION_OF_OXIDATIVE_STRESS_INDUCED_NEURON_DEATH"             
    ##  [66] "GOBP_SYNAPSE_ASSEMBLY"                                                         
    ##  [67] "GOBP_NEUROTRANSMITTER_TRANSPORT"                                               
    ##  [68] "GOBP_AXONAL_TRANSPORT"                                                         
    ##  [69] "GOBP_POSITIVE_REGULATION_OF_CHROMATIN_ORGANIZATION"                            
    ##  [70] "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_AXONOGENESIS"                               
    ##  [71] "GOBP_REGULATION_OF_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY"                          
    ##  [72] "GOBP_INHIBITORY_SYNAPSE_ASSEMBLY"                                              
    ##  [73] "GOBP_RESPONSE_TO_AXON_INJURY"                                                  
    ##  [74] "GOBP_SYNAPTIC_VESICLE_MEMBRANE_ORGANIZATION"                                   
    ##  [75] "GOBP_DNA_METHYLATION_DEPENDENT_HETEROCHROMATIN_ASSEMBLY"                       
    ##  [76] "GOBP_REGULATION_OF_CHROMATIN_ASSEMBLY"                                         
    ##  [77] "GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING"                                   
    ##  [78] "GOBP_SYNAPTIC_VESICLE_TRANSPORT"                                               
    ##  [79] "GOBP_REGULATION_OF_NEURON_PROJECTION_REGENERATION"                             
    ##  [80] "GOBP_CALCIUM_ION_REGULATED_EXOCYTOSIS_OF_NEUROTRANSMITTER"                     
    ##  [81] "GOBP_DOPAMINERGIC_NEURON_DIFFERENTIATION"                                      
    ##  [82] "GOBP_POSITIVE_REGULATION_OF_NEUROTRANSMITTER_SECRETION"                        
    ##  [83] "GOBP_REGULATION_OF_AXONOGENESIS"                                               
    ##  [84] "GOBP_REGULATION_OF_NEURON_DIFFERENTIATION"                                     
    ##  [85] "GOBP_NEUROBLAST_PROLIFERATION"                                                 
    ##  [86] "GOBP_NEURON_CELLULAR_HOMEOSTASIS"                                              
    ##  [87] "GOBP_NEUROTRANSMITTER_UPTAKE"                                                  
    ##  [88] "GOBP_SYNAPTIC_VESICLE_PRIMING"                                                 
    ##  [89] "GOBP_NEUROTRANSMITTER_REUPTAKE"                                                
    ##  [90] "GOBP_GLUTAMATE_SECRETION"                                                      
    ##  [91] "GOBP_CHROMATIN_DISASSEMBLY"                                                    
    ##  [92] "GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_POSTURE"                                
    ##  [93] "GOBP_NEURON_RECOGNITION"                                                       
    ##  [94] "GOBP_NEUROMUSCULAR_SYNAPTIC_TRANSMISSION"                                      
    ##  [95] "GOBP_REGULATION_OF_LONG_TERM_SYNAPTIC_POTENTIATION"                            
    ##  [96] "GOBP_AXONAL_TRANSPORT_OF_MITOCHONDRION"                                        
    ##  [97] "GOBP_CEREBRAL_CORTEX_NEURON_DIFFERENTIATION"                                   
    ##  [98] "GOBP_POSITIVE_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT"                        
    ##  [99] "GOBP_RETINAL_GANGLION_CELL_AXON_GUIDANCE"                                      
    ## [100] "GOBP_ANTEROGRADE_AXONAL_TRANSPORT"                                             
    ## [101] "GOBP_NEURON_MIGRATION"                                                         
    ## [102] "GOBP_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT"                              
    ## [103] "GOBP_REGULATION_OF_NEUROTRANSMITTER_UPTAKE"                                    
    ## [104] "GOBP_NEURON_MATURATION"                                                        
    ## [105] "GOBP_REGULATION_OF_GLUTAMATE_SECRETION"                                        
    ## [106] "GOBP_AXON_DEVELOPMENT"                                                         
    ## [107] "GOBP_RNA_MEDIATED_GENE_SILENCING_BY_INHIBITION_OF_TRANSLATION"                 
    ## [108] "GOBP_POSITIVE_REGULATION_OF_NEURON_DIFFERENTIATION"                            
    ## [109] "GOBP_REGULATION_OF_NEURON_MIGRATION"                                           
    ## [110] "GOBP_REGULATION_OF_PRESYNAPSE_ORGANIZATION"                                    
    ## [111] "GOBP_NEUROMUSCULAR_PROCESS_CONTROLLING_BALANCE"                                
    ## [112] "GOBP_NEURON_PROJECTION_ARBORIZATION"                                           
    ## [113] "GOBP_POSITIVE_REGULATION_OF_NEURON_MIGRATION"                                  
    ## [114] "GOBP_RESPONSE_TO_LEUKEMIA_INHIBITORY_FACTOR"                                   
    ## [115] "GOBP_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT"                                 
    ## [116] "GOBP_NEUROTROPHIN_SIGNALING_PATHWAY"                                           
    ## [117] "GOBP_NEUROTRANSMITTER_RECEPTOR_INTERNALIZATION"                                
    ## [118] "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION"                                  
    ## [119] "GOBP_NEURONAL_STEM_CELL_POPULATION_MAINTENANCE"                                
    ## [120] "GOBP_NEUROTRANSMITTER_RECEPTOR_TRANSPORT"                                      
    ## [121] "GOBP_REGULATION_OF_PROTEIN_LOCALIZATION_TO_SYNAPSE"                            
    ## [122] "GOBP_REGULATION_OF_RECEPTOR_LOCALIZATION_TO_SYNAPSE"                           
    ## [123] "GOBP_FOREBRAIN_NEURON_DEVELOPMENT"                                             
    ## [124] "GOBP_RECEPTOR_LOCALIZATION_TO_SYNAPSE"                                         
    ## [125] "GOBP_SYNAPTIC_VESICLE_RECYCLING"                                               
    ## [126] "GOBP_NEGATIVE_REGULATION_OF_NEURON_DIFFERENTIATION"                            
    ## [127] "GOBP_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_NEUROTRANSMITTER_RECEPTOR_LEVELS"     
    ## [128] "GOBP_AXONEMAL_DYNEIN_COMPLEX_ASSEMBLY"                                         
    ## [129] "GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS"                                              
    ## [130] "GOBP_SYNAPTIC_VESICLE_CYTOSKELETAL_TRANSPORT"                                  
    ## [131] "GOBP_PRESYNAPTIC_ENDOCYTOSIS"                                                  
    ## [132] "KEGG_NEUROTROPHIN_SIGNALING_PATHWAY"                                           
    ## [133] "GOBP_ADENYLATE_CYCLASE_INHIBITING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY"
    ## [134] "GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION"                                    
    ## [135] "GOBP_PROTEIN_LOCALIZATION_TO_POSTSYNAPSE"                                      
    ## [136] "GOBP_AXONEME_ASSEMBLY"                                                         
    ## [137] "GOBP_AXONAL_FASCICULATION"                                                     
    ## [138] "GOBP_NEUROTROPHIN_TRK_RECEPTOR_SIGNALING_PATHWAY"                              
    ## [139] "GOBP_REGULATION_OF_CHROMATIN_BINDING"                                          
    ## [140] "GOBP_PROTEIN_LOCALIZATION_TO_CHROMATIN"                                        
    ## [141] "GOBP_POSTSYNAPTIC_SIGNAL_TRANSDUCTION"                                         
    ## [142] "GOBP_REGULATION_OF_DNA_METHYLATION_DEPENDENT_HETEROCHROMATIN_ASSEMBLY"

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

    ## # A tibble: 19 × 3
    ##    pathway                                                   rare_NES common_NES
    ##    <chr>                                                        <dbl>      <dbl>
    ##  1 GOBP_NEURON_DEATH                                             2.24       1.34
    ##  2 GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS                      2.07       1.47
    ##  3 GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH                      2.01       1.54
    ##  4 GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT     2.38       1.40
    ##  5 GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS                      1.81       1.44
    ##  6 GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE                    1.65       1.45
    ##  7 GOBP_CHROMATIN_REMODELING                                     2.26       1.63
    ##  8 GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS                    1.60       1.39
    ##  9 GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION            1.83       1.51
    ## 10 GOBP_NEURON_PROJECTION_GUIDANCE                               1.57       1.40
    ## 11 KEGG_AXON_GUIDANCE                                            1.65       1.39
    ## 12 GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY                        1.53       1.46
    ## 13 GOBP_SYNAPSE_ORGANIZATION                                     2.14       1.41
    ## 14 GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJ…     1.80       1.83
    ## 15 GOBP_REGULATION_OF_NEUROGENESIS                               2.17       1.37
    ## 16 GOBP_NEUROTRANSMITTER_SECRETION                               1.50       1.50
    ## 17 GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENES…     1.67       1.60
    ## 18 GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_A…     1.47       1.64
    ## 19 KEGG_MTOR_SIGNALING_PATHWAY                                   1.38       1.61

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
    ##  1 GOBP_NEURON_DEATH           2.24       1.34  0.901      0.382           0.259
    ##  2 GOBP_POSITIVE_REGULATI…     2.07       1.47  0.606      0.382           0.259
    ##  3 GOBP_POSITIVE_REGULATI…     2.01       1.54  0.478      0.382           0.259
    ##  4 GOBP_NEGATIVE_REGULATI…     2.38       1.40  0.982      0.382           0.259
    ##  5 GOBP_POSITIVE_REGULATI…     1.81       1.44  0.363      0.382           0.259
    ##  6 GOBP_VESICLE_MEDIATED_…     1.65       1.45  0.194      0.382           0.259
    ##  7 GOBP_CHROMATIN_REMODEL…     2.26       1.63  0.624      0.382           0.259
    ##  8 GOBP_REGULATION_OF_NEU…     1.60       1.39  0.208      0.382           0.259
    ##  9 GOBP_CENTRAL_NERVOUS_S…     1.83       1.51  0.321      0.382           0.259
    ## 10 GOBP_NEURON_PROJECTION…     1.57       1.40  0.175      0.382           0.259
    ## 11 KEGG_AXON_GUIDANCE          1.65       1.39  0.259      0.382           0.259
    ## 12 GOBP_CHROMATIN_ASSEMBL…     1.53       1.46  0.0668     0.382           0.259
    ## 13 GOBP_SYNAPSE_ORGANIZAT…     2.14       1.41  0.732      0.382           0.259
    ## 14 GOBP_NEURON_PROJECTION…     1.80       1.83 -0.0357     0.382           0.259
    ## 15 GOBP_REGULATION_OF_NEU…     2.17       1.37  0.797      0.382           0.259
    ## 16 GOBP_NEUROTRANSMITTER_…     1.50       1.50  0.00154    0.382           0.259
    ## 17 GOBP_CENTRAL_NERVOUS_S…     1.67       1.60  0.0755     0.382           0.259
    ## 18 GOBP_NEGATIVE_REGULATI…     1.47       1.64 -0.174      0.382           0.259
    ## 19 KEGG_MTOR_SIGNALING_PA…     1.38       1.61 -0.230      0.382           0.259

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
    ##   pathway                                                   NES_diff
    ##   <chr>                                                        <dbl>
    ## 1 GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT    0.982

``` r
most_convergent
```

    ## # A tibble: 1 × 2
    ##   pathway                         NES_diff
    ##   <chr>                              <dbl>
    ## 1 GOBP_NEUROTRANSMITTER_SECRETION  0.00154

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

    ## # A tibble: 19 × 7
    ##    pathway                 rare_NES common_NES NES_diff MAD_filt median_NES_diff
    ##    <chr>                      <dbl>      <dbl>    <dbl>    <dbl>           <dbl>
    ##  1 GOBP_NEURON_DEATH           2.24       1.34  0.901      0.382           0.259
    ##  2 GOBP_POSITIVE_REGULATI…     2.07       1.47  0.606      0.382           0.259
    ##  3 GOBP_POSITIVE_REGULATI…     2.01       1.54  0.478      0.382           0.259
    ##  4 GOBP_NEGATIVE_REGULATI…     2.38       1.40  0.982      0.382           0.259
    ##  5 GOBP_POSITIVE_REGULATI…     1.81       1.44  0.363      0.382           0.259
    ##  6 GOBP_VESICLE_MEDIATED_…     1.65       1.45  0.194      0.382           0.259
    ##  7 GOBP_CHROMATIN_REMODEL…     2.26       1.63  0.624      0.382           0.259
    ##  8 GOBP_REGULATION_OF_NEU…     1.60       1.39  0.208      0.382           0.259
    ##  9 GOBP_CENTRAL_NERVOUS_S…     1.83       1.51  0.321      0.382           0.259
    ## 10 GOBP_NEURON_PROJECTION…     1.57       1.40  0.175      0.382           0.259
    ## 11 KEGG_AXON_GUIDANCE          1.65       1.39  0.259      0.382           0.259
    ## 12 GOBP_CHROMATIN_ASSEMBL…     1.53       1.46  0.0668     0.382           0.259
    ## 13 GOBP_SYNAPSE_ORGANIZAT…     2.14       1.41  0.732      0.382           0.259
    ## 14 GOBP_NEURON_PROJECTION…     1.80       1.83 -0.0357     0.382           0.259
    ## 15 GOBP_REGULATION_OF_NEU…     2.17       1.37  0.797      0.382           0.259
    ## 16 GOBP_NEUROTRANSMITTER_…     1.50       1.50  0.00154    0.382           0.259
    ## 17 GOBP_CENTRAL_NERVOUS_S…     1.67       1.60  0.0755     0.382           0.259
    ## 18 GOBP_NEGATIVE_REGULATI…     1.47       1.64 -0.174      0.382           0.259
    ## 19 KEGG_MTOR_SIGNALING_PA…     1.38       1.61 -0.230      0.382           0.259
    ## # ℹ 1 more variable: converge_diverge <chr>

``` r
divergent <- classify_paths_filt %>% filter(converge_diverge == "divergent")
divergent
```

    ## # A tibble: 10 × 7
    ##    pathway                 rare_NES common_NES NES_diff MAD_filt median_NES_diff
    ##    <chr>                      <dbl>      <dbl>    <dbl>    <dbl>           <dbl>
    ##  1 GOBP_NEURON_DEATH           2.24       1.34    0.901    0.382           0.259
    ##  2 GOBP_POSITIVE_REGULATI…     2.07       1.47    0.606    0.382           0.259
    ##  3 GOBP_POSITIVE_REGULATI…     2.01       1.54    0.478    0.382           0.259
    ##  4 GOBP_NEGATIVE_REGULATI…     2.38       1.40    0.982    0.382           0.259
    ##  5 GOBP_POSITIVE_REGULATI…     1.81       1.44    0.363    0.382           0.259
    ##  6 GOBP_CHROMATIN_REMODEL…     2.26       1.63    0.624    0.382           0.259
    ##  7 GOBP_SYNAPSE_ORGANIZAT…     2.14       1.41    0.732    0.382           0.259
    ##  8 GOBP_REGULATION_OF_NEU…     2.17       1.37    0.797    0.382           0.259
    ##  9 GOBP_NEGATIVE_REGULATI…     1.47       1.64   -0.174    0.382           0.259
    ## 10 KEGG_MTOR_SIGNALING_PA…     1.38       1.61   -0.230    0.382           0.259
    ## # ℹ 1 more variable: converge_diverge <chr>

``` r
nrow(divergent) # 8
```

    ## [1] 10

``` r
convergent <- classify_paths_filt %>% filter(converge_diverge == "convergent")
convergent
```

    ## # A tibble: 9 × 7
    ##   pathway rare_NES common_NES NES_diff MAD_filt median_NES_diff converge_diverge
    ##   <chr>      <dbl>      <dbl>    <dbl>    <dbl>           <dbl> <chr>           
    ## 1 GOBP_V…     1.65       1.45  0.194      0.382           0.259 convergent      
    ## 2 GOBP_R…     1.60       1.39  0.208      0.382           0.259 convergent      
    ## 3 GOBP_C…     1.83       1.51  0.321      0.382           0.259 convergent      
    ## 4 GOBP_N…     1.57       1.40  0.175      0.382           0.259 convergent      
    ## 5 KEGG_A…     1.65       1.39  0.259      0.382           0.259 convergent      
    ## 6 GOBP_C…     1.53       1.46  0.0668     0.382           0.259 convergent      
    ## 7 GOBP_N…     1.80       1.83 -0.0357     0.382           0.259 convergent      
    ## 8 GOBP_N…     1.50       1.50  0.00154    0.382           0.259 convergent      
    ## 9 GOBP_C…     1.67       1.60  0.0755     0.382           0.259 convergent

``` r
nrow(convergent) # 7
```

    ## [1] 9

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
    ##  2:                                GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS
    ##  3:                                GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH
    ##  4:               GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT
    ##  5:                                GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS
    ##  6:                              GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE
    ##  7:                                               GOBP_CHROMATIN_REMODELING
    ##  8:                              GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS
    ##  9:                      GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION
    ## 10:                                         GOBP_NEURON_PROJECTION_GUIDANCE
    ## 11:                                                      KEGG_AXON_GUIDANCE
    ## 12:                                  GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY
    ## 13:                                               GOBP_SYNAPSE_ORGANIZATION
    ## 14: GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE
    ## 15:                                         GOBP_REGULATION_OF_NEUROGENESIS
    ## 16:                                         GOBP_NEUROTRANSMITTER_SECRETION
    ## 17:              GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS
    ## 18:    GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE
    ## 19:                                             KEGG_MTOR_SIGNALING_PATHWAY
    ## 20:                                               GOBP_CHROMATIN_REMODELING
    ## 21: GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE
    ## 22:                                               GOBP_SYNAPSE_ORGANIZATION
    ## 23:                                  GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY
    ## 24:                      GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION
    ## 25:                                         GOBP_REGULATION_OF_NEUROGENESIS
    ## 26:                                             KEGG_MTOR_SIGNALING_PATHWAY
    ## 27:                                GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS
    ## 28:                                         GOBP_NEUROTRANSMITTER_SECRETION
    ## 29:                                GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH
    ## 30:    GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE
    ## 31:                              GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE
    ## 32:                              GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS
    ## 33:                                                       GOBP_NEURON_DEATH
    ## 34:                                GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS
    ## 35:                                         GOBP_NEURON_PROJECTION_GUIDANCE
    ## 36:                                                      KEGG_AXON_GUIDANCE
    ## 37:              GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS
    ## 38:               GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT
    ##                                                                     pathway
    ##             pval       padj   log2err        ES      NES  size  leadingEdge
    ##            <num>      <num>     <num>     <num>    <num> <int>       <list>
    ##  1: 1.058477e-03 0.08975834 0.4550599 0.9595688 2.242748   325 6326, 88....
    ##  2: 9.191707e-03 0.08975834 0.3807304 0.8138417 2.072072    66 23394, 1....
    ##  3: 1.877458e-02 0.08975834 0.3524879 0.7955304 2.013173    80 2904, 14....
    ##  4: 2.418156e-02 0.08975834 0.3524879 0.9651832 2.376862   123 8831, 57....
    ##  5: 3.196803e-02 0.08975834 0.2529611 0.7457535 1.805221   199 23394, 8....
    ##  6: 3.296703e-02 0.08975834 0.2489111 0.6767406 1.647688   183 5728, 14....
    ##  7: 3.796204e-02 0.08975834 0.2311267 0.9450217 2.256908   224 57680, 2....
    ##  8: 3.896104e-02 0.08975834 0.2279872 0.6611369 1.600392   199 6529, 93....
    ##  9: 4.256451e-02 0.08986416 0.3217759 0.7493933 1.827047   157 5728, 85....
    ## 10: 4.295704e-02 0.08986416 0.2165428 0.6573949 1.573610   217 1826, 93....
    ## 11: 4.295704e-02 0.08986416 0.2165428 0.6707905 1.652252   121 1808, 27....
    ## 12: 4.895105e-02 0.09184630 0.2020717 0.6311725 1.530832   196 11011, 1....
    ## 13: 4.895105e-02 0.09184630 0.2020717 0.9249977 2.143933   370 8831, 23....
    ## 14: 5.494505e-02 0.09259259 0.1900233 0.6973979 1.798962    33 1826, 65....
    ## 15: 5.494505e-02 0.09259259 0.1900233 0.9269535 2.166674   326 8831, 23....
    ## 16: 7.192807e-02 0.10643016 0.1644058 0.6095070 1.497732   136 9378, 68....
    ## 17: 7.592408e-02 0.10966811 0.1596467 0.6443596 1.673464    25 2047, 47....
    ## 18: 1.328671e-01 0.16677116 0.1167392 0.5653401 1.466556    24 6585, 64....
    ## 19: 1.468531e-01 0.17700181 0.1101223 0.5405839 1.384413    48 208, 253....
    ## 20: 8.317188e-05 0.01530363 0.5384341 0.3704143 1.632745   227 9759, 21....
    ## 21: 6.746814e-04 0.06207069 0.4772708 0.5462061 1.834620    36 7473, 91....
    ## 22: 1.701939e-03 0.10438556 0.4550599 0.3106873 1.411983   384 4137, 11....
    ## 23: 5.025075e-03 0.12283274 0.4070179 0.3368143 1.464058   191 3010, 83....
    ## 24: 5.031741e-03 0.12283274 0.4070179 0.3530124 1.506449   160 4137, 91....
    ## 25: 5.084028e-03 0.12283274 0.4070179 0.3032521 1.369759   338 4137, 74....
    ## 26: 5.866190e-03 0.12283274 0.4070179 0.4540769 1.613946    49 2475, 60....
    ## 27: 6.008123e-03 0.12283274 0.4070179 0.3294199 1.441839   208 4137, 74....
    ## 28: 8.632143e-03 0.14212747 0.3807304 0.3592682 1.496196   137 783, 990....
    ## 29: 9.672977e-03 0.14212747 0.3807304 0.3951873 1.535365    81 4137, 67....
    ## 30: 1.003452e-02 0.14212747 0.3807304 0.5234268 1.640135    26 7473, 10....
    ## 31: 1.012200e-02 0.14212747 0.3807304 0.3350995 1.453844   186 783, 990....
    ## 32: 1.081405e-02 0.14212747 0.3807304 0.3190564 1.392687   201 783, 990....
    ## 33: 1.246729e-02 0.15048755 0.3807304 0.2979639 1.341441   327 4137, 67....
    ## 34: 1.379969e-02 0.15048755 0.3807304 0.3855854 1.465871    71 4137, 74....
    ## 35: 1.475114e-02 0.15048755 0.3807304 0.3180250 1.398423   221 7473, 91....
    ## 36: 1.543152e-02 0.15048755 0.3807304 0.3384993 1.392806   121 5879, 27....
    ## 37: 1.553948e-02 0.15048755 0.3807304 0.5163934 1.597993    25 91584, 1....
    ## 38: 2.325128e-02 0.21391181 0.3524879 0.3374266 1.395005   128 7473, 23....
    ##             pval       padj   log2err        ES      NES  size  leadingEdge
    ##     rare_NES common_NES     NES_diff  MAD_filt median_NES_diff converge_diverge
    ##        <num>      <num>        <num>     <num>           <num>           <char>
    ##  1: 2.242748   1.341441  0.901307416 0.3823771       0.2594461        divergent
    ##  2: 2.072072   1.465871  0.606201356 0.3823771       0.2594461        divergent
    ##  3: 2.013173   1.535365  0.477808229 0.3823771       0.2594461        divergent
    ##  4: 2.376862   1.395005  0.981856834 0.3823771       0.2594461        divergent
    ##  5: 1.805221   1.441839  0.363381688 0.3823771       0.2594461        divergent
    ##  6: 1.647688   1.453844  0.193844052 0.3823771       0.2594461       convergent
    ##  7: 2.256908   1.632745  0.624163869 0.3823771       0.2594461        divergent
    ##  8: 1.600392   1.392687  0.207705384 0.3823771       0.2594461       convergent
    ##  9: 1.827047   1.506449  0.320598400 0.3823771       0.2594461       convergent
    ## 10: 1.573610   1.398423  0.175187000 0.3823771       0.2594461       convergent
    ## 11: 1.652252   1.392806  0.259446099 0.3823771       0.2594461       convergent
    ## 12: 1.530832   1.464058  0.066774651 0.3823771       0.2594461       convergent
    ## 13: 2.143933   1.411983  0.731950357 0.3823771       0.2594461        divergent
    ## 14: 1.798962   1.834620 -0.035658630 0.3823771       0.2594461       convergent
    ## 15: 2.166674   1.369759  0.796914575 0.3823771       0.2594461        divergent
    ## 16: 1.497732   1.496196  0.001536302 0.3823771       0.2594461       convergent
    ## 17: 1.673464   1.597993  0.075471051 0.3823771       0.2594461       convergent
    ## 18: 1.466556   1.640135 -0.173578786 0.3823771       0.2594461        divergent
    ## 19: 1.384413   1.613946 -0.229532793 0.3823771       0.2594461        divergent
    ## 20: 2.256908   1.632745  0.624163869 0.3823771       0.2594461        divergent
    ## 21: 1.798962   1.834620 -0.035658630 0.3823771       0.2594461       convergent
    ## 22: 2.143933   1.411983  0.731950357 0.3823771       0.2594461        divergent
    ## 23: 1.530832   1.464058  0.066774651 0.3823771       0.2594461       convergent
    ## 24: 1.827047   1.506449  0.320598400 0.3823771       0.2594461       convergent
    ## 25: 2.166674   1.369759  0.796914575 0.3823771       0.2594461        divergent
    ## 26: 1.384413   1.613946 -0.229532793 0.3823771       0.2594461        divergent
    ## 27: 1.805221   1.441839  0.363381688 0.3823771       0.2594461        divergent
    ## 28: 1.497732   1.496196  0.001536302 0.3823771       0.2594461       convergent
    ## 29: 2.013173   1.535365  0.477808229 0.3823771       0.2594461        divergent
    ## 30: 1.466556   1.640135 -0.173578786 0.3823771       0.2594461        divergent
    ## 31: 1.647688   1.453844  0.193844052 0.3823771       0.2594461       convergent
    ## 32: 1.600392   1.392687  0.207705384 0.3823771       0.2594461       convergent
    ## 33: 2.242748   1.341441  0.901307416 0.3823771       0.2594461        divergent
    ## 34: 2.072072   1.465871  0.606201356 0.3823771       0.2594461        divergent
    ## 35: 1.573610   1.398423  0.175187000 0.3823771       0.2594461       convergent
    ## 36: 1.652252   1.392806  0.259446099 0.3823771       0.2594461       convergent
    ## 37: 1.673464   1.597993  0.075471051 0.3823771       0.2594461       convergent
    ## 38: 2.376862   1.395005  0.981856834 0.3823771       0.2594461        divergent
    ##     rare_NES common_NES     NES_diff  MAD_filt median_NES_diff converge_diverge

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
