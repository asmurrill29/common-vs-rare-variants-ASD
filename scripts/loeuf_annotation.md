LOEUF Annotation and Gene-Level Constraint Assessment
================
Amalya Murrill
2026-04-22

This script takes the enriched pathways from FGSEA that have been
statistically categorized as convergent or divergent and extracts
gene-level data. The findings are joined with corresponding LOEUF scores
downloaded from gnomAD v2.1.1. Overall constraint for each pathway was
assessed by calculating the proportion and then percentage of
constrained genes (LOEUF \< 0.35 as dictated by literature) within that
pathway. The data was then formatted for Cytoscape output as a
visualization of constraint within a network.

``` r
knitr::opts_chunk$set(warning = FALSE, message = FALSE)
```

## Libraries

``` r
library(tidyverse)
library(org.Hs.eg.db)
```

## Import FGSEA Results

``` r
#import leadingEdge data per pathway per variant
paths_data <- read.csv("../../results/all_path_data.csv", header = TRUE)

paths_data <- paths_data %>% dplyr::select(variant, pathway, padj,leadingEdge)
head(paths_data)
```

    ##   variant                                                   pathway       padj
    ## 1    rare                                         GOBP_NEURON_DEATH 0.08450704
    ## 2    rare                  GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH 0.08450704
    ## 3    rare GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT 0.08450704
    ## 4    rare        GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION 0.08450704
    ## 5    rare                  GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS 0.08450704
    ## 6    rare                  GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS 0.09208973
    ##                                                                                                                                                                      leadingEdge
    ## 1                                                                                                   6326;8831;23394;2904;1499;6812;2561;11151;2047;4763;7068;5609;3326;2932;7314
    ## 2                                                                                                                     2904;1499;4763;5609;2932;6733;2213;468;2353;2896;2901;5621
    ## 3                                                                                                                                                           8831;5728;53335;2670
    ## 4                                         5728;85358;25942;10716;8861;2047;3326;7314;23314;51684;2560;4781;2048;4929;388585;3320;55079;64843;9693;5881;6095;79625;2909;1630;4133
    ## 5                                                                                                                                 23394;1826;22902;5362;5361;5604;6794;3984;5364
    ## 6 23394;85358;1499;1826;2670;22902;5649;5362;3706;23316;2048;5501;5789;5361;4076;5604;6794;6663;2475;3984;5364;23405;10097;7424;7520;3091;5803;3815;7074;10507;27020;23654;81565

``` r
nrow(paths_data) # 38; 19 per variant
```

    ## [1] 38

``` r
#import pathway classification 
classify_paths_results <- read.csv("../../results/path_results.csv", header = TRUE) 
classify_paths_results <- classify_paths_results %>% dplyr::select(pathway, converge_diverge)
nrow(classify_paths_results) # 19
```

    ## [1] 19

``` r
head(classify_paths_results)
```

    ##                                                     pathway converge_diverge
    ## 1                                         GOBP_NEURON_DEATH        divergent
    ## 2                  GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH       convergent
    ## 3 GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT        divergent
    ## 4        GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION       convergent
    ## 5                  GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS       convergent
    ## 6                  GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS       convergent

## Join Gene Data with Convergence Data

``` r
classify_genes <- left_join(paths_data, classify_paths_results, by = "pathway")
head(classify_genes)
```

    ##   variant                                                   pathway       padj
    ## 1    rare                                         GOBP_NEURON_DEATH 0.08450704
    ## 2    rare                  GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH 0.08450704
    ## 3    rare GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT 0.08450704
    ## 4    rare        GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION 0.08450704
    ## 5    rare                  GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS 0.08450704
    ## 6    rare                  GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS 0.09208973
    ##                                                                                                                                                                      leadingEdge
    ## 1                                                                                                   6326;8831;23394;2904;1499;6812;2561;11151;2047;4763;7068;5609;3326;2932;7314
    ## 2                                                                                                                     2904;1499;4763;5609;2932;6733;2213;468;2353;2896;2901;5621
    ## 3                                                                                                                                                           8831;5728;53335;2670
    ## 4                                         5728;85358;25942;10716;8861;2047;3326;7314;23314;51684;2560;4781;2048;4929;388585;3320;55079;64843;9693;5881;6095;79625;2909;1630;4133
    ## 5                                                                                                                                 23394;1826;22902;5362;5361;5604;6794;3984;5364
    ## 6 23394;85358;1499;1826;2670;22902;5649;5362;3706;23316;2048;5501;5789;5361;4076;5604;6794;6663;2475;3984;5364;23405;10097;7424;7520;3091;5803;3815;7074;10507;27020;23654;81565
    ##   converge_diverge
    ## 1        divergent
    ## 2       convergent
    ## 3        divergent
    ## 4       convergent
    ## 5       convergent
    ## 6       convergent

### Convert EntrezIDs to Gene Symbols

``` r
leadingEdgeGenes <- classify_genes %>% dplyr::select(leadingEdge) %>% pull()

#split each list into individual gene characters
leadingEdgeGenes_split <- strsplit(leadingEdgeGenes, ";")

#use lapply across mapping 
symbols <- lapply(leadingEdgeGenes_split, function(r) mapIds(
  org.Hs.eg.db,
  keys = r,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
))

#rename
symbols <- lapply(symbols, unname)
head(symbols)
```

    ## [[1]]
    ##  [1] "SCN2A"    "SYNGAP1"  "ADNP"     "GRIN2B"   "CTNNB1"   "STXBP1"  
    ##  [7] "GABRB2"   "CORO1A"   "EPHB1"    "NF1"      "THRB"     "MAP2K7"  
    ## [13] "HSP90AB1" "GSK3B"    "UBB"     
    ## 
    ## [[2]]
    ##  [1] "GRIN2B" "CTNNB1" "NF1"    "MAP2K7" "GSK3B"  "SRPK2"  "FCGR2B" "ATF4"  
    ##  [9] "FOS"    "GRN"    "GRIK5"  "PRNP"  
    ## 
    ## [[3]]
    ## [1] "SYNGAP1" "PTEN"    "BCL11A"  "GFAP"   
    ## 
    ## [[4]]
    ##  [1] "PTEN"     "SHANK3"   "SIN3A"    "TBR1"     "LDB1"     "EPHB1"   
    ##  [7] "HSP90AB1" "UBB"      "SATB2"    "SUFU"     "GABRB1"   "NFIB"    
    ## [13] "EPHB2"    "NR4A2"    "HES5"     "HSP90AA1" "FEZF2"    "ISL2"    
    ## [19] "RAPGEF2"  "RAC3"     "RORA"     "NDNF"     "ARHGAP35" "DCC"     
    ## [25] "MAP2"    
    ## 
    ## [[5]]
    ## [1] "ADNP"   "DSCAM"  "RUFY3"  "PLXNA2" "PLXNA1" "MAP2K1" "STK11"  "LIMK1" 
    ## [9] "PLXNB1"
    ## 
    ## [[6]]
    ##  [1] "ADNP"    "SHANK3"  "CTNNB1"  "DSCAM"   "GFAP"    "RUFY3"   "RELN"   
    ##  [8] "PLXNA2"  "ITPKA"   "CUX2"    "EPHB2"   "PPP1CC"  "PTPRD"   "PLXNA1" 
    ## [15] "CAPRIN1" "MAP2K1"  "STK11"   "SOX10"   "MTOR"    "LIMK1"   "PLXNB1" 
    ## [22] "DICER1"  "ACTR2"   "VEGFC"   "XRCC5"   "HIF1A"   "PTPRZ1"  "KIT"    
    ## [29] "TIAM1"   "SEMA4D"  "NPTN"    "PLXNB2"  "NDEL1"

``` r
#re-add to table
classify_genes <- classify_genes %>% mutate(leadingEdge = symbols)
head(classify_genes)
```

    ##   variant                                                   pathway       padj
    ## 1    rare                                         GOBP_NEURON_DEATH 0.08450704
    ## 2    rare                  GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH 0.08450704
    ## 3    rare GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT 0.08450704
    ## 4    rare        GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION 0.08450704
    ## 5    rare                  GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS 0.08450704
    ## 6    rare                  GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS 0.09208973
    ##    leadingEdge converge_diverge
    ## 1 SCN2A, S....        divergent
    ## 2 GRIN2B, ....       convergent
    ## 3 SYNGAP1,....        divergent
    ## 4 PTEN, SH....       convergent
    ## 5 ADNP, DS....       convergent
    ## 6 ADNP, SH....       convergent

## Import gnomAD Constraint Scores

``` r
gnomad <- read.csv("../../data/raw/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", header = TRUE, sep = "\t")
head(gnomad)
```

    ##    gene      transcript obs_mis exp_mis  oe_mis     mu_mis possible_mis
    ## 1 MED13 ENST00000397786     871 1117.80 0.77921 5.5598e-05        14195
    ## 2 NIPBL ENST00000282516     846 1441.50 0.58688 7.3808e-05        18540
    ## 3  SMC3 ENST00000361804     178  630.07 0.28251 3.2489e-05         8109
    ## 4 CNOT1 ENST00000317147     561 1295.90 0.43290 6.9116e-05        15670
    ## 5   RLF ENST00000372771     669  972.87 0.68766 4.7052e-05        12682
    ## 6 PCF11 ENST00000298281     574  783.46 0.73265 3.9067e-05        10106
    ##   obs_mis_pphen exp_mis_pphen oe_mis_pphen possible_mis_pphen obs_syn exp_syn
    ## 1           314        529.75      0.59273               6708     422  387.53
    ## 2           158        543.10      0.29092               7135     496  495.01
    ## 3            21        182.52      0.11506               2197     215  203.25
    ## 4            51        290.68      0.17545               3560     470  456.03
    ## 5           107        321.14      0.33319               4151     358  352.62
    ## 6           130        264.36      0.49176               3521     298  272.69
    ##   oe_syn     mu_syn possible_syn obs_lof     mu_lof possible_lof exp_lof pLI
    ## 1 1.0890 1.9097e-05         4248       0 4.9203e-06         1257  98.429   1
    ## 2 1.0020 2.4942e-05         5211       1 9.4214e-06         1781 150.320   1
    ## 3 1.0578 9.8016e-06         2091       0 4.5403e-06          937  79.490   1
    ## 4 1.0306 2.3979e-05         4564       1 6.8100e-06         1440 125.030   1
    ## 5 1.0153 1.6694e-05         3482       0 4.0155e-06         1024  73.222   1
    ## 6 1.0928 1.3262e-05         3026       0 4.6454e-06          862  74.559   1
    ##        pNull       pRec    oe_lof oe_syn_lower oe_syn_upper oe_mis_lower
    ## 1 8.9436e-40 1.8383e-16 0.0000000        1.005        1.180        0.736
    ## 2 2.9773e-59 3.5724e-24 0.0066527        0.930        1.079        0.554
    ## 3 2.7853e-32 2.1914e-13 0.0000000        0.946        1.184        0.249
    ## 4 2.9924e-49 4.5629e-20 0.0079978        0.955        1.112        0.403
    ## 5 8.4055e-30 2.2842e-12 0.0000000        0.930        1.108        0.645
    ## 6 2.4872e-30 1.3855e-12 0.0000000        0.993        1.203        0.683
    ##   oe_mis_upper oe_lof_lower oe_lof_upper constraint_flag     syn_z  mis_z
    ## 1        0.824        0.000        0.030                 -1.376500 2.6232
    ## 2        0.621        0.001        0.032                 -0.035119 5.5737
    ## 3        0.320        0.000        0.037                 -0.647760 6.3999
    ## 4        0.464        0.002        0.038                 -0.514100 7.2546
    ## 5        0.733        0.000        0.040                 -0.225180 3.4620
    ## 6        0.785        0.000        0.040                 -1.205000 2.6592
    ##     lof_z oe_lof_upper_rank oe_lof_upper_bin oe_lof_upper_bin_6 n_sites
    ## 1  9.1935                 0                0                  0       2
    ## 2 11.2860                 1                0                  0       2
    ## 3  8.2618                 2                0                  0       8
    ## 4 10.2790                 3                0                  0       5
    ## 5  7.9294                 4                0                  0       1
    ## 6  8.0014                 5                0                  0       4
    ##   classic_caf     max_af no_lofs obs_het_lof obs_hom_lof defined          p
    ## 1  1.2058e-05 8.0492e-06  124782           3           0  124785 1.2021e-05
    ## 2  1.1943e-05 7.9636e-06  125693           3           0  125696 1.1934e-05
    ## 3  3.1885e-05 3.9986e-06  125731           8           0  125739 3.1812e-05
    ## 4  1.9952e-05 4.0200e-06  125740           4           0  125744 1.5905e-05
    ## 5  3.9961e-06 3.9961e-06  125122           1           0  125123 3.9961e-06
    ## 6  1.6123e-05 4.0534e-06  124625           3           0  124628 1.2036e-05
    ##   exp_hom_lof classic_caf_afr classic_caf_amr classic_caf_asj classic_caf_eas
    ## 1  1.8031e-05      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00
    ## 2  1.7901e-05      0.0000e+00      0.0000e+00      9.9246e-05      0.0000e+00
    ## 3  1.2725e-04      0.0000e+00      0.0000e+00      9.9364e-05      5.4366e-05
    ## 4  3.1811e-05      0.0000e+00      2.8915e-05      0.0000e+00      5.4645e-05
    ## 5  1.9980e-06      6.1866e-05      0.0000e+00      0.0000e+00      0.0000e+00
    ## 6  1.8054e-05      0.0000e+00      0.0000e+00      0.0000e+00      1.1127e-04
    ##   classic_caf_fin classic_caf_nfe classic_caf_oth classic_caf_sas      p_afr
    ## 1      9.2812e-05      8.8571e-06               0      0.0000e+00 0.0000e+00
    ## 2      0.0000e+00      0.0000e+00               0      6.5338e-05 0.0000e+00
    ## 3      0.0000e+00      4.4068e-05               0      3.2673e-05 0.0000e+00
    ## 4      0.0000e+00      2.6416e-05               0      0.0000e+00 0.0000e+00
    ## 5      0.0000e+00      0.0000e+00               0      0.0000e+00 6.1868e-05
    ## 6      0.0000e+00      1.7840e-05               0      0.0000e+00 0.0000e+00
    ##        p_amr      p_asj      p_eas     p_fin      p_nfe p_oth      p_sas
    ## 1 0.0000e+00 0.0000e+00 0.0000e+00 9.276e-05 8.8276e-06     0 0.0000e+00
    ## 2 0.0000e+00 9.9231e-05 0.0000e+00 0.000e+00 0.0000e+00     0 6.5327e-05
    ## 3 0.0000e+00 9.9211e-05 5.4367e-05 0.000e+00 4.3956e-05     0 3.2663e-05
    ## 4 2.8909e-05 0.0000e+00 5.4367e-05 0.000e+00 1.7581e-05     0 0.0000e+00
    ## 5 0.0000e+00 0.0000e+00 0.0000e+00 0.000e+00 0.0000e+00     0 0.0000e+00
    ## 6 0.0000e+00 0.0000e+00 5.5631e-05 0.000e+00 1.7698e-05     0 0.0000e+00
    ##   transcript_type         gene_id transcript_level cds_length num_coding_exons
    ## 1  protein_coding ENSG00000108510                2       6522               30
    ## 2  protein_coding ENSG00000164190                2       8412               46
    ## 3  protein_coding ENSG00000108055                2       3651               29
    ## 4  protein_coding ENSG00000125107                2       7128               48
    ## 5  protein_coding ENSG00000117000                2       5742                8
    ## 6  protein_coding ENSG00000165494                2       4665               16
    ##        gene_type gene_length exac_pLI exac_obs_lof exac_exp_lof exac_oe_lof
    ## 1 protein_coding      122678        1            0       64.393   0.0000000
    ## 2 protein_coding      189655        1            1      110.570   0.0090443
    ## 3 protein_coding       36946        1            0       58.523   0.0000000
    ## 4 protein_coding      109936        1            3       90.130   0.0332850
    ## 5 protein_coding       79549        1            0       43.607   0.0000000
    ## 6 protein_coding       30464        1            1       48.160   0.0207640
    ##   brain_expression chromosome start_position end_position
    ## 1               NA         17       60019966     60142643
    ## 2               NA          5       36876861     37066515
    ## 3               NA         10      112327449    112364394
    ## 4               NA         16       58553855     58663790
    ## 5               NA          1       40627045     40706593
    ## 6               NA         11       82868030     82898493

``` r
colnames(gnomad) #we want "gene" and "oe_lof_upper"
```

    ##  [1] "gene"               "transcript"         "obs_mis"           
    ##  [4] "exp_mis"            "oe_mis"             "mu_mis"            
    ##  [7] "possible_mis"       "obs_mis_pphen"      "exp_mis_pphen"     
    ## [10] "oe_mis_pphen"       "possible_mis_pphen" "obs_syn"           
    ## [13] "exp_syn"            "oe_syn"             "mu_syn"            
    ## [16] "possible_syn"       "obs_lof"            "mu_lof"            
    ## [19] "possible_lof"       "exp_lof"            "pLI"               
    ## [22] "pNull"              "pRec"               "oe_lof"            
    ## [25] "oe_syn_lower"       "oe_syn_upper"       "oe_mis_lower"      
    ## [28] "oe_mis_upper"       "oe_lof_lower"       "oe_lof_upper"      
    ## [31] "constraint_flag"    "syn_z"              "mis_z"             
    ## [34] "lof_z"              "oe_lof_upper_rank"  "oe_lof_upper_bin"  
    ## [37] "oe_lof_upper_bin_6" "n_sites"            "classic_caf"       
    ## [40] "max_af"             "no_lofs"            "obs_het_lof"       
    ## [43] "obs_hom_lof"        "defined"            "p"                 
    ## [46] "exp_hom_lof"        "classic_caf_afr"    "classic_caf_amr"   
    ## [49] "classic_caf_asj"    "classic_caf_eas"    "classic_caf_fin"   
    ## [52] "classic_caf_nfe"    "classic_caf_oth"    "classic_caf_sas"   
    ## [55] "p_afr"              "p_amr"              "p_asj"             
    ## [58] "p_eas"              "p_fin"              "p_nfe"             
    ## [61] "p_oth"              "p_sas"              "transcript_type"   
    ## [64] "gene_id"            "transcript_level"   "cds_length"        
    ## [67] "num_coding_exons"   "gene_type"          "gene_length"       
    ## [70] "exac_pLI"           "exac_obs_lof"       "exac_exp_lof"      
    ## [73] "exac_oe_lof"        "brain_expression"   "chromosome"        
    ## [76] "start_position"     "end_position"

``` r
gnomad_loeuf <- gnomad %>% dplyr::select(gene, oe_lof_upper)
head(gnomad_loeuf)
```

    ##    gene oe_lof_upper
    ## 1 MED13        0.030
    ## 2 NIPBL        0.032
    ## 3  SMC3        0.037
    ## 4 CNOT1        0.038
    ## 5   RLF        0.040
    ## 6 PCF11        0.040

## Annotate Genes with LOEUF Scores

``` r
#make named vector
loeuf <- deframe(gnomad_loeuf[ ,c("gene", "oe_lof_upper")])
head(loeuf)
```

    ## MED13 NIPBL  SMC3 CNOT1   RLF PCF11 
    ## 0.030 0.032 0.037 0.038 0.040 0.040

``` r
#match to symbols
loeuf_scores <- lapply(symbols, function(x)loeuf[x])
length(loeuf_scores) == length(symbols) #matched
```

    ## [1] TRUE

``` r
loeuf_scores
```

    ## [[1]]
    ##    SCN2A  SYNGAP1     ADNP   GRIN2B   CTNNB1   STXBP1   GABRB2   CORO1A 
    ##    0.127    0.052    0.123    0.061    0.129    0.086    0.385    0.320 
    ##    EPHB1      NF1     THRB   MAP2K7 HSP90AB1    GSK3B      UBB 
    ##    0.257    0.290    0.223    0.232    0.213    0.334    1.560 
    ## 
    ## [[2]]
    ## GRIN2B CTNNB1    NF1 MAP2K7  GSK3B  SRPK2 FCGR2B   ATF4    FOS    GRN  GRIK5 
    ##  0.061  0.129  0.290  0.232  0.334  0.269  0.634  1.091  0.645  0.483  0.285 
    ##   PRNP 
    ##  1.299 
    ## 
    ## [[3]]
    ## SYNGAP1    PTEN  BCL11A    GFAP 
    ##   0.052   0.507   0.322   1.026 
    ## 
    ## [[4]]
    ##     PTEN   SHANK3    SIN3A     TBR1     LDB1    EPHB1 HSP90AB1      UBB 
    ##    0.507    0.123    0.070    0.194    0.291    0.257    0.213    1.560 
    ##    SATB2     SUFU   GABRB1     NFIB    EPHB2    NR4A2     HES5 HSP90AA1 
    ##    0.091    0.111    0.306    0.225    0.184    0.133    1.344    0.353 
    ##    FEZF2     ISL2  RAPGEF2     RAC3     RORA     NDNF ARHGAP35      DCC 
    ##    0.222    0.670    0.197    1.266    0.404    0.512    0.057    0.283 
    ##     MAP2 
    ##    0.105 
    ## 
    ## [[5]]
    ##   ADNP  DSCAM  RUFY3 PLXNA2 PLXNA1 MAP2K1  STK11  LIMK1 PLXNB1 
    ##  0.123  0.189  0.458  0.435  0.262  0.377  0.245  0.225  0.385 
    ## 
    ## [[6]]
    ##    ADNP  SHANK3  CTNNB1   DSCAM    GFAP   RUFY3    RELN  PLXNA2   ITPKA    CUX2 
    ##   0.123   0.123   0.129   0.189   1.026   0.458   0.176   0.435   0.472   0.187 
    ##   EPHB2  PPP1CC   PTPRD  PLXNA1 CAPRIN1  MAP2K1   STK11   SOX10    MTOR   LIMK1 
    ##   0.184   0.336   0.112   0.262   0.311   0.377   0.245   0.209   0.185   0.225 
    ##  PLXNB1  DICER1   ACTR2   VEGFC   XRCC5   HIF1A  PTPRZ1     KIT   TIAM1  SEMA4D 
    ##   0.385   0.166   0.195   0.478   0.147   0.309   0.286   0.305   0.304   0.314 
    ##    NPTN  PLXNB2   NDEL1 
    ##   0.247   0.285   0.296 
    ## 
    ## [[7]]
    ##    PTEN  CTNNB1  STXBP1    PCLO   GSK3B PRKAR1B   VAMP2     GAK   AP2M1   AP3D1 
    ##   0.507   0.129   0.086   0.124   0.334   0.538   0.422   0.408   0.190   0.296 
    ##   BRSK1    SYT1 SLC18A2   GRIK5   RIMS1  SNAP91 RAPGEF4    HPCA  SH3GL1     DDC 
    ##   0.205   0.451   0.511   0.285   0.285   0.193   0.288   0.594   0.639   0.915 
    ##  PPFIA3   ADCY1   ITSN1   VPS35   AP3B1      TH    SYT2 SLC17A6    CANX   GSG1L 
    ##   0.121   0.272   0.221   0.313   0.344   0.788   0.287   0.469   0.359   0.464 
    ##   WNT3A   BRSK2    CNR1    STX3   STON2    SV2A   NLGN2    SYT7  STXBP5   ITGB3 
    ##   0.320   0.329   0.602   0.604   0.680   0.414   0.126   0.329   0.256   0.519 
    ##    DVL1   CADPS   PRKCB   SYNJ1  UNC13A   APBA1    <NA>   ITSN2    ACTB  CYFIP1 
    ##   0.733   0.435   0.081   0.330   0.158   0.388      NA   0.491   0.232   0.302 
    ##    DRD3  ABCA13   BTBD9    DNM3    NUMB   LRRK2 
    ##   0.807   1.035   0.866   0.483   0.593   0.640 
    ## 
    ## [[8]]
    ##     CHD8    KDM6B    KDM5B  SMARCC2    SATB1  SMARCA2   INO80D     CTCF 
    ##    0.082    0.140    0.572    0.117    0.293    0.203    0.228    0.148 
    ##    GATA3 SMARCAD1    BAZ2B    SATB2  SMARCE1    KDM4A 
    ##    0.388    0.083    0.288    0.091    0.191    0.156 
    ## 
    ## [[9]]
    ##    DSCAM    NRXN1   DPYSL2     TBR1    SMAD4    EPHB1     TRIO    GATA3 
    ##    0.189    0.254    0.274    0.194    0.222    0.257    0.139    0.388 
    ##     RELN    CNTN1   PLXNA2    FLRT2     NFIB    UNC5D     ECE1    EPHB2 
    ##    0.176    0.373    0.435    0.413    0.225    0.330    0.321    0.184 
    ##   PLXNA1    SLIT1   VSTM2L    FEZF2    EFNB3     ISL2   SEMA3F      APP 
    ##    0.262    0.222    0.841    0.222    0.607    0.670    0.219    0.419 
    ##   PLXNB1   LYPLA2     <NA>    EFNA2   POU4F3 ARHGAP35      DCC  RPS6KA5 
    ##    0.385    0.344       NA    1.192    0.370    0.057    0.283    0.242 
    ##   SEMA4D      FYN     NPTN   PLXNB2    SLIT3    TUBB3    SIAH1     NTN5 
    ##    0.314    0.278    0.247    0.285    0.318    0.315    0.499    1.739 
    ##     KLF7   BCL11B    WNT3A     OTX2    KIF5B 
    ##    0.245    0.282    0.320    0.376    0.285 
    ## 
    ## [[10]]
    ##  SLC6A1   NRXN1  STXBP1    GFAP     NF1    PCLO   GSK3B    ACHE  CAMK2A   VAMP2 
    ##   0.150   0.254   0.086   1.026   0.290   0.124   0.334   0.207   0.251   0.422 
    ##   ASIC1   BRSK1    SYT1 SLC18A2    PER2   GRIK5   RIMS1  PPFIA3   ADCY1      TH 
    ##   0.344   0.205   0.451   0.511   0.347   0.285   0.285   0.121   0.272   0.788 
    ## SLC6A12  SLC1A2    SYT2   GPER1    CNR1    STX3    SV2A    SYT7    HRH3 SLC6A11 
    ##   0.825   0.424   0.287   0.695   0.602   0.604   0.414   0.329   0.693   0.640 
    ##  SLC6A3  SLC6A2  STXBP5   ITGB3  SLC5A7    DVL1   MEF2C   CADPS   PRKCB   SYNJ1 
    ##   0.251   0.476   0.256   0.519   0.562   0.733   0.608   0.435   0.081   0.330 
    ##  UNC13A  PTPRN2   MCTP1 
    ##   0.158   0.568   0.493 
    ## 
    ## [[11]]
    ## SYNGAP1    ADNP  SLC6A1  GRIN2B    PTEN  SHANK3  CTNNB1   DSCAM   SETD5  SHANK2 
    ##   0.052   0.123   0.150   0.061   0.507   0.123   0.129   0.189   0.227   0.458 
    ##   NRXN1  GABRB3  LRRC4C  GABRB2   DIP2A   TANC2   EPHB1 KIRREL3   LRRC4    PCLO 
    ##   0.254   0.341   0.306   0.385   0.354   0.158   0.257   0.318   0.166   0.124 
    ##   LZTS3    ACHE  LRRTM2      C3    RELN    NFIA    ANK3  GABRA1    ABI2 
    ##   0.228   0.207   0.416   0.305   0.176   0.202   0.089   0.367   0.460 
    ## 
    ## [[12]]
    ## SYNGAP1    ADNP    PTEN  SHANK3  CTNNB1   DSCAM     SKI    GFAP     NF1 
    ##   0.052   0.123   0.507   0.123   0.129   0.189   0.194   1.026   0.290 
    ## 
    ## [[13]]
    ##     TLK2   DNMT3A    SIN3A  SMARCC2     PHF2     CTCF  SMARCE1     UBR2 
    ##    0.114    1.581    0.070    0.117    0.177    0.148    0.191    0.213 
    ## MPHOSPH8   CHAF1A    CDAN1    HDAC5  L3MBTL3     RSF1     <NA>    KMT2D 
    ##    0.525    0.118    0.675    0.141    0.269    0.044       NA    0.103 
    ##    SSRP1     CTR9    KAT6A    KAT6B     EZH2     TLK1  SMARCC1     MTA2 
    ##    0.275    0.135    0.069    0.128    0.146    0.165    0.116    0.364 
    ##  SUV39H2    HDAC1    BEND3     <NA>     MBD2  SMARCA5    ARID2     <NA> 
    ##    0.358    0.405    0.347       NA    0.522    0.074    0.096       NA 
    ##   TRIM28     <NA>   CABIN1 
    ##    0.082       NA    0.586 
    ## 
    ## [[14]]
    ##  DSCAM  SLIT1 SEMA3F SEMA4D  SLIT3  WNT3A SEMA3E 
    ##  0.189  0.222  0.219  0.314  0.318  0.320  0.679 
    ## 
    ## [[15]]
    ##  NRXN1 STXBP1    NF1   PCLO  GSK3B CAMK2A  VAMP2  ASIC1  BRSK1   SYT1  GRIK5 
    ##  0.254  0.086  0.290  0.124  0.334  0.251  0.422  0.344  0.205  0.451  0.285 
    ##  RIMS1 PPFIA3  ADCY1   SYT2  GPER1   CNR1   STX3   SV2A   SYT7   HRH3 STXBP5 
    ##  0.285  0.121  0.272  0.287  0.695  0.602  0.604  0.414  0.329  0.693  0.256 
    ##   DVL1  MEF2C  CADPS  PRKCB  SYNJ1 UNC13A PTPRN2  MCTP1  APBA1 
    ##  0.733  0.608  0.435  0.081  0.330  0.158  0.568  0.493  0.388 
    ## 
    ## [[16]]
    ## EPHB1  NFIB EPHB2 NR4A2   DCC DCLK1  TSKU  SZT2 
    ## 0.257 0.225 0.184 0.133 0.283 0.227 0.474 0.405 
    ## 
    ## [[17]]
    ##  SLIT1 SEMA3F SEMA4D  WNT3A SEMA3E SEMA5B SEMA5A SEMA3D SEMA3G SEMA6C SEMA4F 
    ##  0.222  0.219  0.314  0.320  0.679  0.604  0.483  0.657  1.037  0.635  0.919 
    ## SEMA4B SEMA4G 
    ##  0.533  0.756 
    ## 
    ## [[18]]
    ## STXBP1   PCLO  GSK3B  VAMP2   SYT1  GRIK5  RIMS1 PPFIA3  ADCY1   SYT2   CNR1 
    ##  0.086  0.124  0.334  0.422  0.451  0.285  0.285  0.121  0.272  0.287  0.602 
    ##   STX3   SV2A   SYT7 STXBP5   DVL1  CADPS  PRKCB  SYNJ1 UNC13A  APBA1  LRRK2 
    ##  0.604  0.414  0.329  0.256  0.733  0.435  0.081  0.330  0.158  0.388  0.640 
    ## UNC13C   OTOF UNC13B 
    ##  0.537  0.875  0.586 
    ## 
    ## [[19]]
    ##    AKT2  RICTOR   STK11    MTOR   MAPK3 RPS6KA1   VEGFC   HIF1A  PIK3R1   MLST8 
    ##   0.417   0.155   0.245   0.185   0.614   0.515   0.478   0.309   0.228   0.671 
    ##    BRAF  PIK3R2 RPS6KA2  PRKAA1 
    ##   0.209   0.489   0.445   0.478 
    ## 
    ## [[20]]
    ##   HDAC4    BPTF   CHD1L    <NA>   PBRM1    <NA>  JARID2    <NA>    <NA>    <NA> 
    ##   0.052   0.099   1.299      NA   0.129      NA   0.188      NA      NA      NA 
    ##    <NA>  INO80E   HMGA1    <NA>    MCM2    <NA>  YEATS2    RERE  RUVBL2    <NA> 
    ##      NA   0.950   0.509      NA   0.585      NA   0.280   0.120   0.111      NA 
    ##  ANP32E   KAT6A    <NA> SMARCC1  INO80C    <NA>   ARID2   KDM5B  RUVBL1    CHD1 
    ##   0.423   0.069      NA   0.116   1.416      NA   0.096   0.572   0.194   0.162 
    ##   BAZ2B    <NA>    <NA>   FOXA1     MYC    <NA>   KDM4A    <NA>   NFRKB   SYCP3 
    ##   0.288      NA      NA   0.677   0.164      NA   0.156      NA   0.368   0.949 
    ##    CHD7    <NA>   CENPV    DAXX    <NA>   KDM4D  TSPYL6  YEATS4 GATAD2A   HDAC2 
    ##   0.076      NA   0.736   0.496      NA   1.614      NA   1.189   0.202   0.100 
    ##    <NA>    <NA>  CHRAC1    RSF1  MIS18A    <NA>   BAZ1B  NAP1L4    <NA>    <NA> 
    ##      NA      NA   1.492   0.044   0.995      NA   0.109   0.392      NA      NA 
    ##    <NA> SMARCA2    TFPT 
    ##      NA   0.203   1.123 
    ## 
    ## [[21]]
    ##   WNT3 PLXNA4 SEMA4D    RYK SEMA3G  ALCAM SEMA7A SEMA6C  VEGFA  SLIT2 SEMA5A 
    ##  0.404  0.229  0.314  0.289  1.037  0.386  0.441  0.635  0.840  0.160  0.483 
    ## SEMA3F SEMA5B  MEGF8 SEMA4B SEMA4G   NRP2   <NA> SEMA3D SEMA6D SEMA4F 
    ##  0.219  0.604  0.366  0.533  0.756  0.475     NA  0.657  0.323  0.919 
    ## 
    ## [[22]]
    ##     MAPT   CHRNB2   CACNB2     GPC6     SDK1     FZD5    EPHA7   CACNB3 
    ##    0.569    0.822    0.692    0.511    0.398    0.283    0.128    0.644 
    ##     GRM5   GABRB2    NEGR1     DAG1      BSN     PIN1   SEMA4D    EPHB1 
    ##    0.254    0.385    0.253    0.301    0.175    0.540    0.314    0.257 
    ##    GNPAT     INSR     RHOA     REST   CTNND2     NFIA    DCTN1    WASF2 
    ##    0.417    0.453    0.664    0.286    0.098    0.202    0.364    0.311 
    ##   UNC13A    LRRC4     CDH2   CTNNA2   PPFIA3     LRP8   SHISA7   NEURL1 
    ##    0.158    0.166    0.290    0.378    0.121    0.224    0.474    0.739 
    ##   LRRC24   PPFIA1    FARP1    PTPRF     LGI2   AFG3L2     WASL      ARC 
    ##    0.961    0.067    0.509    0.187    0.894    0.700    0.375    0.491 
    ##    CAMKV   ZNF365    LRFN4      RYK     <NA> PAFAH1B1   SHANK1     SNCA 
    ##    0.221    0.757    0.387    0.289       NA    0.107    0.184    0.429 
    ##   NFATC4     <NA>     GDNF  SLITRK3    EPHB2    SETD5     IL10    OBSL1 
    ##    0.358       NA    0.918    0.340    0.184    0.227    1.039    0.878 
    ##    MDGA1    IGSF9   SHANK2     <NA>   PCDH17     NRG1   CX3CR1    CTBP2 
    ##    0.311    0.705    0.458       NA    0.401    0.258    1.010    0.528 
    ##   DTNBP1   LINGO2     CFL1    PICK1    CNTN2    ACTG1  CHCHD10   CACNB4 
    ##    0.984    0.416    0.716    0.671    0.437    0.858    1.910    0.491 
    ##    NTRK2   FBXO45     LGMN  SYNDIG1    CDC20     ABI2   FCGR2B   CTNNB1 
    ##    0.112    0.274    0.685    0.768    0.838    0.460    0.634    0.129 
    ##     PTEN   GRIN2B  SPARCL1    PSEN1    CRIPT     SNCG   SEMA3F   CHRNB1 
    ##    0.507    0.061    0.967    0.318    1.582    1.625    0.219    0.996 
    ##  CNTNAP1  ABHD17C      TNR    EFNA5      VCP     ADD2    C1QL1    YWHAZ 
    ##    0.422    0.302    0.336    0.296    0.127    0.274    1.907    0.357 
    ##     PFN1    TANC2    CBLN3   COL4A1  CTTNBP2     DLG1     NRP2    MTMR2 
    ##    0.686    0.158    1.287    0.128    0.602    0.286    0.475    0.654 
    ##   GABRA2   PCDHB6     DLG5     DLG4     AGRN     CDH8     <NA>     PRNP 
    ##    0.229    1.371    0.262    0.238    0.435    0.295       NA    1.299 
    ##   LILRB2   SNAPIN     KLK8   SHANK3    NRCAM    SNTA1     <NA>    CBLN2 
    ##    1.119    1.437    0.955    0.123    0.383    0.707       NA    1.165 
    ##     TPBG     NEFH    C1QL3     CRKL    NRXN2     BDNF   EIF4G1  NEUROD2 
    ##    1.160    1.064    0.763    0.641    0.258    0.520    0.074    0.332 
    ##   CACNG2 
    ##    0.376 
    ## 
    ## [[23]]
    ##     MAPT     WNT3   PLXNA4     PAX6     HAP1   TMEM98    EPHA7     GRM5 
    ##    0.569    0.404    0.229    0.169    1.637    0.886    0.128    0.254 
    ##     DAG1   TRIM32   SEMA4D     MTOR     REST     PRTG    ROBO1      ID4 
    ##    0.301    0.855    0.314    0.185    0.286    0.562    0.597    1.833 
    ##     RHEB    SHOX2    PLAG1     LRP8    LPAR3   NEURL1     IST1   TRIM46 
    ##    0.247    0.831    0.188    0.224    0.593    0.739    0.727    0.228 
    ##     ELL3   PPP3CA     BMP7     UFL1    SYNJ1     FAIM      MYC   ZNF365 
    ##    1.163    0.158    0.570    0.843    0.330    1.382    0.164    0.757 
    ##   SS18L1      RYK     <NA>      NF1 PAFAH1B1   SEMA3G    KDM4A   PLXNA2 
    ##    0.306    0.289       NA    0.290    0.107    1.037    0.156    0.435 
    ##   NFATC4     DLL4     <NA>     CHD7    RGS14    EPHB2    SORL1    OBSL1 
    ##    0.358    0.188       NA    0.076    0.576    0.184    0.436    0.878 
    ##     E2F1   SEMA7A    PITX3   SEMA6C    HDAC2     HES5   CX3CR1    VEGFA 
    ##    0.225    0.441    0.479    0.635    0.100    1.344    1.010    0.840 
    ##    SOX10    SLIT2      NIN    NTRK2   SEMA5A   FBXO31     ASPM     SPEN 
    ##    0.209    0.160    0.429    0.112    0.483    0.457    0.743    0.071 
    ##   RNF112   CTNNB1     PTEN   PPP1CC    RAB21    DOCK7    PSEN1     HELT 
    ##    0.541    0.129    0.507    0.336    0.671    0.379    0.318    1.177 
    ##      KIT   SEMA3F      ACE      TNR     HES3    EFNA5    SIRT2    FOXG1 
    ##    0.305    0.219    1.081    0.336    1.246    0.296    0.957    0.329 
    ##    ABCC8   SEMA5B    MEGF8   SEMA4B     SNW1     TP53      SKI   SEMA4G 
    ##    0.765    0.604    0.366    0.533    0.217    0.469    0.194    0.756 
    ##     <NA>   SEMA3D     LDLR    WDR62     TWF2   GOLGA4   SHANK3     <NA> 
    ##       NA    0.657    1.128    0.859    0.302    0.515    0.123       NA 
    ##   SEMA6D      ID2     IFNG     THY1   SEMA4F      LYN   TRIM11     BDNF 
    ##    0.323    0.827    0.822    1.084    0.919    0.221    0.375    0.520 
    ##     EGR2    ISLR2     MAP6   MAP2K1      PTN     DLX2 
    ##    0.604    0.583    0.568    0.377    0.982    0.254 
    ## 
    ## [[24]]
    ##   CACNB2     SV2A    PRRT2    STX19   SNCAIP   STXBP3    GSK3B   UNC13A 
    ##    0.692    0.414    0.560    1.677    0.552    0.593    0.334    0.158 
    ##   PPFIA3    DOC2A    ASIC1    CPLX1   CHRNB4    SYNJ1    VPS18    CSPG5 
    ##    0.121    0.690    0.344    0.521    0.559    0.330    0.395    0.654 
    ##     SYT2      NF1   GRIN3A     SNCA    STX1A   PTPRN2   STXBP2   DNAJC5 
    ##    0.287    0.290    0.386    0.429    0.293    0.568    0.842    0.480 
    ##   SNAP23   SNAP47    CTBP2   CAMK2A   DTNBP1    CADPS    SYT10   CHRNA5 
    ##    0.375    0.994    0.528    0.251    0.984    0.435    0.753    1.136 
    ##    PRKCB    STX11    CPLX3     <NA>    PSEN1     SNCG   STXBP5    RPH3A 
    ##    0.081    1.379    0.851       NA    0.318    1.625    0.256    0.474 
    ##    RAB5A RAB3GAP1   CHRNA3    P2RY1   KCNMB4   SNAPIN     NAPB     CNR1 
    ##    0.536    0.535    1.148    0.597    0.338    1.437    0.582    0.602 
    ##     HRH3    ADCY1 
    ##    0.693    0.272 
    ## 
    ## [[25]]
    ##     MAPT   PLXNA4     GBX1   CHRNB2     PAX6    EPHB1   SPOCK1      ID4 
    ##    0.569    0.229    0.668    0.822    0.169    0.257    0.392    1.833 
    ##     LHX5     LHX3     LHX1    DCLK2  BLOC1S2   ZSWIM6     ISL1     RORA 
    ##    0.348    0.742    0.737    0.326    1.178    0.067    0.408    0.404 
    ##    TULP3 PAFAH1B1     DLL4     LHX6    EPHB2    CEND1    ZMIZ1    MDGA1 
    ##    0.879    0.107    0.188    0.326    0.184    1.029    0.252    0.311 
    ##   TTC21B     HES5   HOXD10    SLIT2      NIN    CNTN2     LDB1    NTRK2 
    ##    0.852    1.344    0.762    0.160    0.429    0.437    0.291    0.112 
    ##   FBXO45     PTEN   GIGYF2    GATA2    CDH11     NFIB    PSEN1   PHOX2A 
    ##    0.274    0.507    0.077    0.292    0.156    0.225    0.318    0.796 
    ##     TAL1      TOX    PROX1 ARHGAP35     TBR1    FOXG1     ISL2    FGFR2 
    ##    0.704    0.200    0.221    0.057    0.194    0.329    0.670    0.270 
    ##  B4GALT6     NRP2     CLN8   SHANK3    FAIM2     EMX1     SOX1    TTLL1 
    ##    0.371    0.475    1.174    0.123    0.820    0.409    0.951    0.870 
    ## 
    ## [[26]]
    ##    <NA>    <NA>    <NA>    <NA>    <NA>   AXIN1    <NA>   HMGA1    <NA>    MCM2 
    ##      NA      NA      NA      NA      NA   0.403      NA   0.509      NA   0.585 
    ##    <NA>    <NA>   DNMT1   KAT6A    <NA>   SIRT1 SMARCC1   CTBP1    <NA>   ARID2 
    ##      NA      NA   0.143   0.069      NA   0.503   0.116   0.278      NA   0.096 
    ##    <NA>    <NA>    <NA>  HIRIP3    <NA>    <NA>    UBN2   ZDBF2   CENPV    DAXX 
    ##      NA      NA      NA   0.867      NA      NA   0.213   0.460   0.736   0.496 
    ##    <NA>    M1AP  TSPYL6    <NA>    <NA>    <NA>  CHRAC1    RSF1  MIS18A    <NA> 
    ##      NA   1.405      NA      NA      NA      NA   1.492   0.044   0.995      NA 
    ##  SMCHD1    <NA>  NAP1L4    <NA>    <NA>    <NA>  ZNF445    TAL1    RIF1  PPHLN1 
    ##   0.153      NA   0.392      NA      NA      NA   0.212   0.704   0.113   0.836 
    ##    <NA>   SIRT2    <NA>    NPM1    TP53   BAHD1     SET    HIRA    <NA>   ASF1B 
    ##      NA   0.957      NA   0.162   0.469   0.300   0.183   0.136      NA   1.531 
    ##    <NA>  PIK3CA   POLE3    <NA>    <NA>    MBD3   BAZ1A    <NA>   MAP1S    <NA> 
    ##      NA   0.117   1.284      NA      NA   0.260   0.176      NA   0.436      NA 
    ## SMARCD3   SIN3A   PADI2  NAP1L1    PAF1   HDAC1    NASP    <NA> 
    ##   0.466   0.070   0.891   0.184   0.622   0.405   0.276      NA 
    ## 
    ## [[27]]
    ##    MTOR    RHEB    ULK3    BRAF  EIF4E2   VEGFA RPS6KA1   MAPK3    AKT2    RPS6 
    ##   0.185   0.247   0.584   0.209   0.495   0.840   0.515   0.614   0.417   0.261 
    ##  PRKAA2  PIK3R2  PIK3R3  PIK3CB RPS6KB1   MLST8  PIK3CA   VEGFB  CAB39L   RPTOR 
    ##   0.921   0.489   1.004   0.252   0.317   0.671   0.117   1.873   1.003   0.101 
    ##   DDIT4 RPS6KB2 
    ##   1.278   1.542 
    ## 
    ## [[28]]
    ##   CACNB2     SV2A    PRRT2    STX19     DNM1      BSN     DNM2  RAPGEF4 
    ##    0.692    0.414    0.560    1.677    0.252    0.175    0.180    0.288 
    ##    GSK3B   DNAJC6   UNC13A     CDH2   PPFIA3    DOC2A    CPLX1    SYNJ2 
    ##    0.334    0.222    0.158    0.290    0.121    0.690    0.521    0.691 
    ##    SYNJ1    VPS18    ITSN1      ARC    CSPG5     SYT2   GRIN3A     SNCA 
    ##    0.330    0.395    0.221    0.491    0.654    0.287    0.386    0.429 
    ##    STX1A  SLC17A8    ROCK1   DNAJC5   SNAP23   SNAP47   PCDH17    CTBP2 
    ##    0.293    0.608    0.196    0.480    0.375    0.994    0.401    0.528 
    ##   DTNBP1    CADPS    SYT10   CHRNA5    ACTG1    PRKCB  SYNDIG1     <NA> 
    ##    0.984    0.435    0.753    1.136    0.858    0.081    0.768       NA 
    ##    STX11   CTNNB1     PTEN  PACSIN1    CPLX3    PSEN1     SNCG   SLC2A4 
    ##    1.379    0.129    0.507    0.197    0.851    0.318    1.625    0.777 
    ##   STXBP5     CANX      DDC    AP3B1    RAB5A RAB3GAP1    AP2M1    P2RY1 
    ##    0.256    0.359    0.915    0.344    0.536    0.535    0.190    0.597 
    ##   SNAPIN     NAPB     CNR1    ADCY1 
    ##    1.437    0.582    0.602    0.272 
    ## 
    ## [[29]]
    ##    MAPT   SRPK2  TFAP2B   EPHA7 MAP3K11    RHOA   FBXW7    REST   CASP8   GSK3B 
    ##   0.569   0.269   0.259   0.128   0.479   0.664   0.232   0.286   0.811   0.334 
    ##   ITGA1     NF1    SNCA   PARP1   TGFB2    TLR4    DAXX   PITX3     BAX   CASP2 
    ##   0.739   0.290   0.429   0.479   0.148   0.992   0.496   0.479   0.751   0.658 
    ##  FCGR2B  CTNNB1  GRIN2B 
    ##   0.634   0.129   0.061 
    ## 
    ## [[30]]
    ##     WNT3   PLXNA4     GBX1     PAX6     RAC1     LGI1    EPHA7     DAG1 
    ##    0.404    0.229    0.668    0.169    0.467    0.191    0.128    0.301 
    ##   SEMA4D    EPHB1     ENAH    ROBO1     LHX9    EPHA5     LHX3    EPHA2 
    ##    0.314    0.257    0.173    0.597    0.288    0.521    0.742    0.555 
    ##     LHX1     BMP7     ISL1    CRMP1      RYK     RAC2   SEMA3G    ALCAM 
    ##    0.737    0.570    0.408    0.281    0.289    0.284    1.037    0.386 
    ##   PLXNA2     <NA>     GDNF    EPHB2    USP33   SEMA7A    IGSF9     EXT1 
    ##    0.435       NA    0.918    0.184    0.331    0.441    0.705    0.261 
    ##    EPHA1   SEMA6C    LAMA1    KIF5A    VEGFA    SLIT2    CNTN2   SEMA5A 
    ##    1.045    0.635    0.621    0.228    0.840    0.160    0.437    0.483 
    ##    SIAH1     EDN3     NFIB   SEMA3F ARHGAP35     TBR1    NRXN3      TNR 
    ##    0.499    0.964    0.225    0.219    0.057    0.194    0.112    0.336 
    ##    EFNA5     NEXN    FOXG1   SEMA5B     ISL2    MEGF8     FEZ1   SEMA4B 
    ##    0.296    0.779    0.329    0.604    0.670    0.366    0.414    0.533 
    ##   SEMA4G     NRP2     <NA>   SEMA3D   OR10A4    LRTM1    GATA3     RHOH 
    ##    0.756    0.475       NA    0.657    1.797    1.903    0.388    1.061 
    ##    NRCAM     <NA>   SEMA6D   SEMA4F 
    ##    0.383       NA    0.323    0.919 
    ## 
    ## [[31]]
    ##     MAPT    SRPK2   TFAP2B    CASP7    EPHA7     SOD1   TFAP2D    HLA-F 
    ##    0.569    0.269    0.259    0.873    0.128    0.978    0.252    0.933 
    ##  MAP3K11   GABRB2  KIR3DL2    EPHB1    HYOU1    AGAP2     RHOA    KCNB1 
    ##    0.479    0.385    1.870    0.257    0.331    0.267    0.664    0.125 
    ##    KDM2B    FBXW7     REST    CASP8    GSK3B  ARL6IP5    ITGA1    SIRT1 
    ##    0.178    0.232    0.286    0.811    0.334    1.254    0.739    0.503 
    ##    FOXB1   DHCR24   RILPL1    TRIM2 TNFRSF21     ISL1     <NA>      NF1 
    ##    0.392    0.523    0.470    0.270    0.606    0.408       NA    0.290 
    ##   MAP2K4  TBC1D24     SNCA   NFATC4     <NA>    PARP1    TGFB2   INPP5A 
    ##    0.220    1.123    0.429    0.358       NA    0.479    0.148    0.308 
    ##     GDNF     BAG5     CD34     CHGA     GPX1    BNIP3     TLR4    SORL1 
    ##    0.918    1.159    1.024    0.568    1.892    1.462    0.992    0.436 
    ##     DAXX     IL10     GBE1    ROCK1    PITX3   DNAJC5      BAX    SCN2A 
    ##    0.496    1.039    0.972    0.196    0.479    0.480    0.751    0.127 
    ##   CX3CR1    EGLN2  TMEM259   DTNBP1      EN2    CASP2    NTRK2     LGMN 
    ##    1.010    0.447    0.241    0.984    0.634    0.658    0.112    0.685 
    ##  SLC23A2     CTSZ    CRLF1   FCGR2B   CTNNB1    SIAH1   SLC9A1  ATP13A2 
    ##    0.528    1.919    0.908    0.634    0.129    0.499    0.375    0.584 
    ##    CNTFR     <NA>   GRIN2B   MTNR1B    PSEN1     NAE1     AKT2     SNCG 
    ##    0.339       NA    0.061    1.293    0.318    0.941    0.417    1.625 
    ## 
    ## [[32]]
    ##     MAPT     WNT3   PLXNA4     PAX6     HAP1     GRM5     DAG1   TRIM32 
    ##    0.569    0.404    0.229    0.169    1.637    0.254    0.301    0.855 
    ##   SEMA4D     MTOR    ROBO1     RHEB    SHOX2    PLAG1     LRP8    LPAR3 
    ##    0.314    0.185    0.597    0.247    0.831    0.188    0.224    0.593 
    ##   NEURL1     IST1     ELL3     UFL1    SYNJ1     FAIM      MYC   ZNF365 
    ##    0.739    0.727    1.163    0.843    0.330    1.382    0.164    0.757 
    ##   SS18L1     <NA> PAFAH1B1   PLXNA2     <NA>    RGS14    EPHB2    OBSL1 
    ##    0.306       NA    0.107    0.435       NA    0.576    0.184    0.878 
    ##     E2F1   SEMA7A    HDAC2   CX3CR1    VEGFA    SOX10    SLIT2      NIN 
    ##    0.225    0.441    0.100    1.010    0.840    0.209    0.160    0.429 
    ##    NTRK2   SEMA5A   FBXO31     ASPM     SPEN   RNF112   CTNNB1   PPP1CC 
    ##    0.112    0.483    0.457    0.743    0.071    0.541    0.129    0.336 
    ##    RAB21 
    ##    0.671 
    ## 
    ## [[33]]
    ##   WNT3 SEMA4D    RYK SEMA3G SEMA7A SEMA6C SEMA5A SEMA3F SEMA5B SEMA4B SEMA4G 
    ##  0.404  0.314  0.289  1.037  0.441  0.635  0.483  0.219  0.604  0.533  0.756 
    ##   <NA> SEMA3D SEMA6D SEMA4F 
    ##     NA  0.657  0.323  0.919 
    ## 
    ## [[34]]
    ##   PLXNA4   CHRNB2    EPHB1 PAFAH1B1    EPHB2    SLIT2      NIN   FBXO45 
    ##    0.229    0.822    0.257    0.107    0.184    0.160    0.429    0.274 
    ##    CDH11     NFIB 
    ##    0.156    0.225 
    ## 
    ## [[35]]
    ##     MAPT     WNT3   PLXNA4   SEMA4D    ROBO1    SHOX2    LPAR3     IST1 
    ##    0.569    0.404    0.229    0.314    0.597    0.831    0.593    0.727 
    ## PAFAH1B1   PLXNA2     <NA>   SEMA7A    VEGFA    SLIT2      NIN    NTRK2 
    ##    0.107    0.435       NA    0.441    0.840    0.160    0.429    0.112 
    ##   SEMA5A 
    ##    0.483 
    ## 
    ## [[36]]
    ##   CACNB2     SV2A    DAGLB    PRRT2    STX19   SNCAIP   STXBP3    GSK3B 
    ##    0.692    0.414    0.995    0.560    1.677    0.552    0.593    0.334 
    ##   UNC13A   PPFIA3  SLC22A2    DOC2A    ASIC1    CPLX1   CHRNB4    SYNJ1 
    ##    0.158    0.121    1.154    0.690    0.344    0.521    0.559    0.330 
    ##  SLC6A11    VPS18     NOS1    CSPG5     SYT2    PDE1B   PRIMA1      NF1 
    ##    0.640    0.395    0.190    0.654    0.287    0.430    0.640    0.290 
    ##   GRIN3A     SNCA    STX1A     GDNF   PTPRN2   STXBP2   DNAJC5   SNAP23 
    ##    0.386    0.429    0.293    0.918    0.568    0.842    0.480    0.375 
    ##   SNAP47    CTBP2   CAMK2A   DTNBP1    CADPS    SYT10   CHRNA5     COMT 
    ##    0.994    0.528    0.251    0.984    0.435    0.753    1.136    1.703 
    ##    PRKCB    STX11    CPLX3     <NA>    PSEN1     SNCG   SLC5A7   STXBP5 
    ##    0.081    1.379    0.851       NA    0.318    1.625    0.562    0.256 
    ##    RPH3A      PAH    RAB5A RAB3GAP1   CHRNA3   GABRA2     CLN8    P2RY1 
    ##    0.474    1.502    0.536    0.535    1.148    0.229    1.174    0.597 
    ##   KCNMB4  SLC22A3   SNAPIN     NAPB     CNR1     HRH3    ADCY1 
    ##    0.338    1.034    1.437    0.582    0.602    0.693    0.272 
    ## 
    ## [[37]]
    ##     WNT3     <NA>    EPHA7   SEMA4D     RHOA   SPOCK1   TRIM46    CRMP1 
    ##    0.404       NA    0.128    0.314    0.664    0.392    0.228    0.281 
    ##   ZNF365      RYK     GFI1     <NA>   BCL11A PAFAH1B1   SEMA3G     TBX6 
    ##    0.757    0.289    0.561       NA    0.322    0.107    1.037    0.688 
    ##     BAG5    FKBP4    EPHB2   SEMA7A   SEMA6C    HDAC2   DTNBP1   SEMA5A 
    ##    1.159    0.504    0.184    0.441    0.635    0.100    0.984    0.483 
    ##     PTEN    MYLIP     <NA>    PSEN1   SEMA3F    TRPV4      TNR   SEMA5B 
    ##    0.507    0.536       NA    0.318    0.219    1.057    0.336    0.604 
    ##   SEMA4B    ADCY6     DAB2     NEU4   SEMA4G     <NA>   SEMA3D  RTN4RL2 
    ##    0.533    0.442    0.419    1.606    0.756       NA    0.657    0.423 
    ##   SNAPIN     KLK8   SEMA6D     THY1   SEMA4F 
    ##    1.437    0.955    0.323    1.084    0.919 
    ## 
    ## [[38]]
    ##   CACNB2     SV2A    PRRT2    STX19    GSK3B   UNC13A   PPFIA3    DOC2A 
    ##    0.692    0.414    0.560    1.677    0.334    0.158    0.121    0.690 
    ##    CPLX1    SYNJ1    VPS18    CSPG5     SYT2   GRIN3A     SNCA    STX1A 
    ##    0.521    0.330    0.395    0.654    0.287    0.386    0.429    0.293 
    ##   DNAJC5   SNAP23   SNAP47    CTBP2   DTNBP1    CADPS    SYT10   CHRNA5 
    ##    0.480    0.375    0.994    0.528    0.984    0.435    0.753    1.136 
    ##    PRKCB    STX11    CPLX3    PSEN1   STXBP5    RAB5A RAB3GAP1    P2RY1 
    ##    0.081    1.379    0.851    0.318    0.256    0.536    0.535    0.597 
    ##   SNAPIN     NAPB     CNR1    ADCY1 
    ##    1.437    0.582    0.602    0.272

``` r
#re-add to table
loeuf_scores <- lapply(loeuf_scores, unname)
classify_genes <- classify_genes %>% mutate(LOEUF = loeuf_scores)
head(classify_genes)
```

    ##   variant                                                   pathway       padj
    ## 1    rare                                         GOBP_NEURON_DEATH 0.08450704
    ## 2    rare                  GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH 0.08450704
    ## 3    rare GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT 0.08450704
    ## 4    rare        GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION 0.08450704
    ## 5    rare                  GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS 0.08450704
    ## 6    rare                  GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS 0.09208973
    ##    leadingEdge converge_diverge        LOEUF
    ## 1 SCN2A, S....        divergent 0.127, 0....
    ## 2 GRIN2B, ....       convergent 0.061, 0....
    ## 3 SYNGAP1,....        divergent 0.052, 0....
    ## 4 PTEN, SH....       convergent 0.507, 0....
    ## 5 ADNP, DS....       convergent 0.123, 0....
    ## 6 ADNP, SH....       convergent 0.123, 0....

## Assess Constraint

``` r
#determine constraint by calculating proportion of genes with LOEUF scores below 0.35 threshold
#calculate proportion ignoring NA values
classify_genes <- classify_genes %>% mutate(percent_constrained = sapply(LOEUF, function(x)((sum(x < 0.35, na.rm = TRUE)/sum(!is.na(x)))*100))) %>% # remove NAs
relocate(percent_constrained, .before = LOEUF)
classify_genes
```

    ##    variant
    ## 1     rare
    ## 2     rare
    ## 3     rare
    ## 4     rare
    ## 5     rare
    ## 6     rare
    ## 7     rare
    ## 8     rare
    ## 9     rare
    ## 10    rare
    ## 11    rare
    ## 12    rare
    ## 13    rare
    ## 14    rare
    ## 15    rare
    ## 16    rare
    ## 17    rare
    ## 18    rare
    ## 19    rare
    ## 20  common
    ## 21  common
    ## 22  common
    ## 23  common
    ## 24  common
    ## 25  common
    ## 26  common
    ## 27  common
    ## 28  common
    ## 29  common
    ## 30  common
    ## 31  common
    ## 32  common
    ## 33  common
    ## 34  common
    ## 35  common
    ## 36  common
    ## 37  common
    ## 38  common
    ##                                                                    pathway
    ## 1                                                        GOBP_NEURON_DEATH
    ## 2                                 GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH
    ## 3                GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT
    ## 4                       GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION
    ## 5                                 GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS
    ## 6                                 GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS
    ## 7                               GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE
    ## 8                                                GOBP_CHROMATIN_REMODELING
    ## 9                                          GOBP_NEURON_PROJECTION_GUIDANCE
    ## 10                              GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS
    ## 11                                               GOBP_SYNAPSE_ORGANIZATION
    ## 12                                         GOBP_REGULATION_OF_NEUROGENESIS
    ## 13                                  GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY
    ## 14 GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE
    ## 15                                         GOBP_NEUROTRANSMITTER_SECRETION
    ## 16              GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS
    ## 17    GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE
    ## 18                                        GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS
    ## 19                                             KEGG_MTOR_SIGNALING_PATHWAY
    ## 20                                               GOBP_CHROMATIN_REMODELING
    ## 21 GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE
    ## 22                                               GOBP_SYNAPSE_ORGANIZATION
    ## 23                                         GOBP_REGULATION_OF_NEUROGENESIS
    ## 24                                         GOBP_NEUROTRANSMITTER_SECRETION
    ## 25                      GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION
    ## 26                                  GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY
    ## 27                                             KEGG_MTOR_SIGNALING_PATHWAY
    ## 28                              GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE
    ## 29                                GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH
    ## 30                                         GOBP_NEURON_PROJECTION_GUIDANCE
    ## 31                                                       GOBP_NEURON_DEATH
    ## 32                                GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS
    ## 33    GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE
    ## 34              GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS
    ## 35                                GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS
    ## 36                              GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS
    ## 37               GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT
    ## 38                                        GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS
    ##          padj  leadingEdge converge_diverge percent_constrained        LOEUF
    ## 1  0.08450704 SCN2A, S....        divergent            86.66667 0.127, 0....
    ## 2  0.08450704 GRIN2B, ....       convergent            58.33333 0.061, 0....
    ## 3  0.08450704 SYNGAP1,....        divergent            50.00000 0.052, 0....
    ## 4  0.08450704 PTEN, SH....       convergent            68.00000 0.507, 0....
    ## 5  0.08450704 ADNP, DS....       convergent            55.55556 0.123, 0....
    ## 6  0.09208973 ADNP, SH....       convergent            78.78788 0.123, 0....
    ## 7  0.09557110 PTEN, CT....       convergent            49.09091 0.507, 0....
    ## 8  0.09671180 CHD8, KD....        divergent            85.71429 0.082, 0....
    ## 9  0.09671180 DSCAM, N....       convergent            68.18182 0.189, 0....
    ## 10 0.09671180 SLC6A1, ....       convergent            51.16279 0.15, 0.....
    ## 11 0.09760766 SYNGAP1,....        divergent            75.86207 0.052, 0....
    ## 12 0.10363636 SYNGAP1,....        divergent            77.77778 0.052, 0....
    ## 13 0.10606061 TLK2, DN....       convergent            74.19355 0.114, 1....
    ## 14 0.10606061 DSCAM, S....       convergent            85.71429 0.189, 0....
    ## 15 0.12800000 NRXN1, S....       convergent            58.06452 0.254, 0....
    ## 16 0.12884753 EPHB1, N....       convergent            75.00000 0.257, 0....
    ## 17 0.18808777 SLIT1, S....       convergent            30.76923 0.222, 0....
    ## 18 0.19524100 STXBP1, ....       convergent            52.00000 0.086, 0....
    ## 19 0.19524100 AKT2, RI....        divergent            42.85714 0.417, 0....
    ## 20 0.02119885 HDAC4, B....        divergent            53.84615 0.052, 0....
    ## 21 0.06722425 WNT3, PL....       convergent            30.00000 0.404, 0....
    ## 22 0.10111381 MAPT, CH....        divergent            40.32258 0.569, 0....
    ## 23 0.12821204 MAPT, WN....        divergent            41.50943 0.569, 0....
    ## 24 0.12821204 CACNB2, ....       convergent            28.57143 0.692, 0....
    ## 25 0.12821204 MAPT, PL....       convergent            50.00000 0.569, 0....
    ## 26 0.12978981 H1-6, H3....       convergent            46.51163 NA, NA, ....
    ## 27 0.12978981 MTOR, RH....        divergent            36.36364 0.185, 0....
    ## 28 0.13977450 CACNB2, ....       convergent            37.28814 0.692, 0....
    ## 29 0.13977450 MAPT, SR....       convergent            43.47826 0.569, 0....
    ## 30 0.13977450 WNT3, PL....       convergent            40.00000 0.404, 0....
    ## 31 0.13977450 MAPT, SR....        divergent            32.46753 0.569, 0....
    ## 32 0.14404240 MAPT, WN....       convergent            46.80851 0.569, 0....
    ## 33 0.14404240 WNT3, SE....       convergent            28.57143 0.404, 0....
    ## 34 0.14404240 PLXNA4, ....       convergent            80.00000 0.229, 0....
    ## 35 0.14404240 MAPT, WN....       convergent            31.25000 0.569, 0....
    ## 36 0.14404240 CACNB2, ....       convergent            25.80645 0.692, 0....
    ## 37 0.18748886 WNT3, MI....        divergent            31.70732 0.404, N....
    ## 38 0.19831313 CACNB2, ....       convergent            27.77778 0.692, 0....

``` r
#explore results
classify_genes %>% arrange(desc(percent_constrained)) # arrange by constraint desc
```

    ##    variant
    ## 1     rare
    ## 2     rare
    ## 3     rare
    ## 4   common
    ## 5     rare
    ## 6     rare
    ## 7     rare
    ## 8     rare
    ## 9     rare
    ## 10    rare
    ## 11    rare
    ## 12    rare
    ## 13    rare
    ## 14    rare
    ## 15  common
    ## 16    rare
    ## 17    rare
    ## 18    rare
    ## 19  common
    ## 20    rare
    ## 21  common
    ## 22  common
    ## 23  common
    ## 24    rare
    ## 25  common
    ## 26  common
    ## 27  common
    ## 28  common
    ## 29  common
    ## 30  common
    ## 31  common
    ## 32  common
    ## 33    rare
    ## 34  common
    ## 35  common
    ## 36  common
    ## 37  common
    ## 38  common
    ##                                                                    pathway
    ## 1                                                        GOBP_NEURON_DEATH
    ## 2                                                GOBP_CHROMATIN_REMODELING
    ## 3  GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE
    ## 4               GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS
    ## 5                                 GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS
    ## 6                                          GOBP_REGULATION_OF_NEUROGENESIS
    ## 7                                                GOBP_SYNAPSE_ORGANIZATION
    ## 8               GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS
    ## 9                                   GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY
    ## 10                                         GOBP_NEURON_PROJECTION_GUIDANCE
    ## 11                      GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION
    ## 12                                GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH
    ## 13                                         GOBP_NEUROTRANSMITTER_SECRETION
    ## 14                                GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS
    ## 15                                               GOBP_CHROMATIN_REMODELING
    ## 16                                        GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS
    ## 17                              GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS
    ## 18               GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT
    ## 19                      GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION
    ## 20                              GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE
    ## 21                                GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS
    ## 22                                  GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY
    ## 23                                GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH
    ## 24                                             KEGG_MTOR_SIGNALING_PATHWAY
    ## 25                                         GOBP_REGULATION_OF_NEUROGENESIS
    ## 26                                               GOBP_SYNAPSE_ORGANIZATION
    ## 27                                         GOBP_NEURON_PROJECTION_GUIDANCE
    ## 28                              GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE
    ## 29                                             KEGG_MTOR_SIGNALING_PATHWAY
    ## 30                                                       GOBP_NEURON_DEATH
    ## 31               GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT
    ## 32                                GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS
    ## 33    GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE
    ## 34 GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE
    ## 35                                         GOBP_NEUROTRANSMITTER_SECRETION
    ## 36    GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE
    ## 37                                        GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS
    ## 38                              GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS
    ##          padj  leadingEdge converge_diverge percent_constrained        LOEUF
    ## 1  0.08450704 SCN2A, S....        divergent            86.66667 0.127, 0....
    ## 2  0.09671180 CHD8, KD....        divergent            85.71429 0.082, 0....
    ## 3  0.10606061 DSCAM, S....       convergent            85.71429 0.189, 0....
    ## 4  0.14404240 PLXNA4, ....       convergent            80.00000 0.229, 0....
    ## 5  0.09208973 ADNP, SH....       convergent            78.78788 0.123, 0....
    ## 6  0.10363636 SYNGAP1,....        divergent            77.77778 0.052, 0....
    ## 7  0.09760766 SYNGAP1,....        divergent            75.86207 0.052, 0....
    ## 8  0.12884753 EPHB1, N....       convergent            75.00000 0.257, 0....
    ## 9  0.10606061 TLK2, DN....       convergent            74.19355 0.114, 1....
    ## 10 0.09671180 DSCAM, N....       convergent            68.18182 0.189, 0....
    ## 11 0.08450704 PTEN, SH....       convergent            68.00000 0.507, 0....
    ## 12 0.08450704 GRIN2B, ....       convergent            58.33333 0.061, 0....
    ## 13 0.12800000 NRXN1, S....       convergent            58.06452 0.254, 0....
    ## 14 0.08450704 ADNP, DS....       convergent            55.55556 0.123, 0....
    ## 15 0.02119885 HDAC4, B....        divergent            53.84615 0.052, 0....
    ## 16 0.19524100 STXBP1, ....       convergent            52.00000 0.086, 0....
    ## 17 0.09671180 SLC6A1, ....       convergent            51.16279 0.15, 0.....
    ## 18 0.08450704 SYNGAP1,....        divergent            50.00000 0.052, 0....
    ## 19 0.12821204 MAPT, PL....       convergent            50.00000 0.569, 0....
    ## 20 0.09557110 PTEN, CT....       convergent            49.09091 0.507, 0....
    ## 21 0.14404240 MAPT, WN....       convergent            46.80851 0.569, 0....
    ## 22 0.12978981 H1-6, H3....       convergent            46.51163 NA, NA, ....
    ## 23 0.13977450 MAPT, SR....       convergent            43.47826 0.569, 0....
    ## 24 0.19524100 AKT2, RI....        divergent            42.85714 0.417, 0....
    ## 25 0.12821204 MAPT, WN....        divergent            41.50943 0.569, 0....
    ## 26 0.10111381 MAPT, CH....        divergent            40.32258 0.569, 0....
    ## 27 0.13977450 WNT3, PL....       convergent            40.00000 0.404, 0....
    ## 28 0.13977450 CACNB2, ....       convergent            37.28814 0.692, 0....
    ## 29 0.12978981 MTOR, RH....        divergent            36.36364 0.185, 0....
    ## 30 0.13977450 MAPT, SR....        divergent            32.46753 0.569, 0....
    ## 31 0.18748886 WNT3, MI....        divergent            31.70732 0.404, N....
    ## 32 0.14404240 MAPT, WN....       convergent            31.25000 0.569, 0....
    ## 33 0.18808777 SLIT1, S....       convergent            30.76923 0.222, 0....
    ## 34 0.06722425 WNT3, PL....       convergent            30.00000 0.404, 0....
    ## 35 0.12821204 CACNB2, ....       convergent            28.57143 0.692, 0....
    ## 36 0.14404240 WNT3, SE....       convergent            28.57143 0.404, 0....
    ## 37 0.19831313 CACNB2, ....       convergent            27.77778 0.692, 0....
    ## 38 0.14404240 CACNB2, ....       convergent            25.80645 0.692, 0....

``` r
classify_genes %>% arrange(desc(percent_constrained)) %>% filter(converge_diverge == "convergent")
```

    ##    variant
    ## 1     rare
    ## 2   common
    ## 3     rare
    ## 4     rare
    ## 5     rare
    ## 6     rare
    ## 7     rare
    ## 8     rare
    ## 9     rare
    ## 10    rare
    ## 11    rare
    ## 12    rare
    ## 13  common
    ## 14    rare
    ## 15  common
    ## 16  common
    ## 17  common
    ## 18  common
    ## 19  common
    ## 20  common
    ## 21    rare
    ## 22  common
    ## 23  common
    ## 24  common
    ## 25  common
    ## 26  common
    ##                                                                    pathway
    ## 1  GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE
    ## 2               GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS
    ## 3                                 GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS
    ## 4               GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS
    ## 5                                   GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY
    ## 6                                          GOBP_NEURON_PROJECTION_GUIDANCE
    ## 7                       GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION
    ## 8                                 GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH
    ## 9                                          GOBP_NEUROTRANSMITTER_SECRETION
    ## 10                                GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS
    ## 11                                        GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS
    ## 12                              GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS
    ## 13                      GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION
    ## 14                              GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE
    ## 15                                GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS
    ## 16                                  GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY
    ## 17                                GOBP_POSITIVE_REGULATION_OF_NEURON_DEATH
    ## 18                                         GOBP_NEURON_PROJECTION_GUIDANCE
    ## 19                              GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE
    ## 20                                GOBP_POSITIVE_REGULATION_OF_AXONOGENESIS
    ## 21    GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE
    ## 22 GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE
    ## 23                                         GOBP_NEUROTRANSMITTER_SECRETION
    ## 24    GOBP_NEGATIVE_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE
    ## 25                                        GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS
    ## 26                              GOBP_REGULATION_OF_NEUROTRANSMITTER_LEVELS
    ##          padj  leadingEdge converge_diverge percent_constrained        LOEUF
    ## 1  0.10606061 DSCAM, S....       convergent            85.71429 0.189, 0....
    ## 2  0.14404240 PLXNA4, ....       convergent            80.00000 0.229, 0....
    ## 3  0.09208973 ADNP, SH....       convergent            78.78788 0.123, 0....
    ## 4  0.12884753 EPHB1, N....       convergent            75.00000 0.257, 0....
    ## 5  0.10606061 TLK2, DN....       convergent            74.19355 0.114, 1....
    ## 6  0.09671180 DSCAM, N....       convergent            68.18182 0.189, 0....
    ## 7  0.08450704 PTEN, SH....       convergent            68.00000 0.507, 0....
    ## 8  0.08450704 GRIN2B, ....       convergent            58.33333 0.061, 0....
    ## 9  0.12800000 NRXN1, S....       convergent            58.06452 0.254, 0....
    ## 10 0.08450704 ADNP, DS....       convergent            55.55556 0.123, 0....
    ## 11 0.19524100 STXBP1, ....       convergent            52.00000 0.086, 0....
    ## 12 0.09671180 SLC6A1, ....       convergent            51.16279 0.15, 0.....
    ## 13 0.12821204 MAPT, PL....       convergent            50.00000 0.569, 0....
    ## 14 0.09557110 PTEN, CT....       convergent            49.09091 0.507, 0....
    ## 15 0.14404240 MAPT, WN....       convergent            46.80851 0.569, 0....
    ## 16 0.12978981 H1-6, H3....       convergent            46.51163 NA, NA, ....
    ## 17 0.13977450 MAPT, SR....       convergent            43.47826 0.569, 0....
    ## 18 0.13977450 WNT3, PL....       convergent            40.00000 0.404, 0....
    ## 19 0.13977450 CACNB2, ....       convergent            37.28814 0.692, 0....
    ## 20 0.14404240 MAPT, WN....       convergent            31.25000 0.569, 0....
    ## 21 0.18808777 SLIT1, S....       convergent            30.76923 0.222, 0....
    ## 22 0.06722425 WNT3, PL....       convergent            30.00000 0.404, 0....
    ## 23 0.12821204 CACNB2, ....       convergent            28.57143 0.692, 0....
    ## 24 0.14404240 WNT3, SE....       convergent            28.57143 0.404, 0....
    ## 25 0.19831313 CACNB2, ....       convergent            27.77778 0.692, 0....
    ## 26 0.14404240 CACNB2, ....       convergent            25.80645 0.692, 0....

``` r
classify_genes %>% arrange(desc(percent_constrained)) %>% filter(converge_diverge == "divergent")
```

    ##    variant                                                   pathway       padj
    ## 1     rare                                         GOBP_NEURON_DEATH 0.08450704
    ## 2     rare                                 GOBP_CHROMATIN_REMODELING 0.09671180
    ## 3     rare                           GOBP_REGULATION_OF_NEUROGENESIS 0.10363636
    ## 4     rare                                 GOBP_SYNAPSE_ORGANIZATION 0.09760766
    ## 5   common                                 GOBP_CHROMATIN_REMODELING 0.02119885
    ## 6     rare GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT 0.08450704
    ## 7     rare                               KEGG_MTOR_SIGNALING_PATHWAY 0.19524100
    ## 8   common                           GOBP_REGULATION_OF_NEUROGENESIS 0.12821204
    ## 9   common                                 GOBP_SYNAPSE_ORGANIZATION 0.10111381
    ## 10  common                               KEGG_MTOR_SIGNALING_PATHWAY 0.12978981
    ## 11  common                                         GOBP_NEURON_DEATH 0.13977450
    ## 12  common GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT 0.18748886
    ##     leadingEdge converge_diverge percent_constrained        LOEUF
    ## 1  SCN2A, S....        divergent            86.66667 0.127, 0....
    ## 2  CHD8, KD....        divergent            85.71429 0.082, 0....
    ## 3  SYNGAP1,....        divergent            77.77778 0.052, 0....
    ## 4  SYNGAP1,....        divergent            75.86207 0.052, 0....
    ## 5  HDAC4, B....        divergent            53.84615 0.052, 0....
    ## 6  SYNGAP1,....        divergent            50.00000 0.052, 0....
    ## 7  AKT2, RI....        divergent            42.85714 0.417, 0....
    ## 8  MAPT, WN....        divergent            41.50943 0.569, 0....
    ## 9  MAPT, CH....        divergent            40.32258 0.569, 0....
    ## 10 MTOR, RH....        divergent            36.36364 0.185, 0....
    ## 11 MAPT, SR....        divergent            32.46753 0.569, 0....
    ## 12 WNT3, MI....        divergent            31.70732 0.404, N....

``` r
collapsed_classified_genes <- classify_genes %>%
  dplyr::select(-LOEUF) %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

#write.csv(collapsed_classified_genes, "../../results/loeuf_anno_results.csv", row.names = FALSE)
```

## Cytoscape File Output

``` r
#select top constrained pathways
top <- classify_genes %>% arrange(desc(percent_constrained)) %>% distinct(pathway, .keep_all = TRUE) %>% head(10)

#pathway level metadata
top_node_attributes <- top %>% 
  dplyr::select(pathway, variant, converge_diverge, percent_constrained) %>%
  distinct() %>% 
  mutate(name = paste0(pathway, "_", variant)) #unique ID for pathway

#edge list (pathway-variant → gene)
top_edge_list <- top %>% 
  dplyr::select(pathway, variant, leadingEdge) %>%
  unnest(leadingEdge) %>% 
  mutate(
    source = paste0(pathway, "_", variant),
    target = leadingEdge) %>%
  dplyr::select(source, target)

#pathway nodes
pathway_nodes <- top_edge_list %>%
  distinct(source) %>%
  dplyr::rename(name = source) %>%
  mutate(node_type = "pathway") %>%
  left_join(top_node_attributes, by = "name")

#gene nodes
gene_nodes <- top_edge_list %>%
  distinct(target) %>%
  dplyr::rename(name = target) %>%
  mutate(node_type = "gene")

#combine nodes
top_all_nodes <- dplyr::bind_rows(pathway_nodes, gene_nodes)

#write CSVs
#write.csv(top_all_nodes, "../../results/cyto_top_all_nodes.csv", row.names = FALSE, quote = FALSE)
#write.csv(top_edge_list, "../../results/cyto_top_edge_list.csv", row.names = FALSE, quote = FALSE)
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
    ##  [1] org.Hs.eg.db_3.20.0  AnnotationDbi_1.68.0 IRanges_2.40.1      
    ##  [4] S4Vectors_0.44.0     Biobase_2.66.0       BiocGenerics_0.52.0 
    ##  [7] lubridate_1.9.4      forcats_1.0.0        stringr_1.5.1       
    ## [10] dplyr_1.1.4          purrr_1.0.4          readr_2.1.5         
    ## [13] tidyr_1.3.1          tibble_3.2.1         ggplot2_3.5.1       
    ## [16] tidyverse_2.0.0     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] generics_0.1.4          RSQLite_2.3.9           stringi_1.8.7          
    ##  [4] hms_1.1.3               digest_0.6.39           magrittr_2.0.4         
    ##  [7] evaluate_1.0.5          grid_4.4.3              timechange_0.3.0       
    ## [10] fastmap_1.2.0           blob_1.2.4              jsonlite_2.0.0         
    ## [13] GenomeInfoDb_1.42.3     DBI_1.2.3               httr_1.4.7             
    ## [16] UCSC.utils_1.2.0        scales_1.3.0            Biostrings_2.74.1      
    ## [19] cli_3.6.5               crayon_1.5.3            rlang_1.1.6            
    ## [22] XVector_0.46.0          bit64_4.6.0-1           munsell_0.5.1          
    ## [25] cachem_1.1.0            withr_3.0.2             yaml_2.3.11            
    ## [28] tools_4.4.3             tzdb_0.4.0              memoise_2.0.1          
    ## [31] colorspace_2.1-2        GenomeInfoDbData_1.2.13 png_0.1-8              
    ## [34] vctrs_0.6.5             R6_2.6.1                lifecycle_1.0.4        
    ## [37] zlibbioc_1.52.0         KEGGREST_1.46.0         bit_4.6.0              
    ## [40] pkgconfig_2.0.3         pillar_1.10.1           gtable_0.3.6           
    ## [43] glue_1.8.0              xfun_0.54               tidyselect_1.2.1       
    ## [46] rstudioapi_0.17.1       knitr_1.49              htmltools_0.5.8.1      
    ## [49] rmarkdown_2.29          compiler_4.4.3
