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
    ##          padj
    ## 1  0.08450704
    ## 2  0.08450704
    ## 3  0.08450704
    ## 4  0.08450704
    ## 5  0.08450704
    ## 6  0.09208973
    ## 7  0.09557110
    ## 8  0.09671180
    ## 9  0.09671180
    ## 10 0.09671180
    ## 11 0.09760766
    ## 12 0.10363636
    ## 13 0.10606061
    ## 14 0.10606061
    ## 15 0.12800000
    ## 16 0.12884753
    ## 17 0.18808777
    ## 18 0.19524100
    ## 19 0.19524100
    ## 20 0.02119885
    ## 21 0.06722425
    ## 22 0.10111381
    ## 23 0.12821204
    ## 24 0.12821204
    ## 25 0.12821204
    ## 26 0.12978981
    ## 27 0.12978981
    ## 28 0.13977450
    ## 29 0.13977450
    ## 30 0.13977450
    ## 31 0.13977450
    ## 32 0.14404240
    ## 33 0.14404240
    ## 34 0.14404240
    ## 35 0.14404240
    ## 36 0.14404240
    ## 37 0.18748886
    ## 38 0.19831313
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              leadingEdge
    ## 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           6326;8831;23394;2904;1499;6812;2561;11151;2047;4763;7068;5609;3326;2932;7314
    ## 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             2904;1499;4763;5609;2932;6733;2213;468;2353;2896;2901;5621
    ## 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   8831;5728;53335;2670
    ## 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 5728;85358;25942;10716;8861;2047;3326;7314;23314;51684;2560;4781;2048;4929;388585;3320;55079;64843;9693;5881;6095;79625;2909;1630;4133
    ## 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         23394;1826;22902;5362;5361;5604;6794;3984;5364
    ## 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         23394;85358;1499;1826;2670;22902;5649;5362;3706;23316;2048;5501;5789;5361;4076;5604;6794;6663;2475;3984;5364;23405;10097;7424;7520;3091;5803;3815;7074;10507;27020;23654;81565
    ## 7                                                                                                                                                                                                                                                                                                                                                                                                           5728;1499;6812;27445;2932;5575;6844;2580;1173;8943;84446;6857;6571;2901;22999;9892;11069;3208;6455;1644;8541;107;6453;55737;8546;7054;127833;57084;821;146395;89780;9024;1268;6809;85439;9900;57555;9066;134957;3690;1855;8618;5579;8867;23025;320;84162;50618;60;23191;1814;154664;114781;26052;8650;120892
    ## 8                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          57680;23135;10765;6601;6304;6595;54891;10664;2625;56916;29994;23314;6605;9682
    ## 9                                                                                                                                                                                                                                                                                                                                                                                                                                                                      1826;9378;1808;10716;4089;2047;7204;2625;5649;1272;5362;23768;4781;137970;1889;2048;5361;6585;128434;55079;1949;64843;6405;351;5364;11313;64855;1943;5459;2909;1630;9252;10507;2534;27020;23654;6586;10381;6477;126147;8609;64919;89780;5015;3799
    ## 10                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           6529;9378;6812;2670;4763;27445;2932;43;815;6844;41;84446;6857;6571;8864;2901;22999;8541;107;7054;6539;6506;127833;2852;1268;6809;9900;9066;11255;6538;6531;6530;134957;3690;60482;1855;4208;8618;5579;8867;23025;5799;79772
    ## 11                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              8831;23394;6529;2904;5728;85358;1499;1826;55209;22941;9378;2562;57689;2561;23181;26115;2047;84623;64101;27445;9762;43;26045;718;5649;4774;288;2554;10152
    ## 12                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        8831;23394;5728;85358;1499;1826;6497;2670;4763
    ## 13                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    11011;1788;25942;6601;5253;10664;6605;23304;54737;10036;146059;10014;84456;51773;3020;8085;6749;9646;7994;23522;2146;9874;6599;9219;79723;3065;57673;8970;8932;8467;196528;85236;10155;55196;23523
    ## 14                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  1826;6585;6405;10507;6586;89780;9723
    ## 15                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       9378;6812;4763;27445;2932;815;6844;41;84446;6857;2901;22999;8541;107;127833;2852;1268;6809;9900;9066;11255;134957;1855;4208;8618;5579;8867;23025;5799;79772;320
    ## 16                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             2047;4781;2048;4929;1630;9201;25987;23334
    ## 17                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            6585;6405;10507;89780;9723;54437;9037;223117;56920;10500;10505;10509;57715
    ## 18                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                6812;27445;2932;6844;6857;2901;22999;8541;107;127833;1268;6809;9900;9066;134957;1855;8618;5579;8867;23025;320;120892;440279;9381;10497
    ## 19                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                208;253260;6794;2475;5595;6195;7424;3091;5295;64223;673;5296;6196;5562
    ## 20                                                                                                                                                                                                                                                                                                                                                                9759;2186;9557;3010;55193;8358;3720;8364;440689;8344;8347;283899;3159;8349;4171;8352;55689;473;10856;8365;81611;7994;8970;6599;125476;8351;196528;10765;8607;1105;29994;3020;3005;3169;4609;3006;9682;8339;4798;50511;55636;8354;201161;1616;8343;55693;388951;8089;54815;3066;85236;3009;54108;51773;54069;8290;9031;4676;132243;8362;8366;6595;29844
    ## 21                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    7473;91584;10507;6259;56920;214;8482;10500;7422;9353;9037;6405;54437;1954;10509;57715;8828;7869;223117;80031;10505
    ## 22 4137;1141;783;10082;221935;7855;2045;784;2915;2561;257194;1605;8927;5300;10507;2047;8443;3643;387;5978;1501;4774;1639;10163;23025;64101;1000;1496;8541;7804;729956;9148;441381;8500;10160;5792;55203;10939;8976;23237;79012;22891;78999;6259;8934;5048;50944;6622;4776;4747;2668;22865;2048;55209;3586;23363;266727;57549;22941;23284;27253;3084;1524;1488;84062;158038;1072;9463;6900;71;400916;785;4915;200933;5641;79953;991;10152;2213;1499;5728;2904;8404;5663;9419;6623;6405;1140;8506;58489;7143;1946;7415;119;10882;7534;5216;26115;643866;1282;83992;1739;8828;8898;2555;56130;9231;1742;375790;1006;22859;5621;10288;23557;11202;85358;4897;6640;5818;147381;7162;4744;389941;1399;9379;627;1981;4761;10369
    ## 23                                                                                                 4137;7473;91584;5080;9001;26022;2045;2915;1605;22954;10507;2475;5978;283659;6091;3400;6009;6474;5324;7804;23566;9148;9798;80128;80237;5530;655;23376;8867;55179;4609;22891;26039;6259;644943;4763;5048;56920;9682;5362;4776;54567;4747;55636;10636;2048;6653;23363;1869;8482;5309;10500;3066;388585;1524;7422;6663;9353;51199;4915;9037;79791;259266;23013;7732;1499;5728;5501;23011;85440;5663;391723;3815;6405;1636;7143;390992;1946;22933;2290;6833;54437;1954;10509;22938;7157;6497;57715;7869;223117;3949;284403;11344;2803;85358;58497;80031;3398;3458;7070;10505;4067;81559;627;1959;57611;4135;5604;5764;1746
    ## 24                                                                                                                                                                                                                                                                                                                                                                                                                                      783;9900;112476;415117;9627;6814;2932;23025;8541;8448;41;10815;1143;8867;57617;10675;127833;4763;116443;6622;6804;5799;6813;80331;8773;116841;1488;815;84062;8618;341359;1138;5579;8676;594855;5071;5663;6623;134957;22895;5868;22930;1136;5028;27345;23557;63908;1268;11255;107
    ## 25                                                                                                                                                                                                                                                                                                                                                                                                       4137;91584;2636;1141;5080;2047;6695;3400;64211;8022;3975;166614;282991;57688;3670;6095;7289;5048;54567;26468;2048;51286;57178;266727;79809;388585;3236;9353;51199;6900;8861;4915;200933;5728;26058;2624;1009;4781;5663;401;6886;9760;5629;2909;10716;2290;64843;2263;9331;8828;2055;85358;23017;2016;6656;25809
    ## 26                                                                                                                                                                                                                                                                           3010;8358;8364;440689;8344;8312;8347;3159;8349;4171;8352;8365;1786;7994;8970;23411;6599;1487;8351;196528;3020;3005;3006;8479;8339;8354;254048;57683;201161;1616;8343;130951;388951;55196;85236;3009;54108;51773;54069;8290;23347;23272;4676;132243;8362;8366;353274;6886;55183;51535;8345;22933;8353;4869;7157;22893;6418;7290;8361;55723;8369;5290;54107;83989;8368;53615;11177;3024;55201;8341;6604;25942;11240;4673;54623;3065;4678;3017
    ## 27                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      2475;6009;25989;673;9470;7422;6195;5595;208;6194;5563;5296;8503;5291;6198;64223;5290;7423;81617;57521;54541;6199
    ## 28                                                                                                                                                                                                                                                                                                                                                                               783;9900;112476;415117;1759;8927;1785;11069;2932;9829;23025;1000;8541;8448;10815;8871;8867;57617;6453;23237;10675;127833;116443;6622;6804;246213;6093;80331;8773;116841;27253;1488;84062;8618;341359;1138;71;5579;79953;80852;8676;1499;5728;29993;594855;5663;6623;6517;134957;821;1644;8546;5868;22930;1173;5028;23557;63908;1268;107
    ## 29                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        4137;6733;7021;2045;4296;387;55294;5978;841;2932;3672;4763;6622;142;7042;7099;1616;5309;581;835;2213;1499;2904
    ## 30                                                                                                                                                                                                                                                                                                                                             7473;91584;2636;5080;5879;9211;2045;1605;10507;2047;55740;6091;56956;2044;8022;1969;3975;655;3670;1400;6259;5880;56920;214;5362;729920;2668;2048;23032;8482;57549;2131;2041;10500;284217;3798;7422;9353;6900;9037;6477;1908;4781;6405;2909;10716;9369;7143;1946;91624;2290;54437;64843;1954;9638;10509;57715;8828;7869;223117;283297;57408;2625;399;4897;5818;80031;10505
    ## 31                                                                                                                                                                                                                                                                                            4137;6733;7021;840;2045;6647;83741;3134;4296;2561;3812;2047;10525;116986;387;3745;84678;55294;5978;841;2932;10550;3672;23411;27023;1718;353116;23321;27242;3670;8934;4763;6416;57465;6622;4776;4747;142;7042;3632;2668;9529;947;1113;2876;664;7099;6653;1616;3586;2632;6093;5309;80331;581;6326;1524;112398;91304;84062;2020;835;4915;5641;9962;1522;9244;2213;1499;6477;6548;23400;1271;5071;2904;4544;5663;8883;208;6623
    ## 32                                                                                                                                                                                                                                                                                                                                                                                                                                               4137;7473;91584;5080;9001;2915;1605;22954;10507;2475;6091;6009;6474;5324;7804;23566;9148;9798;80237;23376;8867;55179;4609;22891;26039;644943;5048;5362;4747;10636;2048;23363;1869;8482;3066;1524;7422;6663;9353;51199;4915;9037;79791;259266;23013;7732;1499;5501;23011
    ## 33                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  7473;10507;6259;56920;8482;10500;9037;6405;54437;10509;57715;7869;223117;80031;10505
    ## 34                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 91584;1141;2047;5048;2048;9353;51199;200933;1009;4781
    ## 35                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              4137;7473;91584;10507;6091;6474;23566;9798;5048;5362;4747;8482;7422;9353;51199;4915;9037
    ## 36                                                                                                                                                                                                                                                                                                                                                                783;9900;221955;112476;415117;9627;6814;2932;23025;8541;6582;8448;41;10815;1143;8867;6538;57617;4842;10675;127833;5153;145270;4763;116443;6622;6804;2668;5799;6813;80331;8773;116841;1488;815;84062;8618;341359;1138;1312;5579;8676;594855;5071;5663;6623;60482;134957;22895;5053;5868;22930;1136;2555;2055;5028;27345;6581;23557;63908;1268;11255;107
    ## 37                                                                                                                                                                                                                                                                                                                                                                                                                                                               7473;23251;2045;10507;387;6695;80128;1400;22891;6259;2672;8934;53335;5048;56920;6911;9529;2288;2048;8482;10500;3066;84062;9037;5728;29116;157285;5663;6405;59341;7143;54437;10509;112;1601;129807;57715;7869;223117;349667;23557;11202;80031;7070;10505
    ## 38                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            783;9900;112476;415117;2932;23025;8541;8448;10815;8867;57617;10675;127833;116443;6622;6804;80331;8773;116841;1488;84062;8618;341359;1138;5579;8676;594855;5663;134957;5868;22930;5028;23557;63908;1268;107
    ##    converge_diverge
    ## 1         divergent
    ## 2        convergent
    ## 3         divergent
    ## 4        convergent
    ## 5        convergent
    ## 6        convergent
    ## 7        convergent
    ## 8         divergent
    ## 9        convergent
    ## 10       convergent
    ## 11        divergent
    ## 12        divergent
    ## 13       convergent
    ## 14       convergent
    ## 15       convergent
    ## 16       convergent
    ## 17       convergent
    ## 18       convergent
    ## 19        divergent
    ## 20        divergent
    ## 21       convergent
    ## 22        divergent
    ## 23        divergent
    ## 24       convergent
    ## 25       convergent
    ## 26       convergent
    ## 27        divergent
    ## 28       convergent
    ## 29       convergent
    ## 30       convergent
    ## 31        divergent
    ## 32       convergent
    ## 33       convergent
    ## 34       convergent
    ## 35       convergent
    ## 36       convergent
    ## 37        divergent
    ## 38       convergent

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
loeuf
```

    ##              MED13              NIPBL               SMC3              CNOT1 
    ##              0.030              0.032              0.037              0.038 
    ##                RLF              PCF11             FNDC3B               TAF1 
    ##              0.040              0.040              0.042              0.043 
    ##               RSF1             NCKAP1              KDM2A               BRD4 
    ##              0.044              0.044              0.045              0.048 
    ##               HELZ               FBN1               XPO1              PRR12 
    ##              0.048              0.049              0.051              0.051 
    ##              USP9X              POLA1            SYNGAP1              PRKDC 
    ##              0.051              0.051              0.052              0.052 
    ##              HDAC4               SMG1              ZC3H4             COL5A1 
    ##              0.052              0.053              0.054              0.055 
    ##            SMARCA4              TNPO1               AGO1           ARHGAP35 
    ##              0.055              0.056              0.057              0.057 
    ##               LRP1               TOP1             TRIP12              KMT2E 
    ##              0.058              0.059              0.060              0.060 
    ##              HCFC1               UBTF              HUWE1              KDM3B 
    ##              0.060              0.060              0.060              0.060 
    ##              TRRAP             GRIN2B               USP7             ATP1A3 
    ##              0.060              0.061              0.061              0.062 
    ##              ASH1L            ANKRD17              SMC1A             NPEPPS 
    ##              0.062              0.062              0.062              0.062 
    ##              CHERP             MED13L              RPRD2              KMT2A 
    ##              0.062              0.064              0.064              0.065 
    ##             CREBBP              SF3B1              EIF3B             PPFIA1 
    ##              0.066              0.066              0.067              0.067 
    ##             ZSWIM6              PRPF3              KAT6A              KMT2B 
    ##              0.067              0.068              0.069              0.070 
    ##              SIN3A               CHD2             ARID1A           GLTSCR1L 
    ##              0.070              0.070              0.071              0.071 
    ##              SCN1A              MED12               TJP1               SPEN 
    ##              0.071              0.071              0.071              0.071 
    ##              WDFY3              HDAC6              EHMT1             JMJD1C 
    ##              0.071              0.072              0.072              0.072 
    ##            MAP3K12             GGNBP2               CASK               TLE3 
    ##              0.072              0.073              0.073              0.073 
    ##               CHD6               UBR5             EIF4G1               NBEA 
    ##              0.073              0.073              0.074              0.074 
    ##            SMARCA5               TSC2             SUPT5H             GREB1L 
    ##              0.074              0.074              0.074              0.075 
    ##              PTCH1               CHD7            RAD54L2              CNOT3 
    ##              0.075              0.076              0.076              0.076 
    ##              AHDC1             U2SURP              RBM25             GIGYF2 
    ##              0.076              0.077              0.077              0.077 
    ##              MBTD1              DCAF5           SUV420H1            RABGAP1 
    ##              0.077              0.078              0.078              0.078 
    ##              KHSRP             TM9SF2              MATR3             PHF21A 
    ##              0.078              0.079              0.079              0.080 
    ##               NAV1             HIVEP2                OGT              RBBP6 
    ##              0.080              0.080              0.080              0.080 
    ##             GTPBP4              MYBL1              TEX11            DYNC1H1 
    ##              0.080              0.080              0.080              0.080 
    ##              SBNO1              THOC2              PRKCB              INO80 
    ##              0.081              0.081              0.081              0.081 
    ##               TAF4               FLNA             TRIM28               CHD8 
    ##              0.081              0.082              0.082              0.082 
    ##              RBM33           SNRNP200               XPO6           SMARCAD1 
    ##              0.082              0.082              0.082              0.083 
    ##             FRMPD4              MACF1                WAC              BRWD3 
    ##              0.083              0.084              0.084              0.084 
    ##          C17orf104             ATP2B1             EIF4G2             ZNF462 
    ##              0.084              0.084              0.085              0.085 
    ##              WDR26              CUL4B             TRIM33             STXBP1 
    ##              0.085              0.085              0.086              0.086 
    ##             TOPBP1                NF2               TET3                MN1 
    ##              0.086              0.086              0.086              0.087 
    ##             SUPT6H       RP11-762I7.5                NHS            RPS6KA3 
    ##              0.087              0.087              0.087              0.087 
    ##              MYT1L               UBP1              PLCL2             SPTBN1 
    ##              0.087              0.088              0.088              0.088 
    ##               MYH9             ZDHHC5             LUC7L3              COPB2 
    ##              0.088              0.088              0.089              0.089 
    ##            DNAJC14              RNF17               ANK3                 F8 
    ##              0.089              0.089              0.089              0.089 
    ##             TM9SF3               DLG3                MGA            ANKRD52 
    ##              0.090              0.090              0.090              0.090 
    ##              STAG2              SATB2              USP34             MAP3K2 
    ##              0.090              0.091              0.091              0.091 
    ##              PREX1              MED14             CYFIP2             ZNF407 
    ##              0.091              0.092              0.092              0.092 
    ##              KPNB1               ETV5               CUL2             HNRNPM 
    ##              0.092              0.092              0.093              0.093 
    ##              PIAS1              DAAM2             EFTUD2             YTHDC1 
    ##              0.094              0.094              0.094              0.094 
    ##              NCBP1              TCF20               DLL1               UBR4 
    ##              0.094              0.095              0.095              0.095 
    ##              STAT3            IGF2BP3               NSD1              SRCAP 
    ##              0.095              0.095              0.095              0.096 
    ##             ZNF777             HNRNPK              SCAF8              ARID2 
    ##              0.096              0.096              0.096              0.096 
    ##               CCT3             COL5A2               PHC3            GATAD2B 
    ##              0.097              0.097              0.097              0.097 
    ##             FBXO11             NOTCH1             COL3A1             PPP3CB 
    ##              0.097              0.097              0.098              0.098 
    ##               ILF3            PLEKHA5             ZNF521               TAB2 
    ##              0.098              0.098              0.098              0.098 
    ##             CTNND2               DDX5              ACTN4               UPF2 
    ##              0.098              0.098              0.099              0.099 
    ##              CARM1             PAXIP1               BPTF              EP300 
    ##              0.099              0.099              0.099              0.099 
    ##             FBXL19             CREBRF              HDAC2                SP1 
    ##              0.100              0.100              0.100              0.100 
    ##           MAPK8IP2              AP2A1              ARCN1             SPTAN1 
    ##              0.100              0.101              0.101              0.101 
    ##              RPTOR              CAND1              COPS2              SALL4 
    ##              0.101              0.101              0.101              0.101 
    ##            CACNA1C             ARID1B              PSMD2              PSMD3 
    ##              0.102              0.102              0.102              0.102 
    ##              KMT2D             UBAP2L             ZNF148                BTK 
    ##              0.103              0.103              0.103              0.103 
    ##               RRM1            PPP2R5E               RFX3              BIRC6 
    ##              0.103              0.103              0.104              0.104 
    ##              LMNB2               CLTC               COPA             COL1A1 
    ##              0.104              0.104              0.105              0.105 
    ##              DHX15               DHX9               MAP2              ZMYM3 
    ##              0.105              0.105              0.105              0.106 
    ##               AGO2               ANK2            RALGAPB               IPO7 
    ##              0.106              0.106              0.106              0.106 
    ##             PRPF19              SCAF4            TBL1XR1            ZCCHC11 
    ##              0.106              0.106              0.106              0.107 
    ##             ZCCHC6            ANKRD11           PAFAH1B1               ZEB2 
    ##              0.107              0.107              0.107              0.107 
    ##              AP2A2               SMU1              NAA25               UBA1 
    ##              0.108              0.108              0.108              0.108 
    ##              BAZ1B            HNRNPH1              VPRBP               CTIF 
    ##              0.109              0.109              0.109              0.109 
    ##              ABCB7              UTP15             SETBP1             ARID5B 
    ##              0.109              0.110              0.110              0.110 
    ##              RBM10                TPR              RCOR1               SUFU 
    ##              0.110              0.111              0.111              0.111 
    ##              TRPS1               MTA1             RUVBL2             EIF4G3 
    ##              0.111              0.111              0.111              0.112 
    ##               MDM4              NTRK2             ANKHD1    ANKHD1-EIF4EBP3 
    ##              0.112              0.112              0.112              0.112 
    ##               RORB              TIAL1              NRXN3              PTPRD 
    ##              0.112              0.112              0.112              0.112 
    ##              WAPAL              RBM26               RIF1                NLK 
    ##              0.112              0.112              0.113              0.113 
    ##               PHIP              LPHN2             DPYSL5               USP8 
    ##              0.113              0.113              0.113              0.113 
    ##             TIPARP              ZBTB4               ISY1             HNRNPU 
    ##              0.113              0.113              0.114              0.114 
    ##               TLK2             TRERF1              DHX30             MYCBP2 
    ##              0.114              0.114              0.115              0.115 
    ##              FUBP1             TRIM24             ZC3H7A              UBXN7 
    ##              0.115              0.115              0.115              0.116 
    ##            SLC17A7             NOTCH2            SMARCC1               RAI1 
    ##              0.116              0.116              0.116              0.116 
    ##              CSE1L             EIF4A1             PIK3CA              KIF5C 
    ##              0.117              0.117              0.117              0.117 
    ##              REV3L               AGPS             ADAM10               MED1 
    ##              0.117              0.117              0.117              0.117 
    ##            SMARCC2               MYRF              DOT1L             CHAF1A 
    ##              0.117              0.117              0.118              0.118 
    ##              PSMD1              SRRM2              PIAS4              ACTR3 
    ##              0.118              0.118              0.118              0.118 
    ##              DDX3X               TSC1                SON              WHSC1 
    ##              0.118              0.118              0.118              0.119 
    ##              TCEB3              STAG1               ATRX             VCPIP1 
    ##              0.119              0.119              0.119              0.119 
    ##             EIF4A3               MAU2              PSMC2             DAZAP1 
    ##              0.119              0.119              0.119              0.119 
    ##               POGZ            PRKAR1A              SF3A3               FXR1 
    ##              0.119              0.119              0.119              0.119 
    ##           CDC42BPB               XPO7               RARB               RERE 
    ##              0.120              0.120              0.120              0.120 
    ##            ZMYND11               RPL4            SUPT16H             ZBTB47 
    ##              0.120              0.120              0.120              0.120 
    ##             PRRC2C             PPFIA3             ATP1A1              KIF11 
    ##              0.121              0.121              0.121              0.121 
    ##              NCOA2              HIPK1           KIAA0430             ARMCX4 
    ##              0.121              0.121              0.121              0.122 
    ##              KMT2C            PIK3AP1             SEMA4C              CYTH1 
    ##              0.122              0.122              0.122              0.122 
    ##              RBM17               RTF1               IPO5              NUP98 
    ##              0.122              0.122              0.122              0.122 
    ##              TRAF2              EIF4B               ADNP                RET 
    ##              0.123              0.123              0.123              0.123 
    ##              TEAD1             CCDC22             SHANK3            CACNA1E 
    ##              0.123              0.123              0.123              0.124 
    ##              PRMT5              JADE3               WDR5              CUL4A 
    ##              0.124              0.124              0.124              0.124 
    ##               PCLO             ZNF219             GAPVD1               EPC1 
    ##              0.124              0.124              0.125              0.125 
    ##              KCNB1              NLGN2             PCDH19              DDX17 
    ##              0.125              0.126              0.126              0.126 
    ##                RB1               WNK1             ZC3H18              CDK12 
    ##              0.126              0.126              0.126              0.126 
    ##             HECTD4            RALGPS1             TNRC6C              SCN2A 
    ##              0.127              0.127              0.127              0.127 
    ##               SOBP                VCP               TLN1              BTAF1 
    ##              0.127              0.127              0.128              0.128 
    ##              EPHA7             COL4A1               SIK3              KAT6B 
    ##              0.128              0.128              0.128              0.128 
    ##               KAT7                ZFR               GSE1              GSK3A 
    ##              0.128              0.128              0.128              0.128 
    ##              L1CAM             CTNNB1              PBRM1           PRICKLE1 
    ##              0.129              0.129              0.129              0.129 
    ##               PUM1             AHCTF1              TNPO2               AAK1 
    ##              0.129              0.129              0.129              0.129 
    ##              RASA1              CDC5L             FERMT2               WTAP 
    ##              0.129              0.129              0.130              0.130 
    ##               UBR3              SCN8A              FNIP2                MSN 
    ##              0.130              0.130              0.130              0.130 
    ##           PPP1R12A             RANBP9              USP47            CCDC88A 
    ##              0.130              0.130              0.130              0.130 
    ##              MLLT3              EIF3A             YTHDC2               KSR2 
    ##              0.131              0.131              0.131              0.131 
    ##             ZMYND8              CCNA2              PTPN6              SCML2 
    ##              0.131              0.132              0.132              0.132 
    ##               SNRK              MINK1               IRF6              ZRSR2 
    ##              0.132              0.132              0.132              0.132 
    ##               RFX7              CSTF2              NR4A2              TULP4 
    ##              0.132              0.133              0.133              0.133 
    ##              MAP1B              USP42             ZNF592              NPAS2 
    ##              0.133              0.133              0.133              0.133 
    ##             ZSWIM8             IQSEC2              U2AF2              ITPR1 
    ##              0.133              0.133              0.133              0.134 
    ##             COL2A1             NPLOC4           MAPKAPK2              DDX42 
    ##              0.134              0.134              0.134              0.134 
    ##            CACNA1A              SUZ12              PHF12              DACH1 
    ##              0.134              0.134              0.135              0.135 
    ##              PORCN               CTR9               RBM5              NRBP1 
    ##              0.135              0.135              0.135              0.135 
    ##               TBX5             FRMD4A             PTPN11             DLGAP3 
    ##              0.135              0.135              0.135              0.136 
    ##               HIRA             SETD1A              NAA35            ELMSAN1 
    ##              0.136              0.136              0.136              0.136 
    ##             CACTIN            ARHGEF6            FAM120C               RBM6 
    ##              0.136              0.136              0.136              0.136 
    ##              KIF1A              MECOM               MBD5              MEIS1 
    ##              0.136              0.136              0.136              0.136 
    ##              INTS6              ZFHX3               JPH4            SMARCA1 
    ##              0.136              0.136              0.137              0.137 
    ##             HEATR1               BAI3            SLC12A5               TAF3 
    ##              0.137              0.138              0.138              0.138 
    ##              SHOC2              PSMC1              RBM22              NFKB1 
    ##              0.138              0.138              0.138              0.138 
    ##               AFF4             RB1CC1               SAFB            TRAPPC8 
    ##              0.138              0.138              0.138              0.138 
    ##             ARID4A               TRIO             YTHDF2             PPFIA2 
    ##              0.139              0.139              0.139              0.139 
    ##              KDM6B               SPTB               SOS1               BAI1 
    ##              0.140              0.140              0.140              0.140 
    ##             AMBRA1              UBE2O             ZNF292               BCOR 
    ##              0.140              0.140              0.140              0.141 
    ##              HDAC5             PODXL2              AFTPH             RANBP3 
    ##              0.141              0.141              0.141              0.141 
    ##              SPAG9              FNBP4             AGPAT3             RANBP2 
    ##              0.142              0.142              0.142              0.142 
    ##              TKTL1             ZBTB44            WHSC1L1              CNOT4 
    ##              0.142              0.142              0.142              0.143 
    ##               JAG2            TFCP2L1               MAOA             SAMD4B 
    ##              0.143              0.143              0.143              0.143 
    ##             CSNK1E              DNMT1               PHEX             ATXN2L 
    ##              0.143              0.143              0.143              0.143 
    ##               XPO4            ARFGEF1              COPS5               SPOP 
    ##              0.143              0.143              0.144              0.144 
    ##              CSDE1               API5               CUL1              DDX46 
    ##              0.144              0.144              0.144              0.144 
    ##              ZFHX4               GNB1             HOMER1              ZFPM2 
    ##              0.144              0.145              0.145              0.145 
    ##              AEBP2              NUP85              CBLL1              RNF38 
    ##              0.145              0.145              0.145              0.145 
    ##              BAZ2A               JAG1              ASXL2                SP3 
    ##              0.146              0.146              0.146              0.146 
    ##              SRRM1              MAML1               CCNC              KIF23 
    ##              0.146              0.146              0.146              0.146 
    ##               EZH2               RFX4              XRCC5              USP24 
    ##              0.146              0.147              0.147              0.147 
    ##               CHD3            HNRNPA1             PPP1CB               FLT4 
    ##              0.147              0.147              0.147              0.147 
    ##              UNC79             ZNF217              EFHC2              WBP11 
    ##              0.147              0.147              0.147              0.148 
    ##              ITPKB               CTCF             ATP2B2               MAOB 
    ##              0.148              0.148              0.148              0.148 
    ##             KLHDC3              PRMT1              TGFB2              DIP2C 
    ##              0.148              0.148              0.148              0.149 
    ##               NFIX             FBXL17                WAS               GCM1 
    ##              0.149              0.149              0.149              0.149 
    ##               YBX2               TLE4             SLC6A1               EBF3 
    ##              0.150              0.150              0.150              0.150 
    ##           MAPK8IP3              AHSA1               PAN3               PKN2 
    ##              0.150              0.150              0.150              0.150 
    ##              ABCD1             SETD1B              CLCN4            ARHGEF9 
    ##              0.150              0.151              0.151              0.151 
    ##            B4GALT5               MSL1           SERPINC1              KALRN 
    ##              0.151              0.151              0.151              0.152 
    ##              CASZ1             MARCH6               SFPQ             BCORL1 
    ##              0.152              0.152              0.152              0.152 
    ##               CMIP              PPRC1             DDX26B           C11orf30 
    ##              0.152              0.152              0.152              0.153 
    ##             EFEMP1            MAPKAP1             ATP1B1              ELFN1 
    ##              0.153              0.153              0.153              0.153 
    ##             SMCHD1              FBXW2               GGA1            CAMSAP3 
    ##              0.153              0.153              0.153              0.154 
    ##                AHR                DMD            PPP1R10              HDLBP 
    ##              0.154              0.154              0.154              0.154 
    ##             ZNF609           ARHGAP42               OCRL             RICTOR 
    ##              0.154              0.155              0.155              0.155 
    ##             ZNF541               FASN             ATAD2B              PCIF1 
    ##              0.155              0.155              0.155              0.155 
    ##            DENND5B              HOOK3            ARHGEF2               XPR1 
    ##              0.155              0.155              0.156              0.156 
    ##             INTS10               USP3              MARK2              IKZF1 
    ##              0.156              0.156              0.156              0.156 
    ##              CDH11              KDM4A              RREB1               CHD4 
    ##              0.156              0.156              0.156              0.156 
    ##              HIPK2              NR2E1              KPNA4               CHD5 
    ##              0.156              0.156              0.156              0.157 
    ##             PLXND1               RXRA              LRRC7             ZFC3H1 
    ##              0.157              0.157              0.157              0.157 
    ##              UBE4B              DOCK3              DMXL2              PELP1 
    ##              0.157              0.158              0.158              0.158 
    ##              TXLNG             UNC13A              MAT2A             HECTD1 
    ##              0.158              0.158              0.158              0.158 
    ##              KCNQ2              TANC2             PPP3CA             LARP4B 
    ##              0.158              0.158              0.158              0.158 
    ##               EPC2              ROBO2              MAST1             POLR2A 
    ##              0.159              0.159              0.159              0.159 
    ##            GRIPAP1             PSMD11             TNRC6A              ESPL1 
    ##              0.159              0.159              0.159              0.159 
    ##               CUL5               FBRS             SAP130              ALG13 
    ##              0.159              0.159              0.159              0.160 
    ##              SMAD5             GNB2L1              SLIT2              PPM1A 
    ##              0.160              0.160              0.160              0.160 
    ##             PPP4R2              GSPT1            PIP4K2B                ZFX 
    ##              0.160              0.160              0.160              0.160 
    ##             KANSL3               RORC             ZNF800               DPF2 
    ##              0.160              0.161              0.161              0.161 
    ##                APC              FKBP8              OPHN1              USP48 
    ##              0.161              0.161              0.161              0.161 
    ##              KDM6A             CAPZA2               NPM1             PSMD12 
    ##              0.161              0.161              0.162              0.162 
    ##             ZNF131               GPC3             GRAMD4              KIF2A 
    ##              0.162              0.162              0.162              0.162 
    ##               FGD5              CKAP5               UPF1             UTP14A 
    ##              0.162              0.162              0.162              0.162 
    ##               CHD1              ZBTB2              CPSF6            FAM120A 
    ##              0.162              0.162              0.163              0.163 
    ##              WDR43               CUX1             ZNF865           MAPK8IP1 
    ##              0.163              0.163              0.163              0.163 
    ##              CDK19              NFKB2              KDM5A              TEX10 
    ##              0.163              0.163              0.163              0.163 
    ##           KIAA0947               MDM2            CASKIN1                MYC 
    ##              0.164              0.164              0.164              0.164 
    ##               CLK2               PAX5             ZNF608            DSCAML1 
    ##              0.164              0.165              0.165              0.165 
    ##               CCNJ               CCNK               TLK1              LARP1 
    ##              0.165              0.165              0.165              0.165 
    ##              ARIH1              LRRC4             RBFOX3             ZNF629 
    ##              0.166              0.166              0.166              0.166 
    ##              KDM5C               PHF6             DICER1             GRIN2D 
    ##              0.166              0.166              0.166              0.167 
    ##               RPL5              PHF23            SIPA1L1              MYH10 
    ##              0.167              0.167              0.167              0.167 
    ##             GABBR2              CECR2                 F9               PAK7 
    ##              0.167              0.167              0.167              0.167 
    ##              USP11            PPP1R9B             SFSWAP              HNF1B 
    ##              0.167              0.168              0.168              0.168 
    ##               NOP2            TMEM63B              STAU1              ARPC2 
    ##              0.168              0.168              0.168              0.168 
    ##                CSK               SOX9              NCOR2               ANK1 
    ##              0.168              0.168              0.169              0.169 
    ##             ZNF236              ABCA2           PPP1R16B              RC3H2 
    ##              0.169              0.169              0.169              0.169 
    ##               PAX6              CENPI             ZNF638              NACC1 
    ##              0.169              0.169              0.169              0.169 
    ##             PDGFRA              CADM4             ACVR2A              TRIM9 
    ##              0.169              0.169              0.170              0.170 
    ##              DYRK2              TRPC5            C7orf60               MKL2 
    ##              0.170              0.170              0.170              0.170 
    ##            GPATCH8              ESCO1              GRID2               WEE1 
    ##              0.170              0.170              0.170              0.170 
    ##             MBTPS2             PRPF4B             ZNF654               CHD9 
    ##              0.170              0.170              0.170              0.170 
    ##              KIF1B               RNF2              DAGLA                MNT 
    ##              0.171              0.171              0.171              0.171 
    ##              PRKCD              KDM4B              CELF1            ANKRD28 
    ##              0.171              0.171              0.171              0.171 
    ##               MEN1              SYVN1               FZR1              ZFP91 
    ##              0.171              0.171              0.172              0.172 
    ##              ARMC8             CAMTA1              LIN54             TRIM71 
    ##              0.172              0.172              0.172              0.172 
    ##               NKAP              TFAP4             ACVR1B              PSMD4 
    ##              0.172              0.172              0.172              0.172 
    ##              FBLN5               EEF2         ISY1-RAB43              CFLAR 
    ##              0.172              0.172              0.172              0.173 
    ##              OTUD5             HIVEP1             GTF3C1               ENAH 
    ##              0.173              0.173              0.173              0.173 
    ##              FOXO1             TNRC6B              EEF1G              TENM4 
    ##              0.173              0.173              0.173              0.173 
    ##               KLF4             SEC24C            CTDSPL2              JADE2 
    ##              0.174              0.174              0.174              0.174 
    ##              PCGF1              TTC28              GATA6              SCN3A 
    ##              0.174              0.174              0.174              0.174 
    ##              PSME3               TET1             NUP153              CELF2 
    ##              0.174              0.174              0.174              0.174 
    ##               DDX6               BAG6              EPHA4              MKLN1 
    ##              0.174              0.174              0.175              0.175 
    ##              ZMYM4            PPP2R5D             ANAPC2              FOXP1 
    ##              0.175              0.175              0.175              0.175 
    ##              TDRD1             MAP4K4              MAPK1               NCDN 
    ##              0.175              0.175              0.175              0.175 
    ##                BSN              ADRM1              NCOR1              CNOT7 
    ##              0.175              0.175              0.175              0.176 
    ##               MTM1              BAZ1A               ABL1            JAKMIP2 
    ##              0.176              0.176              0.176              0.176 
    ##             ZNF516              PTBP1             CLASP1                HTT 
    ##              0.176              0.176              0.176              0.176 
    ##               RELN              BRPF1             R3HDM2              MYO9B 
    ##              0.176              0.176              0.176              0.177 
    ##               PHF2             PAPOLA               ZZZ3          EIF4ENIF1 
    ##              0.177              0.177              0.177              0.177 
    ##              PSME4                NCL                GLA               GPHN 
    ##              0.177              0.177              0.177              0.177 
    ##              BRPF3               MSL3              PRPF8              KDM2B 
    ##              0.177              0.178              0.178              0.178 
    ##               AKT3               HIC2              PNISR              DPP10 
    ##              0.178              0.178              0.178              0.178 
    ##              MGAT5               NKRF             ADRBK1               ZER1 
    ##              0.178              0.178              0.178              0.178 
    ##              EXOC5              SF3B3              SYMPK              SPPL3 
    ##              0.178              0.178              0.179              0.179 
    ##               PKD1                TEK              G3BP2             ZBTB18 
    ##              0.179              0.179              0.179              0.179 
    ##              PTPN2           ARHGEF12              LAS1L               PLK2 
    ##              0.179              0.179              0.179              0.180 
    ##            PITPNM2                IDS                DST               DNM2 
    ##              0.180              0.180              0.180              0.180 
    ##              PDS5B              PTPRT               RPL7               TOX3 
    ##              0.180              0.180              0.181              0.181 
    ##              CENPB             SERBP1             ZRANB2               HCN1 
    ##              0.181              0.181              0.181              0.181 
    ##               RELA              COPS3              TFDP2               WDR7 
    ##              0.181              0.181              0.181              0.181 
    ##               TNIK            ZDHHC17               NOL6              NOVA1 
    ##              0.181              0.181              0.182              0.182 
    ##              LMTK3             ATP5A1             SLC7A3             MLLT10 
    ##              0.182              0.182              0.182              0.182 
    ##              MLLT6     RP11-152F13.10             ATP2C1               RNF8 
    ##              0.182              0.182              0.182              0.183 
    ##                SET              UBE3C              SYCP2                CIC 
    ##              0.183              0.183              0.183              0.183 
    ##               ACAN                SF1              EPHB2            CSNK2A2 
    ##              0.183              0.183              0.184              0.184 
    ##               FAT4             NAP1L1             MED12L              ALAS2 
    ##              0.184              0.184              0.184              0.184 
    ##             CEP350              DDX21              MEIS2               MTF2 
    ##              0.184              0.184              0.184              0.184 
    ##             MAP3K4             SHANK1              KLHL2               MTOR 
    ##              0.184              0.184              0.185              0.185 
    ##              PACS1               GLUL               NAV3              NCAM1 
    ##              0.185              0.185              0.185              0.185 
    ##              SARNP               WWC3             LRRC41             CTNND1 
    ##              0.185              0.185              0.185              0.185 
    ##              LPHN3             CNKSR2               FXR2                SYK 
    ##              0.185              0.185              0.185              0.185 
    ##             PTPN12              PAPD7              PARP6             EEF1A2 
    ##              0.186              0.186              0.186              0.186 
    ##           HNRNPUL1            ZC3H11A             PRRC2A               XKR6 
    ##              0.186              0.186              0.186              0.186 
    ##             ARID4B              EIF3D              NDST1              PTPRF 
    ##              0.187              0.187              0.187              0.187 
    ##              LPHN1            ARHGEF7              GPKOW                EED 
    ##              0.187              0.187              0.187              0.187 
    ##             ZNF423               CUX2             NFATC2              STAT1 
    ##              0.187              0.187              0.187              0.187 
    ##            TMEM108              DIDO1            PHF20L1             PRDM16 
    ##              0.187              0.187              0.187              0.187 
    ##               DLL4              SEPT7              CDC73             JARID2 
    ##              0.188              0.188              0.188              0.188 
    ##            SNRNP70             PRRC2B               SOX5             GRIN2A 
    ##              0.188              0.188              0.188              0.188 
    ##              PLAG1              ASAP1               CCT5              CNOT2 
    ##              0.188              0.188              0.188              0.189 
    ##             ELAVL2              SUGP2             RSPRY1              USP37 
    ##              0.189              0.189              0.189              0.189 
    ##               NRF1             NFKBIZ              RBBP5              DSCAM 
    ##              0.189              0.189              0.189              0.189 
    ##              TTC7B              PRKD2             PAXBP1                UNK 
    ##              0.189              0.189              0.190              0.190 
    ##           KIAA2018              AP2M1               NOS1               CBX2 
    ##              0.190              0.190              0.190              0.190 
    ##              PPARD               YAP1               ESR1             THRAP3 
    ##              0.190              0.190              0.190              0.190 
    ##             UBE2Q1              RNF31              NYAP1             SCARF2 
    ##              0.190              0.190              0.190              0.190 
    ##          C20orf112              DOCK9             COL4A5              NR2F1 
    ##              0.191              0.191              0.191              0.191 
    ##               WNK3               DKC1              SENP5               LGI1 
    ##              0.191              0.191              0.191              0.191 
    ##              BMPR2            SMARCE1              FURIN             SPTBN4 
    ##              0.191              0.191              0.191              0.191 
    ##              RPL7A              TENM2              ASIC2             ACTR3B 
    ##              0.191              0.191              0.192              0.192 
    ##           ARHGAP31              SNX27              AJAP1              PTPRA 
    ##              0.192              0.192              0.192              0.192 
    ##              GRIK3              MAP1A             CEP170               GDI1 
    ##              0.192              0.192              0.192              0.192 
    ##               XKR4                HPN               TNS3              MTPAP 
    ##              0.192              0.192              0.192              0.193 
    ##              VDAC3             SNAP91               EBF2              TAF5L 
    ##              0.193              0.193              0.193              0.193 
    ##             ZBTB16               AGO3              PPME1              GRIA3 
    ##              0.193              0.193              0.193              0.193 
    ##             NFKBIA              TENM1               CYBB             SCAF11 
    ##              0.193              0.193              0.193              0.193 
    ##              TBX18             SEC24B            ZNF280C               ATN1 
    ##              0.193              0.194              0.194              0.194 
    ##              EXOC8            BHLHE40            SH3KBP1                SKI 
    ##              0.194              0.194              0.194              0.194 
    ##            DCUN1D5              SOGA1                CHM             RBFOX2 
    ##              0.194              0.194              0.194              0.194 
    ##               RAE1             RUVBL1              SMEK1               TBR1 
    ##              0.194              0.194              0.194              0.194 
    ##             PSMD13               GLI3              SGIP1              FOXP3 
    ##              0.194              0.195              0.195              0.195 
    ##           KIAA0232              WDR44             CACNG3              SIN3B 
    ##              0.195              0.195              0.195              0.195 
    ##               ATRN             KLHL15              ACTR2              GLYR1 
    ##              0.195              0.195              0.195              0.195 
    ##              CLIP2              PDS5A             PIK3CD              ROCK1 
    ##              0.196              0.196              0.196              0.196 
    ##              SALL1            ATXN7L3               FGD1             LINGO1 
    ##              0.196              0.196              0.196              0.196 
    ##            TNFAIP3             NUFIP2              RIC8B           IL1RAPL1 
    ##              0.196              0.196              0.196              0.197 
    ##             POU2F2              COPS4            DNAJC13               XRN1 
    ##              0.197              0.197              0.197              0.197 
    ##            ANKRD12            PACSIN1            RAPGEF2               ETS2 
    ##              0.197              0.197              0.197              0.197 
    ##              HCFC2              PAIP1              LAMC1              SSBP3 
    ##              0.198              0.198              0.198              0.198 
    ##             NEDD4L           KIAA0368              PDGFC           KIAA2022 
    ##              0.198              0.198              0.198              0.198 
    ##              IPO13              PSMA5                YY1              ASXL3 
    ##              0.198              0.198              0.198              0.199 
    ##              PA2G4             DLGAP2             MPPED2              CPSF7 
    ##              0.199              0.199              0.199              0.199 
    ##               VCAN               PAX8              PRKCE             KCTD16 
    ##              0.199              0.200              0.200              0.200 
    ##              KIF4A             DAB2IP            SIPA1L3                TOX 
    ##              0.200              0.200              0.200              0.200 
    ##              NTNG2             CLASP2              HERC2              DCAF8 
    ##              0.201              0.201              0.201              0.201 
    ##              STIP1              ACIN1              HSPA8               RPL6 
    ##              0.201              0.201              0.201              0.201 
    ##             SETDB1               NFIA             RNF220               JPH3 
    ##              0.201              0.202              0.202              0.202 
    ##            CACNA1I              RFWD2             STAT5A            GATAD2A 
    ##              0.202              0.202              0.202              0.202 
    ##             COL1A2              DMXL1              AP1G1             OSBPL6 
    ##              0.202              0.202              0.202              0.202 
    ##              USP19            SMARCA2           ARHGAP21            CAMSAP2 
    ##              0.203              0.203              0.203              0.203 
    ##               PHF1              PRKCG           KIAA1244              PDZD8 
    ##              0.203              0.203              0.203              0.203 
    ##             BTBD11              UNC5A             CTDSP1              CELF5 
    ##              0.204              0.204              0.204              0.204 
    ##              RBM14           TRAPPC10               CYLD              FTSJ1 
    ##              0.204              0.204              0.204              0.204 
    ##              ATP9A              BRSK1              NFAT5               ATF7 
    ##              0.205              0.205              0.205              0.205 
    ##              IL21R              RBM39             EIF4A2              MORC2 
    ##              0.205              0.205              0.205              0.205 
    ##               FBN2            ZC3H12B            CSNK1G3              PDE4D 
    ##              0.205              0.206              0.206              0.206 
    ##               GNAL            PITPNM3              UBE3A            ZNF512B 
    ##              0.206              0.206              0.206              0.206 
    ##            BHLHE41              SETD2            SMARCD1              RBM27 
    ##              0.206              0.206              0.206              0.206 
    ##             SLC9A6             TOMM40              RIMS2             ANKIB1 
    ##              0.206              0.206              0.206              0.207 
    ##               CDK1             FBXO41              MAPK8               GLG1 
    ##              0.207              0.207              0.207              0.207 
    ##              PTPRC            RASGRF1               ACHE             HMGCS1 
    ##              0.207              0.207              0.207              0.207 
    ##              QSER1              ATAD5              SF3B4               DDI2 
    ##              0.208              0.208              0.208              0.208 
    ##              EP400              PAPD5             GCN1L1             CORO1C 
    ##              0.208              0.208              0.208              0.209 
    ##              IL2RG               MBD6           ARHGAP30              PPP5C 
    ##              0.209              0.209              0.209              0.209 
    ##               LMNA             RAB11A            SYNCRIP              SOX10 
    ##              0.209              0.209              0.209              0.209 
    ##               BRAF               SVOP             HNRNPR               RPGR 
    ##              0.209              0.209              0.210              0.210 
    ##               WNK2             MAPRE1              WDR45              MEPCE 
    ##              0.210              0.210              0.210              0.210 
    ##             BRINP2           C11orf84              NPAS3              USP15 
    ##              0.210              0.210              0.211              0.211 
    ##             ZC3H13              DGCR8            CACNA1G               ING3 
    ##              0.211              0.211              0.211              0.211 
    ##              CSMD1             ERCC6L             DNAJA2               TAF5 
    ##              0.211              0.211              0.211              0.211 
    ##           C10orf12             KIF26B       RP11-159G9.5               UBN1 
    ##              0.211              0.211              0.211              0.211 
    ##              MTSS1              ASTN1              NCOA5              OTUD4 
    ##              0.212              0.212              0.212              0.212 
    ##              SORT1              UCHL1               LEF1                FRY 
    ##              0.212              0.212              0.212              0.212 
    ##             MAP3K7                LCK             ZNF445              IGF2R 
    ##              0.212              0.212              0.212              0.212 
    ##             BTBD18              CHRM1              ASH2L              ACTN1 
    ##              0.212              0.213              0.213              0.213 
    ##             ATF7IP            ATP6AP1           HSP90AB1             PDE10A 
    ##              0.213              0.213              0.213              0.213 
    ##               RFX2              COPB1             RNF111           TBC1D10B 
    ##              0.213              0.213              0.213              0.213 
    ##              WDR33               NRG2               UBN2               UBR2 
    ##              0.213              0.213              0.213              0.213 
    ##            FAM193A             SLC6A8               SARS              MLLT4 
    ##              0.213              0.213              0.213              0.213 
    ##             SRCIN1            CACNA1D                ST5             MAP3K1 
    ##              0.213              0.214              0.214              0.214 
    ##              IGSF1               COG3               JAK1              RAB14 
    ##              0.214              0.214              0.214              0.214 
    ##              EVI5L             DYRK1A             NEURL4               DAB1 
    ##              0.214              0.214              0.214              0.215 
    ##              KCNH4           ARHGEF11              CCAR1             THSD7A 
    ##              0.215              0.215              0.215              0.215 
    ##               SOX6              FGFR1              TTYH3              ROCK2 
    ##              0.215              0.215              0.215              0.215 
    ##             ANKS1B              NFASC              WASF1              SENP6 
    ##              0.215              0.215              0.215              0.215 
    ##            CSNK1A1              WDR48               GCLM               RPS8 
    ##              0.215              0.215              0.216              0.216 
    ##               FAF2             HS6ST3               RBM4            CSNK1G1 
    ##              0.216              0.216              0.216              0.216 
    ##              ACACA            COL11A1              U2AF1              ATP7A 
    ##              0.216              0.216              0.216              0.216 
    ##               TPX2              NUMA1              NR2F2               TCF4 
    ##              0.217              0.217              0.217              0.217 
    ##            HTATSF1              TRAF6              RSRC2             SNRPA1 
    ##              0.217              0.217              0.217              0.217 
    ##               SNW1                NXN               CCT4             SAMD4A 
    ##              0.217              0.218              0.218              0.218 
    ##              RGAG1              ZMIZ2             ANP32A              GPM6A 
    ##              0.218              0.218              0.218              0.218 
    ##               ITCH              XRCC6               AGO4             DDX39B 
    ##              0.218              0.218              0.218              0.218 
    ##             ZNF536          HNRNPA2B1              EWSR1              NUMBL 
    ##              0.218              0.218              0.218              0.218 
    ##              STRBP              CPEB1            SUV39H1              FOXP2 
    ##              0.218              0.218              0.218              0.219 
    ##               RXRB                SP2                AQR             SEMA3F 
    ##              0.219              0.219              0.219              0.219 
    ##               SMG7              SRP54              KDM3A               FRS2 
    ##              0.219              0.219              0.219              0.219 
    ##               BRD2            ADAMTS6               GIT1            ATP13A3 
    ##              0.219              0.220              0.220              0.220 
    ##              RAB2A              ASAP2            PRKAR2B            SMARCB1 
    ##              0.220              0.220              0.220              0.220 
    ##             STAT5B              HAUS7              SYCP1             YTHDF3 
    ##              0.220              0.220              0.220              0.220 
    ##           ARHGAP44             MAP2K4               WDR1              NR2C2 
    ##              0.220              0.220              0.220              0.220 
    ##           ARHGAP23             SEMA6A               BCL3              MGAT1 
    ##              0.220              0.221              0.221              0.221 
    ##              SOCS6               AFF3      RP11-1055B8.7             FNDC3A 
    ##              0.221              0.221              0.221              0.221 
    ##             DPYSL3           C17orf85              ITSN1               NOL4 
    ##              0.221              0.221              0.221              0.221 
    ##              PROX1              CAMKV             ZBTB38             UBQLN1 
    ##              0.221              0.221              0.221              0.221 
    ##            ATXN7L1                LYN              SMAD4              BRWD1 
    ##              0.221              0.221              0.222              0.222 
    ##             ELAVL1              SF3B2              CLCN3              SLIT1 
    ##              0.222              0.222              0.222              0.222 
    ##             TUBB2B              PHF19            FAM135B              FCHO2 
    ##              0.222              0.222              0.222              0.222 
    ##              PRDM2             DNAJC6             ATP2B3            RAPGEF1 
    ##              0.222              0.222              0.222              0.222 
    ##               NTN1              FEZF2             FBXO32              PTBP2 
    ##              0.222              0.222              0.222              0.222 
    ##              EFNB2              CDK17              MAML3               THRB 
    ##              0.222              0.222              0.223              0.223 
    ##              MORC3              RAPH1              CREB1              IREB2 
    ##              0.223              0.223              0.223              0.223 
    ##               GDF6              OLFM1           SLC25A14            SLC16A2 
    ##              0.223              0.223              0.223              0.223 
    ##              TRAF3              MAST3              SPAST              USP32 
    ##              0.224              0.224              0.224              0.224 
    ##               LRP8               GJC1             GABBR1             PPP1R8 
    ##              0.224              0.224              0.224              0.224 
    ##               BMP2               P4HB                GDA             IQSEC1 
    ##              0.224              0.224              0.224              0.224 
    ##              KPNA6               LHX2               NFIB                WIZ 
    ##              0.224              0.224              0.225              0.225 
    ##             ZBTB43             BMPR1B              TTC17              PSMA7 
    ##              0.225              0.225              0.225              0.225 
    ##               E2F1              LIMK1               APC2               CSF1 
    ##              0.225              0.225              0.225              0.226 
    ##               VAV1            IGF2BP1              UFD1L              CDKL5 
    ##              0.226              0.226              0.226              0.226 
    ##              MIER3               CCT2            SHROOM4              TBX21 
    ##              0.226              0.226              0.226              0.226 
    ##           HNRNPUL2             CARD11               PCNX             PCNXL3 
    ##              0.227              0.227              0.227              0.227 
    ##             SCHIP1             RAD23B              PRPF4             COL6A1 
    ##              0.227              0.227              0.227              0.227 
    ##              DCLK1              STK40              FANCB              SETD5 
    ##              0.227              0.227              0.227              0.227 
    ##             HNRNPD             TCERG1               GNL1              CPEB3 
    ##              0.227              0.227              0.227              0.227 
    ##              PYGO1               CUL3              LZTS3              KIF5A 
    ##              0.228              0.228              0.228              0.228 
    ##             INO80D              SAFB2              NAA30              PCBP2 
    ##              0.228              0.228              0.228              0.228 
    ##             ZNF384              LATS1              MAML2             CAMK1D 
    ##              0.228              0.228              0.228              0.228 
    ##             TRIM46             PIK3R1              IL6ST            ARHGAP6 
    ##              0.228              0.228              0.228              0.228 
    ##                REL             MAGED1             GABRA2            JAKMIP1 
    ##              0.229              0.229              0.229              0.229 
    ##             PLXNA4           CACNA2D2             SNAP25              CCND2 
    ##              0.229              0.229              0.229              0.229 
    ##            FAM168A              PATZ1            DENND1A              KDM7A 
    ##              0.229              0.229              0.229              0.229 
    ##              SOX30              WDTC1              THOC1              FMNL1 
    ##              0.230              0.230              0.230              0.230 
    ##              MRTO4              STK39               IPO9            HERPUD1 
    ##              0.230              0.230              0.230              0.230 
    ##              REEP1               GPS1              PLRG1              UBE2K 
    ##              0.230              0.231              0.231              0.231 
    ##              PSMA6             CNOT6L              AGAP1             CSNK1D 
    ##              0.231              0.231              0.231              0.231 
    ##              CELF4              RDH10             MAP2K7              GABPA 
    ##              0.231              0.231              0.232              0.232 
    ##             KIF13A            FAM169A               ACTB             NUP205 
    ##              0.232              0.232              0.232              0.232 
    ##              GNL3L             PABPC4           RALGAPA1              FADS2 
    ##              0.232              0.232              0.232              0.232 
    ##              FBXW7           TMEM132D             ATP11C            CASKIN2 
    ##              0.232              0.232              0.232              0.232 
    ##             NUDT21              BACH2              PTPRM             HMBOX1 
    ##              0.232              0.233              0.233              0.233 
    ##              SOCS7            HNRNPA3             MAGED2             PTGES3 
    ##              0.233              0.233              0.233              0.233 
    ##             SRSF11             LSM14B              NCSTN              AIFM1 
    ##              0.233              0.234              0.234              0.234 
    ##             ZNF704              TTBK1               PBX2              YIPF5 
    ##              0.234              0.234              0.234              0.234 
    ##             KBTBD2            SLC30A1           KIAA1468             ZNF827 
    ##              0.234              0.234              0.234              0.234 
    ##             PICALM              PRKG1              YWHAE              CRTC2 
    ##              0.234              0.235              0.235              0.235 
    ##              PSMD6               SMC2              THBS1             MAN1A2 
    ##              0.235              0.235              0.235              0.235 
    ##              GMEB2              SCFD1           TNFRSF1A              SSBP2 
    ##              0.236              0.236              0.236              0.236 
    ##               CDYL               REV1              POLD3              ESRP1 
    ##              0.236              0.236              0.236              0.236 
    ##                FUS               RARG           EPB41L4B               NDC1 
    ##              0.237              0.237              0.237              0.237 
    ##               MSI2         RBM14-RBM4              ZMAT2              AZIN1 
    ##              0.237              0.237              0.237              0.237 
    ##              CRIM1             PPP1R7              RNF40              MTMR1 
    ##              0.238              0.238              0.238              0.238 
    ##             PRPF39                EVL             HEXIM1               DLG4 
    ##              0.238              0.238              0.238              0.238 
    ##             KANSL1             AMIGO1               DAZL               RPL3 
    ##              0.238              0.238              0.239              0.239 
    ##            CAMSAP1              MDGA2               TLE1              LMTK2 
    ##              0.239              0.239              0.239              0.239 
    ##               HCN4              CHMP6              ACTN2              HDAC7 
    ##              0.239              0.239              0.239              0.239 
    ##              HELLS                ENG              SLMAP            C6orf62 
    ##              0.239              0.240              0.240              0.240 
    ##            CACNA1B            KHDRBS1              TRIB2              VEZF1 
    ##              0.240              0.240              0.240              0.240 
    ##              HERC1              CCND3              GSPT2              TXLNA 
    ##              0.240              0.240              0.240              0.240 
    ##              MED15              FBXL5              PCSK2              FOSL2 
    ##              0.241              0.241              0.241              0.241 
    ##             CAMK2G            SH3GLB1             SLC2A1              TAOK2 
    ##              0.241              0.241              0.241              0.241 
    ##              TAOK3              TSHZ1             FBXO33             YTHDF1 
    ##              0.241              0.241              0.241              0.241 
    ##              YWHAB            TMEM259               UBA2               RBPJ 
    ##              0.241              0.241              0.241              0.241 
    ##              ITPK1            RPS6KA5                SHH             PDLIM4 
    ##              0.242              0.242              0.242              0.242 
    ##              SRSF1              FOXJ2              MEMO1              BCL9L 
    ##              0.242              0.242              0.242              0.242 
    ##             TOPORS              IGBP1              LRCH2            TOMM70A 
    ##              0.242              0.242              0.242              0.242 
    ##              NOVA2            SLC38A2              DCHS1             SRGAP3 
    ##              0.242              0.242              0.242              0.243 
    ##               PHC1               DDB1              PCDH1               SCAI 
    ##              0.243              0.243              0.243              0.243 
    ##             PPP4R1               PTK2             PLAGL2              MEF2D 
    ##              0.243              0.243              0.243              0.243 
    ##             RPL18A           TNFSF13B               PIM3               GMPS 
    ##              0.243              0.244              0.244              0.244 
    ##              GRID1             HNRNPL               NONO               ELL2 
    ##              0.244              0.244              0.244              0.244 
    ##            DNAJC11               GRM7            TMEM131              LONP1 
    ##              0.244              0.244              0.244              0.244 
    ##                SRC              HSPA4          KIAA1211L               SPI1 
    ##              0.244              0.244              0.244              0.244 
    ##             TRIM27              CEPT1              RPL19            EPS15L1 
    ##              0.244              0.244              0.244              0.245 
    ##              KPNA1           CACNA2D1              STK11              GABRD 
    ##              0.245              0.245              0.245              0.245 
    ##              EIF5B              UHRF2             CELSR3               BRD1 
    ##              0.245              0.245              0.245              0.245 
    ##               KLF7          RAB11FIP3              RABL6            GRAMD1B 
    ##              0.245              0.245              0.245              0.245 
    ##              PUF60                FBL               TBX2             ZRANB1 
    ##              0.245              0.246              0.246              0.246 
    ##             ATP2A2               FLT1             ACTL6A              HERC4 
    ##              0.246              0.246              0.246              0.246 
    ##              AKAP4              JADE1              TRA2B              BCAR1 
    ##              0.246              0.246              0.246              0.246 
    ##               ADD3              PSMC6             POLR1B            PRPF40B 
    ##              0.246              0.246              0.246              0.246 
    ##              ASTN2              GRIA2              LNPEP               RLIM 
    ##              0.246              0.247              0.247              0.247 
    ##             RASIP1               NAV2               ILF2               NUS1 
    ##              0.247              0.247              0.247              0.247 
    ##               NPTN                WT1              CLOCK              KCNH3 
    ##              0.247              0.247              0.247              0.247 
    ##              TOP2B              DCTN2              NCOA1               RHEB 
    ##              0.247              0.247              0.247              0.247 
    ##           ANKRD34A              HDAC9           ANKRD13C               WWP1 
    ##              0.247              0.248              0.248              0.248 
    ##             ZNF319              WDR18               CBX8              LENG8 
    ##              0.248              0.248              0.248              0.248 
    ##               KAL1             LRRC8B             ZBTB21              MORC4 
    ##              0.248              0.249              0.249              0.249 
    ##           ARHGAP29              NELL2              PDE8B             ZNF281 
    ##              0.249              0.249              0.249              0.249 
    ##            ARFGEF2              PHF13             INPP4A              MUC5B 
    ##              0.249              0.249              0.249              0.249 
    ##                KDR               BAI2             NLGN4X               LRP2 
    ##              0.249              0.249              0.249              0.249 
    ##              CHRM3            PITPNC1             ZNF496               USF2 
    ##              0.249              0.249              0.249              0.249 
    ##              PDIA3               FLNC              MTMR4             TRIP13 
    ##              0.250              0.250              0.250              0.250 
    ##               ACE2             KIRREL             CELSR1                HGF 
    ##              0.250              0.250              0.250              0.250 
    ##               RYR2             TM9SF4           FAM160B1               SYN1 
    ##              0.250              0.250              0.250              0.251 
    ##                CAD               FAF1              FNIP1               ETF1 
    ##              0.251              0.251              0.251              0.251 
    ##              INHBA           SLC30A10              STT3B             MRC1L1 
    ##              0.251              0.251              0.251              0.251 
    ##              MEX3A             CAMK2A              INTS3              RNPS1 
    ##              0.251              0.251              0.251              0.251 
    ##             NFE2L1             SLC6A3             CELSR2             TFAP2D 
    ##              0.251              0.251              0.251              0.252 
    ##             PIK3CB            ZC3H12C             DOPEY1              ZMIZ1 
    ##              0.252              0.252              0.252              0.252 
    ##             ZYG11B             SKIDA1              SEPT9             SCUBE3 
    ##              0.252              0.252              0.252              0.252 
    ##               KLF6             SFMBT2               DNM1              ZNFX1 
    ##              0.252              0.252              0.252              0.252 
    ##              RBBP4              LEMD3              RUSC2              AKAP6 
    ##              0.253              0.253              0.253              0.253 
    ##           TMEM132B             RBM15B             ZNF644              AUTS2 
    ##              0.253              0.253              0.253              0.253 
    ##               SIM1               FAT3            TP53BP1              HECW1 
    ##              0.253              0.253              0.253              0.253 
    ##              NEGR1               XPO5              SVEP1              SCAF1 
    ##              0.253              0.253              0.254              0.254 
    ##             DLGAP1               DLX2              TTYH1            CSNK1G2 
    ##              0.254              0.254              0.254              0.254 
    ##              RAB35              PIAS2              ADCY5               GRM5 
    ##              0.254              0.254              0.254              0.254 
    ##             BOD1L1              PTPRS              EDNRA              NRXN1 
    ##              0.254              0.254              0.254              0.254 
    ##           ARHGAP32              FRMD5              VPS54             UBE2R2 
    ##              0.254              0.254              0.254              0.255 
    ##           AMMECR1L              CCT6A              TCFL5               PBX1 
    ##              0.255              0.255              0.255              0.255 
    ##               TCP1               IRF1             PRKAG2              LRP12 
    ##              0.255              0.256              0.256              0.256 
    ##               EBF1            FAM131B               SOS2              RAD21 
    ##              0.256              0.256              0.256              0.256 
    ##             MOSPD2             STXBP5              KCNJ3             POLR1A 
    ##              0.256              0.256              0.256              0.256 
    ##               TP73              CWC15             NDFIP1              EPHB1 
    ##              0.256              0.256              0.257              0.257 
    ##             PHLPP1              MAGI2            PPP2R2A             ZNF618 
    ##              0.257              0.257              0.257              0.257 
    ##              HECW2               USP1            PHACTR3             FBXL20 
    ##              0.257              0.257              0.257              0.257 
    ##               BRD7              MTMR3              CDK16             HMGXB3 
    ##              0.258              0.258              0.258              0.258 
    ##              MARK4            FAM208A              CALD1             PPP2R4 
    ##              0.258              0.258              0.258              0.258 
    ##               DLC1               ELF2               LIN9               NRG1 
    ##              0.258              0.258              0.258              0.258 
    ##              NRXN2              GRIA1              STMN2               TSKS 
    ##              0.258              0.259              0.259              0.259 
    ##            IRF2BP1             TMEM57              NTRK3             TFAP2B 
    ##              0.259              0.259              0.259              0.259 
    ##      DKFZP761J1410               DEF6               MKL1              SMAD2 
    ##              0.259              0.259              0.259              0.259 
    ##               MBD3               CBX1              HMGCR             PLAGL1 
    ##              0.260              0.260              0.260              0.260 
    ##              TRIM8             SLC8A2                DSP            PPP1R37 
    ##              0.260              0.260              0.260              0.260 
    ##               OSBP              PITX2              ADCY2              DAAM1 
    ##              0.260              0.260              0.261              0.261 
    ##              DOCK2              TTBK2              DUSP6               RPS6 
    ##              0.261              0.261              0.261              0.261 
    ##               EXT1              PHF20              RPS18            STXBP5L 
    ##              0.261              0.261              0.261              0.261 
    ##             TFAP2A               PCNA              YLPM1                ERF 
    ##              0.261              0.261              0.261              0.261 
    ##                CIT              GNAO1              MKRN1              NR5A1 
    ##              0.261              0.261              0.261              0.262 
    ##              ZBTB1              KCNB2              KCNC1              GRIK2 
    ##              0.262              0.262              0.262              0.262 
    ##             PKNOX1               RELB              ARRB1             CAMTA2 
    ##              0.262              0.262              0.262              0.262 
    ##              RAB8A                RAN              ATG4A             PLXNA1 
    ##              0.262              0.262              0.262              0.262 
    ##               CBX6               DLG5               PHF8               DDX4 
    ##              0.262              0.262              0.263              0.263 
    ##             OPN1LW             PPP2CA              PTGS2               TPP2 
    ##              0.263              0.263              0.263              0.263 
    ##               HUNK              AP1B1              AP2B1            CBFA2T2 
    ##              0.263              0.263              0.263              0.263 
    ##            ARHGAP5               DBN1               MSI1              TENM3 
    ##              0.263              0.263              0.264              0.264 
    ##              FCHO1              HSPH1                ST7              TEAD3 
    ##              0.264              0.264              0.264              0.264 
    ##              ITGB8              SMAP1              CCNL1             SLC9A7 
    ##              0.264              0.264              0.264              0.264 
    ##               EYA3               VRTN            CTDNEP1              RRAS2 
    ##              0.264              0.264              0.264              0.264 
    ##             NFATC3             FBXO38            MORF4L1              CHRM4 
    ##              0.264              0.265              0.265              0.265 
    ##             HECTD2               SBF1            SLC44A1               NKD1 
    ##              0.265              0.265              0.265              0.265 
    ##               GNB2            FAM135A              ARPC4             NMNAT2 
    ##              0.265              0.265              0.265              0.265 
    ##               NMT1              ERBB4             SESTD1             ATRNL1 
    ##              0.265              0.265              0.265              0.266 
    ##             PAPOLG              DIP2B               GIT2              NCOA6 
    ##              0.266              0.266              0.266              0.266 
    ##               AMOT              PSIP1              PDHA1               CBLB 
    ##              0.266              0.266              0.266              0.266 
    ##                EDA             SLC4A8            HDGFRP3             RABEP1 
    ##              0.266              0.266              0.266              0.266 
    ##               IWS1              AGAP2               TP63             SLC7A1 
    ##              0.267              0.267              0.267              0.267 
    ##              FOXO3               ZIC2               ROR1              SMAD1 
    ##              0.267              0.267              0.267              0.267 
    ##           PPARGC1A              MMP14             SRGAP1              CNNM2 
    ##              0.267              0.267              0.267              0.267 
    ##              MEGF9              GANAB               RBL2            RPS6KA4 
    ##              0.267              0.267              0.268              0.268 
    ##              LRFN3              APEX2               MBD1              HSPD1 
    ##              0.268              0.268              0.268              0.268 
    ##               MPP1              STX1B               ACRC               SWT1 
    ##              0.268              0.268              0.268              0.268 
    ##             LPCAT3           ATP6V0A1           ATP6V0D1                CPE 
    ##              0.268              0.268              0.268              0.269 
    ##               DPP8           PPP1R13B               RALA              TERF2 
    ##              0.269              0.269              0.269              0.269 
    ##               SAE1              INHBB              WNT5A             TRIM39 
    ##              0.269              0.269              0.269              0.269 
    ##            L3MBTL3            GPRASP2              LSAMP              SRPK2 
    ##              0.269              0.269              0.269              0.269 
    ##             PABPC1              FOXJ3              NABP2              GRHL2 
    ##              0.270              0.270              0.270              0.270 
    ##              GTF2I              FGFR2            HNRNPH3              TRIM2 
    ##              0.270              0.270              0.270              0.270 
    ##              TGFB3             SLC8A1             APPBP2            HSP90B1 
    ##              0.270              0.271              0.271              0.271 
    ##              RAB10             FBXO28              SUGP1             CHAMP1 
    ##              0.271              0.271              0.271              0.271 
    ##               ODF1              LSM12              MCMBP              MYO16 
    ##              0.271              0.271              0.271              0.272 
    ##              ADCY1             DCAF15               FRYL              RPL15 
    ##              0.272              0.272              0.272              0.272 
    ##              COPS6             PKNOX2              RBPMS              RPS10 
    ##              0.272              0.272              0.272              0.273 
    ##             TCF7L2              ATAD2              SPRY2             PLXNC1 
    ##              0.273              0.273              0.273              0.273 
    ##             SLC4A4                OMG               GMDS              TAF7L 
    ##              0.273              0.273              0.273              0.273 
    ##               BANP              CDC27               RFX1             DPYSL2 
    ##              0.273              0.273              0.273              0.274 
    ##              SPIN1             PPP1CA              TAOK1             FBXO42 
    ##              0.274              0.274              0.274              0.274 
    ##             FBXO45              KCNH2              LMX1A               OGDH 
    ##              0.274              0.274              0.274              0.274 
    ##               SLTM               FMN2               ADD2              KCNN2 
    ##              0.274              0.274              0.274              0.274 
    ##               VWC2              ERP44             IQGAP1             QRICH1 
    ##              0.274              0.274              0.275              0.275 
    ##              SSRP1              EIF3G              RQCD1            FAM193B 
    ##              0.275              0.275              0.275              0.275 
    ##              DAPK1               AFF2             TARDBP             ADAM23 
    ##              0.275              0.275              0.275              0.275 
    ##             MTMR12              NR0B1               PTK7             ZNF770 
    ##              0.275              0.275              0.275              0.275 
    ##               FHL1              RPS3A              PSMB1               NASP 
    ##              0.276              0.276              0.276              0.276 
    ##             ARID3A              RBM15              GNA12           ARHGAP26 
    ##              0.276              0.276              0.276              0.276 
    ##            DCUN1D1              KDM1A             KCNMA1               XPOT 
    ##              0.276              0.276              0.276              0.276 
    ##              SETD8              INTS2         AC005358.1              GRHL3 
    ##              0.276              0.277              0.277              0.277 
    ##             PRR14L               BAP1              PSPC1              ARNTL 
    ##              0.277              0.277              0.277              0.277 
    ##              ZFHX2             GTF2A1            CEP170B               HAS2 
    ##              0.277              0.277              0.277              0.277 
    ##              PDZD2           PRICKLE2              RAB7A               SIK2 
    ##              0.277              0.277              0.277              0.277 
    ##            FAM196B              TESK1              CTBP1              AGFG1 
    ##              0.278              0.278              0.278              0.278 
    ##              SART3              PSMA3           SLC39A10            ATP6V0B 
    ##              0.278              0.278              0.278              0.278 
    ##               SSH2              FOXJ1             ZNF710                FYN 
    ##              0.278              0.278              0.278              0.278 
    ##                CNP              PACS2            RAPGEF6              CLCN5 
    ##              0.279              0.279              0.279              0.279 
    ##            PIK3C2B               ARNT               TLN2              NRIP1 
    ##              0.279              0.279              0.279              0.279 
    ##               RYR3              KCND3              LRCH1              PTPN4 
    ##              0.279              0.279              0.279              0.279 
    ##             ZNF182               ASB7             YEATS2              STRN4 
    ##              0.279              0.280              0.280              0.280 
    ##               CUL9              NAA15              DCAF8               EML1 
    ##              0.280              0.280              0.280              0.280 
    ##              RBBP7               FLI1           TMEM185B             ADAM17 
    ##              0.280              0.280              0.280              0.280 
    ##            DENND4A           IVNS1ABP               FIGN             MICAL3 
    ##              0.280              0.281              0.281              0.281 
    ##              CRMP1             AHCYL1              KCNA4               TBX3 
    ##              0.281              0.281              0.281              0.281 
    ##             SLC1A3               EHD1              PLCG2               CDK6 
    ##              0.281              0.281              0.281              0.282 
    ##            COL12A1              USP10             BCL11B               DBF4 
    ##              0.282              0.282              0.282              0.282 
    ##               MFN2              RIBC1               PLS3              ARIH2 
    ##              0.282              0.282              0.283              0.283 
    ##              PYGO2                DCC           CDC42BPA              LASP1 
    ##              0.283              0.283              0.283              0.283 
    ##             DIAPH2                SP4           KIAA1429               FZD5 
    ##              0.283              0.283              0.283              0.283 
    ##                SHE               RND3               ATF2               DDX1 
    ##              0.283              0.283              0.283              0.284 
    ##               LEPR               CLUH             PPP6R3               ST18 
    ##              0.284              0.284              0.284              0.284 
    ##              RPL18             MAP4K5               RAC2              PSMA1 
    ##              0.284              0.284              0.284              0.284 
    ##              PSMA2             PTP4A2             PTCHD1               AMD1 
    ##              0.284              0.284              0.284              0.284 
    ##        AF196779.12             BRINP1               PUM2            ERBB2IP 
    ##              0.284              0.284              0.284              0.285 
    ##             PLXNB2              TSHZ3               PSAP                SP8 
    ##              0.285              0.285              0.285              0.285 
    ##               CBX4           C9orf172               CHUK               GRM4 
    ##              0.285              0.285              0.285              0.285 
    ##               DPF1             GTF3C2          RAB11FIP4               E2F7 
    ##              0.285              0.285              0.285              0.285 
    ##              NISCH             ZCCHC2               EIF5              RIMS1 
    ##              0.285              0.285              0.285              0.285 
    ##              SENP1              GRIK5              IRAK1              KIF5B 
    ##              0.285              0.285              0.285              0.285 
    ##            ZFP36L2            GALNT13                LOX               PAK3 
    ##              0.285              0.285              0.285              0.286 
    ##               IRS2               PHF3             PTPRZ1              IKZF2 
    ##              0.286              0.286              0.286              0.286 
    ##               REST             UBE2D2             ARRDC3               DLG1 
    ##              0.286              0.286              0.286              0.286 
    ##            TMEM201              CTPS1              MEF2A             MYO18A 
    ##              0.286              0.287              0.287              0.287 
    ##              DDX3Y              KCNJ4             DLGAP4              KIF3C 
    ##              0.287              0.287              0.287              0.287 
    ##               SYT2              UBE2I             TNRC18            SIPA1L2 
    ##              0.287              0.287              0.287              0.287 
    ##              RC3H1             PITPNB           ATP6V1B2            RAPGEF4 
    ##              0.287              0.287              0.287              0.288 
    ##              KLF12              TOP2A              SFRP1               LHX9 
    ##              0.288              0.288              0.288              0.288 
    ##              FBXL3              BAZ2B            EPB41L2                 GK 
    ##              0.288              0.288              0.288              0.288 
    ##             FAM60A              FMNL2               ULK1              MLLT1 
    ##              0.288              0.288              0.288              0.288 
    ##              NSDHL                RYK             CACNG8              CDC23 
    ##              0.289              0.289              0.289              0.289 
    ##            CXorf22             HNRNPF           C16orf72            CNTNAP5 
    ##              0.289              0.289              0.289              0.289 
    ##              DHX37              YWHAG             ZBTB10              CMTR1 
    ##              0.289              0.289              0.289              0.289 
    ##                NRK                NF1               TFE3              PDPK1 
    ##              0.289              0.290              0.290              0.290 
    ##               CDH2            DENND4B              IKBKB             UBE2E2 
    ##              0.290              0.290              0.290              0.290 
    ##            CYP26B1             AHCYL2             MAP7D1                FST 
    ##              0.290              0.290              0.290              0.290 
    ##               MCL1              ARNT2               LDB1             DNAJC7 
    ##              0.290              0.291              0.291              0.291 
    ##              ZZEF1             SLC6A6              GMEB1            SLC40A1 
    ##              0.291              0.291              0.291              0.291 
    ##             GTF2H1           RAB3GAP2                 AR              SRSF3 
    ##              0.291              0.291              0.291              0.291 
    ##              RCOR3             ZNF331               GNAZ             IGFBP5 
    ##              0.292              0.292              0.292              0.292 
    ##             FAM50A              RBM47            CREB3L1                PGD 
    ##              0.292              0.292              0.292              0.292 
    ##              UPF3B             DNAJC2              BEND4              GATA2 
    ##              0.292              0.292              0.292              0.292 
    ##             SPATA2               TUBB            TSPAN14              RAI14 
    ##              0.292              0.292              0.292              0.292 
    ##              SCML1              PPM1E               DLX1              STIM2 
    ##              0.292              0.293              0.293              0.293 
    ##              PPEF1             EEF1A1              TRIM3              STX1A 
    ##              0.293              0.293              0.293              0.293 
    ##           ARHGEF17               ERN1           FAM160A2               RPL9 
    ##              0.293              0.293              0.293              0.293 
    ##               SPEG              FBXO3              SATB1            PPP2R2B 
    ##              0.293              0.293              0.293              0.293 
    ##               TERT               PRLR              KAT2B              MLXIP 
    ##              0.293              0.293              0.294              0.294 
    ##             AKAP13             PHF21B           IL1RAPL2               ERC2 
    ##              0.294              0.294              0.294              0.294 
    ##              ZNRF3            ATP13A1              IKZF3              UVRAG 
    ##              0.294              0.294              0.294              0.294 
    ##               OPA1              MYO5A              RBM20               CDH8 
    ##              0.294              0.295              0.295              0.295 
    ##           KIAA1549              LARGE             MAP4K1              RGS19 
    ##              0.295              0.295              0.295              0.295 
    ##               ATG5             GOLGA2             HNRNPC               PURG 
    ##              0.295              0.295              0.295              0.295 
    ##            RPS6KA6              KCNV1            SEC61A1            TSC22D4 
    ##              0.295              0.295              0.295              0.296 
    ##              STAT6             FBXO22               ELF3              NDEL1 
    ##              0.296              0.296              0.296              0.296 
    ##               PIGA           SH3PXD2A               PPAT              AP3D1 
    ##              0.296              0.296              0.296              0.296 
    ##                BCR              PLCB1              KCMF1              RPL11 
    ##              0.296              0.296              0.296              0.296 
    ##               SETX             ZNF532               ADD1              EFNA5 
    ##              0.296              0.296              0.296              0.296 
    ##             NUP188             SCAMP1            SLC38A5            RTN4RL1 
    ##              0.297              0.297              0.297              0.297 
    ##              SYNPO            CSNK2A1              SEC62            C18orf8 
    ##              0.297              0.297              0.297              0.297 
    ##            CYP46A1               RTN1              HDAC8              TBX22 
    ##              0.297              0.297              0.297              0.297 
    ##              LATS2             C2CD2L             PWWP2A              GDF11 
    ##              0.297              0.297              0.297              0.298 
    ##             PSMD14               HSF5             SEPHS1               DPP9 
    ##              0.298              0.298              0.298              0.298 
    ##            C7orf26             H2AFY2              SF3A1             CALCRL 
    ##              0.298              0.298              0.298              0.298 
    ##               ARF4               CNBP             ZBTB34              CSMD3 
    ##              0.298              0.298              0.298              0.299 
    ##              GON4L               SMG6              SRSF7               ETV3 
    ##              0.299              0.299              0.299              0.299 
    ##           KIAA1239               FAR1             SLC4A2              RTN4R 
    ##              0.299              0.299              0.299              0.300 
    ##               NRP1               FOSB              REXO1              RELL1 
    ##              0.300              0.300              0.300              0.300 
    ##             PTGER4             HP1BP3              MALT1              BAHD1 
    ##              0.300              0.300              0.300              0.300 
    ##            PPP2R5B             AKAP11              EOMES              RRAGC 
    ##              0.300              0.300              0.300              0.300 
    ##              PSMB5             ANTXR1                SRF              REPS1 
    ##              0.300              0.300              0.300              0.301 
    ##             TOPAZ1               SND1              ILKAP               DAG1 
    ##              0.301              0.301              0.301              0.301 
    ##              PIAS3              ZNF24               CCT7              CPEB4 
    ##              0.301              0.301              0.301              0.301 
    ##               SKIL              DOCK4                ELL             KLHL29 
    ##              0.301              0.301              0.301              0.301 
    ##       CTC-432M15.3               GNAQ               IRX5           KIAA0922 
    ##              0.301              0.302              0.302              0.302 
    ##               SNX2            HSPA12A             CYFIP1               TWF2 
    ##              0.302              0.302              0.302              0.302 
    ##              MPRIP              FOXF1             MAP3K3             LPGAT1 
    ##              0.302              0.302              0.302              0.302 
    ##           FAM171A1             PRPF31           DYNC1LI2            ABHD17C 
    ##              0.302              0.302              0.302              0.302 
    ##             CHMP4B             NCAPG2              EGLN1                RS1 
    ##              0.302              0.303              0.303              0.303 
    ##                 CS             UBQLN4              PCGF5              CSTF3 
    ##              0.303              0.303              0.303              0.303 
    ##               GGT7              SEL1L              PRMT8             PDZRN3 
    ##              0.303              0.303              0.303              0.303 
    ##               SUCO            PPP2R1A               GAB2              GLRA2 
    ##              0.303              0.303              0.304              0.304 
    ##               MID1              PSMB3               EEA1              MYO9A 
    ##              0.304              0.304              0.304              0.304 
    ##               MYT1             PPP6R1              SPRTN              TIAM1 
    ##              0.304              0.304              0.304              0.304 
    ##              WDFY4            GLTSCR1              PANK4               TLR7 
    ##              0.304              0.304              0.304              0.305 
    ##              MYH11            RAPGEF5            ANKRD50              BICD2 
    ##              0.305              0.305              0.305              0.305 
    ##              FGF17             ANAPC1               RFC1                MET 
    ##              0.305              0.305              0.305              0.305 
    ##             ZNF711              SENP2                KIT                SMS 
    ##              0.305              0.305              0.305              0.305 
    ##              CDC34                 C3             CNOT10             ATXN1L 
    ##              0.305              0.305              0.305              0.305 
    ##              NR6A1             KLHL10              VDAC1             LRRC4C 
    ##              0.305              0.306              0.306              0.306 
    ##              IGF1R              ITGB1            COL11A2             GABRB1 
    ##              0.306              0.306              0.306              0.306 
    ##              ADCY9              KCNQ5             INPP5D              FOXP4 
    ##              0.306              0.306              0.306              0.306 
    ##                RP2              ACSL4              PARD3             SS18L1 
    ##              0.306              0.306              0.306              0.306 
    ##               INF2              KCNK9             SEMA3A             AGPAT1 
    ##              0.307              0.307              0.307              0.307 
    ##             HOXA13              TSHZ2               RPL8              TDRD5 
    ##              0.307              0.307              0.307              0.307 
    ##               DRD1               LYST              UNC80               NAF1 
    ##              0.307              0.307              0.307              0.307 
    ##              PDE4A               MITF              CENPE             INPP5A 
    ##              0.307              0.308              0.308              0.308 
    ##              AHNAK               ATL1           C16orf70              SPDEF 
    ##              0.308              0.308              0.308              0.308 
    ##               NKTR             CLSTN1              SCN5A              GRIN1 
    ##              0.308              0.308              0.308              0.309 
    ##         HSPE1-MOB4               MCM6               SNPH              OTUB1 
    ##              0.309              0.309              0.309              0.309 
    ##             RNF145             MTSS1L               RPS7              NAMPT 
    ##              0.309              0.309              0.309              0.309 
    ##              HIF1A               XAB2              PTPRK            HNRNPH2 
    ##              0.309              0.309              0.309              0.309 
    ##            PCDHGC4              ATP5B              TBL1X             ZBTB17 
    ##              0.309              0.309              0.309              0.309 
    ##               GLI2           DCAF12L1             OSBPL9              CDC40 
    ##              0.309              0.309              0.309              0.310 
    ##             MAGEE1              NODAL              ABHD8             EIF2S1 
    ##              0.310              0.310              0.310              0.310 
    ##              EIF3F              NLGN3               PLK1            SLC12A2 
    ##              0.310              0.310              0.310              0.310 
    ##               NXF1                 XK             FAM49B             SLC9A3 
    ##              0.310              0.310              0.310              0.310 
    ##               TOB1            ADAMTS2               EPS8              TCOF1 
    ##              0.310              0.311              0.311              0.311 
    ##              LAMA5              MDGA1                NNT              FOXC1 
    ##              0.311              0.311              0.311              0.311 
    ##            CAPRIN1           KIAA0355             FBXL16            FAM155A 
    ##              0.311              0.311              0.311              0.311 
    ##              WASF2            ZCCHC14               ZXDB             LRRTM3 
    ##              0.311              0.311              0.312              0.312 
    ##              SCRIB             SRGAP2             ANKFY1             NUCKS1 
    ##              0.312              0.312              0.312              0.312 
    ##            RUNX1T1             MAGEL2             LRRC8D            COL27A1 
    ##              0.312              0.312              0.312              0.312 
    ##              PLCB3               FZD4               EGR3                SCD 
    ##              0.312              0.312              0.312              0.312 
    ##            CXorf56             SCUBE1                HLF              TBCEL 
    ##              0.313              0.313              0.313              0.313 
    ##             LIMCH1            PIP5K1C             CTNNA1              PTK2B 
    ##              0.313              0.313              0.313              0.313 
    ##             HIVEP3              VPS35             ZNF512            ARHGAP4 
    ##              0.313              0.313              0.313              0.313 
    ##              ABCE1              PTPRB             SEMA4D               DGKI 
    ##              0.314              0.314              0.314              0.314 
    ##             FBXO21              RNF41               MCF2             FBXW11 
    ##              0.314              0.314              0.314              0.314 
    ##              CASC5              INTS8               SQLE              CNIH2 
    ##              0.314              0.314              0.314              0.314 
    ##               HSF2               PAN2              TAF1L             TUBA1B 
    ##              0.315              0.315              0.315              0.315 
    ##             FAM65A              TUBB3              TBX20             PRDM11 
    ##              0.315              0.315              0.315              0.315 
    ##              HMGB2              MOB1A               KLC2               G2E3 
    ##              0.315              0.315              0.315              0.315 
    ##             CLSTN3             SPHKAP               AKT1             NKX2-3 
    ##              0.315              0.315              0.315              0.316 
    ##              YPEL2             KBTBD6              RPL13             CAPN15 
    ##              0.316              0.316              0.316              0.316 
    ##              ATXN1               EMX2              MGEA5              PATL1 
    ##              0.316              0.316              0.316              0.316 
    ##               BRD3              IP6K2           RASGEF1A               MDN1 
    ##              0.316              0.316              0.316              0.316 
    ##           RAP1GAP2               USP5             ZBTB20              ATG9A 
    ##              0.316              0.316              0.317              0.317 
    ##            RPS6KB1               ANO8             TUBA1A             VPS26B 
    ##              0.317              0.317              0.317              0.317 
    ##             RNF165              NUP62            SLC32A1               ZEB1 
    ##              0.317              0.317              0.317              0.317 
    ##              TAF12              LRP1B             TRIM41              DMTF1 
    ##              0.317              0.317              0.317              0.317 
    ##            KIRREL3               ETV6              PSEN1              TFDP1 
    ##              0.318              0.318              0.318              0.318 
    ##                JMY              NPAS4             PRKACA            RASGRP1 
    ##              0.318              0.318              0.318              0.318 
    ##              SLIT3               G6PD              STAU2               LRP6 
    ##              0.318              0.318              0.318              0.318 
    ##              KCNN3              HSPA9             ZBTB12              RBM41 
    ##              0.319              0.319              0.319              0.319 
    ##              STAT2             ELOVL6            PPP2R5C             SPATS2 
    ##              0.319              0.319              0.319              0.319 
    ##             ZNF687              USP12              NPAS1              SYNRG 
    ##              0.319              0.319              0.319              0.319 
    ##             SPRED1               MRC2            ATG16L1            SULT4A1 
    ##              0.319              0.319              0.319              0.319 
    ##              NR3C1               BRS3             DIAPH1               MON2 
    ##              0.320              0.320              0.320              0.320 
    ##              WNT3A              HIAT1             GPRIN1               RGL1 
    ##              0.320              0.320              0.320              0.320 
    ##              LRRN1              FBXO5               PAK2             BEGAIN 
    ##              0.320              0.320              0.320              0.320 
    ##              CCNE1           KIAA0907             CORO1A             NBEAL2 
    ##              0.320              0.320              0.320              0.320 
    ##               PPIG               ECE1            RANBP10              APBB1 
    ##              0.320              0.321              0.321              0.321 
    ##              CLCF1              GATA1               MTF1               PIM2 
    ##              0.321              0.321              0.321              0.321 
    ##             BCL11A             FAM46C             NOTCH3                IDS 
    ##              0.322              0.322              0.322              0.322 
    ##                SHB             COL4A6           SERPING1               RPS5 
    ##              0.322              0.322              0.322              0.322 
    ##             EIF2S3               OFD1             NUP155               SIX3 
    ##              0.322              0.322              0.323              0.323 
    ##               RBX1             DYRK1B              CCAR2             ZC3H7B 
    ##              0.323              0.323              0.323              0.323 
    ##               FGF9             PRKRIR              PVRL2               CLTA 
    ##              0.323              0.323              0.323              0.323 
    ##              SH2B1                GLS            ZKSCAN1             SEMA6D 
    ##              0.323              0.323              0.323              0.323 
    ##            IL13RA1             ZCCHC8             MAP3K8             BABAM1 
    ##              0.323              0.323              0.323              0.324 
    ##            TSC22D2              RAB3A             ZNF483                CFH 
    ##              0.324              0.324              0.324              0.324 
    ##             SEC16A             SFMBT1             SLC5A3             SLC5A3 
    ##              0.324              0.324              0.324              0.324 
    ##              CDK13             IGSF9B              ACAP2             MAP2K2 
    ##              0.324              0.324              0.324              0.325 
    ##             PDCD10             DIS3L2             PDGFRB           C19orf68 
    ##              0.325              0.325              0.325              0.325 
    ##             ZFYVE9                UTY            NCKAP1L               BNC1 
    ##              0.325              0.325              0.325              0.326 
    ##             LRRC25             ABLIM1               LHX6            SLC20A2 
    ##              0.326              0.326              0.326              0.326 
    ##             CLPTM1             TRIM67            CDC37L1              MAGT1 
    ##              0.326              0.326              0.326              0.326 
    ##             NKX2-5              DCLK2              CALM1              ELFN2 
    ##              0.326              0.326              0.327              0.327 
    ##            SLC41A2             ACTR1A             GEMIN8              LTBP4 
    ##              0.327              0.327              0.327              0.327 
    ##              ATG13               CCT8              ATXN7             PITPNA 
    ##              0.328              0.328              0.328              0.328 
    ##           DYNC1LI1            TSC22D1             R3HDM1               FJX1 
    ##              0.328              0.328              0.328              0.328 
    ##             FAM65B              RPS11            SKIV2L2              MAK16 
    ##              0.328              0.328              0.328              0.328 
    ##             HOXA11               DEDD             ZBTB7A              BRSK2 
    ##              0.329              0.329              0.329              0.329 
    ##           CDKN2AIP                ERG                CRK               SYT7 
    ##              0.329              0.329              0.329              0.329 
    ##            TMEM164              FOXG1           GTF2IRD1             CACUL1 
    ##              0.329              0.329              0.329              0.330 
    ##              RAB6B               RARA               BCL6               IPMK 
    ##              0.330              0.330              0.330              0.330 
    ##              HNF1A              UNC5D              PPM1G             POFUT1 
    ##              0.330              0.330              0.330              0.330 
    ##             ZNF296               PRCC               VMP1              USP25 
    ##              0.330              0.330              0.330              0.330 
    ##             EIF2B4              SYNJ1               HIP1           ADAMTS10 
    ##              0.330              0.330              0.330              0.330 
    ##               ANO1               CBFB               ARAF               ARF5 
    ##              0.330              0.330              0.330              0.330 
    ##               DSG1               HHIP              USP33             ANP32B 
    ##              0.331              0.331              0.331              0.331 
    ##             DEPDC5              HYOU1              CLIP1              MYH14 
    ##              0.331              0.331              0.331              0.332 
    ##             RGS7BP            FAM208B               EYA1            ZCCHC24 
    ##              0.332              0.332              0.332              0.332 
    ##             VPS13D            SPECC1L              KPNA3              ZMYM2 
    ##              0.332              0.332              0.332              0.332 
    ##              DERL2              SNED1            NEUROD2             OTUD7A 
    ##              0.332              0.332              0.332              0.333 
    ##              PRDM1               SKP2              ATXN2              ELMO1 
    ##              0.333              0.333              0.333              0.333 
    ##            CCDC132               RCC2               BMP4              GSK3B 
    ##              0.333              0.334              0.334              0.334 
    ##              MAGI1             SMURF2              CTLA4              SNX18 
    ##              0.334              0.334              0.334              0.334 
    ##             SEPT11               MSH2                HK1               TOX4 
    ##              0.334              0.334              0.334              0.334 
    ##             CSRNP2                SP9              LTBP3            ALDH1A1 
    ##              0.334              0.335              0.335              0.335 
    ##             RASAL2              AGAP3             CCDC64              ILDR2 
    ##              0.335              0.335              0.335              0.335 
    ##              GRB10              NCOA3               PURA             INCENP 
    ##              0.335              0.335              0.335              0.336 
    ##                TNR           ARHGAP25             PPP1CC             RBFOX1 
    ##              0.336              0.336              0.336              0.336 
    ##              XYLT1           TMEM178B           KIAA2026                PSD 
    ##              0.336              0.336              0.336              0.337 
    ##              KLHL3            SPTY2D1           C1orf226             ANGPT1 
    ##              0.337              0.337              0.337              0.337 
    ##              ATG4B            SLITRK1           CACNA2D3              SART1 
    ##              0.337              0.337              0.337              0.337 
    ##              RBMX2             ANKS1A             KCNMB4              PSMC3 
    ##              0.337              0.338              0.338              0.338 
    ##               CDH6              FBLN1              GFPT1              STX12 
    ##              0.338              0.338              0.338              0.338 
    ##               DLG2            TMEM39A               SIX2              NUTF2 
    ##              0.338              0.338              0.338              0.338 
    ##             RNF126              DDHD1             SORCS1              TIPRL 
    ##              0.338              0.338              0.338              0.338 
    ##                FN1              PCGF3              PLCH1              PTPRU 
    ##              0.338              0.338              0.338              0.338 
    ##              FGF10             IFNGR2              RUSC1             ZNF318 
    ##              0.339              0.339              0.339              0.339 
    ##              SOX13             ZBTB37              BICC1              CNTFR 
    ##              0.339              0.339              0.339              0.339 
    ##            ZMYND19               LDB2               DDA1              LZTS1 
    ##              0.340              0.340              0.340              0.340 
    ##             PHOX2B            SLITRK3               BNC2              GPR64 
    ##              0.340              0.340              0.340              0.340 
    ##              HPRT1              RAB6A             MAP3K5             SLC4A7 
    ##              0.340              0.340              0.340              0.340 
    ##             GABRB3               JAK2               ELK1            SLC35F1 
    ##              0.341              0.341              0.341              0.341 
    ##               LCP2              CSTF1              KIF1C             RALGDS 
    ##              0.341              0.341              0.341              0.341 
    ##               LRP5             POU5F1               SIK1              XYLT2 
    ##              0.341              0.341              0.341              0.341 
    ##         AC006486.9             HOMER3             BCL2L1              DOCK1 
    ##              0.341              0.341              0.341              0.341 
    ##               GAB3             EIF2S2              SYT14            ARHGEF1 
    ##              0.341              0.342              0.342              0.342 
    ##             ZC3H14                EBP              FEZF1           HSD17B10 
    ##              0.342              0.342              0.342              0.342 
    ##              AP3B2               CLPX             CDC14B            EPB41L1 
    ##              0.342              0.342              0.343              0.343 
    ##        RPS10-NUDT3              MOV10               HMBS            HORMAD1 
    ##              0.343              0.343              0.343              0.343 
    ##             POU2F1              PRPS2           RAP1GDS1             PMEPA1 
    ##              0.343              0.343              0.343              0.343 
    ##              RUNX2             ALKBH5              PFDN6               SNX1 
    ##              0.343              0.343              0.343              0.344 
    ##              ASIC1              COPZ1              MTUS2               IRF8 
    ##              0.344              0.344              0.344              0.344 
    ##            PCDHAC2             POLR2B             CADPS2             LYPLA2 
    ##              0.344              0.344              0.344              0.344 
    ##              LRRC1            MAP3K13               OXR1              P4HA1 
    ##              0.344              0.344              0.344              0.344 
    ##               NTN4              AP3B1              KIF3A               FNTA 
    ##              0.344              0.344              0.344              0.345 
    ##              ITGA5              PI4KB              CXXC1               CD86 
    ##              0.345              0.345              0.345              0.345 
    ##               MAF1               VASP            HERPUD2              LARP4 
    ##              0.345              0.345              0.345              0.345 
    ##               LCOR           C6orf136              CCNT2              ATG2A 
    ##              0.345              0.346              0.346              0.346 
    ##              USP49              NUPL1               ASUN             ATP8A1 
    ##              0.346              0.346              0.346              0.346 
    ##             RNF214              EIF3I              SGMS1             TECPR2 
    ##              0.346              0.346              0.346              0.346 
    ##              VPS4A              DHX36             DROSHA               PER2 
    ##              0.347              0.347              0.347              0.347 
    ##               NEO1               AFF1               BZW1               ETV1 
    ##              0.347              0.347              0.347              0.347 
    ##              THBS2              BECN1              BEND3             LRRC8A 
    ##              0.347              0.347              0.347              0.347 
    ##             MEGF10            RNF144A              STAT4              TRAM2 
    ##              0.347              0.347              0.347              0.347 
    ##             STK38L              EXOC3                SLK              DCAF7 
    ##              0.348              0.348              0.348              0.348 
    ##              EHBP1              USP14               LHX5                RDX 
    ##              0.348              0.348              0.348              0.348 
    ##              PDE2A               E2F3              FOXN1               SGTA 
    ##              0.348              0.349              0.349              0.349 
    ##             ZNF503             RASSF8               CD22               FLT3 
    ##              0.349              0.349              0.349              0.349 
    ##               UGCG             ZNF146             FRMD4B               UBR1 
    ##              0.349              0.349              0.349              0.349 
    ##               YBX1            SLC25A5              RPS4X             NUP214 
    ##              0.349              0.349              0.349              0.349 
    ##            HNRNPLL             SLC4A1            RAP1GAP               MAP4 
    ##              0.349              0.349              0.349              0.349 
    ##              EFNB1               EGFR           RASGEF1B                EMD 
    ##              0.349              0.349              0.349              0.350 
    ##              RHOT1          KIDINS220                EN1              MMP16 
    ##              0.350              0.350              0.350              0.350 
    ##               CBX5              APPL1               LEO1            S100PBP 
    ##              0.350              0.350              0.350              0.350 
    ##             ARID5A              MRGBP             LHFPL3              NR3C2 
    ##              0.350              0.350              0.350              0.350 
    ##            SLC38A1              HOXA3             GTPBP1          RAB11FIP2 
    ##              0.350              0.350              0.351              0.351 
    ##              SREK1              MDGA2            DNAJB14            AGTPBP1 
    ##              0.351              0.351              0.351              0.351 
    ##                MAZ               SDK2               MRC1             PCDH10 
    ##              0.351              0.351              0.352              0.352 
    ##            SHROOM3               BMI1            KANSL1L              SYAP1 
    ##              0.352              0.352              0.352              0.352 
    ##              ASNA1              STRN3            KLHDC10                MYB 
    ##              0.352              0.352              0.353              0.353 
    ##             INPP4B            TBC1D30              TRAF4              UBFD1 
    ##              0.353              0.353              0.353              0.353 
    ##           HSP90AA1              ARAP1              PAPPA                TTN 
    ##              0.353              0.353              0.354              0.354 
    ##            TUBGCP3               MIA3            COLEC12              ARL8A 
    ##              0.354              0.354              0.354              0.354 
    ##             ZNF362              EIF4E              LCORL            ST3GAL2 
    ##              0.354              0.354              0.354              0.354 
    ##             RNF19B              DIP2A               RAF1             ELOVL5 
    ##              0.354              0.354              0.354              0.355 
    ##               XRN2            SLC24A3              JMJD6               PLP1 
    ##              0.355              0.355              0.355              0.355 
    ##              CADM3              SPCS2              UBE2H             BMPR1A 
    ##              0.355              0.355              0.355              0.355 
    ##              PSMC5               CD19             FLT3LG              GRIP1 
    ##              0.355              0.355              0.356              0.356 
    ##              WWTR1             ZNF326              KCNQ3              CALM3 
    ##              0.356              0.356              0.356              0.357 
    ##               GID8             SPTBN2               RCC1              EPHA6 
    ##              0.357              0.357              0.357              0.357 
    ##                FER              YWHAZ               GRB2               IPO8 
    ##              0.357              0.357              0.357              0.358 
    ##            SUV39H2               SVIL            OSBPL11              MIER2 
    ##              0.358              0.358              0.358              0.358 
    ##              UBE2Z             NFATC4                LCT              NTNG1 
    ##              0.358              0.358              0.358              0.358 
    ##             FIP1L1            ATP6V1A           TMEM194A             RNF103 
    ##              0.358              0.358              0.358              0.358 
    ##               WSB2               KTN1              FNBP1              PPTC7 
    ##              0.358              0.359              0.359              0.359 
    ##               CANX              KDM5D              CRTC3               SOX4 
    ##              0.359              0.359              0.359              0.359 
    ##            TNFRSF8             MAP2K6              FMNL3             FAM83H 
    ##              0.359              0.359              0.359              0.359 
    ##             GPR125              NXPH2            PIK3C2A             ZNF628 
    ##              0.359              0.359              0.359              0.359 
    ##               EDC4               FAXC              SEPT5            RASGRF2 
    ##              0.360              0.360              0.360              0.360 
    ##             RAVER1                IL7             ADRBK2            PDCD6IP 
    ##              0.360              0.360              0.360              0.360 
    ##             ALYREF              ETNK1              SMC1B              NR4A3 
    ##              0.360              0.360              0.360              0.360 
    ##             PARD6B              RAB5C              CCNB3              DDX23 
    ##              0.360              0.360              0.361              0.361 
    ##            PPP2R3A              TRAM1             MOSPD1               E2F8 
    ##              0.361              0.361              0.361              0.361 
    ##                ABR               ZIC3             BZRAP1                AES 
    ##              0.361              0.361              0.361              0.361 
    ##               POT1             CAPZA1               GAS7             SMURF1 
    ##              0.362              0.362              0.362              0.362 
    ##              WDR13              DHX40               E2F5               TNS1 
    ##              0.362              0.362              0.362              0.363 
    ##              MARK1              FOXN3              PGAP2             RIMBP2 
    ##              0.363              0.363              0.363              0.363 
    ##          MPHOSPH10              SEPT6             TUBB2A             ADAM22 
    ##              0.363              0.363              0.363              0.363 
    ##              SENP3             USP27X              MEX3D              DCTN1 
    ##              0.363              0.363              0.364              0.364 
    ##           PAFAH1B2              FUBP3               SRRT           COL4A3BP 
    ##              0.364              0.364              0.364              0.364 
    ##              SNX13               CDH4               MTA2              ACTA2 
    ##              0.364              0.364              0.364              0.364 
    ##               ZIC1            TBC1D20             CAMK2B              DDX27 
    ##              0.364              0.364              0.364              0.364 
    ##               ORC2              ZAP70             PTPN14             FAM13B 
    ##              0.364              0.364              0.365              0.365 
    ##               MMP2               STRN             TRMT1L               IPPK 
    ##              0.365              0.365              0.365              0.365 
    ##               DMTN            LRRC16A             FCHSD2               MDH1 
    ##              0.365              0.365              0.366              0.366 
    ##              MEGF8               PJA1              EIF3L              PEX26 
    ##              0.366              0.366              0.366              0.366 
    ##               NTF3              PRKG2               LTN1             PRDM10 
    ##              0.366              0.366              0.366              0.366 
    ##              PVRL3               DGKD               CDH9               NPAT 
    ##              0.366              0.366              0.367              0.367 
    ##             ZNF469              WDR20              HIPK3            ARHGDIG 
    ##              0.367              0.367              0.367              0.367 
    ##               MOB4            TRPC4AP             GABRA1            ADAMTS1 
    ##              0.367              0.367              0.367              0.367 
    ##              KCNH7              IP6K1              CYLC1               DGKA 
    ##              0.367              0.367              0.367              0.367 
    ##               MRAS            FAM214A               RGS7                DR1 
    ##              0.368              0.368              0.368              0.368 
    ##              IGSF3              OXSR1              RPS19               DTNA 
    ##              0.368              0.368              0.368              0.368 
    ##               BIN1              CALM2              NFRKB             PRDM12 
    ##              0.368              0.368              0.368              0.368 
    ##               RPSA              SMEK2            NCKAP5L               GPC4 
    ##              0.368              0.368              0.368              0.368 
    ##              SRSF6              LPPR4               DGKB               POGK 
    ##              0.369              0.369              0.369              0.369 
    ##               USO1              DMBX1           SLC25A28             GABPB1 
    ##              0.369              0.369              0.369              0.369 
    ##              PTBP3               SNX9               TLL1               PSD4 
    ##              0.369              0.369              0.370              0.370 
    ##               HMX3                BBX             POU4F3              CDC37 
    ##              0.370              0.370              0.370              0.370 
    ##             RPS27A             SACM1L            SLC13A4             PLXNB3 
    ##              0.370              0.370              0.370              0.371 
    ##             TFAP2C            APBB1IP             FNBP1L              AMER1 
    ##              0.371              0.371              0.371              0.371 
    ##               PIGR             PIEZO2           PPARGC1B              SEZ6L 
    ##              0.371              0.371              0.371              0.371 
    ##              DSTYK            B4GALT6             GABRA5              NLGN1 
    ##              0.371              0.371              0.372              0.372 
    ##             PHLDB1              SYNE2              SYT11              PCDH9 
    ##              0.372              0.372              0.372              0.372 
    ##               HAS3     RP11-793H13.10             MAP4K3             KCNK12 
    ##              0.372              0.372              0.372              0.372 
    ##              TERF1             PTPN23              ERO1L              TCF12 
    ##              0.372              0.372              0.372              0.372 
    ##             ZDHHC8               GNAS              LIN52               STAM 
    ##              0.372              0.372              0.372              0.372 
    ##               GAD2             ANAPC7            CCDC112               NRD1 
    ##              0.372              0.373              0.373              0.373 
    ##             DOCK11              PHKA2              CNTN1               CD47 
    ##              0.373              0.373              0.373              0.373 
    ##              NR1H3              TRPC4              AKAP1             FAM46A 
    ##              0.373              0.373              0.373              0.374 
    ##             KANSL2              RNF43              KAT2A               ANLN 
    ##              0.374              0.374              0.374              0.374 
    ##              BACE1             GPR143            MSANTD2              KCNT2 
    ##              0.374              0.374              0.374              0.374 
    ##             RSBN1L               DGKH              MYO10          SECISBP2L 
    ##              0.374              0.374              0.374              0.374 
    ##              FHOD3               XIAP             ZC3H15              MFAP1 
    ##              0.374              0.374              0.374              0.374 
    ##             DNMT3B             SNAP23            PPP2R2D              NHSL1 
    ##              0.375              0.375              0.375              0.375 
    ##              STRAP            SLC29A1             SLC9A1               WASL 
    ##              0.375              0.375              0.375              0.375 
    ##             TRIM11             DDX19A            SLC4A10              RPL17 
    ##              0.375              0.375              0.375              0.375 
    ##              EPHB3             RMND5A            CHORDC1               DGKZ 
    ##              0.375              0.375              0.375              0.375 
    ##             SLAIN1             CDC25A               ST14               PES1 
    ##              0.375              0.375              0.375              0.376 
    ##             CACNB1             CACNG2               TFEB             VPS37D 
    ##              0.376              0.376              0.376              0.376 
    ##            DENND5A             HS6ST1            AKIRIN2            PPP2R2C 
    ##              0.376              0.376              0.376              0.376 
    ##              RPL26               OTX2            HNRNPA0               GRM3 
    ##              0.376              0.376              0.376              0.376 
    ##               DDR2             DMRTA2              PRPS1                DSE 
    ##              0.376              0.376              0.377              0.377 
    ##            SLITRK5           TNKS1BP1            ATXN7L2            SEC14L1 
    ##              0.377              0.377              0.377              0.377 
    ##              IKZF4             ANGPT2             SCAPER              CCNE2 
    ##              0.377              0.377              0.377              0.377 
    ##            SLC2A13            FAM196A              UBE2F               ACLY 
    ##              0.377              0.377              0.377              0.377 
    ##             MAP2K1               DPP6             ACVR2B             CDK5R1 
    ##              0.377              0.377              0.377              0.378 
    ##             CTNNA2              RAB30              SCN4A               DLX6 
    ##              0.378              0.378              0.378              0.378 
    ##            HDGFRP2             CLASRP               MAP7           SLC22A23 
    ##              0.378              0.378              0.378              0.378 
    ##              STIM1              SRSF4             ZFAND3            PRPF38B 
    ##              0.378              0.378              0.378              0.379 
    ##              UBAP1              NIPA2              ATOH8              GNAI1 
    ##              0.379              0.379              0.379              0.379 
    ##               GID4               WWP2              PASD1              CDC16 
    ##              0.379              0.379              0.379              0.379 
    ##              CAPZB              SENP7              CASC3              DOCK7 
    ##              0.379              0.379              0.379              0.379 
    ##              FSTL1               RPS3             MESDC1             CACYBP 
    ##              0.379              0.380              0.380              0.380 
    ##              CADM2             PCDH20              KCNA2             GOLGA3 
    ##              0.380              0.380              0.380              0.380 
    ##             POU3F2              CHMP7               CALR              HMHA1 
    ##              0.380              0.380              0.380              0.380 
    ##              NR2C1               CD81              BMP2K               MPP6 
    ##              0.381              0.381              0.381              0.381 
    ##              DDX24           FAM171A2               BCL9               KAZN 
    ##              0.381              0.381              0.381              0.381 
    ##            TMEM204            PPP2R5A              EPAS1               UNKL 
    ##              0.381              0.381              0.381              0.381 
    ##              MAST4            NUP210L            CAPRIN2             TSPAN5 
    ##              0.381              0.381              0.381              0.381 
    ##             SREBF2             ADARB2               LGR4              LOXL1 
    ##              0.382              0.382              0.382              0.382 
    ##             ELAVL4              PTPN1               NUDC               IRF4 
    ##              0.382              0.382              0.382              0.382 
    ##               GRK5               XKR7              NRCAM             PPAP2B 
    ##              0.382              0.382              0.383              0.383 
    ##             TGFBR1                PKM               SS18              LTBP2 
    ##              0.383              0.383              0.383              0.383 
    ##             CSNK2B              BCAR3               STIL              DCAF6 
    ##              0.383              0.383              0.383              0.384 
    ##              APBA2               GDI2              PHRF1             PLXNA3 
    ##              0.384              0.384              0.384              0.384 
    ##               ZHX1            CLEC16A               IDH2            ABHD16A 
    ##              0.384              0.384              0.384              0.384 
    ##               BTG3              CD79B            GORASP2             CCSER2 
    ##              0.384              0.384              0.384              0.384 
    ##             MARCH7              MED26              MEOX2             CLDN11 
    ##              0.384              0.384              0.385              0.385 
    ##              IARS2             IGFBP3             GABRB2             GABRG3 
    ##              0.385              0.385              0.385              0.385 
    ##             PLXNB1             RPL27A              WIPF1                SYP 
    ##              0.385              0.385              0.385              0.385 
    ##             KCNAB2               NGEF              PMP22              TMCC1 
    ##              0.385              0.385              0.386              0.386 
    ##               NEK7             GRIN3A              LRFN1              SCRT2 
    ##              0.386              0.386              0.386              0.386 
    ##              ALCAM              PCBP3              DCTN4             NDUFS7 
    ##              0.386              0.386              0.386              0.386 
    ##             EIF1AX              PSMB2              VCAM1              LRFN4 
    ##              0.386              0.386              0.386              0.387 
    ##              GRASP              PCBP1           TMEM132C              ABCC5 
    ##              0.387              0.387              0.387              0.387 
    ##             DOCK10                ALB              MYOCD              KCNH1 
    ##              0.387              0.387              0.387              0.387 
    ##              MAGI3             VANGL2               POMP               PHC2 
    ##              0.387              0.387              0.387              0.387 
    ##              RPL10              GATA3              APBA1               SPIB 
    ##              0.387              0.388              0.388              0.388 
    ##              ACER3             ZBTB33               FLCN              WNT9A 
    ##              0.388              0.388              0.388              0.388 
    ##              UBE4A              MYEF2              BTBD2              HERC3 
    ##              0.388              0.388              0.388              0.388 
    ##              PRKCH               EPN1            PLEKHA2              PNRC1 
    ##              0.388              0.388              0.388              0.389 
    ##              WSCD2              PVRL4            PIKFYVE               ETS1 
    ##              0.389              0.389              0.389              0.389 
    ##               HAT1            PROSER1                ARX              RWDD1 
    ##              0.389              0.389              0.389              0.390 
    ##              ZC4H2             SPTLC2                C1R            TMEM135 
    ##              0.390              0.390              0.390              0.390 
    ##               SRPR              RAP1B            ADAMTS9              RAB18 
    ##              0.390              0.390              0.390              0.390 
    ##            B3GALT2             PIK3R4             UBE2J2             LONRF1 
    ##              0.390              0.390              0.391              0.391 
    ##              AURKA              TCEA1              PARP8              PRKCQ 
    ##              0.391              0.391              0.391              0.391 
    ##               BLNK             PIK3R5              FGF14             NAP1L4 
    ##              0.391              0.391              0.392              0.392 
    ##              WDR47      RP11-603J24.9              ZNRF1              ADNP2 
    ##              0.392              0.392              0.392              0.392 
    ##             SPOCK1              BMP10               PKN1              KCNA3 
    ##              0.392              0.392              0.392              0.392 
    ##              FOXB1               PLEC              CELF3              SNRPE 
    ##              0.392              0.392              0.392              0.392 
    ##             MARCH5             RPL23A              PTAR1             MGAT4B 
    ##              0.392              0.393              0.393              0.393 
    ##              SYDE1             GTF2E1               CDH7              GNAI2 
    ##              0.393              0.393              0.393              0.393 
    ##              APLP2              SPARC              CSF1R             ARGLU1 
    ##              0.393              0.393              0.393              0.394 
    ##              CTPS2            ONECUT2              NR1D1           C18orf25 
    ##              0.394              0.394              0.394              0.394 
    ##              CAPN6              CLCN7              FOXN4               NCLN 
    ##              0.394              0.394              0.395              0.395 
    ##               OPN5              TNIP1              PRKCA              TXNL1 
    ##              0.395              0.395              0.395              0.395 
    ##            SLC12A6             ZNF287              PHF14             PARP14 
    ##              0.395              0.395              0.395              0.395 
    ##              SOGA3              PDE4B              GNA11              VPS18 
    ##              0.395              0.395              0.395              0.395 
    ##             ZNF367              URGCP             HSPA14             FAM58A 
    ##              0.395              0.395              0.395              0.396 
    ##              WDR82              CABP1               PLAA              RPS14 
    ##              0.396              0.396              0.396              0.396 
    ##              AXIN2                VCL             EXOC6B              ERBB2 
    ##              0.396              0.396              0.396              0.397 
    ##              AKAP8              PPP6C              ABCC1       TRIM39-RPP21 
    ##              0.397              0.397              0.397              0.397 
    ##               MAFB               RTN4             ZBTB46              COPG1 
    ##              0.397              0.397              0.397              0.397 
    ##             DNAJB6               EPN2              FOXK1               GCH1 
    ##              0.397              0.397              0.398              0.398 
    ##           KIAA1432             SHISA6              KCNN1               SDK1 
    ##              0.398              0.398              0.398              0.398 
    ##              KNDC1            TMEM63C              USP43               SYT3 
    ##              0.398              0.398              0.398              0.398 
    ##               NCAN              PPIL3                MKX              MXRA5 
    ##              0.398              0.398              0.398              0.398 
    ##              TAF15           KIAA1107            SNRNP40             KLHL36 
    ##              0.398              0.398              0.398              0.398 
    ##             LRRC4B                HLX               RGS6             PI4K2A 
    ##              0.398              0.399              0.399              0.399 
    ##              SF3A2               LCP1             ZNF598             GABRG2 
    ##              0.399              0.399              0.399              0.399 
    ##            JAKMIP3             ADARB1               KRT1             DCAF10 
    ##              0.399              0.399              0.399              0.399 
    ##                ATR            XPNPEP1              PDZD4             ZSWIM5 
    ##              0.399              0.399              0.400              0.400 
    ##              USP31             ZNF106             GPR116             LAPTM5 
    ##              0.400              0.400              0.400              0.400 
    ##              NCAM2              PRKCZ            PRPF40A               CIZ1 
    ##              0.400              0.400              0.400              0.400 
    ##              SMAD3             CLINT1              KCNK3             PCDH17 
    ##              0.400              0.401              0.401              0.401 
    ##               GCLC              AJUBA               STX5             SLAIN2 
    ##              0.401              0.401              0.402              0.402 
    ##              EDRF1               GARS              PPM1H              LUZP1 
    ##              0.402              0.402              0.402              0.402 
    ##               ARF1             SORBS2              TCF25              CLSPN 
    ##              0.402              0.402              0.402              0.402 
    ##               AATK             PCYT1B             HECTD3              AXIN1 
    ##              0.403              0.403              0.403              0.403 
    ##              N4BP1               TANK               RORA               TNKS 
    ##              0.404              0.404              0.404              0.404 
    ##              AP1M1               CDH5              NOP56            TMEM55B 
    ##              0.404              0.404              0.404              0.404 
    ##              LAMP1               WNT3              AKAP9              CHSY1 
    ##              0.404              0.404              0.404              0.405 
    ##              JAZF1            IRF2BPL            METTL16               UBA3 
    ##              0.405              0.405              0.405              0.405 
    ##             SHISA9               SZT2             TSPYL2              ZBED4 
    ##              0.405              0.405              0.405              0.405 
    ##              DISP2              HDAC1           ARHGAP36           ARHGAP33 
    ##              0.405              0.405              0.405              0.406 
    ##              PSMB7              FGF13              ZC3H6             ABHD13 
    ##              0.406              0.406              0.406              0.406 
    ##              RPL14             GABRA4            MAPKBP1             MAPRE2 
    ##              0.406              0.406              0.406              0.406 
    ##              YIPF6              UHMK1                C4B                OTC 
    ##              0.406              0.406              0.406              0.406 
    ##            SERINC3             ATP11B              MECP2                MCU 
    ##              0.406              0.406              0.407              0.407 
    ##              LMOD1             ATP8B2              HSPA5               CAP2 
    ##              0.407              0.407              0.407              0.407 
    ##             ELOVL4             LRRC47               TAF7               HEY2 
    ##              0.407              0.407              0.407              0.407 
    ##              RNGTT              MAPK6               KLC1               MYOG 
    ##              0.407              0.407              0.408              0.408 
    ##              MYSM1              PELI2             MAP3K9              FSCN1 
    ##              0.408              0.408              0.408              0.408 
    ##            FAM171B                GAK           KIAA1551            PRPF38A 
    ##              0.408              0.408              0.408              0.408 
    ##             LPCAT4             SCNN1G               ISL1            RALGPS2 
    ##              0.408              0.408              0.408              0.408 
    ##            TNFSF12    TNFSF12-TNFSF13              GRIK4              CLDN5 
    ##              0.408              0.408              0.409              0.409 
    ##              HOXB3            ARHGAP1             TSGA10             RPL35A 
    ##              0.409              0.409              0.409              0.409 
    ##              GNA13               IRX3              SESN3               JPH1 
    ##              0.409              0.409              0.409              0.409 
    ##              DCP1A            MAP3K10            HEPACAM               EMX1 
    ##              0.409              0.409              0.409              0.409 
    ##            DENND2A              KCNQ4              SRRM4              ZC3H3 
    ##              0.410              0.410              0.410              0.410 
    ##              CCND1             ZNF750           KIAA1109              CDH20 
    ##              0.410              0.410              0.410              0.411 
    ##               AAMP               BMP6              MMP24              PCGF2 
    ##              0.411              0.411              0.411              0.411 
    ##              NACAD              RNF10               LMO7              PLCG1 
    ##              0.411              0.411              0.411              0.411 
    ##             ZNF652              CLIP3             PPP6R2              PTH1R 
    ##              0.411              0.411              0.411              0.411 
    ##              HSPG2              FTSJ3             KCTD10              HDAC3 
    ##              0.412              0.412              0.412              0.412 
    ##              KLHL9              PTOV1               EML5              ITGAL 
    ##              0.412              0.412              0.412              0.412 
    ##              WDR37              MEX3C              KCTD3               MYCN 
    ##              0.412              0.412              0.412              0.412 
    ##            PLEKHG1               SMG5          UHRF1BP1L                FYB 
    ##              0.412              0.412              0.412              0.412 
    ##              NCALD             PABPC5               IBTK              LEMD2 
    ##              0.413              0.413              0.413              0.413 
    ##             ZNF513              FLRT2              SBNO2             CACNG7 
    ##              0.413              0.413              0.413              0.413 
    ##            SLC38A4               MLF2             B3GNT2              TAF4B 
    ##              0.413              0.413              0.413              0.413 
    ##               DTX3              CD2AP               GCC2               SV2A 
    ##              0.413              0.413              0.414              0.414 
    ##              LMNB1            C1orf21               SUB1               DCP2 
    ##              0.414              0.414              0.414              0.414 
    ##               ST13               HCCS                IDE              TMED2 
    ##              0.414              0.414              0.414              0.414 
    ##            FAM205A            ST8SIA3               HSF1               FEZ1 
    ##              0.414              0.414              0.414              0.414 
    ##             ZBTB41               RRM2               NEMF               GRM1 
    ##              0.414              0.415              0.415              0.415 
    ##             SEL1L3              GLIS2              LMX1B              SOX17 
    ##              0.415              0.415              0.415              0.415 
    ##             FAM98A              USP22               BUB3              EIF3J 
    ##              0.415              0.415              0.416              0.416 
    ##              RAB1B             OSBPL8            NEUROD6             HOXA10 
    ##              0.416              0.416              0.416              0.416 
    ##             ZNF787            GPRASP1             ARPC1A               LRP4 
    ##              0.416              0.416              0.416              0.416 
    ##              PDE7B             LRRTM2             LINGO2             FAM53B 
    ##              0.416              0.416              0.416              0.416 
    ##           SLC25A12             ZNF746              NEDD1             HPCAL1 
    ##              0.416              0.416              0.416              0.416 
    ##              ADAM9              CLVS1             ZNF606               ACO2 
    ##              0.417              0.417              0.417              0.417 
    ##                NOG               CCNY               FMR1               RTN2 
    ##              0.417              0.417              0.417              0.417 
    ##            CRAMP1L               AKT2              GNPAT               LRBA 
    ##              0.417              0.417              0.417              0.417 
    ##              CRTC1              FKBP5              INSM2               CDK8 
    ##              0.418              0.418              0.418              0.418 
    ##              PITX1               EI24             GPCPD1              IPO11 
    ##              0.418              0.418              0.418              0.418 
    ##            ANKRD62               ABI1               CCNF               DAB2 
    ##              0.418              0.418              0.419              0.419 
    ##              FNDC4              KITLG              CDH10            FAM214B 
    ##              0.419              0.419              0.419              0.419 
    ##               MYO6              FRMD7                APP              SCYL2 
    ##              0.419              0.419              0.419              0.419 
    ##             GPR107              ELAC1            ZNF280D             NUP160 
    ##              0.419              0.419              0.419              0.419 
    ##             ZNF202              DUSP4             KDELR1           C6orf106 
    ##              0.419              0.419              0.419              0.419 
    ##              MIER1              NDST2            BCL2L11             USP6NL 
    ##              0.419              0.419              0.420              0.420 
    ##               ADAR            DNTTIP1             MBOAT2            RNF113A 
    ##              0.420              0.420              0.420              0.420 
    ##              CXXC4              CXXC5             GOLGB1               IRF2 
    ##              0.420              0.420              0.420              0.420 
    ##              PCDH7               TBK1             CC2D1A             POU6F1 
    ##              0.420              0.420              0.420              0.420 
    ##            GRID2IP             MFAP3L             FAM84A              IFFO2 
    ##              0.420              0.420              0.420              0.420 
    ##               PIM1            CXorf23              INTS5                ERH 
    ##              0.420              0.421              0.421              0.421 
    ##              CYTH3              VAC14              MMS19               TAF2 
    ##              0.421              0.421              0.421              0.421 
    ##           KIAA1217              ITGA9               KLF9           KIAA0100 
    ##              0.421              0.421              0.421              0.421 
    ##              SHPRH             SEZ6L2               PARG             POU4F1 
    ##              0.421              0.421              0.421              0.421 
    ##               EZH1              SCYL3             ELMOD1            CNTNAP1 
    ##              0.421              0.421              0.422              0.422 
    ##      RP11-598P20.5            SLC6A17              EPHB4            PLEKHA6 
    ##              0.422              0.422              0.422              0.422 
    ##               DVL3               CD40              TAGAP            SLITRK4 
    ##              0.422              0.422              0.422              0.422 
    ##              OVOL1                LBR            RSL24D1              VAMP2 
    ##              0.422              0.422              0.422              0.422 
    ##              USP46              ACAP1             INPPL1              RPS13 
    ##              0.423              0.423              0.423              0.423 
    ##            RTN4RL2              ZBED5              SYNE1               TAB3 
    ##              0.423              0.423              0.423              0.423 
    ##            ZDHHC14               TOB2               SNX6               FZD8 
    ##              0.423              0.423              0.423              0.423 
    ##             ANP32E              CDH22             NOS1AP                VIM 
    ##              0.423              0.423              0.423              0.423 
    ##             RBPMS2              PRDM6                MAG            TRABD2B 
    ##              0.423              0.424              0.424              0.424 
    ##             SLC1A2             FBXO15              RUNX3                NEB 
    ##              0.424              0.424              0.424              0.424 
    ##              TASP1             HMGXB4              RPS12                EHF 
    ##              0.424              0.424              0.425              0.425 
    ##                KLB              MBNL1             CDC25B              STK31 
    ##              0.425              0.425              0.425              0.425 
    ##              STOX2            ZDHHC15             PDLIM7                FGR 
    ##              0.425              0.425              0.425              0.425 
    ##           ANKRD13B               DENR              TTLL7               VBP1 
    ##              0.425              0.425              0.426              0.426 
    ##              CERS2               MUC2              LRIG1            EPB41L5 
    ##              0.426              0.426              0.426              0.426 
    ##           TGFBRAP1             AKAP12            TBC1D15              SMPD3 
    ##              0.426              0.426              0.426              0.427 
    ##              RRAGB              LRRN2             SLC6A4               TBX1 
    ##              0.427              0.427              0.427              0.427 
    ##             ATP10A              PEG10             OSGIN2             RNF123 
    ##              0.427              0.427              0.427              0.427 
    ##           TMEM183A            CXorf65             KIF21B             AMOTL2 
    ##              0.427              0.428              0.428              0.428 
    ##             SORCS2              LRTM2           ANKRD13D               NSG2 
    ##              0.428              0.428              0.428              0.428 
    ##               SIX4              PRDM8            IGF2BP2              GRHL1 
    ##              0.428              0.428              0.428              0.428 
    ##             ARFIP2             DUSP10              GRIA4               MPP5 
    ##              0.429              0.429              0.429              0.429 
    ##                NIN             POU3F3             BAIAP2               HM13 
    ##              0.429              0.429              0.429              0.429 
    ##              TANC1             RPL13A               SNCA              KIF14 
    ##              0.429              0.429              0.429              0.429 
    ##             TROVE2            ATP6AP2              C4BPA                CFP 
    ##              0.429              0.429              0.430              0.430 
    ##               IRS4              CCDC6               CTTN               CLK3 
    ##              0.430              0.430              0.430              0.430 
    ##               CDH1              ZFP82               UTRN              MYO1B 
    ##              0.430              0.430              0.430              0.430 
    ##              PDE1B             NFKBIE              UBE2N               MTDH 
    ##              0.430              0.430              0.431              0.431 
    ##           ADAMTS16               KLF3               MYNN              ABHD2 
    ##              0.431              0.431              0.431              0.431 
    ##           ARHGAP20              ZFP62              NOL10               FAT1 
    ##              0.431              0.431              0.431              0.431 
    ##           PPP1R12C              AP4E1              ENPP2               RPA1 
    ##              0.431              0.431              0.431              0.431 
    ##               CCNI               RBMX              TAF10                FGB 
    ##              0.432              0.432              0.432              0.432 
    ##             HMG20A            HNRNPAB               EMC7            SERTAD2 
    ##              0.432              0.432              0.432              0.432 
    ##              BACH1            DCUN1D3               DET1             EIF2B3 
    ##              0.432              0.432              0.432              0.432 
    ##            TP53BP2                CPD            PITPNM1              CHTOP 
    ##              0.433              0.433              0.433              0.433 
    ##             SLC1A6            SLC6A14             GTPBP2              ARAP3 
    ##              0.433              0.433              0.433              0.433 
    ##            TMEM87A              TIAM2               CMAS               WWC1 
    ##              0.433              0.433              0.433              0.433 
    ##                PML               RPS9              SIRPA             TMEM11 
    ##              0.433              0.434              0.434              0.434 
    ##             MAN1C1              ACOT7            EXOC3L2            GUCY1A2 
    ##              0.434              0.434              0.434              0.434 
    ##           KIAA1671            SPATS2L               PAX2             COBLL1 
    ##              0.434              0.434              0.434              0.434 
    ##               FAR2            MACROD2      RP11-343C2.11               PREP 
    ##              0.434              0.434              0.434              0.434 
    ##              ATG2B              BIRC2               AGRN              NEDD9 
    ##              0.435              0.435              0.435              0.435 
    ##              THOC5              ABCD3              H2AFZ              CRADD 
    ##              0.435              0.435              0.435              0.435 
    ##             PLXNA2              CASD1                 PC              CADPS 
    ##              0.435              0.435              0.435              0.435 
    ##             CHAF1B             GPR137              MAP1S             RNF146 
    ##              0.435              0.435              0.436              0.436 
    ##              EFR3B              SORL1              RRAGD              GAREM 
    ##              0.436              0.436              0.436              0.436 
    ##             KLHL22           SUV420H2            PHYHIPL              KLHL4 
    ##              0.436              0.436              0.436              0.436 
    ##              TPCN1              CADM1             INTS12             ZFYVE1 
    ##              0.436              0.436              0.436              0.436 
    ##           TBC1D22A             KCNIP4              SKOR2              TBX15 
    ##              0.436              0.436              0.436              0.437 
    ##              FOXL2              LIMD1               KLF5               ANO4 
    ##              0.437              0.437              0.437              0.437 
    ##              MYBL2              KPNA2              APLP1                QKI 
    ##              0.437              0.437              0.437              0.437 
    ##                AXL              CNTN2              GPBP1              GLP1R 
    ##              0.437              0.437              0.437              0.437 
    ##              NFIL3             MMS22L              RBM42              RUNX1 
    ##              0.437              0.438              0.438              0.438 
    ##            SLC41A1               EMC4             RNF139             TMED10 
    ##              0.438              0.438              0.438              0.438 
    ##               PDXK            TNFSF13       CTD-2545M3.6               KSR1 
    ##              0.438              0.438              0.438              0.438 
    ##               GLCE               DRD2              ADAP1             ARPP21 
    ##              0.438              0.439              0.439              0.439 
    ##               RRP9              SALL3             PSENEN              MTCH1 
    ##              0.439              0.439              0.439              0.439 
    ##            KATNAL1             SEC31A             PHYHIP               DNER 
    ##              0.439              0.439              0.439              0.439 
    ##               GAB1              PRKCI              SOX11                PNN 
    ##              0.439              0.439              0.440              0.440 
    ##             VPS33A              SASH1              APBB2                EZR 
    ##              0.440              0.440              0.440              0.440 
    ##               VEZT             CSF2RB               RHOQ              SRP14 
    ##              0.441              0.441              0.441              0.441 
    ##              RIPK2             SEMA7A             EIF2B5                FAS 
    ##              0.441              0.441              0.441              0.441 
    ##               KRT5               EPT1             AKAP10              PRPF6 
    ##              0.441              0.441              0.441              0.441 
    ##              MCTS1              DOCK5               SSR4             NELFCD 
    ##              0.441              0.441              0.441              0.441 
    ##               TIA1           ARHGAP17              ADCY6               ENY2 
    ##              0.441              0.442              0.442              0.442 
    ##              MYO1E              KCNH5              KCTD1             TUBB4B 
    ##              0.442              0.442              0.442              0.442 
    ##              SULF2              PEAK1              FRMD3             GPR173 
    ##              0.442              0.442              0.442              0.442 
    ##          ARHGEF10L                STS              CEP68            SLCO3A1 
    ##              0.442              0.442              0.442              0.442 
    ##              KRT18                GCK               NSMF                PLG 
    ##              0.443              0.443              0.443              0.443 
    ##              CD163               OTX1             PPP3CC              COX5B 
    ##              0.443              0.443              0.443              0.443 
    ##             ZNF473           TNFRSF14               FLNB              ITGAV 
    ##              0.443              0.443              0.443              0.443 
    ##            RUNDC3A             ATP2B4             RPL10A             DMRTC2 
    ##              0.443              0.444              0.444              0.444 
    ##             ZNF507              LIN7B              RSBN1             IL1RAP 
    ##              0.444              0.444              0.444              0.444 
    ##              ARRB2              FSTL4             OSBPL7             GTF3C4 
    ##              0.444              0.444              0.444              0.444 
    ##              FADS1              BRIX1             MAMLD1               PPIE 
    ##              0.445              0.445              0.445              0.445 
    ##              RAD17               EMC8              SNTG1               PSD3 
    ##              0.445              0.445              0.445              0.445 
    ##               LSM6                TNC              UBXN1           SERPINE2 
    ##              0.445              0.445              0.445              0.445 
    ##            DENND4C            RPS6KA2              SSBP4            EIF2AK3 
    ##              0.445              0.445              0.445              0.445 
    ##             METTL9             CSRNP3               TBX4              CASS4 
    ##              0.445              0.445              0.445              0.446 
    ##              EHMT2            FAM222A            KATNBL1              ASCC2 
    ##              0.446              0.446              0.446              0.446 
    ##              NDRG3           KIAA1210            FAM129B              ATCAY 
    ##              0.446              0.446              0.446              0.446 
    ##             DCBLD2             PFKFB3               UBR7            ALDH1A2 
    ##              0.446              0.447              0.447              0.447 
    ##             TCF7L1             NFATC1               BAG3            RASGRP2 
    ##              0.447              0.447              0.447              0.447 
    ##      RP11-599B13.6               TYMS        FXYD6-FXYD2              EGLN2 
    ##              0.447              0.447              0.447              0.447 
    ##           TMEM200A             FBXO30               HCN2               BMS1 
    ##              0.447              0.448              0.448              0.448 
    ##               SAT1             UBQLN2             YME1L1             NUP133 
    ##              0.448              0.448              0.448              0.448 
    ##             PTDSS1             TSPYL5                ALK             PHLDA1 
    ##              0.448              0.448              0.448              0.448 
    ##             MTHFD1              CPNE6            CACNA1F              CDIP1 
    ##              0.448              0.448              0.448              0.449 
    ##               NMT2            ZKSCAN5              SOGA2              PLCB4 
    ##              0.449              0.449              0.449              0.449 
    ##             ZNF410              DESI2               DLST              HMGA2 
    ##              0.449              0.449              0.449              0.449 
    ##           ARHGAP39           MAPKAPK5              NLRC5             RALBP1 
    ##              0.449              0.449              0.449              0.449 
    ##              HOXB8              CMPK1          KIAA0319L             CARD10 
    ##              0.449              0.449              0.449              0.449 
    ##              NR5A2               RGL2              APH1A             SREBF1 
    ##              0.449              0.450              0.450              0.450 
    ##     TGIF2-C20orf24               LIG1               RPS2               ICMT 
    ##              0.450              0.450              0.450              0.450 
    ##              CYR61              FOXF2               IRF9             BRINP3 
    ##              0.450              0.450              0.450              0.450 
    ##      CTD-2207O23.3           KIAA1033               EYA4              SUMO2 
    ##              0.450              0.450              0.450              0.450 
    ##               CASR               PRR7              CNIH3              PTPRG 
    ##              0.451              0.451              0.451              0.451 
    ##               SYT1            IRF2BP2              KIF3B             ZFAND5 
    ##              0.451              0.451              0.451              0.451 
    ##               IER5             DDX39A              DDX18              LUC7L 
    ##              0.451              0.451              0.452              0.452 
    ##              SRP68            C1orf51             CNOT11            N4BP2L1 
    ##              0.452              0.452              0.452              0.452 
    ##           EPM2AIP1              TRPM6               SOST               EML6 
    ##              0.452              0.452              0.453              0.453 
    ##               INSR             PLXDC2             UBE2G1               ALG5 
    ##              0.453              0.453              0.453              0.453 
    ##              RBCK1             POLR3E              TRPM3            PLEKHF2 
    ##              0.453              0.453              0.453              0.453 
    ##            PACSIN2                TSN             ZNF444            CCDC88C 
    ##              0.454              0.454              0.454              0.454 
    ##               THRA              ACVR1              SNX30              CDH18 
    ##              0.454              0.454              0.454              0.454 
    ##              AP1S2               AKNA              NUP54               ADSS 
    ##              0.454              0.454              0.454              0.454 
    ##             NCKAP5             CCP110              SUGT1              EPS15 
    ##              0.454              0.454              0.455              0.455 
    ##              CAAP1           ARHGEF18              RBMS3               EML4 
    ##              0.455              0.455              0.455              0.455 
    ##              PCMT1             SPIRE1              PHF5A            ZFYVE20 
    ##              0.455              0.455              0.455              0.455 
    ##             PRDM15             KIF21A              USPL1               BMP1 
    ##              0.455              0.456              0.456              0.456 
    ##              RCOR2              TRAF7               BTRC             PCDH18 
    ##              0.456              0.456              0.456              0.456 
    ##            ZFYVE28               EPG5                HGS               RYR1 
    ##              0.456              0.456              0.456              0.456 
    ##              SRSF5             ZNF579               COPE              ITFG1 
    ##              0.457              0.457              0.457              0.457 
    ##              NAT8L            FAM122B              GDPD2            TSPAN12 
    ##              0.457              0.457              0.457              0.457 
    ##               NFX1             FBXO31             FBXO46            TBC1D14 
    ##              0.457              0.457              0.457              0.457 
    ##             LANCL2              RPL21             SHANK2              IKZF5 
    ##              0.457              0.458              0.458              0.458 
    ##             HMG20B             ZNF395               VPS8           KIAA1598 
    ##              0.458              0.458              0.458              0.458 
    ##              RUFY3               SIM2               DLK2              HMGN2 
    ##              0.458              0.458              0.458              0.459 
    ##              CSPG4              SYT13              GULP1            SLC30A5 
    ##              0.459              0.459              0.459              0.459 
    ##              KIFC3              RPS25              GPR61            RASGRP4 
    ##              0.459              0.459              0.459              0.459 
    ##               RPN2               OSR1             SLC9C2               NAB2 
    ##              0.460              0.460              0.460              0.460 
    ##              ZDBF2                UXT               CA10               ABI2 
    ##              0.460              0.460              0.460              0.460 
    ##               JUNB             MYBPC1          C10orf118              PI4KA 
    ##              0.460              0.460              0.461              0.461 
    ##             ANAPC4              PLCE1               RXRG              RCAN2 
    ##              0.461              0.461              0.461              0.461 
    ##            SLC35A1             PRDM13              ABCG1              PRKD3 
    ##              0.461              0.461              0.461              0.461 
    ##               APOB              SULF1           RAPGEFL1             ELAVL3 
    ##              0.461              0.461              0.461              0.462 
    ##             MAN1A1              TRIB1               WSB1              ACOX1 
    ##              0.462              0.462              0.462              0.462 
    ##             MVB12B             SEMA6B               NRG3            SPATA13 
    ##              0.462              0.462              0.462              0.462 
    ##               FGD4             ELOVL1              LRFN5                MAX 
    ##              0.462              0.462              0.462              0.462 
    ##              PTPN9               ESF1              EDEM3            SLC38A9 
    ##              0.462              0.463              0.463              0.463 
    ##              UTP20                TTK              PDGFB             LONRF3 
    ##              0.463              0.463              0.463              0.463 
    ##              CHST1             GAS2L1              SUMO1              GSG1L 
    ##              0.463              0.463              0.463              0.464 
    ##            POLDIP2              APAF1               NRN1            FAM126B 
    ##              0.464              0.464              0.464              0.464 
    ##              IKBKE            SUPV3L1             TRIM26              ABTB2 
    ##              0.464              0.464              0.464              0.464 
    ##               ODF2              RPLP0               FRS3             HS2ST1 
    ##              0.464              0.464              0.464              0.465 
    ##            SLC39A8             SORCS3               GDF5             ABLIM2 
    ##              0.465              0.465              0.465              0.465 
    ##              ATG14              PARVA               PGK1              MGRN1 
    ##              0.465              0.465              0.465              0.465 
    ##      RP11-1220K2.2             ZNF449              RPS23             DMRTB1 
    ##              0.465              0.465              0.465              0.466 
    ##              WDR36            C1GALT1              HELZ2            LRRC16B 
    ##              0.466              0.466              0.466              0.466 
    ##            SMARCD3             MAPK10              SOCS5            TOMM40L 
    ##              0.466              0.466              0.466              0.466 
    ##               PER1              ARMC1           C10orf54              DACT1 
    ##              0.466              0.466              0.467              0.467 
    ##             TBKBP1              WASF3               SBF2               RAC1 
    ##              0.467              0.467              0.467              0.467 
    ##               MEST             SNRPD3             GPR158              PRR14 
    ##              0.467              0.467              0.467              0.467 
    ##            PLEKHM1              BCAS3             ATP11A             TMEFF1 
    ##              0.467              0.467              0.467              0.467 
    ##               SMOX            DCLRE1B               EYA2             GPR126 
    ##              0.467              0.467              0.467              0.468 
    ##             THAP11              MCRS1             AKAP8L                BRE 
    ##              0.468              0.468              0.468              0.468 
    ##             KIF26A            ADIPOR1            SLC17A6                PHB 
    ##              0.468              0.469              0.469              0.469 
    ##             IQSEC3             FBXL14             MGAT5B               TP53 
    ##              0.469              0.469              0.469              0.469 
    ##               COBL             INO80B             ZNF282             CAMKK2 
    ##              0.469              0.469              0.469              0.470 
    ##             ZNF805              NPRL2           PPP1R12B                PGR 
    ##              0.470              0.470              0.470              0.470 
    ##             RILPL1             ZNF691              CCNL2            ZNF385D 
    ##              0.470              0.470              0.470              0.470 
    ##             AMOTL1               SMG8              SOX12            ATF7IP2 
    ##              0.471              0.471              0.471              0.471 
    ##              KCNT1               HID1             RANBP1              CYTH2 
    ##              0.471              0.471              0.471              0.471 
    ##               FZD2               DKK2            TBC1D25             SWAP70 
    ##              0.471              0.471              0.471              0.471 
    ##               ENO2             ZNF526              RAB1A             PIH1D3 
    ##              0.471              0.471              0.471              0.472 
    ##               IMMT              ITPKA              ITPKC              SOX14 
    ##              0.472              0.472              0.472              0.472 
    ##               FSD1              CWC22             ARID3B            AMMECR1 
    ##              0.472              0.472              0.472              0.472 
    ##              RING1               KAT5              ABCB1             ERGIC1 
    ##              0.472              0.472              0.473              0.473 
    ##             ZNF274               BATF            CCDC85C             MAPK14 
    ##              0.473              0.473              0.473              0.473 
    ##               TMX2              CPNE4              NR1H2               LHX8 
    ##              0.473              0.473              0.473              0.473 
    ##              TRPC1              HMGB1             RNF19A             MLXIPL 
    ##              0.473              0.473              0.473              0.473 
    ##              PDIA6            LRRC37B              LTBP1              ITGA4 
    ##              0.473              0.474              0.474              0.474 
    ##              ZFPM1               GAD1              RPH3A              NR4A1 
    ##              0.474              0.474              0.474              0.474 
    ##              CAMK4              PTPRE              PTPRO               TSKU 
    ##              0.474              0.474              0.474              0.474 
    ##               RFX6            SLC26A9             SHISA7            DCAF8L2 
    ##              0.474              0.474              0.474              0.474 
    ##             SPOCK2              ACSL6              FREM2              KANK2 
    ##              0.475              0.475              0.475              0.475 
    ##               RPN1              PNPT1              UBXN4               PAX3 
    ##              0.475              0.475              0.475              0.475 
    ##              WBP1L               CREM               NRP2                FAU 
    ##              0.475              0.475              0.475              0.476 
    ##              CHIC1              GREB1             GOLT1B              INTS1 
    ##              0.476              0.476              0.476              0.476 
    ##           TNFRSF1B             SLC6A2            PPP1R9A               PBX3 
    ##              0.476              0.476              0.476              0.476 
    ##               E4F1              ACSL1              PLCL1         CSGALNACT2 
    ##              0.476              0.476              0.476              0.476 
    ##             FRMPD3             KIF18A           KIAA1462              DOCK8 
    ##              0.476              0.477              0.477              0.477 
    ##              MYLK2            TMEM248               ATF1                DDN 
    ##              0.477              0.477              0.477              0.477 
    ##              FOXA2              GLRX3            EPB41L3              TMTC2 
    ##              0.477              0.477              0.477              0.477 
    ##           KIAA0513              FSD1L              VDAC2             OTUD7B 
    ##              0.477              0.477              0.477              0.478 
    ##             ABCA12             MAP7D2              VPS4B          C1GALT1C1 
    ##              0.478              0.478              0.478              0.478 
    ##              FOXN2              VEGFC             LPCAT1           PRICKLE3 
    ##              0.478              0.478              0.478              0.478 
    ##             PRKAA1              ABCF2              DNAH1                PXK 
    ##              0.478              0.478              0.478              0.479 
    ##              EXOC1              H2AFY              FOXO4               CHKA 
    ##              0.479              0.479              0.479              0.479 
    ##              PITX3              RAP2A              FEM1C              DMRT1 
    ##              0.479              0.479              0.479              0.479 
    ##            FAM155B            MAP3K11           SLC16A12               ATE1 
    ##              0.479              0.479              0.479              0.479 
    ##             PRKACB             ZDHHC9              PARP1             FBXO18 
    ##              0.479              0.479              0.479              0.479 
    ##               KAT8            ZNF804A              PRR19               RND1 
    ##              0.479              0.479              0.480              0.480 
    ##             DNAJC5              RPL27              NPTX1              PTHLH 
    ##              0.480              0.480              0.480              0.480 
    ##              MMP15              ACTC1              CREB5             TSG101 
    ##              0.480              0.480              0.480              0.481 
    ##              USP20               MSL2             OSBPL2              TNKS2 
    ##              0.481              0.481              0.481              0.481 
    ##              NEDD8            PLEKHO1              RPL31              CEP57 
    ##              0.481              0.481              0.481              0.481 
    ##               SMC4             FAM76B               ING5              ABCC9 
    ##              0.482              0.482              0.482              0.482 
    ##              STX16               BLMH              TDRD7            RPS6KC1 
    ##              0.482              0.482              0.482              0.482 
    ##               DTX4               DNM3              HMCN1                GRN 
    ##              0.482              0.483              0.483              0.483 
    ##            PCDHGC3              SGSM1               ING2           C17orf75 
    ##              0.483              0.483              0.483              0.483 
    ##              PDGFA              ZBTB5            SLC39A9              RPS26 
    ##              0.483              0.483              0.483              0.483 
    ##              TMOD1             MARCH9             SEMA5A               NPR1 
    ##              0.483              0.483              0.483              0.483 
    ##               HECA             GIGYF1             VPS13A             VPS26A 
    ##              0.484              0.484              0.484              0.484 
    ##              RBM46            SLC20A1               CDK2             ZNF248 
    ##              0.484              0.484              0.485              0.485 
    ##            FAM126A             B3GNT5                HK2                TAZ 
    ##              0.485              0.485              0.485              0.485 
    ##              CNTN5             PRKAB2               TPT1              ASF1A 
    ##              0.485              0.485              0.485              0.485 
    ##             PAPPA2              PFDN1               GRM2             STARD7 
    ##              0.485              0.485              0.485              0.485 
    ##              RNF34              PDE9A             TARBP2              ABCA1 
    ##              0.485              0.485              0.486              0.486 
    ##              RBMS1             TADA2B            MSANTD1              REPS2 
    ##              0.486              0.486              0.486              0.486 
    ##               NPC1               SMC5            DENND1B               TLR8 
    ##              0.486              0.486              0.486              0.486 
    ##               GAS6              MROH1              COPS8             STK32C 
    ##              0.486              0.486              0.486              0.486 
    ##               GJB1              RPL30                TNF               UCK2 
    ##              0.486              0.486              0.487              0.487 
    ##             RNF121             KLHL11             TBC1D9            AFAP1L2 
    ##              0.487              0.487              0.487              0.487 
    ##            WBSCR17             ZNF624            FAM91A1             COPS7B 
    ##              0.487              0.488              0.488              0.488 
    ##             HAPLN1            GPBP1L1           SH3PXD2B               ENC1 
    ##              0.488              0.488              0.488              0.488 
    ##             GALNT2              NETO2              STT3A           PPP1R13L 
    ##              0.488              0.488              0.488              0.489 
    ##           KIAA1731               CNN3             PIK3R2              PPP4C 
    ##              0.489              0.489              0.489              0.489 
    ##           VKORC1L1              PSMA4              MYCBP              THAP1 
    ##              0.489              0.489              0.489              0.489 
    ##             FAM20C               SMC6             NCAPD3             MAPRE3 
    ##              0.489              0.489              0.489              0.490 
    ##              ESRRA              ALAS1               RNF4               PDP1 
    ##              0.490              0.490              0.490              0.490 
    ##              NPTXR               SYT9               E2F4              PEX5L 
    ##              0.490              0.490              0.490              0.490 
    ##               NFYC              UGGT1             RIMKLB             MCM3AP 
    ##              0.491              0.491              0.491              0.491 
    ##             LRRTM1              ARMC5             CACNB4              NPRL3 
    ##              0.491              0.491              0.491              0.491 
    ##              ITSN2               LTBR             NUDCD3              SPTA1 
    ##              0.491              0.491              0.491              0.491 
    ##                ARC              RENBP            CXorf21              NUDT4 
    ##              0.491              0.491              0.492              0.492 
    ##              MAGOH              ATXN3               WARS              UBE2M 
    ##              0.492              0.492              0.492              0.492 
    ##              ABCC4               ARSE              UNC5C            KHDRBS3 
    ##              0.492              0.492              0.492              0.493 
    ##            DENND6A              MCTP1              PSMC4               TBCD 
    ##              0.493              0.493              0.493              0.493 
    ##             SPECC1              LRRK1              GRINA             ZNF207 
    ##              0.493              0.493              0.493              0.493 
    ##             CRYBG3            CXorf30            C7orf43             CNKSR3 
    ##              0.493              0.493              0.494              0.494 
    ##            DNAJB12              WNT5B            GLTSCR2             HS3ST4 
    ##              0.494              0.494              0.494              0.494 
    ##             GAREML                UST      CTD-2116N17.1             FAM69A 
    ##              0.494              0.494              0.494              0.494 
    ##              SNAI2               CARS             CCNYL1              SRP72 
    ##              0.494              0.494              0.494              0.494 
    ##              TMCC2               YAF2              SNRPB             UBE2J1 
    ##              0.494              0.494              0.494              0.495 
    ##               DHX8             B3GNT7              SRPK3            ZBTB8OS 
    ##              0.495              0.495              0.495              0.495 
    ##             ATP1A2                GGH             PWWP2B            ST3GAL1 
    ##              0.495              0.495              0.495              0.495 
    ##             EIF4E2               NPR2             KLHL14              FRMD6 
    ##              0.495              0.495              0.495              0.495 
    ##            TMEM115              VPS52            AKIRIN1              NETO1 
    ##              0.495              0.495              0.496              0.496 
    ##               DAXX              ITM2C                DTL               MID2 
    ##              0.496              0.496              0.496              0.496 
    ##              GRB14              NUAK1              IL2RB              SEPT2 
    ##              0.496              0.496              0.496              0.496 
    ##               PKD2     RPL17-C18orf32            RHOBTB2               PAX1 
    ##              0.497              0.497              0.497              0.498 
    ##              DUSP8               TCF3             ZNF768               ZHX3 
    ##              0.498              0.498              0.498              0.498 
    ##              TENC1             DCAF11             PGM2L1              CCNT1 
    ##              0.498              0.498              0.498              0.498 
    ##               MSH6             TGFBR3               ANKH             ZNF250 
    ##              0.498              0.499              0.499              0.499 
    ##              KRIT1             UBE2L3              PQBP1             TM4SF2 
    ##              0.499              0.499              0.499              0.499 
    ##           KIAA1467              SIAH1              YWHAQ              TIMP2 
    ##              0.499              0.499              0.499              0.499 
    ##            GPR137C              APLNR               PWP2              PCSK5 
    ##              0.499              0.499              0.499              0.499 
    ##           TBC1D22B            MACROD1              SNIP1            RNF113B 
    ##              0.499              0.499              0.499              0.499 
    ##            CCDC117               TEX2             ZNF689             ELMOD2 
    ##              0.500              0.500              0.500              0.500 
    ##              LPPR1             STARD8             ZNF358                GAN 
    ##              0.500              0.500              0.500              0.500 
    ##              RUFY2              THAP4           ARHGEF40            CBFA2T3 
    ##              0.500              0.500              0.500              0.500 
    ##               LAG3               MNX1                VGF             GARNL3 
    ##              0.500              0.500              0.500              0.501 
    ##              MED23             ZNF566             ZNF335             ADAM19 
    ##              0.501              0.501              0.501              0.501 
    ##              CDC42              ABCA3               GET4              AP2S1 
    ##              0.501              0.501              0.502              0.502 
    ##              RNF44             DUSP16               ATL2              CD79A 
    ##              0.502              0.502              0.502              0.502 
    ##               PCNP              VPS36               NOS3               TMF1 
    ##              0.502              0.502              0.502              0.502 
    ##            GRAMD1A               IRS1               PJA2             AGPAT6 
    ##              0.502              0.502              0.502              0.502 
    ##            C1QTNF5             ANTXR2              DUSP7             GALNT4 
    ##              0.503              0.503              0.503              0.503 
    ##              LZTS2               RFFL              SIRT1          KIAA1549L 
    ##              0.503              0.503              0.503              0.503 
    ##                CBL              YWHAH               GART              MCF2L 
    ##              0.503              0.504              0.504              0.504 
    ##            UNC119B             TMEM72          GABARAPL1              AKTIP 
    ##              0.504              0.504              0.504              0.504 
    ##            CACNA1H             COL8A1              USP21             TRIP11 
    ##              0.504              0.504              0.504              0.504 
    ##              CDCA3             TECPR1              FKBP4            COL15A1 
    ##              0.504              0.504              0.504              0.504 
    ##              DNM1L               CBX7              HNF4A               CSF2 
    ##              0.505              0.505              0.505              0.505 
    ##              HYDIN              MBNL3               RGS1              MORC1 
    ##              0.505              0.505              0.505              0.505 
    ##             PNPLA8               EPRS           TBC1D10C               IL27 
    ##              0.505              0.505              0.505              0.505 
    ##              ZNF81             SH3RF1               NXT2              PDE3B 
    ##              0.505              0.505              0.505              0.506 
    ##             RHBDF2               CHPF              OLIG3             NAP1L2 
    ##              0.506              0.506              0.506              0.506 
    ##               PEG3              NFXL1              IL1R1                PXN 
    ##              0.506              0.506              0.506              0.506 
    ##              SPDYA               MYPN              MYO1F             KCNAB1 
    ##              0.506              0.506              0.507              0.507 
    ##             SLC9A2              SEH1L             FYTTD1           PPP1R15B 
    ##              0.507              0.507              0.507              0.507 
    ##               FAT2              SMAD7               PTEN             SPRYD3 
    ##              0.507              0.507              0.507              0.507 
    ##              RPS15              AADAT              SSU72             ZNF574 
    ##              0.508              0.508              0.508              0.508 
    ##              SETD3              PALLD              ABCB4               NFYA 
    ##              0.508              0.508              0.508              0.509 
    ##               EDF1               FUT8              PDIA4              HMGA1 
    ##              0.509              0.509              0.509              0.509 
    ##            FAM134A               TPM4             KLHL18               FNTB 
    ##              0.509              0.509              0.509              0.509 
    ##              RBM44             ZNF398             MAP2K5              FARP1 
    ##              0.509              0.509              0.509              0.509 
    ##           RABGAP1L            SLC30A3                MAL              PELI1 
    ##              0.509              0.509              0.509              0.509 
    ##              BTBD7             RNF138              RAB9A              ZNF12 
    ##              0.510              0.510              0.510              0.510 
    ##              DHDDS               FGF8            UNC93B1                HCK 
    ##              0.510              0.510              0.510              0.510 
    ##            CTNNBL1              PDE5A              KCNJ6              LPIN2 
    ##              0.510              0.510              0.510              0.510 
    ##               MXD4              ADAD1             KLHL20            SLC18A2 
    ##              0.510              0.511              0.511              0.511 
    ##               VARS               EIF1              CALB1                TFG 
    ##              0.511              0.511              0.511              0.511 
    ##               GPC6            ZC3H12A              CERS6            ALDH1A3 
    ##              0.511              0.511              0.511              0.511 
    ##              SIPA1            ZSCAN20           KIAA1279               DMPK 
    ##              0.511              0.511              0.511              0.512 
    ##             KIF13B             FAM46D               NAPA           SLC25A11 
    ##              0.512              0.512              0.512              0.512 
    ##            PLA2G15               VASN               SOX8            RASGRP3 
    ##              0.512              0.512              0.512              0.512 
    ##            EMILIN1               NDNF        CTB-102L5.4             ZNF112 
    ##              0.512              0.512              0.512              0.512 
    ##              FSIP2              GPR85              CHTF8               STYX 
    ##              0.512              0.513              0.513              0.513 
    ##           TMEM178A              FOXI1              RAB4B             RNF208 
    ##              0.513              0.513              0.513              0.513 
    ##           ATP6V1E1           TMEM184B              NELFA             NDFIP2 
    ##              0.514              0.514              0.514              0.514 
    ##             GPR114            NEUROD1            SEC23IP             SLC7A5 
    ##              0.514              0.514              0.514              0.514 
    ##              CNTN3              OSBP2               FGF7     MSANTD3-TMEFF1 
    ##              0.514              0.514              0.514              0.514 
    ##              HINFP              LARP6               NEFM              PCSK7 
    ##              0.514              0.514              0.514              0.514 
    ##            RPS6KA1              CNNM1               SHOX              USP36 
    ##              0.515              0.515              0.515              0.515 
    ##            EIF2AK2              PDE3A            ARHGEF3               DMWD 
    ##              0.515              0.515              0.515              0.515 
    ##             GOLGA4               VAX1               NFE2             ZNF366 
    ##              0.515              0.515              0.515              0.515 
    ##             ZNF784             ZNF791            GAL3ST1              RNF25 
    ##              0.515              0.515              0.515              0.515 
    ##              TNPO3               NCS1               BRD8              TAF6L 
    ##              0.515              0.516              0.516              0.516 
    ##              GPR98               ALS2              MAPK9            CNTNAP4 
    ##              0.516              0.516              0.516              0.516 
    ##             DNAJC3            FAM179B              UBAP2             FAM57B 
    ##              0.516              0.516              0.516              0.516 
    ##              CXADR             SEMA3C            FAM134C               CPS1 
    ##              0.516              0.516              0.517              0.517 
    ##               MTA3               MAEA               CDV3           TMEM184C 
    ##              0.517              0.517              0.517              0.517 
    ##              TTC13               CHRD               RHOC               PALM 
    ##              0.517              0.517              0.517              0.517 
    ##               HES7              TUBG1       RP11-73M18.2               ZHX2 
    ##              0.517              0.517              0.518              0.518 
    ##            ADIPOR2            DNAJB11            PAK1IP1              CPNE8 
    ##              0.518              0.518              0.518              0.518 
    ##              PEX14              SNX25             ZNF853                AIP 
    ##              0.518              0.518              0.518              0.518 
    ##             CXCL13               ZIC5             TGFBR2              EPB41 
    ##              0.519              0.519              0.519              0.519 
    ##              TRAK1              G3BP1               ALX4             TSPAN9 
    ##              0.519              0.519              0.519              0.519 
    ##           ALDH18A1              BIRC3              CNTN4             GXYLT1 
    ##              0.519              0.519              0.519              0.519 
    ##               TEFM              ETAA1               RBM3              ESRRB 
    ##              0.519              0.519              0.519              0.519 
    ##              ITGB3              NLRP3               PARN             SCARB1 
    ##              0.519              0.519              0.519              0.519 
    ##               BDNF              HTR2A             CDKAL1            SLC44A5 
    ##              0.520              0.520              0.520              0.520 
    ##              FLRT3              FBXL2              NADK2             RNF170 
    ##              0.520              0.520              0.520              0.520 
    ##                TTL              MAPK4              RPL38              TDRKH 
    ##              0.520              0.520              0.520              0.520 
    ##             ZNF668             RAB39B         AC104534.3              TCF21 
    ##              0.520              0.520              0.520              0.520 
    ##              FBLN2                PAM            RACGAP1             ZNF646 
    ##              0.520              0.520              0.520              0.521 
    ##              CPLX1            CACNA1S              ZFP90              NCOA7 
    ##              0.521              0.521              0.521              0.521 
    ##              NRBP2              EPHA5              CKAP4               VAV2 
    ##              0.521              0.521              0.521              0.521 
    ##             PKMYT1              LIMK2            TMEM181              NAA10 
    ##              0.521              0.521              0.521              0.522 
    ##             COL4A3              NALCN              FEM1B                CFB 
    ##              0.522              0.522              0.522              0.522 
    ##               BRD9               MBD2              AKAP5            OSBPL1A 
    ##              0.522              0.522              0.522              0.522 
    ##               SIX1             CHRDL1               RHOV              RNF26 
    ##              0.522              0.522              0.522              0.522 
    ##               SCAP             DHCR24            POLDIP3              MYLPF 
    ##              0.522              0.523              0.523              0.523 
    ##            PRKAR2A               MSTN               PLD5              KCTD9 
    ##              0.523              0.523              0.523              0.523 
    ##               SSH1               PXDN               RPA2             SCAMP5 
    ##              0.523              0.523              0.523              0.523 
    ##             FBXO34              HACE1             INPP5F           SLC25A19 
    ##              0.524              0.524              0.524              0.524 
    ##              RGS10             ZBTB7B             KLHL12             PPP2CB 
    ##              0.524              0.524              0.524              0.524 
    ##               RBM7                SLA               STK4              GFRA1 
    ##              0.524              0.524              0.524              0.524 
    ##              GKAP1              NHSL2             PHLPP2            AKAP17A 
    ##              0.524              0.524              0.524              0.524 
    ##              DHX29               ZW10             CSRNP1             TNFSF8 
    ##              0.525              0.525              0.525              0.525 
    ##               CTSS           MPHOSPH8             IGSF21              SRRM3 
    ##              0.525              0.525              0.525              0.525 
    ##               JAK3              KCNA6             FAM49A              ARPC5 
    ##              0.525              0.526              0.526              0.526 
    ##             SPPL2A               LARS               RNF6              FARSA 
    ##              0.526              0.526              0.526              0.526 
    ##             RAB33A            SLC11A2              FRMD8             TDRD15 
    ##              0.526              0.526              0.526              0.527 
    ##             SH3RF3            ANKRD44              GATA4               GYS1 
    ##              0.527              0.527              0.527              0.527 
    ##              AKAP2               BTF3            ST8SIA1              KDM1B 
    ##              0.527              0.527              0.527              0.527 
    ##        PALM2-AKAP2             SLC4A3             POM121              EPHA3 
    ##              0.527              0.527              0.527              0.527 
    ##             GOLGA5              POLD1            CCDC88B              CYHR1 
    ##              0.527              0.527              0.528              0.528 
    ##              PDE1A            SHROOM2              GPR88              CTBP2 
    ##              0.528              0.528              0.528              0.528 
    ##              RPL24            TBC1D19              AKAP3             ZNF622 
    ##              0.528              0.528              0.528              0.528 
    ##             CORO2B              PAGE1               MYLK                GK5 
    ##              0.528              0.528              0.528              0.528 
    ##            SLC23A2                NGF              TSNAX           SLC25A42 
    ##              0.528              0.528              0.529              0.529 
    ##        IQCJ-SCHIP1             SMNDC1             CDKN1C              PARK7 
    ##              0.529              0.529              0.529              0.529 
    ##              CTDP1             RHBDD2              SPCS3              ACBD5 
    ##              0.529              0.529              0.529              0.529 
    ##              ITGA6              SRSF9               TECR              DUSP9 
    ##              0.529              0.529              0.529              0.529 
    ##             TRIM13               DRG1              CEP85     C7orf55-LUC7L2 
    ##              0.530              0.530              0.530              0.530 
    ##              FBXO8             SLC5A6              ZFP36           KIAA1522 
    ##              0.530              0.530              0.530              0.530 
    ##              LLGL1             LRPPRC               NFS1               GAS1 
    ##              0.530              0.530              0.530              0.530 
    ##       POC1B-GALNT4           ADAMTSL1             ADAM11             MUM1L1 
    ##              0.530              0.531              0.531              0.531 
    ##             METAP1              NSUN2            PLEKHA7           KIAA1211 
    ##              0.531              0.531              0.531              0.531 
    ##              NRROS              ERP29             TRMT12             RHBDL1 
    ##              0.531              0.531              0.531              0.531 
    ##              KCNK5             KLHDC2              STK35               BZW2 
    ##              0.531              0.531              0.532              0.532 
    ##             CGGBP1           KIAA0247               MTPN             FAM81A 
    ##              0.532              0.532              0.532              0.532 
    ##              PDE8A              WDR11              KCNC2              EIF3E 
    ##              0.532              0.532              0.533              0.533 
    ##               KLC4               RTL1             NFKBID              STK25 
    ##              0.533              0.533              0.533              0.533 
    ##            RANGAP1             SEMA4B              TIMP3               AVL9 
    ##              0.533              0.533              0.533              0.533 
    ##              ITGA3               SSR3              SSTR1              VPS29 
    ##              0.534              0.534              0.534              0.534 
    ##               GATM              WDR81             AGPAT4              PPM1L 
    ##              0.534              0.534              0.534              0.535 
    ##              FBXL7            SLC35F3            COL14A1            COL18A1 
    ##              0.535              0.535              0.535              0.535 
    ##             TSPAN7              HVCN1              HYAL2              WIPF2 
    ##              0.535              0.535              0.535              0.535 
    ##           NIPSNAP1              ADCY8           RAB3GAP1               TAC1 
    ##              0.535              0.535              0.535              0.535 
    ##          KIAA1324L               CHP1             RAD23A              UBTD1 
    ##              0.535              0.535              0.536              0.536 
    ##              ITM2B                 F2             METRNL                DCN 
    ##              0.536              0.536              0.536              0.536 
    ##               CLPB                FGG              MYLIP               SACS 
    ##              0.536              0.536              0.536              0.536 
    ##              RAB5A            ZC2HC1A              RTEL1              TTC18 
    ##              0.536              0.537              0.537              0.537 
    ##              RAP1A            PLEKHH1             UNC13C              UBE2A 
    ##              0.537              0.537              0.537              0.537 
    ##           SLC16A10           ADAMTS19            SLC23A1              NUTM1 
    ##              0.537              0.537              0.537              0.537 
    ##              FGF12            FAM188A               SGCE                MAF 
    ##              0.537              0.537              0.537              0.537 
    ##                NES             KIF18B             OSBPL3               UAP1 
    ##              0.537              0.538              0.538              0.538 
    ##             MFSD2A              TOP3B              SRPX2               FGD6 
    ##              0.538              0.538              0.538              0.538 
    ##             PTGFRN               EHD2               MST4              VPS11 
    ##              0.538              0.538              0.538              0.538 
    ##              PAGR1              PAGR1             FKBP1A            PRKAR1B 
    ##              0.538              0.538              0.538              0.538 
    ##              AMER3            LDLRAD4              VGLL1             CAMK2D 
    ##              0.539              0.539              0.539              0.539 
    ##             RBM12B               DIO2              FGF20              KDM4C 
    ##              0.539              0.539              0.539              0.539 
    ##             REPIN1               TLX1              TSSK6             ERLIN1 
    ##              0.539              0.539              0.539              0.539 
    ##              ACBD3               MDC1             SH3BP1              CPSF1 
    ##              0.539              0.539              0.540              0.540 
    ##            FAM115A            ARFGAP2            TWISTNB                OTP 
    ##              0.540              0.540              0.540              0.540 
    ##               PIN1               SHC3             GOLGA1            SLC7A14 
    ##              0.540              0.540              0.541              0.541 
    ##             B3GNT3              OXCT1              TRPC7              RALYL 
    ##              0.541              0.541              0.541              0.541 
    ##             ASPHD2               MUM1               NCK2            PLEKHM2 
    ##              0.541              0.541              0.541              0.541 
    ##             ZBTB11            SHROOM1             RNF112              WDR75 
    ##              0.541              0.541              0.541              0.541 
    ##              USP39               CPOX              ICAM5               NFIC 
    ##              0.541              0.542              0.542              0.542 
    ##            PIP5K1A           TP53INP2             DNAJA1               ODC1 
    ##              0.542              0.542              0.542              0.543 
    ##               ERC1             BARHL2             KIFAP3                SSB 
    ##              0.543              0.543              0.543              0.543 
    ##             CYB561              SEPT3              EIF5A              BEND2 
    ##              0.543              0.543              0.543              0.544 
    ##            FAM178A      RP11-565P22.6                A2M              IFT46 
    ##              0.544              0.544              0.544              0.544 
    ##            MORF4L2            NEUROG2                 FH              TPST1 
    ##              0.544              0.544              0.544              0.544 
    ##            FAM222B               UBL3              NAA50               MLEC 
    ##              0.544              0.544              0.545              0.545 
    ##             CX3CL1              PHTF2                DLD            C12orf5 
    ##              0.545              0.545              0.545              0.545 
    ##                WLS                MVK              UBE3B              SSFA2 
    ##              0.545              0.545              0.545              0.545 
    ##            GDAP1L1             KLHL32              LMAN1              SURF4 
    ##              0.546              0.546              0.546              0.546 
    ##               GMNC             CRNKL1             UQCRC1               PDK3 
    ##              0.546              0.546              0.546              0.546 
    ##             LRRTM4               CAP1              PHKA1            CNTNAP2 
    ##              0.546              0.547              0.547              0.547 
    ##             LMBRD2              LSM11              CPNE5              FOXA3 
    ##              0.547              0.547              0.547              0.547 
    ##               EMC3              S1PR1              CENPC              CNOT6 
    ##              0.547              0.547              0.547              0.548 
    ##              DMAP1          TNFRSF11B              NUP50              HOXD3 
    ##              0.548              0.548              0.548              0.548 
    ##            LAMTOR3              PPIL2              EIF4H           C12orf66 
    ##              0.548              0.548              0.548              0.548 
    ##             ITGA11             PTCHD2               APOO             SCAMP2 
    ##              0.548              0.548              0.548              0.548 
    ##                HDC              SEPT8                SHF               SYT4 
    ##              0.548              0.548              0.548              0.549 
    ##            MAB21L2              TXNIP               RIN2            ZNF385A 
    ##              0.549              0.549              0.549              0.549 
    ##               NOL9              CCNB1              REEP2            SLC10A3 
    ##              0.549              0.549              0.549              0.550 
    ##               NAPG             LSM14A              IDH3A             AMIGO2 
    ##              0.550              0.550              0.550              0.550 
    ##              CXCR5               DND1             RNF169             IGDCC4 
    ##              0.550              0.550              0.550              0.550 
    ##               CHFR             CCSER1            SULT2B1              INTS7 
    ##              0.550              0.550              0.550              0.550 
    ##             STRIP1              ACSL3      BCL2L2-PABPN1              USP54 
    ##              0.550              0.550              0.550              0.550 
    ##              TRHDE            CCDC120             LRRC59               GBF1 
    ##              0.550              0.551              0.551              0.551 
    ##              H3F3B            ZKSCAN2            ADORA2A               GNL2 
    ##              0.551              0.551              0.551              0.551 
    ##               DOK6             WDR45B             RNF216            ALOX12B 
    ##              0.551              0.551              0.551              0.552 
    ##           ANKRD33B             SNCAIP              ADCY7             PDLIM5 
    ##              0.552              0.552              0.552              0.552 
    ##             ZNF655               MCM4              RLTPR      CTD-2132N18.3 
    ##              0.552              0.552              0.552              0.552 
    ##               PDHB              NOP58              ARL4A               DBNL 
    ##              0.553              0.553              0.553              0.553 
    ##            ARHGEF4            C14orf2            HNRNPDL               SSR1 
    ##              0.553              0.553              0.553              0.553 
    ##              SMAP2             ZNF184              KIF2C              SGPL1 
    ##              0.553              0.553              0.554              0.554 
    ##               WWC2             SEC11C           SLC39A14              RIPK1 
    ##              0.554              0.554              0.554              0.554 
    ##               USP2              SNRPC              SNRPG              CPSF4 
    ##              0.555              0.555              0.555              0.555 
    ##             OSBPL5                HDX             PCMTD1               POLG 
    ##              0.555              0.555              0.555              0.555 
    ##           TMEM229A             TOM1L2              EPHA2            NDUFB11 
    ##              0.555              0.555              0.555              0.556 
    ##               FA2H              ELMO2             CREBL2            DPY19L1 
    ##              0.556              0.556              0.556              0.556 
    ##             NFKBIB              HSPE1             UBALD1             HIF1AN 
    ##              0.556              0.556              0.557              0.557 
    ##              PSMD7             CD40LG             POU6F2            LAMTOR1 
    ##              0.557              0.557              0.557              0.557 
    ##              NANOG            SEC61A2               PAX7                AK5 
    ##              0.557              0.558              0.558              0.558 
    ##               PEX3              FARSB               ROR2               DIO1 
    ##              0.558              0.558              0.558              0.558 
    ##             FBXO47      CTD-3193O13.9             ABI3BP              TRPM7 
    ##              0.558              0.558              0.558              0.558 
    ##               TGFA             NCAPH2              PRKRA                SP7 
    ##              0.558              0.558              0.558              0.558 
    ##            SLC44A2              NUP93              AP1AR             TTC39A 
    ##              0.559              0.559              0.559              0.559 
    ##             CHRNB4             RPS15A            C10orf2            SLC35A2 
    ##              0.559              0.559              0.559              0.559 
    ##               VAPB               NID2                LPP              YIPF3 
    ##              0.559              0.559              0.560              0.560 
    ##             NDUFS2              ZC3H8             SNRPD2            PCDH11X 
    ##              0.560              0.560              0.560              0.560 
    ##              PRRT2             ATP8B1              ESRRG              CD200 
    ##              0.560              0.560              0.560              0.560 
    ##              MBNL2               SUN1             BCAP29            PIP4K2A 
    ##              0.560              0.560              0.560              0.560 
    ##              PTPN5            ANKRD40               E2F2                TUB 
    ##              0.560              0.561              0.561              0.561 
    ##               UBA6              ARAP2              TNIP3               GFI1 
    ##              0.561              0.561              0.561              0.561 
    ##               TSR2              WDR52              KCND1              RBM38 
    ##              0.561              0.561              0.561              0.561 
    ##              WDR19           KIAA1045              HNF4G             EFTUD1 
    ##              0.561              0.561              0.561              0.561 
    ##               UPRT            BHLHE22                 F5             SLC5A7 
    ##              0.561              0.562              0.562              0.562 
    ##               MCM5             FERMT3             COL6A3               CDK5 
    ##              0.562              0.562              0.562              0.562 
    ##               EML3             NYNRIN               PRTG              NELFE 
    ##              0.562              0.562              0.562              0.562 
    ##               HEPH              DNAH2            SLC37A1             CHST11 
    ##              0.563              0.563              0.563              0.563 
    ##               PAK1                TKT              TMTC3              DCLK3 
    ##              0.563              0.563              0.563              0.563 
    ##              NXF2B               PKP1              RPL32            TMEM198 
    ##              0.563              0.563              0.563              0.563 
    ##               GBX2               TFPI               TYK2             ZNF205 
    ##              0.563              0.564              0.564              0.564 
    ##            FAM199X              CCPG1             SLC5A5              TAPT1 
    ##              0.564              0.564              0.564              0.564 
    ##               CD74              CPNE2              MTCH2             SNRPD1 
    ##              0.564              0.564              0.564              0.564 
    ##              RBM8A                BMX             IGLON5            PLEKHG5 
    ##              0.564              0.564              0.565              0.565 
    ##               CLPP              MMGT1              MTMR9               GPC2 
    ##              0.565              0.565              0.565              0.565 
    ##              ESYT2            PLEKHM3              CHSY3             MARCH4 
    ##              0.566              0.566              0.566              0.566 
    ##              HMGB3           GPATCH2L              CDH23                DEK 
    ##              0.566              0.566              0.567              0.567 
    ##            PSTPIP2               ABAT              OVOL2             TMEM47 
    ##              0.567              0.567              0.567              0.567 
    ##              CNPY3               CHN1              RPL34                MFF 
    ##              0.567              0.567              0.567              0.567 
    ##            NEURL1B              SYTL4               IL15           ARHGAP27 
    ##              0.567              0.567              0.567              0.567 
    ##               DSEL             PTPRN2               MAP6              PTPN3 
    ##              0.568              0.568              0.568              0.568 
    ##            POU2AF1       RP11-101E3.5               CHGA              CLCC1 
    ##              0.568              0.568              0.568              0.568 
    ##               HES1              MYRIP              ASCC3                 KL 
    ##              0.568              0.568              0.568              0.568 
    ##             PABPN1              ZFP28               AMPH               SOX2 
    ##              0.569              0.569              0.569              0.569 
    ##           TRAPPC11              BTBD3              TSSC1              ITPR3 
    ##              0.569              0.569              0.569              0.569 
    ##             SPOCK3               MAPT               LMO3               PLK4 
    ##              0.569              0.569              0.569              0.569 
    ##            SLC24A1               BMP7              SEC63             NUP107 
    ##              0.569              0.570              0.570              0.570 
    ##             UBE2E1               SNX8              VGLL2                GFY 
    ##              0.570              0.570              0.570              0.570 
    ##             RABEP2             FBRSL1               GRK6              SNRPF 
    ##              0.570              0.571              0.571              0.571 
    ##              BICD1              SMYD5            ST6GAL1               DTX2 
    ##              0.571              0.571              0.571              0.571 
    ##           ANKRD13A                SRM               DEF8             CRYBB2 
    ##              0.571              0.571              0.571              0.571 
    ##             GPR101              KDM5B              PROS1            EIF2AK4 
    ##              0.571              0.572              0.572              0.572 
    ##              THSD4              WFDC2              PVRL1              CEBPZ 
    ##              0.572              0.572              0.572              0.572 
    ##             CACHD1             SH3BP4            SLC24A2               HEG1 
    ##              0.572              0.572              0.572              0.573 
    ##              SASH3              ABCF3           C10orf76           C16orf62 
    ##              0.573              0.573              0.573              0.573 
    ##               PPOX              DNAH6             KCNMB2              MEX3B 
    ##              0.573              0.573              0.574              0.574 
    ##               IL34             LIN28B               RCL1              KCND2 
    ##              0.574              0.574              0.574              0.574 
    ##            SLC43A2             ZSWIM4              CD248              OPCML 
    ##              0.574              0.575              0.575              0.575 
    ##              KCNJ8             LZTFL1               MLH1           C15orf39 
    ##              0.575              0.575              0.575              0.575 
    ##              CECR6             CTAGE5               MYH7             GLCCI1 
    ##              0.575              0.575              0.575              0.575 
    ##              RRP12              NCAPH            DYNC2H1               RER1 
    ##              0.575              0.575              0.576              0.576 
    ##             SMIM14                BOC               BOLL               OLA1 
    ##              0.576              0.576              0.576              0.576 
    ##               MADD             BRMS1L              TUBE1              FOXE1 
    ##              0.576              0.576              0.576              0.576 
    ##             TTC39C              NSRP1              RGS14              PEA15 
    ##              0.576              0.576              0.576              0.576 
    ##               ANHX              PIANP            TGFB1I1              SPOPL 
    ##              0.576              0.576              0.576              0.577 
    ##              CD244               SGK1            ZNF518B               PCM1 
    ##              0.577              0.577              0.577              0.577 
    ##               ICA1             RPRD1A            FAM211A             COL6A2 
    ##              0.577              0.578              0.578              0.578 
    ##               PISD              SGMS2               RTCB               ASB9 
    ##              0.578              0.578              0.578              0.578 
    ##               KDSR              HOXB5           B4GALNT4            C2orf42 
    ##              0.578              0.578              0.578              0.578 
    ##              SNX17              NDRG4            FAM105B             PTGER3 
    ##              0.579              0.579              0.579              0.579 
    ##            SLC30A4             ZBTB22              NOP14              NOTUM 
    ##              0.579              0.579              0.580              0.580 
    ##               HHEX              STK24              NDRG1               ABL2 
    ##              0.580              0.580              0.580              0.580 
    ##                CKB              TECTA              CTCFL             TRIM23 
    ##              0.580              0.580              0.580              0.580 
    ##              BUB1B            NFKBIL1                ZFY              CNOT8 
    ##              0.580              0.580              0.580              0.580 
    ##              FOXK2            IL18RAP               TFRC             ARMCX2 
    ##              0.580              0.580              0.581              0.581 
    ##             PIEZO1              CHIC2              EEPD1             NECAB1 
    ##              0.581              0.581              0.581              0.581 
    ##              WNT9B               NUB1              RPAP3             D2HGDH 
    ##              0.581              0.581              0.581              0.581 
    ##             MBOAT7             NCAPD2             WRNIP1               MCAM 
    ##              0.581              0.582              0.582              0.582 
    ##              CRHBP               NAPB              FASTK               GPS2 
    ##              0.582              0.582              0.583              0.583 
    ##              CDCA7            CDK2AP2              AGGF1              ISLR2 
    ##              0.583              0.583              0.583              0.583 
    ##              CIITA            ATP13A2              BSDC1               YKT6 
    ##              0.583              0.584              0.584              0.584 
    ##               ELF1           SLC25A25            FAM102A             MINOS1 
    ##              0.584              0.584              0.584              0.584 
    ##           SLC22A17              CDH12               DDR1             SEC24A 
    ##              0.584              0.584              0.584              0.584 
    ##             SEC24D              PTPRN               ULK3             KLHL28 
    ##              0.584              0.584              0.584              0.585 
    ##             TSEN34             COL4A2               MCM2             PIWIL4 
    ##              0.585              0.585              0.585              0.585 
    ##              FHDC1               USF1           TMEM170B              PDE7A 
    ##              0.585              0.585              0.585              0.585 
    ##             NKX2-1             CABIN1             SH2D3C              SP140 
    ##              0.586              0.586              0.586              0.586 
    ##              DDX60              SIAH2              ENPP1              NDST3 
    ##              0.586              0.586              0.586              0.586 
    ##              FSTL3              RPS16             UNC13B              ACAP3 
    ##              0.586              0.586              0.586              0.586 
    ##             DNAJB5               SYN3              KCNC3               SEZ6 
    ##              0.587              0.587              0.587              0.587 
    ##               MTTP               TMX3               URB1              FLRT1 
    ##              0.587              0.587              0.587              0.587 
    ##             GALNT1              DDX43           RALGAPA2            PLA2G4A 
    ##              0.587              0.587              0.587              0.587 
    ##             DNAJC8                JUP               MLH3              ADCK5 
    ##              0.587              0.587              0.587              0.588 
    ##             MBTPS1               ULK2             GTF3C3            GALNT16 
    ##              0.588              0.588              0.588              0.588 
    ##              TRMT5            C7orf25             MYBPC3             TMEM9B 
    ##              0.588              0.588              0.588              0.588 
    ##              RAB4A               CDX2             NSFL1C               ARR3 
    ##              0.589              0.589              0.589              0.589 
    ##               CLP1              CPEB2            RBBP8NL              USP9Y 
    ##              0.589              0.589              0.589              0.589 
    ##               VAPA              FNDC1            SLC19A3           TBC1D10A 
    ##              0.589              0.589              0.589              0.589 
    ##            TMEM180             ZNF330               EGR1              EXOC2 
    ##              0.589              0.590              0.590              0.590 
    ##          C14orf166               CRY2               PEX5             DNAJA3 
    ##              0.590              0.590              0.590              0.590 
    ##               HBP1               FIZ1            ZC3HAV1               ESAM 
    ##              0.590              0.590              0.590              0.590 
    ##              SNX12             SOWAHA             ABLIM3             SETDB2 
    ##              0.590              0.590              0.590              0.590 
    ##              RPS24           C11orf73               NAB1              TRA2A 
    ##              0.591              0.591              0.591              0.591 
    ##             ORMDL3              FOXM1              ADRB2             LETMD1 
    ##              0.591              0.591              0.591              0.591 
    ##             SH3BP5              TACC1               WBP2               COG2 
    ##              0.591              0.591              0.591              0.591 
    ##             PHLDB2            HSD17B6             TRIM25               COIL 
    ##              0.591              0.591              0.591              0.592 
    ##            B4GALT2            DCUN1D4              SMOC1              EIF3H 
    ##              0.592              0.592              0.592              0.592 
    ##            ONECUT1              FLAD1             SLC8A3              GPM6B 
    ##              0.592              0.592              0.592              0.592 
    ##           C19orf26             PGRMC2              GPR26              MED19 
    ##              0.592              0.592              0.593              0.593 
    ##              EDNRB              EXOC6            PTPLAD1              LPAR3 
    ##              0.593              0.593              0.593              0.593 
    ##                VWF            TMEM245              BCL6B              FANCM 
    ##              0.593              0.593              0.593              0.593 
    ##             METAP2             STXBP3               NUMB               LBX1 
    ##              0.593              0.593              0.593              0.593 
    ##            C1orf94              CDH13              SPSB1            FAM192A 
    ##              0.593              0.594              0.594              0.594 
    ##            ZFYVE21               HPCA              RIPK4               RBP4 
    ##              0.594              0.594              0.594              0.594 
    ##             CLSTN2             ADAM12            HLA-DRA             SORBS1 
    ##              0.594              0.594              0.594              0.594 
    ##              TNNC1             PDCD11              ITGA2                SP5 
    ##              0.594              0.594              0.594              0.595 
    ##               IRF5      RP11-407N17.3             COL7A1              ASCL1 
    ##              0.595              0.595              0.595              0.595 
    ##               SIX6             PIK3CG               NGFR               LMO4 
    ##              0.595              0.595              0.595              0.595 
    ##               DPF3              NAA60              NAA60              SPNS1 
    ##              0.595              0.596              0.596              0.596 
    ##           RSPH10B2           ARHGAP12                C1S            SH3BP5L 
    ##              0.596              0.596              0.596              0.596 
    ##                DES             NUSAP1               BDP1              DISP1 
    ##              0.596              0.596              0.596              0.596 
    ##             VANGL1               CD72             RAB33B            UBASH3B 
    ##              0.596              0.596              0.597              0.597 
    ##             VWA5B2            FAM175B             TIMM44            ANKRD46 
    ##              0.597              0.597              0.597              0.597 
    ##                NDN             RANBP6             ZBTB42          TNFRSF11A 
    ##              0.597              0.597              0.597              0.597 
    ##              ROBO1              P2RY1              KIFC1              PDE12 
    ##              0.597              0.597              0.598              0.598 
    ##           TNFRSF25            ZFP36L1            TMEM110              CDC45 
    ##              0.598              0.598              0.598              0.598 
    ##           ARHGEF33             CRYBG3              CIRBP              SRPK1 
    ##              0.598              0.598              0.598              0.598 
    ##               NXF3            ARFGAP1             ATP8A2             UBLCP1 
    ##              0.599              0.599              0.599              0.599 
    ##              HOXC6               CLTB              TICRR              HIF3A 
    ##              0.599              0.599              0.599              0.599 
    ##               CD97          C20orf194                CRX              PGAP1 
    ##              0.599              0.600              0.600              0.600 
    ##              CDYL2             PRDM14              AREL1              SPIN3 
    ##              0.600              0.600              0.600              0.600 
    ##        COMMD3-BMI1               LHFP            TBC1D23               BUB1 
    ##              0.600              0.601              0.601              0.601 
    ##             CDK5R2             EFEMP2              ALDOC             PCED1A 
    ##              0.601              0.601              0.601              0.601 
    ##             ENTPD1              MYPOP               EHD3              GDPD5 
    ##              0.602              0.602              0.602              0.602 
    ##            TBC1D12               CNR1              NOLC1               FZD1 
    ##              0.602              0.602              0.602              0.602 
    ##               BRF1            CTTNBP2             CEP192              FAM3A 
    ##              0.602              0.602              0.602              0.602 
    ##            C7orf63            SLC6A15               DOK4            ZFYVE26 
    ##              0.602              0.602              0.602              0.602 
    ##           C20orf24              MON1B              DACT3              NELFB 
    ##              0.602              0.602              0.603              0.603 
    ##              ZNF35                LSS               WNT1             RPRD1B 
    ##              0.603              0.603              0.603              0.603 
    ##              TADA1               TNXB             ZNF143            ADAMTS3 
    ##              0.603              0.603              0.603              0.604 
    ##               MAEL              RSPO3            RHOBTB1              VASH1 
    ##              0.604              0.604              0.604              0.604 
    ##              FGFR3              SGSM2               AMN1             SEMA5B 
    ##              0.604              0.604              0.604              0.604 
    ##               TMC5               EGR2               STX3             MCOLN1 
    ##              0.604              0.604              0.604              0.604 
    ##                NSF           TMEM106B            SLC2A14             ACSBG1 
    ##              0.604              0.604              0.604              0.605 
    ##              HAUS4             PGGT1B             ADRA2C               SGTB 
    ##              0.605              0.605              0.605              0.605 
    ##              TULP1             CCDC93             PIWIL2             COL5A3 
    ##              0.605              0.605              0.605              0.605 
    ##             RHOXF1                GNS           TNFRSF21            TNFRSF9 
    ##              0.606              0.606              0.606              0.606 
    ##              IL23R              N4BP2              STAG3              NCAPG 
    ##              0.606              0.606              0.606              0.606 
    ##              RGS12               MIDN             STK17A              SEC13 
    ##              0.606              0.606              0.606              0.606 
    ##               FZD7             FKBP15               FLII               SMTN 
    ##              0.606              0.607              0.607              0.607 
    ##            METTL14              PRDM4              DHX38             FCER1G 
    ##              0.607              0.607              0.607              0.607 
    ##            ST3GAL4               ALX1                 IK              EFNB3 
    ##              0.607              0.607              0.607              0.607 
    ##              EIF3M             ARPP19               PSD2              DDX50 
    ##              0.607              0.608              0.608              0.608 
    ##              CLCN6             DOLPP1             GTF2E2            NCKIPSD 
    ##              0.608              0.608              0.608              0.608 
    ##             CSTF2T            SLC17A8            SERINC1              MEF2C 
    ##              0.608              0.608              0.608              0.608 
    ##           KIAA0754             ZC3H10              MAPK7              GLIS1 
    ##              0.609              0.609              0.609              0.609 
    ##               NEK6             NDUFA2             KATNB1             ZNF324 
    ##              0.609              0.609              0.609              0.610 
    ##            KHDRBS2               TFAM                IHH             PDZRN4 
    ##              0.610              0.610              0.610              0.610 
    ##               LYL1              VPS53           C19orf47             ZNF583 
    ##              0.610              0.611              0.611              0.611 
    ##           TMEM120B              COX10             SOHLH2     CCDC169-SOHLH2 
    ##              0.611              0.611              0.611              0.611 
    ##              WNT2B              DAPP1              GNAI3           SMARCAL1 
    ##              0.611              0.611              0.611              0.611 
    ##              SNX33             SH3D19             MTFR1L            PLEKHG2 
    ##              0.611              0.611              0.611              0.611 
    ##              PREX2             AKAP14             BTBD10               UXS1 
    ##              0.611              0.612              0.612              0.612 
    ##              ABCF1             TMEM64         ST6GALNAC5             ZNF101 
    ##              0.612              0.612              0.612              0.612 
    ##               NRAS             GPR124              GLUD1              KLF15 
    ##              0.612              0.612              0.612              0.613 
    ##             SLC1A5                ICK                TAT              CD164 
    ##              0.613              0.613              0.613              0.613 
    ##            SLC39A7              ASIC4              HLA-E           C11orf57 
    ##              0.613              0.614              0.614              0.614 
    ##              MAPK3              TGFB1             RNF166           SLC7A6OS 
    ##              0.614              0.614              0.614              0.614 
    ##               PAX9               LCAT              ITPR2             LONRF2 
    ##              0.614              0.614              0.615              0.615 
    ##               SPIC              OLFM2            SLC30A2              HOMEZ 
    ##              0.615              0.615              0.615              0.615 
    ##              KCNJ2              RNF20             TMEM89               GPC1 
    ##              0.615              0.615              0.615              0.616 
    ##              NOP10              PALM3               FECH            ANKRD10 
    ##              0.616              0.616              0.616              0.616 
    ##               DPP4               AHCY               PDCL             COL4A4 
    ##              0.616              0.616              0.616              0.616 
    ##              NXPH1                DCX              RCSD1               RGMA 
    ##              0.616              0.616              0.616              0.616 
    ##            TMEM39B             CLDND1              EPHB6     RP11-1035H13.3 
    ##              0.617              0.617              0.617              0.617 
    ##              OPRL1             LYSMD3             MAN2A1             POLR2H 
    ##              0.617              0.617              0.617              0.617 
    ##               AARS             LHFPL4              DACH2              ZFP92 
    ##              0.617              0.617              0.618              0.618 
    ##           CDC42SE2             PTPDC1            HSPA12B               RGS8 
    ##              0.618              0.618              0.618              0.618 
    ##              ICAM1              KCNJ9              RPL28              ZNF41 
    ##              0.618              0.618              0.618              0.618 
    ##             CCDC92              NOC4L              OTUD1               IER2 
    ##              0.619              0.619              0.619              0.619 
    ##             ZBTB7C             KLHL13             MAN2A2              TMED7 
    ##              0.619              0.619              0.619              0.619 
    ##                TBP               SOX7              WIPI2               MIB2 
    ##              0.619              0.619              0.619              0.620 
    ##               RCE1               FZD3          MAPK1IP1L             ZNF133 
    ##              0.620              0.620              0.620              0.620 
    ##             TUBA4A              JKAMP              CNIH1              TUSC2 
    ##              0.620              0.620              0.620              0.620 
    ##               PRC1              WDR91               BCAN              LAMA1 
    ##              0.621              0.621              0.621              0.621 
    ##              HADHA            PABPC1L              UNC5B              NR1D2 
    ##              0.621              0.621              0.621              0.621 
    ##             HSPA4L           FAM160A1            ST3GAL3              ESYT1 
    ##              0.621              0.621              0.621              0.622 
    ##            FOXRED2              FADS3              PRKD1              SMOC2 
    ##              0.622              0.622              0.622              0.622 
    ##              KRT10            NSMCE4A               PAF1            PLEKHO2 
    ##              0.622              0.622              0.622              0.622 
    ##              GPR50              HOOK1            SLC25A3               GOT2 
    ##              0.623              0.623              0.623              0.623 
    ##            TMEM123             COQ10B             ZNF304             ATP1B2 
    ##              0.623              0.623              0.624              0.624 
    ##              QSOX1               YARS               GJA1              TDRD3 
    ##              0.624              0.624              0.624              0.624 
    ##               RTN3             PPP4R4                CD2            SMARCD2 
    ##              0.624              0.624              0.624              0.625 
    ##               TLR9           C1orf172              ZBED6             ZNF563 
    ##              0.625              0.625              0.625              0.625 
    ##             RNF167            TXNDC16             CKAP2L              LRRN3 
    ##              0.626              0.626              0.626              0.626 
    ##             SPAG17             SERAC1            TMEM41B             EPS8L2 
    ##              0.626              0.626              0.626              0.627 
    ##           SERPINF2             FBXO24              KEAP1                MTR 
    ##              0.627              0.627              0.627              0.627 
    ##              UTP18              CERS3            COL16A1             SUCLA2 
    ##              0.627              0.627              0.627              0.627 
    ##               TOX2            IL22RA1               FDPS               DRP2 
    ##              0.627              0.627              0.628              0.628 
    ##               SV2C              MRVI1              ACTG2              CENPL 
    ##              0.628              0.628              0.628              0.628 
    ##              HOXA1             DOPEY2            MAB21L1                 HR 
    ##              0.628              0.628              0.628              0.628 
    ##             CITED2            RUNDC3B              NR1H4              SHFM1 
    ##              0.628              0.628              0.628              0.628 
    ##             ERGIC2               CTC1            SLC35B2               CHN2 
    ##              0.628              0.628              0.628              0.628 
    ##              NEDD4              NXPE3              MARK3               CDSN 
    ##              0.628              0.628              0.629              0.629 
    ##             LRRC49              GATA5             SYNPO2              RAB24 
    ##              0.629              0.629              0.629              0.629 
    ##               ZFAT             FAM13C               AMFR            TMEM165 
    ##              0.629              0.630              0.630              0.630 
    ##              EDIL3            CABLES1               MYH3             ZFAND6 
    ##              0.630              0.630              0.630              0.630 
    ##              EFNA1              FSTL5            DNAJC10             NKX3-2 
    ##              0.630              0.630              0.631              0.631 
    ##               LDHA             STARD9              MGAT3             SLC51A 
    ##              0.631              0.631              0.631              0.631 
    ##            ATP6V1H              SYCE2               SV2B              PDDC1 
    ##              0.631              0.631              0.631              0.631 
    ##              PPARA               PLAT              PSMD8              ARL4C 
    ##              0.632              0.632              0.632              0.632 
    ##            HEATR5B             IL17RA               MPP3             UTP14C 
    ##              0.632              0.632              0.632              0.632 
    ##              CAB39             SLC6A9              VSNL1               HPS4 
    ##              0.633              0.633              0.633              0.633 
    ##               LMO1           SLC25A37              PTPRJ             GOLGA7 
    ##              0.633              0.633              0.634              0.634 
    ##               VRK1              MKRN3                EN2             FCGR2B 
    ##              0.634              0.634              0.634              0.634 
    ##             CAMKK1             ACVR1C              KPNA5               MPC2 
    ##              0.634              0.634              0.634              0.634 
    ##               USP4             FAM53C              BRCA2               MEI1 
    ##              0.634              0.634              0.635              0.635 
    ##              ITGAM               EDN1              PPARG               PRR5 
    ##              0.635              0.635              0.635              0.635 
    ##               TJP2               UNCX             SEMA6C       RP5-1052I5.2 
    ##              0.635              0.635              0.635              0.635 
    ##             ZNF382             THSD7B             COL9A2             HCRTR2 
    ##              0.635              0.635              0.636              0.636 
    ##            SLC18A3              ALOX5           CALCOCO1           SLC25A23 
    ##              0.636              0.636              0.636              0.636 
    ##            SLC25A4              TMOD3              SPIN4                TEF 
    ##              0.636              0.636              0.636              0.636 
    ##                UBC               DOK1               SIL1              RCHY1 
    ##              0.636              0.636              0.636              0.636 
    ##            TINAGL1                BGN              UBE2S               HDGF 
    ##              0.636              0.637              0.637              0.637 
    ##             PPFIA4             SLC1A4              ATRIP            DENND2D 
    ##              0.637              0.637              0.637              0.637 
    ##              PRR5L              WDR12             VSTM2A            GALNT10 
    ##              0.637              0.637              0.637              0.637 
    ##               GDF7              SOCS2               TRHR              TPST2 
    ##              0.638              0.638              0.638              0.638 
    ##             ABCB11              NR2F6              AMPD2              ANPEP 
    ##              0.638              0.638              0.638              0.638 
    ##                MPZ               GRM8              KCNG3              PLOD2 
    ##              0.638              0.638              0.638              0.638 
    ##             SYNGR1              TNNT2              DHX16               CYC1 
    ##              0.638              0.638              0.638              0.639 
    ##          GTF2IRD2B            MYCBPAP     RP11-1407O15.2             TBCCD1 
    ##              0.639              0.639              0.639              0.639 
    ##               GCDH               LGR6             SH3GL1             LRRC8C 
    ##              0.639              0.639              0.639              0.639 
    ##               UGT8             CHST14              ICA1L              LRRK2 
    ##              0.639              0.640              0.640              0.640 
    ##             PRIMA1             ZDHHC3               POLE             ZBTB26 
    ##              0.640              0.640              0.640              0.640 
    ##               CNN1              NDST4              GFRA3              GGPS1 
    ##              0.640              0.640              0.640              0.640 
    ##            SLC6A11              RPS27              CYTH4              AEBP1 
    ##              0.640              0.640              0.640              0.640 
    ##               RALY             NDUFA9              LAMB1               MTO1 
    ##              0.640              0.640              0.640              0.640 
    ##               CRKL            GALNT18               ALG6              ZBTB6 
    ##              0.641              0.641              0.641              0.641 
    ##             ARFRP1             CHCHD3           C19orf66           KIAA0195 
    ##              0.641              0.641              0.641              0.641 
    ##             MESDC2              STK38                MOG              THOP1 
    ##              0.641              0.641              0.641              0.641 
    ##            STARD13               PKP4             PNPLA6              LAMP2 
    ##              0.641              0.641              0.641              0.641 
    ##              RGS17             SEMA4A            KREMEN1               CDK4 
    ##              0.641              0.641              0.642              0.642 
    ##               ZNF2             PTPN13              GOSR1            B4GALT3 
    ##              0.642              0.642              0.642              0.642 
    ##             IL20RA               STC2             TRIM37           TMEM151B 
    ##              0.642              0.642              0.642              0.642 
    ##            PACSIN3             PDIK1L              KLHL7             ZNF641 
    ##              0.642              0.643              0.643              0.643 
    ##            TNFAIP1             RAB11B              SSTR2               NLE1 
    ##              0.643              0.643              0.643              0.643 
    ##             B3GAT1               IL21              TRIM7           KIAA1958 
    ##              0.643              0.643              0.643              0.643 
    ##               PURB            TNFSF10             TMEM60              PHF10 
    ##              0.643              0.643              0.643              0.644 
    ##              USP38             ZBTB48           SIGLEC10               CHL1 
    ##              0.644              0.644              0.644              0.644 
    ##               CRBN              LRIG3           SERPINI1              TMEM2 
    ##              0.644              0.644              0.644              0.644 
    ##              C2CD3              SUDS3               ELF4              RAD51 
    ##              0.644              0.644              0.644              0.644 
    ##             MGAT4A           TMEM132E             CACNB3               AQP3 
    ##              0.644              0.644              0.644              0.644 
    ##           MPHOSPH9              AP3M1             GALNT7              CELF6 
    ##              0.644              0.644              0.644              0.645 
    ##             TMEM56              MYO1G               RGS3             CUEDC2 
    ##              0.645              0.645              0.645              0.645 
    ##            CLPTM1L           UHRF1BP1             ZBTB14             ACTL6B 
    ##              0.645              0.645              0.645              0.645 
    ##              DZIP3                FOS             CDKN1B               RIT1 
    ##              0.645              0.645              0.645              0.645 
    ##              USP28             LUC7L2             ARPC1B               DIO3 
    ##              0.645              0.645              0.645              0.645 
    ##                 CP               GZF1           C14orf93               RRAD 
    ##              0.646              0.646              0.646              0.646 
    ##             TIMM50         AC074212.3               EDC3            ADAMTS7 
    ##              0.646              0.646              0.646              0.646 
    ##              CENPN              GPR22             TRIP10              PSKH1 
    ##              0.646              0.646              0.646              0.646 
    ##            GALNT11              PSME1             CEP250               CUBN 
    ##              0.646              0.646              0.647              0.647 
    ##               RBP3           PPP1R16A            TCP11L1           ADAMTSL3 
    ##              0.647              0.647              0.647              0.647 
    ##              SGOL2              NLRP6               ZG16               IL1B 
    ##              0.647              0.647              0.647              0.648 
    ##              LAMC2               SUN2              ITGB4               HIC1 
    ##              0.648              0.648              0.648              0.648 
    ##               RRN3              LPPR5            CXorf57              SYTL5 
    ##              0.648              0.648              0.648              0.648 
    ##               SHC1             KISS1R             ZNF341             FBXO10 
    ##              0.648              0.648              0.648              0.648 
    ##               PEX6              ARL8B           BAIAP2L1              TDRD9 
    ##              0.648              0.649              0.649              0.649 
    ##              ITGB7              ENOX1              WIPI1             FAM57A 
    ##              0.649              0.649              0.649              0.649 
    ##              NPHP3            SLC25A1               LCA5               COG1 
    ##              0.649              0.649              0.649              0.649 
    ##              FXYD6            TMEM168             VPS13B              PPWD1 
    ##              0.649              0.649              0.650              0.650 
    ##             BCAP31           CALCOCO2             CCDC80               YES1 
    ##              0.650              0.650              0.650              0.650 
    ##            ABHD17B             KLHL24               RHAG               NID1 
    ##              0.650              0.650              0.650              0.650 
    ##              POLA2              LOXL2             SLC2A3              FGF16 
    ##              0.650              0.650              0.650              0.650 
    ##             FAM98B             COLCA2             KLHL23                NTM 
    ##              0.650              0.651              0.651              0.651 
    ##            ST8SIA2              MEAF6               HTR4              CDCP1 
    ##              0.651              0.651              0.651              0.651 
    ##               WNT4            SLC35E1            PLEKHG3              FOLR1 
    ##              0.651              0.651              0.651              0.652 
    ##               EPOR              TTYH2              EXOC7            TAX1BP1 
    ##              0.652              0.652              0.652              0.652 
    ##             ZC3HC1               NBR1             HSPBP1              MED25 
    ##              0.652              0.652              0.652              0.652 
    ##             NECAP1              UIMC1               PGM3              USP51 
    ##              0.653              0.653              0.653              0.653 
    ##               CD58               CDS1              PRMT2              ABCB9 
    ##              0.653              0.653              0.653              0.653 
    ##                  T             GATAD1               ALG9              SMCR8 
    ##              0.653              0.653              0.653              0.654 
    ##              CSPG5              ACSF2              TRPC6               BMP3 
    ##              0.654              0.654              0.654              0.654 
    ##               REC8            RABGEF1             NBEAL1              MTMR2 
    ##              0.654              0.654              0.654              0.654 
    ##              KCTD7           KIAA1715             NHP2L1             SLC7A6 
    ##              0.654              0.655              0.655              0.655 
    ##             UBE2E3             IPCEF1              NOL11             ATP2A3 
    ##              0.655              0.655              0.655              0.655 
    ##               PHB2               CCIN              MYD88              LAMB2 
    ##              0.655              0.655              0.655              0.655 
    ##               RAI2             PTPN21             FRMPD1             RNF213 
    ##              0.656              0.656              0.656              0.656 
    ##              NUDT3             CIRH1A               TAP1             NKAIN1 
    ##              0.656              0.656              0.656              0.656 
    ##              KLHL5             TCEAL4                CD6             SEMA3D 
    ##              0.656              0.656              0.656              0.657 
    ##            DPY19L3               ECT2               BRDT                JTB 
    ##              0.657              0.657              0.657              0.657 
    ##              PHKG2               DVL2              KCNF1              P4HA2 
    ##              0.657              0.657              0.657              0.657 
    ##            CREB3L2             ACTR10              BEND5           NOTCH2NL 
    ##              0.657              0.658              0.658              0.658 
    ##              FRAS1               PGM5              SESN1              VPS45 
    ##              0.658              0.658              0.658              0.658 
    ##              PELI3              IL2RA              YPEL1              CASP2 
    ##              0.658              0.658              0.658              0.658 
    ##              NTRK1            EHBP1L1             GPR162             TRIM36 
    ##              0.658              0.659              0.659              0.659 
    ##               SCG3             TRIM44              TMED4             HOMER2 
    ##              0.659              0.659              0.659              0.659 
    ##              PRRX1              TRABD              ENPP4               CCM2 
    ##              0.659              0.659              0.660              0.660 
    ##             HOXC13            ANGPTL2              SOCS3              VASH2 
    ##              0.660              0.660              0.660              0.660 
    ##                LAT              GUCD1             TOMM22               COG7 
    ##              0.660              0.660              0.660              0.661 
    ##               SNX5            PPP1R1B               AQP9              GLRX5 
    ##              0.661              0.661              0.661              0.661 
    ##              CPT1A            TBC1D9B              PRUNE             DENND3 
    ##              0.661              0.661              0.661              0.662 
    ##            HSD11B1              PDCD1           C12orf44               MIA2 
    ##              0.662              0.662              0.662              0.662 
    ##              PPM1B              NUP35            FAM117B                GHR 
    ##              0.662              0.662              0.662              0.662 
    ##           CDK5RAP2               GLRB              ZMYM1              BMPER 
    ##              0.663              0.663              0.663              0.663 
    ##             NFE2L2              TCEB1             MAP2K3             COL8A2 
    ##              0.663              0.663              0.663              0.663 
    ##               EIF6              USP30              MKRN2             ZBTB24 
    ##              0.663              0.663              0.664              0.664 
    ##              ITGA8             ZNF346               PEX1               CD44 
    ##              0.664              0.664              0.664              0.664 
    ##              DDX54               RHOA             PPP3R1               SDC3 
    ##              0.664              0.664              0.664              0.664 
    ##              WDR60           ARHGEF15              AFAP1              SAMD1 
    ##              0.665              0.665              0.665              0.665 
    ##               MYCL              PRPH2             UBE2D3             NKX2-2 
    ##              0.665              0.665              0.665              0.665 
    ##              TTPAL             GOLPH3             ALKBH8             PIK3C3 
    ##              0.665              0.665              0.666              0.666 
    ##                HPD              WARS2             BTBD19            PRELID1 
    ##              0.666              0.666              0.666              0.666 
    ##             DYNLL2              KCNG2              TRPV6             STARD3 
    ##              0.666              0.666              0.666              0.666 
    ##              KCNH8               PEMT            SLCO2B1              SIDT2 
    ##              0.666              0.666              0.667              0.667 
    ##            ATP6V0C             CTDSP2               PLS1             FAM73B 
    ##              0.667              0.667              0.667              0.667 
    ##              PRAF2              MFAP4            SLC35F4              GAPDH 
    ##              0.667              0.667              0.667              0.667 
    ##             FLVCR1              DSCC1               CBX3              RIMS4 
    ##              0.668              0.668              0.668              0.668 
    ##             DFNB31              WDR92             TRANK1               GBX1 
    ##              0.668              0.668              0.668              0.668 
    ##              TRIP4               UTP6            SLC13A5             FAM43A 
    ##              0.668              0.668              0.669              0.669 
    ##            TMEM189     TMEM189-UBE2V1            CNEP1R1            PPP1R18 
    ##              0.669              0.669              0.669              0.669 
    ##               MXD1               EID3              UBA52               KRT2 
    ##              0.669              0.669              0.670              0.670 
    ##             ZNF653               ISL2              SKAP2              KCNS1 
    ##              0.670              0.670              0.670              0.670 
    ##                CA8              TRIP6             CLEC2D               EMC2 
    ##              0.670              0.670              0.670              0.670 
    ##            C12orf4             SCAND3              PICK1              RAB21 
    ##              0.670              0.671              0.671              0.671 
    ##              FOXO6              HSPA2              MLST8              TINF2 
    ##              0.671              0.671              0.671              0.671 
    ##               BIVM              DRAM1              UCHL5             RNF128 
    ##              0.672              0.672              0.672              0.672 
    ##             OR10C1             SUPT7L               RECK                MLX 
    ##              0.672              0.672              0.672              0.672 
    ##              SZRD1               SDC2              CRLF3            TUBGCP5 
    ##              0.672              0.672              0.673              0.673 
    ##              THOC6               PNO1              MED30             TCEAL3 
    ##              0.673              0.673              0.673              0.673 
    ##             DNAH11              IFRD1            C5orf42             CPAMD8 
    ##              0.674              0.674              0.674              0.674 
    ##             TUBB4A              DHRS3             KCNK10            TMEM121 
    ##              0.674              0.674              0.674              0.674 
    ##               BRAP            ZSCAN18             CDCA7L            GUCY1B3 
    ##              0.674              0.674              0.675              0.675 
    ##                DHH              DHX57            ERCC6L2              DNAH8 
    ##              0.675              0.675              0.675              0.675 
    ##             MTMR10              TMOD2              CDAN1               TESC 
    ##              0.675              0.675              0.675              0.675 
    ##                CA7            LRRFIP1             INPP5E              NCBP2 
    ##              0.675              0.675              0.675              0.675 
    ##              VTI1B            RABGGTB            GORASP1              ZBTB3 
    ##              0.675              0.675              0.675              0.676 
    ##              HBS1L              TFCP2               LDB3            FAM118B 
    ##              0.676              0.676              0.676              0.676 
    ##               EVX2              LRCH4              RRBP1               PDPR 
    ##              0.676              0.676              0.676              0.676 
    ##      RP11-766F14.2             GPR174               RPF2              CENPP 
    ##              0.677              0.677              0.677              0.677 
    ##            DENND1C               MTX2              MAT2B              FOXA1 
    ##              0.677              0.677              0.677              0.677 
    ##             VPS13C             ALKBH6               GPAM              KCNK2 
    ##              0.677              0.677              0.677              0.677 
    ##              SCYL1              EXTL3            CCDC168              DHX33 
    ##              0.677              0.677              0.678              0.678 
    ##            DYNC1I1              DDIT3              KRT31               GOT1 
    ##              0.678              0.678              0.678              0.678 
    ##             CDK11B             SF3B14            SLC30A9              WSCD1 
    ##              0.678              0.679              0.679              0.679 
    ##              NSMAF            SLC12A4              C2CD2             SEMA3E 
    ##              0.679              0.679              0.679              0.679 
    ##              SP100                MOS               MUC1             SLAMF6 
    ##              0.679              0.679              0.679              0.679 
    ##              SASS6              AGBL5               DLX5              WDR24 
    ##              0.679              0.679              0.680              0.680 
    ##              STUB1               ETFA               PKP3            TSPAN31 
    ##              0.680              0.680              0.680              0.680 
    ##               HTR7              DIEXF              SEPT4              CCNA1 
    ##              0.680              0.680              0.680              0.680 
    ##              STON2           SLC2A4RG              MYO1C             COX4I1 
    ##              0.680              0.680              0.680              0.680 
    ##             FAM83B             HIATL1              LRIF1            ZDHHC16 
    ##              0.680              0.680              0.680              0.680 
    ##               CLK4               IDNK               MLTK              PKHD1 
    ##              0.681              0.681              0.681              0.681 
    ##            TIMM17A             COMMD9             APCDD1              GPR19 
    ##              0.681              0.681              0.681              0.681 
    ##             PCNXL2              RNF14             POU1F1              ADCY3 
    ##              0.681              0.681              0.681              0.681 
    ##              SNAI1               DARS              ACKR3               TLE2 
    ##              0.681              0.682              0.682              0.682 
    ##               LIG3             UBXN2B            TUBGCP2             PTCHD4 
    ##              0.682              0.682              0.682              0.682 
    ##              SCMH1             TRIM66             RAB22A              BARX1 
    ##              0.682              0.682              0.682              0.683 
    ##           ADAMTS13             ZNF625              RAB9B               TMX4 
    ##              0.683              0.683              0.683              0.683 
    ##             ZNF189               PHF7              CDCA8             CLUAP1 
    ##              0.683              0.683              0.683              0.683 
    ##              ARPC3            SLC48A1              TOP3A             ATP1B3 
    ##              0.683              0.683              0.683              0.683 
    ##              RUFY1             ZNF821              GPSM1              VPS39 
    ##              0.684              0.684              0.684              0.684 
    ##                SMO               LGI4               SSR2             ACTR1B 
    ##              0.684              0.684              0.684              0.684 
    ##             NHLRC2           ANKRD34C               E2F6             COPS7A 
    ##              0.684              0.684              0.684              0.685 
    ##               LGMN               PAG1          CTTNBP2NL           DNASE1L1 
    ##              0.685              0.685              0.685              0.685 
    ##              ERBB3             RASSF5               TLE6               WNT2 
    ##              0.685              0.685              0.685              0.685 
    ##           C15orf57           SERPINE1            SLC16A9              TACC2 
    ##              0.685              0.685              0.685              0.685 
    ##             IL10RA               CCR1           KIAA1841             KIF20B 
    ##              0.685              0.685              0.685              0.686 
    ##             HS3ST2             MICAL2               LHX4              CASP3 
    ##              0.686              0.686              0.686              0.686 
    ##           SREK1IP1              PRSS8              IQCA1               GAR1 
    ##              0.686              0.686              0.686              0.686 
    ##               MYH2               TGM3            PPP1R21              HPSE2 
    ##              0.686              0.686              0.686              0.686 
    ##              ZNF79              ZGLP1              OGDHL             CAMK1G 
    ##              0.686              0.686              0.686              0.686 
    ##               PFN1              FOXD3            PPIP5K1              PANK1 
    ##              0.686              0.687              0.687              0.687 
    ##             MB21D2               AZI2             HDAC11             DNAH10 
    ##              0.687              0.687              0.687              0.687 
    ##              FGF18               IL16              BARX2              PODXL 
    ##              0.687              0.687              0.687              0.687 
    ##                TEC              NUP88                IPP             SPRED3 
    ##              0.687              0.687              0.687              0.687 
    ##              TRADD             FAM8A1              SNTB1              OLFM3 
    ##              0.687              0.687              0.688              0.688 
    ##               TBX6              TTC40             IL17RD              ZAR1L 
    ##              0.688              0.688              0.688              0.688 
    ##      RP11-508N12.4               MAFG           ADAMTS12              MYO5B 
    ##              0.688              0.688              0.688              0.688 
    ##             CHMP1B              BACE2            GPATCH1                GPI 
    ##              0.688              0.689              0.689              0.689 
    ##               MLK4             STRADB             LARP1B               DSC2 
    ##              0.689              0.689              0.689              0.689 
    ##             MAP4K2            PHACTR2            PHACTR4              TONSL 
    ##              0.689              0.689              0.689              0.689 
    ##             CCRN4L            ZFYVE16              GDPD1              ATP9B 
    ##              0.689              0.689              0.689              0.689 
    ##             ZNF527              DOC2A               DTX1              ZNF48 
    ##              0.689              0.690              0.690              0.690 
    ##              GFRA2              RNF11             HSPA13             WRAP53 
    ##              0.690              0.690              0.690              0.690 
    ##               LFNG               NDE1            FAM131A              PRIM1 
    ##              0.690              0.690              0.690              0.690 
    ##             ZDHHC2               ZNF8               ZNF8              LAMA3 
    ##              0.690              0.690              0.690              0.691 
    ##              SYNJ2        CHURC1-FNTB             MTMR14              TPCN2 
    ##              0.691              0.691              0.691              0.691 
    ##            C5orf24           C14orf39             SLC9A5               LNX2 
    ##              0.691              0.691              0.691              0.691 
    ##               TPM3              FAM3C               CSF3              FCGRT 
    ##              0.691              0.691              0.691              0.691 
    ##               OTOA             CACNB2              TNNI3             ACVRL1 
    ##              0.692              0.692              0.692              0.692 
    ##              POLD2              CERS1                CFB               HPS1 
    ##              0.692              0.692              0.692              0.692 
    ##            CCDC121              VARS2              DAPK3             PPAP2A 
    ##              0.692              0.692              0.692              0.692 
    ##           KIAA1024              MRPS5              ZMYM6            ARL14EP 
    ##              0.692              0.693              0.693              0.693 
    ##              RMDN3         BIVM-ERCC5               HRH3             ZBTB40 
    ##              0.693              0.693              0.693              0.693 
    ##              SALL2               ALPL            FAM185A              TBPL1 
    ##              0.693              0.693              0.693              0.693 
    ##              SCRN1            SLITRK2              TREX1             TXNDC9 
    ##              0.694              0.694              0.694              0.694 
    ##              EVA1B             TCEAL5              ANKS6              CPT1C 
    ##              0.694              0.694              0.694              0.694 
    ##              GPER1              IFFO1            C9orf41              DNAH5 
    ##              0.695              0.695              0.695              0.695 
    ##              MTUS1                 C5               STC1              ACER2 
    ##              0.695              0.695              0.695              0.695 
    ##              RIOK2             POU3F4              PROSC               ISM1 
    ##              0.695              0.695              0.696              0.696 
    ##              MUC13              NPHS1            OSBPL10             STOML2 
    ##              0.696              0.696              0.696              0.696 
    ##            TMEM246              UBAC2              PTPRQ               ARL1 
    ##              0.696              0.696              0.696              0.696 
    ##              TTLL4              HABP4             ATXN10           TRAPPC13 
    ##              0.696              0.697              0.697              0.697 
    ##              PSME2            PLEKHA3             TBC1D8              KCTD2 
    ##              0.697              0.697              0.697              0.697 
    ##             FAM63B               MUSK                KMO             RUNDC1 
    ##              0.698              0.698              0.698              0.698 
    ##              SIRT6               GNB4               SKP1               NACA 
    ##              0.698              0.698              0.698              0.698 
    ##            SPATA20              CD226              PBDC1              DEAF1 
    ##              0.698              0.698              0.698              0.698 
    ##               ATF6             PRSS56              GLRA3              PTPRR 
    ##              0.698              0.698              0.698              0.698 
    ##               BMP5               GPT2             FAM78A               JAM2 
    ##              0.698              0.698              0.698              0.698 
    ##               CLMN             BTF3L4              MED24              PCBP4 
    ##              0.698              0.698              0.698              0.698 
    ##               SAV1           FAM189A1             FILIP1              PALMD 
    ##              0.698              0.698              0.698              0.698 
    ##               MSH5             GABRA3              HTR1A              UBE2W 
    ##              0.699              0.699              0.699              0.699 
    ##            SLC5A12              PDE1C               HMX1             TBC1D1 
    ##              0.699              0.699              0.699              0.699 
    ##              RRM2B              MFSD4              SAR1A              HOXA2 
    ##              0.699              0.699              0.699              0.699 
    ##              NDC80              INTS4               TPK1             KCNIP2 
    ##              0.699              0.700              0.700              0.700 
    ##               MIOS             IFNGR1             AFG3L2              ATF6B 
    ##              0.700              0.700              0.700              0.700 
    ##              BEND7             TRIM35              GPD1L              STMN4 
    ##              0.700              0.700              0.701              0.701 
    ##              CPSF3              GNA15            TBC1D2B             ENTPD5 
    ##              0.701              0.701              0.701              0.701 
    ##               THPO             KATNA1               BORA           SLC9A3R2 
    ##              0.701              0.701              0.701              0.701 
    ##              FBXO9              GPR63              PCYT2               SBK1 
    ##              0.701              0.701              0.702              0.702 
    ##              WDR96            SLC27A4               GJD2              UBE2B 
    ##              0.702              0.702              0.702              0.702 
    ##              TRPC3               CARF            SIGLEC6              SIKE1 
    ##              0.702              0.703              0.703              0.703 
    ##             TMBIM6              REEP3               DGKQ              NCOA4 
    ##              0.703              0.703              0.703              0.703 
    ##             NT5C3A            ADAMTS5             KCNIP3              MSMO1 
    ##              0.703              0.703              0.703              0.703 
    ##                RP1                SFN               EBF4               TAL1 
    ##              0.704              0.704              0.704              0.704 
    ##           C14orf80           ARHGEF28              MYO1D            B4GALT1 
    ##              0.704              0.704              0.704              0.704 
    ##               PAWR              NT5C2            SLC39A4               GOPC 
    ##              0.704              0.704              0.704              0.704 
    ##             POU3F1             DNAJB9              IGSF9              PAQR5 
    ##              0.705              0.705              0.705              0.705 
    ##              HERC5              C2CD5             ZNF174               NIP7 
    ##              0.705              0.705              0.705              0.705 
    ##             IMPDH1             SCAMP3               PAK6              TRAIP 
    ##              0.705              0.705              0.705              0.705 
    ##            FAM110B              CAPN1              ITGB5             DAZAP2 
    ##              0.705              0.706              0.706              0.706 
    ##             ZNF613               LIFR            ANKRD32            PCDHGB7 
    ##              0.706              0.706              0.706              0.706 
    ##              TPRA1              NYAP2            ST6GAL2            IL12RB2 
    ##              0.706              0.706              0.707              0.707 
    ##              MAST2               ORC3         AC069368.3              SNTA1 
    ##              0.707              0.707              0.707              0.707 
    ##              CPNE9              TACC3            ISG20L2             CEP85L 
    ##              0.707              0.707              0.707              0.707 
    ##             ANGEL2              CYTIP             TOMM34              IFT88 
    ##              0.708              0.708              0.708              0.708 
    ##              SCN9A              PRAM1            ZNF280B               DGKE 
    ##              0.708              0.708              0.708              0.708 
    ##              PGBD5             SLAMF7             GSTT2B               LNX1 
    ##              0.708              0.708              0.708              0.708 
    ##              RAB43    RNASEK-C17orf49              GPR12                CD4 
    ##              0.708              0.708              0.708              0.709 
    ##              SPSB3           ADAMTS14                ID3            CEACAM1 
    ##              0.709              0.709              0.709              0.709 
    ##             POLR2I              KLHL1             ZNF740                B2M 
    ##              0.709              0.709              0.709              0.709 
    ##               KLF2              HTR1B               JDP2              RASA3 
    ##              0.710              0.710              0.710              0.710 
    ##             ZNF850            EIF2AK1              ROBO3            ADAMTS4 
    ##              0.710              0.710              0.710              0.710 
    ##              CROCC                ATM              PROB1              WNT7A 
    ##              0.710              0.710              0.710              0.710 
    ##              TOR4A               CD28               GMIP              ASAP3 
    ##              0.710              0.711              0.711              0.711 
    ##           LEPROTL1              LLGL2              MUC22               RBL1 
    ##              0.711              0.711              0.711              0.711 
    ##              TUFT1               CALU              IDH3G              IL23A 
    ##              0.711              0.711              0.711              0.711 
    ##              TMED8              AWAT2            DCLRE1C                ITK 
    ##              0.711              0.711              0.711              0.712 
    ##            MSANTD4            DNAJC17              EPHA8               BCL2 
    ##              0.712              0.712              0.712              0.712 
    ##               TCF7               ZFP1            SLC36A3               COG4 
    ##              0.712              0.712              0.712              0.713 
    ##               PLK3            THEMIS2              RRAGA               NUF2 
    ##              0.713              0.713              0.713              0.713 
    ##            SLC45A4           NAALADL2             ZNF703             PTDSS2 
    ##              0.713              0.713              0.713              0.713 
    ##              BCL7B              CRHR1             TRIOBP                ILK 
    ##              0.714              0.714              0.714              0.714 
    ##            ALDH4A1              ANXA6               HPS6               PPIA 
    ##              0.714              0.714              0.714              0.714 
    ##              RELL2              LPAR1             DNAJA4              LRFN2 
    ##              0.714              0.714              0.714              0.714 
    ##               DLK1             GPR123              LTA4H               FUT9 
    ##              0.714              0.714              0.715              0.715 
    ##              AVPR2             STK17B              RGSL1              ADPGK 
    ##              0.715              0.715              0.715              0.715 
    ##              AMER2              TMTC1             GRIN2C             MYO18B 
    ##              0.715              0.715              0.715              0.715 
    ##             TAGLN3               ELK3             DYNLT3              ESRP2 
    ##              0.715              0.715              0.715              0.715 
    ##             SAMD14              MGAT2            TUBGCP4              FASLG 
    ##              0.715              0.716              0.716              0.716 
    ##               RHCG               ARF3             ATPAF1               CFL1 
    ##              0.716              0.716              0.716              0.716 
    ##            L3MBTL2            ZKSCAN7             FAM46B             PRUNE2 
    ##              0.716              0.716              0.717              0.717 
    ##              PAIP2           SLC25A27              CILP2              DOCK6 
    ##              0.717              0.717              0.717              0.717 
    ##          KIAA0895L             MAD2L2               OGFR              SPEF2 
    ##              0.717              0.717              0.717              0.717 
    ##       CTC-273B12.7              VSTM4           AASDHPPT            SERTAD4 
    ##              0.717              0.717              0.717              0.717 
    ##               MTG1              ATG4D               ATG3            TSPAN18 
    ##              0.717              0.717              0.718              0.718 
    ##             SOWAHB             ARPC5L              VLDLR             SHKBP1 
    ##              0.718              0.718              0.718              0.718 
    ##               PDHX              HOXB1              ABCD2              CHRM2 
    ##              0.718              0.718              0.718              0.718 
    ##            NDUFAF4           SLC25A36               MTX1               SYBU 
    ##              0.718              0.718              0.718              0.718 
    ##              WDR66              IER5L               GGA3             RAB7L1 
    ##              0.718              0.718              0.719              0.719 
    ##             TOLLIP           SLC16A14            GADD45G             CLNS1A 
    ##              0.719              0.719              0.719              0.719 
    ##            FAM216A              PIBF1           ADAMTS15              MIEN1 
    ##              0.719              0.719              0.719              0.719 
    ##               ENSA            IL17REL           KIAA0586              PFDN2 
    ##              0.719              0.719              0.719              0.719 
    ##              MYLK3              SFTPB            PCDHGB1              SPAG5 
    ##              0.720              0.720              0.720              0.720 
    ##             ZNF572              NCEH1            DNAJC21              IMPG2 
    ##              0.720              0.720              0.720              0.720 
    ##              GRIK1             SRSF12                CRH              ITGAX 
    ##              0.720              0.720              0.720              0.721 
    ##               KRR1              ZNF34            TMPRSS6           C10orf88 
    ##              0.721              0.721              0.721              0.721 
    ##           TNFRSF19              CLCA1              DNMBP           TMEM161B 
    ##              0.721              0.721              0.721              0.721 
    ##             TXNRD1               DSTN              DPCR1            PIP4K2C 
    ##              0.721              0.721              0.722              0.722 
    ##               ARF6              DNAH9             VPS37A               UGDH 
    ##              0.722              0.722              0.722              0.722 
    ##               UBA5            CCDC136              MORN4               OAZ1 
    ##              0.722              0.722              0.722              0.722 
    ##              FBXW4             NT5DC2              OVGP1              RPL35 
    ##              0.723              0.723              0.723              0.723 
    ##              WDR59             FAM83G                FXN               FBP1 
    ##              0.723              0.723              0.723              0.723 
    ##              BCL7C             CCDC13             PARD3B            DENND2C 
    ##              0.723              0.723              0.723              0.723 
    ##               LGR5               ACP5               RYBP            RHOBTB3 
    ##              0.724              0.724              0.724              0.724 
    ##               VAV3             SCARB2             NKX1-1             RNF185 
    ##              0.724              0.724              0.724              0.724 
    ##            SLCO5A1             NDUFA8              WDR61            SLC26A5 
    ##              0.724              0.724              0.724              0.724 
    ##     TMEM110-MUSTN1           KIAA1324             KLHL42               RRS1 
    ##              0.724              0.724              0.724              0.724 
    ##              USP35              DNAI1             FAM69B              POSTN 
    ##              0.724              0.725              0.725              0.725 
    ##              TTC26            FAM172A            TBC1D13              BCL7A 
    ##              0.725              0.725              0.725              0.725 
    ##            FILIP1L              CALN1               XBP1              SAMD8 
    ##              0.725              0.725              0.725              0.726 
    ##                MVP            MTHFD1L              ANKS3              AURKC 
    ##              0.726              0.726              0.726              0.726 
    ##             SCARF1               UROS              SOX21              VPS25 
    ##              0.726              0.726              0.726              0.727 
    ##              FNDC5               IST1              LPIN1            FAM19A2 
    ##              0.727              0.727              0.727              0.727 
    ##              ACSL5               PIN4              PTGES           KIAA1199 
    ##              0.727              0.727              0.727              0.728 
    ##            PLEKHA1              EBAG9             FAM73A              PDGFD 
    ##              0.728              0.728              0.728              0.728 
    ##                CBS           C20orf26             HOXC10               DKK1 
    ##              0.728              0.728              0.728              0.728 
    ##              RIC8A             CCDC79               URI1           C1orf116 
    ##              0.728              0.728              0.728              0.729 
    ##             HGSNAT             MARCH3             NUP210             POLR2M 
    ##              0.729              0.729              0.729              0.729 
    ##            FGFR1OP             ZNF672                OXT               SLA2 
    ##              0.729              0.729              0.729              0.729 
    ##           C1orf112              CNTD1              CHDC2               PPT2 
    ##              0.730              0.730              0.730              0.730 
    ##              RBM19              KCNA1             PLSCR1                ZP2 
    ##              0.730              0.730              0.730              0.730 
    ##              TCEA3              PCSK1             FGFRL1              LAMA4 
    ##              0.730              0.730              0.730              0.730 
    ##               WDR3            TNFAIP2            PLEKHH3              DDOST 
    ##              0.730              0.731              0.731              0.731 
    ##              ZMAT4               CUL7              BCL10              HAUS5 
    ##              0.731              0.731              0.731              0.731 
    ##             GOLIM4               GRB7               SDHD               MAFK 
    ##              0.731              0.731              0.731              0.731 
    ##              NPTX2              AURKB              LAMA2              RABL5 
    ##              0.731              0.732              0.732              0.732 
    ##               RBP7           ADAMTS17             PKD1L1              KCNC4 
    ##              0.732              0.732              0.732              0.732 
    ##              MIEF1            ANKRD26              CENPF               DVL1 
    ##              0.733              0.733              0.733              0.733 
    ##           C1orf106            ZSCAN10               FGD3              MYO7B 
    ##              0.733              0.733              0.733              0.733 
    ##              STAM2               TPI1            TRAPPC9              CCDC9 
    ##              0.733              0.733              0.733              0.733 
    ##             RBMXL2               HPS5              HCLS1              HOXC8 
    ##              0.733              0.733              0.733              0.734 
    ##                KIN              PCSK6              KRT19               DYSF 
    ##              0.734              0.734              0.734              0.734 
    ##               IRX4               MDFI           ALS2CR11              ZFP14 
    ##              0.734              0.734              0.735              0.735 
    ##             POFUT2            ANAPC11              GHITM               PREB 
    ##              0.735              0.735              0.735              0.735 
    ##              ZNF10           ADAMTSL2              HOXC4              HOXC4 
    ##              0.735              0.735              0.735              0.735 
    ##               BBS7            B3GALT4             IKBKAP              SKOR1 
    ##              0.735              0.735              0.735              0.736 
    ##              GLIS3            WBSCR22              HARS2              NPDC1 
    ##              0.736              0.736              0.736              0.736 
    ##              GABRQ              CENPV            SIGMAR1              DDX41 
    ##              0.736              0.736              0.736              0.736 
    ##             ZBTB8A              TYRO3               RAG1       RP4-539M6.19 
    ##              0.736              0.736              0.736              0.736 
    ##            ZNF585A              TNIP2               LHX1             KDELR2 
    ##              0.736              0.736              0.737              0.737 
    ##            FAM19A1               CNST               BTG4               TLR9 
    ##              0.737              0.737              0.737              0.737 
    ##              OSTM1               RCN2           TRAF3IP1              SETD7 
    ##              0.737              0.737              0.737              0.737 
    ##              RASD2              NHEJ1           C17orf70             NDUFS8 
    ##              0.737              0.738              0.738              0.738 
    ##             KCNIP1            EXOSC10               RRP1              TBRG4 
    ##              0.738              0.738              0.738              0.738 
    ##              LETM1              GTSE1                AK9               BTG1 
    ##              0.738              0.738              0.738              0.738 
    ##             HS3ST5               NOA1              RPL23            RNASEH1 
    ##              0.738              0.738              0.738              0.738 
    ##              INTS9               ING4              ZBTB9             MRPL19 
    ##              0.738              0.738              0.738              0.738 
    ##              LAMB3             ETNPPL               GGT1           C11orf58 
    ##              0.738              0.738              0.738              0.739 
    ##              A4GNT               RGMB               RTTN               TUFM 
    ##              0.739              0.739              0.739              0.739 
    ##           RAD51AP1             UNC45A              F13A1             ZNF831 
    ##              0.739              0.739              0.739              0.739 
    ##            GTPBP10               TBL3              ITGA1                TDG 
    ##              0.739              0.739              0.739              0.739 
    ##             TMEFF2             NEURL1               ASNS             RNF157 
    ##              0.739              0.739              0.740              0.740 
    ##                EGF              OR2H1              RASEF               TTC8 
    ##              0.740              0.740              0.740              0.740 
    ##              WNT8B              NEK10              XRCC1             PARP10 
    ##              0.740              0.740              0.740              0.740 
    ##              GFPT2             KIF20A               INTU            MADCAM1 
    ##              0.740              0.740              0.740              0.740 
    ##              CPSF2              FOLH1              SNX29      RP11-724O16.1 
    ##              0.741              0.741              0.741              0.741 
    ##               CERK              DCAF4              RAPSN             TEX264 
    ##              0.741              0.741              0.741              0.741 
    ##             MAGOHB             SEC23A               FPGS             PNPLA2 
    ##              0.741              0.741              0.741              0.741 
    ##             FBXO25               LHX3              CASC4              FLOT2 
    ##              0.742              0.742              0.742              0.742 
    ##            TMEM50A            ZDHHC21             ERLEC1            PLEKHD1 
    ##              0.742              0.742              0.742              0.742 
    ##          ARHGAP11A              USP18             NOTCH4               HCN3 
    ##              0.742              0.742              0.742              0.742 
    ##               WDR6              CHRM5               FBN3             HAPLN4 
    ##              0.742              0.743              0.743              0.743 
    ##               ASPM              RNPC3              PTPLB             CEP104 
    ##              0.743              0.743              0.743              0.743 
    ##            TBC1D8B            SH3GLB2             DIAPH3              LYZL4 
    ##              0.743              0.743              0.743              0.743 
    ##              VOPP1               NCF2              PALM2               PRG4 
    ##              0.743              0.743              0.743              0.743 
    ##             CCDC50      RP1-130H16.18          RAB11FIP5            SLC35G1 
    ##              0.744              0.744              0.744              0.744 
    ##              PSEN2              LRP11            PCDHGC5              VGLL3 
    ##              0.744              0.744              0.744              0.745 
    ##             MOGAT2              KLHL8               ERI3              EVA1C 
    ##              0.745              0.745              0.745              0.745 
    ##            HSD17B4               PARL            ZSCAN25           KIAA0196 
    ##              0.745              0.745              0.745              0.745 
    ##             POLR2C              UBAC1            SLC45A1              VPS28 
    ##              0.745              0.745              0.745              0.745 
    ##           C19orf33             SYNGR3             COL9A3             PPP1R2 
    ##              0.746              0.746              0.746              0.746 
    ##           SERPINH1              REEP4             GATSL3             MROH2A 
    ##              0.746              0.746              0.746              0.747 
    ##               BAG4              ATP4A           C12orf55            ZDHHC20 
    ##              0.747              0.747              0.747              0.747 
    ##              AMACR            DCLRE1A             CYP2U1              PPM1F 
    ##              0.747              0.747              0.747              0.747 
    ##              TRMT6             GEMIN5              MERTK               SLBP 
    ##              0.747              0.748              0.748              0.748 
    ##              CLDN2              FOSL1           ATP6V0A2              MPZL1 
    ##              0.748              0.748              0.748              0.748 
    ##             PNPLA1             GUCY2D               ESX1              RBBP8 
    ##              0.748              0.748              0.749              0.749 
    ##              USP13              USP16               GGA2               NSA2 
    ##              0.749              0.749              0.749              0.749 
    ##              WNT8A       RP11-295K3.1               OTOG               CPA5 
    ##              0.749              0.749              0.749              0.749 
    ##            NDUFAB1               STX6               EVPL              CPLX2 
    ##              0.749              0.749              0.749              0.750 
    ##              GMCL1               IRGQ              NELL1               RFC5 
    ##              0.750              0.750              0.750              0.750 
    ##              GNRHR             DIXDC1             DNAH17             ZNF706 
    ##              0.750              0.750              0.750              0.750 
    ##            FASTKD2             ZNF570             ELOVL7            FAM120B 
    ##              0.750              0.750              0.750              0.750 
    ##               AHRR            COL17A1               TNS4             ANKFN1 
    ##              0.750              0.750              0.751              0.751 
    ##            FAM189B            PPFIBP1            SLC43A1                BAX 
    ##              0.751              0.751              0.751              0.751 
    ##                BLM            C1orf86            ADCYAP1              RAB5B 
    ##              0.751              0.751              0.752              0.752 
    ##              WDFY1                LIF             SLC4A5            ZC3H12D 
    ##              0.752              0.752              0.752              0.752 
    ##              TCAIM             ERGIC3               GLE1             ZNF276 
    ##              0.752              0.752              0.752              0.752 
    ##            MICALL1              HTRA1                IL6            COL13A1 
    ##              0.752              0.752              0.753              0.753 
    ##               ACP2                 TF              NDRG2              SYT10 
    ##              0.753              0.753              0.753              0.753 
    ##               SOD2              HIP1R              KCNG1             EEF1B2 
    ##              0.753              0.753              0.753              0.753 
    ##               TMC6               MCM3            CCDC175             MYBPHL 
    ##              0.753              0.753              0.753              0.753 
    ##               TYW1               JUND            CCDC182             CCDC39 
    ##              0.754              0.754              0.754              0.754 
    ##              ATAD1               GUSB              CFDP1             LRSAM1 
    ##              0.754              0.754              0.754              0.754 
    ##              HLA-A            ZNF354B             ANAPC5         AP000783.1 
    ##              0.754              0.754              0.754              0.754 
    ##              RAB3B              NACC2             OGFRL1              STAC3 
    ##              0.754              0.754              0.754              0.754 
    ##              PTGFR           KIAA2013               SLU7             IFT172 
    ##              0.754              0.755              0.755              0.755 
    ##               RTKN             FAM50B             CCDC71            WFIKKN2 
    ##              0.755              0.755              0.755              0.755 
    ##              ALDOA               NOS2              CLDN3              EFNA3 
    ##              0.755              0.755              0.755              0.756 
    ##               ING1             SEMA4G               OAZ2              MTMR6 
    ##              0.756              0.756              0.756              0.756 
    ##            ST3GAL5              HTRA2         AC079354.1              CCNG2 
    ##              0.756              0.756              0.756              0.756 
    ##               CRYM             ITGA2B              ABTB1             CLDN16 
    ##              0.756              0.756              0.756              0.757 
    ##             BNIP3L            TXNDC15              TTLL5              DMRT2 
    ##              0.757              0.757              0.757              0.757 
    ##              TRAF1             ZNF365              HOXB9              WDR63 
    ##              0.757              0.757              0.757              0.757 
    ##              WDR70              SPC25              LPPR3             FAM21A 
    ##              0.757              0.757              0.758              0.758 
    ##             FAM78B               POP1              ZFP64              CEP76 
    ##              0.758              0.758              0.758              0.758 
    ##              MUC16              ARL5C             DYNLL1              NOBOX 
    ##              0.759              0.759              0.759              0.759 
    ##              TOR1B             GTF2F1              RHPN2              SAP30 
    ##              0.759              0.759              0.759              0.759 
    ##              PAQR8              ACSS1             ZNF317               SCP2 
    ##              0.759              0.759              0.759              0.760 
    ##         AC073610.5             MINPP1               AUP1            SLC7A11 
    ##              0.760              0.760              0.761              0.761 
    ##              RAB28             C4orf3                CFI              EVI2B 
    ##              0.761              0.761              0.761              0.761 
    ##               EVI5            NOXRED1              CIAO1            PSTPIP1 
    ##              0.761              0.761              0.761              0.761 
    ##              SPAG1                SP6              AGBL4       C8orf44-SGK3 
    ##              0.761              0.761              0.761              0.761 
    ##             KCNJ10               RNMT               SGK3              GLRA1 
    ##              0.761              0.761              0.761              0.762 
    ##             PRKCSH              RAB13               HELQ       RP11-794P6.2 
    ##              0.762              0.762              0.762              0.762 
    ##               GPN1              RSPO2             HOXD10            TP53I11 
    ##              0.762              0.762              0.762              0.762 
    ##           SLC39A12               SRGN              C1QL3             CCDC36 
    ##              0.763              0.763              0.763              0.763 
    ##               HARS             ABHD12              RPAP1           C18orf63 
    ##              0.763              0.763              0.763              0.763 
    ##               M6PR                TTR                MPL             ZSWIM3 
    ##              0.763              0.763              0.763              0.763 
    ##             IFNLR1              MKI67       RP11-298I3.5             CCDC25 
    ##              0.763              0.764              0.764              0.764 
    ##            SLC19A2               CNN2               SDF4              POLE2 
    ##              0.764              0.764              0.764              0.764 
    ##               FCF1             STAMBP            CCDC149             CDC14A 
    ##              0.764              0.764              0.765              0.765 
    ##              CABP7               MATK              ABCC8           ATP6V1C1 
    ##              0.765              0.765              0.765              0.765 
    ##              GRAP2               COMP             GPRC5B             SLC9C1 
    ##              0.765              0.765              0.765              0.765 
    ##             TMEM51           PHOSPHO1              SNTB2              PRRT3 
    ##              0.765              0.765              0.765              0.765 
    ##            SERINC5              MED20             EFCAB8              CSF3R 
    ##              0.765              0.766              0.766              0.766 
    ##            SPATA22               EXT2               RARS               ANO3 
    ##              0.766              0.766              0.766              0.766 
    ##               ODF3            PLEKHG4              ZNF18               EAF1 
    ##              0.766              0.766              0.766              0.766 
    ##             BARHL1            FAM219A               LRP3               UBL7 
    ##              0.767              0.767              0.767              0.767 
    ##             IGDCC3              FCGBP               TFR2            SYNJ2BP 
    ##              0.767              0.767              0.767              0.767 
    ##            CNTNAP3              APOBR              LMCD1             SHCBP1 
    ##              0.767              0.767              0.767              0.767 
    ##              ITIH4              NRDE2             CHRNA7               DBX1 
    ##              0.768              0.768              0.768              0.768 
    ##             PARP12               RINL             NKAIN2              DHX34 
    ##              0.768              0.768              0.768              0.768 
    ##               IARS             PIWIL1              YIPF4            SYNDIG1 
    ##              0.768              0.768              0.768              0.768 
    ##               RHCE             MAP3K6              SFXN5              DDHD2 
    ##              0.768              0.768              0.769              0.769 
    ##              CLCN2               PIGS              PAGE4           TMEM150A 
    ##              0.769              0.769              0.769              0.769 
    ##              RBM4B            PLA2G4E               CLK1              ZFPL1 
    ##              0.769              0.769              0.769              0.769 
    ##            DYNC1I2            C5orf22              VAMP4             CRTAC1 
    ##              0.769              0.769              0.770              0.770 
    ##              RPAP2            TMEM30A              YIPF1               LGI3 
    ##              0.770              0.770              0.770              0.770 
    ##              SUSD4              EFR3A             INSIG2             SLC5A1 
    ##              0.770              0.770              0.771              0.771 
    ##              DEDD2              LIMS1               DMC1            CCDC130 
    ##              0.771              0.771              0.771              0.771 
    ##              ATMIN             MARCH1            SLC23A3            TBC1D16 
    ##              0.771              0.771              0.771              0.771 
    ##              TNNI2            MOV10L1              REXO2              PAPD4 
    ##              0.771              0.771              0.771              0.772 
    ##              PRRT1               PLD3             IL17RE               MXI1 
    ##              0.772              0.772              0.772              0.772 
    ##              NUAK2             SOWAHC            ZKSCAN4             ZNF830 
    ##              0.772              0.772              0.772              0.772 
    ##               NET1             TSPAN3            PRPSAP2             NAP1L5 
    ##              0.772              0.772              0.772              0.772 
    ##              CMSS1            FAM204A               PFAS              ACPL2 
    ##              0.773              0.773              0.773              0.773 
    ##              EPB42            FAM102B              GRWD1              KCNA7 
    ##              0.773              0.773              0.773              0.773 
    ##            PCDHGB6             UNC13D            SLC7A10               PIGK 
    ##              0.773              0.773              0.774              0.774 
    ##            APCDD1L              KCNS2               SRPX              BUD13 
    ##              0.774              0.774              0.774              0.774 
    ##           C1orf101             SPACA1              BASP1           FAM160B2 
    ##              0.774              0.774              0.774              0.775 
    ##           B4GALNT3              CENPT                F12            ZNF385B 
    ##              0.775              0.775              0.775              0.775 
    ##           C19orf57             PFKFB2            ANAPC15              SYT16 
    ##              0.775              0.775              0.775              0.775 
    ##             RNF122             VPS37B               AZI1              NLRP4 
    ##              0.776              0.776              0.776              0.776 
    ##            ALDH1L1             GGNBP1               CHGB            L3MBTL1 
    ##              0.776              0.776              0.776              0.777 
    ##              ZGPAT             ANGEL1              OR6C1               NBAS 
    ##              0.777              0.777              0.777              0.777 
    ##              ATP5O            ARL6IP1              EDEM2               FTH1 
    ##              0.777              0.777              0.777              0.777 
    ##             SLC2A4              GTSF1               AVIL              SFXN2 
    ##              0.777              0.778              0.778              0.778 
    ##             RNASEK           HS3ST3B1              HTR1E        AP000304.12 
    ##              0.778              0.778              0.778              0.778 
    ##              TEX15              DHODH               CIPC             GRPEL1 
    ##              0.778              0.778              0.778              0.778 
    ##              PHGDH              MED17               PLEK             ZNF451 
    ##              0.778              0.778              0.778              0.778 
    ##            ZNF354A              RFWD3              RRP1B               CCR7 
    ##              0.778              0.778              0.778              0.778 
    ##               CDON               MARS               NEXN              USP53 
    ##              0.779              0.779              0.779              0.779 
    ##             ZDHHC6              DEGS1             CORO1B              EIF1B 
    ##              0.779              0.779              0.779              0.779 
    ##            SLC25A6              CETN2             LRRC71               ADGB 
    ##              0.779              0.779              0.779              0.780 
    ##               ATF5             GPR176               CINP            SLC46A1 
    ##              0.780              0.780              0.780              0.780 
    ##               SYT6             ZNF33B             MRPL27              ZMAT1 
    ##              0.780              0.780              0.780              0.780 
    ##               IRX2               CTSL                TK1             PSMD10 
    ##              0.780              0.781              0.781              0.781 
    ##            RAB3IL1             RASAL3               NARS               JPH2 
    ##              0.781              0.781              0.781              0.781 
    ##             CC2D2A              PSMB8             ZNF428                TXK 
    ##              0.781              0.781              0.781              0.781 
    ##              FOXD1              TRAT1              TMTC4             PCMTD2 
    ##              0.781              0.782              0.782              0.782 
    ##              CAPN5               CAST            EFCAB12              WDPCP 
    ##              0.782              0.782              0.782              0.782 
    ##              ZBED1                XPC            C2orf44           C16orf45 
    ##              0.782              0.782              0.782              0.782 
    ##                OS9            TSPAN15            CCDC108              PTGS1 
    ##              0.782              0.782              0.783              0.783 
    ##               TPM1               TPP1             BTBD17              XRRA1 
    ##              0.783              0.783              0.783              0.783 
    ##             SAMD9L            CYP4F22              TTC37            DNAJB13 
    ##              0.783              0.783              0.783              0.783 
    ##              DTX3L            DNAJC16              PAQR9              TOR2A 
    ##              0.783              0.783              0.783              0.783 
    ##               PAK4            SLC35D1              SLMO2              TIMP1 
    ##              0.784              0.784              0.784              0.784 
    ##               TAF6               GJA4              TTC27            SDCCAG8 
    ##              0.784              0.784              0.784              0.784 
    ##              SCN7A              DCTN6              BCCIP                GSS 
    ##              0.784              0.784              0.784              0.785 
    ##             MAMDC2           CCDC144A             IGSF11             ADRA1B 
    ##              0.785              0.785              0.785              0.785 
    ##                CKM             ITGA10            TMEM182               BRK1 
    ##              0.785              0.785              0.785              0.785 
    ##              BRIP1              ERCC3             OPN1SW            PRPS1L1 
    ##              0.786              0.786              0.786              0.786 
    ##               YOD1             GPR161              HJURP           KIAA0226 
    ##              0.786              0.786              0.786              0.787 
    ##               BTLA             ZNF414             ATP10D               TAP2 
    ##              0.787              0.787              0.787              0.787 
    ##             L2HGDH             FIGNL1             INPP5K               EMP1 
    ##              0.787              0.787              0.787              0.787 
    ##              GPSM3             TBC1D2              CNGA2               EHD4 
    ##              0.787              0.787              0.787              0.787 
    ##               MED4              FOXC2              NAIF1                 TH 
    ##              0.788              0.788              0.788              0.788 
    ##              BTBD1              SFXN1             ZDHHC7                AK4 
    ##              0.788              0.788              0.788              0.788 
    ##           SLC4A1AP      RP11-108K14.8              ZNF45               SGCA 
    ##              0.788              0.789              0.789              0.789 
    ##             SH3GL2            RNPEPL1              DISC1              CWC25 
    ##              0.789              0.789              0.789              0.789 
    ##                MBP              KCTD5              GPR37                AMT 
    ##              0.789              0.789              0.789              0.789 
    ##              NXPH3               CA11               DSG2               DTNB 
    ##              0.789              0.790              0.790              0.790 
    ##           MIS18BP1               WEE2               PI15              CAPN7 
    ##              0.790              0.790              0.790              0.790 
    ##             TRIM47              KIF15             MAPK11              ARMC2 
    ##              0.790              0.790              0.790              0.790 
    ##               THBD           C17orf99            FAM122A            MAGEA11 
    ##              0.791              0.791              0.791              0.791 
    ##              MKNK2              COX5A             NKAIN3              GFOD1 
    ##              0.791              0.791              0.791              0.791 
    ##              CHID1              NSUN7             KCTD20               RAG2 
    ##              0.791              0.791              0.792              0.792 
    ##             DCBLD1             PTP4A1           TMEM167B             TESPA1 
    ##              0.792              0.792              0.792              0.792 
    ##             ZNF142             PAPSS1             TUBA1C              RPL12 
    ##              0.792              0.792              0.792              0.793 
    ##           B4GALNT1             ZNF501           C1orf222              PSMB4 
    ##              0.793              0.793              0.793              0.793 
    ##             GAS2L2            SAP30BP              AIM1L            FAM149A 
    ##              0.793              0.793              0.793              0.793 
    ##              CDH24              STK10             CGRRF1               MGLL 
    ##              0.794              0.794              0.794              0.794 
    ##               TTC5             PCDH15               OSR2             RCBTB1 
    ##              0.794              0.794              0.794              0.794 
    ##            CCDC155             SSX2IP              CDIPT             METTL3 
    ##              0.794              0.794              0.794              0.794 
    ##              ASXL1             MRPL40             CEP120               CTSD 
    ##              0.794              0.794              0.795              0.795 
    ##              EXOC4       TMED7-TICAM2             SP140L             DEPDC1 
    ##              0.795              0.795              0.795              0.795 
    ##            RAD21L1                AGL            C8orf46             KLHL41 
    ##              0.795              0.795              0.795              0.795 
    ##              DRAP1             IQGAP2                TRO             KCNK13 
    ##              0.795              0.795              0.795              0.795 
    ##             TICAM2              VWC2L              CHST3              POMT2 
    ##              0.795              0.795              0.795              0.795 
    ##              NLRP5            SLC35A3            FAM186A              CLIP4 
    ##              0.795              0.795              0.796              0.796 
    ##               TPRN              RPS28               MFN1               OLR1 
    ##              0.796              0.796              0.796              0.796 
    ##               HNMT             TSEN15             PHOX2A             STRADA 
    ##              0.796              0.796              0.796              0.796 
    ##              EIF2D            ZDHHC19             ZNF839              TEAD2 
    ##              0.796              0.796              0.796              0.796 
    ##               DOK3               VWCE            TMEM104              TGFBI 
    ##              0.796              0.796              0.796              0.796 
    ##               GCM2                PGF              IFT81              MATN1 
    ##              0.796              0.797              0.797              0.797 
    ##        RBAK-RBAKDN             ZNF460             RNF114               TARS 
    ##              0.797              0.797              0.797              0.797 
    ##             WRAP73               CDC7             SCNN1B            TMEM38B 
    ##              0.797              0.797              0.798              0.798 
    ##            VIPAS39               ZPBP              ARL5A               NEK8 
    ##              0.798              0.798              0.798              0.798 
    ##              DNAH3             NLRP11             IMPDH2                ADK 
    ##              0.798              0.798              0.798              0.798 
    ##                MX2           SLC39A13             SLAMF1            TNFSF14 
    ##              0.798              0.799              0.799              0.799 
    ##               SNX4               CTNS              ITGA7             MRPS23 
    ##              0.799              0.799              0.799              0.799 
    ##            SLC52A3              LIN7C               MOB2            SLC39A6 
    ##              0.800              0.800              0.800              0.800 
    ##               GCC1               GFM1              SPNS2              WDHD1 
    ##              0.800              0.800              0.800              0.800 
    ##             CHMP2A            CCDC124              CARKD             MYO15A 
    ##              0.800              0.800              0.800              0.800 
    ##             NDUFC2               ELF5               GPX3              CD274 
    ##              0.800              0.801              0.801              0.801 
    ##            CCNDBP1           KIAA0408              RAP2B               TLR3 
    ##              0.801              0.801              0.801              0.801 
    ##             FAM35A                GBA              KANK4           PCDHGA12 
    ##              0.801              0.801              0.801              0.801 
    ##             B3GNT1             PRKAG1              SCN3B            C2orf81 
    ##              0.801              0.801              0.801              0.801 
    ##             RNF130               SIX5            SLC39A1           TMEM132A 
    ##              0.801              0.802              0.802              0.802 
    ##             SQSTM1             SORBS3              ERCC6              RIOK3 
    ##              0.802              0.802              0.802              0.802 
    ##              EFNA3              PDILT             ALKBH4               TNK2 
    ##              0.802              0.802              0.802              0.802 
    ##               LIG4             ZSCAN2               CT55            ABHD17A 
    ##              0.802              0.802              0.802              0.803 
    ##                RP9             SNRPB2           C12orf50              RNFT2 
    ##              0.803              0.803              0.803              0.803 
    ##             CHST15              TUBG2              SH2B3             SH3GL3 
    ##              0.803              0.803              0.803              0.803 
    ##            TMEM179             CLEC5A            BCL2L13            TMEM214 
    ##              0.803              0.803              0.803              0.804 
    ##           C10orf71               HEY1                 F3             ZNF397 
    ##              0.804              0.804              0.804              0.804 
    ##               SGSH              KLHL6              STAB1            C2orf69 
    ##              0.804              0.805              0.805              0.805 
    ##            ZDHHC18             MTHFD2              RBMS2              ERCC5 
    ##              0.805              0.805              0.805              0.805 
    ##     RP11-977G19.10               TBCE            SLC26A7              PCDH8 
    ##              0.805              0.805              0.805              0.806 
    ##             SPOCD1               FGD2               TAP2              PDAP1 
    ##              0.806              0.806              0.806              0.806 
    ##               BYSL               IGF1            UBE2QL1                GIP 
    ##              0.806              0.806              0.806              0.806 
    ##              ZNF70               QARS              RPL36               DRD3 
    ##              0.806              0.806              0.806              0.807 
    ##              WIPF3            SLC33A1              DNHD1               POLB 
    ##              0.807              0.807              0.807              0.807 
    ##              TESK2               MDH2             RNF115               DPP3 
    ##              0.807              0.807              0.807              0.807 
    ##              CHMP5           CRISPLD1             NDUFS1             NT5DC3 
    ##              0.807              0.808              0.808              0.808 
    ##              PDE4C             GAPDHS              NXNL2            SUPT20H 
    ##              0.808              0.808              0.808              0.808 
    ##             GABRG1                NOV           SLC25A22              ANXA2 
    ##              0.808              0.808              0.808              0.808 
    ##             FBXO43               CCNH               KPTN              MFSD5 
    ##              0.809              0.809              0.809              0.809 
    ##             CCDC91           KIAA1147            TMEM145               AASS 
    ##              0.809              0.809              0.809              0.809 
    ##            SPATA17              INSM1              MAT1A              RBBP9 
    ##              0.809              0.809              0.809              0.810 
    ##                ADO              NOMO3              LIN37                IL2 
    ##              0.810              0.810              0.810              0.810 
    ##             LRRC73             PAPOLB             AMDHD2             ZNF300 
    ##              0.810              0.810              0.810              0.810 
    ##              CHST2             RPP25L              CASP8             EIF4E3 
    ##              0.811              0.811              0.811              0.811 
    ##               GAS2               HRH2              SRSF2             CDC123 
    ##              0.811              0.811              0.811              0.811 
    ##              GSTCD             GTF2H4              TEAD4               PIGO 
    ##              0.811              0.811              0.811              0.812 
    ##              CEP97              KCNU1             GPR133              PPIL1 
    ##              0.812              0.812              0.812              0.812 
    ##             ATP5C1              LRCH3              SPG11               PDK1 
    ##              0.812              0.812              0.812              0.812 
    ##              WNT11           ARHGAP40             CEP135             DIRAS2 
    ##              0.812              0.812              0.813              0.813 
    ##             SOHLH1              BRCC3             OR12D2               GJA5 
    ##              0.813              0.813              0.813              0.813 
    ##               IL6R             SH3TC2             PRPF18             PRSS16 
    ##              0.813              0.813              0.814              0.814 
    ##             COL9A1              PARP9            PIP5K1B             KCTD17 
    ##              0.814              0.814              0.814              0.814 
    ##               CRY1              KCNQ1               DAW1             RAVER2 
    ##              0.814              0.814              0.814              0.814 
    ##              ACAD9              COPG2             UBE2D1              IL12A 
    ##              0.814              0.815              0.815              0.815 
    ##             ATP12A              WDR35              RBM28              ALPK2 
    ##              0.815              0.815              0.815              0.815 
    ##              FKBP9              PLCD3              COX7B         ST6GALNAC6 
    ##              0.815              0.815              0.815              0.815 
    ##                AMN             ENTPD6              GPR34            COLEC11 
    ##              0.816              0.816              0.816              0.816 
    ##            RABGGTA            CTNNAL1               GRK1             IGFBP2 
    ##              0.816              0.816              0.816              0.816 
    ##              ATAT1            KATNAL2               IRGC              PFDN4 
    ##              0.816              0.816              0.816              0.816 
    ##                DBT             ZNF408             BDKRB2            CYP17A1 
    ##              0.816              0.816              0.817              0.817 
    ##           PCDHGA10              RWDD4              DDX31               CD63 
    ##              0.817              0.817              0.817              0.817 
    ##             VSIG10              SYPL1              GREM1              CAMK1 
    ##              0.817              0.818              0.818              0.818 
    ##              VPS16               FDXR             ZNF565             IFT140 
    ##              0.818              0.818              0.818              0.818 
    ##            TXNDC11               NEK2             SPINT1            CYP21A2 
    ##              0.818              0.818              0.818              0.818 
    ##              PLOD3             GPR179           C19orf43               MSX2 
    ##              0.819              0.819              0.819              0.819 
    ##              CBWD2             AMICA1            HSD11B2               TMX1 
    ##              0.819              0.819              0.819              0.819 
    ##             ZBTB32              ABCG4               NPR3              GFOD2 
    ##              0.819              0.819              0.819              0.819 
    ##             CERCAM               TTI1              FAIM2            CCDC101 
    ##              0.819              0.820              0.820              0.820 
    ##             ENOPH1              WNT7B            GABARAP            PPIP5K2 
    ##              0.820              0.820              0.820              0.820 
    ##             TANGO6             DCAF17             LCLAT1               DGKG 
    ##              0.820              0.820              0.820              0.820 
    ##              AGFG2               DSG4               ISPD               ASPH 
    ##              0.821              0.821              0.821              0.821 
    ##            ZNF585B               BEX4            SLC16A1               NMD3 
    ##              0.821              0.821              0.821              0.821 
    ##              THOC7               VWA9             GNPTAB           TMEM161A 
    ##              0.821              0.821              0.821              0.821 
    ##               PCNT              RBM12               TLX3             DNAJC1 
    ##              0.821              0.822              0.822              0.822 
    ##               IFNG            PLEKHA4             CHRNA4             CHRNB2 
    ##              0.822              0.822              0.822              0.822 
    ##             RNF150              H2AFV              MFSD6            SLC38A7 
    ##              0.822              0.822              0.822              0.822 
    ##              EIF3K              UBOX5               AREG            GALNTL6 
    ##              0.822              0.823              0.823              0.823 
    ##              ARVCF           RASGEF1C                GIF             PGRMC1 
    ##              0.823              0.823              0.823              0.823 
    ##           C16orf52               SELT               SGCD             MFHAS1 
    ##              0.823              0.824              0.824              0.824 
    ##              TMED9            PDE4DIP              RAB15            SLC12A1 
    ##              0.824              0.824              0.824              0.824 
    ##              TRAK2              CD247              STRA6              WDR77 
    ##              0.824              0.824              0.824              0.824 
    ##             CEP152             CYP7A1               LIPF     TMEM256-PLSCR3 
    ##              0.824              0.824              0.824              0.824 
    ##           SLC25A33             CACNG4            ALDH3A2             HEATR3 
    ##              0.825              0.825              0.825              0.825 
    ##               NOL8             TFIP11           ANKRD18B            FAM133B 
    ##              0.825              0.825              0.825              0.825 
    ##               GMNN           TRNAU1AP               ETFB               NCK1 
    ##              0.825              0.825              0.825              0.825 
    ##               SDHB              PEX16             KBTBD8            SLC6A12 
    ##              0.825              0.825              0.825              0.825 
    ##               TMC7               HBG2             BCL2L2              CHPT1 
    ##              0.825              0.826              0.826              0.826 
    ##              RINT1              TLDC2               GSX2              LMBR1 
    ##              0.826              0.826              0.826              0.826 
    ##              PAICS               GLTP              KRT85             ATP5G3 
    ##              0.826              0.826              0.826              0.826 
    ##            EMILIN3              AP4B1             POU2F3                ID2 
    ##              0.827              0.827              0.827              0.827 
    ##      RP11-302B13.5               WNT6              SRPRB              JOSD1 
    ##              0.827              0.827              0.827              0.827 
    ##               TPM2             ATAD3A              SNX14              CPT1B 
    ##              0.827              0.827              0.827              0.827 
    ##            SLC16A7             ZNF775             KIF16B               PUS7 
    ##              0.827              0.827              0.827              0.827 
    ##             LILRB3             NECAB3              AASDH              DTHD1 
    ##              0.827              0.827              0.828              0.828 
    ##              ASIC3           C17orf47             SCN11A                SRR 
    ##              0.828              0.828              0.828              0.828 
    ##            CSRP2BP              UFSP2              KNTC1           C15orf26 
    ##              0.828              0.828              0.828              0.828 
    ##            LDLRAD3               APEH             CNTROB               CRB1 
    ##              0.829              0.829              0.829              0.829 
    ##             DNAJB4             KCNJ15              BANK1             FAM13A 
    ##              0.829              0.829              0.829              0.829 
    ##               TUT1            ANKRD54               TLX2              GPR21 
    ##              0.829              0.829              0.829              0.830 
    ##             TMEM66             LRRC39               ALG8              CCDC8 
    ##              0.830              0.830              0.830              0.830 
    ##             MPPED1             C2CD4C            MSANTD3 CSNK2B-LY6G5B-1181 
    ##              0.830              0.830              0.830              0.830 
    ##             MAN1B1              COX11              PRRG4              MIEF2 
    ##              0.830              0.830              0.830              0.831 
    ##                CLU              ECEL1            TMEM117               GLS2 
    ##              0.831              0.831              0.831              0.831 
    ##              KCNK4             FBXL12              SERP2              SHOX2 
    ##              0.831              0.831              0.831              0.831 
    ##              CDR2L            COL26A1               LIX1             HIBADH 
    ##              0.832              0.832              0.832              0.832 
    ##            DENND6B              COPRS              YPEL5              RPLP2 
    ##              0.832              0.832              0.832              0.832 
    ##             SPTLC3                LEP            C2orf78          OVCH1-AS1 
    ##              0.832              0.832              0.832              0.832 
    ##            SLC35F5              TAF11              KIFC2               IL4R 
    ##              0.832              0.832              0.832              0.832 
    ##              CDKN3             IL18R1              GPR65            TSNARE1 
    ##              0.833              0.833              0.833              0.833 
    ##             SLC1A1               SYF2              ESCO2             ATP2A1 
    ##              0.833              0.833              0.833              0.833 
    ##              RAB23           ARHGAP24              GPLD1             ERLIN2 
    ##              0.833              0.833              0.833              0.833 
    ##              PGAP3              SPG21               GPD1              PSMB6 
    ##              0.833              0.834              0.834              0.834 
    ##               DOK5               FZD9             INPP5B              SNRPN 
    ##              0.834              0.834              0.834              0.834 
    ##            TMEM63A           ADAMTS20            SLC13A3             MGAT4C 
    ##              0.834              0.834              0.834              0.834 
    ##             AKT1S1             UQCRC2                EYS             NDUFV2 
    ##              0.834              0.834              0.834              0.834 
    ##             TSPAN8              AARS2              TTC14               LAT2 
    ##              0.835              0.835              0.835              0.835 
    ##               BCAM              IFIT2            TRAPPC3            SLC26A8 
    ##              0.835              0.835              0.835              0.835 
    ##             TMEM65             PDLIM2               SYT5              MTIF2 
    ##              0.835              0.836              0.836              0.836 
    ##            CIAPIN1            ARHGAP9            PCDHGA9              NLRP1 
    ##              0.836              0.836              0.836              0.836 
    ##              RAB12               POLH         AC104809.3              RSRC1 
    ##              0.836              0.836              0.836              0.836 
    ##             MAMSTR            SLC16A6             PPHLN1              SNURF 
    ##              0.836              0.836              0.836              0.836 
    ##               KDM8              PNMA1                GSN             TATDN2 
    ##              0.836              0.836              0.837              0.837 
    ##               PHKB             CEP164              PARP4               VAT1 
    ##              0.837              0.837              0.837              0.837 
    ##              DNAH7           HLA-DQB2               HHAT               LSM7 
    ##              0.837              0.837              0.837              0.837 
    ##              DUSP3              MLYCD               SNX3             ISYNA1 
    ##              0.837              0.837              0.837              0.837 
    ##              PDIA5           C15orf59           KIAA1524              NPY1R 
    ##              0.837              0.837              0.837              0.837 
    ##             LANCL3                MCC             RNF180              SPG20 
    ##              0.837              0.837              0.837              0.838 
    ##            SLCO4C1               TEP1              KRT13             RTFDC1 
    ##              0.838              0.838              0.838              0.838 
    ##              QSOX2             DCTPP1              CDC20              MCM10 
    ##              0.838              0.838              0.838              0.838 
    ##               NEU1               ENAM              LALBA             RAB40C 
    ##              0.838              0.838              0.838              0.838 
    ##              ABCA5              EEF1D             MRE11A               BOD1 
    ##              0.838              0.838              0.839              0.839 
    ##             IMPAD1              PRDM5        ERCC6-PGBD3               ACCS 
    ##              0.839              0.839              0.839              0.839 
    ##              NARG2             CEP128            CWF19L2             ZNF350 
    ##              0.839              0.839              0.839              0.839 
    ##           CDC42BPG                GGN              CARS2            GUCY1A3 
    ##              0.839              0.839              0.840              0.840 
    ##              PLAUR             TM7SF3              VEGFA              VPS41 
    ##              0.840              0.840              0.840              0.840 
    ##            C1orf43               GIN1              MMP11              GDF10 
    ##              0.840              0.840              0.840              0.840 
    ##               MGAM              NUP43              MGME1              TARS2 
    ##              0.840              0.840              0.840              0.841 
    ##             VSTM2L               PER3            FAM184B             ARMCX3 
    ##              0.841              0.841              0.841              0.841 
    ##             TARBP1              FBXO6               MPP2               EBI3 
    ##              0.841              0.841              0.841              0.841 
    ##             PI4K2A              LAMC3              RFTN2             MUC5AC 
    ##              0.841              0.841              0.841              0.842 
    ##            ZCCHC17             TMEM33               BBS9             EFCAB5 
    ##              0.842              0.842              0.842              0.842 
    ##             STXBP2              DERL1           C6orf132            CABLES2 
    ##              0.842              0.842              0.842              0.842 
    ##               FKTN                POR              FOXH1              LARP7 
    ##              0.842              0.842              0.842              0.842 
    ##            PIK3IP1              PUS7L           C11orf68               ZIM2 
    ##              0.842              0.843              0.843              0.843 
    ##            SLC45A3              TROAP              WISP3             PRKAB1 
    ##              0.843              0.843              0.843              0.843 
    ##               UFL1           B3GALNT2               JAM3            FAM104A 
    ##              0.843              0.843              0.843              0.843 
    ##              KRT6A             ADSSL1              MYO5C               CUTC 
    ##              0.843              0.843              0.844              0.844 
    ##      RP11-192H23.4               LIAS           KIAA0930             NIPAL3 
    ##              0.844              0.844              0.844              0.844 
    ##               FHL3                FES             CCDC40              STRA8 
    ##              0.845              0.845              0.845              0.845 
    ##             TMEM19            NR2C2AP               URB2            C9orf84 
    ##              0.845              0.845              0.845              0.846 
    ##           CCDC102A             SLC9B2              KRT14             PDXDC1 
    ##              0.846              0.846              0.846              0.846 
    ##               ETV4              FXYD1             ANKLE2              NUGGC 
    ##              0.846              0.846              0.847              0.847 
    ##            TMEM59L              VSIG8             TADA2A              BTNL3 
    ##              0.847              0.847              0.847              0.847 
    ##              DDX53            SPATA2L              SMYD2               TTC6 
    ##              0.847              0.847              0.847              0.847 
    ##              VPS72              PRRT4               TLR6               RGS9 
    ##              0.847              0.847              0.847              0.847 
    ##               MTL5           LRRC37A3              FRRS1             GTF2F2 
    ##              0.848              0.848              0.848              0.848 
    ##               DDB2              PTGDS               MECR          TNFRSF10B 
    ##              0.848              0.848              0.848              0.848 
    ##              HOXD8             IQGAP3               JRKL              MYO7A 
    ##              0.848              0.848              0.848              0.848 
    ##              PTAFR              UBE2C              OLIG1              MFAP5 
    ##              0.848              0.848              0.849              0.849 
    ##             CCDC12              CDHR2              CCSAP              CDC26 
    ##              0.849              0.849              0.849              0.850 
    ##             FKBP10             CTNNA3             ZNF667            TMEM249 
    ##              0.850              0.850              0.850              0.850 
    ##              PTCH2                FTO              CD101              WDR17 
    ##              0.850              0.850              0.850              0.850 
    ##            L3MBTL4                CA3              GPR45                UNG 
    ##              0.850              0.850              0.850              0.850 
    ##             HHIPL1              CPLX3               LIPE              MTCP1 
    ##              0.851              0.851              0.851              0.851 
    ##              MROH8           KIAA1644             RNF149            TMEM163 
    ##              0.851              0.851              0.851              0.851 
    ##               TPH2              KIF22           ARHGEF25            SLC35G2 
    ##              0.851              0.851              0.851              0.851 
    ##              WDR25            C4orf32              HIPK4             DCAF13 
    ##              0.852              0.852              0.852              0.852 
    ##              SMAD9              NUCB1              BNIP2             CCHCR1 
    ##              0.852              0.852              0.852              0.852 
    ##             FAM47E               IPO4             TTC21B              DDX56 
    ##              0.852              0.852              0.852              0.852 
    ##               MICA               GPN2             CYB5R3                AK7 
    ##              0.852              0.852              0.853              0.853 
    ##               SCG2              GSTP1             MANEAL            C14orf1 
    ##              0.853              0.853              0.853              0.853 
    ##               PLAU             TMEM35            CCDC104              WDR76 
    ##              0.853              0.853              0.853              0.853 
    ##             DNAJB1            TNFSF11             ERO1LB               F13B 
    ##              0.853              0.853              0.853              0.853 
    ##              PLCB2              TCEA2             FERMT1              IFT20 
    ##              0.853              0.854              0.854              0.854 
    ##             MALRD1               CDH3             TRIM32               TWF1 
    ##              0.854              0.854              0.855              0.855 
    ##              LIMA1            ST8SIA5              SHMT2               CYGB 
    ##              0.855              0.855              0.855              0.855 
    ##             ZBTB39               PUS1               SOX3              PROM2 
    ##              0.855              0.855              0.855              0.855 
    ##           ANKRD18A             FAM21C              NICN1              OSTF1 
    ##              0.856              0.856              0.856              0.856 
    ##             CCDC82            SLC43A3               MYH6           STARD3NL 
    ##              0.856              0.856              0.856              0.856 
    ##         AC004381.6               OCA2             MFSD11           ARHGAP28 
    ##              0.856              0.856              0.856              0.856 
    ##              LMOD2              SIRT7              SNX19           CATSPERG 
    ##              0.857              0.857              0.857              0.857 
    ##              NOC3L             RAB27B               GCGR              LARS2 
    ##              0.857              0.857              0.857              0.857 
    ##            STK11IP            TBC1D31              ACACB             GOLT1A 
    ##              0.857              0.857              0.857              0.857 
    ##              CRTAP            CDK2AP1              USH2A                XDH 
    ##              0.857              0.857              0.857              0.857 
    ##            CYP26C1             HAVCR2               TBL2            TMEM45A 
    ##              0.857              0.857              0.858              0.858 
    ##              ACTG1           C16orf47               ARSD              DDAH2 
    ##              0.858              0.858              0.858              0.858 
    ##              VAMP7             IFT122             SNAPC2             PCYT1A 
    ##              0.858              0.858              0.858              0.858 
    ##             DCAF12              WDR62               SDHA           CRISPLD2 
    ##              0.858              0.859              0.859              0.859 
    ##            ZMYND10                RAX              STK33               ARG2 
    ##              0.859              0.859              0.859              0.859 
    ##           ARHGAP18            FAM19A5            CCDC135             KCNJ13 
    ##              0.859              0.859              0.859              0.859 
    ##             ZNF639               GPD2              PEX10              BCAS1 
    ##              0.860              0.860              0.860              0.860 
    ##             VPS37C      NDUFC2-KCTD14               GLB1              PAM16 
    ##              0.860              0.860              0.860              0.860 
    ##              LOXL3               NEK9             PRSS54             MRPL10 
    ##              0.860              0.860              0.860              0.860 
    ##               SELE              VTI1A                CD9               VIL1 
    ##              0.860              0.860              0.860              0.861 
    ##               ORC5              CLVS2               GPR4               NXF2 
    ##              0.861              0.861              0.861              0.861 
    ##             CREBZF              PARK2            FAM117A             SPIRE2 
    ##              0.861              0.861              0.861              0.861 
    ##             SNAP29              STK36               CDS2             NSMCE2 
    ##              0.861              0.862              0.862              0.862 
    ##                DBP             RAD54B           KIAA1875              AP1S1 
    ##              0.862              0.862              0.862              0.862 
    ##            CD163L1             LRRC55               ARL9             PLA2G6 
    ##              0.862              0.862              0.862              0.862 
    ##               NEK1               DLAT             NLGN4Y             LRRC3B 
    ##              0.862              0.863              0.863              0.863 
    ##             PINLYP           SLC38A10             UBE2G2            CCDC174 
    ##              0.863              0.863              0.863              0.863 
    ##              FOXS1             PCOLCE            FAM134B              FZD10 
    ##              0.863              0.863              0.863              0.863 
    ##               FGF1             GALNT3              MDFIC             SH2D4B 
    ##              0.863              0.863              0.863              0.863 
    ##              NARS2              TRPV2               ENO1              GPSM2 
    ##              0.863              0.863              0.864              0.864 
    ##             SPDYE4             ADAM15             RNF212               PROC 
    ##              0.864              0.864              0.864              0.864 
    ##               GALT               HN1L              ABCA6           C11orf82 
    ##              0.864              0.864              0.864              0.864 
    ##             FAM21B       FAM47E-STBD1             TUBA3C              FABP3 
    ##              0.864              0.865              0.865              0.865 
    ##             SCAMP4              MANEA              THAP7             HCRTR1 
    ##              0.865              0.865              0.865              0.865 
    ##             FAM20A              TBRG1             CARNS1              KLF10 
    ##              0.865              0.865              0.865              0.865 
    ##               PIGU              GFI1B             MCIDAS          C14orf105 
    ##              0.865              0.865              0.865              0.866 
    ##              IFT52              ITGAD             TCIRG1             OR52I2 
    ##              0.866              0.866              0.866              0.866 
    ##             LRRC31            ZSCAN12              N4BP3            EFCAB4B 
    ##              0.866              0.866              0.866              0.866 
    ##              BTBD9             RNF133               POLK              PDE6C 
    ##              0.866              0.866              0.866              0.866 
    ##             RNF219          ADCYAP1R1             DNAJC9              HTR2C 
    ##              0.866              0.866              0.867              0.867 
    ##               PIGH             ATP5F1               VWA7              TADA3 
    ##              0.867              0.867              0.867              0.867 
    ##            ANKRD55             IL17RC             TRIM14            C3orf58 
    ##              0.867              0.867              0.867              0.867 
    ##              BCKDK           CATSPERB             ZNF471              CMTM4 
    ##              0.867              0.867              0.867              0.867 
    ##               ZNF7             HIRIP3              BUD31          RAB11FIP1 
    ##              0.867              0.867              0.867              0.867 
    ##           TMEM150C               SYNC             GPR139               ALAD 
    ##              0.867              0.867              0.867              0.868 
    ##              FGFR4            SLC30A6             PLA2R1             THEMIS 
    ##              0.868              0.868              0.868              0.868 
    ##               ST7L              NXPH4                MDK              OPLAH 
    ##              0.868              0.868              0.868              0.868 
    ##             POLR3B            ANKRD31              GPR20              LPIN3 
    ##              0.868              0.868              0.868              0.869 
    ##               LTV1              LPAR2            SLC2A12              DDX55 
    ##              0.869              0.869              0.869              0.869 
    ##            SLC24A4              ACOX2             MRPS27              TTLL1 
    ##              0.870              0.870              0.870              0.870 
    ##                LSR            TMEM240           TRAF3IP3            ARHGDIA 
    ##              0.870              0.870              0.870              0.870 
    ##            NEUROD4              LCMT1               TBCK              ZNF91 
    ##              0.870              0.871              0.871              0.871 
    ##                ELN               STK3             IFRG15              MARS2 
    ##              0.871              0.871              0.871              0.871 
    ##              BANF1             ZNF862             ZNF699             ENTPD7 
    ##              0.872              0.872              0.872              0.872 
    ##              MTURN            COL23A1               FBF1           FGFR1OP2 
    ##              0.872              0.872              0.872              0.872 
    ##               SLX4             TIMM8A              SIDT1         AC037459.4 
    ##              0.872              0.872              0.872              0.872 
    ##               MCM9            SLC35D3              CH25H              TJAP1 
    ##              0.872              0.872              0.873              0.873 
    ##             PDLIM1              OBSCN               ZIC4              CASP7 
    ##              0.873              0.873              0.873              0.873 
    ##               DHPS               DLL3             ZNF529            SLC35C2 
    ##              0.873              0.873              0.873              0.873 
    ##             MRPS31               DPM1              RAD50              HGFAC 
    ##              0.873              0.873              0.873              0.874 
    ##             TRIM59             LEPROT             OR10A7               VSX1 
    ##              0.874              0.874              0.874              0.874 
    ##             MRPL48               PYGM              ALMS1             CCDC47 
    ##              0.874              0.874              0.874              0.875 
    ##              TIGD2               OTOF              MATN2              CEBPG 
    ##              0.875              0.875              0.875              0.875 
    ##               MPND              FHOD1              NPY2R             PARD6G 
    ##              0.875              0.876              0.876              0.876 
    ##             HS3ST6             SUCLG1              TRNT1             TRAFD1 
    ##              0.876              0.876              0.876              0.876 
    ##             POLR2F             KCTD15               HLCS             KBTBD3 
    ##              0.876              0.876              0.876              0.876 
    ##               PLD1             PCDHA1               KRI1              DDX59 
    ##              0.876              0.876              0.876              0.876 
    ##                NDP            RPGRIP1             POLR3A             NMNAT1 
    ##              0.876              0.876              0.876              0.877 
    ##            CYP19A1              IMPA1             LRRC45              PRDX3 
    ##              0.877              0.877              0.877              0.877 
    ##              FOCAD     PTGES3L-AARSD1            COL25A1               ZFP3 
    ##              0.877              0.877              0.877              0.877 
    ##               LMLN             SKIV2L            LAPTM4A              FBXO7 
    ##              0.877              0.877              0.877              0.877 
    ##             ARNTL2             LRRC10              CARD8              HSDL2 
    ##              0.878              0.878              0.878              0.878 
    ##             ZNF436              OBSL1            C8orf88               GRPR 
    ##              0.878              0.878              0.878              0.878 
    ##              PCGF6             GPRC5C            SLC37A4               STRC 
    ##              0.878              0.878              0.878              0.878 
    ##               HFE2             NHLRC1              UBXN6                AGK 
    ##              0.878              0.878              0.878              0.878 
    ##           COLGALT1              FREM3              EIF2A                VDR 
    ##              0.878              0.878              0.878              0.878 
    ##            XPNPEP3                CPM              SKAP1              MOB1B 
    ##              0.879              0.879              0.879              0.879 
    ##             GABPB2                CR2             PCDHA5               CTSA 
    ##              0.879              0.879              0.879              0.879 
    ##           SERPINB9             LRRC8E              TULP3             SDR9C7 
    ##              0.879              0.879              0.879              0.879 
    ##                GCG              NUDT5              IFT74           RNASEH2A 
    ##              0.880              0.880              0.880              0.880 
    ##                ACR               AIDA              CYB5B           C1orf228 
    ##              0.880              0.880              0.880              0.880 
    ##              LIMS2            TSC22D3             FBXL18              SARDH 
    ##              0.880              0.880              0.880              0.880 
    ##             OR51J1            PCDHGB3               A1CF              SYNPR 
    ##              0.880              0.880              0.880              0.880 
    ##            SDR16C5               AQP1              CMTM6             ERRFI1 
    ##              0.880              0.881              0.881              0.881 
    ##          C17orf103               ORC1            CYP11A1            KIRREL2 
    ##              0.881              0.881              0.881              0.881 
    ##              DCST2             STRIP2            ZKSCAN8            RSPH10B 
    ##              0.881              0.882              0.882              0.882 
    ##               MMAB             TSEN54               ELP2              SCRT1 
    ##              0.882              0.882              0.882              0.882 
    ##              STAB2              PARS2             KRBOX4               SNF8 
    ##              0.882              0.882              0.882              0.883 
    ##              LRRD1              THSD1           C14orf37               CHAT 
    ##              0.883              0.883              0.883              0.883 
    ##             ACAD11               CHDH             TTC21A              ITGAE 
    ##              0.883              0.883              0.883              0.884 
    ##               TTF2              DGCR2             KCTD19               RBP2 
    ##              0.884              0.884              0.884              0.884 
    ##                IYD              CKS1B               DHFR             DNASE2 
    ##              0.884              0.884              0.884              0.884 
    ##             SLC6A7              HKDC1             NIPAL4            SLC12A7 
    ##              0.884              0.884              0.884              0.885 
    ##               PPT1               TTC3            FAM184A              KIF24 
    ##              0.885              0.885              0.885              0.885 
    ##               TIE1             TANGO2             CYP2S1             TMEM98 
    ##              0.885              0.885              0.886              0.886 
    ##              CCKAR              ERMP1            ZMYND15            AFAP1L1 
    ##              0.886              0.886              0.886              0.886 
    ##             COMMD1              HHATL              PEBP1               DRG2 
    ##              0.886              0.886              0.886              0.886 
    ##             FANCD2               LRMP           SIGLEC14              GPR56 
    ##              0.886              0.886              0.887              0.887 
    ##               COQ2             MAATS1              HERC6              LHCGR 
    ##              0.887              0.887              0.887              0.887 
    ##             RAB3IP              CHMP3            CCDC85A               RHOG 
    ##              0.887              0.887              0.887              0.887 
    ##              PPIL4       RNF103-CHMP3               CRB2            TMEM169 
    ##              0.887              0.887              0.888              0.888 
    ##              USH1C              DHX32             NT5DC1            SLC15A1 
    ##              0.888              0.888              0.888              0.888 
    ##             SPTSSA      MEF2BNB-MEF2B             CHRDL2              PEAR1 
    ##              0.888              0.888              0.888              0.888 
    ##               NXT1              SCLT1             DPYSL4             GEMIN7 
    ##              0.888              0.888              0.889              0.889 
    ##             ARMCX5            METTL13               RPIA              IL4I1 
    ##              0.889              0.889              0.889              0.889 
    ##               MPDZ            C5orf30              BPIFC               BBS5 
    ##              0.890              0.890              0.890              0.890 
    ##              KHNYN             SAPCD2           C9orf147              ZNF74 
    ##              0.890              0.890              0.890              0.890 
    ##               NSL1             GTF3C5             SPICE1               EID1 
    ##              0.890              0.890              0.890              0.890 
    ##            PCDHGA1             RPL36A             CYBRD1               PYGB 
    ##              0.891              0.891              0.891              0.891 
    ##              ITIH3             ATP5G1            TMEM210              PADI2 
    ##              0.891              0.891              0.891              0.891 
    ##             ZNF697              SH2B2             PFKFB1              FANCE 
    ##              0.891              0.891              0.892              0.892 
    ##             MRPS22               MDM1                 KY              LRP10 
    ##              0.892              0.892              0.892              0.892 
    ##                CEL              WDR27              CES4A             DDX19B 
    ##              0.892              0.892              0.892              0.892 
    ##              CNTRL              ADCY4              ENOX2            FAM71E2 
    ##              0.892              0.892              0.892              0.892 
    ##             MRPL47               UBA7           SLC25A43             ZNF720 
    ##              0.892              0.892              0.892              0.892 
    ##              DPPA3            FAM131C            PPP2R3B              AGBL2 
    ##              0.892              0.892              0.892              0.893 
    ##            POGLUT1              TTC9B             GABRR2             SAMM50 
    ##              0.893              0.893              0.893              0.893 
    ##             MARCH2               ATF3              PANK3           KIAA1683 
    ##              0.893              0.893              0.893              0.893 
    ##              QRSL1               OPN3            DNAJC27               PAX4 
    ##              0.893              0.894              0.894              0.894 
    ##             DUSP12            TMEM55A              BEST4               PCCB 
    ##              0.894              0.894              0.894              0.894 
    ##               LGI2            ATP6V1D               LSG1             CCDC41 
    ##              0.894              0.894              0.894              0.894 
    ##              ZNF32              PGBD1              KIF27              MYH13 
    ##              0.895              0.895              0.895              0.895 
    ##               PRKX              ERCC4               FCN3                NMU 
    ##              0.895              0.896              0.896              0.896 
    ##              RXFP1              SOX15              PARVG              C1QL4 
    ##              0.896              0.896              0.896              0.896 
    ##              ABCC3             ZNF180               LSP1                 TG 
    ##              0.896              0.896              0.896              0.896 
    ##             TRIM56                EPO               CD46              MFGE8 
    ##              0.896              0.896              0.896              0.897 
    ##               TTC4              GDAP2              IQCB1           TMPRSS13 
    ##              0.897              0.897              0.897              0.897 
    ##               LAP3             TM9SF1              PLCD4               UPP2 
    ##              0.897              0.897              0.897              0.897 
    ##               H1F0              MCCC2            CCDC159             POLR2D 
    ##              0.897              0.897              0.897              0.897 
    ##             GPRC5D              ABCD4             ZBTB8B               MYL3 
    ##              0.898              0.898              0.898              0.898 
    ##            PCDHGA3              PROP1             TARSL2             NAP1L3 
    ##              0.898              0.898              0.898              0.898 
    ##              CWC27             CYSTM1            SLC17A2               ASB3 
    ##              0.899              0.899              0.899              0.899 
    ##               MBIP              SMYD3              PMPCB             GALNT6 
    ##              0.899              0.899              0.899              0.899 
    ##             DLGAP5               DOLK         GPR75-ASB3              GCNT7 
    ##              0.899              0.899              0.899              0.900 
    ##             SFT2D2               TBCB               RFC4              OTOGL 
    ##              0.900              0.900              0.900              0.900 
    ##               MLPH             TMEM8A           CYB561A3               PFN2 
    ##              0.900              0.900              0.901              0.901 
    ##              RPP40           CDC42EP4               ESM1            SLC36A1 
    ##              0.901              0.901              0.901              0.901 
    ##            ANKRD24              RNF13             KLHDC1            ST8SIA4 
    ##              0.901              0.901              0.901              0.901 
    ##           HSD17B12           CCDC163P                 F7              GSTA4 
    ##              0.901              0.902              0.902              0.902 
    ##              FUCA2              MCCC1              PRDX2              SGSM3 
    ##              0.902              0.902              0.902              0.902 
    ##               RIN3             POPDC3               VTA1      RP11-196G11.1 
    ##              0.902              0.902              0.902              0.902 
    ##            TRAPPC1              ACSS2              NAGLU            TMEM173 
    ##              0.902              0.903              0.903              0.903 
    ##              UTP23             ZBTB49             MRPL39                NLN 
    ##              0.903              0.903              0.903              0.903 
    ##              RBM24               CD80           SERPINB5             ZNF680 
    ##              0.903              0.903              0.903              0.903 
    ##             LRRC58             IL27RA           ARHGAP15               NOM1 
    ##              0.903              0.903              0.903              0.903 
    ##            N4BP2L2               MYOF            DNTTIP2               PIGQ 
    ##              0.904              0.904              0.904              0.904 
    ##             MCF2L2            GADD45B              SSBP1              CMYA5 
    ##              0.904              0.904              0.904              0.904 
    ##               SRP9                 C7               ATL3           C12orf75 
    ##              0.904              0.904              0.904              0.904 
    ##            TBC1D32             MCMDC2               ICOS              KCNK1 
    ##              0.904              0.904              0.905              0.905 
    ##              POTEE            NDUFAF1               KARS               UFC1 
    ##              0.905              0.905              0.905              0.905 
    ##                DYM              KCNS3               CKS2             PLXDC1 
    ##              0.905              0.905              0.905              0.905 
    ##              PCDP1               CLN3              CPXM1               EDAR 
    ##              0.906              0.906              0.906              0.906 
    ##              ZNF25            TIMM17B              LRIG2              SESN2 
    ##              0.906              0.906              0.906              0.906 
    ##              ARMC3           SLC25A46               PCCA            DPY19L2 
    ##              0.906              0.906              0.906              0.906 
    ##             SLC3A2            ZMYND12               GLMN                MIP 
    ##              0.906              0.906              0.906              0.907 
    ##           C12orf65             DONSON              PALD1              TACR1 
    ##              0.907              0.907              0.907              0.907 
    ##              IL1RN              AMPD3              EPHX4           KIAA0556 
    ##              0.907              0.907              0.907              0.907 
    ##         AC087645.1              TAGLN           TOR1AIP1            SLC44A4 
    ##              0.907              0.907              0.907              0.907 
    ##               TGS1               DUS2              ZCRB1              TIGD5 
    ##              0.908              0.908              0.908              0.908 
    ##              MEGF6              CRLF1               CHML           ANKRD30B 
    ##              0.908              0.908              0.908              0.908 
    ##             ANKRD6            OR52B1P               NME7           PCDHGA11 
    ##              0.908              0.909              0.909              0.909 
    ##             BCL2A1             UNC45B                GSR             HEATR6 
    ##              0.909              0.909              0.909              0.909 
    ##                AFP            ADAMTS8              ROBO4              RNF24 
    ##              0.909              0.909              0.909              0.909 
    ##              FIBIN              CISD2              ABCA4               PELO 
    ##              0.910              0.910              0.910              0.910 
    ##                FUZ              P2RY2             ZNF175            RPS6KL1 
    ##              0.910              0.910              0.910              0.910 
    ##               HSCB            IL12RB1               TNMD              DLEC1 
    ##              0.910              0.910              0.910              0.911 
    ##               ANO9              GCNT4             ZNF696              FYCO1 
    ##              0.911              0.911              0.911              0.911 
    ##             ZNF772              STMN3             ZCCHC7             LILRA6 
    ##              0.911              0.911              0.911              0.911 
    ##              ZNRD1              ALPK3              NAT10             LYSMD2 
    ##              0.911              0.911              0.911              0.911 
    ##              TGIF1             SCNN1A              SFXN4          TMPRSS11E 
    ##              0.911              0.912              0.912              0.912 
    ##              FCRL1               DERA              HLA-B            COL24A1 
    ##              0.912              0.912              0.912              0.912 
    ##              ILVBL              HMGN1            TMEM205              PQLC1 
    ##              0.912              0.912              0.912              0.913 
    ##             FAM20B            CCDC138             KLHL31             KRTDAP 
    ##              0.913              0.913              0.913              0.913 
    ##               COCH              SIMC1             ZNF311             ZNF329 
    ##              0.913              0.913              0.913              0.913 
    ##               SDE2            SLC26A6              NSUN4               CYCS 
    ##              0.913              0.914              0.914              0.914 
    ##                GH1             EPS8L3               IAH1           C18orf64 
    ##              0.914              0.914              0.914              0.915 
    ##               POC5              BRCA1              ZFP30              DARS2 
    ##              0.915              0.915              0.915              0.915 
    ##           ARHGEF38              IFIT5             ZNF322             ZNF813 
    ##              0.915              0.915              0.915              0.915 
    ##                DDC               ELK4              SAP18           ATP6AP1L 
    ##              0.915              0.915              0.915              0.915 
    ##              HSDL1           COLGALT2              STK16             TBC1D4 
    ##              0.915              0.915              0.915              0.915 
    ##            ZCCHC18              P2RX7             MICAL1               CHKB 
    ##              0.915              0.916              0.916              0.916 
    ##              ORAI2               TYMP             MRPL45             EEF1E1 
    ##              0.916              0.916              0.916              0.916 
    ##               HEXB             ACSBG2          TMPRSS11D            B3GALT1 
    ##              0.916              0.916              0.916              0.916 
    ##             WNT10B              PAQR3               FICD               ARL6 
    ##              0.916              0.916              0.916              0.917 
    ##            GLYATL1               CD69            KREMEN2               GDF1 
    ##              0.917              0.917              0.917              0.917 
    ##             ZNF484               SNX7             POLR1E             CC2D1B 
    ##              0.917              0.917              0.917              0.917 
    ##              DUSP1             ZCCHC9               ARG1                CA6 
    ##              0.917              0.917              0.917              0.917 
    ##            PHACTR1               IRF3               TPPP            PLA2G4D 
    ##              0.917              0.917              0.918              0.918 
    ##              VWA3A               GDNF            UBASH3A              GTF2B 
    ##              0.918              0.918              0.918              0.918 
    ##              FBXL4              NLRC4               RBAK              UEVLD 
    ##              0.918              0.918              0.918              0.918 
    ##             LRRC46             SLFN14             ICOSLG               PDK2 
    ##              0.918              0.918              0.919              0.919 
    ##             SEMA4F            ZSCAN29             ATP1B4              GCNT1 
    ##              0.919              0.919              0.919              0.919 
    ##              MEP1A            C2orf43              PRMT6            FASTKD5 
    ##              0.919              0.919              0.919              0.919 
    ##               INVS       RP11-295P9.3           C19orf38             FMR1NB 
    ##              0.919              0.919              0.919              0.919 
    ##              KCNA5      CTD-3074O7.11              IL12B              ADCK2 
    ##              0.919              0.919              0.919              0.920 
    ##             TXNDC2             TRIM43               RMI1              DUSP5 
    ##              0.920              0.920              0.920              0.920 
    ##              TUBD1               CD27             RBMXL1            TMEM50B 
    ##              0.920              0.920              0.920              0.920 
    ##             TMEM70             PNMAL1           C11orf63           CATSPERD 
    ##              0.920              0.920              0.920              0.920 
    ##             COL6A5             TYROBP             ZNF454           ARHGAP22 
    ##              0.920              0.920              0.921              0.921 
    ##             PRKAA2            TP53I13             TDRD10             TM6SF2 
    ##              0.921              0.921              0.921              0.921 
    ##             IGFBP6             LRRC32               FRZB              LIX1L 
    ##              0.921              0.921              0.921              0.922 
    ##            ALDH1L2            COL22A1               KNG1             ADAM32 
    ##              0.922              0.922              0.922              0.922 
    ##             ADAM33              MTHFR             TRMT2A            SLC29A4 
    ##              0.922              0.922              0.922              0.922 
    ##               ECM2              CEP41              BLVRA               CCZ1 
    ##              0.922              0.922              0.922              0.922 
    ##              OTUD3             PBXIP1              THBS4            TMEM203 
    ##              0.922              0.922              0.922              0.922 
    ##                NVL              RARS2             SLC2A2               RIN1 
    ##              0.923              0.923              0.923              0.923 
    ##               CD3E                TXN                MME             KBTBD7 
    ##              0.923              0.923              0.923              0.923 
    ##               KRT9             ZNF626             KLHL40              VPS51 
    ##              0.923              0.923              0.924              0.924 
    ##              IMPA2             CHI3L1              PSMG1             ANKMY2 
    ##              0.924              0.924              0.924              0.924 
    ##           C17orf59              DIS3L              CASP1               GLI1 
    ##              0.924              0.924              0.924              0.924 
    ##               PKLR              THOC3            TUBGCP6            COL19A1 
    ##              0.924              0.925              0.925              0.925 
    ##               MYOT              FLOT1             SH3D21              USP40 
    ##              0.925              0.925              0.925              0.925 
    ##               EXO1            FAM178B              HOOK2                YY2 
    ##              0.925              0.925              0.925              0.925 
    ##            CCDC180               MED9              AQPEP            C3orf67 
    ##              0.925              0.925              0.925              0.925 
    ##               KIF7              ACBD4              SYT17              SARS2 
    ##              0.926              0.926              0.926              0.926 
    ##              ADAM8               LY75         LY75-CD302              KRT80 
    ##              0.926              0.926              0.926              0.926 
    ##              RNFT1           KIAA1328              TUSC3             ZNF227 
    ##              0.926              0.926              0.926              0.926 
    ##               HRAS               GSG1              THBS3             FLVCR2 
    ##              0.926              0.927              0.927              0.927 
    ##              TELO2               DNA2               CDO1               PERP 
    ##              0.927              0.927              0.927              0.927 
    ##                VHL              HMOX1               OAS1             EGFLAM 
    ##              0.927              0.927              0.927              0.927 
    ##               SCD5              CALCR             TRMT13              PDSS1 
    ##              0.927              0.928              0.928              0.928 
    ##               IER3             PIH1D1               RELT              CNGA1 
    ##              0.928              0.928              0.928              0.928 
    ##              TMEM5             PRMT10              XIRP2            C9orf72 
    ##              0.928              0.928              0.928              0.928 
    ##           ADAMTS18              LMAN2             LMBRD1             CCDC30 
    ##              0.928              0.928              0.928              0.928 
    ##             PDZD11              NUP37              PDCD4              NUBPL 
    ##              0.928              0.928              0.928              0.929 
    ##              CENPW              ALPK1            ZNF286A            FTCDNL1 
    ##              0.929              0.929              0.929              0.929 
    ##                HRC             LY6G5C            CCDC42B            LAMTOR2 
    ##              0.929              0.929              0.929              0.929 
    ##             SEC61G             BCLAF1               TCN2             NT5DC4 
    ##              0.929              0.929              0.930              0.930 
    ##               MTG2             ANTXRL               LIPC              GPR52 
    ##              0.930              0.930              0.930              0.930 
    ##               XKR3               GGCX               PTMS              DNPEP 
    ##              0.930              0.930              0.930              0.930 
    ##               POLI             SPINK5               HFM1               TGM1 
    ##              0.930              0.930              0.930              0.930 
    ##             ZNF232              EGLN3              POC1A              GHRHR 
    ##              0.930              0.931              0.931              0.931 
    ##             GPR144               UMPS              GLP2R             ZNF557 
    ##              0.931              0.931              0.931              0.931 
    ##              CHEK1              PRR16              CSPP1           KIAA1430 
    ##              0.931              0.931              0.931              0.931 
    ##              KLF13              CNDP2             HEXIM2              CRLS1 
    ##              0.931              0.931              0.931              0.932 
    ##               ACO1             CYP7B1           ARHGAP10             FAM76A 
    ##              0.932              0.932              0.932              0.932 
    ##                LTA              MGST2            MARCH11             SCARA5 
    ##              0.932              0.932              0.932              0.932 
    ##              P2RX3            SLC10A4            SYNPO2L              SUSD1 
    ##              0.932              0.932              0.932              0.932 
    ##               TAB1                ECD              POC1B             HS6ST2 
    ##              0.932              0.933              0.933              0.933 
    ##              ITM2A             QRICH2               TLR2               NPNT 
    ##              0.933              0.933              0.933              0.933 
    ##             LRRC36              HLA-F           TMEM176A              WDR78 
    ##              0.933              0.933              0.933              0.933 
    ##              RSPH3               TPH1               ESPN               DOK7 
    ##              0.933              0.933              0.934              0.934 
    ##              GIPC3              ACRBP               RHOU              AMHR2 
    ##              0.934              0.934              0.934              0.934 
    ##              POMT1             DUSP22             TAGLN2               FMN1 
    ##              0.934              0.934              0.934              0.934 
    ##             HEATR2             KLHL33             MAD2L1              RPL29 
    ##              0.935              0.935              0.935              0.935 
    ##          AGAP2-AS1              YPEL4             RHBDF1              PEX19 
    ##              0.935              0.935              0.935              0.935 
    ##               ATG7                HGD              ANXA7               AOX1 
    ##              0.936              0.936              0.936              0.936 
    ##              OR2T4              APOC2                ASL           CD200R1L 
    ##              0.936              0.936              0.936              0.936 
    ##             NPC1L1             ANKEF1               NNAT             ZNF259 
    ##              0.936              0.936              0.937              0.937 
    ##            SEC14L2                AK1               GALC             ABCB10 
    ##              0.937              0.937              0.937              0.937 
    ##            ALDH5A1            RTN4IP1             ZNF605                LUM 
    ##              0.937              0.937              0.937              0.937 
    ##             LRRIQ1              PLOD1             ORMDL2             PITRM1 
    ##              0.938              0.938              0.938              0.938 
    ##             ATP1A4            PLEKHH2              ESPNL        CORO7-PAM16 
    ##              0.938              0.938              0.938              0.938 
    ##             SAMHD1  XXbac-BPG181M17.5              DUOX1              CAND2 
    ##              0.938              0.939              0.939              0.939 
    ##              LRGUK            NADSYN1            CCDC106         AC011530.4 
    ##              0.939              0.939              0.939              0.939 
    ##               ASB1             RSPH4A           ATP6V1C2            PCDHGA5 
    ##              0.939              0.939              0.939              0.939 
    ##             FBXL13             C9orf3               LIPN            CLEC17A 
    ##              0.939              0.939              0.940              0.940 
    ##            CEP57L1              IGFN1          KRTAP4-12               GBA2 
    ##              0.940              0.940              0.940              0.940 
    ##               CTSK             PRSS22               PHAX               TTF1 
    ##              0.940              0.940              0.940              0.940 
    ##              VWA3B            PCDHGB4       RP11-178C3.1             ABHD11 
    ##              0.940              0.941              0.941              0.941 
    ##              CNNM4             CYB5R4             ARFIP1               NAE1 
    ##              0.941              0.941              0.941              0.941 
    ##               H6PD             ABCC10              HAND2            GALNT12 
    ##              0.941              0.941              0.941              0.941 
    ##             RPS4Y1            SLC4A11              FCRL4            FAM207A 
    ##              0.942              0.942              0.942              0.942 
    ##              PRDX4             MCOLN3               RFC3            TMPRSS2 
    ##              0.942              0.942              0.942              0.942 
    ##              HBEGF              WDR87           TMEM151A           ATP6V0E1 
    ##              0.943              0.943              0.943              0.943 
    ##              PADI3              GMPPA              EXPH5             POLRMT 
    ##              0.943              0.943              0.943              0.943 
    ##           APOBEC3C              DUS3L               SHC4              CIDEC 
    ##              0.943              0.943              0.943              0.943 
    ##               MPP7              P4HA3             TOMM20              FUCA1 
    ##              0.943              0.943              0.943              0.944 
    ##             TMEM67             PLCXD3             RSPH6A              UPF3A 
    ##              0.944              0.944              0.944              0.944 
    ##             UBE2Q2            C1orf50              LIPT2              STMN1 
    ##              0.944              0.944              0.944              0.944 
    ##               OSMR            ARFGAP3           HLA-DQB1              AADAC 
    ##              0.944              0.944              0.944              0.944 
    ##              MYH7B           SLC22A15              CNPY2             CCDC60 
    ##              0.944              0.944              0.944              0.944 
    ##               BVES             KBTBD4              TPTE2             RCBTB2 
    ##              0.945              0.945              0.945              0.945 
    ##              SP110            ALDH6A1              DNAL1              PDCD6 
    ##              0.945              0.945              0.945              0.945 
    ##            COL10A1                WRN              STX17              TRMT1 
    ##              0.945              0.945              0.945              0.945 
    ##               SDF2              CPXM2               MCM8              MEF2B 
    ##              0.946              0.946              0.946              0.946 
    ##               TNK1             POU4F2              KRT71            TMEM119 
    ##              0.946              0.946              0.946              0.946 
    ##             ZRANB3            ATP13A4              AMELX              RAB17 
    ##              0.947              0.947              0.947              0.947 
    ##               DSC1              WDR46       NT5C1B-RDH14              SOAT1 
    ##              0.947              0.947              0.947              0.947 
    ##             ZNF510            NAPEPLD              NUDT9             ZNF729 
    ##              0.947              0.947              0.947              0.947 
    ##              MASTL               IMP4            FAM133A            PIP5KL1 
    ##              0.947              0.947              0.947              0.947 
    ##           SLC25A44                ADC              SFRP4              DPPA4 
    ##              0.947              0.947              0.948              0.948 
    ##             C6orf1               DPH2             TGOLN2               PIGT 
    ##              0.948              0.948              0.948              0.948 
    ##             CDKN2C              NUDT6             TXNRD3             LEFTY2 
    ##              0.948              0.948              0.948              0.948 
    ##                CGN              SPAG4            CCDC142               TCAP 
    ##              0.948              0.948              0.948              0.948 
    ##              OSER1           ATP6V1B1             SDF2L1              DDX10 
    ##              0.948              0.948              0.948              0.949 
    ##            C1orf27              NKAPL               ANO5             CDKN2A 
    ##              0.949              0.949              0.949              0.949 
    ##             KCTD13              KANK1              SYCP3             ZNF500 
    ##              0.949              0.949              0.949              0.949 
    ##            UQCRFS1             ZCCHC4             ZNF197            FAM110D 
    ##              0.949              0.949              0.949              0.950 
    ##              NGLY1              PPM1M             CHST10             INO80E 
    ##              0.950              0.950              0.950              0.950 
    ##              CDK14              ACADS            ARL6IP4               PKN3 
    ##              0.950              0.950              0.950              0.950 
    ##              SGPP2               NPM2             ZNF649               AHI1 
    ##              0.950              0.950              0.950              0.950 
    ##             TMEM8B             MROH2B              ACSS3             FKBP1B 
    ##              0.950              0.950              0.951              0.951 
    ##               GNMT              INSRR               GEN1               COQ5 
    ##              0.951              0.951              0.951              0.951 
    ##             CCDC43                RHD            RAPGEF3               AOAH 
    ##              0.951              0.951              0.951              0.951 
    ##               SMPX            TNFSF15               SOX1              ITFG2 
    ##              0.951              0.951              0.951              0.951 
    ##             AGTRAP             TRIM29             SCUBE2               OCLN 
    ##              0.951              0.951              0.951              0.951 
    ##               B9D1              DDAH1              SIVA1            TMEM241 
    ##              0.952              0.952              0.952              0.952 
    ##              A2ML1             HOXD13              RAB2B                ZP3 
    ##              0.952              0.952              0.952              0.952 
    ##             LRRC63             ZNF568             SMIM11               RSU1 
    ##              0.952              0.952              0.953              0.953 
    ##             CHCHD4              MYOM1              TFB2M             FER1L6 
    ##              0.953              0.953              0.953              0.953 
    ##             RNF207              PRR11             RASAL1             SHISA2 
    ##              0.953              0.953              0.953              0.953 
    ##              STAC2            IGHMBP2              TEX14             PGPEP1 
    ##              0.953              0.953              0.954              0.954 
    ##             DHRS13             METTL4               CDK9              UBTD2 
    ##              0.954              0.954              0.954              0.954 
    ##              VGLL4               KLK8             CHRNB3              PSMB9 
    ##              0.954              0.955              0.955              0.955 
    ##            PLEKHA8             PAPSS2            EFCAB14               GALT 
    ##              0.955              0.955              0.955              0.955 
    ##            CYP11B1             GNPDA2             SLC9A9                F2R 
    ##              0.955              0.955              0.955              0.955 
    ##             ARMCX1              NAA16            SLC26A3              ECHS1 
    ##              0.955              0.955              0.955              0.956 
    ##               NIFK           RPGRIP1L           C15orf40              PLCH2 
    ##              0.956              0.956              0.956              0.956 
    ##              TMCC3              RGS22            NGFRAP1           SLC25A17 
    ##              0.956              0.956              0.956              0.956 
    ##               SDC1               PAPL              TTC24              ZNRF2 
    ##              0.956              0.957              0.957              0.957 
    ##              SIRT2              ZUFSP              CASC1           C21orf91 
    ##              0.957              0.957              0.957              0.957 
    ##              IBA57             BPIFB1              STX18           CDC42EP3 
    ##              0.957              0.957              0.957              0.957 
    ##              HTR3B               MOGS             PFKFB4               EMCN 
    ##              0.957              0.957              0.957              0.958 
    ##               PGM2             ZNF540             ZNF546            WBSCR16 
    ##              0.958              0.958              0.958              0.958 
    ##              FBLN7             NLRP13         AC018470.1                CTH 
    ##              0.958              0.958              0.958              0.959 
    ##             SLC7A7             DBNDD2              GIPC1             DGCR14 
    ##              0.959              0.959              0.959              0.959 
    ##           C19orf67             SH2D1A                REN            C5orf54 
    ##              0.959              0.959              0.959              0.960 
    ##              CNGA4             ZNF431            PPP2R3C               ANO2 
    ##              0.960              0.960              0.960              0.961 
    ##               ASZ1             PARP15             VPS33B             CTDSPL 
    ##              0.961              0.961              0.961              0.961 
    ##              DDX11             LRRC24            RASL10A            RNASET2 
    ##              0.961              0.961              0.961              0.961 
    ##       TMEM56-RWDD3           KIAA1377               LINS              MYO1H 
    ##              0.961              0.962              0.962              0.962 
    ##           C11orf80            ZSCAN5A             DIRAS1              SETD9 
    ##              0.962              0.962              0.962              0.962 
    ##             SEC16B            CCDC114              LAMP3                OAT 
    ##              0.962              0.962              0.962              0.962 
    ##              BBS12            PPP1R3C            ALOX5AP             NT5C1A 
    ##              0.963              0.963              0.963              0.963 
    ##               SYS1              HAUS6            BLOC1S1             NPIPB4 
    ##              0.963              0.963              0.963              0.963 
    ##               YBX3              FNDC7            SLC35B4             ZNF33A 
    ##              0.963              0.963              0.963              0.963 
    ##               GZMK              CDH17              LACE1               EDN3 
    ##              0.963              0.963              0.964              0.964 
    ##              RPL37                BSG              CERS5                CR1 
    ##              0.964              0.964              0.964              0.964 
    ##                AUH             MYBPC2              PCSK4            SPATA24 
    ##              0.964              0.964              0.964              0.964 
    ##           SLC35E2B             ZNF75D         ST6GALNAC3              FOXQ1 
    ##              0.964              0.964              0.964              0.965 
    ##             LRRCC1             PPAP2C              PPEF2              SFTA3 
    ##              0.965              0.965              0.965              0.965 
    ##             ZCWPW1              DUS1L               COG5              KRT35 
    ##              0.965              0.965              0.965              0.965 
    ##            SLC27A2               PYGL            COLEC10              CAPN3 
    ##              0.965              0.965              0.965              0.965 
    ##            SLC35F2              MANBA               CD99               PASK 
    ##              0.965              0.965              0.965              0.965 
    ##             PITHD1              SYCE1              CDHR4             CC2D2B 
    ##              0.965              0.965              0.966              0.966 
    ##               CPN2               MTBP            PYROXD1               SC5D 
    ##              0.966              0.966              0.966              0.966 
    ##              RFTN1            ASPSCR1             TBC1D5              HMOX2 
    ##              0.966              0.966              0.966              0.966 
    ##              JMJD7               PPCS            STARD10            TRMT61A 
    ##              0.966              0.966              0.966              0.966 
    ##               PDP2            SIGLEC9            TMPRSS9             SLC6A5 
    ##              0.966              0.967              0.967              0.967 
    ##              GLB1L              AIMP2             STEAP2             CCDC18 
    ##              0.967              0.967              0.967              0.967 
    ##            CCDC158              PANK2            SPARCL1             MAN2B1 
    ##              0.967              0.967              0.967              0.967 
    ##            FAM212B              H1FOO             GUCY2F              THADA 
    ##              0.968              0.968              0.968              0.968 
    ##              CXCR3              MMRN2              SEPT1             ABCA10 
    ##              0.968              0.968              0.968              0.968 
    ##              MTRF1            TNFAIP8              KIF12             ZNF263 
    ##              0.969              0.969              0.969              0.969 
    ##               NPM3             MAD1L1              OOSP1               AOC1 
    ##              0.969              0.969              0.969              0.969 
    ##             MRGPRF            ZFAND2B             TM4SF5             RAD54L 
    ##              0.969              0.969              0.970              0.970 
    ##            C3orf17               DLX3           C17orf49               CD96 
    ##              0.970              0.970              0.970              0.970 
    ##             ZNF607             LRP2BP           HLA-DQA1              ELMO3 
    ##              0.970              0.970              0.970              0.970 
    ##              CORO6              KIF19             PI4K2B                F10 
    ##              0.970              0.970              0.970              0.970 
    ##             ZNF569               MUC6               ASB6               CA5A 
    ##              0.970              0.971              0.971              0.971 
    ##              CABP2              SCML4            SLCO1C1              SMYD1 
    ##              0.971              0.971              0.972              0.972 
    ##             ITGBL1              CENPJ               GBE1              NKPD1 
    ##              0.972              0.972              0.972              0.972 
    ##               MRS2              PLVAP                XPA               PFKM 
    ##              0.972              0.972              0.972              0.972 
    ##              RABL3             EHHADH               BIN2             ZNF599 
    ##              0.972              0.972              0.973              0.973 
    ##             CCDC14         NEDD8-MDP1              CHST7              CNNM3 
    ##              0.973              0.973              0.973              0.973 
    ##             CARD16             TTC39B               UACA             CD300E 
    ##              0.973              0.973              0.973              0.973 
    ##               CCR4             SRD5A3               HEXA              APBB3 
    ##              0.973              0.973              0.973              0.974 
    ##               GNB5              CASP9              GSKIP               PSTK 
    ##              0.974              0.974              0.974              0.974 
    ##             GPR155               USP6           KIAA0895              TRPV5 
    ##              0.974              0.974              0.974              0.974 
    ##              WDR53               CLGN                PZP            SLC52A2 
    ##              0.974              0.975              0.975              0.975 
    ##             ZNF394               TTC9            SLC35C1              NPHP1 
    ##              0.975              0.975              0.975              0.975 
    ##               TMC1             NUFIP1           PAFAH1B3               ERI2 
    ##              0.975              0.975              0.975              0.975 
    ##              UROC1              TNNT3             SNAPC4               SNCB 
    ##              0.975              0.975              0.976              0.976 
    ##             ZNF619             GUCA2B           ARHGEF10              TSTD1 
    ##              0.976              0.976              0.976              0.976 
    ##             EIF2B1             EXOSC4          HIST2H2AB             ZYG11A 
    ##              0.976              0.976              0.976              0.976 
    ##               GDF2               GPR6               CCR2             CLDN15 
    ##              0.976              0.976              0.976              0.977 
    ##             DUOXA1           KIAA1009            HLA-DMB               RFC2 
    ##              0.977              0.977              0.977              0.977 
    ##              CAGE1             HRSP12             RECQL5           TIMELESS 
    ##              0.977              0.977              0.977              0.977 
    ##              GPAA1              SRBD1              ABCB6             OR52E8 
    ##              0.977              0.977              0.977              0.977 
    ##         AP000295.9              CHRND             SEC23B               SOD1 
    ##              0.978              0.978              0.978              0.978 
    ##             KLHL26              RFXAP              MASP1             SPTLC1 
    ##              0.978              0.978              0.978              0.978 
    ##           CATSPER3            SLC22A4              DPY30            SLC6A13 
    ##              0.978              0.978              0.978              0.978 
    ##              EMID1             GRIN3B              RSAD1              MROH7 
    ##              0.978              0.978              0.979              0.979 
    ##              BEND6             ZNF783              FAM3D            SHARPIN 
    ##              0.979              0.979              0.979              0.979 
    ##             AGPAT5           TRAF3IP2              PATL2               IQCH 
    ##              0.979              0.979              0.979              0.979 
    ##              TBX19            TMEM174              SFTPD              CELA1 
    ##              0.979              0.979              0.979              0.980 
    ##             IGFBP4             STK32A              CCKBR              MROH9 
    ##              0.980              0.980              0.980              0.980 
    ##             CHRNA9            SLC14A2               TCHH              CBLN1 
    ##              0.980              0.980              0.980              0.980 
    ##             ANXA11             PCDH12              THAP9               BBS1 
    ##              0.980              0.980              0.980              0.981 
    ##              GNAT2             MAP7D3            SIGLEC5             SMTNL2 
    ##              0.981              0.981              0.981              0.981 
    ##             TOM1L1             GALNT5            EXOC3L4              CRYGS 
    ##              0.981              0.981              0.981              0.981 
    ##                GAA             KNSTRN              RRP7A             CLTCL1 
    ##              0.981              0.981              0.981              0.982 
    ##                PTN           C15orf27              CENPU          GOLGA6L18 
    ##              0.982              0.982              0.982              0.982 
    ##              ABHD5             PTGER2              DNAI2             NT5C1B 
    ##              0.982              0.982              0.982              0.982 
    ##              MTMR7           TMEM255A             ZNF679              SVOPL 
    ##              0.982              0.982              0.982              0.982 
    ##              CCDC7              NPHP4              OR9I1            ARHGEF5 
    ##              0.982              0.983              0.983              0.983 
    ##             SEPT14            LRRFIP2               MTRR            SIGLEC5 
    ##              0.983              0.983              0.983              0.983 
    ##             SPINT2            RASL11B              DDX47              KLK11 
    ##              0.983              0.983              0.983              0.983 
    ##             NDUFB8            CYP27B1              VSIG1            HLA-DMA 
    ##              0.983              0.983              0.983              0.983 
    ##               MDS2             PTPN18              BPNT1             PNLDC1 
    ##              0.984              0.984              0.984              0.984 
    ##            C6orf47            DEPDC1B              ABHD6               GLI4 
    ##              0.984              0.984              0.984              0.984 
    ##              MUC19               ASS1              CEBPD             DTNBP1 
    ##              0.984              0.984              0.984              0.984 
    ##                PRX               TEPP             CEP290            CCDC183 
    ##              0.984              0.985              0.985              0.985 
    ##               EOGT               NWD1              TMED3               FMOD 
    ##              0.985              0.985              0.985              0.985 
    ##              ECSIT              RSPH1            METTL7B             TSPAN6 
    ##              0.985              0.985              0.985              0.985 
    ##            HORMAD2             PCDHA7              MUS81              CRYAA 
    ##              0.986              0.986              0.986              0.986 
    ##             SLC2A6             GIMAP8                OGN             ZNF823 
    ##              0.986              0.986              0.986              0.986 
    ##      CTD-2140B24.4            CD300LG             KDELC1             DNAAF3 
    ##              0.986              0.987              0.987              0.987 
    ##              PRMT7              GORAB              CMTM5            ANKRD63 
    ##              0.987              0.987              0.987              0.987 
    ##              AIFM3            TSPAN32              CCNG1             ZNF597 
    ##              0.988              0.988              0.988              0.988 
    ##               PIGG              ACTR8             NKAIN4               AAAS 
    ##              0.988              0.988              0.988              0.988 
    ##          FAM120AOS            SEPSECS            TSPAN33            CEACAM5 
    ##              0.989              0.989              0.989              0.989 
    ##            SLC31A1               MSH4            PLEKHJ1             IL1RL2 
    ##              0.989              0.989              0.989              0.989 
    ##              SAR1B              PNLIP             LGALS9             COMMD6 
    ##              0.989              0.989              0.989              0.990 
    ##             TXNRD2              CHRNE              TCTN2              PARP2 
    ##              0.990              0.990              0.990              0.990 
    ##            SLC22A6              ITGB2             PRSS50              CEP95 
    ##              0.990              0.990              0.990              0.990 
    ##           ARHGEF26             CCDC67             PKD2L2              PKDCC 
    ##              0.990              0.990              0.990              0.990 
    ##             TRIM22           CATSPER1              UBE3D              LACTB 
    ##              0.990              0.990              0.990              0.990 
    ##               AQP4              RAD18           C6orf163               LSM4 
    ##              0.991              0.991              0.991              0.991 
    ##            PPP1R3B               CIR1             MAPK12               GBP2 
    ##              0.991              0.991              0.991              0.991 
    ##             AXDND1            C1QTNF4               CD53              IFT80 
    ##              0.991              0.991              0.991              0.991 
    ##               SCO1                RGN             PET112              TPPP3 
    ##              0.991              0.991              0.992              0.992 
    ##              CHTF8              CD180             ZNF212              DCAKD 
    ##              0.992              0.992              0.992              0.992 
    ##             DCDC2C             NUTM2D               TLR4             PM20D2 
    ##              0.992              0.992              0.992              0.992 
    ##              MYO3B             ZNF264             ZNF268                 C6 
    ##              0.992              0.992              0.992              0.993 
    ##               MUC4               IQUB               BFAR              WDR34 
    ##              0.993              0.993              0.993              0.993 
    ##               CRAT               FBP2            PLEKHB1            SLC15A4 
    ##              0.993              0.993              0.993              0.993 
    ##           KIAA0319            PPP1R1A             LPCAT2              PREPL 
    ##              0.993              0.994              0.994              0.994 
    ##             CPSF3L             CCDC73            GTF2H2C               LPXN 
    ##              0.994              0.994              0.994              0.994 
    ##               NADK             MMADHC             SNAP47             BRI3BP 
    ##              0.994              0.994              0.994              0.994 
    ##              ISCA1              WDFY2              CLCN1             SLC7A8 
    ##              0.994              0.994              0.994              0.994 
    ##             SEC11A               USE1              MARC2              BRAT1 
    ##              0.994              0.994              0.994              0.995 
    ##              DAGLB             MIS18A                NRM              CEP70 
    ##              0.995              0.995              0.995              0.995 
    ##            FAM195A               SBSN               WIBG              CEBPE 
    ##              0.995              0.995              0.995              0.995 
    ##             FBXO39             ZNF879            ALDH3B1           TMEM194B 
    ##              0.995              0.995              0.995              0.995 
    ##            OSGEPL1               CST6               GHDC             CARD14 
    ##              0.995              0.995              0.995              0.995 
    ##             ASPHD1             CHRNB1               WTIP              GCSAM 
    ##              0.996              0.996              0.996              0.996 
    ##              DCTN5              CNIH4             CRYBB1               CD1D 
    ##              0.996              0.996              0.996              0.996 
    ##            OLFML2A              SYDE2               KLK5              SAAL1 
    ##              0.996              0.997              0.997              0.997 
    ##              TMCO6              TRPT1             AKR1B1              ASMTL 
    ##              0.997              0.997              0.997              0.997 
    ##              FBXO2               CCL2                RFK              MMEL1 
    ##              0.997              0.997              0.997              0.997 
    ##              RAB8B            CCDC141            R3HCC1L             LRRC34 
    ##              0.997              0.997              0.997              0.998 
    ##               HADH              CPED1              MRPS9             SYCP2L 
    ##              0.998              0.998              0.998              0.998 
    ##             TMEM68             STK32B             SEPHS2               LIPA 
    ##              0.998              0.998              0.998              0.998 
    ##           KIAA1614             RNF217             CCDC19               MED6 
    ##              0.998              0.998              0.998              0.998 
    ##              RPLP1             TTC30A              DDX51               CHAD 
    ##              0.998              0.999              0.999              0.999 
    ##          TMPRSS11F               EXD3              TRPM5               GJA8 
    ##              0.999              0.999              0.999              0.999 
    ##             CAPN10           B3GALNT1               USB1              TAF1A 
    ##              0.999              0.999              0.999              0.999 
    ##               NRAP              RAB34            FAM71F2             DNAAF1 
    ##              0.999              1.000              1.000              1.000 
    ##              P2RY4            DNAJC18            SLC10A7              OR2A7 
    ##              1.000              1.000              1.000              1.000 
    ##            SLCO4A1              GAP43             TMEM53             ZNF75A 
    ##              1.000              1.000              1.000              1.000 
    ##               ARL2               NAT6               ZNF3              SPAG7 
    ##              1.000              1.000              1.000              1.000 
    ##              SPSB2              ACOT8             ZNF615               BPGM 
    ##              1.001              1.001              1.001              1.001 
    ##               LY6H             SH3BGR           C9orf114              SRP19 
    ##              1.001              1.001              1.001              1.001 
    ##             TYSND1               GSAP            PLEKHG6            ZNF804B 
    ##              1.001              1.001              1.001              1.001 
    ##              ACER1             ASAH2C            SLCO6A1             SIGIRR 
    ##              1.002              1.002              1.002              1.002 
    ##              FBXW9             CORO2A             CRYBA1              CNDP1 
    ##              1.002              1.002              1.002              1.002 
    ##             CXCL10              APOOL               TC2N              RGS20 
    ##              1.002              1.002              1.002              1.002 
    ##             PRSS23           C16orf11             UBALD2               CLMP 
    ##              1.003              1.003              1.003              1.003 
    ##            C2orf16              CDRT1             SMIM23            CCDC147 
    ##              1.003              1.003              1.003              1.003 
    ##             NKX6-1             CAB39L               UTP3             ZNF235 
    ##              1.003              1.003              1.003              1.004 
    ##              ELTD1           ATP6V0A4              CBLN4               CD59 
    ##              1.004              1.004              1.004              1.004 
    ##              RXFP2              PQLC2            FAM118A               LIPK 
    ##              1.004              1.004              1.004              1.004 
    ##            C9orf78           HIST1H3J             PIK3R3              KLKB1 
    ##              1.004              1.004              1.004              1.004 
    ##             DDX60L            C4orf21              HHLA2               PPA1 
    ##              1.005              1.005              1.005              1.005 
    ##               YDJC               TMC4            COL21A1             SLC9B1 
    ##              1.005              1.005              1.005              1.005 
    ##               DPYD              P4HTM            KBTBD11             ZNF416 
    ##              1.005              1.005              1.005              1.005 
    ##              TTC29               ZXDC                DBH                BTC 
    ##              1.005              1.005              1.005              1.005 
    ##               TMPO               CDC6             RSL1D1           C16orf54 
    ##              1.005              1.005              1.005              1.006 
    ##              ACTR5           C6orf211            PPP1R3A          C10orf129 
    ##              1.006              1.006              1.006              1.006 
    ##              PALB2                BMF             NUDCD1             FAM71C 
    ##              1.006              1.006              1.006              1.006 
    ##              CCBE1               AIM1               DIS3            TMEM209 
    ##              1.006              1.006              1.006              1.006 
    ##            PLA2G4C            FAM175A              ADCK3          C20orf196 
    ##              1.006              1.007              1.007              1.007 
    ##              BLZF1              FKBP3              TTC22               ESR2 
    ##              1.007              1.007              1.007              1.007 
    ##             ZNF614              ZPLD1              HACL1            NKIRAS2 
    ##              1.007              1.007              1.007              1.007 
    ##             GEMIN4           C16orf80              MEIS3              BSCL2 
    ##              1.007              1.008              1.008              1.008 
    ##            TPD52L2               FGL2             RASSF2              LAMB4 
    ##              1.008              1.008              1.008              1.008 
    ##               ECE2             KCTD12              MTFMT              KCNN4 
    ##              1.008              1.008              1.008              1.008 
    ##             MRPL33               EMC1              PTF1A               SCEL 
    ##              1.009              1.009              1.009              1.009 
    ##               GSG2             HARBI1              CREB3             CX3CR1 
    ##              1.009              1.010              1.010              1.010 
    ##             NHLRC3               PWP1              RPS20               F11R 
    ##              1.010              1.010              1.010              1.010 
    ##            LEPREL2              SYTL2            SLC45A2              NEIL1 
    ##              1.010              1.010              1.010              1.010 
    ##               MSX1                NBN               SELL            POMGNT2 
    ##              1.010              1.010              1.010              1.010 
    ##              CRIP3              FXYD5             AHNAK2               G6PC 
    ##              1.010              1.010              1.010              1.010 
    ##               COQ9      RP11-514O12.4             MAN2B2            CYP51A1 
    ##              1.011              1.011              1.011              1.011 
    ##     RPL36A-HNRNPH2               CTSH           CACNA2D4        CTC-534A2.2 
    ##              1.011              1.011              1.011              1.011 
    ##     ARHGAP19-SLIT1             ATP8B3            SNRNP27              UNC50 
    ##              1.011              1.011              1.012              1.012 
    ##               GHSR              POTEF               STAC             LOXHD1 
    ##              1.012              1.012              1.012              1.012 
    ##             POLR1D               DNTT              LIMD2              TRIM5 
    ##              1.013              1.013              1.013              1.013 
    ##             COL6A6              GSTM1              GALNS             SLC7A2 
    ##              1.013              1.013              1.013              1.014 
    ##             UBQLN3             SCN10A            TMPRSS4             SLC9A8 
    ##              1.014              1.014              1.014              1.014 
    ##               SFI1            FAM161A       RP11-758M4.1              HMCN2 
    ##              1.014              1.014              1.015              1.015 
    ##           EPB41L4A              PEX13               PMS1            SOSTDC1 
    ##              1.015              1.015              1.015              1.015 
    ##             EFCAB3            ZNF286B            TRMT61B             RIMKLA 
    ##              1.015              1.015              1.015              1.015 
    ##              TLR10             ZNF302             ZNF441               ABT1 
    ##              1.015              1.015              1.015              1.015 
    ##            C2orf82              MED27               TRDN              ANXA4 
    ##              1.016              1.016              1.016              1.016 
    ##           C1orf198               INHA               EDN2             BHLHB9 
    ##              1.016              1.016              1.016              1.016 
    ##             PARD6A               DBR1            C7orf31              KRT75 
    ##              1.016              1.016              1.016              1.016 
    ##     RP11-155D18.14              HYAL1           KIAA0101           APOBEC3H 
    ##              1.016              1.017              1.017              1.017 
    ##             ENTPD4              TM2D2             TM7SF2              NLRP7 
    ##              1.017              1.017              1.017              1.017 
    ##                BOK              CORO7              HLA-C                PGP 
    ##              1.017              1.017              1.017              1.017 
    ##               TCHP              MICU3              SEPN1            C4orf47 
    ##              1.017              1.017              1.017              1.018 
    ##           NFATC2IP              ZMAT3               RALB              CNGB3 
    ##              1.018              1.018              1.018              1.018 
    ##                 HP              CEP78             TRMT2B               SSH3 
    ##              1.018              1.018              1.018              1.018 
    ##       RP11-796G6.2              TEX40              SYT12             ZNF461 
    ##              1.018              1.018              1.018              1.018 
    ##               WBP4              PDE6D           ARHGEF16             ZNF771 
    ##              1.018              1.018              1.018              1.018 
    ##               CSAD              CAMLG              CAPN8           APOBEC3F 
    ##              1.019              1.019              1.019              1.019 
    ##             FUNDC2               HPS3            WBSCR27             CCDC57 
    ##              1.019              1.019              1.019              1.019 
    ##             TM9SF1               RHOJ            FAM110A            SLC35B1 
    ##              1.019              1.019              1.019              1.019 
    ##              HOXB4               PODN            C3orf38              CMTM8 
    ##              1.019              1.019              1.020              1.020 
    ##              AIF1L            ZCCHC10             VKORC1               COG6 
    ##              1.020              1.020              1.020              1.020 
    ##                JUN              MATN4            GPALPP1              NPY5R 
    ##              1.020              1.020              1.020              1.020 
    ##             SNAPC3               TREH            MAP3K19              CALB2 
    ##              1.020              1.020              1.020              1.020 
    ##            CREB3L4            PCDHAC1               COA7              HAUS8 
    ##              1.020              1.020              1.020              1.021 
    ##               CDK7            DYNLRB1             ZWILCH              EARS2 
    ##              1.021              1.021              1.021              1.021 
    ##            TSPAN17                FGA              RGAG4              STAP1 
    ##              1.021              1.021              1.021              1.021 
    ##               EGR4              NOMO1              BCAT1              SYCE3 
    ##              1.021              1.022              1.022              1.022 
    ##           SLC22A11       RP11-47I22.4              SMPD4             GEMIN2 
    ##              1.022              1.022              1.022              1.022 
    ##            GYLTL1B           MARVELD1              HTR5A             CCDC34 
    ##              1.022              1.022              1.023              1.023 
    ##              GABRE               STX7             METTL1             ERMARD 
    ##              1.023              1.023              1.023              1.023 
    ##            MRPS18B               NBL1            SLC30A7              MYL10 
    ##              1.023              1.023              1.023              1.024 
    ##              FANCI              SPHK2              TREM1              ERAL1 
    ##              1.024              1.024              1.024              1.024 
    ##               CD34             CCDC86             THNSL2            VSIG10L 
    ##              1.024              1.024              1.024              1.024 
    ##              CEBPB             ATP10B               ELP3             HDGFL1 
    ##              1.025              1.025              1.025              1.025 
    ##              ENPP3              PLCZ1              TAF1B               ARTN 
    ##              1.025              1.025              1.025              1.025 
    ##               TBCA           PPAPDC1A              SEP15               PIGN 
    ##              1.025              1.025              1.025              1.025 
    ##            TRMT10B              DCHS2             ALKBH1            PPP1R11 
    ##              1.025              1.025              1.026              1.026 
    ##               CDR2              CRIP2              PSMA8               GFAP 
    ##              1.026              1.026              1.026              1.026 
    ##               VAX2             LILRA3               EML2             SH2D3A 
    ##              1.026              1.026              1.026              1.026 
    ##             SPPL2C               CRB3              ABCC6             UBIAD1 
    ##              1.026              1.026              1.026              1.026 
    ##              NAA20            SLC50A1              GNAT3                 C2 
    ##              1.026              1.027              1.027              1.027 
    ##             DNAJB2               PON2            ANKDD1B            CCDC125 
    ##              1.027              1.027              1.027              1.027 
    ##             ZNF485             EPS8L1               AGXT               LMO2 
    ##              1.027              1.027              1.028              1.028 
    ##             LIN28A           C17orf80               ELP4           C12orf49 
    ##              1.028              1.028              1.028              1.028 
    ##     RP11-216L13.17           SECISBP2             ZNF420                PPL 
    ##              1.028              1.028              1.028              1.028 
    ##              SEMG1              CEND1               ERN2             CNRIP1 
    ##              1.028              1.029              1.029              1.029 
    ##              TTLL8               MYH8            TMEM87B            ANKRD45 
    ##              1.029              1.029              1.029              1.030 
    ##              TWSG1              DZIP1                CD5            TMEM38A 
    ##              1.030              1.030              1.030              1.030 
    ##             ZCWPW2              CLCA4             NDUFS3              ANXA1 
    ##              1.030              1.030              1.030              1.030 
    ##             CACNG5              FAIM3             POLR3C              ADAM7 
    ##              1.030              1.030              1.030              1.031 
    ##               GBP5              COASY              G6PC3             STEAP4 
    ##              1.031              1.031              1.031              1.031 
    ##              MED16              HMGN3              SDCBP            TMEM86A 
    ##              1.031              1.031              1.031              1.031 
    ##              FKBP6               HAGH              EXTL2              IFI16 
    ##              1.031              1.031              1.031              1.031 
    ##               PBX4               MST1             SECTM1             ZNF564 
    ##              1.032              1.032              1.032              1.032 
    ##             BCKDHA              NOXA1               ALX3              MKNK1 
    ##              1.032              1.032              1.032              1.032 
    ##               AOC3           CEACAM16              SUGCT              RAB44 
    ##              1.032              1.032              1.032              1.032 
    ##              THG1L          HIST1H2AC              ALG12              PLBD2 
    ##              1.032              1.033              1.033              1.033 
    ##               RTP4              GCNT6              HEMK1            ADPRHL2 
    ##              1.033              1.033              1.033              1.033 
    ##           KIAA1755              S1PR5             GPR156       TRIM6-TRIM34 
    ##              1.033              1.034              1.034              1.034 
    ##              PRIM2            SLC22A3            IFITM10             UBAP1L 
    ##              1.034              1.034              1.034              1.034 
    ##             ADCY10             CYP4V2               APOE               TMC3 
    ##              1.034              1.034              1.034              1.034 
    ##              GOLM1             LRRC14               NOP9              AREGB 
    ##              1.034              1.034              1.034              1.034 
    ##             IL31RA               KLK6             ABCA13              AGBL3 
    ##              1.034              1.034              1.035              1.035 
    ##              MIPEP            GOLGA7B             NKX2-6               PFKL 
    ##              1.035              1.035              1.035              1.035 
    ##             POLR3H             BAIAP3                DCK           SLC25A40 
    ##              1.035              1.035              1.035              1.035 
    ##               CTSC               UFM1              TTC33              ATHL1 
    ##              1.035              1.035              1.035              1.036 
    ##           LRRC37A2             ZNF562              DIMT1              RPL22 
    ##              1.036              1.036              1.036              1.036 
    ##              RBM45             ELOVL3              TSSK3               ARSJ 
    ##              1.036              1.036              1.036              1.037 
    ##           SLC25A53             SEMA3G              PSMD5             TRIM34 
    ##              1.037              1.037              1.037              1.037 
    ##             ARL13B              AAGAB             CD99L2             CD2BP2 
    ##              1.037              1.037              1.037              1.037 
    ##              BFSP2               VRK2             TMEM62            COL28A1 
    ##              1.038              1.038              1.038              1.038 
    ##             TMCO5A            C7orf61            SLC5A10               GFM2 
    ##              1.038              1.038              1.039              1.039 
    ##             P2RY10               ISLR               IL10              MRPL1 
    ##              1.039              1.039              1.039              1.039 
    ##              UQCC1              UQCC2                ZYX             ALOX12 
    ##              1.039              1.039              1.039              1.039 
    ##             SGK494               TDO2              TRPM8               PLTP 
    ##              1.039              1.039              1.039              1.040 
    ##           FAM114A2                 SI              PUS10              CDCA2 
    ##              1.040              1.040              1.040              1.040 
    ##               ABI3              RTKN2              HADHB              NIPA1 
    ##              1.040              1.040              1.040              1.040 
    ##              TCTN3              MYO3A           DCAF12L2               IL7R 
    ##              1.040              1.040              1.040              1.040 
    ##              ZNF23              TUBB6              DDX20               FIBP 
    ##              1.041              1.041              1.041              1.041 
    ##               VNN1               DRC1             LHFPL2             FBXW12 
    ##              1.041              1.041              1.041              1.041 
    ##             ZNF567               DTD1           C1orf127               MYH4 
    ##              1.041              1.041              1.041              1.041 
    ##              NDNL2               LYAR              ZMYM5               CD37 
    ##              1.041              1.041              1.041              1.041 
    ##              INADL               MKS1             LRRC48             ZNF713 
    ##              1.042              1.042              1.042              1.042 
    ##            CCDC181            POMGNT1              FANCC              BIRC5 
    ##              1.042              1.042              1.043              1.043 
    ##           C16orf86               MTX3               TOM1             GTF2H3 
    ##              1.043              1.043              1.043              1.043 
    ##           INS-IGF2             LGALSL               CCR6             TRIM31 
    ##              1.043              1.043              1.043              1.043 
    ##               FGF2            KLHDC8B              TCP11         AC003002.6 
    ##              1.043              1.043              1.043              1.043 
    ##               OSTC                TNN               VSX2              SDAD1 
    ##              1.043              1.043              1.043              1.044 
    ##               SHBG              ETFDH            SLC17A5               C1QB 
    ##              1.044              1.044              1.044              1.044 
    ##               CCL3                AGA              BTNL8            LRRC37A 
    ##              1.044              1.044              1.044              1.044 
    ##             NDUFB9             ADHFE1             SCARA3              EPHA1 
    ##              1.044              1.045              1.045              1.045 
    ##              S1PR2              SIRT3            CLEC18B         AC174470.1 
    ##              1.045              1.045              1.045              1.045 
    ##              OXA1L               TTI2              SPDYC               LLPH 
    ##              1.045              1.045              1.045              1.045 
    ##              ELAC2             FAM65C            POM121C             RHBDL3 
    ##              1.046              1.046              1.046              1.046 
    ##               SGCG             TM6SF1              DESI1              TIPIN 
    ##              1.046              1.047              1.047              1.047 
    ##            PABPC4L              SOCS4             KRT222            ARL6IP6 
    ##              1.047              1.047              1.047              1.047 
    ##             CRELD1              MYOD1              ARL5B           SLC25A16 
    ##              1.048              1.048              1.048              1.048 
    ##              TTC12              MON1A               WNK4       ABHD14A-ACY1 
    ##              1.048              1.049              1.049              1.049 
    ##            ZDHHC12            ZDHHC22            TMEM101               DKK3 
    ##              1.049              1.049              1.049              1.049 
    ##          KRTAP17-1              RBM23              GNB1L           FAM149B1 
    ##              1.049              1.049              1.049              1.049 
    ##              PROCR        AP000275.65             SEC31B              TIGIT 
    ##              1.050              1.050              1.050              1.050 
    ##             ACADVL            METTL25             PNPLA5              PPM1J 
    ##              1.050              1.050              1.050              1.050 
    ##               NOX3            TMPRSS5              VIPR2              LCA5L 
    ##              1.050              1.050              1.050              1.050 
    ##              GPR27              PROM1             CRELD2              ACAT1 
    ##              1.051              1.051              1.051              1.051 
    ##                ZAN           C9orf152              EPDR1            EMILIN2 
    ##              1.051              1.051              1.051              1.051 
    ##              GNRH1              FBXO4              MRPL9           C1orf167 
    ##              1.051              1.052              1.052              1.052 
    ##              RAD9A              ATP4B             LYSMD1              NMUR2 
    ##              1.052              1.052              1.052              1.052 
    ##              IGSF8              CXCR4              IRAK4           KIAA1191 
    ##              1.052              1.052              1.052              1.052 
    ##             ZNF555              ATOH1            UGT2B17              GSDMB 
    ##              1.052              1.052              1.052              1.053 
    ##               PDX1             TRIM72             RPL7L1              ARMC9 
    ##              1.053              1.053              1.053              1.053 
    ##            SLCO2A1           CXorf40A            SLC22A7            SLC16A4 
    ##              1.053              1.053              1.053              1.053 
    ##              KIF17              KANK3             KCNJ12            SLC15A2 
    ##              1.053              1.053              1.053              1.053 
    ##              SOX18             TVP23B             ARMC12               FAAH 
    ##              1.053              1.053              1.053              1.054 
    ##     RP11-426L16.10                CAT              SNX15           C17orf53 
    ##              1.054              1.054              1.054              1.054 
    ##             CNPPD1           ITPRIPL2           C15orf41               TYW3 
    ##              1.054              1.054              1.055              1.055 
    ##              GNG12            C9orf43            C8orf58             LILRB1 
    ##              1.055              1.055              1.055              1.055 
    ##                CGA               DMKN              ATP6C              ACSM3 
    ##              1.055              1.055              1.055              1.055 
    ##               ZBBX           FAM188B2             ZNF213              DCP1B 
    ##              1.055              1.055              1.055              1.055 
    ##             ALOXE3               LCN9              SAMD7              IL5RA 
    ##              1.055              1.055              1.055              1.055 
    ##           SLC25A10               POLL           C11orf85               PIGC 
    ##              1.055              1.056              1.056              1.056 
    ##             ZFAND4               MPP4            ANAPC10              SYTL1 
    ##              1.056              1.056              1.056              1.056 
    ##               QDPR            OLFML2B             TWIST1         AL356356.1 
    ##              1.056              1.056              1.056              1.056 
    ##            FAM153C               FSHB                OSM            PKHD1L1 
    ##              1.056              1.056              1.056              1.057 
    ##             ENTHD1              TRPV4              USP45             CCDC15 
    ##              1.057              1.057              1.057              1.057 
    ##           C19orf12               GHRH             LGALS8              NRBF2 
    ##              1.057              1.057              1.057              1.058 
    ##               PEF1             SMIM19              TOR1A             GPR146 
    ##              1.058              1.058              1.058              1.058 
    ##               BROX              SFRP2           C21orf59             SPDYE3 
    ##              1.058              1.058              1.058              1.058 
    ##            SIGLEC8              BOLA3               GMPR               PKP2 
    ##              1.058              1.059              1.059              1.059 
    ##            PCDHB10              POLG2               CTGF              UGGT2 
    ##              1.059              1.059              1.059              1.059 
    ##           SH3BGRL3              HOXD4            METTL17            C1orf65 
    ##              1.059              1.059              1.059              1.059 
    ##           C9orf142              RIMS3             FGFBP1               RNLS 
    ##              1.060              1.060              1.060              1.060 
    ##               AQP5             TAMM41               KRT4               FMO2 
    ##              1.060              1.060              1.060              1.060 
    ##            TBC1D21                VTN               ELP6              DSCR3 
    ##              1.060              1.060              1.061              1.061 
    ##             GPR183             ABHD10              MYH15              MUC3A 
    ##              1.061              1.061              1.061              1.061 
    ##               GPX8               RHOH         CDKN2AIPNL               DSN1 
    ##              1.061              1.061              1.062              1.062 
    ##              CRYL1               ARSB            SPATA16            C2orf49 
    ##              1.062              1.062              1.062              1.062 
    ##             SPTBN5              PPM1K              FMO6P               ACY1 
    ##              1.062              1.062              1.062              1.063 
    ##              TPGS2             GPR108              ACOX3             DDRGK1 
    ##              1.063              1.063              1.063              1.063 
    ##                ME2              BTBD6               CPN1            SEC14L6 
    ##              1.063              1.063              1.063              1.063 
    ##             CCDC38            C2orf15             MRPL30               PIGZ 
    ##              1.063              1.063              1.063              1.063 
    ##              FCRL5            PPP2R1B               PIGB              NOSIP 
    ##              1.063              1.063              1.063              1.064 
    ##               NOX1               NEFH            GOLPH3L              NARFL 
    ##              1.064              1.064              1.064              1.064 
    ##                SST               GMFB             MYO15B             UGT1A6 
    ##              1.064              1.064              1.064              1.064 
    ##                EVC               EVC2            GPATCH2            TCERG1L 
    ##              1.064              1.064              1.064              1.064 
    ##              DOC2B            APOBEC2               PLD2             SFT2D1 
    ##              1.065              1.065              1.065              1.065 
    ##               AATF              ILDR1              ATOX1             RABAC1 
    ##              1.065              1.065              1.065              1.066 
    ##               BBS2             ADORA3               VWA1             CCDC42 
    ##              1.066              1.066              1.066              1.066 
    ##             CPSF4L             ZNF577              HOXD9            HMGCLL1 
    ##              1.066              1.066              1.066              1.066 
    ##             PCDHA6              MFSD8              KCTD6               DLX4 
    ##              1.067              1.067              1.067              1.067 
    ##            SLC5A11              MIS12               KRT8              ABCA7 
    ##              1.067              1.067              1.067              1.067 
    ##               MED7             TCEAL1              DYRK3             LMBR1L 
    ##              1.067              1.067              1.067              1.067 
    ##              WDR93              FSCN3            NDUFAF3             POLR3F 
    ##              1.067              1.067              1.067              1.067 
    ##              CCER1             TRIM17             ZNF548             DNAAF2 
    ##              1.068              1.068              1.068              1.068 
    ##              BATF2              RGPD8              SOCS1             GALNT9 
    ##              1.068              1.068              1.068              1.068 
    ##           SPATA5L1              TRPM1               IRX1             FAM83C 
    ##              1.068              1.068              1.069              1.069 
    ##           C19orf59              AIPL1               GPN3              GPA33 
    ##              1.069              1.069              1.069              1.069 
    ##             ANKRA2             TCEANC       RP11-371E8.4      RP11-295D22.1 
    ##              1.069              1.069              1.069              1.070 
    ##              CECR1             SPRED2               CD83           ARHGEF37 
    ##              1.070              1.070              1.070              1.070 
    ##              ABCA9            FAM194B              SUMF1            DNAJC25 
    ##              1.070              1.070              1.070              1.070 
    ##              PMPCA               ADSL            TMEM147             TTLL11 
    ##              1.070              1.071              1.071              1.071 
    ##               MINA              MYOZ2              USP44               DOK2 
    ##              1.071              1.071              1.071              1.071 
    ##            CCDC127            CCDC146              RCAN1              ANKAR 
    ##              1.071              1.072              1.072              1.072 
    ##               ART4                SDS              TAAR6              NOMO2 
    ##              1.072              1.072              1.072              1.072 
    ##            PTGES3L              PNMA2           C10orf53              CAPN2 
    ##              1.072              1.072              1.072              1.072 
    ##               CPVL             OCIAD1                ME3              VAT1L 
    ##              1.072              1.072              1.072              1.072 
    ##              POLD4              ABHD4             TM4SF1                CD7 
    ##              1.073              1.073              1.073              1.073 
    ##             NPEPL1              LPAR5              SGPP1               ERI1 
    ##              1.073              1.073              1.073              1.073 
    ##              ACSM1      RP11-644F5.10              ITLN1              ETNK2 
    ##              1.074              1.074              1.074              1.074 
    ##               MICB             PRSS36              SF3B5             UAP1L1 
    ##              1.074              1.074              1.074              1.074 
    ##            DPY19L4             TCHHL1             TTLL12                GCA 
    ##              1.074              1.074              1.074              1.074 
    ##              WDR72              MMP17           TMEM185A               VRK3 
    ##              1.074              1.074              1.074              1.074 
    ##              DDX58           PRICKLE4             CYB5RL             LGALS4 
    ##              1.075              1.075              1.075              1.075 
    ##               PIGW              RAMP2              BSPH1              GINS4 
    ##              1.075              1.075              1.075              1.075 
    ##            MAB21L3             R3HDM4             MBOAT1              PTCRA 
    ##              1.075              1.075              1.075              1.075 
    ##            C5orf34            FAM179A              LACC1              AP1M2 
    ##              1.075              1.075              1.076              1.076 
    ##              GPR82              UQCRB            GRAMD1C           C16orf96 
    ##              1.076              1.076              1.076              1.076 
    ##              MSTO1             POPDC2              SMYD4             ZNF334 
    ##              1.076              1.076              1.076              1.076 
    ##             LRRC38       RP11-529K1.3           ALDH16A1            CXorf38 
    ##              1.076              1.076              1.076              1.076 
    ##              MOB3B        ARPC4-TTLL3               IDI1              NTMT1 
    ##              1.077              1.077              1.077              1.077 
    ##              TTLL6             USHBP1              PANX3               ART5 
    ##              1.077              1.077              1.077              1.077 
    ##              DRAM2              PTGR2              ADAP2           C11orf21 
    ##              1.077              1.077              1.078              1.078 
    ##             ANXA13            TMEM232              PCID2            SLC17A4 
    ##              1.078              1.078              1.078              1.078 
    ##               IDO1              MITD1               TJP3             MAN2C1 
    ##              1.078              1.078              1.078              1.078 
    ##             PNPLA7               DPH6              AP1S3            TMEM136 
    ##              1.078              1.079              1.079              1.079 
    ##                MMD             ZNF782            FAM198B             CHRNA1 
    ##              1.079              1.079              1.079              1.079 
    ##               FZD6            FAM228A               SUOX             GUCY2C 
    ##              1.079              1.079              1.079              1.079 
    ##             SPATA5            ATG16L2              DMBT1               PNMT 
    ##              1.079              1.080              1.080              1.080 
    ##                NPY              RBPJL             TRIM45              TNNC2 
    ##              1.080              1.080              1.080              1.080 
    ##               NGDN              SMUG1             UGT1A1              CEP55 
    ##              1.080              1.080              1.080              1.080 
    ##             PEX11B               DFFB              SYNE3              RLBP1 
    ##              1.080              1.080              1.080              1.080 
    ##            ANKRD35              NEIL3              EEF2K              ABCB8 
    ##              1.081              1.081              1.081              1.081 
    ##              PTCD3            CYP27C1            C5orf15                ACE 
    ##              1.081              1.081              1.081              1.081 
    ##               MYH1           CDK5RAP3            PTTG1IP            TMEM208 
    ##              1.081              1.081              1.081              1.081 
    ##                LOR              EXTL1              NDOR1             FBXO40 
    ##              1.081              1.082              1.082              1.082 
    ##             CLDN19             YAE1D1               NHP2             ALS2CL 
    ##              1.082              1.082              1.082              1.082 
    ##             OR11L1               PTRF             ACOT11             SLFNL1 
    ##              1.082              1.082              1.082              1.082 
    ##             DPAGT1               CDX1              MOCS1              QTRT1 
    ##              1.082              1.082              1.082              1.082 
    ##               MZF1             FCGR1A             CHCHD5              ACTR6 
    ##              1.082              1.083              1.083              1.083 
    ##              AIFM2                RHO               GPC5            LACTBL1 
    ##              1.083              1.083              1.083              1.083 
    ##               AQP1             CCDC62                FRK              TRPV1 
    ##              1.083              1.083              1.084              1.084 
    ##      ZNF559-ZNF177              OSGEP              CCL15               SHPK 
    ##              1.084              1.084              1.084              1.084 
    ##              RGS18              SEPT1               THY1              FUT10 
    ##              1.084              1.084              1.084              1.084 
    ##              LRRC6               DOHH               ARSK                 GC 
    ##              1.084              1.084              1.084              1.084 
    ##               XCR1              IGSF6             ACAD10            PSMC3IP 
    ##              1.084              1.085              1.085              1.085 
    ##              TEKT1               PIRT              SMPD1               RFX8 
    ##              1.085              1.085              1.085              1.085 
    ##            SLC35E4              TAF13             BTN1A1               NIT1 
    ##              1.085              1.085              1.085              1.085 
    ##              SPIDR            C5orf51              ERCC2              CD151 
    ##              1.085              1.086              1.086              1.086 
    ##             MTMR11              CDKL1              CLIC2              RWDD3 
    ##              1.086              1.086              1.086              1.086 
    ##               RHOD              WDR64            TBC1D26            SLC37A3 
    ##              1.086              1.086              1.086              1.086 
    ##               GRK4              NHLH2            CCDC172             ACOT12 
    ##              1.087              1.087              1.087              1.087 
    ##               AQP6              NTSR2              CCDC3              USH1G 
    ##              1.087              1.087              1.087              1.087 
    ##             RPUSD2              REXO4            C6orf89       RP11-58C22.1 
    ##              1.087              1.087              1.087              1.088 
    ##              PDCL2            SLC6A16             CASP12              EDEM1 
    ##              1.088              1.088              1.088              1.088 
    ##               FGF4               AACS            FAM161B              RPS21 
    ##              1.088              1.088              1.088              1.088 
    ##              TACO1             FRRS1L            CYP2A13              CASP6 
    ##              1.088              1.088              1.089              1.089 
    ##             LACTB2              APPL2              TRAF5             ZNF627 
    ##              1.089              1.089              1.089              1.089 
    ##             PTPN22           EBNA1BP2            HSD17B8              FKBP2 
    ##              1.089              1.089              1.089              1.089 
    ##              CEP19              TMOD4              DDX52               DSG3 
    ##              1.090              1.090              1.090              1.090 
    ##              SPDL1             STXBP4            NAALAD2            SPATA18 
    ##              1.090              1.090              1.090              1.091 
    ##              PXDNL            FAM168B             CUEDC1               TTPA 
    ##              1.091              1.091              1.091              1.091 
    ##            HSD17B3            SIGLEC7               VILL            SLC34A2 
    ##              1.091              1.091              1.091              1.091 
    ##              AP5M1             PDCD2L              FANCG               CIB2 
    ##              1.091              1.091              1.091              1.091 
    ##               ATF4             DUSP11             LYPLA1             ZNF621 
    ##              1.091              1.091              1.091              1.092 
    ##              ITIH2              CABP5             ZNF135              GALK1 
    ##              1.092              1.092              1.092              1.092 
    ##                KEL            CCDC71L             EXOSC7               WIF1 
    ##              1.092              1.092              1.092              1.092 
    ##              GMPPB              MED10              PRRX2              KRT27 
    ##              1.092              1.092              1.092              1.092 
    ##               NINL            FAM217A           ADAMTSL4               FMO1 
    ##              1.092              1.093              1.093              1.093 
    ##                CA9            C1QTNF1              SPERT              SNX22 
    ##              1.093              1.093              1.093              1.093 
    ##              HAUS3              XRCC4             CLDN12            NOSTRIN 
    ##              1.093              1.094              1.094              1.094 
    ##             SAP30L               RFT1               RD3L               C1QC 
    ##              1.094              1.094              1.094              1.094 
    ##             ATP8B4             SH2D2A              SH2D5      STON1-GTF2A1L 
    ##              1.094              1.094              1.094              1.094 
    ##              ALDH2              PTPN7              SEPP1            NDUFB10 
    ##              1.094              1.094              1.094              1.094 
    ##            TRMT112             CCDC94               RAD1            C8orf34 
    ##              1.094              1.095              1.095              1.095 
    ##               ULK4            MID1IP1              SNX31              TUBA8 
    ##              1.095              1.095              1.095              1.095 
    ##              RIOK1               MSMP             ZNF620                LPL 
    ##              1.095              1.095              1.095              1.095 
    ##              SCNM1              KRT25               PNCK             EXOSC2 
    ##              1.095              1.096              1.096              1.096 
    ##               UCK1              CLIC6              MUC20            GPIHBP1 
    ##              1.096              1.096              1.096              1.096 
    ##              IFT57              PPM1D              SUMO3       CTC-435M10.3 
    ##              1.096              1.096              1.096              1.096 
    ##              FOLR3            TMPRSS7               ARL3             ARMC10 
    ##              1.096              1.097              1.097              1.097 
    ##               QPRT            GALNT14           BAIAP2L2              ACOXL 
    ##              1.097              1.097              1.097              1.097 
    ##            NDUFA10             CCDC33             BDKRB1             CYP4B1 
    ##              1.097              1.097              1.097              1.097 
    ##            BCL2L15           ARHGAP19            C1orf74             SLC2A8 
    ##              1.097              1.097              1.098              1.098 
    ##              RPL39               GGT5               MED8             SCCPDH 
    ##              1.098              1.098              1.098              1.098 
    ##            ZNF324B               SMG9              HTR1F               ROM1 
    ##              1.098              1.098              1.099              1.099 
    ##             CYP2W1              CCL25             MRPL38               PRR3 
    ##              1.099              1.099              1.099              1.099 
    ##              TRPM4             DEPDC7              VSIG2               ASPA 
    ##              1.099              1.099              1.099              1.099 
    ##               RRAS              PTGIS              AFMID            CCDC90B 
    ##              1.099              1.099              1.099              1.099 
    ##              ALDOB                FAP            CLEC12B             SLC9A4 
    ##              1.100              1.100              1.100              1.100 
    ##               SWI5              KLK15             ADORA1             CDADC1 
    ##              1.100              1.100              1.100              1.100 
    ##             MIPOL1               PIDD             ZNF677              PRODH 
    ##              1.100              1.100              1.100              1.100 
    ##              P2RX6             RPUSD3                NRL              YIF1B 
    ##              1.100              1.100              1.101              1.101 
    ##             ZNF630               HLTF              TTC38             MARCKS 
    ##              1.101              1.101              1.101              1.101 
    ##              KRT34             SMIM20              OLFM4             CCDC83 
    ##              1.102              1.102              1.102              1.102 
    ##             CDC20B               MAL2            MTERFD2              CCM2L 
    ##              1.102              1.102              1.102              1.102 
    ##              DHX58               MRI1               COLQ             GPR149 
    ##              1.102              1.102              1.102              1.102 
    ##            C1orf52             ZNF560              PLCD1              RBM18 
    ##              1.103              1.103              1.103              1.103 
    ##               AARD             MRPS30             BCKDHB              NTAN1 
    ##              1.103              1.103              1.103              1.104 
    ##               PNOC               POLM             CCDC53              CENPM 
    ##              1.104              1.104              1.104              1.104 
    ##          LINC00692           HIST1H1B              TGIF2               RFX5 
    ##              1.104              1.104              1.104              1.104 
    ##               KIF6            ANKRD27           KIAA0753             PNPLA3 
    ##              1.105              1.105              1.105              1.105 
    ##               STX2            FAM92A1              CDH15             KCNK18 
    ##              1.105              1.105              1.105              1.105 
    ##             PMFBP1              RGS11               DAP3              FBXL6 
    ##              1.105              1.105              1.105              1.105 
    ##              ODF3B             AGPAT9              TXLNB              ZFP69 
    ##              1.105              1.106              1.106              1.106 
    ##              MYOM3              WDR88            C2orf71            SLC26A2 
    ##              1.106              1.106              1.106              1.106 
    ##              KLK13              FABP5            SLC35F6               GLDC 
    ##              1.106              1.106              1.106              1.106 
    ##            BCL2L14              IFRD2              CPNE1            TMEM138 
    ##              1.106              1.106              1.106              1.106 
    ##            TSPAN13               TDP1              RMND1              S1PR3 
    ##              1.106              1.107              1.107              1.107 
    ##           SLC25A38              SUSD2             ZNF438             ZNF446 
    ##              1.107              1.107              1.107              1.107 
    ##             GPR180                TPO             CYP4X1             LHFPL1 
    ##              1.107              1.107              1.107              1.107 
    ##              ADPRM              ATG9B        RP4-559A3.7              MTFR2 
    ##              1.107              1.107              1.107              1.107 
    ##            XPNPEP2              ADCK4            FAM129C             RETSAT 
    ##              1.107              1.108              1.108              1.108 
    ##        CTC-360G5.8               HSF4              TRIM6             AARSD1 
    ##              1.108              1.108              1.108              1.108 
    ##              ESYT3             CLEC4G              GDPD4           C15orf60 
    ##              1.109              1.109              1.109              1.109 
    ##           C16orf87            ZFYVE27              AQP10              SNX11 
    ##              1.109              1.109              1.109              1.109 
    ##               AMBN            C2orf74              RAD52            SNRNP25 
    ##              1.109              1.109              1.109              1.109 
    ##              FRAT1             FAM45A             EXOSC9            B4GALT7 
    ##              1.109              1.110              1.110              1.110 
    ##              CDK18            METAP1D           C19orf44              WDR41 
    ##              1.110              1.111              1.111              1.111 
    ##              RPE65             XXYLT1            ZCCHC12              SSC5D 
    ##              1.111              1.111              1.111              1.111 
    ##           FAM189A2              RPP30            CYP3A43                LY9 
    ##              1.112              1.112              1.112              1.112 
    ##               RND2             UGT1A9              VMA21             CELA2B 
    ##              1.112              1.112              1.112              1.112 
    ##               EMC9             CRISP3              SPO11              WDR65 
    ##              1.112              1.112              1.112              1.112 
    ##              AP3M2              TRPM2             LRRC70              HMGCL 
    ##              1.112              1.112              1.112              1.112 
    ##             SH3YL1             MAPK13             CAPN14             HS3ST1 
    ##              1.112              1.112              1.112              1.113 
    ##              TLCD1           METTL11B              DCST1              KRBA1 
    ##              1.113              1.113              1.113              1.113 
    ##               CMBL               TTC1              HMCES             PAFAH2 
    ##              1.113              1.113              1.113              1.113 
    ##             CAPNS1             SLC4A9             TALDO1              TUBB3 
    ##              1.113              1.114              1.114              1.114 
    ##              TIMP4            CCDC151              DUS4L             TRIM10 
    ##              1.114              1.114              1.114              1.114 
    ##              RHPN1             ZNF707              TM2D1               AIRE 
    ##              1.114              1.114              1.115              1.115 
    ##             DUSP27             DYX1C1             LRRC57             IL10RB 
    ##              1.115              1.115              1.115              1.115 
    ##           SERPINB6               PUS3           MAPKAPK3             TRIM38 
    ##              1.115              1.115              1.115              1.115 
    ##              ACMSD             CLCNKA             MRPL37              MORN1 
    ##              1.116              1.116              1.116              1.116 
    ##              TTC19               ZXDA             ZNF558               ROS1 
    ##              1.116              1.116              1.116              1.116 
    ##               LY6E           SERPINE3               WDR4            CLEC11A 
    ##              1.116              1.116              1.116              1.116 
    ##            FAM229B             FBXW10             THAP10             TRMT11 
    ##              1.117              1.117              1.117              1.117 
    ##             ZNF177               VWDE             PCDHB5               FAN1 
    ##              1.117              1.117              1.117              1.117 
    ##              PANX1              FREM1               EMR1             TMBIM1 
    ##              1.117              1.117              1.117              1.117 
    ##      RP11-618P17.4            SERTAD1              ACAA1               HMMR 
    ##              1.118              1.118              1.118              1.118 
    ##             FBXL15               UPP1             CLEC2L            CCDC178 
    ##              1.118              1.118              1.118              1.118 
    ##             IGFBP7            TMEM260              GSTA3              KCTD8 
    ##              1.118              1.118              1.118              1.118 
    ##               UGP2             ZNF396              SNTG2              RBM43 
    ##              1.118              1.118              1.118              1.119 
    ##               GBAS               VNN3              BFSP1             LILRB2 
    ##              1.119              1.119              1.119              1.119 
    ##            OCSTAMP            TMPRSS3              CHST8             CSF2RA 
    ##              1.119              1.119              1.119              1.119 
    ##              ITLN2                NMI              CDH19             PTP4A3 
    ##              1.119              1.119              1.119              1.119 
    ##              FGF19              PXDC1            CCDC64B             INSIG1 
    ##              1.119              1.119              1.119              1.119 
    ##             LRRC42             PCDHB2                DUT              PDE6A 
    ##              1.119              1.119              1.120              1.120 
    ##               TSR1             RFXANK             GLB1L2              MED21 
    ##              1.120              1.120              1.120              1.120 
    ##             MRPL55             DNAH12             SLC7A9            SLC37A2 
    ##              1.120              1.120              1.121              1.121 
    ##             AKR1A1             PPP5D1             AVPR1B              STON1 
    ##              1.121              1.121              1.121              1.121 
    ##              THAP5               EXD2            TMEM102              PP2D1 
    ##              1.121              1.121              1.121              1.121 
    ##              ASAH2             GAS2L3             CCDC65               PHYH 
    ##              1.122              1.122              1.122              1.122 
    ##             FERD3L             UMODL1             GPR128             IFNAR2 
    ##              1.122              1.122              1.122              1.122 
    ##           KIAA1586              FCRL2              WDR74             ZNF283 
    ##              1.122              1.122              1.122              1.122 
    ##              PDCD7              RAMP1              LEKR1             CYP3A4 
    ##              1.122              1.122              1.122              1.122 
    ##           CATSPER4             LRRIQ4             TNFSF9              ACTL8 
    ##              1.123              1.123              1.123              1.123 
    ##              CCT6B               HUS1             IMPACT              PRRG1 
    ##              1.123              1.123              1.123              1.123 
    ##               AQP8             CLEC4M             SIRPB2              LOXL4 
    ##              1.123              1.123              1.123              1.123 
    ##              MRPS7               TFPT              RAB31           DNASE1L3 
    ##              1.123              1.123              1.123              1.123 
    ##             SLC5A2            TBC1D24            PLEKHN1               DPM2 
    ##              1.123              1.123              1.124              1.124 
    ##               TMC8             TRIM21            NANOGNB              GPAT2 
    ##              1.124              1.124              1.124              1.124 
    ##            DCAF4L1              STOX1          C10orf107         AC024580.1 
    ##              1.124              1.124              1.124              1.124 
    ##            CCDC105             ATRAID             ZNF584            CCDC150 
    ##              1.124              1.124              1.124              1.125 
    ##              DFNA5               RGP1         ST6GALNAC1              SURF6 
    ##              1.125              1.125              1.125              1.125 
    ##            C4orf27                SHD          KRTAP16-1              TRUB1 
    ##              1.125              1.125              1.125              1.125 
    ##               GYS2               UROD              CARD9             COQ10A 
    ##              1.125              1.125              1.126              1.126 
    ##              CKAP2               CD38           LOH12CR1             LYPD6B 
    ##              1.126              1.126              1.126              1.126 
    ##              WISP1           SLC25A24           KIAA1919               ECM1 
    ##              1.126              1.126              1.126              1.126 
    ##           CCNB1IP1             CAPN13             AKR1D1            RASL10B 
    ##              1.126              1.126              1.127              1.127 
    ##              MSRB1              IQCF2             CEP112              ABCB5 
    ##              1.127              1.127              1.127              1.127 
    ##               CCR9              EPHX1              CCL19             KLHL21 
    ##              1.127              1.127              1.127              1.127 
    ##              RTDR1            SEC14L4             ARRDC4               KLF1 
    ##              1.127              1.128              1.128              1.128 
    ##               LDLR             GRPEL2             TTC30B              SPCS1 
    ##              1.128              1.128              1.128              1.128 
    ##            TAS2R38                HPX           TMEM229B              ACBD6 
    ##              1.128              1.128              1.129              1.129 
    ##              LCN15              CABYR              ELANE                DAK 
    ##              1.129              1.129              1.129              1.129 
    ##             ANKUB1              ANO10            NPIPB11                GSC 
    ##              1.129              1.129              1.129              1.129 
    ##           DYNC2LI1             UBXN10               MEFV             TRIM16 
    ##              1.129              1.129              1.129              1.129 
    ##             TEX261              NOL12               RCN1             ZBTB45 
    ##              1.130              1.130              1.130              1.130 
    ##               ORC4             NPFFR2             SETMAR               RTCA 
    ##              1.130              1.130              1.130              1.130 
    ##              TSTA3              MXRA8               ECH1                NYX 
    ##              1.130              1.130              1.130              1.131 
    ##             MEGF11              DHX35         AC109829.1             SMIM12 
    ##              1.131              1.131              1.131              1.131 
    ##               FHL2              NUCB2             MVB12A            GADD45A 
    ##              1.131              1.131              1.131              1.131 
    ##               SDC4            FAM188B            TMEM175             ZNF891 
    ##              1.131              1.131              1.131              1.131 
    ##               IGF2               PGS1              ERAP1              CXCR2 
    ##              1.131              1.131              1.132              1.132 
    ##               MYF6               GSX1             TTLL10              BTBD8 
    ##              1.132              1.132              1.132              1.132 
    ##               MFI2               ASTL              OR2W1                FUK 
    ##              1.132              1.132              1.132              1.132 
    ##               OMA1             RNF152             SLC5A4              APOL5 
    ##              1.132              1.133              1.133              1.133 
    ##               NCF1               CLN6             SPINK1        AC007040.11 
    ##              1.133              1.133              1.133              1.133 
    ##             ZBTB25               HBG1             LRRC20            SLC47A2 
    ##              1.133              1.133              1.133              1.133 
    ##             RNASEL              PDSS2               OPA3             UGT1A4 
    ##              1.133              1.133              1.134              1.134 
    ##               TDP2               SHQ1               COG8             GLT8D1 
    ##              1.134              1.134              1.134              1.134 
    ##           C17orf66              PDZD7             CXCL11             SH3BP2 
    ##              1.134              1.134              1.135              1.135 
    ##               TFEC             MRFAP1              BRMS1            CLEC10A 
    ##              1.135              1.135              1.136              1.136 
    ##              PADI1              PADI4            TMEM30C              PVRIG 
    ##              1.136              1.136              1.136              1.136 
    ##             HOXC11              KRT28             CHRNA5               RIC3 
    ##              1.136              1.136              1.136              1.136 
    ##             COMMD3             OR10K1            PIK3C2G           HS3ST3A1 
    ##              1.136              1.136              1.136              1.136 
    ##              PINK1              ACTL9             CCDC58               GUF1 
    ##              1.136              1.137              1.137              1.137 
    ##            B4GALT4            C9orf24                NPL               FUT4 
    ##              1.137              1.137              1.137              1.137 
    ##            CCDC85B              OVOL3               PNKD              POF1B 
    ##              1.137              1.137              1.138              1.138 
    ##              LMOD3                SPR              ASTE1            WFIKKN1 
    ##              1.138              1.138              1.138              1.138 
    ##            ALDH9A1             ZNF502               NOX4             IL17RB 
    ##              1.138              1.138              1.138              1.139 
    ##           ITGB1BP1              CENPH       RAD51L3-RFFL           HSD17B13 
    ##              1.139              1.139              1.139              1.139 
    ##       CTC-512J12.6           C6ORF165               NEBL             ZNF185 
    ##              1.139              1.139              1.139              1.139 
    ##           FAM177A1            TRMT10C             FBXO44             FAM81B 
    ##              1.139              1.139              1.139              1.139 
    ##              HEXDC               ASPG             EXOSC5               LZIC 
    ##              1.140              1.140              1.140              1.140 
    ##              CAPN9              HSPB8            HSPBAP1               RDH8 
    ##              1.140              1.140              1.140              1.140 
    ##          C14orf164              CALR3      CTD-2260A17.2              OSCP1 
    ##              1.140              1.140              1.140              1.140 
    ##              TMA16             FAM64A               NAGS               ANO6 
    ##              1.140              1.141              1.141              1.141 
    ##            PCDHGB2           C11orf70              HTR3A              COX7C 
    ##              1.141              1.141              1.141              1.141 
    ##               HAAO              PAPLN              IL17A               SKA3 
    ##              1.141              1.141              1.141              1.141 
    ##               STOM             SAMD12             SLC2A9           C15orf52 
    ##              1.141              1.141              1.141              1.142 
    ##             UBXN11            PPP1R36              KCNJ5              AGBL1 
    ##              1.142              1.142              1.142              1.142 
    ##             CCDC37               CBLC             GRAMD3            ZFYVE19 
    ##              1.142              1.142              1.142              1.143 
    ##               KIF9              CDCA5               DMP1              ENPEP 
    ##              1.143              1.143              1.143              1.143 
    ##               RBFA              SARM1            TERF2IP             DALRD3 
    ##              1.143              1.143              1.143              1.143 
    ##              CCER2           PPP1R14A               PRG3            C2orf61 
    ##              1.143              1.143              1.143              1.143 
    ##              SCRN3              ACAT2             GPR113               DQX1 
    ##              1.143              1.143              1.143              1.143 
    ##              HDHD1               SCTR              ADAD2              TSEN2 
    ##              1.143              1.144              1.144              1.144 
    ##            C3orf20              CSRP1              LRRC2              SNAI3 
    ##              1.144              1.144              1.144              1.144 
    ##      RP11-385D13.1              TRIT1               HYKK              KRT84 
    ##              1.144              1.144              1.145              1.145 
    ##                MPO            SLC22A8              VAMP1               HAO1 
    ##              1.145              1.145              1.145              1.145 
    ##             ZNF786             TATDN3               CES3           CSNK1A1L 
    ##              1.145              1.145              1.146              1.146 
    ##               IQCK             GPR112              MASP2             TIMM21 
    ##              1.146              1.146              1.146              1.146 
    ##            SLC13A2               CSH2              TRPV3                ESD 
    ##              1.146              1.146              1.146              1.146 
    ##             MRPS24            ANKRD61              TTC16             TRIM69 
    ##              1.146              1.146              1.146              1.146 
    ##             ZNF432               TSR3              NBPF1              HMGN5 
    ##              1.146              1.147              1.147              1.147 
    ##              MRPL2             SH3RF2               ARSG               CIB1 
    ##              1.147              1.147              1.147              1.148 
    ##            KIR2DL4       ATP5J2-PTCD1               TGDS              PTCD1 
    ##              1.148              1.148              1.148              1.148 
    ##              CDKL3             CHRNA3             RILPL2             SPATA7 
    ##              1.148              1.148              1.148              1.149 
    ##             SCPEP1            FAM195B                CA4             DFNB59 
    ##              1.149              1.149              1.149              1.149 
    ##               PLK5               GCKR             SH2D1B              CD276 
    ##              1.149              1.149              1.149              1.149 
    ##              IQCF1             ZNF195              HOXC9             ZNF738 
    ##              1.150              1.150              1.150              1.150 
    ##             FUNDC1               TIFA               TPMT               VWA8 
    ##              1.150              1.150              1.150              1.150 
    ##            CWF19L1              TSSK2                LPO           ADAMTSL5 
    ##              1.150              1.150              1.150              1.151 
    ##              OTOP3             CAPN11              NBPF3            ZSCAN31 
    ##              1.151              1.151              1.151              1.151 
    ##            CCDC154               TMC2           ARHGEF19              MFAP3 
    ##              1.151              1.151              1.151              1.151 
    ##               RPTN               TLL2             GCSAML           TMPRSS15 
    ##              1.151              1.151              1.151              1.151 
    ##               BAG2               PRCP            EXOC3L1                PNP 
    ##              1.151              1.152              1.152              1.152 
    ##             COMMD4              FARS2               PDPN               MUL1 
    ##              1.152              1.152              1.152              1.152 
    ##          SERPINB13             MBOAT4            EDARADD               C1QA 
    ##              1.152              1.153              1.153              1.153 
    ##               PFKP              ACADM             PLA2G3              CNGB1 
    ##              1.153              1.153              1.153              1.153 
    ##               SCLY             CHRNA6             PLSCR4            SLC35A5 
    ##              1.153              1.153              1.153              1.153 
    ##              WHAMM              ARMC4             RNF141              CNTN6 
    ##              1.153              1.153              1.153              1.153 
    ##             ZNF623              MOCOS              YPEL3            SLC22A2 
    ##              1.153              1.154              1.154              1.154 
    ##              SPRY3             TWIST2              NR1I3            FAM183A 
    ##              1.154              1.154              1.154              1.154 
    ##            TMEM127               FIGF            FASTKD3             ZNF165 
    ##              1.154              1.154              1.154              1.154 
    ##              OFCC1               NEK4            LEPREL1            MTHFD2L 
    ##              1.154              1.154              1.155              1.155 
    ##           C18orf54               C1RL           HLA-DRB1               HRH1 
    ##              1.155              1.156              1.156              1.156 
    ##            NCR3LG1           C9orf117              DIRC2               OXSM 
    ##              1.156              1.156              1.156              1.156 
    ##               OIT3             IZUMO2              DCDC2            DCSTAMP 
    ##              1.156              1.156              1.156              1.156 
    ##               ISCU            PCDHA13              LRRC9              TACR2 
    ##              1.156              1.156              1.157              1.157 
    ##               XKRX             KHDC3L           C16orf97              MSRB3 
    ##              1.157              1.157              1.157              1.157 
    ##              OR8D4             NDUFS4              TRPA1           SLC9A3R1 
    ##              1.157              1.157              1.157              1.157 
    ##           C11orf88              ABCC2             PTPMT1             ORMDL1 
    ##              1.158              1.158              1.158              1.158 
    ##             TFAP2E               MSH3             SBSPON             EFCAB6 
    ##              1.158              1.158              1.158              1.158 
    ##               STX4              H3F3A               SIAE              EPCAM 
    ##              1.158              1.158              1.159              1.159 
    ##              CLIC4             ZNF226            SLC12A9             ZNF595 
    ##              1.159              1.159              1.159              1.159 
    ##               MLNR           EIF4EBP2               HPSE               BAG5 
    ##              1.159              1.159              1.159              1.159 
    ##           SLC25A10           SLC25A29            FAM153A              CEP44 
    ##              1.159              1.159              1.159              1.160 
    ##              CMTR2            HSD17B7            FAM169B               TPBG 
    ##              1.160              1.160              1.160              1.160 
    ##                MRO               GYPA              SLX1A             ZNF675 
    ##              1.160              1.160              1.160              1.160 
    ##             ZNF200             CHTF18             ZNHIT6              GCNT2 
    ##              1.160              1.161              1.161              1.161 
    ##               DFFA               RHOB        AC006547.14              TMED6 
    ##              1.161              1.161              1.161              1.161 
    ##             RPUSD4              ANKK1            PCDHA11             ADAM28 
    ##              1.162              1.162              1.162              1.162 
    ##               MCAT            PLEKHF1              GDAP1               SYNM 
    ##              1.162              1.162              1.162              1.162 
    ##              TULP2               PIGX              OCEL1                ADM 
    ##              1.162              1.162              1.162              1.162 
    ##          SERPINB11                GNE               MLC1              ACADL 
    ##              1.162              1.162              1.162              1.162 
    ##            LDLRAP1            C2orf62               ELL3               PPAN 
    ##              1.163              1.163              1.163              1.163 
    ##        PPAN-P2RY11               CA13               BAG1              FUT11 
    ##              1.163              1.163              1.163              1.163 
    ##            SLC24A5            SLC19A1              NPSR1              MEGT1 
    ##              1.163              1.164              1.164              1.164 
    ##             TRIM62              TIGD1                LTF           SLC22A18 
    ##              1.164              1.164              1.164              1.164 
    ##              TPD52             SAMD15           C19orf55            LGALS9B 
    ##              1.164              1.164              1.165              1.165 
    ##             TSPAN2              CRTAM           ZC3HAV1L           GTF2IRD2 
    ##              1.165              1.165              1.165              1.165 
    ##              CBLN2               NOD1              GRHPR            POLR3GL 
    ##              1.165              1.165              1.165              1.165 
    ##            PCDHB17             CYP2E1             CRYZL1             FAM83D 
    ##              1.165              1.165              1.166              1.166 
    ##            PYROXD2             PKDREJ                HK3             AMDHD1 
    ##              1.166              1.166              1.166              1.166 
    ##               IDUA             DMRTA1            SLC47A1             DZANK1 
    ##              1.166              1.166              1.166              1.167 
    ##            GNPNAT1              ERMAP             DUSP14           FRA10AC1 
    ##              1.167              1.167              1.167              1.167 
    ##           SLC25A20              PDCD5             CCDC51              GCFC2 
    ##              1.167              1.167              1.167              1.167 
    ##            PLA2G4F           APOBEC3A             QTRTD1                BLK 
    ##              1.167              1.168              1.168              1.168 
    ##             ZDHHC1              ITFG3              TMEM9              XIRP1 
    ##              1.168              1.168              1.168              1.168 
    ##               LAX1             LRPAP1             IZUMO1              CASP4 
    ##              1.169              1.169              1.169              1.169 
    ##               DPYS               PPID              PLIN3                MPI 
    ##              1.169              1.169              1.169              1.169 
    ##             COX6B1                ME1             SEPT10             KDELC2 
    ##              1.169              1.169              1.169              1.170 
    ##              RRP15             ZNF764             LINGO4               GRM6 
    ##              1.170              1.170              1.170              1.170 
    ##          TNFRSF10D             DUSP18            SLC16A3             MRPL41 
    ##              1.170              1.170              1.170              1.171 
    ##               NEK5            FAM153B            SNRNP48            TMEM41A 
    ##              1.171              1.171              1.171              1.171 
    ##              DBF4B               PMM1            FAM227A               CPT2 
    ##              1.171              1.171              1.171              1.172 
    ##               ABRA              DMRT3               HELB              LTB4R 
    ##              1.172              1.172              1.172              1.172 
    ##             ACADSB              GSTM3              TKTL2             OR10T2 
    ##              1.172              1.172              1.172              1.172 
    ##                INS             LRRC17              WNT16               GDE1 
    ##              1.172              1.172              1.172              1.172 
    ##            SH3BGRL              PHTF1               ADH6              DMGDH 
    ##              1.172              1.172              1.173              1.173 
    ##             FAM53A             IL1RL1            SLC39A3               CROT 
    ##              1.173              1.173              1.173              1.173 
    ##             IFLTD1              RBM11              PLAC8             FAM26E 
    ##              1.173              1.173              1.174              1.174 
    ##             ATAD3C              EGFL8             LINGO3             GLT6D1 
    ##              1.174              1.174              1.174              1.174 
    ##               ADAL              NR1I2               CCNO              CHADL 
    ##              1.174              1.174              1.174              1.174 
    ##              RAB3C                DCT              CLIC5               CLN8 
    ##              1.174              1.174              1.174              1.174 
    ##               SSPN              FCAMR               DSPP             TAS1R1 
    ##              1.174              1.174              1.174              1.174 
    ##             INPP5J             PRSS27            FDXACB1             MRPL42 
    ##              1.175              1.175              1.175              1.175 
    ##              MYO19              PAQR7               DPCD            SLC26A4 
    ##              1.175              1.175              1.175              1.175 
    ##               ALG2              MED22               FMO4             NUDT17 
    ##              1.175              1.175              1.175              1.175 
    ##           GPATCH11              MED29              NPHS2              CD320 
    ##              1.175              1.175              1.175              1.176 
    ##           C12orf10              IGFL3             CYP1B1            CCDC126 
    ##              1.176              1.176              1.176              1.176 
    ##           C6orf203              ADAM2             FRMPD2              CLDN1 
    ##              1.176              1.176              1.177              1.177 
    ##             CHMP1A            RANBP3L              TMLHE             UTP11L 
    ##              1.177              1.177              1.177              1.177 
    ##              MDH1B               HELT             TRIM55               GBP7 
    ##              1.177              1.177              1.177              1.177 
    ##              CEBPA           C16orf78              CCNB2            BLOC1S2 
    ##              1.177              1.177              1.177              1.178 
    ##           C16orf93            CCDC173               LIPG              NOP16 
    ##              1.178              1.178              1.178              1.178 
    ##              CHST9             NDUFA1             NSMCE1             TTLL13 
    ##              1.178              1.178              1.178              1.178 
    ##            SLC46A3             CCDC89              IL17F           RAD51AP2 
    ##              1.179              1.179              1.179              1.179 
    ##              UBL4A               CTBS             ZCCHC3             STXBP6 
    ##              1.179              1.179              1.179              1.179 
    ##              CXCL1             OR52I1             CXCL12            GOLGA6C 
    ##              1.179              1.179              1.179              1.179 
    ##             NANOS3              CLRN3             GTF3C6              ACOT9 
    ##              1.179              1.179              1.179              1.180 
    ##               ISM2              FOXR1           C20orf78               TGM7 
    ##              1.180              1.180              1.180              1.180 
    ##              FANK1           MARVELD2            ZNF354C              COTL1 
    ##              1.180              1.180              1.180              1.180 
    ##              ASAH1              BCMO1             UGT2A1              SQRDL 
    ##              1.180              1.180              1.180              1.180 
    ##              STK19            METTL2B               NT5E             EEFSEC 
    ##              1.181              1.181              1.181              1.181 
    ##            SLC15A5              ABHD3              ENPP6             HEPHL1 
    ##              1.181              1.181              1.181              1.181 
    ##              CUZD1             ZNF404           KRTAP1-3                OAF 
    ##              1.181              1.181              1.181              1.181 
    ##              PTPRH               ASPN            FAM71F1            PCDHGA7 
    ##              1.182              1.182              1.182              1.182 
    ##            IL13RA2             LURAP1              EFCC1              SSUH2 
    ##              1.182              1.182              1.182              1.182 
    ##             EIF1AY             SEL1L2           TMEM179B            TXNDC17 
    ##              1.182              1.182              1.182              1.182 
    ##             EFCAB1           PNLIPRP2            C6orf10              FARP2 
    ##              1.183              1.183              1.183              1.183 
    ##              TAPBP              KLRG1              OR6A2            SRCRB4D 
    ##              1.183              1.183              1.183              1.183 
    ##               CALY              IFIT3            PCDHGA4               BBS4 
    ##              1.184              1.184              1.184              1.184 
    ##            ALOX15B             METTL6              NINJ1              RASD1 
    ##              1.184              1.184              1.184              1.184 
    ##            TRMT10A                C8A            C2orf48           PNLIPRP3 
    ##              1.184              1.184              1.185              1.185 
    ##               EID2             LMAN1L               TGM2             SPINK4 
    ##              1.185              1.185              1.185              1.186 
    ##               MPC1              STYK1             ARMCX6                TK2 
    ##              1.186              1.186              1.186              1.186 
    ##                HFE             EIF1AD               SGCZ             CYP2R1 
    ##              1.186              1.186              1.186              1.186 
    ##              RAP2C             GPR157              HIBCH              ABCA8 
    ##              1.186              1.186              1.186              1.186 
    ##                ACD              PDE6B               TAC3            THUMPD2 
    ##              1.187              1.187              1.187              1.187 
    ##           C11orf48              WDR31              ADAT1            C1orf35 
    ##              1.187              1.187              1.187              1.187 
    ##              CNTLN            C7orf73              OR2B6             CAPZA3 
    ##              1.187              1.187              1.187              1.187 
    ##               PSPH              USP41             CLEC2B               GAS8 
    ##              1.188              1.188              1.188              1.188 
    ##              ATP7B             ZNF554             CYP8B1            TMEM109 
    ##              1.188              1.188              1.188              1.188 
    ##           HSD11B1L               GDF9             MCOLN2            PTPLAD2 
    ##              1.188              1.188              1.188              1.188 
    ##               ALG3              HYAL3               SIT1             PARP11 
    ##              1.188              1.188              1.188              1.189 
    ##           SLC25A30               PRF1             YEATS4            PCDHA12 
    ##              1.189              1.189              1.189              1.189 
    ##              CD109             TMEM40             PRODH2           ATP6V0E2 
    ##              1.189              1.189              1.189              1.189 
    ##              WDR89      CTD-2410N18.5               PPIH            ABHD16B 
    ##              1.189              1.189              1.189              1.190 
    ##                AK8              ACBD7             PRKAG3             CLDN22 
    ##              1.190              1.190              1.190              1.190 
    ##              CMTM7              CENPA              IQCF5              APOA4 
    ##              1.190              1.190              1.190              1.190 
    ##              HDDC2               GANC             SLC2A5             TRDMT1 
    ##              1.190              1.191              1.191              1.191 
    ##           PLA2G12B              SFXN3              MUTYH            ZKSCAN3 
    ##              1.191              1.191              1.191              1.191 
    ##            CYP26A1            PGLYRP3              EPPIN        EPPIN-WFDC6 
    ##              1.191              1.191              1.191              1.191 
    ##              ERAP2              DDX49               APOM              RNF39 
    ##              1.191              1.191              1.191              1.191 
    ##              DGUOK               MAP9               RPF1             EFCAB7 
    ##              1.192              1.192              1.192              1.192 
    ##              EFNA2              AWAT1              OPRD1             BTN3A1 
    ##              1.192              1.192              1.193              1.193 
    ##             IGFBP1             STARD6              TSTD2            ANGPTL6 
    ##              1.193              1.193              1.193              1.193 
    ##           TRAPPC12              SH2D7            FAM163A              RUFY4 
    ##              1.193              1.193              1.193              1.193 
    ##             SLAMF8             ANKRD9            C3orf30            CCDC171 
    ##              1.193              1.194              1.194              1.194 
    ##               CD68             CALHM3              GCOM1             STEAP1 
    ##              1.194              1.195              1.195              1.195 
    ##             SEC22A            ST3GAL6             KCNAB3              KRT86 
    ##              1.195              1.195              1.195              1.195 
    ##              MMRN1              TCF19             TSPEAR              ALG14 
    ##              1.195              1.195              1.195              1.196 
    ##             SLC5A8              ARL15              RASA2              RASA4 
    ##              1.196              1.196              1.196              1.196 
    ##             TCP10L              IL36B          C17orf112           C17orf51 
    ##              1.196              1.196              1.196              1.196 
    ##           C10orf82             MRPL22               CA5B            SLC39A2 
    ##              1.196              1.196              1.196              1.196 
    ##               DSC3             PTCHD3               PMCH             ZNF575 
    ##              1.196              1.196              1.196              1.196 
    ##              CBWD1              DHRSX              ARL4D           SLC25A51 
    ##              1.197              1.197              1.197              1.197 
    ##            NDUFAF6              DYDC1              TMCO3              EPPK1 
    ##              1.197              1.197              1.197              1.198 
    ##             ALKBH2            SLC12A3           C10orf11               GZMM 
    ##              1.198              1.198              1.198              1.198 
    ##               PEX2              FITM1             RNF168             ARRDC1 
    ##              1.198              1.198              1.198              1.198 
    ##               NCF4             ZNF225             SHISA4            ZNF385C 
    ##              1.198              1.198              1.198              1.199 
    ##              FXYD2           C19orf54              MCHR2             FAHD2A 
    ##              1.199              1.199              1.199              1.199 
    ##              CHPF2             CXCL14              ENKUR              MMP21 
    ##              1.199              1.200              1.200              1.200 
    ##              NUBP2             OR10P1               SOX7     CTD-2192J16.24 
    ##              1.200              1.200              1.200              1.200 
    ##              GNG11            SLC39A5             CYP2A6            NDUFAF5 
    ##              1.200              1.200              1.201              1.201 
    ##              BSPRY               CTSF            KLHDC8A             HMGCS2 
    ##              1.201              1.201              1.201              1.201 
    ##             ZNF37A              TNNI1        AC006946.15               PTX3 
    ##              1.201              1.201              1.201              1.202 
    ##           C1orf141             IMMP2L             HEATR4              AMPD1 
    ##              1.202              1.202              1.202              1.202 
    ##               SCOC             LRTOMT             PCDHB8              SYPL2 
    ##              1.202              1.202              1.202              1.202 
    ##               BRF2              BNIPL             MRPL12              RPS29 
    ##              1.203              1.203              1.203              1.203 
    ##               LIPM               MOBP              TCP10              ALG10 
    ##              1.203              1.203              1.203              1.203 
    ##              PLEK2              HABP2            COL20A1            SULT1C4 
    ##              1.203              1.204              1.204              1.204 
    ##           NAALADL1               LDHB            PCDHA10             MYL12B 
    ##              1.204              1.204              1.204              1.204 
    ##             NUDCD2                MAK                RPE              THEGL 
    ##              1.205              1.205              1.205              1.205 
    ##              EFHC1            ANKRD39               GYPC            MYBBP1A 
    ##              1.205              1.205              1.205              1.206 
    ##             RASSF1               NFYB              TMUB1            CEACAM8 
    ##              1.206              1.206              1.206              1.206 
    ##             CLEC9A              P2RX1               TGM4               FUT7 
    ##              1.206              1.206              1.206              1.206 
    ##              HOXD1               LSM1               MMP9            TBC1D17 
    ##              1.206              1.206              1.206              1.206 
    ##             OR51I2               HYPK            C8orf47              MPPE1 
    ##              1.206              1.207              1.207              1.207 
    ##               NME1            SLC27A5          HIST2H2BE             ZSWIM7 
    ##              1.207              1.207              1.207              1.207 
    ##            ZSCAN22              LAIR1             RAB39A               A1BG 
    ##              1.208              1.208              1.208              1.208 
    ##            PCDHGA6             DCDC2B           C10orf90             VSTM2B 
    ##              1.208              1.208              1.208              1.208 
    ##             MFSD6L               LIPH               NTN3               POP4 
    ##              1.209              1.209              1.209              1.209 
    ##             ZNF808               PID1                TRH             KCNJ11 
    ##              1.209              1.209              1.209              1.210 
    ##              KCTD4             RECQL4              TIGD6             MOSPD3 
    ##              1.210              1.210              1.210              1.210 
    ##               PNKP             GLB1L3             ZNF852               GRK7 
    ##              1.210              1.210              1.210              1.211 
    ##              ITIH5                SRL             IFNAR1             NDUFV1 
    ##              1.211              1.211              1.211              1.211 
    ##             OR52N1             MFSD2B              CCL21                MX1 
    ##              1.211              1.211              1.212              1.212 
    ##            FAM180B            SUPT4H1            OR51H1P              MYOZ1 
    ##              1.212              1.212              1.212              1.212 
    ##              RADIL               LCN8              KLRD1              MOB3C 
    ##              1.212              1.213              1.213              1.213 
    ##               NOB1            TMEM231                PBK             EXOSC3 
    ##              1.213              1.213              1.213              1.213 
    ##             RRNAD1             TCEAL8              ARMC7             ANKMY1 
    ##              1.213              1.213              1.214              1.214 
    ##               MXD3             RASL12               ODF4               LIPJ 
    ##              1.214              1.214              1.214              1.214 
    ##             BTN2A2             GUCA2A             CMKLR1              PAAF1 
    ##              1.214              1.214              1.214              1.214 
    ##              DEGS2           SERPINF1              GNPTG          C17orf105 
    ##              1.214              1.214              1.215              1.215 
    ##               POLQ               CDT1               DRGX               FGGY 
    ##              1.215              1.215              1.215              1.215 
    ##             TMEM54             FBXO16              IP6K3                WRB 
    ##              1.215              1.215              1.215              1.215 
    ##              PARM1              SAMD3             CCDC66               FSHR 
    ##              1.215              1.215              1.215              1.215 
    ##               RNH1              TDRD6               BEX1            GOLGA8N 
    ##              1.216              1.216              1.216              1.216 
    ##              FOXL1               NAGK             OR10Z1             CHMP2B 
    ##              1.216              1.216              1.216              1.216 
    ##              AIMP1             ZNF829               GRAP             PLSCR2 
    ##              1.216              1.216              1.217              1.217 
    ##             HS1BP3            CEACAM6             PDLIM3               DPH7 
    ##              1.217              1.217              1.217              1.217 
    ##            PPFIBP2               IRG1             STRA13               GNG3 
    ##              1.217              1.217              1.217              1.218 
    ##       RP11-433C9.2      RP11-477N12.3              DGAT1            TMEM244 
    ##              1.218              1.218              1.218              1.218 
    ##            PLA2G2C               ASB2              C1QL2             RAD51D 
    ##              1.218              1.218              1.218              1.218 
    ##             LRRC23             CCDC81             MAP6D1               CD3G 
    ##              1.218              1.219              1.219              1.219 
    ##                EFS               MLF1              ZNF71              PGAM5 
    ##              1.219              1.219              1.219              1.219 
    ##              TMPPE             LRRIQ3                LTK      RP11-162P23.2 
    ##              1.220              1.220              1.220              1.220 
    ##              SPRY1             ZNF480             UGT2B4             ZNF692 
    ##              1.220              1.220              1.220              1.220 
    ##               SKA1              BCAT2            PCOLCE2              HEMGN 
    ##              1.220              1.220              1.221              1.221 
    ##             MRPL53             ADRA2B               LDHD             CD300A 
    ##              1.221              1.221              1.221              1.221 
    ##             KLHL30             ZNF789              UPK3B              KRT39 
    ##              1.221              1.221              1.221              1.221 
    ##      RP11-1012A1.4            C6orf15              HHLA1            ZDHHC13 
    ##              1.221              1.222              1.222              1.222 
    ##             NUTM2F              VEPH1             NANOS1              PANX2 
    ##              1.222              1.222              1.222              1.222 
    ##           KIAA0391              AKAP7               GAB4           SLC25A39 
    ##              1.222              1.222              1.222              1.222 
    ##              SCFD2              LUZP2              ODF2L              ACTA1 
    ##              1.222              1.222              1.222              1.223 
    ##              OVCH1           CHRFAM7A           C12orf40               MFRP 
    ##              1.223              1.223              1.223              1.223 
    ##              TFPI2              CEP89              NFAM1              SPAG6 
    ##              1.223              1.223              1.223              1.223 
    ##               SCG5            EFCAB11                AK2       CTC-454I21.3 
    ##              1.223              1.223              1.223              1.223 
    ##               GKN2             TMEM79           C12orf45           CATSPER2 
    ##              1.223              1.223              1.224              1.224 
    ##             CLDN10             PIWIL3            SLC13A1             SEC22C 
    ##              1.224              1.224              1.224              1.224 
    ##             GTPBP3             PRSS38             KCNA10             TAPBPL 
    ##              1.224              1.224              1.224              1.224 
    ##           PLEKHG4B              NEK11             GPR115              OBFC1 
    ##              1.224              1.224              1.224              1.224 
    ##              HAND1             DUSP26               GBP1            C1QTNF2 
    ##              1.225              1.225              1.225              1.225 
    ##            SLC12A8              MCUR1             SLC7A4             FAXDC2 
    ##              1.225              1.225              1.225              1.225 
    ##               SBDS               GJC2              SGOL1               CLN5 
    ##              1.225              1.225              1.225              1.226 
    ##                EPX             LEPRE1               MYF5               PLB1 
    ##              1.226              1.226              1.226              1.226 
    ##             ZNF154              NAA11              SYTL3             STARD5 
    ##              1.226              1.226              1.226              1.226 
    ##              ZNF19            CCDC153            PTPRCAP              TTLL3 
    ##              1.226              1.226              1.226              1.226 
    ##              WDR49             GPR110              AANAT             TTC23L 
    ##              1.227              1.227              1.227              1.227 
    ##            SLC27A1             ANKRD7              HOXA5               IQCC 
    ##              1.227              1.227              1.227              1.227 
    ##             SAMSN1           SIGLEC12       RP11-723O4.6                DPT 
    ##              1.227              1.227              1.228              1.228 
    ##              MFAP2             POLR3K             ZNF275             CHCHD7 
    ##              1.228              1.228              1.228              1.228 
    ##               FTCD               KXD1            TNFAIP6            ANKRD42 
    ##              1.228              1.228              1.228              1.229 
    ##           ATP6V1E2              ARL10              KRT38           C15orf62 
    ##              1.229              1.229              1.229              1.229 
    ##             MRPS36               EME1               APOF              MEDAG 
    ##              1.230              1.230              1.230              1.230 
    ##             LANCL1           ARL14EPL              ITIH1            SIGLEC1 
    ##              1.230              1.230              1.230              1.230 
    ##           XRCC6BP1           KIAA1407             ADRA2A               GJE1 
    ##              1.230              1.230              1.231              1.231 
    ##              GRSF1             AKNAD1               LCTL               NAAA 
    ##              1.231              1.231              1.231              1.232 
    ##             LRRC29               NFU1              CORIN             ZNF260 
    ##              1.232              1.232              1.232              1.232 
    ##              SPEF1             ENTPD8                MUT               KRT3 
    ##              1.232              1.232              1.232              1.232 
    ##              TCTN1              CHST5               ECI2           GOLGA6L2 
    ##              1.233              1.233              1.233              1.233 
    ##           PNLIPRP1            KBTBD12              GPR75          TNFAIP8L1 
    ##              1.233              1.233              1.233              1.233 
    ##             OLFML1              HTRA3               CD1C             HSPA1L 
    ##              1.233              1.233              1.233              1.233 
    ##           C9orf156           ATP6V1G3              AP1G2               PTMA 
    ##              1.234              1.234              1.234              1.234 
    ##              ITGB6              ORAI3               CD93            C9orf91 
    ##              1.234              1.234              1.235              1.235 
    ##             CACFD1            ZNF724P               CMC1              LECT1 
    ##              1.235              1.235              1.235              1.235 
    ##           PPP1R15A               NAT9            FAM200A       RP11-542P2.1 
    ##              1.235              1.235              1.235              1.235 
    ##             TXNDC5             ANKZF1              SIRT5              VAMP8 
    ##              1.235              1.235              1.235              1.235 
    ##              DERL3               TYW5               SGK2             PHLDB3 
    ##              1.236              1.236              1.236              1.236 
    ##          C21orf140              CLDN7              MYZAP              ACSM4 
    ##              1.236              1.236              1.236              1.236 
    ##             LRRC43           MARVELD3             ABCC11           SLC25A26 
    ##              1.236              1.236              1.236              1.236 
    ##            CCDC110              TCL1B              BCS1L              PLA1A 
    ##              1.236              1.236              1.237              1.237 
    ##              CWH43               CHIA             PCNXL4              HESX1 
    ##              1.237              1.237              1.237              1.237 
    ##                DDO        C1QTNF1-AS1            C9orf64              CDH26 
    ##              1.237              1.237              1.237              1.237 
    ##               GCAT               CPA6             CYP2B6           SLC25A13 
    ##              1.237              1.237              1.237              1.237 
    ##             MS4A18             SLC8B1             TEX101         PAXIP1-AS2 
    ##              1.238              1.238              1.238              1.238 
    ##           ATP6V1G2               PIGV              GPR83                C1D 
    ##              1.238              1.238              1.238              1.238 
    ##           C17orf96             IGSF10             FAM43B             CLCNKB 
    ##              1.238              1.238              1.238              1.239 
    ##              AP5B1             DUSP19              APBA3             NUDT10 
    ##              1.239              1.239              1.239              1.239 
    ##               MLKL               EFHB               FIG4             GPR89B 
    ##              1.239              1.240              1.240              1.240 
    ##            CCDC116               KRAS              ZMAT5               GBP4 
    ##              1.240              1.240              1.240              1.240 
    ##             ANXA10             HAPLN3            ALDH7A1             PRSS12 
    ##              1.240              1.240              1.240              1.240 
    ##              SPAM1            FAM105A                IGJ           TMEM200C 
    ##              1.240              1.241              1.241              1.241 
    ##               PDK4              ADH1B             ANKRD2              BAMBI 
    ##              1.241              1.241              1.241              1.241 
    ##              ARMC6            LDLRAD2              WFDC5             CYP2F1 
    ##              1.241              1.241              1.241              1.242 
    ##              PARP3              OR2B2               DCPS            FAM186B 
    ##              1.242              1.242              1.242              1.242 
    ##              CGNL1              WDR16               MBD4            SLC28A3 
    ##              1.242              1.242              1.242              1.242 
    ##             GUCA1B           C11orf87              CNBD2               NGRN 
    ##              1.242              1.242              1.243              1.243 
    ##            GPATCH4              SAMD9              PRR15              NLRC3 
    ##              1.243              1.243              1.243              1.243 
    ##            FAM129A             ZNF333             ANKS4B              GPR97 
    ##              1.243              1.244              1.244              1.244 
    ##            TPD52L1              PACRG             PYHIN1              PDZD9 
    ##              1.244              1.244              1.244              1.244 
    ##              IL17B              CCBL1             RASSF9           C22orf15 
    ##              1.244              1.245              1.245              1.245 
    ##             PCDHA3              MMP25             KLHL38              BNIP1 
    ##              1.245              1.245              1.245              1.245 
    ##              INHBC             LRRC56             TMEM88             CDC25C 
    ##              1.245              1.245              1.245              1.245 
    ##             BTBD16             SPRYD7               HES3              PRAME 
    ##              1.246              1.246              1.246              1.246 
    ##              KIF25             POLR3D             SPATC1            TMEM100 
    ##              1.246              1.246              1.246              1.246 
    ##              TRUB2            CCDC169            TMEM237              KNOP1 
    ##              1.246              1.246              1.246              1.246 
    ##         KRTAP10-12            TMEM130             GABRR1             UBE2V2 
    ##              1.246              1.246              1.246              1.246 
    ##              OR2J3              SETD4               CFC1               PGA3 
    ##              1.246              1.246              1.247              1.247 
    ##              DGAT2               EMG1             ZSWIM2               DKK4 
    ##              1.247              1.247              1.247              1.247 
    ##           HIST3H2A             DHTKD1            CCDC122             SUCNR1 
    ##              1.247              1.247              1.247              1.247 
    ##             LY6G6F           RNASEH2B             AMIGO3               PENK 
    ##              1.247              1.247              1.248              1.248 
    ##             KCNK15            SERINC2             PKD2L1               ATIC 
    ##              1.248              1.248              1.248              1.248 
    ##              RAB36            FAM181A           C17orf72            THUMPD3 
    ##              1.249              1.249              1.249              1.249 
    ##             ZNF580              MATN3               OAS2              MFSD1 
    ##              1.249              1.249              1.249              1.249 
    ##              PHKG1               IDO2           C1orf146              PXMP4 
    ##              1.250              1.250              1.250              1.250 
    ##               HPGD            SEC14L5             OGFOD1             APOPT1 
    ##              1.250              1.250              1.250              1.250 
    ##             UGT3A2             OLFML3             ZNF383              GSDMC 
    ##              1.250              1.250              1.250              1.250 
    ##            PCDHB16              RHNO1              KRT40            CLEC14A 
    ##              1.251              1.251              1.251              1.251 
    ##              CHAC1              SPNS3             BPIFB3               CTSV 
    ##              1.251              1.251              1.251              1.251 
    ##              GALK2               PGLS                BPI            SLC28A1 
    ##              1.251              1.252              1.252              1.252 
    ##               REG4               NT5C               STX8              SFTPC 
    ##              1.252              1.252              1.252              1.252 
    ##           ATP6V1G1              FSIP1              SPAG8             ZNF234 
    ##              1.252              1.252              1.252              1.252 
    ##               SUN5             MRPS35              ICAM2              ZNF16 
    ##              1.252              1.253              1.253              1.253 
    ##            METTL10           HEPACAM2              MGST3             ZNF141 
    ##              1.253              1.253              1.253              1.253 
    ##               EMR2               MYOC             IGSF23               UMOD 
    ##              1.253              1.253              1.253              1.253 
    ##             FBLIM1            ARL6IP5             ECHDC2             PCSK1N 
    ##              1.254              1.254              1.254              1.254 
    ##           C11orf96               RNF7               PIGM             TRIM50 
    ##              1.254              1.254              1.254              1.254 
    ##         C14orf166B              TEKT2               AOC2              ZNF14 
    ##              1.255              1.255              1.255              1.255 
    ##               CD82               GATS               ARV1              IRAK2 
    ##              1.255              1.255              1.255              1.255 
    ##               MELK            LEPREL4               UPB1               OPTN 
    ##              1.255              1.256              1.256              1.256 
    ##              CARD6               SDHC             SPATS1              TIGD4 
    ##              1.256              1.256              1.256              1.256 
    ##             MRPL51              KRT73              PYCR1             NUDT11 
    ##              1.256              1.257              1.257              1.257 
    ##               INIP             SOWAHD             HIGD1A             SRFBP1 
    ##              1.257              1.257              1.257              1.258 
    ##              HOXB6              ACSM5              GFRAL         AP001468.1 
    ##              1.258              1.258              1.258              1.258 
    ##               PEPD             LY6G6D           PRAMEF12             SUPT3H 
    ##              1.258              1.258              1.258              1.258 
    ##             CHRNA2             KLHL17              FFAR4              ZPBP2 
    ##              1.259              1.259              1.259              1.259 
    ##               TSHR              WDR73               CDK3              POTEC 
    ##              1.259              1.259              1.259              1.260 
    ##           SLC22A14            SLC6A20               GAMT              PRAP1 
    ##              1.260              1.260              1.260              1.260 
    ##             NDUFB5             TBC1D3              PLIN1      CTD-2287O16.3 
    ##              1.260              1.260              1.260              1.261 
    ##           C16orf92            RNF144B            TCP11L2              ITFG3 
    ##              1.261              1.261              1.261              1.261 
    ##              RNPEP              MTERF            C5orf45              FFAR2 
    ##              1.261              1.261              1.262              1.262 
    ##              CEP63             TMEM17               PON3             METTL5 
    ##              1.262              1.262              1.262              1.262 
    ##             CRYBA2                GEM              HTR3D             GDF5OS 
    ##              1.263              1.263              1.263              1.263 
    ##             HAVCR1              PIPOX              DACT2            SLC36A2 
    ##              1.263              1.263              1.263              1.263 
    ##            GPATCH3              NLRP8           ANKRD30A              MST1R 
    ##              1.263              1.263              1.263              1.263 
    ##            TRIM16L               EXD1              EPHX2              KRT77 
    ##              1.263              1.263              1.263              1.263 
    ##              FRG1B             GPR160            RIPPLY2              YARS2 
    ##              1.263              1.263              1.263              1.264 
    ##               PIGF              IDH3B              GPR25              PLBD1 
    ##              1.264              1.264              1.264              1.264 
    ##              SPEM1              ANXA3             ZFAND1              CLIC1 
    ##              1.264              1.264              1.265              1.265 
    ##              SCN1B             ZNF239              DNPH1               HAO2 
    ##              1.265              1.265              1.265              1.265 
    ##             ZNF474             TAS2R7               GNL3              PDZD3 
    ##              1.265              1.265              1.266              1.266 
    ##            PCDHGA2               RAC3           SERPINA7               PMS2 
    ##              1.266              1.266              1.266              1.266 
    ##             DGCR6L              MS4A1             MLLT11             CYP3A7 
    ##              1.266              1.266              1.266              1.266 
    ##             ZNF467              RMDN1              LPAR6              NMUR1 
    ##              1.266              1.267              1.267              1.267 
    ##              P2RX5               MURC              MEOX1               PRG2 
    ##              1.267              1.267              1.268              1.268 
    ##              DCDC1               DCTD              NXPE2            SCGB1D2 
    ##              1.268              1.268              1.268              1.268 
    ##              GREM2             SFTPA2            TIMMDC1              MUC15 
    ##              1.268              1.268              1.268              1.268 
    ##               OXTR             BPIFB6              BMP8B              FBXW8 
    ##              1.268              1.268              1.269              1.269 
    ##               DPH5             NDUFC1               GBP6            FAM150B 
    ##              1.269              1.269              1.269              1.269 
    ##               FSBP               LMF1            FAM221B            CCDC176 
    ##              1.269              1.269              1.269              1.269 
    ##              VAMP3            GTF2A1L           STAMBPL1              ACAD8 
    ##              1.269              1.269              1.269              1.270 
    ##               OAS3            SLITRK6              DUOX2            SULT1E1 
    ##              1.270              1.270              1.270              1.270 
    ##            PPP1R32              BATF3          KIAA0226L            SLCO1A2 
    ##              1.270              1.270              1.270              1.270 
    ##               FSD2               ADH4          TNFRSF13C             PNPLA4 
    ##              1.270              1.270              1.270              1.271 
    ##            FAM101B               ZFR2              FGF11              PXMP2 
    ##              1.271              1.271              1.271              1.271 
    ##             NUTM2G               BEX2             RHBDD3               UTS2 
    ##              1.271              1.271              1.271              1.271 
    ##               PLLP              ZNF44             YJEFN3         AC008271.1 
    ##              1.271              1.271              1.272              1.272 
    ##             CCDC84               RRP8             VWA5B1            PGLYRP2 
    ##              1.272              1.272              1.272              1.272 
    ##              KRT6B               CPA1            ANGPTL1              GLUD2 
    ##              1.272              1.273              1.273              1.273 
    ##              CDHR3               GLDN               OC90               OC90 
    ##              1.273              1.273              1.273              1.273 
    ##           FLJ27365              CLCA2              REEP5               CD55 
    ##              1.273              1.273              1.273              1.273 
    ##              PSAT1              INPP1              KIF4B             CLDN18 
    ##              1.273              1.273              1.273              1.273 
    ##             TRIM74              APEX1            TMEM206      RP11-383H13.1 
    ##              1.273              1.273              1.274              1.274 
    ##             ABCC12              MARC1             SERGEF            SMPDL3B 
    ##              1.274              1.274              1.274              1.274 
    ##              ZADH2              EFHD2               URM1              MCFD2 
    ##              1.274              1.275              1.275              1.275 
    ##             RPS27L               HTR6               PPA2       RP11-176H8.1 
    ##              1.275              1.275              1.275              1.275 
    ##            C3orf22               SORD             TMEM43                OMD 
    ##              1.275              1.275              1.276              1.276 
    ##             SH3TC1            TMEM160             DHRS12              MED31 
    ##              1.276              1.276              1.276              1.276 
    ##              KRBA2            MTERFD1              GMPR2              DPEP3 
    ##              1.276              1.276              1.276              1.276 
    ##             PTGES2             FAM63A           SLC25A18              GINS1 
    ##              1.276              1.276              1.276              1.277 
    ##            FAM228B               BCO2            FAM132B              BMP15 
    ##              1.277              1.277              1.277              1.277 
    ##            STEAP1B               FADD              CDHR5              MESP2 
    ##              1.277              1.277              1.277              1.277 
    ##           ARHGEF39              DDIT4             TSPYL4               GPX4 
    ##              1.278              1.278              1.278              1.278 
    ##             SLAMF9              PHF11              NTPCR               FEN1 
    ##              1.278              1.278              1.278              1.278 
    ##              PARVB             ZNF488             GOT1L1               SSX3 
    ##              1.278              1.278              1.278              1.278 
    ##            TMEM177              YIPF2              FCRL3              TTC32 
    ##              1.279              1.279              1.279              1.279 
    ##            GALNT15              MEP1B               MCM7     CTD-3105H18.16 
    ##              1.279              1.279              1.279              1.279 
    ##            EFCAB10                HN1              SUMF2              MCHR1 
    ##              1.279              1.279              1.279              1.280 
    ##            C4orf29              DAPK2              CRHR2               GM2A 
    ##              1.280              1.280              1.280              1.280 
    ##            SLC16A8              KRCC1             SAMD11            C1QTNF3 
    ##              1.280              1.280              1.280              1.280 
    ##               LIM2              FGF23               FPGT              NAA40 
    ##              1.280              1.280              1.281              1.281 
    ##              CDH16             ENTHD2            TRABD2A       hsa-mir-1199 
    ##              1.281              1.281              1.281              1.281 
    ##              MYBPH            C2orf47              PRELP           PDCD1LG2 
    ##              1.281              1.281              1.281              1.282 
    ##             NIF3L1            EIF2S3L            C6orf57             CNKSR1 
    ##              1.282              1.282              1.282              1.282 
    ##             EXOSC6              PRDX5              OXLD1              P2RX4 
    ##              1.282              1.282              1.282              1.282 
    ##               RERG             ARRDC2             CHST12      RP11-455G16.1 
    ##              1.282              1.282              1.283              1.283 
    ##              ADPRH               NME5             CCDC69             TMEM8C 
    ##              1.283              1.283              1.283              1.283 
    ##              ADCK1            ANGPTL4              LITAF               IL18 
    ##              1.283              1.283              1.284              1.284 
    ##             ZNF132              POLE3              LDOC1              NUPL2 
    ##              1.284              1.284              1.284              1.284 
    ##             RNF125           C1orf177               COA1              KCNJ1 
    ##              1.284              1.284              1.285              1.285 
    ##              TECTB                 MB               ERMN            ANKDD1A 
    ##              1.285              1.285              1.285              1.285 
    ##           C11orf49               MKKS              GTDC1              TTLL9 
    ##              1.285              1.285              1.285              1.285 
    ##              TTC7A              ASB15             BPIFB2             EPHA10 
    ##              1.285              1.286              1.286              1.286 
    ##             UGT1A5          TMPRSS11A            SLC41A3             ZNF732 
    ##              1.286              1.286              1.286              1.286 
    ##           C17orf98              EGFL6              SERP1              OR5P3 
    ##              1.286              1.286              1.286              1.286 
    ##             TAS1R2               NSG1              TTC9C             CASP10 
    ##              1.286              1.286              1.286              1.287 
    ##            ATP13A5            PRPSAP1               RHOF             UNC5CL 
    ##              1.287              1.287              1.287              1.287 
    ##               PCK1              PEBP4              CBLN3               NARF 
    ##              1.287              1.287              1.287              1.287 
    ##               PMEL            TMEM192           SYNDIG1L              MACC1 
    ##              1.287              1.287              1.287              1.288 
    ##             ZNF442             MBLAC2             HAPLN2             SAMD10 
    ##              1.288              1.288              1.288              1.288 
    ##            FAM19A4               CES2             FCHSD1           ADAMDEC1 
    ##              1.288              1.288              1.288              1.288 
    ##             ASPRV1               PECR            B3GALTL             TVP23A 
    ##              1.288              1.288              1.289              1.289 
    ##             SRD5A1               EMC6           C9orf171              NSUN3 
    ##              1.289              1.289              1.289              1.289 
    ##             TRIML2               MYL6              TAF1C             BPIFB4 
    ##              1.289              1.290              1.290              1.290 
    ##             CASP14           C15orf65               ORM1               MMAA 
    ##              1.290              1.290              1.290              1.290 
    ##             RMND5B             OR51F1               PON1              MICU1 
    ##              1.290              1.290              1.290              1.291 
    ##               DYTN             HSPA1A            FAM200B             CDK11A 
    ##              1.291              1.291              1.291              1.291 
    ##            TSPAN11              RAB32               IQCJ             FAM32A 
    ##              1.291              1.291              1.291              1.291 
    ##             GPR153               MISP               CILP            C1QTNF9 
    ##              1.292              1.292              1.292              1.292 
    ##            PLA2G2A              CD207             CLEC4E             DUSP13 
    ##              1.292              1.292              1.292              1.292 
    ##               NRTN             AVPR1A               KLK7             MRPL35 
    ##              1.292              1.292              1.293              1.293 
    ##              CCZ1B             MTNR1B              CALCB           C20orf27 
    ##              1.293              1.293              1.293              1.293 
    ##             LRRC28               OPN4             MB21D1             ZNF793 
    ##              1.293              1.293              1.293              1.293 
    ##               MEPE             CCDC63           SELENBP1           CCDC102B 
    ##              1.294              1.294              1.294              1.294 
    ##              APMAP              FNDC8             ZNF573            C7orf72 
    ##              1.294              1.294              1.294              1.294 
    ##              SYNE4             LMAN2L              PSRC1              C4BPB 
    ##              1.294              1.295              1.295              1.295 
    ##              GPR62             GPR171              S1PR4             GPR119 
    ##              1.295              1.295              1.295              1.295 
    ##               ALPP             RNF187             ZNF430             LRRC69 
    ##              1.295              1.295              1.295              1.295 
    ##             LRRC72              CENPO              HTR1D              FRMD1 
    ##              1.295              1.296              1.296              1.296 
    ##              ERCC1              ZNF20             PCDHB3              KRT20 
    ##              1.296              1.296              1.296              1.296 
    ##               MAVS              CAPS2            ARHGDIB             DNMT3L 
    ##              1.296              1.297              1.297              1.297 
    ##            ANGPTL5             NDUFS6             DZIP1L              WDR90 
    ##              1.297              1.297              1.297              1.297 
    ##               ZAR1             RNF215             IL20RB              ZRSR1 
    ##              1.297              1.297              1.298              1.298 
    ##           KIAA0141              IL36A             OR52E2             HIATL2 
    ##              1.298              1.298              1.298              1.298 
    ##               IL19              NSUN5              CFHR4               OGG1 
    ##              1.298              1.298              1.299              1.299 
    ##               SSX4         AC068533.7              ZNF76            OR10AD1 
    ##              1.299              1.299              1.299              1.299 
    ##            ALDH1B1               PRNP              CHD1L           SLC22A12 
    ##              1.299              1.299              1.299              1.299 
    ##               RIT2             ZNF285              BTNL2            METTL24 
    ##              1.299              1.299              1.300              1.300 
    ##             WASH4P              TCF15               APTX                PVR 
    ##              1.300              1.300              1.300              1.301 
    ##              RBM48            ANGPTL3             CELA3A              PTGR1 
    ##              1.301              1.301              1.301              1.301 
    ##              HAUS1              MORN3              REG3A               PMP2 
    ##              1.301              1.301              1.301              1.301 
    ##            ZC2HC1C            FAM124B           HIST1H3C            TRIM64B 
    ##              1.301              1.302              1.302              1.302 
    ##           C1orf159               NOL7             GPRC6A              ENDOU 
    ##              1.302              1.302              1.302              1.303 
    ##               NME9               PLD4               XYLB            EIF4E1B 
    ##              1.303              1.303              1.303              1.303 
    ##               NAGA              PPIL6              FGF21             OR8B12 
    ##              1.303              1.303              1.303              1.304 
    ##                AFM           FAM103A1               REM2             OR56A4 
    ##              1.304              1.304              1.304              1.304 
    ##                BTD              OR2D3             DNAJB8              UBE2T 
    ##              1.304              1.304              1.304              1.304 
    ##             CLEC4A             ZSCAN9             TRIM15            C7orf49 
    ##              1.305              1.305              1.305              1.305 
    ##                HAL              TMCO4              ACAA2             MRPL18 
    ##              1.305              1.305              1.306              1.306 
    ##             IGFLR1             IGSF22              AS3MT               NAIP 
    ##              1.306              1.306              1.306              1.306 
    ##             CLECL1              BEST2           PRAMEF11               PRB3 
    ##              1.306              1.306              1.306              1.306 
    ##             SLURP1            PLAC8L1             DIABLO              KLRC1 
    ##              1.306              1.307              1.307              1.307 
    ##              XAGE3             CLEC4F             MS4A6E               CFTR 
    ##              1.307              1.307              1.307              1.307 
    ##            PPP1R3F             CCDC77            FAM212A              LRWD1 
    ##              1.307              1.307              1.307              1.307 
    ##              METRN               SHPK             LRRC3C              OPRM1 
    ##              1.307              1.307              1.307              1.307 
    ##            AKR1B10               FRG1            DNAJC19                DXO 
    ##              1.307              1.308              1.308              1.308 
    ##              PLIN2               RGL4              UCHL3              RECQL 
    ##              1.308              1.308              1.308              1.308 
    ##              OLIG2              CANT1               IL1A              SUSD3 
    ##              1.308              1.309              1.309              1.309 
    ##              CLYBL             ZNF581             YY1AP1              NLRP9 
    ##              1.309              1.309              1.309              1.309 
    ##              FOXD2               CAV1           SERPINB7             OGFOD3 
    ##              1.310              1.310              1.310              1.310 
    ##              CLUL1              CDHR1              DHRS4              KRT37 
    ##              1.310              1.311              1.311              1.311 
    ##             TRIM77               GDF3             ZNF674              CRYGC 
    ##              1.311              1.311              1.311              1.311 
    ##              WDR54             SLC1A7            SULT1B1              CSRP3 
    ##              1.311              1.311              1.311              1.311 
    ##              OPRK1             PCDHB1       RP11-65D24.2             PRSS37 
    ##              1.311              1.311              1.312              1.312 
    ##              TTC31           C17orf62              DPEP1              DDX25 
    ##              1.312              1.312              1.312              1.312 
    ##              AGTR2           LOH12CR2             VPS9D1            TAS2R40 
    ##              1.312              1.312              1.312              1.312 
    ##              FXYD3              ZWINT              AP4M1             MFSD10 
    ##              1.313              1.313              1.313              1.313 
    ##              TPRKB              TECRL               PI16            CCDC137 
    ##              1.313              1.313              1.313              1.313 
    ##               KLF8              OR6S1              RCAN3             NKX6-2 
    ##              1.313              1.313              1.313              1.313 
    ##             ZNF511              TEX30             NFE2L3            C6orf58 
    ##              1.313              1.313              1.314              1.314 
    ##             ADRA1D       RP11-944L7.5                ZP4              MPDU1 
    ##              1.314              1.314              1.314              1.314 
    ##             CARTPT            MTERFD3             MARCH8         AP001362.1 
    ##              1.314              1.314              1.314              1.314 
    ##            SLC35E3               INSC              PTCD2            CYP2C18 
    ##              1.315              1.315              1.315              1.315 
    ##              NHLH1            RANBP17              CMTM1             ALPPL2 
    ##              1.315              1.315              1.315              1.315 
    ##            APOA1BP              FIGLA              CHRNG              SAGE1 
    ##              1.315              1.316              1.316              1.316 
    ##              ACCSL             FAM89A        MINOS1-NBL1             MRPL46 
    ##              1.316              1.316              1.316              1.316 
    ##             TMEM18             CLEC4C              ASCC1                FEV 
    ##              1.317              1.317              1.317              1.317 
    ##             ARL13A              UBE2U            SLC2A11               CD84 
    ##              1.317              1.317              1.317              1.317 
    ##               EPYC                CDA              NAGPA               FKRP 
    ##              1.317              1.318              1.318              1.318 
    ##              ZFP37           C6orf222              DPEP2              KRT12 
    ##              1.318              1.318              1.318              1.319 
    ##               SELP             PRSS21            GPR137B                PGC 
    ##              1.319              1.319              1.319              1.319 
    ##               PGM1           SLC25A32               CPA3               KLK1 
    ##              1.319              1.319              1.319              1.319 
    ##              PRR26              JMJD4               CD33           TRAPPC6B 
    ##              1.319              1.319              1.319              1.320 
    ##            FOXRED1             TRIM42          KRTAP10-7               BOP1 
    ##              1.320              1.320              1.320              1.320 
    ##            SLC25A2               PDXP              PRDM9            BCDIN3D 
    ##              1.320              1.320              1.320              1.321 
    ##               LRAT            TCP10L2             TSPYL1            C4orf48 
    ##              1.321              1.321              1.321              1.321 
    ##                ADA             ZNF157              RAB41              KRT72 
    ##              1.321              1.321              1.321              1.322 
    ##            TNFRSF4            PCYOX1L               PKIG            SLC29A3 
    ##              1.322              1.322              1.322              1.322 
    ##           C15orf61             OTUD6B               GJB3              SNX24 
    ##              1.322              1.323              1.323              1.323 
    ##               CTSG            PPAPDC3              OR6C2             MRPL17 
    ##              1.323              1.323              1.323              1.323 
    ##             ZNF160             ALKBH3      JMJD7-PLA2G4B            PLA2G4B 
    ##              1.323              1.324              1.324              1.324 
    ##               YRDC      RP11-503N18.3             TXNDC8              SMAGP 
    ##              1.324              1.324              1.324              1.324 
    ##             FAM71B            CCDC152           KIAA1143              OR8B4 
    ##              1.324              1.324              1.324              1.324 
    ##            ANKRD29             ANGPT4              MTFP1              KDM4E 
    ##              1.325              1.325              1.325              1.325 
    ##              CCNJL            C9orf96             ZNF134             ZNF790 
    ##              1.325              1.325              1.325              1.325 
    ##              MYO1A               MNS1               PIF1             MAMDC4 
    ##              1.325              1.325              1.325              1.326 
    ##             SMTNL1           TMEM126A              TCEB2      CTD-2267D19.3 
    ##              1.326              1.326              1.326              1.326 
    ##               CA12               NRG4              ACOT2              GSTT1 
    ##              1.327              1.327              1.327              1.327 
    ##              NOC2L           ZMPSTE24             SPRYD4               IDH1 
    ##              1.327              1.327              1.327              1.327 
    ##                BSX             UNC119              F2RL1              COX15 
    ##              1.327              1.327              1.327              1.327 
    ##         AP000721.4             ZNF792               ASB8              L1TD1 
    ##              1.327              1.327              1.327              1.328 
    ##              ALG11              UCKL1            PCDHB12            TXNDC12 
    ##              1.328              1.328              1.328              1.328 
    ##             DHRS11              GP1BB              OOSP2            PRIMPOL 
    ##              1.328              1.328              1.328              1.328 
    ##           MAP1LC3A               COQ6               SHC2              PFDN5 
    ##              1.328              1.328              1.329              1.329 
    ##             NMRAL1              AICDA             NUDT12             ZNF571 
    ##              1.329              1.329              1.329              1.329 
    ##              BCAS2               BCHE               LDHC              GNGT1 
    ##              1.329              1.329              1.329              1.329 
    ##              LCN12              LYPD2                LBP              HOXB2 
    ##              1.330              1.330              1.330              1.330 
    ##             ZNF517             MRPL30               RNF5             PCDHA8 
    ##              1.330              1.330              1.330              1.330 
    ##                ZP1            METTL15               IRX6             LYSMD4 
    ##              1.330              1.331              1.331              1.331 
    ##             RWDD2B             B3GAT3           KRTAP6-3               CTSB 
    ##              1.331              1.331              1.331              1.331 
    ##             UGT1A3             UGT1A8             RHEBL1             HSD3B7 
    ##              1.331              1.331              1.332              1.332 
    ##               CTSO               NOTO             ACTL7A                C8B 
    ##              1.332              1.333              1.333              1.334 
    ##              THAP6            TMEM233         AC040160.1              HTR2B 
    ##              1.334              1.334              1.334              1.334 
    ##             MANSC1              CECR5               GGCT              IQCF6 
    ##              1.334              1.335              1.335              1.335 
    ##             R3HCC1            RHOXF2B              ZNF22               CCR8 
    ##              1.335              1.335              1.335              1.335 
    ##            SULT1C2             METTL8              SPC24              TIMD4 
    ##              1.335              1.335              1.335              1.336 
    ##           ITGB1BP2             LILRB4     C15orf38-AP3S2             CLEC3B 
    ##              1.336              1.336              1.336              1.336 
    ##             ENGASE             ANKRD1               GPX5         AC007952.5 
    ##              1.336              1.336              1.336              1.336 
    ##             UBE2L6              OARD1              ELOF1             PRR15L 
    ##              1.337              1.337              1.337              1.337 
    ##            CLEC18A              KRT23            HRASLS5            MAP3K15 
    ##              1.337              1.337              1.337              1.337 
    ##            KM-PA-2             P2RY14              APOA1              HTRA4 
    ##              1.337              1.338              1.338              1.338 
    ##         CSGALNACT1              LRRN4               CD1B              MAP10 
    ##              1.338              1.338              1.338              1.338 
    ##             ALOX15            FOXD4L1              ENKD1               NME8 
    ##              1.338              1.339              1.339              1.339 
    ##              EDA2R              HTR3C              LIN7A              EPM2A 
    ##              1.339              1.339              1.339              1.339 
    ##             PHLDA3              KRT17              KCNK6              BLCAP 
    ##              1.339              1.339              1.339              1.339 
    ##              CHIT1              STPG1              MARCO              CABS1 
    ##              1.339              1.339              1.340              1.340 
    ##            C4orf33              CCBL2             ZNF683               LAD1 
    ##              1.340              1.340              1.340              1.340 
    ##               ARSA               APOD           KIAA1161              ZFP57 
    ##              1.340              1.340              1.340              1.340 
    ##               FUT1             OR10J5              OR1A2               NEK3 
    ##              1.340              1.340              1.340              1.340 
    ##            GAL3ST4             OCIAD2         ST20-MTHFS              PCSK9 
    ##              1.340              1.340              1.340              1.341 
    ##               CPB2             ZNF273            DNAJC12               GBP3 
    ##              1.341              1.341              1.341              1.341 
    ##            MRPS18A             RAB40B             CYP1A2           KIAA1257 
    ##              1.341              1.341              1.341              1.341 
    ##             COX4I2           KRTAP5-1               PTK6              PYCR2 
    ##              1.341              1.342              1.342              1.342 
    ##              BIRC7               RTP1               NEU3              KPNA7 
    ##              1.342              1.343              1.343              1.343 
    ##               CAV2               CISH       RP11-178L8.4               UTF1 
    ##              1.343              1.343              1.343              1.343 
    ##             SPATA9           CYB561D2            ANKRD37             CCDC11 
    ##              1.343              1.344              1.344              1.344 
    ##               MPST             STEAP3               MCEE            ANKRD49 
    ##              1.344              1.344              1.344              1.344 
    ##               COQ3              EVI2A               HES5            TMEM139 
    ##              1.344              1.344              1.344              1.345 
    ##              LYVE1             CXCL16           SLC25A34           PPP1R14C 
    ##              1.345              1.345              1.345              1.345 
    ##              RDH11               IL26               RGS5              TIGD7 
    ##              1.345              1.345              1.345              1.345 
    ##           C15orf43              TAF1D             SLFN12               AGER 
    ##              1.345              1.346              1.346              1.346 
    ##             GPR182            SLC35G6             FCGR2A              FHAD1 
    ##              1.346              1.346              1.346              1.346 
    ##          SPATA31E1              AGTR1             HSPB11              HSPB6 
    ##              1.347              1.347              1.347              1.347 
    ##              AQP11             FAM98C               ORM2              UQCRH 
    ##              1.347              1.347              1.347              1.347 
    ##              GALR1              GOSR2               AAR2            PCDHB15 
    ##              1.347              1.347              1.348              1.348 
    ##             ZNF215              BMP8A           C9orf131             POMZP3 
    ##              1.348              1.348              1.348              1.348 
    ##               GJC3             ENOSF1              DHCR7             ZNF169 
    ##              1.348              1.348              1.349              1.349 
    ##              CENPQ              BTNL9             WDYHV1              LONP2 
    ##              1.349              1.349              1.349              1.349 
    ##               BPHL                IVD            CCDC170               THEG 
    ##              1.349              1.350              1.350              1.350 
    ##              GINM1              JAGN1              PDCL3       URGCP-MRPS24 
    ##              1.350              1.350              1.350              1.350 
    ##            FLYWCH1            HSD17B1               RGL3             NUDT22 
    ##              1.350              1.351              1.351              1.351 
    ##              YIPF7             ARL2BP           CCDC109B              CDK15 
    ##              1.351              1.351              1.351              1.351 
    ##              TUBB1              GSTK1             LRRC40             B3GAT2 
    ##              1.352              1.352              1.352              1.352 
    ##               PPIF            METTL2A            PCDHGA8             MS4A4A 
    ##              1.352              1.352              1.352              1.353 
    ##               H1FX               PAOX             RFPL3S               ITPA 
    ##              1.353              1.353              1.353              1.353 
    ##             NDUFS5              CEP72               KYNU              FNDC9 
    ##              1.353              1.353              1.353              1.353 
    ##            HLA-DOB             COMTD1            METTL12            SLC11A1 
    ##              1.353              1.353              1.353              1.354 
    ##              H2AFJ              IMPG1            SLCO1B1              MROH6 
    ##              1.354              1.354              1.354              1.354 
    ##             MRPL24             ZNF251               GNG2            TP53TG5 
    ##              1.354              1.354              1.354              1.354 
    ##            METTL22            SDCCAG3              DPPA2            GPR37L1 
    ##              1.354              1.354              1.355              1.355 
    ##              ENDOV               LSM5             NLRP12              LYNX1 
    ##              1.355              1.355              1.355              1.355 
    ##              ERCC8            TMEM144              C5AR2            CXorf36 
    ##              1.355              1.355              1.355              1.355 
    ##              UPK1B             TMEM91           CDK5RAP1              NABP1 
    ##              1.356              1.356              1.356              1.356 
    ##             TRMT44             R3HDML              ATG10              CLHC1 
    ##              1.356              1.356              1.356              1.356 
    ##              RAB38              BBS10              KLF17             GRAMD2 
    ##              1.356              1.357              1.357              1.357 
    ##           C1orf173              CHODL            FAM210B               KLC3 
    ##              1.357              1.358              1.358              1.358 
    ##           APOBEC3B               GGT2               MYL1               BST1 
    ##              1.358              1.358              1.358              1.358 
    ##             CKMT1B             COMMD7              BARD1              THAP3 
    ##              1.358              1.358              1.358              1.358 
    ##            SLC22A5            CCDC103               NME6             SLC51B 
    ##              1.358              1.359              1.359              1.359 
    ##             SLC5A9             WDSUB1              FBLL1              TRIM4 
    ##              1.359              1.359              1.359              1.359 
    ##              SLFN5             GTF2A2              IFI35           C12orf29 
    ##              1.359              1.359              1.359              1.360 
    ##            OPN1MW2             FAHD2B              PCBD2              OR1L8 
    ##              1.360              1.360              1.360              1.360 
    ##            TMEM213               CNFN              MOCS3               FDX1 
    ##              1.360              1.360              1.360              1.361 
    ##      RP11-701P16.2           SLC25A45             CRISP1               PRH1 
    ##              1.361              1.361              1.361              1.361 
    ##             CAPNS2              CCR10             TBXA2R              PVALB 
    ##              1.361              1.361              1.361              1.361 
    ##                LXN              TMCO1              TBX10             TMEM26 
    ##              1.361              1.362              1.362              1.362 
    ##              KCNE4              TEX35               KLK2              SHMT1 
    ##              1.362              1.362              1.362              1.362 
    ##             PRSS53               FN3K           C11orf53           METTL21A 
    ##              1.363              1.363              1.363              1.363 
    ##              RPP14           NDUFA4L2              POTEJ             KHDC1L 
    ##              1.363              1.363              1.363              1.363 
    ##              KLK12               NOD2              SUSD5            NKIRAS1 
    ##              1.363              1.363              1.363              1.363 
    ##           ERVFRD-1             ACSM2A              FATE1              GPR68 
    ##              1.363              1.364              1.364              1.364 
    ##             CYB5R1            PPP1R35             TUBA3E              ADTRP 
    ##              1.364              1.364              1.364              1.365 
    ##            C3orf62             ALG10B         AC117834.1             POLR2J 
    ##              1.365              1.365              1.365              1.365 
    ##              ABCG8               CCR3               LYG1            ANKRD53 
    ##              1.365              1.365              1.365              1.365 
    ##            SLC15A3           MAP1LC3B            CD200R1              FANCA 
    ##              1.366              1.366              1.366              1.366 
    ##            CYP20A1           C17orf74             PSMB10            EFCAB13 
    ##              1.366              1.366              1.366              1.366 
    ##               CRNN           ITPRIPL1           C20orf96               TGM6 
    ##              1.366              1.366              1.367              1.367 
    ##              RIPK3              ASIC5              ATP5I              MYLK4 
    ##              1.367              1.367              1.367              1.367 
    ##             ATP5SL              MMP20             FCER1A              FCER2 
    ##              1.367              1.367              1.368              1.368 
    ##              RAD9B      RP11-105C20.2              PILRB               PIGP 
    ##              1.368              1.368              1.368              1.368 
    ##              LCN10             RHBDL2             WNT10A             TRIML1 
    ##              1.368              1.369              1.369              1.369 
    ##             UBE2V1             ODF3L2            SLC34A3               BAK1 
    ##              1.369              1.369              1.369              1.369 
    ##               DAD1           TMEM184A             TDRD12              XRCC2 
    ##              1.369              1.370              1.370              1.370 
    ##             TSPAN1              FETUB              NTSR1            TMEM216 
    ##              1.370              1.370              1.370              1.370 
    ##            FASTKD1             PCDHB6               ANO7             OR4F15 
    ##              1.371              1.371              1.371              1.371 
    ##             ABRACL              OR5M8              CNTD2               PFN3 
    ##              1.371              1.371              1.371              1.371 
    ##              SAP25             NDUFA5              HDHD3               CCL4 
    ##              1.372              1.372              1.372              1.372 
    ##              PDCD2               DPH1                CA2           C16orf59 
    ##              1.372              1.372              1.373              1.373 
    ##            MICALL2            SLC17A1                CCS              RGPD1 
    ##              1.373              1.373              1.373              1.373 
    ##           SLC26A10             SPDYE1              IKBKG            SHCBP1L 
    ##              1.373              1.373              1.373              1.373 
    ##             ACTR3C              FEM1A               CSN2             KLHL34 
    ##              1.374              1.374              1.374              1.374 
    ##              SPRY4             PCDHB4            ADPRHL1               AMBP 
    ##              1.374              1.374              1.374              1.374 
    ##             ELOVL2            C9orf50              CASQ2               BDH2 
    ##              1.374              1.375              1.375              1.375 
    ##            CYP4F12             CYP4F2              TOMM5             CDKN1A 
    ##              1.375              1.375              1.375              1.375 
    ##            CCDC157             ZNF136              CD209               HAX1 
    ##              1.375              1.375              1.375              1.375 
    ##           SLC25A35           KIAA0020               CTSW              PPCDC 
    ##              1.375              1.376              1.376              1.376 
    ##           SLC22A13             SPINT4            TGFBR3L             MRPS25 
    ##              1.376              1.376              1.376              1.376 
    ##           C21orf33           KRTAP5-5            CYP27A1            CYP39A1 
    ##              1.376              1.376              1.376              1.376 
    ##             ACTRT2               EPN3               MYL4            TMEM251 
    ##              1.377              1.377              1.377              1.377 
    ##            METTL7A             GABRA6            BCL2L12              MSLNL 
    ##              1.378              1.378              1.378              1.378 
    ##            LURAP1L               GNB3              STX11             NUDT13 
    ##              1.378              1.378              1.379              1.379 
    ##               HAS1           SLC25A15           TMEM200B              ZNF99 
    ##              1.379              1.379              1.379              1.380 
    ##              LYPD5              MRPL3              MSS51               HES6 
    ##              1.380              1.380              1.380              1.380 
    ##              KIF2B               TET2             THNSL1              SURF2 
    ##              1.380              1.380              1.380              1.381 
    ##              RAB37              OR2B3              CPNE7             NKX3-1 
    ##              1.381              1.381              1.381              1.381 
    ##              OR6N2               TLR5             ORAOV1             CYP4F3 
    ##              1.381              1.382              1.382              1.382 
    ##             LRRC61               FAIM             MRPL15             SDCBP2 
    ##              1.382              1.382              1.382              1.382 
    ##              DEFA4             CASP16           CTNNBIP1               PRPH 
    ##              1.382              1.383              1.383              1.383 
    ##              OR1D5            SLC27A6               AIG1            RIPPLY3 
    ##              1.383              1.383              1.383              1.383 
    ##             PHYKPL              CDCA4               POMC             RAET1G 
    ##              1.383              1.383              1.384              1.384 
    ##           TP53INP1             UGT2A2              HSPA6              ABCG5 
    ##              1.384              1.384              1.384              1.385 
    ##            ANKRD16               BAAT                MSC               PRM2 
    ##              1.385              1.385              1.385              1.385 
    ##             ZNF230               VWA2              MOCS2               OCM2 
    ##              1.385              1.386              1.386              1.386 
    ##         ST6GALNAC4              RFESD            PPP1R3E             CYP3A5 
    ##              1.387              1.387              1.387              1.387 
    ##               ACP6           C16orf74               GNG7              TBPL2 
    ##              1.387              1.387              1.387              1.387 
    ##           SLC25A31              ATP5H                GH2             LRRC53 
    ##              1.387              1.387              1.387              1.387 
    ##             GALNT8              MILR1             OR2A12             CRYBA4 
    ##              1.388              1.388              1.388              1.388 
    ##               MRM1            METTL18              SNRPA             SPACA3 
    ##              1.388              1.388              1.388              1.388 
    ##               GALE            TMEM235            DNAJC28               NPC2 
    ##              1.388              1.388              1.389              1.389 
    ##            SLC17A9               MLIP              MYOZ3             RNF186 
    ##              1.389              1.389              1.389              1.389 
    ##              IKBIP             CITED1           C16orf58             CYP2J2 
    ##              1.389              1.389              1.390              1.390 
    ##               BDH1               MTAP             NECAB2             ASRGL1 
    ##              1.390              1.390              1.390              1.391 
    ##                AK3              SCRN2            TMEM30B              GSTM2 
    ##              1.391              1.391              1.391              1.392 
    ##             TBXAS1            FOXD4L3            FAM107B             TUBAL3 
    ##              1.392              1.392              1.392              1.392 
    ##            FAM187A             RABEPK              PTH2R             ADRA1A 
    ##              1.392              1.392              1.392              1.393 
    ##            CCDC129               ORC6            MARCH10               PCK2 
    ##              1.393              1.393              1.393              1.393 
    ##            SPINK14            TRAM1L1               ARSH              ATP5S 
    ##              1.394              1.394              1.394              1.394 
    ##              RGS13             CDKN2D           C1orf220              QPCTL 
    ##              1.394              1.394              1.394              1.394 
    ##            SLC2A10              KRT78            ABHD14A              SPA17 
    ##              1.394              1.395              1.395              1.395 
    ##           DNASE1L2              FDCSP               PMVK          KRTAP10-5 
    ##              1.395              1.395              1.395              1.395 
    ##               MND1              OR2V2            SLC35B3           C1orf109 
    ##              1.395              1.395              1.395              1.396 
    ##              KRT36              CRLF2            C21orf2           TOR1AIP2 
    ##              1.396              1.396              1.396              1.396 
    ##             LRRC27             PM20D1           KRTAP3-2           KRTAP4-8 
    ##              1.397              1.397              1.397              1.397 
    ##              MED28               HEYL             ZNF878           IRAK1BP1 
    ##              1.397              1.397              1.397              1.398 
    ##             ZNF426           C12orf39             MS4A15              ZNF17 
    ##              1.398              1.398              1.398              1.398 
    ##             GPANK1               CPA2            GLYATL3              SNX32 
    ##              1.398              1.399              1.399              1.399 
    ##             CTHRC1              RSPO1           TMEM106C             ZNF841 
    ##              1.399              1.399              1.399              1.399 
    ##            C5orf38              TEX26             LTB4R2              TFB1M 
    ##              1.399              1.399              1.400              1.400 
    ##             NDUFA4              ROMO1        AC024592.12            LYPLAL1 
    ##              1.400              1.400              1.400              1.400 
    ##              MYRFL              MS4A3                INA               CPA4 
    ##              1.401              1.401              1.401              1.401 
    ##              KCNH6             POLR1C               BHMT              SCN2B 
    ##              1.401              1.401              1.401              1.401 
    ##             ZNF524             TMEM42              HSBP1             CYP2C9 
    ##              1.401              1.402              1.402              1.402 
    ##              KRT26             ARID3C               ACPP               DPM3 
    ##              1.402              1.402              1.403              1.403 
    ##              P2RX2             PDE11A              OTOP2            CCDC74B 
    ##              1.403              1.404              1.404              1.404 
    ##            FAM101A         KRTAP10-11           MAD2L1BP           ATP6V0D2 
    ##              1.404              1.404              1.404              1.404 
    ##              MYL6B             UBXN2A             CCDC17             ACTL10 
    ##              1.405              1.405              1.405              1.405 
    ##               M1AP             ZNF345               GYG1           C1orf174 
    ##              1.405              1.405              1.405              1.405 
    ##              ATG4C              GABRP             CCDC87             ZNF181 
    ##              1.405              1.405              1.405              1.405 
    ##           C6orf201             SERTM1            C4orf51               PFN4 
    ##              1.406              1.406              1.406              1.406 
    ##             P2RY12            SCGB2A1               MALL               DBX2 
    ##              1.406              1.406              1.406              1.406 
    ##           SLC26A11          KRTAP10-3              MUC21            CXorf58 
    ##              1.407              1.407              1.407              1.407 
    ##              TPPP2              SPHK1             NUDT15              KRT83 
    ##              1.407              1.408              1.408              1.408 
    ##             FAM83A               WBP1              SCIMP               BST2 
    ##              1.408              1.408              1.408              1.409 
    ##              CETN3             ERICH2              ASPDH               SCGN 
    ##              1.409              1.409              1.409              1.409 
    ##             TRIM51             TMEM82             ZNF525           C12orf43 
    ##              1.409              1.409              1.409              1.410 
    ##             CHMP4C              OR1J4            C9orf89              KISS1 
    ##              1.410              1.410              1.410              1.410 
    ##            TMEM14A            IER3IP1              IFT27              PTTG1 
    ##              1.410              1.410              1.411              1.411 
    ##             TSPAN4            SLC29A2             ZNF391             CCDC68 
    ##              1.411              1.411              1.411              1.411 
    ##           KRTAP1-4             SPACA7            ITGB3BP               GPX7 
    ##              1.412              1.412              1.412              1.412 
    ##             RPH3AL            ZDHHC24            C3orf33              TRAP1 
    ##              1.412              1.412              1.412              1.412 
    ##                VIP              VSTM1              MS4A8              FDX1L 
    ##              1.412              1.412              1.413              1.413 
    ##             NIPAL2           SERPINB8              SPSB4             TMEM92 
    ##              1.413              1.413              1.413              1.413 
    ##             STYXL1             LRRC15              MMP10               ARSI 
    ##              1.414              1.414              1.414              1.414 
    ##               ASB4                FAH              DDX28             MRPL16 
    ##              1.414              1.414              1.414              1.414 
    ##              SERF2               AVEN             BTN3A2          C14orf119 
    ##              1.414              1.414              1.415              1.415 
    ##              MYADM         GADD45GIP1       RP11-468E2.6                PRL 
    ##              1.415              1.415              1.415              1.415 
    ##              OR2J1           C18orf42             IMMP1L             RNMTL1 
    ##              1.415              1.415              1.415              1.415 
    ##            DNAJC30             KCNK16             GUCA1A            C3orf84 
    ##              1.416              1.416              1.416              1.416 
    ##              WFDC8             INO80C                GP2               MMP3 
    ##              1.416              1.416              1.416              1.416 
    ##              OR3A1               PEX7              PRDM7             SH2D4A 
    ##              1.416              1.416              1.416              1.416 
    ##              ETV3L            FAM209A       RP11-96O20.4              ZNRF4 
    ##              1.416              1.417              1.417              1.417 
    ##             GLT1D1            CEACAM4               RGS2              MEIOB 
    ##              1.417              1.417              1.417              1.417 
    ##             RANGRF              OR2G2                PDC              ISOC1 
    ##              1.417              1.418              1.418              1.418 
    ##             SAMD13       RP11-89K11.1               ETV2              OR2G3 
    ##              1.418              1.418              1.418              1.419 
    ##           MPHOSPH6              TINAG              GJA10               ICT1 
    ##              1.419              1.419              1.419              1.419 
    ##            TMEM225             AKR1E2             MAGEH1         AP002884.3 
    ##              1.419              1.419              1.419              1.420 
    ##             ZNF266            ANAPC16           TMPRSS12             ACSM2B 
    ##              1.420              1.420              1.420              1.420 
    ##             ZNF514             SPINK7             ZNF662             MFSD12 
    ##              1.421              1.421              1.421              1.421 
    ##              LRIT3             POLR2G            SLC6A19              ROGDI 
    ##              1.421              1.421              1.422              1.422 
    ##              ROPN1            TRAPPC5            C8orf37               POLN 
    ##              1.422              1.422              1.423              1.423 
    ##           FANCD2OS            C2orf40              KRT76               GCSH 
    ##              1.423              1.423              1.423              1.423 
    ##               MMD2             KCTD18            ZDHHC23              KTI12 
    ##              1.423              1.423              1.424              1.424 
    ##              FOLR2              SSTR3               CD8B             STMND1 
    ##              1.424              1.424              1.425              1.425 
    ##              ATP5J               ART3            TMEM143            LDLRAD1 
    ##              1.425              1.425              1.425              1.425 
    ##             LHFPL5                GK2             MRPS15        FPGT-TNNI3K 
    ##              1.425              1.425              1.425              1.425 
    ##            B3GALT6            C4orf17              PYCRL               ACY3 
    ##              1.425              1.426              1.426              1.426 
    ##           C1orf123           LGALS3BP            TMEM14C              CYB5A 
    ##              1.426              1.426              1.426              1.427 
    ##              FCRLA             NAPRT1              AP3S2               GUK1 
    ##              1.427              1.427              1.427              1.427 
    ##             GLTPD1              ACRV1             DEPTOR            CCDC113 
    ##              1.427              1.427              1.427              1.427 
    ##         CTB-96E2.2              GADL1             TNNI3K             NDUFB4 
    ##              1.427              1.427              1.427              1.428 
    ##               FGF6     CTD-2521M24.10              DKKL1             OR51V1 
    ##              1.428              1.428              1.428              1.428 
    ##             LGALS3               LAYN             TVP23C              PAMR1 
    ##              1.428              1.428              1.428              1.429 
    ##         AC104794.4              CRYAB            CHRNA10              PPM1N 
    ##              1.429              1.429              1.429              1.429 
    ##             LILRB5            ZSCAN30             ADAM21            CLEC12A 
    ##              1.429              1.429              1.430              1.430 
    ##            ANKRD36             FAM92B               NOX5              TBATA 
    ##              1.430              1.430              1.430              1.430 
    ##             OR13A1              CDKL2             RNF175              BEST3 
    ##              1.430              1.431              1.431              1.431 
    ##            FAM227B              RHOT2              RBM34              OR1E2 
    ##              1.431              1.431              1.431              1.431 
    ##               CDX4              STAP2              GSDMA               OAZ3 
    ##              1.432              1.432              1.432              1.432 
    ##             TCEAL7              TAAR1      RP11-487E13.1              HYLS1 
    ##              1.432              1.432              1.433              1.433 
    ##               HPDL               EXO5            ZSCAN21             CCDC27 
    ##              1.433              1.433              1.433              1.433 
    ##              ANXA9              GCNT3             FAM83F            SERINC4 
    ##              1.433              1.434              1.434              1.434 
    ##               SRA1             ZNHIT2            FAM124A              SFRP5 
    ##              1.434              1.434              1.434              1.434 
    ##              INHBE              OR6C6         AC007557.1              CCL16 
    ##              1.434              1.434              1.434              1.435 
    ##              RDH12               REM1               LRR1          C10orf128 
    ##              1.435              1.435              1.435              1.435 
    ##           RPA3-AS1             DHFRL1         AC022431.2              OR2C3 
    ##              1.435              1.435              1.435              1.436 
    ##             MRPL54               PNPO               TGM5             ANKLE1 
    ##              1.436              1.436              1.436              1.437 
    ##              OR6Q1              LYRM7             SNAPIN              AGAP5 
    ##              1.437              1.437              1.437              1.438 
    ##               NCR2             CCDC61              OR6B2             RNF135 
    ##              1.438              1.438              1.438              1.438 
    ##             ATPAF2             UGT2B7               HRH4                SAG 
    ##              1.438              1.438              1.438              1.438 
    ##             PRSS33              CMPK2               PMF1         PMF1-BGLAP 
    ##              1.439              1.439              1.439              1.439 
    ##             CLEC6A               RGS4            S100A12           CYP4F31P 
    ##              1.439              1.440              1.440              1.440 
    ##              OR4M1           C21orf88             GPR152            CD300LF 
    ##              1.440              1.441              1.441              1.441 
    ##               HKR1             SHISA5              KLRB1             OR13D1 
    ##              1.441              1.441              1.441              1.441 
    ##           PLA2G12A              MRPL4              MS4A5               LBX2 
    ##              1.441              1.441              1.441              1.441 
    ##              MCPH1             CRABP1      RP11-650K20.3             HSD3B1 
    ##              1.441              1.441              1.441              1.442 
    ##      RP11-831H9.11          C20orf195               NRGN      RP11-201K10.3 
    ##              1.442              1.442              1.442              1.442 
    ##         AL078585.1           C19orf60           C1orf137               KRT7 
    ##              1.442              1.443              1.443              1.443 
    ##         AC018755.1              PSMG2             CRABP2             OR10A5 
    ##              1.443              1.443              1.443              1.444 
    ##             OSGIN1           SLC16A11              OR6T1          C20orf173 
    ##              1.444              1.444              1.444              1.444 
    ##               ERAS              CISD1               AMZ2                F11 
    ##              1.445              1.445              1.445              1.445 
    ##                AGT            C6orf25               SUN3             TMEM44 
    ##              1.445              1.445              1.445              1.445 
    ##                RRH            ZSCAN23              KRT24                KHK 
    ##              1.445              1.445              1.445              1.445 
    ##                C4A             OPALIN             CLEC1B               NENF 
    ##              1.446              1.446              1.446              1.446 
    ##             SCAND1              GSDMD            FAM230A              ANXA5 
    ##              1.446              1.446              1.446              1.446 
    ##              PSMD9             TMEM71             NLRP10             CPPED1 
    ##              1.446              1.446              1.447              1.447 
    ##              SNX21            C1QTNF7             CHI3L2             ZNF773 
    ##              1.447              1.447              1.447              1.447 
    ##             LILRA4               FMO5              QRFPR             TRIM68 
    ##              1.448              1.448              1.448              1.448 
    ##              GPR87             MRPL32              NIM1K              TCTE1 
    ##              1.448              1.448              1.448              1.449 
    ##              LECT2               TPTE           METTL21C               NCR3 
    ##              1.449              1.449              1.449              1.449 
    ##           GLIPR1L1              CTRB2             SFTPA1            FAM210A 
    ##              1.450              1.450              1.450              1.450 
    ##              NUDT2           GLIPR1L2           C19orf77            SLC35A4 
    ##              1.450              1.450              1.450              1.450 
    ##            C4orf50              MYCT1             FAM84B            KRTCAP3 
    ##              1.451              1.451              1.451              1.451 
    ##              RPP14         AC004076.9             NECAP2             RASA4B 
    ##              1.451              1.451              1.451              1.452 
    ##           SLC25A21                AVP            FAM170B              MTFR1 
    ##              1.452              1.452              1.452              1.452 
    ##             PARPBP              CDK10            ZMYM6NB              MOXD1 
    ##              1.453              1.453              1.453              1.453 
    ##               TSFM             PLCXD2            CAMK2N2               FHIT 
    ##              1.453              1.453              1.453              1.453 
    ##            TMEM116             ARRDC5               EAF2               COA6 
    ##              1.454              1.454              1.454              1.454 
    ##               NCR1             PAIP2B            C1orf85            MRPS18C 
    ##              1.454              1.454              1.454              1.454 
    ##            SPATA25           SLC22A24              KLRC4             ZNF669 
    ##              1.454              1.455              1.455              1.455 
    ##       CTD-3088G3.8              GLRA4               CES1             ZFP69B 
    ##              1.455              1.455              1.456              1.456 
    ##            SLC38A6              ETHE1         AC099552.4                IL5 
    ##              1.456              1.456              1.456              1.456 
    ##             S100A9               CAPG            TBC1D3B              ATP5D 
    ##              1.457              1.457              1.457              1.457 
    ##             ROPN1B             NT5C3B               TOE1              HLA-G 
    ##              1.457              1.457              1.457              1.458 
    ##             GLT8D2               IBSP               CD3D            ZSCAN16 
    ##              1.458              1.458              1.458              1.459 
    ##               KLK4              ADAT3             PROCA1             TMIGD1 
    ##              1.459              1.459              1.459              1.459 
    ##               CA14               FTMT            TMEM199               NREP 
    ##              1.459              1.459              1.460              1.460 
    ##            LDHAL6B           C9orf135             RIMBP3              TEKT3 
    ##              1.460              1.460              1.460              1.460 
    ##             OXNAD1             AKR1C2             ZNF320              HAGHL 
    ##              1.460              1.460              1.461              1.461 
    ##              CYLC2             CCDC24             PIK3R6          NME1-NME2 
    ##              1.461              1.461              1.461              1.462 
    ##               NME2            EP400NL            ALDH8A1              WISP2 
    ##              1.462              1.462              1.462              1.462 
    ##                SRI            LGALS7B              BNIP3            KRTCAP2 
    ##              1.462              1.462              1.462              1.462 
    ##           ALS2CR12             CYB5D2             FAM96A              GPR17 
    ##              1.462              1.463              1.463              1.463 
    ##           HSD17B11               KLK3             UGT1A7            PLEKHB2 
    ##              1.463              1.463              1.463              1.463 
    ##              VWA5A              TIMM9              EFNA4             CYP4Z1 
    ##              1.463              1.464              1.464              1.464 
    ##               GLO1            KAZALD1                GP9                GML 
    ##              1.464              1.464              1.464              1.465 
    ##               CHP2            ST7-OT4             TMEM59             ZNF763 
    ##              1.465              1.465              1.465              1.465 
    ##            SLC6A18          LINC01100              KRT79           RFPL4AL1 
    ##              1.465              1.466              1.466              1.466 
    ##               SGCB              NBPF4            FAM206A              LCE3C 
    ##              1.466              1.466              1.466              1.466 
    ##              MMP19               LCN6              C1QBP              PAGE5 
    ##              1.467              1.467              1.467              1.467 
    ##               MDP1           RNASEH2C             KCTD11             CCDC97 
    ##              1.467              1.467              1.467              1.468 
    ##       RP11-181C3.1             ZNF747               MAS1          SERPINB10 
    ##              1.468              1.468              1.468              1.468 
    ##              APOC1             MRPL28             ZNF415             TRIM63 
    ##              1.468              1.469              1.469              1.469 
    ##                DBI             TMEM74           B4GALNT2               MANF 
    ##              1.469              1.469              1.469              1.469 
    ##              AP5Z1            CREB3L3            TMEM219              KLF11 
    ##              1.469              1.469              1.470              1.470 
    ##             ASNSD1              MORN5             OR13C9              ACKR2 
    ##              1.470              1.470              1.470              1.470 
    ##               CAV3             SPATA6            APOBEC4              STPG2 
    ##              1.470              1.470              1.470              1.470 
    ##              OR5B2               XKR8               TAL2             ENTPD2 
    ##              1.471              1.471              1.471              1.471 
    ##            C3orf49            L3HYPDH             NUTM2B              OR4N5 
    ##              1.471              1.471              1.472              1.472 
    ##              RGPD2           C11orf42            C5orf46              GFRA4 
    ##              1.472              1.472              1.472              1.472 
    ##              MSGN1              ENPP5            SLC44A3             OR10H5 
    ##              1.472              1.473              1.473              1.473 
    ##            C2orf76              LCMT2             NDUFA3             PHYHD1 
    ##              1.473              1.473              1.473              1.474 
    ##              CKMT2             UBE2D4               FEZ2             GPRIN3 
    ##              1.474              1.474              1.474              1.474 
    ##      RP1-170O19.20             KCNK17               LST3              OR2T5 
    ##              1.474              1.474              1.475              1.475 
    ##              CETN1              PRRG2              OSCAR              AZGP1 
    ##              1.475              1.475              1.475              1.475 
    ##              RIBC2             OR10D3            TRAPPC4             OR5B12 
    ##              1.475              1.476              1.476              1.476 
    ##               MSLN              BBOX1             ATAD3B              PRRC1 
    ##              1.476              1.477              1.477              1.477 
    ##             MRPL43            FAM151B              HOXA7             LGALS2 
    ##              1.477              1.477              1.477              1.478 
    ##            TSPAN16             HSPA1B             ZNF343             PARP16 
    ##              1.478              1.478              1.478              1.478 
    ##                CPQ             ZNF419                TST              NSUN6 
    ##              1.479              1.479              1.479              1.479 
    ##        AC074091.13              FSCN2             MTRF1L              COX20 
    ##              1.479              1.479              1.480              1.480 
    ##           C16orf90            SLC27A3              HSPB2               NQO2 
    ##              1.480              1.480              1.480              1.480 
    ##              SOAT2             NCCRP1              UVSSA              TNNT1 
    ##              1.480              1.480              1.480              1.480 
    ##             SPDYE2               EBPL              DHRS2              OR1E1 
    ##              1.480              1.481              1.481              1.481 
    ##             MYEOV2              OR4E2              LYPD6           C1orf216 
    ##              1.481              1.481              1.481              1.481 
    ##       TVP23C-CDRT4               GPR1              WDR5B           TMEM106A 
    ##              1.481              1.481              1.481              1.481 
    ##             DCAF16               EMR3             SGK223              EVPLL 
    ##              1.481              1.482              1.482              1.482 
    ##             ZNF684            SLC26A1               ZIM3               COA4 
    ##              1.482              1.482              1.482              1.482 
    ##             PIH1D2             SEC61B             RASSF3            TMEM254 
    ##              1.482              1.483              1.483              1.483 
    ##               ADM2               ALLC              GSTZ1            SLC17A3 
    ##              1.483              1.483              1.483              1.483 
    ##          KIAA1024L              PRLHR              MYOM2               HBE1 
    ##              1.484              1.484              1.484              1.484 
    ##              KLRF1              DLEU1             ZNF774           FAM114A1 
    ##              1.484              1.484              1.484              1.485 
    ##            ADORA2B             PCYOX1             GPR89A             EXOSC1 
    ##              1.485              1.485              1.485              1.485 
    ##             ZNF664            TMEM222             MBLAC1            TAX1BP3 
    ##              1.485              1.485              1.485              1.485 
    ##              ACYP1              CNGA3             UGT3A1              OR4L1 
    ##              1.486              1.486              1.486              1.486 
    ##             PEX11A              PRMT3               ALPI             CALHM1 
    ##              1.486              1.486              1.486              1.486 
    ##                EMB             TPRG1L           C10orf99              NLRP2 
    ##              1.486              1.486              1.487              1.487 
    ##              CISD3             LILRA5             TMBIM4              MUC17 
    ##              1.487              1.487              1.487              1.487 
    ##             POLR2E              F2RL2             SPATA8               DCXR 
    ##              1.487              1.487              1.487              1.488 
    ##              DHRS7               MYL9           SERPINA6            PPP1R42 
    ##              1.488              1.488              1.488              1.488 
    ##              HDHD2           TCTEX1D2               GJA3             TEX13A 
    ##              1.488              1.488              1.488              1.488 
    ##              NMRK1              BEST1             ABHD15           SERPINB1 
    ##              1.489              1.489              1.489              1.489 
    ##               CD48             GPRIN2                MIA            SNRNP35 
    ##              1.489              1.489              1.489              1.489 
    ##              FXYD7             SPINK2          LINC00935            NUP62CL 
    ##              1.489              1.490              1.490              1.490 
    ##            C9orf40               COQ7               IL11             RAD51C 
    ##              1.491              1.491              1.491              1.491 
    ##              SRXN1             SLX4IP            SPATA19            SPATA32 
    ##              1.491              1.491              1.491              1.491 
    ##               MREG             RPL10L           C22orf31             CHRAC1 
    ##              1.491              1.491              1.492              1.492 
    ##              KCNE3              TIGD3              FXYD4               KERA 
    ##              1.492              1.492              1.492              1.492 
    ##              EFHD1               MZB1              POTEM            NDUFA12 
    ##              1.492              1.492              1.492              1.493 
    ##             ZNF490              OR8D2              SNUPN               ENO4 
    ##              1.493              1.493              1.493              1.494 
    ##              NAA38              CD302             RPL37A            FAM181B 
    ##              1.494              1.494              1.494              1.494 
    ##                NMB               LY96               CD8A             TMEM97 
    ##              1.494              1.494              1.495              1.495 
    ##             MIF4GD              OR2J2          C17orf107             DHRS7B 
    ##              1.495              1.495              1.495              1.495 
    ##            FAM151A               GALM             BPIFA2                HRG 
    ##              1.495              1.495              1.495              1.495 
    ##              REEP6             ITPRIP              IFI44              NAT16 
    ##              1.495              1.496              1.496              1.496 
    ##             DUSP15             FCGR3A            SLC14A1            C8orf87 
    ##              1.496              1.496              1.496              1.496 
    ##               TBCC               CTU2           METTL21B              CREG2 
    ##              1.496              1.496              1.497              1.497 
    ##            PRELID2              STX10              FCRLB           SERPIND1 
    ##              1.497              1.497              1.497              1.497 
    ##             MRPL49              OR4D6             DYNLT1             SUCLG2 
    ##              1.497              1.497              1.497              1.497 
    ##               ASB5            SPATA21               LHPP              RPP25 
    ##              1.498              1.498              1.498              1.498 
    ##             HHIPL2             NIPAL1               ALG1              NRSN1 
    ##              1.499              1.499              1.499              1.499 
    ##               CETP              PSMF1             ZNF214             DHRS7C 
    ##              1.499              1.499              1.499              1.500 
    ##              ECT2L              LCE3E           C1QTNF9B            TMEM190 
    ##              1.500              1.500              1.501              1.501 
    ##             RABL2B             GNPDA1              HDDC3              NRARP 
    ##              1.501              1.501              1.501              1.501 
    ##              OR8A1           CYB561D1          SERPINA12             GPR151 
    ##              1.501              1.501              1.501              1.502 
    ##                PAH             GXYLT2              PTGDR            TRIM49B 
    ##              1.502              1.502              1.502              1.502 
    ##             B3GNT4               PPIC           C11orf71               CBY1 
    ##              1.502              1.502              1.502              1.502 
    ##            TMEM217            TMEM74B             ATP5J2             P2RY13 
    ##              1.503              1.503              1.503              1.503 
    ##          LINC00923            GOLGA6A              CRYGD           EIF4EBP1 
    ##              1.504              1.504              1.504              1.505 
    ##              FABP2               PIGL               BTG2            SLC52A1 
    ##              1.505              1.506              1.506              1.506 
    ##      RP11-111M22.2            LGALS9C               NKD2            SLC18B1 
    ##              1.506              1.506              1.506              1.506 
    ##             MAPK15            KLHDC7B             PCDHA2              RPRML 
    ##              1.506              1.507              1.507              1.507 
    ##           C19orf10              GDPD3               CKLF            GOLGA8M 
    ##              1.507              1.507              1.508              1.508 
    ##            C8orf44               NQO1             HDAC10             OR4K15 
    ##              1.508              1.508              1.508              1.508 
    ##              DYRK4              MS4A7             SPATA4             TIMM23 
    ##              1.508              1.508              1.508              1.508 
    ##           TRAPPC3L              GPR35             NAMPTL       CTD-2510F5.6 
    ##              1.509              1.509              1.509              1.509 
    ##              WDR86           TSNAXIP1             ZNF589               DUXA 
    ##              1.510              1.510              1.510              1.510 
    ##               EAPP             ZNF257           C14orf28              ASB17 
    ##              1.510              1.510              1.510              1.511 
    ##             OR6C75              CREG1           C6orf229              MICU2 
    ##              1.511              1.511              1.511              1.511 
    ##             ZNF668            CCDC148       ZHX1-C8ORF76             DNALI1 
    ##              1.512              1.512              1.512              1.512 
    ##              CD160               OTOR               GGT6            TRIM64C 
    ##              1.513              1.513              1.513              1.513 
    ##              IFI30               CTU1             GTPBP8               ACPT 
    ##              1.513              1.513              1.513              1.513 
    ##              HAUS2            TMEM202            DCUN1D2             ZNF611 
    ##              1.513              1.513              1.513              1.513 
    ##                 C9           C15orf38              OR2A2              SMIM7 
    ##              1.514              1.514              1.514              1.515 
    ##             NKX2-8               ODAM            SLC30A8              FOXB2 
    ##              1.515              1.515              1.515              1.516 
    ##              SMCO2             DNAH14           C16orf71              OR9A2 
    ##              1.516              1.517              1.517              1.517 
    ##             EIF2B2            FAM213B              OR2T2               CTRL 
    ##              1.517              1.517              1.517              1.517 
    ##               EVX1              ZNF66             IFI44L             CXCL17 
    ##              1.517              1.517              1.517              1.518 
    ##               RGCC              GINS3              PRR18           C11orf54 
    ##              1.518              1.518              1.518              1.518 
    ##              ITIH6         AP001652.1             MRPL13             PRSS42 
    ##              1.518              1.518              1.518              1.519 
    ##             ZNF682             CHCHD6             ECHDC1               DPH3 
    ##              1.519              1.519              1.519              1.520 
    ##            CEACAM3              MNAT1             OR51G1         AC011897.1 
    ##              1.520              1.520              1.520              1.520 
    ##             CAPN12              ULBP1            FAM162A            CYP4F11 
    ##              1.520              1.520              1.521              1.521 
    ##           C11orf16              IFNL1            IFI27L1               SRMS 
    ##              1.521              1.521              1.521              1.522 
    ##              TSSK4              HOXB7             AKR7A2            SLC16A5 
    ##              1.522              1.522              1.522              1.522 
    ##             RSC1A1             CAMKMT           SLC25A41             KRBOX1 
    ##              1.522              1.523              1.523              1.523 
    ##               MMP7               COA5              PAQR6              IGFL4 
    ##              1.523              1.523              1.523              1.523 
    ##               APLN              TAAR8                LPA     RP11-1167A19.2 
    ##              1.524              1.524              1.524              1.524 
    ##                IL4             OR51S1              MMP13              FANCF 
    ##              1.524              1.524              1.524              1.525 
    ##               MMP1             ZNF763             NLRP14              KRT15 
    ##              1.525              1.525              1.525              1.526 
    ##             ZDHHC4              EVA1A             GPR135             KLHDC9 
    ##              1.526              1.526              1.526              1.526 
    ##           KIAA0825            CEACAM7             DRAXIN           C15orf32 
    ##              1.526              1.527              1.527              1.527 
    ##               CFL2             LDOC1L               FCN2             ZNF418 
    ##              1.527              1.527              1.527              1.527 
    ##            C7orf76            PPP1R3D             SNAPC1              IFIT1 
    ##              1.527              1.527              1.527              1.528 
    ##             OR52H1               ETV7         AC012123.1               CBY3 
    ##              1.528              1.528              1.528              1.528 
    ##              ACOT1               NANP               APLF            SERTAD3 
    ##              1.528              1.528              1.528              1.528 
    ##               RBKS             HOXD12          LINC00955             RNF182 
    ##              1.529              1.529              1.529              1.529 
    ##              FBXW5              WDR83             GPR142              CHEK2 
    ##              1.529              1.529              1.529              1.530 
    ##           SERPINB2            MICALCL              ASGR2              CFHR1 
    ##              1.530              1.530              1.530              1.530 
    ##             UGT2A3             FKBP11             DBNDD1           HSD17B14 
    ##              1.530              1.530              1.531              1.531 
    ##            LGALS14              ASF1B         KRTAP10-10            SLC28A2 
    ##              1.531              1.531              1.531              1.531 
    ##           SERPINI2            TACSTD2            MPV17L2             TAS2R5 
    ##              1.531              1.531              1.531              1.531 
    ##               WWOX               TSPO           KIAA0040              TMUB2 
    ##              1.532              1.532              1.532              1.532 
    ##               KLK9            FAM127A          TNFRSF12A             HOXD11 
    ##              1.532              1.532              1.532              1.533 
    ##             UGT2A1              TSACC               RTP3              ZNF57 
    ##              1.533              1.533              1.533              1.533 
    ##             PRADC1              TUBB8            ALDH3A1               CBR1 
    ##              1.533              1.533              1.533              1.533 
    ##                BAD            C2orf54           C18orf32             AQP12A 
    ##              1.533              1.534              1.534              1.534 
    ##               COA3               QPCT            C9orf69             GPR150 
    ##              1.534              1.534              1.534              1.535 
    ##              MFSD9              PSMG3              HOXA6         AC092687.4 
    ##              1.535              1.535              1.535              1.535 
    ##                MT3              RCCD1             LTB4R2            NDUFA11 
    ##              1.536              1.536              1.536              1.536 
    ##        RP1-66C13.4             OR51D1         AC022532.1           C1orf100 
    ##              1.536              1.536              1.537              1.537 
    ##            TM4SF19              DCTN3            C5orf52              DTYMK 
    ##              1.537              1.537              1.537              1.538 
    ##            LGALS13               DHDH              LYRM4             B3GNT8 
    ##              1.538              1.538              1.538              1.538 
    ##              SEPW1             IL11RA              AGMAT             AGPAT2 
    ##              1.538              1.538              1.539              1.539 
    ##           KRTAP2-2              ZNF77              SMIM4         AP000349.1 
    ##              1.539              1.539              1.539              1.539 
    ##             SLFN13               CMC2              TREX2              CPNE3 
    ##              1.540              1.540              1.540              1.540 
    ##              FAM3B                NGB              CABP4              OR2K2 
    ##              1.540              1.540              1.541              1.541 
    ##            ZDHHC11              HEBP2              RNF32            DNAJC22 
    ##              1.541              1.541              1.541              1.541 
    ##            TMEM252             KRBOX1              SMLR1             OR51M1 
    ##              1.541              1.541              1.541              1.541 
    ##              LETM2              REG1A              COX17              CENPK 
    ##              1.541              1.541              1.541              1.542 
    ##              CERKL              OR4P4           C14orf23             FAM25A 
    ##              1.542              1.542              1.542              1.542 
    ##            CYP24A1             COX7A2            RPS6KB2              TARM1 
    ##              1.542              1.542              1.542              1.542 
    ##           PRAMEF10              ADRB1            CCDC134               GPR3 
    ##              1.542              1.542              1.542              1.543 
    ##           C1orf227            SULT1A1           SLC22A10            PLEKHG7 
    ##              1.543              1.543              1.543              1.543 
    ##              NRIP2             FN3KRP             ZNF596         AC016251.1 
    ##              1.544              1.544              1.544              1.544 
    ##              RAB3D              TEX33             GLTPD2              OR6X1 
    ##              1.544              1.544              1.544              1.544 
    ##             EFCAB2            TPD52L3             NKX2-4              NLRX1 
    ##              1.544              1.545              1.545              1.545 
    ##             C2CD4A            SLCO1B3              GSTM5          SERPINA11 
    ##              1.545              1.545              1.545              1.545 
    ##              IL17D             LGALS1               IL24              OBP2B 
    ##              1.545              1.545              1.545              1.545 
    ##            SLC10A6              GNAT1             S100A7              PDRG1 
    ##              1.546              1.546              1.546              1.547 
    ##             CALML4                TES               SAT2          C14orf159 
    ##              1.547              1.547              1.548              1.548 
    ##              IFIH1              MTHFS             MOGAT1         AC005003.1 
    ##              1.548              1.548              1.548              1.548 
    ##               FIS1              BAALC             TPRX2P              ABCG2 
    ##              1.548              1.549              1.549              1.549 
    ##             OR52E6               ARSF              ISOC2            MEF2BNB 
    ##              1.549              1.549              1.549              1.549 
    ##               UCP2              THTPA              IRAK3             RASSF4 
    ##              1.549              1.549              1.549              1.549 
    ##             GTPBP6              CASP5              BET1L              ATP5L 
    ##              1.549              1.550              1.550              1.550 
    ##           C20orf62            DNAJC5G            TMEM141           C18orf21 
    ##              1.550              1.550              1.551              1.551 
    ##             ZNF470             ROPN1L             TMEM95              SNX20 
    ##              1.551              1.551              1.551              1.552 
    ##           C12orf54              NRIP3          C14orf180            BHLHE23 
    ##              1.552              1.553              1.553              1.553 
    ##           C1orf158               CAMP               PLP2               OASL 
    ##              1.553              1.553              1.553              1.553 
    ##              HOXC5              FITM2              LRIT2              RAB25 
    ##              1.553              1.553              1.554              1.554 
    ##             TIMM13             SPATA3              OR1K1            GOLGA8G 
    ##              1.554              1.554              1.554              1.554 
    ##           PPP1R14B              NINJ2              TEKT4              APOL4 
    ##              1.554              1.554              1.555              1.555 
    ##             ZNF267            C1orf95              MIXL1               CIB4 
    ##              1.555              1.555              1.555              1.555 
    ##              CXCR6            TMEM52B         AL589765.1             WFDC13 
    ##              1.555              1.555              1.555              1.555 
    ##              KLF16                UCN             HRASLS             TXNL4B 
    ##              1.556              1.556              1.556              1.556 
    ##              RMDN2                PLN               LGSN             ZNF277 
    ##              1.556              1.557              1.557              1.557 
    ##              BCAS4          TNFRSF10A             GPR132            CCDC160 
    ##              1.557              1.557              1.558              1.558 
    ##              OR8G5             TMEM25               MFNG               XKR9 
    ##              1.558              1.558              1.559              1.559 
    ##       RP11-934B9.3          ARHGAP11B          C14orf183               EMP3 
    ##              1.559              1.559              1.559              1.560 
    ##           C19orf81            PPP1R17             OR10G8                UBB 
    ##              1.560              1.560              1.560              1.560 
    ##              TOMM6            RPL22L1               TRMU              GNGT2 
    ##              1.561              1.561              1.561              1.561 
    ##             RNF148               GJA9           SIGLEC11              KRT82 
    ##              1.561              1.561              1.561              1.562 
    ##               GPX6               ADH5            FAM194A             SLFN11 
    ##              1.562              1.562              1.562              1.562 
    ##             OR11H6              AGAP6              TPGS1              MSRB2 
    ##              1.562              1.562              1.562              1.563 
    ##             TSGA13             TAS1R3             MRPS21               NME4 
    ##              1.563              1.563              1.564              1.565 
    ##              MZT2B             TUBA3D             TRIM58              LRIT1 
    ##              1.565              1.565              1.565              1.566 
    ##           TCTEX1D1         AC073343.1              TMCO2              NRSN2 
    ##              1.566              1.566              1.566              1.566 
    ##              KLK14             SWSAP1            ZSCAN32            ST8SIA6 
    ##              1.567              1.567              1.567              1.567 
    ##             ZNF799            FAM127B             STOML1            TMEM239 
    ##              1.567              1.567              1.567              1.567 
    ##               VNN2               IRF7            CD300LD               GYG2 
    ##              1.568              1.568              1.568              1.568 
    ##             KCNJ14           HLA-DPA1             ADAM18               EPGN 
    ##              1.568              1.568              1.569              1.569 
    ##               AQP2             GRXCR2               PBLD               NPPC 
    ##              1.569              1.569              1.570              1.570 
    ##             SMIM13             PLA2G5               POMK              ADAT2 
    ##              1.570              1.570              1.570              1.570 
    ##             ACTL7B             NKX6-3             CLEC4D            GOLGA8Q 
    ##              1.571              1.571              1.571              1.571 
    ##              BBIP1             CD300C               EREG             CYB5R2 
    ##              1.571              1.571              1.571              1.572 
    ##           C17orf58             PROKR2              FTSJ2            CYP11B2 
    ##              1.572              1.572              1.572              1.572 
    ##             HENMT1             PCDHA4             N6AMT2             OR5AP2 
    ##              1.572              1.572              1.573              1.573 
    ##               TEX9              RAMP3           C10orf35             HIGD2A 
    ##              1.573              1.574              1.574              1.574 
    ##             OR2A25             OR52A1              FDFT1            C5orf60 
    ##              1.575              1.575              1.575              1.575 
    ##            FAM180A               UBL5            CD164L2               TAF8 
    ##              1.575              1.575              1.575              1.575 
    ##               B9D2             RASSF7              HSPB7        RP1-139D8.6 
    ##              1.575              1.575              1.576              1.576 
    ##            TMEM215              CASQ1               ASIP           CEACAM19 
    ##              1.576              1.576              1.576              1.576 
    ##              AKIP1            C7orf57             OR2AG2           PPP1R14D 
    ##              1.576              1.576              1.576              1.576 
    ##             ZNF544               IQCE           C1orf110                MLN 
    ##              1.577              1.577              1.577              1.577 
    ##            THUMPD1       CTC-260F20.3            C1orf54             SYNGR2 
    ##              1.577              1.577              1.577              1.578 
    ##            FAM71E1             ZNF492            TMEM186             ZNHIT1 
    ##              1.578              1.578              1.578              1.578 
    ##              CES5A            C3orf80              CLDN8           S100A7L2 
    ##              1.578              1.579              1.579              1.579 
    ##               MSMB     CTD-3214H19.16              AUNIP              CSRP2 
    ##              1.579              1.579              1.579              1.579 
    ##           SLC22A25            RPL26L1           C12orf56        AC002310.13 
    ##              1.579              1.579              1.579              1.579 
    ##              OR3A2              LYZL6             IL15RA         AC021218.2 
    ##              1.580              1.580              1.580              1.580 
    ##         AL354898.1              NXPE4               SRRD             DNMT3A 
    ##              1.580              1.580              1.581              1.581 
    ##          GABARAPL2         AL441883.1              DPPA5               COQ4 
    ##              1.581              1.581              1.581              1.581 
    ##             OR2L13            NDUFAF7            ABHD12B               AGR2 
    ##              1.581              1.582              1.582              1.582 
    ##            C2orf70              CRIPT             CHMP4A             COX6A2 
    ##              1.582              1.582              1.582              1.582 
    ##             CYP2A7            DNAJC5B              SMPD2               PRCD 
    ##              1.582              1.583              1.583              1.583 
    ##           C1orf105              KRT74               TSLP            DHRS4L2 
    ##              1.583              1.583              1.584              1.584 
    ##            TOMM20L              CLDN6             PTGER1             TATDN1 
    ##              1.584              1.584              1.584              1.584 
    ##              YIF1A               FMO3               GNG5             TMIGD2 
    ##              1.584              1.584              1.585              1.585 
    ##              KLF14              PGBD2              GPR84            SULT1A2 
    ##              1.585              1.585              1.585              1.585 
    ##              C3AR1              WFDC3             OR2B11             SCNN1D 
    ##              1.585              1.585              1.585              1.586 
    ##            LGALS12              AKR7L             OR5B17               TFF3 
    ##              1.586              1.586              1.586              1.586 
    ##               HMX2              ASB16               PIFO              SIAH3 
    ##              1.587              1.587              1.587              1.587 
    ##           C21orf67                VIT             AKR1C1            FLYWCH2 
    ##              1.588              1.588              1.588              1.588 
    ##              OXGR1              OR4D9              TTC36              UTS2B 
    ##              1.588              1.588              1.588              1.588 
    ##               BBC3              RAB26              WDR38              SSTR5 
    ##              1.589              1.589              1.589              1.589 
    ##          HIST2H2AC              INSL3            APOBEC1             OR2T35 
    ##              1.589              1.590              1.590              1.590 
    ##              GBGT1               RBP5              GIPC2               MEA1 
    ##              1.590              1.590              1.591              1.591 
    ##             OR52B6            UGT2B11             FIBCD1            PCDHB11 
    ##              1.591              1.591              1.591              1.591 
    ##              RPP21              OR6Y1             MRPS16             OR10W1 
    ##              1.591              1.591              1.592              1.592 
    ##               MRRF            SLC18A1               CD14            SLC36A4 
    ##              1.592              1.592              1.592              1.592 
    ##            TMEM187             GLIPR2              GLOD4            C3orf70 
    ##              1.592              1.593              1.593              1.593 
    ##            C4orf45           HIST1H1E              SEMG2              RSPO4 
    ##              1.593              1.593              1.593              1.593 
    ##              RDH16             ZSCAN1               RHBG        RP5-850E9.3 
    ##              1.593              1.593              1.593              1.594 
    ##              KLK10            FAM198A             OR4K14               SPP1 
    ##              1.594              1.594              1.594              1.595 
    ##           C19orf18              NTHL1             OR4D10              PRDX6 
    ##              1.595              1.595              1.595              1.595 
    ##             CLPSL2            TMEM242             SPINK9             ENTPD3 
    ##              1.595              1.595              1.595              1.596 
    ##              APOA5               CST3              CALCA           KRTAP1-5 
    ##              1.596              1.596              1.596              1.596 
    ##             ZNF660              RAB20            SLC35D2               BIN3 
    ##              1.596              1.596              1.596              1.597 
    ##            ZSCAN5C               PPBP              PLIN4              TEX36 
    ##              1.597              1.597              1.597              1.597 
    ##               MAFF               ZBP1           MAP3K7CL             GTF2H2 
    ##              1.597              1.597              1.597              1.597 
    ##              MUC12              IFT43               DARC             CKMT1A 
    ##              1.598              1.598              1.598              1.598 
    ##             SS18L2            KLHDC7A            TM4SF18               RDM1 
    ##              1.598              1.598              1.598              1.598 
    ##              AP3S1            ANXA8L2                HPR              OR1D2 
    ##              1.599              1.599              1.599              1.599 
    ##       RP11-20I23.1           HLA-DQA2              OR1C1              GTF3A 
    ##              1.600              1.600              1.600              1.600 
    ##              DLEU7             OR2A42             WBP2NL               IL25 
    ##              1.600              1.600              1.600              1.600 
    ##               PPIB              SMKR1            SCGB1D4          KRTAP11-1 
    ##              1.600              1.600              1.600              1.601 
    ##              OR4D2            NDUFAF2             EIF5A2              RGS16 
    ##              1.601              1.601              1.601              1.601 
    ##              OR6K6           HLA-DPB1               POP5             KCNMB1 
    ##              1.601              1.602              1.602              1.602 
    ##             GLYCTK               BET1              PAGE2               MC4R 
    ##              1.602              1.602              1.602              1.602 
    ##           PDZK1IP1              PILRA              RPL41              CFHR3 
    ##              1.602              1.603              1.603              1.603 
    ##      RP11-195F19.5             FAM96B          SPATA31D1               GZMA 
    ##              1.603              1.603              1.604              1.604 
    ##              ADH1A            FAM122C       RP11-770J1.4            CLEC19A 
    ##              1.604              1.604              1.604              1.605 
    ##              EMC10            TMEM88B            FAM177B              DTWD2 
    ##              1.605              1.605              1.605              1.605 
    ##             FAM83E             MS4A13               DPP7             SPIN2B 
    ##              1.605              1.605              1.605              1.605 
    ##             ELMOD3         AL627309.1             ZNF676              CCNI2 
    ##              1.606              1.606              1.606              1.606 
    ##             GAGE2D              RITA1               NEU4             OR51L1 
    ##              1.606              1.606              1.606              1.606 
    ##              CLRN2         AC079612.1              RPL3L              SMDT1 
    ##              1.606              1.607              1.607              1.607 
    ##            C1orf56               PGK2             PCDHA9              TSPO2 
    ##              1.608              1.608              1.608              1.608 
    ##            CLEC18C             TBC1D7             TRIAP1                GP6 
    ##              1.608              1.608              1.608              1.609 
    ##              FRG2B               TEN1               SLPI           C1orf186 
    ##              1.609              1.609              1.609              1.609 
    ##             LILRA1            C3orf79            C2orf91            SLC22A9 
    ##              1.609              1.610              1.610              1.610 
    ##             ACTRT1            KIR2DL3           C19orf53            C8orf76 
    ##              1.611              1.611              1.611              1.611 
    ##              ASB11          NIPSNAP3B            RARRES1           C13orf35 
    ##              1.611              1.611              1.611              1.612 
    ##               RDH5              DUSP2              MTIF3              THEM4 
    ##              1.612              1.612              1.612              1.612 
    ##              REG1B              PTGIR               FHL5              NUBP1 
    ##              1.612              1.612              1.613              1.613 
    ##               SAA1              LAMP5              GSTT2               STAR 
    ##              1.613              1.613              1.613              1.613 
    ##            TMEM154              FOXE3              COPZ2              TCF23 
    ##              1.613              1.613              1.614              1.614 
    ##              FANCL              KRT16              GSTO2            PPP1R3G 
    ##              1.614              1.614              1.614              1.614 
    ##              NBPF6             MAGEB3             PRSS45        AP000322.53 
    ##              1.614              1.614              1.614              1.614 
    ##              KDM4D             TMEM52           PHOSPHO2            FAM173A 
    ##              1.614              1.614              1.615              1.615 
    ##              RPP38             C9orf9               GYPE             OR51G2 
    ##              1.615              1.616              1.616              1.616 
    ##            TMEM128             MRPS33              OR5P2            SDR39U1 
    ##              1.616              1.616              1.616              1.616 
    ##               ADH7             ZNF337             ZNF681              GSTM4 
    ##              1.617              1.617              1.617              1.617 
    ##           KIAA0125      CTD-2368P22.1               NOL3            ELSPBP1 
    ##              1.617              1.617              1.617              1.617 
    ##             DUSP28      RP11-347C12.3            ZC2HC1B             TM4SF4 
    ##              1.618              1.618              1.618              1.618 
    ##            SDR42E2              DHRS9              HOXA9             RABL2A 
    ##              1.618              1.618              1.618              1.618 
    ##               CTSE             OR2AT4      RP11-680G10.1             NUDT18 
    ##              1.618              1.618              1.618              1.619 
    ##            ATP6V1F             ZSCAN4              NOXO1             RWDD2A 
    ##              1.619              1.619              1.619              1.619 
    ##               AIF1             HSD3B2              ASB10               CR1L 
    ##              1.619              1.619              1.619              1.619 
    ##              TPBGL       PRR5-ARHGAP8               OLAH             SMIM17 
    ##              1.620              1.620              1.620              1.620 
    ##              OR5M9              GLRX2        CTA-299D3.8            C2orf73 
    ##              1.620              1.621              1.621              1.621 
    ##               UCMA              LYRM9            CD300LB               IFNE 
    ##              1.621              1.621              1.621              1.621 
    ##             ZNF551             OR51I1              LYPD4         AC019171.1 
    ##              1.622              1.622              1.622              1.622 
    ##              RRP36             SLC3A1                GPT            C1QTNF6 
    ##              1.622              1.623              1.623              1.623 
    ##               MYL7         AP000769.1              FOPNL             PODNL1 
    ##              1.623              1.624              1.624              1.624 
    ##              FKBPL              SNX16               SNCG             GATSL2 
    ##              1.624              1.624              1.625              1.625 
    ##              LELP1              LTC4S               FUOM              HPGDS 
    ##              1.626              1.626              1.626              1.626 
    ##             OR10V1               CYBA               RSG1           SERPINA1 
    ##              1.626              1.626              1.626              1.626 
    ##              RSPH9         AC005481.5            HEATR5A             GPRC5A 
    ##              1.626              1.626              1.626              1.627 
    ##             CELA3B             ZNF440              OR5A1             PLCXD1 
    ##              1.627              1.627              1.627              1.627 
    ##            TNFSF18           C16orf46              OR8G1              OR9Q1 
    ##              1.627              1.627              1.627              1.628 
    ##            SLFN12L            LAMTOR4               IL13              HTR3E 
    ##              1.628              1.628              1.628              1.629 
    ##            GOLGA8F              SMCO3            C7orf62              RABIF 
    ##              1.629              1.629              1.629              1.629 
    ##              KCNRG           GOLGA6L1           CDC42SE1              F2RL3 
    ##              1.630              1.630              1.630              1.630 
    ##                LBH           C19orf70             ZNF100              HCAR2 
    ##              1.630              1.630              1.630              1.630 
    ##            SLC46A2              BANF2              OR8B8              NAPSA 
    ##              1.630              1.631              1.631              1.631 
    ##              ACOT6            C2orf80                NMS               EXOG 
    ##              1.631              1.631              1.631              1.631 
    ##               IL22           C12orf79             SELPLG               FCAR 
    ##              1.631              1.631              1.631              1.632 
    ##             OR10G3            TMEM140              H2AFX              HOXA4 
    ##              1.632              1.632              1.632              1.632 
    ##           COX6A1P2                GP5              GRTP1             TMEM69 
    ##              1.632              1.632              1.633              1.633 
    ##              SETD6              GSTA1               CPB1            DNASE2B 
    ##              1.633              1.633              1.633              1.633 
    ##             SIRPB1              TOR3A             ZNF233              PRR22 
    ##              1.633              1.633              1.633              1.633 
    ##              AP4S1           MARCKSL1                DAP            TM4SF20 
    ##              1.633              1.633              1.634              1.634 
    ##             RAD51B             ALKBH7               LY6K              CDRT1 
    ##              1.634              1.634              1.634              1.634 
    ##              MPEG1            GALNTL5         AC027763.2            C4orf40 
    ##              1.634              1.634              1.635              1.635 
    ##             S100A6              KLRC3              DYDC2           TMEM126B 
    ##              1.635              1.635              1.635              1.635 
    ##              UPK1A            IGFBPL1             BPIFA1              MCTP2 
    ##              1.635              1.635              1.636              1.636 
    ##               AIM2            TMEM239               HAP1            UGT2B10 
    ##              1.636              1.636              1.637              1.637 
    ##              SATL1               SBK2             CALHM2              ENDOG 
    ##              1.638              1.638              1.638              1.638 
    ##              GLOD5              NR0B2             OR13C2            WBSCR28 
    ##              1.638              1.638              1.638              1.638 
    ##               GMFG             TIMM22             NUDT14              ZNF98 
    ##              1.638              1.638              1.639              1.639 
    ##            FAM216B               DDI1             MRPS17              H2BFM 
    ##              1.639              1.639              1.639              1.639 
    ##             ERICH1            FAM159B       RP11-211G3.3              BHMT2 
    ##              1.639              1.639              1.640              1.640 
    ##               PTX4            FAM132A            TMEM230            RPSAP58 
    ##              1.640              1.640              1.641              1.641 
    ##             LY6G5B              ASB18            C4orf22             UBQLNL 
    ##              1.641              1.641              1.642              1.642 
    ##           ZDHHC11B               AZU1            LDHAL6A             SRSF10 
    ##              1.642              1.642              1.642              1.642 
    ##              JOSD2              AHSA2            CYP4A11            TMEM159 
    ##              1.642              1.642              1.643              1.643 
    ##             FABP12              IGLL1             MRPL44              P2RY8 
    ##              1.643              1.643              1.644              1.644 
    ##              MXRA7              AMY1C          C20orf197           C17orf78 
    ##              1.644              1.644              1.644              1.644 
    ##             PCP4L1              PROK2            PRR23D1             ENDOD1 
    ##              1.644              1.645              1.645              1.645 
    ##                NTS            SEC14L3            C3orf83         AC112715.2 
    ##              1.645              1.646              1.646              1.646 
    ##               MAFA               CST8              OR8B3             CCT8L2 
    ##              1.646              1.646              1.646              1.647 
    ##              FABP7              OR2C1             MRPL50              TTC34 
    ##              1.647              1.647              1.647              1.647 
    ##            CCDC28A            PLEKHS1             OR6C74              OR6J1 
    ##              1.647              1.647              1.647              1.647 
    ##              PRTN3            RASL11A             FAM89B              CHST4 
    ##              1.647              1.647              1.648              1.648 
    ##            GLYATL2              RGS21             RNF151               IFI6 
    ##              1.648              1.648              1.648              1.648 
    ##               CUTA             ZSWIM1              SH2D6              ASCL3 
    ##              1.649              1.649              1.649              1.649 
    ##              CDK20         AC091801.1              IGFL2             OR10G2 
    ##              1.649              1.649              1.650              1.650 
    ##                MPG               SPG7             TP53RK             FCGR3B 
    ##              1.650              1.651              1.651              1.651 
    ##             SAYSD1               INMT              TEX29            C8orf48 
    ##              1.651              1.651              1.651              1.651 
    ##      RP11-422N16.3            WDR83OS            FAM167B               GFER 
    ##              1.652              1.652              1.652              1.652 
    ##            AADACL2              PRR13              FABP6             RAB27A 
    ##              1.652              1.652              1.652              1.653 
    ##         AL021546.6      RP11-286N22.8            GOLGA8I           C1ORF220 
    ##              1.653              1.654              1.654              1.654 
    ##             NANOS2              CCL27              CCL28              LIME1 
    ##              1.654              1.654              1.654              1.654 
    ##              OR4D1             A4GALT            CCDC166              ACOT4 
    ##              1.654              1.654              1.655              1.655 
    ##           FLJ00273              SNX10           C11orf74            FAM19A3 
    ##              1.655              1.655              1.655              1.655 
    ##             FKBP14             ZNF221              EPHX3            SLC35E2 
    ##              1.655              1.655              1.655              1.656 
    ##             MOGAT3            ALDH3B2              TIRAP             PRSS46 
    ##              1.656              1.656              1.656              1.656 
    ##             TP53I3               PDYN              OBP2A              PLIN5 
    ##              1.656              1.656              1.656              1.656 
    ##            FAM166A           GOLGA6L6            TMEM156             SHISA3 
    ##              1.657              1.657              1.657              1.657 
    ##         AC110619.2         AC011366.3              HCAR1               DTD2 
    ##              1.657              1.657              1.657              1.657 
    ##               SCIN           APOBEC3D              PDE6G      RP11-322L20.1 
    ##              1.657              1.657              1.657              1.658 
    ##              FFAR1               TAC4            ZNF780B               RLN3 
    ##              1.658              1.658              1.658              1.658 
    ##            TMEM247             PRSS35             ZNF788                MGP 
    ##              1.659              1.659              1.659              1.659 
    ##                PI3              LCE1E            PCDH11Y            C6orf52 
    ##              1.659              1.659              1.659              1.659 
    ##             ATP5G2              PTRH1            PTPN20B              GP1BA 
    ##              1.659              1.659              1.659              1.660 
    ##             PEX11G               APOH          KRTAP10-1             OR5AC2 
    ##              1.660              1.660              1.660              1.660 
    ##              ATOH7              CAPSL              GPR42            FAM173B 
    ##              1.661              1.661              1.661              1.662 
    ##               AHSG                ISX         AC096582.1             SMIM18 
    ##              1.662              1.662              1.662              1.662 
    ##             LRRC52             KRT33B               PSCA              OR4C5 
    ##              1.662              1.662              1.663              1.663 
    ##             ZNF253              S100Z             MRPS14                 XG 
    ##              1.663              1.663              1.664              1.664 
    ##           SERPINA4              OR5V1            FAM221A            S100A10 
    ##              1.664              1.664              1.664              1.664 
    ##              WDR55            DGAT2L6             MPLKIP              CDKL4 
    ##              1.665              1.665              1.665              1.665 
    ##               CD5L               DEXI       RP11-62N21.1            TCEANC2 
    ##              1.666              1.666              1.666              1.666 
    ##              KCNK7           CNTNAP3B            C9orf62             NDUFB7 
    ##              1.666              1.666              1.667              1.667 
    ##             BTN3A3             PLSCR5             PABPC3              LCA10 
    ##              1.667              1.667              1.667              1.667 
    ##           APOBEC3G             ZNF778           C1orf195              MMP26 
    ##              1.667              1.667              1.667              1.668 
    ##             NDUFB6        CTC-429P9.4    LL22NC03-63E9.3              OR2M4 
    ##              1.668              1.668              1.668              1.669 
    ##             SLC2A7           SLC25A47              RTBDN              THAP8 
    ##              1.669              1.669              1.669              1.669 
    ##           C6orf118            CAMK2N1             ZNF688            TMEM253 
    ##              1.670              1.670              1.670              1.670 
    ##            RIMBP3B               SOD3              CSDC2             ZNF506 
    ##              1.671              1.671              1.671              1.671 
    ##       RP11-484M3.5      RP11-597K23.2             ZNF224             OR10A2 
    ##              1.671              1.671              1.671              1.672 
    ##                MOK              MED18            S100A13             MRPS26 
    ##              1.672              1.672              1.672              1.673 
    ##             ZNF587               HMSD            C4orf26              POTEH 
    ##              1.673              1.673              1.673              1.673 
    ##             OR52A5             TMSB4X            SLC22A1            PCDHB14 
    ##              1.673              1.673              1.673              1.674 
    ##              PLET1             KCNJ16              TEKT5                DCD 
    ##              1.674              1.674              1.674              1.674 
    ##             TMEM81             ZNF528             FAM26D               PCTP 
    ##              1.674              1.674              1.674              1.675 
    ##              CTXN1             GGTLC1            FAM219B       RP11-352D3.2 
    ##              1.675              1.675              1.675              1.675 
    ##             N6AMT1             TNFSF4              PGAM1             ZNF425 
    ##              1.675              1.675              1.676              1.676 
    ##              DECR1               IL20      CTD-2117L12.1             ZNF211 
    ##              1.676              1.676              1.676              1.676 
    ##              OR6C3              DTWD1             OR51B2              BEAN1 
    ##              1.676              1.677              1.677              1.677 
    ##              GSTA2              STX19              AMY2A               RTP2 
    ##              1.677              1.677              1.678              1.678 
    ##             LRRC26            U2AF1L4              SMCO1                CFD 
    ##              1.678              1.678              1.678              1.678 
    ##          KRTAP10-8              IGSF5         ST6GALNAC2          MAP1LC3B2 
    ##              1.678              1.678              1.678              1.678 
    ##             KCTD21              GPNMB             ZNF837            C1orf87 
    ##              1.678              1.678              1.678              1.679 
    ##              POTEG           C19orf52              CERS4            DEFB123 
    ##              1.679              1.679              1.679              1.679 
    ##            AADACL4             HPCAL4           KRTAP5-3              COX6C 
    ##              1.679              1.679              1.679              1.680 
    ##               AGRP             MRPS11              ASB12              POLE4 
    ##              1.680              1.680              1.680              1.680 
    ##             OR51E1             TOP1MT             ZNF749           C21orf49 
    ##              1.680              1.680              1.680              1.680 
    ##            C7orf55             OR51B4             OR51E2               GNG8 
    ##              1.680              1.681              1.681              1.681 
    ##             ZNF543         RP11-1C1.5            TMEM220             CRISP2 
    ##              1.681              1.681              1.681              1.681 
    ##            C9orf85                CA1               PRND            C3orf35 
    ##              1.682              1.682              1.682              1.683 
    ##            DEFB126             MAGEF1               CNR2     CTD-2207O23.12 
    ##              1.683              1.683              1.683              1.683 
    ##           SLC25A48              ADRB3              DECR2               PRH2 
    ##              1.683              1.683              1.684              1.684 
    ##              VSIG4            RPL36AL                PTS            ZFAND2A 
    ##              1.684              1.684              1.685              1.685 
    ##           C12orf60            PRTFDC1             EXOSC8            FAM203A 
    ##              1.685              1.685              1.685              1.686 
    ##              UQCRQ           C1orf170              CMTM3             UNC93A 
    ##              1.686              1.686              1.686              1.686 
    ##               ENO3               GPX2              GCHFR              RFPL1 
    ##              1.686              1.686              1.686              1.686 
    ##                PIP         AC093157.1             TRIM54              GALR2 
    ##              1.686              1.687              1.687              1.687 
    ##      RP11-520P18.5               ADI1            TBC1D29               SBK3 
    ##              1.687              1.687              1.687              1.687 
    ##              FRAT2            BLOC1S6              USMG5               SAA4 
    ##              1.688              1.688              1.688              1.688 
    ##           KRTAP4-5        AP000350.10              DEFB1               VIMP 
    ##              1.689              1.689              1.689              1.689 
    ##           TMEM255B             KRT33A           C11orf52                MIF 
    ##              1.689              1.689              1.690              1.690 
    ##              IL1R2           C1orf131            HSBP1L1              PYURF 
    ##              1.690              1.690              1.690              1.690 
    ##                HBM             CACNG6              TM2D3            UGT2B28 
    ##              1.690              1.690              1.690              1.691 
    ##         AL450307.1             GIMAP2              ASB14              CXCL9 
    ##              1.691              1.691              1.691              1.691 
    ##             ZNF433               RILP              ISG15            ONECUT3 
    ##              1.691              1.691              1.691              1.692 
    ##              AGXT2         AP000350.4             ZNF138              THRSP 
    ##              1.692              1.692              1.692              1.692 
    ##             RAET1E               CTRC              PLAC9            DEFB127 
    ##              1.692              1.692              1.692              1.692 
    ##            C2orf72             DUSP23              HYAL4              OR4B1 
    ##              1.692              1.692              1.692              1.692 
    ##              FAM9B             ZNF880             OR52W1             ZNF439 
    ##              1.693              1.693              1.693              1.693 
    ##           C1orf192              RERGL               SELO              ICAM4 
    ##              1.693              1.693              1.693              1.693 
    ##               NKG7               NXF5             ZNF593            C3orf18 
    ##              1.693              1.693              1.693              1.694 
    ##                RGR              MRPS6             MRPS17             PACRGL 
    ##              1.694              1.694              1.694              1.694 
    ##              ZNF90              ISG20                PTH              SAMD5 
    ##              1.694              1.694              1.695              1.695 
    ##         AL645608.1           L34079.2            GOLGA6B      RP11-664I21.6 
    ##              1.695              1.695              1.695              1.695 
    ##           C21orf58              PATE3               GATC            ABHD14B 
    ##              1.695              1.695              1.695              1.695 
    ##          SERPINA10              TLCD2             GIMAP1              MMP27 
    ##              1.696              1.696              1.696              1.696 
    ##              FKBP7               PCP2              MUC3A             EIF3CL 
    ##              1.696              1.696              1.696              1.696 
    ##               CBR4           CEACAM21             ZNF561              PAQR4 
    ##              1.696              1.696              1.696              1.697 
    ##              COX18           C19orf25            CYSLTR1             ASAH2B 
    ##              1.697              1.697              1.697              1.697 
    ##                DAO               NTF4         AC040977.1           SLC38A11 
    ##              1.697              1.697              1.697              1.697 
    ##               OOEP            ANKRD60               MT1F              RXFP4 
    ##              1.697              1.698              1.698              1.698 
    ##        CTB-167G5.5               CD52              C5AR1               ELP5 
    ##              1.698              1.698              1.698              1.698 
    ##         AP000889.3            TRIM49C             ZNF140           C1orf168 
    ##              1.699              1.699              1.699              1.699 
    ##           CEACAM18             ZNF648             CCDC59               GALP 
    ##              1.699              1.699              1.699              1.700 
    ##              TMED1              CMTM2            ZNF280A               RBP1 
    ##              1.700              1.700              1.701              1.701 
    ##            AADACL3              OR2F2              NXNL1               FPR2 
    ##              1.701              1.701              1.701              1.701 
    ##           C11orf45               DAOA            ZNF705A               MSR1 
    ##              1.701              1.701              1.701              1.701 
    ##           C16orf89            C5orf48              DNAL4                MVD 
    ##              1.702              1.702              1.702              1.702 
    ##             COX6B2              PSKH2               XAF1              PTPLA 
    ##              1.702              1.702              1.702              1.702 
    ##         AP000758.1             GRXCR1               SDPR             MRPL52 
    ##              1.702              1.702              1.702              1.702 
    ##              PTRH2              MOB3A            C5orf58               COMT 
    ##              1.703              1.703              1.703              1.703 
    ##               CMC4       RP13-512J5.1              CNBD1                PPY 
    ##              1.703              1.703              1.704              1.704 
    ##               GKN1            C5orf49           C22orf34       RP4-758J18.2 
    ##              1.704              1.704              1.704              1.704 
    ##             MTHFSD            C2orf68           C18orf56               AGMO 
    ##              1.705              1.705              1.705              1.705 
    ##               CSN3               PKIB                PYY             NKG2-E 
    ##              1.705              1.705              1.705              1.705 
    ##           TCTEX1D4             APITD1        APITD1-CORT               RCN3 
    ##              1.705              1.706              1.706              1.706 
    ##              THEM6              CPLX4              GPR15             OR2A14 
    ##              1.706              1.707              1.707              1.707 
    ##             NKX1-2               MZT1             ZNF549              GPR55 
    ##              1.707              1.707              1.707              1.707 
    ##             ZNF700             PLA2G7             PRSS57             SPAG16 
    ##              1.708              1.708              1.708              1.708 
    ##              LIPT1     XXcos-LUCA11.5           HIST1H1D             ACTRT3 
    ##              1.708              1.709              1.710              1.710 
    ##            CSNK2A3           C16orf13              GALR3             WFDC11 
    ##              1.710              1.710              1.710              1.711 
    ##               GAST             NEURL2              CIDEA         AC005609.1 
    ##              1.711              1.711              1.711              1.711 
    ##             TSSK1B            EFCAB4A             SPACA5             CRYBB3 
    ##              1.711              1.711              1.712              1.712 
    ##              GNA14               GNG4             RHBDD1             IZUMO4 
    ##              1.712              1.712              1.712              1.713 
    ##         AC132216.1             HOXB13            TMEM256               IQCG 
    ##              1.713              1.713              1.713              1.713 
    ##           SLC16A13             GTF2H5             SMIM15             CSN1S1 
    ##              1.713              1.714              1.714              1.714 
    ##              GLYAT              PDZK1             NDUFA6            DEFB128 
    ##              1.714              1.714              1.714              1.715 
    ##              LRRC3               IL31           C6orf100             PTGDR2 
    ##              1.715              1.715              1.715              1.715 
    ##               TCN1              FADS6              WFDC1           C14orf64 
    ##              1.715              1.715              1.716              1.716 
    ##             TRIM52              OR8H3            RARRES2               ADM5 
    ##              1.716              1.716              1.716              1.716 
    ##            C3orf65               LCN2             ZNF814             TREML2 
    ##              1.716              1.716              1.717              1.717 
    ##            LRRC10B             ZNF117               MBL2            ANKRD22 
    ##              1.717              1.717              1.717              1.718 
    ##            CCDC107              APOA2             SPESP1            FAM213A 
    ##              1.718              1.718              1.718              1.718 
    ##              HINT2           RPS19BP1           TMEM150B              RPEL1 
    ##              1.718              1.719              1.719              1.719 
    ##             EPSTI1               LIPI             AKR1C4         AL050302.1 
    ##              1.719              1.719              1.719              1.719 
    ##              MGST1              TPRG1              CCL23           SERPINA5 
    ##              1.720              1.720              1.720              1.720 
    ##              TSSC4             WFDC12              MORN2              RGPD3 
    ##              1.720              1.720              1.720              1.720 
    ##              SCRG1              XRCC3             LEFTY1            HCFC1R1 
    ##              1.721              1.721              1.721              1.722 
    ##             OR5AK2             MRPS34               ACN9              SSX4B 
    ##              1.722              1.722              1.722              1.722 
    ##             MAGEC3              AMY2B           C21orf62             OR2AE1 
    ##              1.722              1.723              1.723              1.723 
    ##               PRR4               HES2              SLX1B              PROX2 
    ##              1.723              1.723              1.723              1.724 
    ##              RGPD6              SURF1           C12orf71       RP11-47I22.3 
    ##              1.724              1.724              1.724              1.724 
    ##              CBWD6       RP11-944C7.1              APOL1           NUDT16L1 
    ##              1.725              1.725              1.725              1.725 
    ##         AP003774.4            HSD17B2             SSSCA1             TEX13B 
    ##              1.725              1.725              1.725              1.725 
    ##             ZNF121            FAM136A            SPANXN1            LRRC14B 
    ##              1.725              1.725              1.725              1.725 
    ##                RD3               UCN3              S100B              SMR3A 
    ##              1.726              1.726              1.726              1.726 
    ##              VAMP5             DNAJC4             STOML3           C9orf173 
    ##              1.726              1.727              1.727              1.727 
    ##               TXN2           C14orf79              LYZL1              IL17C 
    ##              1.727              1.727              1.727              1.727 
    ##            TMSB15A             SAPCD1              DYNAP             SAC3D1 
    ##              1.728              1.728              1.728              1.729 
    ##               CD1E         AC092675.3         AC026407.1           CDC42EP1 
    ##              1.729              1.729              1.729              1.729 
    ##               IL37             ZNF124             IZUMO3               PROZ 
    ##              1.729              1.730              1.730              1.730 
    ##             MRPS10             CCDC96         SLC22A18AS         AC009365.3 
    ##              1.730              1.730              1.730              1.730 
    ##            SLC7A13             POLR3G              VIPR1             OR13G1 
    ##              1.730              1.730              1.730              1.731 
    ##             PHLDA2             DEPDC4              KLRG2             OR11A1 
    ##              1.731              1.731              1.731              1.732 
    ##       RP11-180C1.1               ART1          C20orf141             CHST13 
    ##              1.732              1.732              1.732              1.732 
    ##           MRFAP1L1                GRP              MS4A2            KIR3DL3 
    ##              1.732              1.733              1.733              1.733 
    ##              MFSD7             SYCE1L                PF4       RP11-661C8.3 
    ##              1.733              1.733              1.733              1.733 
    ##            TMSB15B            UGT2B15              PROL1           C22orf24 
    ##              1.733              1.733              1.733              1.734 
    ##              OR5W2              DHRS1            BCL2L10             RFPL4A 
    ##              1.734              1.735              1.735              1.735 
    ##               PAEP           C12orf77             OR56A1               RPRM 
    ##              1.735              1.735              1.735              1.735 
    ##           C11orf94            PPAPDC2            FAM86B2            PLA2G2E 
    ##              1.735              1.735              1.736              1.736 
    ##              NEIL2              TACR3              PATE2              TUSC1 
    ##              1.736              1.736              1.736              1.736 
    ##            TMEM158         AC092964.1         AC007204.1              ASGR1 
    ##              1.737              1.737              1.737              1.737 
    ##              LACRT              APH1B               UCP3           HIST1H4J 
    ##              1.737              1.737              1.737              1.737 
    ##             ECHDC3             RETNLB              IFI27               MIIP 
    ##              1.737              1.737              1.737              1.737 
    ##          KRTAP10-4             OR52D1              SMR3B              LYPD1 
    ##              1.737              1.737              1.737              1.738 
    ##           PSORS1C2               LSM2              IL36G               OST4 
    ##              1.738              1.738              1.738              1.738 
    ##            SPATA6L             POTEB2             ZNF229               NTN5 
    ##              1.739              1.739              1.739              1.739 
    ##               CRCP      RP11-497E19.2               GJB4               GJB6 
    ##              1.739              1.739              1.739              1.739 
    ##            C2orf50              ATG12       RP11-343C2.9         AC114783.1 
    ##              1.739              1.740              1.740              1.740 
    ##                HYI             APOLD1             FBXO36         AC004076.7 
    ##              1.740              1.740              1.740              1.740 
    ##               HTN1              TMED5              OR1A1       RP11-166B2.1 
    ##              1.740              1.741              1.741              1.741 
    ##             H2BFWT               NPPA              SENP8        AC002472.13 
    ##              1.741              1.741              1.741              1.741 
    ##              OR5C1              OR6N1            SULT1A4           C11orf35 
    ##              1.742              1.742              1.742              1.742 
    ##              LENG1            SLCO1B7            UGT1A10               APIP 
    ##              1.742              1.742              1.742              1.742 
    ##       CTD-2203A3.1          CEBPZ-AS1            PPP1R1C             OR10H2 
    ##              1.742              1.742              1.742              1.742 
    ##               CNTF               ZFP2              SMIM2       RP11-171N4.2 
    ##              1.743              1.743              1.743              1.743 
    ##             ODF3L1          SPATA31A2            C8orf86               CTF1 
    ##              1.743              1.743              1.744              1.744 
    ##            MYADML2         AC002365.1         AC008686.1      RP11-451M19.3 
    ##              1.744              1.744              1.744              1.744 
    ##           ATXN7L3B              GSTO1            PCDHB13      RP11-706O15.1 
    ##              1.744              1.745              1.745              1.745 
    ##               OSTN             GGTLC2         AC093802.1              MROH5 
    ##              1.745              1.745              1.745              1.745 
    ##                CLC             ZNF222            TRAPPC2             TMEM27 
    ##              1.746              1.746              1.746              1.746 
    ##               CD70             CLDN23       RP11-362K2.2            SPANXN3 
    ##              1.746              1.746              1.746              1.746 
    ##          LINC00998             RESP18              SIRPG             GLIPR1 
    ##              1.747              1.747              1.747              1.747 
    ##          TMPRSS11B             TMEM37             LILRA2            C3orf55 
    ##              1.747              1.748              1.748              1.748 
    ##               CCL8              THAP2               TMA7             LGALS7 
    ##              1.748              1.748              1.748              1.748 
    ##              OR2F1               PSPN          KRTAP4-11          KRTAP5-11 
    ##              1.748              1.749              1.749              1.749 
    ##            C9orf57               VMAC             OR6C76             FKBP1C 
    ##              1.749              1.749              1.749              1.749 
    ##               SELK            C5orf66               LNP1         AC106876.2 
    ##              1.749              1.749              1.749              1.749 
    ##            ANGPTL7             ZNF582           C11orf65               LY86 
    ##              1.750              1.750              1.750              1.750 
    ##              TUSC5                GAL             CLDN14              CCL11 
    ##              1.750              1.750              1.751              1.751 
    ##               MYL2              PGBD3      RP11-1021N1.1            ARHGAP8 
    ##              1.751              1.751              1.751              1.751 
    ##              OR5J2              CCL13              NPY4R       CTD-2370N5.3 
    ##              1.751              1.752              1.752              1.752 
    ##              OR6B1           TMEM176B           C1orf115              CCL14 
    ##              1.752              1.753              1.753              1.753 
    ##            FAM166B                SNN               IL33              OR5R1 
    ##              1.753              1.753              1.753              1.753 
    ##              CCL26              TTC23             MMACHC             ZNF284 
    ##              1.754              1.754              1.755              1.755 
    ##              ZNF30              TRNP1           C11orf24              PRRG3 
    ##              1.755              1.755              1.756              1.756 
    ##             RGS9BP               OIP5           C1orf122             ZNF714 
    ##              1.756              1.756              1.756              1.756 
    ##            FAM109A             MYL12A              NCMAP            B3GALT5 
    ##              1.756              1.756              1.756              1.757 
    ##              PINX1             FBXO48            TMEM196             EFCAB9 
    ##              1.757              1.757              1.757              1.757 
    ##               G0S2             COX6A1           C11orf44           C11orf91 
    ##              1.757              1.757              1.757              1.757 
    ##      RP11-343C2.12           C17orf89             CALML6           C6orf223 
    ##              1.757              1.757              1.757              1.757 
    ##              APOC3               IGIP              BLVRB             SHISA8 
    ##              1.757              1.757              1.758              1.758 
    ##            TMEM129              HEBP1         AL590708.2          SPATA31A1 
    ##              1.758              1.758              1.758              1.758 
    ##               PLD6         AC110615.1               GIPR           C22orf46 
    ##              1.758              1.758              1.758              1.758 
    ##              PHGR1         AC002553.1              ZNF43             ZNF254 
    ##              1.759              1.759              1.759              1.759 
    ##              TRIQK             GEMIN6           DEFB105A         AC135048.1 
    ##              1.759              1.759              1.759              1.760 
    ##             FBXO27             ZNF610         AC023590.1               CSH1 
    ##              1.760              1.760              1.760              1.760 
    ##             TREML4             TRIM61             GIMAP6               DDC8 
    ##              1.760              1.760              1.760              1.760 
    ##           SLC25A52              CRIP1             MS4A14              OR2G6 
    ##              1.760              1.761              1.761              1.761 
    ##               MMP8            SPDYE2B               PKIA             OR4F17 
    ##              1.761              1.761              1.761              1.761 
    ##             UQCR10            PGLYRP4             RNF224            PRKRIP1 
    ##              1.761              1.761              1.761              1.761 
    ##               NPVF             ATPIF1             PKD1L2              CNPY4 
    ##              1.762              1.763              1.763              1.764 
    ##            METTL20                LHB              TOMM7             DLEU2L 
    ##              1.764              1.764              1.764              1.764 
    ##            C4orf46              WFDC9               SCO2             LY6G6C 
    ##              1.764              1.764              1.764              1.764 
    ##            A3GALT2               ST20       RP1-228P16.5             TRIM73 
    ##              1.764              1.764              1.765              1.765 
    ##               SNTN           C1orf194          KRTAP29-1               CLNK 
    ##              1.765              1.765              1.765              1.766 
    ##              PGAM2             COMMD8              LPAR4             ZNF443 
    ##              1.766              1.766              1.766              1.766 
    ##            C1orf68               LSM3              TYW1B             BPIFA3 
    ##              1.766              1.766              1.766              1.767 
    ##       RP11-849H4.2              OR2W3           CDC42EP2             LSMEM2 
    ##              1.767              1.767              1.767              1.767 
    ##              CTXN3             MPV17L         AL031590.1         AC018867.2 
    ##              1.767              1.767              1.767              1.768 
    ##           C19orf69           HNRNPCL1           TNFRSF6B             TEDDM1 
    ##              1.768              1.768              1.769              1.769 
    ##            C8orf22            LGALS16             MTNR1A               CLPS 
    ##              1.769              1.769              1.769              1.770 
    ##            CYSLTR2            METTL23              MRP63           C22orf39 
    ##              1.770              1.770              1.770              1.770 
    ##              SPRR4           AURKAIP1             OGFOD2              OR4K5 
    ##              1.770              1.770              1.770              1.771 
    ##              HOGA1            C9orf47            SPINK13            LAMTOR5 
    ##              1.771              1.771              1.771              1.771 
    ##             KCTD14             FBXO17                AMH            NDUFA13 
    ##              1.771              1.771              1.771              1.772 
    ##               NIT2              ACSF3           C22orf26              OR1Q1 
    ##              1.772              1.772              1.772              1.772 
    ##               GZMH              SIRT4             ZNF728               PCP4 
    ##              1.773              1.773              1.773              1.773 
    ##             FAM26F              REG3G          C10orf131             S100A8 
    ##              1.773              1.773              1.774              1.774 
    ##            MAGEB16      RP11-1085N6.3            RIPPLY1              PDHA2 
    ##              1.774              1.774              1.774              1.774 
    ##           C22orf43          TNFAIP8L2           C19orf35       RP11-712L6.5 
    ##              1.775              1.775              1.775              1.775 
    ##             TAS2R3      RP11-187E13.1       RP11-321F6.1           C19orf45 
    ##              1.775              1.776              1.776              1.776 
    ##              BOLA1             CLDND2              ORAI1               TCTA 
    ##              1.777              1.777              1.777              1.777 
    ##         AC136604.1              REP15               MIOX              OR2H2 
    ##              1.777              1.777              1.777              1.777 
    ##             LRRC19              SYT15              CCRL2            CXorf28 
    ##              1.777              1.777              1.777              1.778 
    ##              NRN1L             OR5D18              ASB13             CYP1A1 
    ##              1.778              1.778              1.778              1.778 
    ##        CTB-186H2.3               UCN2               RFNG              INSL5 
    ##              1.778              1.779              1.779              1.779 
    ##              TEX12               KPRP             CCDC70                CPO 
    ##              1.779              1.780              1.780              1.780 
    ##              COX16          KRTAP26-1         AC010336.1            CCDC115 
    ##              1.780              1.780              1.780              1.780 
    ##              AGAP7               BSND             OPN1MW          TNFAIP8L3 
    ##              1.780              1.780              1.780              1.780 
    ##          KRTAP10-9             OR52L1              UPK3A             IFITM1 
    ##              1.780              1.780              1.781              1.781 
    ##               FRG2      RP11-571M6.15       RP11-10J21.3              CCL18 
    ##              1.781              1.781              1.781              1.782 
    ##            KIR3DL1               FGF3               MRAP               TNP1 
    ##              1.782              1.782              1.782              1.782 
    ##               TNP2             NPBWR2               POP7              OR9G1 
    ##              1.782              1.783              1.783              1.783 
    ##             ZNF594               MGMT             ZNF556              TLDC1 
    ##              1.783              1.783              1.783              1.784 
    ##            C8orf12              SCN4B            CCDC179              ZBED2 
    ##              1.784              1.784              1.784              1.784 
    ##           C20orf85         AL049829.1             MRPL14             FAM86A 
    ##              1.784              1.784              1.784              1.785 
    ##              KHDC1              HINT3             PROKR1            PROSER2 
    ##              1.785              1.785              1.785              1.785 
    ##               XCL2             ZNF727              CYYR1             ZNF844 
    ##              1.785              1.785              1.785              1.785 
    ##           C22orf42            ATP5EP2        KB-1507C5.2            PLA2G2D 
    ##              1.785              1.786              1.786              1.786 
    ##              FCRL6          HTR5A-AS1            ZCCHC16               HCST 
    ##              1.786              1.786              1.786              1.786 
    ##         AC103809.2            C7orf50             MS4A12            PGLYRP1 
    ##              1.786              1.787              1.787              1.787 
    ##            C5orf28                LYZ              DAND5              AP5S1 
    ##              1.787              1.787              1.787              1.787 
    ##               UCP1            FAM109B             ZNF479            ZSCAN5D 
    ##              1.788              1.788              1.788              1.788 
    ##              APOL3        RP6-24A23.6             CCDC78            SULT1A3 
    ##              1.788              1.788              1.788              1.789 
    ##               NMBR          ANKRD20A1              LYZL2            C2orf53 
    ##              1.789              1.789              1.789              1.790 
    ##             AQP12B            SMPDL3A              ATP5E              SMIM5 
    ##              1.790              1.790              1.790              1.790 
    ##              LEMD1      CTD-2014B16.3             ZNF586             GTSF1L 
    ##              1.790              1.791              1.791              1.791 
    ##               RETN               CDNF              OR2A5              ADIRF 
    ##              1.791              1.791              1.791              1.791 
    ##           C12orf23            GAL3ST3            TMEM155              SFTA2 
    ##              1.791              1.792              1.792              1.792 
    ##         AC111200.1             B3GNT9              MPV17              TREM2 
    ##              1.792              1.792              1.792              1.793 
    ##              ECSCR           TMEM167A               NAT2                HRK 
    ##              1.793              1.793              1.793              1.793 
    ##            WFDC10B               CGB1               CAPS              CDCP2 
    ##              1.793              1.793              1.793              1.794 
    ##               PRM3             PDGFRL            RNASE10              ENPP7 
    ##              1.794              1.795              1.795              1.795 
    ##              TPRX1             RNF181            B3GNTL1          C17orf100 
    ##              1.795              1.795              1.795              1.795 
    ##            C3orf52              RDH14            ANKRD66             C8orf4 
    ##              1.795              1.795              1.795              1.796 
    ##             AGAP10       RP13-672B3.2              KRT32            FAM174A 
    ##              1.796              1.796              1.796              1.796 
    ##         AC006538.4              NAT14           HIST1H4F              OR5T3 
    ##              1.796              1.796              1.796              1.796 
    ##           C9orf169             OR10A4            TMEM212             FAM72C 
    ##              1.797              1.797              1.797              1.797 
    ##            LAPTM4B   LL22NC03-75H12.2              PAGE3              KLRK1 
    ##              1.797              1.797              1.797              1.797 
    ##             TIMM8B             ZNF665            SLC31A2             C2CD4D 
    ##              1.797              1.797              1.797              1.798 
    ##               CER1             PCDHB7              GPR18              POTED 
    ##              1.798              1.798              1.798              1.798 
    ##             ZNF737          TNFRSF10C         AP002956.1         AC064874.1 
    ##              1.798              1.798              1.798              1.799 
    ##              COX14                PIR              OR4N2           C10orf67 
    ##              1.799              1.799              1.799              1.799 
    ##              MCCD1              INCA1               GYPB              PSG11 
    ##              1.799              1.799              1.800              1.800 
    ##             PRSS58             ZNF487              AAED1            PPP1R27 
    ##              1.801              1.801              1.801              1.801 
    ##              G6PC2              OR5M1              CTXN2              OR1B1 
    ##              1.801              1.801              1.802              1.802 
    ##             OR51T1            S100A11      RP11-113D6.10               ACP1 
    ##              1.802              1.802              1.802              1.802 
    ##                OMP           HIST1H3I              OR6M1             ZNHIT3 
    ##              1.802              1.802              1.802              1.803 
    ##        AP001631.10              ACKR4             POU5F2       RP1-32I10.10 
    ##              1.803              1.803              1.803              1.803 
    ##            C7orf33         AL035588.1            C3orf14              CCL22 
    ##              1.803              1.803              1.804              1.804 
    ##               CBR3         AC005008.2             GAGE2C               TSHB 
    ##              1.804              1.804              1.804              1.805 
    ##           C16orf95            SPATC1L              CLRN1            RNASE12 
    ##              1.805              1.805              1.805              1.805 
    ##              NUDT8             HIGD1C            S100A14            SLC10A5 
    ##              1.805              1.806              1.806              1.806 
    ##              PSMG4            ZNF587B            C6orf99          KRTAP21-1 
    ##              1.806              1.806              1.806              1.806 
    ##              OR5K3              LAGE3            FAM154B              POTEB 
    ##              1.806              1.806              1.806              1.807 
    ##      RP11-315D16.2              KLRC2         CKLF-CMTM1            TMEM258 
    ##              1.807              1.807              1.807              1.807 
    ##            GOLGA8A            PRKCDBP         AL031666.2              PDE6H 
    ##              1.807              1.807              1.807              1.807 
    ##         AC011308.1               RAX2            SPANXN4            TIMM23B 
    ##              1.807              1.807              1.807              1.808 
    ##          NIPSNAP3A              STATH               AMTN          C10orf126 
    ##              1.808              1.808              1.808              1.808 
    ##            BLOC1S5         AC025287.1               SDSL               OPTC 
    ##              1.808              1.808              1.808              1.809 
    ##              ZNF93             OR51Q1              CST11              JMJD8 
    ##              1.809              1.809              1.809              1.809 
    ##              RDH13       RP11-113D6.6           C12orf36          C20orf202 
    ##              1.809              1.809              1.809              1.809 
    ##          KRTAP10-2              HINT1              GPR39            PRAMEF2 
    ##              1.810              1.810              1.810              1.810 
    ##               HAMP             TCEAL2           C19orf40               SFR1 
    ##              1.810              1.810              1.811              1.811 
    ##         AF165138.7            TBC1D28           SLC22A31             HSF2BP 
    ##              1.811              1.811              1.811              1.812 
    ##            TMEM223             SDHAF2          C17orf102         AL354993.1 
    ##              1.812              1.812              1.812              1.812 
    ##               TMIE          PABPC1L2A             ZNF547           KRTAP4-6 
    ##              1.812              1.812              1.813              1.813 
    ##         AC131097.4             COMMD5              GINS2           C2orf27A 
    ##              1.813              1.813              1.813              1.813 
    ##              RSAD2         AC003002.4           KRTAP5-2            TMEM45B 
    ##              1.813              1.813              1.814              1.814 
    ##              OR1N2              TAAR2               BRI3           C21orf37 
    ##              1.814              1.814              1.814              1.814 
    ##              OR5T2         AC106873.4           C1orf185               FCN1 
    ##              1.814              1.815              1.815              1.815 
    ##             MANBAL              OTOL1             LRRC66            DNAJC24 
    ##              1.815              1.815              1.815              1.815 
    ##             ZNF616               CCL1           CDRT15L2       CTD-2600O9.1 
    ##              1.815              1.816              1.816              1.816 
    ##          LINC01098              PQLC3      RP11-167N24.6      RP11-215A19.2 
    ##              1.816              1.816              1.816              1.816 
    ##         AC068987.1             MRPL20          C20ORF135             OR10Q1 
    ##              1.816              1.817              1.817              1.817 
    ##            PLA2G2F            C7orf71          C14orf182             TRIM49 
    ##              1.817              1.817              1.817              1.817 
    ##             OR4K17             S100A3       RP4-576H24.4               CCL5 
    ##              1.817              1.817              1.818              1.818 
    ##              CFHR2       RP11-126K1.2           C12orf68            C8orf59 
    ##              1.818              1.818              1.818              1.818 
    ##              CRYGA            TMEM114             CITED4              GPR78 
    ##              1.818              1.818              1.818              1.818 
    ##              S100P               SSX7         AP000708.1             OR12D3 
    ##              1.818              1.818              1.819              1.819 
    ##             IL18BP            TSPAN19       RP11-17M16.1              ASCL2 
    ##              1.819              1.819              1.819              1.819 
    ##                CEA             CACNG1             FAM69C               PTER 
    ##              1.819              1.819              1.819              1.819 
    ##             COMMD2            FAM187B             ZNF559           C19orf84 
    ##              1.819              1.820              1.820              1.820 
    ##              OR1G1           C16orf98           C1orf162              SMIM1 
    ##              1.820              1.820              1.820              1.820 
    ##             MAGEE2            GAGE12H               WBP5             FAM72A 
    ##              1.820              1.820              1.821              1.821 
    ##             GPR111                VCX             LSMEM1             STARD4 
    ##              1.821              1.821              1.821              1.822 
    ##             TXNL4A         AC017028.1             P2RY11             PET100 
    ##              1.822              1.823              1.823              1.823 
    ##              CHST6      RP11-386G21.2         AL035406.1              ZFP41 
    ##              1.823              1.823              1.823              1.823 
    ##              ZFP41             MS4A10    RP11-1396O13.13              CRYGB 
    ##              1.823              1.823              1.823              1.824 
    ##           C21orf54          C14orf142              NACA2           C1orf234 
    ##              1.824              1.824              1.824              1.825 
    ##     CTD-2583A14.10             OR5AS1            AKR1CL1               GJB5 
    ##              1.825              1.825              1.825              1.825 
    ##           C11orf83          HIST1H2AE             CCL3L1               MT1G 
    ##              1.825              1.825              1.825              1.826 
    ##             MUSTN1           C9orf153            EIF5AL1              CLLU1 
    ##              1.826              1.826              1.826              1.826 
    ##            TMEM218           TP53TG3D              CLDN4             NBPF12 
    ##              1.826              1.826              1.827              1.827 
    ##             ZNF519               CIB3           C19orf82             PET117 
    ##              1.827              1.827              1.827              1.827 
    ##            C7orf69            UBE2Q2L           C1orf233             GAGE10 
    ##              1.827              1.828              1.828              1.828 
    ##             CYB5D1            SULT2A1         AL162407.1            SLC34A1 
    ##              1.828              1.828              1.828              1.828 
    ##             BTN2A1             ZNF670           EIF4EBP3              CFHR5 
    ##              1.829              1.829              1.829              1.829 
    ##             PRSS48              TCTE3           C22orf23       RP11-770J1.5 
    ##              1.829              1.829              1.829              1.829 
    ##               LCN1              OR1J1       RP11-219B4.6             PSMB11 
    ##              1.830              1.830              1.830              1.830 
    ##              CTRB1           SIGLECL1               APCS         AC013468.1 
    ##              1.830              1.830              1.830              1.830 
    ##             CELA2A              OR2A4            FAM174B              ASCL4 
    ##              1.830              1.830              1.830              1.831 
    ##              ZNF28             ZNF576             KCNE1L            HRASLS2 
    ##              1.831              1.831              1.831              1.831 
    ##             CD3EAP           C15orf56             ZNF530         AC007956.1 
    ##              1.831              1.831              1.831              1.831 
    ##              PNRC2               MNDA           C1orf204              SMIM9 
    ##              1.832              1.832              1.832              1.832 
    ##           C16orf91      CTD-2210P24.4             TRIM40      RP11-463J10.2 
    ##              1.832              1.832              1.832              1.832 
    ##                ID4               MT1E              PGBD4         AL136218.1 
    ##              1.833              1.833              1.833              1.833 
    ##         AL138815.1         AL929472.1          HIST1H2BJ             OR5M11 
    ##              1.833              1.833              1.833              1.833 
    ##             RNASE9             CDKN2B             OR51B6             OR13F1 
    ##              1.833              1.833              1.833              1.833 
    ##               LYG2              VCX3B             CARD17                MR1 
    ##              1.833              1.833              1.833              1.833 
    ##       RP11-77K12.1            C2orf66             MMP23B           C10orf10 
    ##              1.833              1.834              1.834              1.834 
    ##             CLDN17               TAF9              APOL2              SSNA1 
    ##              1.834              1.834              1.834              1.834 
    ##              PCBD1            C1orf63           C10orf91             ARL17B 
    ##              1.835              1.835              1.835              1.835 
    ##             NDUFA7             CCL4L1              DGCR6              HHLA3 
    ##              1.835              1.835              1.835              1.835 
    ##              MAGIX         AC112205.1              ACYP2               PSG6 
    ##              1.835              1.835              1.836              1.836 
    ##             OR51F2             RNASE7            C5orf64              COX8A 
    ##              1.836              1.837              1.837              1.837 
    ##              MEIG1              CHAC2             TREML1               MC5R 
    ##              1.837              1.837              1.837              1.837 
    ##               DDTL              OR4S1            C8orf33             OR2T33 
    ##              1.837              1.837              1.837              1.837 
    ##            ANXA8L1               ECI1       RP11-830F9.6         AC003006.7 
    ##              1.837              1.837              1.837              1.838 
    ##              PHPT1              AAMDC             ATP2C2           C17orf64 
    ##              1.838              1.838              1.838              1.838 
    ##              LCE2A             SCP2D1              OR1L1          C14orf177 
    ##              1.838              1.838              1.839              1.839 
    ##             IFITM3              AMY1A               PSG3         AC011500.1 
    ##              1.839              1.839              1.839              1.839 
    ##            SCGB3A1             SPINK8         AC009802.1              OR4F5 
    ##              1.839              1.839              1.840              1.840 
    ##              CCL24              MTMR8             ZNF730              CRYGN 
    ##              1.840              1.840              1.840              1.840 
    ##             GDPGP1               CGB2       RP1-127H14.3         AC015987.2 
    ##              1.840              1.840              1.840              1.840 
    ##         AC015989.1             CLEC2A           DEFB104B      DKFZP779J2370 
    ##              1.840              1.841              1.841              1.841 
    ##       RP11-368I7.4              S100G              RPAIN             VCPKMT 
    ##              1.841              1.841              1.841              1.841 
    ##                UBD       RP11-80A15.1         RP11-9B6.1              ZNF92 
    ##              1.841              1.841              1.841              1.841 
    ##             GPR148            PRAMEF7             OR10R2             IGFALS 
    ##              1.842              1.842              1.842              1.842 
    ##            DYNLRB2              UBL4B             ZNF736              LCE1B 
    ##              1.843              1.843              1.843              1.843 
    ##             TRIM65            FAM86B1         AC096644.1           KRTAP2-3 
    ##              1.843              1.844              1.844              1.844 
    ##             FAM72D           PRAMEF15           C1orf200              LRP5L 
    ##              1.844              1.844              1.844              1.844 
    ##                DDT         AC006486.1             GAGE13            DEFB124 
    ##              1.844              1.844              1.844              1.844 
    ##             RNF183               VMO1              MMP28                ID1 
    ##              1.844              1.845              1.845              1.845 
    ##              HSPB3             UPK3BL              LCE5A             TMEM61 
    ##              1.845              1.845              1.845              1.845 
    ##           SLC22A16             CHURC1              OR4F6               MC2R 
    ##              1.846              1.846              1.846              1.846 
    ##          HIST1H2AI       RP11-10A14.4           C1orf145            CXorf66 
    ##              1.846              1.846              1.847              1.847 
    ##              PEX12           POM121L2             TIMM10            C11orf1 
    ##              1.847              1.847              1.847              1.847 
    ##             ZNF695          TRAPPC2P1             OR10X1         AC003101.1 
    ##              1.847              1.847              1.847              1.847 
    ##              CXCL3       RP11-664D7.4             PMAIP1            OR56B3P 
    ##              1.848              1.848              1.848              1.848 
    ##            TMEM221           C1orf134            CXorf31           DEFB107A 
    ##              1.848              1.848              1.848              1.848 
    ##              ULBP2       CTC-236F12.4              VTCN1             GNG5P2 
    ##              1.848              1.848              1.848              1.849 
    ##       RP11-45H22.3       RP11-738G5.2              ZG16B             ZNF550 
    ##              1.849              1.849              1.849              1.849 
    ##             OR5AN1              OR6P1             CLDN24             FBXL22 
    ##              1.849              1.849              1.849              1.849 
    ##           HLA-DRB5              RPS17               SVIP             ZNF417 
    ##              1.849              1.849              1.849              1.849 
    ##               PSG9            FAM150A              VSTM5           C9orf129 
    ##              1.849              1.849              1.849              1.850 
    ##           C21orf90            HLA-DOA              TPSD1         AC012313.1 
    ##              1.850              1.850              1.850              1.850 
    ##             IL1F10             ZNF114              FOXI2              KRT6C 
    ##              1.850              1.850              1.850              1.850 
    ##               HBA1             NMNAT3         AC104667.3               IQCD 
    ##              1.850              1.850              1.850              1.851 
    ##              ARGFX               TLR1         AC018816.3              NXPE1 
    ##              1.851              1.851              1.851              1.851 
    ##              RFPL2               SSX5              MFSD3             MAGED4 
    ##              1.851              1.851              1.852              1.852 
    ##               HCRT             OR10G7      RP11-131H24.4            FAM154A 
    ##              1.852              1.852              1.852              1.852 
    ##             OR10A3            BHLHA15         AL359736.1           HIST1H3F 
    ##              1.852              1.853              1.853              1.853 
    ##              ULBP3            C1orf61              CYTL1              OR5K2 
    ##              1.853              1.853              1.853              1.853 
    ##             NPIPA1              PATE4            POU5F1B            MAGEB18 
    ##              1.853              1.853              1.853              1.853 
    ##            TMEM262              RCVRN              COX8C             CARD18 
    ##              1.853              1.854              1.854              1.854 
    ##           KRTAP4-4               CSTB               CYS1       FAM47E-STBD1 
    ##              1.854              1.854              1.854              1.854 
    ##            C1orf53            C1orf64          PABPC1L2B              SLMO1 
    ##              1.855              1.855              1.855              1.855 
    ##               NT5M               RLN1           CDC42EP5             ZNF208 
    ##              1.855              1.855              1.855              1.855 
    ##             ZNF776             NDUFB2     RP11-834C11.12              OTUB2 
    ##              1.855              1.855              1.855              1.856 
    ##             ZNF708             ZNF709           SERPINB3            C7orf65 
    ##              1.856              1.856              1.856              1.856 
    ##           C12orf42            ZCCHC13           C15orf48         AC121757.1 
    ##              1.856              1.856              1.856              1.856 
    ##            CLLU1OS             OR13C8             LRRC30             FAM72B 
    ##              1.856              1.856              1.856              1.857 
    ##         AC023469.1             NPIPA5              LUZP4              HMGN4 
    ##              1.857              1.857              1.857              1.857 
    ##               GNLY             TMSB4Y             ZNF347             OR10H1 
    ##              1.857              1.858              1.858              1.858 
    ##              LCE1C            ANAPC13              OR9G4               PRLH 
    ##              1.858              1.858              1.858              1.858 
    ##             RNF222             ZCCHC5      RP11-625H11.1            SPAG11B 
    ##              1.858              1.858              1.858              1.859 
    ##           C1orf213             KLHDC4         AC090616.2               DRD4 
    ##              1.859              1.859              1.860              1.860 
    ##           C19orf71              OR2L8               MT1H            PLA2G16 
    ##              1.860              1.860              1.860              1.860 
    ##         AC120194.1         AC003043.1             SMIM10            IL22RA2 
    ##              1.860              1.860              1.861              1.861 
    ##             ZNF223             AKR1C3            SDR42E1            RNASE13 
    ##              1.861              1.861              1.861              1.861 
    ##                IL8              OR5H6          C10orf113             HILPDA 
    ##              1.861              1.861              1.861              1.861 
    ##           HIST1H1C            FAM110C              OR2T7               NPPB 
    ##              1.861              1.862              1.862              1.862 
    ##              MRPS2      RP11-1026M7.2       CTD-2021H9.3             OR52N5 
    ##              1.862              1.862              1.862              1.862 
    ##             OR7A17      RP11-366L20.2       RP11-404L6.2             CCL3L3 
    ##              1.862              1.862              1.862              1.863 
    ##            TSPAN10           KRTAP9-2             OR2AJ1               TDRP 
    ##              1.863              1.863              1.863              1.863 
    ##              PRSS1            CCDC74A              ERP27           C11orf31 
    ##              1.863              1.863              1.863              1.863 
    ##           C12orf57               ENHO             RAET1L              FABP4 
    ##              1.863              1.864              1.864              1.864 
    ##             ZNF709              MPZL2             NPIPB9               NNMT 
    ##              1.864              1.864              1.864              1.864 
    ##             S100A4              UFSP1      RP11-872D17.8            RAB40AL 
    ##              1.864              1.864              1.864              1.865 
    ##           PSORS1C1             ZNF578         AL592284.1         AC140481.2 
    ##              1.865              1.865              1.865              1.865 
    ##      RP11-1105G2.3              GDF15              OR4N4              CTAG2 
    ##              1.865              1.865              1.865              1.865 
    ##      DNAJC25-GNG10            FAM159A              CXCR1      RP11-826N14.2 
    ##              1.865              1.865              1.865              1.865 
    ##            SCGB3A2           GOLGA6L4               PRM1             NDUFV3 
    ##              1.866              1.866              1.866              1.866 
    ##              ZNF69           TRAPPC6A              XAGE5               XCL1 
    ##              1.866              1.866              1.866              1.866 
    ##            COMMD10         AC009041.2      RP11-293M10.1       CTD-3203P2.2 
    ##              1.866              1.866              1.867              1.867 
    ##              ARL16              OR8I2            SPACA5B             DNAJB7 
    ##              1.867              1.867              1.867              1.867 
    ##             SPRR1A     RP11-1102P16.1         AC110781.3              OR2A1 
    ##              1.867              1.867              1.867              1.867 
    ##           KRTAP5-6             RNF223            CARHSP1             TMEM31 
    ##              1.867              1.867              1.867              1.867 
    ##          C10orf115             S100A1         AC017104.2         AL162431.1 
    ##              1.868              1.868              1.868              1.868 
    ##      RP4-583P15.14             ZNF766              LYRM2              IFNL2 
    ##              1.868              1.868              1.868              1.868 
    ##               SMN2            C7orf34               IL9R      RP11-683L23.1 
    ##              1.868              1.868              1.868              1.868 
    ##          HIST1H2BA              LEAP2         AL390778.1              THEM5 
    ##              1.869              1.869              1.869              1.869 
    ##             KCNMB3              RXFP3            ZSCAN5B              TSPY1 
    ##              1.869              1.869              1.869              1.869 
    ##         AC005082.1              VCX3A           KRTAP8-1              PRR24 
    ##              1.869              1.870              1.870              1.870 
    ##           SIGLEC15            C2ORF15              USP50         AP001024.2 
    ##              1.870              1.870              1.870              1.870 
    ##             OR10J3             OR10S1            TMEM211            KIR3DL2 
    ##              1.870              1.870              1.870              1.870 
    ##          KRTAP13-4              OR4D5            KBTBD13               SELM 
    ##              1.870              1.870              1.870              1.870 
    ##            CYP2C19         AL365202.1            CYP4A22               LY6D 
    ##              1.871              1.871              1.871              1.871 
    ##            ZNF321P              GRAPL            FAM163B            PABPN1L 
    ##              1.871              1.871              1.871              1.872 
    ##              TTLL2         AC073569.1              HSPB1         AC005549.3 
    ##              1.872              1.872              1.872              1.872 
    ##              OR4F4             CGREF1         AC104057.1            BLOC1S4 
    ##              1.872              1.872              1.872              1.872 
    ##            TMEM234            GOLGA8B         CTB-54O9.9             GAGE2A 
    ##              1.872              1.872              1.873              1.873 
    ##             MRPL21         AL136376.1              OR6V1              ZNF85 
    ##              1.873              1.873              1.873              1.873 
    ##              VEGFB               GAPT              TCF24              RAB19 
    ##              1.873              1.873              1.873              1.873 
    ##             SFT2D3               DNLZ           C15orf54              OR5B3 
    ##              1.873              1.873              1.874              1.874 
    ##              OR5H1            ZNF780A            TMEM171               PGA4 
    ##              1.874              1.874              1.874              1.874 
    ##             PLGRKT             RNASE1            FAM90A1              OR1S2 
    ##              1.874              1.874              1.874              1.874 
    ##         AL590714.1           HIST1H3H            ANKRD33              CCL20 
    ##              1.875              1.875              1.875              1.875 
    ##             PAGE2B              P2RY6             MRPS28             MRPL11 
    ##              1.875              1.875              1.875              1.875 
    ##               CGB5             CT45A6           PRAMEF18               CMA1 
    ##              1.875              1.875              1.876              1.876 
    ##          PAPPA-AS1           SLC39A11            ZNF705B         AC108868.1 
    ##              1.876              1.876              1.876              1.876 
    ##             CLEC7A               FGF5           C6orf226             OR2T29 
    ##              1.876              1.876              1.876              1.876 
    ##             PLGLB2          C20orf201            TMEM257             DDIT4L 
    ##              1.876              1.876              1.876              1.876 
    ##      CTD-2547L24.3            TMEM134           KRTAP4-1               EQTN 
    ##              1.876              1.877              1.877              1.877 
    ##         AL163636.6              FFAR3              LAIR2              MRAP2 
    ##              1.877              1.877              1.877              1.878 
    ##                NPB              OR5H2                AEN            KIR2DL1 
    ##              1.878              1.878              1.878              1.878 
    ##            ANKRD23        CTB-60B18.6             GUCA1C            GOLGA8J 
    ##              1.878              1.878              1.878              1.878 
    ##                CGB            PRAMEF1               SYCN              OR8H2 
    ##              1.879              1.879              1.879              1.879 
    ##           SERPINB4           KRTAP9-7              SMIM8             POLR2K 
    ##              1.879              1.879              1.879              1.879 
    ##               PMM2           C12orf73            FAM111B          GOLGA6L10 
    ##              1.879              1.879              1.879              1.879 
    ##           TNFRSF18             NBPF24              LCE2B              OR8S1 
    ##              1.880              1.880              1.880              1.880 
    ##             ZNF256             CT45A1            ZBED6CL             SERF1B 
    ##              1.880              1.880              1.880              1.880 
    ##             PRSS55             OR56B4             OR4A16             ZNF671 
    ##              1.880              1.881              1.881              1.881 
    ##              MOAP1              IL3RA              DUPD1             SEPT12 
    ##              1.881              1.881              1.881              1.881 
    ##            FAM162B               PIGY             MANSC4             SPIN2A 
    ##              1.881              1.881              1.881              1.881 
    ##            DEFB118            DEFB119            DEFB129              DAPL1 
    ##              1.882              1.882              1.882              1.882 
    ##            C2orf57              LCNL1           PRAMEF16           CXorf40B 
    ##              1.882              1.882              1.882              1.882 
    ##             OR11G2         AC187652.1             KBTBD4             OR10A6 
    ##              1.882              1.883              1.883              1.883 
    ##             VPREB1           HIST1H3D             CYP2C8             GIMAP5 
    ##              1.883              1.883              1.883              1.883 
    ##              OR6K2            DEFB121         AP000688.1             ZNF716 
    ##              1.883              1.883              1.883              1.883 
    ##               LST1              OTOP1            PP13004             CCDC23 
    ##              1.883              1.883              1.883              1.883 
    ##            TMEM261               MT1B             SYNGR4              OR2D2 
    ##              1.883              1.884              1.884              1.884 
    ##           C9orf141         AC090427.1           C1orf143             RBMY1F 
    ##              1.884              1.884              1.884              1.884 
    ##              CDRT4              APOL6           TNFRSF17       CTB-134H23.2 
    ##              1.884              1.884              1.884              1.884 
    ##               EMP2            DEFB130              OR1S1          SPATA31A3 
    ##              1.884              1.885              1.885              1.885 
    ##              TSPY2             NBPF11             BHLHA9            C9orf66 
    ##              1.885              1.885              1.885              1.885 
    ##               LRG1         AC112721.1             B3GNT6             MBD3L4 
    ##              1.885              1.885              1.885              1.885 
    ##            GAL3ST2               HTN3             RHOXF2             DMRTC1 
    ##              1.885              1.886              1.886              1.886 
    ##              ARL11            NEUROG1               CGB7          HIST1H2AJ 
    ##              1.886              1.886              1.886              1.886 
    ##      RP11-210M15.2            C6orf48           KRTAP5-8              OR2T6 
    ##              1.886              1.886              1.886              1.886 
    ##               NAT1         AL161784.1              WFDC6           SERPINA9 
    ##              1.886              1.886              1.886              1.886 
    ##           C1orf111             ZNF721            SPG20OS      RP11-386G21.1 
    ##              1.886              1.886              1.886              1.886 
    ##              SIRPD              TEX38             ZNF726             VPREB3 
    ##              1.886              1.887              1.887              1.887 
    ##               SPRN           C10orf32            C4orf19              CBWD3 
    ##              1.887              1.887              1.887              1.887 
    ##         AC022400.2              TSPY8           C1orf180            TMEM238 
    ##              1.887              1.887              1.887              1.887 
    ##               CST9            C8orf31            DEFB113           ANKRD34B 
    ##              1.887              1.887              1.888              1.888 
    ##              TRIB3           KRTAP3-1           DEFB104A         AC009403.2 
    ##              1.888              1.888              1.888              1.888 
    ##         AC012360.1         AL161450.1             ZNF429            DNAJC15 
    ##              1.888              1.888              1.888              1.889 
    ##             RPUSD1              OR1J2              OR4C6           PRAMEF13 
    ##              1.889              1.889              1.889              1.889 
    ##               TFF2             CDRT15           DEFB105B              TYRP1 
    ##              1.889              1.889              1.889              1.889 
    ##               URAD              OVCA2              CLIC3             IFITM2 
    ##              1.889              1.889              1.889              1.889 
    ##             MBD3L3      CTD-2215E18.1             ZNF534             PPP3R2 
    ##              1.889              1.889              1.889              1.889 
    ##      RP11-998D10.1           PRAMEF14              OR9Q2              HCG27 
    ##              1.889              1.890              1.890              1.890 
    ##             MS4A4E            SCGB2B2            S100A16           TRIM49D1 
    ##              1.890              1.890              1.890              1.890 
    ##              OR6K3            RBMY1A1             RBMY1E            C7orf66 
    ##              1.890              1.890              1.890              1.891 
    ##             HIGD1B           KRTAP9-9              MGARP              OR4K2 
    ##              1.891              1.891              1.891              1.891 
    ##             ZNF846            DEFB136               KLLN         AL121963.1 
    ##              1.891              1.891              1.891              1.891 
    ##             OR5AU1              ARL14               GPX1             MRPS12 
    ##              1.891              1.892              1.892              1.892 
    ##            DEFB114               DAZ2               SPP2             SNAPC5 
    ##              1.892              1.892              1.892              1.893 
    ##                ANG              KAAG1              KCNG4             OR4C12 
    ##              1.893              1.893              1.893              1.893 
    ##               CST5                SLN              AVPI1              SPON2 
    ##              1.893              1.893              1.893              1.893 
    ##             FGFBP3              OR8K3              LCE1F             RBMY1J 
    ##              1.893              1.893              1.893              1.893 
    ##           FAM90A26           KRTAP5-7            SPANXN2            PRAMEF3 
    ##              1.894              1.894              1.894              1.894 
    ##             ZNF486              DEFA3       RP11-763F8.1         AC026369.1 
    ##              1.894              1.894              1.894              1.894 
    ##             OR10H4             COX7A1             PYCARD            RNASE11 
    ##              1.895              1.895              1.895              1.895 
    ##             OTUD6A             CLDN20         AC009065.1             CHCHD1 
    ##              1.895              1.895              1.895              1.895 
    ##             IL36RN         AC011239.1            HIST4H4               NPFF 
    ##              1.895              1.895              1.895              1.896 
    ##               IFNK           C12orf74           KRTAP9-6             TSPY10 
    ##              1.896              1.896              1.896              1.896 
    ##              FOXD4               SKA2            C8orf82           C11orf86 
    ##              1.896              1.896              1.897              1.897 
    ##            CENPBD1          LINC00632         AC019206.1         AL355390.1 
    ##              1.897              1.897              1.897              1.897 
    ##            DEFB133             NDUFB1           C17orf67             SPACA4 
    ##              1.897              1.897              1.897              1.897 
    ##            COX7A2L             KLHL35              TPSB2          C10orf120 
    ##              1.897              1.898              1.898              1.898 
    ##             CCL4L2      RP11-404P21.8             FAM25C               PSG5 
    ##              1.898              1.898              1.898              1.898 
    ##            SLC38A8         AL591025.1              OR6B3             OR13J1 
    ##              1.898              1.898              1.898              1.899 
    ##           KRTAP6-2           C10orf95              CLDN9              OR8K5 
    ##              1.899              1.899              1.899              1.899 
    ##             LRRC18       RP11-379H8.1             FAM24B              OR4X1 
    ##              1.899              1.899              1.899              1.899 
    ##         AC026202.1             ZNF658             OR10G9           C17orf50 
    ##              1.899              1.900              1.900              1.900 
    ##              PRSS3       RP11-410N8.4      RP11-867G23.8             ZNF845 
    ##              1.900              1.900              1.900              1.900 
    ##             NUDT16              PROK1             SDHAF1             ZNF107 
    ##              1.901              1.901              1.901              1.901 
    ##            C4orf36            FAM182B          HIST2H2BF             OR52N2 
    ##              1.901              1.901              1.901              1.902 
    ##             FCGR1B            TMEM107            CXorf64               MT1X 
    ##              1.902              1.902              1.902              1.902 
    ##              RGPD4           SERPINA3           C11orf40           KRTAP1-1 
    ##              1.902              1.902              1.902              1.903 
    ##               GJD3           DNAH10OS            TMEM243            HTATIP2 
    ##              1.903              1.903              1.903              1.903 
    ##               UPK2         AC008060.7            FAM209B              LRTM1 
    ##              1.903              1.903              1.903              1.903 
    ##            NEUROG3          KRTAP13-3             CALML3                PDF 
    ##              1.903              1.904              1.904              1.904 
    ##         AC114494.1                CCK              LCE2D              LCE3D 
    ##              1.904              1.904              1.904              1.904 
    ##             CHCHD2            DEFB112             SPRR2F            SPANXN5 
    ##              1.904              1.904              1.904              1.904 
    ##       RP11-408E5.4             CASC10               RLN2            GOLGA8O 
    ##              1.904              1.905              1.905              1.905 
    ##         AC138647.1             RNASE4            C5orf63              OR8J1 
    ##              1.905              1.905              1.905              1.905 
    ##           MAP1LC3C       RP11-536G4.1         RP11-6L6.2             OR11H4 
    ##              1.905              1.905              1.905              1.905 
    ##              OR1L3               MT1A            SCGB2A2             ZNF493 
    ##              1.905              1.905              1.905              1.905 
    ##              PRDX1               GSC2             CLEC3A                HBZ 
    ##              1.905              1.905              1.905              1.905 
    ##         AL391421.1              GPR31          LINC01119              KRT81 
    ##              1.905              1.905              1.906              1.906 
    ##         DNAH17-AS1             ADIPOQ          C10orf105               CCL7 
    ##              1.906              1.906              1.906              1.906 
    ##              OR1I1              OR2V1               TFF1             PRED62 
    ##              1.906              1.906              1.906              1.906 
    ##           C12orf76               GJB7             LY6G6E             RIIAD1 
    ##              1.906              1.906              1.906              1.906 
    ##              OR3A3             OR4C11           C19orf73              C1QL1 
    ##              1.907              1.907              1.907              1.907 
    ##            OR10AG1      RP11-552I14.1              SLIRP              OR1L4 
    ##              1.907              1.907              1.907              1.907 
    ##              CSHL1              FAM9A              LSMD1       RP11-127H5.1 
    ##              1.907              1.907              1.907              1.907 
    ##                NPW              FABP1           C1orf210              LCE1A 
    ##              1.907              1.907              1.907              1.907 
    ##               CT83               MT1M              HMGB4              COX19 
    ##              1.907              1.907              1.907              1.908 
    ##              MUCL1         AC008394.1      RP11-195B21.3             ZNF552 
    ##              1.908              1.908              1.908              1.908 
    ##          C14orf144              TEX37             IFIT1B             SPTSSB 
    ##              1.908              1.908              1.908              1.908 
    ##             RNASE6             OR52R1           SH3BGRL2              OR2S2 
    ##              1.909              1.909              1.909              1.909 
    ##              ZBED3        KB-1980E6.3               HOPX            NPIPB15 
    ##              1.909              1.910              1.910              1.910 
    ##         AC011294.3              GSTA5           KRTAP2-4           C1orf189 
    ##              1.910              1.910              1.910              1.910 
    ##            CHCHD10              OR1L6              LYPD3            FAM229A 
    ##              1.910              1.910              1.910              1.911 
    ##          C20orf187           DEFB106B             OR2AP1              NUPR1 
    ##              1.911              1.911              1.911              1.911 
    ##            C5orf56             CCDC54           PRAMEF22               HBQ1 
    ##              1.911              1.911              1.911              1.912 
    ##              TAF9B           KRTAP4-7              POTEI      RP11-187E13.2 
    ##              1.912              1.912              1.912              1.912 
    ##              OR8K1         AC016885.1         AP005482.1              PF4V1 
    ##              1.912              1.912              1.912              1.913 
    ##           SLC9B1P1              ABHD1           KRTAP5-4              LSM10 
    ##              1.913              1.913              1.913              1.913 
    ##             ZNF717             TAS2R4      RP11-160N1.10            BLOC1S3 
    ##              1.913              1.913              1.913              1.913 
    ##         AC137932.1               RMI2              OR8U1            TMEM86B 
    ##              1.913              1.913              1.913              1.913 
    ##               RBL1             MAGEC1              LCE1D             OR52N4 
    ##              1.913              1.913              1.913              1.913 
    ##              IFNL3              ISCA2                BIK              CXCL5 
    ##              1.914              1.914              1.914              1.914 
    ##          HNRNPA1L2      RP11-553A10.1             MALSU1              FBXL8 
    ##              1.914              1.914              1.914              1.914 
    ##             OR52J3             KLHL25          GOLGA6L19            TIMM10B 
    ##              1.915              1.915              1.915              1.915 
    ##              NUDT7       RP11-93B14.6               CSTA              CSTL1 
    ##              1.915              1.915              1.915              1.915 
    ##          HIST1H2BI       RP11-279O9.4      RP11-116D17.1            SCGB1A1 
    ##              1.915              1.915              1.915              1.915 
    ##           KRTAP2-1            C9orf16             OR51A4              VENTX 
    ##              1.915              1.916              1.916              1.916 
    ##               ZIK1              JSRP1          HIST1H2AL            FAM107A 
    ##              1.916              1.916              1.916              1.916 
    ##              OR2T1              LYRM1                LTB               DEC1 
    ##              1.916              1.917              1.917              1.917 
    ##                BID              INSL6               MSRA         RBM12B-AS1 
    ##              1.917              1.917              1.917              1.917 
    ##             FGFBP2              DSCR4               SELV               IL32 
    ##              1.917              1.917              1.917              1.917 
    ##             SPANXD                IL9              GNRH2             TMSB10 
    ##              1.917              1.917              1.918              1.918 
    ##             NPFFR1            DEFB110            C5orf27          HIST1H2BC 
    ##              1.918              1.918              1.918              1.918 
    ##              FOXR2               YBEY         AC073063.1             OR7A10 
    ##              1.918              1.918              1.918              1.918 
    ##            GOLGA8K              TCL1A              TDGF1               CTSZ 
    ##              1.919              1.919              1.919              1.919 
    ##             OR4D11               FLG2          HIST1H2BB              GGACT 
    ##              1.919              1.919              1.919              1.919 
    ##             BOD1L2              EGFL7         AC010646.3              OR5M3 
    ##              1.919              1.919              1.919              1.919 
    ##           C10orf62              ICAM3             NPIPB8             RPS4Y2 
    ##              1.920              1.920              1.920              1.920 
    ##             OR13H1           C6orf164                NPS          SPATA31A6 
    ##              1.920              1.920              1.920              1.920 
    ##             ZNF816              CCL17            FAM104B           C1orf147 
    ##              1.920              1.921              1.921              1.921 
    ##            FAM115C           C2orf27B      RP11-812E19.9                CRP 
    ##              1.921              1.921              1.921              1.921 
    ##        CTB-78H18.1         AC011997.1         HIST2H3PS2            TMEM236 
    ##              1.921              1.921              1.921              1.922 
    ##            RARRES3             SPINT3      RP11-817J15.3         AL022328.1 
    ##              1.922              1.922              1.922              1.922 
    ##         AL627171.1               FUT2              IGFL1             PLGLB1 
    ##              1.922              1.923              1.923              1.923 
    ##      RP11-521M14.2           C15orf53               PRB1            CCDC28B 
    ##              1.923              1.923              1.923              1.923 
    ##               CST7              CSAG1              CIDEB             SPANXC 
    ##              1.923              1.923              1.923              1.924 
    ##               ZACN           HIST2H3D         AP002348.1             DUOXA2 
    ##              1.924              1.924              1.924              1.924 
    ##              TEX22              THYN1              OR8J3              OR9A4 
    ##              1.924              1.924              1.924              1.924 
    ##                MT4            GOLGA8T               CT62              PRR25 
    ##              1.924              1.925              1.925              1.925 
    ##         AC018867.1       RP1-241P17.4              LCE2C              OR6F1 
    ##              1.925              1.925              1.925              1.925 
    ##              UTS2R             TCEAL6       RP11-15E18.4           FLJ20306 
    ##              1.925              1.925              1.926              1.926 
    ##          C20orf166           KIAA0087           C9orf116              NMRK2 
    ##              1.926              1.926              1.926              1.926 
    ##           U82695.9              EIF3C              LCE3B              OR5K4 
    ##              1.926              1.926              1.926              1.926 
    ##              OR8B2            SULT1C3              KLRF2            POLR2J2 
    ##              1.926              1.926              1.926              1.926 
    ##              KCNE1              LENEP            SPATA33               CST4 
    ##              1.927              1.927              1.927              1.927 
    ##              FGF22              MYEOV             NPIPB7             ZNF155 
    ##              1.927              1.927              1.927              1.927 
    ##             OR56A3              OR5L2             CXXC11          KRTAP24-1 
    ##              1.927              1.927              1.927              1.927 
    ##               PGA5      RP11-139J15.7           PPAPDC1B      RP11-321M21.3 
    ##              1.928              1.928              1.928              1.928 
    ##            PLA2G10               CRYZ         AC011475.1            AKR1B15 
    ##              1.928              1.928              1.928              1.928 
    ##               NEU2             AKR7A3            TAS2R39             OR2T10 
    ##              1.928              1.928              1.928              1.928 
    ##             CLEC1A             CLPSL1        AE000662.92             ZNF678 
    ##              1.929              1.929              1.929              1.929 
    ##               CGB8           HIST1H4I              CNPY1                C8G 
    ##              1.929              1.929              1.929              1.929 
    ##              CDPF1            C5orf47               SMN1             SPRR2D 
    ##              1.929              1.929              1.929              1.929 
    ##              PRAC1             OR10G4             OR10J1            PTPN20A 
    ##              1.929              1.929              1.929              1.929 
    ##            SCGB1C1         AL161645.2               WFS1             DEFB4B 
    ##              1.929              1.930              1.930              1.930 
    ##             SMIM21             SPRR1B               OTOS               CCR5 
    ##              1.930              1.930              1.930              1.930 
    ##              UBBP4            LRRN4CL             MRPL23      RP11-322E11.6 
    ##              1.930              1.931              1.931              1.931 
    ##              OR1F1          ANKRD20A3             NBPF15            ZNF705G 
    ##              1.931              1.931              1.931              1.931 
    ##               MC3R              OR2Z1             OR4C16             OR4K13 
    ##              1.931              1.931              1.931              1.931 
    ##              OR4Q3               IAPP              OR7D4             UBE2NL 
    ##              1.931              1.931              1.931              1.931 
    ##              INSL4             RASSF6             OR13C4            PRAMEF4 
    ##              1.932              1.932              1.932              1.932 
    ##                IL3             FAM24A             NPIPB6           TMEM170A 
    ##              1.932              1.932              1.932              1.933 
    ##               GLRX      RP11-324D17.1              AMY1B         AL139333.1 
    ##              1.933              1.933              1.933              1.933 
    ##         AL646016.1            IFI27L2               GJD4              MED11 
    ##              1.933              1.933              1.933              1.933 
    ##             OR51B5             OR52K2          KRTAP27-1           KRTAP4-9 
    ##              1.933              1.933              1.933              1.933 
    ##           KRTAP5-9             CLDN25       RP11-276H1.3              OVCH2 
    ##              1.933              1.934              1.934              1.934 
    ##             MBD3L5             HOXC12              TAAR5               ASMT 
    ##              1.934              1.934              1.934              1.934 
    ##               HBA2                HBD               FGL1           PRAMEF20 
    ##              1.934              1.934              1.934              1.934 
    ##             OR52B2             OR52E4            GOLGA6D         AC062017.1 
    ##              1.934              1.934              1.934              1.934 
    ##              OR5K1          HIST1H2AB               MUC7              OR9K2 
    ##              1.934              1.935              1.935              1.935 
    ##             UQCR11               GZMB              RAB42              APOC4 
    ##              1.935              1.935              1.935              1.935 
    ##               PRR9              GPHA2              FAAH2             OR5B21 
    ##              1.935              1.935              1.935              1.935 
    ##                SPN           C10orf55          HIST1H2AK             OR2T11 
    ##              1.936              1.936              1.936              1.936 
    ##             ACOT13               DRD5            FAM127C             OR52B4 
    ##              1.936              1.936              1.936              1.936 
    ##             NPBWR1           HIST1H4B              OR5I1            WFDC10A 
    ##              1.936              1.936              1.936              1.936 
    ##              FOLR4           C15ORF37               IDI2               SYT8 
    ##              1.936              1.936              1.936              1.936 
    ##            DEFB125             NUTM2A          HIST1H2AG             OR13C3 
    ##              1.936              1.936              1.937              1.937 
    ##              PPDPF       RP11-219B4.5              KCNE2             GPR141 
    ##              1.937              1.937              1.937              1.937 
    ##            TAS2R41         AL135998.1             ZNF468            TMEM105 
    ##              1.937              1.937              1.937              1.937 
    ##             OR4A47           KRTAP3-3         AC005544.1             KDELR3 
    ##              1.938              1.938              1.938              1.938 
    ##             OR6C65             OR6C70              FAM9C               CST2 
    ##              1.938              1.938              1.939              1.939 
    ##            C2orf88               EME2             NBPF14             OR5D14 
    ##              1.939              1.939              1.939              1.939 
    ##             FAM27B               AMZ1             NPIPA8         AL137026.1 
    ##              1.939              1.939              1.939              1.940 
    ##           HIST1H4E            HIST3H3             OR5AR1             OR6C68 
    ##              1.940              1.940              1.940              1.940 
    ##            DEFB131             SPRR2E              ARMS2       RP11-297N6.4 
    ##              1.940              1.940              1.940              1.940 
    ##              HCAR3             NPIPA7         AC007401.2         AL033381.1 
    ##              1.940              1.940              1.940              1.940 
    ##              RFPL3         AL353354.1           HIST1H4G           PRAMEF19 
    ##              1.940              1.941              1.941              1.941 
    ##     RP11-481A20.11             PTRHD1               CORT            CCDC167 
    ##              1.941              1.941              1.941              1.941 
    ##             NDUFB3             OR10K2             NUPR1L             ZNF785 
    ##              1.941              1.941              1.941              1.941 
    ##              GAGE1               PTH2          LINC01118              NUDT1 
    ##              1.941              1.942              1.942              1.942 
    ##          HIST1H2AH               SAA2         AC005841.1   XXbac-BPG32J3.20 
    ##              1.942              1.942              1.942              1.942 
    ##      CTD-2054N24.2          HIST1H2BN             OR5H14               PSG4 
    ##              1.942              1.942              1.942              1.943 
    ##               BEX5             C2CD4B             SSMEM1               VCX2 
    ##              1.943              1.943              1.943              1.943 
    ##            C8orf74          SPATA31A7          KRTAP10-6              OR4K1 
    ##              1.943              1.943              1.943              1.943 
    ##              ALG1L        RP11-67H2.1               PXT1            SLC10A2 
    ##              1.944              1.944              1.944              1.944 
    ##          ANKRD20A4           C19orf80              MESP1             NUTM2E 
    ##              1.944              1.944              1.944              1.944 
    ##       RP11-247C2.2               LMF2         AC010441.1            S100A7A 
    ##              1.944              1.944              1.944              1.944 
    ##                OCM             OR51A7           TXNRD3NB             CYP2D6 
    ##              1.944              1.944              1.944              1.944 
    ##           HIST1H3B              TIFAB               DPRX           KIAA1456 
    ##              1.944              1.945              1.945              1.945 
    ##              OR7C2          C14orf178         AC008372.1            TMEM207 
    ##              1.945              1.945              1.945              1.945 
    ##           C9orf139         AC016757.3         AL355490.1          HIST1H2BO 
    ##              1.945              1.945              1.945              1.945 
    ##              AMELY                FTL             POLR2L           KRTAP4-2 
    ##              1.945              1.946              1.946              1.946 
    ##              OR8H1               CD1A        hsa-mir-150             ACTBL2 
    ##              1.946              1.946              1.946              1.946 
    ##              GNG13              CBWD5              IQCF3          HIST1H2AD 
    ##              1.946              1.946              1.947              1.947 
    ##             NPIPB3           KRTAP6-1            TBC1D3F            TMEM14B 
    ##              1.947              1.947              1.947              1.947 
    ##               HES4          KRTAP13-2             ZNF812         AC006435.1 
    ##              1.947              1.947              1.947              1.947 
    ##           C10orf25            FAM211B           TP53AIP1             DNASE1 
    ##              1.948              1.948              1.948              1.948 
    ##              MPZL3           C13orf45              CXCL2               HRNR 
    ##              1.948              1.948              1.948              1.948 
    ##       RP11-527L4.2             GAGE2B           TRAPPC2L          C21orf128 
    ##              1.948              1.948              1.948              1.949 
    ##             IFITM5             SPINK6          C20orf144            FAM111A 
    ##              1.949              1.949              1.949              1.949 
    ##            FAM86C1               NME3              DIRC3            PGPEP1L 
    ##              1.949              1.949              1.950              1.950 
    ##             ZNF701              OR4A5             OR10H3             BRICD5 
    ##              1.950              1.950              1.950              1.950 
    ##             MS4A6A       RP11-514P8.7             TPSAB1             OR2T27 
    ##              1.950              1.950              1.950              1.950 
    ##               MT2A             OR5M10             GPBAR1           PRAMEF17 
    ##              1.950              1.951              1.951              1.951 
    ##           PRAMEF26            C3orf72             MRPL34            SPAG11A 
    ##              1.951              1.951              1.951              1.951 
    ##               RPA3             SERHL2        AC009892.10            DEFB115 
    ##              1.951              1.952              1.952              1.952 
    ##              MLANA             CT47B1            GAGE12J                CPZ 
    ##              1.952              1.953              1.953              1.953 
    ##             SMIM22        AP003068.23               AQP7           HIST1H4D 
    ##              1.953              1.953              1.953              1.953 
    ##           C19orf24              BGLAP              RP1L1           C19orf83 
    ##              1.954              1.954              1.954              1.954 
    ##              PATE1            DEFB116               FUT3              CRCT1 
    ##              1.954              1.954              1.954              1.954 
    ##           HIST1H4C           HIST1H4L               PSG8             OR11H1 
    ##              1.954              1.954              1.955              1.955 
    ##              OXCT2            SULT6B1              OR5F1             ZNF836 
    ##              1.955              1.955              1.955              1.955 
    ##          ANKRD20A2               AHSP              HMHB1             OR56B1 
    ##              1.956              1.956              1.956              1.956 
    ##             S100A2              FABP9              OR4X2             GIMAP4 
    ##              1.956              1.956              1.956              1.957 
    ##              DSCR8             DUSP21            DEFB132            DEFB134 
    ##              1.957              1.957              1.957              1.957 
    ##            C9orf92         AL133481.1          CCDC144NL      RP11-272B17.2 
    ##              1.957              1.957              1.957              1.957 
    ##             S100A5             OR4C15               MYL5      CH17-132F21.1 
    ##              1.957              1.957              1.958              1.958 
    ##               APRT       RP11-351M8.1          TNFRSF13B              OR4M2 
    ##              1.958              1.958              1.958              1.958 
    ##           HIST1H3E          HIST3H2BB             OR5H15              OR6C4 
    ##              1.959              1.959              1.959              1.959 
    ##             ZNF765             OR14J1           DEFB106A              GNG10 
    ##              1.959              1.959              1.959              1.959 
    ##         AL645922.1           HIST1H4K          SERPINB12                FLG 
    ##              1.960              1.960              1.960              1.960 
    ##              LEUTX             OR52M1               CST1           DEFB108B 
    ##              1.960              1.960              1.960              1.960 
    ##               KNCN            ANKRD65              MZT2A      RP11-144F15.1 
    ##              1.961              1.961              1.961              1.962 
    ##             ATP5L2         AC004466.1           C9orf170              TBL1Y 
    ##              1.962              1.962              1.962              1.963 
    ##            TBC1D3C              PDIA2          HIST1H2AA             TMEM80 
    ##              1.963              1.963              1.963              1.963 
    ##           HIST1H3G           HIST1H4H              LYRM5              SSTR4 
    ##              1.963              1.963              1.964              1.964 
    ##                TYR            C2orf83            PLA2G1B      RP11-318A15.7 
    ##              1.964              1.964              1.965              1.965 
    ##            TMEM125               IMP3              TPSG1         AC005493.1 
    ##              1.965              1.966              1.966              1.966 
    ##               GHRL             DEFB4A              NBPF8              NBPF9 
    ##              1.966              1.967              1.967              1.967 
    ##      RP11-467N20.5              SMCO4               MC1R              FAHD1 
    ##              1.967              1.967              1.967              1.967 
    ##           HIST1H4A                SCT              IGLL5        TMPRSS11BNL 
    ##              1.967              1.967              1.968              1.968 
    ##              DEFA5       RP11-297M9.1            SCGB1D1           C17orf97 
    ##              1.968              1.968              1.968              1.969 
    ##              OR4S2             OR52K1          HIST1H2BM            FAM167A 
    ##              1.969              1.969              1.969              1.970 
    ##              CXCL6          HIST1H2BL             NBPF20              DEFA6 
    ##              1.970              1.970              1.971              1.971 
    ##          HIST1H2BD          HIST1H2BE          HIST1H2BK            C1QTNF8 
    ##              1.971              1.971              1.971              1.972 
    ##          HIST1H2BH           HIST1H3A       CTD-2228K2.5               MIB1 
    ##              1.972              1.972              1.972              1.973 
    ##              SMIM3          HIST1H2AM             CT45A5              KCNV2 
    ##              1.973              1.973              1.973              1.973 
    ##          HIST1H2BG               PRB2              OXER1           HIST1H1A 
    ##              1.974              1.974              1.975              1.975 
    ##              SMAD6            SLC10A1          KRTAP5-10             OR4C13 
    ##              1.975              1.975              1.975              1.975 
    ##             NUDT19              CST9L            GOLGA8S             NBPF16 
    ##              1.975              1.975              1.976              1.976 
    ##              OR4C3               AGR3            GOLGA8H              PUSL1 
    ##              1.976              1.976              1.977              1.977 
    ##             OR2AG1               PRB4             C4orf6            PRAMEF6 
    ##              1.977              1.977              1.978              1.978 
    ##                HBB             OR4C46            DEFB135            GOLGA8R 
    ##              1.978              1.980              1.980              1.981 
    ##               PSG2          HIST1H2BF               GJB2             TRIM48 
    ##              1.981              1.981              1.984              1.985 
    ##               PSG1             NBPF10              LZTR1               CD36 
    ##              1.986              1.988              1.994              1.995 
    ##               SSX1             COX7B2             DUX4L4              EBLN1 
    ##              1.996                 NA                 NA                 NA 
    ##              EBLN2              H1FNT             H2AFB1              H3F3C 
    ##                 NA                 NA                 NA                 NA 
    ##              HUS1B               IRGM             MT1HL1              CKS1B 
    ##                 NA                 NA                 NA                 NA 
    ##             RFPL4B            FAM183B              VN1R1              VN1R2 
    ##                 NA                 NA                 NA                 NA 
    ##              VN1R4              HEPN1             PRKACG      RP11-382J12.1 
    ##                 NA                 NA                 NA                 NA 
    ##         RP11-3B7.1      RP11-429E11.3      RP11-480I12.4      RP11-863K10.7 
    ##                 NA                 NA                 NA                 NA 
    ##       RP11-89N17.1               SMCP            TMEM14E              AGAP4 
    ##                 NA                 NA                 NA                 NA 
    ##              AGAP8             ANP32D            CCDC140           C19orf48 
    ##                 NA                 NA                 NA                 NA 
    ##              EID2B            C8orf49            C8orf56           C9orf163 
    ##                 NA                 NA                 NA                 NA 
    ##             FAM71A               NAT8             TSPYL6               TSRM 
    ##                 NA                 NA                 NA                 NA 
    ##            CCDC177          KRTAP12-1          KRTAP12-2          KRTAP12-3 
    ##                 NA                 NA                 NA                 NA 
    ##          KRTAP12-4          KRTAP13-1          KRTAP15-1          KRTAP19-1 
    ##                 NA                 NA                 NA                 NA 
    ##          KRTAP19-2          KRTAP19-3          KRTAP19-4          KRTAP19-5 
    ##                 NA                 NA                 NA                 NA 
    ##          KRTAP19-6         BX088651.1         BX255923.1         AF131215.5 
    ##                 NA                 NA                 NA                 NA 
    ##              OR2Y1             OR4A15             OR51A2              HSPB9 
    ##                 NA                 NA                 NA                 NA 
    ##             PCED1B             FKSG52             ZNF843             ZNF860 
    ##                 NA                 NA                 NA                 NA 
    ##              GVQW1             NCBP2L              MAS1L          POM121L12 
    ##                 NA                 NA                 NA                 NA 
    ##              PTTG2              ZFP42           C12orf61          C14orf132 
    ##                 NA                 NA                 NA                 NA 
    ##           FLJ00104           FLJ00418           FLJ20373           FLJ45079 
    ##                 NA                 NA                 NA                 NA 
    ##             MAGEA1            MAGEA10            MAGEA12             MAGEA3 
    ##                 NA                 NA                 NA                 NA 
    ##             MAGEA4             MAGEA6             MAGEA8            MAGEA9B 
    ##                 NA                 NA                 NA                 NA 
    ##             MAGEB1            MAGEB10            MAGEB17             MAGEB2 
    ##                 NA                 NA                 NA                 NA 
    ##             MAGEB4             MAGEB5             MAGEB6             MAGEC2 
    ##                 NA                 NA                 NA                 NA 
    ##              PBOV1              PLAC1              SPHAR          GABARAPL3 
    ##                 NA                 NA                 NA                 NA 
    ##          SSBP3-AS1             ANXA2R         AP000679.2         AP000695.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AP000867.1         AP000974.1         AP001350.1         AP001816.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AP001885.1         AP003062.1         AP003733.1             FTHL17 
    ##                 NA                 NA                 NA                 NA 
    ##              HRCT1             MRGPRD             MRGPRE             MRGPRG 
    ##                 NA                 NA                 NA                 NA 
    ##            MRGPRX1            MRGPRX2            MRGPRX3            MRGPRX4 
    ##                 NA                 NA                 NA                 NA 
    ##            OR11H12             OR13C5            OR14A16             OR14A2 
    ##                 NA                 NA                 NA                 NA 
    ##            OR14C36             OR14I1             OR14K1              OR1M1 
    ##                 NA                 NA                 NA                 NA 
    ##              OR1N1                STH             RAB40A              RAB6C 
    ##                 NA                 NA                 NA                 NA 
    ##              SDIM1              SUMO4              SRRM5                SRY 
    ##                 NA                 NA                 NA                 NA 
    ##              TEX19             ZNF788              ZNF80             TRIM60 
    ##                 NA                 NA                 NA                 NA 
    ##           DEFB107B             ADAM20             ADAM29             ADAM30 
    ##                 NA                 NA                 NA                 NA 
    ##              GPR32                IVL               SBP1            C9orf37 
    ##                 NA                 NA                 NA                 NA 
    ##            C9orf38            C9orf53            FOXD4L5            FOXD4L6 
    ##                 NA                 NA                 NA                 NA 
    ##           MTRNR2L1          MTRNR2L10          MTRNR2L12          MTRNR2L13 
    ##                 NA                 NA                 NA                 NA 
    ##           MTRNR2L2           MTRNR2L3           MTRNR2L4           MTRNR2L5 
    ##                 NA                 NA                 NA                 NA 
    ##           MTRNR2L6           MTRNR2L7           MTRNR2L8           MTRNR2L9 
    ##                 NA                 NA                 NA                 NA 
    ##             SPRR2A             SPRR2B             SPRR2G              SPRR3 
    ##                 NA                 NA                 NA                 NA 
    ##            TMEM133             TLX1NB              TPRXL            C3orf27 
    ##                 NA                 NA                 NA                 NA 
    ##              CEMP1          LINC00346          LINC00696          LINC01101 
    ##                 NA                 NA                 NA                 NA 
    ##          LINC01124              PNMA3              PNMA5             PNMAL2 
    ##                 NA                 NA                 NA                 NA 
    ##           C15orf37               VHLL           C17orf77           C17orf82 
    ##                 NA                 NA                 NA                 NA 
    ##             CALML5             FCGR2C            FP15737               FPR1 
    ##                 NA                 NA                 NA                 NA 
    ##             GIMAP7              PGAM4             PSAPL1               RPA4 
    ##                 NA                 NA                 NA                 NA 
    ##           USP17L10           USP17L13           USP17L15           USP17L17 
    ##                 NA                 NA                 NA                 NA 
    ##           USP17L18            USP17L2            USP17L5         AC010536.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AC011484.1         AC012215.1         AC012493.2         AC013449.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AC015989.2         AC016586.1         AC016745.1         AC019294.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AC020907.1         AC020922.1         AC021860.1         AC022498.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AC023632.1         AC024940.1         AC025278.1         AC026310.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AC026703.1         AC037199.1         AC061992.1         AC073342.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AL109927.1         AL133373.1         AL136115.1         AL136531.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AL138764.1         AL139099.1         AL139147.1         AL158091.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AL161915.1         AL162424.1         AL353354.2         AL353698.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AL357673.1         AL358113.1         AL359195.1         AL359878.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AL360004.1         AL391152.1         AL445989.1         AL583828.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AL590822.1         AL590822.2         AL603965.1         AL627171.2 
    ##                 NA                 NA                 NA                 NA 
    ##         AL645728.1         AL645730.2         AL807752.1             ERV3-1 
    ##                 NA                 NA                 NA                 NA 
    ##         ERVMER34-1             ERVV-1             ERVV-2             ERVW-1 
    ##                 NA                 NA                 NA                 NA 
    ##               F8A1              LCE3A              LCE4A              LCE6A 
    ##                 NA                 NA                 NA                 NA 
    ##              OR5A2             OR5D13             OR5D16              OR5L1 
    ##                 NA                 NA                 NA                 NA 
    ##              OR5T1              OR7A5             RBMXL3            SLC35G3 
    ##                 NA                 NA                 NA                 NA 
    ##            SLC35G5      WI2-3308P17.2           C6orf120           C6orf141 
    ##                 NA                 NA                 NA                 NA 
    ##           C6orf195            C7orf13            C8orf17               CDR1 
    ##                 NA                 NA                 NA                 NA 
    ##             COLCA1             CRIPAK            FAM106A               FPR3 
    ##                 NA                 NA                 NA                 NA 
    ##               FUT5               FUT6             ATXN3L          KRTAP19-7 
    ##                 NA                 NA                 NA                 NA 
    ##          KRTAP19-8          KRTAP20-1          KRTAP20-2          KRTAP20-3 
    ##                 NA                 NA                 NA                 NA 
    ##          KRTAP20-4          KRTAP21-2          KRTAP21-3          KRTAP22-1 
    ##                 NA                 NA                 NA                 NA 
    ##          KRTAP22-2          KRTAP23-1          KRTAP25-1         KRTAP4-16P 
    ##                 NA                 NA                 NA                 NA 
    ##           KRTAP4-3           KRTAP9-1           KRTAP9-3           KRTAP9-4 
    ##                 NA                 NA                 NA                 NA 
    ##           KRTAP9-8               LSP1             OR2AK2              OR2L2 
    ##                 NA                 NA                 NA                 NA 
    ##              OR2L3              OR2L5              OR2M2              OR2M3 
    ##                 NA                 NA                 NA                 NA 
    ##              OR2M5              OR2M7             OR2T12              OR2T3 
    ##                 NA                 NA                 NA                 NA 
    ##             OR2T34              OR2T8            PPP1R26              PRR21 
    ##                 NA                 NA                 NA                 NA 
    ##             PRR23A             PRR23B             PRR23C            SPATA12 
    ##                 NA                 NA                 NA                 NA 
    ##              LYPD8            TBC1D3G              TFDP3              USP26 
    ##                 NA                 NA                 NA                 NA 
    ##              USP29               SPZ1         AP006621.5          C10orf111 
    ##                 NA                 NA                 NA                 NA 
    ##           C10orf40           C10orf85             TCEB3B             TCEB3C 
    ##                 NA                 NA                 NA                 NA 
    ##            TCEB3CL           TCEB3CL2           C11orf72           C11orf89 
    ##                 NA                 NA                 NA                 NA 
    ##            C3orf36            C3orf56       EPB41L4A-AS2              TIAF1 
    ##                 NA                 NA                 NA                 NA 
    ##             TICAM1             ZNF781            FAM215A            FAM217B 
    ##                 NA                 NA                 NA                 NA 
    ##            FAM218A            FAM220A              ZNF83             ZNF835 
    ##                 NA                 NA                 NA                 NA 
    ##             CTAGE1            CTAGE15             CTAGE4             CTAGE6 
    ##                 NA                 NA                 NA                 NA 
    ##             CTAGE9        CTB-58E17.5        CTC-241N9.1             HIGD2B 
    ##                 NA                 NA                 NA                 NA 
    ##           HIST1H1T      CTD-2144E22.5         AC124890.1         AC127496.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AC129492.6         AC131971.1         AC132186.1         AC135178.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AC139426.2         AC142381.1         AC145676.2         AC226150.4 
    ##                 NA                 NA                 NA                 NA 
    ##              LENG9            ARIH2OS             MRPL36            PP13439 
    ##                 NA                 NA                 NA                 NA 
    ##      RP11-204N11.1      RP11-268J15.5       RP11-315O6.2      RP11-342M21.2 
    ##                 NA                 NA                 NA                 NA 
    ##      RP11-347C12.1            TGIF2LX            TGIF2LY              PYDC1 
    ##                 NA                 NA                 NA                 NA 
    ##              PYDC2            RNASE11             RNASE2             RNASE3 
    ##                 NA                 NA                 NA                 NA 
    ##             RNASE8      RP11-1055B8.6     RP11-1070N10.3       RP11-156E8.1 
    ##                 NA                 NA                 NA                 NA 
    ##            DCAF4L2            DCAF8L1               FSCB             TMEM75 
    ##                 NA                 NA                 NA                 NA 
    ##             TMEM78         AC003102.1         AC004899.1         AC005477.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AC006372.1         AC007421.1        AC008132.13         AC008267.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AC010327.2         AL020996.1         AL031320.1         AL049840.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AL079342.1              LUZP6               OCLM             RPL39L 
    ##                 NA                 NA                 NA                 NA 
    ##             TAS2R1            TAS2R10            TAS2R13            TAS2R14 
    ##                 NA                 NA                 NA                 NA 
    ##            TAS2R16            TAS2R19            TAS2R20            TAS2R30 
    ##                 NA                 NA                 NA                 NA 
    ##            TAS2R31            TAS2R42            TAS2R43            TAS2R46 
    ##                 NA                 NA                 NA                 NA 
    ##            TAS2R50            TAS2R60             TAS2R8             TAS2R9 
    ##                 NA                 NA                 NA                 NA 
    ##         AC074389.6         AC079210.1         AC079341.1         AC079354.2 
    ##                 NA                 NA                 NA                 NA 
    ##         AC083862.1         AC090186.1         AC092811.1         AC093677.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AC097381.1         AC103801.2         AC104472.1         AC105020.1 
    ##                 NA                 NA                 NA                 NA 
    ##         AC107021.1         AC112693.2         AC114546.1         AC115618.1 
    ##                 NA                 NA                 NA                 NA 
    ##           ARHGEF35           C1orf229           C22orf29            CXorf24 
    ##                 NA                 NA                 NA                 NA 
    ##            CXorf27            CXorf67             EDDM3A             EDDM3B 
    ##                 NA                 NA                 NA                 NA 
    ##              IFNA1             IFNA10             IFNA13             IFNA14 
    ##                 NA                 NA                 NA                 NA 
    ##             IFNA16             IFNA17              IFNA2             IFNA21 
    ##                 NA                 NA                 NA                 NA 
    ##              IFNA4              IFNA5              IFNA6              IFNA7 
    ##                 NA                 NA                 NA                 NA 
    ##              IFNA8              IFNB1              IFNW1             MBD3L1 
    ##                 NA                 NA                 NA                 NA 
    ##             MBD3L2           MGC10955             NAP1L6            PPIAL4A 
    ##                 NA                 NA                 NA                 NA 
    ##            PPIAL4G               PRNT               QRFP           REXO1L1P 
    ##                 NA                 NA                 NA                 NA 
    ##              SMIM6       ADAMTSL4-AS1              BIRC8              BLACE 
    ##                 NA                 NA                 NA                 NA 
    ##               BLID            C16orf3          C21orf119             DIRAS3 
    ##                 NA                 NA                 NA                 NA 
    ##              DIRC1      DKFZP434O1614      DKFZP667F0711            FAM231B 
    ##                 NA                 NA                 NA                 NA 
    ##            FAM231D            FAM27E1            FAM27E2             FAM47A 
    ##                 NA                 NA                 NA                 NA 
    ##             FAM47B             FAM47C             NHLRC4              NPAP1 
    ##                 NA                 NA                 NA                 NA 
    ##       SLC25A21-AS1             TMEM99           Z97053.1             ZNF491 
    ##                 NA                 NA                 NA                 NA 
    ##             ZNF497            ZNF705D            C5orf20            C5orf55 
    ##                 NA                 NA                 NA                 NA 
    ##            C6ORF50             CPXCR1          DHRS4-AS1          GRIK1-AS2 
    ##                 NA                 NA                 NA                 NA 
    ##          IBA57-AS1             NPIPA2             NPIPB5              OR7C1 
    ##                 NA                 NA                 NA                 NA 
    ##              OR7D2             OR7E24              OR7G1              OR7G2 
    ##                 NA                 NA                 NA                 NA 
    ##              OR7G3              OR8D1             ZNF600             ZNF645 
    ##                 NA                 NA                 NA                 NA

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
