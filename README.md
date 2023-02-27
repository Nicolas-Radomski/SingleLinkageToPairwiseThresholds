# Usage
The R scripts SingleLinkage2PairwiseThresholds.R (detailed algorithm with Rstudio) and SingleLinkageToPairwiseThresholds.R (automatic algorithm with Rscript) aim at identifying relationship between thresholds of pairwise differences and single linkage in a context of outbreak investigation based on microbial mutations (e.g. cg/wgMLST, genes, SNPs, InDels, kmers) (Profiles.csv).
# Objectives
The terms sensitivity and specificity are defined as the proportions of true positive and true negative samples, respectively. Regarding the main goal to identify pairwise difference thresholds according to single linkage thresholds, the successive objectives aim at:
- Objective 1: creating big single linkage clusters at four different single linkage thresholds (SL),
- Objective 2: splitting randomly samples from each single linkage clusters,
- Objective 3: assessing correlation between intra cluster pairwise differences (PW) and single linkage thresholds (SL),
- Objective 4: optimizing sensibility and specificity of pairwise differences between a selected single linkage cluster and the other.
# Approaches
After the creation of single linkage cluster from mutation profiles (i.e. arguments argSLA, argSLB, argSLC and argSLD), the clusters presenting many samples are retained (i.e. argument argMNS) to avoid small clusters which would not allow estimation of accurate sensitivity. Then, these big single linkage clusters are randomly splitted to create two collections (i.e. argument argPSP). One dataset is dedicated to a correlation approach, while the other dataset is dedicated to a optimization approach, both with the objective to identify relationship between single linkage thresholds and pairwise difference thresholds.
## Correlation approach
The correlation approach aims successively at:
- extracting pairwise differences inside each single linkage cluster (i.e. intra cluster tag in the R script) of each single linkage thresholds,
- estimating variability of intra cluster pairwise differences for each single linkage cluster of each single linkage thresholds,
- and screening for correlation between intra cluster pairwise differences and single linkage thresholds, or between maximal values of intra cluster pairwise differences and single linkage thresholds.
## Optimization approach
The optimization approach aims successively at:
- extracting pairwise differences between all single linkage clusters (i.e. extra cluster tag in the R script),
- identifying the single linkage cluster presenting the minimal extra cluster pairwise difference in order to optimize downstream estimations of sensitivity and specificity of pairwise differences,
- randomly selecting one cluster in case of multiple clusters presenting identical minimal extra cluster pairwise differences,
- extracting pairwise differences between the single linkage cluster of interest and itself (i.e. intra cluster tag in the R script), and between the single linkage cluster of interest and the other single linkage clusters for each single linkage thresholds (i.e. extra cluster tag in the R script),
- calculating sensitivity and specificity of pairwise differences through Receiver Operating Characteristic (ROC) analyses for each single linkage thresholds according to true positive (i.e. intra cluster tag in the R script) and true negative (i.e. extra cluster tag in the R script) samples,
- and estimating the pairwise difference thresholds presenting 100% of sensitivity and 100% of specificity for each single linkage thresholds.
# Input
## Explainations (i.e. Profiles.csv)
- S stands for sample
- L stands for locus
- G stands for genotype
- 0 stands for missing data
- The string "sample" must appear in the first line of the first column
## Profiles of microbial mutations (i.e. Profiles.csv)
```
sample  L1  L2  L3  L4  L5  L6  L7  L8   L9  L10 L11 L12 L13 L14 L15 L16 L17 L18 L19 L20 L21 L22 L23 L24 L25 L26 L27 L28 L29 L30 L31 L32 L33 L34 L35
    S1 G20 G15 G55 G12 G30 G11 G24 G66  G12  G55 G66  G5 G86 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
    S2 G20 G15 G55 G12 G30 G11 G24 G66  G12  G55 G66  G2 G87 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
    S3 G20 G15 G55 G12 G30 G11 G24 G66  G12  G55 G66  G2 G86 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
    S4 G20 G31 G55 G30 G30 G11 G55 G66  G55  G55 G66  G5 G87 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
    S5 G20 G31 G55 G30 G30 G11 G55 G66  G55  G55 G66  G5 G98 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
    S6 G10 G15 G10 G12 G30 G10 G24 G10  G12  G55 G66  G5 G98 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
    S7 G41 G22 G41 G22 G22 G41 G27 G41  G27  G27 G66  G9 G86 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
    S8 G41 G22 G41 G22 G22 G41 G27 G41  G27  G27 G66  G9 G86 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
    S9 G41 G15 G41 G12 G30 G41 G24 G41  G12  G55 G66  G8 G97 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
   S10 G50 G22 G50 G22 G55 G51 G27 G50  G27  G66 G66  G8 G97 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
   S11 G10 G54 G15 G41 G65 G88 G75 G89 G420 G998 G66  G5 G86 G11 G10  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
   S12 G10 G54 G98 G41 G65 G88 G75 G89 G420 G998 G66  G8 G86 G14  G1  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
   S13 G20 G15 G55 G12 G30 G11 G24 G66  G12  G55 G66  G2 G87 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
   S14 G41 G15 G41 G12 G30 G41 G24 G41  G12  G55 G66  G8 G97 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
   S15 G20 G15 G55 G12 G30 G11 G24 G66  G12  G55 G66  G5 G86 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
   S16 G10 G54 G15 G41 G65 G88 G75 G89 G420 G998 G66  G5 G86 G11 G10  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
   S17 G20 G31 G55 G30 G30 G11 G55 G66  G55  G55 G66  G5 G87 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
   S18 G41 G22 G41 G22 G22 G41 G27 G41  G27  G27 G66  G9 G86 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
   S19 G41 G22 G41 G22 G22 G41 G27 G41  G27  G27 G66  G9 G86 G54 G47  G5  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19  G6  G6  G6  G6  G6
   S20 G10 G54 G98 G41 G65 G88 G75 G89 G420 G998 G66  G8 G86 G14  G1  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19 G20  G6  G6  G6  G6  G6
   S21 G10 G15 G10 G12 G30 G10 G24 G10  G12  G55 G66  G5 G98 G54 G47  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19 G20  G6  G6  G6  G6  G6
   S22 G41 G22 G41 G22 G22 G41 G27 G41  G27  G27 G66  G9 G86 G54 G47  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19 G20  G6  G6  G6  G6  G6
   S23 G41 G22 G41 G22 G22 G41 G27 G41  G27  G27 G66  G9 G86 G54 G47  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19 G20  G6  G6  G6  G6  G6
   S24 G41 G15 G41 G12 G30 G41 G24 G41  G12  G55 G66  G8 G97 G54 G47  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19 G20  G6  G6  G6  G6  G6
   S25 G50 G22 G50 G22 G55 G51 G27 G50  G27  G66 G66  G8 G97 G54 G47  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19 G20  G6  G6   0  G6  G6
   S26 G10 G54 G15 G41 G65 G88 G75 G89 G420 G998 G66  G5 G86 G11 G10  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19 G20  G6  G6  G6  G6  G6
   S27 G10 G54 G98 G41 G65 G88 G75 G89 G420 G998 G66  G8 G86 G14  G1  G6  G7  G8  G9 G10 G11 G12 G13 G14 G15 G16 G17 G18 G19 G20  G6  G6  G6  G6  G6
```
# Arguments
The arguments of this R script have been defined according to the objectives.
- argDFP: dataframe of profiles (i.e. Profiles.csv)
- argSLA: single linkage threshold A (i.e. 4)
- argSLB: single linkage threshold B (i.e. 7)
- argSLC: single linkage threshold C (i.e. 10)
- argSLD: single linkage threshold D (i.e. 13)
- argMNS: minimum number of samples to retain inside each single linkage cluster (i.e. 6)
- argPSP: proportion of spiting (%) defining the optimization dataframes  (i.e. 50)
# Output
## Objective 1:
- Objective1-Dendrogram.pdf
- Objective1-PairwiseMatrix.csv
- Objective1-SingleLinkageAllClusters.csv
- Objective1-SingleLinkageAllClustersStatPerCluster.csv
- Objective1-SingleLinkageAllClustersStatPerSL.csv
- Objective1-SingleLinkageBigClusters.csv
- Objective1-SingleLinkageBigClustersStatPerCluster.csv
- Objective1-SingleLinkageBigClustersStatPerSL.csv
## Objective 2:
- Objective2-SingleLinkageBigClustersCorrelation.csv
- Objective2-SingleLinkageBigClustersCorrelationStatPerCluster.csv
- Objective2-SingleLinkageBigClustersCorrelationStatPerSL.csv
- Objective2-SingleLinkageBigClustersOptimization.csv
- Objective2-SingleLinkageBigClustersOptimizationStatPerCluster.csv
- Objective2-SingleLinkageBigClustersOptimizationStatPerSL.csv
## Objective 3:
- Objective3-PairwiseSingleLinkageCorrelationBoxplot.pdf
- Objective3-PairwiseSingleLinkageCorrelationBoxplot.tiff
- Objective3-PairwiseSingleLinkageCorrelationDataframe.csv
- Objective3-PairwiseSingleLinkageCorrelationLinearMax.pdf
- Objective3-PairwiseSingleLinkageCorrelationLinearMax.tiff
- Objective3-PairwiseSingleLinkageCorrelationLinear.pdf
- Objective3-PairwiseSingleLinkageCorrelationLinear.tiff
- Objective3-PairwiseSingleLinkageCorrelationStat.csv
## Objective 4:
- Objective4-PairwiseSingleLinkageOptimizationDataframe.csv
- Objective4-PairwiseSingleLinkageOptimizationSLAThresholds.csv
- Objective4-PairwiseSingleLinkageOptimizationSLBThresholds.csv
- Objective4-PairwiseSingleLinkageOptimizationSLCThresholds.csv
- Objective4-PairwiseSingleLinkageOptimizationSLDThresholds.csv
- Objective4-PairwiseSingleLinkageOptimizationThresholdsBest.txt

# Install R dependencies and launch with R
## 1/ Install dependencies from the R console
The R scripts SingleLinkage2PairwiseThresholds.R (detailed algorithm with Rstudio) and SingleLinkageToPairwiseThresholds.R (automatic algorithm with Rscript) were prepared and tested with R version 4.1.2 and RStudio 2021.09.1.
```
install.packages("data.table")
install.packages("spaa")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("forcats")
install.packages("ggpmisc")
install.packages("pROC")
```
## 2/ Launch each command from Rstudio (i.e. SingleLinkage2PairwiseThresholds.R detailed algorithm)
```
git clone https://github.com/Nicolas-Radomski/SingleLinkageToPairwiseThresholds.git
cd SingleLinkageToPairwiseThresholds
rstudio SingleLinkage2PairwiseThresholds.R
```
## 3/ Launch the whole script from Rscript (i.e. SingleLinkageToPairwiseThresholds.R automatic algorithm)
```
git clone https://github.com/Nicolas-Radomski/SingleLinkageToPairwiseThresholds.git
cd SingleLinkageToPairwiseThresholds
Rscript SingleLinkageToPairwiseThresholds.R Profiles.csv 4 7 10 13 6 50
```
# Illustration
![Non-parametric figure](https://github.com/Nicolas-Radomski/SingleLinkageToPairwiseThresholds/blob/main/illustration.png)
# Acknowledgment
The GENPAT-IZASM Staff for our discussions aiming at managing arguments
# Author
Nicolas Radomski
