####################################################################################
###### Semi-automatic SingleLinkageToPairwiseThresholds.R script with RStudio ######
################################## R version 4.1.2 #################################
######################## RStudio version 2021.09.1 Build 372 #######################
####################################################################################

# INPUT: 
## Profiles.csv (fictitious input cgMLST profiles)
### a "sample" string must appear in the first row of the first column
### S stands for sample (i.e. rows) (e.g. n = 12)
### L stands for locus (i.e. columns) (e.g. n = 15)
### A stands for allele (i.e. data) (e.g. n= 180)
### 0 stands for missing data
# ARGUMENT SETTING
## argSLA: single linkage threshold A (Objective 1)
## argSLB: single linkage threshold B (Objective 1)
## argSLC: single linkage threshold C (Objective 1)
## argSLD: single linkage threshold D (Objective 1)
## argMNS: minimum number of samples to retain inside each single linkage cluster (Objective 1)
## argPSP: proportion of splitting (%) defining the optimization dataframes  (Objective 2)
# OUTPUT objective 1: Big single linkage clusters at four different single linkage thresholds (SL)
## Objective1-Dendrogram.pdf
## Objective1-PairwiseMatrix.csv
## Objective1-SingleLinkageAllClusters.csv
## Objective1-SingleLinkageAllClustersStatPerCluster.csv
## Objective1-SingleLinkageAllClustersStatPerSL.csv
## Objective1-SingleLinkageBigClusters.csv
## Objective1-SingleLinkageBigClustersStatPerCluster.csv
## Objective1-SingleLinkageBigClustersStatPerSL.csv
# OUTPUT objective 2: Objective 2: Random splitting of samples from each single linkage clusters
## Objective2-SingleLinkageBigClustersCorrelation.csv
## Objective2-SingleLinkageBigClustersCorrelationStatPerCluster.csv
## Objective2-SingleLinkageBigClustersCorrelationStatPerSL.csv
## Objective2-SingleLinkageBigClustersOptimization.csv
## Objective2-SingleLinkageBigClustersOptimizationStatPerCluster.csv
## Objective2-SingleLinkageBigClustersOptimizationStatPerSL.csv
# OUTPUT objective 3: Correlation between intra cluster pairwise differences (PW) and single linkage thresholds (SL)
## Objective3-PairwiseSingleLinkageCorrelationBoxplot.pdf
## Objective3-PairwiseSingleLinkageCorrelationBoxplot.tiff
## Objective3-PairwiseSingleLinkageCorrelationDataframe.csv
## Objective3-PairwiseSingleLinkageCorrelationLinearMax.pdf
## Objective3-PairwiseSingleLinkageCorrelationLinearMax.tiff
## Objective3-PairwiseSingleLinkageCorrelationLinear.pdf
## Objective3-PairwiseSingleLinkageCorrelationLinear.tiff
## Objective3-PairwiseSingleLinkageCorrelationStat.csv
# OUTPUT objective 4: Optimization of sensibility and specificity of pairwise differences between a single linkage clusters and the other
## Objective4-PairwiseSingleLinkageOptimizationDataframe.csv
## Objective4-PairwiseSingleLinkageOptimizationSLAThresholds.csv
## Objective4-PairwiseSingleLinkageOptimizationSLBThresholds.csv
## Objective4-PairwiseSingleLinkageOptimizationSLCThresholds.csv
## Objective4-PairwiseSingleLinkageOptimizationSLDThresholds.csv
## Objective4-PairwiseSingleLinkageoptimizationThresholdsBest.txt

####################################################################################
################################## Before to start #################################
####################################################################################

# install packages
install.packages("data.table")
install.packages("spaa")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("forcats")
install.packages("ggpmisc")
install.packages("pROC")

# set working directory for Linux and Mac
setwd("/home/IZSNT/n.radomski/Downloads")

# call library
library(data.table)
library(spaa)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggpmisc)
library(pROC)

# clean environment
rm(list=ls())
# reset graph devices
graphics.off()

# define arguments
## argSLA: single linkage threshold A (Objective 1)
argSLA <- 4
## argSLB: single linkage threshold B (Objective 1)
argSLB <- 7
## argSLC: single linkage threshold C (Objective 1)
argSLC <- 10
## argSLD: single linkage threshold D (Objective 1)
argSLD <- 13
## argMNS: minimum number of samples to retain inside each single linkage cluster (Objective 1)
argMNS <- 6
## argPSP: proportion of splitting (%) defining the optimization dataframes (Objective 2)
argPSP <- 50

####################################################################################
#################################### Objective 1 ###################################
####################################################################################

# prepare profiles
## read dataframe of profiles with missing data (i.e. Profiles.csv)
dffp <- read.table("Profiles.csv", dec = ".", header=TRUE, sep = ",", quote = "")
## make sure that each variable of the dataframe is a character
dffp <- data.frame(lapply(dffp, as.character))
## check nature of variables (must be character for each variable)
str(dffp)
## check dimension (i.e. [1] 40 36)
dim(dffp)
## replace missing data encoded 0 with NA
dfpNA <- dffp
dfpNA[ dfpNA == "0" ] <- NA
## remove columns harboring NAs (i.e. locus with missing data)
dfp <- dfpNA[ , colSums(is.na(dfpNA)) == 0]
## check dimension (i.e. [1] 40 31)
dim(dfp)

# transform profiles into matrix of pairwise differences
## transpose dataframe
tdfp <- transpose(dfp, keep.names = "locus", make.names = "sample")
## check dimension (i.e. [1] 30 41)
dim(tdfp)

## calculate pairwise differences between all vectors (i.e. independently of the sample amount)
### multiple for-loops adding results in an empty vector
v <- integer() # create empty vector v
for (i in tdfp[, 2:ncol(tdfp)]) { # first for loop from the column 2 of the dataframe
  for (j in tdfp[, 2:ncol(tdfp)]){ # second for loop from the column 2 of the dataframe
    #print(i) # for i checking
    #print(j) # for j checking
    output <- sum(!i == j) # count the number of FALSE (i.e. reverse (!) of paired vectors i versus j)
    v <- c(v, output) # over right the output vector into the empty vector v
  }
}
### check the output vector v
print(v)
### check size of the vector v (i.e. [1] 1600)
length(v)
### create a matrix adding the vector v by row independently of the sample amount
mat <- matrix(v, nrow = (ncol(tdfp)-1) , ncol = (ncol(tdfp)-1), byrow = TRUE)
### add row and column names of the matrix
#### keep the names of sample identifiers into a vector
samples = dfp$sample
#### add sample names to matrix rows
rownames(mat) <- samples
#### add sample names to matrix columns
colnames(mat) <- samples
### check class (i.e. [1] "matrix" "array")
class(mat)
### check pairwise differences
mat
### export the matrix into a csv file
write.csv(mat,file="Objective1-PairwiseMatrix.csv", row.names=TRUE, quote=FALSE)

# transform the matrix of pairwise differences into distobject and dataframe
## transform the matrix into a dist object
distobj = as.dist(mat)
## check class (i.e. [1] "dist")
class(distobj)
## check the dist object
distobj
## transform the dist object into a long dataframe of pairwise differences
dfl = dist2list(distobj)

# define 4 different single linkage thresholds (SL) (see arguments argSLA, argSLB, argSLC and argSLD in the script header)
# extract single linkage clusters
slc <- hclust(distobj, method = 'single')
## plot Dendrogram forcing the height of clusters (i.e. single linkage threshold h)
### open a pdf file
pdf("Objective1-Dendrogram.pdf", width=8, height=5)
### plot Dendrogram
plot(as.dendrogram(slc), 
     main=paste("Clusters of single linkage thresholds","\n",argSLA,"(A in blue),",argSLB,"(B in green),",argSLC,"(C in red) and",argSLD,"(D in violet)"),
     xlab= "Samples (S)",
     ylab="Single linkage thresholds (SL)",
     cex.lab=1,
     cex.axis=1,
     cex.main=1,)
rect.hclust(slc, h = argSLA, border = "blue")
rect.hclust(slc, h = argSLB, border = "green")
rect.hclust(slc, h = argSLC, border = "red")
rect.hclust(slc, h = argSLD, border = "violet")
### close the pdf file
dev.off()

# create clusters of many samples according to different single linkage threshold h (i.e. height of clusters)
## retrieve all clusters of samples for different single linkage threshold h
clusters <- cutree(slc, h = c(argSLA,argSLB,argSLC,argSLD))
## check class (i.e. [1] "matrix" "array")
class(clusters)
## remove clusters <= n samples
### transform as dataframe
dfc <- as.data.frame(clusters)
### rename all variable
  colnames(dfc) <- c("A","B","C","D")
  ### add prefix to each variable
  colnames(dfc) <- paste("SL",colnames(dfc),sep="")
  #### add a variable of sample
  dfc$sample <- samples
### add C character before each cluster number of each cutoff
dfc$SLA <- sub("^", "C", dfc$SLA)
dfc$SLB <- sub("^", "C", dfc$SLB)
dfc$SLC <- sub("^", "C", dfc$SLC)
dfc$SLD <- sub("^", "C", dfc$SLD)
#### make sure all variable are characters
dfc <- data.frame(lapply(dfc, as.character))

# measure and identify big cluster
## count samples of each clusters and rename variables
### for single linkage threshold A
SLA <- count(dfc, (dfc)[1])
colnames(SLA) <- c("cluster","sample")
### for single linkage threshold B
SLB <- count(dfc, (dfc)[2])
colnames(SLB) <- c("cluster","sample")
### for single linkage threshold C
SLC <- count(dfc, (dfc)[3])
colnames(SLC) <- c("cluster","sample")
### for single linkage threshold D
SLD <- count(dfc, (dfc)[4])
colnames(SLD) <- c("cluster","sample")
## define minimum number of samples (MNS) to retain inside each single linkage cluster (e.g. >= MNS samples) (see argument argMNS in the script header)
## subset clusters >= MNS
SLAhigher <- subset(SLA, sample >= argMNS)
SLBhigher <- subset(SLB, sample >= argMNS)
SLChigher <- subset(SLC, sample >= argMNS)
SLDhigher <- subset(SLD, sample >= argMNS)
## flag clusters with higher (or equal) or lower (strictly) than MNS into the whole dataframe of clusters
dfc$SLAsize <- ifelse(dfc$SLA %in% SLAhigher$cluster, "higher", "lower")
dfc$SLBsize <- ifelse(dfc$SLB %in% SLBhigher$cluster, "higher", "lower")
dfc$SLCsize <- ifelse(dfc$SLC %in% SLChigher$cluster, "higher", "lower")
dfc$SLDsize <- ifelse(dfc$SLD %in% SLDhigher$cluster, "higher", "lower")

# retrieve all clusters for each single linkage threshold, rename and add variables
## for single linkage threshold A
SLAunfiltered <- subset(dfc,dfc$SLAsize %in% c("higher","lower"))
SLAunfiltered <- subset(SLAunfiltered, select = c(SLA, sample))
names(SLAunfiltered)[names(SLAunfiltered) == "SLA"] <- "cluster"
SLAunfiltered$SL <- "A"
## for single linkage threshold B
SLBunfiltered <- subset(dfc,dfc$SLBsize %in% c("higher","lower"))
SLBunfiltered <- subset(SLBunfiltered, select = c(SLB, sample))
names(SLBunfiltered)[names(SLBunfiltered) == "SLB"] <- "cluster"
SLBunfiltered$SL <- "B"
## for single linkage threshold C
SLCunfiltered <- subset(dfc,dfc$SLCsize %in% c("higher","lower"))
SLCunfiltered <- subset(SLCunfiltered, select = c(SLC, sample))
names(SLCunfiltered)[names(SLCunfiltered) == "SLC"] <- "cluster"
SLCunfiltered$SL <- "C"
## for single linkage threshold D
SLDunfiltered <- subset(dfc,dfc$SLDsize %in% c("higher","lower"))
SLDunfiltered <- subset(SLDunfiltered, select = c(SLD, sample))
names(SLDunfiltered)[names(SLDunfiltered) == "SLD"] <- "cluster"
SLDunfiltered$SL <- "D"
## combine and export dataframes into a csv files
SLall <- rbind(SLAunfiltered, SLBunfiltered, SLCunfiltered, SLDunfiltered)
write.csv(SLall,file="Objective1-SingleLinkageAllClusters.csv", row.names=FALSE, quote=FALSE)
## synthesis number of samples and export in a csv file
### per single linkage threshold (SL)
SLallStatPerSL <- SLall %>% group_by(SL) %>% tally()
colnames(SLallStatPerSL) <- c("SL","sample")
write.csv(SLallStatPerSL,file="Objective1-SingleLinkageAllClustersStatPerSL.csv", row.names=FALSE, quote=FALSE)
### per cluster
SLallStatPerCluster <- SLall %>% group_by(cluster) %>% tally()
colnames(SLallStatPerCluster) <- c("cluster","sample")
write.csv(SLallStatPerCluster,file="Objective1-SingleLinkageAllClustersStatPerCluster.csv", row.names=FALSE, quote=FALSE)

# retrieve big clusters (>= MNS) for each single linkage threshold, rename and add variables
## for single linkage threshold A
SLAfiltered <- subset(dfc,dfc$SLAsize %in% c("higher"))
SLAfiltered <- subset(SLAfiltered, select = c(SLA, sample))
names(SLAfiltered)[names(SLAfiltered) == "SLA"] <- "cluster"
SLAfiltered$SL <- "A"
## for single linkage threshold B
SLBfiltered <- subset(dfc,dfc$SLBsize %in% c("higher"))
SLBfiltered <- subset(SLBfiltered, select = c(SLB, sample))
names(SLBfiltered)[names(SLBfiltered) == "SLB"] <- "cluster"
SLBfiltered$SL <- "B"
## for single linkage threshold C
SLCfiltered <- subset(dfc,dfc$SLCsize %in% c("higher"))
SLCfiltered <- subset(SLCfiltered, select = c(SLC, sample))
names(SLCfiltered)[names(SLCfiltered) == "SLC"] <- "cluster"
SLCfiltered$SL <- "C"
## for single linkage threshold D
SLDfiltered <- subset(dfc,dfc$SLDsize %in% c("higher"))
SLDfiltered <- subset(SLDfiltered, select = c(SLD, sample))
names(SLDfiltered)[names(SLDfiltered) == "SLD"] <- "cluster"
SLDfiltered$SL <- "D"
## combine and export dataframes into a csv files
SLbig <- rbind(SLAfiltered, SLBfiltered, SLCfiltered, SLDfiltered)
write.csv(SLbig,file="Objective1-SingleLinkageBigClusters.csv", row.names=FALSE, quote=FALSE)
## synthesis number of samples and export in a csv file
### per single linkage threshold (SL)
SLbigStatPerSL <- SLbig %>% group_by(SL) %>% tally()
colnames(SLbigStatPerSL) <- c("SL","sample")
write.csv(SLbigStatPerSL,file="Objective1-SingleLinkageBigClustersStatPerSL.csv", row.names=FALSE, quote=FALSE)
### per cluster
SLbigStatPerCluster <- SLbig %>% group_by(cluster) %>% tally()
colnames(SLbigStatPerCluster) <- c("cluster","sample")
write.csv(SLbigStatPerCluster,file="Objective1-SingleLinkageBigClustersStatPerCluster.csv", row.names=FALSE, quote=FALSE)

# checking
## number of samples including all single linkage thresholds
### original ([1] 160)
#### first check
dim(SLall)
#### second check
sum(SLallStatPerSL$sample)
### after filtering ([1] 122)
#### first check
dim(SLbig)
#### second check
sum(SLbigStatPerSL$sample)
## number of clusters including all single linkage thresholds
### original ([1] 28)
#### first check
length(SLA$cluster) + length(SLB$cluster) + length(SLC$cluster) + length(SLD$cluster)
#### second check
dA <- subset(SLall, SLall$SL %in% c("A"))
dB <- subset(SLall, SLall$SL %in% c("B"))
dC <- subset(SLall, SLall$SL %in% c("C"))
dD <- subset(SLall, SLall$SL %in% c("D"))
ddA <- count(dA, dA$cluster)
ddB <- count(dB, dB$cluster)
ddC <- count(dC, dC$cluster)
ddD <- count(dC, dD$cluster)
dim(ddA) + dim(ddB) + dim(ddC) + dim(ddD)
### after filtering ([1] 11)
#### first check
dA <- count(SLAfiltered, SLAfiltered$cluster)
dB <- count(SLBfiltered, SLBfiltered$cluster)
dC <- count(SLCfiltered, SLCfiltered$cluster)
dD <- count(SLDfiltered, SLDfiltered$cluster)
dim(dA) + dim(dB) + dim(dC) + dim(dD)
#### second check
dA <- subset(SLbig, SLbig$SL %in% c("A"))
dB <- subset(SLbig, SLbig$SL %in% c("B"))
dC <- subset(SLbig, SLbig$SL %in% c("C"))
dD <- subset(SLbig, SLbig$SL %in% c("D"))
ddA <- count(dA, dA$cluster)
ddB <- count(dB, dB$cluster)
ddC <- count(dC, dC$cluster)
ddD <- count(dD, dD$cluster)
dim(ddA) + dim(ddB) + dim(ddC) + dim(ddD)

####################################################################################
#################################### Objective 2 ###################################
####################################################################################

# proportion based random splitting of samples in each clusters
## define proportion of splitting (%) defining the optimization dataframes (see argument argPSP in the script header)
## define proportion of splitting (number between 0 and 1)
psn <- argPSP/100
## split randomly by cluster
### for single linkage threshold A
SLAoptimization <- SLAfiltered %>% group_by(cluster) %>% sample_frac(psn)
SLAoptimization <- as.data.frame(SLAoptimization)
### for single linkage threshold B
SLBoptimization <- SLBfiltered %>% group_by(cluster) %>% sample_frac(psn)
SLBoptimization <- as.data.frame(SLBoptimization)
### for single linkage threshold C
SLCoptimization <- SLCfiltered %>% group_by(cluster) %>% sample_frac(psn)
SLCoptimization <- as.data.frame(SLCoptimization)
### for single linkage threshold D
SLDoptimization <- SLDfiltered %>% group_by(cluster) %>% sample_frac(psn)
SLDoptimization <- as.data.frame(SLDoptimization)
## combine and export dataframes into a csv files
SLoptimization <- rbind(SLAoptimization, SLBoptimization, SLCoptimization, SLDoptimization)
write.csv(SLoptimization,file="Objective2-SingleLinkageBigClustersOptimization.csv", row.names=FALSE, quote=FALSE)
## synthesis number of samples and export in a csv file
### per single linkage threshold (SL)
SLoptimizationStatPerSL <- SLoptimization %>% group_by(SL) %>% tally()
colnames(SLoptimizationStatPerSL) <- c("SL","sample")
write.csv(SLoptimizationStatPerSL,file="Objective2-SingleLinkageBigClustersOptimizationStatPerSL.csv", row.names=FALSE, quote=FALSE)
### per cluster
SLoptimizationStatPerCluster <- SLoptimization %>% group_by(cluster) %>% tally()
colnames(SLoptimizationStatPerCluster) <- c("cluster","sample")
write.csv(SLoptimizationStatPerCluster,file="Objective2-SingleLinkageBigClustersOptimizationStatPerCluster.csv", row.names=FALSE, quote=FALSE)

# retrieve remaining samples defining the correlation dataframes
## perform differences between dataframes
### for single linkage threshold A
SLAcorrelation <- setdiff(SLAfiltered, SLAoptimization)
### for single linkage threshold B
SLBcorrelation <- setdiff(SLBfiltered, SLBoptimization)
### for single linkage threshold C
SLCcorrelation <- setdiff(SLCfiltered, SLCoptimization)
### for single linkage threshold D
SLDcorrelation <- setdiff(SLDfiltered, SLDoptimization)
## combine and export dataframes into a csv files
SLcorrelation <- rbind(SLAcorrelation, SLBcorrelation, SLCcorrelation, SLDcorrelation)
write.csv(SLcorrelation,file="Objective2-SingleLinkageBigClustersCorrelation.csv", row.names=FALSE, quote=FALSE)
## synthesis number of samples and export in a csv file
### per single linkage threshold (SL)
SLcorrelationStatPerSL <- SLcorrelation %>% group_by(SL) %>% tally()
colnames(SLcorrelationStatPerSL) <- c("SL","sample")
write.csv(SLcorrelationStatPerSL,file="Objective2-SingleLinkageBigClustersCorrelationStatPerSL.csv", row.names=FALSE, quote=FALSE)
### per cluster
SLcorrelationStatPerCluster <- SLcorrelation %>% group_by(cluster) %>% tally()
colnames(SLcorrelationStatPerCluster) <- c("cluster","sample")
write.csv(SLcorrelationStatPerCluster,file="Objective2-SingleLinkageBigClustersCorrelationStatPerCluster.csv", row.names=FALSE, quote=FALSE)

# checking
## checking of number of samples including all single linkage thresholds
### optimization dataframe ([1] 47)
#### first check
sum(SLoptimizationStatPerSL$sample)
#### second check
dim(SLoptimization)
### correlation dataframe ([1] 47)
#### first check
sum(SLcorrelationStatPerSL$sample)
#### second check
dim(SLcorrelation)

## checking of cluster of samples including all single linkage thresholds
### optimization dataframe ([1] 11)
#### first check
dA <- count(SLAoptimization, SLAoptimization$cluster)
dB <- count(SLBoptimization, SLBoptimization$cluster)
dC <- count(SLCoptimization, SLCoptimization$cluster)
dD <- count(SLDoptimization, SLDoptimization$cluster)
dim(dA) + dim(dB) + dim(dC) + dim(dD)
#### second check
dA <- subset(SLoptimization, SLoptimization$SL %in% c("A"))
dB <- subset(SLoptimization, SLoptimization$SL %in% c("B"))
dC <- subset(SLoptimization, SLoptimization$SL %in% c("C"))
dD <- subset(SLoptimization, SLoptimization$SL %in% c("D"))
ddA <- count(dA, dA$cluster)
ddB <- count(dB, dB$cluster)
ddC <- count(dC, dC$cluster)
ddD <- count(dD, dD$cluster)
dim(ddA) + dim(ddB) + dim(ddC) + dim(ddD)
### correlation dataframe ([1] 11)
#### first check
dA <- count(SLAcorrelation, SLAcorrelation$cluster)
dB <- count(SLBcorrelation, SLBcorrelation$cluster)
dC <- count(SLCcorrelation, SLCcorrelation$cluster)
dD <- count(SLDcorrelation, SLDcorrelation$cluster)
dim(dA) + dim(dB) + dim(dC) + dim(dD)
#### second check
dA <- subset(SLcorrelation, SLcorrelation$SL %in% c("A"))
dB <- subset(SLcorrelation, SLcorrelation$SL %in% c("B"))
dC <- subset(SLcorrelation, SLcorrelation$SL %in% c("C"))
dD <- subset(SLcorrelation, SLcorrelation$SL %in% c("D"))
ddA <- count(dA, dA$cluster)
ddB <- count(dB, dB$cluster)
ddC <- count(dC, dC$cluster)
ddD <- count(dD, dD$cluster)
dim(ddA) + dim(ddB) + dim(ddC) + dim(ddD)

####################################################################################
#################################### Objective 3 ###################################
####################################################################################

# transform the dist object into a long dataframe of pairwise differences
## pass from dist object to long dataframe
dfl = dist2list(distobj)
## check class (i.e. [1] "data.frame")
class(dfl)
## check dimension (i.e. [1] 1600    3)
dim(dfl)
## check names of variables
colnames(dfl)
## rename variables
names(dfl)[names(dfl) == "col"] <- "firstsample"
names(dfl)[names(dfl) == "row"] <- "secondsample"
names(dfl)[names(dfl) == "value"] <- "differences"

# retrieve pairwise differences belonging to correlation dataframe
## for single linkage threshold A
### add single linkage clusters into the correlation dataframe
#### rename variable
names(SLAcorrelation)[names(SLAcorrelation) == "sample"] <- "firstsample"
#### joint left correlation dataframe and dataframe of single linkage clusters
dflSLAcorrelation <- left_join(dfl, SLAcorrelation, by = "firstsample", keep = FALSE)
#### remove rows harboring NAs (i.e. pairwise not involved in single linkage clusters)
dflSLAcorrelation = dflSLAcorrelation %>% drop_na()
#### rename variable
names(dflSLAcorrelation)[names(dflSLAcorrelation) == "cluster"] <- "firstcluster"
names(dflSLAcorrelation)[names(dflSLAcorrelation) == "SL"] <- "firstSL"
#### reiterate with the second samples
names(SLAcorrelation)[names(SLAcorrelation) == "firstsample"] <- "secondsample"
dflSLAcorrelation <- left_join(dflSLAcorrelation, SLAcorrelation, by = "secondsample", keep = FALSE)
dflSLAcorrelation = dflSLAcorrelation %>% drop_na()
names(dflSLAcorrelation)[names(dflSLAcorrelation) == "cluster"] <- "secondcluster"
names(dflSLAcorrelation)[names(dflSLAcorrelation) == "SL"] <- "secondSL"
## for single linkage threshold B
names(SLBcorrelation)[names(SLBcorrelation) == "sample"] <- "firstsample"
dflSLBcorrelation <- left_join(dfl, SLBcorrelation, by = "firstsample", keep = FALSE)
dflSLBcorrelation = dflSLBcorrelation %>% drop_na()
names(dflSLBcorrelation)[names(dflSLBcorrelation) == "cluster"] <- "firstcluster"
names(dflSLBcorrelation)[names(dflSLBcorrelation) == "SL"] <- "firstSL"
names(SLBcorrelation)[names(SLBcorrelation) == "firstsample"] <- "secondsample"
dflSLBcorrelation <- left_join(dflSLBcorrelation, SLBcorrelation, by = "secondsample", keep = FALSE)
dflSLBcorrelation = dflSLBcorrelation %>% drop_na()
names(dflSLBcorrelation)[names(dflSLBcorrelation) == "cluster"] <- "secondcluster"
names(dflSLBcorrelation)[names(dflSLBcorrelation) == "SL"] <- "secondSL"
## for single linkage threshold C
names(SLCcorrelation)[names(SLCcorrelation) == "sample"] <- "firstsample"
dflSLCcorrelation <- left_join(dfl, SLCcorrelation, by = "firstsample", keep = FALSE)
dflSLCcorrelation = dflSLCcorrelation %>% drop_na()
names(dflSLCcorrelation)[names(dflSLCcorrelation) == "cluster"] <- "firstcluster"
names(dflSLCcorrelation)[names(dflSLCcorrelation) == "SL"] <- "firstSL"
names(SLCcorrelation)[names(SLCcorrelation) == "firstsample"] <- "secondsample"
dflSLCcorrelation <- left_join(dflSLCcorrelation, SLCcorrelation, by = "secondsample", keep = FALSE)
dflSLCcorrelation = dflSLCcorrelation %>% drop_na()
names(dflSLCcorrelation)[names(dflSLCcorrelation) == "cluster"] <- "secondcluster"
names(dflSLCcorrelation)[names(dflSLCcorrelation) == "SL"] <- "secondSL"
## for single linkage threshold D
names(SLDcorrelation)[names(SLDcorrelation) == "sample"] <- "firstsample"
dflSLDcorrelation <- left_join(dfl, SLDcorrelation, by = "firstsample", keep = FALSE)
dflSLDcorrelation = dflSLDcorrelation %>% drop_na()
names(dflSLDcorrelation)[names(dflSLDcorrelation) == "cluster"] <- "firstcluster"
names(dflSLDcorrelation)[names(dflSLDcorrelation) == "SL"] <- "firstSL"
names(SLDcorrelation)[names(SLDcorrelation) == "firstsample"] <- "secondsample"
dflSLDcorrelation <- left_join(dflSLDcorrelation, SLDcorrelation, by = "secondsample", keep = FALSE)
dflSLDcorrelation = dflSLDcorrelation %>% drop_na()
names(dflSLDcorrelation)[names(dflSLDcorrelation) == "cluster"] <- "secondcluster"
names(dflSLDcorrelation)[names(dflSLDcorrelation) == "SL"] <- "secondSL"

# flag the pairwise differences from intra single linkage clusters
dflSLAcorrelation$cluster <- ifelse((dflSLAcorrelation$firstcluster == dflSLAcorrelation$secondcluster), "intra", "extra")
dflSLBcorrelation$cluster <- ifelse((dflSLBcorrelation$firstcluster == dflSLBcorrelation$secondcluster), "intra", "extra")
dflSLCcorrelation$cluster <- ifelse((dflSLCcorrelation$firstcluster == dflSLCcorrelation$secondcluster), "intra", "extra")
dflSLDcorrelation$cluster <- ifelse((dflSLDcorrelation$firstcluster == dflSLDcorrelation$secondcluster), "intra", "extra")

# subset intra single linkage clusters
SLAPWcorrelationIntra <- subset(dflSLAcorrelation,dflSLAcorrelation$cluster %in% c("intra"))
SLBPWcorrelationIntra <- subset(dflSLBcorrelation,dflSLBcorrelation$cluster %in% c("intra"))
SLCPWcorrelationIntra <- subset(dflSLCcorrelation,dflSLCcorrelation$cluster %in% c("intra"))
SLDPWcorrelationIntra <- subset(dflSLDcorrelation,dflSLDcorrelation$cluster %in% c("intra"))

# estimate variability of pairwise differences for each single linkage cluster
## for single linkage threshold A
### minimum
SLAmin <- aggregate(x = SLAPWcorrelationIntra$differences, by= list(SLAPWcorrelationIntra$firstcluster), FUN = min)
names(SLAmin)[names(SLAmin) == "Group.1"] <- "SL_cluster"
names(SLAmin)[names(SLAmin) == "x"] <- "PW_min"
### mean
SLAmean <- aggregate(x = SLAPWcorrelationIntra$differences, by= list(SLAPWcorrelationIntra$firstcluster), FUN = mean)
names(SLAmean)[names(SLAmean) == "Group.1"] <- "SL_cluster"
names(SLAmean)[names(SLAmean) == "x"] <- "PW_mean"
### standard deviation
SLAsd <- aggregate(x = SLAPWcorrelationIntra$differences, by= list(SLAPWcorrelationIntra$firstcluster), FUN = sd)
names(SLAsd)[names(SLAsd) == "Group.1"] <- "SL_cluster"
names(SLAsd)[names(SLAsd) == "x"] <- "PW_sd"
### maximum
SLAmax <- aggregate(x = SLAPWcorrelationIntra$differences, by= list(SLAPWcorrelationIntra$firstcluster), FUN = max)
names(SLAmax)[names(SLAmax) == "Group.1"] <- "SL_cluster"
names(SLAmax)[names(SLAmax) == "x"] <- "PW_max"
### joint dataframes
SLAminmean <- inner_join(SLAmin, SLAmean, by = "SL_cluster")
SLAminmeansd <- inner_join(SLAminmean, SLAsd, by = "SL_cluster")
SLAstat <- inner_join(SLAminmeansd, SLAmax, by = "SL_cluster")
## for single linkage threshold B
SLBmin <- aggregate(x = SLBPWcorrelationIntra$differences, by= list(SLBPWcorrelationIntra$firstcluster), FUN = min)
names(SLBmin)[names(SLBmin) == "Group.1"] <- "SL_cluster"
names(SLBmin)[names(SLBmin) == "x"] <- "PW_min"
SLBmean <- aggregate(x = SLBPWcorrelationIntra$differences, by= list(SLBPWcorrelationIntra$firstcluster), FUN = mean)
names(SLBmean)[names(SLBmean) == "Group.1"] <- "SL_cluster"
names(SLBmean)[names(SLBmean) == "x"] <- "PW_mean"
SLBsd <- aggregate(x = SLBPWcorrelationIntra$differences, by= list(SLBPWcorrelationIntra$firstcluster), FUN = sd)
names(SLBsd)[names(SLBsd) == "Group.1"] <- "SL_cluster"
names(SLBsd)[names(SLBsd) == "x"] <- "PW_sd"
SLBmax <- aggregate(x = SLBPWcorrelationIntra$differences, by= list(SLBPWcorrelationIntra$firstcluster), FUN = max)
names(SLBmax)[names(SLBmax) == "Group.1"] <- "SL_cluster"
names(SLBmax)[names(SLBmax) == "x"] <- "PW_max"
SLBminmean <- inner_join(SLBmin, SLBmean, by = "SL_cluster")
SLBminmeansd <- inner_join(SLBminmean, SLBsd, by = "SL_cluster")
SLBstat <- inner_join(SLBminmeansd, SLBmax, by = "SL_cluster")
## for single linkage threshold C
SLCmin <- aggregate(x = SLCPWcorrelationIntra$differences, by= list(SLCPWcorrelationIntra$firstcluster), FUN = min)
names(SLCmin)[names(SLCmin) == "Group.1"] <- "SL_cluster"
names(SLCmin)[names(SLCmin) == "x"] <- "PW_min"
SLCmean <- aggregate(x = SLCPWcorrelationIntra$differences, by= list(SLCPWcorrelationIntra$firstcluster), FUN = mean)
names(SLCmean)[names(SLCmean) == "Group.1"] <- "SL_cluster"
names(SLCmean)[names(SLCmean) == "x"] <- "PW_mean"
SLCsd <- aggregate(x = SLCPWcorrelationIntra$differences, by= list(SLCPWcorrelationIntra$firstcluster), FUN = sd)
names(SLCsd)[names(SLCsd) == "Group.1"] <- "SL_cluster"
names(SLCsd)[names(SLCsd) == "x"] <- "PW_sd"
SLCmax <- aggregate(x = SLCPWcorrelationIntra$differences, by= list(SLCPWcorrelationIntra$firstcluster), FUN = max)
names(SLCmax)[names(SLCmax) == "Group.1"] <- "SL_cluster"
names(SLCmax)[names(SLCmax) == "x"] <- "PW_max"
SLCminmean <- inner_join(SLCmin, SLCmean, by = "SL_cluster")
SLCminmeansd <- inner_join(SLCminmean, SLCsd, by = "SL_cluster")
SLCstat <- inner_join(SLCminmeansd, SLCmax, by = "SL_cluster")
## for single linkage threshold D
SLDmin <- aggregate(x = SLDPWcorrelationIntra$differences, by= list(SLDPWcorrelationIntra$firstcluster), FUN = min)
names(SLDmin)[names(SLDmin) == "Group.1"] <- "SL_cluster"
names(SLDmin)[names(SLDmin) == "x"] <- "PW_min"
SLDmean <- aggregate(x = SLDPWcorrelationIntra$differences, by= list(SLDPWcorrelationIntra$firstcluster), FUN = mean)
names(SLDmean)[names(SLDmean) == "Group.1"] <- "SL_cluster"
names(SLDmean)[names(SLDmean) == "x"] <- "PW_mean"
SLDsd <- aggregate(x = SLDPWcorrelationIntra$differences, by= list(SLDPWcorrelationIntra$firstcluster), FUN = sd)
names(SLDsd)[names(SLDsd) == "Group.1"] <- "SL_cluster"
names(SLDsd)[names(SLDsd) == "x"] <- "PW_sd"
SLDmax <- aggregate(x = SLDPWcorrelationIntra$differences, by= list(SLDPWcorrelationIntra$firstcluster), FUN = max)
names(SLDmax)[names(SLDmax) == "Group.1"] <- "SL_cluster"
names(SLDmax)[names(SLDmax) == "x"] <- "PW_max"
SLDminmean <- inner_join(SLDmin, SLDmean, by = "SL_cluster")
SLDminmeansd <- inner_join(SLDminmean, SLDsd, by = "SL_cluster")
SLDstat <- inner_join(SLDminmeansd, SLDmax, by = "SL_cluster")

## add variable, combine and export the dataframe into a csv file
SLAstat$SL <- "A"
SLBstat$SL <- "B"
SLCstat$SL <- "C"
SLDstat$SL <- "D"
SLPWcorrelationIntrastat <- rbind(SLAstat, SLBstat, SLCstat, SLDstat)
write.csv(SLPWcorrelationIntrastat,file="Objective3-PairwiseSingleLinkageCorrelationStat.csv", row.names=FALSE, quote=FALSE)

# graphical representation of pairwise differences for each single linkage cluster
## compile data into a long dataframe
### make copies of dataframes
SLAPWcorrelationIntra <- SLAPWcorrelationIntra
SLBPWcorrelationIntra <- SLBPWcorrelationIntra
SLCPWcorrelationIntra <- SLCPWcorrelationIntra
SLDPWcorrelationIntra <- SLDPWcorrelationIntra
### rename variables
names(SLAPWcorrelationIntra)[names(SLAPWcorrelationIntra) == "firstSL"] <- "SLthreshold"
names(SLBPWcorrelationIntra)[names(SLBPWcorrelationIntra) == "firstSL"] <- "SLthreshold"
names(SLCPWcorrelationIntra)[names(SLCPWcorrelationIntra) == "firstSL"] <- "SLthreshold"
names(SLDPWcorrelationIntra)[names(SLDPWcorrelationIntra) == "firstSL"] <- "SLthreshold"
names(SLAPWcorrelationIntra)[names(SLAPWcorrelationIntra) == "firstcluster"] <- "SLcluster"
names(SLBPWcorrelationIntra)[names(SLBPWcorrelationIntra) == "firstcluster"] <- "SLcluster"
names(SLCPWcorrelationIntra)[names(SLCPWcorrelationIntra) == "firstcluster"] <- "SLcluster"
names(SLDPWcorrelationIntra)[names(SLDPWcorrelationIntra) == "firstcluster"] <- "SLcluster"
names(SLAPWcorrelationIntra)[names(SLAPWcorrelationIntra) == "differences"] <- "PWdifferences"
names(SLBPWcorrelationIntra)[names(SLBPWcorrelationIntra) == "differences"] <- "PWdifferences"
names(SLCPWcorrelationIntra)[names(SLCPWcorrelationIntra) == "differences"] <- "PWdifferences"
names(SLDPWcorrelationIntra)[names(SLDPWcorrelationIntra) == "differences"] <- "PWdifferences"
### remove useless variables
SLAPWcorrelationIntra$secondcluster <- NULL
SLAPWcorrelationIntra$secondSL <- NULL
SLAPWcorrelationIntra$cluster <- NULL
SLBPWcorrelationIntra$secondcluster <- NULL
SLBPWcorrelationIntra$secondSL <- NULL
SLBPWcorrelationIntra$cluster <- NULL
SLCPWcorrelationIntra$secondcluster <- NULL
SLCPWcorrelationIntra$secondSL <- NULL
SLCPWcorrelationIntra$cluster <- NULL
SLDPWcorrelationIntra$secondcluster <- NULL
SLDPWcorrelationIntra$secondSL <- NULL
SLDPWcorrelationIntra$cluster <- NULL
### combine dataframes by rows
SLPWcorrelationIntra <- rbind (SLAPWcorrelationIntra, SLBPWcorrelationIntra, SLCPWcorrelationIntra, SLDPWcorrelationIntra)
### double the variable SLthreshold
SLPWcorrelationIntra$chrSLthreshold <- SLPWcorrelationIntra$SLthreshold
SLPWcorrelationIntra$numSLthreshold <- SLPWcorrelationIntra$SLthreshold
SLPWcorrelationIntra$SLthreshold <- NULL
### rename levels of variables
SLPWcorrelationIntra$numSLthreshold <- factor(SLPWcorrelationIntra$numSLthreshold)
levels(SLPWcorrelationIntra$numSLthreshold)
levels(SLPWcorrelationIntra$numSLthreshold) <- c(argSLA, argSLB, argSLC, argSLD)
levels(SLPWcorrelationIntra$numSLthreshold)
### transform numerics as characters
SLPWcorrelationIntra$numSLthreshold <- as.character(SLPWcorrelationIntra$numSLthreshold)
## export the dataframe into csv files
write.csv(SLPWcorrelationIntra,file="Objective3-PairwiseSingleLinkageCorrelationDataframe.csv", row.names=FALSE, quote=FALSE)

## plot boxplots
### reorganize levels of variables
SLPWcorrelationIntra$chrSLthreshold <- factor(SLPWcorrelationIntra$chrSLthreshold)
levels(SLPWcorrelationIntra$chrSLthreshold)
SLPWcorrelationIntra$chrSLthreshold <- factor(SLPWcorrelationIntra$chrSLthreshold, levels=c("A", "B", "C", "D"))
### reorganize levels of variables of character by frequency
SLPWcorrelationIntra$SLcluster <- factor(SLPWcorrelationIntra$SLcluster)
levels(SLPWcorrelationIntra$SLcluster)
SLPWcorrelationIntra$SLcluster <- fct_rev(fct_infreq(SLPWcorrelationIntra$SLcluster))
levels(SLPWcorrelationIntra$SLcluster)
### plot
p = ggplot(data = SLPWcorrelationIntra, aes(x = SLcluster, y = PWdifferences)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(height = 0.1), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  ggtitle("Boxplot of pairwise differences intra single linkage clusters \n from single linkage thresholds") +
  scale_y_continuous(name = "Pairwise differences of single linkage clusters") +
  scale_x_discrete(name = paste("Single linkage clusters (e.g. C1) of single linkage thresholds \n " ,argSLA,"(i.e. A argument),",argSLB,"(i.e. B argument),",argSLC,"(i.e. C argument) and",argSLD,"(i.e. D argument)")) +
  theme(plot.title = element_text(color="black", size=12, face="bold.italic", hjust = 0.5),
        axis.text.x = element_text (color = "#000000", size = 8, angle = 90, vjust = 0.5),
        axis.title.x = element_text (size = 10),
        axis.title.y = element_text (size = 10),
        strip.text.x = element_text(size=16, face = "bold"),
        strip.text.y = element_text(size=16, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  facet_grid(~ chrSLthreshold)
p
plot(p)
ggsave("Objective3-PairwiseSingleLinkageCorrelationBoxplot.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("Objective3-PairwiseSingleLinkageCorrelationBoxplot.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

## plot correlation
### transform factor as characters, then character as numerics 
class(SLPWcorrelationIntra$numSLthreshold)
SLPWcorrelationIntra$numSLthreshold <- as.numeric(SLPWcorrelationIntra$numSLthreshold)
### plot
my.formula <- y ~ x
p = ggplot(SLPWcorrelationIntra, aes(x = numSLthreshold, y = PWdifferences)) +
  theme_light(base_size = 16) +
  geom_point(shape=18, color = "#000000") +
  geom_smooth(size = 0.5, color="#000000", fill = "#A9A9A9", method=lm, linetype="dashed", se=TRUE, formula = my.formula) +
  ggtitle("Correlation between pairwise differences of single linkage clusters \n and single linkage thresholds \n ") +
  scale_y_continuous(name = "Pairwise differences intra single linkage clusters") +
  scale_x_continuous(name = paste("Single linkage clusters of single linkage thresholds \n " ,argSLA,"(i.e. A argument),",argSLB,"(i.e. B argument),",argSLC,"(i.e. C argument) and",argSLD,"(i.e. D argument)"), breaks = c(argSLA,argSLB,argSLC,argSLD)) +
  theme(plot.title = element_text(color="black", size=12, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text (size = 10),
        axis.title.y = element_text (size = 10)) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,label.y = "top", label.x = "right", rr.digits = 3, coef.digits = 3, f.digits = 6)
p
plot(p)
ggsave("Objective3-PairwiseSingleLinkageCorrelationLinear.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("Objective3-PairwiseSingleLinkageCorrelationLinear.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

## plot correlation for max values of pairwise differences
### retrieve max values of pairwise differences for each single linkage threshold
SLPWcorrelationIntramax <- SLPWcorrelationIntra
SLPWcorrelationIntramax <- SLPWcorrelationIntramax %>% 
  group_by(chrSLthreshold) %>% 
  filter(PWdifferences==max(PWdifferences))
### plot
my.formula <- y ~ x
p = ggplot(SLPWcorrelationIntramax, aes(x = numSLthreshold, y = PWdifferences)) +
  theme_light(base_size = 16) +
  geom_point(shape=18, color = "#000000") +
  geom_smooth(size = 0.5, color="#000000", fill = "#A9A9A9", method=lm, linetype="dashed", se=TRUE, formula = my.formula) +
  ggtitle("Correlation between maximum pairwise differences of single linkage clusters \n and single linkage thresholds \n ") +
  scale_y_continuous(name = "Maximum pairwise differences intra single linkage clusters") +
  scale_x_continuous(name = paste("Single linkage clusters of single linkage thresholds \n " ,argSLA,"(i.e. A argument),",argSLB,"(i.e. B argument),",argSLC,"(i.e. C argument) and",argSLD,"(i.e. D argument)"), breaks = c(argSLA,argSLB,argSLC,argSLD)) +
  theme(plot.title = element_text(color="black", size=12, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text (size = 10),
        axis.title.y = element_text (size = 10)) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,label.y = "top", label.x = "right", rr.digits = 3, coef.digits = 3, f.digits = 6)
p
plot(p)
ggsave("Objective3-PairwiseSingleLinkageCorrelationLinearMax.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("Objective3-PairwiseSingleLinkageCorrelationLinearMax.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

####################################################################################
#################################### Objective 4 ###################################
####################################################################################

# make copy of optimization dataframes (from random splitting)
SLAoptimizationBIS <- SLAoptimization
SLBoptimizationBIS <- SLBoptimization
SLCoptimizationBIS <- SLCoptimization
SLDoptimizationBIS <- SLDoptimization

# retrieve pairwise differences belonging to optimizaton dataframe
## for single linkage threshold A
names(SLAoptimizationBIS)[names(SLAoptimizationBIS) == "sample"] <- "firstsample"
dflSLAoptimization <- left_join(dfl, SLAoptimizationBIS, by = "firstsample", keep = FALSE)
dflSLAoptimization = dflSLAoptimization %>% drop_na()
names(dflSLAoptimization)[names(dflSLAoptimization) == "cluster"] <- "firstcluster"
names(dflSLAoptimization)[names(dflSLAoptimization) == "SL"] <- "firstSL"
names(SLAoptimizationBIS)[names(SLAoptimizationBIS) == "firstsample"] <- "secondsample"
dflSLAoptimization <- left_join(dflSLAoptimization, SLAoptimizationBIS, by = "secondsample", keep = FALSE)
dflSLAoptimization = dflSLAoptimization %>% drop_na()
names(dflSLAoptimization)[names(dflSLAoptimization) == "cluster"] <- "secondcluster"
names(dflSLAoptimization)[names(dflSLAoptimization) == "SL"] <- "secondSL"
## for single linkage threshold B
names(SLBoptimizationBIS)[names(SLBoptimizationBIS) == "sample"] <- "firstsample"
dflSLBoptimization <- left_join(dfl, SLBoptimizationBIS, by = "firstsample", keep = FALSE)
dflSLBoptimization = dflSLBoptimization %>% drop_na()
names(dflSLBoptimization)[names(dflSLBoptimization) == "cluster"] <- "firstcluster"
names(dflSLBoptimization)[names(dflSLBoptimization) == "SL"] <- "firstSL"
names(SLBoptimizationBIS)[names(SLBoptimizationBIS) == "firstsample"] <- "secondsample"
dflSLBoptimization <- left_join(dflSLBoptimization, SLBoptimizationBIS, by = "secondsample", keep = FALSE)
dflSLBoptimization = dflSLBoptimization %>% drop_na()
names(dflSLBoptimization)[names(dflSLBoptimization) == "cluster"] <- "secondcluster"
names(dflSLBoptimization)[names(dflSLBoptimization) == "SL"] <- "secondSL"
## for single linkage threshold C
names(SLCoptimizationBIS)[names(SLCoptimizationBIS) == "sample"] <- "firstsample"
dflSLCoptimization <- left_join(dfl, SLCoptimizationBIS, by = "firstsample", keep = FALSE)
dflSLCoptimization = dflSLCoptimization %>% drop_na()
names(dflSLCoptimization)[names(dflSLCoptimization) == "cluster"] <- "firstcluster"
names(dflSLCoptimization)[names(dflSLCoptimization) == "SL"] <- "firstSL"
names(SLCoptimizationBIS)[names(SLCoptimizationBIS) == "firstsample"] <- "secondsample"
dflSLCoptimization <- left_join(dflSLCoptimization, SLCoptimizationBIS, by = "secondsample", keep = FALSE)
dflSLCoptimization = dflSLCoptimization %>% drop_na()
names(dflSLCoptimization)[names(dflSLCoptimization) == "cluster"] <- "secondcluster"
names(dflSLCoptimization)[names(dflSLCoptimization) == "SL"] <- "secondSL"
## for single linkage threshold D
names(SLDoptimizationBIS)[names(SLDoptimizationBIS) == "sample"] <- "firstsample"
dflSLDoptimization <- left_join(dfl, SLDoptimizationBIS, by = "firstsample", keep = FALSE)
dflSLDoptimization = dflSLDoptimization %>% drop_na()
names(dflSLDoptimization)[names(dflSLDoptimization) == "cluster"] <- "firstcluster"
names(dflSLDoptimization)[names(dflSLDoptimization) == "SL"] <- "firstSL"
names(SLDoptimizationBIS)[names(SLDoptimizationBIS) == "firstsample"] <- "secondsample"
dflSLDoptimization <- left_join(dflSLDoptimization, SLDoptimizationBIS, by = "secondsample", keep = FALSE)
dflSLDoptimization = dflSLDoptimization %>% drop_na()
names(dflSLDoptimization)[names(dflSLDoptimization) == "cluster"] <- "secondcluster"
names(dflSLDoptimization)[names(dflSLDoptimization) == "SL"] <- "secondSL"

# remove useless variables
dflSLAoptimization$secondSL <- NULL
dflSLBoptimization$secondSL <- NULL
dflSLCoptimization$secondSL <- NULL
dflSLDoptimization$secondSL <- NULL

# prepare dataframes
## for single linkage threshold A
dflSLAoptimizationExtra <- dflSLAoptimization
#### flag the pairwise differences from extra single linkage clusters
dflSLAoptimizationExtra$cluster <- ifelse((dflSLAoptimizationExtra$firstcluster == dflSLAoptimizationExtra$secondcluster), "intra", "extra")
#### subset extra single linkage clusters
SLAPWoptimizationExtra <- subset(dflSLAoptimizationExtra,dflSLAoptimizationExtra$cluster %in% c("extra"))
#### retrieve clusters and smallest extra pairwise differences
SLAPWoptimizationExtramax <- aggregate(x = SLAPWoptimizationExtra$differences, by= list(SLAPWoptimizationExtra$firstcluster), FUN = min)
names(SLAPWoptimizationExtramax)[names(SLAPWoptimizationExtramax) == "Group.1"] <- "SL_cluster"
names(SLAPWoptimizationExtramax)[names(SLAPWoptimizationExtramax) == "x"] <- "PW_min"
#### retrieve clusters presenting the smallest extra pairwise differences
SLAselected <- subset(SLAPWoptimizationExtramax, PW_min == min(PW_min))
#### select randomly one cluster in case of multiple clusters
clusterSLA <- sample (SLAselected$SL_cluster, size=1, replace = F)
### retrieve pairwise differences of interest (i.e. selected cluster versus all clusters)
dflSLAoptimizationClusterInterest <- subset(dflSLAoptimization, firstcluster == clusterSLA)
### flag the pairwise differences from positive (i.e. belonging to the cluster of interest) and negative (i.e. not belonging to the cluster of interest) controls
dflSLAoptimizationClusterInterest$controls <- ifelse((dflSLAoptimizationClusterInterest$firstcluster == dflSLAoptimizationClusterInterest$secondcluster), "positive", "negative")
### rename variables
names(dflSLAoptimizationClusterInterest)[names(dflSLAoptimizationClusterInterest) == "firstcluster"] <- "interestcluster"
names(dflSLAoptimizationClusterInterest)[names(dflSLAoptimizationClusterInterest) == "secondcluster"] <- "othercluster"
names(dflSLAoptimizationClusterInterest)[names(dflSLAoptimizationClusterInterest) == "differences"] <- "PWdifferences"
names(dflSLAoptimizationClusterInterest)[names(dflSLAoptimizationClusterInterest) == "firstSL"] <- "SLthreshold"
## for single linkage threshold B
dflSLBoptimizationExtra <- dflSLBoptimization
dflSLBoptimizationExtra$cluster <- ifelse((dflSLBoptimizationExtra$firstcluster == dflSLBoptimizationExtra$secondcluster), "intra", "extra")
SLBPWoptimizationExtra <- subset(dflSLBoptimizationExtra,dflSLBoptimizationExtra$cluster %in% c("extra"))
SLBPWoptimizationExtramax <- aggregate(x = SLBPWoptimizationExtra$differences, by= list(SLBPWoptimizationExtra$firstcluster), FUN = min)
names(SLBPWoptimizationExtramax)[names(SLBPWoptimizationExtramax) == "Group.1"] <- "SL_cluster"
names(SLBPWoptimizationExtramax)[names(SLBPWoptimizationExtramax) == "x"] <- "PW_min"
SLBselected <- subset(SLBPWoptimizationExtramax, PW_min == min(PW_min))
clusterSLB <- sample (SLBselected$SL_cluster, size=1, replace = F)
dflSLBoptimizationClusterInterest <- subset(dflSLBoptimization, firstcluster == clusterSLB)
dflSLBoptimizationClusterInterest$controls <- ifelse((dflSLBoptimizationClusterInterest$firstcluster == dflSLBoptimizationClusterInterest$secondcluster), "positive", "negative")
names(dflSLBoptimizationClusterInterest)[names(dflSLBoptimizationClusterInterest) == "firstcluster"] <- "interestcluster"
names(dflSLBoptimizationClusterInterest)[names(dflSLBoptimizationClusterInterest) == "secondcluster"] <- "othercluster"
names(dflSLBoptimizationClusterInterest)[names(dflSLBoptimizationClusterInterest) == "differences"] <- "PWdifferences"
names(dflSLBoptimizationClusterInterest)[names(dflSLBoptimizationClusterInterest) == "firstSL"] <- "SLthreshold"
## for single linkage threshold C
dflSLCoptimizationExtra <- dflSLCoptimization
dflSLCoptimizationExtra$cluster <- ifelse((dflSLCoptimizationExtra$firstcluster == dflSLCoptimizationExtra$secondcluster), "intra", "extra")
SLCPWoptimizationExtra <- subset(dflSLCoptimizationExtra,dflSLCoptimizationExtra$cluster %in% c("extra"))
SLCPWoptimizationExtramax <- aggregate(x = SLCPWoptimizationExtra$differences, by= list(SLCPWoptimizationExtra$firstcluster), FUN = min)
names(SLCPWoptimizationExtramax)[names(SLCPWoptimizationExtramax) == "Group.1"] <- "SL_cluster"
names(SLCPWoptimizationExtramax)[names(SLCPWoptimizationExtramax) == "x"] <- "PW_min"
SLCselected <- subset(SLCPWoptimizationExtramax, PW_min == min(PW_min))
clusterSLC <- sample (SLCselected$SL_cluster, size=1, replace = F)
dflSLCoptimizationClusterInterest <- subset(dflSLCoptimization, firstcluster == clusterSLC)
dflSLCoptimizationClusterInterest$controls <- ifelse((dflSLCoptimizationClusterInterest$firstcluster == dflSLCoptimizationClusterInterest$secondcluster), "positive", "negative")
names(dflSLCoptimizationClusterInterest)[names(dflSLCoptimizationClusterInterest) == "firstcluster"] <- "interestcluster"
names(dflSLCoptimizationClusterInterest)[names(dflSLCoptimizationClusterInterest) == "secondcluster"] <- "othercluster"
names(dflSLCoptimizationClusterInterest)[names(dflSLCoptimizationClusterInterest) == "differences"] <- "PWdifferences"
names(dflSLCoptimizationClusterInterest)[names(dflSLCoptimizationClusterInterest) == "firstSL"] <- "SLthreshold"
## for single linkage threshold D
dflSLDoptimizationExtra <- dflSLDoptimization
dflSLDoptimizationExtra$cluster <- ifelse((dflSLDoptimizationExtra$firstcluster == dflSLDoptimizationExtra$secondcluster), "intra", "extra")
SLDPWoptimizationExtra <- subset(dflSLDoptimizationExtra,dflSLDoptimizationExtra$cluster %in% c("extra"))
SLDPWoptimizationExtramax <- aggregate(x = SLDPWoptimizationExtra$differences, by= list(SLDPWoptimizationExtra$firstcluster), FUN = min)
names(SLDPWoptimizationExtramax)[names(SLDPWoptimizationExtramax) == "Group.1"] <- "SL_cluster"
names(SLDPWoptimizationExtramax)[names(SLDPWoptimizationExtramax) == "x"] <- "PW_min"
SLDselected <- subset(SLDPWoptimizationExtramax, PW_min == min(PW_min))
clusterSLD <- sample (SLDselected$SL_cluster, size=1, replace = F)
dflSLDoptimizationClusterInterest <- subset(dflSLDoptimization, firstcluster == clusterSLD)
dflSLDoptimizationClusterInterest$controls <- ifelse((dflSLDoptimizationClusterInterest$firstcluster == dflSLDoptimizationClusterInterest$secondcluster), "positive", "negative")
names(dflSLDoptimizationClusterInterest)[names(dflSLDoptimizationClusterInterest) == "firstcluster"] <- "interestcluster"
names(dflSLDoptimizationClusterInterest)[names(dflSLDoptimizationClusterInterest) == "secondcluster"] <- "othercluster"
names(dflSLDoptimizationClusterInterest)[names(dflSLDoptimizationClusterInterest) == "differences"] <- "PWdifferences"
names(dflSLDoptimizationClusterInterest)[names(dflSLDoptimizationClusterInterest) == "firstSL"] <- "SLthreshold"
## combine and export the dataframe into csv files
dflSLoptimizationClusterInterestAll <- rbind(dflSLAoptimizationClusterInterest, dflSLBoptimizationClusterInterest, dflSLCoptimizationClusterInterest, dflSLDoptimizationClusterInterest)
write.csv(dflSLoptimizationClusterInterestAll,file="Objective4-PairwiseSingleLinkageOptimizationDataframe.csv", row.names=FALSE, quote=FALSE)

# perform Receiver Operating Characteristic (ROC) analysis
## for single linkage threshold A
### run a analysis
SLAROC <- roc(dflSLAoptimizationClusterInterest$controls, dflSLAoptimizationClusterInterest$PWdifferences, percent=TRUE, quiet = TRUE, levels = c("negative", "positive"), direction = ">")
## extract sensitivity and specificity of all thresholds of pairwise differences
SLAoptimizationPWthresholds <- coords(SLAROC, ret=c("threshold", "sensitivity", "specificity"))
## rename variables
colnames(SLAoptimizationPWthresholds) <- c("PWthreshold","sensitivity","specificity")
## export the dataframe into csv files
write.csv(SLAoptimizationPWthresholds,file="Objective4-PairwiseSingleLinkageOptimizationSLAThresholds.csv", row.names=FALSE, quote=FALSE)
## calculate the best pairwise threshold
SLAbestPW <- coords(SLAROC, "best", ret=c("threshold", "sensitivity", "specificity"))
## rename variables
colnames(SLAbestPW) <- c("PWthreshold","sensitivity","specificity")
## for single linkage threshold B
SLBROC <- roc(dflSLBoptimizationClusterInterest$controls, dflSLBoptimizationClusterInterest$PWdifferences, percent=TRUE, quiet = TRUE, levels = c("negative", "positive"), direction = ">")
SLBoptimizationPWthresholds <- coords(SLBROC, ret=c("threshold", "sensitivity", "specificity"))
colnames(SLBoptimizationPWthresholds) <- c("PWthreshold","sensitivity","specificity")
write.csv(SLBoptimizationPWthresholds,file="Objective4-PairwiseSingleLinkageOptimizationSLBThresholds.csv", row.names=FALSE, quote=FALSE)
SLBbestPW <- coords(SLBROC, "best", ret=c("threshold", "sensitivity", "specificity"))
colnames(SLBbestPW) <- c("PWthreshold","sensitivity","specificity")
## for single linkage threshold C
SLCROC <- roc(dflSLCoptimizationClusterInterest$controls, dflSLCoptimizationClusterInterest$PWdifferences, percent=TRUE, quiet = TRUE, levels = c("negative", "positive"), direction = ">")
SLCoptimizationPWthresholds <- coords(SLCROC, ret=c("threshold", "sensitivity", "specificity"))
colnames(SLCoptimizationPWthresholds) <- c("PWthreshold","sensitivity","specificity")
write.csv(SLCoptimizationPWthresholds,file="Objective4-PairwiseSingleLinkageOptimizationSLCThresholds.csv", row.names=FALSE, quote=FALSE)
SLCbestPW <- coords(SLCROC, "best", ret=c("threshold", "sensitivity", "specificity"))
colnames(SLCbestPW) <- c("PWthreshold","sensitivity","specificity")
## for single linkage threshold D
SLDROC <- roc(dflSLDoptimizationClusterInterest$controls, dflSLDoptimizationClusterInterest$PWdifferences, percent=TRUE, quiet = TRUE, levels = c("negative", "positive"), direction = ">")
SLDoptimizationPWthresholds <- coords(SLDROC, ret=c("threshold", "sensitivity", "specificity"))
colnames(SLDoptimizationPWthresholds) <- c("PWthreshold","sensitivity","specificity")
write.csv(SLDoptimizationPWthresholds,file="Objective4-PairwiseSingleLinkageOptimizationSLDThresholds.csv", row.names=FALSE, quote=FALSE)
SLDbestPW <- coords(SLDROC, "best", ret=c("threshold", "sensitivity", "specificity"))
colnames(SLDbestPW) <- c("PWthreshold","sensitivity","specificity")

# export setting, selected samples and best thresholds
## retrieve best thresholds
valueSLAbestPW <- SLAbestPW[1,1]
valueSLBbestPW <- SLBbestPW[1,1]
valueSLCbestPW <- SLCbestPW[1,1]
valueSLDbestPW <- SLDbestPW[1,1]
## export in a txt file
sink("Objective4-PairwiseSingleLinkageOptimizationThresholdsBest.txt")
cat("The SingleLinkagePairwise.R script was executed with the arguments:","\n")
cat("- argSLA: single linkage threshold A =", argSLA,"\n")
cat("- argSLB: single linkage threshold B =", argSLB,"\n")
cat("- argSLC: single linkage threshold C =", argSLC,"\n")
cat("- argSLD: single linkage threshold D =", argSLD,"\n")
cat("- argMNS: minimum number of samples to retain inside each single linkage cluster =", argMNS,"\n")
cat("- argPSP: proportion of splitting (%) defining the optimization dataframe =", argPSP,"\n")
cat("Sensitivity and specificity of pairwise differences (PW) were estimated based on:","\n")
cat("- the single linkage cluster", clusterSLA, ", concerning the single linkage threashold (SL)" ,argSLA, "(i.e. A argument)","\n")
cat("- the single linkage cluster", clusterSLB, ", concerning the single linkage threashold (SL)" ,argSLB, "(i.e. B argument)","\n")
cat("- the single linkage cluster", clusterSLC, ", concerning the single linkage threashold (SL)" ,argSLC, "(i.e. C argument)","\n")
cat("- the single linkage cluster", clusterSLD, ", concerning the single linkage threashold (SL)" ,argSLD, "(i.e. D argument)","\n")
cat("Receiver Operating Characteristic (ROC) analyses allowed estimation of:","\n")
cat("- the best pairwise difference threshold (PW)", valueSLAbestPW, ", concerning the single linkage threashold (SL)" ,argSLA, "(i.e. A argument)","\n")
cat("- the best pairwise difference threshold (PW)", valueSLBbestPW, ", concerning the single linkage threashold (SL)" ,argSLB, "(i.e. B argument)","\n")
cat("- the best pairwise difference threshold (PW)", valueSLCbestPW, ", concerning the single linkage threashold (SL)" ,argSLC, "(i.e. C argument)","\n")
cat("- the best pairwise difference threshold (PW)", valueSLDbestPW, ", concerning the single linkage threashold (SL)" ,argSLD, "(i.e. D argument)","\n")
sink()
