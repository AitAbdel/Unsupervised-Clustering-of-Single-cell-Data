# Université Libre De Bruxelles
# Name : Ait Oujkal
# First Name : Abdellatif
# MASTER 2 - COMPUTER SCIENCE

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#####################################################
# TO INSTALL PACKAGES FROM GITHUB
#####################################################

install.packages("devtools")
library(devtools)

#####################################################
# CRAN PACKAGES
#####################################################

install.packages(c("reshape2", "fields", "ggbeeswarm", "gridExtra",
                   "dynamicTreeCut", "dendextend", "RColorBrewer",
                   "locfit", "KernSmooth","mnormt","mclust"))


#####################################################
# CUSTERING METHODS AND SIMULATION PACKAGES
#####################################################

BiocManager::install(c("Biobase", "BiocGenerics", "BiocParallel","ComplexHeatmap","scater","splatter","scran",
						"SingleCellExperiment", "GenomeInfoDb", "GenomeInfoDbData","IMB-Computational-Genomics-Lab/ascend",
						"VCCRI/CIDR","Linnorm","monocle","pcaMethods","JustinaZ/pcaReduce","RaceID","SC3","SIMLR","Seurat",
						"sincell","TSCAN"))

#####################################################
# DATA MANIPULATION & PLOTTING
#####################################################

install.packages("ggplot2")
install.packages("ggfortify")
install.packages("GGally")
install.packages("dplyr")
install.packages("evaluate")
install.packages("mclust")
install.packages("tidyr")
install.packages("WriteXLS")