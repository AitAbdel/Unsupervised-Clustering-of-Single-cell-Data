# Université Libre De Bruxelles
# Name : Ait Oujkal
# First Name : Abdellatif
# MASTER 2 - COMPUTER SCIENCE

# import needed packages
library(scater)
library(ggplot2)

#system("mkdir datasets_figures")

# ALL DATASETS
DATASET1_CATEGORY1 = "simulated_data/Category_1/sim_c1000_g1000_gr4_balanced_de0.1_drop0.Rdata"
DATASET2_CATEGORY1 = "simulated_data/Category_1/sim_c1000_g1000_gr4_balanced_de0.5_drop0.Rdata"
DATASET3_CATEGORY1 = "simulated_data/Category_1/sim_c1000_g1000_gr4_balanced_de0.9_drop0.Rdata"

DATASET1_CATEGORY2 = "simulated_data/Category_2/sim_c1000_g1000_gr4_balanced_de0.5_drop0.Rdata"
DATASET2_CATEGORY2 = "simulated_data/Category_2/sim_c1000_g1000_gr4_balanced_de0.5_drop2.Rdata"
DATASET3_CATEGORY2 = "simulated_data/Category_2/sim_c1000_g1000_gr4_balanced_de0.5_drop4.Rdata"
DATASET4_CATEGORY2 = "simulated_data/Category_2/sim_c1000_g1000_gr4_balanced_de0.5_drop6.Rdata"

DATASET1_CATEGORY3 = "simulated_data/Category_3/sim_c500_g1000_gr4_balanced_de0.5_drop0.Rdata"
DATASET2_CATEGORY3 = "simulated_data/Category_3/sim_c1000_g1000_gr8_unbalanced_de0.5_drop0.Rdata"
DATASET3_CATEGORY3 = "simulated_data/Category_3/sim_c5000_g1000_gr16_unbalanced_de0.5_drop0.Rdata"

#################################################################################################
######################################### CATEGORY 1 ############################################
#################################################################################################

load(DATASET1_CATEGORY1)
sce <- logNormCounts(sce)
sce <- runPCA(sce)
plotPCA(sce, colour_by = "Group")
ggsave("datasets_figures/DATASET1_CATEGORY1.png")

load(DATASET2_CATEGORY1)
sce <- logNormCounts(sce)
sce <- runPCA(sce)
plotPCA(sce, colour_by = "Group")
ggsave("datasets_figures/DATASET2_CATEGORY1.png")

load(DATASET3_CATEGORY1)
sce <- logNormCounts(sce)
sce <- runPCA(sce)
plotPCA(sce, colour_by = "Group")
ggsave("datasets_figures/DATASET3_CATEGORY1.png")

#################################################################################################
######################################### CATEGORY 2 ############################################
#################################################################################################

load(DATASET1_CATEGORY2)
sce <- logNormCounts(sce)
sce <- runPCA(sce)
plotPCA(sce, colour_by = "Group")
ggsave("datasets_figures/DATASET1_CATEGORY2.png")

load(DATASET2_CATEGORY2)
sce <- logNormCounts(sce)
sce <- runPCA(sce)
plotPCA(sce, colour_by = "Group")
ggsave("datasets_figures/DATASET2_CATEGORY2.png")

load(DATASET3_CATEGORY2)
sce <- logNormCounts(sce)
sce <- runPCA(sce)
plotPCA(sce, colour_by = "Group")
ggsave("datasets_figures/DATASET3_CATEGORY2.png")

load(DATASET4_CATEGORY2)
sce <- logNormCounts(sce)
sce <- runPCA(sce)
plotPCA(sce, colour_by = "Group")
ggsave("datasets_figures/DATASET4_CATEGORY2.png")


#################################################################################################
######################################### CATEGORY 3 ############################################
#################################################################################################

load(DATASET1_CATEGORY3)
sce <- runTSNE(sce, exprs_values = "logcounts", perplexity = 10)
plotTSNE(sce, colour_by = "Group")
ggsave("datasets_figures/DATASET1_CATEGORY3.png")

load(DATASET2_CATEGORY3)
sce <- runTSNE(sce, exprs_values = "logcounts", perplexity = 10)
plotTSNE(sce, colour_by = "Group")
ggsave("datasets_figures/DATASET2_CATEGORY3.png")

load(DATASET3_CATEGORY3)
sce <- runTSNE(sce, exprs_values = "logcounts", perplexity = 10)
plotTSNE(sce, colour_by = "Group")
ggsave("datasets_figures/DATASET3_CATEGORY3.png")