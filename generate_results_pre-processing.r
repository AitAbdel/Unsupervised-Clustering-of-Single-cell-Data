# Université Libre De Bruxelles
# Name : Ait Oujkal
# First Name : Abdellatif
# MASTER 2 - COMPUTER SCIENCE

# import needed packages
library(ggplot2)
library(mclust)
library(splatter)
library(scater)

LINNORM_WITH_PRE = "Results_simulated_data/Res_Category_3/sim_c5000_g1000_gr16_balanced_de0.5_drop0/Linnorm_method_specific_PCA_TRUE_3_fixed_estimate_NA.Rdata"
LINNORM_WITHOUT_PRE = "Results_simulated_data/Res_Category_3/sim_c5000_g1000_gr16_balanced_de0.5_drop0/Linnorm_none_PCA_TRUE_3_fixed_estimate_NA.Rdata"

SC3_WITH_PRE = "Results_simulated_data/Res_Category_3/sim_c5000_g1000_gr16_balanced_de0.5_drop0/SC3_method_specific_internal_internal_NA_fixed_estimate_NA.Rdata"
SC3_WITHOUT_PRE = "Results_simulated_data/Res_Category_3/sim_c5000_g1000_gr16_balanced_de0.5_drop0/SC3_none_internal_internal_NA_fixed_estimate_NA.Rdata"

x <- c()
load(LINNORM_WITHOUT_PRE)
x <- c(x,adjustedRandIndex(true_labels, Result))
load(SC3_WITHOUT_PRE)
x <- c(x,adjustedRandIndex(true_labels, Result))

y <- c()
load(LINNORM_WITH_PRE)
y <- c(y,adjustedRandIndex(true_labels, Result))
load(SC3_WITH_PRE)
y <- c(y,adjustedRandIndex(true_labels, Result))

value_matrix = matrix(, nrow = 2, ncol = 2)

value_matrix[1,] = x
value_matrix[2,] = y

png(paste0("results_figures/pre-processing.png"))

par(mar = c(5,4,4,10))

barplot(main = "Pre-processing vs No Pre-processing",
        ylab = "ARI",
        xlab= "Clustering methods",
        value_matrix, 
        names.arg = c("Linnorm","SC3"), 
        beside = TRUE, 
        col = c("peachpuff", "skyblue"),
        ylim=c(0,1))
        
legend(x = "topright", legend = c("without pre-process", "with pre-process"), fill = c("peachpuff", "skyblue"),inset=c(-0.5,0), xpd = TRUE)

dev.off()