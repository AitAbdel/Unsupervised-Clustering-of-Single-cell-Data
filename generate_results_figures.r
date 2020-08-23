# Université Libre De Bruxelles
# Name : Ait Oujkal
# First Name : Abdellatif
# MASTER 2 - COMPUTER SCIENCE

# import needed packages
library(scater)
library(mclust)
library(ggplot2)
library(RColorBrewer)

#system("mkdir results_figures")       # creates a new directory where all the plots will be stored

# Plots the ARI performance of the different methods for a given dataset
plot_bar_ARI <- function(dataset,title,ylabel,methods,ARI,coul){

    png(paste0("results_figures/",dataset,".png"))

    par(mar = c(5,4,4,7))

    barplot(ARI,
    main = title,
    ylab = ylabel,
    ylim = c(0,1),
    names.arg = methods,
    col = coul,
    horiz = FALSE,
    las=2)

    legend(x = "topright", legend = unique(methods), fill = unique(coul),inset=c(-0.3,0), xpd = TRUE)

    dev.off()
}

# Plots the Computational time of the different methods for a given dataset
plot_bar_time <- function(dataset,title,ylabel,methods,ARI,coul){

    png(paste0("results_figures/",dataset,".png"))

    par(mar = c(5,4,4,7))

    barplot(ARI,
    main = title,
    ylab = ylabel,
    names.arg = methods,
    col = coul,
    horiz = FALSE,
    las=2)

    legend(x = "topright", legend = unique(methods), fill = unique(coul),inset=c(-0.3,0), xpd = TRUE)

    dev.off()
}

# Computes the ARI of each method for a given dataset
compute_ARI <- function(dataset,category,methods){

    ARI <- list()
    for(i in 1:length(methods)){

        f = paste0("Results_simulated_data/Res_Category_",category,"/",dataset,"/",methods[i])

        if (file.exists(f)) {
            data=load(f)

            ARI[i]=adjustedRandIndex(true_labels, Result)
        }else
        {
            ARI[i]=0
        }
    }
    return(ARI)
}

# Computes the time it took to cluster the given dataset for each method
compute_TIME <- function(dataset,category,methods){

    TIME <- list()
    for(i in 1:length(methods)){
        data=load(paste0("Results_simulated_data/Res_Category_",category,"/",dataset,"/",methods[i]))
        t = Run_time/60
        TIME[i]=log(as.double(t)+1)
    }
    return(TIME)
}

# ALL SIMULATED DATASETS
DATASET1_CATEGORY3 = "sim_c500_g1000_gr16_balanced_de0.5_drop0"
DATASET2_CATEGORY3 = "sim_c1000_g1000_gr16_balanced_de0.5_drop0"
DATASET3_CATEGORY3 = "sim_c5000_g1000_gr16_balanced_de0.5_drop0"
DATASET4_CATEGORY3 = "sim_c1000_g1000_gr8_balanced_de0.5_drop0"
DATASET5_CATEGORY3 = "sim_c1000_g1000_gr8_unbalanced_de0.5_drop0"

DATASET1_CATEGORY1 = "sim_c1000_g1000_gr4_balanced_de0.1_drop0"
DATASET2_CATEGORY1 = "sim_c1000_g1000_gr4_balanced_de0.5_drop0"
DATASET3_CATEGORY1 = "sim_c1000_g1000_gr4_balanced_de0.9_drop0"

DATASET1_CATEGORY2 = "sim_c1000_g1000_gr4_balanced_de0.5_drop0"
DATASET2_CATEGORY2 = "sim_c1000_g1000_gr4_balanced_de0.5_drop2"
DATASET3_CATEGORY2 = "sim_c1000_g1000_gr4_balanced_de0.5_drop4"
DATASET4_CATEGORY2 = "sim_c1000_g1000_gr4_balanced_de0.5_drop6"

# CREATE ARRAY OF COLORS FOR THE PLOT FUNCTION
coul <- brewer.pal(5, "Set2")

colors = c(coul[1],coul[1],coul[1],coul[1],
            coul[2],coul[2],coul[2],coul[2],coul[2],coul[2],coul[2],coul[2],coul[2],coul[2],coul[2],coul[2],coul[2],coul[2],coul[2],coul[2],
            coul[3],coul[3],coul[3],coul[3],coul[3],coul[3],coul[3],coul[3],
            coul[4],coul[4],coul[4],coul[4],
            coul[5],coul[5])

######################################################################################################
##################################### METHODS PERFORMANCE ############################################
######################################################################################################

methods = list.files(paste0("Results_simulated_data/Res_Category_2/",DATASET1_CATEGORY2))

################# CATEGORY 1 #################

ARI = compute_ARI(DATASET1_CATEGORY1, "1", methods)

plot_bar_ARI(DATASET1_CATEGORY1,"de probability = 0.1","ARI",sapply(strsplit(methods,"_"), `[`, 1), unlist(ARI),colors)

ARI = compute_ARI(DATASET2_CATEGORY1, "1", methods)

plot_bar_ARI(paste0(DATASET3_CATEGORY1,"_"),"de probability = 0.5","ARI",sapply(strsplit(methods,"_"), `[`, 1), unlist(ARI),colors)

ARI = compute_ARI(DATASET3_CATEGORY1, "1", methods)

plot_bar_ARI(DATASET3_CATEGORY1,"de probability = 0.9","ARI",sapply(strsplit(methods,"_"), `[`, 1), unlist(ARI),colors)

################# CATEGORY 3 #################

ARI = compute_ARI(DATASET1_CATEGORY3, "3", methods)

plot_bar_ARI(DATASET1_CATEGORY3,"Number of cells = 500","ARI",sapply(strsplit(methods,"_"), `[`, 1), unlist(ARI),colors)

ARI = compute_ARI(DATASET2_CATEGORY3, "3", methods)

plot_bar_ARI(DATASET2_CATEGORY3,"Number of cells = 1000","ARI",sapply(strsplit(methods,"_"), `[`, 1), unlist(ARI),colors)

ARI = compute_ARI(DATASET3_CATEGORY3, "3", methods)

plot_bar_ARI(DATASET3_CATEGORY3,"Number of cells = 5000","ARI",sapply(strsplit(methods,"_"), `[`, 1), unlist(ARI),colors)

ARI = compute_ARI(DATASET4_CATEGORY3, "3", methods)

plot_bar_ARI(DATASET4_CATEGORY3,"balanced","ARI",sapply(strsplit(methods,"_"), `[`, 1), unlist(ARI),colors)

ARI = compute_ARI(DATASET5_CATEGORY3, "3", methods)

plot_bar_ARI(DATASET5_CATEGORY3,"unbalanced","ARI",sapply(strsplit(methods,"_"), `[`, 1), unlist(ARI),colors)

################# CATEGORY 2 #################

ARI = compute_ARI(DATASET1_CATEGORY2, "2", methods)

plot_bar_ARI(DATASET1_CATEGORY2,"dropout = 24%","ARI",sapply(strsplit(methods,"_"), `[`, 1), unlist(ARI),colors)

ARI = compute_ARI(DATASET2_CATEGORY2, "2", methods)

plot_bar_ARI(DATASET2_CATEGORY2,"dropout = 45%","ARI",sapply(strsplit(methods,"_"), `[`, 1), unlist(ARI),colors)

ARI = compute_ARI(DATASET3_CATEGORY2, "2", methods)

plot_bar_ARI(DATASET3_CATEGORY2,"dropout = 74%","ARI",sapply(strsplit(methods,"_"), `[`, 1), unlist(ARI),colors)

ARI = compute_ARI(DATASET4_CATEGORY2, "2", methods)

plot_bar_ARI(DATASET4_CATEGORY2,"dropout = 91%","ARI",sapply(strsplit(methods,"_"), `[`, 1), unlist(ARI),colors)

###################################################################################################### 
##################################### COMPUTATIONAL TIME ############################################# 
###################################################################################################### 

# DATASETS
DATASET1 = "sim_c500_g1000_gr16_balanced_de0.5_drop0"
DATASET2 = "sim_c1000_g1000_gr16_balanced_de0.5_drop0"
DATASET3 = "sim_c5000_g1000_gr16_balanced_de0.5_drop0"

methods = list.files(paste0("Results_simulated_data/Res_Category_3/",DATASET1))

TIME = compute_TIME(DATASET1, "3", methods)

plot_bar_time(paste0(DATASET1,"_time"),"Number of cells = 500","TIME log (run time in minutes)",sapply(strsplit(methods,"_"), `[`, 1), unlist(TIME),colors)

methods = list.files(paste0("Results_simulated_data/Res_Category_3/",DATASET2))

TIME = compute_TIME(DATASET2, "3", methods)

plot_bar_time(paste0(DATASET2,"_time"),"Number of cells = 1000","TIME log (run time in minutes)",sapply(strsplit(methods,"_"), `[`, 1), unlist(TIME),colors)

methods = list.files(paste0("Results_simulated_data/Res_Category_3/",DATASET3))

TIME = compute_TIME(DATASET3, "3", methods)

plot_bar_time(paste0(DATASET3,"_time"),"Number of cells = 5000","TIME log (run time in minutes)",sapply(strsplit(methods,"_"), `[`, 1), unlist(TIME),colors)