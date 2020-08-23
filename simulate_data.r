# Université Libre De Bruxelles
# Name : Ait Oujkal
# First Name : Abdellatif
# MASTER 2 - COMPUTER SCIENCE

# import needed packages
library(splatter)

# Pre-process the given dataset and stores the pre-processed dataset
pre_process <- function(dataset,output_file){
    
    load(dataset)
    
    # **** First preprocessing: gene filtering ****
    # Filter genes with average expression count (adjusted by library size) equal to 0
    rowData(sce)$ave.counts <- calculateAverage(sce, exprs_values = "counts", use_size_factors=FALSE)
    to.keep <- rowData(sce)$ave.counts > 0
    sce <- sce[to.keep,]

    # **** Second preprocessing: normalization after filtering ****
    # Normalize dataset
    set.seed(100)
    if(dim(sce)[2]<500){
     min.size=50
    } else {
     min.size=200   
    }

    clusters <- quickCluster(sce, min.size=min.size, min.mean=0.1, method="igraph")
    sce <- computeSumFactors(sce, cluster=clusters, min.mean=0.1)
    sce <- scater::normalize(sce, return_log = FALSE)
    counts(sce) <- normcounts(sce) #assign normalized counts to slot counts

    # Save normalized, filtered and qc dataset 
    save(sce, file=output_file)
}

# Simulates all the datasets 
simulateSplatter <- function(nCells, nGenes, group.prob, type, de.prob, dropout.type, dropout.mid, category_n){
    
    #system(paste0("mkdir simulated_data/Category_", category_n))
    #system(paste0("mkdir pre-processed_data/Category_", category_n))
    
    for(i in 1:length(nCells)){
        for(j in 1:length(nGenes)){
            for(k in 1:length(group.prob)){
                for(l in 1:length(de.prob)){
                    for(m in 1:length(dropout.mid)){
                        
                        params <- newSplatParams()
                        params <- setParams(params, batchCells=nCells[i], nGenes=nGenes[j], group.prob = unlist(group.prob[k]), de.prob=de.prob[l], dropout.type=dropout.type)
                        sce <- splatSimulate(params, method = "groups", verbose = FALSE)
                        rowData(sce)$feature_symbol <- rownames(sce)
                        logcounts(sce) <- log2(counts(sce)+1)
                        sce <- splatter:::splatSimDropout(sce, setParam(params, "dropout.mid", dropout.mid[m]))
                        
                        # save simulated data
                        fname=paste0("sim_c", nCells[i] , "_g", nGenes[j], "_gr", length(unlist(group.prob[k])), "_", type[k], "_de", de.prob[l], "_drop", dropout.mid[m])
                        save(sce, file=paste0("simulated_data/Category_", category_n, "/", fname, ".Rdata"))
                        
                        # save pre-processed data
                        #pre_process(paste0("simulated_data/Category_", category_n, "/", fname, ".Rdata"), paste0("pre-processed_data/Category_", category_n, "/", fname, ".Rdata"))
                    }
                }
            }
        }
    }
}

#####################################################################################################
######################################### DATASETS GENERATION #######################################
#####################################################################################################

#system(paste0("mkdir simulated_data/"))
#system(paste0("mkdir pre-processed_data/"))  

######################################### DATASETS CATEGORY 3 #######################################
nCells=c(500,1000,5000)
nGenes=1000         # fixed
group.prob=list(rep(1/4,4), 
                rep(1/8,8), 
                rep(1/16,16),
                c(0.05, 0.05, 0.4, 0.5), 
                c(0.01, 0.01, 0.04 , 0.04, 0.2, 0.2, 0.2, 0.3),
                c(0.01, 0.01, 0.04 , 0.04, 0.1, 0.1, 0.1, 0.1, 0.01, 0.01, 0.04 , 0.04, 0.1, 0.1, 0.1, 0.1))

type=c("balanced", 
        "balanced", 
        "balanced",
        "unbalanced",
        "unbalanced",
        "unbalanced")
de.prob=0.5
dropout.type ="experiment"
dropout.mid=0
category_n = "3"

# simulation
simulateSplatter(nCells, nGenes, group.prob, type, de.prob, dropout.type, dropout.mid, category_n)

######################################### DATASETS CATEGORY 1 #######################################
nCells=1000
nGenes=1000
group.prob=list(c(0.25, 0.25, 0.25, 0.25))
type="balanced"
de.prob=c(0.1,0.5,0.9)
dropout.type ="experiment"
dropout.mid=0
category_n = "1"

# simulation
simulateSplatter(nCells, nGenes, group.prob, type, de.prob, dropout.type, dropout.mid, category_n)

######################################### DATASETS CATEGORY 2 #######################################
nCells=1000
nGenes=1000
group.prob=list(c(0.25, 0.25, 0.25, 0.25))
type="balanced"
de.prob=0.5
dropout.type = "experiment"
dropout.mid=c(0,2,4,6)
category_n = "2"

# simulation
simulateSplatter(nCells, nGenes, group.prob, type, de.prob, dropout.type, dropout.mid, category_n)