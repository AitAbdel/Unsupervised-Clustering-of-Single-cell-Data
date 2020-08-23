run_Seurat <- function(data, preproc, dim_red, use_dims_input ,dims_input, clust_tech, clust, k_input){
    library(Seurat)
  

    cso <- Seurat::CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)  #filter
    cso <- subset(cso, subset = nFeature_RNA > 40)
    cso <- Seurat::NormalizeData(cso, normalization.method = "LogNormalize", scale.factor = 10000) #normalize
    cso <- FindVariableFeatures(cso, selection.method = "vst", nfeatures = 1000)

    # Data clustering
    all.genes <- rownames(cso)
    cso <- Seurat::ScaleData(cso, features = all.genes) #scaling data is a must
    
    cso <- Seurat::RunPCA(cso,  features = VariableFeatures(object = cso))
    cso <- Seurat::FindNeighbors(cso, dims = 1:10)
    cso <- Seurat::FindClusters(cso, resolution = 0.1)

    clusters <- as.numeric(Idents(cso))
    detach(package:Seurat)
    return(clusters)
}
