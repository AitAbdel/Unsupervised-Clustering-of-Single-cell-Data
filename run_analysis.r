# Université Libre De Bruxelles
# Name : Ait Oujkal
# First Name : Abdellatif
# MASTER 2 - COMPUTER SCIENCE

# import needed packages
library(scater)
library(tidyr)

source("run_methods/run_methods.R")

output_directory = "Results_simulated_data"

log_file <-file(paste0(output_directory, "/log.txt"))
sink(log_file)

P=read.csv(file="Valid_parameter_settings.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
parameter_combinations=P
for(i in 1:dim(P)[2]){
    parameter_combinations=separate_rows(parameter_combinations, names(parameter_combinations)[i], sep = "/", convert = TRUE)
}

parameter_settings = parameter_combinations
dimensions_input = 3

category_1 = list.files("simulated_data/Category_1", full.names = TRUE)
df=data.frame(t(sapply(strsplit(category_1, "_"),c)))
df=df[order(as.numeric(sub("c", "", df$X4)), decreasing=FALSE),]
category_1 <- paste(df$X1, df$X2, df$X3, df$X4, df$X5, df$X6, df$X7, df$X8, df$X9, sep='_')
#datasets=list("Category_2" = list.files("simulated_data/Category_2", full.names = TRUE), "Category_3" = list.files("simulated_data/Category_3", full.names = TRUE), "Category_1" = category_1)
datasets=list("Category_3" = list.files("simulated_data/Category_3", full.names = TRUE))

for(h in 1:length(datasets)){
   #system(paste0("mkdir ", output_directory, "/Res_", names(datasets)[h]))
    
    for(i in 1:length(datasets[[h]])){
        
        data=load(datasets[[h]][i])
        data.name=sub('\\.Rdata$', '', unlist(strsplit(datasets[[h]][i], "/"))[length(unlist(strsplit(datasets[[h]][i], "/")))])

        #system(paste0("mkdir ", output_directory, "/Res_", names(datasets)[h], "/", data.name))

        for(j in 1:dim(parameter_settings)[1]){

            preproc=parameter_settings[j,]$preproc
            method=parameter_settings[j,]$method
            dim_red=parameter_settings[j,]$dim_red
            use_dims_input=parameter_settings[j,]$use_dims_input
            clust_tech=parameter_settings[j,]$clust_tech
            clust=parameter_settings[j,]$clust

            if(clust=="set"){
                k_input=length(unique(sce$Group))
            }else if(clust=="estimate"){
                k_input=NA
            }

            if(use_dims_input=="TRUE"){
                dims_input=dimensions_input
            }else if(use_dims_input=="internal" | use_dims_input=="FALSE"){
                dims_input=NA
            }

            file_name=paste0(method,"_", preproc,"_",dim_red,"_",use_dims_input,"_", dims_input, "_", clust_tech,"_", clust, "_", k_input)
            print(paste0("Running..." , file_name, "on ", datasets[[h]][i]))

            #Running benchmark
            start_time <- Sys.time()
            Result0=tryCatch(
                run_methods(data=counts(sce), method=method, preproc=preproc, dim_red=dim_red, use_dims_input=use_dims_input, dims_input=dims_input, clust_tech=clust_tech, clust=clust, k_input=k_input),
                error = function(e) { write(as.character(e), file=log_file, append = TRUE, sep = "\n") }
                )

            end_time <- Sys.time()
            Run_time=difftime(end_time, start_time, units = "secs")

            #Saving results
            if(is.integer(Result0)==TRUE | is.numeric(Result0)==TRUE){
                Result=Result0
                true_labels=sce$Group
            }else if(is.null(Result0)){
                Result=Result0
                true_labels=sce$Group
            }else{
                Result=Result0$clusters
                true_labels=sce$Group[-Result0$missing_ids]}

            file_path=paste0("Res_",names(datasets)[h], "/", data.name ,"/", file_name , ".Rdata")

            save(Result, Run_time, true_labels, file=paste0(output_directory,"/",file_path))
        }
    }
}
sink()