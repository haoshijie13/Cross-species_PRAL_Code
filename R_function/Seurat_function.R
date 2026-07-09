save_seurat_key <- function(obj,path){
    save_key_list <- list('reductions'=obj@reductions,'meta.data'=obj@meta.data,'graphs'=obj@graphs,'meta.features'=obj@assays$RNA@meta.features)
    saveRDS(save_key_list,path)
}
read_seurat_key <- function(obj,path){
    save_key_list <- readRDS(path)
    obj@reductions <- save_key_list$reduction
    obj@meta.data <- save_key_list$meta.data
    obj@graphs <- save_key_list$graphs
    obj@assays$RNA@meta.features <- save_key_list$meta.features
    return(obj)
}

