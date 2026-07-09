library(RANN)
#initial_df and query is the reference and query coordiante dataframe, include(x, y) or (PC1, PC2, PC3...)
#sm_vector is the reference vector which will map to query
#round is the average round
#knn is the n nearest point from reference
smooth_kNN <- function(initial_df,query_df,sm_vector,round=100,knn=30){
    result <- RANN::nn2(initial_df,query_df,k=knn)
    if(!identical(as.matrix(initial_df), as.matrix(query_df))){
        sm_vector <- rowMeans(matrix(sm_vector[result$nn.idx],nrow = nrow(result$nn.idx), byrow = FALSE))
        return(sm_vector)
        stop()
    }
    for(i in 1:round){
        sm_vector <- rowMeans(matrix(sm_vector[result$nn.idx],nrow = nrow(result$nn.idx), byrow = FALSE))
    }
    return(sm_vector)
}

smooth_kNN2 <- function(initial_df,query_df,sm_vector,round=100,knn=30){
    result <- RANN::nn2(initial_df,query_df,k=knn)
    query_Nx <- nrow(query_df)
    initial_Nx <- nrow(initial_df)
    knnIdx <- result$nn.idx
    W <- Matrix::sparseMatrix(rep(seq_len(query_Nx), knn), c(knnIdx), x=1, dims = c(query_Nx, initial_Nx))
    W <- W/rowSums(W)
    if(!identical(as.matrix(initial_df), as.matrix(query_df))){
        sm_vector <- as.matrix(W %*% as.matrix(sm_vector))
        return(sm_vector)
        stop()
    }
    for(i in 1:round){
        sm_vector <- as.matrix(W %*% as.matrix(sm_vector))
    }
    return(sm_vector)
}


#initial_df and query is the reference and query coordiante dataframe, include(x, y) or (PC1, PC2, PC3...)
#sm_vector is the reference vector which will map to query
#knn is the n nearest point from reference
#If the value counts exceed the threshold, specific actions apply; otherwise, use the nearest point as the query value.
winner_kNN <- function(initial_df,query_df,sm_vector,knn=5,threhold=1){
    result <- RANN::nn2(initial_df,query_df,k=knn)
    sm_vector <- matrix(sm_vector[result$nn.idx],nrow = nrow(result$nn.idx), byrow = FALSE)
    sm_vector <- sapply(1:nrow(sm_vector),function(x){
        tmp_table <- table(sm_vector[x,])
        name <- names(sort(tmp_table,decreasing = T))[1]
        if(tmp_table[name]>threhold){return(name)}else{return(sm_vector[x,1])}
    })
    return(sm_vector)
}
magic_knn <- function(smooth_df,sm_vector,round=3,knn=15,ka=3,epsilon=1){
    Nx <- nrow(smooth_df)
    result <- RANN::nn2(smooth_df,smooth_df,k=knn)
    knnIdx <- result$nn.idx
    knnDist <- result$nn.dists
    
    knnDist <- knnDist / knnDist[,ka]
    #Markov compute
    W <- Matrix::sparseMatrix(rep(seq_len(Nx), knn), c(knnIdx), x=c(knnDist), dims = c(Nx, Nx))
    W <- W + Matrix::t(W)
    #Compute Kernel
    W@x <- exp(-(W@x / epsilon^2))
    #Markov normalization
    W <- W / Matrix::rowSums(W)
    #Initialize Matrix
    Wt <- W
    #Computing Diffusion Matrix
    for(i in seq_len(round)){
      Wt <- Wt %*% W
    }
    sm_vector_sm <- as.matrix(Wt %*% sm_vector)
    return(sm_vector_sm)
}
ROC_kNN <- function (initial_df, vector,normalize=FALSE, knn = 9) 
{
    result <- RANN::nn2(initial_df, initial_df, k = knn)
    vector_n <- matrix(vector[result$nn.idx], nrow = nrow(result$nn.idx), 
        byrow = FALSE)
    x_n <- matrix(initial_df[, 1][result$nn.idx], nrow = nrow(result$nn.idx), 
        byrow = FALSE)
    y_n <- matrix(initial_df[, 2][result$nn.idx], nrow = nrow(result$nn.idx), 
        byrow = FALSE)
    vector_max <- matrixStats::rowMaxs(vector_n)
    vector_min <- matrixStats::rowMins(vector_n)
    gradient_top <- abs(vector_max - vector_min)
    gradient_bottom <- sapply(1:nrow(vector_n), function(i) {
        x_max <- x_n[i, ][vector_n[i, ] == vector_max[i]][1]
        y_max <- y_n[i, ][vector_n[i, ] == vector_max[i]][1]
        x_min <- x_n[i, ][vector_n[i, ] == vector_min[i]][1]
        y_min <- y_n[i, ][vector_n[i, ] == vector_min[i]][1]
        distance <- ((x_max - x_min)^2 + (y_max - y_min)^2)^0.5
        return(distance)
    })
    result <- gradient_top/gradient_bottom
    result[is.na(result)] <- 0
    result[is.infinite(result)] <- 0
    if(normalize){
        result <- (result-min(result))/(max(result)-min(result))
        result[is.na(result)] <- 0}
    return(result)
}

smooth_obj_data <- function(obj,knn=16,round=2,smooth_function='spatial',assays='RNA',slot='data',reduction='pca',run_PCA=TRUE){
    
    if(smooth_function=='spatial'){
        data_matrix <- as.matrix(slot(obj@assays[[assays]],slot))
        pb <- txtProgressBar(style=3)
        n <- 1
        all_n <- nrow(data_matrix)
        for(i in rownames(data_matrix)){
            data_matrix[i,] <- smooth_kNN(obj@meta.data[,c('x','y')],obj@meta.data[,c('x','y')],data_matrix[i,],knn=knn,round=round)
            setTxtProgressBar(pb, n/all_n)
            n <- n+1}
        close(pb)
        
    
    }else if(smooth_function=='spatial2'){
        data_matrix <- Matrix::t((slot(obj@assays[[assays]],slot)))
        data_matrix <- smooth_kNN2(obj@meta.data[,c('x','y')],obj@meta.data[,c('x','y')],data_matrix,knn=knn,round=round)
        data_matrix <- t(data_matrix)
        colnames(data_matrix) <- colnames(obj)
    
    
    }else if(smooth_function=='magic'){
        if(run_PCA){
            obj <- FindVariableFeatures(obj,nfeatures = 2000)
            obj <- ScaleData(obj)
            obj <- RunPCA(obj,verbose=FALSE)}

        data_matrix <- magic_knn(obj[[reduction]]@cell.embeddings,
                         Matrix::t(slot(obj@assays[[assays]],slot)),
                         knn=knn,
                         round=round)
        data_matrix <- t(data_matrix)
        colnames(data_matrix) <- colnames(obj)
    
    }
    data_matrix <- CreateAssayObject(data=data_matrix)
    slot(obj@assays[[assays]],slot) <- data_matrix@data
    return(obj)
}