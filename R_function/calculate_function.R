#Calculate function
matrix_groupby_raw <- function(matrix_input,group.by,type='row',cal='sum'){
    if(type=='col'){matrix_input <- t(matrix_input)}
    combined_matrix <- data.frame('group'=group.by,'raw'=rownames(matrix_input),'identify'=1)
    combined_matrix <- reshape2::acast(combined_matrix,group~raw,fill = 0,value.var = 'identify')
    combined_matrix <- combined_matrix[,rownames(matrix_input)]
    if(cal=='sum'){matrix_output <- combined_matrix%*%matrix_input}else if(cal=='mean'){
        combined_matrix <- combined_matrix/rowSums(combined_matrix)
        matrix_output <- combined_matrix%*%matrix_input
    }
    if(type=='col'){matrix_output <- t(matrix_output)}
    return(matrix_output)}



matrix_groupby <- function(matrix_input, group.by, type = 'row', cal = 'sum') {
    if (type == 'col') {
    matrix_input <- Matrix::t(matrix_input)
    }
    
    f_group <- factor(group.by)
    row_idx <- as.numeric(f_group)
    group_levels <- levels(f_group)
    rownames(matrix_input) <- as.character(1:nrow(matrix_input))
    col_idx <- as.numeric(factor(rownames(matrix_input)))
    
    combined_matrix <- Matrix::sparseMatrix(
    i = row_idx, j = col_idx, x = 1,
    dimnames = list(group_levels, levels(factor(rownames(matrix_input))))
    )
    combined_matrix <- combined_matrix[, rownames(matrix_input),drop=FALSE]
    
    if (cal == 'sum') {
        matrix_output <- combined_matrix %*% matrix_input
    } else if (cal == 'mean') {
        combined_matrix <- combined_matrix / Matrix::rowSums(combined_matrix)
        matrix_output <- combined_matrix %*% matrix_input
    } else if(cal=='max'){
        group_rows <- split(seq_len(nrow(matrix_input)), row_idx)
        matrix_output <- sapply(group_rows, function(rows) apply(matrix_input[rows, , drop=FALSE], 2, max))
        matrix_output <- t(matrix_output)
    } else if(cal=='min'){
        group_rows <- split(seq_len(nrow(matrix_input)), row_idx)
        matrix_output <- sapply(group_rows, function(rows) apply(matrix_input[rows, , drop=FALSE], 2, min))
        matrix_output <- t(matrix_output)
    } else {
        stop(paste("cal just sum/mean/min/max", cal))
    }
    
    if (type == 'col') {
    matrix_output <- Matrix::t(matrix_output)
    }
    
    matrix_output <- as(matrix_output, "dgCMatrix")
    rownames(matrix_output) <- if (type == 'row') group_levels else colnames(matrix_input)
    colnames(matrix_output) <- if (type == 'row') colnames(matrix_input) else group_levels
    
    return(matrix_output)
}

matrix_group_scale <- function(matrix_input, group.by, type = "row") {
      scaled_matrix <- as.matrix(matrix_input)
      
      if (type == "row") {
        for (g in unique(group.by)) {
          row_idx <- which(group.by == g)
          scaled_matrix[row_idx, ] <- scale(scaled_matrix[row_idx, ])
        }
      } else if(type == "col"){ # 
        for (g in unique(group.by)) {
          col_idx <- which(group.by == g)
          scaled_matrix[, col_idx] <- Matrix::t(scale(Matrix::t(scaled_matrix[, col_idx])))
        }
      }
      scaled_matrix[is.na(scaled_matrix)] <- 0 
  return(scaled_matrix)
}

numeric_vector_group <- function(vec, n_groups,round_num = NULL) {
    n_groups <- as.integer(n_groups)
    if (length(vec) < n_groups) stop("")
    
    quantiles <- quantile(vec, probs = seq(0, 1, length.out = n_groups + 1), na.rm = TRUE)
    groups <- cut(vec, breaks = quantiles, include.lowest = TRUE, labels = FALSE)
    group_means <- tapply(vec, groups, mean, na.rm = TRUE)
    if(length(round_num)==1){group_means <- round(group_means, round_num)}
    group_means <- as.numeric(group_means[groups])
    names(group_means) <- groups
  return(group_means)
}
#Calculate pearson C.C. after filter NA
cor_myself <- function(matrix_input,return.p=T){
    matrix_output <- matrix(0,nrow = ncol(matrix_input),ncol=ncol(matrix_input))
    rownames(matrix_output) <- colnames(matrix_input)
    colnames(matrix_output) <- colnames(matrix_input)
    if(return.p){matrix.p <- matrix_output}
    for(i in colnames(matrix_input)){
        for(j in colnames(matrix_input)){
            use_index <- !(is.na(matrix_input[,i])|is.na(matrix_input[,j]))
            result <- cor.test(matrix_input[use_index,i],matrix_input[use_index,j])
            matrix_output[i,j] <- result$estimate
            if(return.p){matrix.p[i,j] <- result$p.value}
            }
    }
    if(return.p){return(list('R'=matrix_output,'p.value'=matrix.p))}else(return(matrix_output))
}
cor_twomatrix <- function(matrix_1,matrix_2,return.p=F){
    matrix_output <- matrix(0,nrow = ncol(matrix_1),ncol=ncol(matrix_2))
    rownames(matrix_output) <- colnames(matrix_1)
    colnames(matrix_output) <- colnames(matrix_2)
    matrix.p <- matrix_output
    n <- 1
    all_n <- ncol(matrix_1)*ncol(matrix_2)
    pb <- txtProgressBar(style=3)
    for(i in colnames(matrix_1)){
        for(j in colnames(matrix_2)){
            use_index <- !(is.na(matrix_1[,i])|is.na(matrix_2[,j]))
            result <- cor.test(matrix_1[use_index,i],matrix_2[use_index,j])
            matrix_output[i,j] <- result$estimate
            if(return.p){matrix.p[i,j] <- result$p.value}
            setTxtProgressBar(pb, n/all_n)
            n <- n+1
            }
    }
    close(pb)
    if(return.p){return(list('R'=matrix_output,'p.value'=matrix.p))}else(return(matrix_output))
}
Mean_normalize <- function(input_vector_raw){
    input_vector <- input_vector_raw
    up_vector <- input_vector[input_vector_raw>=mean(input_vector_raw)]
    down_vector <- input_vector[input_vector_raw<mean(input_vector_raw)]
    input_vector[input_vector_raw>=mean(input_vector_raw)] <- ((up_vector-min(up_vector))/(max(up_vector)-min(up_vector)))+1
    input_vector[input_vector_raw<mean(input_vector_raw)] <- (down_vector-min(down_vector))/(max(down_vector)-min(down_vector))
    return(input_vector)
}
Manual_mean_normalize <- function(input_vector_raw,manual_mean=0){
    input_vector <- input_vector_raw
    up_vector <- input_vector[input_vector_raw>=manual_mean]
    down_vector <- input_vector[input_vector_raw<manual_mean]
    input_vector[input_vector_raw>=manual_mean] <- ((up_vector-min(up_vector))/(max(up_vector)-min(up_vector)))+1
    input_vector[input_vector_raw<manual_mean] <- (down_vector-min(down_vector))/(max(down_vector)-min(down_vector))
    return(input_vector)
}
normalize_columns <- function(mat) {
  min_vals <- apply(mat, 2, min)
  max_vals <- apply(mat, 2, max)
  normalized_mat <- t((t(mat) - min_vals) / (max_vals - min_vals))
  return(normalized_mat)
}
normalize_row <- function(mat) {
  min_vals <- apply(mat, 1, min)
  max_vals <- apply(mat, 1, max)
  normalized_mat <- (mat - min_vals) / (max_vals - min_vals)
  return(normalized_mat)
}
order_matrix <- function(input_matrix,row_order){
    input_matrix <- scale(input_matrix)
    input_matrix[is.na(input_matrix)] <- 0
    return_df <- dplyr::bind_rows(lapply(colnames(input_matrix),function(x){
        max_value <- max(input_matrix[,x])
        max_label <- rownames(input_matrix)[input_matrix[,x]==max_value][1]
        df <- data.frame(var=x,label=max_label,value=max_value)
        return(df)}))
    rownames(return_df) <- return_df$var
    return_df <- dplyr::bind_rows(lapply(row_order,function(x){
        tmp_df <- return_df[return_df$label==x,]
        tmp_df <- tmp_df[order(tmp_df$value,decreasing =T),]
        return(tmp_df)
    }))
    return(return_df)
}
upset2inter <- function(upset_data){
    intersection <- lapply(1:nrow(upset_data),function(X){
        tmp <- upset_data[X,]
        intersection <- sort(colnames(tmp)[tmp!=0])
        intersection <- paste(intersection,collapse="/")
        return(intersection)
                                               })
    intersection <- as.data.frame(cbind(rownames(upset_data),intersection))
    intersection$V1 <- as.character(intersection$V1)
    intersection$intersection <- as.character(intersection$intersection)
    intersection$Number <- rowSums(upset_data)
    return(intersection)
}

upset2inter_s <- function(upset_data){
    intersection <- upset_data
    for(i in colnames(intersection)){
        intersection[intersection[,i]!=0,i] <- i
        intersection[intersection[,i]==0,i] <- ''
    }
    intersection <- apply(X = intersection,MARGIN = 1,FUN = function(row) paste(row[row != ""], collapse = "/"))
                                                   
    intersection <- as.data.frame(cbind(rownames(upset_data),intersection))
    intersection$V1 <- as.character(intersection$V1)
    intersection$intersection <- as.character(intersection$intersection)
    intersection$Number <- rowSums(upset_data)
    return(intersection)
}
                          
ztest <- function(vector,mu){
    vector_add <- scale(c(vector,mu))
    return(1-pnorm(vector_add,mean=0,sd=1)[length(vector_add)])
}