plot_spatial_feature <- function(obj,features,assays='RNA',slot='data',
                                 coor_use=c('x','y'),
                                 max.cutoff=1,min.cutoff=0,vmid=NULL,limits=NULL,
                                 smooth_function='spatial',
                                 smooth=FALSE,knn=25,round=2,
                                 border=FALSE,border_use=Border,
                                 size=3,height=NULL,width=NULL,
                                 raster='scatter',
                                 color=c('gray50','gray','gray97','red','darkred')){
    obj@meta.data[,c('x','y')] <- obj@meta.data[,coor_use]
    plot.list <- lapply(features,function(gene){
        if(smooth){
            if(smooth_function=='spatial'){
                exp <- smooth_kNN(obj@meta.data[,c('x','y')],
                                  obj@meta.data[,c('x','y')],
                                  slot(obj@assays[[assays]],slot)[gene,],
                                  round=round,knn=knn)
            }else if(smooth_function=='magic'){
                exp <- magic_knn(obj@reductions$pca@cell.embeddings,
                                 t(as.matrix(slot(obj@assays[[assays]],slot)[gene,,drop=FALSE])),
                                 knn=knn,
                                 round=round)
            }
        }else{exp <- slot(obj@assays[[assays]],slot)[gene,]}

        if(min.cutoff!=0 | max.cutoff!=1){
            max.cutoff <- quantile(exp[!is.na(exp)],max.cutoff)
            min.cutoff <- quantile(exp[!is.na(exp)],min.cutoff)
            exp[exp>max.cutoff&(!is.na(exp))] <- max.cutoff
            exp[exp<min.cutoff&(!is.na(exp))] <- min.cutoff}
        obj$exp <- exp

        xmid <- (max(obj$x)+min(obj$x))/2
        ymid <- (max(obj$y)+min(obj$y))/2
        if(is.null(height)){height <- max(obj$y)-min(obj$y)+10}
        if(is.null(width)){width <- max(obj$x)-min(obj$x)+10}
        
        p1 <- ggplot()
        if(raster=='scatter'){
            p1 <- p1+scattermore::geom_scattermore(data=obj@meta.data,aes(x=x,y=y,color=exp),pointsize = size)+
            scale_color_gradientn(colours = color,name = gene,na.value = "gray50",limits=limits)
        }else if(raster=='tile'){
            p1 <- p1+geom_raster(data=obj@meta.data,aes(x=x,y=y,fill=exp))+
            scale_fill_gradientn(colours = color,name = gene,na.value = "gray50",limits=limits)
        }
        p1 <- p1+
        theme_void()+
        xlim(xmid-width/2,xmid+width/2)+
        ylim(ymid-height/2,ymid+height/2)+
        theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'))+
        coord_fixed()
        
        if(!is.null(vmid) & raster=='scatter'){
        p1 <- p1+scale_color_gradientn(colours = color,name = gene,
                             values=scales::rescale(c(min(exp),quantile(exp,vmid),max(exp))),na.value = "gray",limits=limits)}
        if(!is.null(vmid) & raster=='tile'){
        p1 <- p1+scale_fill_gradientn(colours = color,name = gene,
                             values=scales::rescale(c(min(exp),quantile(exp,vmid),max(exp))),na.value = "gray",limits=limits)}
        
        if(border){p1 <- p1+geom_segment(data=border_use,aes(x=X,y=Y,xend=X1,yend=Y1),lwd=0.8)}
        return(p1)
    })
    p <- plot.list[[1]]
    if(length(features)>1){
    for(i in 2:length(features)){p <- p+plot.list[[i]]}
    return(p)
        }else{
    return(p)}
                        }

plot_spatial_col <- function(obj,col_name,coor_use=c('x','y'),
                             border=FALSE,border_use=Border,
                             smooth=FALSE,knn=25,round=2,
                             max.cutoff=1,min.cutoff=0,vmid=NULL,limits=NULL,
                             size=3,height=NULL,width=NULL,
                             raster='scatter',
                             color=rev(RColorBrewer::brewer.pal(11,'Spectral'))){
    obj[,c('x','y')] <- obj[,coor_use]
    xmid <- (max(obj$x)+min(obj$x))/2
    ymid <- (max(obj$y)+min(obj$y))/2
    exp <- obj[,col_name]
    if(min.cutoff!=0 | max.cutoff!=1){
        max.cutoff <- quantile(exp[!is.na(exp)],max.cutoff)
        min.cutoff <- quantile(exp[!is.na(exp)],min.cutoff)
        exp[exp>max.cutoff&(!is.na(exp))] <- max.cutoff
        exp[exp<min.cutoff&(!is.na(exp))] <- min.cutoff}
    if(smooth){
        exp <- smooth_kNN(obj[,c('x','y')],
                          obj[,c('x','y')],
                          exp,
                          round=round,knn=knn)}
    if(is.null(height)){height <- max(obj$y)-min(obj$y)+10}
    if(is.null(width)){width <- max(obj$x)-min(obj$x)+10}
    obj$exp <- exp
    p1 <- ggplot()
    if(raster=='scatter'){
        p1 <- p1+scattermore::geom_scattermore(data=obj,aes(x=x,y=y,color=exp),pointsize = size)+
        scale_color_gradientn(colours = color,name = col_name,na.value = "gray50",limits=limits)+
        xlim(xmid-width/2,xmid+width/2)+
        ylim(ymid-height/2,ymid+height/2)+
        theme_void()+
        theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'))+
        coord_fixed()
    }else if(raster=='tile'){
        p1 <- p1+geom_raster(data=obj,aes(x=x,y=y,fill=exp))+
        scale_fill_gradientn(colours = color,name = col_name,na.value = "gray",limits=limits)+
        xlim(xmid-width/2,xmid+width/2)+
        ylim(ymid-height/2,ymid+height/2)+
        theme_void()+
        theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'))+
        coord_fixed()
    }
    
    if(!is.null(vmid) & raster=='scatter'){
    p1 <- p1+scale_color_gradientn(colours = color,name = col_name,
                         values=scales::rescale(c(min(exp),quantile(exp,vmid),max(exp))),na.value = "gray50",limits=limits)}
    if(!is.null(vmid) & raster=='tile'){
    p1 <- p1+scale_fill_gradientn(colours = color,name = col_name,
                         values=scales::rescale(c(min(exp),quantile(exp,vmid),max(exp))),na.value = "gray50",limits=limits)}
    
    if(border){p1 <- p1+geom_segment(data=border_use,aes(x=X,y=Y,xend=X1,yend=Y1),lwd=0.8)}
    return(p1)
}

save_png_plot <- function(plot,dir_path,name,height=200,width=500,save='png'){
    if(!dir.exists(dir_path)){dir.create(dir_path)}
    if(save=='ggsave'){ggsave(paste0(dir_path, "/", name, ".png"),plot + NoLegend(), height = height, width = width,bg = "transparent", units = 'px')
                      }else if(save=='png'){
        png(paste0(dir_path, "/", name, ".png"), height = height, width = width,bg = "transparent")
        print(plot + NoLegend())
        dev.off()}
    legend <- cowplot::get_legend(plot)
    pdf(paste0(dir_path,'/',name,'.legend.pdf'))
    print(grid::grid.draw(legend))
    dev.off()
}

save_real_plot <- function(save_plot,dir_path,name,plot_x_vector,plot_y_vector,edge_add=50){
    if(!dir.exists(dir_path)){dir.create(dir_path)}
    
    x_range <- c(min(plot_x_vector)-edge_add, max(plot_x_vector)+edge_add)
    y_range <- c(min(plot_y_vector)-edge_add, max(plot_y_vector)+edge_add)
    x_length <- diff(x_range)
    y_length <- diff(y_range)
    x_offset <- min(x_range)
    y_offset <- min(y_range)

    resolution <- 1
    width_px <- x_length
    height_px <- y_length
    dpi <- 96

    save_plot <- save_plot +
    coord_cartesian(xlim = x_range, ylim = y_range, expand = FALSE)+
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0),legend.position = "none")

    ggsave(
      filename = paste0(dir_path, "/", name,'_',round(x_offset,2),'_',round(y_offset,2),".tif"),
      plot = save_plot,
      width = width_px / dpi,  
      height = height_px / dpi,
      dpi = dpi,               
      device = "tiff",         
      units = "in",            
      bg = "white",            
      compression = "lzw"      
    ) 
}

read_real_mask_plot <- function(path,x_offset,y_offset){
    tif_array <- tiff::readTIFF(path, native = FALSE)
    tif_array <- tif_array*255
    tif_array <- tif_array[nrow(tif_array):1, ]
    tif_df <- reshape2::melt(tif_array)
    colnames(tif_df) <- c("y", "x", "mask")
    tif_df <- tif_df[tif_df$mask!=0,]
    tif_df$x <- tif_df$x + x_offset
    tif_df$y <- tif_df$y + y_offset
    return(tif_df)
}

df_add_mask <- function(df_raw,df_mask,raw_coord=c('x','y'),mask_coord=c('x','y'),mask_col='mask',imputation=TRUE){
    rownames(df_mask) <- paste0(round(df_mask[,mask_coord[1]],0),'_',round(df_mask[,mask_coord[2]],0))
    raw_index <- paste0(round(df_raw[,raw_coord[1]],0),'_',round(df_raw[,raw_coord[2]],0))
    df_mask <- df_mask[rownames(df_mask)%in%raw_index ,]
    df_raw[,mask_col] <- df_mask[raw_index,mask_col]
    if(imputation&(TRUE%in%is.na(df_raw[,mask_col]))){
        df_raw[is.na(df_raw[,mask_col]),mask_col] <- winner_kNN(initial_df = df_raw[!is.na(df_raw[,mask_col]),raw_coord],
                                                        query_df = df_raw[is.na(df_raw[,mask_col]),raw_coord],
                                                        sm_vector = df_raw[!is.na(df_raw[,mask_col]),mask_col])}
    return(df_raw)
}