plot_quasirandom <- function(meta.data,X,Y,group_down=TRUE,group_size=500,y_dashed=c(0),raster=FALSE,split=FALSE){
    set.seed(123)
    df_tmp <- meta.data
    if(split){
        y_split <- sapply(unique(df_tmp$`Pr-Al-Type`),function(x){
            PrAl_type_df <- df_tmp[df_tmp$`Pr-Al-Type`==x,]
            return(min(PrAl_type_df$`Pr-Al-Index`))})
        y_split <- unique(unlist(y_split))
        y_split <- c(y_split,1)}
    
    if(group_down){
        df_tmp <- dplyr::bind_rows(lapply(unique(df_tmp[,X]),function(x){
            group_df <- df_tmp[df_tmp[,X]==x,]
            return(group_df[sample(rownames(group_df),size = group_size),])}))}
    df_tmp$group <- df_tmp[,X]
    df_tmp$value <- df_tmp[,Y]

    p <- ggplot(data=df_tmp)
    if(raster){
        p <- p+geom_quasirandom_rast( aes(x = group, y = value,color=value),dodge.width = 0.9,size=0.5)
    }else{
        p <- p+geom_quasirandom( aes(x = group, y = value,color=value),dodge.width = 0.9,size=0.5)
    }
    p <- p+
    scale_color_gradientn(colours  = c('#2873B3','#2873B3','gray','#A14462','#A14462'),limits=c(-1,1),name = X)+
    theme_classic()+
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust = 1))+
    NoLegend()+
    ylim(-1,1)+
    ylab(Y)+
    xlab(X)
    p <- p+geom_hline(yintercept = y_dashed,linetype='dashed')
    return(p)
}
plot_cor_raster <- function(meta.data,X,Y,var,fitting=FALSE,color=rev(RColorBrewer::brewer.pal(11,"Spectral")[-6])){
    meta.data$X <- meta.data[,X]
    meta.data$Y <- meta.data[,Y]
    meta.data$var <- meta.data[,var]
    p <- ggplot(data=meta.data,aes(x=X,y=Y))+
    scattermore::geom_scattermore(aes(color=var),
                     pointsize = 3)+
    scale_color_gradientn(colours = color,name=var)+
    xlab(X)+
    ylab(Y)+
    theme_classic()+
    theme(aspect.ratio=1,text=element_text(size=16))+
    ggtitle(paste0('R=',round(cor(meta.data$X,meta.data$Y),2),' p=',signif(cor.test(meta.data$X,meta.data$Y)$p.value,3)))
    if(fitting){p <- p+geom_smooth(method = 'lm',lwd=2,se = T,color='gray')}
    return(p)
}
plot_cor_dot <- function(meta.data,X,Y,var,fitting=FALSE,method='pearson',color=rev(RColorBrewer::brewer.pal(11,"Spectral")[-6])){
    meta.data$X <- meta.data[,X]
    meta.data$Y <- meta.data[,Y]
    meta.data$var <- meta.data[,var]
    p <- ggplot(data=meta.data,aes(x=meta.data$X,y=meta.data$Y))+
    geom_point(aes(color=meta.data$var),show.legend = T, shape=21, size = 3, stroke = 1.5)+
    scale_color_gradientn(colours = color,name=var)+
    xlab(X)+
    ylab(Y)+
    theme_classic()+
    theme(aspect.ratio=1,text=element_text(size=16))+
    ggtitle(paste0('R=',round(cor(meta.data$X,meta.data$Y,method=method),2),' ',
                   'p=',signif(cor.test(meta.data$X,meta.data$Y,method=method)$p.value,3)))
    if(fitting){p <- p+geom_smooth(method = 'lm',lwd=2,se = T,color='gray')}
    return(p)
}

plot_cor_raster_discrete <- function(meta.data,X,Y,var,show_dot=TRUE,show_text=TRUE,fitting=FALSE,color=ggsci::pal_igv()(50)){
    meta.data$X <- meta.data[,X]
    meta.data$Y <- meta.data[,Y]
    meta.data$var <- meta.data[,var]
    agg_meta.data <- aggregate(list('X'=meta.data$X,'Y'=meta.data$Y),by=list('var'=meta.data$var),mean)
    p <- ggplot(data=meta.data,aes(x=X,y=Y))+
    scattermore::geom_scattermore(aes(color=var),pointsize = 3)+
    scale_color_manual(values = color,name=var)+
    xlab(X)+
    ylab(Y)+
    theme_classic()+
    theme(aspect.ratio=1,text=element_text(size=16))+
    ggtitle(paste0('R=',round(cor(meta.data$X,meta.data$Y),2),' p=',signif(cor.test(meta.data$X,meta.data$Y)$p.value,3)))
    if(fitting){p <- p+geom_smooth(method = 'lm',lwd=2,se = T,color='gray')}
    if(show_dot){p <- p+geom_point(data=agg_meta.data,aes(x=X,y=Y,color=var))}
    if(show_text){p <- p+geom_text(data=agg_meta.data,aes(x=X,y=Y,label=var))}
    return(p)
}
