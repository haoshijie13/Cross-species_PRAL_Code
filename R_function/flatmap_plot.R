#flat_mar <- read_sf('data/MRI_shp_file/Marmoset.shp')
#flat_mac <- read_sf('data/MRI_shp_file/Macaque.shp')
#flat_hum <- read_sf('data/MRI_shp_file/Human.shp')
#flat_mou <- read_sf('data/MRI_shp_file/Mouse.shp')

plot_mar_flatmap_col <- function(obj.meta,flatmap,max.cutoff=1,min.cutoff=0,plot_names,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){
    exp <- obj.meta[flatmap$region,plot_names]
    if(min.cutoff!=0 | max.cutoff!=1){
        max.cutoff <- quantile(exp[!is.na(exp)],max.cutoff)
        min.cutoff <- quantile(exp[!is.na(exp)],min.cutoff)
        exp[exp>max.cutoff&(!is.na(exp))] <- max.cutoff
        exp[exp<min.cutoff&(!is.na(exp))] <- min.cutoff}
    p1 <- ggplot(flatmap)+
    geom_sf(aes(fill=exp),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=plot_names)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
    return(p1)
}
plot_mac_flatmap_col <- function(obj.meta,flatmap,max.cutoff=1,min.cutoff=0,plot_names,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){
    tmp <- obj.meta
    tmp[tmp$region%in%c('Iam','Iapm'),'region'] <- 'Iam-Iapm'
    tmp[tmp$region%in%c('PEc','PEci'),'region'] <- 'PEc-PEci'
    tmp[tmp$region%in%c('Ia','Id'),'region'] <- 'Ia-Id'
    tmp[tmp$region%in%c('PG','Opt','DP'),'region'] <- 'PG-Opt-DP'
    tmp[tmp$region%in%c('36r','36p'),'region'] <- '36r-36p'
    tmp[tmp$region%in%c('V4','V4t'),'region'] <- 'V4-V4t'
    tmp[tmp$region%in%c('PF','PFG'),'region'] <- 'PF-PFG'
    tmp[tmp$region%in%c('13a','13b'),'region'] <- '13a-13b'
    tmp[tmp$region%in%c('EO','EI','ELr','ELc'),'region'] <- 'EO-EI-ELr-ELc'
    tmp <- aggregate(list('feature'=tmp[,plot_names]),by=list('region'=tmp$region),mean)
    rownames(tmp) <- tmp$region

    exp <- tmp[flatmap$region,'feature']
    if(min.cutoff!=0 | max.cutoff!=1){
        max.cutoff <- quantile(exp[!is.na(exp)],max.cutoff)
        min.cutoff <- quantile(exp[!is.na(exp)],min.cutoff)
        exp[exp>max.cutoff&(!is.na(exp))] <- max.cutoff
        exp[exp<min.cutoff&(!is.na(exp))] <- min.cutoff}
    
    p1 <- ggplot(flatmap) +
    geom_sf(aes(fill=exp),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=plot_names)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
    return(p1)
}

plot_hum_flatmap_col <- function(obj.meta,flatmap,plot_names,max.cutoff=1,min.cutoff=0,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){
    exp <- obj.meta[flatmap$region,plot_names]
    if(min.cutoff!=0 | max.cutoff!=1){
        max.cutoff <- quantile(exp[!is.na(exp)],max.cutoff)
        min.cutoff <- quantile(exp[!is.na(exp)],min.cutoff)
        exp[exp>max.cutoff&(!is.na(exp))] <- max.cutoff
        exp[exp<min.cutoff&(!is.na(exp))] <- min.cutoff}

    p1 <- ggplot(flatmap) +
    geom_sf(aes(fill=exp),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=plot_names)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
    return(p1)
}

plot_mou_flatmap_col <- function(obj.meta,flatmap,plot_names,max.cutoff=1,min.cutoff=0,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){
    exp <- obj.meta[flatmap$region,plot_names]
    if(min.cutoff!=0 | max.cutoff!=1){
        max.cutoff <- quantile(exp[!is.na(exp)],max.cutoff)
        min.cutoff <- quantile(exp[!is.na(exp)],min.cutoff)
        exp[exp>max.cutoff&(!is.na(exp))] <- max.cutoff
        exp[exp<min.cutoff&(!is.na(exp))] <- min.cutoff}

    p1 <- ggplot(flatmap) +
    geom_sf(aes(fill=exp),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=plot_names)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
return(p1)
}

#plot_features
plot_mar_flatmap_feature <- function(obj,flatmap,feature,assays='RNA',slot='scale.data',max.cutoff=1,min.cutoff=0,vmid=FALSE,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){ 
    exp <- slot(obj@assays[[assays]],slot)[feature,]
    if(min.cutoff!=0 | max.cutoff!=1){
            max.cutoff <- quantile(exp[!is.na(exp)],max.cutoff)
            min.cutoff <- quantile(exp[!is.na(exp)],min.cutoff)
            exp[exp>max.cutoff&(!is.na(exp))] <- max.cutoff
            exp[exp<min.cutoff&(!is.na(exp))] <- min.cutoff}
    obj$exp <- exp
    p1 <- ggplot(flatmap) +
    geom_sf(aes(fill=obj@meta.data[flatmap$region,'exp']),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=feature)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
    
    if(vmid!=FALSE){
        p1 <- p1+scale_fill_gradientn(colours = color,name = gene,
                             values=scales::rescale(c(min(exp),quantile(exp,vmid),max(exp))),na.value = "gray")}
    
    return(p1)
}
plot_mac_flatmap_feature <- function(obj,flatmap,feature,assays='RNA',slot='scale.data',max.cutoff=1,min.cutoff=0,vmid=FALSE,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){
    exp <- slot(obj@assays[[assays]],slot)[feature,]
    if(min.cutoff!=0 | max.cutoff!=1){
            max.cutoff <- quantile(exp[!is.na(exp)],max.cutoff)
            min.cutoff <- quantile(exp[!is.na(exp)],min.cutoff)
            exp[exp>max.cutoff&(!is.na(exp))] <- max.cutoff
            exp[exp<min.cutoff&(!is.na(exp))] <- min.cutoff}
    obj$exp <- exp
    tmp <- obj@meta.data
    tmp[tmp$region%in%c('Iam','Iapm'),'region'] <- 'Iam-Iapm'
    tmp[tmp$region%in%c('PEc','PEci'),'region'] <- 'PEc-PEci'
    tmp[tmp$region%in%c('Ia','Id'),'region'] <- 'Ia-Id'
    tmp[tmp$region%in%c('PG','Opt','DP'),'region'] <- 'PG-Opt-DP'
    tmp[tmp$region%in%c('36r','36p'),'region'] <- '36r-36p'
    tmp[tmp$region%in%c('V4','V4t'),'region'] <- 'V4-V4t'
    tmp[tmp$region%in%c('PF','PFG'),'region'] <- 'PF-PFG'
    tmp[tmp$region%in%c('13a','13b'),'region'] <- '13a-13b'
    tmp[tmp$region%in%c('EO','EI','ELr','ELc'),'region'] <- 'EO-EI-ELr-ELc'
    tmp <- aggregate(list('feature'=tmp[,'exp']),by=list('region'=tmp$region),mean)
    rownames(tmp) <- tmp$region
    
    p1 <- ggplot(flatmap) +
    geom_sf(aes(fill=tmp[flatmap$region,]$feature),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=feature)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
    if(vmid!=FALSE){
        p1 <- p1+scale_fill_gradientn(colours = color,name = gene,
                             values=scales::rescale(c(min(exp),quantile(exp,vmid),max(exp))),na.value = "gray")}
    return(p1)
}

plot_hum_flatmap_feature <- function(obj,flatmap,feature,assays='RNA',slot='scale.data',max.cutoff=1,min.cutoff=0,vmid=FALSE,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){
    exp <- slot(obj@assays[[assays]],slot)[feature,]
    if(min.cutoff!=0 | max.cutoff!=1){
            max.cutoff <- quantile(exp[!is.na(exp)],max.cutoff)
            min.cutoff <- quantile(exp[!is.na(exp)],min.cutoff)
            exp[exp>max.cutoff&(!is.na(exp))] <- max.cutoff
            exp[exp<min.cutoff&(!is.na(exp))] <- min.cutoff}
    obj$exp <- exp
    tmp <- obj@meta.data
    p1 <- ggplot(flatmap) +
    geom_sf(aes(fill=tmp[flatmap$region,'exp']),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=feature)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
    if(vmid!=FALSE){
        p1 <- p1+scale_fill_gradientn(colours = color,name = gene,
                             values=scales::rescale(c(min(exp),quantile(exp,vmid),max(exp))),na.value = "gray")}
    return(p1)
}

plot_mou_flatmap_feature <- function(obj,flatmap,feature,assays='RNA',slot='scale.data',max.cutoff=1,min.cutoff=0,vmid=FALSE,color=rev(RColorBrewer::brewer.pal(11,'Spectral')[-6])){
    exp <- slot(obj@assays[[assays]],slot)[feature,]
    if(min.cutoff!=0 | max.cutoff!=1){
            max.cutoff <- quantile(exp[!is.na(exp)],max.cutoff)
            min.cutoff <- quantile(exp[!is.na(exp)],min.cutoff)
            exp[exp>max.cutoff&(!is.na(exp))] <- max.cutoff
            exp[exp<min.cutoff&(!is.na(exp))] <- min.cutoff}
    obj$exp <- exp
    tmp <- obj@meta.data
    p1 <- ggplot(flatmap) +
    geom_sf(aes(fill=tmp[flatmap$region,'exp']),color=NA,show.legend = T)+
    scale_fill_gradientn(colours = color,name=feature)+
    theme_void()+
    theme(text=element_text(size=16))+
    theme(legend.position = "bottom")
    if(vmid!=FALSE){
        p1 <- p1+scale_fill_gradientn(colours = color,name = gene,
                             values=scales::rescale(c(min(exp),quantile(exp,vmid),max(exp))),na.value = "gray")}
    return(p1)
}
