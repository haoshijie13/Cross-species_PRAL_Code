Run_GO <- function(gene_list,species,change_human=TRUE,ont='BP'){
    if(change_human){
        gene_symbols <- change_gene(gene_list,species,'human')
        gene_symbols[is.na(gene_symbols)] <- gene_list[is.na(gene_symbols)]
                    }else(
        gene_symbols <- gene_list
                    )
    
    gene_entrez <- bitr(gene_symbols,
                        fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)
    gene <- gene_entrez$ENTREZID
    go_enrich_bp <- enrichGO(gene = gene,
                             OrgDb = org.Hs.eg.db,
                             ont = ont,
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.05)
}

Run_KEGG <- function(gene_list,species,change_human=TRUE){
    if(change_human){
        gene_symbols <- change_gene(gene_list,species,'human')
        gene_symbols[is.na(gene_symbols)] <- gene_list[is.na(gene_symbols)]
                    }else(
        gene_symbols <- gene_list
                    )
    gene_entrez <- bitr(gene_symbols,
                        fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)
    gene <- gene_entrez$ENTREZID
    go_enrich_bp <- enrichKEGG(gene = gene,
                               organism = 'hsa',
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.05,
                               use_internal_data =T)
}

get_Description_gene <- function(plot_df,Description){
    gene_entrez <- strsplit(plot_df[ plot_df$Description==Description,'geneID'],'/')[[1]]
    gene_symbols <- bitr(gene_entrez,
                        fromType = "ENTREZID",
                        toType = "SYMBOL",
                        OrgDb = org.Hs.eg.db)
    return(gene_symbols$SYMBOL)
}

get_ENTREZID_gene <- function(gene_list){
    gene_entrez <- strsplit(gene_list,'/')[[1]]
    gene_symbols <- bitr(gene_entrez,
                        fromType = "ENTREZID",
                        toType = "SYMBOL",
                        OrgDb = org.Hs.eg.db)
    return(gene_symbols$SYMBOL)
}

plot_GO <- function(plot_df,Description,p.adjust=FALSE){
    plot_df <- plot_df[plot_df$Description%in%Description,]
    if(class(plot_df$GeneRatio)=='character'){
    plot_df$GeneRatio <- sapply(plot_df$GeneRatio, function(x){eval(parse(text=x))})}
    plot_df$Description <- factor(plot_df$Description,levels = plot_df[order(plot_df$GeneRatio),'Description'])
    if(p.adjust){plot_df$pvalue <- plot_df$p.adjust}
    p1 <- ggplot()+
    geom_point(data=plot_df,aes(x=GeneRatio,y=Description,color=-log10(pvalue),size=Count))+
    scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdYlGn")))+
    scale_size(range = c(5,8))+
    theme_classic()+
    theme(text=element_text(size=16))
    return(p1)
}
plot_GO_bar <- function(plot_df,Description,p.adjust=FALSE){
    plot_df <- plot_df[plot_df$Description%in%Description,]
    if(class(plot_df$GeneRatio)=='character'){
    plot_df$GeneRatio <- sapply(plot_df$GeneRatio, function(x){eval(parse(text=x))})}
    plot_df$Description <- factor(plot_df$Description,levels = plot_df[order(plot_df$GeneRatio),'Description'])
    if(p.adjust){plot_df$pvalue <- plot_df$p.adjust}
    p1 <- ggplot()+
    geom_bar(data=plot_df,aes(x=GeneRatio,y=Description,fill=-log10(pvalue)),stat='identity')+
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdYlGn")))+
    theme_classic()+
    theme(text=element_text(size=16))
    return(p1)
}
plot_GO_double <- function(plot_df_left,plot_df_right,Description,order=TRUE,p.adjust=FALSE,return_order=FALSE){
    
    
    plot_df_left <- plot_df_left[plot_df_left$Description%in%Description,]
    plot_df_left$show <- 'left'
    if(class(plot_df_left$GeneRatio)=='character'){
        plot_df_left$GeneRatio <- sapply(plot_df_left$GeneRatio, function(x){eval(parse(text=x))})
    }else{
        plot_df_left$GeneRatio <- plot_df_left$GeneRatio
    }
                                     
    plot_df_right <- plot_df_right[plot_df_right$Description%in%Description,]
    plot_df_right$show <- 'right'
    if(class(plot_df_right$GeneRatio)=='character'){
    plot_df_right$GeneRatio <- sapply(plot_df_right$GeneRatio, function(x){eval(parse(text=x))})}

    plot_df <- rbind(plot_df_left,plot_df_right)
    if(p.adjust){plot_df$pvalue <- plot_df$p.adjust}
    plot_df <- plot_df[plot_df$pvalue<0.05,]
    plot_df$change_p <- -log10(plot_df$pvalue)
    plot_df[plot_df$show=='left','GeneRatio'] <- -plot_df[plot_df$show=='left','GeneRatio']
    
    if(order){
        Description_value <- sapply(Description,function(x){
            tmp_df <- plot_df[plot_df$Description==x,]
            tmp_max <- max(abs(tmp_df$GeneRatio))
            return(tmp_df[abs(tmp_df$GeneRatio)==tmp_max,'GeneRatio'][1])})
        plot_df$Description <- factor(plot_df$Description,levels = Description[order(Description_value)])
        
    }else{
        plot_df$Description <- factor(plot_df$Description,levels = Description)
    }
    
                  
    p1 <- ggplot()+
    geom_segment(data=plot_df,aes(x=GeneRatio,y=Description,xend=0,yend=Description,color=change_p),lwd=1)+
    geom_point(data=plot_df,aes(x=GeneRatio,y=Description,color=change_p,size=Count))+
    geom_vline(xintercept = 0,lwd=0.5)+
    #scale_color_manual(values=c('#2873B3','#A14462'),breaks=c('left','right'))+
    scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdYlGn")),'-log10(pvalue)')+
    scale_size(range = c(2,8))+
    theme_classic()+
    theme(text=element_text(size=16))+
    xlab('GeneRatio')+
    xlim(c(-max(plot_df$GeneRatio),max(plot_df$GeneRatio)))
    if(return_order){return(Description[order(Description_value)])}else{return(p1)} 
}

plot_GO_double_color <- function(plot_df_left,plot_df_right,Description,order=TRUE,p.adjust=FALSE,return_order=FALSE){
    
    
    plot_df_left <- plot_df_left[plot_df_left$Description%in%Description,]
    plot_df_left$show <- 'left'
    if(class(plot_df_left$GeneRatio)=='character'){
        plot_df_left$GeneRatio <- sapply(plot_df_left$GeneRatio, function(x){eval(parse(text=x))})
    }else{
        plot_df_left$GeneRatio <- plot_df_left$GeneRatio
    }
                                     
    plot_df_right <- plot_df_right[plot_df_right$Description%in%Description,]
    plot_df_right$show <- 'right'
    if(class(plot_df_right$GeneRatio)=='character'){
    plot_df_right$GeneRatio <- sapply(plot_df_right$GeneRatio, function(x){eval(parse(text=x))})}

    plot_df <- rbind(plot_df_left,plot_df_right)
    if(p.adjust){plot_df$pvalue <- plot_df$p.adjust}
    plot_df <- plot_df[plot_df$pvalue<0.05,]
    plot_df$change_p <- -log10(plot_df$pvalue)
    plot_df[plot_df$show=='left','change_p'] <- -plot_df[plot_df$show=='left','change_p']
    if(order){
        Description_value <- sapply(Description,function(x){
            tmp_df <- plot_df[plot_df$Description==x,]
            tmp_max <- max(abs(tmp_df$change_p))
            return(tmp_df[abs(tmp_df$change_p)==tmp_max,'change_p'][1])})
        plot_df$Description <- factor(plot_df$Description,levels = Description[order(Description_value)])
        
    }else{
        plot_df$Description <- factor(plot_df$Description,levels = Description)
    }
    
                  
    p1 <- ggplot()+
    geom_segment(data=plot_df,aes(x=change_p,y=Description,xend=0,yend=Description,color=show),lwd=1)+
    geom_point(data=plot_df,aes(x=change_p,y=Description,color=show,size=Count))+
    geom_vline(xintercept = 0,lwd=0.5)+
    scale_color_manual(values=c('#2873B3','#A14462'),breaks=c('left','right'))+
    #scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(10,"RdYlGn")),'-log10(pvalue)')+
    scale_size(range = c(2,8))+
    theme_classic()+
    theme(text=element_text(size=16))+
    xlab('-log10(pvalue)')+
    xlim(c(-max(plot_df$change_p),max(plot_df$change_p)))
    if(return_order){return(Description[order(Description_value)])}else{return(p1)} 
}