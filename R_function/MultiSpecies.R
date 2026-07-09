PlotMutiSpeciesSTAssays <- function(obj,
                                    feature_name,
                                    max.cutoff=0.97,
                                    min.cutoff=0,
                                    limits = NULL,
                                    assays = 'RNA',
                                    slot = 'data',
                                    plot_id = c(1,2,3),
                                    species = NULL,
                                    smooth = FALSE,
                                    knn = 30,
                                    round = 1,
                                    use_color=rev(c("#ff6b35","#f7c59f","#efefd0","#004e89","#1a659e"))){

    obj$var <- slot(obj@assays[[assays]], slot)[feature_name, ]
    plot_df <- obj@meta.data

    plot_config <- list(
        `1` = list(
            marm = list(
                list(sample = 'marm_T454', hemisphere = NULL, size = 3.5)
            ),
            mous = list(
                list(sample = 'mous_T313', hemisphere = NULL, size = 3.5)
            ),
            bird = list(
                list(sample = 'bird_a1', hemisphere = NULL, size = 3.5),
                list(sample = 'bird_d4', hemisphere = NULL, size = 3.5),
                list(sample = 'bird_a5', hemisphere = NULL, size = 3.5)
            )
        ),
        `2` = list(
            turt = list(
                list(sample = 'turt_a2', hemisphere = 'left',  size = 4),
                list(sample = 'turt_a2', hemisphere = 'right', size = 4),
                list(sample = 'turt_f5', hemisphere = 'left',  size = 4),
                list(sample = 'turt_f5', hemisphere = 'right', size = 4)
            ),
            axol = list(
                list(sample = 'axol_all', hemisphere = 'left',  size = 6),
                list(sample = 'axol_all', hemisphere = 'right', size = 4)
            )
        ),
        `3` = list(
            lung = list(
                list(sample = 'lung_all', hemisphere = 'left_C',  size = 4),
                list(sample = 'lung_all', hemisphere = 'right_C', size = 4)
            ),
            fish = list(
                list(sample = 'fish_D2', hemisphere = NULL, size = 10),
                list(sample = 'fish_D4', hemisphere = NULL, size = 10)
            ),
            lamp = list(
                list(sample = 'lamp_4',  hemisphere = NULL, size = 6),
                list(sample = 'lamp_5',  hemisphere = NULL, size = 6),
                list(sample = 'lamp_34', hemisphere = NULL, size = 6)
            )
        )
    )

    all_species <- c("marm", "mous", "bird", "turt", "axol", "lung", "fish", "lamp")

    if (is.null(species) || "all" %in% species) {
        species <- all_species
    }

    species <- intersect(all_species, species)
    if (length(species) == 0) {
        stop("species must be one or more of: marm, mous, bird, turt, axol, lung, fish, lamp")
    }

    plot_list <- list()

    for (pid in plot_id) {
        pid_chr <- as.character(pid)
        cfg_sub <- plot_config[[pid_chr]]
        species_in_panel <- intersect(names(cfg_sub), species)

        if (length(species_in_panel) == 0) {
            next
        }

        p_list <- list()

        for (sp in species_in_panel) {
            for (cfg in cfg_sub[[sp]]) {
                df_sub <- plot_df[plot_df$sample == cfg$sample, ]
                if (!is.null(cfg$hemisphere)) {
                    df_sub <- df_sub[df_sub$hemisphere == cfg$hemisphere, ]
                }

                p_list[[length(p_list) + 1]] <- plot_spatial_col(
                    df_sub,
                    col_name = 'var',
                    color = use_color,
                    smooth = smooth,
                    knn = knn,
                    round = round,
                    max.cutoff = max.cutoff,
                    min.cutoff = min.cutoff,
                    size = cfg$size,
                    limits = limits
                )
            }
        }

        if (length(p_list) > 0) {
            if (pid == 3) {
                widths_full <- c(1,1,0.7,0.7,0.7,0.7,0.7)
                plot_list[[pid_chr]] <- Reduce(`+`, p_list) +
                    plot_layout(ncol = length(p_list), widths = widths_full[seq_along(p_list)])
            } else {
                plot_list[[pid_chr]] <- Reduce(`+`, p_list) +
                    plot_layout(ncol = length(p_list))
            }
        }
    }

    if (length(plot_list) == 0) {
        stop("No plots available for the selected species and plot_id.")
    }

    options(repr.plot.width = 22, repr.plot.height = 10)
    return(Reduce(`/`, plot_list))
}

PlotMutiSpeciesSTMeta <- function(meta.data,
                                  col_name,
                                  max.cutoff=0.97,
                                  min.cutoff=0,
                                  limits = NULL,
                                  plot_id = c(1,2,3),
                                  species = NULL,
                                  smooth = FALSE,
                                  knn = 30,
                                  round = 1,
                                  use_color=rev(c("#ff6b35","#f7c59f","#efefd0","#004e89","#1a659e"))){

    plot_df <- meta.data
    plot_df$var <- plot_df[, col_name]

    plot_config <- list(
        `1` = list(
            marm = list(
                list(sample = 'marm_T454', hemisphere = NULL, size = 3.5)
            ),
            mous = list(
                list(sample = 'mous_T313', hemisphere = NULL, size = 3.5)
            ),
            bird = list(
                list(sample = 'bird_a1', hemisphere = NULL, size = 3.5),
                list(sample = 'bird_d4', hemisphere = NULL, size = 3.5),
                list(sample = 'bird_a5', hemisphere = NULL, size = 3.5)
            )
        ),
        `2` = list(
            turt = list(
                list(sample = 'turt_a2', hemisphere = 'left',  size = 4),
                list(sample = 'turt_a2', hemisphere = 'right', size = 4),
                list(sample = 'turt_f5', hemisphere = 'left',  size = 4),
                list(sample = 'turt_f5', hemisphere = 'right', size = 4)
            ),
            axol = list(
                list(sample = 'axol_all', hemisphere = 'left',  size = 6),
                list(sample = 'axol_all', hemisphere = 'right', size = 4)
            )
        ),
        `3` = list(
            lung = list(
                list(sample = 'lung_all', hemisphere = 'left_C',  size = 4),
                list(sample = 'lung_all', hemisphere = 'right_C', size = 4)
            ),
            fish = list(
                list(sample = 'fish_D2', hemisphere = NULL, size = 10),
                list(sample = 'fish_D4', hemisphere = NULL, size = 10)
            ),
            lamp = list(
                list(sample = 'lamp_4',  hemisphere = NULL, size = 6),
                list(sample = 'lamp_5',  hemisphere = NULL, size = 6),
                list(sample = 'lamp_34', hemisphere = NULL, size = 6)
            )
        )
    )

    all_species <- c("marm", "mous", "bird", "turt", "axol", "lung", "fish", "lamp")

    if (is.null(species) || "all" %in% species) {
        species <- all_species
    }

    species <- intersect(all_species, species)
    if (length(species) == 0) {
        stop("species must be one or more of: marm, mous, bird, turt, axol, lung, fish, lamp")
    }

    plot_list <- list()

    for (pid in plot_id) {
        pid_chr <- as.character(pid)
        cfg_sub <- plot_config[[pid_chr]]
        species_in_panel <- intersect(names(cfg_sub), species)

        if (length(species_in_panel) == 0) {
            next
        }

        p_list <- list()

        for (sp in species_in_panel) {
            for (cfg in cfg_sub[[sp]]) {
                df_sub <- plot_df[plot_df$sample == cfg$sample, ]
                if (!is.null(cfg$hemisphere)) {
                    df_sub <- df_sub[df_sub$hemisphere == cfg$hemisphere, ]
                }

                p_list[[length(p_list) + 1]] <- plot_spatial_col(
                    df_sub,
                    col_name = 'var',
                    color = use_color,
                    smooth = smooth,
                    knn = knn,
                    round = round,
                    max.cutoff = max.cutoff,
                    min.cutoff = min.cutoff,
                    size = cfg$size,
                    limits = limits
                )
            }
        }

        if (length(p_list) > 0) {
            if (pid == 3) {
                widths_full <- c(1,1,0.7,0.7,0.7,0.7,0.7)
                plot_list[[pid_chr]] <- Reduce(`+`, p_list) +
                    plot_layout(ncol = length(p_list), widths = widths_full[seq_along(p_list)])
            } else {
                plot_list[[pid_chr]] <- Reduce(`+`, p_list) +
                    plot_layout(ncol = length(p_list))
            }
        }
    }

    if (length(plot_list) == 0) {
        stop("No plots available for the selected species and plot_id.")
    }

    options(repr.plot.width = 22, repr.plot.height = 10)
    return(Reduce(`/`, plot_list))
}

Enrichment_myself <- function(input_gene, database_term_list, database_gene_list, background = NULL) {

  db_df <- data.frame(
    Term = as.character(database_term_list),
    OG_list = as.character(database_gene_list),
    stringsAsFactors = FALSE
  ) %>% distinct(Term, .keep_all = TRUE)
  
  # 将逗号分隔的OG拆分为向量
  db_df$OG_vec <- lapply(db_df$OG_list, function(x) {
    unique(na.omit(strsplit(x, split = ",")[[1]]))
  })
  # 剔除无OG的term
  db_df <- db_df[sapply(db_df$OG_vec, length) > 0, ]
  if(nrow(db_df) == 0) stop("数据库中无有效term-OG对应关系！")
  
  # ===================== 2. 定义背景基因集(OG) =====================
  # 自动获取背景：数据库中所有唯一OG（可手动传入background覆盖）
  if(is.null(background)) {
    all_bg_og <- unique(unlist(db_df$OG_vec))
  } else {
    all_bg_og <- unique(na.omit(as.character(background)))
  }
  bg_total <- length(all_bg_og) # 背景总OG数
  input_total <- length(input_gene) # 输入OG总数
  
  # ===================== 3. 遍历每个term，执行超几何检验 =====================
  enrich_res <- lapply(1:nrow(db_df), function(i) {
    term_now <- db_df$Term[i]
    og_in_term <- db_df$OG_vec[[i]] # 当前term下的所有OG
    term_total <- length(og_in_term) # 当前term的OG总数
    
    # 计算：输入OG中落在当前term的数目（交集）
    hit_og <- intersect(input_gene, og_in_term)
    hit_num <- length(hit_og)
    if(hit_num == 0) return(NULL) # 无交集，跳过
    
    # 超几何分布计算P值（核心：phyper函数，富集分析标准算法）
    # phyper(交集数-1, 通路OG数, 背景-通路OG数, 输入OG数, lower.tail=F)
    p_val <- phyper(
      q = hit_num - 1,
      m = term_total,
      n = bg_total - term_total,
      k = input_total,
      lower.tail = FALSE
    )
    
    # 计算比例（和clusterProfiler一致）
    GeneRatio <- paste0(hit_num, "/", input_total)
    BgRatio <- paste0(term_total, "/", bg_total)
    
    # 返回单行结果
    return(data.frame(
      Term = term_now,
      Count = hit_num,
      GeneRatio = GeneRatio,
      BgRatio = BgRatio,
      pvalue = p_val,
      Hit_OG = paste(hit_og, collapse = ","), # 命中的OG列表
      stringsAsFactors = FALSE
    ))
  })
  
  # ===================== 4. 结果合并与BH校正 =====================
  enrich_res <- do.call(rbind, enrich_res)
  if(is.null(enrich_res) || nrow(enrich_res) == 0) {
    warning("无显著富集结果！")
    return(NULL)
  }
  
  # BH多重检验校正
  enrich_res$p.adjust <- p.adjust(enrich_res$pvalue, method = "BH")
  # 按p值升序排序
  enrich_res <- enrich_res %>% arrange(pvalue) %>% dplyr::relocate(Term, Count, GeneRatio, BgRatio, pvalue, p.adjust, Hit_OG)
  
  return(enrich_res)
}

