prepare_violin_data_colors <-
  function(my_object,
           genes_to_investigate,
           colors) {
    n_genes <- length(genes_to_investigate)
    
    n_clusters <-
      my_object$final_cluster_labels %>%
      unique() %>%
      length()
    
    ordered_colors <- rev(colors)
    
    new_levels <-
      Seurat::Idents(object = my_object) %>%
      levels() %>%
      rev()
    
    ordered_colors <-
      purrr::set_names(
        ordered_colors,
        as.character(new_levels)
      )
    
    Seurat::Idents(object = my_object) <-
      Seurat::Idents(object = my_object) %>%
      factor(levels = new_levels)
    
    cell_clusters <-
      tibble::tibble(
        cluster_id = Seurat::Idents(my_object),
        cell_id = names(Seurat::Idents(my_object))
      )
    
    fetch_data <-
      Seurat::FetchData(
        object = my_object,
        vars = genes_to_investigate,
        slot = "data"
      ) %>%
      tibble::rownames_to_column("cell_id")
    
    # Extracted the noise code from `Seurat:::SingleExIPlot`
    # (cf. lines 18-31). By setting the random number generator seed,
    # Seurat adds a tiny bit of noise to the expression values in the
    # plot. The same noise vector gets added to each gene. I've adapted
    # the strategy below to achieve the same plots as Seurat. The noise
    # changes how the violins get rendered.
    seed.use <- 42
    
    set.seed(seed = seed.use)
    
    noise <- rnorm(n = nrow(fetch_data)) / 1e+05
    
    seurat_data <-
      fetch_data %>%
      tidyr::gather(gene_id, expression_value, -cell_id) %>%
      # Maintain the ordering of the genes as requested in the vector.
      # Otherwise, the genes will be ordered alphabetically.
      dplyr::mutate(gene_id = factor(gene_id, levels = genes_to_investigate)) %>%
      dplyr::inner_join(cell_clusters) %>%
      # Add the noise vector generated above.
      dplyr::mutate(expression_value_noise = expression_value + rep(noise, n_genes))
    
    list(
      seurat_data = seurat_data,
      ordered_colors = ordered_colors
    )
  }


construct_violin_plot <- function(my_object,
                                  genes_to_investigate,
                                  colors,
                                  use_noise = TRUE,
                                  scale = "free_x") {
  
  object_data <-
    prepare_violin_data_colors(
      my_object = my_object,
      genes_to_investigate = genes_to_investigate,
      colors = colors
    )
  
  p_cluster_expression <-
    object_data$seurat_data %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = cluster_id,
        y = if(use_noise) expression_value_noise else expression_value
      )
    ) +
    ggplot2::geom_violin(
      ggplot2::aes(fill = cluster_id),
      scale = "width",
      adjust = 1
    ) +
    ggplot2::facet_grid(~ gene_id, scale = scale) +
    ggplot2::scale_fill_manual(
      values = object_data$ordered_colors
    ) +
    ggplot2::coord_flip()
  
  p_cluster_expression
}


theme_violins <- function() {
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = rel(1.5)),
    legend.position = "none"
  )
}

compare_regional_expression <- function(dge_list,
                                        celltype_vector,
                                        genes_of_interest){
  
  celltype_dge_list <- vector(length = length(dge_list), mode = "list")
  for(i in 1:length(dge_list)){
    index_of_celltype_dge <- which(names(dge_list[[i]]) == celltype_vector[i])
    df <- dge_list[[i]][[index_of_celltype_dge]]
    df <- df %>%
      filter(gene %in% genes_of_interest) %>%
      select("gene", "p_val", "avg_logFC",
             "pct.1", "pct.2", "p_val_adj", 
             "delta.pct", "exp.group1", "exp.group2",
             "logFC.pseudobulk") %>%
      mutate(dataset = names(dge_list)[i],
             celltype = celltype_vector[i]) 
    
    celltype_dge_list[[i]] <- df
  }
  
  return_tidy_data <- do.call(rbind, celltype_dge_list)
  return(return_tidy_data)
  
}



regional_figure <- function(object, 
                            cell_population,
                            celltype_vector, 
                            gene_of_interest,
                            logFC_variable,
                            dge_list,
                            dataset_factor_levels,
                            dataset_colors,
                            upper_ylim){
  
  Idents(object) <- "celltype"
  my_subset <- subset(object, idents = cell_population)
  my_subset[["region"]] <- factor(my_subset[["region"]][,"region"], 
                                  levels = c("perifovea", "fovea"))
  Idents(my_subset) <- "region"
  
  vln_plot <- VlnPlot(my_subset, 
                      features = gene_of_interest, 
                      pt.size = 0, 
                      cols = c("grey10", "grey60")) + 
    coord_flip() + 
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank())
  
  tidy_data <- compare_regional_expression(dge_list = dge_list, 
                                           celltype_vector = celltype_vector, 
                                           genes_of_interest = gene_of_interest)
  tidy_data <- tidy_data %>% 
    mutate(dataset = factor(dataset, levels = dataset_factor_levels))
  
  large_logFCs <- tidy_data %>% 
    filter(abs(!!rlang::sym(logFC_variable)) > upper_ylim) %>%
    mutate(saved_logFC = !!rlang::sym(logFC_variable)) %>%
    mutate(plot_logFC = ifelse(!!rlang::sym(logFC_variable) < 0, -1*upper_ylim, upper_ylim))
  
  
  my_ggplot <- ggplot(tidy_data, aes(x = dataset, y = !!rlang::sym(logFC_variable))) +
    geom_bar(stat = "identity", 
             aes(fill = dataset)) +
    scale_fill_manual(values = dataset_colors, 
                      breaks = dataset_factor_levels,
                      drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    ylim(c(upper_ylim*-1, upper_ylim)) +
    cowplot::theme_half_open() +
    theme(#axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none") + 
    geom_hline(yintercept = 0, size = 1.75)
  
  my_ggplot <- my_ggplot +
    geom_bar(data = large_logFCs, 
             mapping = aes(x = dataset, y = plot_logFC), 
             stat = "identity",
             fill = "black") +
    geom_text(data = large_logFCs, aes(x = dataset, y = -1*(plot_logFC), label = round(saved_logFC, 1)), angle = 30, hjust = 0)
  
  my_ggplot <- my_ggplot + 
    theme(axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank())
  
  final <- (vln_plot | my_ggplot)
  plot(final)
  return(final)
  
}