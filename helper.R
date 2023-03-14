# tidyverse libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(purrr)
library(tibble)
library(readr)

# plotting libraries
library(cowplot)
library(ggpubr)
library(viridis)
library(gridExtra)
library(grid)
library(ggridges)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(ggforce)
library(scattermore)
library(ggrastr)

# analysis libraries
library(Seurat)
library(SeuratObject)

# helper functions -----

set_idents <- function(so, ident) {
  Idents(so) <- ident
  return(so)
}

set_assay <- function(so, assay) {
  DefaultAssay(so) <- assay
  return(so)
}

# named_group_split
# credit goes to romainfrancois on https://github.com/tidyverse/dplyr/issues/4223
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))

  grouped %>%
    group_split() %>%
    rlang::set_names(names)
}

plot_grid_shared_legend <- function(plotlist,
                                    ...,
                                    legend.position = "bottom",
                                    rel_heights = c(4, 1),
                                    rel_widths = c(4, 1),
                                    align = "none",
                                    axis = "none") {
  legend <-
    cowplot::get_legend(plotlist[[1]] + theme(legend.position = legend.position))
  plotlist <-
    purrr::map(plotlist, ~ ..1+theme(legend.position = "none"))
  plots <- cowplot::plot_grid(plotlist = plotlist, ...)
  if (legend.position == "bottom") {
    cowplot::plot_grid(
      plots,
      legend,
      ncol = 1,
      rel_heights = rel_heights,
      align = align,
      axis = axis
    )
  } else if (legend.position == "right") {
    cowplot::plot_grid(
      plots,
      legend,
      ncol = 2,
      rel_widths = rel_widths,
      align = align,
      axis = axis
    )
  } else {
    stop("Invalid input to legend.position. ")
  }
}

plot_gsea_curve <- function(leg_df,
                            group_var,
                            pt_size = 2,
                            pt_alpha = 0.8,
                            line_size = 1) {
  group_var <- enquo(group_var)
  curve_plot <- leg_df %>%
    filter(in_gs) %>%
    ggplot() +
    aes(rank, running_es) +
    geom_point(aes(color = !!group_var),
               size = pt_size, alpha = pt_alpha) +
    geom_hline(yintercept = 0, size = line_size) +
    xlim(min(leg_df$rank), max(leg_df$rank)) +
    ylab("Running ES") +
    scale_x_continuous(expand = c(0, 0)) +
    theme(
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "top"
    )

  return(curve_plot)

}

plot_gsea_ticks <- function(leg_df,
                            group_var,
                            pt_size = 3,
                            pt_stroke = 1) {
  group_var <- enquo(group_var)

  tick_plot <- leg_df %>%
    filter(in_gs) %>%
    ggplot() +
    aes(x = rank, y = !!group_var) +
    geom_point(
      shape = '|',
      alpha = 0.5,
      size = pt_size,
      stroke = pt_stroke
    ) +
    xlim(min(leg_df$rank), max(leg_df$rank)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )

  n_groups <-
    leg_df[[rlang::as_name(group_var)]] %>% unique() %>% length()
  if (n_groups == 1) {
    tick_plot <- tick_plot +
      theme(axis.title.y = element_blank())
  } else {
    tick_plot <- tick_plot +
      theme(axis.title.y = element_text(
        hjust = 1,
        vjust = 0.5,
        angle = 0
      ))
  }

  return(tick_plot)

}

plot_gsea_curve_pretty <- function(leg_df,
                                   group_var,
                                   palette = NA,
                                   rel_heights = c(4, 1),
                                   show.legend = TRUE,
                                   pt_size = 3,
                                   text_size = 1,
                                   line_size = 1,
                                   tick_size = 3,
                                   tick_stroke = 0.1,
                                   annotate_df = NA,
                                   annotate_coords = c(0.7, 0.8),
                                   annotate_df_FDR_col = q.val,
                                   annotate_df_sscore_col = sscore,
                                   annotate_df_filters = NULL) {
  # e.g. list(Cluster = "prog.")
  group_var <- enquo(group_var)
  annotate_df_FDR_col <- enquo(annotate_df_FDR_col)
  annotate_df_sscore_col <- enquo(annotate_df_sscore_col)

  p_leg <- plot_gsea_curve(
    leg_df,
    group_var = !!group_var,
    pt_size = pt_size,
    line_size = line_size
  ) + theme(plot.margin = margin(r = 5, l = 5))
  p_leg_ticks <- plot_gsea_ticks(
    leg_df,
    group_var = !!group_var,
    pt_size = tick_size,
    pt_stroke = tick_stroke
  ) + theme(plot.margin = margin(r = 5, l = 5))

  if (!is.na(palette_cluster)) {
    p_leg <- p_leg +
      scale_color_manual(values = palette_cluster)
  }

  if (!show.legend) {
    p_leg <- p_leg +
      theme(legend.position = "none")
  }

  if (!is.na(annotate_df)) {
    if (!is.null(annotate_df_filters)) {
      for (name in names(annotate_df_filters)) {
        annotate_df <- annotate_df %>%
          filter(!!sym(name) == annotate_df_filters[[name]])
      }
    }

    p_leg <- p_leg +
      geom_text(
        aes(
          x = max(leg_df$rank)*annotate_coords[[1]],
          y = min(leg_df$running_es) + annotate_coords[[2]]*diff(range(leg_df$running_es)),
          label = paste0(
            "FDR=",
            format(q.val, digits = 2, scientific = TRUE),
            "\n",
            "sscore=",
            format(sscore, digits = 2)
          )
        ),
        data = annotate_df,
        size = text_size,
        hjust = 0
      )
  }

  plot_grid(
    p_leg,
    p_leg_ticks,
    ncol = 1,
    align = "v",
    axis = "tblr",
    rel_heights = rel_heights
  )
}

plot_gsea_curve_grid <- function(leg_list,
                                 gsea_df,
                                 rel_plot_size = 2,
                                 annotate_coords = c(0.05, 0.8),
                                 text_size = TEXT_SIZE * GGPLOT_TEXT_SCALE_FACTOR) {
  leg_list %>%
    imap( ~ {
      contrast_str <- str_extract(..2, ".* / ") %>% str_sub(end = -4)
      cluster_str <-
        str_extract(..2, " / .*") %>% str_sub(start = 4)
      plot_gsea_curve_pretty(
        ..1,
        Clusters,
        palette_cluster,
        pt_size = 0.1,
        line_size = LINE_WIDTH,
        tick_stroke = 0.1,
        rel_heights = c(rel_plot_size, 1),
        show.legend = FALSE,
        annotate_df = gsea_df,
        annotate_coords = annotate_coords,
        annotate_df_filters = list(Clusters = cluster_str,
                                   contrast = contrast_str),
        text_size = text_size
      )
    }) %>%
    plot_grid(plotlist = ., ncol = 4)
}


get_leading_edge_df <- function(rank_vec, gs_vec) {
  tibble(name = names(rank_vec), value = rank_vec) %>%
    arrange(desc(value)) %>%
    mutate(in_gs = name %in% gs_vec) %>%
    mutate(rank = row_number()) %>%
    mutate(score = ifelse(in_gs, 1 / sum(in_gs),-1 / (max(rank) - sum(in_gs)))) %>%
    mutate(running_es = cumsum(score)) %>%
    mutate(
      max_es = max(running_es),
      min_es = min(running_es),
      max_abs_es = max(abs(running_es))
    ) %>%
    mutate(direction = ifelse(max_abs_es == max_es, "up", "dn")) %>%
    mutate(
      inflection_rank = case_when(
        direction == "up" & running_es == max_es ~ rank,
        direction == "dn" & running_es == min_es ~ rank,
        TRUE ~ NA_integer_
      )
    ) %>%
    tidyr::fill(inflection_rank, .direction = "downup") %>%
    mutate(
      is_leading_edge = case_when(
        direction == "up" ~ rank <= inflection_rank,
        direction == "dn" ~ rank >= inflection_rank
      )
    ) %>%
    select(name,
           value,
           in_gs,
           rank,
           running_es,
           direction,
           is_leading_edge)
}

make_galaxy_plot <- function(metadata_df,
                             sample_n = 800,
                             umap1 = "UMAP_1", umap2 = "UMAP_2") {
  if (sample_n != 0)
    data_sampled <- slice_sample(metadata_df, n = sample_n)
  retplot <- ggplot(metadata_df) +
    aes_string(umap1, umap2) +
    stat_density_2d(aes(fill = ..density..), geom = 'raster', contour = FALSE) +
    scale_fill_viridis(option = "magma") +
    coord_cartesian(expand = FALSE, xlim = c(min(metadata_df[[umap1]]), max(metadata_df[[umap1]])),
                    ylim = c(min(metadata_df[[umap2]]), max(metadata_df[[umap2]])))
  if (sample_n != 0)
    retplot <- retplot + geom_point(shape = '.', col = 'white', data = data_sampled)
  return(retplot)
}


calculate_proportion <- function(.data, var, norm_var) {
  var <- enquo(var)
  norm_var <- enquo(norm_var)
  .data %>%
    group_by(!!var, !!norm_var) %>%
    summarise(count = n()) %>%
    group_by(!!norm_var) %>%
    mutate(total = sum(count), prop = count / total) %>%
    group_by(!!var) %>%
    mutate(total_prop = sum(prop), norm_prop = prop / total_prop)
}


# pseudobulk functions ------


make_pseudobulk <-
  function (so, split_var = NULL, sample_var, vars) {
    sce <- seurat_to_sce(so, split_var, sample_var, vars)
    if (is.null(split_var)) {
      groups <- SummarizedExperiment::colData(sce)[, c(sample_var,
                                                       vars)]
      pb <-
        Matrix.utils::aggregate.Matrix(
          Matrix::t(SingleCellExperiment::counts(sce)),
          groupings = groups,
          fun = "sum"
        )
      return(pb)
    }
    else {
      groups <- map(sce, ~ SummarizedExperiment::colData(..1)[,
                                                              c(sample_var, vars)])
      pb_list <-
        map2(
          sce,
          groups,
          ~ Matrix.utils::aggregate.Matrix(
            Matrix::t(SingleCellExperiment::counts(..1)),
            groupings = ..2,
            fun = "sum"
          )
        ) %>% set_names(names(sce))
      return(pb_list)
    }
  }

pb_to_deseq <-
  function (pb,
            sample_var,
            vars,
            formula_str = paste0("~", paste0(vars,
                                             collapse = "+"))) {
    counts_mtx <- as.matrix(pb)
    counts_metadata <- tibble::tibble(id = rownames(pb)) %>%
      tidyr::separate(id, into = c(sample_var, vars), sep = "_")
    de <-
      DESeqDataSetFromMatrix(Matrix::t(counts_mtx),
                             colData = counts_metadata,
                             design = as.formula(formula_str))
    de <- DESeq2::DESeq(de)
  }

# SCP export functions -------

scp_export_clustering_file <-
  function (so,
            str_subset_regex_group = "RNA_snn_res.",
            str_subset_regex_numeric = "_RNA",
            reduction_str = "umap",
            output_dir,
            filename = "clustering_scp") {
    umap_coords <- so@reductions[[reduction_str]]@cell.embeddings %>%
      as_tibble(rownames = "cell")
    metadata_df <- so@meta.data %>% as_tibble(rownames = "cell") %>%
      left_join(umap_coords, by = "cell")
    other_cols_to_select_group <- str_subset(colnames(metadata_df),
                                             str_subset_regex_group)
    other_cols_to_select_num <- str_subset(colnames(metadata_df),
                                           str_subset_regex_numeric)
    clustering_cols <-
      c("NAME",
        "X",
        "Y",
        other_cols_to_select_group,
        other_cols_to_select_num)
    clustering_types_vec <-
      c("TYPE",
        "numeric",
        "numeric",
        rep("group",
            (length(
              other_cols_to_select_group
            ))),
        rep("numeric",
            (length(
              other_cols_to_select_num
            )))) %>% set_names(clustering_cols)
    clustering_types_line <-
      paste0(clustering_types_vec, collapse = "\t")
    clustering_scp <-
      metadata_df %>% mutate(NAME = cell, X = UMAP_1,
                             Y = UMAP_2) %>% select(names(clustering_types_vec))
    write_tsv(clustering_scp, file.path(output_dir, "tmp.txt"))
    f_clustering_raw <- file(file.path(output_dir, "tmp.txt"),
                             open = "r")
    lines <- readLines(f_clustering_raw)
    close(f_clustering_raw)
    new_lines <-
      c(lines[1], clustering_types_line, lines[2:length(lines)])
    f_clustering <- file(file.path(output_dir, paste0(filename,
                                                      ".txt")), open = "w")
    writeLines(new_lines, f_clustering)
    close(f_clustering)
    file.remove(file.path(output_dir, "tmp.txt"))
    print("done!")
  }


scp_export_metadata_file <- function(so,
                                     output_dir,
                                     filename = "scp_metadata",
                                     biosample_id = orig.ident,
                                     donor_id = orig.ident,
                                     species,
                                     species__ontology_label,
                                     disease,
                                     disease__ontology_label,
                                     organ,
                                     organ__ontology_label,
                                     library_preparation_protocol,
                                     library_preparation_protocol__ontology_label,
                                     sex = "unknown",
                                     str_subset_regex_group = NULL,
                                     str_subset_regex_numeric = NULL) {
  biosample_id <- enquo(biosample_id)
  donor_id <- enquo(donor_id)

  metadata_df <- so@meta.data %>%
    as_tibble(rownames = "NAME")

  required_metadata_df <- metadata_df %>%
    dplyr::rename(biosample_id = !!biosample_id) %>%
    mutate(donor_id = !!donor_id) %>%
    mutate(species = species, species__ontology_label = species__ontology_label) %>%
    mutate(disease = disease, disease__ontology_label = disease__ontology_label) %>%
    mutate(organ = organ, organ__ontology_label = organ__ontology_label) %>%
    mutate(
      library_preparation_protocol = library_preparation_protocol,
      library_preparation_protocol__ontology_label = library_preparation_protocol__ontology_label
    ) %>%
    mutate(sex = sex) %>%
    select(
      NAME,
      biosample_id,
      donor_id,
      species,
      species__ontology_label,
      disease,
      disease__ontology_label,
      organ,
      organ__ontology_label,
      library_preparation_protocol,
      library_preparation_protocol__ontology_label,
      sex
    )

  type_vec <- rep("group", ncol(required_metadata_df) - 1)

  if (length(str_subset_regex_group) != 0) {
    message("Looking for matching group columns...")
    group_cols <-
      str_subset(colnames(metadata_df), str_subset_regex_group)
    group_metadata_df <- metadata_df %>%
      dplyr::select(NAME, all_of(group_cols))
    message(paste0("Found column: ", group_cols, "\n"))
    type_vec <-
      c(rep("group", ncol(group_metadata_df) - 1), type_vec)
  }
  if (length(str_subset_regex_numeric) != 0) {
    message("Looking for matching numeric columns...")
    numeric_cols <-
      str_subset(colnames(metadata_df), str_subset_regex_numeric)
    numeric_metadata_df <- metadata_df %>%
      dplyr::select(NAME, all_of(numeric_cols))
    message(paste0("Found column: ", numeric_cols, "\n"))
    type_vec <-
      c(rep("numeric", ncol(numeric_metadata_df) - 1), type_vec)
  }

  type_vec <- c("TYPE", type_vec)
  types_line <- paste0(type_vec, collapse = "\t")

  scp_metadata_df <- numeric_metadata_df %>%
    left_join(group_metadata_df, by = "NAME") %>%
    left_join(required_metadata_df, by = "NAME")
  colnames(scp_metadata_df) <-
    str_replace_all(colnames(scp_metadata_df), "\\.", "\\_")

  write_tsv(scp_metadata_df, file.path(output_dir, "tmp.txt"))
  f_raw <- file(file.path(output_dir, "tmp.txt"), open = "r")
  lines <- readLines(f_raw)
  close(f_raw)
  new_lines <- c(lines[1], types_line, lines[2:length(lines)])
  f_new <-
    file(file.path(output_dir, paste0(filename, ".txt")), open = "w")
  writeLines(new_lines, f_new)
  close(f_new)
  file.remove(file.path(output_dir, "tmp.txt"))
  message("done!")

}

# Seurat helper functions -------


#' @name get_metadata_from_so
#'
#' @param so Seurat object
#' @param assay
#' @param slot one of "counts", "data", or "scale.data"
#' @param genes a character vector of genes or NULL
#' @param metadata a character vector of columns from metadata, "all", or NULL
#' @param reduction a string specifying the reduction.
#'
#' @return data table

get_metadata_from_so <- function(so,
                                 assay = "RNA",
                                 genes = NULL,
                                 slot = "data",
                                 metadata = "all",
                                 reduction = "umap") {
  if (is.null(metadata)) {
    out_df <- so@meta.data %>%
      as_tibble(rownames = "cell") %>%
      dplyr::select(cell)
  } else if (length(metadata) > 1) {
    out_df <- so@meta.data %>%
      as_tibble(rownames = "cell") %>%
      dplyr::select(cell, all_of(metadata))
  } else if (metadata == "all") {
    out_df <- so@meta.data %>%
      as_tibble(rownames = "cell")
  } else {
    out_df <- so@meta.data %>%
      as_tibble(rownames = "cell") %>%
      dplyr::select(cell, all_of(metadata))
  }

  if (!is.null(genes)) {
    # get index
    gene_i <- match(genes, rownames(so))

    # check that genes exists
    if (any(is.na(gene_i))) {
      stop(paste0(genes[which(is.na(gene_i))], " not found in data."))
    }

    # get gene expression
    if (length(genes) == 1) {
      if (slot == "data") {
        gene_df <- so@assays[[assay]]@data[genes, ] %>%
          as_tibble(rownames = "cell")
      } else if (slot == "scale.data") {
        gene_df <- so@assays[[assay]]@scale.data[genes, ] %>%
          as_tibble(rownames = "cell")
      } else {
        stop("Oops, haven't implemented that slot yet. Try again. ")
      }
      gene_df[[genes]] <- gene_df$value
      gene_df <- gene_df %>%
        dplyr::select(-value)
    } else {
      if (slot == "data") {
        gene_df <- so@assays[[assay]]@data[genes, ] %>%
          as.matrix() %>% t() %>%
          as_tibble(rownames = "cell")
      } else if (slot == "scale.data") {
        gene_df <- so@assays[[assay]]@scale.data[genes, ] %>%
          as.matrix() %>% t() %>%
          as_tibble(rownames = "cell")
      } else {
        stop("Oops, haven't implemented that slot yet. Try again. ")
      }
    }

    out_df <- out_df %>% left_join(gene_df, by = "cell")
  }

  if (!is.null(reduction)) {
    if (length(reduction) > 1) {
      stop(paste0("reduction must be of length 1"))
    } else if (!(reduction %in% names(so@reductions))) {
      stop(paste0(reduction, " not found in reductions. "))
    }

    reduc_df <- so@reductions[["umap"]]@cell.embeddings %>%
      as_tibble(rownames = "cell")
    out_df <- out_df %>% left_join(reduc_df, by = "cell")
  }

  return(out_df)
}



#' @name plot_violin
#'
#' @param .tbl a data frame or tibble with x and y coordinates
#' @param var variable to plot as symbol
#' @param group group to plot as symbol
#'
#' @return ggplot object

plot_violin <- function(.tbl, group, var, pt.size = 0, pt.stroke = 1) {
  var <- enquo(var)
  group <- enquo(group)
  p <- .tbl %>%
    ggplot() +
    aes(!!group, !!var, fill = !!group) +
    geom_violin(scale = "width") +
    scale_y_continuous(expand = c(0, 0)) +
    remove_x_spine() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  if (pt.size != 0) {
    p <- p +
      geom_jitter(width = 0.25, size = pt.size, stroke = pt.stroke)
  }
  return(p)
}


#' @name remove_x_spine
#'
#' @return a ggplot theme to remove the x spines

remove_x_spine <- function(...) {
  remove.x.spine.theme <- theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    validate = TRUE,
    ...
  )
  return(remove.x.spine.theme)
}

#' @name remove_spines
#'
#' @return a ggplot theme to remove the spines

remove_spines <- function(...) {
  remove.spine.theme <- theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    validate = TRUE,
    ...
  )
  return(remove.spine.theme)
}

