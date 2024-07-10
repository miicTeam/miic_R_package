#*******************************************************************************
# Filename   : interactmiic.utils.R                Creation date: 13 june 2024
#
# Description: Utility functions for interactMIIC
#
# Author     : Louise DUPUIS
#*******************************************************************************

#-------------------------------------------------------------------------------
# interactmiic.mosaic
#-------------------------------------------------------------------------------
#' Create Mosaic datatable for interactMIIC
#'
#' @description Function that takes an seurat_object raw counts,
#'  senders and receivers information metadata and create a mosaic
#'  datable that contains first senders network
#'  variables then receivers network variables
#'
#'
#' @param seurat_object [a Seurat object]
#' A Single-Cell transcriptomics Seurat object
#'
#' @param assay_name [a string]. Gives the name of the assay to take
#' the transcriptomics raw counts from (usually 'RNA')
#'
#' @param interact_ident [a string]. Gives the name of the metadata
#' containing the celltypes population
#'
#' @param senders_name [a string]. Gives the name of the senders
#'
#' @param receivers_name [a string]. Gives the name of the receivers
#'
#'
#' @param genes_senders [a vector] A list of selected genes (strings) for the
#' senders cells. Overlap with genes_receivers is not an issue
#'
#' @param genes_receivers [a vector] A list of selected genes (strings) for the
#' receivers cells. Overlap with genes_senders is not an issue
#'
#' @param metadata_senders [a vector] A list of metadata (strings) that will
#' be added to the genes network. They can be relevant to senders
#' and/or receivers cells.
#'
#' @param metadata_receivers [a vector] A list of metadata (strings) that will
#' be added to the genes network. They can be relevant to senders
#' and/or receivers cells.
#'
#' @return A mosaic datatable that can be used for interactMIIC
#' network reconstruction.
#'
#' @export
#'
#-------------------------------------------------------------------------------

interactmiic.mosaic <- function(seurat_object,
                                  assay_name = "RNA",
                                  interact_ident,
                                  senders_name,
                                  receivers_name,
                                  genes_senders,
                                  genes_receivers,
                                  metadata_senders = NULL,
                                  metadata_receivers = NULL) {

  # Validate inputs
  if (!inherits(seurat_object, "Seurat")) {
    stop("The 'seurat_object' must be a Seurat object.")
  }
  if (!is.character(assay_name) || length(assay_name) != 1) {
    stop("The 'assay_name' must be a single character string.")
  }
  if (!(assay_name %in% names(seurat_object@assays))) {
    stop(paste("The assay", assay_name, "is not present in the seurat_object object."))
  }
  if (!is.character(interact_ident) || length(interact_ident) != 1) {
    stop("The 'interact_ident' must be a single character string.")
  }
  if (!interact_ident %in% colnames(seurat_object@meta.data)) {
    stop(paste("Metadata", interact_ident, "is not found in the Seurat object"))
  }
  if (!is.character(senders_name) || length(senders_name) != 1) {
    stop("The 'senders_name' must be a single character string.")
  }
  if (!is.character(receivers_name) || length(receivers_name) != 1) {
    stop("The 'receivers_name' must be a single character string.")
  }
  if (!is.vector(genes_senders) || !is.character(genes_senders)) {
    stop("The 'genes_senders' must be a character vector.")
  }
  if (!is.vector(genes_receivers) || !is.character(genes_receivers)) {
    stop("The 'genes_receivers' must be a character vector.")
  }
  if (!is.null(metadata_senders) && 
        (!is.vector(metadata_senders) || !is.character(metadata_senders))) {
    stop("The 'metadata_senders' must be a character vector or NULL.")
  }

  if (!is.null(metadata_receivers) && 
        (!is.vector(metadata_receivers) || !is.character(metadata_receivers))) {
    stop("The 'metadata_receivers' must be a character vector or NULL.")
  }

  #Start preparation
  Seurat::Idents(seurat_object) <- interact_ident


  # We need to detect genes that are present in both senders and
  # receivers variables lists to duplicate them
  # (to have one node in each network)
  duplicated_genes <- intersect(genes_senders, genes_receivers)
  duplicated_meta <- intersect(metadata_senders, metadata_receivers)

  # SENDERS preparation
  print(paste("Preparing", senders_name, "cells table ..."))
  so_senders <- subset(seurat_object, idents = senders_name)
  sub_matrix_senders <- as.data.table(t(as.matrix(Seurat::GetAssayData(so_senders, assay = assay_name, slot = "counts"))))
  shuff_senders <- sample(nrow(sub_matrix_senders))
  sub_matrix_senders <- sub_matrix_senders[shuff_senders,]
  tmp_senders <- sub_matrix_senders[, ..genes_senders]

  if (!is.null(metadata_senders)) {
    for (onemeta in metadata_senders) {
      if (!onemeta %in% colnames(so_senders@meta.data)) {
        stop(paste("Metadata", onemeta, "is not found in the Seurat object"))
      }
      labels <- so_senders@meta.data[[onemeta]][shuff_senders]
      if (all(is.na(labels))) {
        warning(paste("in the senders cells", senders_name, ", metadata",
        onemeta, "contains only NA values and will be excluded from the senders network"))
        next  # Skip to the next iteration of the loop
      }
      tmp_senders[, onemeta] <- labels
    }
  }

  rm(so_senders)
  rm(sub_matrix_senders)
  gc()

  # RECEIVERS preparation
  print(paste("Preparing", receivers_name, "cells table ..."))
  so_receivers <- subset(seurat_object, idents = receivers_name)
  sub_matrix_receivers <- as.data.table(t(as.matrix(Seurat::GetAssayData(so_receivers, assay = assay_name, slot = "counts"))))
  shuff_receivers <- sample(nrow(sub_matrix_receivers))
  sub_matrix_receivers <- sub_matrix_receivers[shuff_receivers,]
  tmp_receivers <- sub_matrix_receivers[, ..genes_receivers]

  if (!is.null(metadata_receivers)) {
    for (onemeta in metadata_receivers) {
      if (!onemeta %in% colnames(so_receivers@meta.data)) {
        stop(paste("Metadata", onemeta, "not found in receivers' metadata."))
      }
      labels <- so_receivers@meta.data[[onemeta]][shuff_receivers]
      if (all(is.na(labels))) {
        warning(paste("in the receivers cells", receivers_name, ", metadata",
        onemeta, "contains only NA values and will be excluded from the receivers network")) 
        next  # Skip to the next iteration of the loop
      }
      tmp_receivers[, onemeta] <- labels
    }
  }

  rm(so_receivers)
  rm(sub_matrix_receivers)
  gc()


  #MERGE - Create the mosaic dataset
  print("Merge into the mosaic datable for interactMIIC")
  
  if (length(duplicated_genes) >0) {
    # Rename duplicated genes
    data.table::setnames(tmp_receivers, old = duplicated_genes, new = paste0(duplicated_genes, "_receivers"))
    data.table::setnames(tmp_senders, old = duplicated_genes, new = paste0(duplicated_genes, "_senders"), skip_absent = TRUE)
  }
  if (length(duplicated_meta) >0) {
    # Rename duplicated metadata in senders and receivers
      data.table::setnames(tmp_receivers, old = duplicated_meta, new = paste0(duplicated_meta, "_receivers"))
      data.table::setnames(tmp_senders, old = duplicated_meta, new = paste0(duplicated_meta, "_senders"), skip_absent = TRUE)
  }


  new_genes_senders <- colnames(tmp_receivers)
  new_genes_receivers <- colnames(tmp_senders)

  tmp_receivers[, new_genes_receivers] <- NA
  tmp_senders[, new_genes_senders] <- NA

  tmp <- plyr::rbind.fill(tmp_receivers, tmp_senders)

  return(tmp)
}


#-------------------------------------------------------------------------------
# interactmiic.state_order
#-------------------------------------------------------------------------------
#' Create a state order file for interactMIIC
#'
#' @description The state order is an mandatory file that allows you to input
#' optional information for the computation (such as contextual variables,
#' types of variables (otherwise detected automatically)... ) and information
#' for the display (groups of nodes, levels ordering of categorical data ...).
#'
#'
#' @param mosaic_data_table [a data frame]
#' An output of interactmiic.mosaic()
#'
#' @param genes_senders [a vector] A list of selected genes (strings) for the
#' senders cells. Overlap with genes_receivers is not an issue
#'
#' @param genes_receivers [a vector] A list of selected genes (strings) for the
#' receivers cells. Overlap with genes_senders is not an issue
#'
#' @param ligands [a vector] A list of selected CCC genes (strings) for the
#' senders cells. Overlap with genes_senders is not an issue
#'
#' @param receptors [a vector] A list of selected CCC genes (strings) for the
#' receivers cells. Overlap with genes_receivers is not an issue
#'
#' @param metadata_senders [a vector] A named list of metadata (strings) with
#' levels as items. See exemple
#' 
#' @param metadata_receivers [a vector] A named list of metadata (strings) with
#' levels as items. See exemple
#'
#' @examples
#' \donttest{
#' metadata <- list(
#'  treatment = c("Control", "Treated"),
#'  another_meta = c("Level1", "Level2", "Level3")
#')
#' }
#'
#' @return A state order that can be used for interactMIIC
#' network reconstruction.
#'
#' @export
#'
#-------------------------------------------------------------------------------

interactmiic.state_order <- function(mosaic_data_table,
                                     metadata_senders,
                                     metadata_receivers,
                                     genes_senders,
                                     ligands,
                                     receptors,
                                     genes_receivers) {

  # Input checks
  if (!is.data.frame(mosaic_data_table)) {
    stop("mosaic_data_table must be a data frame.")
  }
  if (!is.list(metadata_senders)) {
    stop("metadata_senders is not a named list")
  }
  if (!is.list(metadata_receivers)) {
    stop("metadata_receivers is not a named list")
  }
  if (!is.vector(genes_senders) || !is.character(genes_senders)) {
    stop("genes_senders must be a character vector.")
  }
  if (!is.vector(genes_receivers) || !is.character(genes_receivers)) {
    stop("genes_receivers must be a character vector.")
  }
  if (!is.vector(ligands) || !is.character(ligands)) {
    stop("ligands must be a character vector.")
  }
  if (!is.vector(receptors) || !is.character(receptors)) {
    stop("receptors must be a character vector.")
  }

  # Initialize state_order data frame
  state_order <- data.frame(var_names = character(),
                            var_type = numeric(),
                            levels_increasing_order = character(),
                            is_contextual = numeric(),
                            group = character(),
                            group_color = character())

  # Add metadata information
  duplicated_meta <- intersect(names(metadata_senders), names(metadata_receivers))
  if (length(metadata_senders) >0) {
    for (onemeta in names(metadata_senders)) {
        levels_increasing_order <- paste(metadata_senders[[onemeta]], collapse = ",")
        if (ifelse(onemeta %in% duplicated_meta, paste0(onemeta, "_senders"), onemeta) %in% colnames(mosaic_data_table)) {
          state_order[nrow(state_order) + 1,] <- list(ifelse(onemeta %in% duplicated_meta, paste0(onemeta, "_senders"), onemeta), 0, levels_increasing_order, 0, "metadata", "FFD050")
        }
      }
  }

  if (length(metadata_receivers) >0) {
    for (onemeta in names(metadata_receivers)) {
        levels_increasing_order <- paste(metadata_receivers[[onemeta]], collapse = ",")
        if (ifelse(onemeta %in% duplicated_meta, paste0(onemeta, "_receivers"), onemeta) %in% colnames(mosaic_data_table)) {
          state_order[nrow(state_order) + 1,] <- list(ifelse(onemeta %in% duplicated_meta, paste0(onemeta, "_receivers"), onemeta), 0, levels_increasing_order, 0, "metadata", "FFD050")
        }
      }
  }

  # Identify duplicated genes (special case)
  dups <- c(paste0(intersect(genes_senders, genes_receivers), "_senders"),paste0(intersect(genes_senders, genes_receivers), "_receivers"))

  interact_meta <- as.vector(state_order$var_names)
  
  # Assign the genes
  for(col in setdiff(colnames(mosaic_data_table), interact_meta)) {

    # Initialize variables
    group <- ""
    color <- ""
    
    # Test if gene is duplicated and assign group
    if (col %in% dups) {
      gene <- stringr::str_split(col, "_")[[1]][1]
      family <- stringr::str_split(col, "_")[[1]][2]
    } else {
      gene <- col
      family <- ifelse(gene %in% genes_senders, "senders", "receivers")
    }

    if (family == "senders") {
        group <- "senders genes"
        color <- "EBF5D8"
        if (gene %in% ligands) {
          group <- "ligands"
          color <- "28827A"
        }
      } else {
        group <- "receivers genes"
        color <- "DBF5D8"
        if (gene %in% receptors) {
          group <- "receptors"
        color <- "D37538"
        }
      }

    # Add to state_order
    if (unique(dplyr::n_distinct(mosaic_data_table[,col]) <= 5)) {
      state_order[nrow(state_order) + 1,] <- list(col, 0, "", 0, group, color)
    } else {
      state_order[nrow(state_order) + 1,] <- list(col, 1, "", 0, group, color)
    }
  }
  return(state_order)
}


#-------------------------------------------------------------------------------
# interactmiic.MIselection
#-------------------------------------------------------------------------------
#' Select informative genes for the network reconstruction
#'
#' @description This function offers an unsupervized feature selection tool
#' based on fast pairwise mutual information computation. If the user do not
#' gives genes_senders/genes_receivers, this function will select relevant
#' genes for the interactMIIC network
#'
#' @param seurat_object [a Seurat object]
#' A Single-Cell transcriptomics Seurat object
#'
#' @param assay_name [a string]. Gives the name of the assay to take
#' the transcriptomics raw counts from (usually 'RNA')
#'
#' @param interact_ident [a string]. Gives the name of the metadata
#' containing the celltypes population
#'
#' @param oneinteract [a string]. Either the senders tag or the receivers tag
#'
#' @param goi [a vector] A list of genes of interest in your dataset.
#' Must be of length <15. These genes will be individually used as pivot
#'  variable to look for mutual information in your dataset.
#'
#' @param metadata_list [a vector] A list of metadata names of interest in
#' your dataset. Must be of length <15. These metadata will be individually
#' used as pivot variable to look for mutual information in your dataset.
#'
#' @param n_genes [a numeric] The number of genes to keep after pairwise
#' mutual information ranking. Default is 30
#'
#' @param plot [a boolean] A boolean to specify if you want a heatmap plot
#' of the highest mutual information found with the features of interest,
#' default is FALSE. Will be saved in output_dir
#'
#' @param save [a boolean] A boolean to specify if you want to save the full
#' mutual information table, default is FALSE. Will be saved in output_dir
#'
#' @param output_dir [a string] A string to specify an output directory to
#' export optional outputs to, default is "MI_tables"
#'
#' @return A vector of genes of length n_genes sharing the most mutual
#' information with the given features of interest
#'
#' @export
#'
#-------------------------------------------------------------------------------

interactmiic.MIselection <- function(seurat_object,
                                     assay_name = "RNA",
                                     interact_ident,
                                     oneinteract,
                                     goi,
                                     metadata_list,
                                     n_genes = min(length(goi)*5,100),
                                     plot = FALSE,
                                     save = FALSE,
                                     output_dir = "MI_tables") {

  # Gather features of interest
  foi <- unique(c(goi, metadata_list))

  # Input checks
  if (!inherits(seurat_object, "Seurat")) {
    stop("The 'seurat_object' must be a Seurat object.")
  }
  if (!is.character(assay_name) || length(assay_name) != 1) {
    stop("The 'assay_name' must be a single character string.")
  }
  if (!(assay_name %in% names(seurat_object@assays))) {
    stop(paste("The assay", assay_name, "is not present in the seurat_object object."))
  }
  if (!is.character(interact_ident) || length(interact_ident) != 1) {
    stop("The 'interact_ident' must be a single character string.")
  }
  if (!is.character(foi) || length(foi) > 30) {
    stop("the features of interest goi+metadata must be a character vector with length less than or equal to 30.")
  }
  if (!is.numeric(n_genes) || length(n_genes) != 1 || n_genes <= 0 || n_genes > dim(seurat_object)[1]) { 
    stop("The 'n_genes' must be a positive integer less than or equal to the number of unique features of interest.")
  }
  if (!is.logical(plot) || length(plot) != 1) {
    stop("The 'plot' argument must be a boolean.")
  }
  if (!is.logical(save) || length(save) != 1) {
    stop("The 'save' argument must be a boolean.")
  }
  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("The 'output_dir' must be a single character string.")
  }

  # Initialization
  Seurat::Idents(seurat_object) <- interact_ident
  MI_oneinteract <- data.table(variables = character(), foi = character(), value = numeric())

  # Prepare dataset
  message("Performing feature selection for ", oneinteract, " cells ...")
  tmp <- subset(seurat_object, idents = oneinteract)
  seurat_df <- as.data.frame(t(as.matrix(Seurat::GetAssayData(tmp, assay = assay_name, slot = "counts"))))

  for (onemeta in metadata_list) {
    if (!onemeta %in% colnames(tmp@meta.data)) {
      stop(paste("Metadata", onemeta, "is not found in the Seurat object"))
    }
    if (all(is.na(tmp@meta.data[[onemeta]]))) {
      warning(paste("in the", oneinteract, "cells, metadata",
                    onemeta, "contains only NA values and will be excluded from the features of interest"))
      foi <- foi[foi != onemeta]
      next  # Skip to the next iteration of the loop
    }
    seurat_df[, onemeta] <- as.character(tmp@meta.data[[onemeta]])
  }

  # Calculate pairwise MI for each of the features of interest
  for (onefoi in foi) {
    message("Computing pairwise mutual information with ", onefoi, " in ", oneinteract, " population ...")
    if (!(onefoi %in% colnames(seurat_df))) {
      stop(onefoi, " not found in data.")
    }

    values <- seurat_df[, onefoi]
    MI_table <- data.frame("to_be_renamed" = rep(NA, ncol(seurat_df)), stringsAsFactors = FALSE)
    colnames(MI_table) <- onefoi
    rownames(MI_table) <- colnames(seurat_df)
    pairwise_variables <- rownames(MI_table)[is.na(MI_table[, onefoi]) & (rownames(MI_table) != onefoi)]
    
    cat("Computing MI with ", onefoi, " for ", length(pairwise_variables), " variables")
    mat <- seurat_df[, pairwise_variables]
    col_low <- 1
    loop_range <- 100
    while (col_low <= ncol(mat)) {
      col_max <- min(ncol(mat), col_low + loop_range - 1)
      
      df_loop <- mat[, col_low:col_max]
      df_st <- data.frame(
        "var_names" = c(colnames(df_loop), "Metadata"),
        "var_type" = c(rep(1, ncol(df_loop)), 1),
        "is_contextual" = c(rep(0, ncol(df_loop)), 0),
        "is_consequence" = c(rep(1, ncol(df_loop)), 0),
        stringsAsFactors = FALSE
      )
      df_loop$Metadata <- values
      for (col in colnames(df_loop)) {
        if (unique(n_distinct(df_loop[, col]) <= 5)) {
          df_st[df_st$var_names == col, ]$var_type <- 0
        }
      }
      
      miic_res <- suppressMessages({
        suppressWarnings({
          miic(
            input_data = df_loop,
            state_order = df_st,
            orientation = FALSE,
            latent = "no", 
            n_threads = 8
          )
        })
      })
      
      miic_res <- miic_res$all.edges.summary
      rownames(miic_res) <- NULL
      rownames(miic_res)[miic_res$x != "Metadata"] <- miic_res[miic_res$x != "Metadata", "x"]
      rownames(miic_res)[miic_res$y != "Metadata"] <- miic_res[miic_res$y != "Metadata", "y"]
      if (col_low == 1) {
        df_miic <- miic_res[, "info_shifted", FALSE]
      } else {
        df_miic <- rbind(df_miic, miic_res[, "info_shifted", FALSE])
      }
      
      if ((col_low-1) %% 2500 == 0) {
        cat( col_low, "VARIABLES DONE OUT OF", length(pairwise_variables))
      }
      col_low <- col_low + loop_range
    }
    
    MI_table[rownames(df_miic), onefoi] <- round(df_miic[, "info_shifted"], 6)
    MI_table[is.na(MI_table[, onefoi]), onefoi] <- 0
    MI_table[, onefoi] <- as.data.frame(as.numeric(MI_table[, onefoi]))

    # Add results to master table
    MI_oneinteract <- rbind(MI_oneinteract, data.table(variables = row.names(MI_table), foi = onefoi, value = MI_table[, 1]))
    rm(MI_table)
  }

  rm(seurat_df)
  MI_oneinteract$foi <- factor(MI_oneinteract$foi, levels = unique((MI_oneinteract$foi)[order(MI_oneinteract$value, decreasing = FALSE)]))

  max_values <- MI_oneinteract %>%
    dplyr::group_by(variables) %>%
    dplyr::slice(which.max(value))
  ordered_variables <- as.character(max_values$variables[order(max_values$value, decreasing = TRUE)])

  # Writing table
  if (save) {
    # Create the output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    message("Saving results in ", output_dir,  " ...")
    write.table(MI_oneinteract,
                file = file.path(output_dir, paste0(oneinteract, "_MI_table.csv")),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE)
  }

  # Plotting
  if (plot) {
    message("Generating heatmap of top 150 MI values for each feature of interest ...")

    # Create the output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    plot <- MI_oneinteract %>%
      dplyr::filter(variables %in% ordered_variables[1:150]) %>%
      ggplot2::ggplot(ggplot2::aes(x = reorder(variables, -value),
                                   y = foi,
                                   fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(low = "floralwhite", high = "deeppink4") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                     axis.text.y = ggplot2::element_text(size = 10),
                     axis.title = ggplot2::element_blank())

    plot_filename <- file.path(output_dir,
                               paste0(oneinteract, "_MI_heatmap_filtered.png"))
    message("Saving heatmap in ", output_dir,  " ...")
    ggplot2::ggsave(plot, filename = plot_filename,
                    dpi = 300, units = "px",
                    width = 8000, height = 500 + 500 * length(foi))

    # Remove plot object from memory
    rm(plot)
  }

  return(ordered_variables[1:n_genes])
}




#-------------------------------------------------------------------------------
# interactmiic.layout
#-------------------------------------------------------------------------------
#' Create a layout file for interactMIIC
#'
#' @description This function generates an automatic json file to input into
#' your interactMIIC job so your network benefits from a clear interaction
#' display. This initial layout can be modified with your mouse cursor or
#' with the "Relax" button on the MIIC server for a more organic display.
#'
#'
#' @param state_order A data frame containing the state order information.
#' output of interactmiic.state_order()
#'
#' @param labelSize The size of the labels in the layout. Default is 16.
#'
#' @param threshold The threshold for log confidence on the edges. Default is 0.
#'
#' @param edgeLength The default edge length for the network. Default is 120.
#'
#' @param network_height The height of the  network grids for gene senders and
#' receivers in terms of number of nodes. Default is automatically calculated 
#' from the number of genes in the network.
#'
#' @param CCC_sizefactor The multiplier factor for the ligands and receptors
#' grids height (network_height*CCC_sizefactor). Default is 1.3.
#'
#' @param node_spacing The spacing between nodes in the layout. Default is 100.
#'
#' @return A json file with positions for each node of the network. See how to
#' include it in your interactMIIC job on the MIIC webserver.
#'
#' @export
#'
#-------------------------------------------------------------------------------


interactmiic.layout <- function(state_order,
                                labelSize = 16,
                                threshold = 0,
                                edgeLength = 120,
                                network_height = NULL,
                                CCC_sizefactor = 1.3,
                                node_spacing = 100) {


# Input checks
  if (!is.data.frame(state_order)) {
    stop("state_order must be a data frame.")
  }
  if (!"var_names" %in% names(state_order)) {
    stop("Input data must contain 'var_names' column.")
  }
  if (!"group" %in% names(state_order)) {
    stop("Input data must contain 'group' column.")
  }
  if (!is.numeric(labelSize) || length(labelSize) != 1 || labelSize <= 0 || labelSize != floor(labelSize)) {
    stop("labelSize must be a positive number.")
  }
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold < 0) {
    stop("threshold must be a non-negative number.")
  }
  if (!is.numeric(edgeLength) || length(edgeLength) != 1 || edgeLength <= 0) {
    stop("edgeLength must be a positive number.")
  }
  if (!is.null(network_height) & (!is.numeric(network_height) || length(network_height) != 1 || network_height <= 0 || network_height != floor(network_height))) {
    stop("network_height must be a positive integer.")
  }
  if (!is.numeric(CCC_sizefactor) || length(CCC_sizefactor) != 1 || CCC_sizefactor <= 0) {
    stop("CCC_sizefactor must be a positive number.")
  }
  if (!is.numeric(node_spacing) || length(node_spacing) != 1 || node_spacing <= 0) {
    stop("node_spacing must be a positive number.")
  }

  # Generate empty layout json file
  network_layout <- list(labelSize = c(labelSize),
                         threshold = c(threshold),
                         edgeLength = c(edgeLength),
                         savedPositions = list())

  # Generate empty list of coordinates
  tmp_list <- vector(mode = "list", length = nrow(state_order))
  names(tmp_list) <- state_order$var_names

  # Filter each group
  genes_senders <- filter(state_order, group == "senders genes")
  ligands <- filter(state_order, group == "ligands")
  receptors <- filter(state_order, group == "receptors")
  genes_receivers <- filter(state_order, group == "receivers genes")

  # Define the height of the layout so there is approximatively 2 columns of ligands/receptors 
  network_height <- floor(max(nrow(ligands), nrow(receptors))/(2*CCC_sizefactor))

  # Calculate widths for each group
  LR_height <- floor(CCC_sizefactor * network_height)
  width_senders <- ceiling(nrow(genes_senders) / network_height)
  width_ligands <- ceiling(nrow(ligands) / LR_height)
  width_receptors <- ceiling(nrow(receptors) / LR_height)
  width_receivers <- ceiling(nrow(genes_receivers) / network_height) 


  # Initial x_start positions
  x_start_senders <- -1000
  x_start_ligands <- x_start_senders + width_senders * node_spacing + node_spacing / 2
  x_start_receptors <- x_start_ligands + width_ligands * node_spacing + node_spacing * 3
  x_start_receivers <- x_start_receptors + width_receptors * node_spacing + node_spacing / 2

  # Genes Senders layout (rectangle shape)
  y_start_senders <- - (network_height * node_spacing) / 2
  k <- 1
  for (i in 1:width_senders) {
    for (j in 1:network_height) {
      if (k <= nrow(genes_senders)) {
        tmp_list[[genes_senders$var_names[k]]] <- c(x_start_senders + (i-1) * node_spacing, y_start_senders + (j-1) * node_spacing)
        k <- k + 1
      }
    }
  }
  
  # Ligands layout (line shape)
  y_start_ligands <- - (LR_height * node_spacing) / 2
  k <- 1
  for (i in 1:width_ligands) {
    for (j in 1:LR_height) {
      if (k <= nrow(ligands)) {
        tmp_list[[ligands$var_names[k]]] <- c(x_start_ligands + (i-1) * node_spacing, y_start_ligands + (j-1) * node_spacing)
        k <- k + 1
      }
    }
  }
  
  # Receptors layout (line shape)
  y_start_receptors <- - (LR_height * node_spacing) / 2
  k <- 1
  for (i in 1:width_receptors) {
    for (j in 1:LR_height) {
      if (k <= nrow(receptors)) {
        tmp_list[[receptors$var_names[k]]] <- c(x_start_receptors + (i-1) * node_spacing, y_start_receptors + (j-1) * node_spacing)
        k <- k + 1
      }
    }
  }
  
  # Genes Receivers layout (rectangle shape)
  y_start_receivers <- - (network_height * node_spacing) / 2
  k <- 1
  for (i in 1:width_receivers) {
    for (j in 1:network_height) {
      if (k <= nrow(genes_receivers)) {
        tmp_list[[genes_receivers$var_names[k]]] <- c(x_start_receivers + (i-1) * node_spacing, y_start_receivers + (j-1) * node_spacing)
        k <- k + 1
      }
    }
  }
  
  # Metadata: Position above the corresponding squares in a horizontal line
  if ("metadata" %in% unique(state_order$group)) {
    metadata <- filter(state_order, group == "metadata")
    sender_metadata_count <- 0
    receiver_metadata_count <- 0
    unknown_metadata_count <- 0
    metadata_offset <- node_spacing * 2
    for (i in 1:nrow(metadata)) {
      var_name <- metadata$var_names[i]
      if (grepl("_senders$", var_name)) {
        x_pos <- x_start_senders + sender_metadata_count * node_spacing
        y_pos <- y_start_senders + network_height * node_spacing + metadata_offset
        sender_metadata_count <- sender_metadata_count + 1
      } else if (grepl("_receivers$", var_name)) {
        x_pos <- x_start_receivers + receiver_metadata_count * node_spacing
        y_pos <- y_start_receivers + network_height * node_spacing + metadata_offset
        receiver_metadata_count <- receiver_metadata_count + 1
      } else {
        x_pos <- x_start_ligands + unknown_metadata_count * node_spacing
        y_pos <- y_start_ligands + LR_height * node_spacing + metadata_offset
        unknown_metadata_count <- unknown_metadata_count + 1
      }
      tmp_list[[var_name]] <- c(x_pos, y_pos)
    }
  }

  network_layout[["savedPositions"]] <- tmp_list

  json <- rjson::toJSON(network_layout)

  return(json)

}



#-------------------------------------------------------------------------------
# interactmiic.CCClinks
#-------------------------------------------------------------------------------
#' This function identifies cell-cell communication links to include in the 
#' interactMIIC network using the LIANA package.
#'
#' @description This function generates an LIANA output of senders and receivers
#' significant interactions.
#'
#'
#' @param seurat_object A Seurat object containing single-cell data.
#'
#' @param species A character string indicating the species. Either "human" 
#' (default) or "mouse".
#'
#' @param assay_name A character string specifying the assay name in the 
#' Seurat object. Default is "RNA".
#'
#' @param senders A non-empty character vector of sender cell types.
#'
#' @param receivers A non-empty character vector of receiver cell types.
#'
#' @param interact_ident A character string specifying the identity column in 
#' the Seurat object where senders and receivers can be found.
#' 
#' @param CCC_method A character string specifying the CCC method to use, see 
#' show_methods() of LIANA
#' 
#' @param n_CCClinks A numeric specifiying how many L-R to keep. Default is 60.
#'
#' @return A data frame containing the top cell-cell communication links.
#'
#' @examples
#' \dontrun{
#' result <- interactmiic.CCClinks(seurat_obj, species = "human", senders = c("T_cells"), 
#' receivers = c("B_cells"), interact_ident = "celltype")
#' }
#'
#' @export
#'
#-------------------------------------------------------------------------------


interactmiic.CCClinks <- function(seurat_object,
                                  species = "human",
                                  assay_name = 'RNA',
                                  CCC_method = "cellphonedb",
                                  n_CCClinks = 30,
                                  senders,
                                  receivers,
                                  interact_ident) {

  # List of callable methods score and p_values
  CCC_score = list(
    connectome = c("weight_sc",NA),
    logfc = c("logfc_comb",NA),
    natmi = c("prod_weight","edge_specificity"),
    sca = c("LRscore", NA),
    cellphonedb = c("lr.mean","pvalue")
  )
  # Input checks
  if (!inherits(seurat_object, "Seurat")) {
    stop("seurat_object must be a Seurat object.")
  }
  if (!is.character(species) || !(species %in% c("mouse", "human"))) {
    stop("species must be either 'mouse' or 'human'.")
  }
  if (!is.character(assay_name) || assay_name == "") {
    stop("assay_name must be a non-empty string.")
  }
  if (!is.character(interact_ident) || interact_ident == "") {
    stop("interact_ident must be a non-empty string.")
  }
  if (!is.character(senders) || length(senders) == 0) {
    stop("senders must be a non-empty character vector.")
  }
  if (!is.character(receivers) || length(receivers) == 0) {
    stop("receivers must be a non-empty character vector.")
  }
  if ((length(CCC_method) !=1) || !(CCC_method %in%  names(CCC_score))) {
    stop(paste("CCC_method must be an unique method selected in ", paste( names(CCC_score), collapse = ","), "as LIANA did not include more methods yet."))
  }
  if (!is.numeric(n_CCClinks) || (n_CCClinks<= 0) || (n_CCClinks != floor(n_CCClinks))) {
    stop("n_CCClinks must be a positive integer.")
  }

  
 if (species == "mouse") {
    # Convert LIANA's Consensus resource to murine symbols
    op_resource <- select_resource("Consensus")[[1]]
    
    # Generate orthologous resource
    ortholog_resource <- generate_homologs(op_resource = op_resource,
                                           target_organism = 10090) # mouse
    
    # Run LIANA with the orthologous resource
    liana_res <- liana_wrap(seurat_object,
                            method = CCC_method,
                            idents_col = interact_ident,
                            assay = assay_name,
                            resource = 'custom', # resource has to be set to 'custom' to work with external resources
                            external_resource = ortholog_resource) # provide orthologous resource  
    
  } else {
    liana_res <- liana_wrap(seurat_object,
                            method = CCC_method,
                            idents_col = interact_ident,
                            assay = assay_name)
  }
  liana_links <- liana_res %>%
     dplyr::filter((source %in% senders) & (target %in% receivers)) 
  liana_top <- liana_links[order(liana_links[,CCC_score[[CCC_method]][1]], decreasing = T),]
  
  # If a measure of specificity is present, filter only significant edges
  if (!(is.na(CCC_score[[CCC_method]][2]))) {
    liana_top <- liana_top[liana_top[CCC_score[[CCC_method]][2]] <=0.05,] 
  }
  
  return(liana_top[1:min(nrow(liana_links),n_CCClinks),])
  
}