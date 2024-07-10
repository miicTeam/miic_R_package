#*******************************************************************************
# Filename   : interactmiic.wrapper.R                Creation date: 13 june 2024
#
# Description: Wrapper function for interactMIIC
#
# Author     : Louise DUPUIS
#*******************************************************************************

#-------------------------------------------------------------------------------
# interactmiic.wrapper
#-------------------------------------------------------------------------------
#' Run the preprocessing of an Seurat object for interactMIIC
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

interactmiic.wrapper <- function(seurat_object,
                                assay_name = "RNA",
                                interact_ident = NULL,
                                senders_name = NULL,
                                receivers_name = NULL,
                                goi_senders = NULL,
                                goi_receivers = NULL,
                                genes_senders = NULL,
                                genes_receivers = NULL,
                                species = 'human',
                                do_CCC = TRUE,
                                CCC_method = "cellphonedb",
                                n_CCClinks = 20,
                                ligands = NULL,
                                receptors = NULL,
                                metadata_senders = NULL,
                                metadata_receivers = NULL,
                                wd_path = getwd()) {

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
  if (!is.null(senders_name) & (!is.character(senders_name) || length(senders_name) == 0)) {
    stop("senders_name must be a non-empty character vector.")
  }
  if (!is.null(receivers_name) & (!is.character(receivers_name) || length(receivers_name) == 0)) {
    stop("receivers_name must be a non-empty character vector.")
  }
  if (!is.numeric(n_CCClinks) || (n_CCClinks<= 0) || (n_CCClinks != floor(n_CCClinks))) {
    stop("n_CCClinks must be a positive integer.")
  }

  #Automatic detection when they did not give parameters

  if (is.null(interact_ident)) {
    print("No metadata name was given in interact_ident to find the interacting populations, the active.ident will be used and called 'interactmiic_ident'...")
    seurat_object$interactmiic_ident <- seurat_object@active.ident
    interact_ident <- "interactmiic_ident"
    interact_levels <- unique(SeuratObject@meta.data[[interact_ident]])
    if (length(interact_levels) >=2 || (!(is.null(senders_name)) & !is.null(receivers_name))) {
      if (is.null(senders_name)) {
        senders_name <- interact_levels[1]
        print(paste("No senders_name was given, the senders population automatically selected is the", senders_name, "from the active.ident"))
      } else if (!(senders_name %in% interact_levels)) {
        stop(paste("The senders population", senders_name, "is not in the active.ident levels. Did you forget to specify the metadata column ?"))
      }
      if (is.null(receivers_name)) {
        receivers_name <- interact_levels[2]
        print(paste("No receivers_name was given, the receivers population automatically selected is the", receivers_name, "from the active.ident"))
      } else if (!(receivers_name %in% interact_levels)) {
        stop(paste("The receivers population", receivers_name, "is not in the active.ident levels. Did you forget to specify the metadata column ?"))
      }
    } else {
      stop("There is less than 2 levels in the active.ident we can't automatically select senders and receivers populations. Did you forget to specify the metadata column ?")
    }
  } else {
    interact_levels <- unique(SeuratObject@meta.data[[interact_ident]])
    if (length(interact_levels) >=2 || (!(is.null(senders_name)) & !is.null(receivers_name))) {
        if (is.null(senders_name)) {
          senders_name <- interact_levels[1]
          print(paste("No senders_name was given, the senders population automatically selected is the", senders_name, "from the",interact_ident))
        } else if (!(senders_name %in% interact_levels)) {
          stop(paste("The senders population", senders_name, "is not in", interact_ident))
        }
        if (is.null(receivers_name)) {
          receivers_name <- interact_levels[2]
          print(paste("No receivers_name was given, the receivers population automatically selected is the", receivers_name, "from the",interact_ident))
        } else if (!(receivers_name %in% interact_levels)) {
          stop(paste("The receivers population", receivers_name, "is not in", interact_ident))
        }
      } else {
        stop("There is less than 2 levels in",interact_ident,", we can't automatically select senders and receivers populations. Did you the right metadata column ?")
      }
  }

  cat(paste("Wrapping your Seurat object to run interactMIIC:\n",
  senders_name, "are the senders population\n",
  receivers_name, "are the receivers population\n",
  "and are present in the", interact_ident, "column."))
  Sys.sleep(2)

  Idents(seurat_object) <- interact_ident
  seurat_object <- subset(seurat_object, idents = c(senders_name, receivers_name))

  if (do_CCC) {
    #get the LIANA significant found interactions from senders to receivers
    cat(paste("Searching the top", n_CCClinks, "L-R interactions using LIANA", CCC_method, "pipeline ..."))
  result <- interactmiic.CCClinks(seurat_object,
                                  species = species,
                                  assay_name = assay_name,
                                  interact_ident = interact_ident,
                                  senders = senders_name,
                                  receivers = receivers_name,
                                  CCC_method = CCC_method,
                                  n_CCClinks = n_CCClinks)

    write.table(result, file = file.path(wd_path,paste0("pbmc_liana_", CCC_method, ".csv")), quote = F, sep = ",")
    interact_edges <- data.frame(ligands = character(), receptors = character())

    for (i in 1:nrow(result)) {
      oneligand <- result$ligand.complex[i]
      onereceptor <- result$receptor.complex[i]
      onereceptor <- unique(unlist(str_split(onereceptor, "_")))
      for (onerecp in onereceptor) {
        interact_edges[nrow(interact_edges) +1,] <- c(oneligand, onerecp)
      }
    }

    interact_edges <- interact_edges[!(duplicated(interact_edges)),]
    #final genes lists
    ligands_CCC <- unique(interact_edges$ligand)
    receptors_CCC <- unique(interact_edges$receptors)
  } else {
    ligands <- intersect(ligands, row.names(seurat_object))
    receptors <- intersect(receptors, row.names(seurat_object))
    cat("LIANA search for L-R interactions will not be done as do_CCC is set as False.")
    ligands_CCC <- NULL
    receptors_CCC <- NULL
  }

  ligands <- unique(c(ligands_CCC,ligands))
  receptors <- unique(c(receptors_CCC,receptors))





}