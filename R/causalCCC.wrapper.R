#*******************************************************************************
# Filename   : causalCCC.wrapper.R                Creation date: 13 june 2024
#
# Description: Wrapper function for causalCCC
#
# Author     : Louise DUPUIS
#*******************************************************************************

#-------------------------------------------------------------------------------
# causalCCC.wrapper
#-------------------------------------------------------------------------------
#' Run the preprocessing of an Seurat object for causalCCC
#'
#' @description Function that takes an data_input raw counts,
#'  senders and receivers information metadata and create the
#' output files necessary to run causalCCC
#'
#' @import Seurat
#' @import data.table
#'
#' @param data_input [a dataframe or a Seurat object]
#' A Single-Cell transcriptomics object, dataframe or Seurat. If
#' dataframe must contains genes and metadata as variables and
#' cells as observations.
#' 
#' @param assay_name [a string]. Gives the name of the assay to take
#' the transcriptomics raw counts from (usually 'RNA')
#'
#' @param interact_ident [a string]. Gives the name of the metadata
#' containing the celltypes population
#'
#' @param senders_name [a string]. Gives the name of the senders population.
#'
#' @param receivers_name [a string]. Gives the name of the receivers population.
#'
#' @param do_MIselect A boolean. If TRUE, perform Mutual Information-based feature selection.
#' 
#' @param save A boolean. If TRUE, save the MI table of MIselection to wd_path.
#' 
#' @param plot A boolean. If TRUE, generate MI heatmaps plots and save them in wd_path.
#' 
#' @param n_genes_senders [a numeric] The number of genes to keep after pairwise
#' mutual information ranking.  Default depends on the length of goi
#' (between 15 and 100).
#' 
#' @param n_genes_receivers [a numeric] The number of genes to keep after pairwise
#' mutual information ranking.  Default depends on the length of goi
#' (between 15 and 100).
#' 
#' @param goi_senders A vector. A list of genes of interest for the senders cells.
#' 
#' @param goi_receivers A vector. A list of genes of interest for the receivers cells.
#'
#' @param genes_senders [a vector] An optional list of genes to add to the senders network.
#'
#' @param genes_receivers [a vector] An optional list of genes to add to the receivers network.
#'
#' @param metadata_senders [a named list] A named list of metadata with their levels
#'  that will be added to the genes network. They are relevant to senders. See exemple
#'
#' @param metadata_receivers [a named list] A named list of metadata with their levels
#'  that will be added to the genes network. They are relevant to receivers. See exemple
#' 
#' @param species A string. The species ('mouse' or 'human').
#' 
#' @param do_CCC A boolean. If TRUE, perform Cell-Cell Communication analysis.
#' 
#' @param do_layout A boolean. If TRUE, create a network layout.json file for causalCCC web server.
#' 
#' 
#' @param CCC_method A character string specifying the CCC method to use, see 
#' show_methods() of LIANA
#' 
#' @param n_CCClinks A numeric specifiying how many L-R to keep. Default is 60.
#' 
#' @param ligands A vector. An optional list of ligand genes to add to the network.
#' 
#' @param receptors A vector. An optional list of receptor genes to add to the network.
#' 
#' @param wd_path A string. The working directory path for saving results.
#' 
#' @examples
#' \donttest{
#' metadata <- list(
#'  treatment = c("Control", "Treated"),
#'  another_meta = c("Level1", "Level2", "Level3"),
#'  continusous_meta = NULL
#')
#' }
#'
#' @return A list containing mosaic datatable, state order, interact edges, 
#' and network layout for causalCCC network reconstruction.
#'
#' @export
#'
#-------------------------------------------------------------------------------

causalCCC.wrapper <- function(data_input,
                                assay_name = "RNA",
                                interact_ident = NULL,
                                senders_name = NULL,
                                receivers_name = NULL,
                                do_CCC = FALSE,
                                CCC_method = "cellphonedb",
                                n_CCClinks = 20,
                                interact_edges = NULL,
                                do_MIselect = FALSE,
                                n_genes_senders = NULL,
                                n_genes_receivers = NULL,
                                save = TRUE,
                                plot = TRUE,
                                goi_senders = NULL,
                                goi_receivers = NULL,
                                metadata_senders = list(),
                                metadata_receivers = list(),
                                species = 'human',
                                genes_senders = NULL,
                                genes_receivers = NULL,
                                do_layout = FALSE,
                                wd_path = getwd()) {

  # Input checks
  if (!is.character(species) || !(species %in% c("mouse", "human"))) {
    stop("species must be either 'mouse' or 'human'.")
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
  if (!is.null(n_genes_senders) & (!is.numeric(n_genes_senders) || (n_genes_senders<= 0) || (n_genes_senders != floor(n_genes_senders)))) {
    stop("n_genes_senders must be NULL or a positive integer.")
  }
  if (!is.null(n_genes_senders) && n_genes_senders > 500) {
    stop("n_genes_senders must be ")
  }
  if (!is.null(n_genes_receivers) & (!is.numeric(n_genes_receivers) || (n_genes_receivers<= 0) || (n_genes_receivers != floor(n_genes_receivers)))) {
    stop("n_genes_receivers must be NULL or a positive integer.")
  }
  if (!is.null(n_genes_receivers) && n_genes_receivers > 500) {
    stop("n_genes_receivers must be ")
  }
  if (!is.logical(do_MIselect)) {
    stop("do_MIselect must be a boolean value.")
  }
  if (!is.logical(save)) {
    stop("save must be a boolean value.")
  }
  if (!is.logical(plot)) {
    stop("plot must be a boolean value.")
  }
  if (!is.logical(do_CCC)) {
    stop("do_CCC must be a boolean value.")
  }
  if ((do_CCC == FALSE)&(is.null(interact_edges))) {
    stop("do_CCC is set as false but you did not give the interact_edges file. It must be the cell-cell communication links in a form of a dataframe with ligands and receptors as column names.")
  }
  if (!is.character(wd_path)) {
    stop("wd_path must be a non-empty string.")
  }
  
  if (!is.null(interact_edges) & (!is.data.frame(interact_edges) || !(all(c("ligands", "receptors") %in% colnames(interact_edges))))) {
    stop("interact_edges must be a dataframe containing two columns 'ligands' and 'receptors'.")
  }

   #Seurat input
  if (inherits(data_input, "Seurat")) {
    is_seurat <- TRUE
    if (!is.character(assay_name) || length(assay_name) != 1) {
    stop("The 'assay_name' must be a single character string.")
    }
    if (!(assay_name %in% names(data_input@assays))) {
      stop(paste("The assay", assay_name, "is not present in the data_input object."))
    }
  #Dataframe input
  } else if (is.data.frame(data_input)) {
    is_seurat <- FALSE
  } else {
    stop("The 'data_input' must be either a Seurat object or a data frame of raw counts and metadata.")
  }


  #Automatic detection when they did not give parameters
  if(is_seurat) {
    if (is.null(interact_ident)) {
      print("No metadata name was given in interact_ident to find the interacting populations, the active.ident will be used and called 'causalCCC_ident'...")
      data_input$causalCCC_ident <- data_input@active.ident
      interact_ident <- "causalCCC_ident"
      interact_levels <- unique(data_input@meta.data[[interact_ident]])
      if (length(interact_levels) >=2 || (!(is.null(senders_name)) & !is.null(receivers_name))) {
        if (is.null(senders_name)) {
          senders_name <- as.character(interact_levels[1])
          print(paste("No senders_name was given, the senders population automatically selected is the", senders_name, "from the active.ident"))
        } else if (!(senders_name %in% interact_levels)) {
          stop(paste("The senders population", senders_name, "is not in the active.ident levels. Did you forget to specify the metadata column ?"))
        }
        if (is.null(receivers_name)) {
          receivers_name <- as.character(interact_levels[2])
          print(paste("No receivers_name was given, the receivers population automatically selected is the", receivers_name, "from the active.ident"))
        } else if (!(receivers_name %in% interact_levels)) {
          stop(paste("The receivers population", receivers_name, "is not in the active.ident levels. Did you forget to specify the metadata column ?"))
        }
      } else {
        stop("There is less than 2 levels in the active.ident we can't automatically select senders and receivers populations. Did you forget to specify the metadata column ?")
      }
    } else {
      interact_levels <- unique(data_input@meta.data[[interact_ident]])
      if (length(interact_levels) >=2 || (!(is.null(senders_name)) & !is.null(receivers_name))) {
          if (is.null(senders_name)) {
            senders_name <- as.character(interact_levels[1])
            print(paste("No senders_name was given, the senders population automatically selected is the", senders_name, "from the",interact_ident))
          } else if (!(senders_name %in% interact_levels)) {
            stop(paste("The senders population", senders_name, "is not in", interact_ident))
          }
          if (is.null(receivers_name)) {
            receivers_name <- as.character(interact_levels[2])
            print(paste("No receivers_name was given, the receivers population automatically selected is the", receivers_name, "from the",interact_ident))
          } else if (!(receivers_name %in% interact_levels)) {
            stop(paste("The receivers population", receivers_name, "is not in", interact_ident))
          }
        } else {
          stop("There is less than 2 levels in",interact_ident,", we can't automatically select senders and receivers populations. Did you the right metadata column ?")
        }
    }
  } else {
     if (is.null(interact_ident)) {
      stop("You did not precise which column to use to automatically select senders and receivers populations. Please precise the name in interact_ident. ")
     } 
     if (!(interact_ident %in% colnames(data_input))){
       stop(paste(interact_ident, "is not found in the data_input columns."))
     }
     if (is.null(senders_name)) {
      stop("You did not precise which cells populations are the senders.")
    } else if (!(senders_name %in% unique(data_input[,interact_ident]))) {
      stop(paste("The senders population", senders_name, "is not in", interact_ident))
    }
     if (is.null(receivers_name)) {
      stop("You did not precise which cells populations are the receivers.")
        } else if (!(receivers_name %in% unique(data_input[,interact_ident]))) {
      stop(paste("The receivers population", receivers_name, "is not in", interact_ident))
    }
  }

  # Print to make sure the wrapping is accurate
  cat(paste("Wrapping your single-cell object to run causalCCC:\n",
  senders_name, "are the senders population\n",
  receivers_name, "are the receivers population\n",
  "and are present in the", interact_ident, "column.\n"))
  Sys.sleep(2)

 
  #Initialization
  if (is_seurat) {
     Idents(data_input) <- interact_ident
    data_input <- subset(data_input, idents = c(senders_name, receivers_name))
  } else {
    data_input <- data_input[data_input[,interact_ident] %in% c(senders_name, receivers_name),]
  }
  if (!is.null(interact_edges)) {
    ligands <- intersect(unique(interact_edges$ligands), row.names(data_input))
    receptors <- intersect(unique(interact_edges$receptors), row.names(data_input))
  }

  # CCC analysis - optional
  if (do_CCC) {
    #If they do not know the L-R links they can not know the upstream and downstream genes so we will seelct them.
    do_MIselect <- TRUE
    if (!(is_seurat)) {
      stop("Our embedded cell-cell communication analysis can only be performed with Seurat objects. Please give a Seurat object or perform your own cell-cell communication analysis")
    }

    #get the LIANA significant found interactions from senders to receivers
    cat("\n----------------------------------------------\n")
    cat(paste("Ligands-Receptors selection using LIANA", CCC_method, "pipeline ... \n"))
    cat("----------------------------------------------\n")

    result <- causalCCC.links(data_input,
                                    species = species,
                                    assay_name = assay_name,
                                    interact_ident = interact_ident,
                                    senders = senders_name,
                                    receivers = receivers_name,
                                    CCC_method = CCC_method,
                                    n_CCClinks = n_CCClinks)
    write.table(result, file = file.path(wd_path,paste0("CCCanalysis_LIANA_", CCC_method, ".csv")), quote = F, sep = ",")
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
    ligands <- c()
    receptors <- c()
  } else {
    cat("LIANA search for L-R interactions will not be done as do_CCC is set as FALSE.\n")
    ligands_CCC <- c()
    receptors_CCC <- c()
  }

  ligands <- unique(c(ligands_CCC,ligands))
  ligands <- ligands[!is.na(ligands)]
  receptors <- unique(c(receptors_CCC,receptors))
  receptors <- receptors[!is.na(receptors)]
  cat(paste("For the senders population", senders_name,
    "the ligands are ",
    paste(unique(c(ligands)), collapse = ","), "\n"))
  cat(paste("For the receivers population", senders_name,
    "the receptors are ",
    paste(unique(c(receptors)), collapse = ","), "\n"))

  cat("\n")
  if (do_MIselect) {
    #Find genes that share the most information with your biological question :
    
    cat("\n----------------------------------------------\n")
    cat(paste("Mutual Information-based feature selection ... \n"))
    cat("----------------------------------------------\n")

    goi_senders <- unique(c(goi_senders, ligands[1:min(length(ligands), 8)]))
    goi_senders <- goi_senders[!is.na(goi_senders)]
    if (length(unique(c(names(metadata_senders),goi_senders)))==0) {
      stop("There are no features of interest for the MI feature selection. Meaning no ligands, 
      no metadata and no genes of interest were detected. If you set do_CCC = TRUE maybe the selected method
      did not find any ligands receptors pairs.")
    }
    cat(paste("For the senders population", senders_name, 
    "the features of interest are ",
    paste(unique(c(names(metadata_senders),goi_senders)), collapse = ","),"\n",
    "The selection should take ~", length(unique(c(names(metadata_senders),goi_senders))),
    "minutes.\n"))
    Sys.sleep(2)
    if (!is.null(n_genes_senders)) {
      MI_senders_genes <- causalCCC.MIselection(data_input = data_input,
                                    assay_name = assay_name,
                                    interact_ident = interact_ident,
                                    oneinteract = senders_name,
                                    goi = goi_senders,
                                    metadata_list = names(metadata_senders),
                                    n_genes = n_genes_senders,
                                    save = save,
                                    output_dir = file.path(wd_path, "MI_tables"),
                                    plot = plot,
                                    color_heatmap = "darkgreen")
    } else {
      MI_senders_genes <- causalCCC.MIselection(data_input = data_input,
                                          assay_name = assay_name,
                                          interact_ident = interact_ident,
                                          oneinteract = senders_name,
                                          goi = goi_senders,
                                          metadata_list = names(metadata_senders),
                                          save = save,
                                          output_dir = file.path(wd_path, "MI_tables"),
                                          plot = plot,
                                          color_heatmap = "darkgreen")
    }
    

    goi_receivers <- unique(c(goi_receivers, receptors[1:min(length(receptors), 8)]))
    goi_receivers <- goi_receivers[!is.na(goi_receivers)]
    if (length(unique(c(names(metadata_receivers),goi_receivers)))==0) {
      stop("There are no features of interest for the MI feature selection. Meaning no receptors, 
      no metadata and no genes of interest were detected. If you set do_CCC = TRUE maybe the selected method
      did not find any ligands receptors pairs.")
    }
    cat(paste("For the receivers population", receivers_name,
    "the features of interest are ",
    paste(unique(c(names(metadata_receivers),goi_receivers)), collapse = ","),"\n",
    "The selection should take ~", length(unique(c(names(metadata_receivers),goi_receivers))),
    "minutes.\n"))
    if (!is.null(n_genes_receivers)) {
      MI_receivers_genes <- causalCCC.MIselection(data_input = data_input,
                                    assay_name = assay_name,
                                    interact_ident = interact_ident,
                                    oneinteract = receivers_name,
                                    goi = goi_receivers,
                                    metadata_list = names(metadata_receivers),
                                    n_genes = n_genes_receivers,
                                    save = save,
                                    output_dir = file.path(wd_path, "MI_tables"),
                                    plot = plot,
                                    color_heatmap = "darkorange")
    } else {
    MI_receivers_genes <- causalCCC.MIselection(data_input = data_input,
                                    assay_name = assay_name,
                                    interact_ident = interact_ident,
                                    oneinteract = receivers_name,
                                    goi = goi_receivers,
                                    metadata_list = names(metadata_receivers),
                                    save = save,
                                    output_dir = file.path(wd_path, "MI_tables"),
                                    plot = plot,
                                    color_heatmap = "darkorange")
    }
  } else {
    cat("Mutual Information-based feature selection will not be done as do_MIselect is set as FALSE.\n")
    MI_senders_genes <- c()
    MI_receivers_genes <- c()
    if (is.null(genes_senders)) {
      stop("do_MIselect is set as FALSE but you didn't give any genes_senders")
    }
    if (is.null(genes_receivers)) {
      stop("do_MIselect is set as FALSE but you didn't give any genes_receivers")
    }
  }

  #Create the necessary causalCCC files to reconstruct a network.

  ## Create final variables list
  genes_senders <- unique(c(ligands, MI_senders_genes, goi_senders, genes_senders))
  message(paste(c("Here are the genes used to reconstruct the senders MIIC network:",paste0(genes_senders, collapse = ", ")), collapse = "\n"))
   message(paste(c("Here are the ligands in the senders MIIC network:",paste0(ligands, collapse = ", ")), collapse = "\n"))
  
  if (length(metadata_senders) >0) {
    message(paste(c("Here are the metadata used to reconstruct the senders MIIC network:",paste0(c(names(metadata_senders)), collapse = ", ")), collapse = "\n"))
  }

  #RECEIVER
  genes_receivers <- unique(c(receptors, MI_receivers_genes, goi_receivers, genes_receivers ))
  message(paste(c("Here are the genes used to reconstruct the receivers MIIC network:",paste0(genes_receivers, collapse = ", ")), collapse = "\n"))
   message(paste(c("Here are the receptors in the receivers MIIC network:",paste0(receptors, collapse = ", ")), collapse = "\n"))
   
  if (length(metadata_receivers) >0) {
    message(paste(c("Here are the metadata used to reconstruct the receivers MIIC network:",paste0(c(names(metadata_receivers)), collapse = ", ")), collapse = "\n"))
  }

  ## Create the input mosaic matrix
  causalCCC_df <- causalCCC.mosaic(data_input = data_input,
                                  assay_name = assay_name,
                                  interact_ident = interact_ident,
                                  senders_name = senders_name,
                                  receivers_name = receivers_name,
                                  genes_senders= genes_senders,
                                  genes_receivers = genes_receivers,
                                  metadata_senders = names(metadata_senders),
                                  metadata_receivers = names(metadata_receivers))

  ## Create the state order
  causalCCC_st <- causalCCC.state_order(mosaic_data_table = causalCCC_df,
                                              genes_senders= genes_senders,
                                              genes_receivers = genes_receivers,
                                              ligands = ligands,
                                              receptors = receptors,
                                              metadata_senders = metadata_senders,
                                              metadata_receivers = metadata_receivers)


  ## Save files
  cat(paste("Saving the causalCCC files in", file.path(wd_path,"CausalCCC_files/")))

  ## Check existence of directory
  dir.create(file.path(wd_path,"CausalCCC_files/"))
  
  ## Create the network layout
  if (do_layout) {
    network_layout <- causalCCC.layout(causalCCC_st, network_height = 8)
    file <- file(file.path(wd_path,"CausalCCC_files/causalCCC_layout.json"))
    writeLines(network_layout, file)
    close(file)
  } else {
    network_layout <- NULL
  }


  write.table(causalCCC_df, file = file.path(wd_path,"CausalCCC_files/causalCCC_df.csv"), quote = F, sep = ",", row.names = F)
  write.table(causalCCC_st, file = file.path(wd_path, "CausalCCC_files/causalCCC_st.tsv"), quote = F, sep = "\t", row.names = F)

  #Fix ligands and receptors names when duplicated
  duplicated_ligands <- intersect(ligands,genes_receivers)
  duplicated_receptors <- intersect(receptors, genes_senders)

  interact_edges$ligands <- ifelse(interact_edges$ligands %in% duplicated_ligands, paste0(interact_edges$ligands, "_senders"), interact_edges$ligands)
  interact_edges$receptors <- ifelse(interact_edges$receptors %in% duplicated_receptors, paste0(interact_edges$receptors, "_receivers"), interact_edges$receptors)

  write.table(interact_edges, file.path(wd_path,"CausalCCC_files/causalCCC_interactEdges.tsv"), row.names=F, quote=F, sep="\t")



  output_files <- list(
    mosaic_table = causalCCC_df,
    state_order = causalCCC_st,
    interact_edges = interact_edges,
    network_layout = network_layout
  )

  return(output_files)




}