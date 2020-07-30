#' MIIC, causal network learning algorithm including latent variables
#'
#' @description MIIC (Multivariate Information based Inductive Causation) combines
#' constraint-based and information-theoretic approaches to disentangle direct
#' from indirect effects amongst correlated variables, including cause-effect
#' relationships and the effect of unobserved latent causes.
#'
#' @details for non temporal series, the method, starting from a complete graph, 
#' iteratively removes dispensable edges, by uncovering significant information 
#' contributions from indirect paths, and assesses edge-specific confidences 
#' from randomization of available data. The remaining edges are then oriented 
#' based on the signature of causality in observational data.
#' 
#' For temporal series, miic reorganizes the dataset using the \emph{tau}
#' and \emph{delta_tau} parameters to transform the timesteps 
#' into lagged samples. As starting point, a lagged graph is created with 
#' only edges having at least one node laying on the last timestep. 
#' Then, miic standard algorithm is applied to remove dispensable edges. 
#' The remaining edges are then oriented by using the temporality and the 
#' signature of causality in observational data.
#'
#' @references
#' \itemize{
#' \item Verny et al., \emph{PLoS Comp. Bio. 2017.}
#' }
#'
#' @param inputData [a data frame]\cr
#' A data frame that contains the observational data. Each column corresponds 
#' to one variable and each row is a sample that gives the values for all 
#' the observed variables. The column names correspond to the names of 
#' the observed variables.\cr
#' Numeric columns will be treated as continuous values, factors 
#' and character as categorical.\cr
#' For temporal series, (when \emph{tau} parameter is >= 1), in addition 
#' to the variables, the first column must provide the timestep information 
#' in ascending order for each time series.
#'
#' @param blackBox [a data frame]
#' An optional data frame containing the pairs of variables that should be
#' considered as independent. Each row contains one column for each of 
#' the two variables. The variable name must correspond to the one in 
#' the \emph{inputData} data frame.
#'
#' @param neff [a positive integer]
#' The N samples given in the \emph{inputdata} data frame are expected
#' to be independent. In case of correlated samples such as in 
#' Monte Carlo sampling approaches, the effective number of independent samples
#' \emph{neff} can be estimated using the decay of the autocorrelation function
#' (Verny \emph{et al.}, PLoS Comp. Bio. 2017). This \emph{effective} number \emph{neff}
#' of \emph{independent} samples can be provided using this parameter.
#'
#' @param cplx [a string; \emph{c("nml", "mdl")}]
#' In practice, the finite size of the input
#' dataset requires that the 2-point and 3-point information measures should be
#' \emph{shifted} by a \emph{complexity} term. The finite size corrections can be
#' based on the Minimal Description Length (MDL) criterion (set the option with "mdl").
#' In practice, the MDL complexity criterion tends to underestimate the relevance of
#' edges connecting variables with many different categories, leading to the removal of
#' false negative edges. To avoid such biases with finite datasets, the (universal)
#' Normalized Maximum Likelihood (NML) criterion can be used (set the option with "nml").
#' The default is "nml" (see Affeldt \emph{et al.}, UAI 2015).
#'
#' @param latent [a string; \emph{c("no", "yes", "ort")}]
#' When set to "yes", the network reconstruction is taking into account hidden (latent)
#' variables. When set to "ort", latent variables are not considered during the skeleton
#' reconstruction but allows bi-directed edges during the orientation. Dependence
#' between two observed variables due to a latent variable is indicated with a '6' in
#' the adjacency matrix and in the network edges.summary and by a bi-directed edge
#' in the (partially) oriented graph.
#'
#' @param orientation [a boolean value]
#' The miic network skeleton can be partially directed by orienting and 
#' propagating edge directions. The orientation is based on the sign and 
#' magnitude of the conditional 3-point information of unshielded triples
#' and, for temporal graphs, on time lags.
#' The propagation procedure relyes on probabilities; for more details, 
#' see Verny \emph{et al.}, PLoS Comp. Bio. 2017).
#' If set to FALSE the orientation step is not performed.
#'
#' @param propagation [a boolean value]
#' If set to FALSE, the skeleton is partially oriented with only the
#' v-structure orientations, plus time for temporal graphs. 
#' Otherwise, the v-structure orientations are propagated to downstream 
#' undirected edges in unshielded triples following the orientation method.
#'
#' @param categoryOrder [a data frame] An optional data frame giving information
#' about how to order the various states of categorical variables. It will be
#' used to compute the signs of the edges (using partial correlation coefficient)
#' by sorting each variable's levels accordingly to the given category order.
#'
#' @param trueEdges [a data frame]  An optional data frame containing all the
#' true edges of the graph. Each line corresponds to one edge.
#'
#' @param edges [a data frame] The miic$edges object returned by an execution
#' of the miic function. It represents the result of the skeleton step. If this
#' object is provided, the skeleton step will not be done, and the required
#' orientation will be performed using this edges data frame.
#'
#' @param confidenceShuffle [a positive integer] The number of shufflings of
#' the original dataset in order to evaluate the edge specific confidence
#' ratio of all inferred edges.
#'
#' @param confidenceThreshold [a positive floating point] The threshold used
#' to filter the less probable edges following the skeleton step. See Verny
#' \emph{et al.}, PLoS Comp. Bio. 2017.
#'
#' @param confList [a data frame] An optional data frame containing the
#' confFile data frame returned by a miic execution. It is useful when a
#' second run of the same input data set has to be performed with a different
#' confidence threshold and the same confidenceShuffle value. In
#' this way the computations based on the randomized dataset do not need to
#' be performed again, and the values in this data frame are used instead.
#'
#' @param sampleWeights [a numeric vector]
#' An optional vector containing the weight of each observation.
#'
#' @param testMAR [a boolean value]
#' If set to TRUE, distributions with missing values will be tested with Kullback-Leibler
#' divergence : conditioning variables for the given link \eqn{X\arrow Y}\eqn{Z} will be
#' considered only if the divergence between the full distribution and the non-missing
#' distribution \eqn{KL(P(X,Y) | P(X,Y)_{!NA})} is low enough (with \eqn{P(X,Y)_{!NA}} as
#' the joint distribution of \eqn{X} and \eqn{Y} on samples which are not missing on Z.
#' This is a way to ensure that data are missing at random for the considered
#' interaction and to avoid selection bias. Set to TRUE by default
#'
#' @param consistent [a string; \emph{c("no", "orientation", "skeleton")}]
#' if "orientation": iterate over skeleton and orientation steps to ensure
#' consistency of the network;
#' if "skeleton": iterate over skeleton step to get a consistent skeleton, then
#' orient edges and discard inconsistent orientations to ensure consistency of
#' the network.
#'
#' @param nThreads [a positive integer]
#' When set greater than 1, nThreads parallel threads will be used for computation. Make sure
#' your compiler is compatible with openmp if you wish to use multithreading.
#'
#' @param tau [an integer] Optional, -1 by default.\cr
#' Max lag used for temporal series. If \emph{tau} is supplied (integer >= 1), 
#' miic switches to temporal mode: it contructs a lagged graph over 
#' \emph{tau} / \emph{delta_tau} periods of time, and looks both for temporal 
#' and contemporaneous edges.\cr
#' Note that if \emph{delta_tau} is also supplied, \emph{tau} must be a
#' multiple of \emph{delta_tau}.
#' 
#' @param movavg [an integer] Optional, -1 by default.\cr
#' Used only in temporal mode, If \emph{movavg} is supplied (integer > 1), 
#' a moving average operation is applied to each time series.\cr
#' 
#' @param delta_tau [an integer] Optional, 1 by default.\cr
#' When \emph{delta_tau} is supplied (integer > 1), the samples will be 
#' construted using 1 timestep every \emph{delta_tau} timesteps starting 
#' from the last.\cr 
#' i.e.: on 1000 timesteps with  \emph{tau} = 14 and \emph{delta_tau} = 7, 
#' the timesteps kept for the samples conversion will be 1000, 993, 986 
#' for the first sample, the next sample will use 999, 992, 985 and so on.\cr
#'  
#' @param verbose [a boolean value] If TRUE, debugging output is printed.
#'
#' @return A \emph{miic-like} object that contains:
#' \itemize{
#'  \item{all.edges.summary:}{ a data frame with information about the relationship between
#'  each pair of variables
#'  \itemize{
#'  \item \emph{x:} X node
#'  \item \emph{y:} Y node
#'  \item \emph{type:} contains 'N' if the edge has
#'  been removed or 'P' for retained edges. If a true edges file is given,
#'  'P' becomes 'TP' (True Positive) or 'FP' (False Positive), while
#'  'N' becomes 'TN' (True Negative) or 'FN' (False Negative).
#'  \item \emph{ai:} the contributing nodes found by the method which participate in
#'  the mutual information between \emph{x} and \emph{y}, and possibly separate them.
#'  \item \emph{info:} provides the final mutual information times \emph{Nxy_ai} for
#'  the pair (\emph{x}, \emph{y}) when conditioned on the collected nodes \emph{ai}.
#'  \item \emph{cplx:} gives the computed complexity between the (\emph{x}, \emph{y})
#'  variables taking into account the contributing nodes \emph{ai}.
#'  \item \emph{Nxy_ai:} gives the number of samples on which the information and the
#'  complexity have been computed. If the input dataset has no missing value, the
#'  number of samples is the same for all pairs and corresponds to the total
#'  number of samples.
#'  \item \emph{log_confidence:} represents the \emph{info} - \emph{cplx} value.
#'  It is a way to quantify the strength of the edge (\emph{x}, \emph{y}).
#'  \item \emph{confidenceRatio:} this column is present if the confidence cut
#'  is > 0 and it represents the ratio between the probability to reject
#'  the edge (\emph{x}, \emph{y}) in the dataset versus the mean probability
#'  to do the same in multiple (user defined) number of randomized datasets.
#'  \item \emph{infOrt:} the orientation of the edge (\emph{x}, \emph{y}). It is
#'  the same value as in the adjacency matrix at row \emph{x} and column \emph{y}.
#'  \item \emph{trueOrt:} the orientation of the edge (\emph{x}, \emph{y}) present
#'  in the true edges file (if true edges file is provided).
#'  \item \emph{isOrtOk:} information about the consistency of the inferred graph's
#'  orientations with a reference graph is given (i.e. if true edges file is provided).
#'  Y: the orientation is consistent; N: the orientation is not consistent with
#'  the PAG derived from the given true graph.
#'  \item \emph{sign:} the sign of the partial correlation between variables
#'  \emph{x} and \emph{y}, conditioned on the contributing nodes \emph{ai}.
#'  \item \emph{partial_correlation:} value of the partial correlation for the
#'  edge (\emph{x}, \emph{y}) conditioned on the contributing nodes \emph{ai}.
#'  }
#'  }
#'
#'  \item{retained.edges.summary:} {a data frame in the format of all.edges.summary containing only the inferred edges.}
#'
#'  \item{orientations.prob:} {this data frame lists the orientation probabilities of the two edges of all unshielded triples
#'  of the reconstructed network with the structure: node1 -- mid-node -- node2:}
#'  \itemize{
#'  \item node1: node at the end of the unshielded triplet
#'  \item p1: probability of the arrowhead node1 <- mid-node
#'  \item p2: probability of the arrowhead node1 -> mid-node
#'  \item mid-node: node at the center of the unshielded triplet
#'  \item p3: probability of the arrowhead mid-node <- node2
#'  \item p4: probability of the arrowhead mid-node -> node2
#'  \item node2: node at the end of the unshielded triplet
#'  \item NI3: 3 point (conditional) mutual information * N
#'  }
#'
#'
#'  \item {AdjMatrix:} the adjacency matrix is a square matrix used to represent
#'  the inferred graph. The entries of the matrix indicate whether pairs of
#'  vertices are adjacent or not in the graph. The matrix can be read as a
#'  (row, column) set of couples where the row represents the source node and
#'  the column the target node. Since miic can reconstruct mixed networks
#'  (including directed, undirected and bidirected edges), we will have a
#'  different digit for each case:
#'  \itemize{
#'  \item 1: (\emph{x}, \emph{y}) edge is undirected
#'  \item 2: (\emph{x}, \emph{y}) edge is directed as \emph{x} -> \emph{y}
#'  \item -2: (\emph{x}, \emph{y}) edge is directed as \emph{x} <- \emph{y}
#'  \item 6: (\emph{x}, \emph{y}) edge is bidirected
#'  }
#' }
#' @export
#' @useDynLib miic
#'
#' @examples
#' library(miic)
#'
#' # EXAMPLE HEMATOPOIESIS
#' data(hematoData)
#'
#' # execute MIIC (reconstruct graph)
#' miic.res <- miic(
#'   inputData = hematoData, latent = TRUE,
#'   confidenceShuffle = 10, confidenceThreshold = 0.001
#' )
#'
#' # plot graph
#' miic.plot(miic.res)
#' \dontrun{
#'
#' # write graph to graphml format. Note that to correctly visualize
#' # the network we created the miic style for Cytoscape (http://www.cytoscape.org/).
#'
#' miic.write.network.cytoscape(g = miic.res, file = file.path(tempdir(), "temp"))
#'
#' # EXAMPLE CANCER
#' data(cosmicCancer)
#' data(cosmicCancer_stateOrder)
#' # execute MIIC (reconstruct graph)
#' miic.res <- miic(
#'   inputData = cosmicCancer, categoryOrder = cosmicCancer_stateOrder, latent = "yes",
#'   confidenceShuffle = 100, confidenceThreshold = 0.001
#' )
#'
#' # plot graph
#' miic.plot(miic.res, igraphLayout = igraph::layout_on_grid)
#'
#' # write graph to graphml format. Note that to correctly visualize
#' # the network we created the miic style for Cytoscape (http://www.cytoscape.org/).
#' miic.write.network.cytoscape(g = miic.res, file = file.path(tempdir(), "temp"))
#'
#' # EXAMPLE OHNOLOGS
#' data(ohno)
#' data(ohno_stateOrder)
#' # execute MIIC (reconstruct graph)
#' miic.res <- miic(
#'   inputData = ohno, latent = "yes", categoryOrder = ohno_stateOrder,
#'   confidenceShuffle = 100, confidenceThreshold = 0.001
#' )
#'
#' # plot graph
#' miic.plot(miic.res)
#'
#' # write graph to graphml format. Note that to correctly visualize
#' # the network we created the miic style for Cytoscape (http://www.cytoscape.org/).
#' miic.write.network.cytoscape(g = miic.res, file = file.path(tempdir(), "temp"))
#' }
#'
miic <- function(inputData,
                 categoryOrder = NULL,
                 trueEdges = NULL,
                 blackBox = NULL,
                 nThreads = 1,
                 cplx = c("nml", "mdl"),
                 orientation = TRUE,
                 propagation = TRUE,
                 latent = c("no", "yes", "ort"),
                 neff = -1,
                 edges = NULL,
                 confidenceShuffle = 0,
                 confidenceThreshold = 0,
                 confList = NULL,
                 sampleWeights = NULL,
                 testMAR = TRUE,
                 consistent = c("no", "orientation", "skeleton"),
                 tau = -1,
                 movavg = -1,
                 delta_tau = 1,
                 verbose = FALSE
                 ) {
  res <- NULL
  skeleton <- TRUE

  #### Check the input arguments
  if (is.null(inputData)) {
    stop("The input data file is required")
  }
  if (!is.data.frame(inputData)) 
    {
    stop("The input data is not a dataframe")
    }
  
  if (tau > 0)
    {
    # If we use temporal version of miic, convert history into lagged nodes and samples
    #
    cat ("Using temporal version of miic\n")
    struct_ret <- tmiic.transform_data_for_miic (inputData, tau, 
        categoryOrder=categoryOrder, movavg=movavg, delta_tau=delta_tau)
    inputData <- struct_ret$inputData
    categoryOrder <- struct_ret$categoryOrder
    }

  effnAnalysis <- miic.evaluate.effn(inputData, plot = F)
  if (effnAnalysis$neff < 0.5 * nrow(inputData)) {
    if (effnAnalysis$exponential_decay) {
      warning(
        paste0(
          "Your samples in the datasets seem to be correlated! We ",
          "suggest to re run the method specifying ",
          effnAnalysis$neff,
          " in the neff parameter. See the ",
          "autocorrelation plot for more details."
        )
      )
    } else {
      warning(
        paste0(
          "Your samples in the datasets seem to be correlated but ",
          "the correlation decay is not exponential. Are your ",
          "samples correlated in some way? See the autocorrelation ",
          "plot for more details."
        )
      )
    }
  }

  # if (!is.null(trueEdges)) {
  #   if (length(which(!c(as.vector(trueEdges[,1]), as.vector(trueEdges[,2])) %in% colnames(inputData))) > 0){
  #     stop("True edges file does not correspond to the input data matrix. Please check node names.")
  #   }
  # }

  cplx <- tryCatch(
    {
      match.arg(cplx)
    },
    error = function(e) {
      if (grepl("object .* not found", e$message)) {
        message(e, "")
        return("")
      }
      return(toString(cplx))
    }
  )
  cplx <- match.arg(cplx)

  latent <- tryCatch(
    {
      match.arg(latent)
    },
    error = function(e) {
      if (grepl("object .* not found", e$message)) {
        message(e, "")
        return("")
      }
      return(toString(latent))
    }
  )
  latent <- match.arg(latent)

  consistent <- tryCatch(
    {
      match.arg(consistent)
    },
    error = function(e) {
      if (grepl("object .* not found", e$message)) {
        message(e, "")
        return("")
      }
      return(toString(consistent))
    }
  )
  consistent <- match.arg(consistent)


  if (neff > nrow(inputData)) {
    stop(
      paste0(
        "The number of effective samples cannot be greater than the ",
        "number of samples."
      )
    )
  }
  # edges
  if (!is.null(edges)) {
    skeleton <- FALSE
    if (length(colnames(edges)) != 10) {
      stop(
        paste0(
          "The edges data frame is not correct. The required data ",
          "frame is the $edges output of the miic method"
        )
      )
    }
  }

  # propagation
  if (propagation != TRUE && propagation != FALSE) {
    stop("The propagation type is not correct. Allowed types are TRUE or FALSE")
  }

  # orientation
  if (orientation != TRUE && orientation != FALSE) {
    stop("The orientation type is not correct. Allowed types are TRUE or FALSE")
  }

  if (verbose) {
    cat("START miic...\n")
  }

  # continuous or discrete?
  cntVar <- sapply(inputData, is.numeric)
  for (col in colnames(inputData)) {
    unique_values <- length(unique(inputData[[col]][!is.na(inputData[[col]])]))
    if (cntVar[[col]] &&
      (unique_values <= 40) && (nrow(inputData) > 40)) {
      if (cntVar[[col]] &&
        (unique_values <= 2) && (nrow(inputData) > 40)) {
        stop(
          paste0(
            "Numerical variable ",
            col,
            " only has ",
            unique_values,
            " non-NA unique values. Is this a factor ?"
          )
        )
      }
      warning(
        paste0(
          "Numerical variable ",
          col,
          " is treated as continuous but only has ",
          unique_values,
          " unique values."
        )
      )
    }
    if ((!cntVar[[col]]) &&
      (unique_values >= 40) && (nrow(inputData) > 40)) {
      warning(paste0(
        col,
        " is treated as discrete but has many levels (",
        unique_values,
        ")."
      ))
    }
  }
  typeOfData <- 0 # Assume all discrete
  if (any(cntVar)) {
    typeOfData <- 2 # Mixed if any are continuous
    if (all(cntVar)) {
      typeOfData <- 1 # All continuous
    }
  }

    # STATE ORDER FILE
    if (!is.null(categoryOrder)) {
      err_code <- checkStateOrder(categoryOrder, inputData)
      if (err_code != "0") {
        print(errorCodeToString(err_code))
        print("WARNING: Category order file will be ignored!")
        categoryOrder <- NULL
      }
    }

  err_code <- checkInput(inputData, "miic")
  if (err_code != "0") {
    print(errorCodeToString(err_code))
  } else {
    if (verbose) {
      cat("\t# -> START reconstruction...\n")
    }
    res <-
      miic.reconstruct(
        inputData = inputData,
        stateOrder = categoryOrder,
        nThreads = nThreads,
        cplx = cplx,
        latent = latent,
        effN = neff,
        blackBox = blackBox,
        confidenceShuffle = confidenceShuffle,
        edges = edges,
        orientation = orientation,
        propagation = propagation,
        confidenceThreshold = confidenceThreshold,
        verbose = verbose,
        cntVar = cntVar,
        typeOfData = typeOfData,
        sampleWeights = sampleWeights,
        testMAR = testMAR,
        consistent = consistent,
        tau = tau
      )
    if (res$interrupted) {
      stop("Interupted by user")
    }
    time <- res$time
    if (verbose) {
      cat("\t# -> END reconstruction...\n\t# --------\n")
    }

    if (!is.null(trueEdges)) {
      err_code <- checkTrueEdges(trueEdges)
      if (err_code != "0") {
        print(errorCodeToString(err_code))
        print("WARNING: True edges file will be ignored!")
        trueEdges <- NULL
      }
    }

    # Summarize the results
    # --------
    res$all.edges.summary <- summarizeResults(
      observations = inputData,
      edges = res$edges,
      true_edges = trueEdges,
      state_order = categoryOrder,
      adj_matrix = res$adjMatrix,
      orientation_probabilities = res$orientations.prob,
      verbose = verbose
    )
  }

  res
}
