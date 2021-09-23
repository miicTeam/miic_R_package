#' MIIC, causal network learning algorithm including latent variables
#'
#' @description MIIC (Multivariate Information based Inductive Causation) combines
#' constraint-based and information-theoretic approaches to disentangle direct
#' from indirect effects amongst correlated variables, including cause-effect
#' relationships and the effect of unobserved latent causes.
#'
#' @details Starting from a complete graph, the method iteratively removes
#' dispensable edges, by uncovering significant information contributions from
#' indirect paths, and assesses edge-specific confidences from randomization of
#' available data. The remaining edges are then oriented based on the signature
#' of causality in observational data.
#'
#' The method relies on an information theoretic based (conditional) independence
#' test which is described in (Verny \emph{et al.}, PLoS Comp. Bio. 2017),
#' (Cabeli \emph{et al.}, PLoS Comp. Bio. 2020). It deals with both categorical
#' and continuous variables by performing optimal context-dependent discretization.
#' As such, the input data frame may contain both numerical columns which will be
#' treated as continuous, or character / factor columns which will be treated
#' as categorical. For further details on the optimal discretization method and
#' the conditional independence test, see the function discretizeMutual.
#' The user may also choose to run miic with scheme presented in
#' (Li \emph{et al.}, NeurIPS 2019) to improve the end result's interpretability
#' by ensuring consistent separating set during the skeleton iterations.
#'
#' @seealso \code{\link{discretizeMutual}} for optimal discretization and
#' (conditional) independence test.
#'
#' @references
#' \itemize{
#' \item Verny et al., \emph{PLoS Comp. Bio. 2017.}  https://doi.org/10.1371/journal.pcbi.1005662
#' \item Cabeli et al., \emph{PLoS Comp. Bio. 2020.}  https://doi.org/10.1371/journal.pcbi.1007866
#' \item Li et al., \emph{NeurIPS 2019} http://papers.nips.cc/paper/9573-constraint-based-causal-structure-learning-with-consistent-separating-sets.pdf
#' }
#'
#' @param input_data [a data frame]
#' A n*d data frame (n samples, d variables) that contains the observational data.
#' Each column corresponds to one variable and each row is a sample that gives the
#' values for all the observed variables. The column names correspond to the
#' names of the observed variables. Numeric columns will be treated as continuous
#' values, factors and character as categorical.
#'
#' @param black_box [a data frame]
#' An optional E*2 data frame containing E pairs of variables that will be considered
#' as independent during the network reconstruction. In practice, these edges will not
#' be included in the skeleton initialization and cannot be part of the final result.
#' Variable names must correspond to the \emph{input_data} data frame.
#'
#' @param n_eff [a positive integer]
#' The n samples given in the \emph{input_data} data frame are expected
#' to be independent. In case of correlated samples such as in time series or
#' Monte Carlo sampling approaches, the effective number of independent samples
#' \emph{n_eff} can be estimated using the decay of the autocorrelation function
#' (Verny \emph{et al.}, PLoS Comp. Bio. 2017). This \emph{effective} number \emph{n_eff}
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
#' @param latent [a string; \emph{c("no", "yes", "orientation")}]
#' When set to "yes", the network reconstruction is taking into account hidden (latent)
#' variables. When set to "orientation", latent variables are not considered during the skeleton
#' reconstruction but allows bi-directed edges during the orientation. Dependence
#' between two observed variables due to a latent variable is indicated with a '6' in
#' the adjacency matrix and in the network edges.summary and by a bi-directed edge
#' in the (partially) oriented graph.
#'
#' @param orientation [a boolean value]
#' The miic network skeleton can be partially directed
#' by orienting and propagating edge directions, based on the sign and magnitude
#' of the conditional 3-point information of unshielded triples. The propagation
#' procedure relyes on probabilities; for more details, see Verny \emph{et al.}, PLoS Comp. Bio. 2017).
#' If set to FALSE the orientation step is not performed.
#'
#' @param ori_proba_ratio [a floating point between 0 and 1] The threshold when
#' deducing the type of an edge tip (head/tail) from the probability of
#' orientation. For a given edge tip, denote by p the probability of it being a
#' head, the orientation is accepted if (1 - p) / p < ori_proba_ratio. 0 means
#' reject all orientations, 1 means accept all orientations.
#'
#' @param ori_consensus_ratio [a floating point between 0 and 1] The threshold
#' when deducing the type of an consensus edge tip (head/tail) from the average
#' probability of orientation. For a given edge tip, denote by p the probability
#' of it being a head, the orientation is accepted if (1 - p) / p <
#' ori_consensus_ratio. 0 means reject all orientations, 1 means accept all
#' orientations.
#'
#' @param propagation [a boolean value]
#' If set to FALSE, the skeleton is partially oriented with only the
#' v-structure orientations. Otherwise, the v-structure orientations are
#' propagated to downstream undirected edges in unshielded triples following
#' the orientation method
#'
#' @param state_order [a data frame] An optional data frame providing extra
#' information for variables. It must have d rows where d is the number of input
#' variables, and the following structure (named columns):
#'
#' "var_names" (required) contains the name of each variable as specified
#' by colnames(input_data).
#'
#' "var_type" (optional) contains a binary value that specifies if each
#' variable is to be considered as discrete (0) or continuous (1).
#'
#' "levels_increasing_order" (optional) contains a single character string
#' with all of the unique levels of the ordinal variable in increasing order,
#' delimited by comma ','. It will be used during the post-processing to compute
#' the sign of an edge using Spearman's rank correlation. If a variable is
#' continuous or is categorical but not ordinal, this column should be NA.
#'
#' "is_contextual" (optional) contains a binary value that specifies if a
#' variable is to be considered as a contextual variable (1) or not (0).
#' Contextual variables cannot be the child node of any other variable (cannot
#' have edge with arrowhead pointing to them).
#'
#' @param true_edges [a data frame]
#' An optional E*2 data frame containing the E edges of the true graph for
#' computing performance after the run.
#'
#' @param n_shuffles [a positive integer] The number of shufflings of
#' the original dataset in order to evaluate the edge specific confidence
#' of all retained edges.
#'
#' @param conf_threshold [a positive floating point] The threshold used
#' to filter the less probable edges following the skeleton step. See Verny
#' \emph{et al.}, PLoS Comp. Bio. 2017.
#'
#' @param sample_weights [a numeric vector]
#' An optional vector containing the weight of each observation.
#'
#' @param test_mar [a boolean value]
#' If set to TRUE, distributions with missing values will be tested with Kullback-Leibler
#' divergence : conditioning variables for the given link \eqn{X\rightarrow Y}\eqn{Z} will be
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
#' the network. See (Li \emph{et al.}, NeurIPS 2019) for details.
#'
#' @param max_iteration [a positive integer] When the \emph{consistent} parameter
#' is set to "skeleton" or "orientation", the maximum number of iterations
#' allowed when trying to find a consistent graph. Set to 100 by default.
#'
#' @param consensus_threshold [a floating point between 0.5 and 1.0] When the
#' \emph{consistent} parameter is set to "skeleton" or "orientation", and when
#' the result graph is inconsistent, or is a union of more than one inconsistent
#' graphs, a consensus graph will be produced based on a pool of graphs. If the
#' result graph is inconsistent, then the pool is made of [max_iteration] graphs
#' from the iterations, otherwise it is made of those graphs in the union. In
#' the consensus graph, an edge is present when the proportion of non-zero
#' status in the pool is above the threshold. For example, if the pool contains
#' [A, B, B, 0, 0], where "A", "B" are different status of the edge and "0"
#' indicates the absence of the edge. Then the edge is set to connected ("1") if
#' the proportion of non-zero status (0.6 in the example) is equal to or higher
#' than [consensus_threshold]. (When set to connected, the orientation of the
#' edge will be further determined by the average probability of orientation.)
#' Set to 0.8 by default.
#'
#' @param verbose [a boolean value] If TRUE, debugging output is printed.
#'
#' @param n_threads [a positive integer] When set greater than 1, n_threads
#' parallel threads will be used for computation. Make sure your compiler is
#' compatible with openmp if you wish to use multithreading.
#'
#' @param negative_info [a boolean value] For test purpose only. FALSE by
#' default. If TRUE, negative shifted mutual information is allowed during the
#' computation when mutual information is inferior to the complexity term. For
#' small dateset with complicated structures, e.g., discrete variables with many
#' levels, allowing for negative shifted mutual information may help identifying
#' weak v-structures related to those discrete variables, as the negative
#' three-point information in those cases will come from the difference between
#' two negative shifted mutual information terms (expected to be negative due to
#' the small sample size). However, under this setting, a v-structure (X -> Z <-
#' Y) in the final graph does not necessarily imply that X is dependent on Y
#' conditioning on Z, As a consequence, the interpretability of the final graph
#' is hindered. In practice, it's advised to keep this parameter as FALSE.
#'
#' @return A \emph{miic-like} object that contains:
#' \itemize{
#'  \item{all.edges.summary:} {a data frame with information about the
#'  relationship between each pair of variables
#'  \itemize{
#'  \item \emph{x:} X node name
#'  \item \emph{y:} Y node name
#'  \item \emph{type:} 'N': the edge is removed,'P': the edge is retained.
#'  If the true_edges file is given,
#'  'N' is replaced by 'TN' (True Negative) or 'FN' (False Negative),
#'  'P' is replaced by 'TP' (True Positive) or 'FP' (False Positive).
#'  \item \emph{ai:} the set of conditioning nodes used in the attempt to
#'  separate \emph{x} and \emph{y}.
#'  \item \emph{info:} the mutual information \emph{I(x;y)} times the number of
#'  complete samples \emph{Nxy} used in the computation.
#'  \item \emph{info_cond:} the conditional mutual information \emph{I(x;y|ai)}
#'  times the number of complete samples \emph{Nxy_ai} used in the computation,
#'  equal to \emph{info} when \emph{ai} is an empty set.
#'  \item \emph{cplx:} the complexity term for the pair (\emph{x, y})
#'  taking into account the separating set \emph{ai}. The edge between \emph{x}
#'  and \emph{y} is retained if \emph{info_cond} is greater than \emph{cplx}.
#'  \item \emph{Nxy_ai:} the number of complete samples taking into account
#'  \emph{x}, \emph{y} and all nodes in \emph{ai}, based on which the
#'  information and the complexity terms are computed. Without missing value, it
#'  is the same for all pairs and is equal to the total number of samples.
#'  \item \emph{info_shifted:} value equal to \emph{info} - \emph{cplx}.
#'  Used to decide whether the edge is removed (positive) or retained (zero,
#'  possibly negative when the parameter \emph{negative_info} is set to TRUE).
#'  \item \emph{ort_inferred:} the inferred orientation of edge (\emph{x, y}).
#'  0: edge removed, 1: undirected, 2: directed from X to Y, -2: directed from Y
#'  to X, 6: bidirected. When the \emph{consistent} option is turned on and
#'  there is more than one graph in the consistent cycle, this is the inferred
#'  ororientation of the edge in the last graph in the cycle.
#'  \item \emph{ort_ground_truth:} the orientation of the edge (\emph{x, y}) as
#'  specified in the true_edges file (if provided).
#'  \item \emph{is_inference_correct:} TRUE: the inferred orientation agrees
#'  with the provided ground truth, FALSE: the inferred orientation disagrees
#'  with the provided ground truth.
#'  \item \emph{is_causal:} boolean value indicating the causal nature of the
#'  arrow tips of a directed edge, measured based on the probabilities given in
#'  the columns \emph{p_y2x} and \emph{p_x2y}.
#'  TRUE: both the head and the tail are set with high probability (threshold
#'  set by the parameter \emph{ori_proba_ratio}), implying a likely causation.
#'  FALSE: only the head is set with high probability, while the tail
#'  probability is not extreme enough to be either head or tail. In this case,
#'  the edge is undecided between a directed edge (implying causation) and a
#'  bidirected edge (implying possible latent confounders).
#'  \item \emph{ort_consensus:} when the \emph{consistent} option is turned on
#'  and there is more than one graph in the consistent cycle, the consensus of
#'  all graphs in the cycle on the orientation of the edge (\emph{x, y}).
#'  \item \emph{is_causal_consensus:} same as \emph{is_causal} but for
#'  \emph{ort_consensus}.
#'  \item \emph{edge_stats} when the \emph{consistent} option is turned on and
#'  there is more than one graph in the consistent cycle, a histogram of all
#'  \emph{ort_inferred} values present in the cycle for the edge (\emph{x, y}),
#'  in the format [percentage(orientation)], separated by ";".
#'  \item \emph{sign:} the sign of the partial correlation between variables
#'  \emph{x} and \emph{y}, conditioned on the contributing nodes \emph{ai}.
#'  \item \emph{partial_correlation:} value of the partial correlation for the
#'  edge (\emph{x, y}) conditioned on the contributing nodes \emph{ai}.
#'  \item \emph{p_y2x:} probability of the arrowhead from \emph{y} to \emph{x}, of
#'  the inferred orientation, derived from the three-point mutual information
#'  (cf Affeldt & Isambert, UAI 2015 proceedings).
#'  \item \emph{p_x2y:} probability of the arrowhead from \emph{x} to \emph{y}, of
#'  the inferred orientation, derived from the three-point mutual information
#'  (cf Affeldt & Isambert, UAI 2015 proceedings).
#'  \item \emph{confidence:} when the parameter \emph{n_shuffles} is positive,
#'  the ratio between the term \emph{exp(-info_shifted(x;y|ai))} of the original
#'  dataset and that of the randomized dataset (averaged over \emph{n_shuffles}
#'  randomizations, as a measure of the strength of the retained edge.
#'  }
#'  }
#'
#'  \item{orientations.prob:} {a data frame of the orientation probabilities of
#'  the two edges of all unshielded triples (node1 -- mid-node -- node2) of the
#'  reconstructed network:}
#'  \itemize{
#'  \item node1: left node of the unshielded triple
#'  \item p1: probability of the arrowhead node1 <- mid-node
#'  \item p2: probability of the arrowhead node1 -> mid-node
#'  \item mid-node: middle node of the unshielded triple
#'  \item p3: probability of the arrowhead mid-node <- node2
#'  \item p4: probability of the arrowhead mid-node -> node2
#'  \item node2: right node of the unshielded triple
#'  \item NI3: 3-point information * \emph{N} where \emph{N} is the number of
#'  complete samples taking into account node1, mid-node, node2 and \emph{ai} of
#'  the pair (node1, node2).
#'  }
#'
#'
#'  \item {adj_matrix:} the adjacency matrix of the inferred graph.
#'  For an element at position (row, column), the row is regarded as the source
#'  node \emph{x} and the column the target node \emph{y}, with the value
#'  representing the status of the tip of edge from \emph{x} to \emph{y}:
#'  \itemize{
#'  \item 0: the edge is removed
#'  \item 1: the tip is a tail \emph{x} -- \emph{y}
#'  \item 2: the tip is a head \emph{x} -> \emph{y}
#'  }
#' }
#' @export
#' @useDynLib miic
#' @import Rcpp
#'
#' @examples
#' library(miic)
#'
#' # EXAMPLE HEMATOPOIESIS
#' data(hematoData)
#'
#' # execute MIIC (reconstruct graph)
#' miic.res <- miic(
#'   input_data = hematoData[1:1000,], latent = "yes",
#'   n_shuffles = 10, conf_threshold = 0.001
#' )
#'
#' # plot graph
#' if(require(igraph)) {
#'  plot(miic.res, method="igraph")
#' }
#'
#' \donttest{
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
#'   input_data = cosmicCancer, state_order = cosmicCancer_stateOrder, latent = "yes",
#'   n_shuffles = 100, conf_threshold = 0.001
#' )
#'
#' # plot graph
#' if(require(igraph)) {
#'  plot(miic.res)
#' }
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
#'   input_data = ohno, latent = "yes", state_order = ohno_stateOrder,
#'   n_shuffles = 100, conf_threshold = 0.001
#' )
#'
#' # plot graph
#' if(require(igraph)) {
#'  plot(miic.res)
#' }
#'
#' # write graph to graphml format. Note that to correctly visualize
#' # the network we created the miic style for Cytoscape (http://www.cytoscape.org/).
#' miic.write.network.cytoscape(g = miic.res, file = file.path(tempdir(), "temp"))
#' }
#'
miic <- function(input_data,
                 state_order = NULL,
                 true_edges = NULL,
                 black_box = NULL,
                 n_threads = 1,
                 cplx = c("nml", "mdl"),
                 orientation = TRUE,
                 ori_proba_ratio = 1,
                 ori_consensus_ratio = NULL,
                 propagation = TRUE,
                 latent = c("no", "yes", "orientation"),
                 n_eff = -1,
                 n_shuffles = 0,
                 conf_threshold = 0,
                 sample_weights = NULL,
                 test_mar = TRUE,
                 consistent = c("no", "orientation", "skeleton"),
                 max_iteration = 100,
                 consensus_threshold = 0.8,
                 negative_info = FALSE,
                 verbose = FALSE) {
  res <- NULL

  if (is.null(input_data)) {
    stop("The input data file is required")
  }

  if (!is.data.frame(input_data)) {
    stop("The input data is not a dataframe")
  }
  # Remove rows with only NAs
  input_data <- input_data[rowSums(is.na(input_data)) != ncol(input_data), ]
  if (length(input_data) == 0) {
    stop("The input data is empty or contains only NAs")
  }

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


  if (n_eff > nrow(input_data)) {
    stop(
      paste0(
        "The number of effective samples cannot be greater than the ",
        "number of samples."
      )
    )
  }

  if (propagation != TRUE && propagation != FALSE) {
    stop("The propagation type is not correct. Allowed types are TRUE or FALSE")
  }

  if (orientation != TRUE && orientation != FALSE) {
    stop("The orientation type is not correct. Allowed types are TRUE or FALSE")
  }

  if (is.null(ori_consensus_ratio)) {
    ori_consensus_ratio <- ori_proba_ratio
  }

  if (verbose) {
    cat("START miic...\n")
  }

  is_contextual <- rep(0, ncol(input_data))
  names(is_contextual) <- colnames(input_data)
  is_continuous <- sapply(input_data, is.numeric)
  # Parse "state_order" file
  if (!is.null(state_order)) {
    err_code <- checkStateOrder(state_order, input_data)
    if (err_code != "0") {
      warning(paste(errorCodeToString(err_code),
          "state_order file will be ignored.", sep=", "), call.=FALSE)
      state_order <- NULL
    } else if (is.null(state_order$var_names)) {
      warning(paste("Named column var_names is required for state_order,",
                    "state_order file will be ignored."), call.=FALSE)
      state_order <- NULL
    } else {
      mismatch <- is.na(match(state_order$var_names, colnames(input_data)))
      if (any(mismatch)) {
        var_str <- paste(state_order$var_names[mismatch], collapse=", ")
        warning(paste("Variable(s)", var_str, "specified in state_order file",
            "will not match any name in input_data, and will be ignored."),
            call.=FALSE)
        state_order <- state_order[!mismatch, ]
      }
      not_found <- is.na(match(colnames(input_data), state_order$var_names))
      if (any(not_found)) {
        var_str <- paste(colnames(input_data)[not_found], collapse=", ")
        warning(paste("Variable(s)", var_str, "in input_data not found",
            "in state_order file."), call.=FALSE)
      }
      for (row in 1:nrow(state_order)) {
        col <- as.character(state_order[row, "var_names"])
        if (!is.null(state_order$var_type)) {
          if (state_order[row, "var_type"] == 0) {
            input_data[, col] <- factor(input_data[, col])
            is_continuous[[col]] <- FALSE
          }
        }
        if (!is.null(state_order$is_contextual)) {
          if (state_order[row, "is_contextual"] == 1) {
            is_contextual[[col]] <- 1
          }
        }
        if (!is.null(state_order$levels_increasing_order)) {
          order_string <- state_order[row, "levels_increasing_order"]
          if (!is.na(order_string)) {
            if (is_continuous[[col]] == TRUE) {
              warning(paste(col, "is considered as a continuous variable,",
                  "the provided orders will be ignored."), call.=FALSE)
            } else {
              orders <- unlist(strsplit(as.character(order_string), ","))
              values <- unique(input_data[[col]][!is.na(input_data[[col]])])
              absent <- is.na(match(values, orders))
              if (any(absent)) {
                var_str <- ""
                if (length(values[absent]) > 10) {
                  # Only show the first 3 if the list is too long
                  var_str <- paste(paste(values[absent][1:3], collapse=", "),
                      "and", length(values[absent]) - 3, "more")
                } else {
                  var_str <- paste(values[absent], collapse=", ")
                }
                warning(paste("Variable", col, "has value(s)", var_str,
                    "that will not match the provided orders,",
                    "the provided orders will be ignored."), call.=FALSE)
              }
            }
          }
        }
      }
    }
  }
  # Check the number of unique values of continuous and discrete variables
  for (col in colnames(input_data)) {
    unique_values <- length(unique(input_data[[col]][!is.na(input_data[[col]])]))
    if (is_continuous[[col]] &&
      (unique_values <= 40) && (nrow(input_data) > 40)) {
      if (is_continuous[[col]] &&
        (unique_values <= 2) && (nrow(input_data) > 40)) {
        # Less than 3 unique values does not make sense for a continuous variable
        stop(
          paste0(
            "Numerical variable ",
            col,
            " only has ",
            unique_values,
            " non-NA unique values. Is this a factor?"
          )
        )
      }
      # Less than 40 unique variables can be discretized but may not be truly
      # continuous
      warning(
        paste0(
          "Numerical variable ",
          col,
          " is treated as continuous but only has ",
          unique_values,
          " unique values."
        ), call.=FALSE
      )
    }
    if ((!is_continuous[[col]]) &&
      (unique_values >= 40) && (nrow(input_data) > 40)) {
      warning(paste0(
        col,
        " is treated as discrete but has many levels (",
        unique_values,
        ")."
      ), call.=FALSE)
    }
  }

  err_code <- checkInput(input_data, "miic")
  if (err_code != "0") {
    warning(errorCodeToString(err_code), call.=FALSE)
  } else {
    if (verbose) {
      cat("\t# -> START reconstruction...\n")
    }
    res <-
      miic.reconstruct(
        input_data = input_data,
        n_threads = n_threads,
        cplx = cplx,
        latent = latent,
        n_eff = n_eff,
        black_box = black_box,
        n_shuffles = n_shuffles,
        orientation = orientation,
        ori_proba_ratio = ori_proba_ratio,
        propagation = propagation,
        conf_threshold = conf_threshold,
        verbose = verbose,
        is_contextual = is_contextual,
        is_continuous = is_continuous,
        sample_weights = sample_weights,
        test_mar = test_mar,
        consistent = consistent,
        max_iteration = max_iteration,
        negative_info = negative_info
      )
    if (res$interrupted) {
      stop("Interupted by user")
    }
    if (verbose) {
      cat("\t# -> END reconstruction...\n\t# --------\n")
    }

    if (!is.null(true_edges)) {
      err_code <- checkTrueEdges(true_edges)
      if (err_code != "0") {
        warning(paste(errorCodeToString(err_code),
            "true_edges file will be ignored.", sep=", "), call.=FALSE)
        true_edges <- NULL
      }
    }

    res$all.edges.summary <- summarizeResults(
      observations = input_data,
      results = res,
      true_edges = true_edges,
      state_order = state_order,
      consensus_threshold = consensus_threshold,
      ori_consensus_ratio = ori_consensus_ratio,
      latent = latent != "no",
      propagation = propagation,
      verbose = verbose
    )
  }

  class(res) <- "miic"
  return(res)
}


#' Basic plot function of a miic network inference result
#'
#' @description This function calls \code{\link{miic.export}} to build a
#' plottable object from the result returned by \code{\link{miic}} and plot it.
#'
#' @details See the documentation of \code{\link{miic.export}} for further
#' details.
#'
#' @param x [a miic graph object]
#' The graph object returned by \code{\link{miic}}.
#' @param method A string representing the plotting method. Default to "igraph".
#' Currently only "igraph" is supported.
#' @param \dots Additional plotting parameters. See the corresponding plot function
#' for the complete list.
#' For igraph, see \code{\link[igraph]{igraph.plotting}}.
#'
#' @export
#'
#' @seealso \code{\link{miic.export}} for generic exports,
#' \code{\link{getIgraph}} for igraph export,
#' \code{\link[igraph]{igraph.plotting}}
#'
plot.miic = function(x, method = 'igraph', ...) {
  if (class(x) != "miic"){
    stop("Not a miic object.")
  }
  if (method == 'igraph'){
    if (base::requireNamespace("igraph", quietly = TRUE)) {
      igraph::plot.igraph(miic.export(x, 'igraph'), ...)
    } else {
      stop("Package 'igraph' is required.")
    }
  } else {
    stop("Method not supported. See ?miic.export for supported methods.")
  }
}
