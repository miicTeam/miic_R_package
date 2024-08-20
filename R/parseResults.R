#-------------------------------------------------------------------------------
# summarizeResults
#-------------------------------------------------------------------------------
# Summarize the list of edges, the list will contain:
# - edges that exist in the miic reconstruction (oriented or not)
# - edges that were conditioned away with a non-empty separating set
# - if ground truth is known, edges present in ground truth but not in the
#   2 previous categories, it corresponds to the true edges removed without
#   conditioning
# The summary is sorted by info_shifted (log likelihood), decreasing.
#
# Inputs:
# - observations: a data frame, mandatory, the input data
# - results: the miic c++ part output, mandatory.
# - true_edges: a 2 columns data frame, optional, NULL by default.
#   The ground truth, if known
# - state_order: the state order data frame used, optional, NULL by default.
# - consensus_threshold: a float, optional, 0.8 by default. Used when
#   consistency is activated to construct the consensus graph skeleton.
# - ori_consensus_ratio: a float, optional, 0.1 by default. Used when
#   consistency is activated for the consensus graph orientation
# - latent: a boolean, optional, TRUE by default. Indicates if latent
#   variables discovery was activated during the network reconstruction.
# - propagation: a boolean, optional, FALSE by default. Indicates if orientation
#   propagation was activated during the network reconstruction.
# Return:
# - a data frame: the summary data frame.
#   * x: 1st node (in alphanumerical order) of the edge
#   * y: 2nd node (in alphanumerical order) of the edge
#   * type: the miic prediction : "P"(ositive) or "N"(egative) for respectively
#     presence or absence of an edge, without considering orientation.
#     If ground truth is known, edges are classified as True or False
#     Positives/Negatives ("TP", "FP", "TN", "FN").
#   * ai: list containing the conditioning nodes
#     (corresponds to ai.vect in the C++ output)
#   * raw_contributions: raw contributions of each ai to the conditional
#     independence, measured by I'(x;y;ai|{aj}) / I'(x;y),
#     where {aj} is the separating set before adding ai.
#   * contributions : contributions of each ai to the reduction of conditional
#     mutual information measured by I'(x;y;ai|{aj}) / I'(x;y|{aj}),
#     where {aj} is the separating set before adding ai.
#   * info: mutual information (corresponds to Ixy in the C++ output)
#   * info_cond: conditional mutual information  (corresponds to Ixy_ai in the
#     C++ output)
#   * cplx: the NML complexity (used for independence testing)
#   * Nxy_ai: the number of samples without missing values used for this edge
#   * info_shifted: the difference between conditional MI and cplx
#   * infOrt: the inferred edge orientation
#   * trueOrt: is the true edge orientation if known. NA if truth is unknown
#   * isOrtOk: indicate if the inferred edge is correctly oriented.
#     NA if truth is unknown, "Y" or "N" if truth is known
#   * sign: sign of partial correlation between x and y conditioned on "ai"s
#   * partial_correlation: coefficient of partial correlation between x and y
#     conditioned on "ai"s
#   * is_causal: indicates if edge is genuine causal.
#     Note that the genuine causality is deducible only when latent variables
#     are allowed and propagation is not allowed
#   * proba: contains the orientation likelihoods as computed by miic
#     (cf Affeldt & Isambert, UAI 2015 proceedings) : the probabilities of
#     both orientations separated by a semi-colon.
#   * confidence: the ratio of info_shifted between randomized
#     and normal samples
#-------------------------------------------------------------------------------
summarizeResults = function (observations, results,
                             true_edges = NULL, state_order = NULL,
                             consensus_threshold = 0.8, ori_consensus_ratio = 0.1,
                             latent = TRUE, propagation = FALSE)
  {
  # Keep only edges remaining and edges removed using conditioning
  #
  summary <- results$edges [ (results$edges$category == 1)
                           | (! is.na (results$edges$contributions) ), , drop=F]
  #
  # If true edges is supplied, check, and add if needed, the False Negative
  # (the edges removed without conditioning are not present in the summary
  # but need to be there to be marked as FN)
  #
  if ( ! is.null(true_edges) )
    {
    # To match easily and quickly edges between data frames, define for each edge
    # the nodes names ordered alphanumerically as rowname
    #
    rownames(summary) <- apply ( summary, 1, FUN=function (x) {
      ifelse (x[["x"]] < x[["y"]],
              paste0 (x[["x"]], "-", x[["y"]]),
              paste0 (x[["y"]], "-", x[["x"]])) } )
    rownames(true_edges) <- apply ( true_edges, 1, FUN=function (x) {
      ifelse (x[[1]] < x[[2]],
              paste0 (x[[1]], "-", x[[2]]),
              paste0 (x[[2]], "-", x[[1]])) } )
    #
    # Detect missing edges present in ground truth but not in summary
    #
    missing_fn <- rownames(true_edges)[
      ! ( rownames(true_edges) %in% rownames(summary) ) ]
    if (length (missing_fn) > 0)
      {
      # Add missing edges coming from the ground truth
      #
      summary[ ( nrow (summary) + 1 ) :
               ( nrow (summary) + length (missing_fn) ), ] = NA
      rownames (summary) [ (nrow (summary) - length (missing_fn) + 1) :
                           nrow(summary) ] = missing_fn
      summary[ missing_fn, c("x","y") ] = true_edges [ missing_fn, ]
      #
      # Pick values from 'edges' data frame for these missing FN edges
      # To avoid iterations over a possibly huge 'edges' data frame
      # apply a pre-filtering to reduce search
      # (normally, there should not be a lot of missing FN)
      #
      pre_filter = unique (unlist (true_edges [ missing_fn, ]) )
      edges_4_fn = results$edges[ (results$edges$x %in% pre_filter)
                                | (results$edges$y %in% pre_filter), , drop=F]
      cols_2_pick = colnames (results$edges)
      cols_2_pick = cols_2_pick[ (cols_2_pick != "x") & (cols_2_pick != "y") ]
      for (i in (nrow (summary) - length (missing_fn) + 1) : nrow(summary) )
        {
        one_edge = edges_4_fn[ ( (edges_4_fn$x == summary[i, "x"])
                               & (edges_4_fn$y == summary[i, "y"]) )
                             | ( (edges_4_fn$x == summary[i, "y"])
                               & (edges_4_fn$y == summary[i, "x"]) ), , drop=F]
        summary[i, cols_2_pick] = one_edge[1, cols_2_pick]
        }
      }
    }
  #
  # If no edge in the summary, returns directly an empty data frame
  #
  if (nrow (summary) == 0)
    return (data._frame (x = character(0), y = character(0),
                         type = character(0), ai = character(0),
                         raw_contributions = character(0),
                         contributions = character(0),
                         info = numeric(0), info_cond = numeric(0),
                         cplx = numeric(0), Nxy_ai = numeric(0),
                         info_shifted = numeric(0), infOrt = numeric(0),
                         trueOrt = numeric(0), isOrtOk = character(0),
                         sign = character(0), partial_correlation = numeric(0),
                         is_causal = character(0), proba = character(0),
                         confidence = character(0), stringsAsFactors = FALSE) )
  #
  # Edge ordering (A<-B or B->A) is given by alphanumerical order
  #
  summary [, c("x","y")] <- t (apply (summary[,c("x","y")], 1,
                                      function (row) { sort(row) } ) )
  #
  # Edge 'type' corresponds to the miic prediction : P(ositive) or N(egative)
  # for respectively presence or absence of an edge, without considering
  # orientation. If ground truth is known, edges are classified as True or
  # False Positives/Negatives (TP, FP, TN, FN).
  #
  if ( is.null(true_edges) )
    summary$type <- ifelse (summary$category == 1, "P", "N")
  else
    {
    summary$type[ (summary$category == 1)
                & (rownames(summary) %in% rownames(true_edges)) ] <- "TP"
    summary$type[ (summary$category == 1)
                & (! (rownames(summary) %in% rownames(true_edges))) ] <- "FP"
    summary$type[ (summary$category == 0)
                & (! (rownames(summary) %in% rownames(true_edges))) ] <- "TN"
    summary$type[ (summary$category == 0)
                & (rownames(summary) %in% rownames(true_edges)) ] <- "FN"
    }
  #
  # ai is a list containing the conditioning nodes
  #
  colnames(summary)[ which (colnames(summary) == "ai.vect") ] <- "ai"
  #
  # info and info_cond contain the (conditional) mutual information values
  #
  colnames(summary)[ which (colnames(summary) == "Ixy")] <- "info"
  colnames(summary)[ which (colnames(summary) == "Ixy_ai")] <- "info_cond"
  #
  # info_shifted is the difference between MI and cplx
  #
  summary$info_shifted <- summary$info_cond - summary$cplx
  #
  # confidence is the ratio of info_shifted between randomized
  # and normal samples
  #
  summary$confidence [summary$confidence == -1] <- NA
  #
  # infOrt is the inferred edge orientation
  #
  adj_matrix <- results$adj_matrix
  summary$infOrt <- apply (summary, 1, function (row, adj_matrix) {
                              adj_matrix[row[1], row[2]] }, adj_matrix)
  #
  # trueOrt is the true edge orientation (if known)
  #
  var_names <- colnames(adj_matrix)
  if ( is.null(true_edges) )
    summary$trueOrt <- NA
  else
    {
    true_adj_matrix <- matrix (0, ncol = ncol(adj_matrix),
                               nrow = nrow(adj_matrix),
                               dimnames = list(var_names,var_names) )
    for ( i in 1:nrow (true_edges) )
      {
      true_edge <- as.character (unlist (true_edges[i, ]) )
      true_adj_matrix[true_edge[1], true_edge[2]] <- 2
      true_adj_matrix[true_edge[2], true_edge[1]] <- -2
      }
    summary$trueOrt <- apply (summary, 1, function (row, true_adj_matrix) {
      true_adj_matrix[row[1], row[2]] }, true_adj_matrix)
    }
  #
  # isOrtOk indicates if the inferred edge is correctly oriented.
  # NA if truth is unknown, "Y" or "N" if truth is known
  #
  summary$isOrtOk <- ifelse (summary$infOrt == summary$trueOrt, "Y", "N")
  #
  # Sign and coefficient of partial correlation between x and y conditioned
  # on "ai"s.
  #
  summary [, c("sign", "partial_correlation")] <- compute_partial_correlation (
    summary, observations, state_order)
  summary$partial_correlation <- as.numeric (summary$partial_correlation)
  #
  # Probabilities of both orientations separated by a semi-colon.
  #
  orientation_probabilities <- results$orientations.prob
  summary$proba <- sapply (1:nrow(summary), function(i)
    {
    row <- summary[i, ]
    id_x <- match(row$x, var_names)
    id_y <- match(row$y, var_names)
    proba_adj <- results$proba_adj_matrix
    if (!is.null(results$adj_matrices) && ncol(results$adj_matrices) > 1)
      {
      proba_adj <- results$proba_adj_average
      }
    return (paste (proba_adj[id_y, id_x], proba_adj[id_x, id_y], sep = ";") )
    })
  #
  # Genuine causality is deducible only when latent variables are allowed and
  # propagation is not allowed
  #
  causality_deducible <- latent && (!propagation)
  #
  # If consistent parameter is turned on and the result graph is a union of
  # more than one inconsistent graphs, get the possible orientations of each
  # edge with the corresponding frequencies and the consensus status.
  #
  summary$is_causal = NA
  if ( ( ! is.null (results$adj_matrices) )
     && (ncol(results$adj_matrices) > 1) )
    {
    # use split to turn summary data frame to list only to be able to use sapply
    # with simplify = FALSE, otherwise apply() will force simplification, which
    # is annoying and surprising when the return value of the function is the
    # same for all rows. Theoretically it can never happen here, but this is R!
    # Always be careful.
    #
    edge_stats_table <- sapply (split (summary, seq (nrow (summary) ) ),
                                get_edge_stats_table,
                                var_names,
                                results$adj_matrices,
                                simplify = FALSE)
    target <- which(names(summary) == "infOrt")[1]
    summary <- cbind(
      summary[, 1:target, drop = FALSE],
      consensus = sapply(edge_stats_table,
        get_consensus_status,
        consensus_threshold
      ),
      edge_stats = sapply(edge_stats_table, get_edge_stats_str),
      is_causal_consensus = NA,
      summary[, (target + 1):length(summary), drop = FALSE]
    )
    #
    # Set consensus edge status according to the average probabilities
    #
    for (i in 1:nrow(summary))
      {
      row <- summary[i, ]
      if (causality_deducible)
        {
        # Set initial values if deducible
        if (row$infOrt != 0)
          summary[i, ]$is_causal <- "N"
        if (row$consensus != 0)
          summary[i, ]$is_causal_consensus <- "N"
        }
      if (row$consensus == 0) next

      id_x <- match(row$x, var_names)
      id_y <- match(row$y, var_names)
      # probability of an edge tip being a head (<,>), * means head or tail (-)
      proba_x2y <- results$proba_adj_average[id_x, id_y]  # proba of x *-> y
      proba_y2x <- results$proba_adj_average[id_y, id_x]  # proba of x <-* y
      ratio_x2y <- (1 - proba_x2y) / proba_x2y
      ratio_y2x <- (1 - proba_y2x) / proba_y2x

      if (  (ratio_x2y < ori_consensus_ratio)
         && (ratio_y2x < ori_consensus_ratio) )
        {
        summary[i, ]$consensus <- 6
        }
      else if (  (ratio_x2y < ori_consensus_ratio)
              && (ratio_y2x >= ori_consensus_ratio) )
        {
        summary[i, ]$consensus <- 2
        if (1 / ratio_y2x < ori_consensus_ratio && causality_deducible)
          {
          summary[i, ]$is_causal_consensus <- "Y"
          if (row$infOrt == 2)
            {
            summary[i, ]$is_causal <- "Y"
            }
          }
        }
      else if (ratio_y2x < ori_consensus_ratio &&
                 ratio_x2y >= ori_consensus_ratio)
        {
        summary[i, ]$consensus <- -2
        if (1 / ratio_x2y < ori_consensus_ratio && causality_deducible)
          {
          summary[i, ]$is_causal_consensus <- "Y"
          if (row$infOrt == -2)
            {
            summary[i, ]$is_causal <- "Y"
            }
          }
        }
      else
        {
        summary[i, ]$consensus <- 1
        }
      }
    }
  else if (causality_deducible)
    {
    # if (is.null(results$adj_matrices) || ncol(results$adj_matrices) <= 1)
    # set is_causal by results$proba_adj_matrix
    for (i in 1:nrow(summary))
      {
      row <- summary[i, ]
      if (row$infOrt == 0)
        next
      summary[i, ]$is_causal <- "N"

      id_x <- match(row$x, var_names)
      id_y <- match(row$y, var_names)
      # probability of an edge tip being a head (<,>), * means head or tail (-)
      proba_x2y <- results$proba_adj_matrix[id_x, id_y]  # proba of x *-> y
      proba_y2x <- results$proba_adj_matrix[id_y, id_x]  # proba of x <-* y
      ratio_x2y <- (1 - proba_x2y) / proba_x2y
      ratio_y2x <- (1 - proba_y2x) / proba_y2x
      if (   (row$infOrt == 2)
          && (ratio_x2y < ori_consensus_ratio)
          && (1 / ratio_y2x < ori_consensus_ratio) )
        {
        summary[i, ]$is_causal <- "Y"
        }
      if (  (row$infOrt == -2)
         && (ratio_y2x < ori_consensus_ratio)
         && (1 / ratio_x2y < ori_consensus_ratio) )
        {
        summary[i, ]$is_causal <- "Y"
        }
      }
    }
  #
  # Sort summary by log likelihood, keep only some cols and returns
  #
  columns_kept = c ("x", "y", "type", "ai", "raw_contributions",
                    "contributions", "info", "info_cond", "cplx",
                    "Nxy_ai", "info_shifted", "infOrt",
                    "consensus", "edge_stats", "is_causal_consensus",
                    "trueOrt", "isOrtOk", "sign", "partial_correlation",
                    "is_causal", "proba", "confidence")
  columns_kept = columns_kept[columns_kept %in% colnames(summary)]
  summary <- summary[order(summary$info_shifted, decreasing = TRUE),
                     columns_kept, drop=F]
  rownames(summary) <- c()
  return (summary)
  }

#===============================================================================
# TO_BE_DELETED AFTER TESTS
#===============================================================================
# summarizeResults_old
#-------------------------------------------------------------------------------
summarizeResults_old <- function(observations = NULL, results = NULL,
                             true_edges = NULL, state_order = NULL,
                             consensus_threshold = 0.8,
                             ori_consensus_ratio = 0.1, latent = TRUE,
                             propagation = FALSE, verbose = FALSE) {
  # Reduced list of edges that will be summarized. There are 3 categories:
  # - Edges that exist in the miic reconstruction (oriented or not)
  # - Edges that were conditioned away with a non-empty separating set
  # - If ground truth is known, any other positive edge
  summarized_edges <- matrix(ncol = 2, nrow = 0)
  adj_matrix <- results$adj_matrix
  var_names <- colnames(adj_matrix)

  # List of edges found by miic
  half_adj_matrix = adj_matrix
  half_adj_matrix[ lower.tri(adj_matrix, diag = TRUE) ] <- 0
  predicted_edges <- which(half_adj_matrix != 0, arr.ind = TRUE, useNames = FALSE)
  predicted_edges <- apply(predicted_edges, 2, function(x) { var_names[x] })
  # Add to summarized edges list
  if(length(predicted_edges > 0)) summarized_edges <- predicted_edges

  edges <- results$edges
  # List of negative edges with non null conditioning set
  conditioned_edges <- as.matrix(edges[(!is.na(edges$ai.vect)) &
    edges$category != 1, c("x", "y")])
  # List of negative edges with null conditioning set
  indep_null_cond_set_edges <- as.matrix(edges[(is.na(edges$ai.vect)) &
    edges$category != 1, c("x", "y")])
  # Add to summarized edges list
  summarized_edges <- rbind(summarized_edges, conditioned_edges)
                            #indep_null_cond_set_edges)

  if (!is.null(true_edges)) {
    # List of False Negative edges
    false_negative_edges <- data.frame(
      x = character(), y = character(),
      stringsAsFactors = FALSE
    )
    for (i in 1:nrow(true_edges)) {
      true_edge <- as.character(unlist(true_edges[i, ]))
      if (!any(apply(summarized_edges, 1, function(row, edge) {
        all(row %in% edge)
      }, true_edge))) {
        # If this edge is not already in summarized edges
        false_negative_edges[nrow(false_negative_edges) + 1, ] <- true_edge
      }
    }
    # Add to summarized edges list
    summarized_edges <- rbind(summarized_edges, false_negative_edges,
      stringsAsFactors = FALSE
    )
  }

  n <- nrow(summarized_edges)
  summary <- data.frame(
    x = character(n), y = character(n), type = character(n), ai = character(n),
    raw_contributions = character(n), contributions = character(n),
    info = numeric(n), info_cond = numeric(n), cplx = numeric(n),
    Nxy_ai = numeric(n), info_shifted = numeric(n), infOrt = numeric(n),
    trueOrt = numeric(n), isOrtOk = character(n), sign = character(n),
    partial_correlation = numeric(n), is_causal = character(n), proba = character(n),
    confidence = character(n), stringsAsFactors = FALSE
  )
  if(n == 0) return(summary)

  #Initialize is_causal column as NA
  summary$is_causal = NA

  # Edge ordering (A<-B or B->A) is given by lexicographical sort
  summary[,c('x','y')] = t(apply(as.data.frame(summarized_edges)[,c(1,2)], 1, function(row){sort(row)}))

  # Edge 'type' corresponds to the miic prediction : P(ositive) or N(egative)
  # for respectively presence or absence of an edge, without considering
  # orientation. If ground truth is known, edges are classified as True or
  # False Positives/Negatives (TP, FP, TN, FN).
  if (is.null(true_edges)) {
    summary$type <- apply(summary, 1, function(row, adj_matrix) {
      ifelse(adj_matrix[row[1], row[2]] == 0, "N", "P")
    }, adj_matrix)
  } else {
    type <- character(n)
    for (i in 1:nrow(summary)) {
      row <- summary[i, ]
      row_is_true_edge <- any(apply(
        true_edges, 1,
        function(true_edge, edge) {
          all(true_edge %in% edge)
        },
        c(row$x, row$y)
      ))
      if (adj_matrix[row$x, row$y] != 0) {
        # Inferred positive : either True Positive or False Positive
        type[i] <- ifelse(row_is_true_edge, "TP", "FP")
      }
      else {
        # Inferred negative : either True Negative or False Negative
        type[i] <- ifelse(row_is_true_edge, "FN", "TN")
      }
    }
    summary$type <- type
  }

  # Ai is a list containing the conditioning nodes
  summary$ai <- fill_summary_column(summary, edges, "x", "y", "ai.vect")

  # Raw contributions of each ai to the conditional independence, measured by
  # I'(x;y;ai|{aj}) / I'(x;y), where {aj} is the separating set before adding ai.
  summary$raw_contributions <- fill_summary_column(
      summary, edges, "x", "y", "raw_contributions")

  # Contributions of each ai to the reduction of conditional mutual information
  # measured by I'(x;y;ai|{aj}) / I'(x;y|{aj}), where {aj} is the separating set
  # before adding ai.
  summary$contributions <- fill_summary_column(
      summary, edges, "x", "y", "contributions")

  # info and info_cond contain the (conditional) mutual information values
  summary$info <- fill_summary_column(summary, edges, "x", "y", "Ixy")
  summary$info_cond <- fill_summary_column(summary, edges, "x", "y", "Ixy_ai")

  # cplx is the NML complexity (used for independence testing)
  summary$cplx <- fill_summary_column(summary, edges, "x", "y", "cplx")

  # Nxy_ai is the number of samples without missing values used for this edge
  summary$Nxy_ai <- fill_summary_column(summary, edges, "x", "y", "Nxy_ai")

  # info_shifted is the difference between MI and cplx
  summary$info_shifted <- summary$info_cond - summary$cplx

  # confidence is the ratio of info_shifted between randomized and normal sample
  summary$confidence <- fill_summary_column(
      summary, edges, "x", "y", "confidence")
  summary$confidence[summary$confidence == -1] <- NA

  # infOrt is the inferred edge orientation
  summary$infOrt <- apply(summary, 1, function(row, adj_matrix) {
    adj_matrix[row[1], row[2]]
  }, adj_matrix)

  # trueOrt is the true edge orientation (if known)
  if (is.null(true_edges)) {
    summary$trueOrt <- rep(NA, n)
  } else {
    true_adj_matrix <- matrix(0,
      ncol = dim(adj_matrix)[1],
      nrow = dim(adj_matrix)[1]
    )
    colnames(true_adj_matrix) <- var_names
    rownames(true_adj_matrix) <- var_names
    for (i in 1:nrow(true_edges)) {
      true_edge <- as.character(unlist(true_edges[i, ]))
      true_adj_matrix[true_edge[1], true_edge[2]] <- 2
      true_adj_matrix[true_edge[2], true_edge[1]] <- -2
    }

    summary$trueOrt <- apply(summary, 1, function(row, true_adj_matrix) {
      true_adj_matrix[row[1], row[2]]
    }, true_adj_matrix)
  }

  # isOrtOk
  summary$isOrtOk <- ifelse(summary$infOrt == summary$trueOrt, "Y", "N")

  # Sign and coefficient of partial correlation between x and y conditioned
  # on "ai"s.
  summary[, c("sign", "partial_correlation")] <- compute_partial_correlation(
    summary, observations, state_order
  )
  summary$partial_correlation <- as.numeric(summary$partial_correlation)

  # proba contains the orientation likelihoods as computed by miic (cf
  # Affeldt & Isambert, UAI 2015 proceedings) : the probabilities of
  # both orientations separated by a semi-colon.
  orientation_probabilities <- results$orientations.prob
  summary$proba <- sapply(1:nrow(summary), function(i) {
    row <- summary[i, ]
    id_x <- match(row$x, var_names)
    id_y <- match(row$y, var_names)
    proba_adj <- results$proba_adj_matrix
    if (!is.null(results$adj_matrices) && ncol(results$adj_matrices) > 1) {
      proba_adj <- results$proba_adj_average
    }
    return(paste(proba_adj[id_y, id_x], proba_adj[id_x, id_y], sep = ";"))
  })

  # Genuine causality is deducible only when latent variables are allowed and
  # propagation is not allowed
  causality_deducible <- latent && !propagation
  # If consistent parameter is turned on and the result graph is a union of
  # more than one inconsistent graphs, get the possible orientations of each
  # edge with the correponding frequencies and the consensus status.
  if (!is.null(results$adj_matrices) && ncol(results$adj_matrices) > 1) {
    # use split to turn summary data frame to list only to be able to use sapply
    # with simplify = FALSE, otherwise apply() will force simplification, which
    # is annoying and surprising when the return value of the function is the
    # same for all rows. Theoretically it can never happen here, but this is R!
    # Always be careful.
    edge_stats_table <- sapply(split(summary, seq(nrow(summary))),
      get_edge_stats_table,
      var_names,
      results$adj_matrices,
      simplify = FALSE
    )
    target <- which(names(summary) == "infOrt")[1]
    summary <- cbind(
      summary[, 1:target, drop = FALSE],
      consensus = sapply(edge_stats_table,
        get_consensus_status,
        consensus_threshold
      ),
      edge_stats = sapply(edge_stats_table, get_edge_stats_str),
      is_causal_consensus = NA,
      summary[, (target + 1):length(summary), drop = FALSE]
    )
    # Set consensus edge status according to the average probabilities
    for (i in 1:nrow(summary)) {
      row <- summary[i, ]
      if (causality_deducible) {
        # Set initial values if deducible
        if (row$infOrt != 0)
          summary[i, ]$is_causal <- "N"
        if (row$consensus != 0)
          summary[i, ]$is_causal_consensus <- "N"
      }
      if (row$consensus == 0) next

      id_x <- match(row$x, var_names)
      id_y <- match(row$y, var_names)
      # probability of an edge tip being a head (<,>), * means head or tail (-)
      proba_x2y <- results$proba_adj_average[id_x, id_y]  # proba of x *-> y
      proba_y2x <- results$proba_adj_average[id_y, id_x]  # proba of x <-* y
      ratio_x2y <- (1 - proba_x2y) / proba_x2y
      ratio_y2x <- (1 - proba_y2x) / proba_y2x

      if (ratio_x2y < ori_consensus_ratio && ratio_y2x < ori_consensus_ratio) {
        summary[i, ]$consensus <- 6
      } else if (ratio_x2y < ori_consensus_ratio &&
                 ratio_y2x >= ori_consensus_ratio) {
        summary[i, ]$consensus <- 2
        if (1 / ratio_y2x < ori_consensus_ratio && causality_deducible) {
          summary[i, ]$is_causal_consensus <- "Y"
          if (row$infOrt == 2) {
            summary[i, ]$is_causal <- "Y"
          }
        }
      } else if (ratio_y2x < ori_consensus_ratio &&
                 ratio_x2y >= ori_consensus_ratio) {
        summary[i, ]$consensus <- -2
        if (1 / ratio_x2y < ori_consensus_ratio && causality_deducible) {
          summary[i, ]$is_causal_consensus <- "Y"
          if (row$infOrt == -2) {
            summary[i, ]$is_causal <- "Y"
          }
        }
      } else {
        summary[i, ]$consensus <- 1
      }
    }
  } else if (causality_deducible) {
    # if (is.null(results$adj_matrices) || ncol(results$adj_matrices) <= 1)
    # set is_causal by results$proba_adj_matrix
    for (i in 1:nrow(summary)) {
      row <- summary[i, ]
      if (row$infOrt == 0) {
        next
      }
      summary[i, ]$is_causal <- "N"

      id_x <- match(row$x, var_names)
      id_y <- match(row$y, var_names)
      # probability of an edge tip being a head (<,>), * means head or tail (-)
      proba_x2y <- results$proba_adj_matrix[id_x, id_y]  # proba of x *-> y
      proba_y2x <- results$proba_adj_matrix[id_y, id_x]  # proba of x <-* y
      ratio_x2y <- (1 - proba_x2y) / proba_x2y
      ratio_y2x <- (1 - proba_y2x) / proba_y2x
      if (row$infOrt == 2 &&
          ratio_x2y < ori_consensus_ratio &&
          1 / ratio_y2x < ori_consensus_ratio) {
        summary[i, ]$is_causal <- "Y"
      }
      if (row$infOrt == -2 &&
          ratio_y2x < ori_consensus_ratio &&
          1 / ratio_x2y < ori_consensus_ratio) {
        summary[i, ]$is_causal <- "Y"
      }
    }
  }

  # Sort summary by log confidence and return it
  summary <- summary[order(summary$info_shifted, decreasing = TRUE), ]
  rownames(summary) <- c()
  return(summary)
}

matrix_from_3_columns <- function(df, rows, columns, values) {
  x = df[[rows]]
  y = df[[columns]]
  nodes = unique(c(x,y))
  with(df, {
  out <- matrix(0, nrow=length(nodes), ncol=length(nodes),
                dimnames=list(nodes, nodes))
  out[cbind(x, y)] <- df[[values]]
  out[cbind(y, x)] <- df[[values]]
  out
  })
}

fill_summary_column <- function(summary, matrix, rows, columns, values) {
  # This function uses the information from matrix to return a vector of
  # values (from the column named `values`).
  # It creates a 2D matrix indexed with the `rows` and `columns` values for
  # faster lookup, and returns the values for all the combinations of rows and
  # columns observed in the `summary` matrix.

  wide_matrix <- matrix_from_3_columns(matrix, rows, columns, values)

  apply(summary, 1, function(row, wide_matrix) {
    wide_matrix[row[1], row[2]]
  }, wide_matrix)
}
#===============================================================================
# END TO_BE_DELETED AFTER TESTS
#===============================================================================

compute_partial_correlation <- function(summary, observations, state_order) {
  ppcor_results <- data.frame(
    sign = character(nrow(summary)),
    partial_correlation = numeric(nrow(summary)),
    stringsAsFactors = FALSE
  )

  for (j in 1:ncol(observations)) {
    col <- colnames(observations)[j]
    # row index in state_order, NA if not found/available
    row <- match(col, state_order$var_names)
    order_string <- NULL
    if (!is.null(state_order$levels_increasing_order) && !is.na(row)) {
      order_string <- state_order[row, "levels_increasing_order"]
    }
    # If the variable is described in the state order file, transform to a
    # factor with the right category order.
    if (!is.null(order_string) && !is.na(order_string) &&
        (is.null(state_order$var_type) || state_order[row, "var_type"] == 0)) {
      # Convert ordered categorical features to integers
      observations[, col] <- factor(observations[, col])
      ordered_levels <- unlist(strsplit(as.character(order_string), ","))
      ordered_levels <- tryCatch(
        as.numeric(ordered_levels),
        warning = function(w) {return(ordered_levels)})
      # levels(observations[,col]) = ordered_levels
      observations[, col] <- ordered(observations[, col], ordered_levels)
      observations[, col] <- as.numeric(observations[, col])
    } else if (is.factor(observations[,col]) &&
      suppressWarnings(all(!is.na(as.numeric(levels(observations[,col])))))) {
      # If the variable is not described but categorical and numerical, assume
      # its order from its values.
      observations[, col] <- as.numeric(observations[, col])
    }
  }

  for (i in which(summary$type %in% c('P', 'TP', 'FP'))) {
    x <- summary[i, "x"]
    y <- summary[i, "y"]
    ai <- summary[i, "ai"]  # String with Ais separated by comma, or NA
    ppcor_results[i, ] <- c("NA", NA)

    if (is.na(ai)) {
      ai <- NULL
    } else {
      ai <- unlist(strsplit(ai, ","))  # Character vector
    }

    if (!all(sapply(observations[, c(x, y, ai)], is.numeric))) next

    if (is.null(ai)) {
      OK <- stats::complete.cases(observations[, c(x, y)])
      if (sum(OK) < 2) next

      edge_res <- stats::cor.test(
        observations[OK, x],
        observations[OK, y],
        method = "spearman",
        exact = FALSE
      )
    } else {
      OK <- stats::complete.cases(observations[, c(x, y, ai)])
      if (sum(OK) < 2) next

      edge_res <- ppcor::pcor.test(
        observations[OK, x],
        observations[OK, y],
        observations[OK, ai],
        method = "spearman"
      )
    }

    # Save sign and coef
    ppcor_results[i, ] <- c(
      ifelse(edge_res$estimate >= 0, "+", "-"),
      edge_res$estimate
    )
    # Sign is either positive or negative... maybe put a threshold for very
    # low absolute values ?
  }

  return(ppcor_results)
}

#-------------------------------------------------------------------------------
# get_edge_stats_table
#-------------------------------------------------------------------------------
get_edge_stats_table <- function(row, var_names, adj_matrices) {
  # row[1]: variable name of x
  # row[2]: variable name of y
  n_var <- length(var_names)
  # adj_matrices is of dimension (n_var * n_var, n_cycle), i.e., each
  # column is a 1-d adjacency matrix
  index_1d <- n_var * (match(row[1], var_names) - 1) + match(row[2], var_names)
  n_cycle <- dim(adj_matrices)[2]
  # edge stats table, count replaced by frequency (percentage)
  t <- table(adj_matrices[index_1d,]) / n_cycle
  t <- t[order(t, decreasing = TRUE), drop = FALSE]
  return(t)
}

#-------------------------------------------------------------------------------
# get_edge_stats_str
#-------------------------------------------------------------------------------
get_edge_stats_str <- function(stats_table) {
  t <- sapply(stats_table, scales::percent_format())
  # return a ";" separated string of format "percentage(orientation)"
  return(paste(t, "(", names(t), ")", sep = "", collapse = ";"))
}

#-------------------------------------------------------------------------------
# get_consensus_status
#-------------------------------------------------------------------------------
# 0: unconnected, 1: connected
#-------------------------------------------------------------------------------
get_consensus_status <- function(stats_table, consensus_threshold) {
  if (length(stats_table) < 1)
    return(NA)
  freq_no_edge <- unname(stats_table["0"])
  # "0" doesn't exist or the percentage of non-"0" is above the threshold
  if (is.na(freq_no_edge) || freq_no_edge < 1 - consensus_threshold)
    return(1)
  return(0)
}
