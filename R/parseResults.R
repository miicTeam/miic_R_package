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
# - ort_consensus_ratio: a float, optional, 0.1 by default. Used to determine
#   if oriented edges are genuine causal and, when consistency is activated,
#   to determine the consensus graph orientations.
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
#   * raw_contributions: raw contributions of each ai to the conditional
#     independence, measured by I'(x;y;ai|{aj}) / I'(x;y),
#     where {aj} is the separating set before adding ai.
#   * contributions : contributions of each ai to the reduction of conditional
#     mutual information measured by I'(x;y;ai|{aj}) / I'(x;y|{aj}),
#     where {aj} is the separating set before adding ai.
#   * info: mutual information (corresponds to i_xy in the C++ output)
#   * n_xy: the number of samples without missing values for the pair of
#     variables
#   * info_cond: conditional mutual information  (corresponds to i_xy_ai in the
#     C++ output)
#   * cplx: the NML complexity (used for independence testing)
#   * n_xy_ai: the number of samples without missing values for the pair of
#     variables and the contributors
#   * info_shifted: the difference between conditional MI and cplx
#   * ort_inferred: the inferred edge orientation
#   * ort_ground_truth: is the true edge orientation if known.
#     NA if truth is unknown
#   * is_inference_correct: indicate if the inferred edge is correctly oriented.
#     NA if truth is unknown, TRUE or FALSE if truth is known
#   * is_causal: indicates if edge is genuine causal.
#     Note that the genuine causality is deducible only when latent variables
#     are allowed and propagation is not allowed
#   * ort_consensus: if consistency is activated, indicates the consensus
#     orientation of the edge, possible values are 0: not connected,
#     1: not oriented, -2 or 2: oriented or 6: bi-directional (latent variable)
#     NA if consistency is not activated.
#   * is_causal_consensus: if consistency is activated, indicates if the
#     consensus orientation is genuine causal. NA if consistency is not
#     activated.
#   * edge_stats: if consistency is activated, contains the orientation
#     frequencies of each orientation present in the cycle of graphs.
#     NA if consistency is not activated.
#   * sign: sign of partial correlation between x and y conditioned on "ai"s
#   * partial_correlation: coefficient of partial correlation between x and y
#     conditioned on "ai"s
#   * p_y2x: probability of the arrowhead from y to x, NA for removed edges.
#   * p_x2y: probability of the arrowhead from x to y, NA for removed edges.
#   * confidence: the ratio of info_shifted between randomized
#     and normal samples
#-------------------------------------------------------------------------------
summarizeResults = function (observations, results,
                             true_edges = NULL, state_order = NULL,
                             consensus_threshold = 0.8, ort_consensus_ratio = 0.1,
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
    return (data.frame (x = character(0), y = character(0),
      type = character(0), ai = character(0),
      raw_contributions = character(0), contributions = character(0),
      info = numeric(0), n_xy = numeric(0), info_cond = numeric(0),
      cplx = numeric(0), n_xy_ai = numeric(0), info_shifted = numeric(0),
      ort_inferred = integer(0), ort_ground_truth = integer(0),
      is_inference_correct = logical(0), is_causal = logical(0),
      ort_consensus = integer(0), is_causal_consensus = logical(0),
      edge_stats = character(0), sign = character(0),
      partial_correlation = numeric(0), p_y2x = numeric(0), p_x2y = numeric(0),
      confidence = numeric(0), stringsAsFactors = FALSE) )
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
  # info and info_cond contain the (conditional) mutual information values
  #
  colnames(summary)[ which (colnames(summary) == "i_xy")] <- "info"
  colnames(summary)[ which (colnames(summary) == "i_xy_ai")] <- "info_cond"
  #
  # info_shifted is the difference between MI and cplx
  #
  summary$info_shifted <- summary$info_cond - summary$cplx
  #
  # confidence is the ratio of info_shifted between randomized
  # and normal samples
  #
  summary$confidence [summary$confidence == -1] <- NA_real_
  #
  # ort_inferred is the inferred edge orientation
  #
  summary$ort_inferred <- apply (summary, 1, function (row, adj_mat) {
      adj_mat[row[1], row[2]] }, results$adj_matrix)
  #
  # ort_ground_truth is the true edge orientation (if known)
  #
  var_names <- colnames (results$adj_matrix)
  n_vars <- length (var_names)
  if ( is.null(true_edges) )
    summary$ort_ground_truth <- NA_integer_
  else
    {
    true_adj_matrix <- matrix (0, ncol=n_vars, nrow=n_vars,
                               dimnames=list(var_names,var_names) )
    for ( i in 1:nrow (true_edges) )
      {
      true_edge <- as.character (unlist (true_edges[i, ]) )
      true_adj_matrix[true_edge[1], true_edge[2]] <- 2
      true_adj_matrix[true_edge[2], true_edge[1]] <- -2
      }
    summary$ort_ground_truth <- apply (summary, 1, function (row, true_adj_matrix) {
      true_adj_matrix[row[1], row[2]] }, true_adj_matrix)
    }
  #
  # is_inference_correct indicates if the inferred edge is correctly oriented.
  # NA if truth is unknown, TRUE or FALSE if truth is known
  #
  summary$is_inference_correct <- ifelse (summary$ort_inferred == summary$ort_ground_truth,
                             TRUE, FALSE)
  #
  # Sign and coefficient of partial correlation between x and y conditioned
  # on "ai"s.
  #
  summary [, c("sign", "partial_correlation")] <- compute_partial_correlation (
    summary, observations, state_order)
  summary$partial_correlation <- as.numeric (summary$partial_correlation)
  #
  # Probabilities of orientations
  #
  if (!is.null(results$adj_matrices) && length(results$adj_matrices) > 1)
    tmp_proba_adj <- results$proba_adj_average
  else
    tmp_proba_adj <- results$proba_adj_matrix
  summary$p_y2x <- unlist (lapply (1:nrow(summary), function(i) {
    proba_of_edge <- tmp_proba_adj[ summary[i, "y"], summary[i, "x"] ]
    ifelse (proba_of_edge == -1, NA_real_, proba_of_edge)
    } ) )
  summary$p_x2y <- unlist (lapply (1:nrow(summary), function(i) {
    proba_of_edge <- tmp_proba_adj[ summary[i, "x"], summary[i, "y"] ]
    ifelse (proba_of_edge == -1, NA_real_, proba_of_edge)
    } ) )
  #
  # Genuine causality is deducible only when latent variables are allowed and
  # propagation is not allowed
  #
  causality_deducible <- latent && (!propagation)
  summary$is_causal = as.logical (NA)
  summary$ort_consensus = NA_integer_
  summary$is_causal_consensus = as.logical (NA)
  summary$edge_stats = NA_character_

  if (  ( ! is.null (results$adj_matrices) )
     && (length (results$adj_matrices) > 1) )
    {
    # If consistent parameter is turned on and the result graph is a union of
    # more than one inconsistent graphs, get the possible orientations of each
    # edge with the corresponding frequencies and the consensus status.
    #
    n_cycles = length (results$adj_matrices)
    edge_stats_table <- lapply (1:nrow(summary), function(i) {
      list_adj <- unlist (lapply (results$adj_matrices,
        FUN=function (z) { z[ summary[i, "x"], summary[i, "y"] ] }) )
      t <- table(list_adj) / n_cycles
      t <- t[order(t, decreasing = TRUE), drop = FALSE]
      return (t)
      })

    summary$ort_consensus <- unlist (lapply (edge_stats_table,
                                             get_consensus_status,
                                             consensus_threshold) )
    summary$edge_stats <- unlist (lapply (edge_stats_table,
                                          get_edge_stats_str) )
    #
    # Set consensus edge status according to the average probabilities
    #
    for (i in 1:nrow(summary))
      {
      row <- summary[i, ]
      if (causality_deducible)
        {
        # Set initial values if deducible
        if (row$ort_inferred != 0)
          summary[i, "is_causal"] <- FALSE
        if (row$ort_consensus != 0)
          summary[i, "is_causal_consensus"] <- FALSE
        }
      if (row$ort_consensus == 0)
        next

      # probability of an edge tip being a head (<,>), * means head or tail (-)
      proba_x2y <- results$proba_adj_average[row$x, row$y]  # proba of x *-> y
      proba_y2x <- results$proba_adj_average[row$y, row$x]  # proba of x <-* y
      ratio_x2y <- (1 - proba_x2y) / proba_x2y
      ratio_y2x <- (1 - proba_y2x) / proba_y2x

      if (  (ratio_x2y < ort_consensus_ratio)
         && (ratio_y2x < ort_consensus_ratio) )
        summary[i, "ort_consensus"] <- 6
      else if (  (ratio_x2y < ort_consensus_ratio)
              && (ratio_y2x >= ort_consensus_ratio) )
        {
        summary[i, "ort_consensus"] <- 2
        if (1 / ratio_y2x < ort_consensus_ratio && causality_deducible)
          {
          summary[i, "is_causal_consensus"] <- TRUE
          if (row$ort_inferred == 2)
            summary[i, "is_causal"] <- TRUE
          }
        }
      else if (  (ratio_y2x < ort_consensus_ratio)
              && (ratio_x2y >= ort_consensus_ratio) )
        {
        summary[i, "ort_consensus"] <- -2
        if (1 / ratio_x2y < ort_consensus_ratio && causality_deducible)
          {
          summary[i, "is_causal_consensus"] <- TRUE
          if (row$ort_inferred == -2)
            summary[i, "is_causal"] <- TRUE
          }
        }
      else
        summary[i, "ort_consensus"] <- 1
      }
    summary$ort_consensus = as.integer (summary$ort_consensus)
    }
  else if (causality_deducible)
    {
    # Consistency not activated or only one graph in the cycle
    #
    for (i in 1:nrow(summary))
      {
      row <- summary[i, ]
      if (row$ort_inferred == 0)
        next
      summary[i, "is_causal"] <- FALSE

      # probability of an edge tip being a head (<,>), * means head or tail (-)
      proba_x2y <- results$proba_adj_matrix[row$x, row$y]  # proba of x *-> y
      proba_y2x <- results$proba_adj_matrix[row$y, row$x]  # proba of x <-* y
      ratio_x2y <- (1 - proba_x2y) / proba_x2y
      ratio_y2x <- (1 - proba_y2x) / proba_y2x
      if (  (row$ort_inferred == 2)
         && (ratio_x2y < ort_consensus_ratio)
         && (1 / ratio_y2x < ort_consensus_ratio) )
        summary[i, "is_causal"] <- TRUE
      if (  (row$ort_inferred == -2)
         && (ratio_y2x < ort_consensus_ratio)
         && (1 / ratio_x2y < ort_consensus_ratio) )
        summary[i, "is_causal"] <- TRUE
      }
    }
  #
  # Sort summary by log likelihood, keep only some cols and returns
  #
  columns_kept = c ("x", "y", "type", "ai", "raw_contributions", "contributions",
    "info", "n_xy", "info_cond", "cplx", "n_xy_ai", "info_shifted",
    "ort_inferred", "ort_ground_truth", "is_inference_correct", "is_causal",
    "ort_consensus", "is_causal_consensus", "edge_stats",
    "sign", "partial_correlation", "p_y2x", "p_x2y", "confidence")
  columns_kept = columns_kept[columns_kept %in% colnames(summary)]
  summary <- summary[order(summary$info_shifted, decreasing = TRUE),
                     columns_kept, drop=F]
  rownames(summary) <- c()
  return (summary)
  }

#-------------------------------------------------------------------------------
# compute_partial_correlation
#-------------------------------------------------------------------------------
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
