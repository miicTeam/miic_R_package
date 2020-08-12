summarizeResults <- function(observations = NULL, results = NULL,
                             true_edges = NULL, state_order = NULL,
                             consensus_threshold = 0.8, verbose = FALSE) {
  # Reduced list of edges that will be summarized. There are 3 categories:
  # - Edges that exist in the miic reconstruction (oriented or not)
  # - Edges that were conditioned away with a non-empty separating set
  # - If ground truth is known, any other positive edge
  summarized_edges <- matrix(ncol = 2)
  adj_matrix <- results$adj_matrix
  var_names <- colnames(adj_matrix)

  # List of edges found by miic
  half_adj_matrix = adj_matrix
  half_adj_matrix[ lower.tri(adj_matrix, diag = TRUE) ] <- 0
  predicted_edges <- which(half_adj_matrix != 0, arr.ind = T, useNames = F)
  predicted_edges <- apply(predicted_edges, 2, function(x) { var_names[x] })
  # Add to summarized edges list
  summarized_edges <- predicted_edges

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
      stringsAsFactors = F
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
      stringsAsFactors = F
    )
  }


  n <- nrow(summarized_edges)
  summary <- data.frame(
    x = character(n), y = character(n), type = character(n), ai = character(n),
    info = numeric(n), info_cond = numeric(n), cplx = numeric(n),
    Nxy_ai = numeric(n), log_confidence = numeric(n), infOrt = numeric(n),
    trueOrt = numeric(n), isOrtOk = character(n), sign = character(n),
    partial_correlation = numeric(n), isCausal = character(n),
    proba = character(n),
    stringsAsFactors = FALSE
  )

  # Edge ordering (A<-B or B->A) is given by lexicographical sort
  summary[,c('x','y')] = t(apply(summarized_edges[,c(1,2)], 1, function(row){sort(row)}))

  # Edge 'type' correponds to the miic prediction : P(ositive) or N(egative)
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

  # info and info_cond contain the (conditional) mutual information values
  summary$info <- fill_summary_column(summary, edges, "x", "y", "Ixy")
  summary$info_cond <- fill_summary_column(summary, edges, "x", "y", "Ixy_ai")

  # cplx is the NML complexity (used for independence testing)
  summary$cplx <- fill_summary_column(summary, edges, "x", "y", "cplx")

  # Nxy_ai is the number of samples without missing values used for this edge
  summary$Nxy_ai <- fill_summary_column(summary, edges, "x", "y", "Nxy_ai")

  # log_confidence is the difference between MI and cplx
  summary$log_confidence <- summary$info_cond - summary$cplx

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

  orientation_probabilities <- results$orientations.prob
  # isCausal considers whether the source node of an oriented edge is also
  # at the tip of a V-structure. If so, then the source node has a stronger
  # indication of being causal.
  summary$isCausal <- "NA"
  if ((!is.null(nrow(orientation_probabilities))) &&
      nrow(orientation_probabilities) > 0 ) {
    summary$isCausal <- is_causal(summary, orientation_probabilities)
  }

  # proba contains the orientation likelihoods as computed by miic (cf
  # Affeldt & Isambert, UAI 2015 proceedings) : the probabilities of
  # both orientations separated by a semi-colon.
  summary$proba <- sapply(1:nrow(summary), function(i) {
    paste(findProba(summary, orientation_probabilities, i), collapse = ";")
  })

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
      summary[, (target + 1):length(summary), drop = FALSE]
    )
  }

  # Sort summary by log confidence and return it
  summary <- summary[order(summary$log_confidence, decreasing = T), ]
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

compute_partial_correlation <- function(summary, observations, state_order) {
  ppcor_results <- data.frame(
    sign = character(nrow(summary)),
    partial_correlation = numeric(nrow(summary)),
    stringsAsFactors = FALSE
  )

  for (j in 1:ncol(observations)) {
    col <- colnames(observations)[j]
    # If the variable is described in the state order file, transform to a
    # factor with the right category order.
    if (!is.null(state_order) &&
      state_order[j, "var_type"] == 0 &&
      "levels_increasing_order" %in% colnames(state_order) &&
      !is.na(state_order[state_order$var_names == col, "levels_increasing_order"])) {
      rownames(state_order) <- state_order$var_names
      # Convert ordered categorical features to integers
      observations[, col] <- factor(observations[, col])
      ordered_levels <- as.character(
        state_order[state_order$var_names == col, "levels_increasing_order"])
      ordered_levels <- unlist(strsplit(ordered_levels, ","))
      # levels(observations[,col]) = ordered_levels
      observations[, col] <- ordered(observations[, col], ordered_levels)
      observations[, col] <- as.numeric(observations[, col])
    } else if (is.factor(observations[,col]) &&
      suppressWarnings(all(!is.na(as.numeric(levels(observations[,col])))))) {
      # If the variable is not described but numerical, assume its order from
      # the numerical categories.
      observations[, col] <- as.numeric(observations[, col])
    }
  }

  for (i in which(summary$type=="P")) {
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
        method = "spearman"
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

is_causal <- function(summary, probas) {
  is_causal_results <- rep("N", (nrow(summary)))
  v_structs <- probas$NI3 < 0
  if (length(which(v_structs)) == 0) {
    return(is_causal_results)
  }
  probas_V_structs = probas[v_structs,]

  vstruct_matrix <- matrix_from_3_columns(
    probas_V_structs, "source1",
    "source2", "target"
  )
  error_matrix <- matrix_from_3_columns(
    probas_V_structs, "source1",
    "source2", "Error"
  )

  for (i in 1:nrow(summary)) {
    row <- summary[i, ]
    if (row$type %in% c("TN", "FN", "N") || # Negative edge
      row$infOrt == 1) { # Non oriented edge
      # If any of these conditions is true, edge cannot be causal
      next
    }
    if (row$infOrt == 6) {
      is_causal_results[i] <- "Y"
      next
    }

    if (row$infOrt == 2) {
      source <- row$x
      target <- row$y
    } else {
      source <- row$y
      target <- row$x
    }

    if (!any(probas_V_structs[, "target"] == source)) {
      # Propagated orientation
      next
    }

    # If there is another V-structure pointing to the source node, which
    # is not an error.
    if (any(error_matrix[vstruct_matrix == source] == 0)) {
      is_causal_results[i] <- "Y"
    }
  }

  return(is_causal_results)
}

findProba <- function(outputSummary.df, proba, index) {
  # I'm sorry... I didn't rewrite that one

  # find probability of this edge in the proba file
  if (outputSummary.df[index, "infOrt"] == 2) {
    posProba <- which((proba[, 4] == outputSummary.df[index, 2] &
      proba[, 7] == outputSummary.df[index, 1]))
    # source2 -> target
    if (length(posProba) > 0) {
      posProba <- posProba[1]
      return(c(proba[posProba, "p3"], proba[posProba, "p4"]))
    }

    posProba <- which((proba[, 1] == outputSummary.df[index, 1] &
      proba[, 4] == outputSummary.df[index, 2]))
    if (length(posProba) > 0) {
      posProba <- posProba[1]
      return(c(proba[posProba, "p2"], proba[posProba, "p1"]))
    }

    posProba <- which((proba[, 4] == outputSummary.df[index, 1] &
      proba[, 1] == outputSummary.df[index, 2]))
    if (length(posProba) > 0) {
      posProba <- posProba[1]
      return(c(proba[posProba, "p1"], proba[posProba, "p2"]))
    }

    posProba <- which((proba[, 4] == outputSummary.df[index, 1] &
      proba[, 7] == outputSummary.df[index, 2]))
    if (length(posProba) > 0) {
      posProba <- posProba[1]
      return(c(proba[posProba, "p4"], proba[posProba, "p3"]))
    }
  } else if (outputSummary.df[index, "infOrt"] == -2) {
    posProba <- which((proba[, 4] == outputSummary.df[index, 1] &
      proba[, 7] == outputSummary.df[index, 2]))
    if (length(posProba) > 0) {
      posProba <- posProba[1]
      return(c(proba[posProba, "p3"], proba[posProba, "p4"]))
    }

    posProba <- which((proba[, 1] == outputSummary.df[index, 2] &
      proba[, 4] == outputSummary.df[index, 1]))
    if (length(posProba) > 0) {
      posProba <- posProba[1]
      return(c(proba[posProba, "p2"], proba[posProba, "p1"]))
    }

    posProba <- which((proba[, 4] == outputSummary.df[index, 2] &
      proba[, 1] == outputSummary.df[index, 1]))
    if (length(posProba) > 0) {
      posProba <- posProba[1]
      return(c(proba[posProba, "p1"], proba[posProba, "p2"]))
    }

    posProba <- which((proba[, 4] == outputSummary.df[index, 2] &
      proba[, 7] == outputSummary.df[index, 1]))
    if (length(posProba) > 0) {
      posProba <- posProba[1]
      return(c(proba[posProba, "p4"], proba[posProba, "p3"]))
    }
  } else if (outputSummary.df[index, "infOrt"] == 6) {
    posProba <- which((proba[, 4] == outputSummary.df[index, 1] &
      proba[, 7] == outputSummary.df[index, 2]))
    if (length(posProba) > 0) {
      posProba <- posProba[1]
      return(c(proba[posProba, "p4"], proba[posProba, "p3"]))
    }

    posProba <- which((proba[, 7] == outputSummary.df[index, 1] &
      proba[, 4] == outputSummary.df[index, 2]))
    if (length(posProba) > 1) {
      posProba <- posProba[1]
      return(c(proba[posProba, "p3"], proba[posProba, "p4"]))
    }

    posProba <- which((proba[, 1] == outputSummary.df[index, 1] &
      proba[, 4] == outputSummary.df[index, 2]))
    if (length(posProba) > 0) {
      posProba <- posProba[1]
      return(c(proba[posProba, "p2"], proba[posProba, "p1"]))
    }

    posProba <- which((proba[, 4] == outputSummary.df[index, 1] &
      proba[, 1] == outputSummary.df[index, 2]))
    if (length(posProba) > 0) {
      posProba <- posProba[1]
      return(c(proba[posProba, "p1"], proba[posProba, "p2"]))
    }
  }
  return(NA)
}

get_edge_stats_table <- function(row, var_names, adj_matrices) {
  # row[1]: variable name of x
  # row[2]: variable name of y
  n_var <- length(var_names)
  # adj_matrices is of dimention (n_var * n_var, n_cycle), i.e., each
  # column is a 1-d adjacency matrix
  index_1d <- n_var * (match(row[1], var_names) - 1) + match(row[2], var_names)
  n_cycle <- dim(adj_matrices)[2]
  # edge stats table, count replaced by frequency (percentage)
  t <- table(adj_matrices[index_1d,]) / n_cycle
  t <- t[order(t, decreasing = TRUE), drop = FALSE]
  return(t)
}

get_edge_stats_str <- function(stats_table) {
  t <- sapply(stats_table, scales::percent_format())
  # return a ";" separated string of format "percentage(orientation)"
  return(paste(t, "(", names(t), ")", sep = "", collapse = ";"))
}

get_consensus_status <- function(stats_table, consensus_threshold) {
  if (length(stats_table) < 1)
    return(NA)
  # stats_table is sorted, stats_table[[1]] is the highest frequency
  if (stats_table[[1]] < consensus_threshold)
    return(1)  # undirected
  return(as.numeric(names(stats_table)[[1]]))
}
