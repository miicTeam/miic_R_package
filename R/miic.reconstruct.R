miic.reconstruct <- function(input_data = NULL,
                             is_contextual = NULL,
                             is_consequence = NULL,
                             is_continuous = NULL,
                             black_box = NULL,
                             n_threads = 1,
                             n_eff = -1,
                             cplx = "nml",
                             eta = 1,
                             latent = "orientation",
                             n_shuffles = 0,
                             orientation = TRUE,
                             ort_proba_ratio = 1,
                             propagation = FALSE,
                             conf_threshold = 0,
                             verbose = FALSE,
                             sample_weights = NULL,
                             test_mar = TRUE,
                             consistent = "no",
                             max_iteration = NULL,
                             mode = "S",
                             n_layers = NULL,
                             delta_t = NULL,
                             negative_info = FALSE
                             ) {
  n_samples <- nrow(input_data)
  n_nodes <- ncol(input_data)
  # Numeric factor matrix, level starts from 0, NA mapped to -1
  input_factor <- apply(input_data, 2, function(x)
                        (as.numeric(factor(x, levels = unique(x))) - 1))
  input_factor[is.na(input_factor)] <- -1
  max_level_list <- as.numeric(apply(input_factor, 2, max)) + 1
  input_factor <- as.vector(as.matrix(input_factor))
  # Data list, numeric for continuous columns, -1 for discrete columns
  input_double <- matrix(nrow = n_samples, ncol = n_nodes)
  # Order list, order(column) for continuous columns (index starting from 0, NA
  # mapped to -1), -1 for discrete columns
  input_order <- matrix(nrow = n_samples, ncol = n_nodes)
  for (i in c(1:ncol(input_data))) {
    if (is_continuous[i]) {
      input_double[, i] <- as.numeric(input_data[, i])
      n_NAs <- sum(is.na(input_data[, i]))
      input_order[, i] <- c(order(input_data[, i], na.last=NA) - 1,
                            rep_len(-1, n_NAs))
    } else {
      input_double[, i] <- rep_len(-1, n_samples)
      input_order[, i] <- rep_len(-1, n_samples)
    }
  }
  input_order <- as.vector(input_order)
  input_double <- as.vector(input_double)

  var_names <- colnames(input_data)
  n_vars <- length (var_names)

  arg_list <- list(
    "conf_threshold" = conf_threshold,
    "consistent" = consistent,
    "cplx" = cplx,
    "degenerate" = FALSE,
    "eta" = eta,
    "half_v_structure" = 0,
    "is_continuous" = as.numeric(is_continuous),
    "is_k23" = TRUE,
    "latent" = latent,
    "levels" = max_level_list,
    "max_iteration" = max_iteration,
    "n_eff" = n_eff,
    "n_nodes" = n_nodes,
    "n_samples" = n_samples,
    "n_shuffles" = n_shuffles,
    "n_threads" = n_threads,
    "no_init_eta" = FALSE,
    "orientation" = orientation,
    "ort_proba_ratio" = ort_proba_ratio,
    "propagation" = propagation,
    "test_mar" = test_mar,
    "mode" = mode,
    "negative_info" = negative_info,
    "max_bins" = min(50, n_samples),
    "var_names" = var_names,
    "verbose" = verbose
  )
  if (!is.null(black_box)) {
    # transform var names to var indices
    black_box[] <- sapply(black_box, function(x) {
      match(as.character(x), colnames(input_data)) - 1 } )
    black_box[] <- black_box[stats::complete.cases(black_box),]
    arg_list[["black_box"]] <- as.vector(as.matrix(t(black_box)))
  }
  if (!is.null(sample_weights))
    arg_list[["sample_weights"]] <- sample_weights
  if (!is.null(is_contextual))
    arg_list[["is_contextual"]] <- is_contextual
  if (!is.null(is_consequence))
    arg_list[["is_consequence"]] <- is_consequence
  if (!is.null(n_layers))
    arg_list[["n_layers"]] <- n_layers
  if (!is.null(delta_t))
    arg_list[["delta_t"]] <- delta_t

  cpp_input <- list("factor" = input_factor, "double" = input_double,
                    "order" = input_order)
  # Call C++ function
  res <- reconstruct(cpp_input, arg_list)
  if (res$interrupted)
    return(list(interrupted = TRUE))

  # R-formalize returned object
  #
  # Table of edges information
  #
  n_row <- length(res$edges) - 1
  header <- unlist(res$edges[1])
  df <- data.frame(matrix(unlist(res$edges[2:(n_row + 1)]), nrow = n_row,
                          byrow = TRUE), stringsAsFactors = FALSE)
  colnames(df) <- header
  df[df == "NA"] <- NA
  df$i_xy <- as.numeric(df$i_xy)
  df$i_xy_ai <- as.numeric(df$i_xy_ai)
  df$cplx <- as.numeric(df$cplx)
  df$r_xyz_ai <- as.numeric(df$r_xyz_ai)
  df$n_xy_ai <- as.numeric(df$n_xy_ai)
  df$confidence <- as.numeric(df$confidence)
  res$edges <- df
  #
  # Reshape items returned as list into their correct shape,
  # add column/row names
  #
  res$adj_matrix <- matrix (unlist (res$adj_matrix),
                            ncol=n_vars, nrow=n_vars, byrow=TRUE,
                            dimnames=list (var_names, var_names) )

  res$proba_adj_matrix <- matrix (unlist(res$proba_adj_matrix),
                                  ncol=n_vars, nrow=n_vars, byrow=TRUE,
                                  dimnames=list (var_names, var_names) )
  #
  # Same : reshape items returned when consistent parameter is turned on
  #
  if (length (res$adj_matrices) > 0)
    {
    tmp_reshape = list()
    for (i in 1:length (res$adj_matrices) )
      tmp_reshape[[i]] = matrix (unlist (res$adj_matrices[[i]]),
                                 ncol=n_vars, nrow=n_vars, byrow=TRUE,
                                 dimnames=list (var_names, var_names) )
    res$adj_matrices = tmp_reshape
    }

  if (length (res$proba_adj_matrices) > 0)
    {
    # First reshape with n_vars * n_vars rows, n_cycles columns to compute mean
    #
    tmp_reshape <- matrix (unlist (res$proba_adj_matrices),
                           ncol=length(res$proba_adj_matrices))
    tmp_reshape[tmp_reshape == -1] <- NA
    adj_average <- rowMeans (tmp_reshape, na.rm = TRUE)
    res$proba_adj_average <- matrix (unlist (adj_average),
                                     ncol=n_vars, nrow=n_vars, byrow=TRUE,
                                     dimnames=list (var_names, var_names) )
    #
    # Final reshape into a list of n_cycles items, each item is n_vars rows,
    # n_vars columns matrix
    #
    res$proba_adj_matrices = list()
    for (i in 1:ncol (tmp_reshape) )
      res$proba_adj_matrices[[i]] = matrix (unlist (tmp_reshape[,i]),
                                            ncol=n_vars, nrow=n_vars, byrow=TRUE,
                                            dimnames=list (var_names, var_names) )
    }

  # save time
  time <- strsplit(as.character(res$time), " ")
  time[which(time == 0)] <- NA

  res$time <- stats::setNames(
    as.numeric(time),
    c("init", "iter", "cut", "ort", "cpp")
  )

  # create the data frame of the structures after orientation
  orientations_prob <- res$triples

  if (length(res$triples) > 0) {
    a <- length(orientations_prob[[1]])
    b <- length(unlist(orientations_prob))
    tmp <- unlist(res$triples)[1:a]
    res1 <- unlist(res$triples)[(a + 1):b]
    orientations_prob <- data.frame(matrix(
      res1,
      nrow = length(orientations_prob) - 1,
      byrow = TRUE
    ),
    stringsAsFactors = FALSE
    )
    colnames(orientations_prob) <- tmp

    orientations_prob[, c(2:3)] <- sapply(orientations_prob[, c(2:3)], as.numeric)
    orientations_prob[, c(5:6)] <- sapply(orientations_prob[, c(5:6)], as.numeric)
    orientations_prob[, c(8:9)] <- sapply(orientations_prob[, c(8:9)], as.numeric)
  }
  # update the returned matrix
  res$triples <- orientations_prob

  res$interrupted <- FALSE

  res
}
