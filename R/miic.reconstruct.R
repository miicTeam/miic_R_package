miic.reconstruct <- function(input_data = NULL,
                             is_continuous = NULL,
                             black_box = NULL,
                             n_threads = 1,
                             n_eff = -1,
                             cplx = "nml",
                             eta = 1,
                             latent = "no",
                             n_shuffles = 0,
                             orientation = TRUE,
                             ori_proba_ratio = 1,
                             propagation = TRUE,
                             conf_threshold = 0,
                             verbose = FALSE,
                             sample_weights = NULL,
                             test_mar = TRUE,
                             consistent = "no",
                             max_iteration = NULL
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
    "ori_proba_ratio" = ori_proba_ratio,
    "propagation" = propagation,
    "test_mar" = test_mar,
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
  if (!is.null(sample_weights)) {
    arg_list[["sample_weights"]] <- sample_weights
  }

  cpp_input <- list("factor" = input_factor, "double" = input_double,
                    "order" = input_order)
  # Call C++ function
  res <- reconstruct(cpp_input, arg_list)
  if (res$interrupted) {
    return(list(interrupted = TRUE))
  }

  # R-formalize returned object
  # table of edges infomation
  n_row <- length(res$edges) - 1
  header <- unlist(res$edges[1])
  df <- data.frame(matrix(unlist(res$edges[2:(n_row + 1)]), nrow = n_row,
                          byrow = TRUE), stringsAsFactors = FALSE)
  colnames(df) <- header
  df[df == "NA"] <- NA
  df$Ixy <- as.numeric(df$Ixy)
  df$Ixy_ai <- as.numeric(df$Ixy_ai)
  df$cplx <- as.numeric(df$cplx)
  df$Rxyz_ai <- as.numeric(df$Rxyz_ai)
  df$Nxy_ai <- as.numeric(df$Nxy_ai)
  df$confidence <- as.numeric(df$confidence)
  res$edges <- df

  #  adj_matrix
  res$adj_matrix <- matrix(unlist(res$adj_matrix), nrow = length(input_data), byrow = TRUE)
  colnames(res$adj_matrix) <- var_names
  rownames(res$adj_matrix) <- var_names

  # adj_matrices (when consistent parameter is turned on)
  if (length(res$adj_matrices) > 0) {
    res$adj_matrices <- matrix(unlist(res$adj_matrices),
                               ncol = length(res$adj_matrices))
  }

  # save time
  time <- strsplit(as.character(res$time), " ")
  time[which(time == 0)] <- NA

  res$time <- stats::setNames(
    as.numeric(time),
    c("init", "iter", "cut", "ori", "total")
  )

  # create the data frame of the structures after orientation
  orientations_prob <- res$orientations.prob

  if (length(res$orientations.prob) > 0) {
    a <- length(orientations_prob[[1]])
    b <- length(unlist(orientations_prob))
    tmp <- unlist(res$orientations.prob)[1:a]
    res1 <- unlist(res$orientations.prob)[(a + 1):b]
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
  res$orientations.prob <- orientations_prob

  res$interrupted <- FALSE

  res
}
