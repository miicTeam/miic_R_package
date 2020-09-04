miic.reconstruct <- function(input_data = NULL,
                             is_continuous = NULL,
                             black_box = NULL,
                             n_threads = n_threads,
                             n_eff = -1,
                             cplx = c("nml", "mdl"),
                             eta = 1,
                             latent = c("no", "yes", "orientation"),
                             n_shuffles = 0,
                             orientation = TRUE,
                             ori_proba_ratio = 1,
                             propagation = TRUE,
                             conf_threshold = 0,
                             verbose = FALSE,
                             sample_weights = NULL,
                             test_mar = TRUE,
                             consistent = c(
                               "no",
                               "orientation",
                               "skeleton"
                             ),
                             max_iteration = max_iteration
                             ) {
  # Numeric factor matrix, level starts from 0, NA mapped to -1
  input_factor <- apply(input_data, 2, function(x)
                        (as.numeric(factor(x, levels = unique(x))) - 1))
  input_factor[is.na(input_factor)] <- -1
  max_level_list <- as.numeric(apply(input_factor, 2, max)) + 1
  input_factor <- data.frame(t(input_factor))
  # Data list, numeric for continuous columns, empty for others
  input_double <- list()
  # Order list, order(column) for continuous columns (index starting from 0, NA
  # mapped to -1), empty for others
  input_order <- list()
  for (i in c(1:length(input_data))) {
    if (is_continuous[i]) {
      input_double[[i]] <- as.numeric(input_data[, i])
      n_NAs <- sum(is.na(input_data[, i]))
      input_order[[i]] <- c(order(input_data[, i], na.last=NA) - 1,
                            rep_len(-1, n_NAs))
    } else {
      input_double[[i]] <- numeric(0)
      input_order[[i]] <- numeric(0)
    }
  }

  var_names <- colnames(input_data)
  nameToIndex <- function(name, var_indices) {
    index <- var_indices[name]
    if (is.na(index))
      stop(paste0(name, " in the black box does not match names in the input"))
    else
      return(index)
  }
  if (!is.null(black_box)) {
    # transform var names to var indices
    var_indices <- c(0: (length(input_data) - 1))
    names(var_indices) <- var_names
    black_box[] <- lapply(black_box, nameToIndex, var_indices)
    black_box <- data.frame(t(black_box))
  } else {
    black_box <- list()  # pass to cpp as empty vector
  }

  if (is.null(sample_weights)) {
    sample_weights <- numeric(0)  # pass to cpp as empty vector<double>
  }

  arg_list <- list(
    "black_box" = black_box,
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
    "n_shuffles" = n_shuffles,
    "n_threads" = n_threads,
    "no_init_eta" = FALSE,
    "orientation" = orientation,
    "ori_proba_ratio" = ori_proba_ratio,
    "propagation" = propagation,
    "sample_weights" = sample_weights,
    "test_mar" = test_mar,
    "max_bins" = 50,
    "var_names" = var_names,
    "verbose" = verbose
  )
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
