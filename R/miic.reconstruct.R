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
  var_names <- colnames(input_data)

  if (!is.null(black_box)) {
    black_box <- as.vector(as.character(t(as.matrix(black_box))))
  } else {
    black_box <- character()  # pass to cpp as empty vector<string>
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
    "max_iteration" = max_iteration,
    "n_eff" = n_eff,
    "n_shuffles" = n_shuffles,
    "n_threads" = n_threads,
    "no_init_eta" = FALSE,
    "orientation" = orientation,
    "propagation" = propagation,
    "sample_weights" = sample_weights,
    "test_mar" = test_mar,
    "var_names" = var_names,
    "verbose" = verbose
  )
  if (base::requireNamespace("Rcpp", quietly = TRUE)) {
    cpp_input <- data.frame(t(sapply(input_data, as.character)),
                            stringsAsFactors = FALSE)
    # Call C++ function
    res <- reconstruct(cpp_input, arg_list)
    if (res$interrupted) {
      return(list(interrupted = TRUE))
    }
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

  if (n_shuffles > 0) {
    # create the data frame for the confidence file
    confData <- res$confData
    a <- (length(confData[[1]]) + 1)
    b <- length(unlist(confData))
    tmp <- unlist(confData)[1:length(confData[[1]])]
    res1 <- unlist(confData)[a:b]
    df <-
      data.frame(matrix(res1, nrow = length(confData) - 1, byrow = T),
        stringsAsFactors = FALSE
      )
    colnames(df) <- tmp
    df[, "confidence_ratio"] <- as.numeric(df[, "confidence_ratio"])

    df <- df[order(df[, "confidence_ratio"]), ]

    res$confData <- df
  }

  # save time
  time <- strsplit(as.character(res$time), " ")
  time[which(time == 0)] <- NA

  res$time <- stats::setNames(
    as.numeric(time),
    c("init", "iter", "cut", "skeleton")
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
      byrow = T
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
