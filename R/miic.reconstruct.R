miic.reconstruct <- function (list_in)
  {
  print (head (list_in$input_data))
  n_samples <- nrow (list_in$input_data)
  n_nodes <- ncol (list_in$input_data)
  print (n_samples)
  print (n_nodes)
  #
  # Convert discrete vars as factors
  #
  for ( i in 1:nrow(list_in$state_order) )
    if (list_in$state_order[i, "var_type"] == 0)
      list_in$input_data[, i] <- factor (list_in$input_data[, i])
  print ("convers factor done")
  #
  # Numeric factor matrix, level starts from 0, NA mapped to -1
  #
  input_factor <- apply(list_in$input_data, 2, function(x)
                        (as.numeric(factor(x, levels = unique(x))) - 1))
  print ("input factor base done")
  input_factor[is.na(input_factor)] <- -1
  max_level_list <- as.numeric(apply(input_factor, 2, max)) + 1
  input_factor <- as.vector(as.matrix(input_factor))
  #
  # Data list, numeric for continuous columns, -1 for discrete columns
  #
  input_double <- matrix(nrow = n_samples, ncol = n_nodes)
  #
  # Order list, order(column) for continuous columns (index starting from 0,
  # NA mapped to -1), -1 for discrete columns
  #
  input_order <- matrix(nrow = n_samples, ncol = n_nodes)
  for (i in c(1:ncol(list_in$input_data)))
    {
    if (list_in$state_order[i,"var_type"] == 1)
      {
      input_double[, i] <- as.numeric(list_in$input_data[, i])
      n_NAs <- sum(is.na(list_in$input_data[, i]))
      input_order[, i] <- c(order(list_in$input_data[, i], na.last=NA) - 1,
                            rep_len(-1, n_NAs))
      }
    else
      {
      input_double[, i] <- rep_len(-1, n_samples)
      input_order[, i] <- rep_len(-1, n_samples)
      }
    }
  input_order <- as.vector(input_order)
  input_double <- as.vector(input_double)

  var_names <- colnames(list_in$input_data)
  n_vars <- length (var_names)
  #
  # Parameters always supplied to C++ reconstruct
  #
  arg_list <- list(
    #
    # Parameters coming from miic parameters (can be defined by the user)
    #
    "n_threads" = list_in$params$n_threads,
    "cplx" = list_in$params$cplx,
    "orientation" = list_in$params$orientation,
    "ort_proba_ratio" = list_in$params$ort_proba_ratio,
    "propagation" = list_in$params$propagation,
    "latent" = list_in$params$latent,
    "n_eff" = list_in$params$n_eff,
    "n_shuffles" = list_in$params$n_shuffles,
    "conf_threshold" = list_in$params$conf_threshold,
    "test_mar" = list_in$params$test_mar,
    "consistent" = list_in$params$consistent,
    "max_iteration" = list_in$params$max_iteration,
    "negative_info" = list_in$params$negative_info,
    "mode" = list_in$params$mode,
    "verbose" = list_in$params$verbose,
    #
    # Parameters part of miic state order (can be defined by the user)
    #
    "is_continuous" = as.numeric (list_in$state_order$var_type),
    #
    # Parameters deduced from data
    #
    "var_names" = var_names,
    "n_nodes" = n_nodes,
    "n_samples" = n_samples,
    "levels" = max_level_list,
    "max_bins" = min(50, n_samples),
    #
    # Parameters fixed to a default value, not supplied by the user to miic
    #
    "degenerate" = FALSE,
    "half_v_structure" = 0,
    "no_init_eta" = FALSE
    )
  #
  # Optional parameters
  #
  if (!is.null(list_in$black_box))
    {
    # transform var names to var indices
    black_box <- list_in$black_box
    black_box[] <- sapply(black_box, function(x) {
      match(as.character(x), colnames(list_in$input_data)) - 1 } )
    black_box[] <- black_box[stats::complete.cases(black_box),]
    arg_list[["black_box"]] <- as.vector(as.matrix(t(black_box)))
    }
  if (!is.null(list_in$params$sample_weights))
    arg_list[["sample_weights"]] <- list_in$params$sample_weights
  if (!is.null(list_in$state_order$is_contextual))
    arg_list[["is_contextual"]] <- list_in$state_order$is_contextual
  if (!is.null(list_in$state_order$is_consequence))
    arg_list[["is_consequence"]] <- list_in$state_order$is_consequence
  if (!is.null(list_in$state_order$n_layers))
    arg_list[["n_layers"]] <- list_in$state_order$n_layers
  if (!is.null(list_in$state_order$delta_t))
    arg_list[["delta_t"]] <- list_in$state_order$delta_t
  if (!is.null(list_in$state_order$layers))
    # # -1 for C++ indices starting at 0
    arg_list[["nodes_layers"]] <- list_in$state_order$layers
  if (!is.null(list_in$layers))
    # {
    # -1 for C++ indices starting at 0
    # list_in$layers$contributors = unlist (lapply (list_in$layers$contributors,
    #   FUN=function(x) { ifelse (is.na (x),
    #     NA_character_,
    #     paste (as.integer (strsplit (x, ";", fixed=T)[[1]]) - 1, collapse=";") )
    #   } ) )
    # FALSE/TRUE transformed into 0/1 for C++
    # list_in$layers$preoriented = unlist (lapply (list_in$layers$preoriented,
    #   FUN=function(x) { ifelse (is.na (x),
    #     NA_character_,
    #     paste (as.integer (as.logical (strsplit (x, ";", fixed=T)[[1]]) ), collapse=";") )
    #   } ) )
    arg_list[["layers"]] <- list_in$layers
    # }

  cpp_input <- list("factor" = input_factor, "double" = input_double,
                    "order" = input_order)
  # Call C++ function
  res <- reconstruct(cpp_input, arg_list)
  if (res$interrupted)
    return(list(interrupted = TRUE))
  #
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
  #
  # Execution time
  #
  time <- strsplit(as.character(res$time), " ")
  time[which(time == 0)] <- NA
  res$time <- stats::setNames(
    as.numeric(time),
    c("init", "iter", "cut", "ort", "cpp")
  )
  #
  # Create the data frame of the structures after orientation
  #
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
  return (res)
  }
