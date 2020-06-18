miic.reconstruct <- function(inputData = NULL,
                             cntVar = NULL,
                             blackBox = NULL,
                             stateOrder = NULL,
                             nThreads = nThreads,
                             effN = -1,
                             cplx = c("nml", "mdl"),
                             eta = 1,
                             latent = c("no", "yes", "ort"),
                             confidenceShuffle = 0,
                             edges = NULL,
                             orientation = TRUE,
                             propagation = TRUE,
                             confidenceThreshold = 0,
                             verbose = FALSE,
                             sampleWeights = NULL,
                             testMAR = TRUE,
                             consistent = c(
                               "no",
                               "orientation",
                               "skeleton"
                             )) {
  isTplReuse <- TRUE
  isK23 <- TRUE
  isDegeneracy <- FALSE
  isNoInitEta <- FALSE

  n_node <- length(inputData)

  inData <- c(
    colnames(inputData),
    as.vector(as.character(t(
      as.matrix(inputData)
    )))
  )
  if (!is.null(blackBox)) {
    bB <- as.vector(as.character(t(as.matrix(blackBox))))
  } else {
    bB <- c("")
  }

  if (!is.null(edges)) {
    edges <- as.vector(as.character(t(as.matrix(blackBox))))
  } else {
    edges <- c("")
  }

  if (!is.null(stateOrder)) {
    stateOrder <- as.vector(as.character(t(as.matrix(stateOrder))))
  } else {
    stateOrder <- c("")
  }

  if (is.null(sampleWeights)) {
    sampleWeights <- c(-1, rep(0, nrow(inputData) - 1))
  }
  hvs <- 0
  cntVar <- as.numeric(cntVar)
  if (base::requireNamespace("Rcpp", quietly = TRUE)) {
    res <- .Call(
      "reconstruct",
      inData,
      cntVar,
      n_node,
      nThreads,
      edges,
      bB,
      effN,
      cplx,
      eta,
      hvs,
      latent,
      isTplReuse,
      isK23,
      isDegeneracy,
      orientation,
      propagation,
      isNoInitEta,
      confidenceShuffle,
      confidenceThreshold,
      sampleWeights,
      consistent,
      testMAR,
      verbose,
      PACKAGE = "miic"
    )
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
  res$adj_matrix <- matrix(unlist(res$adj_matrix), nrow = n_node, byrow = TRUE)
  colnames(res$adj_matrix) <- colnames(inputData)
  rownames(res$adj_matrix) <- colnames(inputData)

  # adj_matrices (when consistent parameter is turned on)
  if (length(res$adj_matrices) > 0) {
    res$adj_matrices <- matrix(unlist(res$adj_matrices),
                               ncol = length(res$adj_matrices))
  }

  if (confidenceShuffle > 0) {
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
