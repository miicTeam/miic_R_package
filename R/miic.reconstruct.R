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
                             typeOfData = NULL,
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

  numNodes <- length(inputData)

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
      typeOfData,
      cntVar,
      numNodes,
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

  # if(shuffle == 0) {
  # create the data frame of the edges
  tmp <- unlist(res$edges)[1:10]
  res1 <- unlist(res$edges)[11:length(unlist(res$edges))]
  df <-
    data.frame(matrix(res1, nrow = length(res$edges) - 1, byrow = T),
      stringsAsFactors = FALSE
    )
  row.names(df) <- df[, 1]
  df <- df[, -1]
  colnames(df) <- tmp
  df[df == "NA"] <- NA
  df[, c(6:10)] <- sapply(df[, c(6:10)], as.numeric)
  # update the returned object
  res$edges <- df
  # create the data frame of the adj matrix
  a <- (length(res$adjMatrix[[1]]) + 1)
  b <- length(unlist(res$adjMatrix))
  tmp <- unlist(res$adjMatrix)[1:length(res$adjMatrix[[1]])]
  res1 <- unlist(res$adjMatrix)[a:b]
  df <- data.frame(matrix(res1,
    nrow = length(res$adjMatrix) - 1,
    byrow = T
  ),
  stringsAsFactors = FALSE
  )
  df <- df[, -1]
  colnames(df) <- tmp
  df <- sapply(df, as.numeric)
  row.names(df) <- tmp
  # update the returned adj matrix
  res$adjMatrix <- df
  # } else {
  #   # create the data frame of the adj matrix
  #   tmp <- unlist(res$edges.shuffle)[1:3]
  #   res1 <- unlist(res$edges.shuffle)[(length(res$edges.shuffle[[1]])+1):length(unlist(res$edges.shuffle))]
  #   df <- data.frame(matrix(res1, nrow=length(res$edges.shuffle)-1, byrow=T),stringsAsFactors=FALSE)
  #   colnames(df) <-tmp

  #   # update the returned adj matrix
  #   res$edges.shuffle <- df

  #   tmp <- unlist(res$edges.mean)[1:2]
  #   res1 <- unlist(res$edges.mean)[(length(res$edges.mean[[1]])+1):length(unlist(res$edges.mean))]
  #   df <- data.frame(matrix(res1, nrow=length(res$edges.mean)-1, byrow=T),stringsAsFactors=FALSE)
  #   colnames(df) <-tmp

  #   # update the returned adj matrix
  #   res$edges.mean <- df
  # }

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
    c("init", "iter", "initIter", "cut")
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
