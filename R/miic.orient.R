miic.orient <- function(inputData = NULL, cntVar = NULL, stateOrder = NULL,
                        edges = NULL, effN = -1, cplx = "nml", eta = 1,
                        latent = FALSE, propagation = TRUE, typeOfData = NULL,
                        sampleWeights = NULL, hvs = FALSE, continuous = 0,
                        verbose = FALSE) {
  isK23 = TRUE
  isDegeneracy = FALSE

  if (hvs == FALSE)
    hvs = 0
  else if (hvs == TRUE)
    hvs = 1
  else
    stop("Half v-structures must be TRUE or FALSE")

  if (is.null(inputData))
    stop("The input data file is required")

  if (is.null(edges))
    stop("The edge file path is required")

  if (cplx != "nml" && cplx != "mdl")
    stop("The complexity method is not allowed")

  # Set the script directory
  # rootDir = file.path( "~")

  numNodes <- length(inputData)
  edges <- as.vector(as.character(t(as.matrix(edges))))

  inData <- c(colnames(inputData),
              as.vector(as.character(t(as.matrix(inputData)))))


  if (!is.null(stateOrder)) {
    stateOrder <- as.vector(as.character(t(as.matrix(stateOrder))))
  } else {
    stateOrder = c("")
  }

  if (is.null(sampleWeights)) {
    sampleWeights = c(-1, rep(0, nrow(inputData) - 1))
  }
  cntVar = as.numeric(cntVar)
  res <- .Call('orientationProbability', inData, typeOfData, cntVar, numNodes,
               edges, effN, cplx, eta, hvs, latent, isK23, isDegeneracy,
               propagation, sampleWeights, verbose)

  #create the data frame of the structures after orientation
  df = res$orientations.prob
  if (length(res$orientations.prob) > 0) {
    tmp <- unlist(res$orientations.prob)[1:length(res$orientations.prob[[1]])]
    res1 <- unlist(res$orientations.prob)[(length(res$orientations.prob[[1]]) + 1):length(unlist(res$orientations.prob))]
    df <- data.frame(matrix(res1, nrow = length(res$orientations.prob) - 1, byrow = T),
                     stringsAsFactors = FALSE)
    colnames(df) <-tmp

    df[, c(2:3)] = sapply(df[, c(2:3)], as.numeric)
    df[, c(5:6)] = sapply(df[, c(5:6)], as.numeric)
    df[, c(8:9)] = sapply(df[, c(8:9)], as.numeric)
  }

  # update the returned matrix
  res$orientations.prob <- df

  # create the data frame of the adj matrix
  tmp <- unlist(res$adjMatrix)[1:length(res$adjMatrix[[1]])]
  res1 <- unlist(res$adjMatrix)[(length(res$adjMatrix[[1]]) + 1):length(unlist(res$adjMatrix))]
  df <- data.frame(matrix(res1, nrow = length(res$adjMatrix) - 1, byrow = T),
                   stringsAsFactors = FALSE)
  df <- df[, -1]
  colnames(df) <-tmp
  df = sapply(df, as.numeric)
  row.names(df) <- tmp

  # update the returned adj matrix
  res$adjMatrix <- df

  res
}
