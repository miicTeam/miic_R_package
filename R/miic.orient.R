
miic.orient <- function(inputData = NULL, method = c("probabilistic", "logic"), stateOrder = NULL, edges = NULL,
                             effN = -1, cplx = c("nml", "mdl"), eta = 1,
                             latent = FALSE, propagation = TRUE, hvs = FALSE, continuous = 0, verbose = FALSE)
{
  isK23 = TRUE
  isDegeneracy = FALSE

  if(hvs == FALSE)
    hvs = 0
  else if(hvs == TRUE)
    hvs = 1
  else
    stop("Half v-structures must be TRUE or FALSE")

  if( is.null( inputData ) )
  { stop("The input data file is required") }

  if( is.null( edges ) )
  { stop("The edge file path is required") }

  if(length(cplx) == 2)
    cplx = "nml"

  if(cplx != "nml" && cplx != "mdl")
  { stop("The complexity method is not allowed") }

  if(length(method) == 2)
    method = "probabilistic"

  if(method != "probabilistic" && method != "logic")
  { stop("This orientation method is not allowed") }

  # Set the script directory
  # rootDir = file.path( "~")

  numNodes <- length(inputData)
  edges <- as.vector(as.character(t(as.matrix(edges))))

  inData <- c(colnames(inputData), as.vector(as.character(t(as.matrix(inputData)))))


  if(!is.null(stateOrder)){
    stateOrder <- as.vector(as.character(t(as.matrix(stateOrder))))
  } else {
    stateOrder = c("")
  }

  if(method=="logic" ){

    res <- .Call('orientation',inData, numNodes, edges, effN, cplx, eta, latent,
                isK23, isDegeneracy, propagation, continuous, verbose, PACKAGE = "miic")

    # create the data frame of the structures before orientation
    tmp <- unlist(res$tableOfOrientationsBeforePropagation)[1:length(res$tableOfOrientationsBeforePropagation[[1]])]
    res1 <- unlist(res$tableOfOrientationsBeforePropagation)[(length(res$tableOfOrientationsBeforePropagation[[1]])+1):length(unlist(res$tableOfOrientationsBeforePropagation))]
    df <- data.frame(matrix(res1, nrow=length(res$tableOfOrientationsBeforePropagation)-1, byrow=T),stringsAsFactors=FALSE)
    df <- df[,-1]
    colnames(df) <-tmp
    df[ df == "NA" ] = NA
    df[,c(2)] = sapply(df[,c(2)], as.numeric)
    df[,c(4)] = sapply(df[,c(4)], as.numeric)
    df[,c(6:13)] = sapply(df[,c(6:13)], as.numeric)
    # update the returned matrix
    res$tableOfOrientationsBeforePropagation <- df


    #create the data frame of the structures after orientation
    tmp <- unlist(res$tableOfOrientationsAfterPropagation)[1:length(res$tableOfOrientationsAfterPropagation[[1]])]
    res1 <- unlist(res$tableOfOrientationsAfterPropagation)[(length(res$tableOfOrientationsAfterPropagation[[1]])+1):length(unlist(res$tableOfOrientationsAfterPropagation))]
    df <- data.frame(matrix(res1, nrow=length(res$tableOfOrientationsAfterPropagation)-1, byrow=T),stringsAsFactors=FALSE)
    df <- df[,-1]
    colnames(df) <-tmp
    df[ df == "NA" ] = NA
    df[,c(2)] = sapply(df[,c(2)], as.numeric)
    df[,c(4)] = sapply(df[,c(4)], as.numeric)
    df[,c(6:13)] = sapply(df[,c(6:13)], as.numeric)

    # update the returned matrix
    res$tableOfOrientationsAfterPropagation <- df


  } else if( method=="probabilistic" ){

    res <- .Call('orientationProbability',inData, numNodes, edges, effN, cplx, eta, hvs, latent,
                 isK23, isDegeneracy, propagation, continuous, verbose)

    #create the data frame of the structures after orientation
    df = res$orientations.prob
    if(length(res$orientations.prob) > 0)
    {
      tmp <- unlist(res$orientations.prob)[1:length(res$orientations.prob[[1]])]
      res1 <- unlist(res$orientations.prob)[(length(res$orientations.prob[[1]])+1):length(unlist(res$orientations.prob))]
      df <- data.frame(matrix(res1, nrow=length(res$orientations.prob)-1, byrow=T),stringsAsFactors=FALSE)
      colnames(df) <-tmp

      df[,c(2:3)] = sapply(df[,c(2:3)], as.numeric)
      df[,c(5:6)] = sapply(df[,c(5:6)], as.numeric)
      df[,c(8:9)] = sapply(df[,c(8:9)], as.numeric)
    }

    # update the returned matrix
    res$orientations.prob <- df
  }

  # create the data frame of the adj matrix

  tmp <- unlist(res$adjMatrix)[1:length(res$adjMatrix[[1]])]
  res1 <- unlist(res$adjMatrix)[(length(res$adjMatrix[[1]])+1):length(unlist(res$adjMatrix))]
  df <- data.frame(matrix(res1, nrow=length(res$adjMatrix)-1, byrow=T),stringsAsFactors=FALSE)
  df <- df[,-1]
  colnames(df) <-tmp
  df = sapply(df, as.numeric)
  row.names(df) <- tmp

  # update the returned adj matrix
  res$adjMatrix <- df

  res
}
