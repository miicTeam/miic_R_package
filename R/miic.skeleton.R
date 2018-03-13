
miic.skeleton <- function(inputData = NULL, blackBox = NULL, stateOrder = NULL, nThreads= nThreads, effN = -1,
                               cplx = c("nml", "mdl"), eta = 1, latent = FALSE, continuous = 0,
                               confidenceShuffle = 0, confidenceThreshold = 0, verbose = FALSE
){

  isTplReuse = TRUE
  isK23 = TRUE
  isDegeneracy = FALSE
  isNoInitEta = FALSE

  if( is.null( inputData ) )
  { stop("The input data file is required") }

  if(length(cplx) == 2)
    cplx = "nml"

  if(cplx != "nml" && cplx != "mdl")
  { stop("The complexity method is not allowed") }

  numNodes <- length(inputData)

  inData <- c(colnames(inputData), as.vector(as.character(t(as.matrix(inputData)))))
  if(!is.null(blackBox)){
    bB <- as.vector(as.character(t(as.matrix(blackBox))))
  } else {
    bB = c("")
  }

  if(!is.null(stateOrder)){
    stateOrder <- as.vector(as.character(t(as.matrix(stateOrder))))
  } else {
    stateOrder = c("")
  }

  if (base::requireNamespace("Rcpp", quietly = TRUE)) {
      res <- .Call('skeleton', inData, numNodes, nThreads, bB, effN, cplx, eta, latent, isTplReuse,
               isK23, isDegeneracy, isNoInitEta, continuous, confidenceShuffle, confidenceThreshold, verbose, PACKAGE = "miic")
  }

  # if(shuffle == 0) {
    # create the data frame of the edges
    tmp <- unlist(res$edges)[1:10]
    res1 <- unlist(res$edges)[11:length(unlist(res$edges))]
    df <- data.frame(matrix(res1, nrow=length(res$edges)-1, byrow=T),stringsAsFactors=FALSE)
    row.names(df) <- df[,1]
    df <- df[,-1]
    colnames(df) <-tmp
    df[ df == "NA" ] = NA
    df[,c(6:10)] = sapply(df[,c(6:10)], as.numeric)

    # update the returned object
    res$edges <- df

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

  if(confidenceShuffle > 0){
    # create the data frame for the confidence file
    tmp <- unlist(res$confData)[1:length(res$confData[[1]])]
    res1 <- unlist(res$confData)[(length(res$confData[[1]])+1):length(unlist(res$confData))]
    df <- data.frame(matrix(res1, nrow=length(res$confData)-1, byrow=T),stringsAsFactors=FALSE)
    colnames(df) <-tmp
    df[,"confidence_ratio"] = as.numeric(df[,"confidence_ratio"])

    df <-df[order(df[,"confidence_ratio"]),]


    res$confData <- df

  }

  # save time
  time = strsplit(as.character(res$time)," ")
  time[which(time == 0)]=NA

  res$time <- stats::setNames(as.numeric(time),c("init", "iter", "initIter", "cut"))



  res
}
