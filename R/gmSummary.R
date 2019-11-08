gmSummary <- function(inputData=NULL, edges = NULL, trueEdges = NULL,
                      stateOrder = NULL, adjMatrix = NULL, verbose = FALSE)
{

  if( is.null( inputData ) )
  { stop("The input data file is required") }

  if( is.null( edges ) )
  { stop("The inferred edges data frame is required") }

  if( is.null( adjMatrix ) )
  { stop("The adjacency matrix data frame is required") }

  startTime.whole <- proc.time()

  correlationToEvaluate = TRUE
  if(length(stateOrder) == 0 ){
    if(length(which(apply(inputData, 2, is.numeric) == FALSE)) != 0){
      correlationToEvaluate = FALSE
    }
  }


  if(!is.null(trueEdges))
    colnames(trueEdges)=c("x","y")

  #### Set an env variable
  gV <- new.env()

  #### Load fist the adjacency matrix to get all the variables name
  adjMat <- adjMatrix
  gV$allProperties = colnames(adjMat)

  #### Load also the  inferred edges with their mutual information
  allMutInfo.df <- edges

  #### Check if the orientations have been learned with assuming latent variables
  gV$isLatent = FALSE
  if( max(adjMat) > 2 | min(adjMat) < -2 ){ gV$isLatent = TRUE }

  #### Load the list of the true edges
  true.edgesList = NULL
  if( length(trueEdges) > 0 )
  {
    if( !is.null(trueEdges))
    { true.edgesList <- trueEdges }

    #### Set as rowname the key of each edge
    rownames( true.edgesList ) <- apply( true.edgesList, MARGIN = 1 , function(myRow){ binVectToStr( myRow[1:2], gV ) } )
  }

  #### Convert the adj mat it to upper triangular matrix
  adjMat[ lower.tri( adjMat, diag = TRUE ) ] = 0

  #### Get all inferred edges and gather their idx and val into one list
  edgesNbr = 0
  edges.list = list( nOrt = list( idx = data.frame( nrow = numeric(0), ncol = numeric(0) ), val = numeric(0), key = character(0) )
                     , ort = list( idx = data.frame( nrow = numeric(0), ncol = numeric(0) ), val = numeric(0), key = character(0) )
                     , ph = list( idx = data.frame( nrow = numeric(0), ncol = numeric(0) ), val = numeric(0), key = character(0) ) )

  #### ----
  # -- not oriented
  edges.list[["nOrt"]][["idx"]] <- which( adjMat == 1, arr.ind = TRUE )
  edges.list[["nOrt"]][["val"]] <- rep( 1, nrow(edges.list[["nOrt"]][["idx"]]) )
  edgesNbr = edgesNbr + length( edges.list[["nOrt"]][["val"]] )

  # -- oriented
  #### get the number of oriented edges
  ortN = ( nrow( which( abs( adjMat ) >= 2, arr.ind = TRUE ) ) )

  if( ortN > 0 )
  {
    #### Initialize the dimension array of index with the nbr of oriented edges
    edges.list[["ort"]][["idx"]] <- data.frame( nrow = numeric(ortN), ncol = numeric(ortN) )
    edges.list[["ort"]][["val"]] <- rep(NA, ortN)

    countRow = 0
    for( iOrt in c(2,-2,4,-4,6) )
    {
      #### Get the array of index
      tmp.ortArr.idx <- which( adjMat == iOrt, arr.ind = TRUE )

      #### Get the number of edges
      tmp.ortN = nrow( tmp.ortArr.idx )
      tmp.ortArr.val <- rep( iOrt, (tmp.ortN) )

      #### Insert in the initalised matrix
      if( tmp.ortN > 0 )
      {
        edges.list[["ort"]][["idx"]][(countRow+1):(countRow+tmp.ortN),] <- tmp.ortArr.idx
        edges.list[["ort"]][["val"]][(countRow+1):(countRow+tmp.ortN)] <- tmp.ortArr.val

        countRow = (countRow + tmp.ortN)
      }
    }
  }
  edgesNbr = edgesNbr + length( edges.list[["ort"]][["val"]] )

  # -- phantom
  edges.list[["ph"]][["idx"]] <- which( adjMat == 0, arr.ind = TRUE )

  if( nrow( edges.list[["ph"]][["idx"]] ) > 0 )
  {
    edges.list[["ph"]][["idx"]] <- edges.list[["ph"]][["idx"]][which( edges.list[["ph"]][["idx"]][,2] > edges.list[["ph"]][["idx"]][,1] ),, drop=FALSE]
    edges.list[["ph"]][["val"]] <- rep( 0, nrow( edges.list[["ph"]][["idx"]] ) )
  }
  edgesNbr = edgesNbr + length( edges.list[["ph"]][["val"]] )

  #### Initialize the output data frame
  ### LV // adding 2 columns: signs & partial_correlation
  outputSummary.df <- data.frame( x = character(edgesNbr), y = character(edgesNbr)
                                  , type = character(edgesNbr), ai = character(edgesNbr)
                                  , info = numeric(edgesNbr) , info_cond = numeric(edgesNbr)
                                  , cplx = numeric(edgesNbr), Nxy_ai = numeric(edgesNbr)
                                  , log_confidence = numeric(edgesNbr)
                                  , infOrt = numeric(edgesNbr), trueOrt = numeric(edgesNbr)
                                  , isOrt = character(edgesNbr), isOrtOk = character(edgesNbr)
                                  , essential = character(edgesNbr)
                                  , sign = character(edgesNbr)
                                  , partial_correlation = character(edgesNbr)
                                  , is_causal = character(edgesNbr)
                                  , proba = character(edgesNbr)
                                  , stringsAsFactors = FALSE )

  #### ----
  countRow <- 0
  for( iEdgeCategory in names( edges.list ) )
  {
    if( nrow( edges.list[[iEdgeCategory]][["idx"]] ) > 0 )
    {
      #### Replace numbers by names
      edges.list[[iEdgeCategory]][["idx"]] = t(apply( edges.list[[iEdgeCategory]][["idx"]], MARGIN = c(1), FUN = function(x) { x = gV$allProperties[x] } ))

      #### Compute the key of each edge
      edges.list[[iEdgeCategory]][["key"]] = t(apply( edges.list[[iEdgeCategory]][["idx"]], MARGIN = c(1), FUN = function(x) { x = binVectToStr( x, gV ) } ))

      #### Set the bounds for insertion
      currentN = nrow(edges.list[[iEdgeCategory]][["idx"]])
      lastRow = ( countRow + currentN )

      #### Insert the new data
      outputSummary.df[c((countRow+1):lastRow), ] <- data.frame(
        x = as.character( as.vector( edges.list[[iEdgeCategory]][["idx"]][,1] ) )
        , y = as.character( as.vector( edges.list[[iEdgeCategory]][["idx"]][,2] ) )
        , type = rep( NA, currentN )
        , ai = rep( NA, currentN )
        , info = rep( NA, currentN )
        , cplx = rep( NA, currentN )
        , Nxy_ai = rep( NA, currentN )
        , log_confidence = rep( NA, currentN )
        , infOrt = edges.list[[iEdgeCategory]][["val"]]
        , trueOrt = rep( NA, currentN )
        , isOrt = rep( NA, currentN ), isOrtOk = rep( NA, currentN )
        , essential = rep( NA, currentN )
        , sign = rep( NA, currentN )
        , partial_correlation = rep( NA, currentN )
        , is_causal = rep( NA, currentN )
        , proba = rep( NA, currentN )
        , stringsAsFactors = FALSE )

      #### Update the total row count
      countRow = lastRow
    }
  }

  #### Set the key of each edge as rowname
  rownames( outputSummary.df ) <- do.call(c, list(edges.list$nOrt$key, edges.list$ort$key, edges.list$ph$key))

  #### Get the true links and check the TP, FP, FN, TN
  if( !is.null(true.edgesList) )
  {
    #### Set the name of the columns for the true edges list
    colnames( true.edgesList )[1:2] = c( "x", "y" )
    if( ncol(true.edgesList) == 3 ){ colnames( true.edgesList )[3] = "essential" }

    #### Distinguish the TP from the FP among the non phantom links
    allEdge.idx <- which( abs( outputSummary.df[,"infOrt"] ) > 0 )

    # [TP]
    TP.type.idx = c()
    if( length(allEdge.idx) > 0 )
    {
      #### First, initialise all inferred edge as FP
      outputSummary.df[allEdge.idx, "type"] <- 'FP'

      #### Then, find the TP
      #### ie., among the inferred edges, which are in the true edges list
      TP.type.idx <- which( abs( outputSummary.df[,"infOrt"] ) > 0 & ( rownames( outputSummary.df ) %in% rownames( true.edgesList ) ) )

      if( length(TP.type.idx) > 0 ) { outputSummary.df[TP.type.idx, "type"] <- 'TP' }
    }
    #### ----

    #### Distinguish the TN from the FN among the phantom links
    allPh.idx <- which( outputSummary.df[,"infOrt"] == 0 )

    # [FN]
    FN.type.idx = c()
    if( length(allPh.idx) > 0 )
    {
      #### First, initialise all inferred phantom as TN
      outputSummary.df[allPh.idx, "type"] <- 'TN'

      #### Then, find the FN
      #### ie., among the phantom, which have a key found in the true edges list
      FN.type.idx <- which( rownames( true.edgesList ) %in% rownames(outputSummary.df)[allPh.idx] )

      if( length(FN.type.idx) > 0 ) { outputSummary.df[rownames( true.edgesList )[FN.type.idx], "type"] <- 'FN' }
    }
    #### ----

    #### Add the true orientations into the summary
    #### ----
    if( length(FN.type.idx) > 0 )
    {
      #### For the FN (forward = 2, backward = -2)
      FN.ort <- sapply( rownames( true.edgesList )[FN.type.idx]
                        , FUN=function(myKey)
                        {
                          ifelse( isTRUE(all.equal(outputSummary.df[myKey, c("x", "y")], true.edgesList[myKey, c("x", "y")] )), 2, -2 )
                        } )
      if( length( names( FN.ort ) ) > 0)
      {
        outputSummary.df[names(FN.ort), "trueOrt"] = FN.ort
      }
    }

    if( length(TP.type.idx) > 0 )
    {
      #### For the TP (forward = 2, backward = -2)
      TP.ort <- sapply( rownames( outputSummary.df )[TP.type.idx]
                        , FUN=function(myKey)
                        {
                          # print("\n")
                          # print(outputSummary.df[myKey, c("x", "y")])
                          # print(true.edgesList[myKey, c("x", "y")] )
                          # print(all.equal(outputSummary.df[myKey, c("x", "y")], true.edgesList[myKey, c("x", "y")] ))
                          #ifelse( isTRUE(all.equal(outputSummary.df[myKey, c("x", "y")], true.edgesList[myKey, c("x", "y")] )), 2, -2 )
                          ifelse((outputSummary.df[myKey, "x"] == true.edgesList[myKey,"x"] & outputSummary.df[myKey, "y"] == true.edgesList[myKey,"y"] ), 2, -2)
                        } )
      if( length( names( TP.ort ) ) > 0)
      {
        outputSummary.df[names(TP.ort), "trueOrt"] = TP.ort
      }
    }
    #### ----

    #### Get the information about the correct orientation
    #### ----
    # -- oriented TP edge
    ortEdge.idx = which( ( abs( outputSummary.df[, "infOrt"] ) > 1 ) & ( outputSummary.df[, "type"] == 'TP' ) )
    if( length( ortEdge.idx ) > 0 ) { outputSummary.df[ortEdge.idx, "isOrt"] = 'Y' }

    # -- compatible orientation
    compatibleOrt.idx = which( ( outputSummary.df[, "type"] == 'TP') &
                                 ( ( sign( outputSummary.df[, "infOrt"] ) == sign( outputSummary.df[, "trueOrt"] ) ) |
                                     ( outputSummary.df[, "infOrt"] %in% c(1,6) ) ) )
    if( length( compatibleOrt.idx ) > 0 ) { outputSummary.df[compatibleOrt.idx, "isOrtOk"] = 'Y' }

    # -- not compatible orientation
    notCompatibleOrt.idx = which( ( outputSummary.df[, "type"] == 'TP' ) & ( is.na( outputSummary.df[, "isOrtOk"] ) ) )
    if( length( notCompatibleOrt.idx ) > 0 ) { outputSummary.df[notCompatibleOrt.idx, "isOrtOk"] = 'N' }
    #### ----

    #### Add the essentiality info if any
    if( "essential" %in% colnames(true.edgesList) )
    { outputSummary.df[rownames(true.edgesList), "essential"] <- true.edgesList[,"essential"] }

  } else {

    #### First, initialise all inferred edge as N
    outputSummary.df[, "type"] <- 'N'

    #### Then, get the index of oriented edges and set them to P
    remainingEdges.idx <- which( abs(outputSummary.df[, "infOrt"]) > 0 )
    if( length(remainingEdges.idx) > 0 )
    { outputSummary.df[rownames(outputSummary.df[remainingEdges.idx,]), "type"] <- 'P' }
  }

  #### Use the mutual information values
  #### ----
  #### Add the mutual information values, the cplx values and the difference (ie, confidence) if the cplx exisits...
  if( "cplx" %in% colnames( allMutInfo.df ) )
  {
    outputSummary.df[rownames(outputSummary.df), c( "info", "cplx")] <- allMutInfo.df[rownames(outputSummary.df), c("Ixy_ai", "cplx")]
    outputSummary.df[, "log_confidence"] = as.numeric(as.numeric(outputSummary.df[, "info"]) - as.numeric(outputSummary.df[, "cplx"]) )

  } else {
    outputSummary.df[rownames(outputSummary.df), c( "info")] <- allMutInfo.df[rownames(outputSummary.df), c("Ixy_ai")]
  }

  #### Add the Nxy_ui if it exists
  if( "Nxy_ai" %in% colnames( allMutInfo.df ) )
  { outputSummary.df[rownames(outputSummary.df), "Nxy_ai"] <- allMutInfo.df[rownames(outputSummary.df), "Nxy_ai"] }

  #### Add the {ai}
  if( "ai.vect" %in% colnames( allMutInfo.df ) )
  { outputSummary.df[rownames(outputSummary.df), "ai"] <- allMutInfo.df[rownames(outputSummary.df), "ai.vect"] }

  #### Order by decreasing value of confidence
  outputSummary.df <- outputSummary.df[order(as.numeric(outputSummary.df[,"log_confidence"]), decreasing = TRUE),]


  ### AJOUT LV ###
  ### Compute the sign of each edge and fill the two last columns (sign and partial correlation)
  if( verbose == TRUE )
    cat("\t# -> Computing the sign of the edges ---- \n")
  inputData.df  <- inputData
  #### Remove the lines that are all 'NA'
  allNAs.idx = which( rowSums( is.na( inputData.df ) ) == ncol( inputData.df ) )
  if( length( allNAs.idx ) > 0 ) { inputData.df = inputData.df[-allNAs.idx,] }

  gV$data = inputData.df
  rm(inputData.df)
  stateOrderRet <- dataToStateOrder( gV, stateOrder, inputData )

  if( verbose == TRUE ){
    cat("\tSigns are calculated with the PARTIAL correlation method\n")
  }

  if(correlationToEvaluate)
    outputSummary.df[,c("sign","partial_correlation")] <- computeSign.pcor(outputSummary.df, gV)

  res = list()

  if(!is.null(trueEdges)){

    #### Compute complementary stats
    #### Open the final summary and count the number of TP, FP, TN, FN, TPnort, precision, recall, fscore
    retListSkeleton = list()
    retListSkeleton = computeSingleSampleStat( myFinalSummary = outputSummary.df )

    #### SHD2 for ORIENTATIONS
    infMethod = "miic"
    retListOrient = list()


    retListOrient = computeSHD2SampleStat( myTrueGraph = trueEdges, myFinalSummary = outputSummary.df, skeletonStats = retListSkeleton, myInfMethod = infMethod, myAllProperties = gV$allProperties )

    outputSummary.df = retListOrient$summary
    res$statistics = retListOrient$stats
  }

  spentTime.whole <- (proc.time() - startTime.whole)


  outputSummary.df <- outputSummary.df[,!colnames(outputSummary.df) %in% c("essential")]
  rownames(outputSummary.df) <- NULL
  res$all.edges.summary <- outputSummary.df


  outputSummary.df <-outputSummary.df[which(outputSummary.df$type %in% c("P","TP","FP")),]

  res$retained.edges.summary <- outputSummary.df
  res
}

