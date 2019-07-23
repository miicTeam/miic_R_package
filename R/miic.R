#' MIIC, causal network learning algorithm including latent variables
#'
#' @description MIIC (Multivariate Information based Inductive Causation) combines
#' constraint-based and information-theoretic approaches to disentangle direct
#' from indirect effects amongst correlated variables, including cause-effect
#' relationships and the effect of unobserved latent causes.
#'
#' @details Starting from a complete graph, the method iteratively removes
#' dispensable edges, by uncovering significant information contributions from
#' indirect paths, and assesses edge-specific confidences from randomization of
#' available data. The remaining edges are then oriented based on the signature
#' of causality in observational data.
#'
#' @references
#' \itemize{
#' \item Verny et al., \emph{PLoS Comp. Bio. 2017.}
#' }
#'
#' @param inputData [a data frame]
#' A data frame that contains the observational data. Each
#' column corresponds to one variable and each row is a sample that gives the
#' values for all the observed variables. The column names correspond to the
#' names of the observed variables. Numeric columns will be treated as continuous 
#' values, factors and character as categorical.
#'
#' @param blackBox [a data frame]
#' An optional data frame containing the
#' pairs of variables that should be considered as independent. Each row contains
#' one column for each of the two variables. The variable name must correspond to
#' the one in the \emph{inputData} data frame.
#'
#' @param neff [a positive integer]
#' The N samples given in the \emph{inputdata} data frame are expected
#' to be independent. In case of correlated samples such as in time series or
#' Monte Carlo sampling approaches, the effective number of independent samples
#' \emph{neff} can be estimated using the decay of the autocorrelation function
#' (Verny \emph{et al.}, PLoS Comp. Bio. 2017). This \emph{effective} number \emph{neff}
#' of \emph{independent} samples can be provided using this parameter.
#'
#' @param cplx [a string; \emph{c("nml", "mdl")}]
#' In practice, the finite size of the input
#' dataset requires that the 2-point and 3-point information measures should be
#' \emph{shifted} by a \emph{complexity} term. The finite size corrections can be
#' based on the Minimal Description Length (MDL) criterion (set the option with "mdl").
#' In practice, the MDL complexity criterion tends to underestimate the relevance of
#' edges connecting variables with many different categories, leading to the removal of
#' false negative edges. To avoid such biases with finite datasets, the (universal)
#' Normalized Maximum Likelihood (NML) criterion can be used (set the option with "nml").
#' The default is "nml" (see Affeldt \emph{et al.}, UAI 2015).
#'
#' @param latent [a boolean value]
#' When set to TRUE, the network reconstruction is taking into account
#' hidden (latent) variables. Dependence between two observed variables due to
#' a latent variable is indicated with a '6' in the adjacency matrix and in the
#' network edges.summary and by a bi-directed edge in the (partially) oriented graph.
#'
#' @param orientation [a boolean value]
#' The miic network skeleton can be partially directed
#' by orienting and propagating edge directions, based on the sign and magnitude
#' of the conditional 3-point information of unshielded triples. The propagation
#' procedure relyes on probabilities; for more details, see Verny \emph{et al.}, PLoS Comp. Bio. 2017).
#' If set to FALSE the orientation step is not performed.
#'
#' @param propagation [a boolean value]
#' If set to FALSE, the skeleton is partially oriented with only the
#' v-structure orientations. Otherwise, the v-structure orientations are
#' propagated to downstream undirected edges in unshielded triples following
#' the orientation method
#'
#' @param categoryOrder [a data frame] An optional data frame giving information
#' about how to order the various states of categorical variables. It will be
#' used to compute the signs of the edges (using partial correlation coefficient)
#' by sorting each variable’s levels accordingly to the given category order.
#'
#' @param trueEdges [a data frame]  An optional data frame containing all the
#' true edges of the graph. Each line corresponds to one edge.
#'
#' @param edges [a data frame] The miic$edges object returned by an execution
#' of the miic function. It represents the result of the skeleton step. If this
#' object is provided, the skelethon step will not be done, and the required
#' orientation will be performed using this edges data frame.
#'
#' @param confidenceShuffle [a positive integer] The number of shufflings of
#' the original dataset in order to evaluate the edge specific confidence
#' ratio of all inferred edges.
#'
#' @param confidenceThreshold [a positive floating point] The threshold used
#' to filter the less probable edges following the skeleton step. See Verny
#' \emph{et al.}, PLoS Comp. Bio. 2017.
#'
#' @param confList [a data frame] An optional data frame containing the
#' confFile data frame returned by a miic execution. It is useful when a
#' second run of the same input data set has to be performed with a different
#' confidence threshold and the same confidenceShuffle value. In
#' this way the computations based on the randomized dataset do not need to
#' be performed again, and the values in this data frame are used instead.
#'
#' @param sampleWeights [a numeric vector]
#' An optional vector containing the weight of each observation.
#'
#' @param testMAR [a boolean value]
#' If set to TRUE, distributions with missing values will be tested with Kullback-Leibler 
#' divergence : conditioning variables for the given link \eqn{X\arrow Y}\eqn{Z} will be 
#' considered only if the divergence between the full distribution and the non-missing 
#' distribution \eqn{KL(P(X,Y) | P(X,Y)_{!NA})} is low enough (with \eqn{P(X,Y)_{!NA}} as 
#' the joint distribution of \eqn{X} and \eqn{Y} on samples which are not missing on Z. 
#' This is a way to ensure that data are missing at random for the considered
#' interaction and to avoid selection bias. Set to TRUE by default
#'
#' @param consistent [a boolean value] If TRUE, runs the consistent version.
#'
#' @param verbose [a boolean value] If TRUE, debugging output is printed.
#'
#' @param nThreads [a positive integer]
#' When set greater than 1, nThreads parallel threads will be used for computation. Make sure 
#' your compiler is compatible with openmp if you wish to use multithreading. 
#'
#' @param doConsensus [a positive integer] If doConsensus is larger than 0 (default
#' value), it will create bootstraping skeletons for the estimation of a consensus skeleton
#' made by the associations that were common among the bootstraping skeletons. The positive
#' integer indicates the % of the bootstraping skeletons that the edge must exist in order to
#' be in the consensus skeleton. A value of zero indicates that no bootstraping should
#' be performed.
#'
#' @param nSkeletons [a positive integer] The number of bootstraping skeletons that will
#' be inferred for the estimation of a consensus skeleton. The last skeleton will be used for
#' the next steps of the MIIC algorithm.
#' 
#' @return A \emph{miic-like} object that contains:
#' \itemize{
#'  \item{all.edges.summary:}{ a data frame with information about the relationship between
#'  each pair of variables
#'  \itemize{
  #'  \item \emph{x:} X node
  #'  \item \emph{y:} Y node
  #'  \item \emph{type:} contains 'N' if the edge has
  #'  been removed or 'P' for retained edges. If a true edges file is given,
  #'  'P' becomes 'TP' (True Positive) or 'FP' (False Positive), while
  #'  'N' becomes 'TN' (True Negative) or 'FN' (False Negative).
  #'  \item \emph{ai:} the contributing nodes found by the method which participate in
  #'  the mutual information between \emph{x} and \emph{y}, and possibly separate them.
  #'  \item \emph{info:} provides the final mutual information times \emph{Nxy_ai} for
  #'  the pair (\emph{x}, \emph{y}) when conditioned on the collected nodes \emph{ai}.
  #'  \item \emph{cplx:} gives the computed complexity between the (\emph{x}, \emph{y})
  #'  variables taking into account the contributing nodes \emph{ai}.
  #'  \item \emph{Nxy_ai:} gives the number of samples on which the information and the
  #'  complexity have been computed. If the input dataset has no missing value, the
  #'  number of samples is the same for all pairs and corresponds to the total
  #'  number of samples.
  #'  \item \emph{log_confidence:} represents the \emph{info} - \emph{cplx} value.
  #'  It is a way to quantify the strength of the edge (\emph{x}, \emph{y}).
  #'  \item \emph{confidenceRatio:} this column is present if the confidence cut
  #'  is > 0 and it represents the ratio between the probability to reject
  #'  the edge (\emph{x}, \emph{y}) in the dataset versus the mean probability
  #'  to do the same in multiple (user defined) number of randomized datasets.
  #'  \item \emph{infOrt:} the orientation of the edge (\emph{x}, \emph{y}). It is
  #'  the same value as in the adjacency matrix at row \emph{x} and column \emph{y}.
  #'  \item \emph{trueOrt:} the orientation of the edge (\emph{x}, \emph{y}) present
  #'  in the true edges file (if true edges file is provided).
  #'  \item \emph{isOrtOk:} information about the consistency of the inferred graph’s
  #'  orientations with a reference graph is given (i.e. if true edges file is provided).
  #'  Y: the orientation is consistent; N: the orientation is not consistent with
  #'  the PAG derived from the given true graph.
  #'  \item \emph{sign:} the sign of the partial correlation between variables
  #'  \emph{x} and \emph{y}, conditioned on the contributing nodes \emph{ai}.
  #'  \item \emph{partial_correlation:} value of the partial correlation for the
  #'  edge (\emph{x}, \emph{y}) conditioned on the contributing nodes \emph{ai}.
  #'  }
#'  }
#'
#'  \item{retained.edges.summary:} {a data frame in the format of all.edges.summary containing only the inferred edges.}
#'
#'  \item{orientations.prob:} {this data frame lists the orientation probabilities of the two edges of all unshielded triples
#'  of the reconstructed network with the structure: node1 -- mid-node -- node2:}
#'  \itemize{
#'  \item node1: node at the end of the unshielded triplet
#'  \item p1: probability of the arrowhead node1 <- mid-node
#'  \item p2: probability of the arrowhead node1 -> mid-node
#'  \item mid-node: node at the center of the unshielded triplet
#'  \item p3: probability of the arrowhead mid-node <- node2
#'  \item p4: probability of the arrowhead mid-node -> node2
#'  \item node2: node at the end of the unshielded triplet
#'  \item NI3: 3 point (conditional) mutual information * N
#'  }
#'
#'
#'  \item {AdjMatrix:} the adjacency matrix is a square matrix used to represent
#'  the inferred graph. The entries of the matrix indicate whether pairs of
#'  vertices are adjacent or not in the graph. The matrix can be read as a
#'  (row, column) set of couples where the row represents the source node and
#'  the column the target node. Since miic can reconstruct mixed networks
#'  (including directed, undirected and bidirected edges), we will have a
#'  different digit for each case:
  #'  \itemize{
    #'  \item 1: (\emph{x}, \emph{y}) edge is undirected
    #'  \item 2: (\emph{x}, \emph{y}) edge is directed as \emph{x} -> \emph{y}
    #'  \item -2: (\emph{x}, \emph{y}) edge is directed as \emph{x} <- \emph{y}
    #'  \item 6: (\emph{x}, \emph{y}) edge is bidirected
  #'  }
  #'}
#' @export
#' @useDynLib miic
#'
#' @examples
#' library(miic)
#'
#' # EXAMPLE HEMATOPOIESIS
#' data(hematoData)
#'
#' # execute MIIC (reconstruct graph)
#' miic.res = miic(inputData = hematoData, latent = TRUE,
#' confidenceShuffle = 10, confidenceThreshold = 0.001)
#'
#' # plot graph
#' miic.plot(miic.res)
#'\dontrun{
#'
#' # write graph to graphml format. Note that to correctly visualize
#' # the network we created the miic style for Cytoscape (http://www.cytoscape.org/).
#'
#' miic.write.network.cytoscape(g = miic.res, file = file.path(tempdir(),"temp"))
#'
#' # EXAMPLE CANCER
#' data(cosmicCancer)
#' data(cosmicCancer_stateOrder)
#' # execute MIIC (reconstruct graph)
#' miic.res = miic(inputData = cosmicCancer, categoryOrder = cosmicCancer_stateOrder, latent = TRUE,
#' confidenceShuffle = 100, confidenceThreshold = 0.001)
#'
#' # plot graph
#' miic.plot(miic.res, igraphLayout=igraph::layout_on_grid)
#'
#' # write graph to graphml format. Note that to correctly visualize
#' # the network we created the miic style for Cytoscape (http://www.cytoscape.org/).
#' miic.write.network.cytoscape(g = miic.res, file = file.path(tempdir(),"temp"))
#'
#' # EXAMPLE OHNOLOGS
#' data(ohno)
#' data(ohno_stateOrder)
#' # execute MIIC (reconstruct graph)
#' miic.res = miic(inputData = ohno, latent = TRUE, categoryOrder = ohno_stateOrder,
#' confidenceShuffle = 100, confidenceThreshold = 0.001)
#'
#' # plot graph
#' miic.plot(miic.res)
#'
#' # write graph to graphml format. Note that to correctly visualize
#' # the network we created the miic style for Cytoscape (http://www.cytoscape.org/).
#' miic.write.network.cytoscape(g = miic.res, file = file.path(tempdir(),"temp"))
#'}

miic <- function(inputData, categoryOrder= NULL, trueEdges = NULL, blackBox = NULL, nThreads=1, 
                 cplx = c("nml", "mdl"), orientation = TRUE, propagation = TRUE, latent = FALSE,
                 neff = -1, edges=NULL, confidenceShuffle = 0, confidenceThreshold = 0, 
                 confList = NULL, sampleWeights = NULL, testMAR = TRUE, consistent = FALSE, 
                 verbose = FALSE, doConsensus=0, nSkeletons=0)
{
  res = NULL
  skeleton = TRUE

  #### Check the input arguments
  if( is.null( inputData ) )
  { stop("The input data file is required") }

  if( !is.data.frame(inputData))
  { stop("The input data is not a dataframe") }

  effnAnalysis = miic.evaluate.effn(inputData, plot=F)
  if(effnAnalysis$neff < 0.5 * nrow(inputData) ){
    if(effnAnalysis$exponential_decay){
        print(paste("Warning! Your samples in the datasets seem to be correlated! We suggest to re run the method specifying",
          effnAnalysis$neff, " in the neff parameter. See the autocorrelation plot for more details."))
      } else {
        print("Warning! Your samples in the datasets seem to be correlated but the correlation decay is not exponential. Are your samples correlated in some way? See the autocorrelation plot for more details.")
      }
  }

  # if(!is.null(trueEdges)){
  #   if(length(which(!c(as.vector(trueEdges[,1]), as.vector(trueEdges[,2])) %in% colnames(inputData))) > 0){
  #     stop("True edges file does not correspond to the input data matrix. Please check node names.")
  #   }
  # }

  #cplx
  if(length(cplx) == 2)
    cplx = "nml"

  if(cplx != "nml" && cplx != "mdl")
  { stop("The complexity method is not allowed") }

  if(neff > nrow(inputData))
  { stop("The number of effective samples can't be greater than the number of samples.") }

  #edges
  if(!is.null(edges)){
    skeleton = FALSE
    if(length(colnames(edges)) != 10)
      stop("The edges data frame is not correct. The required data frame is the $edges output of the miic method")
  }

  #propagation
  if(propagation != TRUE && propagation != FALSE)
    stop("The propagation type is not correct. Allowed types are TRUE or FALSE")

  #orientation
  if(orientation != TRUE && orientation != FALSE)
    stop("The orientation type is not correct. Allowed types are TRUE or FALSE")

  #latent
  if(latent != TRUE && latent != FALSE)
    stop("The latent type is not correct. Allowed types are TRUE or FALSE")
  
  #Bootstraping
  if(doConsensus == 0 && nSkeletons != 0)
    stop(paste0("It's useless to set the nSkeletons parameter to different ",
                "than zero if you did not set the doConsensus parameter."))
  if(doConsensus != 0 && nSkeletons == 0)
    stop(paste0("You must set nSkeletons if you want to build a consensus",
                " skeleton."))
  if(doConsensus > 100 || doConsensus < 0)
    stop("doConsensus can not be smaller than 0 or greater than 100.")
  if(nSkeletons < 0)
    stop("nSkeletons can not be smaller than 0.")
  
  if(verbose)
    cat("START miic...\n")

  # continuous or discrete?
  cntVar = sapply(inputData, is.numeric)
  for(col in colnames(inputData)){
    if(cntVar[[col]] && (length(unique(inputData[[col]])) <= 40) && (nrow(inputData) > 40 )){
      warning(paste0("Variable ", col, " is treated as continuous but only has ", length(unique(inputData[[col]])), " unique values."))
    }
    if((!cntVar[[col]]) && (length(unique(inputData[[col]])) >= 40) && (nrow(inputData) > 40 )){
      warning(paste0("Variable ", col, " is treated as discrete but has many levels (", length(unique(inputData[[col]])), ")."))
    }
  }
  typeOfData = 0 # Assume all discrete
  if(any(cntVar)){ 
    typeOfData = 2 # Mixed if any are continuous
    if(all(cntVar)){
      typeOfData = 1 # All continuous
    }
  }

  err_code = checkInput(inputData, "miic")
  if(err_code != "0"){
    print(errorCodeToString(err_code))
  } else {

    if(skeleton){
      if(verbose)
        cat("\t# -> START skeleton...\n")
      if (doConsensus == 0) {
        res <- miic.skeleton(inputData = inputData, stateOrder= categoryOrder, nThreads= nThreads, cplx = cplx, latent = latent,
                             effN = neff, blackBox = blackBox, confidenceShuffle = confidenceShuffle,
                             confidenceThreshold= confidenceThreshold, verbose= verbose, cntVar = cntVar, typeOfData = typeOfData,
                             sampleWeights = sampleWeights, testMAR = testMAR, consistent = consistent)
      } else {
        skeletons <- NULL
        # All inferred skeletons will be used to compute the consensus skeleton. The last
        # skeleton inferred will also be used for the next steps of the algorithm.
        for (i in seq(nSkeletons)) {
          # Bootstrap: sample with replacement
          input <- inputData[sample(nrow(inputData), replace=TRUE), ]
          # If it's the last skeleton, do it from the not-bootstrapped sample
          if (i == nSkeletons) {
            res <- miic.skeleton(inputData = inputData, stateOrder= categoryOrder, nThreads= nThreads, cplx = cplx, latent = latent,
                                 effN = neff, blackBox = blackBox, confidenceShuffle = confidenceShuffle,
                                 confidenceThreshold= confidenceThreshold, verbose= verbose, cntVar = cntVar, typeOfData = typeOfData,
                                 sampleWeights = sampleWeights, testMAR = testMAR, consistent = consistent)
          } else {
            # Infer skeleton from a bootstraping sample
            # Try to make miic.sckeleton smarter for this case, like not saving anythinf for orientation
            res <- miic.skeleton(inputData = input, stateOrder= categoryOrder, nThreads= nThreads, cplx = cplx, latent = latent,
                                 effN = neff, blackBox = blackBox, confidenceShuffle = confidenceShuffle,
                                 confidenceThreshold= confidenceThreshold, verbose= verbose, cntVar = cntVar, typeOfData = typeOfData,
                                 sampleWeights = sampleWeights, testMAR = testMAR, consistent = consistent)
          }
          # These lines and the global variables that change their scope are temporary and
          # will be removed at the end of the development of this feature
          # (consensus skeleton).
          x <- res$edges[res$edges$category == 3, 'x']
          y <- res$edges[res$edges$category == 3, 'y']
          I <- res$edges[res$edges$category == 3, 'Ixy_ai']
          ai_vect <- res$edges[res$edges$category == 3, 'ai.vect']
          ai_vect_n <- sapply(ai_vect, function (x) ifelse(is.na(x),
                                                           0,
                                                           length(unlist(strsplit(x, ',')))))
          
          skeletons <- rbind(skeletons, cbind(x, y, I, ai_vect_n))
          rownames(skeletons) <- NULL
        }
        skeletons <<- skeletons
        res <<- res
        # Computing consensus
        
      }
      if(res$interrupted){
        warning("Interupted by user")
        return(NULL)
      }
      # print(res)
      edges = res$edges
      confData = res$confData
      time = res$time
      if(verbose)
        cat("\t# -> END skeleton...\n\t# --------\n")
    }

    if( confidenceShuffle < 0 | confidenceThreshold < 0 ){
      cat("Warning! ConfidenceShuffle and confidenceThreshold must be greater than 0, the confidence cut step will not be performed.")
      confidenceShuffle =0
    }


    timeOrt=0
    if(orientation){

      if(verbose)
        cat("\tSTART orientation...")
      ptm <- proc.time()
      res = miic.orient(inputData= inputData, stateOrder = categoryOrder, edges = edges, effN = neff,
                        cplx = cplx,  latent = latent, propagation = propagation, cntVar = cntVar, 
                        typeOfData = typeOfData, verbose = FALSE)

      timeOrt=(proc.time() - ptm)["elapsed"]
      timeInitIterOrt = time["initIter"]+timeOrt

      res$edges <- edges
      res$confData <- confData
      if(verbose)
        cat("\tEND orientation...")
    }

    # Summarize the results
    # --------

    if( !is.null(trueEdges ) ){
      err_code = checkTrueEdges(trueEdges)
      if(err_code != "0"){
        print(errorCodeToString(err_code))
        print("WARNING: True edges file will be ignored!")
        trueEdges = NULL
      }
    }


    #STATE ORDER FILE
    if(  !is.null( categoryOrder ) ){
      err_code = checkStateOrder(categoryOrder, inputData)
      if(err_code != "0"){
        print(errorCodeToString(err_code))
        print("WARNING: Cathegory order file will be ignored!")
        categoryOrder = NULL
      }
    }

    # Call the function

    ptm <- proc.time()
    resGmSummary <- gmSummary(inputData = inputData, edges = res$edges,
                              adjMatrix= res$adjMatrix, trueEdges = trueEdges, stateOrder = categoryOrder, verbose = verbose)

    timeInitIterOrt = timeOrt + time[4]
    timeSum=(proc.time() - ptm)["elapsed"]
    timeTotal = timeInitIterOrt+timeSum
    timeVec = c(time, timeOrt, timeInitIterOrt, timeSum, timeTotal)
    timeVec[which(timeVec==0)]=NA
    res$time = stats::setNames(timeVec,
                c("initialization", "iteration", "confidenceCut", "skeleton", "orientation", "skeleton+Orientation", "summary", "total"))

    res$all.edges.summary <- resGmSummary$all.edges.summary
    res$retained.edges.summary <- resGmSummary$retained.edges.summary
    res$statistics = resGmSummary$statistics

    rm(resGmSummary)

    if( confidenceShuffle > 0 & confidenceThreshold > 0 )
    {
      # Insert the confidence ratio
      tmp_sum = res$all.edges.summary
      conf_col = rep(1, nrow( tmp_sum ) )
      isCut = rep(NA, nrow(tmp_sum))
      tmp_sum = cbind(tmp_sum,conf_col,isCut)

      tmp_sum = cbind(tmp_sum,conf_col)

      tmp_pval = res$confData
      for(r in 1:nrow(tmp_pval))
      {
        tmp_sum[which(tmp_sum[,"x"] == tmp_pval[r,"x"] & tmp_sum[,"y"] == tmp_pval[r,"y"]),'confidence_ratio'] =tmp_pval[r,"confidence_ratio"]
        if(tmp_pval[r,"confidence_ratio"] < confidenceThreshold){
          tmp_sum[which(tmp_sum[,"x"]==tmp_pval[r,"x"] & tmp_sum[,"y"]==tmp_pval[r,"y"]),'isCut']='N'
        }
        else{
          tmp_sum[which(tmp_sum[,"x"]==tmp_pval[r,"x"] & tmp_sum[,"y"]==tmp_pval[r,"y"]),'isCut']='Y'
        }
      }
      tmp_sum = tmp_sum[,c('x','y','type','ai','info','cplx','Nxy_ai','log_confidence','confidence_ratio','infOrt','trueOrt', 'isOrt', 'isOrtOk', 'sign','partial_correlation', 'isCut')]

      res$all.edges.summary= tmp_sum

      res$retained.edges.summary= tmp_sum[which(tmp_sum$type %in% c("P","TP","FP")),]
    }
    if(verbose)
      cat("END miic")
  }
  res
}
