#' Discretize two distributions with
#' @description This function discretizes two distributions by maximizing their mutual information, with MDL.
#'
#' @param myDist1 [a vector]
#' A vector that contains the observational data of the first variable.
#' @param myDist2 [a vector]
#' A vector that contains the observational data of the second variable.
#' @param matrixU [a numeric matrix]
#' The matrix containing the values of the conditioning variables with as many columns as variables.
#' @param maxbins [an int] 
#' The maximum number of bins desired in the discretization.
#' @param initbins [an int]
#' The number of bins of the equal
#' @param cplx [a string]
#' The complexity used in the dynamic programming. Either "mdl" for Minimum description Length or "nml" for
#' Normalized Maximum Likelihood, which is less punitive in the finite sample case and will create more bins than mdl.
#'
#' @return A list with the two vectors containing the cutpoints of the best discretization for both variables.
#' @export
#' @useDynLib miic

discretizeMutual <- function(myDist1 = NULL, myDist2 = NULL, matrixU=NULL, maxbins=NULL, initbins=NULL, cplx="mdl", pxy=1, plot=T)
{
  result = list()
  #### Check the input arguments
  if( !is.vector(myDist1) || !is.vector(myDist2) ) {
    stop("Please provide the two samples myDist1 and myDist2 as vectors.")
  }

  if( length(myDist1) != length(myDist2) ){
    stop(paste("The two samples must have the same number of observation (", length(myDist1), "vs",length(myDist2), ")."))
  }

  if( (!is.null(matrixU) && !is.matrix(matrixU)) || (!is.null(matrixU) && nrow(matrixU) != length(myDist1)) ) {
    stop("matrixU is not a matrix or its number of rows differs from the number of observations.")
  }

  if((maxbins > length(myDist1)) || is.null(maxbins))
    maxbins=length(myDist1)

  if((initbins > length(myDist1)) || is.null(initbins))
    initbins=round(length(myDist1)**(1/3))

  myDist1[is.na(myDist1)] = -1
  myDist2[is.na(myDist2)] = -1


  if(is.null(matrixU)){
    nbrU = 0
    flatU = c(0)
  } else{
    nbrU = dim(matrixU)[2]
    flatU <- as.vector(as.matrix(matrixU))
  }

  if(cplx=="mdl"){
    intcplx = 0
  } else if (cplx=="nml"){
    intcplx = 1
  } else {
    print("cplx parameter not understood, please specify either \'mdl\' or \'nml\'. Running with the default option (mdl).")
    intcplx = 0
  }

  if (base::requireNamespace("Rcpp", quietly = TRUE)) {
    rescpp <- .Call('mydiscretizeMutual', myDist1, myDist2, flatU, nbrU, maxbins, initbins, intcplx, pxy, PACKAGE = "miic")
  }
  niterations = length(rescpp$cutpoints1)/maxbins

  result$niterations = niterations
  for(i in 0:(niterations-1)){
      result[[paste0("iteration",i+1)]] = list()
    for(l in 1:(nbrU+2)){
      clean_cutpoints = rescpp$cutpointsmatrix[,l][(maxbins*i)+(1:maxbins)]
      clean_cutpoints = clean_cutpoints[clean_cutpoints != -1]
      if(l==1) clean_cutpoints = unique(sort(myDist1)[c(1, clean_cutpoints+1, length(myDist1))])
      else if(l==2) clean_cutpoints = unique(sort(myDist2)[c(1, clean_cutpoints+1, length(myDist1))])
      else clean_cutpoints = unique(sort(matrixU[,l-2])[c(1, clean_cutpoints+1, length(myDist1))])
      result[[paste0("iteration",i+1)]][[paste0("cutpoints",l)]] = clean_cutpoints
    }
  }
  for(l in 1:(nbrU+2)){
    result[[paste0("cutpoints",l)]] = result[[paste0("iteration", niterations)]][[paste0("cutpoints",l)]]
  }

  library(infotheo)
  result$info = infotheo::mutinformation(cut(myDist1, result$cutpoints1), cut(myDist2, result$cutpoints2))
  if(nbrU>0){
    result$info = infotheo::condinformation(X = cut(myDist1, result$cutpoints1),
                                            Y = cut(myDist2, result$cutpoints2), 
                                            S = cut(matrixU[,1], result$cutpoints3))
  }
  result$infobits = infotheo::natstobits(result$info)

  if(plot) {
    require(ggplot2)
    require(gridExtra)
    jointplot = jointplot_hist(myDist1, myDist2, result)
    result$plot = jointplot
    jointplot
  }

  result
}

#####
# Plot functions

axisprint <- function(x) sprintf("%6s", x)

theme_side_hist <- function () {
  theme_classic() %+replace%
    theme(title = element_text(family = "", face = "plain",
                              color = NA, size = theme_classic()$text$size),
          axis.text = element_text(color = NA, size = theme_classic()$axis.text$size),
          axis.line.x = element_line(colour = NA),
          axis.line.y = element_line(colour = NA),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white", colour = "white")
          )
}


jointplot_hist <- function(myDist1, myDist2, result, title="Joint histogram"){

  library(ggplot2)
  cut_points1 = result$cutpoints1
  cut_points2 = result$cutpoints2
  info = result$info

  # Custom density matrix for colour filling relative to 2D bin area in the 2d histogram
  bin_count = table(cut(myDist1, cut_points1), cut(myDist2, cut_points2))
  bin_areas =
     (cut_points1[-1] - cut_points1[1:(length(cut_points1)-1)]) %*%
    t(cut_points2[-1] - cut_points2[1:(length(cut_points2)-1)])
  fill_density = bin_count / bin_areas
  fill_density = fill_density / sum(fill_density)
  fill_density_flat = data.frame(xstart=numeric(), xend=numeric(),
                                 ystart=numeric(), yend=numeric(), density=numeric())
  for(j in 1:(ncol(fill_density))) {
    for(i in 1:(nrow(fill_density))) {
      fill_density_flat[(j-1)*nrow(fill_density)+i,] = c(cut_points1[i], cut_points1[i+1],
                                                         cut_points2[j], cut_points2[j+1],
                                                         fill_density[i,j])
    }
  }
  fill_density_flat[fill_density_flat$density==0, "density"] = NA

  hist2d = ggplot(fill_density_flat) + 
    geom_rect(aes(xmin=xstart, xmax=xend, ymin=ystart, ymax=yend, fill=density), na.rm = T, show.legend = F) + 
    scale_fill_gradient(low = "#e1e3f2", high = "#0013a3", position = "left", na.value = "white", limits=c(0,max(fill_density_flat$density))) +
    geom_vline(xintercept=cut_points1, linetype="dashed", color="grey") +
    geom_hline(yintercept=cut_points2, linetype="dashed", color="grey") +
    geom_point(data = data.frame(myDist1, myDist2), aes(x=myDist1, y=myDist2), shape=21, alpha=.7, fill="#ffef77", size=2) +
    theme_classic()

  g = ggplot_build(hist2d)
  labels = g$layout$panel_ranges[[1]]$y.labels


  side_hist_top = ggplot(data.frame(myDist1), aes(x=myDist1)) +
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   breaks = cut_points1,
                   colour="black", fill="white") +
    geom_density(adjust=0.5, alpha=.5, fill="#c1c6ee") +  # Overlay with transparent density plot
    theme_side_hist() %+replace% theme(plot.margin = margin(5.5,5.5,-25,5.5,"pt")) +
    scale_y_continuous(labels=labels, breaks=seq(0,1,length.out = length(labels)), expand=c(0,0)) #Pass hist2d's labels to align both X axes

  side_hist_bot = ggplot(data.frame(myDist2), aes(x=myDist2)) +
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   breaks = cut_points2,
                   colour="black", fill="white") +
    geom_density(adjust=0.5, alpha=.5, fill="#c1c6ee") +  # Overlay with transparent density plot
    theme_side_hist() %+replace% theme(plot.margin = margin(5.5,5.5,5.5,-30,"pt")) +
    scale_y_continuous(expand=c(0,0)) +
    coord_flip()

  I2 = info
  N = length(myDist1)
  rp = length(cut_points1)
  rq = length(cut_points2)
  bic = 1/2*(rp-1)*(rq-1)*log(N)
  cpl = log(choose(N-1, rp-1))+ log(choose(N-1, rq-1))
  I2p = I2 - cpl/N

  empty <- ggplot()+geom_point(aes(1,1), colour="white")+
    geom_text(aes(x=1, y=0.5, label=paste("I(X;Y) =", round(I2,3)))) +
    geom_text(aes(x=1, y=0, label=paste("I(X;Y)-kBIC =", round(I2p,3)))) +
    theme(axis.ticks=element_blank(),
          panel.background=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank())


  return(gridExtra::grid.arrange(side_hist_top, empty, hist2d, side_hist_bot, ncol=2, nrow=2,
                                 widths=c(4.2, 1), heights=c(1, 4.2), bottom=title))
}
