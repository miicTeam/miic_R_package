#' Discretize two distributions with
#' @description This function discretizes two distributions by maximizing their mutual information, with MDL.
#'
#' @param myDist1 [a vector]
#' A vector that contains the observational data of the first variable.
#' @param myDist2 [a vector]
#' A vector that contains the observational data of the second variable.
#' @param maxbins [an int] The maximum number of bins to test for.
#' @return A list with the two vectors containing the cutpoints of the best discretization for both variables.
#' @export
#' @useDynLib miic

discretizeMutual <- function(myDist1 = NULL, myDist2 = NULL, maxbins=50, plot=T)
{
  result = list()
  #### Check the input arguments
  if( is.null( myDist1 ) || is.null(myDist2) )
  { stop("The input data file is required") }

  if(maxbins > length(myDist1))
    maxbins=length(myDist1)

  myDist1[is.na(myDist1)] = -1
  myDist2[is.na(myDist2)] = -1

  if (base::requireNamespace("Rcpp", quietly = TRUE)) {
    rescpp <- .Call('mydiscretizeMutual', myDist1, myDist2, maxbins, PACKAGE = "miic")
  }
  niterations = length(rescpp$cutpoints1)/maxbins

  result$niterations = niterations
  for(i in 0:(niterations-1)){
    clean_cutpoints1 = rescpp$cutpoints1[(maxbins*i)+(1:maxbins)]
    clean_cutpoints1 = clean_cutpoints1[clean_cutpoints1 != -1]
    clean_cutpoints1 = unique(sort(myDist1)[c(1, clean_cutpoints1+1, length(myDist1))])
    clean_cutpoints2 = rescpp$cutpoints2[(maxbins*i)+(1:maxbins)]
    clean_cutpoints2 = clean_cutpoints2[clean_cutpoints2 != -1]
    clean_cutpoints2 = unique(sort(myDist2)[c(1, clean_cutpoints2+1, length(myDist2))])
    result[[paste0("iteration",i+1)]] = list(cutpoints1=clean_cutpoints1, cutpoints2=clean_cutpoints2)
  }
  result$cutpoints1 = result[[paste0("iteration", niterations)]]$cutpoints1
  result$cutpoints2 = result[[paste0("iteration", niterations)]]$cutpoints2

  library(infotheo)
  result$info = infotheo::mutinformation(cut(myDist1, result$cutpoints1), cut(myDist2, result$cutpoints2))
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

axisprint <- function(x) sprintf("%.2f", x)

theme_side_hist <- function () {
  theme_classic() %+replace%
    theme(title = element_text(family = "", face = "plain",
                              color = NA, size = theme_classic()$text$size,
                              hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                              margin = margin(), debug = FALSE),
          text = element_text(family = "", face = "plain",
                              color = NA, size = theme_classic()$text$size,
                              hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                              margin = margin(), debug = FALSE),
          axis.text = element_text(family = "", face = "plain",
                              color = NA, size = theme_classic()$text$size,
                              hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                              margin = margin(), debug = FALSE),
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

  hist1 = ggplot(data.frame(myDist1), aes(x=myDist1)) +
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   breaks = cut_points1,
                   colour="black", fill="white") +
    geom_density(adjust=0.5, alpha=.5, fill="#c1c6ee") +  # Overlay with transparent density plot
    theme_side_hist() %+replace% theme(plot.margin = margin(5.5,5.5,-25,5.5,"pt")) +
    scale_y_continuous(labels=axisprint)

  hist2 = ggplot(data.frame(myDist2), aes(x=myDist2)) +
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   breaks = cut_points2,
                   colour="black", fill="white") +
    geom_density(adjust=0.5, alpha=.5, fill="#c1c6ee") +  # Overlay with transparent density plot
    theme_side_hist() %+replace% theme(plot.margin = margin(5.5,5.5,5.5,-25,"pt")) +
    coord_flip()

  hist2d = ggplot(data.frame(myDist1, myDist2), aes(x=myDist1, y=myDist2)) +
    stat_bin2d(aes(fill = ..density..), breaks=list(x=cut_points1, y=cut_points2), show.legend = F) +
    scale_fill_gradient(low = "#e1e3f2", high = "#0013a3") +
    geom_vline(xintercept=cut_points1, linetype="dashed", color="grey") +
    geom_hline(yintercept=cut_points2, linetype="dashed", color="grey") +
    geom_point(shape=21, alpha=.7, fill="#ffef77", size=2) +
    theme_classic() +
    scale_y_continuous(labels=axisprint) +
    scale_x_continuous(labels=axisprint)

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


  return(gridExtra::grid.arrange(hist1, empty, hist2d, hist2, ncol=2, nrow=2,
                                 widths=c(4, 1), heights=c(1, 4), bottom=title))
}
