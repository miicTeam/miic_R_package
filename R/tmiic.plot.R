#*****************************************************************************
# Filename   : tmiic.plot.R                     Creation date: 24 march 2020
#
# Description: Function to plot lagged graphs of temporal miic (tmiic)
#
# Author     : Franck SIMON (fsimon.informaticien@wanadoo.fr)
#
# Changes history:
# - 24 march 2020 : initial version
#*****************************************************************************

#-----------------------------------------------------------------------------
# tmiic.plot
#-----------------------------------------------------------------------------
#' tmiic.plot
#'
#' Igraph plotting functions for tmiic
#' 
#' @description This function plots the temporal network with the given layout (if specified) 
#' using the igraph package.
#'
#' @details The plot reports the partial correlation or the log_confidence as strength of the edges.
#'
#' @param g [a miic graph object]
#' The graph object returned by the miic execution.
#' @param method [a string; \emph{c("pcor", "log_confidence")}]
#' The column used to plot the strength of the edges. Default: pcor.
#' @param userLayout [a data frame]
#' An optional data frame reporting the position of nodes. Each line corresponds to the \emph{(x,y)}
#' coordinates of each vertex. This data frame must have three columns, the first containing the
#' name of the vertex as indicated in the colnames of the input data frame, the two others reporting the x and y positions.
#' @param igraphLayout [an igraph layout object]
#' When set it is used to plot the network. See the igraph manual for more information.
#'  Default: \emph{layout_with_kk}
#' @param title [a string] Optional, NULL by default. The title of the plot.
#' @param filename [a string] Optional, NULL by default. If supplied, plot will be saved here.
#' @param curve_edges [a boolean value] Optional, TRUE by default. If TRUE, edges are curved. 
#' @param verbose [a boolean value] If TRUE, debugging output is printed.
#--------------------------------------------------------------------------------
tmiic.plot <- function (g, method="log_confidence", igraphLayout=NULL, userLayout=NULL,
                        title=NULL, filename=NULL, curve_edges=TRUE, verbose=FALSE) 
  {
  DEBUG <- FALSE
  
  if (is.null(g$adjMatrix)) 
    {
    message ("The learnt graphical model adjacency matrix does not exist")
    return()
    }
  if (DEBUG)
    {
    print ("tmiic.plot:")
    print (paste ("input title=", title, sep="") )
    print (paste ("input filename=", filename, sep="") )
    print (paste ("curve_edges=", curve_edges, sep="") )
    print ("input adjacency matrix:")
    print (g$adjMatrix)
    }
  #
  # Set a layout if none is supplied
  #
  if ( is.null(userLayout) ) 
    {
    list_nodes <- colnames(g$adjMatrix)
    #
    # Get the max lag from nodes names
    #
    max_lag <- 0
    for (one_node in list_nodes)
      {
      pos_lag_x <- str_locate(one_node, "_lag")
      lag <- 0
      if ( !is.na(pos_lag_x[1]) )
        {
        lag <- str_remove(one_node, ".*_lag")
        lag <- strtoi (lag)
        }
      if (lag > max_lag)
        max_lag <- lag
      } 
    max_lag_plus1 <- max_lag + 1
    n_nodes <- length(list_nodes) %/% max_lag_plus1
    userLayout = data.frame (list_nodes, rep(1:max_lag_plus1, each=n_nodes), 
                                         rep(1:n_nodes, times=max_lag_plus1) )    

    if (DEBUG)
      {
      print ("no user layout given:")
      print ("found these nodes=")
      print (list_nodes)
      print (paste ("found max lag=", max_lag, sep="") )
      print (paste ("found n_nodes=", n_nodes, sep="") )
      print ("layout defined=")
      print (userLayout)
      }
    } 
  #
  # Plot the graph
  #
  if (! is.null(filename) )
    {
    png (filename=filename)
    }
  miic.plot (g=g, method=method, igraphLayout=igraphLayout, 
             userLayout=userLayout, curve_edges=curve_edges, 
             title=title, verbose=verbose, call_from_tmiic=TRUE) 
  if (! is.null(filename) )
    {
    dev.off()
    }
  }


