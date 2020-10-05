#*****************************************************************************
# Filename   : tmiic.plot.R                     Creation date: 24 march 2020
#
# Description: Plotting for temporal miic (tmiic)
#
# Author     : Franck SIMON (fsimon.informaticien@wanadoo.fr)
#
# Changes history:
# - 24 march 2020 : initial version
# - 11 may 2020   : modification to handle condensed networks
# - 11 june 2020 : initial version
#   The plot function is a rewrite from miic.plot.R and gmPlot.lib.R. 
#   + add capabiliy to plot multiple edges between the same nodes
# - 11 sept 2020 : rewrite to be aligned with 1.5.0 miic plotting 
#*****************************************************************************

#-----------------------------------------------------------------------------
# tmiic.export
#-----------------------------------------------------------------------------
#' Export temporal miic (tmiic) result to different plotting methods
#'
#' @description This function creates an object built from the result returned
#' by \code{\link{miic}} executed in temporal mode that is ready to be fed to 
#' different plotting methods.
#'
#' @details See the details of specific function for each method.
#' For igraph, see \code{\link{tmiic.getIgraph}}.
#'
#' @param tmiic.res [a tmiic graph object]
#' The graph object returned by the \code{\link{miic}} execution in temporal mode.
#' 
#' @param flatten_mode [a string]. Optional, default value "compact".
#' Possible values are \emph{"none"}, \emph{"compact"}, \emph{"combine"}, 
#' \emph{"unique"}, \emph{"drop"}:
#' \itemize{
#' \item When \emph{flatten_mode} = \emph{"none"}, the export function will
#'   use the tmiic graph object as it, leading to the return of a lagged
#'   graph. 
#' \item When \emph{flatten_mode} = \emph{"compact"}, the default, nodes 
#'   and edges are converted into a flattened version to produce a compact 
#'   view of the temporal network whilst still presenting all the information
#'   in the export.\cr
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 become respectively X->Y lag=1, 
#'   X<-Y lag=2.  
#' \item When \emph{flatten_mode} = \emph{"combine"}, prior to the export,
#'   a preprocessing will be applied to kept only one edge
#'   per couple of nodes. The log_confidence will be the highest one 
#'   of the summarized edges whilst the lag and orientation of the
#'   summarized edge will be an aggregation.\cr 
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 will become X<->Y lag=1,2 with
#'   the log_confidence of X_lag1->Y_lag0 if log_confidence of 
#'   X_lag1->Y_lag0 > X_lag2<-Y_lag0.
#' \item When \emph{flatten_mode} = \emph{"unique"}, prior to the export,
#'   a preprocessing will be applied to kept only the edges having the
#'   highest log_confidence for a couple of nodes. 
#'   If several edges between the sames nodes have the same
#'   log_confidence, then the edge kept is the one with the minimum lag.\cr 
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 with log_confidence of 
#'   X_lag1->Y_lag0 > X_lag2<-Y_lag0 become X->Y lag=1.
#' \item When \emph{flatten_mode} = \emph{"drop"}, prior to the export,
#'   a preprocessing will be applied to kept only the edges having the 
#'   highest log_confidence for a couple of nodes. 
#'   If several edges between the sames nodes have the same
#'   log_confidence, then the edge kept is the one with the minimum lag.\cr
#'   i.e. :  X_lag1->Y_lag0, X_lag2<-Y_lag0 with log_confidence of 
#'   X_lag1->Y_lag0 > X_lag2<-Y_lag0 become X->Y. 
#'   The lag information is dropped during the preprocessing and 
#'   will not be exported.
#' }
#' 
#' @param show_self_loops [a boolean] Optional, TRUE by default.
#' When TRUE, the edges like X_lag0-X_lag1 are exported.
#' When FALSE, only edges having different nodes are exported.
#' 
#' @param method A string representing the plotting method.
#' Currently only "igraph" is supported.
#'
#' @return A graph object adapted to the method.
#'
#' @export
#'
#' @seealso
#' \code{\link{tmiic.getIgraph}} for details on the igraph exported object.
#'
#' @examples
#' \donttest{
#' library(miic)
#' data(covidCases)
#' # execute MIIC (reconstruct graph in temporal mode)
#' tmiic.res <- miic(input_data = covidCases, tau = 2, movavg = 14)
#'
#' # Plot default compact temporal network Using igraph
#' if(require(igraph)) {
#' g = tmiic.export(tmiic.res, method="igraph")
#' plot(g) # Default visualisation, calls igraph::plot.igraph()
#'
#' # Plot full temporal network Using igraph
#' if(require(igraph)) {
#' g = tmiic.export(tmiic.res, flatten_mode="none", method="igraph")
#' plot(g) # Default visualisation, calls igraph::plot.igraph()
#'
#' # Specifying layout (see ?igraph::layout_)
#' l <- layout_on_grid(g, width = 5, height = 3, dim = 2)
#' plot(g, layout=l)
#' 
#' # Override some graphical parameters
#' plot(g, edge.arrow.size = 0.75)
#' plot(g, vertex.shape="none", edge.color="gray85", vertex.label.color="gray10")
#' 
#' # For compact graphs, please be aware that the rendering of
#' # igraph::plot.igraph() is not optimal when the graph contains
#' # multiple edges between the same nodes.
#' # So the recommend way to plot a compact graph is to use tmiic plotting:  
#' plot(tmiic.res)
#' }
#'
#' }
#-----------------------------------------------------------------------------
tmiic.export <- function (tmiic.res, flatten_mode="compact", 
                          show_self_loops=TRUE, method="igraph") {
  if (is.null(tmiic.res$all.edges.summary)) {
    stop("The inferred network does not exist")
  }
  if (is.null(method)) {
    stop("Plotting method is required")
  }
  if (method != "igraph") {
    stop("Method not supported")
  }
  return(tmiic.getIgraph(tmiic.res, flatten_mode=flatten_mode, 
                         show_self_loops=show_self_loops))
}

#-----------------------------------------------------------------------------
# tmiic.getIgraph
#-----------------------------------------------------------------------------
#' Igraph plotting function for tmiic (temporal mode of miic)
#'
#' @description This functions returns an igraph object built from the result
#' returned by:
#' \code{\link{miic}} executed in temporal temporal mode.
#' \code{\link{tmiic.flatten_network}}
#' 
#' @details
#' Edges attributes are passed to the igraph graph and can be accessed with
#' e.g. \code{E(g)$partial_correlation}. See \code{\link{miic}} for more
#' details on edge parameters. By default, edges are colored according to the
#' partial correlation between two nodes conditioned on the conditioning set
#' (negative is blue, null is gray and positive is red) and their width is
#' based on the conditional mutual information minus the complexity cost.
#'
#' @param tmiic.res [a tmiic graph object]
#' The graph object returned by the \code{\link{miic}} execution in temporal mode
#'
#' @param flatten_mode [a string]. Optional, default value "compact".
#' Possible values are \emph{"none"}, \emph{"compact"}, \emph{"combine"}, 
#' \emph{"unique"}, \emph{"drop"}:
#' \itemize{
#' \item When \emph{flatten_mode} = \emph{"none"}, the function will
#'   use the tmiic graph object as it, leading to the return of a lagged
#'   graph.  
#' \item When \emph{flatten_mode} = \emph{"compact"}, the default, nodes 
#'   and edges are converted into a flattened version to produce a compact 
#'   view of the temporal network whilst still presenting all the information.\cr
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 become respectively X->Y lag=1, 
#'   X<-Y lag=2.  
#' \item When \emph{flatten_mode} = \emph{"combine"}, 
#'   a preprocessing will be applied to kept only one edge
#'   per couple of nodes. The log_confidence will be the highest one 
#'   of the summarized edges whilst the lag and orientation of the
#'   summarized edge will be an aggregation.\cr 
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 will become X<->Y lag=1,2 with
#'   the log_confidence of X_lag1->Y_lag0 if log_confidence of 
#'   X_lag1->Y_lag0 > X_lag2<-Y_lag0.
#' \item When \emph{flatten_mode} = \emph{"unique"}, 
#'   a preprocessing will be applied to kept only the edges having the
#'   highest log_confidence for a couple of nodes. 
#'   If several edges between the sames nodes have the same
#'   log_confidence, then the edge kept is the one with the minimum lag.\cr 
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 with log_confidence of 
#'   X_lag1->Y_lag0 > X_lag2<-Y_lag0 become X->Y lag=1.
#' \item When \emph{flatten_mode} = \emph{"drop"}, prior to the plotting,
#'   a preprocessing will be applied to kept only the edges having the 
#'   highest log_confidence for a couple of nodes. 
#'   If several edges between the sames nodes have the same
#'   log_confidence, then the edge kept is the one with the minimum lag.\cr
#'   i.e. :  X_lag1->Y_lag0, X_lag2<-Y_lag0 with log_confidence of 
#'   X_lag1->Y_lag0 > X_lag2<-Y_lag0 become X->Y. 
#'   The lag information is dropped during the preprocessing.
#' }
#' 
#' @param show_self_loops [a boolean] Optional, TRUE by default.
#' When TRUE, the edges like X_lag0-X_lag1 are included in the iGraph object. 
#' When FALSE, only edges having different nodes are present in the iGraph 
#' object.
#' 
#' @return An igraph graph object.
#' 
#' @export
#'
#' @seealso
#' \code{\link{miic}} for details on edge parameters in the returned object,
#' \code{\link[igraph]{igraph.plotting}} for the detailed description of the
#' plotting parameters and \code{\link[igraph]{layout}} for different layouts.
#-----------------------------------------------------------------------------
tmiic.getIgraph <- function (tmiic.res, flatten_mode="compact", show_self_loops=TRUE){
  if (class(tmiic.res) != "tmiic") {
    stop("Not a tmiic object.")
  }
  
  if (flatten_mode != "none")
    tmiic.res <- tmiic.flatten_network(tmiic.res, flatten_mode=flatten_mode,
                                       keep_edges_on_same_node=show_self_loops)
  
  if (tmiic.res$tmiic_specific[["graph_type"]] != "lagged") {
    list_nodes <- tmiic.res$tmiic_specific[["nodes_not_lagged"]]
    tmiic.res$adj_matrix <- matrix(NA, nrow=0, ncol=length (list_nodes)) 
    colnames(tmiic.res$adj_matrix) <- list_nodes
  }
  graph <- getIgraph(tmiic.res)
  
  if (tmiic.res$tmiic_specific[["graph_type"]] == "lagged") {
    igraph::V(graph)$label.dist = 1
    igraph::V(graph)$label.degree = pi/2
    igraph::E(graph)$curved = TRUE
  }
  else {
    igraph::E(graph)$curved = FALSE
    if ( "lag" %in% colnames(tmiic.res$all.edges.summary) )
      igraph::E(graph)$label <- tmiic.res$all.edges.summary$lag
  }
  return(graph)
}
 
#-----------------------------------------------------------------------------
# tmiic.prepareEdgesForPlotting
#-----------------------------------------------------------------------------
# Prepare the edges for plotting
#
# @description This function firstly filters the edges in the summary to keep
# only the ones detected by miic and adds to every edge an id constructed
# using the couple of nodes ordered alphabetically ("node1" < "node2") 
# 
# @param  [a tmiic graph object] The graph object returned by the miic 
# execution in temporal mode, eventually flattened
# 
# @return tmiic.res [a tmiic object] The modified tmiic object
#-----------------------------------------------------------------------------
tmiic.prepareEdgesForPlotting <- function (tmiic.res) {
  df_edges <- tmiic.res$all.edges.summary[tmiic.res$all.edges.summary$type %in% c('P', 'TP', 'FP'), ]
  #
  # Ensure all edges have an id xy where x < y
  #
  df_edges$xy = NULL
  for (edge_idx in 1:nrow(df_edges)) {
    one_edge <- df_edges[edge_idx,]
    if (one_edge$x < one_edge$y) 
      df_edges[edge_idx, "xy"] <- paste (one_edge$x, "-", one_edge$y, sep="")
    else
      df_edges[edge_idx, "xy"] <- paste (one_edge$y, "-", one_edge$x, sep="")
  }
  #
  # Order the edges so that all orientations goes from x to y
  #
  for(row in 1:nrow(df_edges)){
    if(df_edges[row, "infOrt"] == -2){
      df_edges[row, c("x","y")] = df_edges[row, c("y","x")]
      df_edges[row, "infOrt"] = 2
      if(!is.na(df_edges[row, "proba"])){
        df_edges[row, "proba"] = paste0(rev(
          strsplit(df_edges[row, "proba"], ";")[[1]]), collapse=";")
      }
      if(!is.na(df_edges[row, "trueOrt"])){
        df_edges[row, "trueOrt"] = 2
      }
    }
  }
  
  tmiic.res$all.edges.summary <- df_edges
  return (tmiic.res)
} 

#-----------------------------------------------------------------------------
# tmiic.getMultipleEdgesForPlotting
#-----------------------------------------------------------------------------
# Look for mutiple edges (that needs specific plotting)
#
# @description This function identifies the couple of nodes having mutiples 
# edges
# 
# @param  [a tmiic graph object]
# The graph object returned by the miic execution in temporal mode and
# flattened (if the tmiic object is not flattened, the function does nothing) 
# 
# @return df_mult [a dataframe] The dataframe containing the multiple edges
#-----------------------------------------------------------------------------
tmiic.getMultipleEdgesForPlotting <- function (tmiic.res) {
  df_mult <- tmiic.res$all.edges.summary
  df_mult$count <- 1
  df_mult <- stats::aggregate(data.frame(count = df_mult$count), 
                              by = list(xy = df_mult$xy), sum)
  df_mult <- df_mult[df_mult$count > 1,]
  return (df_mult)
} 
  
#-----------------------------------------------------------------------------
# plot.tmiic
#-----------------------------------------------------------------------------
#' Basic plot function of a tmiic network inference result
#'
#' @description This function calls \code{\link{tmiic.export}} to build a
#' plottable object from the result returned by \code{\link{miic}} in 
#' temporal mode and plot it.
#'
#' @details See the documentation of \code{\link{tmiic.export}} for further
#' details.
#'
#' @param x [a tmiic graph object]
#' The graph object returned by \code{\link{miic}} in temporal mode
#' 
#' @param flatten_mode [a string]. Optional, default value "compact".
#' Possible values are \emph{"none"}, \emph{"compact"}, \emph{"combine"}, 
#' \emph{"unique"}, \emph{"drop"}:
#' \itemize{
#' \item When \emph{flatten_mode} = \emph{"none"}, the plot function will
#'   use the tmiic graph object as it, leading to the display of a lagged
#'   graph. Unless a specific layout is specified, nodes will be positioned 
#'   on a grid. 
#' \item When \emph{flatten_mode} = \emph{"compact"}, the default, nodes 
#'   and edges are converted into a flattened version to produce a compact 
#'   view of the temporal network whilst still presenting all the information
#'   in the plotting.\cr
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 become respectively X->Y lag=1, 
#'   X<-Y lag=2.  
#' \item When \emph{flatten_mode} = \emph{"combine"}, prior to the plotting,
#'   a preprocessing will be applied to kept only one edge
#'   per couple of nodes. The log_confidence will be the highest one 
#'   of the summarized edges whilst the lag and orientation of the
#'   summarized edge will be an aggregation.\cr 
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 will become X<->Y lag=1,2 with
#'   the log_confidence of X_lag1->Y_lag0 if log_confidence of 
#'   X_lag1->Y_lag0 > X_lag2<-Y_lag0.
#' \item When \emph{flatten_mode} = \emph{"unique"}, prior to the plotting,
#'   a preprocessing will be applied to kept only the edges having the
#'   highest log_confidence for a couple of nodes. 
#'   If several edges between the sames nodes have the same
#'   log_confidence, then the edge kept is the one with the minimum lag.\cr 
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 with log_confidence of 
#'   X_lag1->Y_lag0 > X_lag2<-Y_lag0 become X->Y lag=1.
#' \item When \emph{flatten_mode} = \emph{"drop"}, prior to the plotting,
#'   a preprocessing will be applied to kept only the edges having the 
#'   highest log_confidence for a couple of nodes. 
#'   If several edges between the sames nodes have the same
#'   log_confidence, then the edge kept is the one with the minimum lag.\cr
#'   i.e. :  X_lag1->Y_lag0, X_lag2<-Y_lag0 with log_confidence of 
#'   X_lag1->Y_lag0 > X_lag2<-Y_lag0 become X->Y. 
#'   The lag information is dropped during the preprocessing and 
#'   will not be displayed on the final plotting.
#' }
#' 
#' @param show_self_loops [a boolean] Optional, TRUE by default.
#' When TRUE, the edges like X_lag0-X_lag1 are included in the iGraph object. 
#' When FALSE, only edges having different nodes are present in the iGraph 
#' object.
#' 
#' @param method A string representing the plotting method. Default to "igraph".
#' Currently only "igraph" is supported.
#' 
#' @param \dots Additional plotting parameters. See the corresponding plot function
#' for the complete list.
#' 
#' For igraph, see \code{\link[igraph]{igraph.plotting}}.
#'
#' @export
#'
#' @seealso \code{\link{tmiic.export}} for generic exports,
#' \code{\link{tmiic.getIgraph}} for igraph export,
#' \code{\link[igraph]{igraph.plotting}}
#' 
#' @examples
#' library(miic)
#'
#' #' # EXAMPLE COVID CASES (timeseries demo)
#' data(covidCases)
#' # execute MIIC (reconstruct graph in temporal mode)
#' tmiic.res <- miic(input_data = covidCases, tau = 2, movavg = 14)
#'
#' # to plot the default compact graph
#' if(require(igraph)) {
#'  plot(tmiic.res)
#'  
#' # to plot the full temporal graph
#' if(require(igraph)) {
#'  plot(tmiic.res, flatten_mode="none")
#' }
#'
#' }
#-----------------------------------------------------------------------------
plot.tmiic = function(x, flatten_mode="compact", show_self_loops=TRUE,
                      method = 'igraph', ...) {
  
  if (class(x) != "tmiic")
    stop("Not a tmiic object.")
  if (method != 'igraph')
    stop("Method not supported. See ?tmiic.export for supported methods.")
  if ( !base::requireNamespace("igraph", quietly = TRUE) ) 
    stop("Package 'igraph' is required.")
  if ( is.null (x$adj_matrix) ) 
    stop ("The learnt graphical model adjacency matrix does not exist")
  
  if ( (x$tmiic_specific[["graph_type"]] == "lagged") & (flatten_mode != "none") )
    x <- tmiic.flatten_network(x, flatten_mode=flatten_mode, 
                               keep_edges_on_same_node=show_self_loops)
  
  x <- tmiic.prepareEdgesForPlotting(x)
  
  df_mult <- data.frame(count=integer(), stringsAsFactors = FALSE)
  if (x$tmiic_specific[["graph_type"]] == "compact")
    df_mult <- tmiic.getMultipleEdgesForPlotting(x)
  #
  # Set a layout if none supplied by user : grid for lagged, 
  # layout_with_kk for flatten
  #
  layout <- NULL
  if ( ! ("layout" %in% names(list(...))) ) {
    if (x$tmiic_specific[["graph_type"]] == "lagged") {
      n_nodes_not_lagged <- length(x$tmiic_specific[["nodes_not_lagged"]])
      list_nodes <- colnames (x$adj_matrix)
      max_lag_plus1 <- length(list_nodes) %/% n_nodes_not_lagged
      layout <- data.frame (rep(1:max_lag_plus1, each=n_nodes_not_lagged),
                            rep(1:n_nodes_not_lagged, times=max_lag_plus1) )
      layout <- as.matrix (layout)
    }
    else {
      layout <- igraph::layout_with_kk
    }
  }
  #
  # Export the graph to a graphical object and plot
  #
  graph <- tmiic.export (x, flatten_mode="none", method=method)
  if (nrow (df_mult) <= 0) {
    # No multiple edges between the same nodes, we draw in one go
    #
    if ( is.null(layout) )
      igraph::plot.igraph (graph, ...)
    else
      igraph::plot.igraph (graph, layout=layout, ...)
  }
  else {
    # Multiple edges between the same nodes exist, draw iteratively
    #
    df_edges <- x$all.edges.summary
    edges_colors_iter <- igraph::E(graph)$color
    edges_labels_iter <- igraph::E(graph)$label
    #
    # On a first step, we will draw all the graph except multiple edges
    # The multiple edges will be drawn with invisible color "#FF000000"
    # and with no labels 
    #
    for ( edge_idx in 1:nrow(df_edges) ) {
      one_edge <- df_edges[edge_idx,]
      if (one_edge$xy %in% df_mult$xy) {
        edges_colors_iter[[edge_idx]] <- "#FF000000"
        edges_labels_iter[[edge_idx]] <- NA
      }
    }
    igraph::plot.igraph (graph, layout=layout,
                         edge.color=edges_colors_iter, 
                         edge.label=edges_labels_iter, ...)
    #
    # Draw each group of multiple edges
    #
    for ( mult_idx in 1:nrow(df_mult) ) {
      one_mult <- df_mult[mult_idx,]
      nodes_of_mult <- strsplit (one_mult$xy, "-")[[1]]
      if (nodes_of_mult[[1]] == nodes_of_mult[[2]]) {
        #
        # for self loop, we will go over 2*pi around the node
        #
        step_pos <- 0
        step_inc <- (2 * pi) / one_mult$count
      }
      else {
        #
        # otherelse, we will curve edges from -0.5 to +0.5
        #
        if (one_mult$count > 4) {
          #
          # if more than 4 edges, curve more 
          #
          step_pos <- -1
          step_inc <- 2.0 / (one_mult$count - 1)
        }
        else {
          step_pos <- -0.5
          step_inc <- 1.0 / (one_mult$count - 1)
        }
      }
      #
      # Draw multiple egdes one by one
      #
      list_to_draw = which(df_edges[, "xy"] == one_mult$xy)
      for (idx_to_draw in 1:length(list_to_draw) ) {
        edge_to_draw = list_to_draw[[idx_to_draw]]
        #
        # We hide all edges except one
        #
        edges_colors_iter <- rep ("#FF000000", nrow (df_edges) )
        edges_labels_iter <- rep (NA, nrow (df_edges) )
        edges_colors_iter[[edge_to_draw]] <- igraph::E(graph)[[edge_to_draw]]$color
        edges_labels_iter[[edge_to_draw]] <- igraph::E(graph)[[edge_to_draw]]$label

        if (nodes_of_mult[[1]] == nodes_of_mult[[2]])
          igraph::plot.igraph (graph, layout=layout, add=TRUE,
                               edge.color=edges_colors_iter, 
                               edge.label=edges_labels_iter, 
                               edge.loop.angle=step_pos, ...)
        else {
          if (df_edges[edge_to_draw,]$x < df_edges[edge_to_draw,]$y)
            igraph::plot.igraph (graph, layout=layout, add=TRUE,
                                 edge.color=edges_colors_iter, 
                                 edge.label=edges_labels_iter, 
                                 edge.curved=step_pos, ...)
          else # y > x we need to curve the edge on the opposite way
            igraph::plot.igraph (graph, layout=layout, add=TRUE,
                                 edge.color=edges_colors_iter, 
                                 edge.label=edges_labels_iter, 
                                 edge.curved=-step_pos, ...)
        }
        #
        # Update position for next edge
        #
        step_pos <- step_pos + step_inc
      }
    }
  }
}


