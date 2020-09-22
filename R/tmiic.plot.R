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
# - 05 oct 2020 : include the flattening as a parameter 
# - 06 oct 2020 : addition of "lagged" graph (edges duplicated over history)
# - 10 fev 2021 : plotting modified to manage contextual variables
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
#' @param display [a string]. Optional, default value "compact".
#' Possible values are \emph{"raw"}, \emph{"lagged"}, \emph{"compact"}, 
#' \emph{"combine"}, \emph{"unique"}, \emph{"drop"}:
#' \itemize{
#' \item When \emph{display} = \emph{"raw"}, the export function will
#'   use the tmiic graph object as it, leading to the return of a lagged
#'   graph. 
#' \item When \emph{display} = \emph{"lagged"}, the export function will
#'   repeat the edges over history assuming stationarity and return a lagged
#'   graph. 
#' \item When \emph{display} = \emph{"compact"}, the default, nodes 
#'   and edges are converted into a flattened version to produce a compact 
#'   view of the temporal network whilst still presenting all the information
#'   in the export.\cr
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 become respectively X->Y lag=1, 
#'   X<-Y lag=2.  
#' \item When \emph{display} = \emph{"combine"}, prior to the export,
#'   a preprocessing will be applied to kept only one edge
#'   per couple of nodes. The info_shifted will be the highest one 
#'   of the summarized edges whilst the lag and orientation of the
#'   summarized edge will be an aggregation.\cr 
#'   i.e.: X_lag2->Y_lag0, X_lag0<-Y_lag1 will become X<->Y lag=1-2 with
#'   the info_shifted of X_lag2->Y_lag0 if info_shifted of 
#'   X_lag2->Y_lag0 > X_lag0<-Y_lag1.
#' \item When \emph{display} = \emph{"unique"}, prior to the export,
#'   a preprocessing will be applied to kept only the edges having the
#'   highest info_shifted for a couple of nodes. 
#'   If several edges between the sames nodes have the same
#'   info_shifted, then the edge kept is the one with the minimum lag.\cr 
#'   i.e.: X_lag1->Y_lag0, X_lag0<-Y_lag2 with info_shifted of 
#'   X_lag1->Y_lag0 > X_lag0<-Y_lag2 become X->Y lag=1.
#' \item When \emph{display} = \emph{"drop"}, prior to the export,
#'   a preprocessing will be applied to kept only the edges having the 
#'   highest info_shifted for a couple of nodes. 
#'   If several edges between the sames nodes have the same
#'   info_shifted, then the edge kept is the one with the minimum lag.\cr
#'   i.e. :  X_lag1->Y_lag0, X_lag0<-Y_lag2 with info_shifted of 
#'   X_lag1->Y_lag0 > X_lag0<-Y_lag2 become X->Y. 
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
#' # Plot raw temporal network Using igraph
#' g = tmiic.export(tmiic.res, display="raw", method="igraph")
#' plot(g) # Default visualisation, calls igraph::plot.igraph()
#'
#' # Plot full temporal network Using igraph
#' g = tmiic.export(tmiic.res, display="lagged", method="igraph")
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
tmiic.export <- function (tmiic.res, display="compact", 
                          show_self_loops=TRUE, method="igraph") {
  if (is.null(tmiic.res$all.edges.summary)) {
    stop("Error: The inferred network does not exist")
  }
  if (is.null(method)) {
    stop("Error: Plotting method is required")
  }
  if (method != "igraph") {
    stop("Error: Method not supported")
  }
  return(tmiic.getIgraph(tmiic.res, display=display, 
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
#' @param display [a string]. Optional, default value "compact".
#' Possible values are \emph{"raw"}, \emph{"lagged"}, \emph{"compact"}, 
#' \emph{"combine"}, \emph{"unique"}, \emph{"drop"}:
#' \itemize{
#' \item When \emph{display} = \emph{"raw"}, the function will
#'   use the tmiic graph object as it, leading to the return of a lagged
#'   graph.  
#' \item When \emph{display} = \emph{"lagged"}, the function will
#'   repeat the edges over history assuming stationarity and return a lagged
#'   graph. 
#' \item When \emph{display} = \emph{"compact"}, the default, nodes 
#'   and edges are converted into a flattened version to produce a compact 
#'   view of the temporal network whilst still presenting all the information.\cr
#'   i.e.: X_lag1->Y_lag0, X_lag0<-Y_lag2 become respectively X->Y lag=1, 
#'   X<-Y lag=2.  
#' \item When \emph{display} = \emph{"combine"}, 
#'   a preprocessing will be applied to kept only one edge
#'   per couple of nodes. The info_shifted will be the highest one 
#'   of the summarized edges whilst the lag and orientation of the
#'   summarized edge will be an aggregation.\cr 
#'   i.e.: X_lag2->Y_lag0, X_lag0<-Y_lag1 will become X<->Y lag=1,2 with
#'   the info_shifted of X_lag2->Y_lag0 if info_shifted of 
#'   X_lag2->Y_lag0 > X_lag0<-Y_lag1.
#' \item When \emph{display} = \emph{"unique"}, 
#'   a preprocessing will be applied to kept only the edges having the
#'   highest info_shifted for a couple of nodes. 
#'   If several edges between the sames nodes have the same
#'   info_shifted, then the edge kept is the one with the minimum lag.\cr 
#'   i.e.: X_lag1->Y_lag0, X_lag0<-Y_lag2 with info_shifted of 
#'   X_lag1->Y_lag0 > X_lag0<-Y_lag2 become X->Y lag=1.
#' \item When \emph{display} = \emph{"drop"}, prior to the plotting,
#'   a preprocessing will be applied to kept only the edges having the 
#'   highest info_shifted for a couple of nodes. 
#'   If several edges between the sames nodes have the same
#'   info_shifted, then the edge kept is the one with the minimum lag.\cr
#'   i.e. :  X_lag1->Y_lag0, X_lag0<-Y_lag2 with info_shifted of 
#'   X_lag1->Y_lag0 > X_lag0<-Y_lag2 become X->Y. 
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
tmiic.getIgraph <- function (tmiic.res, display="compact", show_self_loops=TRUE){
  if (class(tmiic.res) != "tmiic") {
    stop("Error: Not a tmiic object.")
  }
  
  if (display == "lagged") {
    tmiic.res <- tmiic.repeat_edges_over_history(tmiic.res)
  }
  else {
    if (display != "raw")
      tmiic.res <- tmiic.flatten_network(tmiic.res, flatten_mode=display,
                                       keep_edges_on_same_node=show_self_loops)
  }
  if ( ! tmiic.res$tmiic_specific[["graph_type"]] %in% c("raw", "lagged") ) {
    list_nodes <- tmiic.res$tmiic_specific[["nodes_not_lagged"]]
    tmiic.res$adj_matrix <- matrix(NA, nrow=0, ncol=length (list_nodes)) 
    colnames(tmiic.res$adj_matrix) <- list_nodes
  }
  graph <- getIgraph(tmiic.res)
  
  if ( tmiic.res$tmiic_specific[["graph_type"]] %in% c("raw", "lagged") ) {
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
  if (nrow(df_edges) <= 0) {
    df_edges$xy = character(0)
  }
  else {
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
  if (nrow(df_mult) <= 0)
    df_mult$count <- numeric(0)
  else {
    df_mult$count <- 1
    df_mult <- stats::aggregate(data.frame(count = df_mult$count), 
                                by = list(xy = df_mult$xy), sum)
  }
  df_mult <- df_mult[df_mult$count > 1,]
  return (df_mult)
} 

#-----------------------------------------------------------------------------
# tmiic.computeRowLayoutGreedyBase
#-----------------------------------------------------------------------------
# Internal function to compute an optimized raw layout used to construct
# a pseudo grid layout for the display of raw and lagged graphs
#
# @description The function starts by choosing two nodes: the node with the 
# maximum degree and the one sharing the most of edges with the first. 
# Then, the other nodes are placed recurvely using these two nodes.
# If some nodes are still not positioned after  the recursion, 
# the whole process is done over until all nodes are placed.
# 
# @param  list_nodes [a list] The list of nodes to be positioned
# @param  df_edges [a dataframe] The list of edges, with the edges nodes
# stored in columns x and y and count columns (>1 when edge exists with 
# multiple lags between the two nodes)
# 
# @return [a list] The list of nodes ordered to avoid crossing edges when 
# the network is displayed as raw or lagged graphs
#-----------------------------------------------------------------------------
tmiic.computeRowLayoutGreedyBase <- function (list_nodes, df_edges)
  {
  if (length (list_nodes) <= 0)
    return ( list() )
  if (nrow(df_edges) <= 0)
    return (list_nodes)
  #
  # Count nodes degrees
  #
  df_nodes <- data.frame (nodes=unlist(list_nodes) )
  df_nodes$degree <- 0
  for (i in 1:nrow(df_nodes))
    df_nodes[i,2] <- sum (  (df_edges$x == df_nodes[i,1])
                          | (df_edges$y == df_nodes[i,1]) )
  #
  # Select first node with the max degree
  #
  max_degree <- max (df_nodes$degree)
  node_chosen <- df_nodes[ (df_nodes$degree == max_degree), 1]
  node_chosen <- node_chosen[[1]]
  #
  # Select node having the maximum number of edges with the max degree one
  #
  cond_rel_with_chosen <- ( (df_edges$x == node_chosen) 
                          | (df_edges$y == node_chosen) )
  df_rel_with_chosen <- df_edges[cond_rel_with_chosen,]
  max_rel_with_chosen <- max (df_rel_with_chosen$count)
  edges_max_rel_with_chosen <- df_rel_with_chosen[ df_rel_with_chosen$count == max_rel_with_chosen,]
  edge_max_rel_with_chose <- edges_max_rel_with_chosen[1,] 
  if (edge_max_rel_with_chose$x == node_chosen)
    other_node = edge_max_rel_with_chose$y
  if (edge_max_rel_with_chose$y == node_chosen)
    other_node = edge_max_rel_with_chose$x
  #
  # Remove the two selected nodes from the lists of nodes and edges
  #
  cond_edge_chosen_other = (  (df_edges$x == node_chosen) & (df_edges$y == other_node)
                            | (df_edges$y == node_chosen) & (df_edges$x == other_node) )
  df_edges <- df_edges[(!cond_edge_chosen_other),] 
  
  cond_node_chosen_other <- (  (df_nodes$nodes == node_chosen)
                             | (df_nodes$nodes == other_node) )
  df_nodes <- df_nodes[ (!cond_node_chosen_other), ]
  #
  # Compute recursively the positions of nodes in regard of the two selected
  #
  ret <- tmiic.tmiic.computeRowLayoutGreedyRecurs (node_chosen, other_node, 
                                      df_nodes$nodes, df_nodes, df_edges)  
  #
  # If some nodes are still not positioned, do separate graph(s) beside
  #
  while (! is.null (ret$nodes_no_care) )
    {
    #
    # Remove all edgees havbing their two nodes positioned
    #
    cond <- ( (df_edges$x %in% ret$nodes_positioned) 
            & (df_edges$y %in% ret$nodes_positioned) )
    df_edges <- df_edges[ (!cond), ]
    #
    # Construct a separate layout with remaining edges
    #
    ret_next <- tmiic.computeRowLayoutGreedyBase (ret$nodes_no_care, df_edges)  
    ret$nodes_positioned <- c(ret$nodes_positioned, ret_next$nodes_positioned)
    ret$nodes_no_care <- ret_next$nodes_no_care
    }
  return (ret)
  }

#-----------------------------------------------------------------------------
# tmiic.tmiic.computeRowLayoutGreedyRecurs
#-----------------------------------------------------------------------------
# Internal recursive function to compute an optimized raw layout used to 
# construct a pseudo grid layout for the display of raw and lagged graphs
#
# @description The function starts by using two nodes as separators:
# the other nodes are placed in sets depending on how they are placed 
# regarding the two separators: left, center and rigth. These sets are 
# processed in a recursive way until becoming empty, then the backtrack 
# generates the lit of nodes representing the raw layout with minimal
# crossing
# 
# @param node_left [a string] The left separator node
# @param node_right [a string] The right separator node
# @param list_nodes_to_affect[a list] The list is nodes to be postioned
# @param df_nodes [a dataframe] The list of nodes with columns nodes and degree
# @param df_edges [a dataframe] The list of edges with columns x, y containing
# the nodes and count columns (>1 when edge exists with multiple lags between 
# the two nodes)
# 
# @return [a list] The list of nodes ordered to avoid crossing edges
#-----------------------------------------------------------------------------
tmiic.tmiic.computeRowLayoutGreedyRecurs <- function (node_left, node_right, 
                                    list_nodes_to_affect, df_nodes, df_edges)
  {
  #
  # Remove from nodes and edges the right and left nodes that we chose to 
  # position the others 
  #
  cond_node_right_left <- (  (list_nodes_to_affect == node_left)
                           | (list_nodes_to_affect == node_right) )
  list_nodes_to_affect <- list_nodes_to_affect[ (!cond_node_right_left)]
  
  cond_edge_right_left = (  (df_edges$x == node_left) & (df_edges$y == node_right)
                          | (df_edges$y == node_left) & (df_edges$x == node_right) )
  df_edges <- df_edges[(!cond_edge_right_left),] 
  
  cond_node_right_left <- (  (df_nodes$nodes == node_left)
                           | (df_nodes$nodes == node_right) )
  df_nodes <- df_nodes[ (!cond_node_right_left), ]
  #
  # Position the other nodes compared with the right and left chosen nodes
  #
  nodes_left <- list()
  nodes_center <- list()
  nodes_right <- list()
  nodes_no_care <- list()
  for (one_node in list_nodes_to_affect)
    {
    cond_rel_with_left <- any ( ( (df_edges$x == one_node) | (df_edges$y == one_node) )
                              & ( (df_edges$x == node_left) | (df_edges$y == node_left) ) )
    cond_rel_with_right <- any ( ( (df_edges$x == one_node) | (df_edges$y == one_node) )
                               & ( (df_edges$x == node_right) | (df_edges$y == node_right) ) )
    cond_rel_with_both <- cond_rel_with_left & cond_rel_with_right
    
    if (cond_rel_with_both)
      {
      nodes_center[[(length(nodes_center)+1)]] <- one_node
      next
      }
    if (cond_rel_with_left)
      {
      nodes_left[[length(nodes_left)+1]] <- one_node
      next
      }
    if (cond_rel_with_right)
      {
      nodes_right[[length(nodes_right)+1]]<- one_node
      next
      }
    nodes_no_care[[length(nodes_no_care)+1]] <- one_node
    }
  #
  # If there is no interest to position some nodes, end recursion
  #
  if ( sum(length(nodes_left), length(nodes_center), length(nodes_right)) <= 0)
    {
    ret <- list (nodes_positioned=unlist(c(node_left, node_right)),
                 nodes_no_care=unlist(nodes_no_care) ) 
    return (ret)
    }
  #
  # There is some  interest to position some nodes
  #
  find_node_max_degre <- function (list_possible_nodes, df_nodes)
    {
    df_nodes <- df_nodes[(df_nodes$nodes %in% list_possible_nodes),]
    max_edges <- max(df_nodes$degree)
    ret_node <- df_nodes[(df_nodes$degree == max_edges),1]
    ret_node <- ret_node[[1]]
    return (ret_node)
    }
  
  nodes_positioned_left <- list()
  nodes_positioned_center_left <- list()
  nodes_positioned_center_right <- list()
  nodes_positioned_right <- list()
  if (length(nodes_left) > 0)
    {
    new_node_sep <- find_node_max_degre (nodes_left, df_nodes)
    ret <- tmiic.tmiic.computeRowLayoutGreedyRecurs (new_node_sep, node_left, 
                                          append (nodes_left, nodes_no_care), 
                                          df_nodes, df_edges)
    nodes_positioned_left  <- ret$nodes_positioned
    nodes_no_care <- ret$nodes_no_care
    df_nodes <- df_nodes[ (!df_nodes$nodes %in% nodes_positioned_left), ]
    }
  if (length(nodes_center) > 0)
    {
    new_node_sep <- find_node_max_degre (nodes_center, df_nodes)
    ret <- tmiic.tmiic.computeRowLayoutGreedyRecurs (node_left, new_node_sep, 
                                        append (nodes_center, nodes_no_care),
                                        df_nodes, df_edges)
    nodes_positioned_center_left <- ret$nodes_positioned
    df_nodes <- df_nodes[ (!df_nodes$nodes %in% nodes_positioned_center_left), ]
    
    ret <- tmiic.tmiic.computeRowLayoutGreedyRecurs (new_node_sep, node_right, 
                                                     ret$nodes_no_care,
                                                     df_nodes, df_edges)
    nodes_positioned_center_right  <- ret$nodes_positioned
    nodes_no_care <- ret$nodes_no_care
    df_nodes <- df_nodes[ (!df_nodes$nodes %in% nodes_positioned_center_right), ]
    }
  if (length(nodes_right) > 0)
    {
    new_node_sep <- find_node_max_degre (nodes_right, df_nodes)
    ret <- tmiic.tmiic.computeRowLayoutGreedyRecurs (node_right, new_node_sep, 
                                          append (nodes_right, nodes_no_care),
                                          df_nodes, df_edges)
    nodes_positioned_right <- ret$nodes_positioned
    nodes_no_care <- ret$nodes_no_care
    }
  #
  # Concat nodes that have been positioned the and return
  #
  nodes_pos_all <- c(nodes_positioned_left, node_left, nodes_positioned_center_left,
    nodes_positioned_center_right, node_right, nodes_positioned_right)
  ret <- list (nodes_positioned=unlist(nodes_pos_all),
               nodes_no_care=unlist(nodes_no_care) ) 
  return (ret)
  }

#-----------------------------------------------------------------------------
# tmiic.computeRowLayoutGreedy
#-----------------------------------------------------------------------------
# Internal function to compute an optimized pseudo grid layout for the display
# of raw and lagged graphs
#
# @description 
# The function counts edges per couple of nodes whatever their lags are and 
# exclude self loop. Then it call tmiic.computeRowLayoutGreedyBase
# with the nodes having at least one edge to compute a layer 0 layout.
# The layout is completed with nodes without edges to produce the final 
# layer 0 layout.
#
# @param  tmiic.res [a tmiic graph object] The graph object returned by the 
# miic execution in temporal mode, "raw" graph_type
# 
# @return [a list] The position along an axis for each node
#-----------------------------------------------------------------------------
tmiic.computeRowLayoutGreedy <- function (tmiic.res)
  {
  list_nodes_not_lagged <- tmiic.res$tmiic_specific[["nodes_not_lagged"]]
  #
  # Filter out self edges, count and summarize edges regardless their lags
  #
  tmiic.flat <- miic:::tmiic.flatten_network(tmiic.res)
  df_edges <- tmiic.flat$all.edges.summary
  df_edges <- df_edges[(df_edges$x != df_edges$y),]
  for (edge_idx in 1:nrow(df_edges)) 
    {
    one_edge <- df_edges[edge_idx,]
    if (one_edge$x >= one_edge$y) 
      df_edges[edge_idx, c("x","y")] <- c(one_edge$y, one_edge$x)
    }
  df_edges$count <- 1
  df_edges <- stats::aggregate(data.frame(count = df_edges$count), 
                                  by = list(x=df_edges$x, y=df_edges$y), sum)
  #
  # Find nodes not part of an edges or at least part of one edge
  #
  list_nodes_no_edges <- list()
  for (one_node in list_nodes_not_lagged)
    if ( (! one_node %in% df_edges$x) & (! one_node %in% df_edges$y) )
      list_nodes_no_edges[(length(list_nodes_no_edges)+1)] <- one_node
  
  list_nodes_with_edge <- list_nodes_not_lagged[ (!list_nodes_not_lagged %in% list_nodes_no_edges) ]
  #
  # Compute layer 0 layout (without nodes not part of an edges)
  #
  ret_recurs <- miic:::tmiic.computeRowLayoutGreedyBase (list_nodes_with_edge, df_edges)
  layout_unique_nodes = unique (ret_recurs$nodes_positioned)
  
  layout_row <- list()
  max_p1 <- length(layout_unique_nodes) + 1
  for (one_node in list_nodes_not_lagged)
    {
    layout_pos <- which (layout_unique_nodes == one_node)
    if (length(layout_pos) <= 0)
      {
      layout_row[[ (length(layout_row)+1) ]] <- max_p1
      max_p1 <- max_p1 + 1
      }
    else
      layout_row[[ (length(layout_row)+1) ]] = layout_pos[[1]]
    }
  return ( unlist (layout_row) )
  }

#-----------------------------------------------------------------------------
# tmiic.computeRowLayoutLayers
#-----------------------------------------------------------------------------
# Internal function to precompute a layout suited for the display of raw and 
# lagged graphs
#
# @description This function computes the layout so that the less layers 
# has a node, the more to the exteriors it will be placed.
#
# @param  tmiic.res [a tmiic graph object] The graph object returned by the 
# miic execution in temporal mode ("raw" graph_type)
# 
# @return [a list] The position along an axis for each node
#-----------------------------------------------------------------------------
tmiic.computeRowLayoutLayers <- function (tmiic.res)
  {
  n_nodes_not_lagged <- length(tmiic.res$tmiic_specific[["nodes_not_lagged"]])
  list_taus <- tmiic.res$tmiic_specific[["tau"]]
  list_taus [which (list_taus < 0)] <- 0
  tau_max <- max (list_taus)
  #
  # Precompute the rows on the grid, putting nodes with the less lags 
  # on the exteriors while the nodes having the most lags are in the center
  #
  list_pos_of_nodes <- list()
  idx_top <- 1
  idx_end <- n_nodes_not_lagged
  for (tau_ix in 0:tau_max) {
    list_nodes_idx_for_tau <- which(list_taus == tau_ix)
    if (length (list_nodes_idx_for_tau) > 0) {
      nb_top <- (length (list_nodes_idx_for_tau) + 1) %/% 2
      nb_end <- length (list_nodes_idx_for_tau) - nb_top
      i <- 1
      while (nb_top > 0) {
        node_idx <- list_nodes_idx_for_tau[[i]]
        list_pos_of_nodes[[node_idx]] <- idx_top
        idx_top <- idx_top + 1
        i <- i + 1
        nb_top <- nb_top - 1
      }
      if (nb_end > 0) {
        i <- length(list_nodes_idx_for_tau)
        while (nb_end > 0) {
          node_idx <- list_nodes_idx_for_tau[[i]]
          list_pos_of_nodes[[node_idx]] <- idx_end
          idx_end <- idx_end - 1
          i <- i - 1
          nb_end <- nb_end - 1
        } 
      }
    } 
  } 
  return (unlist (list_pos_of_nodes) )
  } 

#-----------------------------------------------------------------------------
# tmiic.computeRowLayoutSugiyama
#-----------------------------------------------------------------------------
# Internal function to precompute a layout suited for the display of raw and 
# lagged graphs
#
# @description This function computes the layout using Sugiyama algorihtm to
# minimize crossing edges
#
# @param  tmiic.res [a tmiic graph object] The graph object returned by the 
# miic execution in temporal mode ("raw" graph_type)
# 
# @return [a list] The position along an axis for each node
#-----------------------------------------------------------------------------
tmiic.computeRowLayoutSugiyama <- function (tmiic.res)
  {
  list_nodes_not_lagged <- tmiic.res$tmiic_specific[["nodes_not_lagged"]]
  n_nodes_not_lagged <- length(list_nodes_not_lagged)
  #
  # Filter out self edges, count and summarize edges regardless their lags
  #
  tmiic.flat <- miic:::tmiic.flatten_network(tmiic.res)
  df_edges <- tmiic.flat$all.edges.summary
  df_edges <- df_edges[(df_edges$x != df_edges$y),]
  for (edge_idx in 1:nrow(df_edges)) 
    {
    one_edge <- df_edges[edge_idx,]
    if (one_edge$x > one_edge$y) 
      df_edges[edge_idx, c("x","y")] <- c(one_edge$y, one_edge$x)
    }
  df_edges$count <- 1
  df_edges <- stats::aggregate(data.frame(count = df_edges$count), 
                                  by = list(x=df_edges$x, y=df_edges$y), sum)
  #
  # Create a dummy graph and apply Sugiyama algotrithm to get the layout
  #
  g_tmp <- igraph::graph_from_data_frame (df_edges, vertices=list_nodes_not_lagged)
  nodes_layers <- rep(1,n_nodes_not_lagged)
  edges_weight <- df_edges$count
  ret_sugiyama <- igraph::layout_with_sugiyama (g_tmp, layers=nodes_layers, 
                                     weights=edges_weight, attributes="all")
  list_pos_of_nodes <- ret_sugiyama$layout[,1]
  list_pos_of_nodes <- list_pos_of_nodes + 1
  return (list_pos_of_nodes)
  }
  
#-----------------------------------------------------------------------------
# tmiic.computeGridLayout
#-----------------------------------------------------------------------------
# Internal function to compute a pseudo grid layout to display raw and lagged 
# graphs
#
# @description Depending on the algorithm parameter, a layer 0 layout will
# be produced and be expanded over time.
# 
# @param  tmiic.res [a tmiic graph object]
# The graph object returned by the miic execution in temporal mode.
# The tmiic object can be used as it ("raw" graph_type) or modified into a 
# "lagged" graph by calling tmiic.repeat_edges_over_history.
# 
# @param display [a string]. Optional, default value "raw".
# Possible values are "raw" and "lagged".
# 
# @param  positioning [a string] Optional, default:"greedy". 
# The method used to position nodes. 
# Possible values are "none", "alphabetical", "layers", 
# "greedy" and "sugiyama":
# When positioning = "none": 
# - Th{ nodes are positioned as they appear in the miic result
# When positioning = "alphabetical":
# - The nodes are positioned alphabeticaly in ascending order
# When positioning = "layers":
# - The nodes with the less lags wil be placed on the exteriors 
#   while the nodes having the most lags are in the center
# When positioning = "greedy":
# - A greedy algorithm will be used to placed the nodes in a way minimizing
#   the crossing edges
# When positioning = "sugiyama":
# - The sugiyama algorithm will be used to placed the nodes in a way 
#   minimizing the crossing edges
# 
# @param  orientation [a character] Optional, default:"L". 
# The orientation of the draw. 
# Possible values are landscape:"L" or portrait: "P".
# 
# @return layout [a matrix] The layout to use for drawing
#-----------------------------------------------------------------------------
tmiic.computeGridLayout <- function (tmiic.res, display="raw", 
                                     positioning="greedy", orientation="L")
  {
  if (! display %in% c("raw", "lagged") )
    stop ("Error: Invalid display parameter")
  if (! positioning %in% c("none", "alphabetical", "layers", "greedy", "sugiyama") )
    stop ("Error: Invalid positioning parameter")
  if (! orientation %in% c("L", "P") )
    stop ("Error: Invalid orientation parameter")
  
  nodes_not_lagged <- tmiic.res$tmiic_specific[["nodes_not_lagged"]]
  n_nodes_not_lagged <- length(nodes_not_lagged)
  #
  # Precompute the layer 0 layout
  #
  list_pos_of_nodes <- NULL
  if (positioning == "none")
    list_pos_of_nodes = 1:n_nodes_not_lagged
  if (positioning == "alphabetical")
    {
    list_pos_of_nodes <- list()
    list_sorted <- sort (nodes_not_lagged)
    for (one_node in nodes_not_lagged)
      list_pos_of_nodes[[ (length(list_pos_of_nodes)+1) ]] <- which (list_sorted == one_node)[[1]]
    list_pos_of_nodes <- unlist (list_pos_of_nodes)
    }
  if (positioning == "layers")
    list_pos_of_nodes <- tmiic.computeRowLayoutLayers (tmiic.res)
  if (positioning == "greedy")
    list_pos_of_nodes <- tmiic.computeRowLayoutGreedy (tmiic.res)
  if (positioning == "sugiyama")
    list_pos_of_nodes <- tmiic.computeRowLayoutSugiyama (tmiic.res)
  if ( is.null (list_pos_of_nodes) )
    stop ("Error: Layout can not be infered")
  #
  # As contextual nodes are placed in an extra column/row when display is "raw",
  # here we update the nodes positions to maintain a "nice" display
  #
  is_contextual <- tmiic.res$tmiic_specific[["is_contextual"]]
  if (is.null(is_contextual))
    is_contextual <- rep (0, n_nodes_not_lagged)
  else
    is_contextual <- is_contextual[1:n_nodes_not_lagged]
    
  if ( (display == "raw") & (sum(is_contextual) > 0) )
    {
    list_pos_upd <- list_pos_of_nodes
    #
    # Identify contextual nodes
    #
    list_ctx_idx <- which (is_contextual != 0)
    n_ctx <- length (list_ctx_idx)
    #
    # Identify the order we need to follow to update postions
    #
    list_ctx_pos <- list_pos_of_nodes[list_ctx_idx]
    list_ctx_pos_order <- sort (list_ctx_pos)
    #
    # Distance between contextual nodes
    #
    max_pos <- n_nodes_not_lagged - n_ctx
    ctx_pos_shift <- max(1, max_pos / (n_ctx + 1) )
    #
    # Update the positions of the contextual node and shift the others
    #
    for (i in 1:n_ctx)
      {
      one_pos <- list_ctx_pos_order[[i]]
      node_idx <- which (list_pos_of_nodes == one_pos)
      list_pos_upd[node_idx] <- round(i * ctx_pos_shift, 0)
      #
      # Shift higher positions of non contextual nodes
      #
      node_shift <- i - 1
      pos_to_upd <- which ((is_contextual == 0) & (list_pos_upd >= one_pos - node_shift))
      list_pos_upd[pos_to_upd] <- list_pos_upd[pos_to_upd] - 1
      }
    list_pos_of_nodes <- list_pos_upd
    }
  #
  # In iGraph, drawing starts from bottom to top 
  # => reverse nodes order to display from top to bottom
  #
  max_node_pos <- max(list_pos_of_nodes)
  list_pos_of_nodes <- -list_pos_of_nodes + (max_node_pos + 1) 
  #
  # Place contextual and lag0 nodes 
  #
  list_taus <- tmiic.res$tmiic_specific[["tau"]]
  list_taus [which (list_taus < 0)] <- 0
  tau_max <- max (list_taus)
  list_delta_taus <- tmiic.res$tmiic_specific[["delta_tau"]]
  max_lags <- max (list_taus * list_delta_taus)
    
  df_layout <- data.frame ( col=integer(), row=integer() )
  for (i in 1:n_nodes_not_lagged)
    {
    if (is_contextual[[i]])
      {
      if (display == "raw")
        col_display <- max_lags + (max_lags / max(list_taus))
      else
        col_display <- 0
      }
    else
      col_display <- max_lags
    df_layout [i,] <- c(col_display, list_pos_of_nodes[[i]])
    }
  #
  # Place each lagged node using its lag (tau * delta_tau)
  #
  for (tau_idx in 1:tau_max)
    for (node_idx in 1:n_nodes_not_lagged) 
      if (tau_idx <= list_taus[[node_idx]]) 
        {
        col_display <- max_lags - tau_idx * list_delta_taus[[node_idx]]
        df_layout [nrow(df_layout)+1,] <- c (col_display, 
                                             list_pos_of_nodes[[node_idx]])
        }        
  #
  # If laout orientation is portrait
  #
  if (orientation == "P")
    {
    df_layout <- df_layout[,c(2,1)]
    max_pos <- max(df_layout[,1])
    df_layout[,1] <- -df_layout[,1] + (max_pos+1)
    max_pos <- max(df_layout[,2])
    df_layout[,2] <- -df_layout[,2] + (max_pos+1)
    }
  
  layout = as.matrix (df_layout) 
  return (layout)
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
#' @param display [a string]. Optional, default value "compact".
#' Possible values are \emph{"raw"}, \emph{"lagged"}, \emph{"compact"}, 
#' \emph{"combine"}, \emph{"unique"}, \emph{"drop"}:
#' \itemize{
#' \item When \emph{display} = \emph{"raw"}, the plot function will
#'   use the tmiic graph object as it, leading to the display of a lagged
#'   graph. Unless a specific layout is specified, nodes will be positioned 
#'   on a grid. 
#' \item When \emph{display} = \emph{"lagged"}, the function will
#'   repeat the edges over history assuming stationarity and plot a lagged
#'   graph. Unless a specific layout is specified, nodes will be positioned 
#'   on a grid.
#' \item When \emph{display} = \emph{"compact"}, the default, nodes 
#'   and edges are converted into a flattened version to produce a compact 
#'   view of the temporal network whilst still presenting all the information
#'   in the plotting.\cr
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 become respectively X->Y lag=1, 
#'   X<-Y lag=2.  
#' \item When \emph{display} = \emph{"combine"}, prior to the plotting,
#'   a preprocessing will be applied to kept only one edge
#'   per couple of nodes. The info_shifted will be the highest one 
#'   of the summarized edges whilst the lag and orientation of the
#'   summarized edge will be an aggregation.\cr 
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 will become X<->Y lag=1,2 with
#'   the info_shifted of X_lag1->Y_lag0 if info_shifted of 
#'   X_lag1->Y_lag0 > X_lag2<-Y_lag0.
#' \item When \emph{display} = \emph{"unique"}, prior to the plotting,
#'   a preprocessing will be applied to kept only the edges having the
#'   highest info_shifted for a couple of nodes. 
#'   If several edges between the sames nodes have the same
#'   info_shifted, then the edge kept is the one with the minimum lag.\cr 
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 with info_shifted of 
#'   X_lag1->Y_lag0 > X_lag2<-Y_lag0 become X->Y lag=1.
#' \item When \emph{display} = \emph{"drop"}, prior to the plotting,
#'   a preprocessing will be applied to kept only the edges having the 
#'   highest info_shifted for a couple of nodes. 
#'   If several edges between the sames nodes have the same
#'   info_shifted, then the edge kept is the one with the minimum lag.\cr
#'   i.e. :  X_lag1->Y_lag0, X_lag2<-Y_lag0 with info_shifted of 
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
#' @positioning_for_grid [a string] Optional, "greedy" by default.
#' Used only when the display is "raw" or "lagged and no layout is supplied.
#' Possible values are \emph{"none"}, \emph{"alphabetical"}, \emph{"layers"}
#' \emph{"greedy"}, \emph{"sugiyama"}
#' \itemize{
#' \item When \emph{positioning_for_grid} = \emph{"none"}
#'  The nodes are positioned as they appear in the miic result
#' \item When \emph{positioning_for_grid} = \emph{"alphabetical"}
#'  The nodes are positioned alphabeticaly in ascending order
#' \item When \emph{positioning_for_grid} = \emph{"layers"}
#'  The nodes with the less lags wil be placed on the exteriors 
#'  while the nodes having the most lags are in the center
#' \item When \emph{positioning_for_grid} = \emph{"greedy"}
#'  A greedy algorithm will be used to placed the nodes in a way minimizing
#'  the crossing edges
#' \item When \emph{positioning_for_grid} = \emph{"sugiyama"}
#'  The sugiyama algorithm will be used to placed the nodes in a way 
#'  minimizing the crossing edges
#' }
#' 
#' @orientation_for_grid [a string] Optional, "L" by default.
#' Used only when the display is "raw" or "lagged and no layout is supplied.
#' Indicates the orientation of the draw, possible values are landscape: "L"
#' or portrait: "P".
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
#' \donttest{
#' library(miic)
#'
#' #' # EXAMPLE COVID CASES (timeseries demo)
#' data(covidCases)
#' # execute MIIC (reconstruct graph in temporal mode)
#' tmiic.res <- miic(input_data = covidCases, tau = 2, movavg = 14)
#'
#' # to plot the default compact graph
#' if(require(igraph)) {
#'   plot(tmiic.res)
#' }
#'  
#' # to plot the raw temporal network Using igraph
#' if(require(igraph)) {
#'   plot(tmiic.res, display="raw")
#' }
#'
#' # to plot the full temporal network Using igraph
#' if(require(igraph)) {
#'   plot(tmiic.res, display="lagged")
#' }
#'
#' }
#-----------------------------------------------------------------------------
plot.tmiic = function(x, display="compact", show_self_loops=TRUE, 
                      positioning_for_grid="greedy", orientation_for_grid="L", 
                      method = 'igraph', ...) {
  
  if (class(x) != "tmiic")
    stop("Error: Not a tmiic object.")
  if (method != 'igraph')
    stop("Error: Method not supported. See ?tmiic.export for supported methods.")
  if ( !base::requireNamespace("igraph", quietly = TRUE) ) 
    stop("Error: Package 'igraph' is required.")
  if ( is.null (x$adj_matrix) ) 
    stop ("Error: The learnt graphical model adjacency matrix does not exist")
  #
  # Set a layout if none supplied by user : grid like for lagged, 
  # layout_with_kk for flatten displays
  #
  if ( ! ("layout" %in% names(list(...))) ) 
    {
    if (display %in% c("raw", "lagged") )
      layout <- tmiic.computeGridLayout (x, display=display,
              positioning=positioning_for_grid, orientation=orientation_for_grid)
    else 
      layout <- igraph::layout_with_kk
    }
  #
  # Prepare data for the required display
  #
  if ( (display == "lagged") & (x$tmiic_specific[["graph_type"]] == "raw") )
    x <- tmiic.repeat_edges_over_history (x)
  if ( (! display %in% c("raw", "lagged") ) & (x$tmiic_specific[["graph_type"]] == "raw") )
    x <- tmiic.flatten_network(x, flatten_mode=display, 
                               keep_edges_on_same_node=show_self_loops)
  
  x <- tmiic.prepareEdgesForPlotting(x)
  
  df_mult <- data.frame(count=integer(), stringsAsFactors = FALSE)
  if (x$tmiic_specific[["graph_type"]] == "compact")
    df_mult <- tmiic.getMultipleEdgesForPlotting(x)
  #
  # Export the graph to a graphical object and plot
  #
  graph <- tmiic.export (x, display="raw", method=method)
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


