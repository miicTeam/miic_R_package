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
#' The graph object returned by the miic execution in temporal mode, 
#' eventually flattened.
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
#' # Plot temporal network Using igraph
#' if(require(igraph)) {
#' g = tmiic.export(tmiic.res, "igraph")
#' plot(g) # Default visualisation, calls igraph::plot.igraph()
#'
#' # Specifying layout (see ?igraph::layout_)
#' l <-layout_with_kk(g)
#' plot(g, layout=l)
#' 
#' # Override some graphical parameters
#' plot(g, edge.arrow.size = 0.75)
#' plot(g, vertex.shape="none", edge.color="gray85", vertex.label.color="gray10")
#' 
#' # For compact graph, please be aware that the rendering of
#' # igraph::plot.igraph() is not optimal when the graph contains
#' # multiple edges between the same nodes.
#' # So the recommend way to plot a compact graph is to use tmiic:  
#' flatten.res <- tmiic.flatten_network(tmiic.res)
#' plot(flatten.res)
#' }
#'
#' }
#-----------------------------------------------------------------------------
tmiic.export <- function (tmiic.res, method = "igraph") 
  {
  if (is.null(tmiic.res$all.edges.summary)) 
    {
    stop("The inferred network does not exist")
    }
  if (is.null(method)) 
    {
    stop("Plotting method is required")
    }
  if (method != "igraph") 
    {
    stop("Method not supported")
    }
  return(tmiic.getIgraph(tmiic.res))
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
#' The graph object returned by the miic execution in temporal mode, 
#' eventually flattened
#'
#' @return An igraph graph object.
#'
#' @seealso
#' \code{\link{miic}} for details on edge parameters in the returned object,
#' \code{\link[igraph]{igraph.plotting}} for the detailed description of the
#' plotting parameters and \code{\link[igraph]{layout}} for different layouts.
#-----------------------------------------------------------------------------
tmiic.getIgraph <- function (tmiic.res) 
  {
  if (class(tmiic.res) != "tmiic")
    {
    stop("Not a tmiic object.")
    }
  
  class(tmiic.res) <- "miic"
  graph <- getIgraph(tmiic.res)
  class(tmiic.res) <- "tmiic"
  
  if ( tmiic.isLaggeg (colnames(tmiic.res$adj_matrix)) )
    {
    igraph::V(graph)$label.dist = 1
    igraph::V(graph)$label.degree = pi/2
    igraph::E(graph)$curved = TRUE
    }
  else
    {
    igraph::E(graph)$label <- tmiic.res$all.edges.summary$lag
    }
  return(graph)
  }
 
#-----------------------------------------------------------------------------
# tmiic.isLaggeg
#-----------------------------------------------------------------------------
#' Test if the list of nodes corresponds to a lagged or non lagged graph
#'
#' @param list_nodes [a list of string]
#'
#' @return A boolean value, TRUE if the network is lagged, FALSE otherwise
#-----------------------------------------------------------------------------
tmiic.isLaggeg <- function (list_nodes) 
  {
  is_graph_lagged = TRUE
  for (one_node in list_nodes)
    {
    regex_found = grepl ("_lag.*$", one_node)
    if (!regex_found)
      {
      is_graph_lagged = FALSE
      break
      }
    }
  return (is_graph_lagged)
  }

#-----------------------------------------------------------------------------
# tmiic.getNbNodesNonLagged
#-----------------------------------------------------------------------------
#' Get the number of non lagged nodes in a list of nodes
#'
#' @description This functions iterates on the list of nodes until if
#' find a node not ending with lag0 and returns the number of iterations done. 
#' 
#' @param list_nodes [a list of string] The expected list of nodes is 
#' a lagged one: all nodes finishing by "lagX"
#'
#' @return A number
#-----------------------------------------------------------------------------
tmiic.getNbNodesNonLagged <- function (list_nodes) 
  {
  n_nodes <- 0
  for (one_node in list_nodes)
    {
    pos_lag_x <- stringr::str_locate(one_node, "_lag")
    lag <- 0
    if ( !is.na(pos_lag_x[1]) )
      {
      lag <- stringr::str_remove(one_node, ".*_lag")
      lag <- strtoi (lag)
      }
    if (lag > 0)
      break
    n_nodes <- n_nodes + 1
    } 
  return (n_nodes)
  } 
  
#-----------------------------------------------------------------------------
# tmiic.prepareEdgesForPlotting
#-----------------------------------------------------------------------------
#' Prepare the edges for plotting
#'
#' @description This function firstly filters the edges in the summary to keep
#' only the ones detected by miic (type = "P") then makes sure that, for each  
#' edge, node name X >= node name Y and adds the couple of nodes of the edge 
#' as id 
#' 
#' @param  [a tmiic graph object] The graph object returned by the miic 
#' execution in temporal mode, eventually flattened
#' 
#' @return tmiic.res [a tmiic object] The modified tmiic object
#-----------------------------------------------------------------------------
tmiic.prepareEdgesForPlotting <- function (tmiic.res) 
  {
  DEBUG <- FALSE
  
  df_edges <- 
    tmiic.res$all.edges.summary[tmiic.res$all.edges.summary$type == 'P', ]
  #
  # Ensure all edges have node x < node y and add an id xy 
  #
  df_edges$xy = NULL
  for (edge_idx in 1:nrow(df_edges)) 
    {
    one_edge <- df_edges[edge_idx,]
    if (one_edge$x < one_edge$y)
      {
      df_edges[edge_idx,"x"] <- one_edge$y
      df_edges[edge_idx,"y"] <- one_edge$x
      if (abs(one_edge$infOrt) == 2) 
        df_edges[edge_idx, "infOrt"] <- -one_edge$infOrt
      if ( !is.na(one_edge$trueOrt) ) 
        if (abs(one_edge$trueOrt) == 2) 
          df_edges[edge_idx, "trueOrt"] <- -one_edge$trueOrt
      }
    df_edges[edge_idx, "xy"] <- paste (df_edges[edge_idx,"x"], "-",
                                       df_edges[edge_idx,"y"], sep="")
    }
  
  tmiic.res$all.edges.summary <- df_edges
  if (DEBUG) 
    {
    print ("df_edges:")
    if ("lag" %in% tmiic.res$all.edges.summary)
      print (tmiic.res$all.edges.summary %>% dplyr::select(xy,x,y,lag,type,infOrt,trueOrt,sign,log_confidence) )
    }
  return (tmiic.res)
  } 

#-----------------------------------------------------------------------------
# tmiic.getMultipleEdgesForPlotting
#-----------------------------------------------------------------------------
#' Look for mutiple edges that needs specific plotting
#'
#' @description This function identifies the couple of nodes having mutiples 
#' edges
#' 
#' @param  [a tmiic graph object]
#' The graph object returned by the miic execution in temporal mode and
#' flattened (if the tmiic object is not flattened, the function does nothing) 
#' 
#' @return df_mult [a dataframe] The dataframe containing the multiple edges
#'
#' @importFrom magrittr "%>%"                             
#-----------------------------------------------------------------------------
tmiic.getMultipleEdgesForPlotting <- function (tmiic.res) 
  {
  DEBUG <- FALSE
  #
  # Find couples of nodes having mutiple edges
  #
  df_mult <- dplyr::group_by (tmiic.res$all.edges.summary, xy, x, y)
  df_mult <- dplyr::summarise (df_mult, count=dplyr::n())
  df_mult <- df_mult[df_mult$count > 1,]
  df_mult <- as.data.frame (df_mult)
  
  if (DEBUG) 
    {
    print ("df_mult:")
    print (df_mult)
    }
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
#' @param tmiic.res [a tmiic graph object]
#' The graph object returned by \code{\link{miic}} in temporal mode, 
#' eventually flatten
#' 
#' @param method A string representing the plotting method. Default to "igraph".
#' Currently only "igraph" is supported.
#' 
#' @param \dots Additional plotting parameters. See the corresponding plot function
#' for the complete list.
#' 
#' For igraph, see \code{\link[igraph]{igraph.plotting}}.
#'
#' @importFrom magrittr "%>%"                             
#' 
#' @export
#'
#' @seealso \code{\link{tmiic.export}} for generic exports,
#' \code{\link{tmiic.getIgraph}} for igraph export,
#' \code{\link[igraph]{igraph.plotting}}
#-----------------------------------------------------------------------------
plot.tmiic = function(tmiic.res, method = 'igraph', ...) 
  {
  DEBUG <- TRUE
  if (class(tmiic.res) != "tmiic")
    {
    stop("Not a tmiic object.")
    }
  if (method != 'igraph')
    {
    stop("Method not supported. See ?tmiic.export for supported methods.")
    }
  if ( !base::requireNamespace("igraph", quietly = TRUE) ) 
    {
    stop("Package 'igraph' is required.")
    }
  if ( is.null (tmiic.res$adj_matrix) ) 
    {
    stop ("The learnt graphical model adjacency matrix does not exist")
    }
  
  tmiic.res <- tmiic.prepareEdgesForPlotting(tmiic.res)
  df_mult <- tmiic.getMultipleEdgesForPlotting(tmiic.res)
  list_nodes <- colnames (tmiic.res$adj_matrix)
  is_graph_lagged = tmiic.isLaggeg (list_nodes)
  #
  # Set a layout if none supplied by user : grid for lagged, 
  # layout_with_kk for flatten
  #
  layout <- NULL
  if ( ! ("layout" %in% names(list(...))) )
    {
    if (is_graph_lagged)
      {
      n_nodes <- tmiic.getNbNodesNonLagged (list_nodes)
      max_lag_plus1 <- length(list_nodes) %/% n_nodes
      layout <- data.frame (rep(1:max_lag_plus1, each=n_nodes),
                            rep(1:n_nodes, times=max_lag_plus1) )
      layout <- as.matrix (layout)
      }
    else
      {
      layout <- igraph::layout_with_kk
      }
    }
  #
  # Export the graph to a graphical objet and plot
  #
  graph <- tmiic.export (tmiic.res, method)
  df_edges <- tmiic.res$all.edges.summary
  if (nrow (df_mult) <= 0)
    {
    # No multiple edges between the same nodes, we draw in one go
    #
    if ( is.null(layout) )
      igraph::plot.igraph (graph, ...)
    else
      igraph::plot.igraph (graph, layout=layout, ...)
    }
  else
    {
    # Multiple edges between the same nodes exist, draw iteratively
    #
    edges_colors_iter <- igraph::E(graph)$color
    edges_labels_iter <- igraph::E(graph)$label
    #
    # On a first step, we will draw all the graph except multiple edges
    # The multiple edges will be drawn with invisible color "#FF000000"
    # and with no labels 
    #
    for ( edge_idx in 1:nrow(df_edges) )
      {
      one_edge <- df_edges[edge_idx,]
      if (one_edge$xy %in% df_mult$xy)
        {
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
    for ( mult_idx in 1:nrow(df_mult) )
      {
      one_mult <- df_mult[mult_idx,]
      if (one_mult$x == one_mult$y)
        {
        # for self loop, we will go over 2*pi around the node
        #
        step_pos <- 0
        step_inc <- (2 * pi) / one_mult$count
        }
      else
        {
        # otherelse, we will curve edges from -0.5 to +0.5
        #
        if (one_mult$count > 4) # if more than 4 edges, curve more
          {
          step_pos <- -1
          step_inc <- 2.0 / (one_mult$count - 1)
          }
        else
          {
          step_pos <- -0.5
          step_inc <- 1.0 / (one_mult$count - 1)
          }
        }
      #
      # Draw mutliple egdes one by one
      #
      list_to_draw = which(df_edges[, "xy"] == one_mult$xy)
      for (idx_to_draw in 1:length(list_to_draw) )
        {
        edge_to_draw = list_to_draw[[idx_to_draw]]
        #
        # We hide all edges except one
        #
        edges_colors_iter <- rep ("#FF000000", nrow (df_edges) )
        edges_labels_iter <- rep (NA, nrow (df_edges) )
        edges_colors_iter[[edge_to_draw]] <- igraph::E(graph)[[edge_to_draw]]$color
        edges_labels_iter[[edge_to_draw]] <- igraph::E(graph)[[edge_to_draw]]$label

        if (one_mult$x == one_mult$y)
          igraph::plot.igraph (graph, layout=layout, add=TRUE,
                               edge.color=edges_colors_iter, 
                               edge.label=edges_labels_iter, 
                               edge.loop.angle=step_pos, ...)
        else
          igraph::plot.igraph (graph, layout=layout, add=TRUE,
                               edge.color=edges_colors_iter, 
                               edge.label=edges_labels_iter, 
                               edge.curved=step_pos, ...)
        #
        # Update position for next edge
        #
        step_pos <- step_pos + step_inc
        }
      }
    }
  }

  

