#*****************************************************************************
# Filename   : tmiic.plot.R                     Creation date: 24 march 2020
#
# Description: Function to plot lagged graphs of temporal miic (tmiic)
#
# Author     : Franck SIMON (fsimon.informaticien@wanadoo.fr)
#
# Changes history:
# - 24 march 2020 : initial version
# - 11 may 2020   : modification to handle condensed networks
#*****************************************************************************

#-----------------------------------------------------------------------------
# tmiic.plot
#-----------------------------------------------------------------------------
#' Igraph plotting function for tmiic (temporal miic)
#' 
#' @description This function plots the temporal network with the given layout
#' (if specified) using the igraph package.
#'
#' @details The plot reports the partial correlation or the log_confidence 
#' as strength of the edges. The function process differently condensed and
#' lagged networks.
#' 
#' For lagged networks (the ones with nodes including '_lagX' at their names 
#' ends), the function will define an user layout if none is supplied 
#' and use curved edges.
#' 
#' For condensed networks (the network obtained after removing the '_lagX' 
#' at the end of the nodes names), the function will plot the lag 
#' information on the edges.
#' 
#' @param g [a miic graph object]
#' The graph object returned by the miic execution.
#' @param method [a string; \emph{c("pcor", "log_confidence")}]
#' Optional, "log_confidence" by default. The column used to plot the 
#' strength of the edges. 
#' @param igraph_layout [an igraph layout object]
#' Optional, \emph{layout_with_kk} by default. When set, it is used to plot the network. 
#' See the igraph manual for more information.
#' @param user_layout [a data frame]
#' Optional, NULL by default. A data frame reporting the position of nodes. 
#' Each line corresponds to the \emph{(x,y)} coordinates of each vertex.
#' When the data frame has two columns, the first one is assocated with
#' \emph{x} and the second to \emph{y}. When the data frame has more than 
#' two columns, columns two and three are extracted and mapped, 
#' respectively, to \emph{x} and \emph{y} positions.
#' @param miic_defaults [a boolean] Optional, TRUE by default. 
#' When TRUE, several graphic parameters are intialised with default values
#' suited for most of miic results:
#' \itemize{
#' \item For nodes, font family: "Helvetica", color: "lightblue", size: 10, 
#' label.cex: 0.6
#' \item For edges, arrow.size: 0.2, arrow.width: 3, width: 3, curved: FALSE
#' }
#' Note that if the user supplies specific values for these parameters
#' (ie a specific node size), the user choices will overwrite the miic 
#' defaults.
#' @param filename [a string] Optional, NULL by default. 
#' If supplied, plot will be saved here.
#' @param file_figsize [a list; \emph{c(length, height)}] Optional, NULL by default. 
#' When plots are drawn to a file, the size of the draw in pixels can be specified
#' by a couple of values. The first is the lentgh, the second the height.
#' @param font_family [a string] Optional, NULL by default. The font to use.\cr
#' Note that the font may not apply in all the graphical parts, depending
#' on how the plotted objects manage the font parameter.
#' @param title [a string] Optional, NULL by default. The title of the plot. 
#' @param title_cex [a number] Optional, 1.5 by default. The size of the title. 
#' @param draw_legend [a boolean] Optional, TRUE by default. 
#' When TRUE, a legend is drawn on the right side of the plot.
#' @param graph_margin [a list of floats: c(bottom,left,top,right)]
#' Optional, NULL by default. Margin applied around the graph. 
#' When NULL, the value is determined by the presence or absence of 
#' self loop(s) in the graph: 
#' \itemize{
#' \item no self loop = no margin 
#' \item at least one self loop = margin of 0.25 on the left and right sides 
#' }
#' @param nodes_label_dist [a number] Optional, NULL by default. 
#' The distance of the labels from the nodes. 
#' @param nodes_label_degree [a number] Optional, NULL by default.
#' The degre where to draw the labels from the nodes. 
#' @param nodes_label_cex [a number] Optional, NULL by default.
#' The size of the nodes labels. 
#' @param nodes_shape [a string] Optional, NULL by default.
#' The shape of the nodes. See the igraph manual for more information.
#' @param nodes_colors [a string or a vector] Optional, NULL by default.
#' The colors of the nodes. See the igraph manual for more information.
#' @param nodes_sizes [a number or a vector] Optional, NULL by default.
#' The sizes of the nodes. See the igraph manual for more information.
#' @param edges_label_cex [a number] Optional, NULL by default.
#' The size of the edges labels. 
#' @param edges_width [a number] Optional, NULL by default.
#' The width of the edges.  
#' @param edges_curved [a boolean or a vector] Optional, NULL by default.
#' The curvatures to apply to the edges:\cr
#' \itemize{
#' \item when TRUE, edges are curved\cr
#' \item when FALSE, edges are straitgh\cr
#' \item when NULL, the function will choose the best curvature depending 
#' on the type of network: curved for a lagged graphs and straitgh for 
#' condensed graphs.
#' }
#' Note that mutiple edges between the same nodes will always be curved
#' whatever value has edges_curved. 
#' See the igraph manual for more information. 
#' @param edges_arrow_size [a number] Optional, NULL by default.
#' The size of the edges arrows.  
#' @param edges_arrow_width [a number] Optional, NULL by default.
#' The width of the edges arrows.  
#' @param verbose [a boolean value] Optional, FALSE by default. If TRUE, 
#' debugging output is printed.
#' 
#' @export
#' @useDynLib miic
#--------------------------------------------------------------------------------
tmiic.plot <- function (g, method="log_confidence", igraph_layout=NULL, 
                        user_layout=NULL, miic_defaults=TRUE, 
                        filename=NULL, file_figsize=NULL, 
                        font_family=NULL, title=NULL, title_cex=1.5, 
                        draw_legend=TRUE, graph_margin=NULL, 
                        nodes_label_dist=NULL, nodes_label_degree=NULL, 
                        nodes_label_cex=NULL, nodes_shape=NULL, 
                        nodes_colors=NULL, nodes_sizes=NULL, 
                        edges_label_cex=NULL, edges_width=NULL, edges_curved=NULL, 
                        edges_arrow_size=NULL, edges_arrow_width=NULL,
                        verbose=FALSE) 
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
    print (paste ("edges_curved=", edges_curved, sep="") )
    print ("input adjacency matrix:")
    print (g$adjMatrix)
    print ("input summary:")
    print (g$all.edges.summary [ g$all.edges.summary[["type"]] == "P" ])
    }
  #
  # Check if the network is the lagged or the flatten one:
  # If any node does not end with "_lag*", the network is not lagged
  #
  is_graph_lagged = TRUE
  list_nodes <- colnames (g$adjMatrix)
  for (one_node in list_nodes)
    {
    regex_found = grepl ("_lag.*$", one_node)
    if (!regex_found)
      {
      is_graph_lagged = FALSE
      break
      }
    }
  #
  # Set some options depending of the type of graph (condensed or lagged)
  #
  edges_labels = NULL
  if (! is_graph_lagged)
    {
    # If the graph is condensed, we display lag values on edges
    #
    edges_labels = g$all.edges.summary$lag
    }
  else 
    { 
    # For lagged graph, Set some params if not supplied
    #
    if (is.null(nodes_label_dist) )
      nodes_label_dist = 1
    if (is.null(nodes_label_degree) )
      nodes_label_degree = pi/2
    if (is.null(edges_curved) )
      edges_curved = TRUE
    #
    # Set a layout if none is supplied
    #
    if ( is.null(user_layout) )
      {
      #
      # Get the number of non lagged nodes (the first not lag0)
      #
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
      max_lag_plus1 <- length(list_nodes) %/% n_nodes
      user_layout = data.frame (list_nodes, rep(1:max_lag_plus1, each=n_nodes), 
                                           rep(1:n_nodes, times=max_lag_plus1) )    
  
      if (DEBUG)
        {
        print ("no user layout given:")
        print ("found these nodes=")
        print (list_nodes)
        print (paste ("found max lag=", max_lag_plus1 - 1, sep="") )
        print (paste ("found n_nodes=", n_nodes, sep="") )
        print ("layout defined=")
        print (user_layout)
        }
      }
    } 
  if (DEBUG)
    {
    print (paste ("is_graph_lagged =", is_graph_lagged) )
    print (paste ("edges_curved=", edges_curved, sep="") )
    }
  #
  # Plot the graph
  #
  miic.plot (g, method=method, igraph_layout=igraph_layout, 
             user_layout=user_layout, miic_defaults=miic_defaults, 
             filename=filename, file_figsize=file_figsize, 
             font_family=font_family, title=title, title_cex=title_cex,
             draw_legend=draw_legend, graph_margin=graph_margin,
             nodes_label_dist=nodes_label_dist, nodes_label_degree=nodes_label_degree, 
             nodes_label_cex=nodes_label_cex, nodes_shape=nodes_shape, 
             nodes_colors=nodes_colors, nodes_sizes=nodes_sizes, 
             edges_labels=edges_labels, edges_label_cex=edges_label_cex, 
             edges_width=edges_width, edges_curved=edges_curved, 
             edges_arrow_size=edges_arrow_size, edges_arrow_width=edges_arrow_width,
             verbose=verbose) 
  }
