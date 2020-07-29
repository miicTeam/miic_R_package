#*****************************************************************************
# Filename   : miic.plot.R                     Creation date: 11 june 2020
#
# Description: iGraph plotting for miic
#
# Author     : Franck SIMON (fsimon.informaticien@wanadoo.fr)
#
# Changes history:
# - 11 june 2020 : initial version
#   This version is a rewrite from miic.plot.R and gmPlot.lib.R. 
#   miic.plot is now able to plot multiple edges between the same nodes
#
# TODO:
# - Transfer each of the 3 parts (data preparation, graph, plot) into a sub
#   function
#*****************************************************************************

#-----------------------------------------------------------------------------
# miic.plot
#-----------------------------------------------------------------------------
#' Igraph plotting function for miic
#' 
#' @description This functions plots the network with the given layout 
#' (if specified) using the igraph package.
#'
#' @details The plot reports the partial correlation or the log_confidence
#' as strength of the edges.
#'
#' @param g [a miic graph object]
#' The graph object returned by the miic execution.
#' @param method [a string; \emph{c("pcor", "log_confidence")}]
#' Optional, "log_confidence" by default. The column used to plot the 
#' strength of the edges. 
#' @param igraph_layout [an igraph layout object]
#' Optional, \emph{layout_with_kk} by default. When set, it is used to plot 
#' the network. See the igraph manual for more information.
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
#' @param file_figsize [a list: \emph{c(length, height)}] Optional, NULL by 
#' default. When plots are drawn to a file, the size of the draw in pixels 
#' can be specified by a couple of values. The first is the lentgh, the 
#' second the height.
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
#' @param edges_labels [a vector] Optional, NULL by default.
#' The labels to display on the edges. 
#' @param edges_label_cex [a number] Optional, NULL by default.
#' The size of the edges labels. 
#' @param edges_width [a number] Optional, NULL by default.
#' The width of the edges.  
#' @param edges_curved [a boolean or a vector] Optional, NULL by default.
#' The curvatures to apply to the edges. When TRUE, edges are curved 
#' and when FALSE, edges are straitgh.\cr
#' Note that mutiple edges between the same nodes will always be curved
#' whatever value has \emph{edges_curved}. 
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
#-----------------------------------------------------------------------------
miic.plot <- function (g, method = "log_confidence", 
                       igraph_layout = NULL, user_layout = NULL, 
                       miic_defaults=TRUE, filename=NULL, file_figsize=NULL,
                       font_family=NULL, title=NULL, title_cex=1.5, 
                       draw_legend=TRUE, graph_margin=NULL, 
                       nodes_label_dist=NULL, nodes_label_degree=NULL, 
                       nodes_label_cex=NULL, nodes_shape=NULL, 
                       nodes_colors=NULL, nodes_sizes=NULL, 
                       edges_labels=NULL, edges_label_cex=NULL, 
                       edges_width=NULL, edges_curved=NULL, 
                       edges_arrow_size=NULL, edges_arrow_width=NULL,
                       verbose = FALSE) 
  {
  DEBUG <- FALSE
  if ( (verbose) | (DEBUG) )
      cat ("# --------\n# -> START Plot...\n")
  #
  # Check inputs
  #
  if ( (method != "pcor") & (method != "log_confidence") )
    stop ("incorrect method argument: must be \"log_confidence\" or \"pcor\" ")
  if ( is.null (g$adjMatrix) ) 
    stop ("The learnt graphical model adjacency matrix does not exist")
  if ( is.null (g$all.edges.summary) ) 
    stop ("The learnt graphical model summary does not exist")
  if (method == "pcor") 
    {
    if (is.na (g$all.edges.summary$partial_correlation[1])) 
      {
      cat ("Impossible to plot correlation without data in summary\n")
      cat ("# -> END Plot...\n# --------\n")
      return()
      } 
    } 
  else
    {
    if (is.na (g$all.edges.summary$log_confidence[1]) ) 
      {
      cat ("Impossible to plot log_confidence without data in summary\n")
      cat ("# -> END Plot...\n# --------\n")
      return()
      } 
    } 
  #
  ############################################################################
  # PREPARATION OF DATA
  ############################################################################
  #
  # Extract nodes and edges. 
  #
  list_nodes <- colnames (g$adjMatrix)
  if (DEBUG) 
    {
    print ("list_nodes:")
    print (list_nodes)
    }
  #
  # Edge are filtered to remove "TN", "N", "FN"
  #
  df_edges <- g$all.edges.summary
  cond_filter <- ( (df_edges[["type"]] != "TN") 
                 & (df_edges[["type"]] != "N") 
                 & (df_edges[["type"]] != "FN") ) 
  df_edges <- df_edges[cond_filter,]
  #
  # As we have a bug with curving of multiple edges between same nodes
  # We will apply a specific treatment to deal with these edges. 
  # At first, we ensure X string < Y string and we create an extra
  # column containing x-y so we have an id to perform a duplicate check
  #
  df_edges$xy = NULL
  for (edge_idx in 1:nrow(df_edges)) 
    {
    one_edge <- df_edges[edge_idx,]
    if (one_edge$x < one_edge$y)
      {
      df_edges[edge_idx,"x"] <- one_edge$y
      df_edges[edge_idx,"y"] <- one_edge$x
      # inverse edge oriention if oriented and not bidirectional
      # not sure -4,4 orientations still exists, it was present in the old code
      if ( (abs(one_edge$infOrt) == 2) | (abs(one_edge$infOrt) == 4) )
        df_edges[edge_idx, "infOrt"] <- -one_edge$infOrt
      if ( !is.na(one_edge$trueOrt) ) 
        if ( (abs(one_edge$trueOrt) == 2) | (abs(one_edge$trueOrt) == 4) )
          df_edges[edge_idx, "trueOrt"] <- -one_edge$trueOrt
      }
    df_edges[edge_idx, "xy"] <- paste (df_edges[edge_idx,"x"], "-",
                                       df_edges[edge_idx,"y"], sep="")
    }
  if (DEBUG) 
    {
    print ("df_edges:")
    if ( "lag" %in% colnames (df_edges) )
      print (df_edges %>% select(xy,x,y,lag,type,infOrt,trueOrt,sign,log_confidence) )
    else
      print (df_edges %>% select(xy,x,y,type,infOrt,trueOrt,sign,log_confidence) )
    }
  #
  # Handle the layout : firstly, take user layout if supplied, 
  # otherwise, use specific iGraph layout if supplied
  # otherwise, use igraph::layout_with_kk
  #
  if ( !is.null (user_layout) ) 
    {
    layout <- user_layout
    if (ncol (layout) > 2) 
      # keep only posX and posY
      layout <- layout[, 2:3]
    layout <- as.matrix (layout)
    if ( (verbose) | (DEBUG) )
      cat ("Load the positions of the vertices\n")
    }
  else 
    {
    if ( is.null (igraph_layout) ) 
      {
      layout <- igraph::layout_with_kk
      if ( (verbose) | (DEBUG) )
        cat ("Use igraph::layout_with_kk layout as none supplied\n")
      }
    else 
      {
      layout <- igraph_layout
      if ( (verbose) | (DEBUG) )
        cat ("Use supplied igraph::layout\n")
      }
    }
  #
  # Define the color gradients
  #
  blue.gradient <- grDevices::rainbow (100, start = 3 / 6, end = 4 / 6)
  red.gradient <- rev (grDevices::rainbow (100, start = 0, end = 0.16) )
  #
  # Prepare the orientation and color of each edge in two vectors
  #
  edges_orients <- rep (NA, nrow (df_edges) )
  edges_colors <- rep (NA, nrow (df_edges) )
  for (edge_idx in 1:nrow(df_edges)) 
    {
    # Orientation part
    #
    one_edge_orient <- df_edges[edge_idx, "infOrt"]
    if ( (one_edge_orient == 2) | (one_edge_orient == 4) ) # forward oriented
      edges_orients[edge_idx] <- 2
    else if ( (one_edge_orient == -2) | (one_edge_orient == -4) ) # backward oriented
      edges_orients[edge_idx] <- 1
    else if (one_edge_orient == 6) # bidirectional
      edges_orients[edge_idx] <- 3
    else  # not oriented
      edges_orients[edge_idx] <- 0
    #
    # Color part depending on the method
    #
    edge_color_idx <- NULL
    if (method == "pcor") 
      {
      edge_color_idx <- df_edges[edge_idx, "partial_correlation"]
      edge_color_idx <- round ( abs (edge_color_idx) * 100 )
      if (edge_color_idx == 0) 
        edge_color_idx <- 1
      }
    else # log_confidence
      {
      edge_color_idx <- df_edges[edge_idx, "log_confidence"]
      if (edge_color_idx < 1) # min confidence
        edge_color_idx <- 1
      if (edge_color_idx > 100) # max confidence
        edge_color_idx <- 100
      }
    #
    # Get the sign of the link to look at the correct color gradient
    #
    if ( is.na (df_edges[edge_idx, "sign"]) )
      {
      edges_colors[edge_idx] <- "grey88"
      next
      }
    edge_sign <- df_edges[edge_idx, "sign"]
    if (edge_sign == "+") 
      edges_colors[edge_idx] <- red.gradient [edge_color_idx]
    else
      edges_colors[edge_idx] <- blue.gradient [edge_color_idx]
    }
  #
  # To later identify multiple edges between same nodes 
  #
  df_mult <- dplyr::group_by (df_edges, xy, x, y)
  df_mult <- dplyr::summarise (df_mult, count=dplyr::n())
  if ( is.null (graph_margin) )
    {
    # Define a margin (if none supplied) on left and right sides
    # of the graph when we have self loop(s) to avoid self loops 
    # to be drawn outside of the plotting area
    #
    graph_margin <- c(0,0,0,0)
    if (sum (df_mult[["x"]] == df_mult[["y"]]) > 0)
      graph_margin <- c(0,0.25,0,0.25)
    }
  df_mult <- df_mult[df_mult$count > 1,]
  df_mult <- as.data.frame (df_mult)
  #
  ############################################################################
  # GRAPH PART
  ############################################################################
  #
  # Edges colors and orientations are ready, we can create our graph
  #
  graph <- igraph::graph.data.frame (df_edges[,c("x","y")], vertices=list_nodes, directed=TRUE)
  #
  # Set vertices options
  #
  if (miic_defaults)
    {
    igraph::V(graph)$label.family <- "Helvetica"
    igraph::V(graph)$color <- "lightblue"
    igraph::V(graph)$size <- 10
    igraph::V(graph)$label.cex <- 0.6
    }
  
  if ( !is.null (nodes_shape) )
    igraph::V(graph)$shape <- nodes_shape
  if ( !is.null (font_family) )
    igraph::V(graph)$label.family <- font_family
  if ( !is.null (nodes_colors) )
    igraph::V(graph)$color <- nodes_colors
  if ( !is.null (nodes_sizes) )
    igraph::V(graph)$size <- nodes_sizes
  if ( !is.null (nodes_label_dist) )
    igraph::V(graph)$label.dist <- nodes_label_dist
  if ( !is.null (nodes_label_degree) )
    igraph::V(graph)$label.degree <- nodes_label_degree
  if ( !is.null (nodes_label_cex) )
    igraph::V(graph)$label.cex <- nodes_label_cex

  igraph::V(graph)$label <- list_nodes
  #
  # Set the edges options
  #
  if (miic_defaults)
    {
    igraph::E(graph)$arrow.size <- 0.2
    igraph::E(graph)$arrow.width <- 3
    igraph::E(graph)$width <- 3
    igraph::E(graph)$curved <- FALSE
    }
  
  if ( !is.null (edges_labels) )
    igraph::E(graph)$label <- edges_labels
  # Generate an error : graphical parameter "family" has the wrong length 
  # (iGraph github issue #37)
  # if ( !is.null (font_family) )
  #   igraph::E(graph)$label.family <- font_family
  if ( !is.null (edges_label_cex) )
    igraph::E(graph)$label.cex <- edges_label_cex
  if ( !is.null (edges_width) )
    igraph::E(graph)$width <- edges_width
  if ( !is.null(edges_curved) )
    igraph::E(graph)$curved <- edges_curved
  if ( !is.null (edges_arrow_size) )
    igraph::E(graph)$arrow.size <- edges_arrow_size
  if ( !is.null (edges_arrow_width) )
    igraph::E(graph)$arrow.width <- edges_arrow_width

  igraph::E(graph)$color <- edges_colors
  igraph::E(graph)$arrow.mode <- edges_orients
  #
  ############################################################################
  # PLOT PART
  ############################################################################
  #
  # Save the graphic settings (the ones that can be modified)
  #
  sav_config <- par (no.readonly=TRUE) 
  if ( !is.null(font_family) )
    par (family=font_family)
  #
  # If a filename is supplied, redirect the plot output
  #
  if (! is.null(filename) )
    {
    if (! is.null(file_figsize) )    
      png (filename=filename, width=file_figsize[[1]], 
           height=file_figsize[[2]], units="px")
    else
      png (filename=filename)
    }
  #
  # If we draw the legend, divide the layout for plot and legend
  # and set margins 
  #
  if (draw_legend)
    {
    margins <- par("mar")
    if (method == "pcor")
      {
      graphics::layout (t(1:2), widths=c(5, 1))
      graphics::par (mar=c(.5, .5, .5, .5),
                     oma=c(0.5, .5, margins[[2]], .5) ) 
      }
    else # log_confidence
      {
      graphics::layout (t(1:3), widths = c(5, 1, 1))
      graphics::par (mar=c(0, .5, .5, .5),
                     oma=c(0, .5, (margins[[2]]+1), .5) ) 
      }
    }
  #
  # Plotting of the network
  #
  if (nrow (df_mult) <= 0)
    {
    # No multiple edges between the same nodes, we draw in one go
    #
    graphics::plot (graph, layout=layout, margin=graph_margin)
    }
  else
    {
    # Multiple edges between the same nodes exist, draw iteratively
    #
    edges_colors_iter <- edges_colors
    edges_labels_iter <- edges_labels
    #
    # On a first step, we will draw all the graph except multiple edges
    # The multiple edges will be drawn with invisible color "#FF000000"
    # and with no labels (if any)
    #
    for ( edge_idx in 1:nrow(df_edges) )
      {
      one_edge <- df_edges[edge_idx,]
      if (one_edge$xy %in% df_mult$xy)
        {
        edges_colors_iter[[edge_idx]] <- "#FF000000"
        if (! is.null (edges_labels_iter) )
          edges_labels_iter[[edge_idx]] <- NA
        }
      }
    graphics::plot (graph, layout=layout, margin=graph_margin, 
                    edge.color=edges_colors_iter, edge.label=edges_labels_iter)
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
        edges_colors_iter[[edge_to_draw]] <- edges_colors[[edge_to_draw]]
        edges_labels_iter[[edge_to_draw]] <- edges_labels[[edge_to_draw]]
        
        if (one_mult$x == one_mult$y)
          graphics::plot (graph, layout=layout, margin=graph_margin, add=TRUE, 
                          edge.color=edges_colors_iter, edge.label=edges_labels_iter,
                          edge.loop.angle=step_pos)
        else
          graphics::plot (graph, layout=layout, margin=graph_margin, add=TRUE, 
                          edge.color=edges_colors_iter, edge.label=edges_labels_iter,
                          edge.curved=step_pos)
        #
        # Update position for next edge
        #
        step_pos <- step_pos + step_inc
        }
      }
    }
  #
  # Finition of the plot, add title
  #
  if (! is.null(title) )
    {
    if (draw_legend)
      graphics::mtext (title, cex=title_cex, side=3, outer=TRUE, line=-1)
    else  
      graphics::title (title, cex=title_cex)
    }
  #
  # Add legend
  #
  if (draw_legend)
    {
    if (method == "pcor")
      {
      leg_colors <- c(red.gradient, blue.gradient)
      legend_image <- grDevices::as.raster (matrix (leg_colors, ncol=1) )
      graphics::plot (x=c(0, 5), y=c(0, 0.8), type="n", axes=F, xlab="", ylab="")
      # graphics::par (adj=0.5)
      graphics::title ("Partial\ncorrelation", cex.main=1, line=-2)
      graphics::text (x=1.5, y=seq(0, 0.75, l=5), 
                      labels=seq(-1, 1, l=5), cex=1)
      graphics::rasterImage(legend_image, 2.5, 0, 3.5, 0.75)
      }
    else # log_confidence
      {
      # Positive correlations
      # 
      legend_image <- grDevices::as.raster (matrix (red.gradient, ncol=1))
      graphics::plot (x=c(0, 4), y=c(0.2, 0.8), type="n", axes=F, xlab="", ylab="")
      graphics::par (adj=1)
      graphics::title ("Confidence\npcor+", cex.main=1, line=-4)
      graphics::rasterImage (legend_image, 3.3, 0.25, 3.8, 0.75)
      #
      # Negative correlations
      # 
      legend_image <- grDevices::as.raster (matrix (rev(blue.gradient), ncol=1) )
      graphics::plot (x=c(0, 4), y=c(0.2, 0.8), type="n", axes=F, xlab="", ylab="")
      graphics::par (adj=0.5)
      graphics::title ("(NI' = -log Pxy)\npcor-", cex.main=1, line=-4)
      graphics::rasterImage(legend_image, 1.5, 0.25, 2, 0.75)
      
      graphics::text (x=rep(.5, 3), y=c(0.26, 0.5, 0.74), 
                      labels=c("< 1", "50", "> 100"), cex=1)
      }
    }
  
  if (! is.null(filename) )
    dev.off()
  par (sav_config, new=FALSE)
  if ( (verbose) | (DEBUG) )
      cat ("# -> END Plot...\n# --------\n")
  }

