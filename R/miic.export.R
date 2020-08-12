#' Export miic result to different plotting methods
#'
#' @description This function returns an object built from the results
#' of the miic main function that is ready to be fed to different
#' plotting methods.
#'
#' @details See the details of specific function for each method.
#' 
#' @param miic.res [a miic graph object]
#' The graph object returned by the miic execution.
#' @param method A string representing the plotting method. Currently only "igraph" is supported
#' 
#' @export
#' 
#' @return A graph object adapted to the method.
#' 
#' @seealso
#' \code{\link{getIgraph}} for returning an igraph object built from the results.
#'
#' @examples
#' \dontrun{
#' library(miic)
#'
#' # Using igraph
#' library(igraph)
#' g = miic.export(miic.res, "igraph")
#' plot(g) # Default visualisation, calls igraph::plot.igraph()
#'
#' }
#'

miic.export <- function(miic.res, method = NULL) {
  if (is.null(miic.res$all.edges.summary)) {
    stop("The inferred network does not exist")
  }
  if (is.null(method)) {
    stop("Plotting method is required")
  } 
  if (method == "igraph") {
    return(getIgraph(miic.res))
  } else {
    stop("Method not supported")
  }
}


#' Igraph plotting function for miic
#' @description This functions returns an igraph object built from the results
#' of the miic main function.
#'
#' @details
#' Edges attributes are passed to the igraph graph and can be accessed with
#' e.g. \code{E(g)$partial_correlation}. See \code{\link{miic}} for more
#' details on edge parameters. By default, edges are colored according to the
#' partial correlation between two nodes conditioned on the conditioning set
#' (negative is blue, null is gray and positive is red) and their width is
#' based on the conditional mutual information minus the complexity cost.
#'
#' @param miic.res [a miic graph object]
#' The graph object returned by the miic execution.
#' 
#' @export
#'
#' @return An igraph graph object.
#'
#' @seealso
#' \code{\link{miic}} for details on edge parameters in the returned object,
#' \code{\link[igraph]{igraph.plotting}} for the detailed description of the
#' plotting parameters and \code{\link[igraph]{layout}} for different layouts.
#'
#' @examples
#' \dontrun{
#' library(miic)
#' library(igraph)
#'
#' g = getIgraph(miic.res)
#' plot(g) # Default visualisation, calls igraph::plot.igraph()
#'
#' # Specifying layout (see ?igraph::layout_)
#' l <-layout_with_kk(g)
#' plot(g, layout=l)
#'
#' # Override some graphical parameters
#' plot(g, edge.curved = .2)
#' plot(g, vertex.shape="none", edge.color="gray85", vertex.label.color="gray10")
#' }
#'

getIgraph <- function(miic.res) {
  if (is.null(miic.res$all.edges.summary)) {
    stop("The inferred network does not exist")
  }
  if (!base::requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package is required to use this function")
  }

  #### summary
  # plots for all technics on partial correlation
  mySummary <-
    miic.res$all.edges.summary[miic.res$all.edges.summary$type == 'P', ]
  ig_graph = igraph::graph_from_data_frame(mySummary)

  # Set correct orientations
  igraph::E(ig_graph)$arrow.mode = rep(0, igraph::gsize(ig_graph))
  igraph::E(ig_graph)$arrow.mode[igraph::E(ig_graph)$infOrt == 2]  = 2
  igraph::E(ig_graph)$arrow.mode[igraph::E(ig_graph)$infOrt == -2] = 1
  igraph::E(ig_graph)$arrow.mode[igraph::E(ig_graph)$infOrt == 6]  = 3

  # Set visuals
  igraph::V(ig_graph)$color <- "lightblue"
  igraph::V(ig_graph)$label.family <- "Helvetica"
  igraph::V(ig_graph)$label.cex <- 0.8
  igraph::V(ig_graph)$size <- 12

  igraph::E(ig_graph)$width <-
    log10(igraph::E(ig_graph)$info_cond - igraph::E(ig_graph)$cplx)
  
  # Negative pcors are blue, null is dark grey and positive are red
  igraph::E(ig_graph)$color <- "darkgray"
  pcor_palette = grDevices::colorRampPalette(c("blue", "darkgrey", "red"))
  edge_colors_indices = sapply(igraph::E(ig_graph)$partial_correlation,
                               function(pcor) {
                                 ifelse(is.na(pcor), 100, abs(round(pcor * 100)) + 100 * (pcor > 0))
                               })
  igraph::E(ig_graph)$color <-
    pcor_palette(200)[edge_colors_indices]

  return(ig_graph)
}
