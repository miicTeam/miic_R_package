#-------------------------------------------------------------------------------
# export
#-------------------------------------------------------------------------------
#' Export miic result for plotting (with igraph)
#'
#' @description This function creates an object built from the result returned
#' by \code{\link{miic}} that is ready to be fed to the plotting method.
#'
#' @details The behavior depends on the method used for the export.
#'
#' For igraph, edges attributes are passed to the igraph graph
#' and can be accessed with e.g. \code{E(g)$partial_correlation}.
#' See \code{\link{miic}} for more details on edge parameters.
#' By default, edges are colored according to the partial correlation
#' between two nodes conditioned on the conditioning set
#' (negative is blue, null is gray and positive is red)
#' and their width is based on the conditional mutual information
#' minus the complexity cost.
#'
#' @param mo [a miic object, required]
#'
#' The object returned by the \code{\link{miic}} execution.
#'
#' @param method  [a string, optional, default value "igraph"]
#'
#' The plotting method, currently only "igraph" is supported.
#'
#' @param pcor_palette [a color palette, optional, default value
#' grDevices::colorRampPalette(c("blue", "darkgrey", "red")]
#'
#' Used to represent the partial correlations (the color of the edges).
#' The palette must be able to handle 201 shades to cover the correlation range
#' from -100 to +100.
#'
#' @param display [a string, optional, default value "compact"]
#'
#' Used only when exporting object returned by miic in temporal mode.
#' It allows different representations of the temporal graph.
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
#'   e.g. X_lag1->Y_lag0, X_lag2<-Y_lag0 become respectively X->Y lag=1,
#'   X<-Y lag=2.
#' \item When \emph{display} = \emph{"combine"}, prior to the export,
#'   a pre-processing will be applied to kept only one edge
#'   per couple of nodes. The info_shifted will be the highest one
#'   of the summarized edges whilst the lag and orientation of the
#'   summarized edge will be an aggregation.\cr
#'   e.g. X_lag2->Y_lag0, X_lag0<-Y_lag1 will become X<->Y lag=1-2 with
#'   the info_shifted of X_lag2->Y_lag0 if info_shifted of
#'   X_lag2->Y_lag0 > X_lag0<-Y_lag1.
#' \item When \emph{display} = \emph{"unique"}, prior to the export,
#'   a pre-processing will be applied to kept only the edges having the
#'   highest info_shifted for a couple of nodes.
#'   If several edges between the sames nodes have the same
#'   info_shifted, then the edge kept is the one with the minimum lag.\cr
#'   e.g. X_lag1->Y_lag0, X_lag0<-Y_lag2 with info_shifted of
#'   X_lag1->Y_lag0 > X_lag0<-Y_lag2 become X->Y lag=1.
#' \item When \emph{display} = \emph{"drop"}, prior to the export,
#'   a pre-processing will be applied to kept only the edges having the
#'   highest info_shifted for a couple of nodes.
#'   If several edges between the sames nodes have the same
#'   info_shifted, then the edge kept is the one with the minimum lag.\cr
#'   e.g. X_lag1->Y_lag0, X_lag0<-Y_lag2 with info_shifted of
#'   X_lag1->Y_lag0 > X_lag0<-Y_lag2 become X->Y.
#'   The lag information is dropped during the preprocessing and
#'   will not be exported.
#' }
#'
#' @param show_self_loops [a boolean, optional, TRUE by default]
#'
#' Used only when exporting object returned by miic in temporal mode.
#' When TRUE, the lagged edges starting and ending on the same node
#' are included in the igraph  object.
#' When FALSE, only edges having different nodes are present in the igraph
#' object.
#'
#' @export
#'
#' @return A graph object adapted to the method.
#'
#' @examples
#' \donttest{
#' library(miic)
#' data(hematoData)
#'
#' # execute MIIC (reconstruct graph)
#' miic_obj <- miic(
#'   input_data = hematoData, latent = "yes",
#'   n_shuffles = 10, conf_threshold = 0.001
#' )
#'
#' # Using igraph
#' if(require(igraph)) {
#' g = export(miic_obj, "igraph")
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
#' # In temporal mode, execute MIIC
#' data(covidCases)
#' tmiic_obj <- miic(input_data = covidCases, mode = "TS", n_layers = 3, delta_t = 1, movavg = 14)
#'
#' # Plot by default the compact display of the temporal network using igraph
#' if(require(igraph)) {
#' g = export (tmiic_obj)
#' plot(g)
#'
#' # Plot the raw temporal network using igraph
#' g = export(tmiic_obj, display="raw")
#' plot(g)
#'
#' # Plot the complete temporal network using igraph (completed by stationarity)
#' g = export(tmiic_obj, display="lagged")
#' plot(g)
#'
#' # Specifying layout (see ?igraph::layout_)
#' l <- layout_on_grid(g, width = 5, height = 3, dim = 2)
#' plot(g, layout=l)
#'
#' # For compact temporal display, please be aware that the rendering of
#' # igraph::plot.igraph() is not optimal when the graph contains
#' # multiple edges between the same nodes.
#' # So, the recommend way to plot a compact graph is to use tmiic plotting:
#' plot(tmiic_obj)
#' }
#'
#' }
#-------------------------------------------------------------------------------
export <- function (mo, method="igraph", pcor_palette=NULL,
                    display="compact", show_self_loops=TRUE)
  {
  if ( is.null(mo$summary) )
    stop("The inferred network does not exist")
  if ( (!is.null(method)) && (method != "igraph") )
    stop("Method not supported")

  if ( is.null(mo$tmiic) )
    return (getIgraph(mo, pcor_palette=pcor_palette))
  else
    return (tmiic_getIgraph (mo, pcor_palette=pcor_palette,
                             display=display, show_self_loops=show_self_loops))
  }

#-------------------------------------------------------------------------------
#' Igraph export function for miic
#'
#' @description This functions returns an igraph object built from the result
#' returned by \code{\link{miic}}.
#'
#' @details
#' Edges attributes are passed to the igraph graph and can be accessed with
#' e.g. \code{E(g)$partial_correlation}. See \code{\link{miic}} for more
#' details on edge parameters. By default, edges are colored according to the
#' partial correlation between two nodes conditioned on the conditioning set
#' (negative is blue, null is gray and positive is red) and their width is
#' based on the conditional mutual information minus the complexity cost.
#'
#' @param mo [a miic object]
#' The object returned by the \code{\link{miic}} execution.
#' @param pcor_palette The color palette used to represent the partial correlations
#' (the color of the edges). The palette must be able to handle 201 shades
#' to cover the correlation range from -100 to +100. The default palette is
#' grDevices::colorRampPalette(c("blue", "darkgrey", "red").
#'
#' @return An igraph graph object.
#'
#' @noRd
#'
#' @seealso
#' \code{\link{miic}} for details on edge parameters in the returned object,
#' \code{\link[igraph]{igraph.plotting}} for the detailed description of the
#' plotting parameters and \code{\link[igraph]{layout}} for different layouts.
#-------------------------------------------------------------------------------
getIgraph <- function(mo, pcor_palette = NULL) {
  if (is.null(mo$summary)) {
    stop("The inferred network does not exist.")
  }
  if (!base::requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required.")
  }

  summary = mo$summary[mo$summary$type %in% c('P', 'TP', 'FP'), ]
  if (nrow(summary) > 0) {
    # Re-order summary so that all edges go from "x" to "y"
    for(row in 1:nrow(summary)){
      if(summary[row, "ort_inferred"] == -2){
        summary[row, c("x","y")] = summary[row, c("y","x")]
        summary[row, "ort_inferred"] = 2
        if (  (!is.na(summary[row, "p_y2x"]))
           && (!is.na(summary[row, "p_x2y"])) ) {
          temp <- summary[row, "p_y2x"]
          summary[row, "p_y2x"] <- summary[row, "p_x2y"]
          summary[row, "p_x2y"] <- temp
        }
        if(!is.na(summary[row, "ort_ground_truth"])){
          summary[row, "ort_ground_truth"] = 2
        }
      }
    }
  }

  # Create igraph object from summary
  ig_graph = igraph::graph_from_data_frame(summary,
                                           vertices=colnames(mo$adj_matrix))

  # Set nodes visuals
  igraph::V(ig_graph)$color <- "lightblue"
  igraph::V(ig_graph)$label.cex <- 0.8
  igraph::V(ig_graph)$size <- 12

  if (nrow(summary) == 0) {
    # When no edge, returns immediately the graph, do not define edges visuals
    return (ig_graph)
  }

  # Set correct orientations
  igraph::E(ig_graph)$arrow.mode = rep(0, igraph::gsize(ig_graph))
  igraph::E(ig_graph)$arrow.mode[igraph::E(ig_graph)$ort_inferred == 2]  = 2
  igraph::E(ig_graph)$arrow.mode[igraph::E(ig_graph)$ort_inferred == -2] = 1
  igraph::E(ig_graph)$arrow.mode[igraph::E(ig_graph)$ort_inferred == 6]  = 3

  # Set edges visuals
  min_width = 0.2
  igraph::E(ig_graph)$width <-
    pmax(log10(igraph::E(ig_graph)$info_cond - igraph::E(ig_graph)$cplx),
         min_width)
  igraph::E(ig_graph)$arrow.size <- scales::rescale(igraph::E(ig_graph)$width, to=c(0.2,1))

  # By default, negative pcors are blue, null is dark grey and positive are red
  if ( is.null(pcor_palette) )
    pcor_palette = grDevices::colorRampPalette(c("blue", "darkgrey", "red"))

  igraph::E(ig_graph)$color <- "darkgray"
  edge_colors_indices = sapply(
    igraph::E(ig_graph)$partial_correlation,
    function(pcor) {
      ifelse (is.na (pcor), 101, round (pcor * 100) + 101)
    }
  )
  igraph::E(ig_graph)$color <- pcor_palette(201)[edge_colors_indices]

  return(ig_graph)
}

#' Basic plot function of a miic network inference result
#'
#' @description This function calls \code{\link{export}} to build a
#' plottable object from the result returned by \code{\link{miic}} and plot it.
#'
#' @details See the documentation of \code{\link{export}} for further
#' details.
#'
#' @param x [a miic object, required]
#'
#' The object returned by \code{\link{miic}} execution.
#'
#' @param method  [a string, optional, default value "igraph"]
#'
#' The plotting method, currently only "igraph" is supported.
#'
#' @param pcor_palette [a color palette, optional, default value
#' grDevices::colorRampPalette(c("blue", "darkgrey", "red")]
#'
#' Used to represent the partial correlations (the color of the edges).
#' The palette must be able to handle 201 shades to cover the correlation range
#' from -100 to +100.
#'
#' @param \dots Additional plotting parameters. See the corresponding plot
#' function for the complete list.
#'
#' For igraph, see \code{\link[igraph]{igraph.plotting}}.
#'
#' @export
#'
#' @seealso \code{\link{export}} for graphical exports,
#' \code{\link[igraph]{igraph.plotting}}
#'
plot.miic = function(x, method = 'igraph', pcor_palette = NULL, ...) {
  if (method == 'igraph'){
    if (base::requireNamespace("igraph", quietly = TRUE)) {
      igraph_obj = export (x, 'igraph', pcor_palette = pcor_palette)
      igraph::plot.igraph (igraph_obj, ...)
    } else {
      stop("Package 'igraph' is required.")
    }
  } else {
    stop("Method not supported. See ?export for supported methods.")
  }
}
