% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmiic.plot.R
\name{plot.tmiic}
\alias{plot.tmiic}
\title{Basic plot function of a temporal miic (tmiic) network inference result}
\usage{
\method{plot}{tmiic}(
  x,
  display = "compact",
  show_self_loops = TRUE,
  positioning_for_grid = "greedy",
  orientation_for_grid = "L",
  method = "igraph",
  pcor_palette = NULL,
  ...
)
}
\arguments{
\item{x}{[a tmiic graph object]
The graph object returned by \code{\link{miic}} in temporal mode}

\item{display}{[a string]. Optional, default value "compact".
Possible values are \emph{"raw"}, \emph{"lagged"}, \emph{"compact"},
\emph{"combine"}, \emph{"unique"}, \emph{"drop"}:
\itemize{
\item When \emph{display} = \emph{"raw"}, the plot function will
  use the tmiic graph object as it, leading to the display of a lagged
  graph. Unless a specific layout is specified, nodes will be positioned
  on a grid.
\item When \emph{display} = \emph{"lagged"}, the function will
  repeat the edges over history assuming stationarity and plot a lagged
  graph. Unless a specific layout is specified, nodes will be positioned
  on a grid.
\item When \emph{display} = \emph{"compact"}, the default, nodes
  and edges are converted into a flattened version to produce a compact
  view of the temporal network whilst still presenting all the information
  in the plotting.\cr
  i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 become respectively X->Y lag=1,
  X<-Y lag=2.
\item When \emph{display} = \emph{"combine"}, prior to the plotting,
  a preprocessing will be applied to kept only one edge
  per couple of nodes. The info_shifted will be the highest one
  of the summarized edges whilst the lag and orientation of the
  summarized edge will be an aggregation.\cr
  i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 will become X<->Y lag=1,2 with
  the info_shifted of X_lag1->Y_lag0 if info_shifted of
  X_lag1->Y_lag0 > X_lag2<-Y_lag0.
\item When \emph{display} = \emph{"unique"}, prior to the plotting,
  a preprocessing will be applied to kept only the edges having the
  highest info_shifted for a couple of nodes.
  If several edges between the sames nodes have the same
  info_shifted, then the edge kept is the one with the minimum lag.\cr
  i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 with info_shifted of
  X_lag1->Y_lag0 > X_lag2<-Y_lag0 become X->Y lag=1.
\item When \emph{display} = \emph{"drop"}, prior to the plotting,
  a preprocessing will be applied to kept only the edges having the
  highest info_shifted for a couple of nodes.
  If several edges between the sames nodes have the same
  info_shifted, then the edge kept is the one with the minimum lag.\cr
  i.e. :  X_lag1->Y_lag0, X_lag2<-Y_lag0 with info_shifted of
  X_lag1->Y_lag0 > X_lag2<-Y_lag0 become X->Y.
  The lag information is dropped during the preprocessing and
  will not be displayed on the final plotting.
}}

\item{show_self_loops}{[a boolean] Optional, TRUE by default.
When TRUE, the edges like X_lag0-X_lag1 are included in the iGraph object.
When FALSE, only edges having different nodes are present in the iGraph
object.}

\item{positioning_for_grid}{[a string] Optional, "greedy" by default.
Used only when the display is "raw" or "lagged and no layout is supplied.
Possible values are \emph{"none"}, \emph{"alphabetical"}, \emph{"layers"}
\emph{"greedy"} and \emph{"sugiyama"}
\itemize{
\item When \emph{positioning_for_grid} = \emph{"none"}
 The nodes are positioned as they appear in the miic result
\item When \emph{positioning_for_grid} = \emph{"alphabetical"}
 The nodes are positioned alphabeticaly in ascending order
\item When \emph{positioning_for_grid} = \emph{"layers"}
 The nodes with the less lags wil be placed on the exteriors
 while the nodes having the most lags are in the center
\item When \emph{positioning_for_grid} = \emph{"greedy"}
 A greedy algorithm will be used to placed the nodes in a way minimizing
 the crossing edges
\item When \emph{positioning_for_grid} = \emph{"sugiyama"}
 The sugiyama algorithm will be used to placed the nodes in a way
 minimizing the crossing edges
}}

\item{orientation_for_grid}{[a string] Optional, "L" by default.
Used only when the display is "raw" or "lagged and no layout is supplied.
Indicates the orientation of the draw, possible values are landscape: "L"
or portrait: "P".}

\item{method}{A string representing the plotting method. Default to "igraph".
Currently only "igraph" is supported.}

\item{pcor_palette}{Optional. The color palette used to represent the partial
correlations (the color of the edges). See \code{\link{getIgraph}} for details.}

\item{\dots}{Additional plotting parameters. See the corresponding plot function
for the complete list.

For igraph, see \code{\link[igraph]{igraph.plotting}}.}
}
\description{
This function calls \code{\link{tmiic.export}} to build a
plottable object from the result returned by \code{\link{miic}} in
temporal mode and plot it.
}
\details{
See the documentation of \code{\link{tmiic.export}} for further
details.
}
\examples{
\donttest{
library(miic)

#' # EXAMPLE COVID CASES (time series demo)
data(covidCases)
# execute MIIC (reconstruct graph in temporal mode)
tmiic_res <- miic(input_data = covidCases, mode = "TS", n_layers = 3, delta_t = 1, movavg = 14)

# to plot the default compact graph
if(require(igraph)) {
  plot(tmiic_res)
}

# to plot the raw temporal network Using igraph
if(require(igraph)) {
  plot(tmiic_res, display="raw")
}

# to plot the full temporal network Using igraph
if(require(igraph)) {
  plot(tmiic_res, display="lagged")
}

}
}
\seealso{
\code{\link{tmiic.export}} for generic exports,
\code{\link[igraph]{igraph.plotting}}
}
