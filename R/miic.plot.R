#' Igraph plotting function for miic
#' @description This functions plots the network with the given layout (if specified) using the igraph package.
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
#' @param verbose [a boolean value] If TRUE, debugging output is printed.
#' @export
#' @useDynLib miic

miic.plot <-
  function(g,
           method = "log_confidence",
           igraphLayout = NULL,
           userLayout = NULL,
           verbose = F) {
    if (is.null(g$all.edges.summary)) {
      stop("The learnt graphical model summary does not exist")
    }

    if (verbose) {
      if (is.null(userLayout)) {
        cat("\t# --W--> No file path for the nodes user layout given.\n")
      }
    }

    if (verbose) {
      cat("# --------\n# -> START Plot...\n")
    }

    #### Set Global variables
    myVariables <- colnames(g$adjMatrix)

    #### Read the first line of the inputdata file to get all the properties
    #### saved in the correct order
    gV <- new.env()
    #### ----

    #### Load the vertices positions
    if (verbose == TRUE) {
      cat("# Load the positions of the vertices\n")
    }

    gV$myVerticesPos <- NULL
    if (!is.null(userLayout)) {
      gV$myVerticesPos <- userLayout
      if (ncol(gV$myVerticesPos) > 2) {
        gV$myVerticesPos <- gV$myVerticesPos[, -1]
      }
      gV$myVerticesPos <- as.matrix(gV$myVerticesPos)
    }

    #### summary
    # plots for all technics on partial correlation
    mySummary <- plot.loadSummary(g$all.edges.summary)
    # create the color vector based on partial correlation values


    # handle the layout
    if (!is.null(userLayout)) {
      myLayout <- gV$myVerticesPos
    } else {
      if (is.null(igraphLayout)) {
        myLayout <- igraph::layout_with_kk
        # if (length(myVariables) >=40){
        #  myLayout = igraph::layout.circle
        # }
      } else {
        myLayout <- igraphLayout
      }
    }

    if (method == "pcor") {
      if (!is.na(g$all.edges.summary$partial_correlation[1])) {
        myColors <- pCor.edgeCol(mySummary, myVariables)
        # create the graph object
        myGraph <- modif.Graph(mySummary, myVariables, myColors)
        # plot the Partial Correlation Graph
        graphics::layout(t(1:2), widths = c(5, 1))
        # Set margins and turn all axis labels horizontally (with `las=1`)
        graphics::par(
          mar = rep(.5, 4),
          oma = c(3, 3, 3, 1),
          las = 1
        )
        graphics::plot(myGraph, layout = myLayout)
        blue.gradient <- grDevices::rainbow(100, start = 3 / 6, end = 4 / 6)
        red.gradient <- grDevices::rainbow(100, start = 0, end = 0.16)
        leg_colors <- c(red.gradient, blue.gradient)
        legend_image <-
          grDevices::as.raster(matrix(leg_colors, ncol = 1))
        graphics::plot(
          c(0, 5),
          c(-1, 1),
          type = "n",
          axes = F,
          xlab = "",
          ylab = "",
          main = "Partial correlation",
          cex.main = 0.7
        )
        graphics::text(
          x = 1.5,
          y = seq(-1, 1, l = 5),
          labels = seq(-1, 1, l = 5),
          cex = 0.7
        )
        graphics::rasterImage(legend_image, 2.5, -1, 3.5, 1)
      } else {
        print(
          paste(
            "It is not possible to plot correlation for categorical ",
            "variables without an order"
          )
        )
      }
      # plot the Confidence Graph  only for miic on log_confidence
    } else if (method == "log_confidence") {
      myColors <- conf.edgeCol(mySummary, myVariables)
      myGraph <- modif.Graph(mySummary, myVariables, myColors)
      graphics::layout(t(1:3), widths = c(5, 1, 1))

      # Set margins and turn all axis labels horizontally (with `las=1`)
      graphics::par(
        mar = c(.5, .1, .5, .1),
        oma = c(3, 3, 3, 1),
        las = 1
      )
      graphics::plot(myGraph, layout = myLayout)
      blue.gradient <- grDevices::rainbow(100, start = 3 / 6, end = 4 / 6)
      red.gradient <- grDevices::rainbow(100, start = 0, end = 0.16)
      # Legend
      # Positive correlations

      legend_image <-
        grDevices::as.raster(matrix(red.gradient, ncol = 1))
      graphics::plot(
        c(0, 4),
        c(0.2, 0.8),
        type = "n",
        axes = F,
        xlab = "",
        ylab = ""
      )
      graphics::par(adj = 1)


      if (!is.na(g$all.edges.summary$partial_correlation[1])) {
        leg_colors <- c(red.gradient, blue.gradient)

        graphics::title("Confidence\npcor+",
          cex.main = 1,
          line = -4
        )
        graphics::rasterImage(legend_image, 3.3, 0.25, 3.8, 0.75)
        legend_image <-
          grDevices::as.raster(matrix(rev(blue.gradient), ncol = 1))
        graphics::plot(
          c(0, 4),
          c(0.2, 0.8),
          type = "n",
          axes = F,
          xlab = "",
          ylab = ""
        )
        graphics::par(adj = 0.5)
        graphics::title("(NI' = -log Pxy)\npcor-",
          cex.main = 1,
          line = -4
        )
        graphics::text(
          x = rep(.5, 3),
          y = c(0.26, 0.5, 0.74),
          labels = c("< 1", "50", "> 100"),
          cex = 1
        )
        graphics::rasterImage(legend_image, 1.5, 0.25, 2, 0.75)
      } else {
        leg_colors <- c(red.gradient)
        graphics::title("Confidence   \n(NI' = -log Pxy)",
          cex.main = 1,
          line = -3
        )
        graphics::rasterImage(legend_image, 2.5, 0.25, 3, 0.75)
        graphics::plot(
          c(0, 4),
          c(0.2, 0.8),
          type = "n",
          axes = F,
          xlab = "",
          ylab = ""
        )
        graphics::par(adj = 0.3)
        graphics::text(
          x = rep(.5, 3),
          y = c(0.26, 0.5, 0.74),
          labels = c("< 1", "50", "> 100"),
          cex = 1
        )
        # graphics::rasterImage(legend_image, 1.5, 0.25, 2, 0.75)

        # legend_image <- grDevices::as.raster(matrix(leg_colors, ncol=1))
        # graphics::plot(c(0,5),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "Confidence (NI' = -log Pxy)", cex.main = 0.7)
        # graphics::text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5), cex = 0.7)
        # graphics::rasterImage(legend_image, 3.3, 0, 3.8, 1)
      }
    } else {
      stop("incorrect method argument: must be \"log_confidence\" or \"pcor\" ")
    }
    if (verbose) {
      cat("\t# --------\n# -> END Plot...\n")
    }
  }
