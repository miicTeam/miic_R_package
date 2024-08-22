fromStringToNumberArrowType <- function(val) {
  ret <- 0
  if (val == "arrow") {
    ret <- 6
  } else if (val == "TRUE") {
    ret <- 15
  }
  return(ret)
}

#' GraphML converting function for miic graph
#'
#' @description Convert miic graph to [GraphML format](http://graphml.graphdrawing.org/).
#' @param g The graph object returned by [miic][miic()].
#' @param file A string. Path to the output file containing file name without
#' extension (.graphml will be appended).
#' @param layout An optional data frame of 2 (or 3) columns containing the
#' coordinate `x` and `y` for each node. The optional first column can contain
#' node names. If node names is not given, the order of the input file will be
#' assigned to the list of positions.
#' @export
#' @useDynLib miic
#' @md

miic.write.network.cytoscape <- function(g, file, layout = NULL) {
  ##################################### NETWORK IN GRAPHML
  if (missing(file)) {
    stop("The file path is necessary")
  }

  if (is.null(g$all.edges.summary)) {
    stop("The result of the miic execution is required")
  }

  summary <- g$all.edges.summary
  adj_matrix <- g$adj_matrix

  if (is.null(layout)) {
    line <- "<graphml>\n"

    # attributes part nodes
    line <- paste(
      line,
      "\t<key id=\"weight\" for=\"node\" attr.name=\"weight\" attr.type=\"double\">\n",
      sep = ""
    )
    line <- paste(line, "\t\t<default>0.2</default>\n", sep = "")
    line <- paste(line, "\t</key>\n", sep = "")
    line <- paste(
      line,
      "\t<key id=\"label\" for=\"node\" attr.name=\"label\" attr.type=\"string\"/>\n",
      sep = ""
    )

    # attributes part edges
    line <- paste(
      line,
      "\t<key id=\"weight\" for=\"edge\" attr.name=\"weight\" attr.type=\"double\"/>\n",
      sep = ""
    )
    line <- paste(
      line,
      "\t<key id=\"label\" for=\"edge\" attr.name=\"label\" attr.type=\"string\"/>\n",
      sep = ""
    )
    line <- paste(
      line,
      "\t<key id=\"sourceArrowShape\" for=\"edge\" attr.name=\"sourceArrowShape\" attr.type=\"string\"/>\n",
      sep = ""
    )
    line <- paste(
      line,
      "\t<key id=\"targetArrowShape\" for=\"edge\" attr.name=\"targetArrowShape\" attr.type=\"string\"/>\n",
      sep = ""
    )
    line <- paste(
      line,
      "\t<key id=\"upstream\" for=\"edge\" attr.name=\"upstream\" attr.type=\"string\"/>\n",
      sep = ""
    )
    line <- paste(
      line,
      "\t<key id=\"complexity\" for=\"edge\" attr.name=\"complexity\" attr.type=\"double\"/>\n",
      sep = ""
    )
    line <- paste(
      line,
      "\t<key id=\"nSamples\" for=\"edge\" attr.name=\"nSamples\" attr.type=\"int\"/>\n",
      sep = ""
    )
    line <- paste(
      line,
      "\t<key id=\"info_shifted\" for=\"edge\" attr.name=\"info_shifted\" attr.type=\"double\"/>\n",
      sep = ""
    )
    line <- paste(
      line,
      "\t<key id=\"confidenceRatio\" for=\"edge\" attr.name=\"confidenceRatio\" attr.type=\"double\"/>\n",
      sep = ""
    )
    line <- paste(
      line,
      "\t<key id=\"sign\" for=\"edge\" attr.name=\"sign\" attr.type=\"string\"/>\n",
      sep = ""
    )
    line <- paste(
      line,
      "\t<key id=\"partialCorrelation\" for=\"edge\" attr.name=\"partialCorrelation\" attr.type=\"double\"/>\n",
      sep = ""
    )
    line <- paste(
      line,
      "\t<key id=\"edgeType\" for=\"edge\" attr.name=\"edgeType\" attr.type=\"int\"/>\n",
      sep = ""
    )
    line <- paste(line, "\n", sep = "")
    line <- paste(line, "\t<graph edgedefault=\"directed\">\n", sep = "")

    # cicle on nodes
    for (node in colnames(adj_matrix)) {
      line <- paste(line, "\t\t<node id=\"", node, "\">\n", sep = "")
      line <- paste(line,
        "\t\t\t<data key=\"label\">",
        node,
        "</data>\n",
        sep = ""
      )
      line <- paste(line, "\t\t\t<data key=\"weight\">0.5</data>\n",
        sep =
          ""
      )
      line <- paste(line, "\t\t</node>\n", sep = "")
    }

    line <- paste(line, "\n", sep = "")
    indexes <- which(summary["type"] == "P" |
      summary["type"] == "TP" | summary["type"] == "FP")

    # cicle on edges
    for (index in indexes) {
      if (!is.na(summary[index, "info_shifted"])) {
        weigth <- summary[index, "info_shifted"]
      } else {
        weigth <- (summary[index, "partial_correlation"])
      }
      if (summary[index, "ort_inferred"] == 1) {
        line <- paste(
          line,
          "\t\t<edge target=\"",
          summary[index, 2],
          "\" source=\"",
          summary[index, 1],
          "\" directed=\"false\">\n",
          sep = ""
        )
        line <- paste(
          line,
          "\t\t\t<data key=\"label\">",
          summary[index, 2],
          "---",
          summary[index, 1],
          "</data>\n",
          sep = ""
        )
        line <- paste(line, "\t\t\t<data key=\"edgeType\">1</data>\n",
          sep =
            ""
        )
      } else if (summary[index, "ort_inferred"] == 2) {
        line <- paste(
          line,
          "\t\t<edge target=\"",
          summary[index, 2],
          "\" source=\"",
          summary[index, 1],
          "\" directed=\"true\">\n",
          sep = ""
        )
        line <- paste(line,
          "\t\t\t<data key=\"sourceArrowShape\">none</data>\n",
          sep = ""
        )
        if (is.na(summary[index, "partial_correlation"])) {
          line <- paste(line,
            "\t\t\t<data key=\"targetArrowShape\">arrow</data>\n",
            sep = ""
          )
          line <- paste(
            line,
            "\t\t\t<data key=\"label\">",
            summary[index, 2],
            "--&gt;",
            summary[index, 1],
            "</data>\n",
            sep = ""
          )
        } else {
          if (summary[index, "partial_correlation"] > 0) {
            line <- paste(line,
              "\t\t\t<data key=\"targetArrowShape\">arrow</data>\n",
              sep = ""
            )
            line <- paste(
              line,
              "\t\t\t<data key=\"label\">",
              summary[index, 2],
              "--&gt;",
              summary[index, 1],
              "</data>\n",
              sep = ""
            )
          } else {
            line <- paste(line,
              "\t\t\t<data key=\"targetArrowShape\">T</data>\n",
              sep = ""
            )
            line <- paste(
              line,
              "\t\t\t<data key=\"label\">",
              summary[index, 2],
              "--|",
              summary[index, 1],
              "</data>\n",
              sep = ""
            )
          }
        }
        line <- paste(line, "\t\t\t<data key=\"edgeType\">2</data>\n",
          sep =
            ""
        )
      } else if (summary[index, "ort_inferred"] == -2) {
        line <- paste(
          line,
          "\t\t<edge target=\"",
          summary[index, 1],
          "\" source=\"",
          summary[index, 2],
          "\" directed=\"true\">\n",
          sep = ""
        )
        line <- paste(line,
          "\t\t\t<data key=\"sourceArrowShape\">none</data>\n",
          sep = ""
        )
        if (is.na(summary[index, "partial_correlation"])) {
          line <- paste(line,
            "\t\t\t<data key=\"targetArrowShape\">arrow</data>\n",
            sep = ""
          )
          line <- paste(
            line,
            "\t\t\t<data key=\"label\">",
            summary[index, 1],
            "--&gt;",
            summary[index, 2],
            "</data>\n",
            sep = ""
          )
        } else {
          if (summary[index, "partial_correlation"] > 0) {
            line <- paste(line,
              "\t\t\t<data key=\"targetArrowShape\">arrow</data>\n",
              sep = ""
            )
            line <- paste(
              line,
              "\t\t\t<data key=\"label\">",
              summary[index, 1],
              "--&gt;",
              summary[index, 2],
              "</data>\n",
              sep = ""
            )
          } else {
            line <- paste(line,
              "\t\t\t<data key=\"targetArrowShape\">T</data>\n",
              sep = ""
            )
            line <- paste(
              line,
              "\t\t\t<data key=\"label\">",
              summary[index, 1],
              "--|",
              summary[index, 2],
              "</data>\n",
              sep = ""
            )
          }
        }
        line <- paste(line, "\t\t\t<data key=\"edgeType\">2</data>\n",
          sep =
            ""
        )
      } else if (summary[index, "ort_inferred"] == 6) {
        line <- paste(
          line,
          "\t\t<edge target=\"",
          summary[index, 2],
          "\" source=\"",
          summary[index, 1],
          "\" directed=\"true\">\n",
          sep = ""
        )
        if (is.na(summary[index, "partial_correlation"])) {
          line <- paste(line,
            "\t\t\t<data key=\"sourceArrowShape\">arrow</data>\n",
            sep = ""
          )
          line <- paste(line,
            "\t\t\t<data key=\"targetArrowShape\">arrow</data>\n",
            sep = ""
          )
          line <- paste(
            line,
            "\t\t\t<data key=\"label\">",
            summary[index, 2],
            "&lt;-&gt;",
            summary[index, 1],
            "</data>\n",
            sep = ""
          )
        } else {
          if (summary[index, "partial_correlation"] > 0) {
            line <- paste(line,
              "\t\t\t<data key=\"sourceArrowShape\">arrow</data>\n",
              sep = ""
            )
            line <- paste(line,
              "\t\t\t<data key=\"targetArrowShape\">arrow</data>\n",
              sep = ""
            )
            line <- paste(
              line,
              "\t\t\t<data key=\"label\">",
              summary[index, 2],
              "&lt;-&gt;",
              summary[index, 1],
              "</data>\n",
              sep = ""
            )
          } else {
            line <- paste(line,
              "\t\t\t<data key=\"sourceArrowShape\">T</data>\n",
              sep = ""
            )
            line <- paste(line,
              "\t\t\t<data key=\"targetArrowShape\">T</data>\n",
              sep = ""
            )
            line <- paste(
              line,
              "\t\t\t<data key=\"label\">",
              summary[index, 2],
              "|-|",
              summary[index, 1],
              "</data>\n",
              sep = ""
            )
          }
        }
        line <- paste(line, "\t\t\t<data key=\"edgeType\">6</data>\n",
          sep =
            ""
        )
      }

      if (!all(is.na(summary[, "info_shifted"]))) {
        if (summary[index, "info_shifted"] <= 1) {
          value <- 1
        } else if (summary[index, "info_shifted"] >= 20) {
          value <- 8
        } else {
          value <- summary[index, "info_shifted"] * 8 / 20
        }
      } else {
        value <- (abs(summary[index, "partial_correlation"]) + 1) * 4
      }
      line <- paste(line,
        "\t\t\t<data key=\"weight\">",
        value,
        "</data>\n",
        sep = ""
      )
      line <- paste(line,
        "\t\t\t<data key=\"upstream\">",
        summary[index, "ai"],
        "</data>\n",
        sep = ""
      )
      line <- paste(line,
        "\t\t\t<data key=\"complexity\">",
        summary[index, "cplx"],
        "</data>\n",
        sep = ""
      )
      line <- paste(line,
        "\t\t\t<data key=\"nSamples\">",
        summary[index, "Nxy_ui"],
        "</data>\n",
        sep = ""
      )
      line <- paste(line,
        "\t\t\t<data key=\"confidenceRatio\">",
        summary[index, "confidence_ratio"],
        "</data>\n",
        sep = ""
      )
      line <- paste(line,
        "\t\t\t<data key=\"info_shifted\">",
        summary[index, "info_shifted"],
        "</data>\n",
        sep = ""
      )
      line <- paste(line,
        "\t\t\t<data key=\"sign\">",
        summary[index, "sign"],
        "</data>\n",
        sep = ""
      )
      line <- paste(line,
        "\t\t\t<data key=\"partialCorrelation\">",
        summary[index, "partial_correlation"],
        "</data>\n",
        sep = ""
      )
      line <- paste(line, "\t\t</edge>\n", sep = "")
    }
    line <- paste(line, "\t</graph>\n", sep = "")
    line <- paste(line, "</graphml>\n", sep = "")
    writeLines(line, paste(file, ".graphml", sep = ""))
  } else {
    if (ncol(layout) == 2) {
      xcol <- 1
      ycol <- 2
      rownames(layout) <- colnames(adj_matrix)
    } else {
      xcol <- 2
      ycol <- 3
      rownames(layout) <- layout[, 1]
    }

    line <- "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
    line <- paste(line, "<graph label=\"graph\"", sep = "")
    line <- paste(line,
      " xmlns:dc=\"http://purl.org/dc/elements/1.1/\"",
      sep = ""
    )
    line <- paste(line,
      " xmlns:xlink=\"http://www.w3.org/1999/xlink\"",
      sep = ""
    )
    line <- paste(line,
      " xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"",
      sep = ""
    )
    line <- paste(line, " xmlns:cy=\"http://www.cytoscape.org\"", sep = "")
    line <- paste(line, " xmlns=\"http://www.cs.rpi.edu/XGMML\"", sep = "")
    line <- paste(line, " directed=\"1\">\n", sep = "")
    # cicle on nodes
    for (node in colnames(adj_matrix)) {
      line <- paste(line,
        "\t\t<node label=\"",
        node,
        "\" id=\"",
        node,
        "\">\n",
        sep = ""
      )
      line <- paste(line,
        "\t\t\t<att name=\"size\" type=\"integer\" value=\"32\"/>\n",
        sep = ""
      )
      x <- layout[node, xcol] * 10
      y <- -layout[node, ycol] * 10
      line <- paste(
        line,
        "\t\t\t<graphics fill=\"#f5f5f5\" x=\"",
        x,
        "\" y=\"",
        y,
        "\" cy:nodeLabelFont=\"Arial-0-11\" labelanchor=\"c\" type=\"ELLIPSE\" cy:nodeTransparency=\"0.8\" h=\"32\" width=\"1\" outline=\"#666666\" w=\"32\"/>\n",
        sep = ""
      )
      line <- paste(line, "\t\t</node>\n", sep = "")
    }

    line <- paste(line, "\n", sep = "")

    indexes <- which(summary["type"] == "P" |
      summary["type"] == "TP" |
      summary["type"] == "FP")

    # cycle on edges
    for (index in indexes) {
      sourceArrowNum <- 0
      targetArrowNum <- 0
      if (summary[index, "ort_inferred"] == 1) {
        line <- paste(
          line,
          "\t\t<edge label=\"",
          summary[index, 2],
          "---",
          summary[index, 1],
          "\" target=\"",
          summary[index, 2],
          "\" source=\"",
          summary[index, 1],
          "\">\n",
          sep = ""
        )
        line <- paste(line,
          "\t\t\t<att name=\"edgeType\" type=\"integer\" value=\"1\"/>\n",
          sep = ""
        )
      } else if (summary[index, "ort_inferred"] == 2) {
        if (is.na(summary[index, "partial_correlation"])) {
          value <- "arrow"
          varchar <- intToUtf8(187)
          label <- paste(summary[index, 1], "--&gt;", summary[index, 2],
            sep =
              ""
          )
        } else {
          if (summary[index, "partial_correlation"] > 0) {
            value <- "arrow"
            varchar <- intToUtf8(187)
            label <- paste(summary[index, 1], "--&gt;", summary[index, 2],
              sep =
                ""
            )
          } else {
            value <- "T"
            varchar <- "|"
            label <- paste(summary[index, 1], "--|", summary[index, 2],
              sep =
                ""
            )
          }
        }
        line <- paste(
          line,
          "\t\t<edge label=\"",
          label,
          "\" target=\"",
          summary[index, 2],
          "\" source=\"",
          summary[index, 1],
          "\">\n",
          sep = ""
        )
        line <- paste(
          line,
          "\t\t\t<att name=\"targetArrowShape\" type=\"string\" value=\"",
          value,
          "\"/>\n",
          sep = ""
        )
        line <- paste(
          line,
          "\t\t\t<att name=\"sourceArrowShape\" type=\"string\" value=\"none\"/>\n",
          sep = ""
        )
        line <- paste(line,
          "\t\t\t<att name=\"edgeType\" type=\"integer\" value=\"2\"/>\n",
          sep = ""
        )
        sourceArrowNum <- 0
        targetArrowNum <- fromStringToNumberArrowType(value)
      } else if (summary[index, "ort_inferred"] == -2) {
        if (is.na(summary[index, "partial_correlation"])) {
          value <- "arrow"
          varchar <- intToUtf8(187)
          label <- paste(summary[index, 2], "--&gt;", summary[index, 1],
            sep =
              ""
          )
        } else {
          if (summary[index, "partial_correlation"] > 0) {
            value <- "arrow"
            varchar <- intToUtf8(187)
            label <- paste(summary[index, 2], "--&gt;", summary[index, 1],
              sep =
                ""
            )
          } else {
            value <- "T"
            varchar <- "|"
            label <- paste(summary[index, 2], "--|", summary[index, 1],
              sep =
                ""
            )
          }
        }
        line <- paste(
          line,
          "\t\t<edge label=\"",
          label,
          "\" target=\"",
          summary[index, 1],
          "\" source=\"",
          summary[index, 2],
          "\">\n",
          sep = ""
        )
        line <- paste(
          line,
          "\t\t\t<att name=\"targetArrowShape\" type=\"string\" value=\"",
          value,
          "\"/>\n",
          sep = ""
        )
        line <- paste(
          line,
          "\t\t\t<att name=\"sourceArrowShape\" type=\"string\" value=\"none\"/>\n",
          sep = ""
        )
        line <- paste(line,
          "\t\t\t<att name=\"edgeType\" type=\"integer\" value=\"2\"/>\n",
          sep = ""
        )
        sourceArrowNum <- 0
        targetArrowNum <- fromStringToNumberArrowType(value)
      } else if (summary[index, "ort_inferred"] == 6) {
        if (is.na(summary[index, "partial_correlation"])) {
          value <- "arrow"
          varchar <- intToUtf8(187)
          label <- paste(summary[index, 2], "&lt;-&gt;", summary[index, 1],
            sep =
              ""
          )
        } else {
          if (summary[index, "partial_correlation"] > 0) {
            value <- "arrow"
            varchar <- intToUtf8(187)
            label <- paste(summary[index, 2], "&lt;-&gt;", summary[index, 1],
              sep =
                ""
            )
          } else {
            value <- "T"
            varchar <- "|"
            label <- paste(summary[index, 2], "|-|", summary[index, 1],
              sep =
                ""
            )
          }
        }
        line <- paste(
          line,
          "\t\t<edge label=\"",
          label,
          "\" target=\"",
          summary[index, 1],
          "\" source=\"",
          summary[index, 2],
          "\">\n",
          sep = ""
        )
        line <- paste(
          line,
          "\t\t\t<att name=\"targetArrowShape\" type=\"string\" value=\"",
          value,
          "\"/>\n",
          sep = ""
        )
        line <- paste(
          line,
          "\t\t\t<att name=\"sourceArrowShape\" type=\"string\" value=\"",
          value,
          "\"/>\n",
          sep = ""
        )
        line <- paste(line,
          "\t\t\t<att name=\"edgeType\" type=\"integer\" value=\"6\"/>\n",
          sep = ""
        )

        sourceArrowNum <- fromStringToNumberArrowType(value)
        targetArrowNum <- fromStringToNumberArrowType(value)
      }

      if (summary[index, "info_shifted"] <= 1) {
        value <- 1
      } else if (summary[index, "info_shifted"] >= 20) {
        value <- 8
      } else {
        value <- summary[index, "info_shifted"] * 8 / 20
      }

      line <- paste(
        line,
        "\t\t\t<att name=\"weight\" type=\"double\" value=\"",
        value,
        "\"/>\n",
        sep = ""
      )
      line <- paste(
        line,
        "\t\t\t<att name=\"upstream\" type=\"string\" value=\"",
        summary[index, "ai"],
        "\"/>\n",
        sep = ""
      )
      line <- paste(
        line,
        "\t\t\t<att name=\"complexity\" type=\"double\" value=\"",
        summary[index, "cplx"],
        "\"/>\n",
        sep = ""
      )
      line <- paste(
        line,
        "\t\t\t<att name=\"nSamples\" type=\"integer\" value=\"",
        summary[index, "Nxy_ai"],
        "\"/>\n",
        sep = ""
      )
      line <- paste(
        line,
        "\t\t\t<att name=\"confidenceRatio\" type=\"double\" value=\"",
        summary[index, "confidence_ratio"],
        "\"/>\n",
        sep = ""
      )
      line <- paste(
        line,
        "\t\t\t<att name=\"info_shifted\" type=\"double\" value=\"",
        summary[index, "info_shifted"],
        "\"/>\n",
        sep = ""
      )
      line <- paste(
        line,
        "\t\t\t<att name=\"sign\" type=\"string\" value=\"",
        summary[index, "sign"],
        "\"/>\n",
        sep = ""
      )
      line <- paste(
        line,
        "\t\t\t<att name=\"partialCorrelation\" type=\"double\" value=\"",
        summary[index, "partial_correlation"],
        "\"/>\n",
        sep = ""
      )

      if (summary[index, "sign"] == "+") {
        fillColor <- "#ff3300"
      } else if (summary[index, "sign"] == "-") {
        fillColor <- "#aaaaff"
      } else {
        fillColor <- "#808080"
      }
      line <- paste(
        line,
        "\t\t\t<graphics cy:sourceArrowColor=\"#000000\" cy:targetArrowColor=\"#000000\" width=\"",
        value,
        "\" cy:edgeLineType=\"SOLID\" cy:targetArrow=\"",
        targetArrowNum,
        "\" cy:sourceArrow=\"",
        sourceArrowNum,
        "\" fill=\"",
        fillColor,
        "\"/>\n",
        sep = ""
      )
      line <- paste(line, "\t\t</edge>\n", sep = "")
    }
    line <- paste(line, "\t</graph>\n", sep = "")
    # name = basename(file_path_sans_ext(outDirPath))
    writeLines(line, paste(file, ".graphml", sep = ""))
  }
}
