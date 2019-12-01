fromStringToNumberArrowType <- function(val) {
  ret = 0
  if (val == "arrow") {
    ret = 6
  } else if (val == "T") {
    ret = 15
  }
  return(ret)
}

dataToStateOrder <- function(myGv, stateOrder, inputDataFile) {
  # Check if a description file of the state order is specified
  if (!is.null(stateOrder)) {
    myStateOrder = stateOrder
    myStateOrder = as.matrix(myStateOrder)

    if (!"var_names" %in% colnames(myStateOrder))
      stop(
        paste(
          "State order must have column names as 'var_names'",
          "and 'levels_increasing_order'"
        )
      )

    rownames(myStateOrder) = myStateOrder[, "var_names"]
    myVariables = rownames(myStateOrder)

    for (myVar in myVariables) {
      myCatStr = unlist(strsplit(myStateOrder[myVar,
                                              "levels_increasing_order"],
                                 ", "))
      myGv$data[, myVar] = as.character(myGv$data[, myVar])
      ### Security part:
      ### check if myCatStr corresponds to unique(myData[,myVar])
      mySafety = sort(unique(myGv$data[which(!is.na(myGv$data[, myVar])),
                                       myVar]))
      if (identical(mySafety, sort(myCatStr) == FALSE))
        stop(
          paste(
            "## --w--> WARNING !! The state order file does not ",
            "correspond to the dataset for ",
            myVar,
            "\n"
          )
        )
      #print(myGv$data[1:5,myVar])
      ### Set the non-NA datas as integer corresponding to the correct levels
      myGv$data[which(!is.na(myGv$data[, myVar])), myVar] =
        as.integer(factor(myGv$data[which(!is.na(myGv$data[, myVar])), myVar],
                          levels = myCatStr))
      ### Convert the NA Values as well

      myGv$data[, myVar] = as.integer(myGv$data[, myVar])
      #print(myGv$data[1:5,myVar])
    }
  }	else {
    ### Compute the stateOrder.tsv file
    myStateOrder = matrix(
      NA,
      nrow = ncol(myGv$data),
      ncol = 2,
      dimnames = list(
        colnames(myGv$data),
        c("var_names",
          "levels_increasing_order")
      )
    )
    #inputDataDir = unlist(strsplit(inputDataFile, "/"))
    #inputDataDir =
    # paste(inputDataDir[1:length(inputDataDir)-1], collapse = "/")
    ### Fill the table with the "unique" vectors
    for (var in colnames(myGv$data)) {
      ### Get the vector of unique values
      tmp.unique = unique(stats::na.omit(myGv$data[, var]))
      myStateOrder[var, "var_names"] = var
      myStateOrder[var, "levels_increasing_order"] = paste(tmp.unique,
                                                           collapse = ", ")
    }


    myGv$data[, colnames(myGv$data)] <-
      as.data.frame(lapply(myGv$data[, colnames(myGv$data)], factor))
    myGv$data[, colnames(myGv$data)] <-
      as.data.frame(lapply(myGv$data[, colnames(myGv$data)], as.numeric))
    myGv$data <- as.matrix(myGv$data)
  }
  myStateOrder
}

computeSign <- function(edgeList, myGv) {
  edgeList.signed <- cbind(edgeList,
                           c(rep(NA, nrow(edgeList))),
                           c(rep(NA, nrow(edgeList))))
  colnames(edgeList.signed) <- c(colnames(edgeList), "Sign", "Diff")
  for (edge in 1:nrow(edgeList)) {
    if (edgeList[edge, "infOrt"] != 1 &
        edgeList[edge, "infOrt"] != 0 &
        edgeList[edge, "infOrt"] != 6) {
      ### Deal with the "reversed" edges
      if (edgeList[edge, "infOrt"] < 0) {
        x <- edgeList[edge, "y"]
        y <- edgeList[edge, "x"]
      } else {
        x <- edgeList[edge, "x"]
        y <- edgeList[edge, "y"]
      }

      ### Get the different max/min for the two considered variables
      xMax <- max(myGv$data[, x], na.rm = T)
      xMin <- min(myGv$data[, x], na.rm = T)
      yMax <- max(myGv$data[, y], na.rm = T)
      yMin <- min(myGv$data[, y], na.rm = T)

      ### Gather the different counts needed to compute the probabilities
      y.xMax <- myGv$data[which(myGv$data[, x] == xMax), y]
      y.xMin <- myGv$data[which(myGv$data[, x] == xMin), y]
      yMax.xMax <- length(which(y.xMax == yMax))
      yMax.xMin <- length(which(y.xMin == yMax))

      if (sign(yMax.xMax / length(y.xMax) - yMax.xMin / length(y.xMin)) == 1) {
        edgeList.signed[edge, "Sign"] <- "+"
      } else {
        edgeList.signed[edge, "Sign"] <- "-"
      }
      edgeList.signed[edge, "Diff"] <-
        (yMax.xMax / length(y.xMax)) - (yMax.xMin / length(y.xMin))
    } else {
      edgeList[edge, "Sign"] <- NA
      edgeList[edge, "Diff"] <- 0
    }
  }
  return(edgeList.signed[, c(ncol(edgeList.signed) - 1,
                             ncol(edgeList.signed))])
}


# To compute sign with partial correlation, based on the previously
# determined ai set
computeSign.pcor <- function(edgeList, myGv) {
  edgeList.signed <- cbind(edgeList,
                           c(rep(NA, nrow(edgeList))),
                           c(rep(NA, nrow(edgeList))))
  colnames(edgeList.signed) <- c(colnames(edgeList), "Sign", "Diff")
  tmp.core = NA
  for (edge in 1:nrow(edgeList)) {
    if (edgeList[edge, "infOrt"] != 0) {
      ### get the Ai set of separation if it exists
      x <- myGv$data[, edgeList[edge, "x"]]
      y <- myGv$data[, edgeList[edge, "y"]]
      x.name <- edgeList[edge, "x"]
      y.name <- edgeList[edge, "y"]
      if (!is.na(edgeList[edge, "ai"])) {
        # partial correlation only if ui exists
        if (!is.na(edgeList[edge, "ai"])) {
          ai = myGv$data[, unlist(strsplit(edgeList[edge, "ai"], ","))]
          ai.name = unlist(strsplit(edgeList[edge, "ai"], ","))
          # Remove NA samples on x,y and ai
          NaVal =  sort(unique(which(is.na(
            cbind(x, y, ai)
          ), arr.ind = T)[, 1]))

          if (length(NaVal) != 0) {
            if (length(ai.name) > 1) {
              ai = ai[-NaVal, ]
            } else {
              ai = ai[-NaVal]
            }

            tmp.pcor = try(ppcor::pcor.test(x[-NaVal],
                                            y[-NaVal],
                                            ai,
                                            method = "spearman")$"estimate",
                           silent = T)
            #--- Some of the ai do not have a consistent pcor
            if (class(tmp.pcor) == "try-error") {
              ai.tokeep = c()
              for (n in ai.name) {
                if (class(try(ppcor::pcor.test(x[-NaVal],
                                               y[-NaVal],
                                               ai[, n],
                                               method = "spearman")$"estimate",
                              silent = T)) != "try-error") {
                  ai.tokeep = c(ai.tokeep, n)
                }
              }
              ai = ai[, ai.tokeep]
              tmp.pcor = ppcor::pcor.test(x[-NaVal],
                                          y[-NaVal],
                                          ai,
                                          method = "spearman")$"estimate"
            }
          } else {
            tmp.pcor = ppcor::pcor.test(x, y, ai, method = "spearman")$"estimate"
          }
          # if no ai, we simply compute the spearman's correlation
        } else {
          NaVal = sort(unique(which(is.na(
            cbind(x, y)
          ), arr.ind = T)[, 1]))
          if (length(NaVal) != 0) {
            tmp.pcor = stats::cor(x[-NaVal], y[-NaVal], method = "spearman")
          } else {
            tmp.pcor = stats::cor(x, y, method = "spearman")
          }
        }
        # if no ai, we simply compute the spearman's correlation
      } else {
        NaVal = sort(unique(which(is.na(cbind(
          x, y
        )), arr.ind = T)[, 1]))
        if (length(NaVal) != 0) {
          tmp.pcor = stats::cor(x[-NaVal],
                                y[-NaVal],
                                method = "spearman")
        } else {
          tmp.pcor = stats::cor(x, y, method = "spearman")
        }
      }

      if (sign(tmp.pcor) == 1) {
        edgeList.signed[edge, "Sign"] <- "+"
      } else {
        edgeList.signed[edge, "Sign"] <- "-"
      }
      edgeList.signed[edge, "Diff"] <- (tmp.pcor)
    } else {
      edgeList[edge, "Sign"] <- NA
      edgeList[edge, "Diff"] <- 0
    }
  }
  return(edgeList.signed[, c(ncol(edgeList.signed) - 1, ncol(edgeList.signed))])
}


# Function to compute the signs according to the total correlation (without
# conditionning on the ui set)
computeSign.Cor <- function(edgeList, myGv) {
  edgeList.signed <- cbind(edgeList,
                           c(rep(NA, nrow(edgeList))),
                           c(rep(NA, nrow(edgeList))))
  colnames(edgeList.signed) <- c(colnames(edgeList), "Sign", "Diff")
  for (edge in 1:nrow(edgeList)) {
    if (edgeList[edge, "infOrt"] != 0) {
      x <- myGv$data[, edgeList[edge, "x"]]
      y <- myGv$data[, edgeList[edge, "y"]]
      x.name <- edgeList[edge, "x"]
      y.name <- edgeList[edge, "y"]

      NaVal = sort(unique(which(is.na(cbind(
        x, y
      )), arr.ind = T)[, 1]))
      if (length(NaVal) != 0) {
        tmp.pcor = stats::cor(x[-NaVal], y[-NaVal], method = "spearman")
      } else {
        tmp.pcor = stats::cor(x, y, method = "spearman")
      }

      if (sign(tmp.pcor) == 1) {
        edgeList.signed[edge, "Sign"] <- "+"
      } else {
        edgeList.signed[edge, "Sign"] <- "-"
      }
      edgeList.signed[edge, "Diff"] <- (tmp.pcor)
    } else {
      edgeList[edge, "Sign"] <- NA
      edgeList[edge, "Diff"] <- 0
    }
  }
  return(edgeList.signed[, c(ncol(edgeList.signed) - 1, ncol(edgeList.signed))])
}
### FIN AJOUT LV ####
