isContinuous <- function(tmpArray) {
  nbLevels <- unique(tmpArray)
  nbSamples <- length(tmpArray)
  if (length(nbLevels) >= 0.8 * nbSamples |
    length(nbLevels) >= 100) {
    print(paste(
      "VARIABLE CONTINUOUS LIKE:",
      length(nbLevels),
      "different levels"
    ))
    return(1)
  } else {
    return(0)
  }
}

checkInput <- function(dataFile, method) {
  errCode <- "0"
  isCnt_test <- 0
  isContinuousArg <- 0

  str <- unlist(paste(as.character(dataFile), collapse = ""))
  str_names <- unlist(strsplit(str, "\n"))[1]
  if (grepl("#", str) | grepl("&", str) | grepl("'", str)) {
    errCode <- "111"
  } else {
    isCnt_test <- 0 # sum(apply( dataFile, 2, isContinuous ))

    if (method == "miic" &
      length(unique(dataFile[1, ])) != ncol(dataFile)) {
      errCode <- "117"
    } else if (method == "miic" &
      isCnt_test > 0.1 * ncol(dataFile) &
      isContinuousArg == 0) {
      errCode <- "118"
    } else {
      if (method == "miic" & isCnt_test > 0 & isContinuousArg == 0) {
        errCode <- "018"
      }
    }
  }
  return(errCode)
}


checkTrueEdges <- function(edgesFile) {
  errCode <- "0"
  str <- unlist(paste(as.character(edgesFile), collapse = ""))
  if (grepl("#", str) | grepl("&", str) | grepl("'", str)) {
    errCode <- "121"
  } else {
    if (ncol(edgesFile) != 2 & ncol(edgesFile) != 3) {
      errCode <- "023"
    }
  }
  return(errCode)
}

checkLayout <- function(layoutFile) {
  errCode <- "0"

  str <- unlist(paste(as.character(layoutFile), collapse = ""))
  if (grepl("#", str) | grepl("&", str) | grepl("'", str)) {
    errCode <- "131"
  } else {
    if (ncol(layoutFile) != 2 & ncol(layoutFile) != 3) {
      errCode <- "033"
    } else {
      if (ncol(layoutFile) == 2) {
        if (!is.numeric(layoutFile[, 1]) | !is.numeric(layoutFile[, 2])) {
          errCode <- "034"
        }
      } else if (ncol(layoutFile) == 3) {
        if (!is.numeric(layoutFile[, 2]) | !is.numeric(layoutFile[, 3])) {
          errCode <- "038"
        }
      }
    }
  }
  return(errCode)
}

checkStateOrder <- function(stateOrderFile, dataFile) {
  errCode <- "0"

  str <- unlist(paste(as.character(stateOrderFile), collapse = ""))
  str_names <- unlist(strsplit(str, "\n"))[1]
  if (grepl("#", str) | grepl("&", str) | grepl("'", str)) {
    errCode <- "141"
  } else {
    rownames(stateOrderFile) <- stateOrderFile[, "var_names"]
    myVariables <- rownames(stateOrderFile)

    if (ncol(stateOrderFile) != 2 & ncol(stateOrderFile) != 3) {
      errCode <- "043"
    }
  }
  return(errCode)
}

errorCodeToString <- function(error_code) {
  errorList1 <- list(
    "0" = "Warning:",
    "1" = "Fatal error:"
  )
  errorList2 <- list(
    "0" = "Unknown Error",
    "1" = "input data frame",
    "2" = "trueEdge data frame",
    "3" = "layout data frame",
    "4" = "cathegory order data frame"
  )

  errorList3 <- list(
    "0" = "does not exist",
    "1" = paste(
      "is not readable, check the file format. Special characters ",
      "like #, &, ' are not allowed."
    ),
    "2" = "should not have rownames",
    "3" = "should have exactly two columns",
    "4" = "should be numerical",
    "5" = "should not have column names",
    "6" = "problem as states in stateOrderFile are not consistent with dataset",
    "7" = "has duplicated column names",
    "8" = paste(
      "contains one or several variables with too many different ",
      " levels (might be continuous data)"
    ),
    "9" = "occured."
  )
  error_string <- unlist(strsplit(error_code, ""))
  return(paste(errorList1[[error_string[1]]], errorList2[[error_string[2]]],
    errorList3[[error_string[3]]],
    sep = " "
  ))
}

fromStringToNumberArrowType <- function(val) {
  ret <- 0
  if (val == "arrow") {
    ret <- 6
  } else if (val == "T") {
    ret <- 15
  }

  return(ret)
}
