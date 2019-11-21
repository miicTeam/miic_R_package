#!/usr/bin/Rscript

computeSingleSampleStat <- function(myFinalSummary) {
	myKeys = c("raw__TP", "raw__TPnort", "raw__FP", "raw__TN", "raw__FN",
	           "raw__TPess", "raw__FNess",  "precision__", "recall__", "fscore__",
	           "fscore__0_5")
	myRetList = vector("list", length(myKeys))
	names(myRetList) <- myKeys
	myRetList = lapply(myRetList, function(x) x = NA)

	#### Count the number of TP, FP, TN, FN
	myRetList[["raw__TP"]] = length(which(myFinalSummary[, "type"] == "TP"))
	myRetList[["raw__TPnort"]] = length(which(myFinalSummary[, "isOrt"] == "N"))
	myRetList[["raw__FP"]] = length(which(myFinalSummary[, "type"] == "FP"))
	myRetList[["raw__TN"]] = length(which(myFinalSummary[, "type"] == "TN"))
	myRetList[["raw__FN"]] = length(which(myFinalSummary[, "type"] == "FN"))
	myRetList[["raw__TPess"]] = length(which(myFinalSummary[, "type"] == "TP" &
	                                         myFinalSummary[, "essential"] == "Y"))
	myRetList[["raw__FNess"]] = length(which(myFinalSummary[, "type"] == "FN" &
	                                         myFinalSummary[, "essential"] == "Y"))
				
	#### Compute precision, recall, fscore
	myRetList[["precision__"]] = ((myRetList[["raw__TP"]]) /
                                (myRetList[["raw__TP"]] +
                                 myRetList[["raw__FP"]]))
	myRetList[["recall__"]] = ((myRetList[["raw__TP"]]) /
                             (myRetList[["raw__TP"]] +
                              myRetList[["raw__FN"]]))
	myRetList[["fscore__"]] = (2 * myRetList[["precision__"]] *
                             myRetList[["recall__"]]) /
                             (myRetList[["precision__"]] +
                             myRetList[["recall__"]])
	myRetList[["fscore__0_5"]] = ((1 + 1 / 4) * myRetList[["precision__"]] *
                                myRetList[["recall__"]]) / ((1 / 4) *
                                myRetList[["precision__"]] +
                                myRetList[["recall__"]])

	return(myRetList)
}


