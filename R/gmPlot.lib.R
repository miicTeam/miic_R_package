plot.loadSummary <- function( mySummary )
{
    #### Load the summary of the edges
    rownames(mySummary) = c()

    #### Ignore the TN edges
    myTypesToIgnore = c()
    myTypesToIgnore = c("TN","N","FN")
    myLinesToIgnore = which( mySummary[,"type"] %in% myTypesToIgnore )
    if( length( myLinesToIgnore ) > 0 ) { mySummary = mySummary[-myLinesToIgnore,]}

   return(mySummary)
}

plot.createDefaultGraph <- function( mySummary, myAllGenes )
{
    #### Replace names by numbers to create the graph
    inf.edgesList.nbr <- apply( mySummary[, c("x","y")], MARGIN = c(1,2), function(x) { x <- which( myAllGenes == x ) } )

    #### Create an unoriented igraph with all the nodes
    inf.graph <- igraph::graph( t( inf.edgesList.nbr ), length( myAllGenes ), directed = TRUE )

    #### Set the vertices options
    igraph::V(inf.graph)$label <- myAllGenes
    igraph::V(inf.graph)$shape <- "circle"
    igraph::V(inf.graph)$color <- "lightblue"
    igraph::V(inf.graph)$label.family <- "Helvetica"
    igraph::V(inf.graph)$label.cex <- 0.6
    igraph::V(inf.graph)$size <- 10

    #### Set the general edges options
    igraph::E(inf.graph)$arrow.size <- 0.5
    igraph::E(inf.graph)$arrow.width <- 3
    igraph::E(inf.graph)$width <- 3
    igraph::E(inf.graph)$curved <- FALSE
    igraph::E(inf.graph)$color <- "red2"
    igraph::E(inf.graph)$lty <- "solid"
    igraph::E(inf.graph)$arrow.mode <- 0
    return(inf.graph)
}

plot.setOrientation <- function( mySummary, inf.graph )
{
    #### Set the options specific for forward oriented
    ort.fwd.idx <- which( ( mySummary[, "infOrt"] %in% c(2,4) ) | ( mySummary[, "type"] == 'FN' & mySummary[, "trueOrt"] == 2 ) )
    if( length( ort.fwd.idx ) > 0 ){ igraph::E(inf.graph)[ort.fwd.idx]$arrow.mode <- 2 }

    #### Set the options specific for backward oriented
    ort.bck.idx <- which( ( mySummary[, "infOrt"] %in% c(-2,-4) ) | ( mySummary[, "type"] == 'FN' & mySummary[, "trueOrt"] == (-2) ) )
    if( length( ort.bck.idx ) > 0 ) { igraph::E(inf.graph)[ort.bck.idx]$arrow.mode <- 1 }

    #### Set the options specific for bidirectional orientations
    bidir.idx <- which( mySummary[, "infOrt"] == 6 )
    if( length( bidir.idx ) > 0 )
    {
        igraph::E(inf.graph)[bidir.idx]$arrow.mode <- 3
    }
    return(inf.graph)

}

littlefunc <- function(edge1,edge2)
{
    return(paste(edge1, collapse=",") == paste(rev(edge2), collapse=","))
}

# ---- Function to plot graphes with edges matching to their partial correlation
pCor.edgeCol <- function(summary, features)
{
    # Define the color gradients
    blue.gradient = grDevices::rainbow(100, start = 3/6, end=4/6)
    red.gradient = rev(grDevices::rainbow(100, start=0, end=0.16))

    myEdgesColor = rep(NA, nrow(summary)) # Set the color vector for the edges
    max.pcor.neg = suppressWarnings(min(summary[which(summary[,"sign"] == "-"),"partial_correlation"])) # get the maximum negative pcor
    max.pcor.pos = suppressWarnings(max(summary[which(summary[,"sign"] == "+"),"partial_correlation"])) # get the maximum positive pcor
    for(edge in 1:nrow(summary)) # loop on all the edges present in the network
    {
        # Set the correct tmp.max.pcor
        if(sign(summary[edge, "partial_correlation"]) == -1){tmp.max.pcor = max.pcor.neg}
        else {tmp.max.pcor = max.pcor.pos}
        # Compute the ratio between the tmp.max.pcor and the edge's pcor, and use it as an index to get a color
        edge.pCor.ind = abs(summary[edge, "partial_correlation"])
        edge.colIndex = round(edge.pCor.ind * 100)
        if( edge.colIndex == 0 ) { edge.colIndex = 1}
        ### Get the sign of the link to look at the correct color gradient
        if(! is.na(summary[edge, "sign"]) )
        {
            if(summary[edge, "sign"] == "+")
            {
                myEdgesColor[edge] = red.gradient[edge.colIndex]
            }
            else
            {
                myEdgesColor[edge] = blue.gradient[edge.colIndex]
            }
        }
        else { myEdgesColor[edge] = "grey88" }
    }

    return(myEdgesColor)
}

# ---- Function to plot graphes with edges matching to their mutual information (confidence column in summary)
conf.edgeCol <- function(summary, features)
{
    # Define the color gradients
    blue.gradient = grDevices::rainbow(100, start = 3/6, end=4/6)
    red.gradient =  rev(grDevices::rainbow(100, start=0, end=0.16))
    myEdgesColor = rep(NA, nrow(summary)) # Set the color vector for the edges
    max.conf = 100 # get the maximum
    min.conf = 1 # get the minimum

    for(edge in 1:nrow(summary)) # loop on all the edges present in the network
    {
        # Set the correct tmp.max.pcor
        # Compute the ratio between the tmp.max.pcor and the edge's pcor, and use it as an index to get a color
        edge.colIndex = round(summary[edge, "log_confidence"])
        if( edge.colIndex < min.conf  ) { edge.colIndex = 1 }
        else if( edge.colIndex > max.conf ){ edge.colIndex = 100 }
        ### Get the sign of the link to look at the correct color gradient
        if(! is.na(summary[edge, "sign"]) )
        {
            if(summary[edge, "sign"] == "+")
            {
                myEdgesColor[edge] = red.gradient[edge.colIndex]
            }
            else
            {
                myEdgesColor[edge] = blue.gradient[edge.colIndex]
            }
        }
        else { myEdgesColor[edge] = red.gradient[edge.colIndex] }
    }

    return(myEdgesColor)
}

# ---- Function which returns a graph object from edges colors and node sizes eventually
modif.Graph <- function(summary, features, edgeColors, nodeSizes = 10, nodeColors = 'lightblue')
{
    mygraph = plot.createDefaultGraph(summary, features)
    igraph::E(mygraph)$color = edgeColors
    igraph::E(mygraph)$arrow.size = 0.2
    igraph::V(mygraph)$color = nodeColors
    mygraph = plot.setOrientation(summary, mygraph)
    igraph::V(mygraph)$size = nodeSizes
    return(mygraph)
}
