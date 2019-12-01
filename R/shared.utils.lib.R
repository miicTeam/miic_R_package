#### Convert the 'binary' vector corresponding to the properties into a string
#### Exp:
#### Prop: A B C D E F
#### Vect: D A C
#### Bin : 1 0 1 1 0 0
#### Str : "101100"
binVectToStr <- function(myEns, myGv) {
  myBinVect = (myGv$allProperties %in% myEns)
  return(paste(as.numeric(myBinVect), collapse = ''))
}
