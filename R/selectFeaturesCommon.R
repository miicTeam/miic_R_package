#*******************************************************************************
# Filename   : selectFeaturesCommon.R           Creation date: 17 October 2024
#
# Description: Common functions to selectFeatures and selectFeaturesPath
#
# Author     : Franck SIMON
#*******************************************************************************

#===============================================================================
# FUNCTIONS (internal)
#===============================================================================
# sf_check_input_data
#-------------------------------------------------------------------------------
# Check input data
# Params:
# - input_data: a data frame or a matrix (including sparse matrix)
#   Expected layout is samples as rows and variables as columns.
#   Column names must contain the names of the variables.
#-------------------------------------------------------------------------------
sf_check_input_data <- function (input_data)
  {
  if (  ( ! is.data.frame (input_data) )
     && ( ! is.matrix(input_data) )
     && ( ! inherits(input_data, "Matrix") ) )
    miic_error ("parameters",
      "the input data must be a data frame or a matrix.")
  if ( (ncol (input_data) <= 0) || (nrow (input_data) <= 0) )
    miic_error ("parameters", "the input data is empty.")
  if ( is.data.frame (input_data) )
    # Ensure we have a true data frame, e.g. not a tibble
    # (but let matrices unchanged to avoid warnings on large memory allocation)
    # TODO evaluate run time impact on very large data frames
    input_data <- as.data.frame (input_data)
  if ( is.null (colnames (input_data) ) )
    miic_error ("parameters", "the input data must have column names.")
  return (input_data)
  }

