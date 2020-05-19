#*****************************************************************************
# Filename   : tmiic_datasets.R                   Creation date: 24 march 2020
#
# Description: time series data sets generation
#
# Author     : Franck SIMON (fsimon.informaticien@wanadoo.fr)
#
# Changes history:
# - 24 march 2020 : initial version
#*****************************************************************************
#
# 5 predefined models are defined :
#
# Model 0 (debug)
#
# X1(t) = sample * 10000 + node * 1000 + timestep
# X2(t) = sample * 10000 + node * 1000 + timestep
# X3(t) = sample * 10000 + node * 1000 + timestep
# X4(t) = sample * 10000 + node * 1000 + timestep
# X5(t) = sample * 10000 + node * 1000 + timestep
#
# Model 1
#
# X1(t) = 0.20 * X1(t-1) + noise1(t)
# X2(t) = 0.60 * X2(t-1) - 0.287 * f( X5(t-1) ) + 0.287 * f( X4(t-1) ) + noise2(t)
# X3(t) = noise3(t) 
# X4(t) = 0.40 * X4(t-1) - 0.287 * f( X1(t-1) ) + 0.287 * f( X5(t-1) ) + noise4(t)
# X5(t) = 0.60 * X5(t-1) - 0.287 * f( X2(t-1) ) + noise5(t)
#
# Model 2
#
# X1(t) <- 0.40 * X1(t-1) + 0.287 * f( X5(t-1) ) + 0.287 * f( X4(t-1) ) + noise1[t]
# X2(t) <- 0.20 * X2(t-1) - 0.287 * f( X4(t-1) ) + noise2[t]
# X3(t) <- noise3[t]
# X4(t) <- 0.80 * X4(t-1) + 0.287 * f( X5(t-1) ) - 0.287 * f( X3(t-2) ) + noise4[t]
# X5(t) <- 0.60 * X5(t-1)  + noise5[t]
#
# Model 3
#
# X1(t) <- - 0.287 * f( X4(t-1) ) + noise1[t]
# X2(t) <-   0.40 * X2(t-1) - 0.287 * f( X4(t-1) ) + 0.287 * f( X1(t-1) ) +  noise2[t]
# X3(t) <-   0.90 * X3(t-1) - 0.287 * f( X1(t-2) ) + noise3[t]
# X4(t) <-   0.90 * X4(t-1) + noise4[t]
# X5(t) <-   0.80 * X5(t-1) - 0.287 * f( X3(t-2) )  + noise5[t]
#
# Model 4
#
# X1(t) <- 0.40 * X1(t-1) + noise1(t)
# X2(t) <- 0.20 * X2(t-1) - 0.287 * f( X1(t-1) ) + noise2(t)
# X3(t) <- noise3(t)
# X4(t) <- 0.80 * X4(t-1) + 0.287 * f( X1(t-1) ) - 0.287 * f( X2(t-1) ) + noise4(t)
# X5(t) <- 0.60 * X5(t-1) + 0.287 * f( X1(t-2) ) - 0.287 * f( X4(t) ) + noise5(t)
#
# Model 5 : test isolated edges orientation using time
#
# X1(t) <- 0.287 * X2(t-1) + noise1(t)
# X2(t) <- noise2(t)
# X3(t) <- 0.287 * X3(t-1) + noise3(t)
#
# Model 6 : orientation of non lagged edges with 4 points
#
# X1(t) <- noise1(t)
# X2(t) <- - 0.287 * X1(t-1) + noise2(t)
# X3(t) <- + 0.287 * X2(t)   - 0.287 * X4(t-1) + noise3(t)
# X4(t) <- noise4(t)
#
#
# Model 7 : latent variable on X3 (erased by new white noise after generation)
#
# X1(t) <- 0.40 * X1(t-1) + noise1[t]
# X2(t) <- 0.30 * X2(t-1) - 0.287 * f( X1(t-1) ) + 0.287 * f( X3(t-1) ) + noise2[t]
# X3(t) <- noise3[t]
# X4(t) <- 0.20 * X4(t-1) + 0.287 * f( X5(t-1) ) - 0.287 * f( X3(t-1) ) + noise4[t]
# X5(t) <- 0.60 * X5(t-1)  + noise5[t]
#
# noises ~ N(0, 1)
# f1 (x) = x
# f2 (x) = ( 1 - 4 * e[-(x^2)/2] ) * x
# f3 (x) = ( 1 - 4 * x^3 * e[-(x^2)/2] ) * x
#
#
# True graphs
#
tmiic.DF_TRUE_EDGES_MODEL_1 = data.frame (
  orig     = c(1   , 2   , 4   , 5   , 1     , 2     , 5     , 4     , 5     ), 
  dest     = c(1   , 2   , 4   , 5   , 4     , 5     , 2     , 2     , 4     ),
  lag      = c(1   , 1   , 1   , 1   , 1     , 1     , 1     , 1     , 1     ),
  strength = c(+0.2, +0.6, +0.4, +0.6, -0.287, -0.287, -0.287, +0.287, +0.287)  ,
  stringsAsFactors=FALSE)

tmiic.DF_TRUE_EDGES_MODEL_2 = data.frame (
  orig     = c(1   , 2   , 4   , 5   , 3     , 4     , 4     , 5     , 5    ), 
  dest     = c(1   , 2   , 4   , 5   , 4     , 1     , 2     , 1     , 4    ),
  lag      = c(1   , 1   , 1   , 1   , 2     , 1     , 1     , 1     , 1    ),
  strength = c(+0.4, +0.2, +0.8, +0.6, -0.287, +0.287, -0.287, +0.287,+0.287),
  stringsAsFactors=FALSE)

tmiic.DF_TRUE_EDGES_MODEL_3 = data.frame (
  orig     = c(2,    3,    4,    5,    1     , 1     , 3     , 4     , 4     ), 
  dest     = c(2,    3,    4,    5,    2     , 3     , 5     , 1     , 2     ),
  lag      = c(1   , 1   , 1   , 1   , 1     , 2     , 2     , 1     , 1     ),
  strength = c(+0.4, +0.9, +0.9, +0.8, +0.287, -0.287, -0.287, -0.287, -0.287),
  stringsAsFactors=FALSE)

tmiic.DF_TRUE_EDGES_MODEL_4 = data.frame (
  orig     = c(1   , 2,    4,    5,    1  ,    1  ,    2     , 1     , 4     ), 
  dest     = c(1   , 2,    4,    5,    2  ,    4  ,    4     , 5     , 5     ),
  lag      = c(1   , 1   , 1   , 1   , 1     , 1     , 1     , 2     , 0     ),
  strength = c(+0.4, +0.2, +0.8, +0.6, -0.287, +0.287, -0.287, +0.287, -0.287),
  stringsAsFactors=FALSE)

tmiic.DF_TRUE_EDGES_MODEL_5 = data.frame (
  orig     = c(2     , 3     ), 
  dest     = c(1     , 3     ),
  lag      = c(1     , 1     ),
  strength = c(+0.287, +0.287),
  stringsAsFactors=FALSE)

tmiic.DF_TRUE_EDGES_MODEL_6 = data.frame (
  orig     = c(1     , 2     , 4     ), 
  dest     = c(2     , 3     , 3     ),
  lag      = c(1     , 0     , 1     ),
  strength = c(-0.287, +0.287, -0.287),
  stringsAsFactors=FALSE)

tmiic.DF_TRUE_EDGES_MODEL_7 = data.frame (
  orig     = c(1   , 2   , 1     , 3     , 4   , 3     , 5     , 5   ), 
  dest     = c(1   , 2   , 2     , 2     , 4   , 4     , 4     , 5   ),
  lag      = c(1   , 1   , 1     , 1     , 1   , 1     , 1     , 1   ),
  strength = c(+0.4, +0.3, -0.287, +0.287, +0.2, -0.287, +0.287, +0.6),
  stringsAsFactors=FALSE)

#-----------------------------------------------------------------------------
# tmiic.f1
#-----------------------------------------------------------------------------
#' tmiic.f1
#'
#' @description 
#' basic linear function used for sample generation f(x) = x
#'
#' @param x [a number]
#'
#' @return the number received as argument 
#-----------------------------------------------------------------------------
tmiic.f1 <- function (x) 
  {
  return (x)
  }

#-----------------------------------------------------------------------------
# tmiic.f2
#-----------------------------------------------------------------------------
#' tmiic.f2
#'
#' @description 
#' non linear function used for sample generation 
#' f(x) = (1 - 4 * exp[-(x^2)/2]) * x 
#'
#' @param x [a number]
#'
#' @return the number computed by f(x)
#-----------------------------------------------------------------------------
tmiic.f2 <- function (x) 
  {
  return ( ( 1 - 4 * exp(-(x^2)/2) ) * x )
  }

#-----------------------------------------------------------------------------
# tmiic.f3
#-----------------------------------------------------------------------------
#' tmiic.f3
#'
#' @description 
#' non linear function used for sample generation 
#' f(x) = ( 1 - 4 * x^3 * exp[-(x^2)/2] ) * x
#'
#' @param x [a number]
#'
#' @return the number computed by f(x)
#-----------------------------------------------------------------------------
tmiic.f3 <- function (x) 
  {
  return ( ( 1 - 4 * x^3 * exp(-(x^2)/2) ) * x )
  }

#-----------------------------------------------------------------------------
# Constant to iterate over possible models and functions
#-----------------------------------------------------------------------------
tmiic.LIST_MODELS = list (model1 = tmiic.DF_TRUE_EDGES_MODEL_1, 
                          model2 = tmiic.DF_TRUE_EDGES_MODEL_2, 
                          model3 = tmiic.DF_TRUE_EDGES_MODEL_3, 
                          model4 = tmiic.DF_TRUE_EDGES_MODEL_4, 
                          model5 = tmiic.DF_TRUE_EDGES_MODEL_5, 
                          model6 = tmiic.DF_TRUE_EDGES_MODEL_6, 
                          model7 = tmiic.DF_TRUE_EDGES_MODEL_7)
tmiic.LIST_FUNCTS = list (f1 = tmiic.f1, f2 = tmiic.f2, f3 = tmiic.f3)

#-----------------------------------------------------------------------------
# tmiic.call_funct
#-----------------------------------------------------------------------------
#' tmiic.call_funct
#'
#' @description 
#' Utility function to call function f1, f2 or f3 supplied as parameter
#'
#' @param funct [the function to use for the generation : f1, f2 or f3]
#'
#' @return the number computed by f(x)
#-----------------------------------------------------------------------------
tmiic.call_funct <- function(funct, x)
  { 
  funct(x)
  }

#-----------------------------------------------------------------------------
# tmiic.generate_predefined_dataset
#-----------------------------------------------------------------------------
#' tmiic.generate_predefined_dataset
#'
#' @description 
#' Generate n_samples samples for x nodes over n_time timesteps
#' using a predefined model
#'
#' @param model_idx [a number betwen 0 and 4] The model for the data generation.
#' Model 0 is for debug (intialazing nodes with easily identiable values)
#' Model 1 uses only t-1 information. 
#' Models 2 and 3 uses both t-1 and t-2.
#' Model 4 used t-1, t-2 and contemporanous edge
#' @param funct [a function] the function usse for the generation : 
  #' tmiic.f1, f2 or f3 or any other function computing a float from a float 
#' @param n_samples [the number of samples to generate]
#' @param n_time [the number of timesteps of the time series generated]
#'
#' @return a list with two entries :
#' - samples as an array of dimensions n_samples * n_nodes * n_time
#' - a dataframe with the true edes of the network
#-----------------------------------------------------------------------------
tmiic.generate_predefined_dataset <- function (model_idx, funct, n_samples, n_time) 
  {
  DEBUG <- FALSE
  if (DEBUG)
    {
    print ("tmiic.generate_predefined_dataset:")
    print (paste ("model_idx=", model_idx, sep="") )
    print (paste ("n_samples=", n_samples, sep="") )
    print (paste ("n_time=", n_time, sep="") )
    }
  
  if ( (model_idx < 0) | (model_idx > length(tmiic.LIST_MODELS)) )
    {
    stop (paste ("Predefined model", model_idx, "doesn't exist") )
    }
  #
  # All the predefined models have 5 nodes
  #
  list_nodes = list ("X1", "X2", "X3", "X4", "X5")
  n_nodes = length(list_nodes)
  #
  # Model 0 is debug
  #
  if (model_idx == 0)
    {
    true_edges <- data.frame ( orig=character(), dest=character(),
                               lag=integer(), strength=double(), 
                               stringsAsFactors = FALSE)
    data_tab <- array ( data=NA, dim=c(n_samples, n_nodes, n_time),
                        dimnames=list(seq(1,n_samples), list_nodes, seq(1,n_time )) )
    print (dim(data_tab))
    for(sample_idx in 1:n_samples)
      {
      for(node_idx in 1:n_nodes)
        {
        for(time_idx in 1:n_time)
          {
          data_tab [sample_idx, node_idx, time_idx] <- sample_idx * 10000 + node_idx * 1000 + time_idx
          }
        }
      }
    }
  else
    {
    #
    # Model >= 1 are predefined models
    #
    true_edges <- tmiic.LIST_MODELS[[model_idx]]
    data_tab <- tmiic.generate_dataset (true_edges, funct, list_nodes, n_samples, n_time) 
    #
    # For model 7, erase node 3 data which is the latent node and modify the 
    # true edges dataframe to remove edges starting from node3 and add a latent
    # edge between node 2 and 4 (indicated by strengh = 0)
    #
    if (model_idx == 7)
      {
      data_tab[,3,] <- rnorm(n_samples * n_time)
      true_edges <- true_edges[true_edges$orig != 3,]
      true_edges [nrow(true_edges) + 1,] = c(2, 4, 0, 0)
      }
    }
  if (DEBUG)
    {
    print ("")
    print ("Returned true edges:")
    print (true_edges)
    print ("Returned array dimension:")
    print (dim(data_tab))
    print ("Returned data for first sample=1, 5 first timestep:")
    print (data_tab[1,,1:5])
    print (paste ("Returned array on last sample=", n_samples, ", 5 last timesteps", sep="") )
    print (data_tab[n_samples,,(n_time-5):n_time])
    }
  return ( list(data_tab, true_edges) )
  }

#-----------------------------------------------------------------------------
# tmiic.generate_dataset
#-----------------------------------------------------------------------------
#' tmiic.generate_dataset
#'
#' @description 
#' Generate n_samples samples for n_nodes nodes over n_time timesteps
#'
#' @details 
#' The function uses the funct and the true_edge parameters to generate temporal 
#' samples. ie: if in the true edges, the node1 is connected to itself at t-1
#' with strengh 0.5 and to the node2 at t-2 with strength 0.2, the generation 
#' will compute for each timestep of node 1 :
#' node1(t) <- 0.5 * funct ( node1[t-1] ) + 0.2 * funct ( node2[t-2] ) + white noise
#' generate_temporal_dataset discards the first iterations so the tempral data
#' generated are not affected by the intial values.
#' 
#' @param true_edges [a dataframe] The true_edges dataframe must at least contains
#' the columns orig, dest, lag and strengh. The orig and dest columns are the 
#' indexes of the edges's nodes. Lag is a postive or 0 integer. Strengh is 
#' a real number positive or negative indicating the impact of the orig node to
#' the dest one.
#' CAUTION: if a 0 lag is used, the order of the edges in the imput dataframe
#' is important. The values of orig node(s) of the lag 0 edge must be completly
#' computed before the dest node of the lag 0 edge is computed. 
#' i.e. : if we have edges X3(t-1) -> X2(t) and X2(t) -> X1(t)
#' the edge X3(t-1) -> X2(t) must appear first in the true_edges parameter.
#' @param funct [the function to apply when an edge exist between two nodes]
#' @param list_nodes [the list of nodes to generate]
#' @param n_samples [the number of samples to generate]
#' @param n_time [the number of timesteps of the time series generated]
#'
#' @return an array of dimensions n_samples * n_nodes * n_time
#-----------------------------------------------------------------------------
tmiic.generate_dataset <- function (true_edges, funct, list_nodes, n_samples, n_time) 
  {
  DEBUG <- FALSE
  if (DEBUG)
    {
    print ("tmiic.generate_dataset:")
    print ("true edges=")
    print (true_edges)
    print ("list_nodes=")
    print (list_nodes)
    print (paste ("n_samples=", n_samples, sep="") )
    print (paste ("n_time=", n_time, sep="") )
    }
  n_nodes = length (list_nodes)
  max_lag = max(true_edges$lag)
  n_time_generation = (n_time + max_lag) * 2
  #
  # Init the array with white noise, dimensions n_samples * n_nodes * n_time_generation
  #
  data_tab <- array ( data=rnorm(n_samples * n_nodes * n_time_generation),
                      dim=c(n_samples, n_nodes, n_time_generation),
                      dimnames=list ( seq(1,n_samples), list_nodes, seq(1,n_time_generation) ) )
  #
  # Generate n_samples time series
  #
  for(sample_idx in 1:n_samples)
    {
    #
    # Generate the time serie
    #
    for (time_idx in (max_lag+1):n_time_generation)
      {
      #
      # For each true edge, compute new value of nodes
      #
      if ( (DEBUG) & (sample_idx == 1) & (time_idx == n_time_generation - n_time + 1) )
        {
        print (paste ("time_idx=", time_idx, sep="") ) 
        print ("t-2=") 
        print (data_tab[sample_idx, , time_idx-2] ) 
        print ("t-1=") 
        print (data_tab[sample_idx, , time_idx-1] ) 
        print ("old val (noise) at t=") 
        print (data_tab[sample_idx, , time_idx] ) 
        }
      
      for (edge_idx in 1:nrow(true_edges) )
        {
        one_edge <- true_edges[edge_idx,]
        #
        # Function is applied only when edge is from a different node
        #
        if (one_edge$orig == one_edge$dest)
          {
          fct_res <- data_tab[sample_idx, one_edge$orig, time_idx - one_edge$lag]
          }
        else
          {
          fct_res <- tmiic.call_funct (funct, data_tab[sample_idx, one_edge$orig, time_idx - one_edge$lag])
          }
        
        if ( (DEBUG) & (sample_idx == 1) & (time_idx == n_time_generation - n_time + 1) )
          {
          print ("edge=") 
          print (one_edge) 
          print (paste ("old val of dest node ", list_nodes[[one_edge$dest]], "=",
                        data_tab[sample_idx, one_edge$dest, time_idx], sep="") )
          print (paste ("use val of orig node ", list_nodes[[one_edge$orig]], 
                        " lag ", one_edge$lag, "=",
                        data_tab[sample_idx, one_edge$orig, time_idx - one_edge$lag], sep="") )
          print (paste ("function result=", fct_res, sep="") )
          print (paste ("function result * strength=", one_edge$strength * fct_res, sep="") )
          }
        
        data_tab[sample_idx, one_edge$dest, time_idx] <- data_tab[sample_idx, one_edge$dest, time_idx] + one_edge$strength * fct_res

        if ( (DEBUG) & (sample_idx == 1) & (time_idx == n_time_generation - n_time + 1) )
          {
          print (paste ("new val at t of dest node ", list_nodes[[one_edge$dest]], "=",
                        data_tab[sample_idx, one_edge$dest, time_idx], sep="") ) 
          }
        }

      if ( (DEBUG) & (sample_idx == 1) & (time_idx == n_time_generation - n_time + 1) )
        {
        print ("new val at t=") 
        print (data_tab[sample_idx, , time_idx] ) 
        }
      }
    }
  #
  # Keep only the last n_time timesteps so the tempral data
  # generated are not affected by the intial values.
  #
  data_tab <- data_tab [, , (n_time_generation-n_time+1):n_time_generation]
  dim_data <- dim(data_tab)
  if ( length (dim_data) <= 2)
    {
    #
    # CAUTION : If we have only one sample, the dimension is reduced to 2
    # => change back to 3
    #
    tmp_dim_names <- dimnames(data_tab)
    dim(data_tab) <- c(1, dim_data)
    dimnames(data_tab) <- c(c(1), tmp_dim_names) 
    dim_data <- dim(data_tab)
    }
  dimnames(data_tab)[[3]] <- seq(1, n_time)
  
  if (DEBUG)
    {
    print ("Returned array dimension")
    print ( dim(data_tab) )
    print ("Returned array on first sample=1, 5 first timesteps")
    print (data_tab [1,,1:5])
    print (paste ("Returned array on last sample=", n_samples, ", 5 last timesteps", sep="") )
    print (data_tab [n_samples,,(n_time-5):n_time])
    }
  return (data_tab)
  }


#-----------------------------------------------------------------------------
# tmiic.plot_one_sample
#-----------------------------------------------------------------------------
#' tmiic.plot_one_sample
#'
#' @description 
#' plot one sample of a temporal serie
#'
#' @param data_tab [an array of dimensions n_samples * n_nodes * n_time]
#' @param sample_idx [an integer] the sample index to plot
#' @param title [a string] Optional, NULL by default. The title of the plot
#' @param filename [a string] Optional, NULL by default. If supplied, the plot
#' is saved in this file
#'
#' @return None
#-----------------------------------------------------------------------------
tmiic.plot_one_sample <- function (data_tab, sample_idx, title=NULL, filename=NULL)
  {
  DEBUG <- FALSE
  
  if (! is.null(filename) )
    {
    png (filename=filename)
    }
  data_plot <- data_tab[sample_idx,,]
  data_plot <- t(data_plot)
  df_plot <- as.data.frame ( ts(data_plot) )
  colnames(df_plot) <- dimnames(data_tab)[[2]]
  
  if (DEBUG)
    {
    print ("head df to plot")
    print (head (df_plot))
    }
  plot.ts (df_plot)
  if (! is.null(title) )
    {
    title (title)
    }
  if (! is.null(filename) )
    {
    dev.off()
    }
  }

