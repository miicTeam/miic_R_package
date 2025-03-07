#*******************************************************************************
# Filename   : selectFeatures.R                 Creation date: 17 October 2024
#
# Description: Features selection based on Mutual Information
#
# Author     : Franck SIMON
#*******************************************************************************

#===============================================================================
# FUNCTIONS (internal)
#===============================================================================
# plot_top_features
#-------------------------------------------------------------------------------
plot_top_features = function (df_mis, var_of_interest_name, n_plots=25,
  x_lab=NULL, y_lab=NULL, values=T, annotate=T, font_size=11,
  box_fill="steelblue", box_text="white")
  {
  # Check the parameters that can have been supplied by the user
  # (passed by the ... extra params of selectFeatures function)
  #
  n_plots <- miic:::check_param_int (n_plots, "number of features to plot",
                                     default=25, min=1)
  if ( is.null (x_lab) )
    x_lab <- paste0 ("Top features for ", var_of_interest_name)
  else
    x_lab <- as.character (x_lab)
  if ( is.null (y_lab) )
    y_lab <- "MI (bits)"
  else
    y_lab <- as.character (y_lab)
  values <- miic:::check_param_logical (values, "plotting of values", default=T)
  annotate <- miic:::check_param_logical (values, "plotting of annotatation", default=T)
  font_size <- miic:::check_param_int (font_size, "font size",
                                       default=11, min=1)
  # TODO: add a function for color checking
  if ( is.null (box_fill) )
    box_fill <- "steelblue"
  if ( is.null (box_text) )
    box_text <- "white"
  #
  # Extract and order desc the MIs for the var of interest to plot
  #
  v_mis = df_mis[, var_of_interest_name, drop=T]
  v_mis = v_mis[ (!is.na (v_mis)) & (v_mis > 0) ]
  v_mis = sort (v_mis, decreasing=T)
  df = data.frame ("Features" = names(v_mis), "MI" = v_mis)
  if (nrow (df) > n_plots)
    df = df[1:n_plots, , drop=F]
  feat_remaining = length(v_mis) - n_plots
  mi_min = min(df[,"MI"])

  p <- ggplot(data=df, aes(x=Features, y=MI)) +
    geom_bar (stat="identity", fill=box_fill) +
    scale_x_discrete (limits = df$Features ) +
    theme_classic() +
    theme ( text=element_text(size=font_size) ) +
    theme ( axis.text.x=element_text (angle=30, vjust =1, hjust=1) ) +
    xlab (x_lab) +
    ylab (y_lab)

  if (values)
    {
    if (mi_min < 1)
      p = p + geom_text(aes(label=round(MI, 2)), vjust=1.6, color=box_text, size=.pt * 0.75)
    else if (mi_min < 10)
      p = p + geom_text(aes(label=round(MI, 1)), vjust=1.6, color=box_text, size=.pt * 0.75)
    else
      p = p + geom_text(aes(label=round(MI, 0)), vjust=1.6, color=box_text, size=.pt * 0.75)
    }

  if ( annotate && (feat_remaining > 0) )
    {
    annotation <- data.frame (x = nrow(df), y = max(df$MI),
       label = paste0 ("... ", feat_remaining , " more with MI > 0") )
    p = p + geom_text (data=annotation, aes( x=x, y=y, label=label), size=.pt,
                       hjust=1, vjust=1)
    }
  return (p)
  }

#===============================================================================
# FUNCTIONS (exported)
#===============================================================================
# selectFeatures
#-------------------------------------------------------------------------------
#' Features selection based on Mutual Information (MI)
#'
#' @description Select features by selecting the set of variables having
#' the highest MIs with the variable(s) of interest. When multiple variables
#' of interest are supplied, the function will try to select an equal
#' number of features per variable of interest.
#'
#' Variables not sharing information at all (MI = 0) can not be selected
#' as features, so the number of features returned globally or for a variable
#' of interest can be lower than requested. This can also lead to
#' an unbalanced number of features per variable of interest.
#'
#' @param input_data [a data frame or a matrix, required]
#'
#' Expected layout is samples as rows and variables as columns. Column names
#' must contain the names of the variables.
#'
#' @param n_features [an integer, required]
#'
#' The total number of features to select, must be in the range [0, number
#' of variables in \emph{input_data}]. If several variables of interest
#' are supplied, the function will try to select an equal number of features
#' per variable.
#'
#' When \emph{n_features} is 0, no feature selection will be performed
#' but the MIs between the variables in \emph{input_data} and the variable(s)
#' of interest will be computed and returned.
#' This can be used to compute once for all the MIs and perform several
#' features selection at a later stage, supplying MIs matrix in the
#' \emph{precomputed_mis} parameter.
#'
#' @param var_of_interest_names [a string or vector of strings, optional,
#' NULL by default]
#'
#' For the variable(s) of interest that are part of the \emph{input_data},
#' you should supply their names here.
#'
#' @param var_of_interest_values [a data frame, optional, NULL by default]
#'
#' For the variables of interest that are not in \emph{input_data},
#' a data frame can be supplied. The column names are the names
#' of the variables of interest and rows are the samples
#' ordered in the same way as the \emph{input_data}.
#' Typically, such variables are metadata associated to samples but not
#' stored in \emph{input_data}, e.g. a "Treatment" vs "Control" variable
#' in an experiment and a count matrix with the expression of genes
#' in \emph{input_data}.
#'
#' @param skip_cheks [a boolean, optional, FALSE by default]
#'
#' Before computing MI between the variable of interest and the features,
#' \emph{input_data} is checked to filter out constant features and rows full
#' of NAs. When the \emph{input_data} does not need such filtering,
#' these checks can be skipped to speed up the process.
#'
#' @param precomputed_mis [a matrix, optional, NULL by default]
#'
#' if MIs has been previously computed between some variables
#' in the \emph{input_data} and variable(s) of interest,
#' supplying these precomputed MIs will speed up the process as the existing
#' MIs (values different from NA) will not be recomputed.
#' This matrix must have variables names from the \emph{input_data}
#' as row names and variables of interest names as column names
#' (the layout is the same as the \emph{mis} matrix returned).
#'
#' @param n_threads [a positive integer, optional, 1 by default]
#'
#' When set greater than 1, \emph{n_threads} parallel threads will be used for
#' computation. Make sure your compiler is compatible with openmp
#' if you wish to use multithreading.
#'
#' @param verbose [an integer, optional, 3 by default]
#'
#' Level of verbosity: 0=no display, 1=summary, 2=progress per variable of
#' interest, 3=same as 2 with display of estimated time remaining.
#'
#' @param plot [a boolean, optional, FALSE by default]
#'
#' If set to TRUE, a plot with the top features for each variable of
#' interest will be generated.
#'
#' Extra parameters can be used to customize the plot rendering:
#'
#' \itemize{
#' \item \emph{n_plots}: an integer, optional, 25 by default, the number of
#'   features to plot.
#' \item \emph{x_lab}: a string, optional, "Top features for
#'   variable_of_interest_name" by default, the X axis label.
#' \item \emph{y_lab}: a string, optional, "MI (bits)" by default,
#'   the Y axis label.
#' \item \emph{values}: a boolean, optional, TRUE by default,
#'   displays the MI of each feature when activated.
#' \item \emph{annotate}: a boolean, optional, TRUE by default,
#'   displays the remaining number of features with MI > 0 when activated.
#' \item \emph{font_size}: an integer, optional, 11 by default, the font size.
#' \item \emph{box_fill}: a string, optional, "steelblue" by default,
#'   the filling color of the bar plot.
#' \item \emph{box_text}: a string, optional, "white" by default,
#'   the text color of MI values in the bar plot.
#' }
#'
#' @return A named list with three items:
#'
#' \itemize{
#' \item \emph{features}: a vector with the features selected.
#' \item \emph{mis}: a matrix with the MIs (in bits) between the variables in
#'   \emph{input_data} (as rows) and the variable(s) of interest (as columns).
#'   Row and column names are sorted alphabetically.\cr
#'   If a precomputed MIs matrix was supplied, new values computed
#'   are added to the existing matrix.
#' \item \emph{plots}: when the \emph{plot} parameter is turned to TRUE,
#'   a list of plots with the top features for each variable of interest.
#'   The list is empty if the \emph{plot} parameter is set to FALSE.
#' }
#'
#' @examples
#' library(miic)
#'
#' # Get the 10 top features related to "tp53" (gene mutation)
#' # and "TP53" (gene expression) in the dataset
#' ret <- selectFeatures (cosmicCancer,
#'                        var_of_interest_names=c("TP53", "tp53"),
#'                        n_features=10)
#' message ("Tops 10 features = ", paste (ret$features, collapse=", ") )
#'
#' # Features selection with reuse of above computed MIs
#' ret <- selectFeatures (cosmicCancer,
#'                        var_of_interest_names=c("TP53", "tp53"),
#'                        n_features=20,
#'                        precomputed_mis=ret$mis)
#' message ("Tops 20 features = ", paste (ret$features, collapse=", ") )
#'
#' # Features selection using an external metadata "Ploidy"
#' # (simulation by extracting "Ploidy" out of the dataset)
#' df_external_meta <- data.frame ("Ploidy" = cosmicCancer$Ploidy)
#' df_data <- cosmicCancer[ , ! (colnames(cosmicCancer) == "Ploidy") ]
#' ret <- selectFeatures (df_data,
#'                        var_of_interest_values=df_external_meta,
#'                        n_features=10)
#' message ("Tops 10 features = ", paste (ret$features, collapse=", ") )
#'
#' # Same features selection with plot of the result
#' if ( require(ggplot2) ) {
#'   ret <- selectFeatures (df_data,
#'                          var_of_interest_values=df_external_meta,
#'                          n_features=10,
#'                          plot=TRUE)
#'   print (ret$plots[[1]])
#' }
#'
#' # Features selection using variables of interest coming both from
#' # the dataset ("tp53", "TP53") and an external metadata ("Ploidy")
#' # (simulation by extracting "Ploidy" out of the dataset)
#' df_external_meta <- data.frame ("Ploidy" = cosmicCancer$Ploidy)
#' df_data <- cosmicCancer[ , ! (colnames(cosmicCancer) == "Ploidy") ]
#' ret <- selectFeatures (df_data,
#'                        var_of_interest_names=c("TP53", "tp53"),
#'                        var_of_interest_values=df_external_meta,
#'                        n_features=30)
#' message ("Tops 30 features = ", paste (ret$features, collapse=", ") )
#'
#' # Same features selection with plotting using a customized rendering
#' if ( require(ggplot2) && require(gridExtra) ) {
#'   ret <- selectFeatures (df_data,
#'                          var_of_interest_names=c("TP53", "tp53"),
#'                          var_of_interest_values=df_external_meta,
#'                          n_features=30,
#'                          plot=TRUE,
#'                          n_plots=10,
#'                          box_fill="lightgreen",
#'                          box_text="black")
#'   do.call("grid.arrange", c(ret$plots, ncol=2))
#' }
#'
#' @export
#-------------------------------------------------------------------------------
selectFeatures <- function (input_data, n_features,
  var_of_interest_names=NULL, var_of_interest_values=NULL,
  skip_cheks=F, precomputed_mis=NULL, n_threads=1, plot=F, verbose=3, ...)
  {
  n_features <- miic:::check_param_int (n_features, "number of features",
                                        default=0, min=0)
  #
  # Check other params and compute MIs
  #
  mat_mis = miic:::compute_mi_batch (input_data=input_data,
    var_of_interest_names=var_of_interest_names,
    var_of_interest_values=var_of_interest_values,
    skip_cheks=skip_cheks, precomputed_mis=precomputed_mis,
    n_threads=n_threads, verbose=verbose)
  #
  # The mat_mis can contain more rows than the features to select
  # e.g. we computed the MI for all genes and now we want select only
  # the TFs. In this case, in input_data, the variables are only the TFs
  # while mat_mis would contain all genes. Same for the columns as
  # we can have precomputed more variables of interest than the ones we use now
  #
  all_voi_names = unique (c (var_of_interest_names,
                             colnames(var_of_interest_values) ) )
  n_vois = length (all_voi_names)
  mat_mis_filt = mat_mis[rownames(mat_mis) %in% colnames(input_data),
                         colnames(mat_mis) %in% all_voi_names,
                         drop=F]
  #
  # Plot if requested
  #
  list_plots = list ()
  if (plot)
    for ( one_col in colnames (mat_mis_filt) )
      list_plots[[one_col]] = miic:::plot_top_features (
        df_mis=mat_mis_filt, var_of_interest_name=one_col, ...)
  #
  # If no feature selection, can end here
  #
  if (n_features <= 0)
    return (list ("features"=c(), "mis"=mat_mis, "plots"=list_plots) )
  #
  # For each voi requested, order mis decreasing and filter on MI' > 0
  #
  list_mis_sorted = list()
  uniq_poss_feats = c()
  for (one_voi in all_voi_names)
    {
    list_mis_sup_0 = mat_mis_filt[, one_voi] [mat_mis_filt[, one_voi] > 0]
    list_mis_sorted[[one_voi]] = sort (list_mis_sup_0, decreasing=T)
    uniq_poss_feats = unique (c (uniq_poss_feats,
                                 names(list_mis_sup_0) ) )
    }
  #
  # Check if too much features to select
  #
  if (n_features >= length(uniq_poss_feats) )
    {
    if (n_features > length(uniq_poss_feats) )
      miic:::miic_warning ("features selection",
        "too few variables share information with the variable(s)",
        " of interest to select the top ", n_features,
        ", returning all the ", length(uniq_poss_feats),
        " variables sharing information.")
    return (list ("features"=uniq_poss_feats, "mis"=mat_mis, "plots"=list_plots) )
    }
  #
  # Get top variables for each voi
  #
  list_tops_prec = c()
  list_tops = c()
  #
  # Select top features until we get enough or a bit too much
  #
  list_tops = c()
  list_tops_prec = c()
  n_tops_to_sel = n_features %/% n_vois
  while (T)
    {
    list_tops = unique (unlist (lapply (list_mis_sorted, FUN=function(x) {
      names(x)[1:n_tops_to_sel] } ) ) )
    n_missing = n_features - length (list_tops)
    if (n_missing <= 0)
      break
    list_tops_prec = list_tops
    n_tops_to_sel = n_tops_to_sel + max (1, n_missing %/% n_vois)
    }
  #
  # If a bit too much, select in the last features added those with highest MI
  #
  if (n_missing < 0)
    {
    list_tops_to_test = unlist (lapply (list_mis_sorted,
      FUN=function(x) { names(x)[n_tops_to_sel] } ) )
    list_mis_to_test = unlist (lapply (list_mis_sorted,
      FUN=function(x) { x[n_tops_to_sel] } ) )
    list_tops_to_test = list_tops_to_test[order (list_mis_to_test,
                                                  decreasing=T)]
    list_tops_to_test = unique (list_tops_to_test)
    list_tops_to_test = list_tops_to_test[
      ! (list_tops_to_test %in% list_tops_prec) ]

    n_missing = n_features - length (list_tops_prec)
    list_tops = c (list_tops_prec, list_tops_to_test[1:n_missing])
    }
  if (verbose >= 1)
    cat (paste0 (length(list_tops), " features selected.\n") )
  return (list ("features"=list_tops, "mis"=mat_mis, "plots"=list_plots) )
  }

