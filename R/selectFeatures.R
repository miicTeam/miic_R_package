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
  box_fill="#1F78B4", box_text="white")
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
  annotate <- miic:::check_param_logical (annotate, "plotting of annotatation", default=T)
  font_size <- miic:::check_param_int (font_size, "font size",
                                       default=11, min=1)
  # TODO: add a function for color checking
  if ( is.null (box_fill) )
    box_fill <- "#1F78B4"
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

  if ( annotate && (length(v_mis) > n_plots) )
    {
    annotation <- data.frame (x = nrow(df), y = max(df$MI),
       label = paste0 ("Top ", n_plots, " shown of ",
                       length(v_mis), " features with MI > 0") )
    p = p + geom_text (data=annotation, aes( x=x, y=y, label=label), size=.pt,
                       hjust=1, vjust=1)
    }
  return (p)
  }

#-------------------------------------------------------------------------------
# plot_top_features_dpi
#-------------------------------------------------------------------------------
plot_top_features_dpi = function (dpis, couple, n_plots=20,
  x_lab=NULL, y_lab=NULL, values=T, annotate=T, threshold=T, font_size=11,
  box_fill=c("#A6CEE3", "#1F78B4"), box_text="white")
  {
  # Check the parameters that can have been supplied by the user
  # (passed by the ... extra params of selectFeaturesPath function)
  #
  n_plots <- miic:::check_param_int (n_plots, "number of features to plot",
                                     default=20, min=1)
  couple_name = paste0 (couple[1,1], "-", couple[1,2])
  if ( is.null (x_lab) )
    x_lab <- paste0 ("Top features for ", couple_name)
  else
    x_lab <- as.character (x_lab)
  if ( is.null (y_lab) )
    y_lab <- "MI (bits)"
  else
    y_lab <- as.character (y_lab)
  values <- miic:::check_param_logical (values, "plotting of values", default=T)
  annotate <- miic:::check_param_logical (annotate, "plotting of annotatation", default=T)
  threshold <- miic:::check_param_logical (threshold, "plotting of the MI threshold", default=T)
  font_size <- miic:::check_param_int (font_size, "font size",
                                       default=11, min=1)
  # TODO ? add a function for color checking
  if ( is.null (box_fill) )
    box_fill=c("#A6CEE3", "#1F78B4")
  if ( is.null (box_text) )
    box_text <- "white"
  #
  # Extract and order desc the MIs for the var of interest to plot
  #
  feat_x_2 = unlist (lapply (rownames(dpis), FUN=function(x) { c(x,x) } ) )
  df = data.frame ("Features" = as.factor (feat_x_2),
   "Voi" = as.factor (rep ( c (couple[1,1], couple[1,2]), nrow(dpis)) ),
   "MI" = NA_real_, stringsAsFactors=F)
  df$MI[ df$Voi == couple[1,1] ] = dpis[, couple[1,1] ]
  df$MI[ df$Voi == couple[1,2] ] = dpis[, couple[1,2] ]
  if (nrow (df) > n_plots * 2)
    df = df[1:(n_plots*2), , drop=F]
  mi_min = min(df[,"MI"])

  p <- ggplot(data=df, aes(x=Features, y=MI, fill=Voi)) +
    geom_bar ( stat="identity", position=position_dodge(width=1.8) ) +
    scale_x_discrete (limits = df$Features) +
    scale_fill_manual("Variables\nof\ninterest", values = box_fill) +
    theme_classic() +
    theme ( text=element_text(size=font_size),
            axis.text.x=element_text (angle=30, vjust =1, hjust=1),
            legend.title = element_text (hjust = 0.5) ) +
    xlab (x_lab) +
    ylab (y_lab)

  if (values)
    {
    if (mi_min < 1)
      p = p + geom_text(aes(label=round(MI, 2)), position=position_dodge(width=1.8),
                        vjust=1.6, color=box_text, size=.pt * 0.75)
    else if (mi_min < 10)
      p = p + geom_text(aes(label=round(MI, 1)), position=position_dodge(width=1.8),
                        vjust=1.6, color=box_text, size=.pt * 0.75)
    else
      p = p + geom_text(aes(label=round(MI, 0)), position=position_dodge(width=1.8),
                        vjust=1.6, color=box_text, size=.pt * 0.75)
    }

  if ( annotate && (nrow(dpis) > n_plots) )
    {
    annot_text = paste0 ("Top ", n_plots, " shown of ", nrow(dpis), " features passing the dpi test")
    p <- p + annotate ("text", label = annot_text, x = nrow(df), y = max(df$MI),
                       hjust=1, vjust=1, size=.pt)
    }

  if (threshold)
    {
    df_legend2 = data.frame (
      label= paste0 (round(couple[1,"mi"],4), " bits"),
      value=couple[1,"mi"])
    p <- p + geom_hline (data = df_legend2, linetype="dashed", linewidth=0.6,
                         aes(yintercept = value, color = label) ) +
      scale_color_manual("MI thresold",
                         values="black",
                         guide = guide_legend(override.aes = list(fill = NA)))
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
#' \item \emph{box_fill}: a string, optional, "#1F78B4" by default
#'   (darker blue of the Paired brewer palette), the bar plot filling color.
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
#'   do.call ("grid.arrange", c(ret$plots, ncol=2,
#'     top="Top features on cosmicCancer dataset" ) )
#' }
#'
#' @export
#-------------------------------------------------------------------------------
selectFeatures <- function (input_data, n_features,
  var_of_interest_names=NULL, var_of_interest_values=NULL,
  skip_cheks=F, precomputed_mis=NULL, n_threads=1, verbose=3, plot=F, ...)
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

#-------------------------------------------------------------------------------
# selectFeaturesPath
#-------------------------------------------------------------------------------
#' Features selection on the path between two sets of variables of interest
#' based on Mutual Information (MI)
#'
#' @description Select features by using data propagation inequality (DPI)
#' to select the set of variables that are the most likely to be in the path
#' between two set of variable(s) of interest.
#' When multiple variables of interest are supplied in one or both sides of
#' the path, the function will try to select an equal number of features
#' per couple (one of each side) of variables of interest .
#'
#' Variables not sharing information at all (MI = 0) or not in the path
#' between the set of variables considering the DPI can not be selected
#' as features. So the number of features returned globally or linked to a
#' set of variables of interest can be lower than requested. This can also
#' lead to unbalanced number of features per couple of variables of interest.
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
#' @param var_of_interest_names_side1 [a string or vector of strings, optional,
#' NULL by default]
#'
#' For the variable(s) of interest that are part of the \emph{input_data}
#' on the first side of the path, you should supply their names here.
#'
#' @param var_of_interest_values_side1 [a data frame, optional, NULL by default]
#'
#' For the variables of interest that are not in \emph{input_data}
#' and on the first side of the path, a data frame can be supplied.
#' The column names are the names of the variables of interest
#' and rows are the samples ordered in the same way as the \emph{input_data}.
#' Typically, such variables are metadata associated to samples but not
#' stored in \emph{input_data}, e.g. a "Treatment" vs "Control" variable
#' in an experiment and a count matrix with the expression of genes
#' in \emph{input_data}.
#'
#' @param var_of_interest_names_side2 [a string or vector of strings, optional,
#' NULL by default]
#'
#' Same as \emph{var_of_interest_names_side1} for the second side of the path.
#'
#' @param var_of_interest_values_side2 [a data frame, optional, NULL by default]
#'
#' Same as \emph{var_of_interest_values_side1} for the second side of the path.
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
#' @param filter_v_struct [a boolean, optional, FALSE by default]
#'
#' If set to TRUE, the features passing the DPI test are filtered
#' out if they form a v-structure with the couple of variables of interest.
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
#' If set to TRUE, a plot with the top features for each couple of variables of
#' interest will be generated.
#'
#' Extra parameters can be used to customize the plot rendering:
#'
#' \itemize{
#' \item \emph{n_plots}: an integer, optional, 25 by default, the number of
#'   features to plot.
#' \item \emph{x_lab}: a string, optional, "Top features for
#'   variable_of_interest_1-variable_of_interest_2" by default,
#'   the X axis label.
#' \item \emph{y_lab}: a string, optional, "MI (bits)" by default,
#'   the Y axis label.
#' \item \emph{values}: a boolean, optional, TRUE by default,
#'   displays the MI values of each feature when activated.
#' \item \emph{annotate}: a boolean, optional, TRUE by default,
#'   displays the total number of features passing the DPI test
#'   when activated.
#' \item \emph{threshold}: a boolean, optional, TRUE by default,
#'   displays the MI threshold of the DPI test when activated.
#' \item \emph{font_size}: an integer, optional, 11 by default, the font size.
#' \item \emph{box_fill}: a vector, optional, c("#A6CEE3", "#1F78B4") by
#'   default (paired blues of brewer palette), the bar plot filling colors.
#' \item \emph{box_text}: a string, optional, "white" by default,
#'   the text color of MI values in the bar plot.
#' }
#'
#' @return A named list with six items:
#'
#' \itemize{
#' \item \emph{features}: a vector with the features selected.
#' \item \emph{mis}: a matrix with the MIs (in bits) between the variables in
#'   \emph{input_data} (as rows) and the variable(s) of interest (as columns).
#'   Row and column names are sorted alphabetically.\cr
#'   If a precomputed MIs matrix was supplied, new values computed
#'   are added to the existing matrix.
#' \item \emph{couples}: a data frames with the MIs of each couple of variables
#'   of interest used as threshold for the DPI test.
#' \item \emph{dpis}: a list if data frames with the MIs (in bits) used for the
#'   DPI tests. Each item corresponds to a couple of variables of interest.
#'   The row names of the data frames are the features passing the DPI test
#'   and data frames are ordered by sum of MIs with the variables of interest
#'   descending.
#' \item \emph{plots}: when the \emph{plot} parameter is turned to TRUE,
#'   a list of plots with the top features for each couple of variables of
#'   interest. The list is empty if the \emph{plot} parameter is set to FALSE.
#' \item \emph{ais}: when the \emph{filter_v_struct} parameter is turned to
#'   TRUE, a list with the separating set for each couple of variables of
#'   interest.
#'   The list is empty if the \emph{filter_v_struct} parameter is set to FALSE.
#' \item \emph{v_struct_excluded}: the \emph{filter_v_struct} parameter is
#'   turned to TRUE, a list with the features excluded (because forming a
#'   v-structure) for each couple of variables of interest.
#'   The list is empty if the \emph{filter_v_struct} parameter is set to FALSE.
#' }
#'
#' @examples
#' library(miic)
#'
#' # Features selection on the path between an external metadata "Ploidy"
#' # (simulation by extracting "Ploidy" out of the dataset)
#' # and the gene expression of TP53
#' df_external_meta <- data.frame ("Ploidy" = cosmicCancer$Ploidy)
#' df_data <- cosmicCancer[ , ! (colnames(cosmicCancer) == "Ploidy") ]
#' ret <- selectFeaturesPath (df_data,
#'                            var_of_interest_values_side1=df_external_meta,
#'                            var_of_interest_names_side2="TP53",
#'                            n_features=10)
#' message ("Tops 10 features = ", paste (ret$features, collapse=", ") )
#'
#' # Same features selection with reuse of the MIs computed above and plot
#' if ( require(ggplot2) ) {
#'   ret <- selectFeaturesPath (df_data,
#'                              var_of_interest_values_side1=df_external_meta,
#'                              var_of_interest_names_side2="TP53",
#'                              precomputed_mis=ret$mis,
#'                              n_features=10,
#'                              plot=TRUE)
#'   print (ret$plots[[1]])
#' }
#'
#' # Features selection using multiple variables on each side
#' ret <- selectFeaturesPath (cosmicCancer,
#'   var_of_interest_names_side1=c("TP53", "Ploidy"),
#'   var_of_interest_names_side2=c("FOXM1", "AURKA"),
#'   n_features=10)
#' message ("Tops 10 features = ", paste (ret$features, collapse=", ") )
#'
#' # Same features selection with plotting using a customized rendering
#' if ( require(ggplot2) && require(gridExtra) ) {
#'   ret <- selectFeaturesPath (cosmicCancer,
#'                              var_of_interest_names_side1=c("TP53", "Ploidy"),
#'                              var_of_interest_names_side2=c("FOXM1", "AURKA"),
#'                              n_features=10,
#'                              plot=TRUE,
#'                              n_plots=10,
#'                              box_fill=c("lawngreen", "limegreen"),
#'                              box_text="black")
#'   do.call ("grid.arrange", c(ret$plots, ncol=2,
#'     top="Top features on cosmicCancer dataset") )
#' }
#'
#' @export
#-------------------------------------------------------------------------------
selectFeaturesPath <- function (input_data, n_features,
  var_of_interest_names_side1=NULL, var_of_interest_values_side1=NULL,
  var_of_interest_names_side2=NULL, var_of_interest_values_side2=NULL,
  skip_cheks=F, precomputed_mis=NULL, filter_v_struct=F,
  n_threads=1, verbose=3, plot=F, ...)
  {
  print ("TODO ? exlude voi from returned values ?")
  print ("TODO ? warning discrete number of levels (as miic) ?")
  # As we need to check both var_of_interest_names_sideX and/or
  # var_of_interest_values_sideX, so we can not check only n_features
  # as selectFeatures and let compute_mi_batch perform to other checks.
  # => do all the checks as compute_mi_batch + selectFeatures + some specific
  #
  # Check input_data
  #
  if (  ( ! is.data.frame (input_data) )
     && ( ! is.matrix(input_data) )
     && ( ! inherits(input_data, "Matrix") ) )
    miic:::miic_error  ("parameters",
      "the input data must be a data frame or a matrix.")
  if ( (ncol (input_data) <= 0) || (nrow (input_data) <= 0) )
    miic:::miic_error  ("parameters", "the input data is empty.")
  if ( is.data.frame (input_data) )
    # Ensure we have a true data frame, e.g. not a tibble
    # TODO evaluate run time impact on very large data frames
    # (let matrices unchanged to avoid warnings on large memory allocation)
    input_data <- as.data.frame (input_data)
  if ( is.null (colnames (input_data) ) )
    miic:::miic_error  ("parameters", "the input data must have column names.")
  if ( length(unique(colnames (input_data))) != ncol(input_data) )
    miic:::miic_error  ("parameters", "the input data have some column names duplicated.")
  #
  # Check VOI of each side
  #
  i = 2
  extra_voi_names_side1 = c()
  extra_voi_names_side2= c()
  for (i in 1:2)
    {
    if (i == 1)
      {
      var_of_interest_names = var_of_interest_names_side1
      var_of_interest_values = var_of_interest_values_side1
      }
    else
      {
      var_of_interest_names = var_of_interest_names_side2
      var_of_interest_values = var_of_interest_values_side2
      }
    #
    # Same kind of tests as compute_mi_batch
    #
    if ( is.null (var_of_interest_names) && is.null (var_of_interest_values) )
      miic:::miic_error  ("parameters", "the name of the variable(s) of interest",
        " or a data frame with the variable(s) of interest values must be supplied",
        " for side ", i, ".")

    if (is.null (var_of_interest_names) )
      var_of_interest_names <- c()
    else
      {
      for (one_var_name in var_of_interest_names)
        if ( miic:::test_param_wrong_string (one_var_name, colnames(input_data) ) )
          miic:::miic_error ("parameters",  "Some of the variable of interest",
            " names for side ", i, " are incorrect or not in the input_data.")
      }

    extra_voi_names <- c()
    if ( ! is.null (var_of_interest_values) )
      {
      if ( ! is.data.frame (var_of_interest_values) )
        miic:::miic_error  ("parameters",
          "the var_of_interest_values for side ", i, " must be a data frame.")
      # Ensure we have a true data frame, i.e. not a tibble
      var_of_interest_values <- as.data.frame (var_of_interest_values)
      if (ncol (var_of_interest_values) <= 0)
        {
        miic:::miic_warning  ("parameters",
          "the var_of_interest_values data frame for side ", i,
          " has been supplied but is empty.")
        var_of_interest_values <- NULL
        }
      else if ( nrow (var_of_interest_values) != nrow (input_data) )
        miic:::miic_error  ("parameters",
          "the variable of interest values for side ", i,
          " does not match the number of samples.")
      else
        {
        # Data frame OK, checks variables names not in data
        #
        extra_voi_names <- colnames (var_of_interest_values)
        #
        # Error or warning if voi requested as external in var_of_interest_values
        # are present in the input data
        #
        poss_wrong_idx = which ( extra_voi_names %in% colnames(input_data) )
        if (length (poss_wrong_idx) >= 1)
          {
          for (one_var_name in extra_voi_names[poss_wrong_idx])
            {
            one_var_voi_vals = var_of_interest_values[,one_var_name]
            one_var_input_vals = input_data[,one_var_name]

            if (any ( ( is.na (one_var_voi_vals) != is.na (one_var_input_vals) )
                    | (one_var_voi_vals[ !is.na(one_var_voi_vals) ] != one_var_input_vals[ !is.na(one_var_input_vals) ]) ) )
              miic:::miic_error  ("parameters",
                "the variable ", one_var_name, " is supplied both in input data",
                " and in variables of interest values of side ", i, ".")
            #
            # Supplied in both in input data and in variables of interest values
            # and with identical values => just a warning, use input_data
            # and ignore variables of interest values
            #
            miic:::miic_warning  ("parameters",
              "the variable ", one_var_name, " is supplied both in input data",
              " and in variables of interest values of side ", i, ".")
            var_of_interest_names = unique (c (var_of_interest_names, one_var_name) )
            var_of_interest_values[,one_var_name] = NULL
            }
          extra_voi_names <- colnames (var_of_interest_values)
          }
        }
      }

    if (i == 1)
      {
      var_of_interest_names_side1 = var_of_interest_names
      var_of_interest_values_side1 = var_of_interest_values
      extra_voi_names_side1 = extra_voi_names
      }
    else
      {
      var_of_interest_names_side2 = var_of_interest_names
      var_of_interest_values_side2 = var_of_interest_values
      extra_voi_names_side2= extra_voi_names
      }
    }
  #
  # Test on each side done, now test one side against the other
  #
  all_voi_names = c (var_of_interest_names_side1, var_of_interest_names_side2,
                     extra_voi_names_side1, extra_voi_names_side2)
  are_duplicated =  duplicated (all_voi_names)
  if (any (are_duplicated))
    miic:::miic_error  ("parameters",
      "Some variable(s) have been supplied in both side: ",
      miic:::list_to_str (all_voi_names[are_duplicated], n_max=10), ".")
  #
  # Create variables that store the vois of both sides
  #
  var_of_interest_names_all_sides = unique (c(var_of_interest_names_side1,
                                              var_of_interest_names_side2))
  if (  is.null (var_of_interest_values_side1)
     && is.null (var_of_interest_values_side2) )
    var_of_interest_values_all_sides = NULL
  else if (  is.null (var_of_interest_values_side1) )
    var_of_interest_values_all_sides = var_of_interest_values_side2
  else if (  is.null (var_of_interest_values_side2) )
    var_of_interest_values_all_sides = var_of_interest_values_side1
  else
    var_of_interest_values_all_sides = cbind (var_of_interest_values_side1,
                                              var_of_interest_values_side2)
  extra_voi_names_all_sides = unique (c (extra_voi_names_side1,
                                         extra_voi_names_side2) )
  #
  # Check parameters that will not be checked by compute_mi_batch
  #
  n_features <- miic:::check_param_int (n_features, "number of features",
                                        default=0, min=0)
  #
  # Init extra returned value to empty list for now
  #
  list_ret_ais = list()
  list_ret_filter_out = list()
  list_plots = list()
  #
  # Check other params and compite MIs
  #
  mat_mis = miic:::compute_mi_batch (input_data=input_data,
    var_of_interest_names=var_of_interest_names_all_sides,
    var_of_interest_values=var_of_interest_values_all_sides,
    skip_cheks=skip_cheks, precomputed_mis=precomputed_mis,
    n_threads=n_threads, verbose=verbose)
  n_vois = length (all_voi_names)
  mat_mis_filt = mat_mis[rownames(mat_mis) %in% colnames(input_data),
                         colnames(mat_mis) %in% all_voi_names,
                         drop=F]
  #
  # For each voi requested, filter on MI' > 0
  #
  list_mis_sup_0 = list()
  for (one_voi in all_voi_names)
    {
    mis_4_one_voi = mat_mis_filt[, one_voi]
    mis_4_one_voi = mis_4_one_voi[ ( ! is.na (mis_4_one_voi) )
                                 & ( mis_4_one_voi > 0) ]
    list_mis_sup_0[[one_voi]] = mis_4_one_voi
    }
  #
  # For each couple of voi (one from each side), search the MI of the
  # 2 voi
  #
  all_voi_names_side1 = c (var_of_interest_names_side1,
                           colnames(var_of_interest_values_side1) )
  all_voi_names_side2 = c (var_of_interest_names_side2,
                           colnames(var_of_interest_values_side2) )
  couples <- expand.grid (all_voi_names_side1, all_voi_names_side2,
                          stringsAsFactors=F)
  couples$mi = unlist (apply (couples, MARGIN=1, function (x) {
    if ( ! (x[[2]] %in% extra_voi_names_side2) )
      return (mat_mis_filt[ x[[2]], x[[1]] ])
    if ( ! (x[[1]] %in% extra_voi_names_side1) )
      return (mat_mis_filt[ x[[1]], x[[2]] ])
    #
    # 2 variables given as metadata => the MI has not been computed
    #
    list_vois_vals = list ("voi1" = var_of_interest_values_side1[ , x[[1]] ],
                           "voi2" = var_of_interest_values_side2[ , x[[2]] ])
    are_continuous = unlist (lapply (list_vois_vals, FUN=function(y) {
      return (  is.numeric (y)
             && (length (unique (y[!is.na(y)]) ) >= miic:::MIIC_CONTINUOUS_TRESHOLD) )
      } ) )
    ret = miic::computeMutualInfo (list_vois_vals[[1]], list_vois_vals[[2]],
                                   is_continuous=are_continuous, plot=F)
    return (ret$infok)
    } ) )
  couples_mi_0_test = (couples$mi <= 0)
  if ( any(couples_mi_0_test) )
    {
    couples_mi_0 = couples[couples_mi_0_test, , drop=F]
    couples = couples[!couples_mi_0_test, , drop=F]
    miic:::miic_warning ("path feature selection", "MI = 0 for ",
      paste ( apply (couples_mi_0, MARGIN=1, FUN=function(x) {
                paste0 (x[[1]], "-", x[[2]]) } ), collapse=", "),
      ", no feature selection possible on these couple(s)." )
    if (nrow (couples) <= 0)
      return (list ("features"=c(), "mis"=mat_mis, "couples"=couples,
        "dpis"=list(), "plots"=list_plots,
        "ais"= list_ret_ais, "v_struct_excluded"=list_ret_filter_out) )
    }
  #
  # Compute MI for DPI test on each feature and order by MI desc
  #
  list_dpi_sorted = list()
  i = 1
  for ( i in 1:nrow(couples) )
    {
    feat_voi1 = list_mis_sup_0[[ couples[i, 1] ]]
    feat_voi2 = list_mis_sup_0[[ couples[i, 2] ]]
    mi_threshold = couples[i, "mi"]
    feat_voi1 = feat_voi1[feat_voi1 >= mi_threshold]
    feat_voi2 = feat_voi2[feat_voi2 >= mi_threshold]
    feat_voi1 = feat_voi1[names(feat_voi1) %in% names(feat_voi2)]
    feat_voi2 = feat_voi2[names(feat_voi2) %in% names(feat_voi1)]
    dpi_sorted = feat_voi1 [names(feat_voi1)] + feat_voi2 [names(feat_voi1)]
    dpi_sorted = sort (dpi_sorted, decreasing=T)
    df_dpi_sorted = as.data.frame (dpi_sorted)
    colnames(df_dpi_sorted) = "sum"
    df_dpi_sorted[, couples[i, 1] ] = feat_voi1 [rownames(df_dpi_sorted)]
    df_dpi_sorted[, couples[i, 2] ] = feat_voi2 [rownames(df_dpi_sorted)]
    list_dpi_sorted[[i]] = df_dpi_sorted
    }
  #
  # If filter v-structures activated
  #
  if (filter_v_struct)
    {
    if (verbose == 1)
      miic:::miic_msg ("Filtering features using v-stucture...")
    i = 1
    for ( i in 1:nrow(couples) )
      {
      if (verbose >= 2)
        miic:::miic_msg ("Filtering with v-stuctures for couple ",
                  couples[i, 1], "-", couples[i, 2], "...")
      #
      # For each couple of voi, prepare an extract of data with all variables
      # verifying dpi
      #
      if (couples[i, 1] %in% var_of_interest_names_side1)
        data_1 = input_data[, couples[i, 1], drop=F]
      else
        data_1 = var_of_interest_values_side1[, couples[i, 1], drop=F]
      if (couples[i, 2] %in% var_of_interest_names_side2)
        data_2 = input_data[, couples[i, 2], drop=F]
      else
        data_2 = var_of_interest_values_side2[, couples[i, 2], drop=F]
      mis_4_couple = list_dpi_sorted[[i]]
      if (nrow (mis_4_couple) <= 0)
        next
      if (nrow (mis_4_couple) > 100)
        {
        miic:::miic_warning ("v-structure filtering",
          "number of features passsing the DPI test > 100 for couple ",
          couples[i, 1], "-" , couples[i, 2],
          ", only the 100 tops will be used to look for a separating set.")
        mis_4_couple = mis_4_couple[1:100, , drop=F]
        }
      data_3 = input_data[, rownames(mis_4_couple) ]
      if ( couples[i, 1] %in% colnames(data_3) )
        stop ("Bug : voi was in dpi")
      if ( couples[i, 2] %in% colnames(data_3) )
        stop ("Bug : voi was in dpi")
      df = data.frame (data_1, data_2, data_3, stringsAsFactors=F)
      #
      # Prepare a black box removing all edges expect between 2 voi
      #
      bb <- expand.grid (colnames(df), colnames(df),
                         stringsAsFactors=F)
      bb <- bb[ bb[,1] < bb[,2], , drop=F ]
      bb <- bb[ ! ( ( (bb[,1] == couples[i, 1])
                    & (bb[,2] == couples[i, 2]) )
                  | ( (bb[,1] == couples[i, 2])
                    & (bb[,2] == couples[i, 1]) ) ), , drop=F ]
      #
      # Run miic with latent on skeleton to find a contributors set
      #
      if (verbose >= 3)
        miic:::miic_msg ("Search for separating set...")
      res = miic (input_data=df, black_box=bb, latent="yes", orientation=F,
                  n_threads=n_threads, verbose=0)
      if (nrow (res$summary) != 1)
        miic:::miic_err ("v-structure filtering",
          " number of rows (", nrow(res$summary), ") incorrect in summary.")
      sum_row = res$summary[1, ]
      if (sum_row$type == "P")
        {
        miic:::miic_warning ("v-structure filtering", "no separating set",
          " found for ", couples[i, 1], "-", couples[i, 2], ".")
        next
        }
      if (is.na (sum_row$ai) || (sum_row$ai == "") )
        {
        miic:::miic_error ("v-structure filtering", "separting set empty for ",
                      couples[i, 1], "-", couples[i, 2], ".")
        next
        }
      #
      # Extra data of the separating set
      #
      list_ais = unlist (strsplit (sum_row$ai, ",", fixed=T) )
      list_ret_ais[[ paste0 (couples[i, 1], "-", couples[i, 2])]] = list_ais
      data_sep = input_data[ , list_ais, drop=F]
      #
      # Prepare are_continuous constant values
      #
      df_tmp = data.frame (data_1, data_2, data_sep, stringsAsFactors=F)
      are_continuous = unlist (lapply (df_tmp, FUN=function(one_col) {
        (  is.numeric (one_col)
        && (length (unique (one_col)) >= miic:::MIIC_CONTINUOUS_TRESHOLD) ) } ) )
      are_continuous[[length(are_continuous)+1]] = NA_integer_
      #
      # Test each feature to see if MI increases
      #
      if (verbose >= 3)
        miic:::miic_msg ("Filtering using v-structures...")
      feat_to_remove = rep ( F, nrow(list_dpi_sorted[[i]]) )
      feat_idx = 1
      for ( feat_idx in 1:nrow(list_dpi_sorted[[i]]) )
        {
        one_feat = rownames (list_dpi_sorted[[i]]) [[feat_idx]]
        if (one_feat %in% list_ais)
          next
        one_feat_vals = input_data[, one_feat]
        df_tmp = data.frame (data_sep, one_feat_vals, stringsAsFactors=F)
        are_continuous[[length(are_continuous)]] = ( ( is.numeric (one_feat_vals) )
          && (length (unique (one_feat_vals[!is.na(one_feat_vals)])) >= miic:::MIIC_CONTINUOUS_TRESHOLD) )

        res = miic::computeMutualInfo (x=data_1, y=data_2,
          df_conditioning=df_tmp, is_continuous=are_continuous)
        if (res$infok > 0)
          feat_to_remove[[feat_idx]] = T
        }
      if (any (feat_to_remove))
        {
        list_ret_filter_out[[ paste0 (couples[i, 1], "-", couples[i, 2])]] =
          rownames (list_dpi_sorted[[i]]) [feat_to_remove]
        miic:::miic_msg (sum (feat_to_remove),
          " features filtered out using v-structures.")
        list_dpi_sorted[[i]] = list_dpi_sorted[[i]] [!feat_to_remove, , drop=F]
        }
      else
        {
        list_ret_filter_out[[ paste0 (couples[i, 1], "-", couples[i, 2])]] = list()
        miic:::miic_msg ("No feature filtered out using v-structures.")
        }
      }
    }
  #
  # Store coumles names the returned list_dpi_sorted
  #
  names(list_dpi_sorted) = paste0 (couples[,1], "-", couples[,2])
  #
  # Plot if requested
  #
  i = 1
  if (plot)
    for ( i in 1:nrow(couples) )
      {
      couple_name = paste0 (couples[i,1], "-", couples[i,2])
      list_plots[[couple_name]] = plot_top_features_dpi (
        dpis=list_dpi_sorted[[i]], couple=couples[i, ,drop=F], ...)
      }
  #
  # Check if too much features to select
  #
  i = 1
  uniq_poss_feats = c()
  for ( i in 1:nrow(couples) )
    uniq_poss_feats = c (uniq_poss_feats, rownames(list_dpi_sorted[[i]]) )
  uniq_poss_feats = unique (uniq_poss_feats)

  if (n_features >= length(uniq_poss_feats) )
    {
    if (n_features > length(uniq_poss_feats) )
      miic:::miic_warning ("features selection",
        "too few variables share information with the variable(s)",
        " of interest to select the top ", n_features,
        ", returning all the ", length(uniq_poss_feats),
        " variables sharing information.")
    return (list ("features"=uniq_poss_feats, "mis"=mat_mis, "couples"=couples,
      "dpis"=list_dpi_sorted, "plots"=list_plots,
      "ais"= list_ret_ais, "v_struct_excluded"=list_ret_filter_out) )
    }
  #
  # Get top variables for each voi
  #
  list_tops_prec = c()
  list_tops = c()
  #
  # Select top features until we get enough or a bit too much
  #
  n_couples = nrow (couples)
  list_tops = c()
  list_tops_prec = c()
  n_tops_to_sel = n_features %/% n_couples
  while (T)
    {
    list_tops = unique (unlist (lapply (list_dpi_sorted, FUN=function(x) {
      rownames(x)[1:n_tops_to_sel] } ) ) )
    n_missing = n_features - length (list_tops)
    if (n_missing <= 0)
      break
    list_tops_prec = list_tops
    n_tops_to_sel = n_tops_to_sel + max (1, n_missing %/% n_couples)
    }
  #
  # If a bit too much, select in the last features added those with highest MI
  #
  if (n_missing < 0)
    {
    list_tops_to_test = unlist (lapply (list_dpi_sorted,
      FUN=function(x) { rownames(x)[n_tops_to_sel] } ) )
    list_dpis_to_test = unlist (lapply (list_dpi_sorted,
      FUN=function(x) { x[n_tops_to_sel] } ) )
    list_tops_to_test = list_tops_to_test[order (list_dpis_to_test,
                                                  decreasing=T)]
    list_tops_to_test = unique (list_tops_to_test)
    list_tops_to_test = list_tops_to_test[
      ! (list_tops_to_test %in% list_tops_prec) ]

    n_missing = n_features - length (list_tops_prec)
    list_tops = c (list_tops_prec, list_tops_to_test[1:n_missing])
    }
  if (verbose >= 1)
    miic:::miic_msg (length(list_tops), " features selected.")
  return (list ("features"=list_tops, "mis"=mat_mis, "couples"=couples,
    "dpis"=list_dpi_sorted, "plots"=list_plots,
    "ais"= list_ret_ais, "v_struct_excluded"=list_ret_filter_out) )
  }

