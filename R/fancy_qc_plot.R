#' @title Plot Quality Control Metrics.
#'
#' @description Plot for visualizing QC metrics and allowing for grouping by different metadata columns.
#'
#' @details This function is readily available for visualizing a variety of quality control metrics. 
#' To get started, the user can easily overview all the available metrics with `return_plotdata = TRUE`.
#' When this parameter is set to TRUE, a vector of characters will be returned detailing all the, for this plot, available metrics. 
#' After deciding what metric to plot, simply give the metric of choice to the `plot_data` parameter.
#' This function expects the user to provide a data frame with collated results (see `collate_results` in GAMBLR) given to the `collated_results` parameter.
#' `fill_by` allows the user to control on what factor the plot will be filled by.
#' In addition, the generated plot can also be returned as an interactive HTML rendering, 
#' allowing the user to easily hover over any of the points in the plot and get expanded information on each data point. 
#' To toggle this function, set the `interactive` parameter to TRUE.
#' If an interactive plot is generated, it is also possible to dictate what information should be available in the plotted data points. 
#' Default for this parameter is sample ID and cohort.
#' Sometimes it can also be useful to see how a subset of samples compares to another group;
#' to do this one could call the function with a vector of additional sample IDs given to the `comparison_samples` parameter.
#' Lastly, the plot can also be configured with custom plot title and axis labels (`plot_title` and `y_axis_lab`). 
#' For more information, see examples and parameter descriptions.
#'
#' @param collated_results Required parameter. A data frame with collated results for the sample IDs of interest. Preferably, the return from `collate_results` in GAMBLR.
#' @param fill_by Parameter for specifying fill variable for grouped bar plot. Can be any factor from incoming metadata, e.g pathology, cohort, etc.
#' @param labels If HTML plot version is rendered, you can specify what labels should be visible when hovering over the dots. Default is sample id and cohort. This parameter expects a vector of characters.
#' @param interactive Boolean parameter for generating interactive plot (HTML). Default is FALSE.
#' @param comparison_samples Optional parameter, give the function a vector of sample IDs to be compared against the main plotting group. Pathology is default.
#' @param plot_data Plotting parameter, define the data type to be plotted.
#' @param plot_title Plotting parameter, plot title.
#' @param y_axis_lab Plotting parameter, label of y-axis.
#' @param return_plotdata Optional parameter, if set to TRUE a vector of acceptable data types for plotting will be returned, and nothing else.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @rawNamespace import(plotly, except = c("last_plot", "add_heatmap", "export"))
#' @import dplyr ggplot2 cowplot ggbeeswarm
#' @export
#'
#' @examples
#' #get data
#' collated_fl = collate_results(join_with_full_metadata = TRUE) %>% dplyr::filter(pathology == "FL")
#' 
#' #build plot
#' fancy_qc_plot(collated_results = collated_fl, plot_data = "AverageBaseQuality")
#'
fancy_qc_plot = function(collated_results,
                         comparison_samples,
                         plot_data,
                         fill_by = "pathology",
                         labels = c("sample_id", "cohort"),
                         interactive = FALSE,
                         plot_title = "",
                         y_axis_lab = "",
                         return_plotdata = FALSE){

  #return a list of acceptable data types for plotting
  if(return_plotdata){
    plotting_variables = c("AverageBaseQuality", "AverageInsertSize", "AverageReadLength",
                           "PairsOnDiffCHR", "TotalReads", "TotalUniquelyMapped",
                           "TotalUnmappedreads", "TotalDuplicatedreads", "ProportionReadsDuplicated",
                           "ProportionReadsMapped", "MeanCorrectedCoverage", "ProportionCoverage10x", "ProportionCoverage30x")

    return(plotting_variables)
  }

  #check if required parameter is provided
  if(missing(collated_results)){
    stop("Please provide a data frame with collated results to the `collated_reesults` parameter...")
  }

  #aggregate sample list with metadata columns
  qc_meta = collated_results %>%
    mutate_if(is.integer, as.factor) %>%
    mutate_if(is.character, as.factor)

  qc_meta$group = "main_sample"

  if(!missing(comparison_samples)){
    comp_meta = comparison_samples
    comp_meta = mutate_if(comp_meta, is.integer, as.factor)
    comp_meta = mutate_if(comp_meta, is.character, as.factor)
    comp_meta$group = "comparison_sample"
    qc_meta = rbind(qc_meta, comp_meta)
  }

  #get gambl colours for selected fill and subset to levels in selected factor
  col_gambl = get_gambl_colours(fill_by) %>%
    as.data.frame()

  col_gambl$factors = rownames(col_gambl)
  colnames(col_gambl)[1] = "hex"
  row.names(col_gambl) <- NULL

  levels_fill = levels(qc_meta[[fill_by]]) %>%
    as.data.frame()

  colnames(levels_fill)[1] = "factors"
  sub_cols = dplyr::left_join(levels_fill, col_gambl, by = "factors")
  list_col = as.list(sub_cols$hex)

  #plotting
  p = ggplot(qc_meta) +
    {if(interactive)aes_string(x = paste0("group"), y = plot_data, fill = fill_by, label1 = labels[1], label2 = labels[2])} +
    {if(!interactive)aes_string(x = paste0("group"), y = plot_data, fill = fill_by)} +
    geom_boxplot(mapping = aes(x = group), outlier.shape = NA) +
    geom_quasirandom() +
    labs(title = plot_title, x = "", y = y_axis_lab) +
    theme_cowplot() +
    scale_fill_manual(values = c(list_col)) +
    theme(legend.position = "right", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank())

  #make plot interactive (html) with plotly.
  if(interactive){
    p = ggplotly(p)
  }
  return(p)
}
