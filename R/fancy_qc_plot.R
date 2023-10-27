#' @title Plot Quality Control Metrics.
#'
#' @description Plot for visualizing QC metrics and allowing for grouping by different metadata columns.
#'
#' @details This function is readily available for visualizing a variety of quality control metrics. To get started, the user can easily overview all the available metrics with `return_plotdata = TRUE`.
#' When this parameter is set to TRUE, a vector of characters will be returned detailing all the, for this plot, available metrics. After deciding what metric to plot, simply give the metric of choice to the `plot_data` parameter.
#' This function also lets the user provide a data frame with sample IDs to be included in the plot. Optionally, the user can also provide an already filtered metadata table with sample IDs of interest to the `these_samples_metadata`.
#' If none of the two parameters are supplied, the user can easily restrict the plot to any cohort and/or pathology of their liking. This is done by calling `keep_cohort` and `keep_pathology`.
#' If these parameters are used, the function will retrieve metadata for all available GAMBL sample IDs and then subset to the specified cohort or pathology.
#' The layout of the returned plot can also be further customized with `sort_by`. This parameter controls the order in which samples would appear. Similarly, `fill_by` allows the user to control on what factor the plot will be filled by.
#' In addition, the generated plot can also be returned as an interactive HTML rendering, allowing the user to easily hover over any of the points in the plot and get expanded information on each data point. To toggle this function, set the `interactive` parameter to TRUE.
#' If an interactive plot is generated, it is also possible to dictate what information should be available in the plotted data points. Default for this parameter is sample ID and cohort.
#' Sometimes it can also be useful to see how a subset of samples compares to another group; to do this one could call the function with a vector of additional sample IDs given to the `comparison_samples` parameter (see examples for more information).
#' lastly, the plot can also be configured with custom plot title and axis labels (`plot_title` and `y_axis_lab`). For more information, see examples and parameter descriptions.
#'
#' @param these_sample_ids Data frame with sample IDs (to be plotted) in the first column (has to be named sample_id).
#' @param keep_cohort Optional parameter to be used when these_sample is NULL. Returns metadata and filters on the cohort supplied in this parameter.
#' @param keep_pathology Optional parameter to be used when these_sample is NULL. Returns metadata and filters on the pathology supplied in this parameter.
#' @param seq_type Selected seq type for incoming QC metrics.
#' @param metadata Optional, user can provide a metadata df to subset sample IDs from.
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process.
#' @param plot_data Plotting parameter, define the data type to be plotted.
#' @param fill_by Parameter for specifying fill variable for grouped bar plot. Can be any factor from incoming metadata, e.g pathology, cohort, etc.
#' @param labels If HTML plot version is rendered, you can specify what labels should be visible when hovering over the dots. Default is sample id and cohort. This parameter expects a vector of charachters.
#' @param interactive Boolean parameter for generating interactive plot (HTML). Default is FALSE.
#' @param comparison_samples Optional parameter, give the function a vector of sample IDs to be compared against the main plotting group. Pathology is default.
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
#' #load packages
#' library(dplyr)
#' 
#' #get sample IDs for available genome samples
#' genome_collated = collate_results(seq_type_filter = "genome") %>% 
#'   pull(sample_id)
#'
#' #subset the collated samples on BL samples
#' my_samples = get_gambl_metadata() %>% 
#'   dplyr::filter(sample_id %in% genome_collated) %>% 
#'   dplyr::filter(pathology == "BL") %>% pull(sample_id)
#'
#' fancy_qc_plot(these_sample_ids = my_samples, plot_data = "AverageBaseQuality")
#'
fancy_qc_plot = function(these_sample_ids,
                         keep_cohort,
                         keep_pathology,
                         seq_type = "genome",
                         metadata,
                         these_samples_metadata,
                         plot_data,
                         fill_by = "pathology",
                         labels = c("sample_id", "cohort"),
                         interactive = FALSE,
                         comparison_samples,
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

  #get gambl metadata (if not supplied)
  if(missing(metadata)){
    this_meta = GAMBLR.helpers::handle_metadata(this_seq_type = seq_type)
  }else{
    this_meta = metadata
  }

  if(!missing(these_samples_metadata)){
    these_sample_ids = dplyr::select(these_samples_metadata, sample_id) %>%
      as.data.frame(strings.as.factors = FALSE)
  }

  #filter metadata on selected cohort/pathology
  if(missing(these_sample_ids) && missing(these_samples_metadata)){
    if(!missing(keep_cohort) && missing(keep_pathology)){
      these_sample_ids = dplyr::filter(this_meta, cohort == keep_cohort)
    }

    if(!missing(keep_pathology) && missing(keep_cohort)){
      these_sample_ids = dplyr::filter(this_meta, pathology == keep_pathology)
    }

    if(!missing(keep_cohort) && !missing(keep_pathology)){
      these_sample_ids = dplyr::filter(this_meta, pathology == keep_pathology, cohort == keep_cohort)
    }

    if(missing(keep_cohort) && missing(keep_pathology)){
      these_sample_ids = dplyr::select(this_meta, sample_id)
    }
  }

  #get QC data for selected samples
  qc_metrics = collate_results(sample_table = these_sample_ids, seq_type_filter = seq_type)
  message(paste0("QC Metric successfully retreived for ", nrow(qc_metrics), " samples out of a total of ", nrow(these_sample_ids), " samples in input sample table."))

  #aggregate sample list with metadata columns
  qc_meta = qc_metrics %>%
    inner_join(this_meta) %>%
    mutate_if(is.integer, as.factor) %>%
    mutate_if(is.character, as.factor)

  qc_meta$group = "main_sample"

  #Retrieve QC metrics for comparison samples, if provided.
  if(!missing(comparison_samples)){
    comp_data = collate_results(sample_table = comparison_samples, seq_type_filter = seq_type)

    #aggregate sample list with metadata columns
    comp_meta = comp_data %>% inner_join(this_meta)
    comp_meta = mutate_if(comp_meta, is.integer, as.factor)
    comp_meta = mutate_if(comp_meta, is.character, as.factor)
    comp_meta$group = "comparison_sample"
    qc_meta = rbind(qc_meta, comp_meta)
  }

  #get gambl colours for selected fill and subset to levels in selected factor
  col_gambl = GAMBLR.helpers::get_gambl_colours(fill_by) %>%
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
