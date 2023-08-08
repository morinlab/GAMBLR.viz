#' @title Proportional Metrics Plot.
#'
#' @details Visualize proportional metrics for selected samples.
#'
#' @description This function takes all the available proportional quality control metrics and a vector of sample IDs (plotted on the x-ais) and plots the vlaues along the y-axis.
#' This function provides straightforward subsetting parameters allowing for a straightforward execution.
#' Either provide a data frame with sample IDs in the first column to the `these_samples_ids` parameter.
#' Or, call one of the optional parameters for using an already subset metadata table (subset to the sample IDs of interest).
#' If `these_samples_ids` and `these_samples_metadata` is not provided, the user can subset al GAMBL samples on the fly with `keep_cohort` and/or `keep_pathology`.
#' This function also provides parameters for easy customization of the plot aesthetics. For example, the subtitle of the plot can easily be controlled with the `plot_subtitle` parameter.
#' For usage examples and more information, refer to the parameter descriptions and examples in the fancy vignette.
#'
#' @param these_sample_ids Data frame with sample IDs (to be plotted) in the first column.
#' @param metadata Optional, user can provide a metadata df to subset sample IDs from.
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process.
#' @param keep_cohort Optional parameter to be used when these_sample is NULL. Returns metadata and filters on the cohort supplied in this parameter.
#' @param keep_pathology Optional parameter to be used when these_sample is NULL. Returns metadata and filters on the pathology supplied in this parameter.
#' @param seq_type Selected seq type for incoming QC metrics.
#' @param plot_subtitle Plotting parameter, subtitle of generated plot.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import dplyr ggplot2 cowplot
#' @export
#'
#' @examples
#' #Example 1 - using these_sample_ids parameter
#' #subset on FL cases with QC metrics available and plot
#' metadata = GAMBLR.data::gambl_metadata
#' kridel_fl = dplyr::filter(metadata, pathology == "FL",
#'                cohort == "FL_Kridel")
#' kridel_fl_samples = dplyr::select(kridel_fl, sample_id)
#'
#' fancy_proportions_plot(these_sample_ids = kridel_fl_samples)
#'
#' #Example 2 - using already filtered metadata (these_samples_metadata)
#' fancy_proportions_plot(these_samples_metadata = kridel_fl)
#'
#' #Example 3 - using in-house metadata filtering options
#' fancy_proportions_plot(keep_cohort = "FL_Kridel",
#'                    keep_pathology = "FL")
#'
fancy_proportions_plot = function(these_sample_ids,
                                  metadata,
                                  these_samples_metadata,
                                  keep_cohort,
                                  keep_pathology,
                                  seq_type = "genome",
                                  plot_subtitle = ""){

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

  #data wrangling
  qc_sub = dplyr::select(qc_metrics, sample_id, ProportionReadsDuplicated, ProportionReadsMapped, ProportionCoverage10x, ProportionCoverage30x) %>%
    gather(Type, Value, -sample_id)

  qc_sub$Type = as.factor(qc_sub$Type)

  levels(qc_sub$Type)[levels(qc_sub$Type)=="ProportionReadsDuplicated"] = "Duplicated Reads"
  levels(qc_sub$Type)[levels(qc_sub$Type)=="ProportionReadsMapped"] = "Mapped Reads"
  levels(qc_sub$Type)[levels(qc_sub$Type)=="ProportionCoverage10x"] = "10X"
  levels(qc_sub$Type)[levels(qc_sub$Type)=="ProportionCoverage30x"] = "30X"

  #get means for each metric
  mean_df = data.frame(Metric = c("Duplicated Reads", "Mapped Reads", "10X", "30X"),
                       Value = c(mean(qc_sub$Value[qc_sub$Type == "Duplicated Reads"]),
                                 mean(qc_sub$Value[qc_sub$Type == "Mapped Reads"]),
                                 mean(qc_sub$Value[qc_sub$Type == "10X"]),
                                 mean(qc_sub$Value[qc_sub$Type == "30X"])))

  #set colors
  these_colors = c("Duplicated Reads" = "#B8794D",
                   "Mapped Reads" = "#366B32",
                   "10X" = "#DDA15E",
                   "30X" = "#433D6B")

  #plotting
  p = ggplot(qc_sub, aes(x = sample_id, y = Value, group = Type)) +
    geom_bar(aes(fill = Type), position = "dodge", stat = "identity", width = 0.7) +
    geom_hline(mean_df, mapping = aes(yintercept = Value, color = Metric), size = 0.5, linetype = "dashed") +
    labs(title = "QC Metrics As Proportions (of All Reads)", subtitle = plot_subtitle, x = "", y = "Proportion") +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, by = 0.1)) +
    scale_fill_manual(values = these_colors) +
    scale_color_manual(values = these_colors) +
    labs(fill = "Metrics", color = "Mean") +
    theme(legend.position = "right", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.background = element_blank(),
          axis.line.x = element_line(color = "black", size = 0.5, linetype = 1),
          axis.line.y = element_line(color = "black", size = 0.5, linetype = 1),
          axis.ticks.y = element_line(color = "black", size = 0.5, linetype = 1))

  return(p)
}
