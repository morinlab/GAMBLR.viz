#' @title Plot Alignment Metrics
#'
#' @description Visualize (stacked barplot) genomic read-subsets (metrics) across a selection of samples.
#'
#' @details This function is available for plotting relevant alignment metrics (read-subsets) for a selection of samples. Per default, this plot returns the following read-metrics;
#' total n reads, total n uniquely mapped reads, total n duplicated reads. This plot can also be superimposed with read metrics from additional samples,
#' allowing for easy comparisons between different sample populations. To run this function, simply specify the sample IDs you are interested in with `these_sample_ids`.
#' This parameter expects a data frame with sample IDs in the first column. Optionally, the user can also provide an already subset (with the sample IDS of interest)
#' metadata table with `these_samples_metadata`. For adding a comparison group to the returned plot, simply give another cohort/set of samples to the `comparison_group` parameter.
#' Similarly to `these_sample_ids`, this parameter also expects a data frame with sample IDs in the first column. In addition, this plot can also add additional read-metrics such as
#' mean values for all plotted metrics and corrected coverage. To enable these features, simply set `add_mean` and `add_corrected_coverage` to TRUE (default).
#'
#' @param these_sample_ids Data frame with sample IDs (to be plotted) in the first column.
#' @param metadata Optional argument, used to derive sample IDs if sample_table is Null.
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process.
#' @param comparison_group Optional argument for plotting mean alignment metrics. Default is plotting the mean for samples provided. This parameter takes a vector of sample IDs.
#' @param seq_type Subset qc metrics to a specific seq_type, default is genome.
#' @param add_mean Set to TRUE to superimpose mean values of plotted variables. Default is TRUE.
#' @param add_corrected_coverage Set to TRUE to add corrected coverage for selected samples.
#' @param keep_cohort If no df with sample IDs is supplied (these_sample_ids = NULL) the function calls get_gambl_metadata and subsets on selected cohorts.
#' @param keep_pathology If no df with sample IDs is supplied (these_sample_ids = NULL) the function calls get_gambl_metadata and subsets on selected pathology.
#' @param this_color_palette Optional parameter that holds the selected colours for the plotted bars.
#' @param plot_sub Optional parameter, add a subtitle to the alignment metric plot.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import ggplot2 cowplot dplyr
#' @export
#'
#' @examples
#' #Example 1 - using these_sample_ids parameter
#' #subset on FL cases with QC metrics available and plot
#' metadata = get_gambl_metadata()
#' kridel_fl = dplyr::filter(metadata, pathology == "FL",
#'                cohort == "FL_Kridel")
#'
#' kridel_fl_samples = dplyr::select(kridel_fl, sample_id)
#'
#' fancy_alignment_plot(these_sample_ids = kridel_fl_samples)
#'
#' #Example 2 - using already filtered metadata (these_samples_metadata)
#' fancy_alignment_plot(these_samples_metadata = kridel_fl)
#'
#' #Example 3 - using in-house metadata filtering options
#' fancy_alignment_plot(keep_cohort = "FL_Kridel",
#'                      keep_pathology = "FL")
#'
fancy_alignment_plot = function(these_sample_ids,
                                metadata,
                                these_samples_metadata,
                                comparison_group,
                                seq_type = "genome",
                                add_mean = TRUE,
                                add_corrected_coverage = TRUE,
                                keep_cohort,
                                keep_pathology,
                                this_color_palette = c("TotalReads" = "#3D405B",
                                                       "TotalUniquelyMapped" = "#81B29A",
                                                       "TotalDuplicatedreads" = "#E07A5F"),
                                plot_sub = ""){

  #get gambl metadata (if not supplied)
  if(missing(metadata)){
    this_meta = get_gambl_metadata(seq_type_filter = seq_type)
  }else{
    this_meta = metadata
  }

  if(!missing(these_samples_metadata)){
    these_sample_ids = dplyr::select(these_samples_metadata, sample_id) %>%
      as.data.frame(strings.as.factors = FALSE)
  }


  #filter metadata on selected cohort/pathology
  if(missing(these_sample_ids)){
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

  #get qc data for selected samples
  qc_metrics = collate_results(sample_table = these_sample_ids, seq_type_filter = seq_type)

  message(paste0("QC Metric successfully retreived for ", nrow(qc_metrics),
                 " samples out of a total of ", nrow(these_sample_ids), " samples in input sample table."))

  #subset alignment metrics
  melt_align = dplyr::select(qc_metrics, c(sample_id, TotalReads, TotalUniquelyMapped, TotalDuplicatedreads)) %>%
    as.data.table() %>% 
    melt(id.var = "sample_id") %>%
    arrange(sample_id)

  mean_cov_df = data.frame(Metric = c("TotalReads", "TotalUniquelyMapped", "TotalDuplicatedreads"),
                           Value = c (mean(melt_align$value[melt_align$variable == "TotalReads"]),
                                      mean(melt_align$value[melt_align$variable == "TotalUniquelyMapped"]),
                                      mean(melt_align$value[melt_align$variable == "TotalDuplicatedreads"])))

  if(!missing(comparison_group)){
    comp_data = collate_results(sample_table = comparison_group, seq_type_filter = seq_type) %>%
      dplyr::select(sample_id, TotalReads, TotalUniquelyMapped, TotalDuplicatedreads) %>%
      as.data.table() %>% 
      melt(id.var = "sample_id") %>%
      arrange(sample_id)

    mean_cov_df_comp = data.frame(Metric = c("TotalReads", "TotalUniquelyMapped", "TotalDuplicatedreads"),
                                  Value = c (mean(comp_data$value[comp_data$variable == "TotalReads"]),
                                             mean(comp_data$value[comp_data$variable == "TotalUniquelyMapped"]),
                                             mean(comp_data$value[comp_data$variable == "TotalDuplicatedreads"])))
  }

  #corrected mean coverage
  melt_cov = dplyr::select(qc_metrics, c(sample_id, MeanCorrectedCoverage)) %>%
    as.data.table() %>% 
    melt(id.var = "sample_id") %>%
    arrange(sample_id)

  #plot alignment data
  p = ggplot() +
    geom_bar(melt_align, mapping = aes(x = sample_id, y = value, fill = variable), position = "dodge", stat = "identity") +
    {if(add_mean)geom_hline(mean_cov_df, mapping = aes(yintercept = Value, color = Metric))} +
    {if(!missing(comparison_group)) geom_hline(mean_cov_df_comp, mapping = aes(yintercept = Value, linetype = Metric), color = "#54555E")} +
    {if(add_corrected_coverage)geom_point(melt_cov, mapping = aes(x = sample_id, y = value * 25000000, shape = variable), fill = "#A892B3", color = "#5C4966", size = 3, position = "dodge", stat = "identity")} +
    {if(add_corrected_coverage)scale_y_continuous(expand = c(0, 0), breaks = seq(0, max(melt_align$value), by = 3e+08), sec.axis = sec_axis(~ . / 2500000000, name = "", labels = function(b){paste0(round(b * 100, 0), "X")}))} +
    {if(!add_corrected_coverage)scale_y_continuous(expand = c(0, 0), breaks = seq(0, max(melt_align$value), by = 3e+08))} +
    labs(title = "Alignment Summary", subtitle = plot_sub, x = "", y = "Reads (n)") +
    scale_fill_manual(values = this_color_palette) +
    scale_linetype_manual(values = c("TotalReads" = "solid", "TotalUniquelyMapped" = "dashed", "TotalDuplicatedreads" = "dotted")) +
    scale_shape_manual(values = c("MeanCorrectedCoverage" = 21)) +
    scale_color_manual(values = c(this_color_palette)) +
    theme_cowplot() +
    labs(linetype = "Comparison Group", shape = "Corrected Coverage (right y-axis)", fill = "Alignment Metrics", color = "Alignment Metrics (Mean)") +
    theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank())

  return(p)
}
