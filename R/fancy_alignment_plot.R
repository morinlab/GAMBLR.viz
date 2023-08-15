#' @title Plot Alignment Metrics
#'
#' @description Visualize (stacked barplot) genomic read-subsets (metrics) across a selection of samples.
#'
#' @details This function is available for plotting relevant alignment metrics (read-subsets) for a selection of samples. Per default, this plot returns the following read-metrics;
#' total n reads, total n uniquely mapped reads, total n duplicated reads. This plot can also be superimposed with read metrics from additional samples,
#' allowing for easy comparisons between different sample populations. To run this function, provide a data frame with collated results [GAMBLR::colalte_results] subset to the sample IDs of interest.
#' Similarly, if a comparison group is to be superimposed to the returned plot, provide another collated results for these sample IDs (with [GAMBLR::collate_results]).
#' In addition, this plot can also add additional read-metrics such as mean values for all plotted metrics and corrected coverage. 
#' To enable these features, simply set `add_mean` and `add_corrected_coverage` to TRUE (default).
#'
#' @param collated_results Required parameter. A data frame with collated results for the sample IDs of interest. Preferably, the return from [GAMBLR::collate_results].
#' @param comparison_group Optional argument for plotting mean alignment metrics from another group off sample IDs. This should be a data frame with collated results [GAMBLR::collate_results] with sample IDs to be compared. 
#' @param add_mean Set to TRUE to superimpose mean values of plotted variables. Default is TRUE.
#' @param add_corrected_coverage Set to TRUE to add corrected coverage for selected samples. Default is TRUE.
#' @param plot_sub Optional parameter, add a subtitle to the alignment metric plot.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import ggplot2 cowplot dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' dlbcl_a_collated = all_collated %>% dplyr::filter(BL_subgroup == "DLBCL-A")
#' dlbcl_b_collated = all_collated %>% dplyr::filter(BL_subgroup == "DLBCL-B")
#' 
#' fancy_alignment_plot(collated_results = dlbcl_a_collated, 
#'                      comparison_group = dlbcl_b_collated, 
#'                      add_mean = TRUE, 
#'                      add_corrected_coverage = TRUE, 
#'                      plot_sub = "Example Plot")
#' }
#'
fancy_alignment_plot = function(collated_results,
                                comparison_group,
                                add_mean = TRUE,
                                add_corrected_coverage = TRUE,
                                plot_sub = ""){
  #subset alignment metrics
  melt_align = dplyr::select(collated_results, c(sample_id, TotalReads, TotalUniquelyMapped, TotalDuplicatedreads)) %>%
    as.data.table() %>% 
    melt(id.var = "sample_id") %>%
    arrange(sample_id)
  
  mean_cov_df = data.frame(Metric = c("TotalReads", "TotalUniquelyMapped", "TotalDuplicatedreads"),
                           Value = c (mean(melt_align$value[melt_align$variable == "TotalReads"]),
                                      mean(melt_align$value[melt_align$variable == "TotalUniquelyMapped"]),
                                      mean(melt_align$value[melt_align$variable == "TotalDuplicatedreads"])))
  
  if(!missing(comparison_group)){
    comp_data = comparison_group %>%
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
  melt_cov = dplyr::select(collated_results, c(sample_id, MeanCorrectedCoverage)) %>%
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
    scale_fill_manual(values = c("TotalReads" = "#3D405B", "TotalUniquelyMapped" = "#81B29A", "TotalDuplicatedreads" = "#E07A5F")) +
    scale_linetype_manual(values = c("TotalReads" = "solid", "TotalUniquelyMapped" = "dashed", "TotalDuplicatedreads" = "dotted")) +
    scale_shape_manual(values = c("MeanCorrectedCoverage" = 21)) +
    scale_color_manual(values = c("TotalReads" = "#3D405B", "TotalUniquelyMapped" = "#81B29A", "TotalDuplicatedreads" = "#E07A5F")) +
    theme_cowplot() +
    labs(linetype = "Comparison Group", shape = "Corrected Coverage (right y-axis)", fill = "Alignment Metrics", color = "Alignment Metrics (Mean)") +
    theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank())
                                           
  return(p)  
}

  
