#' @title Proportional Metrics Plot.
#'
#' @details Visualize proportional metrics for selected samples.
#'
#' @description This function also provides parameters for easy customization of the plot aesthetics. For example, the subtitle of the plot can easily be controlled with the `plot_subtitle` parameter.
#' For usage examples and more information, refer to the parameter descriptions and examples in the fancy vignette.
#'
#' @param collated_results Required parameter. A data frame with collated results for the sample IDs of interest. Preferably, the return from [GAMBLR::collate_results].
#' @param plot_subtitle Plotting parameter, subtitle of generated plot.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import dplyr ggplot2 cowplot
#' @export
#'
#' @examples
#' #get data
#' dlbcl_a_collated = collate_results() %>% dplyr::filter(BL_subgroup == "DLBCL-A")
#' 
#' #build plot
#' fancy_proportions_plot(collated_results = dlbcl_a_collated, plot_subtitle = "Example Plot")  
#'
fancy_proportions_plot = function(collated_results,
                                  plot_subtitle = ""){

  #check if required parameter is provided
  if(missing(collated_results)){
    stop("Please provide a data frame with collated results to the `collated_reesults` parameter...")
  }

  #data wrangling
  qc_sub = dplyr::select(collated_results, sample_id, ProportionReadsDuplicated, ProportionReadsMapped, ProportionCoverage10x, ProportionCoverage30x) %>%
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
