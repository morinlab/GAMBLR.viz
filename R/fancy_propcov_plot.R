#' @title Proportional Coverage Plot.
#'
#' @description Visualize proportional coverage (10X and 30X) for selected samples and add comparison group (optional).
#'
#' @details Create a highly customizable plot visualizing proportional alignment metrics i.e what proportions of the aligned reads that show 10x and 30x coverage.
#' This function can also plot the same results for a comparison group of interest, i.e another table with sample IDs. This can be useful for visualizing how certain cohorts/pathologies compares to each other.
#' To do this one would call the function with `comparison_samples` parameter. For more info on how to use this function, please refer to examples, vignettes (fancy_vignette) and parameter descriptions.
#'
#' @param collated_results Required parameter. A data frame with collated results for the sample IDs of interest. Preferably, the return from `collate_results` in GAMBLR.
#' @param comparison_samples Optional parameter, give the function a vector of sample IDs to be compared against the main plotting group.
#' @param seq_type Selected seq type for incoming QC metrics.
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
#' dlbcl_b_collated = collate_results() %>% dplyr::filter(BL_subgroup == "DLBCL-B")
#' 
#' #build plot
#' fancy_propcov_plot(collated_results = dlbcl_a_collated, comparison_samples = dlbcl_b_collated)
#'
fancy_propcov_plot = function(collated_results,
                              comparison_samples,
                              plot_subtitle = ""){

  #load QC data for selected samples
  if(missing(collated_results)){
     stop("Please provide a data frame with collated results to the `collated_reesults` parameter...")
  }
  
  #retrieve data for comparison, provided as a df with sample IDs in the first column (subset from gambl metadata)
  if(!missing(comparison_samples)){
    comp_data = comparison_samples %>%
      dplyr::select(ProportionCoverage10x, ProportionCoverage30x) %>%
      gather(Type, Value)
    
    comp_data$Type = as.factor(comp_data$Type)
    
    levels(comp_data$Type)[levels(comp_data$Type)=="ProportionCoverage10x"] = "comparison_group_10X"
    levels(comp_data$Type)[levels(comp_data$Type)=="ProportionCoverage30x"] = "comparison_group_30X"
  }

  #data wrangling steps
  sub_metrics = dplyr::select(collated_results, ProportionCoverage10x, ProportionCoverage30x) %>%
    gather(Type, Value)

  sub_metrics$Type = as.factor(sub_metrics$Type)

  levels(sub_metrics$Type)[levels(sub_metrics$Type)=="ProportionCoverage10x"] = "selected_samples_10X"
  levels(sub_metrics$Type)[levels(sub_metrics$Type)=="ProportionCoverage30x"] = "selected_samples_30X"

  #combine comparison data with sample data
  if(!missing(comparison_samples)){
    sub_metrics = rbind(sub_metrics, comp_data) %>%
      mutate(Type = factor(Type, levels = c("selected_samples_10X", "comparison_group_10X", "selected_samples_30X", "comparison_group_30X")))
  }

  #plotting
  p = ggplot(data = sub_metrics, aes(x = Type, y = Value, fill = Type)) +
    geom_violin(trim = FALSE, scale = "width") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white") +
    ylim(0, 1) +
    labs(title = "Proportion Coverage", subtitle = plot_subtitle, x = "", y = "Fraction") +
    theme_cowplot() +
    scale_fill_manual(values = c("#dda15e", "#606c38", "#433D6B", "#6B3254")) +
    theme(legend.position = "right", legend.title = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), panel.background = element_blank())

  return(p)
}
