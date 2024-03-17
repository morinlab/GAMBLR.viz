#' @title Proportional Coverage Plot.
#'
#' @description Visualize proportional coverage (10X and 30X) for selected samples and add comparison group (optional).
#'
#' @details Create a highly customizable plot visualizing proportional alignment metrics i.e what proportions of the aligned reads that show 10x and 30x coverage.
#' This function provides straightforward subsetting parameters allowing for a straightforward execution. Either provide a data frame with sample IDs in the first column to the `these_samples_ids` parameter.
#' Or, call one of the optional parameters for using an already subset metadata table (subset to the sample IDs of interest).
#' If `these_samples_ids` and `these_samples_metadata` is not provided, the user can subset al GAMBL samples on the fly with `keep_cohort` and/or `keep_pathology`.
#' This function can also plot the same results for a comparison group of interest, i.e another table with sample IDs. This can be useful for visualising how certain cohorts/pathologies compares to each other.
#' To do this one would call the function with `comparison_samples` parameter. For more info on how to use this function, please refer to examples, vignettes (fancy_vignette) and parameter descriptions.
#'
#' @param these_sample_ids Data frame with sample IDs (to be plotted) in the first column.
#' @param metadata Optional, user can provide a metadata df to subset sample IDs from.
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process.
#' @param keep_cohort Optional parameter to be used when these_sample is NULL. Obtains metadata and filters on the cohort supplied in this parameter.
#' @param keep_pathology Optional parameter to be used when these_sample is NULL. Obtains metadata and filters on the pathology supplied in this parameter.
#' @param comparison_samples Optional parameter, give the function a vector of sample IDs to be compared against the main plotting group.
#' @param seq_type Selected seq type for incoming QC metrics.
#' @param plot_subtitle Plotting parameter, subtitle of generated plot.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import dplyr ggplot2 GAMBLR.helpers
#' @export
#'
#' @examples
#' #load packages
#' library(dplyr)
#' library(GAMBLR.data)
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
#' fancy_propcov_plot(these_sample_ids = my_samples)
#'
fancy_propcov_plot = function(these_sample_ids,
                              metadata,
                              these_samples_metadata,
                              keep_cohort,
                              keep_pathology,
                              comparison_samples,
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

  #retrieve data for comparison, provided as a df with sample IDs in the first column (subset from gambl metadata)
  if(!missing(comparison_samples)){
    comp_data = collate_results(sample_table = comparison_samples, seq_type_filter = seq_type) %>%
      dplyr::select(ProportionCoverage10x, ProportionCoverage30x) %>%
      gather(Type, Value)

    comp_data$Type = as.factor(comp_data$Type)

    levels(comp_data$Type)[levels(comp_data$Type)=="ProportionCoverage10x"] = "comparison_group_10X"
    levels(comp_data$Type)[levels(comp_data$Type)=="ProportionCoverage30x"] = "comparison_group_30X"
  }

  #get QC data for selected samples
  qc_metrics = collate_results(sample_table = these_sample_ids, seq_type_filter = seq_type)
  message(paste0("QC Metric successfully retreived for ", nrow(qc_metrics), " samples out of a total of ", nrow(these_sample_ids), " samples in input sample table."))

  #data wrangling steps
  sub_metrics = dplyr::select(qc_metrics, ProportionCoverage10x, ProportionCoverage30x) %>%
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
    theme_Morons() +
    scale_fill_manual(values = c("#dda15e", "#606c38", "#433D6B", "#6B3254"))

  return(p)
}
