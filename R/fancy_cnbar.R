#' @title Copy Number states barplot
#'
#' @description Plot sample-specific CN states and affected bases for each segment
#'
#' @details `fancy_cnbar` visualizes copy number (CN) states on sample-level. Similarly to other fancy_x_plots this function
#' accepts either a sample ID, for which the function will get copy number states with [GAMBLR::get_sample_cn_segments]. The function
#' can also accept an already loaded seg file given to the `seg_data` parameter. It can also load a seg file with the `seg_path`
#' parameter. If the user calls either `seg_data` or `seg_path`, there are a collection of parameters available for specifying
#' the relevant columns in the given data frame (`chrom_col`, `starat_col`, `end_col`, `cn_col`). It is also possible to
#' restrict the returned plot to any given chromosome. This is done with the `chr_select` parameter (default is all autosomes).
#' For further control of the returned plot, it is also possible to set the threshold for maximum CN states to be returned (default is 15).
#' With `include_cn2` (Boolean) the user can control if CN segments = 2 should be added to the plot, default is TRUE.
#' The user can also control the annotations of the returned plot with `plot_title` and `plot_subtitle`. Lastly,
#' This function also computes the number of affected bases for each copy number state and plots these values on a secondary y-axis (right),
#' useful for overviewing the extent of each copy number state, in the context of the full genome.
#'
#' @param this_sample_id Sample to be plotted.
#' @param seg_data Optional parameter with copy number df already loaded into R.
#' @param seg_path Optional parameter with path to external cn file.
#' @param chrom_col Index of column with chromosome annotations (to be used with either maf_data or maf_path).
#' @param start_col Index of column with copy number start coordinates (to be used with either maf_data or maf_path).
#' @param end_col Index of column with copy number end coordinates (to be used with either maf_data or maf_path).
#' @param cn_col Index of column holding copy number information (to be used with either maf_data or maf_path).
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Subtitle for created plot.
#' @param chr_select Vector of chromosomes to be included in plot, defaults to autosomes.
#' @param cutoff Set threshold for maximum CN state to be retrieved.
#' @param include_cn2 Optional boolean statement for including CN = 2 states in plot.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import ggplot2 dplyr GAMBLR.helpers
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#'
#' #Return a plot for one sample, with default parameters.
#' fancy_cnbar(this_sample_id = "DOHH-2")
#'
fancy_cnbar = function(this_sample_id,
                       seg_data,
                       seg_path = NULL,
                       chrom_col = 2,
                       start_col = 3,
                       end_col = 4,
                       cn_col = 7,
                       plot_title = paste0(this_sample_id),
                       plot_subtitle = "n CNV Segments (barplots, left y-axis), n Affected bases for each CN state",
                       chr_select = paste0("chr", c(1:22)),
                       cutoff = 15,
                       include_cn2 = TRUE,
                       this_seq_type = "genome") {

  if(!missing(seg_data)){
    seg = seg_data
    seg = as.data.frame(seg)
    colnames(seg)[chrom_col] = "chrom"
    colnames(seg)[start_col] = "start"
    colnames(seg)[end_col] = "end"
    colnames(seg)[cn_col] = "CN"

  }else if (!is.null(seg_path)){
    seg = read.table(seg_path, sep = "\t", header = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
    seg = as.data.frame(seg)
    colnames(seg)[chrom_col] = "chrom"
    colnames(seg)[start_col] = "start"
    colnames(seg)[end_col] = "end"
    colnames(seg)[cn_col] = "CN"
  }

  #get seg data for a specific sample.
  if(missing(seg_data) && is.null(seg_path)){
    seg = get_sample_cn_segments(these_sample_ids = this_sample_id,
                                 streamlined = FALSE,
                                 this_seq_type = this_seq_type
    )
  }

  #add chr prefix if missing
  if(!str_detect(seg$chrom, "chr")[2]){
    seg = mutate(seg, chrom = paste0("chr", chrom))
  }

  #read maf into R and select relevant variables and transformt to factor.
  seg_df = dplyr::select(seg, chrom, start, end, CN) %>%
    mutate_at(vars(chrom), list(factor))

  #subsetting maf based on user-defined parameters
  seg_df = seg_df[seg_df$chrom %in% chr_select, ]

  #transform data type
  seg_df$CN = as.factor(seg_df$CN)

  #count levels of factor
  if(include_cn2){
    cn_states = c(0:cutoff)
    cns_count = dplyr::filter(seg_df, CN %in% cn_states) %>%
      group_by(CN) %>%
      summarize(count = n())
  }else{
    cn_states = c(0:1, 3:cutoff)
    cns_count = dplyr::filter(seg_df, CN %in% cn_states) %>%
      group_by(CN) %>%
      summarize(count = n())
  }

  cns_count$Type = paste0(cns_count$CN)
  cns_count = dplyr::select(cns_count, count, Type)

  #compute lenght of cn segments and transform for plotting
  l_cn_seg = seg
  l_cn_seg$lenght = l_cn_seg$end - l_cn_seg$start
  l_cn_seg$CN = as.factor(l_cn_seg$CN)
  l_cn_seg = dplyr::filter(l_cn_seg, CN %in% cn_states)

  if(!include_cn2){
    l_cn_seg = dplyr::filter(l_cn_seg, CN != 2)
  }

  cn_seg_lenghts = aggregate(l_cn_seg$lenght, list(l_cn_seg$CN), sum)
  colnames(cn_seg_lenghts) = c("CN", "lenght")

  joined_cn = cbind(cns_count, cn_seg_lenghts) %>%
    as.data.frame() %>%
    dplyr::select(CN, count, lenght)

  #get levels of cn states for plotting
  cn_levels = cns_count$Type

  #plot
  p = ggplot(joined_cn, aes(x = CN)) +
    geom_segment(aes(y = 1, yend = lenght/500000, x = CN, xend = CN)) +
    geom_point(aes(y = lenght/500000), colour = "#E6B315", size = 3, group = 2) +
    geom_bar(aes(y = count, fill = CN), position = "stack", stat = "identity") +
    scale_y_log10(limits = c(1, max(joined_cn$count) + 5000), sec.axis = sec_axis(~.*500000, name = "Nucleotides (n)")) +
    labs(title = plot_title, subtitle = plot_subtitle, x = "CN States", y = "CN Segments (n)", fill = "Legend") +
    scale_fill_manual(values = GAMBLR.helpers::get_gambl_colours("copy_number")) +
    scale_x_discrete(limits = cn_levels) +
    geom_text(aes(x = CN, y = count, label = count), colour = "#000000", size = 5, position = position_stack(vjust = 0.5)) +
    theme_Morons() +
    theme(legend.position = "none")

  return(p)
}
