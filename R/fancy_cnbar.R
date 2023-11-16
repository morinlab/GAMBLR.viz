#' @title Copy Number states barplot
#'
#' @description Generate a bar plot visualizing sample-specific copy number states and affected bases for each CN segment.
#'
#' @details `fancy_cnbar` visualizes copy number (CN) states on sample-level. Similarly to other fancy_x_plots this function
#' accepts either a sample ID, for which the function will get copy number states with [GAMBLR::get_sample_cn_segments]. The function
#' can also accept an already loaded seq file given to the `seq_data` parameter. It can also load a seq file with the `seq_path`
#' parameter. If the user calls either `seq_data` or `seq_path`, there are a collection of parameters available for specifying
#' the relevant columns in the given data frame (`chrom_col`, `starat_col`, `end_col`, `cn_col`). It is also possible to
#' restrict the returned plot to any given chromosome. This is done with the `chr_select` parameter (default is all autosomes).
#' For further control of the returned plot, it is also possible to set the threshold for maximum CN states to be returned (default is 15).
#' With `include_cn2` (Boolean) the user can control if CN segments = 2 should be added to the plot, default is TRUE.
#' The user can also control the annotations of the returned plot with `plot_title` and `plot_subtitle`. Lastly,
#' This function also computes the number of affected bases for each copy number state and plots these values on a secondary y-axis (right),
#' useful for overviewing the extent of each copy number state, in the context of the full genome.
#'
#' @param this_sample_id Sample to be plotted.
#' @param seq_data Optional parameter with copy number df already loaded into R.
#' @param seq_path Optional parameter with path to external cn file.
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
#' @import ggplot2 dplyr cowplot
#' @export
#'
#' @examples
#' #Return a plot for one sample, with default parameters.
#' \dontrun{
#' fancy_cnbar(this_sample_id = "DOHH-2")
#' }
#'
fancy_cnbar = function(this_sample_id,
                       seq_data,
                       seq_path = NULL,
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

  if(!missing(seq_data)){
    seq = seq_data
    seq = as.data.frame(seq)
    colnames(seq)[chrom_col] = "chrom"
    colnames(seq)[start_col] = "start"
    colnames(seq)[end_col] = "end"
    colnames(seq)[cn_col] = "CN"

  }else if (!is.null(seq_path)){
    seq = read.table(seq_path, sep = "\t", header = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
    seq = as.data.frame(seq)
    colnames(seq)[chrom_col] = "chrom"
    colnames(seq)[start_col] = "start"
    colnames(seq)[end_col] = "end"
    colnames(seq)[cn_col] = "CN"
  }

  #get maf data for a specific sample.
  if(missing(seq_data) && is.null(seq_path)){
    seq = get_sample_cn_segments(these_sample_ids = this_sample_id,
                                 streamlined = FALSE,
                                 this_seq_type = this_seq_type
    )
  }

  #add chr prefix if missing
  if(!str_detect(seq$chrom, "chr")[2]){
    seq = mutate(seq, chrom = paste0("chr", chrom))
  }

  #read maf into R and select relevant variables and transformt to factor.
  seq_df = dplyr::select(seq, chrom, start, end, CN) %>%
    mutate_at(vars(chrom), list(factor))

  #subsetting maf based on user-defined parameters
  seq_df = seq_df[seq_df$chrom %in% chr_select, ]

  #transform data type
  seq_df$CN = as.factor(seq_df$CN)

  #count levels of factor
  if(include_cn2){
    cn_states = c(0:cutoff)
    cns_count = dplyr::filter(seq_df, CN %in% cn_states) %>%
      group_by(CN) %>%
      summarize(count = n())
  }else{
    cn_states = c(0:1, 3:cutoff)
    cns_count = dplyr::filter(seq_df, CN %in% cn_states) %>%
      group_by(CN) %>%
      summarize(count = n())
  }

  cns_count$Type = paste0(cns_count$CN)
  cns_count = dplyr::select(cns_count, count, Type)

  #compute lenght of cn segments and transform for plotting
  l_cn_seg = seq
  l_cn_seg$lenght = l_cn_seg$end - l_cn_seg$start
  l_cn_seg$CN = as.factor(l_cn_seg$CN)
  l_cn_seg = dplyr::filter(l_cn_seg, CN %in% cn_states)

  if(!include_cn2){
    l_cn_seg = dplyr::filter(l_cn_seg, CN != 2)
  }

  cn_seq_lenghts = aggregate(l_cn_seg$lenght, list(l_cn_seg$CN), sum)
  colnames(cn_seq_lenghts) = c("CN", "lenght")

  joined_cn = cbind(cns_count, cn_seq_lenghts) %>%
    as.data.frame() %>%
    dplyr::select(CN, count, lenght)

  #get levels of cn states for plotting
  cn_levels = cns_count$Type

  #plot
  p = ggplot(joined_cn, aes(x = CN)) +
    geom_segment(aes(y = 1, yend = lenght/500000, x = CN, xend = CN)) +
    geom_point(aes(y = lenght/500000), colour = "#E6B315", size = 3, group = 2) +
    geom_bar(aes(y = count, fill = CN, label = count), position = "stack", stat = "identity") +
    scale_y_log10(limits = c(1, max(joined_cn$count) + 5000), sec.axis = sec_axis(~.*500000, name = "Nucleotides (n)")) +
    labs(title = plot_title, subtitle = plot_subtitle, x = "CN States", y = "CN Segments (n)", fill = "Legend") +
    scale_fill_manual(values = GAMBLR.helpers::get_gambl_colours("copy_number")) +
    scale_x_discrete(limits = cn_levels) +
    geom_text(aes(x = CN, y = count, label = count), colour = "#000000", size = 5, position = position_stack(vjust = 0.5)) +
    theme_cowplot() +
    theme(legend.position = "none")

  return(p)
}
