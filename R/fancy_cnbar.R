#' @title Copy Number states barplot
#'
#' @description Generate a bar plot visualizing sample-specific copy number states and affected bases for each CN segment.
#'
#' @details `fancy_cnbar` visualizes copy number (CN) states on sample-level. This function expects an already loaded SEG file with the `this_seg` parameter.
#' Such a seg file can be retrieved with [GAMBLR.results::get_sample_cn_segments]. It is also possible to specify the absolute path to such a SEG file with `this_seg_path`.
#' The returned plot can also be restricted to any given chromosome. This is done with the `chr_select` parameter (default is all autosomes).
#' For further control of the returned plot, it is also possible to set the threshold for maximum CN states to be returned (default is 15).
#' With `include_cn2` (Boolean) the user can control if CN segments = 2 should be added to the plot, default is TRUE.
#' The user can also control the annotations of the returned plot with `plot_title` and `plot_subtitle`. Lastly,
#' This function also computes the number of affected bases for each copy number state and plots these values on a secondary y-axis (right),
#' useful for overview the extent of each copy number state, in the context of the full genome.
#'
#' @param this_seg An already loaded SEG file. This function is compatible with the return of [GAMBLR.results::get_sample_cn_segments].
#' @param this_seg_path An absolute path to a SEG file. This function is compatible with the return of [GAMBLR.results::get_sample_cn_segments].
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Subtitle for created plot.
#' @param chr_select Vector of chromosomes to be included in plot, defaults to autosomes.
#' @param cutoff Set threshold for maximum CN state to be retrieved.
#' @param include_cn2 Optional Boolean statement for including CN = 2 states in plot.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import ggplot2 dplyr cowplot
#' @export
#'
#' @examples
#' #get data
#' dohh2_seg = GAMBLR.data::sample_data$grch37$seg %>% dplyr::filter(ID == "DOHH-2")
#' 
#' #build plot
#' fancy_cnbar(this_seg = dohh2_seg, include_cn2 = FALSE)
#'
fancy_cnbar = function(this_seg,
                       this_seg_path = NULL,
                       plot_title = paste0("Copy Number Segments"),
                       plot_subtitle = "n CNV Segments (barplots, left y-axis), n Affected bases for each CN state",
                       chr_select = paste0("chr", c(1:22)),
                       cutoff = 15,
                       include_cn2 = TRUE){

  if(!missing(this_seg)){
    seg = this_seg

  }else if (!is.null(this_seg_path)){
    seg = read.table(this_seg_path, sep = "\t", header = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  }

  if(missing(this_seg) && is.null(this_seg_path)){
    stop("Please provide an already loaded SEG file with the `seg_data` parameter. 
  If you have GAMBLR.results installed (and GSC access), you can call `get_sample_cn_segments` to retrieve such a file.
  As an alternative, an absolute path to a SEG file is also accepted by `seg_path`.")
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
  cns_count = dplyr::select(cns_count, count, Type, CN)

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

  joined_cn = merge(cns_count, cn_seg_lenghts) %>%
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
    scale_fill_manual(values = get_gambl_colours("copy_number")) +
    scale_x_discrete(limits = cn_levels) +
    geom_text(aes(x = CN, y = count, label = count), colour = "#000000", size = 5, position = position_stack(vjust = 0.5)) +
    theme_cowplot() +
    theme(legend.position = "none")

  return(p)
}