#' @title Pretty Cumulative CN Alteration Plot
#'
#' @description Genome-wide frequency of copy number gain/loss across a
#' cohort -- the per-bin "Gain"/"loss" bars [pretty_CN_heatmap()] draws as a
#' heatmap annotation, pulled out into its own standalone plot. This is the
#' complement to [pretty_CN_cohort_heatmap()]: same matrix input (e.g. from
#' [segments_to_exact_cn_matrix()] or [segmented_data_to_cn_matrix()]), and
#' the same column-index x-axis and chromosome strip, so the two line up if
#' stacked.
#'
#' @details For each bin (column), the fraction (or count, see
#' `as_fraction`) of samples with `CN > gain_threshold` is drawn as a red
#' area above zero, and the fraction/count with `CN < loss_threshold` as a
#' blue area below zero -- the standard genome-wide CN alteration frequency
#' plot (as produced by, e.g., GISTIC).
#'
#' `gain_threshold`/`loss_threshold` default to 2.5/1.5, not the exact `> 2`/
#' `< 2` [pretty_CN_heatmap()] uses. A matrix of continuous (unrounded) CN
#' values -- the default output of [segments_to_exact_cn_matrix()] -- will
#' have plenty of samples sitting just barely off 2 from measurement noise
#' alone (more so still on a ploidy-adjusted matrix, where "2" means a
#' sample's own recentred baseline rather than a literal integer copy
#' number); counting every such marginal crossing makes this plot much
#' noisier than it should be and buries the loci with real, confident
#' alterations. Narrow the thresholds back toward 2/2 to match
#' [pretty_CN_heatmap()]'s convention exactly, or pass a `rounded = TRUE`
#' matrix instead, if that's what you want.
#'
#' Bins are plotted by column index, not genomic coordinate -- the same
#' simplification [pretty_CN_cohort_heatmap()] makes -- so that stacking the
#' two keeps every bin aligned on the x-axis regardless of its true width.
#'
#' @param cn_state_matrix A sample x region CN matrix, e.g. from
#' [segments_to_exact_cn_matrix()] or [segmented_data_to_cn_matrix()]. Column
#' names must be `"chrom:start-end"`, in genomic order.
#' @param gain_threshold,loss_threshold A sample counts as gained at a bin if
#' its CN there is `> gain_threshold`, and lost if `< loss_threshold`.
#' Defaults 2.5/1.5 -- wider than a strict `> 2`/`< 2` call, so that samples
#' sitting only marginally off 2 (noise, not a confident alteration) don't
#' contribute. See Details.
#' @param as_fraction Set to FALSE to plot (and return) raw sample counts
#' instead of the fraction of the cohort. Default TRUE.
#' @param chrom_list Optional vector of chromosomes to restrict the plot to.
#' Default: every chromosome present in `cn_state_matrix`.
#' @param gain_colour,loss_colour Fill colours for the gain/loss areas.
#' Defaults `"red"`/`"blue"`.
#' @param title Plot title. Default `"Cumulative Copy Number Alterations"`.
#' @param chromosome_strip Set to FALSE to omit the chromosome boundary strip
#' below the plot. Default TRUE.
#' @param return_data Set to TRUE to get the underlying per-bin data frame
#' back instead of a plot.
#'
#' @return A ggplot/cowplot object, or a data frame if `return_data = TRUE`.
#'
#' @import ggplot2 dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' cn_mat = segments_to_exact_cn_matrix(seg_data = all_segments,
#'                                      these_samples_metadata = meta)
#' pretty_CN_cumulative_plot(cn_mat)
#' pretty_CN_cumulative_plot(cn_mat, chrom_list = c("chr17", "chr18"))
#' }
pretty_CN_cumulative_plot = function(cn_state_matrix,
                                     gain_threshold = 2.5,
                                     loss_threshold = 1.5,
                                     as_fraction = TRUE,
                                     chrom_list,
                                     gain_colour = "red",
                                     loss_colour = "blue",
                                     title = "Cumulative Copy Number Alterations",
                                     chromosome_strip = TRUE,
                                     return_data = FALSE){

  column_chromosome = sub(":.*", "", colnames(cn_state_matrix))

  if(!missing(chrom_list)){
    keep_cols = which(column_chromosome %in% chrom_list)
    if(length(keep_cols) == 0){
      stop("None of chrom_list was found among cn_state_matrix's columns.")
    }
    cn_state_matrix = cn_state_matrix[, keep_cols, drop = FALSE]
    column_chromosome = column_chromosome[keep_cols]
  }

  n_samples = nrow(cn_state_matrix)
  n_bins = ncol(cn_state_matrix)
  cumulative_gain = colSums(cn_state_matrix > gain_threshold, na.rm = TRUE)
  cumulative_loss = colSums(cn_state_matrix < loss_threshold, na.rm = TRUE)

  if(as_fraction){
    cumulative_gain = cumulative_gain / n_samples
    cumulative_loss = cumulative_loss / n_samples
  }

  plot_df = data.frame(
    bin_index = seq_len(n_bins),
    bin = colnames(cn_state_matrix),
    chromosome = column_chromosome,
    cumulative_gain = unname(cumulative_gain),
    cumulative_loss = unname(cumulative_loss)
  )

  if(return_data){
    return(plot_df)
  }

  y_lab = if(as_fraction) "Fraction of cohort" else "Number of samples"

  main_plot = ggplot2::ggplot(plot_df, ggplot2::aes(x = bin_index)) +
    ggplot2::geom_area(ggplot2::aes(y = cumulative_gain), fill = gain_colour, colour = gain_colour) +
    ggplot2::geom_area(ggplot2::aes(y = -cumulative_loss), fill = loss_colour, colour = loss_colour) +
    ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0.5, n_bins + 0.5)) +
    ggplot2::labs(title = title, x = NULL, y = y_lab) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(b = 0)
    )

  if(!chromosome_strip){
    return(main_plot)
  }

  chrom_col = c("chr1"="#555FAB", "chr2"="#CE3D31", "chr3"="#749B58",
                "chr4"="#F0E584", "chr5"="#476A85", "chr6"="#BA6338",
                "chr7"="#5CB1DD", "chr8"="#7F2368", "chr9"="#77C269",
                "chr10"="#D595A6", "chr11"="#934823", "chr12"="#857B8E",
                "chr13"="#C85328", "chr14"="#D58F5B", "chr15"="#7A65A7",
                "chr16"="#E3AF69", "chr17"="#3C1C54", "chr18"="#CEDEB7",
                "chr19"="#612B79", "chr20"="#AF2064", "chr21"="#E6C66F",
                "chr22"="#5B665D", "chrX"="#CA992C", "chrY"="white")
  if(!grepl("chr", column_chromosome[1])){
    names(chrom_col) = gsub("chr", "", names(chrom_col))
  }

  chrom_rle = rle(column_chromosome)
  chrom_end = cumsum(chrom_rle$lengths)
  chrom_start = c(1, utils::head(chrom_end, -1) + 1)
  strip_df = data.frame(
    xmin = chrom_start - 0.5,
    xmax = chrom_end + 0.5,
    xmid = (chrom_start + chrom_end) / 2,
    chromosome = chrom_rle$values
  )

  strip_plot = ggplot2::ggplot(strip_df) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = chromosome),
      colour = NA, show.legend = FALSE
    ) +
    ggplot2::geom_text(ggplot2::aes(x = xmid, y = 0.5, label = chromosome), size = 2) +
    ggplot2::scale_fill_manual(values = chrom_col) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0.5, n_bins + 0.5)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 0))

  cowplot::plot_grid(
    main_plot, strip_plot, ncol = 1, align = "v", rel_heights = c(8, 1)
  )
}
