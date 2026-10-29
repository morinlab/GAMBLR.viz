#' @title Plot a sample's Battenberg allele-specific copy number profile
#'
#' @description A quick visual assessment of one sample's Battenberg fit:
#' the allele-specific copy number profile with LOH annotated, the observed
#' segment BAF against the BAF expected for each called state, and a state
#' map placing every segment by its observed BAF and estimated total copy
#' number over the positions expected for integer allele states at the
#' fitted purity. Shows all autosomes or one chromosome.
#'
#' @details Takes the output of [GAMBLR.utils::summarize_battenberg], whose
#' `segments` table holds each segment's observed BAF (major allele
#' frequency), total copy number estimated from the segment's LogR at the
#' fitted purity and ploidy (`ntot`; not rounded to an integer state), clonal allele state and
#' LOH class, and whose `summary` holds purity and ploidy. LOH classes are
#' relative to the sample's modal total copy number, as in
#' `summarize_battenberg`: copy-neutral LOH (`loh_neutral`), LOH with copy
#' loss (`loh_loss`) or gain (`loh_gain`), and homozygous deletion.
#'
#' The expected BAF of a state with `a` copies of the major allele and `b` of
#' the minor one, at purity `p`, is `(p * a + (1 - p)) / (p * (a + b) + 2 *
#' (1 - p))`. Segments far from the expected BAF of their state, or far from
#' every grid point in the state map, suggest a poor fit, subclonality or a
#' different purity/ploidy solution.
#'
#' @param battenberg Output of [GAMBLR.utils::summarize_battenberg] (a list
#' with `segments` and `summary`), covering at least this sample.
#' @param these_samples_metadata Optional one-row metadata data frame; its
#' sample_id selects the sample. Alternatively use `sample`.
#' @param sample Sample name as in `battenberg$summary$sample`. Not needed if
#' `battenberg` holds a single sample or `these_samples_metadata` is given.
#' @param chromosome Optional chromosome(s) to show, e.g. "1", "chr1" or
#' `c("1", "6", "17")`. One chromosome is drawn with a Mb axis; several are
#' laid out side by side, as for the whole genome. Default all autosomes.
#' @param max_cn Copy number at which the y axes are capped. Default 6.
#' @param min_segment_size Segments smaller than this (bp) are left out of
#' the state map (they still appear in the tracks). Default 1e5.
#' @param panels Which panels to draw: any of "profile", "baf" and
#' "state_map". Default all three.
#' @param return_plots Set to TRUE to return the individual ggplot objects
#' (a named list) instead of the combined figure.
#'
#' @return A combined plot (cowplot), or with `return_plots = TRUE` a named
#' list of ggplot objects.
#'
#' @import ggplot2 cowplot
#' @export
#'
#' @examples
#' \dontrun{
#' meta_row = dplyr::filter(get_gambl_metadata(), sample_id == "SP59324")
#' bb = GAMBLR.utils::summarize_battenberg(GAMBLR.results::get_battenberg(meta_row))
#' plot_sample_battenberg(bb, meta_row)
#' plot_sample_battenberg(bb, meta_row, chromosome = "1")
#'
#' # local files
#' bb = GAMBLR.utils::summarize_battenberg(battenberg_dir = "battenberg_check")
#' plot_sample_battenberg(bb, sample = "06-16716T", chromosome = "1")
#' }
plot_sample_battenberg = function(battenberg,
                           these_samples_metadata = NULL,
                           sample = NULL,
                           chromosome = NULL,
                           max_cn = 6,
                           min_segment_size = 1e5,
                           panels = c("profile", "baf", "state_map"),
                           return_plots = FALSE){
  if(!is.list(battenberg) || !all(c("segments", "summary") %in% names(battenberg)) ||
     is.null(battenberg$segments)){
    stop("battenberg must be the output of GAMBLR.utils::summarize_battenberg() (with a segments table)")
  }
  panels = match.arg(panels, c("profile", "baf", "state_map"), several.ok = TRUE)
  if(!is.null(these_samples_metadata)){
    if(nrow(these_samples_metadata) != 1){
      stop("these_samples_metadata must have one row (one sample)")
    }
    sample = as.character(these_samples_metadata$sample_id)
  }
  available = unique(battenberg$summary$sample)
  if(is.null(sample)){
    if(length(available) != 1){
      stop("battenberg holds several samples; choose one with sample or these_samples_metadata")
    }
    sample = available
  }
  if(!sample %in% available){
    stop("sample ", sample, " is not in battenberg")
  }
  info = battenberg$summary[battenberg$summary$sample == sample, ][1, ]
  seg = battenberg$segments[battenberg$segments$sample == sample, ]
  if(!is.null(chromosome)){
    chromosome = unique(gsub("^chr", "", as.character(chromosome)))
    absent = setdiff(chromosome, seg$chrom)
    if(length(absent) == length(chromosome)){
      stop("no segments on chromosome ", paste(chromosome, collapse = ", "))
    }
    if(length(absent) > 0){
      warning("no segments on chromosome ", paste(absent, collapse = ", "))
    }
    seg = seg[seg$chrom %in% chromosome, ]
  }
  single_chromosome = length(unique(seg$chrom)) == 1 && !is.null(chromosome)
  purity = info$purity
  if(is.na(purity)){
    warning("no purity for ", sample, "; expected BAFs assume purity 1")
    purity = 1
  }

  # x positions: Mb on one chromosome, cumulative across several
  if(single_chromosome){
    seg$x0 = seg$start / 1e6
    seg$x1 = seg$end / 1e6
    x_lab = paste0("chr", seg$chrom[1], " (Mb)")
    chrom_breaks = NULL
  }else{
    chroms = as.character(sort(unique(as.integer(seg$chrom))))
    chrom_len = tapply(seg$end, seg$chrom, max)[chroms]
    offset = c(0, cumsum(as.numeric(chrom_len)))[seq_along(chroms)]
    names(offset) = chroms
    seg$x0 = (seg$start + offset[seg$chrom]) / 1e6
    seg$x1 = (seg$end + offset[seg$chrom]) / 1e6
    x_lab = NULL
    chrom_breaks = data.frame(chrom = chroms,
                              start = offset / 1e6,
                              mid = (offset + chrom_len / 2) / 1e6)
  }

  state_colours = c(loh_neutral = "#7B3294", loh_loss = "#2166AC",
                    loh_gain = "#E66101", homdel = "black",
                    imbalance_gain = "#F4A582", imbalance_neutral = "#B8B8B8",
                    imbalance_loss = "#92C5DE", balanced = "grey70")
  loh_states = c("loh_neutral", "loh_loss", "loh_gain", "homdel")
  seg$state = factor(seg$state, levels = names(state_colours))
  expected_baf = function(a, b){
    (purity * a + (1 - purity)) / (purity * (a + b) + 2 * (1 - purity))
  }
  seg$baf_major = pmax(seg$BAF, 1 - seg$BAF)
  seg$baf_expected = ifelse(seg$total > 0, expected_baf(seg$major, seg$minor), NA)
  cap = function(v) pmin(v, max_cn)

  genome_axis = function(p){
    if(is.null(chrom_breaks)){
      return(p + labs(x = x_lab))
    }
    p + geom_vline(xintercept = chrom_breaks$start[-1], colour = "grey85",
                   linewidth = 0.3) +
      scale_x_continuous(breaks = chrom_breaks$mid, labels = chrom_breaks$chrom,
                         expand = c(0.005, 0)) +
      labs(x = NULL)
  }
  loh_shading = geom_rect(data = seg[seg$state %in% loh_states, ],
                          aes(xmin = x0, xmax = x1, ymin = -Inf, ymax = Inf,
                              fill = state),
                          alpha = 0.18, inherit.aes = FALSE)
  plots = list()

  if("profile" %in% panels){
    plots$profile = genome_axis(
      ggplot(seg) + loh_shading +
        geom_segment(aes(x = x0, xend = x1, y = cap(ntot), yend = cap(ntot)),
                     colour = "grey55", linewidth = 0.6) +
        geom_segment(aes(x = x0, xend = x1, y = cap(major) + 0.06,
                         yend = cap(major) + 0.06), colour = "#D55E00",
                     linewidth = 1.1) +
        geom_segment(aes(x = x0, xend = x1, y = cap(minor) - 0.06,
                         yend = cap(minor) - 0.06), colour = "#0072B2",
                     linewidth = 1.1) +
        scale_fill_manual(values = state_colours, drop = TRUE, name = "LOH") +
        scale_y_continuous(limits = c(-0.2, max_cn + 0.3), breaks = 0:max_cn) +
        labs(y = "copy number",
             subtitle = "major (orange) and minor (blue) allele copies; total copy number estimated from LogR (grey)") +
        theme_classic() +
        theme(legend.position = "top", plot.subtitle = element_text(size = 8)))
  }
  if("baf" %in% panels){
    plots$baf = genome_axis(
      ggplot(seg) + loh_shading +
        geom_segment(aes(x = x0, xend = x1, y = baf_expected, yend = baf_expected),
                     colour = "#E08214", linetype = "22", linewidth = 0.6,
                     na.rm = TRUE) +
        geom_segment(aes(x = x0, xend = x1, y = baf_major, yend = baf_major),
                     colour = "black", linewidth = 0.8, na.rm = TRUE) +
        scale_fill_manual(values = state_colours, drop = TRUE, guide = "none") +
        scale_y_continuous(limits = c(0.5, 1)) +
        labs(y = "BAF (major allele)",
             subtitle = "observed (black) vs expected (orange)") +
        theme_classic() +
        theme(plot.subtitle = element_text(size = 8)))
  }
  if("state_map" %in% panels){
    grid = expand.grid(a = 0:max_cn, b = 0:max_cn)
    grid = grid[grid$a >= grid$b & grid$a + grid$b <= max_cn &
                  grid$a + grid$b > 0, ]
    grid$baf = expected_baf(grid$a, grid$b)
    grid$total = grid$a + grid$b
    grid$label = paste0(grid$a, "+", grid$b)
    map_seg = seg[seg$end - seg$start >= min_segment_size & !is.na(seg$ntot), ]
    plots$state_map = ggplot() +
      geom_point(data = grid, aes(x = baf, y = total), shape = 3,
                 colour = "grey60", size = 2) +
      geom_text(data = grid, aes(x = baf, y = total, label = label),
                colour = "grey50", size = 2.4, vjust = -0.9) +
      geom_point(data = map_seg,
                 aes(x = baf_major, y = cap(ntot), colour = state,
                     size = (end - start) / 1e6),
                 alpha = 0.6) +
      scale_colour_manual(values = state_colours, drop = TRUE, name = "state") +
      scale_size_area(max_size = 8, name = "Mb") +
      scale_x_continuous(limits = c(0.45, 1.02)) +
      scale_y_continuous(limits = c(-0.2, max_cn + 0.5), breaks = 0:max_cn) +
      labs(x = "BAF (major allele)", y = "total copy number (from LogR)",
           subtitle = "segments vs integer states (major+minor) at this purity") +
      theme_classic() +
      theme(plot.subtitle = element_text(size = 8))
  }

  if(return_plots){
    return(plots)
  }
  title = paste0(sample,
                 if(!is.null(chromosome))
                   paste0(", chr", paste(sort(as.integer(unique(seg$chrom))), collapse = ", ")) else "",
                 "   purity ", format(round(purity, 2), nsmall = 2),
                 ", ploidy ", format(round(info$ploidy, 2), nsmall = 2),
                 if(isTRUE(info$doubled)) " (doubled)" else "",
                 "   copy-neutral LOH ", round(100 * info$loh_neutral_fraction, 1),
                 "%, LOH with loss ", round(100 * info$loh_loss_fraction, 1), "%")
  tracks = plots[intersect(c("profile", "baf"), names(plots))]
  left = if(length(tracks) > 0){
    cowplot::plot_grid(plotlist = tracks, ncol = 1, align = "v", axis = "lr",
                       rel_heights = if(length(tracks) == 2) c(1.4, 1) else 1)
  }else NULL
  body = if(!is.null(left) && !is.null(plots$state_map)){
    cowplot::plot_grid(left, plots$state_map, nrow = 1, rel_widths = c(2.2, 1))
  }else if(!is.null(left)) left else plots$state_map
  cowplot::plot_grid(cowplot::ggdraw() +
                       cowplot::draw_label(title, x = 0.01, hjust = 0, size = 11),
                     body, ncol = 1, rel_heights = c(0.06, 1))
}
