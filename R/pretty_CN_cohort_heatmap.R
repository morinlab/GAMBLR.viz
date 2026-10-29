#' @title Pretty Cohort CN Heatmap
#'
#' @description Fast-to-redraw cohort-level copy number overview from a
#' precomputed matrix (e.g. from [segments_to_exact_cn_matrix()]). This is
#' the complement to [pretty_CN_heatmap()] for when you want the heatmap
#' driven by actual segment values (via a breakpoint-resolution matrix)
#' rather than a fixed bin count, at a scale (hundreds to low thousands of
#' samples) where the slow part -- turning segments into a matrix -- should
#' happen once, not on every plot tweak.
#'
#' @details This function draws the heatmap body with raw `grid` graphics
#' (`grid::grid.raster()` on a precomputed colour matrix) instead of
#' `ComplexHeatmap::Heatmap()`. That's a deliberate, benchmarked choice, not
#' a style preference: at cohort scale (tested at 1500 samples x ~60,000
#' bins) `Heatmap(..., use_raster = TRUE)` still visits every cell via
#' `grid`-level drawing calls internally to build its raster and took
#' 115-177s to `draw()` in testing, whereas building a colour matrix and
#' handing it to `grid.raster()` directly (which treats it as a literal
#' pixel array) took ~7s end to end for the same data -- `use_raster`
#' saves output file size/vector complexity, it does not avoid visiting
#' every cell at render time. Colour mapping itself also matters at this
#' scale: `circlize::colorRamp2()` does LAB/HCL interpolation per call
#' (~39s on that same matrix), so this function precomputes a small
#' (`n_color_levels`) palette once and buckets values into it with
#' `findInterval()` instead (~11x faster; colours differ by at most 1-2 RGB
#' units, imperceptible at heatmap scale). The legend is still built with
#' `ComplexHeatmap::Legend()`, which is unaffected by matrix size.
#'
#' Row (sample) and column (region) order are not altered -- this function
#' assumes your matrix is already in the order you want (genomic order for
#' columns is guaranteed if it came from [segments_to_exact_cn_matrix()] or
#' [segmented_data_to_cn_matrix()]; see `sample_order` for rows). There is
#' no clustering option here: hierarchical clustering at this scale is its
#' own potentially slow step, deliberately decoupled so it doesn't have to
#' be recomputed on every redraw either -- compute an order yourself (e.g.
#' from a subsample of columns, or from [pretty_CN_heatmap()]'s own
#' clustering on a coarser matrix of the same samples) and pass it in.
#'
#' @param cn_state_matrix A sample x region CN matrix, e.g. from
#' [segments_to_exact_cn_matrix()] (or [segmented_data_to_cn_matrix()]).
#' Column names must be `"chrom:start-end"`, in genomic order.
#' @param these_samples_metadata Optional metadata table (with `sample_id`)
#' for the left annotation strip. If missing, no annotation is drawn.
#' @param metadataColumns Column names from `these_samples_metadata` to
#' show as row annotation colour strips. Default `c("pathology")`.
#' @param sample_order Optional vector of `sample_id`s (a permutation of
#' `rownames(cn_state_matrix)`) giving the row order to plot. Default: the
#' matrix's existing row order.
#' @param chrom_list Optional vector of chromosomes to restrict the plot to
#' (must match `colnames(cn_state_matrix)`'s own `chr`-prefix convention
#' exactly, e.g. `c("chr8", "chr18")` -- no normalization is attempted).
#' Filtering happens before everything else here, so it also narrows the
#' `cumulative_plot` panel, gene-label matching and the chromosome strip for
#' free. Default: every chromosome present in `cn_state_matrix`. Named to
#' match [pretty_CN_cumulative_plot()]'s and [plot_sample_circos()]'s
#' argument of the same name, rather than [pretty_CN_heatmap()]'s
#' `keep_these_chromosomes`.
#' @param col_breaks,col_colors `circlize::colorRamp2()`-style breakpoints
#' and colours for the CN colour scale. Defaults `c(1, 2, 5)` /
#' `c("blue", "white", "red")`, matching [pretty_CN_heatmap()].
#' @param n_color_levels Size of the precomputed colour lookup table. Higher
#' is smoother but slower to build. Default 256 (far more than visually
#' distinguishable).
#' @param cumulative_plot Set to `"top"` or `"bottom"` to draw a genome-wide
#' gain/loss frequency panel (as [pretty_CN_cumulative_plot()] computes)
#' directly above/below the heatmap body, column-aligned with it -- the
#' same pairing [pretty_CN_heatmap()] draws as a built-in bottom annotation.
#' [pretty_CN_cumulative_plot()] is called internally for its data only
#' (`return_data = TRUE`); the panel itself is drawn with the same raw
#' `grid` primitives as the rest of this function (not its ggplot output)
#' so it lines up exactly with the body and stays cheap to redraw -- an
#' embedded ggplot would reserve its own axis/margin space and throw off
#' column alignment. Default `"none"`. For a fully-labelled version with a
#' real y-axis, call [pretty_CN_cumulative_plot()] directly instead; this
#' embedded panel only labels its max value, not a full axis (see
#' `cumulative_plot_height`).
#' @param cumulative_plot_height Height in cm of the `cumulative_plot` panel.
#' Default 3.
#' @param cumulative_gain_threshold,cumulative_loss_threshold Passed to
#' [pretty_CN_cumulative_plot()]. Defaults 2.5/1.5 (see its documentation
#' for why not a strict 2/2).
#' @param cumulative_as_fraction Passed to [pretty_CN_cumulative_plot()] as
#' `as_fraction`. Default TRUE.
#' @param cumulative_gain_colour,cumulative_loss_colour Fill colours for the
#' `cumulative_plot` panel. Defaults `"red"`/`"blue"`.
#' @param labelGainGenes,labelLossGenes Vectors of Hugo gene symbols to label
#' directly on the `cumulative_plot` panel instead of the heatmap body --
#' `labelGainGenes` anchored to the gain (red) curve, `labelLossGenes` to
#' the loss (blue) curve, each with a connector line running from that
#' gene's actual bin value on the curve up/down to its label. This makes it
#' easy to see at a glance whether a known gene sits on a local peak/trough
#' or somewhere unremarkable, which a label on the heatmap body can't show.
#' Requires `cumulative_plot != "none"` (there's no curve to anchor to
#' otherwise); an error asks for it if you forget. Resolved the same way as
#' `labelTheseGenes` (see there for matching/warning behaviour), and shares
#' its `genome_build`/`bin_label_fontsize`/`bin_label_rotation`.
#' @param cumulative_label_headroom When `labelGainGenes`/`labelLossGenes`
#' are used, the `cumulative_plot` panel's y-range is expanded by this
#' factor beyond the data's own max, to leave room for labels above the
#' gain curve/below the loss curve without them overlapping it. Default 1.4.
#' Ignored with no cumulative gene labels.
#' @param cumulative_n_gain_peaks,cumulative_n_loss_peaks Set to label the
#' top N local peaks of the gain/loss curve automatically, instead of (or
#' alongside) naming specific genes via `labelGainGenes`/`labelLossGenes` --
#' useful for a first look at a new cohort, when you don't yet know which
#' genes are worth naming explicitly. Peaks are found greedily (tallest
#' first, skipping any candidate within `cumulative_peak_min_separation` of
#' one already picked, so one broad event doesn't get relabelled several
#' times) and must clear `cumulative_peak_min_height`. Each peak's bin is
#' reverse-looked-up against `GAMBLR.data::grch37_gene_coordinates` (or
#' `hg38_...`, per `genome_build`) for overlapping genes -- multiple genes
#' are joined (capped by `cumulative_peak_max_genes`), and a peak with no
#' overlapping gene falls back to a coordinate label (e.g. `"chr8:128.7Mb"`)
#' rather than being silently dropped. An auto-detected peak that lands on a
#' bin already labelled via `labelGainGenes`/`labelLossGenes` is skipped in
#' favour of the explicit one. Both default 0 (disabled); requires
#' `cumulative_plot != "none"` like the other cumulative-panel label
#' arguments. Unlike the rest of this function's labels, auto-detected ones
#' are not included in `return_data`'s output (they depend on the
#' cumulative panel's own data, which `return_data = TRUE` skips computing
#' to stay fast).
#' @param cumulative_peak_min_height Minimum peak height to be eligible for
#' auto-labelling, as a fraction (0-1) of that curve's own max -- e.g. the
#' default 0.1 means a gain peak must reach at least 10% of the tallest gain
#' peak's height to be considered, filtering out noise without needing to
#' know the data's absolute scale ahead of time.
#' @param cumulative_peak_min_separation Minimum distance between
#' auto-labelled peaks, as a fraction (0-1) of the total bin count. Default
#' 0.02 (2% of bins apart), auto-scaling with however many bins your matrix
#' has rather than assuming a fixed bin count.
#' @param cumulative_peak_max_genes Cap on how many gene symbols to join in
#' one auto-detected peak's label before collapsing the rest to `"+N"`.
#' Default 2.
#' @param labelTheseGenes A vector of Hugo gene symbols to label above the
#' heatmap body, as [pretty_CN_heatmap()]'s argument of the same name does:
#' each gene is mapped to the bin (column) it overlaps via `gene_to_region()`
#' and `genome_build`, and drawn as a rotated label with a connector line
#' down to that column. A gene with no match in `cn_state_matrix` (no
#' overlapping bin, or `gene_to_region()` doesn't recognize it) is skipped
#' with a `warning()` rather than erroring. Unlike [pretty_CN_heatmap()]
#' (which overlays labels on top of the body via `decorate_heatmap_body()`),
#' this draws them in their own row, directly above the body -- see Details
#' in [pretty_CN_cohort_heatmap()]'s own source for why an embedded ggplot
#' isn't used for panels like this; the same column-alignment reasoning
#' applies here, implemented natively in `grid` instead.
#' @param bin_labels Optional named vector/list (`bin name = label text`) to
#' label bins directly, skipping the `gene_to_region()` lookup -- as
#' [pretty_CN_heatmap()]'s argument of the same name. Bin names must match
#' `colnames(cn_state_matrix)` exactly. Merged with (and overridden by, on a
#' collision) whatever `labelTheseGenes` resolves to.
#' @param genome_build Genome build passed to `gene_to_region()` when
#' resolving `labelTheseGenes`. Default `"grch37"`.
#' @param bin_label_fontsize,bin_label_rotation Font size and rotation (in
#' degrees) for gene/bin labels. Defaults 5 and 45, matching
#' [pretty_CN_heatmap()].
#' @param bin_label_nudge How far up (as a fraction, 0-1, of the label row's
#' own height) a label's text sits above its connector line's base. Default
#' 0.15. Note this isn't directly comparable to [pretty_CN_heatmap()]'s
#' `bin_label_nudge`: that one is an npc offset above the heatmap body's own
#' viewport (values a little over 1), since labels there share the body's
#' coordinate space; here, labels get their own dedicated row (see
#' `gene_label_height_cm`), so the nudge is relative to that row instead.
#' @param gene_label_height_cm Height in cm of the gene/bin label row (drawn
#' only if `labelTheseGenes` or `bin_labels` is supplied). Default 2.
#' @param width,height Heatmap body size in cm. Default 20 x 10.
#' @param legend_title Title for the CN colour legend. Default
#' `"Copy\\nNumber"`.
#' @param show_legend Set to FALSE to omit the legend. Default TRUE.
#' @param verbose Default TRUE. Prints a timing breakdown for each step, so
#' you can see exactly where time goes while judging feasibility.
#' @param return_data Set to TRUE to get the precomputed colour matrix back
#' instead of drawing (e.g. to drop into your own custom grid layout).
#' @param newpage Set to FALSE to draw into the current device/viewport
#' instead of starting a new page (e.g. to compose with other grobs).
#'
#' @return Nothing (draws the plot), or a list with `colour_matrix`,
#' `chromosome_columns`, `data`, `labels`, `gain_labels` and `loss_labels`
#' if `return_data = TRUE`.
#'
#' @import ComplexHeatmap circlize dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' dlbcl_genome_meta = get_gambl_metadata() %>%
#'   dplyr::filter(pathology == "DLBCL", seq_type == "genome")
#' all_segments = get_cn_segments(these_samples_metadata = dlbcl_genome_meta)
#'
#' # slow step, run once
#' cn_mat = segments_to_exact_cn_matrix(
#'   seg_data = all_segments,
#'   these_samples_metadata = dlbcl_genome_meta,
#'   merge_tolerance = 100000
#' )
#'
#' # fast step, tweak freely
#' pretty_CN_cohort_heatmap(cn_mat, these_samples_metadata = dlbcl_genome_meta)
#' pretty_CN_cohort_heatmap(cn_mat, these_samples_metadata = dlbcl_genome_meta,
#'                          metadataColumns = c("pathology", "lymphgen"))
#' }
pretty_CN_cohort_heatmap = function(cn_state_matrix,
                                    these_samples_metadata,
                                    metadataColumns = c("pathology"),
                                    sample_order,
                                    chrom_list,
                                    col_breaks = c(1, 2, 5),
                                    col_colors = c("blue", "white", "red"),
                                    n_color_levels = 256,
                                    cumulative_plot = "none",
                                    cumulative_plot_height = 3,
                                    cumulative_gain_threshold = 2.5,
                                    cumulative_loss_threshold = 1.5,
                                    cumulative_as_fraction = TRUE,
                                    cumulative_gain_colour = "red",
                                    cumulative_loss_colour = "blue",
                                    labelGainGenes,
                                    labelLossGenes,
                                    cumulative_label_headroom = 1.4,
                                    cumulative_n_gain_peaks = 0,
                                    cumulative_n_loss_peaks = 0,
                                    cumulative_peak_min_height = 0.1,
                                    cumulative_peak_min_separation = 0.02,
                                    cumulative_peak_max_genes = 2,
                                    labelTheseGenes,
                                    bin_labels,
                                    genome_build = "grch37",
                                    bin_label_fontsize = 5,
                                    bin_label_rotation = 45,
                                    bin_label_nudge = 0.15,
                                    gene_label_height_cm = 2,
                                    width = 20,
                                    height = 10,
                                    legend_title = "Copy\nNumber",
                                    show_legend = TRUE,
                                    verbose = TRUE,
                                    return_data = FALSE,
                                    newpage = TRUE){
  t_start = Sys.time()
  tick = function(label, from){
    if(verbose){
      message(sprintf(
        "[pretty_CN_cohort_heatmap] %s: %.2fs",
        label, as.numeric(difftime(Sys.time(), from, units = "secs"))
      ))
    }
  }

  cumulative_plot = match.arg(cumulative_plot, c("none", "top", "bottom"))

  wants_cumulative_labels = !missing(labelGainGenes) || !missing(labelLossGenes) ||
    cumulative_n_gain_peaks > 0 || cumulative_n_loss_peaks > 0
  if(wants_cumulative_labels && cumulative_plot == "none"){
    stop(paste(
      "labelGainGenes/labelLossGenes/cumulative_n_gain_peaks/cumulative_n_loss_peaks",
      "are drawn on the cumulative_plot panel, so cumulative_plot must be 'top' or",
      "'bottom' to use them."
    ))
  }

  if(!missing(sample_order)){
    cn_state_matrix = cn_state_matrix[sample_order, , drop = FALSE]
  }

  if(!missing(chrom_list)){
    keep_cols = which(sub(":.*", "", colnames(cn_state_matrix)) %in% chrom_list)
    if(length(keep_cols) == 0){
      stop("None of chrom_list was found among cn_state_matrix's columns.")
    }
    cn_state_matrix = cn_state_matrix[, keep_cols, drop = FALSE]
  }

  column_chromosome = sub(":.*", "", colnames(cn_state_matrix))
  chrom_rle = rle(column_chromosome)
  chrom_end = cumsum(chrom_rle$lengths)
  chrom_start = c(1, utils::head(chrom_end, -1) + 1)
  chrom_mid = (chrom_start + chrom_end) / 2
  n_col = ncol(cn_state_matrix)

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

  t0 = Sys.time()
  cn_col_fun = circlize::colorRamp2(col_breaks, col_colors)
  lut_breaks = seq(min(col_breaks), max(col_breaks), length.out = n_color_levels)
  lut = cn_col_fun(lut_breaks)
  col_mat = matrix(
    lut[findInterval(cn_state_matrix, lut_breaks, all.inside = TRUE)],
    nrow = nrow(cn_state_matrix)
  )
  tick(sprintf("build colour matrix (%d x %d)", nrow(col_mat), ncol(col_mat)), t0)

  has_meta = !missing(these_samples_metadata)
  if(has_meta){
    t0 = Sys.time()
    row_df = these_samples_metadata %>%
      dplyr::filter(sample_id %in% rownames(cn_state_matrix)) %>%
      dplyr::select(sample_id, dplyr::all_of(metadataColumns)) %>%
      tibble::column_to_rownames("sample_id")
    row_df = row_df[rownames(cn_state_matrix), , drop = FALSE]
    meta_colours = map_metadata_to_colours(
      metadataColumns = metadataColumns,
      these_samples_metadata = these_samples_metadata,
      verbose = FALSE
    )
    meta_col_mat = matrix(
      unlist(lapply(metadataColumns, function(cn){
        unname(meta_colours[[cn]][as.character(row_df[[cn]])])
      })),
      nrow = nrow(row_df)
    )
    tick("build metadata colour strip", t0)
  }

  needs_gene_lookup = !missing(labelTheseGenes) || !missing(labelGainGenes) || !missing(labelLossGenes)
  resolved_bin_labels = list()
  resolved_gain_labels = list()
  resolved_loss_labels = list()
  if(needs_gene_lookup || !missing(bin_labels)){
    t0 = Sys.time()
    if(!missing(bin_labels)){
      resolved_bin_labels = as.list(bin_labels)
    }
    if(needs_gene_lookup){
      # vectorized bin parsing/matching, not a per-bin loop + region_to_chunks()
      # call per bin like pretty_CN_heatmap's map_bin_to_bin() -- fine at that
      # function's bin counts, but this one can have tens of thousands of bins
      col_start = as.integer(sub("^[^:]+:(\\d+)-.*", "\\1", colnames(cn_state_matrix)))
      col_end = as.integer(sub(".*-(\\d+)$", "\\1", colnames(cn_state_matrix)))
      norm_col_chrom = gsub("^chr", "", column_chromosome, ignore.case = TRUE)

      resolve_gene_bins = function(genes, existing, label_kind){
        result = existing
        for(g in genes){
          gene_region = tryCatch(
            suppressMessages(gene_to_region(g, projection = genome_build)),
            error = function(e) NA_character_
          )
          if(length(gene_region) == 0 || is.na(gene_region)){
            warning(sprintf("no region found for gene '%s' -- skipping its %s", g, label_kind))
            next
          }
          g_chrom = gsub("^chr", "", sub(":.*", "", gene_region), ignore.case = TRUE)
          g_start = as.integer(sub("^[^:]+:(\\d+)-.*", "\\1", gene_region))
          g_end = as.integer(sub(".*-(\\d+)$", "\\1", gene_region))
          # standard inclusive interval overlap -- any bin the gene touches
          # at all, including exactly on a boundary
          match_idx = which(
            norm_col_chrom == g_chrom & g_start <= col_end & g_end >= col_start
          )
          if(length(match_idx) == 0){
            warning(sprintf(
              "gene '%s' (%s) does not overlap any bin in cn_state_matrix -- skipping its %s",
              g, gene_region, label_kind
            ))
            next
          }
          # of the bins the gene overlaps, anchor to whichever one's own
          # centre is closest to the gene's centre (not just the first/
          # start-side match) -- matters for a gene that spans several bins
          g_mid = (g_start + g_end) / 2
          bin_mid = (col_start[match_idx] + col_end[match_idx]) / 2
          this_bin = colnames(cn_state_matrix)[match_idx[which.min(abs(bin_mid - g_mid))]]
          if(!is.null(result[[this_bin]])){
            message(sprintf(
              "Gene %s and gene %s share a bin. I will only show the last one!",
              g, result[[this_bin]]
            ))
          }
          result[[this_bin]] = g
        }
        result
      }

      if(!missing(labelTheseGenes)){
        resolved_bin_labels = resolve_gene_bins(labelTheseGenes, resolved_bin_labels, "label")
      }
      if(!missing(labelGainGenes)){
        resolved_gain_labels = resolve_gene_bins(labelGainGenes, list(), "gain label")
      }
      if(!missing(labelLossGenes)){
        resolved_loss_labels = resolve_gene_bins(labelLossGenes, list(), "loss label")
      }
    }
    tick("resolve gene/bin labels", t0)
  }
  show_gene_labels = length(resolved_bin_labels) > 0
  show_cumulative_gene_labels = length(resolved_gain_labels) > 0 || length(resolved_loss_labels) > 0

  if(return_data){
    tick("TOTAL", t_start)
    return(list(
      colour_matrix = col_mat,
      chromosome_columns = column_chromosome,
      data = cn_state_matrix,
      labels = resolved_bin_labels,
      gain_labels = resolved_gain_labels,
      loss_labels = resolved_loss_labels
    ))
  }

  show_cumulative = cumulative_plot != "none"
  if(show_cumulative){
    t0 = Sys.time()
    cum_df = pretty_CN_cumulative_plot(
      cn_state_matrix,
      gain_threshold = cumulative_gain_threshold,
      loss_threshold = cumulative_loss_threshold,
      as_fraction = cumulative_as_fraction,
      return_data = TRUE
    )
    tick("compute cumulative gain/loss", t0)

    if(cumulative_n_gain_peaks > 0 || cumulative_n_loss_peaks > 0){
      t0 = Sys.time()
      n_bins_total = nrow(cum_df)
      min_separation = max(1, round(n_bins_total * cumulative_peak_min_separation))

      # tallest-first greedy selection: take the highest remaining candidate,
      # then skip any candidate within min_separation of something already
      # picked -- a simple, robust stand-in for proper peak-finding that
      # naturally collapses a broad plateau/event into a single label
      find_top_peaks = function(values, n_peaks, min_height){
        if(n_peaks <= 0) return(integer(0))
        candidate_idx = which(values >= min_height)
        if(length(candidate_idx) == 0) return(integer(0))
        ord = candidate_idx[order(-values[candidate_idx])]
        selected = integer(0)
        for(i in ord){
          if(length(selected) >= n_peaks) break
          if(length(selected) == 0 || all(abs(i - selected) >= min_separation)){
            selected = c(selected, i)
          }
        }
        selected
      }

      gene_coords = if(genome_build == "grch37"){
        GAMBLR.data::grch37_gene_coordinates
      }else{
        GAMBLR.data::hg38_gene_coordinates
      }
      # restrict the reverse (position -> gene) lookup to genes with an
      # approved HGNC symbol -- there's no biotype column in this table to
      # filter on directly, but in practice a clone/contig-based placeholder
      # name (e.g. "RP11-299P2.1") essentially never has one, while a
      # properly named gene essentially always does, so this cleanly drops
      # placeholders without a coalesce() fallback to gene_name (which is
      # what was letting them through before). It does NOT catch a named
      # but non-protein-coding gene that has its own HGNC symbol (e.g. the
      # pseudogene "RPL7L1P1", the snoRNA "SNORD38A") -- there's no safe,
      # general pattern for that distinction here: e.g. excluding anything
      # ending "P<digit>" would also wrongly drop real genes like the AKAP
      # family (AKAP1, AKAP9, ...).
      gene_coords = gene_coords[!is.na(gene_coords$hugo_symbol), ]
      gc_chrom = gsub("^chr", "", gene_coords$chromosome, ignore.case = TRUE)
      gc_symbol = gene_coords$hugo_symbol

      genes_at_bin = function(bin_name){
        chrom = sub(":.*", "", bin_name)
        rest = sub("^[^:]+:", "", bin_name)
        b_start = as.integer(sub("-.*", "", rest))
        b_end = as.integer(sub(".*-", "", rest))
        norm_chrom = gsub("^chr", "", chrom, ignore.case = TRUE)
        hit = which(gc_chrom == norm_chrom & gene_coords$start <= b_end & gene_coords$end >= b_start)
        syms = unique(stats::na.omit(gc_symbol[hit]))
        if(length(syms) == 0){
          return(sprintf("%s:%.1fMb", chrom, round((b_start + b_end) / 2) / 1e6))
        }
        if(length(syms) > cumulative_peak_max_genes){
          syms = c(syms[seq_len(cumulative_peak_max_genes)],
                   sprintf("+%d", length(syms) - cumulative_peak_max_genes))
        }
        paste(syms, collapse = "/")
      }

      if(cumulative_n_gain_peaks > 0){
        peak_idx = find_top_peaks(
          cum_df$cumulative_gain, cumulative_n_gain_peaks,
          cumulative_peak_min_height * max(cum_df$cumulative_gain, 0, na.rm = TRUE)
        )
        for(i in peak_idx){
          bin_name = cum_df$bin[i]
          if(is.null(resolved_gain_labels[[bin_name]])){
            resolved_gain_labels[[bin_name]] = genes_at_bin(bin_name)
          }
        }
      }
      if(cumulative_n_loss_peaks > 0){
        peak_idx = find_top_peaks(
          cum_df$cumulative_loss, cumulative_n_loss_peaks,
          cumulative_peak_min_height * max(cum_df$cumulative_loss, 0, na.rm = TRUE)
        )
        for(i in peak_idx){
          bin_name = cum_df$bin[i]
          if(is.null(resolved_loss_labels[[bin_name]])){
            resolved_loss_labels[[bin_name]] = genes_at_bin(bin_name)
          }
        }
      }
      show_cumulative_gene_labels = length(resolved_gain_labels) > 0 || length(resolved_loss_labels) > 0
      tick("auto-detect peak genes", t0)
    }
  }

  t0 = Sys.time()
  if(newpage) grid::grid.newpage()

  meta_width_cm = if(has_meta) 0.4 * length(metadataColumns) else 0.001
  chrom_height_cm = 0.4
  legend_height_cm = if(show_legend) 1.5 else 0.001

  # row layout is built up as a list of (type, height) so the optional rows
  # (gene labels, cumulative plot) slot in without hardcoding row indices
  # for everything else. Gene labels always sit directly above the body
  # (their connector lines terminate right at its top edge); a top
  # cumulative panel goes above that, so it never gets crossed by a
  # connector line.
  row_types = "body"
  row_heights = height
  if(show_gene_labels){
    row_types = c("gene_labels", row_types)
    row_heights = c(gene_label_height_cm, row_heights)
  }
  if(show_cumulative && cumulative_plot == "top"){
    row_types = c("cumulative", row_types)
    row_heights = c(cumulative_plot_height, row_heights)
  }
  if(show_cumulative && cumulative_plot == "bottom"){
    body_pos = which(row_types == "body")
    row_types = append(row_types, "cumulative", after = body_pos)
    row_heights = append(row_heights, cumulative_plot_height, after = body_pos)
  }
  row_types = c(row_types, "chrom", "legend")
  row_heights = c(row_heights, chrom_height_cm, legend_height_cm)

  body_row = which(row_types == "body")
  cumulative_row = which(row_types == "cumulative")
  gene_labels_row = which(row_types == "gene_labels")
  chrom_row = which(row_types == "chrom")
  legend_row = which(row_types == "legend")

  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    nrow = length(row_types), ncol = 2,
    heights = grid::unit(row_heights, "cm"),
    widths = grid::unit(c(meta_width_cm, width), "cm")
  )))

  # heatmap body
  grid::pushViewport(grid::viewport(layout.pos.row = body_row, layout.pos.col = 2))
  grid::grid.raster(col_mat, width = 1, height = 1, interpolate = FALSE)
  grid::popViewport()

  # gene/bin labels -- own row directly above the body (see row_types
  # comment above), one short connector line per label plus a rotated text
  # label, both positioned on the same column-index x-axis as the body
  if(show_gene_labels){
    grid::pushViewport(grid::viewport(
      layout.pos.row = gene_labels_row, layout.pos.col = 2,
      xscale = c(0.5, n_col + 0.5), yscale = c(0, 1)
    ))
    for(bin_name in names(resolved_bin_labels)){
      col_idx = which(colnames(cn_state_matrix) == bin_name)
      if(length(col_idx) == 0) next
      label = resolved_bin_labels[[bin_name]]
      grid::grid.lines(
        x = grid::unit(c(col_idx, col_idx), "native"),
        y = grid::unit(c(0, bin_label_nudge), "npc"),
        gp = grid::gpar(lwd = 1)
      )
      grid::grid.text(
        label,
        x = grid::unit(col_idx, "native"), y = grid::unit(bin_label_nudge, "npc"),
        rot = bin_label_rotation, just = c("left", "bottom"),
        gp = grid::gpar(fontsize = bin_label_fontsize)
      )
    }
    grid::popViewport()
  }

  # cumulative gain/loss panel -- drawn as a single grid.polygon (one
  # graphics call regardless of bin count) on the same column-index x-axis
  # as everything else here, so it aligns exactly with the body
  if(show_cumulative){
    max_extent = max(cum_df$cumulative_gain, cum_df$cumulative_loss, 0, na.rm = TRUE)
    if(max_extent == 0) max_extent = 1
    n_bins = nrow(cum_df)
    y_extent = if(show_cumulative_gene_labels) max_extent * cumulative_label_headroom else max_extent

    grid::pushViewport(grid::viewport(
      layout.pos.row = cumulative_row, layout.pos.col = 2,
      xscale = c(0.5, n_bins + 0.5), yscale = c(-y_extent, y_extent)
    ))
    grid::grid.rect(gp = grid::gpar(col = NA, fill = "white"))
    grid::grid.polygon(
      x = grid::unit(c(cum_df$bin_index, rev(cum_df$bin_index)), "native"),
      y = grid::unit(c(cum_df$cumulative_gain, rep(0, n_bins)), "native"),
      gp = grid::gpar(fill = cumulative_gain_colour, col = cumulative_gain_colour)
    )
    grid::grid.polygon(
      x = grid::unit(c(cum_df$bin_index, rev(cum_df$bin_index)), "native"),
      y = grid::unit(c(-cum_df$cumulative_loss, rep(0, n_bins)), "native"),
      gp = grid::gpar(fill = cumulative_loss_colour, col = cumulative_loss_colour)
    )
    grid::grid.lines(
      x = grid::unit(c(0.5, n_bins + 0.5), "native"), y = grid::unit(c(0, 0), "native"),
      gp = grid::gpar(col = "black", lwd = 0.5)
    )

    # gain/loss gene labels -- connector runs from the gene's actual bin
    # value on the curve up/down to its label near the expanded headroom
    # edge, so it's visually obvious whether the gene sits on a local
    # peak/trough or somewhere unremarkable on the curve
    if(length(resolved_gain_labels) > 0){
      label_y = y_extent * 0.85
      for(bin_name in names(resolved_gain_labels)){
        row_idx = which(cum_df$bin == bin_name)
        if(length(row_idx) == 0) next
        grid::grid.lines(
          x = grid::unit(rep(cum_df$bin_index[row_idx], 2), "native"),
          y = grid::unit(c(cum_df$cumulative_gain[row_idx], label_y), "native"),
          gp = grid::gpar(lwd = 0.75, col = "black")
        )
        grid::grid.text(
          resolved_gain_labels[[bin_name]],
          x = grid::unit(cum_df$bin_index[row_idx], "native"), y = grid::unit(label_y, "native"),
          rot = bin_label_rotation, just = c("left", "bottom"),
          gp = grid::gpar(fontsize = bin_label_fontsize, col = cumulative_gain_colour)
        )
      }
    }
    if(length(resolved_loss_labels) > 0){
      label_y = -y_extent * 0.85
      for(bin_name in names(resolved_loss_labels)){
        row_idx = which(cum_df$bin == bin_name)
        if(length(row_idx) == 0) next
        grid::grid.lines(
          x = grid::unit(rep(cum_df$bin_index[row_idx], 2), "native"),
          y = grid::unit(c(-cum_df$cumulative_loss[row_idx], label_y), "native"),
          gp = grid::gpar(lwd = 0.75, col = "black")
        )
        grid::grid.text(
          resolved_loss_labels[[bin_name]],
          x = grid::unit(cum_df$bin_index[row_idx], "native"), y = grid::unit(label_y, "native"),
          rot = bin_label_rotation, just = c("left", "top"),
          gp = grid::gpar(fontsize = bin_label_fontsize, col = cumulative_loss_colour)
        )
      }
    }
    grid::popViewport()

    # max-value labels only -- a full numeric axis would need its own
    # reserved gutter width, which would throw off column alignment with
    # the body; see pretty_CN_cumulative_plot() for a fully-labelled version
    grid::pushViewport(grid::viewport(layout.pos.row = cumulative_row, layout.pos.col = 1))
    label_mult = if(cumulative_as_fraction) 100 else 1
    label_suffix = if(cumulative_as_fraction) "%" else ""
    grid::grid.text(
      sprintf("%.0f%s", max(cum_df$cumulative_gain) * label_mult, label_suffix),
      x = grid::unit(1, "npc"), y = grid::unit(0.95, "npc"),
      just = c("right", "top"), gp = grid::gpar(fontsize = 6, col = cumulative_gain_colour)
    )
    grid::grid.text(
      sprintf("%.0f%s", max(cum_df$cumulative_loss) * label_mult, label_suffix),
      x = grid::unit(1, "npc"), y = grid::unit(0.05, "npc"),
      just = c("right", "bottom"), gp = grid::gpar(fontsize = 6, col = cumulative_loss_colour)
    )
    grid::popViewport()
  }

  # chromosome boundary strip
  grid::pushViewport(grid::viewport(
    layout.pos.row = chrom_row, layout.pos.col = 2, xscale = c(0, n_col)
  ))
  grid::grid.rect(
    x = grid::unit(chrom_start - 1, "native"),
    width = grid::unit(chrom_rle$lengths, "native"),
    just = "left",
    gp = grid::gpar(fill = chrom_col[chrom_rle$values], col = NA)
  )
  grid::grid.text(
    chrom_rle$values,
    x = grid::unit(chrom_mid, "native"), y = grid::unit(0.5, "npc"),
    gp = grid::gpar(fontsize = 6)
  )
  grid::popViewport()

  # metadata annotation strip
  if(has_meta){
    grid::pushViewport(grid::viewport(layout.pos.row = body_row, layout.pos.col = 1))
    grid::grid.raster(meta_col_mat, width = 1, height = 1, interpolate = FALSE)
    grid::popViewport()
  }

  # legend
  if(show_legend){
    grid::pushViewport(grid::viewport(layout.pos.row = legend_row, layout.pos.col = 1:2))
    lgd = ComplexHeatmap::Legend(
      col_fun = cn_col_fun, title = legend_title, direction = "horizontal"
    )
    ComplexHeatmap::draw(lgd, x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"))
    grid::popViewport()
  }

  grid::popViewport()
  tick("draw", t0)
  tick("TOTAL", t_start)
}
