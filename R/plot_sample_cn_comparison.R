#' @title Compare copy number segmentations of one sample
#'
#' @description Draw one sample's copy number segments from several sources
#' (e.g. pipelines) as aligned tracks, genome-wide or in a region, to see
#' where the segmentations agree and where they don't.
#'
#' @details By default (`value = "auto"`) each track is drawn on the scale
#' its copy number features are called on ([GAMBLR.utils::cnv_calling_settings];
#' a `cn_type` column, e.g. from [GAMBLR.utils::calibrate_relative_cn],
#' overrides a pipeline's type):
#' pipelines with purity- and ploidy-corrected copy number (e.g. battenberg,
#' PureCN, sequenza; see [GAMBLR.utils::cn_pipeline_types]) as CN, with a line
#' at the sample's dominant copy number; the rest (e.g. CNVkit and gatk_cnv,
#' whose CN is only a rescaled log ratio) as log ratio, centred on the track's
#' baseline (its length-weighted median log ratio over autosomes 1-22, whatever
#' `chromosomes` or `region` shows). Each track has its own y axis. Dashed and
#' dotted lines mark where a gain or loss reaches level 1 and level 2 (the
#' feature set's thresholds, or `cnv_args`), so what is drawn is what is
#' thresholded. `value = "log_ratio"` or `"CN"` draws every track on one scale;
#' threshold lines are then drawn only on tracks whose features are called on
#' that scale. Segments beyond `limit` are drawn at it.
#'
#' With `cnv_features`, the features (e.g. a classifier's arm and focal gains
#' and losses) are called separately from each source's genome-wide segments
#' with [GAMBLR.utils::segmented_data_to_cnv_status], using each source's
#' settings from [GAMBLR.utils::cnv_calling_settings] (and `cnv_args`). Features called by at least one source and within the view
#' are shaded across all tracks and labelled in a "features" panel on top,
#' and each source that calls one gets a marker on its track: an upward
#' triangle for a gain, downward for a loss, filled for a level 2 (high
#' level) event. A feature is not called from single segments: each feature
#' is reduced to one copy number change (for an arm, the level a set fraction
#' of the arm, e.g. 75% for "LymphGen_model", reaches; for a gene, its
#' overlapping segments, or with LymphGen's gene rule the one segment that
#' decides it, see `gene_rule` in [GAMBLR.utils::segmented_data_to_cnv_status])
#' and that is compared with the thresholds. On each
#' track drawn on its calling scale, a short bar across each shaded feature
#' marks that value, so a feature is called when its bar passes the line:
#' black when called, grey when not, with an orange halo when the value is
#' within 0.1 copies of a threshold on either side (a borderline call). Gene features limited to short events
#' (e.g. LymphGen's 30 Mb) aren't called from longer events even when the bar
#' passes. The calls are returned as the plot's `feature_calls` attribute
#' (`attr(p, "feature_calls")`), with each feature's value, change (signed,
#' positive in the feature's direction), margin past the nearer threshold,
#' borderline and status ("called", "below threshold", "event too long",
#' "outweighed by opposite change" (LymphGen's gene rule: a larger or
#' higher-level change in the other direction decides the gene, and the bar
#' shows that change) or "not assessed"). Feature calling supports grch37 only.
#'
#' @param seg_data Segments with ID, chrom, start, end, log.ratio (and CN for
#' `value = "CN"`) and the column named in `by`, e.g. from
#' `GAMBLR.results::get_cn_segments(..., pipeline = "all")`. It may hold
#' other samples too; only the one in `these_samples_metadata` is drawn.
#' @param these_samples_metadata One row of metadata naming the sample to
#' show (its sample_id).
#' @param by Column with the source of each segment, one track per value.
#' Default "pipeline".
#' @param value "auto" (default; each track on the scale its features are
#' called on), "log_ratio" or "CN". See Details.
#' @param center Centre log ratio tracks on their baseline (see Details).
#' Default TRUE.
#' @param region Optional region ("chr:start-end") to zoom in on; default the
#' genome (`chromosomes`).
#' @param chromosomes Chromosomes shown genome-wide. Default 1-22. Each
#' track's baseline always comes from autosomes 1-22.
#' @param projection Genome build of the coordinates ("grch37" or "hg38").
#' Default: from `seg_data` (as set by `get_cn_segments`), else "grch37".
#' @param limit Values beyond this distance from the baseline are drawn at
#' it: one number for every track, or a vector named log_ratio and CN.
#' Default `c(log_ratio = 2, CN = 4)`.
#' @param highlight Optional source(s) to mark in the track labels, e.g. the
#' pipeline a sample's combined-merge segments come from.
#' @param selection Optional table of segmentation scores from
#' [GAMBLR.utils::compare_segmentations], to show the choice of the best
#' segmentation. Score the whole cohort or project with it, not just this
#' sample: its over-segmentation penalty compares each segmentation with the
#' median over every sample of the same seq_type and genome version. The
#' selected source's track is shaded green and tagged "best", each scored
#' track's label shows its QC score (lower is better) and rank, and the
#' subtitle names the choice. Default NULL: no selection.
#' @param order Optional order of the tracks (top to bottom); default
#' alphabetical.
#' @param cnv_features Optional features to overlay: a feature set name for
#' [GAMBLR.utils::get_cnv_features] (e.g. "DLBClass", or "LymphGen_model" for
#' the features LymphGen scores, with its rules, as GAMBLR.predict's
#' lymphgen_input() calls them; "LymphGen" is a broader set of
#' LymphGen-associated features with looser rules), or a feature table as it
#' returns. Default NULL: no overlay. See Details.
#' @param cnv_args Optional named list of further arguments for
#' [GAMBLR.utils::segmented_data_to_cnv_status] (e.g. `blacklist`); its
#' `value_scale`, `centering`, `single_thresh` and `high_thresh` override
#' the per-pipeline settings for every track.
#' @param thresholds Draw the level 1 and level 2 threshold lines. Default
#' TRUE.
#' @param base_size Base font size. Default 9.
#'
#' @return A ggplot object; each track's scale, baseline and thresholds are
#' in `attr(p, "track_settings")`; with `selection`, the scores are in
#' `attr(p, "selection")`, and with `cnv_features` the feature calls are in
#' `attr(p, "feature_calls")`.
#'
#' @import ggplot2 dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' meta = get_gambl_metadata() %>% dplyr::filter(sample_id == "DLBCL-RICOVER_299-Tumor")
#' segs = GAMBLR.results::get_cn_segments(meta, pipeline = "all")
#' plot_sample_cn_comparison(segs, meta)
#' plot_sample_cn_comparison(segs, meta, region = "chr9:19000000-24000000")
#'
#' # mark the pipeline the combined merge uses for this sample
#' av = GAMBLR.results::cn_pipeline_availability(meta)
#' plot_sample_cn_comparison(segs, meta, highlight = av$combined_module)
#'
#' # which sources call the DLBClass copy number features
#' p = plot_sample_cn_comparison(segs, meta, cnv_features = "DLBClass")
#' attr(p, "feature_calls")
#'
#' # the best segmentation: score the whole cohort once, then show any sample
#' cohort_meta = get_gambl_metadata() %>% dplyr::filter(cohort == "dlbcl_chapuy")
#' cohort_segs = GAMBLR.results::get_cn_segments(cohort_meta, pipeline = "all")
#' scores = GAMBLR.utils::compare_segmentations(cohort_segs, cohort_meta)
#' plot_sample_cn_comparison(cohort_segs, dplyr::filter(cohort_meta, sample_id == "DLBCL10541T"),
#'                           selection = scores)
#' }
plot_sample_cn_comparison = function(seg_data,
                                     these_samples_metadata,
                                     by = "pipeline",
                                     value = "auto",
                                     center = TRUE,
                                     region = NULL,
                                     chromosomes = as.character(1:22),
                                     projection = NULL,
                                     limit = NULL,
                                     highlight = NULL,
                                     order = NULL,
                                     cnv_features = NULL,
                                     cnv_args = list(),
                                     thresholds = TRUE,
                                     selection = NULL,
                                     base_size = 9){
  value = match.arg(value, c("auto", "log_ratio", "CN"))
  if(is.null(projection)){
    projection = tryCatch(GAMBLR.utils::get_genome_build(seg_data), error = function(e) NULL)
    if(is.null(projection) || !projection %in% c("grch37", "hg38")) projection = "grch37"
  }
  segs = as.data.frame(seg_data)
  if(!by %in% colnames(segs)) stop("seg_data has no '", by, "' column")
  if(missing(these_samples_metadata) || nrow(these_samples_metadata) != 1){
    stop("these_samples_metadata must be one row, naming the sample to show")
  }
  sample_id = as.character(these_samples_metadata$sample_id)
  segs = segs[segs$ID == sample_id, , drop = FALSE]
  if(nrow(segs) == 0) stop("no segments for ", sample_id)
  segs$chrom = gsub("^chr", "", as.character(segs$chrom))
  segs$track = as.character(segs[[by]])
  limits = c(log_ratio = 2, CN = 4)
  if(!is.null(limit)){
    if(is.null(names(limit))) limits[] = limit else limits[names(limit)] = limit
  }

  # how each track's features are called: its pipeline (the source, or its
  # one module, e.g. for a combined merge) decides the scale and centring;
  # the feature set the thresholds
  feature_set = if(is.character(cnv_features)) cnv_features[1] else
    if(is.data.frame(cnv_features) && "feature_set" %in% colnames(cnv_features)) unique(cnv_features$feature_set)[1] else NULL
  settings = .track_settings(segs, feature_set, cnv_args)
  # the scale each track is drawn on
  settings$display = if(value == "auto") settings$value_scale else value
  missing_cn = settings$display == "CN" & !settings$has_cn
  if(any(missing_cn)){
    if(value == "CN") message("no CN for ", paste(settings$track[missing_cn], collapse = ", "), "; drawn as log ratio")
    settings$display[missing_cn] = "log_ratio"
  }
  # each track's baseline on its drawn scale, from autosomes 1-22 genome-wide
  # (never just the chromosomes or region shown), as the feature calls use it
  settings$baseline = vapply(seq_len(nrow(settings)), function(i){
    d = segs[segs$track == settings$track[i], , drop = FALSE]
    scale = settings$display[i]
    centering = if(scale == settings$value_scale[i]) settings$centering[i] else if(scale == "CN") "mode" else "median"
    b = tryCatch(GAMBLR.utils::cnv_sample_baseline(d, scale, centering), error = function(e) numeric(0))
    if(length(b) == 0) if(scale == "CN") 2 else 0 else unname(b[1])
  }, numeric(1))
  # drawn values: centred log ratio (0 = baseline) or CN (baseline as is)
  settings$shift = ifelse(settings$display == "log_ratio" & center, settings$baseline, 0)
  settings$mid = settings$baseline - settings$shift
  settings$limit = limits[settings$display]
  settings$lo = ifelse(settings$display == "CN", pmax(0, settings$mid - settings$limit), settings$mid - settings$limit)
  settings$hi = settings$mid + settings$limit
  at = match(segs$track, settings$track)
  segs$y = ifelse(settings$display[at] == "CN", segs$CN, segs$log.ratio)
  no_value = is.na(segs$y)
  if(any(no_value)){
    message(sum(no_value), " segment(s) without a value on their track's scale left out (",
            paste(unique(segs$track[no_value]), collapse = ", "), ")")
    segs = segs[!no_value, , drop = FALSE]
    at = match(segs$track, settings$track)
  }
  segs$y = segs$y - settings$shift[at]
  segs$shown = pmin(pmax(segs$y, settings$lo[at]), settings$hi[at])
  # colour: distance from the baseline as a fraction of the track's limit
  segs$tone = pmax(-1, pmin(1, (segs$y - settings$mid[at]) / settings$limit[at]))

  calls = NULL
  if(!is.null(cnv_features)){
    calls = .call_features_by_track(seg_data, sample_id, by, cnv_features, settings, cnv_args, projection)
  }

  arms = if(projection == "hg38") GAMBLR.data::chromosome_arms_hg38 else GAMBLR.data::chromosome_arms_grch37
  arms$chromosome = gsub("^chr", "", as.character(arms$chromosome))
  if(is.null(region)){
    chrom_len = vapply(chromosomes, function(ch) max(arms$end[arms$chromosome == ch]), numeric(1))
    offset = c(0, cumsum(as.numeric(chrom_len)))[seq_along(chromosomes)]
    names(offset) = chromosomes
    segs = segs[segs$chrom %in% chromosomes, , drop = FALSE]
    segs$x0 = offset[segs$chrom] + segs$start
    segs$x1 = offset[segs$chrom] + segs$end
    x_breaks = offset + chrom_len / 2
    x_labels = chromosomes
    bounds = offset[-1]
    x_title = NULL
    if(!is.null(calls)){
      calls$chrom = gsub("^chr", "", as.character(calls$chrom))
      calls = calls[calls$chrom %in% chromosomes, , drop = FALSE]
      calls$x0 = offset[calls$chrom] + calls$start
      calls$x1 = offset[calls$chrom] + calls$end
    }
  } else {
    clean = gsub(",", "", region)
    parts = regmatches(clean, regexec("^(chr)?([0-9XY]+):([0-9.e+]+)-([0-9.e+]+)$", clean))[[1]]
    if(length(parts) == 0) stop("region must look like \"chr9:19000000-24000000\"")
    r_chrom = parts[3]; r_start = as.numeric(parts[4]); r_end = as.numeric(parts[5])
    segs = segs[segs$chrom == r_chrom & segs$end > r_start & segs$start < r_end, , drop = FALSE]
    segs$x0 = pmax(segs$start, r_start) / 1e6
    segs$x1 = pmin(segs$end, r_end) / 1e6
    x_breaks = ggplot2::waiver(); x_labels = ggplot2::waiver(); bounds = NULL
    x_title = paste0("chr", r_chrom, " (Mb)")
    if(!is.null(calls)){
      calls$chrom = gsub("^chr", "", as.character(calls$chrom))
      calls = calls[calls$chrom == r_chrom & calls$end > r_start & calls$start < r_end, , drop = FALSE]
      calls$x0 = pmax(calls$start, r_start) / 1e6
      calls$x1 = pmin(calls$end, r_end) / 1e6
    }
  }

  # track labels: source, number of segments (in the view), highlight mark
  tracks = if(is.null(order)) sort(unique(as.character(seg_data[[by]][seg_data$ID == sample_id]))) else order
  counts = table(factor(segs$track, levels = tracks))
  # selection: this sample's QC scores per source
  scores = NULL
  if(!is.null(selection)){
    if(!is.data.frame(selection)) stop("selection must be a table from GAMBLR.utils::compare_segmentations")
    scores = as.data.frame(selection)
    scores = scores[scores$ID == sample_id, , drop = FALSE]
    if(nrow(scores) == 0) message("selection has no rows for ", sample_id)
  }
  best = if(!is.null(scores) && nrow(scores) > 0) as.character(scores[[by]][scores$selected %in% TRUE]) else character(0)
  score_text = if(!is.null(scores) && nrow(scores) > 0){
    at = match(tracks, scores[[by]])
    ifelse(is.na(at), "", paste0("\nscore ", format(scores$score[at], nsmall = 2), " (#", scores$rank[at], ")",
                                 ifelse(tracks %in% best, " best", "")))
  } else rep("", length(tracks))
  scale_text = vapply(tracks, function(t){
    i = match(t, settings$track)
    if(is.na(i)) return("")
    if(settings$display[i] == "CN") sprintf("\nCN, baseline %.1f", settings$baseline[i]) else "\nlog ratio"
  }, "")
  labels = paste0(tracks, ifelse(tracks %in% highlight, " *", ""), "\n", as.integer(counts), " segs",
                  scale_text, score_text)
  show_features = !is.null(calls) && any(calls$level > 0)
  panel_levels = if(show_features) c("features", labels) else labels
  segs$track = factor(segs$track, levels = tracks, labels = labels)
  segs$track = factor(as.character(segs$track), levels = panel_levels)

  shown_scales = unique(settings$display[settings$track %in% tracks])
  y_title = if(identical(shown_scales, "CN")) "copy number" else
    if(identical(shown_scales, "log_ratio")) if(center) "log ratio (centred)" else "log ratio" else
      "copy number or log ratio (see each track)"
  track_labels = stats::setNames(labels, tracks)
  panel = function(track) factor(track_labels[track], levels = panel_levels)
  per_track = settings[settings$track %in% tracks, , drop = FALSE]
  # each panel's y range (free scales), and the features panel's
  ranges = data.frame(track = factor(c(rep(as.character(panel(per_track$track)), each = 2),
                                       if(show_features) rep("features", 2)), levels = panel_levels),
                      y = c(rbind(per_track$lo, per_track$hi), if(show_features) c(-1, 1)))
  p = ggplot2::ggplot(segs) +
    ggplot2::geom_blank(data = ranges, ggplot2::aes(y = y), inherit.aes = FALSE)
  # the automatically selected source's panel is shaded, under everything else
  if(length(best) > 0){
    p = p + ggplot2::geom_rect(data = data.frame(track = factor(labels[tracks %in% best], levels = panel_levels)),
                               xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                               fill = "#E3F2E1", inherit.aes = FALSE)
  }
  p = p +
    ggplot2::geom_hline(data = data.frame(track = panel(per_track$track), mid = per_track$mid),
                        ggplot2::aes(yintercept = mid), colour = "grey70", linewidth = 0.3)
  if(thresholds){
    lines = .threshold_lines(per_track)
    if(!is.null(lines) && nrow(lines) > 0){
      lines$track = panel(lines$track)
      p = p + ggplot2::geom_hline(data = lines, ggplot2::aes(yintercept = y),
                                  colour = ifelse(lines$direction == "gain", "#C0392B", "#2E6FB7"),
                                  linetype = ifelse(lines$level == 1, "dashed", "dotted"),
                                  linewidth = 0.3, alpha = 0.8)
    }
  }
  if(!is.null(bounds)){
    p = p + ggplot2::geom_vline(xintercept = bounds, colour = "grey85", linewidth = 0.3)
  }
  bars = NULL
  if(show_features){
    called = calls[calls$level > 0, , drop = FALSE]
    # every feature any source calls: a band in all panels and a label on top
    spans = unique(called[, c("name", "label", "direction", "x0", "x1")])
    p = p +
      ggplot2::geom_rect(data = spans, ggplot2::aes(xmin = x0, xmax = x1, ymin = -Inf, ymax = Inf,
                                                    fill = direction),
                         inherit.aes = FALSE, alpha = 0.15) +
      ggplot2::scale_fill_manual(values = c(gain = "#C0392B", loss = "#2E6FB7"), guide = "none") +
      ggrepel::geom_text_repel(data = transform(spans, track = factor("features", levels = panel_levels),
                                                x = (x0 + x1) / 2, y = 0),
                               ggplot2::aes(x = x, y = y, label = label),
                               inherit.aes = FALSE, size = base_size / 3.5, direction = "y",
                               segment.size = 0.2, min.segment.length = 0, max.overlaps = Inf)
    # a marker on each track whose source calls the feature
    called$track = panel(called$source)
    called$x = (called$x0 + called$x1) / 2
    # one row at the top of each track (deep losses, drawn at the bottom, stay visible);
    # the triangle's direction gives the event's
    st = settings[match(called$source, settings$track), , drop = FALSE]
    called$y = st$lo + 0.94 * (st$hi - st$lo)
    p = p + ggplot2::geom_point(data = called, ggplot2::aes(x = x, y = y),
                                shape = ifelse(called$direction == "gain", 24, 25),
                                fill = ifelse(called$level >= 2,
                                              ifelse(called$direction == "gain", "#C0392B", "#2E6FB7"), "white"),
                                colour = ifelse(called$direction == "gain", "#C0392B", "#2E6FB7"),
                                size = 2, stroke = 0.6, inherit.aes = FALSE)
    # the value each shaded feature is judged on, on each track drawn on its
    # calling scale: a feature is called when its bar passes the line
    bars = calls[calls$name %in% spans$name & !is.na(calls$value), , drop = FALSE]
    st = settings[match(bars$source, settings$track), , drop = FALSE]
    bars = bars[!is.na(st$display) & st$display == st$value_scale, , drop = FALSE]
    st = settings[match(bars$source, settings$track), , drop = FALSE]
    if(nrow(bars) > 0){
      bars$y = ifelse(st$display == "CN", st$mid + bars$value,
                      st$mid + log2(pmax(2 + bars$value, 1e-3) / 2))
      bars$y = pmin(pmax(bars$y, st$lo), st$hi)
      bars$track = panel(bars$source)
      bars$bar_colour = ifelse(bars$level > 0, "black", "grey55")
    }else bars = NULL
  }
  p = p +
    ggplot2::geom_segment(ggplot2::aes(x = x0, xend = x1, y = shown, yend = shown, colour = tone),
                          linewidth = 1.6, lineend = "butt") +
    # borderline values get an orange halo behind their bar; the bar's
    # colour always gives the call
    (if(show_features && !is.null(bars) && any(bars$borderline)) ggplot2::geom_segment(
      data = bars[bars$borderline, , drop = FALSE], ggplot2::aes(x = x0, xend = x1, y = y, yend = y),
      colour = "#F5A04A", linewidth = 2.6, inherit.aes = FALSE)) +
    (if(show_features && !is.null(bars)) ggplot2::geom_segment(
      data = bars, ggplot2::aes(x = x0, xend = x1, y = y, yend = y), colour = bars$bar_colour,
      linewidth = 0.7, inherit.aes = FALSE)) +
    ggplot2::scale_colour_gradient2(low = "#2E6FB7", mid = "grey55", high = "#C0392B", midpoint = 0,
                                    limits = c(-1, 1), guide = "none") +
    ggplot2::facet_grid(track ~ ., drop = FALSE, scales = "free_y") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.04)) +
    ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.005, 0)) +
    ggplot2::labs(x = x_title, y = y_title,
                  caption = .null_if_empty(paste(c(
                    if(thresholds) "Lines: dashed = level 1, dotted = level 2 threshold, from each track's baseline (red gain, blue loss).",
                    if(show_features) paste0("Bars: the value each shaded feature is judged on (e.g. the level a set fraction of an arm ",
                                             "reaches); called when it passes the line.\n",
                                             "Black = called, grey = not called; an orange halo = within 0.1 copies of a ",
                                             "threshold (borderline, either side).")),
                    collapse = "\n")),
                  title = sample_id,
                  subtitle = .null_if_empty(paste(c(
                    if(!is.null(highlight) && any(tracks %in% highlight)) "* used in the combined merge",
                    if(length(best) > 0) paste0("best by QC score (shaded): ", paste(best, collapse = ", "))),
                    collapse = "    "))) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_text(angle = 0, hjust = 0),
                   strip.background = ggplot2::element_rect(fill = "grey95"),
                   plot.caption = ggplot2::element_text(hjust = 0, colour = "grey30"))
  attr(p, "track_settings") = settings[, c("track", "pipeline", "cn_type", "value_scale", "centering",
                                           "single_thresh", "high_thresh", "display", "baseline")]
  if(!is.null(calls)) attr(p, "feature_calls") = calls[, setdiff(colnames(calls), c("x0", "x1")), drop = FALSE]
  if(!is.null(scores)) attr(p, "selection") = scores
  p
}

.null_if_empty = function(x) if(length(x) == 0 || identical(x, "")) NULL else x

# Calling settings per track (GAMBLR.utils::cnv_calling_settings): the
# pipeline is the track's name, or its one module (e.g. a combined merge's
# segments); a CN-scale track without CN falls back to log ratios
.track_settings = function(segs, feature_set, cnv_args){
  tracks = sort(unique(segs$track))
  known = GAMBLR.utils::cn_pipeline_types()$pipeline
  pipeline = vapply(tracks, function(t){
    if(t %in% known || !"module" %in% colnames(segs)) return(t)
    modules = unique(stats::na.omit(as.character(segs$module[segs$track == t])))
    if(length(modules) == 1) modules else t
  }, "")
  has_cn = vapply(tracks, function(t) "CN" %in% colnames(segs) && any(!is.na(segs$CN[segs$track == t])), TRUE)
  # a cn_type column (e.g. from GAMBLR.utils::calibrate_relative_cn) overrides
  # the pipeline's type
  cn_type = vapply(tracks, function(t){
    if(!"cn_type" %in% colnames(segs)) return(NA_character_)
    types = unique(stats::na.omit(as.character(segs$cn_type[segs$track == t])))
    if(length(types) == 1) types else NA_character_
  }, "")
  known_type = GAMBLR.utils::cn_pipeline_types()$cn_type[match(pipeline, GAMBLR.utils::cn_pipeline_types()$pipeline)]
  cn_type = ifelse(is.na(cn_type), known_type, cn_type)
  settings = GAMBLR.utils::cnv_calling_settings(pipeline, feature_set, cnv_args, cn_type = cn_type)
  if(any(settings$value_scale == "CN" & !has_cn)){
    fallback = GAMBLR.utils::cnv_calling_settings("relative", feature_set,
                                                  utils::modifyList(cnv_args, list(value_scale = "log_ratio")))
    swap = which(settings$value_scale == "CN" & !has_cn)
    settings[swap, c("value_scale", "centering", "single_thresh", "high_thresh")] =
      fallback[rep(1, length(swap)), c("value_scale", "centering", "single_thresh", "high_thresh")]
  }
  data.frame(track = tracks, settings, has_cn = has_cn, row.names = NULL, stringsAsFactors = FALSE)
}

# Threshold lines for tracks drawn on the scale their features are called on:
# level 1 and 2 gains and losses as copy changes from the baseline, which on
# the log ratio scale are log2((2 + change) / 2) from the baseline
.threshold_lines = function(settings){
  st = settings[settings$display == settings$value_scale, , drop = FALSE]
  if(nrow(st) == 0) return(NULL)
  out = do.call(rbind, lapply(seq_len(nrow(st)), function(i){
    change = c(st$single_thresh[i], st$high_thresh[i])
    gain = if(st$display[i] == "CN") st$mid[i] + change else st$mid[i] + log2((2 + change) / 2)
    loss = if(st$display[i] == "CN") st$mid[i] - change else
      ifelse(change < 2, st$mid[i] + log2(pmax(2 - change, 1e-9) / 2), NA)
    data.frame(track = st$track[i], level = c(1, 2, 1, 2), direction = rep(c("gain", "loss"), each = 2),
               y = c(gain, loss), lo = st$lo[i], hi = st$hi[i], stringsAsFactors = FALSE)
  }))
  out[!is.na(out$y) & out$y >= out$lo & out$y <= out$hi, , drop = FALSE]
}

# Feature calls (long format, one row per source and feature) from each
# source's genome-wide segments for one sample, with the features' direction
# and labels; each source called with its own settings
.call_features_by_track = function(seg_data, sample_id, by, cnv_features, settings, cnv_args, projection){
  if(!projection %in% c("grch37", "hg19")){
    stop("cnv_features: feature calling supports grch37 only")
  }
  features = if(is.character(cnv_features)) {
    GAMBLR.utils::get_cnv_features(projection, feature_set = cnv_features)
  } else cnv_features
  segs = as.data.frame(seg_data)
  segs = segs[segs$ID == sample_id, , drop = FALSE]
  meta = data.frame(sample_id = sample_id, stringsAsFactors = FALSE)
  out = lapply(sort(unique(as.character(segs[[by]]))), function(src){
    s = segs[segs[[by]] == src, , drop = FALSE]
    st = settings[match(src, settings$track), , drop = FALSE]
    own = if(nrow(st) == 1 && !is.na(st$value_scale)) list(value_scale = st$value_scale, centering = st$centering,
                                                           single_thresh = st$single_thresh,
                                                           high_thresh = st$high_thresh,
                                                           gene_rule = st$gene_rule,
                                                           homdel_max_cn = st$homdel_max_cn) else list()
    args = utils::modifyList(utils::modifyList(cnv_args, own), list(seg_data = s, cnv_features = features,
                                            these_samples_metadata = meta, genome_build = projection,
                                            output_format = "long", encoding = "level",
                                            gene_labels = TRUE))
    called = tryCatch(suppressMessages(do.call(GAMBLR.utils::segmented_data_to_cnv_status, args)),
                      error = function(e){
                        message("no feature calls for ", src, ": ", conditionMessage(e))
                        NULL
                      })
    if(is.null(called) || nrow(called) == 0) return(NULL)
    called = as.data.frame(called)
    # the copy number change each call is judged on (e.g. the value 80% of
    # an arm is at or beyond)
    args$encoding = "value"
    value = tryCatch(suppressMessages(do.call(GAMBLR.utils::segmented_data_to_cnv_status, args)),
                     error = function(e) NULL)
    called$value = if(is.null(value)) NA_real_ else as.data.frame(value)$value[match(called$name, value$name)]
    # LymphGen's gene rule says why each of its features is or isn't called
    if(!"rule_status" %in% colnames(called)) called$rule_status = NA_character_
    called$single_thresh = if(nrow(st) == 1) st$single_thresh else NA_real_
    called$high_thresh = if(nrow(st) == 1) st$high_thresh else NA_real_
    called$source = src
    called
  })
  calls = do.call(rbind, out)
  if(is.null(calls)) return(NULL)
  if(!"label" %in% colnames(calls)) calls$label = calls$name
  calls$direction = as.character(features$direction[match(calls$name, features$name)])
  calls$level[is.na(calls$level)] = 0
  # signed change (positive = in the feature's direction), its margin past
  # the nearer threshold, and why a call was or wasn't made
  calls$change = ifelse(calls$direction == "gain", calls$value, -calls$value)
  calls$margin = ifelse(calls$change > calls$high_thresh, calls$change - calls$high_thresh,
                        calls$change - calls$single_thresh)
  near = pmin(abs(calls$change - calls$single_thresh), abs(calls$change - calls$high_thresh))
  calls$borderline = !is.na(near) & near < 0.1
  calls$status = ifelse(calls$level > 0, "called",
                        ifelse(is.na(calls$value), "not assessed",
                               ifelse(calls$change > calls$single_thresh, "event too long", "below threshold")))
  calls$status = ifelse(!is.na(calls$rule_status) & !is.na(calls$value), calls$rule_status, calls$status)
  calls$rule_status = NULL
  calls
}
