#' @title Plot pairwise co-occurrence and mutual exclusivity between events
#'
#' @description Visualise the output of
#' [GAMBLR.utils::pairwise_event_association] as an event-by-event matrix, a
#' volcano plot, or a comparison of two analyses of the same events (e.g.
#' unadjusted Fisher tests against burden-adjusted logistic regression).
#'
#' @details
#' \strong{"matrix"}: one cell per pair of events. Pairs with q below
#' `q_threshold` get a circle coloured by log odds ratio (red for
#' co-occurring, blue for mutually exclusive) and sized by -log10 q; other
#' tested pairs get a faint tile in the same colour scale (hide them with
#' `show_nonsignificant = FALSE`). Pairs of events on the same chromosome
#' that weren't tested (e.g. same-arm pairs, removed by default) are
#' grey. When the table has event coordinates (it does when
#' `pairwise_event_association` was given `feature_info`), events are
#' ordered along the genome, grouped by chromosome, and marked with labelled
#' chromosome strips; events without coordinates (e.g. mutated genes) follow
#' in their own group.
#'
#' \strong{Clustering} (`cluster = TRUE`): events are grouped by their
#' association profiles, i.e. events with the same partners and the same
#' exclusions end up together, which can reveal modules of events spread
#' across chromosomes. Each event's profile is its row of pairwise scores
#' (`cluster_on`), and the distance is 1 - the correlation of two profiles,
#' with average linkage. "signed_p" is the default because rare events give
#' extreme odds ratios from a handful of samples; weighting by significance
#' keeps those from driving the tree. Pairs without a value (e.g. same
#' chromosome pairs that weren't tested) count as 0, so events aren't
#' grouped just for sharing a chromosome. The chromosome strip then labels
#' each event's chromosome. Clustering unadjusted results tends to recover
#' a genomic instability cluster; burden-adjusted results show the
#' remaining structure. This summarises pairwise statistics: groups it
#' suggests are hypotheses to test at the sample level.
#'
#' \strong{"volcano"}: log odds ratio against -log10 q for every tested pair,
#' with the most significant pairs labelled.
#'
#' \strong{"compare"}: the log odds ratio of each pair in `compare_to` (x)
#' against `assoc` (y). Pairs near the diagonal agree; pairs pulled towards
#' zero in `assoc` lost their association, e.g. because it was explained by
#' genomic instability when `assoc` is adjusted for burden.
#'
#' In all styles, log odds ratios beyond `max_log_odds_ratio` (including the
#' infinite values of pairs that never, or always, co-occur) are drawn at
#' that limit.
#'
#' @param assoc Output of [GAMBLR.utils::pairwise_event_association].
#' @param style "matrix" (default), "volcano" or "compare".
#' @param compare_to For `style = "compare"`, a second
#' [GAMBLR.utils::pairwise_event_association] table for the same events,
#' plotted on the x axis.
#' @param q_threshold q value below which a pair is significant. Default 0.05.
#' @param max_log_odds_ratio Limit of the log odds ratio scale. Default 3.
#' @param events Optional character vector of events to show (matrix style).
#' @param drop_insignificant_events For the matrix style, show only events in
#' at least one significant pair. Default TRUE.
#' @param triangle For the matrix style, show only the lower triangle.
#' Default NULL: TRUE in genome order, FALSE with `cluster = TRUE` (a full
#' matrix with dendrograms reads better for clusters).
#' @param cluster For the matrix style, cluster events by their association
#' profiles instead of ordering them along the genome. See Details. Default
#' FALSE.
#' @param cluster_on What to cluster on: "signed_p" (default; -log10 p with
#' the sign of the log odds ratio, capped at 10) or "log_odds_ratio" (capped
#' at `max_log_odds_ratio`).
#' @param cluster_k Optional number of clusters to cut the dendrogram into;
#' the matrix is split into blocks marked by numbered, coloured bars, and
#' each event's cluster is returned.
#' @param max_neg_log_q For the matrix style, the -log10 q at which circles
#' reach full size. Default 10 (q = 1e-10); smaller q values are drawn at
#' full size. Without a cap, a few extremely significant pairs in a large
#' cohort would shrink every other circle.
#' @param show_nonsignificant For the matrix style, draw faint tiles for
#' tested pairs that aren't significant. Default TRUE.
#' @param cell_size Size of each matrix cell in mm. Default 4. The matrix
#' is drawn at this fixed size, and a message gives the device size needed if
#' the figure doesn't fit.
#' @param font_size Font size of event names (matrix) or pair labels
#' (volcano, compare). Default 7.
#' @param chromosome_label_fontsize Font size of the chromosome strip labels.
#' Default 6.
#' @param label_top Number of significant pairs to label in the volcano and
#' compare styles, by q value. Default 10.
#' @param axis_labels For `style = "compare"`, the x and y axis titles.
#' Default `c("log odds ratio (compare_to)", "log odds ratio (assoc)")`.
#'
#' @return For "matrix", invisibly, a list with the drawn ComplexHeatmap
#' object (heatmap) and a data frame (events) of the events in display order
#' with their chromosome and, with `cluster_k`, their cluster (numbered as
#' in the plot, from the top); for
#' "volcano" and "compare", a ggplot object.
#'
#' @import ComplexHeatmap ggplot2 ggrepel dplyr grid
#' @export
#'
#' @examples
#' \dontrun{
#' feature_info = GAMBLR.utils::get_cnv_features("grch37")
#' pga = GAMBLR.utils::calculate_pga(segs, centering = "median", cutoff = 0.1,
#'                                   per_chromosome = TRUE)
#' adjusted = GAMBLR.utils::pairwise_event_association(
#'   cnv_status, dlbcl_meta, feature_info = feature_info, burden = pga,
#'   method = "logistic")
#' unadjusted = GAMBLR.utils::pairwise_event_association(
#'   cnv_status, dlbcl_meta, feature_info = feature_info)
#'
#' pretty_event_association(adjusted)
#' pretty_event_association(adjusted, style = "volcano")
#'
#' # groups of events with similar association profiles
#' pretty_event_association(unadjusted, cluster = TRUE, cluster_k = 4)
#' pretty_event_association(adjusted, style = "compare",
#'                          compare_to = unadjusted,
#'                          axis_labels = c("Fisher (unadjusted)",
#'                                          "Logistic (PGA-adjusted)"))
#' }
pretty_event_association = function(assoc,
                                     style = "matrix",
                                     compare_to = NULL,
                                     q_threshold = 0.05,
                                     max_log_odds_ratio = 3,
                                     max_neg_log_q = 10,
                                     events = NULL,
                                     drop_insignificant_events = TRUE,
                                     triangle = NULL,
                                     cluster = FALSE,
                                     cluster_on = "signed_p",
                                     cluster_k = NULL,
                                     show_nonsignificant = TRUE,
                                     cell_size = 4,
                                     font_size = 7,
                                     chromosome_label_fontsize = 6,
                                     label_top = 10,
                                     axis_labels = c("log odds ratio (compare_to)",
                                                     "log odds ratio (assoc)")){
  if(!style %in% c("matrix", "volcano", "compare")){
    stop("style must be \"matrix\", \"volcano\" or \"compare\"")
  }
  .check_assoc_table = function(x, arg){
    required = c("event_a", "event_b", "log_odds_ratio", "q_value")
    missing_cols = required[!required %in% colnames(x)]
    if(length(missing_cols) > 0){
      stop(paste0(arg, " is missing columns: ",
                  paste(missing_cols, collapse = ", "),
                  ". Use the output of GAMBLR.utils::pairwise_event_association."))
    }
    as.data.frame(x)
  }
  assoc = .check_assoc_table(assoc, "assoc")
  cap = function(x) pmax(pmin(x, max_log_odds_ratio), -max_log_odds_ratio)
  neg_log_q = function(q) -log10(pmax(q, .Machine$double.xmin))
  pair_label = function(a, b) paste(a, "+", b)
  assoc_colours = c("co-occurring" = "#C0392B",
                    "mutually exclusive" = "#2471A3",
                    "not significant" = "grey70")

  if(style == "volcano"){
    df = assoc %>%
      dplyr::filter(!is.na(log_odds_ratio), !is.na(q_value)) %>%
      dplyr::mutate(x = cap(log_odds_ratio),
                    y = neg_log_q(q_value),
                    capped = abs(log_odds_ratio) > max_log_odds_ratio,
                    group = dplyr::case_when(
                      q_value >= q_threshold ~ "not significant",
                      log_odds_ratio > 0 ~ "co-occurring",
                      TRUE ~ "mutually exclusive"),
                    label = pair_label(event_a, event_b))
    labelled = df %>%
      dplyr::filter(q_value < q_threshold) %>%
      dplyr::arrange(q_value) %>%
      utils::head(label_top)
    p = ggplot(df, aes(x = x, y = y, colour = group)) +
      geom_hline(yintercept = neg_log_q(q_threshold), linetype = "dashed",
                 colour = "grey50") +
      geom_vline(xintercept = 0, colour = "grey80") +
      geom_point(aes(shape = capped), size = 1.8, alpha = 0.8) +
      ggrepel::geom_text_repel(data = labelled, aes(label = label),
                               size = font_size / ggplot2::.pt,
                               colour = "black", min.segment.length = 0,
                               max.overlaps = Inf, show.legend = FALSE) +
      scale_colour_manual(values = assoc_colours, name = NULL) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                         labels = c(`FALSE` = "within scale",
                                    `TRUE` = paste0("beyond +/-",
                                                    max_log_odds_ratio,
                                                    " (drawn at the limit)")),
                         name = NULL) +
      scale_x_continuous(limits = c(-max_log_odds_ratio, max_log_odds_ratio)) +
      labs(x = "log odds ratio", y = expression(-log[10]~q)) +
      theme_classic()
    return(p)
  }

  if(style == "compare"){
    if(is.null(compare_to)){
      stop("style = \"compare\" needs a second table in compare_to")
    }
    compare_to = .check_assoc_table(compare_to, "compare_to")
    pair_key = function(a, b) ifelse(a < b, paste(a, b, sep = "\r"),
                                     paste(b, a, sep = "\r"))
    x_df = data.frame(key = pair_key(compare_to$event_a, compare_to$event_b),
                      x_lor = compare_to$log_odds_ratio,
                      x_q = compare_to$q_value)
    df = assoc %>%
      dplyr::mutate(key = pair_key(event_a, event_b)) %>%
      dplyr::inner_join(x_df, by = "key") %>%
      dplyr::filter(!is.na(log_odds_ratio), !is.na(x_lor))
    if(nrow(df) == 0){
      stop("assoc and compare_to have no pairs in common")
    }
    df = df %>%
      dplyr::mutate(x = cap(x_lor),
                    y = cap(log_odds_ratio),
                    sig_x = !is.na(x_q) & x_q < q_threshold,
                    sig_y = !is.na(q_value) & q_value < q_threshold,
                    group = dplyr::case_when(
                      sig_x & sig_y ~ "significant in both",
                      sig_y ~ paste("significant in", axis_labels[2], "only"),
                      sig_x ~ paste("significant in", axis_labels[1], "only"),
                      TRUE ~ "not significant"),
                    best_q = pmin(q_value, x_q, na.rm = TRUE),
                    label = pair_label(event_a, event_b))
    compare_colours = stats::setNames(
      c("#7D3C98", "#C0392B", "#E59866", "grey75"),
      c("significant in both",
        paste("significant in", axis_labels[2], "only"),
        paste("significant in", axis_labels[1], "only"),
        "not significant"))
    labelled = df %>%
      dplyr::filter(sig_x | sig_y) %>%
      dplyr::arrange(best_q) %>%
      utils::head(label_top)
    df = df %>% dplyr::arrange(dplyr::desc(group == "not significant"))
    p = ggplot(df, aes(x = x, y = y, colour = group)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  colour = "grey50") +
      geom_hline(yintercept = 0, colour = "grey85") +
      geom_vline(xintercept = 0, colour = "grey85") +
      geom_point(size = 1.8, alpha = 0.8) +
      ggrepel::geom_text_repel(data = labelled, aes(label = label),
                               size = font_size / ggplot2::.pt,
                               colour = "black", min.segment.length = 0,
                               max.overlaps = Inf, show.legend = FALSE) +
      scale_colour_manual(values = compare_colours, name = NULL,
                          breaks = names(compare_colours)) +
      coord_equal(xlim = c(-max_log_odds_ratio, max_log_odds_ratio),
                  ylim = c(-max_log_odds_ratio, max_log_odds_ratio)) +
      labs(x = axis_labels[1], y = axis_labels[2]) +
      theme_classic()
    return(p)
  }

  # matrix style
  assoc = assoc[!is.na(assoc$log_odds_ratio) & !is.na(assoc$q_value), ,
                drop = FALSE]
  all_events = unique(c(assoc$event_a, assoc$event_b))
  if(!is.null(events)){
    unknown = events[!events %in% all_events]
    if(length(unknown) > 0){
      warning(paste("events not in assoc:", paste(unknown, collapse = ", ")))
    }
    all_events = all_events[all_events %in% events]
  }
  if(drop_insignificant_events){
    significant = assoc[assoc$q_value < q_threshold, , drop = FALSE]
    all_events = all_events[all_events %in%
                              c(significant$event_a, significant$event_b)]
  }
  if(length(all_events) < 2){
    stop("fewer than two events to show; try drop_insignificant_events = FALSE or a higher q_threshold")
  }
  assoc = assoc[assoc$event_a %in% all_events & assoc$event_b %in% all_events, ,
                drop = FALSE]

  # event coordinates, when pairwise_event_association had feature_info
  has_coords = all(c("chrom_a", "start_a", "chrom_b", "start_b") %in%
                     colnames(assoc))
  event_chrom = stats::setNames(rep(NA_character_, length(all_events)),
                                all_events)
  event_start = stats::setNames(rep(NA_real_, length(all_events)), all_events)
  if(has_coords){
    coords = rbind(
      data.frame(event = assoc$event_a, chrom = as.character(assoc$chrom_a),
                 start = as.numeric(assoc$start_a)),
      data.frame(event = assoc$event_b, chrom = as.character(assoc$chrom_b),
                 start = as.numeric(assoc$start_b)))
    coords = coords[!duplicated(coords$event), , drop = FALSE]
    event_chrom[coords$event] = gsub("^chr", "", coords$chrom)
    event_start[coords$event] = coords$start
  }
  chrom_levels = c(as.character(1:22), "X", "Y")
  chrom_rank = match(event_chrom, chrom_levels)
  chrom_rank[is.na(chrom_rank) & !is.na(event_chrom)] = length(chrom_levels) + 1
  chrom_rank[is.na(event_chrom)] = length(chrom_levels) + 2
  ordered = all_events[order(chrom_rank, event_start, all_events)]
  n = length(ordered)
  event_chrom = event_chrom[ordered]
  has_any_coords = any(!is.na(event_chrom))

  lor = matrix(NA_real_, n, n, dimnames = list(ordered, ordered))
  qv = lor
  ia = match(assoc$event_a, ordered)
  ib = match(assoc$event_b, ordered)
  lor[cbind(ia, ib)] = assoc$log_odds_ratio
  lor[cbind(ib, ia)] = assoc$log_odds_ratio
  qv[cbind(ia, ib)] = assoc$q_value
  qv[cbind(ib, ia)] = assoc$q_value
  same_chrom = outer(event_chrom, event_chrom, "==")
  same_chrom[is.na(same_chrom)] = FALSE
  lor_capped = cap(lor)

  # optional clustering on association profiles
  dendrogram = NULL
  display_pos = seq_len(n)
  if(cluster){
    if(!cluster_on %in% c("signed_p", "log_odds_ratio")){
      stop("cluster_on must be \"signed_p\" or \"log_odds_ratio\"")
    }
    if(cluster_on == "signed_p"){
      p_col = if("p_value" %in% colnames(assoc)) "p_value" else "q_value"
      score = pmin(-log10(pmax(assoc[[p_col]], 1e-300)), 10) *
        sign(assoc$log_odds_ratio)
      profile = matrix(0, n, n, dimnames = list(ordered, ordered))
      profile[cbind(ia, ib)] = score
      profile[cbind(ib, ia)] = score
    }else{
      profile = lor_capped
      profile[is.na(profile)] = 0
    }
    diag(profile) = 0
    similarity = suppressWarnings(stats::cor(profile))
    similarity[is.na(similarity)] = 0
    dendrogram = stats::hclust(stats::as.dist(1 - similarity),
                               method = "average")
    display_pos = match(seq_len(n), dendrogram$order)
  }
  if(is.null(triangle)){
    triangle = !cluster
  }

  col_fun = circlize::colorRamp2(
    c(-max_log_odds_ratio, 0, max_log_odds_ratio),
    c(assoc_colours[["mutually exclusive"]], "white",
      assoc_colours[["co-occurring"]]))
  sig_q = qv[!is.na(qv) & qv < q_threshold]
  size_max = min(max_neg_log_q,
                 max(2, ceiling(max(c(neg_log_q(sig_q), 0)))))
  size_min = neg_log_q(q_threshold)
  circle_fraction = function(q){
    0.35 + 0.65 * (pmin(neg_log_q(q), size_max) - size_min) /
      max(size_max - size_min, 1e-9)
  }
  untested_fill = "grey85"
  cell_fun = function(j, i, x, y, width, height, fill){
    if(i == j || (triangle && display_pos[j] > display_pos[i])){
      return(invisible(NULL))
    }
    q = qv[i, j]
    if(is.na(q)){
      if(same_chrom[i, j]){
        grid.rect(x, y, width, height,
                  gp = gpar(fill = untested_fill, col = "white", lwd = 0.5))
      }
      return(invisible(NULL))
    }
    value = lor_capped[i, j]
    grid.rect(x, y, width, height, gp = gpar(fill = NA, col = "grey92",
                                             lwd = 0.5))
    if(q < q_threshold){
      grid.circle(x, y, r = min(unit.c(width, height)) * circle_fraction(q) / 2,
                  gp = gpar(fill = col_fun(value), col = "grey30", lwd = 0.4))
    }else if(show_nonsignificant){
      grid.rect(x, y, width, height,
                gp = gpar(fill = grDevices::adjustcolor(col_fun(value),
                                                        alpha.f = 0.45),
                          col = "white", lwd = 0.5))
    }
  }

  # chromosome strips and grouping
  split = NULL
  top_anno = NULL
  left_anno = NULL
  if(has_any_coords){
    group = ifelse(is.na(event_chrom), "other", event_chrom)
    group_levels = unique(group)
    if(!cluster){
      split = factor(group, levels = group_levels)
    }
    chrom_cols = GAMBLR.helpers::get_gambl_colours("chromosome")
    names(chrom_cols) = gsub("^chr", "", names(chrom_cols))
    chrom_cols = c(chrom_cols, Y = "#7A6F5E", other = "white")
    missing_cols = setdiff(group_levels, names(chrom_cols))
    chrom_cols = c(chrom_cols,
                   stats::setNames(rep("grey60", length(missing_cols)),
                                   missing_cols))
    chrom_cols = chrom_cols[group_levels]
    runs = rle(group)
    run_end = cumsum(runs$lengths)
    run_mid = floor((run_end - runs$lengths + 1 + run_end) / 2)
    label_at = rep(NA_character_, n)
    if(cluster){
      # events are no longer grouped by chromosome, so label each one
      label_at = ifelse(group == "other", NA, group)
    }else{
      label_at[run_mid] = ifelse(runs$values == "other", NA, runs$values)
    }
    fill_rgb = grDevices::col2rgb(chrom_cols[group])
    luminance = (0.299 * fill_rgb[1, ] + 0.587 * fill_rgb[2, ] +
                   0.114 * fill_rgb[3, ]) / 255
    label_col = ifelse(luminance > 0.6, "black", "white")
    strip = function(which){
      anno_simple(group, col = chrom_cols, pch = label_at,
                  pt_gp = gpar(fontsize = chromosome_label_fontsize,
                               col = label_col),
                  which = which,
                  height = if(which == "column") unit(cell_size, "mm") else NULL,
                  width = if(which == "row") unit(cell_size, "mm") else NULL)
    }
  }

  # clusters, numbered in display order (top to bottom), as ComplexHeatmap
  # numbers the slices it cuts from the dendrogram
  event_cluster = NULL
  if(cluster && !is.null(cluster_k)){
    if(!is.numeric(cluster_k) || length(cluster_k) != 1 || cluster_k < 2 ||
       cluster_k > n){
      stop("cluster_k must be a single number between 2 and the number of events shown")
    }
    cluster_k = as.integer(cluster_k)
    raw_cluster = stats::cutree(dendrogram, k = cluster_k)
    first_seen = unique(raw_cluster[dendrogram$order])
    event_cluster = match(raw_cluster, first_seen)
    split = cluster_k
    cluster_cols = if(cluster_k <= 8){
      RColorBrewer::brewer.pal(max(3, cluster_k), "Dark2")[seq_len(cluster_k)]
    }else{
      grDevices::hcl.colors(cluster_k, "Dark 3")
    }
    block = function(which){
      anno_block(gp = gpar(fill = cluster_cols, col = NA),
                 labels = paste0("C", seq_len(cluster_k)),
                 labels_gp = gpar(col = "white", fontsize = font_size + 1,
                                  fontface = "bold"),
                 which = which,
                 height = if(which == "column") unit(4, "mm") else NULL,
                 width = if(which == "row") unit(4, "mm") else NULL)
    }
  }
  column_annos = list()
  row_annos = list()
  if(!is.null(event_cluster)){
    column_annos$cluster = block("column")
    row_annos$cluster = block("row")
  }
  if(has_any_coords){
    column_annos$chromosome = strip("column")
    row_annos$chromosome = strip("row")
  }
  if(length(column_annos) > 0){
    top_anno = do.call(HeatmapAnnotation,
                       c(column_annos, list(show_annotation_name = FALSE,
                                            gap = unit(0.5, "mm"))))
    left_anno = do.call(rowAnnotation,
                        c(row_annos, list(show_annotation_name = FALSE,
                                          gap = unit(0.5, "mm"))))
  }
  slice_gap = if(is.null(event_cluster)) unit(0.6, "mm") else unit(1.5, "mm")
  cluster_by = if(cluster) dendrogram else FALSE
  ht = Heatmap(lor_capped,
               name = "log_odds_ratio",
               col = col_fun,
               rect_gp = gpar(type = "none"),
               cell_fun = cell_fun,
               cluster_rows = cluster_by,
               cluster_columns = cluster_by,
               row_dend_width = unit(15, "mm"),
               column_dend_height = unit(15, "mm"),
               row_split = split,
               column_split = split,
               row_gap = slice_gap,
               column_gap = slice_gap,
               row_title = NULL,
               column_title = NULL,
               row_names_side = "left",
               column_names_side = "bottom",
               row_names_gp = gpar(fontsize = font_size),
               column_names_gp = gpar(fontsize = font_size),
               top_annotation = top_anno,
               left_annotation = left_anno,
               width = unit(cell_size * n, "mm"),
               height = unit(cell_size * n, "mm"),
               show_heatmap_legend = FALSE)

  legend_steps = pretty(c(size_min, size_max), n = 3)
  legend_q = unique(c(q_threshold,
                      10^-legend_steps[legend_steps > size_min &
                                         legend_steps <= size_max]))
  legend_q = sort(legend_q, decreasing = TRUE)
  legend_labels = format(legend_q, digits = 2)
  if(any(neg_log_q(sig_q) > size_max)){
    # circles for smaller q values are drawn at the largest size
    legend_labels[length(legend_labels)] = paste0(
      "<= ", legend_labels[length(legend_labels)])
  }
  legends = list(
    Legend(col_fun = col_fun, title = "log odds ratio",
           at = c(-max_log_odds_ratio, 0, max_log_odds_ratio),
           labels = c(paste0("-", max_log_odds_ratio, " (exclusive)"),
                      "0",
                      paste0(max_log_odds_ratio, " (co-occurring)"))),
    Legend(labels = legend_labels,
           title = "q value",
           type = "points", pch = 21,
           size = unit(cell_size * circle_fraction(legend_q), "mm"),
           legend_gp = gpar(fill = "grey60", col = "grey30"),
           grid_width = unit(cell_size, "mm"),
           grid_height = unit(cell_size, "mm"),
           background = NA))
  if(any(same_chrom & is.na(qv) & row(qv) != col(qv))){
    legends = c(legends, list(
      Legend(labels = "same chromosome,\nnot tested",
             legend_gp = gpar(fill = untested_fill), title = NULL)))
  }
  drawn = draw(ht, heatmap_legend_list = legends,
               heatmap_legend_side = "right")
  # cells have a fixed size, so the figure can be larger than the device
  needed = c(convertWidth(ComplexHeatmap:::width(drawn), "in", valueOnly = TRUE),
             convertHeight(ComplexHeatmap:::height(drawn), "in",
                           valueOnly = TRUE))
  available = grDevices::dev.size("in")
  if(any(needed > available + 0.01)){
    message(sprintf(paste0("The plot needs about %.1f x %.1f in but the ",
                           "device is %.1f x %.1f in, so it may be clipped. ",
                           "Use a larger device or a smaller cell_size."),
                    needed[1], needed[2], available[1], available[2]))
  }
  events_out = data.frame(event = ordered, chromosome = event_chrom,
                          stringsAsFactors = FALSE)
  if(!is.null(event_cluster)){
    events_out$cluster = event_cluster
  }
  display_order = if(cluster) dendrogram$order else seq_len(n)
  events_out = events_out[display_order, , drop = FALSE]
  rownames(events_out) = NULL
  invisible(list(heatmap = drawn, events = events_out))
}
