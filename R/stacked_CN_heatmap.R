#' Stacked Copy Number Heatmap
#'
#' @description Draw one [pretty_CN_heatmap] per group of samples (e.g. per
#' LymphGen class), stacked vertically with the genomic bins aligned, so that
#' each group gets its own mean CN and cumulative gain/loss tracks.
#'
#' @details The groups share everything needed to compare them:
#' \itemize{
#'   \item the same bins, chosen once from the whole cohort with
#'     `drop_bin_if_sd_below`;
#'   \item gain/loss tracks shown as the fraction of each group's samples
#'     (`cumulative_as_fraction = TRUE`) on one y-axis scale;
#'   \item the same colours for every sample annotation, including PGA.
#' }
#' Gene and feature labels are drawn once, above the top group; with
#' `label_associations`, features significantly enriched in a group are also
#' labelled on that group's panel. Within each
#' group, rows are ordered as [pretty_CN_heatmap] would order them (e.g.
#' clustered, or sorted with the `sortBy*` options passed through `...`).
#'
#' With `metadata_summary = TRUE` (the default), a small bar is drawn above
#' each sample annotation column of every group, summarizing that group:
#' for categorical columns, a stacked bar of the proportion of samples with
#' each value (in the annotation colours, missing values light grey, the
#' first value at the top); for numeric columns (e.g. PGA), a box plot on a
#' scale shared by all groups (range printed above the top group's box
#' plot). The annotation names then move below the last group.
#'
#' Not supported: `rotate`, column (bin) clustering, `geneBoxPlot`,
#' `splitByBinState` and `label_style = "classic"`.
#'
#' **Correcting for ploidy.** Build `cn_state_matrix` with
#' `GAMBLR.utils::segmented_data_to_cn_matrix(..., adjust_for_ploidy = TRUE)`
#' and don't pass `scale_by_sample = TRUE`. Scaling after the matrix is made
#' turns bins with no data (filled with CN 2) into apparent deletions in
#' high-ploidy samples, and combining both approaches shifts samples twice.
#' See the Details of [pretty_CN_heatmap].
#'
#' @param cn_state_matrix Copy number matrix with samples as rows and bins as
#' columns, e.g. from [GAMBLR.utils::segmented_data_to_cn_matrix]. Only
#' samples also in `these_samples_metadata` are shown. Build it
#' with `adjust_for_ploidy = TRUE` to correct for sample ploidy (see
#' Details).
#' @param these_samples_metadata Metadata with one row per sample_id
#' (see [GAMBLR.helpers::check_and_clean_metadata]) and the `split_by` column.
#' @param split_by Name of the metadata column that defines the groups.
#' Samples with a missing value are left out.
#' @param split_order Optional vector giving the groups to show, top to
#' bottom. Defaults to the factor levels of `split_by`, or its sorted values.
#' @param metadataColumns,numericMetadataColumns Sample annotations to show,
#' as in [pretty_CN_heatmap].
#' @param drop_bin_if_sd_below Bins whose standard deviation across the whole
#' cohort is at or below this value are dropped from every group. Default 0.
#' @param cumulative_as_fraction Show the gain/loss tracks as the fraction of
#' samples in each group (default TRUE) rather than counts.
#' @param cumulative_position Where each group's gain/loss tracks go:
#' "top" (default, above the heatmap, so labels point at the peaks) or
#' "bottom", as in [pretty_CN_heatmap].
#' @param show_mean_cn Show each group's mean CN strip. Default FALSE.
#' @param cumulative_axis Axis on the gain/loss tracks: "minimal" (default;
#' only 0 and the top of the shared scale), "full" or "none".
#' @param chromosome_strip Show the chromosome strip under the "last" group
#' only (default) or under "all" groups.
#' @param group_title_fontsize,group_title_rot Font size (default 9) and
#' rotation (default 90; use 0 for horizontal) of the group titles.
#' @param labelTheseGenes,labelTheseFeatures Labels drawn above the top group,
#' as in [pretty_CN_heatmap].
#' @param label_associations Optional output of
#' [GAMBLR.utils::test_feature_associations] comparing the groups (its
#' `comparison_column` values must be the `split_by` values). Each CNV
#' feature with `q_value` below `association_q` is labelled on the panel of
#' the group where it is more frequent (from the `n_mutated_<group>` and
#' `n_total_<group>` columns), above that group's gain/loss track, pointing at
#' the bin containing the feature's midpoint. Needs `cnv_features`.
#' @param cnv_features Feature table with the coordinates of the features in
#' `label_associations`, e.g. `GAMBLR.utils::get_cnv_features("grch37")`.
#' Features are matched by name, including gene-labelled names such as
#' "2p16.1~REL_gain".
#' @param association_q q value threshold for `label_associations`. Default
#' 0.05.
#' @param metadata_summary Draw a per-group summary bar above each sample
#' annotation column (see Details). Default TRUE.
#' @param metadata_summary_height Height of the summary bars in mm. Default
#' 12.
#' @param legend_position Where to put the legends. Default "bottom".
#' @param show_legends Set to FALSE to omit all legends. Default TRUE. When
#' drawing the returned `heatmap_list` yourself, pass
#' `show_heatmap_legend = FALSE, show_annotation_legend = FALSE` to
#' `ComplexHeatmap::draw()`.
#' @param group_gap Space between groups in mm. Default 5. When drawing the
#' returned `heatmap_list` yourself, pass `ht_gap = unit(5, "mm")` to
#' `ComplexHeatmap::draw()` for the same spacing.
#' @param return_data Set to TRUE to return the combined heatmap list and the
#' per-group output of [pretty_CN_heatmap] instead of drawing.
#' @param ... Other arguments passed to [pretty_CN_heatmap] for every group.
#'
#' @return Draws the plot. With `return_data = TRUE`, a list with
#' `heatmap_list` (draw it with `ComplexHeatmap::draw()`), `groups` and
#' `group_data` (the per-group [pretty_CN_heatmap] output).
#'
#' @import ComplexHeatmap dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' meta = check_and_clean_metadata(dlbcl_meta, duplicate_action = "keep_first")
#' # adjust for ploidy when building the matrix, not with scale_by_sample
#' cn_mat = GAMBLR.utils::segmented_data_to_cn_matrix(segs,
#'                                                    these_samples_metadata = meta,
#'                                                    adjust_for_ploidy = TRUE,
#'                                                    genome_build = "grch37")
#' stacked_CN_heatmap(cn_mat, meta,
#'                    split_by = "COO_consensus",
#'                    split_order = c("GCB", "ABC"),
#'                    metadataColumns = c("COO_consensus", "lymphgen"),
#'                    labelTheseFeatures = dplyr::filter(
#'                      GAMBLR.utils::get_cnv_features("grch37"), type != "arm"))
#'
#' # label CNV features enriched in each group (q < 0.05)
#' meta$Tumor_Sample_Barcode = meta$sample_id
#' assoc = GAMBLR.utils::test_feature_associations(
#'   metadata = meta, comparison_column = "COO_consensus",
#'   feature_matrix = cnv_status, sample_id_column = NULL)
#' stacked_CN_heatmap(cn_mat, meta, split_by = "COO_consensus",
#'                    label_associations = assoc,
#'                    cnv_features = GAMBLR.utils::get_cnv_features("grch37"))
#' }
stacked_CN_heatmap = function(cn_state_matrix,
                              these_samples_metadata,
                              split_by,
                              split_order,
                              metadataColumns = c("pathology"),
                              numericMetadataColumns,
                              drop_bin_if_sd_below = 0,
                              cumulative_as_fraction = TRUE,
                              cumulative_position = "top",
                              show_mean_cn = FALSE,
                              cumulative_axis = "minimal",
                              chromosome_strip = "last",
                              group_title_fontsize = 9,
                              group_title_rot = 90,
                              labelTheseGenes,
                              labelTheseFeatures,
                              label_associations = NULL,
                              cnv_features = NULL,
                              association_q = 0.05,
                              metadata_summary = TRUE,
                              metadata_summary_height = 12,
                              legend_position = "bottom",
                              group_gap = 5,
                              show_legends = TRUE,
                              return_data = FALSE,
                              ...){
  extra_args = list(...)
  unsupported = c("rotate", "geneBoxPlot", "splitByBinState", "return_data",
                  "heatmap_name", "heatmap_row_title", "cumulative_ylim")
  if(any(names(extra_args) %in% unsupported)){
    stop(paste("stacked_CN_heatmap does not support:",
               paste(intersect(names(extra_args), unsupported), collapse = ", ")))
  }
  if(isTRUE(extra_args$cluster_columns) || isTRUE(extra_args$cluster_regions)){
    stop("column (bin) clustering would misalign the groups; it isn't supported")
  }
  if(identical(extra_args$label_style, "classic")){
    stop("stacked_CN_heatmap only supports label_style = \"mark\"")
  }
  if(missing(split_by) || length(split_by) != 1 ||
     !split_by %in% colnames(these_samples_metadata)){
    stop("split_by must name one column of these_samples_metadata")
  }
  if(!chromosome_strip %in% c("last", "all")){
    stop("chromosome_strip must be \"last\" or \"all\"")
  }
  if(any(duplicated(these_samples_metadata$sample_id))){
    stop(paste("these_samples_metadata has duplicate sample_id values;",
               "see check_and_clean_metadata()"))
  }

  not_in_metadata = setdiff(rownames(cn_state_matrix), these_samples_metadata$sample_id)
  if(length(not_in_metadata) > 0){
    message(paste(length(not_in_metadata),
                  "sample(s) in cn_state_matrix but not in these_samples_metadata are left out"))
  }
  meta = dplyr::filter(these_samples_metadata,
                       sample_id %in% rownames(cn_state_matrix))
  no_group = is.na(meta[[split_by]])
  if(any(no_group)){
    message(paste(sum(no_group), "sample(s) with no", split_by,
                  "value are left out"))
    meta = meta[!no_group, , drop = FALSE]
  }
  present = unique(as.character(meta[[split_by]]))
  if(!missing(split_order)){
    absent = setdiff(split_order, present)
    if(length(absent) > 0){
      message(paste("no samples for:", paste(absent, collapse = ", ")))
    }
    groups = intersect(split_order, present)
  }else if(is.factor(meta[[split_by]])){
    groups = intersect(levels(meta[[split_by]]), present)
  }else{
    groups = sort(present)
  }
  if(length(groups) == 0){
    stop("no groups with samples to plot")
  }
  meta = meta[as.character(meta[[split_by]]) %in% groups, , drop = FALSE]

  # bins are chosen once for the whole cohort so they line up across groups
  cn_state_matrix = cn_state_matrix[meta$sample_id, , drop = FALSE]
  if(is.finite(drop_bin_if_sd_below)){
    bin_sd = apply(cn_state_matrix, 2, function(x) stats::sd(x, na.rm = TRUE))
    cn_state_matrix = cn_state_matrix[, which(bin_sd > drop_bin_if_sd_below),
                                      drop = FALSE]
  }

  has_numeric = !missing(numericMetadataColumns)
  # categorical colours from the whole cohort, so each value keeps its colour
  if(has_numeric){
    categorical = setdiff(metadataColumns, numericMetadataColumns)
  }else{
    categorical = metadataColumns
  }
  annotation_colours = map_metadata_to_colours(metadataColumns = categorical,
                                               these_samples_metadata = meta)
  if(!is.null(extra_args$annotation_colours)){
    annotation_colours = utils::modifyList(annotation_colours,
                                           extra_args$annotation_colours)
    extra_args$annotation_colours = NULL
  }

  # features enriched in each group, from a test_feature_associations table
  group_features = NULL
  if(!is.null(label_associations)){
    group_features = .enriched_features_by_group(label_associations, cnv_features,
                                                 groups, association_q)
  }

  group_args = function(g){
    m = meta[as.character(meta[[split_by]]) == g, , drop = FALSE]
    args = list(cn_state_matrix = cn_state_matrix[m$sample_id, , drop = FALSE],
                these_samples_metadata = m,
                metadataColumns = metadataColumns,
                drop_bin_if_sd_below = -Inf,
                cumulative_as_fraction = cumulative_as_fraction,
                cumulative_position = cumulative_position,
                show_mean_cn = show_mean_cn,
                cumulative_axis = cumulative_axis,
                heatmap_row_title_fontsize = group_title_fontsize,
                heatmap_row_title_rot = group_title_rot,
                legend_position = legend_position,
                return_data = TRUE)
    if(has_numeric){
      args$numericMetadataColumns = numericMetadataColumns
    }
    args = utils::modifyList(extra_args, args)
    # a single sample can't be clustered
    if(nrow(m) < 2){
      args$cluster_rows = FALSE
    }
    args
  }

  # first pass (no clustering) to find shared scales: the tallest gain/loss
  # bar and the PGA range across groups
  bar_max = 0
  pga_max = 0
  for(g in groups){
    args = group_args(g)
    args$cluster_rows = FALSE
    out = suppressMessages(do.call(pretty_CN_heatmap, args))
    bars = c(out$cumulative_gain, out$cumulative_loss)
    if(cumulative_as_fraction){
      # same denominator as pretty_CN_heatmap: samples with data in each bin
      n_with_data = pmax(colSums(!is.na(out$data)), 1)
      bars = bars / c(n_with_data, n_with_data)
    }
    bar_max = max(bar_max, bars, na.rm = TRUE)
    pga_max = max(pga_max, out$row_anno$PGA, na.rm = TRUE)
  }
  if(!"PGA" %in% names(annotation_colours)){
    annotation_colours$PGA = circlize::colorRamp2(c(0, max(pga_max, 1e-6)),
                                                  c("white", "darkorchid4"))
  }

  heatmaps = list()
  group_data = list()
  for(i in seq_along(groups)){
    g = groups[i]
    args = group_args(g)
    args$heatmap_name = paste0("CN_", i)
    # count on its own line so titles fit beside small groups
    args$heatmap_row_title = paste0(g, "\n(n=", nrow(args$cn_state_matrix), ")")
    args$show_chromosome_strip = chromosome_strip == "all" ||
      i == length(groups)
    args$cumulative_ylim = c(0, bar_max)
    args$annotation_colours = annotation_colours
    # labels once, above the top group
    if(i == 1){
      if(!missing(labelTheseGenes)){
        args$labelTheseGenes = labelTheseGenes
      }
      if(!missing(labelTheseFeatures)){
        args$labelTheseFeatures = labelTheseFeatures
      }
    }
    # features enriched in this group, above its gain/loss track
    if(!is.null(group_features[[g]]) && nrow(group_features[[g]]) > 0){
      args$labelTheseFeatures = unique(rbind(args$labelTheseFeatures,
                                             group_features[[g]]))
    }
    # with summary bars above the annotations, name them once, below
    if(metadata_summary){
      args$left_annotation_name_side = "bottom"
    }
    out = do.call(pretty_CN_heatmap, args)
    ho = out$heatmap_object
    if(i > 1){
      ho@heatmap_param$show_heatmap_legend = FALSE
    }
    if(metadata_summary && i < length(groups) && !is.null(ho@left_annotation)){
      for(n in names(ho@left_annotation@anno_list)){
        ho@left_annotation@anno_list[[n]]@name_param$show = FALSE
      }
      ho@left_annotation@extended = unit(c(0, 0, 0, 0), "mm")
    }
    heatmaps[[i]] = ho
    group_data[[g]] = out
  }
  heatmap_list = Reduce(`%v%`, heatmaps)

  if(return_data){
    return(list(heatmap_list = heatmap_list,
                groups = groups,
                group_data = group_data))
  }
  draw(heatmap_list,
       ht_gap = unit(group_gap, "mm"),
       merge_legends = TRUE,
       show_heatmap_legend = show_legends,
       show_annotation_legend = show_legends,
       heatmap_legend_side = legend_position,
       annotation_legend_side = legend_position)
  if(metadata_summary){
    .draw_metadata_summary(heatmaps, group_data, metadata_summary_height)
  }
  invisible(NULL)
}

# After drawing, add a summary of each group's sample annotations above each
# annotation column: stacked proportions (categorical) or a box plot on a
# shared scale (numeric).
.draw_metadata_summary = function(heatmaps, group_data, height, gap = 1){
  annos = lapply(heatmaps, function(ho) ho@left_annotation)
  if(all(vapply(annos, is.null, logical(1)))){
    return(invisible(NULL))
  }
  values = function(i, n){
    v = group_data[[i]]$row_anno[[n]]
    if(is.null(v)) NULL else v
  }
  # shared scales for numeric annotations
  ranges = list()
  for(i in seq_along(heatmaps)){
    if(is.null(annos[[i]])) next
    for(n in names(annos[[i]]@anno_list)){
      v = values(i, n)
      if(is.numeric(v)){
        ranges[[n]] = range(c(ranges[[n]], v), na.rm = TRUE)
      }
    }
  }
  for(i in seq_along(heatmaps)){
    if(is.null(annos[[i]])) next
    for(n in names(annos[[i]]@anno_list)){
      v = values(i, n)
      cm = annos[[i]]@anno_list[[n]]@color_mapping
      if(is.null(v) || is.null(cm)) next
      ok = tryCatch({
        grid::seekViewport(paste0("heatmap_", heatmaps[[i]]@name))
        grid::downViewport(paste0("annotation_", n, "_1"))
        TRUE
      }, error = function(e) FALSE)
      if(!ok){
        grid::upViewport(0)
        next
      }
      grid::pushViewport(grid::viewport(y = unit(1, "npc") + unit(gap, "mm"),
                                        height = unit(height, "mm"),
                                        just = "bottom", clip = "off"))
      if(is.numeric(v)){
        r = ranges[[n]]
        if(diff(r) == 0) r = r + c(-0.5, 0.5)
        grid::pushViewport(grid::viewport(yscale = r))
        v = v[!is.na(v)]
        if(length(v) > 0){
          grid::grid.rect(gp = gpar(col = "grey80", fill = NA, lwd = 0.5))
          ComplexHeatmap::grid.boxplot(v, pos = 0.5, direction = "vertical",
                                       box_width = 0.7, outline = FALSE,
                                       gp = gpar(fill = map_to_colors(cm, stats::median(v)),
                                                 lwd = 0.5))
        }
        if(i == 1){
          grid::grid.text(paste0(signif(r[1], 2), "-", signif(r[2], 2)),
                          y = unit(1, "npc") + unit(0.5, "mm"), rot = 90,
                          just = "left", gp = gpar(fontsize = 5))
        }
        grid::popViewport()
      }else{
        lev = cm@levels
        cols = cm@colors[lev]
        counts = table(factor(as.character(v), levels = lev))
        n_na = sum(is.na(v) | !as.character(v) %in% lev)
        prop = c(as.numeric(counts), n_na) / max(length(v), 1)
        fill = c(unname(cols), "grey90")
        keep = prop > 0
        top = 1 - c(0, cumsum(prop))[seq_along(prop)]
        grid::grid.rect(y = top[keep], height = prop[keep], just = "top",
                        gp = gpar(fill = fill[keep], col = NA))
        grid::grid.rect(gp = gpar(col = "grey40", fill = NA, lwd = 0.5))
      }
      grid::upViewport(0)
    }
  }
  invisible(NULL)
}

# For a test_feature_associations() table, the cnv_features rows enriched in
# each group (q < association_q; the group with the higher frequency).
.enriched_features_by_group = function(associations, cnv_features, groups,
                                       association_q){
  if(is.null(cnv_features)){
    stop("label_associations needs cnv_features (e.g. GAMBLR.utils::get_cnv_features(\"grch37\"))")
  }
  associations = as.data.frame(associations)
  if(!all(c("feature", "q_value") %in% colnames(associations))){
    stop("label_associations must be a test_feature_associations() table (with feature and q_value columns)")
  }
  freq_cols = paste0("n_mutated_", groups)
  total_cols = paste0("n_total_", groups)
  missing_cols = setdiff(c(freq_cols, total_cols), colnames(associations))
  if(length(missing_cols) > 0){
    stop(paste0("label_associations has no counts for the groups shown (missing ",
                paste(missing_cols, collapse = ", "),
                "); its comparison_column values must be the split_by values"))
  }
  significant = associations[!is.na(associations$q_value) &
                               associations$q_value < association_q, , drop = FALSE]
  cnv_features = as.data.frame(cnv_features)
  # match names, allowing gene-labelled ones (region~GENE_direction)
  plain = sub("~[^~]*(_(gain|loss))$", "\\1", significant$feature)
  row = match(significant$feature, cnv_features$name)
  row[is.na(row)] = match(plain[is.na(row)], cnv_features$name)
  if(anyNA(row)){
    message(sum(is.na(row)), " significant feature(s) are not in cnv_features and aren't labelled: ",
            paste(utils::head(significant$feature[is.na(row)], 10), collapse = ", "))
  }
  significant = significant[!is.na(row), , drop = FALSE]
  row = row[!is.na(row)]
  if(nrow(significant) == 0){
    return(stats::setNames(vector("list", length(groups)), groups))
  }
  freq = as.matrix(significant[, freq_cols, drop = FALSE]) /
    pmax(as.matrix(significant[, total_cols, drop = FALSE]), 1)
  enriched_in = groups[max.col(freq, ties.method = "first")]
  # a feature tested in several comparisons is labelled once per group
  keep = !duplicated(paste(row, enriched_in))
  lapply(stats::setNames(groups, groups), function(g){
    cnv_features[unique(row[keep & enriched_in == g]), , drop = FALSE]
  })
}
