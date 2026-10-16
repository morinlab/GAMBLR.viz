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
#' Gene and feature labels are drawn once, above the top group. Within each
#' group, rows are ordered as [pretty_CN_heatmap] would order them (e.g.
#' clustered, or sorted with the `sortBy*` options passed through `...`).
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
    out = do.call(pretty_CN_heatmap, args)
    ho = out$heatmap_object
    if(i > 1){
      ho@heatmap_param$show_heatmap_legend = FALSE
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
  invisible(NULL)
}
