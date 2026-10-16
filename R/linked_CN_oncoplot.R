#' Copy Number Heatmap Linked to an Oncoplot
#'
#' @description Draw a [pretty_CN_heatmap] and a [prettyOncoplot] of the
#' same samples together, with the samples in the same order in both.
#'
#' @details Two layouts are available:
#' \itemize{
#'   \item `layout = "stacked"`: the CN heatmap is rotated (samples as
#'     columns, bins as rows) and the oncoplot sits below it (or above, with
#'     `oncoplot_side = "top"`). Sample metadata tracks are shown on the
#'     oncoplot (`metadataColumns`); the CN heatmap only shows PGA.
#'   \item `layout = "side"`: the CN heatmap is drawn normally (samples as
#'     rows) and the oncoplot is drawn sideways to its right (or left), with
#'     genes as columns and a per-gene alteration bar on top. Its cells use
#'     the same glyphs and colours as [prettyOncoplot]; the oncoplot's other
#'     decorations (percentages, its own sample tracks) aren't shown.
#'     Sample metadata tracks are shown on the CN heatmap.
#' }
#' Sample order comes from `order_by`: "cn" (default) clusters the samples
#' by copy number, "oncoplot" uses the oncoplot's own mutation-based order.
#' `sort_by` sorts samples by any mix of gene mutation status and CN features
#' (e.g. `c("TP53", "17p_loss", "MYD88")`), with `order_by` breaking ties.
#' `cluster_by` instead clusters samples on mutation and CN feature status
#' together.
#'
#' With `split_by`, samples are split into groups and ordered within each
#' group. `split_style` controls how:
#' \itemize{
#'   \item "slices" (default): one plot split into slices, with one set of
#'     CN gain/loss and gene alteration summaries for all samples.
#'   \item "panels" (stacked layout only): one linked plot per group, side by
#'     side, each with its own CN gain/loss summary (as the fraction of the
#'     group's samples, on one shared scale) and per-gene alteration
#'     fractions (also on a shared scale). Rows line up across panels, the
#'     cells have the same width in every panel, genes keep the order given
#'     in `genes`, and one legend is drawn below.
#' }
#'
#' Only samples present in both `cn_state_matrix` and
#' `these_samples_metadata` are used. To correct for ploidy, build
#' `cn_state_matrix` with `adjust_for_ploidy = TRUE` (see the Details of
#' [pretty_CN_heatmap]).
#'
#' @param cn_state_matrix Copy number matrix with samples as rows and bins as
#' columns, e.g. from [GAMBLR.utils::segmented_data_to_cn_matrix].
#' @param these_samples_metadata Metadata with one row per sample_id.
#' @param maf_df Mutations for the oncoplot, as for [prettyOncoplot].
#' @param genes Genes to show in the oncoplot.
#' @param layout "stacked" (default) or "side". See Details.
#' @param oncoplot_side Where the oncoplot goes: "bottom" (default) or "top"
#' for the stacked layout, "right" (default) or "left" for the side layout.
#' @param order_by What orders the samples: "cn" (default) or "oncoplot".
#' With `sort_by`, it only breaks ties.
#' @param sort_by Optional vector of genes and/or CN feature names to sort
#' samples by, in priority order. Genes shown in the oncoplot sort by mutation
#' status (mutated first); other names must be columns of `cnv_status` and
#' sort by event level (strongest first; see `sortByFeatures` in
#' [pretty_CN_heatmap]). Remaining ties are ordered by `order_by`. With
#' `split_by`, samples are sorted within each group.
#' @param cluster_by Optional: cluster samples on mutation and CN feature
#' status together instead of on the CN bins. Either "all" (every gene in
#' the oncoplot plus every column of `cnv_status`) or a vector of genes
#' and/or CN feature names. By default each feature counts as present or
#' absent and samples are clustered with binary (Jaccard) distance, so
#' samples are similar when they share alterations. The resulting dendrogram
#' is shown on the CN heatmap. Can't be combined with `sort_by`; `order_by`
#' is ignored. With `split_by`, each group is clustered separately.
#' @param cluster_levels Set to TRUE to cluster on CN event levels (0/1/2,
#' mutations 0/1) with Euclidean distance instead, so high-level events
#' weigh more. Default FALSE.
#' @param cnv_status Sample-by-feature CN calls with sample IDs as row names,
#' e.g. from `GAMBLR.utils::segmented_data_to_cnv_status()`. Needed when
#' `sort_by` includes CN features.
#' @param split_by Optional metadata column used to split the samples into
#' groups.
#' @param split_order Optional vector giving the groups to show, in order.
#' Defaults to the factor levels of `split_by`, or its sorted values.
#' @param split_style "slices" (default) or "panels". See Details.
#' @param metadataColumns,numericMetadataColumns Sample annotations to show,
#' once (see Details).
#' @param labelTheseFeatures,labelTheseGenes Labels for the CN heatmap, as in
#' [pretty_CN_heatmap]. With `split_style = "panels"` they are shown on the
#' first panel.
#' @param cn_size,oncoplot_size Size in cm of each panel along the gene/bin
#' axis: heights in the stacked layout, widths in the side layout. NULL
#' leaves the panel's own default.
#' @param width Width in cm of both panels in the stacked layout, i.e. the
#' length of the sample axis. Default NULL fills the device. With
#' `split_style = "panels"` it is the total width of the heatmaps across all
#' panels (default 25), shared out in proportion to group size.
#' @param height Height in cm of both panels in the side layout, i.e. the
#' length of the sample axis. Default NULL fills the device.
#' @param cn_args,oncoplot_args Optional named lists of further arguments for
#' [pretty_CN_heatmap] and [prettyOncoplot].
#' @param legend_position Where to put the legends. Default "bottom".
#' @param show_legends Set to FALSE to omit all legends. Default TRUE.
#' @param return_data Set to TRUE to return the plot components instead of
#' drawing them.
#'
#' @return Draws the plot. With `return_data = TRUE`, a list with
#' `heatmap_list` and `draw_args`; draw it with
#' `do.call(ComplexHeatmap::draw, c(list(x$heatmap_list), x$draw_args))`.
#' With `split_style = "panels"`, a list with `panels` (one such
#' `heatmap_list`/`draw_args` pair per group), `panel_widths_cm` and
#' `legend`.
#'
#' @import ComplexHeatmap dplyr grid
#' @export
#'
#' @examples
#' \dontrun{
#' meta = check_and_clean_metadata(dlbcl_meta, duplicate_action = "keep_first")
#' cn_mat = GAMBLR.utils::segmented_data_to_cn_matrix(segs,
#'                                                    these_samples_metadata = meta,
#'                                                    adjust_for_ploidy = TRUE)
#' maf = get_coding_ssm(these_samples_metadata = meta)
#' features = dplyr::filter(GAMBLR.utils::get_cnv_features("grch37"),
#'                          type != "arm")
#'
#' # CN heatmap rotated, oncoplot below, samples clustered by copy number
#' linked_CN_oncoplot(cn_mat, meta, maf,
#'                    genes = c("KMT2D", "CREBBP", "BCL2", "MYD88", "TP53"),
#'                    metadataColumns = c("COO_consensus", "lymphgen"),
#'                    labelTheseFeatures = features)
#'
#' # one linked plot per COO group, side by side, each with its own summaries
#' linked_CN_oncoplot(cn_mat, meta, maf,
#'                    genes = c("KMT2D", "CREBBP", "BCL2", "MYD88", "TP53"),
#'                    split_by = "COO_consensus", split_style = "panels",
#'                    metadataColumns = "lymphgen",
#'                    labelTheseFeatures = features,
#'                    cn_args = list(hide_these_chromosomes = c("X", "Y")))
#'
#' # CN heatmap upright with a sideways oncoplot on the right, split by COO
#' linked_CN_oncoplot(cn_mat, meta, maf,
#'                    genes = c("KMT2D", "CREBBP", "BCL2", "MYD88", "TP53"),
#'                    layout = "side", split_by = "COO_consensus",
#'                    metadataColumns = c("COO_consensus", "lymphgen"))
#' }
linked_CN_oncoplot = function(cn_state_matrix,
                              these_samples_metadata,
                              maf_df,
                              genes,
                              layout = "stacked",
                              oncoplot_side,
                              order_by = "cn",
                              sort_by = NULL,
                              cluster_by = NULL,
                              cluster_levels = FALSE,
                              cnv_status = NULL,
                              split_by,
                              split_order,
                              split_style = "slices",
                              metadataColumns = c("pathology"),
                              numericMetadataColumns,
                              labelTheseFeatures,
                              labelTheseGenes,
                              cn_size = NULL,
                              oncoplot_size = NULL,
                              width = NULL,
                              height = NULL,
                              cn_args = list(),
                              oncoplot_args = list(),
                              legend_position = "bottom",
                              show_legends = TRUE,
                              return_data = FALSE){
  if(!layout %in% c("stacked", "side")){
    stop("layout must be \"stacked\" or \"side\"")
  }
  if(missing(oncoplot_side)){
    oncoplot_side = if(layout == "stacked") "bottom" else "right"
  }
  allowed_sides = if(layout == "stacked") c("bottom", "top") else c("right", "left")
  if(!oncoplot_side %in% allowed_sides){
    stop(paste0("oncoplot_side must be \"", allowed_sides[1], "\" or \"",
                allowed_sides[2], "\" for layout = \"", layout, "\""))
  }
  if(!order_by %in% c("cn", "oncoplot")){
    stop("order_by must be \"cn\" or \"oncoplot\"")
  }
  if(!is.null(sort_by) && !is.null(cluster_by)){
    stop("use either sort_by or cluster_by, not both")
  }
  if(!split_style %in% c("slices", "panels")){
    stop("split_style must be \"slices\" or \"panels\"")
  }
  if(any(duplicated(these_samples_metadata$sample_id))){
    stop(paste("these_samples_metadata has duplicate sample_id values;",
               "see check_and_clean_metadata()"))
  }
  has_split = !missing(split_by) && !is.null(split_by)
  if(has_split && !split_by %in% colnames(these_samples_metadata)){
    stop("split_by must name a column of these_samples_metadata")
  }
  if(split_style == "panels" && (!has_split || layout != "stacked")){
    stop("split_style = \"panels\" needs split_by and layout = \"stacked\"")
  }

  not_in_metadata = setdiff(rownames(cn_state_matrix), these_samples_metadata$sample_id)
  if(length(not_in_metadata) > 0){
    message(paste(length(not_in_metadata),
                  "sample(s) in cn_state_matrix but not in these_samples_metadata are left out"))
  }
  # samples present in both the CN matrix and the metadata
  meta = dplyr::filter(these_samples_metadata,
                       sample_id %in% rownames(cn_state_matrix))
  if(has_split){
    no_group = is.na(meta[[split_by]])
    if(any(no_group)){
      message(paste(sum(no_group), "sample(s) with no", split_by,
                    "value are left out"))
      meta = meta[!no_group, , drop = FALSE]
    }
  }
  dropped = setdiff(these_samples_metadata$sample_id, meta$sample_id)
  if(length(dropped) > 0){
    message(paste(length(dropped), "sample(s) without copy number data",
                  "(or a split_by value) are left out"))
  }
  if(nrow(meta) == 0){
    stop("no samples have both copy number data and metadata")
  }

  build_args = list(
    cn_state_matrix = cn_state_matrix, maf_df = maf_df, genes = genes,
    layout = layout, oncoplot_side = oncoplot_side, order_by = order_by,
    sort_by = sort_by, cluster_by = cluster_by,
    cluster_levels = cluster_levels, cnv_status = cnv_status,
    metadataColumns = metadataColumns,
    numericMetadataColumns = if(missing(numericMetadataColumns)) NULL
                             else numericMetadataColumns,
    labelTheseFeatures = if(missing(labelTheseFeatures)) NULL
                         else labelTheseFeatures,
    labelTheseGenes = if(missing(labelTheseGenes)) NULL else labelTheseGenes,
    cn_size = cn_size, oncoplot_size = oncoplot_size,
    width = width, height = height,
    cn_args = cn_args, oncoplot_args = oncoplot_args
  )

  if(split_style == "panels"){
    present = unique(as.character(meta[[split_by]]))
    if(!missing(split_order)){
      groups = intersect(split_order, present)
    }else if(is.factor(meta[[split_by]])){
      groups = intersect(levels(meta[[split_by]]), present)
    }else{
      groups = sort(present)
    }
    return(.linked_CN_oncoplot_panels(meta, split_by, groups, build_args,
                                      show_legends, return_data))
  }

  if(has_split && !missing(split_order)){
    meta = meta[as.character(meta[[split_by]]) %in% split_order, , drop = FALSE]
  }
  built = do.call(.build_linked_CN_oncoplot, c(list(meta = meta), build_args))
  draw_args = list(
    main_heatmap = built$main_heatmap,
    merge_legends = TRUE,
    show_heatmap_legend = show_legends,
    show_annotation_legend = show_legends,
    heatmap_legend_side = legend_position,
    annotation_legend_side = legend_position
  )
  if(has_split){
    split_levels = if(!missing(split_order)) split_order
                   else if(is.factor(meta[[split_by]])) levels(meta[[split_by]])
                   else sort(unique(as.character(meta[[split_by]])))
    groups = factor(as.character(built$meta[[split_by]]),
                    levels = intersect(split_levels, built$meta[[split_by]]))
    if(layout == "stacked"){
      draw_args$column_split = groups
    }else{
      draw_args$row_split = groups
    }
  }
  if(return_data){
    return(list(heatmap_list = built$heatmap_list, draw_args = draw_args))
  }
  do.call(draw, c(list(built$heatmap_list), draw_args))
  invisible(NULL)
}

# Build one linked CN heatmap + oncoplot (as a ComplexHeatmap list). `onco`
# can be a precomputed prettyOncoplot(return_inputs = TRUE) result, and
# `panel` carries the settings shared across split_style = "panels" groups.
.build_linked_CN_oncoplot = function(cn_state_matrix, meta, maf_df, genes,
                                     layout, oncoplot_side, order_by,
                                     sort_by, cluster_by, cluster_levels,
                                     cnv_status,
                                     metadataColumns, numericMetadataColumns,
                                     labelTheseFeatures, labelTheseGenes,
                                     cn_size, oncoplot_size, width, height,
                                     cn_args, oncoplot_args,
                                     onco = NULL, panel = NULL){
  # oncoplot first: it fixes the sample set and, optionally, the order
  if(is.null(onco)){
    onco = .linked_prettyOncoplot(maf_df, genes, meta, layout,
                                  metadataColumns, numericMetadataColumns,
                                  oncoplot_args)
  }
  onco_heatmap = onco$Heatmap
  samples = onco_heatmap@column_names_param$labels
  sort_status = NULL
  if(!is.null(sort_by)){
    sort_status = .linked_sort_status(onco, samples, sort_by, cnv_status)
  }
  meta = meta[match(samples, meta$sample_id), , drop = FALSE]
  cn_state_matrix = cn_state_matrix[samples, , drop = FALSE]

  # copy number heatmap for the same samples, in the same input order
  cn_call = utils::modifyList(cn_args, list(
    cn_state_matrix = cn_state_matrix,
    these_samples_metadata = meta,
    keep_sample_order = TRUE,
    return_data = TRUE,
    heatmap_name = "CN"
  ))
  if(!is.null(labelTheseFeatures)){
    cn_call$labelTheseFeatures = labelTheseFeatures
  }
  if(!is.null(labelTheseGenes)){
    cn_call$labelTheseGenes = labelTheseGenes
  }
  if(layout == "stacked"){
    cn_call$rotate = TRUE
    # the rotated CN heatmap has a fixed default width; set it explicitly
    # (NULL = fill the device) so both panels share the sample axis
    cn_call["width"] = list(width)
    # metadata tracks are shown once, on the oncoplot
    cn_call$metadataColumns = character(0)
  }else{
    cn_call$metadataColumns = metadataColumns
    if(!is.null(numericMetadataColumns)){
      cn_call$numericMetadataColumns = numericMetadataColumns
    }
  }
  combined_clustering = NULL
  if(!is.null(cluster_by)){
    combined_clustering = .linked_combined_clustering(onco, samples,
                                                      cluster_by, cnv_status,
                                                      cluster_levels)
  }
  if(order_by == "oncoplot" || !is.null(sort_status)){
    # the sample order comes from the oncoplot or from sort_by; nothing else
    # reorders the samples
    cn_call$cluster_rows = FALSE
  }
  if(!is.null(cluster_by)){
    # the CN heatmap draws the combined clustering (its data stay in the
    # oncoplot's sample order, which the clustering was computed on)
    cn_call$cluster_rows = if(is.null(combined_clustering)) FALSE
                           else combined_clustering
  }
  cn_out = do.call(pretty_CN_heatmap, cn_call)
  cn_heatmap = cn_out$heatmap_object
  if(!is.null(sort_status)){
    # Heatmaps combined in one list are matched by position, so both keep
    # their data in the same sample order and the sort is applied as the CN
    # heatmap's display order (it is the main heatmap; the oncoplot follows)
    if(!identical(rownames(cn_out$data), samples)){
      stop("the CN heatmap dropped or reordered samples, so sort_by can't be applied")
    }
    tie_order = if(order_by == "oncoplot"){
      onco_heatmap@column_order
    }else{
      hclust(dist(.neutral_for_clustering(cn_out$data)))$order
    }
    keys = lapply(sort_status[tie_order, , drop = FALSE], function(x){
      x = abs(as.numeric(x))
      x[is.na(x)] = 0
      -x
    })
    display_order = tie_order[do.call(order, c(keys, list(seq_along(tie_order))))]
    if(layout == "stacked"){
      cn_heatmap@column_order = display_order
    }else{
      cn_heatmap@row_order = display_order
    }
  }

  if(layout == "stacked"){
    if(!is.null(cn_size)){
      cn_heatmap@heatmap_param$height = unit(cn_size, "cm")
    }
    if(!is.null(panel)){
      # rebuild the oncoplot with shared-scale per-gene fractions, and gene
      # names only where asked, from the arguments prettyOncoplot used
      oa = onco$oncoprint_args
      oa$right_annotation = rowAnnotation(
        altered = anno_barplot(panel$gene_fractions,
                               gp = gpar(fill = panel$gene_colours, col = NA),
                               ylim = c(0, panel$gene_max),
                               border = FALSE,
                               # axis on the side away from the CN heatmap
                               axis_param = list(
                                 side = if(oncoplot_side == "top") "top" else "bottom"
                               ),
                               width = unit(1.5, "cm")),
        show_annotation_name = FALSE
      )
      oa$show_row_names = panel$show_gene_names
      oa$show_heatmap_legend = FALSE
      oa$width = unit(width, "cm")
      oa$height = unit(if(is.null(oncoplot_size)) 0.6 * length(genes)
                       else oncoplot_size, "cm")
      onco_panel = do.call(oncoPrint, oa)
    }else{
      if(!is.null(oncoplot_size)){
        onco_heatmap@heatmap_param$height = unit(oncoplot_size, "cm")
      }
      if(!is.null(width)){
        onco_heatmap@heatmap_param$width = unit(width, "cm")
      }
      onco_panel = onco_heatmap
    }
    heatmap_list = if(oncoplot_side == "bottom"){
      cn_heatmap %v% onco_panel
    }else{
      onco_panel %v% cn_heatmap
    }
  }else{
    # redraw the oncoplot's cells sideways: samples as rows, genes as columns
    oa = onco$oncoprint_args
    sideways_mat = if(is.list(oa$mat) && !is.matrix(oa$mat)){
      lapply(oa$mat, t)
    }else{
      t(oa$mat)
    }
    n_genes = length(onco$gene_order)
    onco_panel = oncoPrint(
      sideways_mat,
      alter_fun = oa$alter_fun,
      col = oa$col,
      column_order = onco$gene_order,
      # with order_by = "oncoplot", keep its mutation-based sample order
      row_order = if(order_by == "oncoplot" && is.null(sort_by) &&
                     is.null(cluster_by))
                    onco_heatmap@column_order else NULL,
      show_pct = FALSE,
      right_annotation = NULL,
      top_annotation = HeatmapAnnotation(
        genes_altered = anno_oncoprint_barplot(which = "column",
                                               height = unit(1.2, "cm"))
      ),
      show_row_names = FALSE,
      show_column_names = TRUE,
      column_names_side = "top",
      column_names_rot = 90,
      width = unit(if(is.null(oncoplot_size)) 0.5 * n_genes else oncoplot_size,
                   "cm"),
      heatmap_legend_param = list(title = "Alterations"),
      name = "oncoplot"
    )
    if(!is.null(cn_size)){
      cn_heatmap@heatmap_param$width = unit(cn_size, "cm")
    }
    if(!is.null(height)){
      cn_heatmap@heatmap_param$height = unit(height, "cm")
      onco_panel@heatmap_param$height = unit(height, "cm")
    }
    heatmap_list = if(oncoplot_side == "right"){
      cn_heatmap + onco_panel
    }else{
      onco_panel + cn_heatmap
    }
  }
  list(heatmap_list = heatmap_list,
       main_heatmap = if(order_by == "cn" || !is.null(sort_by) ||
                         !is.null(cluster_by)) "CN"
                      else onco_panel@name,
       meta = meta,
       onco = onco)
}

# Missing bins (e.g. fill_missing_with = "nothing") count as neutral when
# computing distances for ordering.
.neutral_for_clustering = function(m){
  m = as.matrix(m)
  m[is.na(m)] = 2
  m
}

# Clustering of samples on mutation and CN feature status together, for
# cluster_by. NULL when there are too few samples to cluster.
.linked_combined_clustering = function(onco, samples, cluster_by, cnv_status,
                                       cluster_levels){
  if(length(samples) < 2){
    return(NULL)
  }
  mat = onco$oncoprint_args$mat
  plotted_genes = if(is.list(mat) && !is.matrix(mat)) rownames(mat[[1]]) else rownames(mat)
  if(identical(cluster_by, "all")){
    cluster_by = c(plotted_genes, colnames(cnv_status))
    if(is.null(cnv_status)){
      message("cluster_by = \"all\" without cnv_status uses mutation status only")
    }
  }
  status = .linked_sort_status(onco, samples, cluster_by, cnv_status,
                               arg_name = "cluster_by")
  values = as.matrix(as.data.frame(lapply(status, function(x){
    x = abs(as.numeric(x))
    x[is.na(x)] = 0
    if(cluster_levels) x else as.numeric(x > 0)
  })))
  rownames(values) = samples
  distances = if(cluster_levels) dist(values) else dist(values, method = "binary")
  hclust(distances)
}

# Sorting table for sort_by: mutation status (0/1) for genes in the oncoplot,
# CN feature levels from cnv_status for everything else (samples as rows).
.linked_sort_status = function(onco, samples, sort_by, cnv_status,
                               arg_name = "sort_by"){
  mat = onco$oncoprint_args$mat
  plotted_genes = if(is.list(mat) && !is.matrix(mat)) rownames(mat[[1]]) else rownames(mat)
  status = data.frame(row.names = samples)
  for(key in sort_by){
    if(key %in% plotted_genes){
      altered = if(is.list(mat) && !is.matrix(mat)){
        Reduce(`|`, lapply(mat, function(m) m[key, samples] != 0 & !is.na(m[key, samples])))
      }else{
        !is.na(mat[key, samples]) & mat[key, samples] != ""
      }
      status[[key]] = as.integer(altered)
    }else{
      if(is.null(cnv_status) || !key %in% colnames(cnv_status)){
        stop(paste0(arg_name, " entry \"", key, "\" is neither a gene in the ",
                    "oncoplot nor a column of cnv_status"))
      }
      status[[key]] = as.data.frame(cnv_status)[samples, key]
    }
  }
  status
}

.linked_prettyOncoplot = function(maf_df, genes, meta, layout,
                                  metadataColumns, numericMetadataColumns,
                                  oncoplot_args){
  onco_call = utils::modifyList(oncoplot_args, list(
    maf_df = maf_df,
    genes = genes,
    these_samples_metadata = meta,
    return_inputs = TRUE
  ))
  if(layout == "stacked"){
    onco_call$metadataColumns = metadataColumns
    if(!is.null(numericMetadataColumns)){
      onco_call$numericMetadataColumns = numericMetadataColumns
    }
  }
  do.call(prettyOncoplot, onco_call)
}

# Fraction of samples with each alteration type, per gene (genes x types),
# from the alteration matrix prettyOncoplot gave to oncoPrint.
.linked_gene_fractions = function(oa){
  mat = oa$mat
  types = if(!is.null(oa$pct_include)) oa$pct_include else names(oa$col)
  if(is.list(mat) && !is.matrix(mat)){
    types = intersect(types, names(mat))
    n = ncol(mat[[1]])
    fractions = sapply(types, function(t) rowSums(mat[[t]] != 0 & !is.na(mat[[t]])) / n)
  }else{
    n = ncol(mat)
    cells = strsplit(ifelse(is.na(mat), "", mat), ";")
    types = intersect(types, unique(unlist(cells)))
    fractions = sapply(types, function(t){
      hit = vapply(cells, function(x) t %in% x, logical(1))
      rowSums(matrix(hit, nrow = nrow(mat))) / n
    })
  }
  fractions = matrix(fractions, nrow = nrow(mat),
                     dimnames = list(rownames(mat), types))
  fractions
}

.linked_CN_oncoplot_panels = function(meta, split_by, groups, build_args,
                                      show_legends, return_data){
  cn_args = build_args$cn_args
  oncoplot_args = build_args$oncoplot_args
  # genes keep the order given, so gene rows line up across panels
  oncoplot_args$keepGeneOrder = TRUE

  # first pass: each group's oncoplot, and the scales the panels will share
  group_meta = list()
  oncos = list()
  gene_fractions = list()
  bar_max = 0
  pga_max = 0
  for(g in groups){
    m = meta[as.character(meta[[split_by]]) == g, , drop = FALSE]
    onco = .linked_prettyOncoplot(build_args$maf_df, build_args$genes, m,
                                  "stacked", build_args$metadataColumns,
                                  build_args$numericMetadataColumns,
                                  oncoplot_args)
    samples = onco$Heatmap@column_names_param$labels
    m = m[match(samples, m$sample_id), , drop = FALSE]
    group_meta[[g]] = m
    oncos[[g]] = onco
    gene_fractions[[g]] = .linked_gene_fractions(onco$oncoprint_args)
    probe = utils::modifyList(cn_args, list(
      cn_state_matrix = build_args$cn_state_matrix[samples, , drop = FALSE],
      these_samples_metadata = m,
      metadataColumns = character(0),
      rotate = TRUE, cluster_rows = FALSE, return_data = TRUE
    ))
    out = suppressMessages(do.call(pretty_CN_heatmap, probe))
    n_with_data = pmax(colSums(!is.na(out$data)), 1)
    bar_max = max(bar_max, out$cumulative_gain / n_with_data,
                  out$cumulative_loss / n_with_data, na.rm = TRUE)
    pga_max = max(pga_max, out$row_anno$PGA, na.rm = TRUE)
  }
  gene_max = max(unlist(lapply(gene_fractions, rowSums)), 1e-6)
  all_types = unique(unlist(lapply(gene_fractions, colnames)))
  type_colours = oncos[[1]]$oncoprint_args$col[all_types]
  pga_colours = circlize::colorRamp2(c(0, max(pga_max, 1e-6)),
                                     c("white", "darkorchid4"))

  # heatmap cells get the same width in every panel
  n_total = sum(vapply(group_meta, nrow, integer(1)))
  total_width = if(is.null(build_args$width)) 25 else build_args$width
  cell_width = total_width / n_total

  panels = list()
  for(i in seq_along(groups)){
    g = groups[i]
    m = group_meta[[g]]
    fractions = gene_fractions[[g]]
    # every panel stacks the same alteration types in the same order
    missing_types = setdiff(all_types, colnames(fractions))
    if(length(missing_types) > 0){
      fractions = cbind(fractions, matrix(0, nrow(fractions), length(missing_types),
                                          dimnames = list(NULL, missing_types)))
    }
    fractions = fractions[, all_types, drop = FALSE]
    panel_cn_args = utils::modifyList(list(
      cumulative_axis = "minimal",
      show_chromosome_strip = i == length(groups)
    ), cn_args)
    panel_cn_args = utils::modifyList(panel_cn_args, list(
      cumulative_as_fraction = TRUE,
      cumulative_ylim = c(0, bar_max),
      annotation_colours = utils::modifyList(
        if(is.null(cn_args$annotation_colours)) list() else cn_args$annotation_colours,
        list(PGA = pga_colours))
    ))
    args = build_args
    args$cn_args = panel_cn_args
    args$oncoplot_args = oncoplot_args
    args$width = nrow(m) * cell_width
    if(i > 1){
      args["labelTheseFeatures"] = list(NULL)
      args["labelTheseGenes"] = list(NULL)
    }
    built = do.call(.build_linked_CN_oncoplot, c(list(
      meta = m, onco = oncos[[g]],
      panel = list(gene_fractions = fractions,
                   gene_colours = unname(type_colours),
                   gene_max = gene_max,
                   show_gene_names = i == length(groups))
    ), args))
    panels[[g]] = list(
      heatmap_list = built$heatmap_list,
      draw_args = list(main_heatmap = built$main_heatmap,
                       show_heatmap_legend = FALSE,
                       show_annotation_legend = FALSE,
                       column_title = paste0(g, " (n=", nrow(m), ")"))
    )
  }

  # one legend for all panels
  legend = NULL
  if(show_legends){
    legends = list(
      Legend(col_fun = circlize::colorRamp2(c(1, 2, 5), c("blue", "white", "red")),
             title = "Copy Number", direction = "horizontal"),
      Legend(col_fun = pga_colours, title = "PGA", direction = "horizontal"),
      Legend(labels = all_types, legend_gp = gpar(fill = unname(type_colours)),
             title = "Alterations", nrow = 2)
    )
    categorical = setdiff(build_args$metadataColumns,
                          build_args$numericMetadataColumns)
    if(length(categorical) > 0){
      metadata_colours = map_metadata_to_colours(metadataColumns = categorical,
                                                 these_samples_metadata = meta)
      for(column in categorical){
        values = sort(unique(as.character(stats::na.omit(meta[[column]]))))
        cols = metadata_colours[[column]]
        if(is.null(cols) || !all(values %in% names(cols))){
          next
        }
        legends = c(legends, list(Legend(labels = values,
                                         legend_gp = gpar(fill = cols[values]),
                                         title = column, nrow = 2)))
      }
    }
    legend = do.call(packLegend, c(legends, list(direction = "horizontal")))
  }

  # measure each panel so the figure is sized to fit them exactly
  target_device = grDevices::dev.cur()
  sizes = vapply(panels, function(p){
    grDevices::pdf(NULL)
    # closing the scratch device must hand drawing back to the caller's device
    on.exit({
      grDevices::dev.off()
      if(target_device > 1){
        grDevices::dev.set(target_device)
      }
    })
    drawn = suppressMessages(do.call(draw, c(list(p$heatmap_list), p$draw_args)))
    c(convertWidth(ComplexHeatmap:::width(drawn), "cm", valueOnly = TRUE),
      convertHeight(ComplexHeatmap:::height(drawn), "cm", valueOnly = TRUE))
  }, numeric(2))
  widths = sizes[1, ]

  if(return_data){
    return(list(panels = panels, panel_widths_cm = widths, legend = legend))
  }
  grid.newpage()
  panel_height = max(sizes[2, ])
  legend_height = if(is.null(legend)) 0
                  else convertHeight(ComplexHeatmap:::height(legend), "cm",
                                     valueOnly = TRUE) + 0.5
  # the whole figure, centred on the page, with the legend just below
  pushViewport(viewport(
    width = unit(sum(widths), "cm"),
    height = unit(panel_height + legend_height, "cm"),
    layout = grid.layout(2, length(panels),
                         widths = unit(widths, "cm"),
                         heights = unit(c(panel_height, legend_height), "cm"))
  ))
  for(i in seq_along(panels)){
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i))
    do.call(draw, c(list(panels[[i]]$heatmap_list), panels[[i]]$draw_args,
                    list(newpage = FALSE)))
    upViewport()
  }
  if(!is.null(legend)){
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = seq_along(panels)))
    draw(legend)
    upViewport()
  }
  upViewport()
  invisible(NULL)
}
