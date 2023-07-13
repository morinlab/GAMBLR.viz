#' @title Mutation Frequency Heatmap
#'
#' @description Count hypermutated bins and generate heatmap/cluster the data.
#'
#' @details This function takes a metadata table with `these_samples_metadata` parameter and internally calls [GAMBLR::calc_mutation_frequency_sliding_windows] (that internally calls [GAMBLR::get_ssm_by_regions])
#' to retrieve mutations for plotting. This plotting function has a variety of useful parameters, providing many customizable plotting options. For more details on how these parameters can be used,
#' and extended usage examples, refer to the SSM tutorial vignette section 1.4.9.
#'
#' @param regions Vector of regions in the format "chr:start-end".
#' @param regions_df Data frame of regions with four columns (chrom,start,end,gene_name).
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process (or full metadata).
#' @param region_padding How many bases will be added on the left and right of the regions to ensure any small regions are sufficiently covered by bins. Default is  1000.
#' @param seq_type The seq_type you want back, default is genome.
#' @param metadataColumns What metadata will be shown in the visualization.
#' @param sortByColumns Which of the metadata to sort on for the heatmap.
#' @param expressionColumns Optional variable for retrieving expression values for a specific gene(s).
#' @param orientation Specify the sample orientation, default is sample_rows.
#' @param skip_regions Regions to be filtered out from the regions data frame. Only applies if `regions_df` is not provided. Default is MYC, BCL2 and IGLL5.
#' @param customColour Optional named list of named vectors for specifying all colours for metadata. Can be generated with map_metadata_to_colours. Default is NULL.
#' @param slide_by How far to shift before starting the next window.
#' @param window_size The width of your sliding window.
#' @param min_count_per_bin Minimum counts per bin, default is 3.
#' @param min_bin_recurrence How many samples a bin must be mutated in to retain in the visualization.
#' @param min_bin_patient How many bins must a patient mutated in to retain in the visualization.
#' @param region_fontsize Font size of regions in plot, default is 8ppt.
#' @param cluster_rows_heatmap Optional parameter to enable/disable clustering of each dimension of the heatmap. Default is FALSE.
#' @param cluster_cols_heatmap Optional parameter to enable/disable clustering of each dimension of the heatmap. Default is FALSE.
#' @param show_gene_colours Optional logical argument indicating whether regions should have associated colours plotted as annotation track of heatmap.
#' @param legend_row Fiddle with these to widen or narrow your legend.
#' @param legend_col Fiddle with these to widen or narrow your legend.
#' @param legend_direction Accepts one of "horizontal" (default) or "vertical" to indicate in which direction the legend will be drawn.
#' @param legendFontSize Font size of legend in plot, default is 10ppt.
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flat files (only works for streamlined data, not full MAF details).
#' @param mode Only works with indexed flat files. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#'
#' @return Nothing
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr tidyr circlize ComplexHeatmap tibble
#' @export
#'
#' @examples
#' #load metadata.
#' metadata = get_gambl_metadata()
#' dlbcl_bl_meta = dplyr::filter(metadata, pathology %in% c("DLBCL", "BL"))
#'
#' #bring together all derived sample-level results from many GAMBL pipelines.
#' dlbcl_bl_meta = collate_results(join_with_full_metadata = TRUE,
#'                                 these_samples_metadata = dlbcl_bl_meta)
#'
#' #get ashm regions
#' some_regions = grch37_ashm_regions
#'
#' get_mutation_frequency_bin_matrix(these_samples_metadata = dlbcl_bl_meta,
#'                                   regions_df = some_regions)
#'
get_mutation_frequency_bin_matrix = function(regions,
                                             regions_df,
                                             these_samples_metadata,
                                             seq_type = "genome",
                                             region_padding = 1000,
                                             metadataColumns = c("pathology"),
                                             sortByColumns = c("pathology"),
                                             expressionColumns = c(),
                                             orientation = "sample_rows",
                                             skip_regions = c("MYC", "BCL2", "IGLL5"),
                                             customColour = NULL,
                                             slide_by = 100,
                                             window_size = 500,
                                             min_count_per_bin = 3,
                                             min_bin_recurrence = 5,
                                             min_bin_patient = 0,
                                             region_fontsize = 8,
                                             cluster_rows_heatmap = FALSE,
                                             cluster_cols_heatmap = FALSE,
                                             show_gene_colours = FALSE,
                                             legend_row = 3,
                                             legend_col = 3,
                                             legend_direction = "horizontal",
                                             legendFontSize = 10,
                                             from_indexed_flatfile = TRUE,
                                             mode = "slms-3"){

  if(missing(regions)){
    if(missing(regions_df)){
      regions_df = grch37_ashm_regions #drop MYC and BCL2
      regions_df = grch37_ashm_regions %>%
        dplyr::filter(!gene %in% skip_regions)
    }
    regions = unlist(apply(regions_df, 1, function(x){paste0(x[1], ":", as.numeric(x[2]) - region_padding, "-", as.numeric(x[3]) + region_padding)})) #add some buffer around each
  }
  dfs = lapply(regions, function(x){calc_mutation_frequency_sliding_windows(
    this_region = x, drop_unmutated = TRUE,seq_type=seq_type,
    slide_by = slide_by, plot_type = "none", window_size = window_size,
    min_count_per_bin = min_count_per_bin, return_count = TRUE,
    metadata = these_samples_metadata,
    from_indexed_flatfile = from_indexed_flatfile, mode = mode)})

  all= do.call("rbind", dfs)

  #add a fake bin for one gene and make every patient not mutated in it (to fill gaps)
  fake = these_samples_metadata %>%
    dplyr::select(sample_id) %>%
    mutate(bin = "1_chrN") %>%
    mutate(mutated = 0)

  all = bind_rows(all, fake)
  completed = complete(all, sample_id, bin, fill = list(mutated = 0))
  widened = pivot_wider(completed, names_from = sample_id, values_from = mutated)
  widened_df = column_to_rownames(widened, var = "bin")

  if(length(expressionColumns)>0){
    these_samples_metadata = these_samples_metadata %>%
      mutate(across(all_of(expressionColumns), ~ trim_scale_expression(.x)))
  }
  meta_show = these_samples_metadata %>%
    select(sample_id, all_of(metadataColumns)) %>%
    arrange(across(all_of(sortByColumns))) %>%
    dplyr::filter(sample_id %in% colnames(widened_df)) %>%
    column_to_rownames(var = "sample_id")

  message(paste("starting with", length(colnames(widened_df)), "patients"))
  patients_show = colnames(widened_df)[which(colSums(widened_df)>= min_bin_patient)]
  message(paste("returning matrix with", length(patients_show), "patients"))
  meta_show = dplyr::filter(meta_show, rownames(meta_show) %in% patients_show)
  to_show = widened_df[which(rowSums(widened_df) > min_bin_recurrence), patients_show]

  bin_col_fun = colorRamp2(c(0, 3, 6, 9), c("white", "orange", "red", "purple"))
  to_show_t = t(to_show)
  meta_show_t = meta_show[rownames(to_show_t),]
  lg_cols = get_gambl_colours("lymphgen")
  path_fun = function(x){
    path_cols = get_gambl_colours("pathology")
    lg_cols = get_gambl_colours("lymphgen")

    return(unname(path_cols[x]))
  }
  path_cols = get_gambl_colours("pathology")

  #assign bins back to regions for better annotation
  assign_bins_to_region = function(bin_names, rdf){
    bin_df = data.frame(bin_name = bin_names)

    separated = bin_df %>%
      separate(bin_name, into = c("start", "chrom")) %>%
      mutate(start = as.integer(start)) %>%
      mutate(end = start + 1)

    separated$bin_name = bin_names
    colnames(rdf)[c(1:3)] = c("chrom", "start", "end")
    rdf = mutate(rdf, start = start - 1500) %>%
      mutate(end = end + 1500)

    regions.dt = as.data.table(rdf)

    setkey(regions.dt, chrom, start, end)
    bin.dt = as.data.table(separated)
    setkey(bin.dt, chrom, start, end)
    bin_overlapped = foverlaps(bin.dt, regions.dt, mult = "first") %>%
      as.data.frame() %>%
      select(bin_name, gene) %>%
      column_to_rownames(var = "bin_name")

    return(bin_overlapped)
  }

  #regions_df = grch37_ashm_regions
  if(is.null(customColour)){
    meta_cols = map_metadata_to_colours(metadataColumns, these_samples_metadata = meta_show, as_vector = F)

  }else{
    meta_cols = customColour
  }
  col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  for(exp in expressionColumns){
    meta_cols[[exp]] = col_fun
  }
  bin_annot = assign_bins_to_region(bin_names = colnames(to_show_t), rdf = regions_df)
  heatmap_legend_param = list(title = "Bin value", nrow = legend_row, ncol = legend_row, legend_direction = legend_direction, labels_gp = gpar(fontsize = legendFontSize))

  annotation_legend_param = list(nrow = legend_row, ncol = legend_col, direction = legend_direction, labels_gp = gpar(fontsize = legendFontSize))

  if(orientation == "sample_rows"){
    row_annot = HeatmapAnnotation(df = meta_show, show_legend = T, which = 'row', col = meta_cols, annotation_legend_param = annotation_legend_param)
    if(show_gene_colours){
      col_annot = HeatmapAnnotation(df = bin_annot, show_legend = F, which = 'col', annotation_legend_param = annotation_legend_param)
    }else{
      col_annot = HeatmapAnnotation(value = anno_empty(border = FALSE))
    }
    Heatmap(to_show_t[rownames(meta_show), rownames(bin_annot)],
            cluster_columns = cluster_cols_heatmap,
            cluster_rows = cluster_rows_heatmap,
            col = bin_col_fun,
            bottom_annotation = col_annot,
            left_annotation = row_annot,
            show_row_names = F,
            show_column_names = F,
            column_split = factor(bin_annot$gene),
            column_title_gp = gpar(fontsize = region_fontsize),
            column_title_rot = 90,
            row_title_gp = gpar(fontsize = 10),
            heatmap_legend_param = heatmap_legend_param)
  }else{
    col_annot = HeatmapAnnotation(df = meta_show, show_legend = T, which = 'col', col = meta_cols, annotation_legend_param = annotation_legend_param)
    if(show_gene_colours){
      row_annot = HeatmapAnnotation(df = bin_annot,show_legend = F, which = 'row', annotation_legend_param = annotation_legend_param)
    }else{
      row_annot = rowAnnotation(value = anno_empty(border = FALSE))
    }
    Heatmap(to_show[rownames(bin_annot),rownames(meta_show)],
            show_heatmap_legend = F,
            cluster_columns = cluster_rows_heatmap,
            cluster_rows = cluster_cols_heatmap,
            col = bin_col_fun,
            bottom_annotation = col_annot,
            left_annotation = row_annot,
            show_row_names = F,
            show_column_names = F,
            row_split = factor(bin_annot$gene),
            row_title_gp = gpar(fontsize = region_fontsize),
            row_title_rot = 0,
            column_title_gp = gpar(fontsize = 8),
            heatmap_legend_param = heatmap_legend_param)
  }
}
