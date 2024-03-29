#' @title Heatmap of mutation counts across sliding windows for multiple regions.
#'
#' @description Obtain a heatmap of mutation counts across sliding windows for multiple regions.
#'
#' @details This function takes a metadata table with `these_samples_metadata` parameter and internally calls [GAMBLR::calc_mutation_frequency_bin_region] (that internally calls [GAMBLR::get_ssm_by_regions]).
#' to retrieve mutation counts for sliding windows across one or more regions and generate a heatmap. May optionally provide any combination of a maf data frame, existing metadata, or a regions data frame or named vector.
#'
#' @param regions_list Named vector of regions in the format c(name1 = "chr:start-end", name2 = "chr:start-end"). Only one (or none) between `regions_list` and `regions_bed` arguments should be provided. If neither regions_list nor regions_bed is specified, the function will use GAMBLR aSHM region information.
#' @param regions_bed Data frame of regions with four columns (chrom, start, end, name). Only one (or nome) between `regions_list` and `regions_bed` arguments should be provided.
#' @param these_samples_metadata Metadata with at least sample_id column. If not providing a maf data frame, seq_type is also required.
#' @param these_sample_ids Vector of sample IDs. Metadata will be subset to sample IDs present in this vector.
#' @param this_seq_type Optional vector of seq_types to include in heatmap. Default c("genome", "capture"). Uses default seq_type priority for samples with >1 seq_type.
#' @param maf_data Optional maf data frame. Will be subset to rows where Tumor_Sample_Barcode matches provided sample IDs or metadata table. If not provided, maf data will be obtained with get_ssm_by_regions().
#' @param mut_freq_matrix Optional matrix of binned mutation frequencies generated outside of this function, usually by [GAMBLR::calc_mutation_frequency_bin_regions].
#' @param projection Genome build the function will operate in. Ensure this matches your provided regions and maf data for correct chr prefix handling. Default grch37.
#' @param region_padding Amount to pad the start and end coordinates by. Default 1000
#' @param drop_unmutated Whether to drop bins with 0 mutations. If returning a matrix format, this will only drop bins with no mutations in any samples.
#' @param skip_regions Optional character vector of genes to exclude from the default aSHM regions.
#' @param only_regions Optional character vector of genes to include from the default aSHM regions.
#' @param slide_by Slide size for sliding window. Default 100.
#' @param window_size Size of sliding window. Default 500.
#' @param metadataColumns Mandatory character vector of metadata columns to use in heatmap annotation. Default c("pathology").
#' @param sortByColumns Optional character vector of metadata columns to order annotations by. Will be ordered by factor levels and sorted in the order specified. Default NULL.
#' @param expressionColumns Optional character vector of numeric metadata columns, usually gene expression, for heatmap annotation.
#' @param orientation Specify whether heatmap should have samples in rows ("sample_rows") or in columns ("sample_cols"). Default sample_rows.
#' @param customColours Optional list of character vectors specifying colours for heatmap annotation with metadataColumns, e.g. list(pathology = c(DLBCL = "green", BL = "purple")). If left blank, the function will attempt to match heatmap annotations with existing colours from [GAMBLR::get_gambl_colours], or will default to the Blood colour palette.
#' @param naColour Colour to use for NA values in metadata/expression. Default "white".
#' @param backgroundColour Optionally specify the colour for heatmap bins with 0 mutations. Default grey90.
#' @param min_count_per_bin Specify the minimum number of mutations per bin to be included in the heatmap. Only bins with all samples falling below this threshold will be dropped. Default 0.
#' @param min_bin_recurrence Specify how many samples a bin must be mutated in to be displayed. Default 5.
#' @param min_mut_tumour Specify how many bins a tumour must be mutated in to be displayed. Default 0.
#' @param region_fontsize Fontsize of region labels on the heatmap. Default 8.
#' @param cluster_rows_heatmap Boolean. Default FALSE.
#' @param cluster_cols_heatmap Boolean.  Default FALSE.
#' @param show_gene_colours Boolean. Whether to add heatmap annotation colours for each region. Default FALSE.
#' @param label_regions_by Specify which feature of the regions to label the heatmap with. Heatmap will be split according to this value, and ordered by factor levels if the specified column is a factor. Default name.
#' @param label_regions_rotate Specify degree by which the label in the previous parameter will be rotated. Default 0 (no rotation). The accepted values are 0, 90, 270.
#' @param legend_row Control aesthetics of the heatmap legend. Default 3.
#' @param legend_col Control aesthetics of the heatmap legend. Default 3.
#' @param legend_direction Control aesthetics of the heatmap legend. Default "horizontal".
#' @param legendFontSize Control aesthetics of the heatmap legend. Default 10.
#' @param legend_side Control aesthetics of the heatmap legend. Default "bottom".
#' @param return_heatmap_obj Boolean. FALSE will plot the heatmap automatically. TRUE will return a heatmap object to allow further tweaking with the draw() function. Default FALSE.
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flat files (only works for streamlined data, not full MAF details).
#' @param mode Only works with indexed flat files. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#'
#'
#' @return A table of mutation counts for sliding windows across one or more regions. May be long or wide.
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose", "melt", "dcast"))
#' @import dplyr tidyr tibble ComplexHeatmap circlize grid
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#' library(dplyr)
#'
#' # get meta data
#' my_meta <- get_gambl_metadata() %>%
#'   filter(sample_id %in% c("DOHH-2", "OCI-Ly10", "OCI-Ly3", "SU-DHL-10", "SU-DHL-4"))
#'
#' # get ashm regions of a set of genes.
#' my_regions <- GAMBLR.data::somatic_hypermutation_locations_GRCh37_v_latest %>%
#'   rename("chrom" = "chr_name", "start" = "hg19_start", "end" = "hg19_end", "name" = "gene") %>%
#'   mutate(chrom = stringr::str_remove(chrom, "chr"))
#'
#' # create heatmap of mutation counts for the specified regions
#' meta_columns <- c("pathology", "lymphgen", "COO_consensus", "DHITsig_consensus")
#' heatmap_mutation_frequency_bin(
#'   regions_bed = my_regions,
#'   these_samples_metadata = my_meta,
#'   metadataColumns = meta_columns,
#'   sortByColumns = meta_columns
#' )
#'
heatmap_mutation_frequency_bin <- function(regions_list = NULL,
                                           regions_bed = NULL,
                                           these_samples_metadata = NULL,
                                           these_sample_ids = NULL,
                                           this_seq_type = c("genome", "capture"),
                                           maf_data,
                                           mut_freq_matrix,
                                           projection = "grch37",
                                           region_padding = 1000,
                                           drop_unmutated = FALSE,
                                           metadataColumns = c("pathology"),
                                           sortByColumns = NULL,
                                           expressionColumns = NULL,
                                           orientation = "sample_rows",
                                           skip_regions,
                                           only_regions,
                                           customColours = NULL,
                                           naColour = "white",
                                           backgroundColour = "grey90",
                                           slide_by = 100,
                                           window_size = 500,
                                           min_count_per_bin = 0,
                                           min_bin_recurrence = 5,
                                           min_mut_tumour = 0,
                                           region_fontsize = 8,
                                           cluster_rows_heatmap = FALSE,
                                           cluster_cols_heatmap = FALSE,
                                           show_gene_colours = FALSE,
                                           label_regions_by = "name",
                                           label_regions_rotate = 0,
                                           legend_row = 3,
                                           legend_col = 3,
                                           legend_direction = "horizontal",
                                           legendFontSize = 10,
                                           legend_side = "bottom",
                                           return_heatmap_obj = FALSE,
                                           from_indexed_flatfile = TRUE,
                                           mode = "slms-3") {
  # check arguments
  stopifnot(
    "Only one (or none) between regions_list and regions_bed arguments should be provided." =
      any(c(is.null(regions_list), is.null(regions_bed)))
  )

  # Get region specifications
  if (missing(skip_regions)) {
    skip_regions <- NULL
  }
  if (missing(only_regions)) {
    only_regions <- NULL
  }
  regions <- process_regions(
    regions_list = regions_list,
    regions_bed = regions_bed,
    region_padding = region_padding,
    skip_regions,
    only_regions
  )
  regions_bed <- regions$regions_bed
  regions <- regions$regions_list

  # Harmonize metadata and sample IDs
  metadata <- id_ease(
    these_samples_metadata,
    these_sample_ids,
    this_seq_type
  )

  these_sample_ids <- metadata$sample_id

  # Ensure all requested metadata columns are present in the metadata
  allMetaCols <- unique(c(metadataColumns, sortByColumns, expressionColumns))

  if (!min(allMetaCols %in% colnames(metadata))) {
    stop("Not all requested columns are present in the metadata table. ")
  }

  if (!missing(mut_freq_matrix)) {
    samples_in_colnames <- max(colnames(mut_freq_matrix) %in% these_sample_ids)
    samples_in_rownames <- max(rownames(mut_freq_matrix) %in% these_sample_ids)
    if (samples_in_colnames == 1) {
      all_matrix <- mut_freq_matrix
      matrix_samps <- colnames(mut_freq_matrix)
    } else if (samples_in_rownames == 1) {
      all_matrix <- t(mut_freq_matrix)
      rownames(all_matrix) <- colnames(mut_freq_matrix)
      colnames(all_matrix) <- rownames(mut_freq_matrix)
    } else {
      stop("Error: sample IDs are not present in either rownames or colnames of provided matrix. ")
    }
  } else {
    # Obtain sliding window mutation frequencies for all regions
    if (missing(maf_data)) {
      maf_data <- NULL
    }

    all_wide <- calc_mutation_frequency_bin_regions(
      maf_data = maf_data,
      regions_bed = regions_bed,
      these_samples_metadata = metadata,
      projection = projection,
      slide_by = slide_by,
      window_size = window_size,
      drop_unmutated = drop_unmutated,
      return_format = "wide",
      from_indexed_flatfile = from_indexed_flatfile,
      mode = mode
    )

    # Convert to a matrix with samples in colnames and bins in rownames
    all_matrix <- data.frame(t(select(all_wide, -sample_id)))
    colnames(all_matrix) <- all_wide$sample_id
    rownames(all_matrix) <- colnames(all_wide)[2:length(colnames(all_wide))]
  }


  # Normalize the expression columns
  if (!is.null(expressionColumns)) {
    these_samples_metadata <- these_samples_metadata %>%
      mutate(across(all_of(expressionColumns), ~ trim_scale_expression(.x)))
  }

  # Subset metadata to specified display columns

  meta_show <- metadata %>%
    select(sample_id, all_of(allMetaCols)) %>%
    mutate(across(all_of(allMetaCols), ~ factor(.x))) %>%
    dplyr::filter(sample_id %in% colnames(all_matrix)) %>%
    column_to_rownames(var = "sample_id")

  # Sort metadata by columns specified in sortByColumn
  if (!is.null(sortByColumns)) {
    meta_show <- arrange(
      meta_show,
      across(all_of(sortByColumns))
    )
  }

  message(paste("starting with", length(colnames(all_matrix)), "samples"))
  samples_show <- colnames(all_matrix)[which(colSums(all_matrix) >= min_mut_tumour)]
  message(paste("returning matrix with", length(samples_show), "samples"))
  if (length(samples_show) < 2) {
    stop("Insufficient samples remaining after filtering. This function requries at least two samples. ")
  }
  meta_show <- meta_show[rownames(meta_show) %in% samples_show, , drop = FALSE]
  matrix_show <- all_matrix[which(rowSums(all_matrix) > min_bin_recurrence), rownames(meta_show)]
  stopifnot(identical(rownames(meta_show), colnames(matrix_show)))

  # Set heatmap colour function
  bin_col_fun <- colorRamp2(
    c(0, 3, 6, 9),
    c(backgroundColour, "orange", "red", "purple")
  )

  # Handle custom colours
  annoColumns <- unique(c(metadataColumns, sortByColumns))
  # Get columns with no custom colours specified
  needsColour <- annoColumns[!annoColumns %in% names(customColours)]

  gamblColours <- NULL
  if (length(needsColour) > 0) {
    gamblColours <- lapply(needsColour, function(x) {
      colours <- GAMBLR.helpers::get_gambl_colours()[levels(meta_show[[x]])]
      colours <- colours[unique(names(colours))][!is.na(names(colours))]
    })
    names(gamblColours) <- needsColour
  }
  annoColoursTmp <- append(gamblColours, customColours)

  # Check that there are enough colour values for each annotation
  annoColours <- lapply(annoColumns, function(x) {
    if (length(levels(meta_show[[x]])[
      !levels(meta_show[[x]]) %in% names(annoColoursTmp[[x]])
    ] > 1)) {
      message(paste(
        "Warning: Insufficient values available for annotation ",
        x,
        "- using default colours. "
      ))
      colours <- GAMBLR.helpers::get_gambl_colours("blood")[1:length(levels(meta_show[[x]]))]
      names(colours) <- levels(meta_show[[x]])
    } else {
      return(annoColoursTmp[[x]])
    }
    return(colours)
  })
  names(annoColours) <- annoColumns

  # Add colour functions for expression columns
  if (!is.null(expressionColumns)) {
    col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
    exprColours <- lapply(expressionColumns, function(x) {
      return(col_fun)
    })
    names(exprColours) <- expressionColumns
    annoColours <- append(annoColours, exprColours)
  }

  # assign bins back to regions for better annotation
  assign_bins_to_region <- function(bin_names, rdf, label_by = "name") {
    bin_df <- data.frame(bin_name = bin_names) %>%
      separate(bin_name, into = c("chrom", "start")) %>%
      mutate(start = as.integer(start)) %>%
      mutate(end = start + 1) %>%
      mutate(bin_name = bin_names)

    regions.dt <- rdf %>%
      mutate(
        start = start - window_size,
        end = end + window_size
      ) %>%
      as.data.table()
    setkey(regions.dt, chrom, start, end)

    bin.dt <- as.data.table(bin_df)
    setkey(bin.dt, chrom, start, end)
    bin_overlapped <- foverlaps(bin.dt, regions.dt) %>%
      as.data.frame() %>%
      rename(label = !!sym(label_by)) %>%
      arrange(label, start) %>%
      select(bin_name, label) %>%
      distinct(bin_name, .keep_all = TRUE) %>%
      column_to_rownames(var = "bin_name") %>%
      drop_na()

    return(bin_overlapped)
  }

  bin_annot <- assign_bins_to_region(
    bin_names = rownames(matrix_show),
    rdf = regions_bed,
    label_by = label_regions_by
  )
  heatmap_legend_param <- list(
    title = "Mutation count",
    at = c(0, 2, 4, 6, 8, 10),
    nrow = legend_row,
    ncol = legend_col,
    legend_direction = legend_direction,
    labels_gp = gpar(fontsize = legendFontSize)
  )

  annotation_legend_param <- list(
    nrow = legend_row,
    ncol = legend_col,
    direction = legend_direction,
    labels_gp = gpar(fontsize = legendFontSize)
  )

  if (orientation == "sample_rows") {
    to_show_t <- t(matrix_show[rownames(bin_annot), ])

    row_annot <- HeatmapAnnotation(
      df = meta_show,
      show_legend = T,
      which = "row",
      col = annoColours,
      na_col = naColour,
      annotation_legend_param = annotation_legend_param
    )
    if (show_gene_colours) {
      col_annot <- HeatmapAnnotation(
        df = bin_annot,
        show_legend = F,
        which = "col",
        annotation_legend_param = annotation_legend_param
      )
    } else {
      col_annot <- HeatmapAnnotation(value = anno_empty(border = FALSE))
    }
    ht <- Heatmap(
      to_show_t[rownames(meta_show), rownames(bin_annot)],
      cluster_columns = cluster_cols_heatmap,
      cluster_rows = cluster_rows_heatmap,
      col = bin_col_fun,
      bottom_annotation = col_annot,
      left_annotation = row_annot,
      show_row_names = F,
      show_column_names = F,
      column_split = factor(bin_annot$label),
      column_title_gp = gpar(fontsize = region_fontsize),
      column_title_rot = label_regions_rotate,
      row_title_gp = gpar(fontsize = 10),
      heatmap_legend_param = heatmap_legend_param
    )
  } else {
    col_annot <- HeatmapAnnotation(
      df = meta_show,
      show_legend = T,
      which = "col",
      col = annoColours,
      na_col = naColour,
      annotation_legend_param = annotation_legend_param
    )
    if (show_gene_colours) {
      row_annot <- HeatmapAnnotation(
        df = bin_annot,
        show_legend = F,
        which = "row",
        annotation_legend_param = annotation_legend_param
      )
    } else {
      row_annot <- rowAnnotation(value = anno_empty(border = FALSE))
    }
    ht <- Heatmap(
      as.matrix(matrix_show)[rownames(bin_annot), rownames(meta_show)],
      show_heatmap_legend = F,
      cluster_columns = cluster_rows_heatmap,
      cluster_rows = cluster_cols_heatmap,
      col = bin_col_fun,
      bottom_annotation = col_annot,
      left_annotation = row_annot,
      show_row_names = F,
      show_column_names = F,
      row_split = factor(bin_annot$label),
      row_title_gp = gpar(fontsize = region_fontsize),
      row_title_rot = label_regions_rotate,
      column_title_gp = gpar(fontsize = 8),
      heatmap_legend_param = heatmap_legend_param
    )
  }

  if (return_heatmap_obj) {
    return(ht)
  } else {
    draw(ht, heatmap_legend_side = legend_side)
  }
}
