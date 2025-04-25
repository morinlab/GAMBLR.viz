#' @title PrettyOncoplot
#'
#' @description Create a highly customizable oncoplot.
#'
#' @details Generatss an oncoplot (a.k.a. oncoprint) that is pretty
#' using ComplexHeatmap. The metadata is expected to follow the structure
#' and column naming used in GAMBL.
#' If you provide your own non-GAMBL samples and metadata,
#' you must include at least the following columns with these names.
#' The first one should match the Tumor_Sample_Barcode in the MAF
#' object or onco_matrix you provide.
#' sample_id, pathology
#'
#' @param maf_df A maf as data frame containing the mutations you
#' want to plot.
#' @param gene_cnv_df An optional data frame of CN status for
#' genes you want included (rows = sample_id,
#' columns = Hugo_Symbol)
#' See [GAMBLR.results::get_cnv_and_ssm_status] for more information.
#' @param binned_cnv_df An optional data frame with the genome-wide CN
#' status of your samples in genomic bins
#' see [GAMBLR.utils::segmented_data_to_cn_matrix] for more information.
#' @param genes An optional vector of genes to restrict
#' your plot to.
#' @param include_noncoding List of non-coding regions to
#' be included, default is NULL. Specify like this:
#' include_noncoding=list("NFKBIZ" = c("3'UTR"), "HNRNPH1" = "Splice_Region")
#' @param keepGeneOrder Set to TRUE if you want to preserve
#' the gene order specified.
#' @param keepSampleOrder Set to TRUE if you want to preserve the sample order
#' specified. The default value is FALSE and respects all of the specified
#' ordering.
#' @param highlightHotspots Set to TRUE to highlight hot spots.
#' Default is FALSE.
#' @param these_samples_metadata Data frame containing metadata for your
#' samples.
#' @param genes_CN_thresh A named vector specifying the genes whose copy
#' number status should be incorporated.
#' The names must be the gene symbols and the values should be integers that
#' indicate the maximum or minimum CN states to consider for that gene.
#' For example:
#' 'REL'=4 would show CN 4 or higher; 'TP53'=1 would show
#' heterozygous and homozygous deletions
#' 'BCL2'=3 would show single-copy gains or higher
#' @param metadataColumns A vector containing the categorical
#' column names you want to plot below.
#' @param numericMetadataColumns A vector containing the numeric columns you
#' want to plot below.
#' @param expressionColumns Optional variable for retreiving expression
#' values for a specific gene(s).
#' @param numericMetadataMax A numeric vector of cutoffs to apply to
#' numeric columns above.
#' @param sortByColumns A vector containing the column names you want
#' to sort columns (patients) on.
#' @param arrange_descending A Boolean parameter. Set to TRUE to sort
#' metadata in descending fashion. Default is FALSE.
#' @param removeNonMutated Set to TRUE to drop unmutated cases.
#' @param minMutationPercent Only genes mutated in more than
#' minMutationPercent \% patients will be included.
#' @param fontSizeGene Font size for gene labels (default 6).
#' @param annoAlpha Optional alpha to apply to annotation colours.
#' @param mutAlpha Optional alpha to apply to mutation colours.
#' @param recycleOncomatrix Set to TRUE most of the time to
#' reuse the oncomatrix saved by maftools.
#' @param box_col Colour of boxes for outlining mutations (can be
#' problematic with larger oncoprints).
#' @param metadataBarHeight Optional argument to adjust the height of
#' bar with annotations. The default is 1.5.
#' @param metadataBarFontsize Optional argument to control for the
#' font size of metadata annotations. The default is 5.
#' @param hideTopBarplot Optional argument for removing top bar plot.
#' Default value is TRUE.
#' @param tally_all_mutations Optional argument. Set to TRUE to tally
#' all mutations. Default is FALSE.
#' @param tally_all_mutations_max Optional argument. Default is 1000.
#' @param hideSideBarplot Optional argument for removing side bar plot.
#' Default value is FALSE.
#' @param splitColumnName Optional argument to indicate which metadata
#' column to split on.
#' @param splitGeneGroups Split genes into groups for better seperation
#' (between different gene-groups) in prettyOncoplot.
#' @param legend_row Fiddle with these to widen or narrow your legend.
#' @param legend_col Fiddle with these to widen or narrow your legend.
#' @param showTumorSampleBarcode Optional argument for showing tumor
#' barcode. Default is FALSE.
#' @param groupNames optional vector of group names to be displayed
#' above heatmap. Should be the same length as the number of groups
#' that will be shown. Default is NULL (no labels).
#' @param verbose Set to TRUE to enable verbose mode (debugging messages.
#' @param hide_annotations Hide annotations for specifc ashms. argument
#' takes a list with annotations.
#' @param hide_annotations_tracks When hide_annotations is supplied
#' with a list of columns, this parameter can optionally also not display
#' those columns as the annotation track. Accepts TRUE and FALSE (default).
#' @param annotate_specific_genes Optional argument, specifying whether
#' the features should be labelled according to their significance in
#' one of the pathologies. Default is FALSE (no annotation).
#' @param this_forest_object If annotate_specific_genes is specified,
#' this arguments takes the output of GAMBLR::prettyForestPlot directly
#' to determine the annotations.
#' @param custom_colours Provide named vector (or named list of vectors)
#' containing custom annotation colours if you do not want to use
#' standardized pallette.
#' @param legend_direction Direction of legend, default is "horizontal".
#' @param ylim Limit for y-axis.
#' @param legend_position Position of legend, default is "bottom".
#' @param annotation_row Row for annotations, default is 2.
#' @param annotation_col Column for annotations, default is 1.
#' @param legendFontSize Font size for legend, default is 10.
#' @param cluster_rows Force clustering of genes with correlated
#' mutation patterns
#' @param cluster_cols Force clustering of patients with correlated
#' mutation patterns
#' @param clustering_distance_rows Distance metric used for clustering
#' when cluster_rows = TRUE
#' @param clustering_distance_cols Distance metric used for clustering
#' when cluster_cols = TRUE
#' @param split_rows_kmeans K value for k-means clustering on rows
#' @param split_columns_kmeans K value for k-means clustering on columns
#' @param dry_run Set to TRUE to more efficiently view the clustering
#' result while debugging cluster_rows/clustering_distance_rows
#' @param simplify_annotation Collapse/group the variant effect categories
#' to only 3 options. This is a much faster option for when many
#' patients/genes are included.
#' @param simplify_bg_colour When simplify_annotation is called, adjust
#' the color of the background by passign a value to this argument.
#' Default is NA.
#' @param gap Size of gap between columns represented as a proportion
#' of the full width of the column. Default 0 (no gap).
#' @param return_inputs Optional flag to return the plot and various
#' other internal objects such as the underlying mutation matrix.
#' @param use_raster Whether to rasterize image
#' @param show_pct TRUE by default. Set to FALSE to hide percentage.
#' @param hide_annotation_name Default: FALSE
#' @param stacked Deprecated. See [GAMBLR.viz::prettyStackedOncoplot]
#' for this functionality.
#' @param cluster_numeric_rows Deprecated.
#' See [GAMBLR.viz::prettyStackedOncoplot] for this functionality.
#' @param cluster_numeric_cols Deprecated.
#' See [GAMBLR.viz::prettyStackedOncoplot] for this functionality.
#' @param numeric_heatmap_type Deprecated.
#' See [GAMBLR.viz::prettyStackedOncoplot] for this functionality.
#' @param numeric_heatmap_location Deprecated.
#' See [GAMBLR.viz::prettyStackedOncoplot] for this functionality.
#' @param cnv_df Deprecated. See gene_cnv_df. 
#' @param summarizeByColumns Optional vector of metadata column
#' names that will be summarized on the left as stacked bar plots.
#' @param longest_label Optional character vector specifying the
#' longest gene name. This is only needed when you are vertically
#' combining multiple prettyOncoplots after you create them
#' @param right_anno_column *Experimental*. Specify one metadata
#' column name that will replace the default bar plot. The incidence
#' of mutations in patients will be shown as a grouped bar chart
#' with a bar for each different value of right_anno_column e.g. "pathology"
#' @param coloured_genes Specify any genes whose labels you
#' want to be coloured red instead of the default black.
#' 
#' @return By default, nothing unless return_inputs is specified,
#' in which case it returns a named list that contains different
#' things depending on how the function was run
#' At the very least, it will contain the Heatmap object a logical
#' matrix indicating the mutation status of each gene and
#' patient shown in the output.
#'
#' @import tidyr dplyr circlize ComplexHeatmap ggplot2 GAMBLR.helpers tibble
#' @export
#'
#' @examples
#'
#' suppressMessages(
#'   suppressWarnings({
#'     # load packages
#'     library(grid)
#'     library(dplyr)
#'     library(GAMBLR.open)
#'
#'     # Using GAMBLR.open
#'     maf_metadata <- GAMBLR.open::get_gambl_metadata(
#'       seq_type_filter = "genome"
#'     ) %>%
#'       dplyr::filter(
#'         pathology %in% c("FL", "DLBCL"),
#'         study == "FL_Dreval"
#'       )
#'
#'
#'     maf_data <- get_coding_ssm(
#'       these_samples_metadata = maf_metadata
#'     )
#'   })
#' )
#'
#' # define some genes of interest
#' fl_genes <- GAMBLR.data::lymphoma_genes %>%
#'   dplyr::filter(FL_Tier == 1) %>%
#'   pull(Gene)
#'
#' dlbcl_genes <- GAMBLR.data::lymphoma_genes %>%
#'   dplyr::filter(DLBCL_Tier == 1, !Gene %in% fl_genes) %>%
#'   pull(Gene)
#'
#' genes <- c(fl_genes, dlbcl_genes)
#'
#' # For splitting into gene sets
#' split_genes <- c(
#'   rep("FL", length(fl_genes)),
#'   rep("DLBCL", length(dlbcl_genes))
#' )
#' names(split_genes) <- genes
#'
#' suppressMessages(
#'   suppressWarnings({
#'     prettyOncoplot(
#'       maf_df = maf_data, genes = genes,
#'       these_samples_metadata = maf_metadata,
#'       splitGeneGroups = split_genes,
#'       minMutationPercent = 5
#'     )
#'   })
#' )
#'
#' # Was that too slow for you? Enable the simplify_annotation
#' # parameter for a quicker result.
#'
#' suppressMessages(
#'   suppressWarnings({
#'     prettyOncoplot(
#'       maf_df = maf_data, genes = genes,
#'       these_samples_metadata = maf_metadata,
#'       splitGeneGroups = split_genes,
#'       minMutationPercent = 5,
#'       simplify_annotation = TRUE
#'     )
#'   })
#' )
#'
#' # Want to include copy number? You have two options.
#' # Option 1:
#' # Incorporate CN status of specific genes into your oncoplot
#' # along with mutations.
#' # There are two ways to go about this.
#' # The original way involves using the helper function get_cnv_and_ssm_status
#'
#' gene_regions <- data.frame(
#'   gene_id = c(
#'     "REL", "CDKN2A",
#'     "MIR17HG", "TP53", "ATM", "FAS", "SMARCA4", "B2M", "TNFRSF14",
#'     "TMEM30A", "TNFAIP3", "BCL2"
#'   ),
#'   cn_thresh = c(4, 1, 4, 1, 1, 1, 1, 1, 1, 1, 1, 3)
#' )
#'
#' # this data frame specifies the threshold and directionality for
#' # each gene's copy number state to display on the oncoplot.
#' # Amplifications will be shown for REL and MIR17HG, gains
#' # for BCL2, deletions for the rest
#' print(gene_regions)
#'
#' suppressMessages(
#'   suppressWarnings({
#'     gene_cnv <- GAMBLR.results::get_cnv_and_ssm_status(
#'       only_cnv = "all",
#'       these_samples_metadata = get_gambl_metadata(),
#'       genes_and_cn_threshs = gene_regions
#'     )
#'
#'
#'
#'
#'     prettyOncoplot(
#'       maf_df = maf_data, genes = c(
#'         "CREBBP", "EZH2", "MYD88",
#'         "TCF3", "BCL2", "BCL7A",
#'         "MEF2B", "POU2F2", "POU2AF1",
#'         "ID3", "MYC",
#'         "RRAGC", "TCL1A", "KMT2D",
#'         "PIM1", "CD79B", "TMSB4X",
#'         "TMEM30A", "TNFAIP3"
#'       ),
#'       these_samples_metadata = maf_metadata,
#'       cluster_rows = TRUE,
#'       metadataColumns = c(
#'         "pathology",
#'         "lymphgen",
#'         "seq_type",
#'         "ffpe_or_frozen"
#'       ),
#'       cluster_cols = FALSE,
#'       simplify_annotation = TRUE,
#'       cnv_df = gene_cnv,
#'       sortByColumns = c("pathology", "lymphgen")
#'     )
#'   })
#' )
#'
#'
#' # Option 2:
#' # The second way to incorporate copy number relies
#' # instead on a binned copy number matrix
#' # If you already have one on hand, this is clearly the preferred approach!
#' # First let's make one with the help of segmented_data_to_cn_matrix
#' \dontrun{
#' all_segments <- get_cn_segments(these_samples_metadata = maf_metadata)
#' all_states_binned <- segmented_data_to_cn_matrix(
#'   seg_data = all_segments,
#'   strategy = "auto_split",
#'   n_bins_split = 1000,
#'   fill_missing_with = "avg_ploidy",
#'   adjust_for_ploidy = TRUE,
#'   these_samples_metadata = maf_metadata
#' )
#'
#'
#' # Note: adjust_for_ploidy = TRUE ensures the relative
#' # CN status is used for high-ploidy cases
#'
#' # as before, we need to specify which genes we want
#' # CN events shown for and what direction (gain or loss)
#' # This is done a bit more easily with the genes_CN_thresh option.
#'
#' CN_thresh <- c(
#'   "REL" = 4,
#'   "CDKN2A" = 1,
#'   "MIR17HG" = 4,
#'   "TP53" = 1,
#'   "TNFRSF14" = 1,
#'   "TNFAIP3" = 1
#' )
#'
#'
#' prettyOncoplot(
#'   maf_df = maf_data,
#'   binned_cnv_df = all_states_binned,
#'   genes_CN_thresh = CN_thresh,
#'   genes = head(genes, 25),
#'   these_samples_metadata = maf_metadata,
#'   cluster_rows = TRUE,
#'   metadataColumns = c(
#'     "pathology",
#'     "genetic_subgroup",
#'     "seq_type",
#'     "ffpe_or_frozen"
#'   ),
#'   cluster_cols = FALSE,
#'   simplify_annotation = TRUE,
#'   sortByColumns = c("pathology", "genetic_subgroup"),
#'   minMutationPercent = 0
#' )
#' }
prettyOncoplot <- function(maf_df, # nolint: object_name_linter.
                           gene_cnv_df,
                           binned_cnv_df,
                           genes,
                           include_noncoding = NULL,
                           keepGeneOrder = FALSE,
                           keepSampleOrder = FALSE,
                           highlightHotspots = FALSE,
                           these_samples_metadata,
                           genes_CN_thresh,
                           metadataColumns,
                           numericMetadataColumns,
                           expressionColumns = c(),
                           numericMetadataMax,
                           sortByColumns,
                           sortByGenes,
                           arrange_descending = FALSE,
                           removeNonMutated = FALSE,
                           minMutationPercent = 0,
                           mutAlpha = 1,
                           recycleOncomatrix = FALSE,
                           splitColumnName,
                           splitGeneGroups,
                           showTumorSampleBarcode = FALSE,
                           groupNames,
                           hide_annotations,
                           hide_annotations_tracks = FALSE,
                           annotate_specific_genes = FALSE,
                           this_forest_object = NULL,
                           custom_colours = NULL,
                           hideTopBarplot = TRUE,
                           tally_all_mutations = FALSE,
                           tally_all_mutations_max = 1000,
                           hideSideBarplot = FALSE,
                           box_col = NA,
                           annoAlpha = 1,
                           legend_direction = "horizontal",
                           ylim = NULL,
                           legend_position = "bottom",
                           legend_row = 3,
                           legend_col = 3,
                           metadataBarHeight = 1.5,
                           metadataBarFontsize = 5,
                           metadataSide = "bottom",
                           legendFontSize = 10,
                           fontSizeGene = 6,
                           annotation_row = 2,
                           annotation_col = 1,
                           verbose = FALSE,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           clustering_distance_rows = "binary",
                           clustering_distance_cols = "binary",
                           split_rows_kmeans,
                           split_columns_kmeans,
                           dry_run = FALSE,
                           simplify_annotation = FALSE,
                           simplify_bg_colour = "transparent",
                           stacked = FALSE,
                           return_inputs = FALSE,
                           gap = 0,
                           use_raster = NULL,
                           plot_width = NULL,
                           plot_height = NULL,
                           show_any_legend = TRUE,
                           pct_side = "left",
                           pctFontSize = 6, # nolint: object_name_linter.
                           row_names_side = "right",
                           show_pct = TRUE,
                           hide_annotation_name = FALSE,
                           cnv_df,
                           summarizeByColumns,
                           longest_label,
                           right_anno_column,
                           hide_left_annotations = FALSE,
                           coloured_genes,
                           annotation_name_side = "bottom",
                           axis_font_size = 5,
                           right_anno_width = 2) {
  all_args <- mget(names(formals()),             # the formal names
                   sys.frame(sys.nframe()),      # the current call frame
                   inherits = FALSE)             # don’t look up the stack
  if (stacked) {
    stop("stacked functionality has migrated to prettyStackedOncoplot!")
  }
  if (!missing(cnv_df)) {
    warning("cnv_df is replaced by gene_cnv_df. Please update your code.")
    gene_cnv_df <- cnv_df
  }
  if ("maf_data" %in% class(maf_df)) {
    maf_df <- strip_genomic_classes(maf_df)
  }
  if (any(these_samples_metadata$seq_type %in% c("mrna", "promethION"))) {
    message("metadata contains rows with incompatible seq_type")
    message("Dropping these for you.")
    these_samples_metadata <- dplyr::filter(
      these_samples_metadata,
      !seq_type %in% c("mrna", "promethION")
    )
  }
  if (!missing(split_columns_kmeans) && !missing(splitColumnName)) {
    stop("split_columns_kmeans and splitColumnName are incompatible. Use one or the other")
  }
  if (missing(maf_df)) {
    stop("You must provide maf data frame.")
  }
  if (!missing(genes_CN_thresh)) {
    if (missing(binned_cnv_df)) {
      error <- paste(
        "genes_CN_thresh provided without binned_cnv_df.",
        "Both arguments must be provided together."
      )
      stop(error)
    } else {
      cn_thresh <- data.frame(
        gene_id = names(genes_CN_thresh),
        cn_thresh = genes_CN_thresh
      )
      # this will not work for users of GAMBLR.open until this function is added
      cnv_df <- get_cnv_and_ssm_status(
        genes_and_cn_threshs = cn_thresh,
        cn_matrix = binned_cnv_df,
        these_samples_metadata = these_samples_metadata,
        maf_df = maf_df,
        only_cnv = "all"
      )
    }
  }
  onco_matrix_coding <- GAMBLR.helpers::coding_class[
    !GAMBLR.helpers::coding_class %in% c(
      "Silent",
      "Splice_Region",
      "Targeted_Region"
    )
  ]
  if (length(include_noncoding) > 0) {
    if (verbose) {
      print("Generating maf_df_noncoding")
    }
    nc <- data.frame(t(data.frame(include_noncoding))) %>%
      rownames_to_column("Hugo_Symbol") %>%
      pivot_longer(-Hugo_Symbol) %>%
      dplyr::select(
        Hugo_Symbol,
        Variant_Classification = value
      )
    maf_df_noncoding <- inner_join(maf_df, nc)
    if (verbose) {
      print(nc)
      print(maf_df_noncoding)
    }
  } else {
    maf_df_noncoding <- maf_df %>%
      dplyr::filter(!Tumor_Sample_Barcode %in% maf_df$Tumor_Sample_Barcode)
  }
  # Ensure Tumor_Sample_Barcode is always present because different
  # versions of metadata may generate different output with inconsistent
  # columns present in the resulting metadata
  if (!"Tumor_Sample_Barcode" %in% colnames(these_samples_metadata)) {
    these_samples_metadata <- these_samples_metadata %>%
      dplyr::mutate(Tumor_Sample_Barcode = sample_id)
  }
  patients <- dplyr::pull(these_samples_metadata, sample_id) %>% unique()
  if (verbose) {
    print(
      paste(
        "There are",
        length(patients),
        "unique patients in the provided metadata."
      )
    )
  }
  # Modify maf_df to include only coding mutations then
  # add specified non-coding mutations
  maf_df <- maf_df %>%
    dplyr::filter(
      Tumor_Sample_Barcode %in% patients,
      Variant_Classification %in% onco_matrix_coding
    ) %>%
    bind_rows(maf_df_noncoding)
  # Make sure that the N of patients matches between metadata and maf
  # so the displayed %ages and counts are correct
  # This will account for cases where the sample is in metadata but has 0
  # mutations in maf
  maf_patients <- unique(as.character(maf_df$Tumor_Sample_Barcode))
  if (any(!maf_patients %in% patients)) {
    extra <- maf_patients[which(!maf_patients %in% patients)]
    patients <- maf_patients[which(maf_patients %in% patients)]
    n_drop <- length(extra)
    message(paste(
      n_drop,
      "patients are not in your metadata,",
      "will drop them from the data before displaying"
    ))
    maf_df <- filter(maf_df, Tumor_Sample_Barcode %in% patients)
    maf_patients <- unique(as.character(maf_df$Tumor_Sample_Barcode))
  }
  # Now do opposite check and make sure that all patients from metadata are
  # also present in maf
  if (any(!patients %in% maf_patients)) {
    maf_df <- supplement_maf(
      incoming_maf = maf_df,
      these_samples_metadata = these_samples_metadata
    )
    maf_df <- maf_df %>%
      mutate(
        Hugo_Symbol = ifelse(
          Hugo_Symbol == "GARBAGE",
          paste0(Hugo_Symbol, row_number()),
          Hugo_Symbol
        )
      )
  }
  ##### Handling of genes to be displayed

  # When no genes provided, display only genes above the specified min
  # mutation frequency threshold. Here we will define those genes
  if (missing(genes)) {
    if (verbose) {
      print("finding genes to include...")
    }
    gene_summary <- maf_df %>%
      distinct(
        Tumor_Sample_Barcode,
        Hugo_Symbol, Variant_Classification, Start_Position,
        End_Position
      ) %>%
      filter(
        Tumor_Sample_Barcode %in% patients,
        Variant_Classification %in% onco_matrix_coding
      ) %>%
      distinct(Tumor_Sample_Barcode, Hugo_Symbol) %>%
      group_by(Hugo_Symbol) %>%
      summarize(MutatedSamples = n(), .groups = "drop") %>%
      arrange(desc(MutatedSamples))
    genes <- gene_summary
    if (verbose) {
      print(paste("numgenes:", length(genes[, 1])))
    }
    colnames(genes)[2] <- "mutload"
    genes$fractMutated <- genes$mutload / length(maf_patients)
    genes <- genes %>%
      dplyr::filter(fractMutated * 100 >= minMutationPercent) %>%
      dplyr::pull(Hugo_Symbol)
    lg <- length(genes)
    if (!recycleOncomatrix) {
      message(paste("creating oncomatrix with", lg, "genes"))
      mat_origin <- GAMBLR.helpers::create_onco_matrix(maf_df, genes)
      mat_origin <- mat_origin[, !colSums(mat_origin == "") ==
        nrow(mat_origin), drop = FALSE]
      tsbs <- maf_df %>%
        dplyr::distinct(
          Tumor_Sample_Barcode,
          Hugo_Symbol,
          Variant_Classification,
          Start_Position,
          End_Position
        ) %>%
        dplyr::filter(
          Tumor_Sample_Barcode %in% patients,
          Variant_Classification %in% onco_matrix_coding
        ) %>%
        dplyr::group_by(Tumor_Sample_Barcode) %>%
        dplyr::summarize(total = n(), .groups = "drop") %>%
        dplyr::arrange(desc(total)) %>%
        dplyr::pull(Tumor_Sample_Barcode)
      if (!removeNonMutated) {
        tsb.include <- matrix( # nolint: object_name_linter.
          data = 0, nrow = nrow(mat_origin),
          ncol = length(tsbs[!tsbs %in% colnames(mat_origin)])
        )
        colnames(tsb.include) <- tsbs[!tsbs %in% colnames(mat_origin)] # nolint: object_name_linter.
        rownames(tsb.include) <- rownames(mat_origin) # nolint: object_name_linter.
        mat_origin <- cbind(mat_origin, tsb.include)
      }
      if (verbose) {
        print(paste("numcases:", length(tsbs)))
      }
    }
  } else {
    if (any(duplicated(genes))) {
      stop(paste(
        "There are duplicated elements in the provided gene list.",
        "Please ensure only unique entries are present in this list."
      ))
    }
    if (verbose) {
      print("summarizing genes...")
    }
    gene_summary <- maf_df %>%
      distinct(
        Tumor_Sample_Barcode,
        Hugo_Symbol, Variant_Classification, Start_Position,
        End_Position
      ) %>%
      filter(
        Hugo_Symbol %in% genes,
        Tumor_Sample_Barcode %in% patients, Variant_Classification %in%
          onco_matrix_coding
      ) %>%
      distinct(
        Tumor_Sample_Barcode,
        Hugo_Symbol
      ) %>%
      group_by(Hugo_Symbol) %>%
      summarize(
        MutatedSamples = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(MutatedSamples))
    gene_summary$fractMutated <- gene_summary$MutatedSamples /
      length(maf_patients)
    gene_summary <- gene_summary %>% filter(fractMutated *
      100 >= minMutationPercent)
    if (verbose) {
      print(gene_summary)
      print(genes)
    }
    if (!recycleOncomatrix) {
      if (verbose) {
        print("creating oncomatrix ")
      }
      npat <- length(patients)
      if (verbose) {
        print(paste("patients:", npat))
      }
      mat_origin <- GAMBLR.helpers::create_onco_matrix(maf_df, genes)
      mat_origin <- mat_origin[, !colSums(mat_origin == "") ==
        nrow(mat_origin), drop = FALSE]
      tsbs <- maf_df %>%
        distinct(
          Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification,
          Start_Position, End_Position,
          .keep_all = TRUE
        )

      tsbs <- tsbs %>%
        filter(
          Tumor_Sample_Barcode %in% patients,
          Variant_Classification %in% onco_matrix_coding
        ) %>%
        ungroup() %>%
        group_by(Tumor_Sample_Barcode) %>%
        summarize(total = n())

      tsbs <- tsbs %>%
        arrange(desc(total)) %>%
        pull(Tumor_Sample_Barcode)
      if (verbose) {
        print(paste("numcases:", length(tsbs)))
        print(paste("numgenes:", length(mat_origin[, 1])))
      }
    }
    if (!removeNonMutated && !recycleOncomatrix) {
      tsb.include <- matrix( # nolint: object_name_linter.
        data = 0, nrow = nrow(mat_origin),
        ncol = length(tsbs[!tsbs %in% colnames(mat_origin)])
      )
      colnames(tsb.include) <- tsbs[!tsbs %in% colnames(mat_origin)] # nolint: object_name_linter.
      rownames(tsb.include) <- rownames(mat_origin) # nolint: object_name_linter.
      mat_origin <- cbind(mat_origin, tsb.include)
      if (verbose) {
        print(head(tsbs))
      }
    } else if (length(include_noncoding) > 0) {
      if (verbose) {
        print("Will include noncoding mutations, remove unmutated patients")
      }
      these_have_noncoding <- maf_df %>%
        filter(
          Tumor_Sample_Barcode %in%
            patients, Hugo_Symbol %in% names(include_noncoding),
          Variant_Classification %in% unlist(unname(include_noncoding))
        ) %>%
        distinct(
          Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification,
          Start_Position, End_Position
        ) %>%
        pull(Tumor_Sample_Barcode)
      tsb.include <- matrix( # nolint: object_name_linter.
        data = 0, nrow = nrow(mat_origin),
        ncol = length(these_have_noncoding[!these_have_noncoding %in%
          colnames(mat_origin)])
      )
      colnames(tsb.include) <- these_have_noncoding[!these_have_noncoding %in% # nolint: object_name_linter.
        colnames(mat_origin)]
      rownames(tsb.include) <- rownames(mat_origin) # nolint: object_name_linter.
      mat_origin <- cbind(mat_origin, tsb.include)
    }
  }
  mat <- mat_origin
  mat[mat == 0] <- ""
  if (!missing(gene_cnv_df)) {
    if (verbose) {
      print("BEFORE:")
      print(dim(cnv_df))
    }
    cnv_df <- gene_cnv_df[rownames(gene_cnv_df) %in% patients, ,
      drop = FALSE
    ]
    if (verbose) {
      print("AFTER")
      print(dim(cnv_df))
      print(colSums(cnv_df))
    }
  }
  ##### The part below will be handling the metadata and colors
  if (missing(metadataColumns)) {
    message("Defaulting to show pathology")
    metadataColumns <- c("pathology") # nolint: object_name_linter.
  }
  if (!missing(numericMetadataColumns)) {
    message("Regular and numeric metadata will be plotted as numeric")
    metadataColumns <- metadataColumns[!metadataColumns %in%
      numericMetadataColumns]
  }
  if (missing(genes)) {
    genes <- rownames(mat)
  }
  col <- GAMBLR.helpers::get_gambl_colours("mutation", alpha = mutAlpha)
  patients_dropped <- patients[which(!patients %in% colnames(mat))]
  if (verbose) {
    message("====DROPPED=====")
    message(patients_dropped)
  }
  patients_kept <- patients[which(patients %in% colnames(mat))]
  genes_kept <- genes[which(genes %in% rownames(mat))]
  if (!missing(cnv_df)) {
    # ensure any genes in this input are added to genes_kept
    cnv_genes <- colnames(cnv_df)
    cnv_patients <- rownames(cnv_df)[which(rownames(cnv_df) %in%
      these_samples_metadata$sample_id)]
    cnv_df <- cnv_df[cnv_patients, , drop = FALSE]
    extra_cnv_genes <- cnv_genes[!cnv_genes %in% genes_kept]
    extra_cnv_patients <- cnv_patients[!cnv_patients %in%
      patients_kept]
  }
  if (recycleOncomatrix) {
    gene_summary_list <- apply(mat, 1, function(x) {
      sum(!is.na(x))
    })
    gene_summary <- data.frame(
      Hugo_Symbol = names(gene_summary_list),
      MutatedSamples = as.numeric(unname(gene_summary_list))
    )
    gene_summary <- arrange(gene_summary, desc(MutatedSamples))
  }

  if (!missing(minMutationPercent) && minMutationPercent > 0) {
    if (recycleOncomatrix) {
      stop(paste(
        "mintMutationPercent not available with custom oncomatrix.",
        "Feel free to implement this if you need it"
      ))
    }
    mutation_counts <- gene_summary %>% select(
      Hugo_Symbol,
      MutatedSamples
    )
    numpat <- length(patients_kept)
    mutation_counts <- mutate(mutation_counts, percent_mutated = 100 *
      MutatedSamples / numpat)
    genes_keep <- mutation_counts %>%
      dplyr::filter(percent_mutated >= minMutationPercent) %>%
      pull(Hugo_Symbol)
    genes_kept <- genes[genes %in% genes_keep]
  }

  mat <- mat[, patients_kept, drop = FALSE]
  mat <- mat[which(rownames(mat) %in% genes_kept), , drop = FALSE]

  if (!missing(cnv_df)) {
    genes_kept <- c(extra_cnv_genes, genes_kept)
    patients_kept <- c(patients_kept, extra_cnv_patients)
    # add missing genes from CNV data
    if (verbose) {
      print("KEPT GENES")
      print(genes_kept)
      print("CNV MAT:")
      print(dim(mat))
    }
    mat1 <- matrix(nrow = length(extra_cnv_genes), ncol = ncol(mat))
    colnames(mat1) <- colnames(mat)
    rownames(mat1) <- extra_cnv_genes
    mat1[] <- 0
    mat <- rbind(mat1, mat)
    for (pat in extra_cnv_patients) {
      mat[, pat] <- 0
    }
    if (verbose) {
      print("MAT:")
      print(dim(mat))
    }
  }
  spacing <- 0
  height_scaling <- 1
  if (simplify_annotation) {
    # make oncomatrix individually from the MAF
    summarize_mutation_by_class <- function(mutation_set) {
      if ("hot_spot" %in% mutation_set) {
        snv_maf <- filter(maf_df, hot_spot == TRUE) %>%
          filter(Hugo_Symbol %in% rownames(mat)) %>%
          select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
          unique() %>%
          mutate(Mutated = TRUE)
      } else {
        snv_maf <- filter(
          maf_df,
          Variant_Classification %in% mutation_set
        ) %>%
          filter(Hugo_Symbol %in% rownames(mat)) %>%
          select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
          unique() %>%
          mutate(Mutated = TRUE)
      }
      if (nrow(snv_maf) == 0) {
        stop(paste(
          "There are no genes left.",
          "Try reducing minMutationPercent"
        ))
      }
      snv_wide <- pivot_wider(snv_maf,
        names_from = "Tumor_Sample_Barcode",
        values_from = "Mutated", values_fill = FALSE
      ) %>%
        column_to_rownames("Hugo_Symbol")
      missing_g <- rownames(mat)[which(!rownames(mat) %in%
        rownames(snv_wide))]
      missing_s <- colnames(mat)[which(!colnames(mat) %in%
        colnames(snv_wide))]
      # add missing rows
      for (g in missing_g) {
        snv_wide[g, ] <- FALSE
      }
      # add missing cols
      for (s in missing_s) {
        snv_wide[, s] <- FALSE
      }
      if (verbose) {
        print(nrow(snv_wide))
        print(ncol(snv_wide))
      }
      return(snv_wide[rownames(mat), colnames(mat)])
    }
    col["Missense"] <- col["Missense_Mutation"]
    col["Truncating"] <- col["Nonsense_Mutation"]
    col["CNV"] <- "purple"
    col["HotSpot"] <- "magenta"
    # heights <- c(
    #  "Missense", "Truncating", "Splice_Site",
    #  "HotSpot", "CNV"
    # )
    path_cols = get_gambl_colours("pathology")
    dilution = 40
    width_adj <- 1 - gap
    alter_fun <- list(
      background = function(x, y, w, h) {
        grid.rect(x,
          y, w, h,
          gp = gpar(fill = NA, col = NA)
        )
      },
      DLBCL = function(x, y, w, h) {
        grid.rect(x,
          y, w, h,
          gp = gpar(fill = paste0(path_cols["DLBCL"],dilution), col = NA)
        )
      },
      BL = function(x, y, w, h) {
        grid.rect(x,
          y, w, h,
          gp = gpar(fill = paste0(path_cols["BL"],dilution), col = NA)
        )
      },
      FL = function(x, y, w, h) {
        grid.rect(x,
          y, w, h,
          gp = gpar(fill = paste0(path_cols["FL"],dilution), col = NA)
        )
      },
      Missense = function(x, y, w, h) {
        grid.rect(x, y, w * width_adj, h * 0.9, gp = gpar(
          fill = col["Missense"],
          col = NA
        ))
      }, Truncating = function(x, y, w, h) {
        grid.rect(x, y, w * width_adj, h * 0.9, gp = gpar(
          fill = col["Truncating"],
          col = NA
        ))
      }, Splice_Site = function(x, y, w, h) {
        grid.rect(x, y, w * width_adj, h * 0.9, gp = gpar(
          fill = col["Splice_Site"],
          col = NA
        ))
      }, CNV = function(x, y, w, h) {
        grid.rect(x, y, w * width_adj, h * 0.6, gp = gpar(
          fill = col["CNV"],
          col = NA
        ))
      }, HotSpot = function(x, y, w, h) {
        grid.rect(x, y, w * width_adj, h * 0.9, gp = gpar(
          fill = col["Missense"],
          col = NA
        ))
        grid.rect(x, y, w * width_adj, h * 0.5, gp = gpar(
          fill = col["HotSpot"],
          col = NA
        ))
      }, col = col
    )
    snv_df <- summarize_mutation_by_class(mutation_set = c(
      "Missense_Mutation",
      "In_Frame_Del", "In_Frame_Ins", "Translation_Start_Site"
    ))



    if (verbose) {
      print("SNV")
      print(dim(snv_df))
    }
    trunc_df <- summarize_mutation_by_class(mutation_set = c(
      "Nonsense_Mutation",
      "Frame_Shift_Del", "Frame_Shift_Ins", "Nonstop_Mutation"
    ))
    splice_df <- summarize_mutation_by_class(mutation_set = "Splice_Site")
    if (highlightHotspots) {
      if (verbose) {
        print("summarize hot spots")
      }
      hotspot_df <- summarize_mutation_by_class(mutation_set = "hot_spot")
      if (verbose) {
        print(dim(hotspot_df))
      }
      snv_df[trunc_df == TRUE | splice_df == TRUE] <- FALSE
      snv_df[hotspot_df == TRUE] <- FALSE
      trunc_df[hotspot_df == TRUE] <- FALSE
      splice_df[hotspot_df == TRUE] <- FALSE
    } else {
      snv_df[trunc_df == TRUE | splice_df == TRUE] <- FALSE
      splice_df[trunc_df == TRUE] <- FALSE
      snv_df[trunc_df == TRUE | splice_df == TRUE] <- FALSE
      splice_df[trunc_df == TRUE] <- FALSE
    }
  
    if(simplify_bg_colour == "pathology"){
      #background colour is based on metadata
      bg_df = matrix(nrow=nrow(snv_df),ncol=ncol(snv_df))
      rownames(bg_df) = rownames(snv_df)
      colnames(bg_df) = colnames(snv_df)
      for(path in unique(these_samples_metadata$pathology)){
        ids = filter(these_samples_metadata,pathology==path) %>% pull(sample_id)
        ids_in = ids[ids %in% colnames(bg_df)]
        bg_df[,ids_in] = path
      }
      #print(head(bg_df[c("EZH2","KMT2D","CREBBP","MYD88"),c(1:10)]))
      #print(table(bg_df["EZH2",]))
      #stop()
    }
    if (!missing(cnv_df)) {
      transposed_cnv_df <- t(cnv_df)
      cn_df <- matrix(
        nrow = nrow(transposed_cnv_df),
        ncol = ncol(transposed_cnv_df)
      )
      rownames(cn_df) <- colnames(cnv_df)
      colnames(cn_df) <- rownames(cnv_df)
      cn_df[] <- FALSE
      cn_df[transposed_cnv_df == 1] <- TRUE
      missing_snv_genes <- rownames(mat)[which(!rownames(mat) %in%
        rownames(cn_df))]
      if (verbose) {
        print("MISSING CNVs:")
        print(missing_snv_genes)
        print("====")
      }
      mat1 <- matrix(nrow = length(missing_snv_genes), ncol = ncol(cn_df))
      colnames(mat1) <- colnames(cn_df)
      rownames(mat1) <- missing_snv_genes
      mat1[] <- FALSE
      cn_df <- rbind(cn_df, mat1)
      cn_df <- as.data.frame(cn_df)
      missing_snv_patients <- colnames(mat)[which(!colnames(mat) %in%
        colnames(cn_df))]
      for (pat in missing_snv_patients) {
        cn_df[missing_snv_genes, pat] <- FALSE
      }
      cn_df <- cn_df[rownames(mat), patients_kept]
      # check for NAs
      na_cols <- names(which(is.na(colSums(cn_df))))
      if (verbose) {
        print("NA COLS:")
        print(na_cols)
      }
      cn_df[, na_cols] <- FALSE
    }
    if (!missing(cnv_df)) {
      # Disabled. Unsure why these were ever here.
      # cn_df[splice_df == TRUE] = FALSE
      # cn_df[trunc_df == TRUE] = FALSE
      if (highlightHotspots) {
        cn_df[hotspot_df == TRUE] <- FALSE
      }
      # cn_df[snv_df==TRUE] = FALSE
      if (verbose) {
        print("ROWNAMES cn_df:")
        print(rownames(cn_df))
        print(table(unlist(cn_df[1, ])))
      }
    }
  } else { # end simplifyAnnotation
    alter_fun <- list(background = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(fill = "#e6e6e6", col = box_col)
      )
    }, RNA = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(fill = "#F2ED36", col = box_col)
      )
    }, `3'UTR` = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(fill = "#F2ED36", col = box_col)
      )
    }, `5'UTR` = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(fill = col["5'UTR"], col = box_col)
      )
    }, Intron = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), 0.75 * h *
        height_scaling, gp = gpar(
        fill = col["Intron"],
        col = box_col
      ))
    }, Nonsense_Mutation = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(fill = col["Nonsense_Mutation"], col = box_col)
      )
    }, Splice_Site = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(fill = col["Splice_Site"], col = box_col)
      )
    }, Splice_Region = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(fill = col["Splice_Region"], col = box_col)
      )
    }, Nonstop_Mutation = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(fill = col["Nonstop_Mutation"], col = box_col)
      )
    }, Translation_Start_Site = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(
          fill = col["Translation_Start_Site"],
          col = box_col
        )
      )
    }, In_Frame_Ins = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(fill = col["In_Frame_Ins"], col = box_col)
      )
    }, In_Frame_Del = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(fill = col["In_Frame_Del"], col = box_col)
      )
    }, Frame_Shift_Del = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(fill = col["Frame_Shift_Del"], col = box_col)
      )
    }, Frame_Shift_Ins = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(fill = col["Frame_Shift_Ins"], col = box_col)
      )
    }, Multi_Hit = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(fill = col["Multi_Hit"], col = box_col)
      )
    }, Missense_Mutation = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(fill = col["Missense_Mutation"], col = box_col)
      )
    }, hot_spot = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), (height_scaling / 5) *
        h, gp = gpar(fill = "magenta", col = box_col))
    }, Silent = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
        gp = gpar(fill = col["Silent"], col = box_col)
      )
    })
  }
  if (highlightHotspots) {
    # Annotate hotspots if necessary
    if (!"hot_spot" %in% colnames(maf_df)) {
      stop(paste(
        "maf_df requires a hot_spot column.",
        "Did you forget to run annotate_hotspots?"
      ))
    }
    if (verbose) {
      print("annotating hot spots")
    }
    hot_samples <- dplyr::filter(maf_df, hot_spot == TRUE &
      Hugo_Symbol %in% genes_kept) %>%
      dplyr::select(
        Hugo_Symbol,
        Tumor_Sample_Barcode
      ) %>%
      mutate(mutated = "hot_spot") %>%
      unique()

    all_genes_df <- data.frame(Hugo_Symbol = rownames(mat))
    all_samples_df <- data.frame(Tumor_Sample_Barcode = colnames(mat))
    hs <- left_join(all_samples_df, hot_samples)
    hot_mat <- hs %>%
      pivot_wider(
        names_from = "Tumor_Sample_Barcode",
        values_from = "mutated"
      ) %>%
      left_join(all_genes_df, .) %>%
      column_to_rownames("Hugo_Symbol") %>%
      as.matrix()

    for (i in colnames(mat)) {
      mat[genes_kept, i][!is.na(hot_mat[genes_kept, i])] <-
        paste0(
          mat[genes_kept, i][!is.na(hot_mat[genes_kept, i])],
          ";",
          hot_mat[genes_kept, i][!is.na(hot_mat[genes_kept, i])]
        )
    }

    if (verbose) {
      print("done annotating hot spots")
    }
  }
  if (verbose) {
    print(colours)
    print("KEPT:")
    print(length(patients_kept))
  }
  if (length(expressionColumns)) {
    if (missing(numericMetadataColumns)) {
      numericMetadataColumns <- expressionColumns # nolint: object_name_linter.
    } else {
      numericMetadataColumns <- c( # nolint: object_name_linter.
        numericMetadataColumns,
        expressionColumns
      )
    }
  }
  if (!missing(numericMetadataColumns)) {
    metadata_df <- dplyr::filter(these_samples_metadata, sample_id %in%
      patients_kept) %>%
      column_to_rownames("sample_id") %>%
      dplyr::select(all_of(c(
        metadataColumns, numericMetadataColumns,
        expressionColumns
      )))
    # Warn the user and deal with NA values, if present
    na_columns <- metadata_df %>%
      summarise(across(all_of(numericMetadataColumns), ~ any(is.na(.)))) %>%
      pivot_longer(everything(),
        names_to = "column",
        values_to = "has_na"
      ) %>%
      filter(has_na) %>%
      pull(column)

    if (length(na_columns) > 0) {
      warning(paste(
        "The following numeric columns contain NA values: ",
        paste(na_columns, collapse = ", "),
        "\nThese will be replaced with 0\n",
        "If NAs are unexpected, you should investigate"
      ))
    }
    metadata_df <- metadata_df %>%
      mutate(across(
        na_columns,
        ~ ifelse(is.na(.x), 0, .x)
      ))
    if (!missing(numericMetadataMax)) {
      max_list <- setNames(numericMetadataMax, numericMetadataColumns)
      metadata_df <- metadata_df %>%
        mutate(across(
          names(max_list),
          ~ ifelse(.x > max_list[[cur_column()]], max_list[[cur_column()]],
            .x
          )
        ))
    }
  } else {
    metadata_df <- dplyr::filter(these_samples_metadata, sample_id %in%
      patients_kept) %>%
      column_to_rownames("sample_id") %>%
      dplyr::select(all_of(c(metadataColumns, expressionColumns)))
  }
  metadata_df <- metadata_df %>% mutate(across(
    all_of(expressionColumns),
    ~ trim_scale_expression(.x)
  ))
  # check for missing colours
  colours <- map_metadata_to_colours(
    metadataColumns = metadataColumns,
    these_samples_metadata = these_samples_metadata, annoAlpha = annoAlpha,
    verbose = verbose
  )
  if (highlightHotspots) {
    colours[["hot_spot"]] <- c(hot_spot = "magenta")
    colours[["HotSpot"]] <- "magenta"
  }
  if (!is.null(custom_colours)) {
    for (colname in names(custom_colours)) {
      colours[[colname]] <- custom_colours[[colname]]
      if (verbose) {
        print("adding:")
        print(colname)
        print(colours[[colname]])
      }
    }
    if (verbose) {
      print(names(colours))
    }
  }
  col_fun <- circlize::colorRamp2(c(1, 2, 5), c(
    "blue", "white",
    "red"
  ))
  for (exp in expressionColumns) {
    colours[[exp]] <- col_fun
  }
  for (annotation in intersect(colnames(metadata_df), names(colours))) {
    if (annotation %in% c("HotSpot", "hot_spot")) {
      next
    }
    if (!missing(numericMetadataColumns)) {
      if (annotation %in% numericMetadataColumns) {
        next
      }
    }
    all_names <- filter(metadata_df, !is.na(annotation)) %>%
      pull(annotation) %>%
      unique() %>%
      as.character()
    all_names <- all_names[!is.na(all_names)]
    if (any(!all_names %in% names(colours[[annotation]]))) {
      print("-=-=-=")
      print(all_names[!all_names %in% names(colours[[annotation]])])
      print("-=-=-=")
      print(unique(all_names))
      print("-=-=-=")
      print(names(colours[[annotation]]))
      print("-=-=-=")
      stop(paste(
        "No colour assigned for one or more values in annotation:",
        annotation,
        "Remove the offending rows or provide custom colours",
        annotation
      ))
    }
  }
  if (verbose) {
    print("COLOURS:")
    print(colours)
  }
  if (!missing(sortByColumns)) {
    if (arrange_descending) {
      metadata_df <- arrange(
        metadata_df,
        across(
          sortByColumns,
          desc
        )
      )
    } else {
      metadata_df <- arrange(metadata_df, across(sortByColumns))
    }
    if(verbose){
      print("done sorting samples")
    }
    patients_kept <- rownames(metadata_df)
  }
  if (verbose) {
    print(genes_kept)
  }
  if (missing(splitGeneGroups)) {
    row_split <- rep("", length(genes_kept))
  } else {
    if (is.factor(splitGeneGroups)) {
      row_split <- splitGeneGroups[genes_kept]
    } else {
      (row_split <- factor(splitGeneGroups[genes_kept],
        levels = unique(splitGeneGroups[genes_kept])
      ))
    }
  }
  if (!missing(groupNames)) {
    column_title <- groupNames
  } else {
    column_title <- NULL
  }
  if (keepGeneOrder) {
    gene_order <- intersect(genes, genes_kept)
  } else {
    gene_order <- NULL
  }
  ##### Now handle the annotations around the oncoplot
  if (missing(hide_annotations)) {
    show_legend <- rep(TRUE, length(colnames(metadata_df)))
  } else {
    show_legend <- rep(TRUE, length(colnames(metadata_df)))
    names(show_legend) <- colnames(metadata_df)
    show_legend[hide_annotations] <- FALSE
  }
  if (hideTopBarplot) {
    top_annotation <- NULL
  } else {
    tally_mutations <- maf_df %>%
      dplyr::filter(Tumor_Sample_Barcode %in% patients_kept) %>%
      group_by(Tumor_Sample_Barcode) %>%
      summarize(n_mutations = n()) %>%
      ungroup() %>%
      arrange(match(
        Tumor_Sample_Barcode,
        patients_kept
      )) %>%
      select(n_mutations) %>%
      mutate(n_mutations = ifelse(n_mutations >
        tally_all_mutations_max, tally_all_mutations_max,
      n_mutations
      ))
    if (is.null(ylim) && !tally_all_mutations) {
      top_annotation <- HeatmapAnnotation(cbar = anno_oncoprint_barplot())
    } else if (!is.null(ylim) && !tally_all_mutations) {
      top_annotation <- HeatmapAnnotation(cbar = anno_oncoprint_barplot(ylim = ylim))
    } else if (is.null(ylim) && tally_all_mutations) {
      top_annotation <- columnAnnotation(` ` = anno_barplot(tally_mutations))
    } else if (!is.null(ylim) && tally_all_mutations) {
      top_annotation <- columnAnnotation(` ` = anno_barplot(tally_mutations,
        ylim = ylim
      ))
    }
  }
  # Handle right annotation for specific genes
  if (annotate_specific_genes && is.null(this_forest_object)) {
    message(paste(
      "WARNING: You requested right annotation",
      "but neglected to provide output of GAMBLR::prettyForestPlot",
      "No right annotation will be drawn."
    ))
    right_annotation <- NULL
  } else if (annotate_specific_genes) {
    these_comparisons <- this_forest_object$mutmat$comparison %>%
      levels()
    enrichment_label <- mat[
      intersect(genes, genes_kept),
      patients_kept
    ] %>%
      rownames_to_column("gene") %>%
      select(gene) %>%
      left_join(this_forest_object$fisher %>%
        select(gene, estimate, q.value)) %>%
      mutate(`Enriched in` = case_when(
        estimate ==
          "Inf" & q.value <= 0.1 ~ these_comparisons[1], estimate ==
          "-Inf" & q.value <= 0.1 ~ these_comparisons[2], is.na(estimate) ~
          "NA", estimate <= 1 & q.value <= 0.1 ~ these_comparisons[2],
        estimate > 1 & q.value <= 0.1 ~ these_comparisons[1],
        TRUE ~ "Neither"
      )) %>%
      pull("Enriched in")
    right_annotation <- rowAnnotation(
      ` ` = enrichment_label,
      col = list(` ` = c(GAMBLR.helpers::get_gambl_colours()[these_comparisons],
        Neither = "#ACADAF", `NA` = "#000000"
      )), simple_anno_size = unit(
        metadataBarHeight,
        "mm"
      ), annotation_legend_param = list(
        title = "Enriched in",
        nrow = legend_row,
        ncol = legend_col,
        direction = legend_direction,
        labels_gp = gpar(fontsize = legendFontSize)
      )
    )
  } else {
    if (hideSideBarplot) {
      right_annotation <- NULL
    } else {
      right_annotation <- rowAnnotation(rbar = anno_oncoprint_barplot())
    }
  }
  if (keepSampleOrder) {
    patients_kept <- patients_kept[order(match(
      patients_kept,
      these_samples_metadata %>% filter(sample_id %in%
        patients_kept) %>%
        pull(sample_id)
    ))]
    metadata_df <- metadata_df[order(match(
      rownames(metadata_df),
      patients_kept
    )), , drop = FALSE]
  }
  if (simplify_annotation) {
    if (!missing(cnv_df)) {
      cn_df <- cn_df[genes_kept, patients_kept]
    }
    trunc_df <- trunc_df[genes_kept, patients_kept]
    snv_df <- snv_df[genes_kept, patients_kept]
    splice_df <- splice_df[genes_kept, patients_kept]
    if (verbose) {
      print("Simplify Annotations")
    }

    any_hit <- trunc_df
    all_hit <- trunc_df
    all_hit[] <- ""
    all_hit[trunc_df == TRUE] <- "Truncating"
    any_hit[snv_df == TRUE] <- TRUE
    all_hit[snv_df == TRUE] <- "Missense"
    any_hit[splice_df == TRUE] <- TRUE
    all_hit[splice_df == TRUE] <- "Splice"
    if (!missing(cnv_df)) {
      any_hit[cn_df == TRUE] <- TRUE
    }
    if (highlightHotspots) {
      any_hit[hotspot_df == TRUE] <- TRUE
    }
    if (!missing(sortByGenes) && !keepSampleOrder) {
      tany_hit <- t(any_hit) %>% as.data.frame()
      tany_hit <- arrange(tany_hit, desc(across(sortByGenes)))
      any_hit <- t(tany_hit) %>% as.data.frame()
      patients_kept <- colnames(any_hit)
      col_order <- patients_kept
    } else if (keepSampleOrder) {
      patients_kept <- colnames(any_hit)
    } else if (!missing(sortByColumns)) {
      if(verbose){
        print("sorting samples again")
        print("n patients:")
        print(length(patients_kept))
        print("size of metadata:")
        print(dim(metadata_df))
        print(table(patients_kept %in% rownames(metadata_df)))

      }
      #metadata_df <- metadata_df[patients_kept, ]
      metadata_df <- metadata_df[order(match(
        rownames(metadata_df),
        patients_kept
      )), , drop = FALSE]
      col_order <- patients_kept
      if (verbose) {
        print("done!")
      }
    } else {
      if (verbose) {
        print(dim(any_hit))
      }
      # sort genes by frequency first
      gene_sums <- rowSums(any_hit)
      genes_sorted <- names(rev(sort(gene_sums)))
      if (verbose) {
        print("ORDER:")
        print(genes_sorted)
      }
      any_hit <- any_hit[genes_sorted, ]
      tany_hit <- t(any_hit) %>% as.data.frame()
      tany_hit <- arrange(tany_hit, desc(across(genes_sorted)))
      any_hit <- t(tany_hit) %>% as.data.frame()
      patients_kept <- colnames(any_hit)
      col_order <- patients_kept
    }
    if (removeNonMutated) {
      patients_kept <- names(which(colSums(any_hit) > 0))
      if (!missing(cnv_df)) {
        cn_df <- cn_df[genes_kept, patients_kept]
      }
      if (verbose) {
        print("KEPT:")
        print(length(patients_kept))
        print(dim(metadata_df))
      }
      trunc_df <- trunc_df[genes_kept, patients_kept]
      snv_df <- snv_df[genes_kept, patients_kept]
      splice_df <- splice_df[genes_kept, patients_kept]
      metadata_df <- metadata_df[patients_kept, , drop = FALSE]
      if (highlightHotspots) {
        hotspot_df <- hotspot_df[genes_kept, patients_kept]
      }
    }
    col_order <- patients_kept
    if(verbose){
        print("sorting samples yet again")
        print("n patients:")
        print(length(patients_kept))
        print("size of metadata:")
        print(head(metadata_df))
        print(dim(metadata_df))
        print(table(patients_kept %in% rownames(metadata_df)))

    }
    metadata_df <- metadata_df[order(match(
      rownames(metadata_df),
      patients_kept
    )), , drop = FALSE]
    if (verbose) {
      
      print(paste("proceeding with these genes:"))
      print(genes_kept)
      print(head(metadata_df))
    }
    mat_list <- list(Missense = as.matrix(snv_df[
      genes_kept,
      patients_kept
    ]), Truncating = as.matrix(trunc_df[
      genes_kept,
      patients_kept
    ]), Splice_Site = as.matrix(splice_df[
      genes_kept,
      patients_kept
    ]))
    if (!missing(cnv_df)) {
      mat_list[["CNV"]] <- as.matrix(cn_df[genes_kept, patients_kept])
    }
    if (highlightHotspots) {
      mat_list[["HotSpot"]] <- as.matrix(hotspot_df[
        genes_kept,
        patients_kept
      ])
    }
    if(simplify_bg_colour == "pathology"){
      path_mat <- as.matrix(bg_df[genes_kept,patients_kept])
      mat_list[["BL"]] = path_mat
      mat_list[["BL"]][path_mat=="BL"] = TRUE
      mat_list[["BL"]][path_mat!="BL"] = FALSE
      mat_list[["FL"]] = path_mat
      mat_list[["FL"]][path_mat=="FL"] = TRUE
      mat_list[["FL"]][path_mat!="FL"] = FALSE
      mat_list[["DLBCL"]] = path_mat
      mat_list[["DLBCL"]][path_mat=="DLBCL"] = TRUE
      mat_list[["DLBCL"]][path_mat!="DLBCL"] = FALSE
    }
    if (verbose) {
      message("done")
    }
    mat_input <- mat_list
  } else { # end simplify
    if (missing(sortByColumns) && keepSampleOrder == FALSE) {
      if (verbose) {
        print("col_order will be NULL")
      }
      col_order <- NULL
    } else {
      col_order <- patients_kept
    }
    mat_input <- mat[intersect(genes, genes_kept), patients_kept, drop = FALSE]
    
    any_hit <- mat_input
    any_hit[] <- 0
    any_hit[mat_input != ""] <- 1
  }
  if (cluster_rows || cluster_cols) {
    if (!cluster_cols) {
      message("clustering mutation rows but not columns")
      if (verbose) {
        print(paste(
          "clustering for row ordering using this distance metric:",
          clustering_distance_rows
        ))
      }
      h_obj <- pheatmap(any_hit,
        clustering_distance_rows = clustering_distance_rows,
        cluster_rows = TRUE, cluster_cols = FALSE, fontsize_row = 6,
        show_colnames = showTumorSampleBarcode, use_raster = use_raster
      )
      if (verbose) {
        print("Done!")
      }
    } else if (!cluster_rows) {
      message("clustering mutation columns but not rows")
      h_obj <- pheatmap(any_hit,
        clustering_distance_cols = clustering_distance_cols,
        cluster_cols = TRUE, cluster_rows = FALSE, fontsize_row = 6,
        show_colnames = showTumorSampleBarcode, use_raster = use_raster
      )
    } else {
      message("clustering mutation rows and columns")
      h_obj <- pheatmap(any_hit,
        clustering_distance_rows = clustering_distance_rows,
        clustering_distance_cols = clustering_distance_cols,
        fontsize_row = 6, show_colnames = showTumorSampleBarcode,
        use_raster = use_raster
      )
    }
    if (dry_run) {
      print(h_obj)
      return(h_obj)
    }
    row_dend <- row_dend(h_obj)
    col_dend <- column_dend(h_obj)
  } else {
    row_dend <- NULL
    col_dend <- NULL
  }
  if (missing(splitColumnName)) {
    column_split <- rep("", length(patients_kept))
  } else {
    if (is.factor(metadata_df[, splitColumnName])) {
      if (verbose) {
        print("FACTOR!")
      }
      metadata_df <- arrange(metadata_df, splitColumnName)
      column_split <- pull(metadata_df, splitColumnName)
      patients_kept <- rownames(metadata_df)
    } else {
      if (verbose) {
        print("NOT FACTOR")
      }
      metadata_df[
        is.na(metadata_df[, splitColumnName]),
        splitColumnName
      ] <- "NA"
      column_split <- factor(metadata_df[
        patients_kept,
        splitColumnName
      ])
    }
  }
  if (missing(hide_annotations)) {
    metadata_df <- metadata_df
  } else if (hide_annotations_tracks) {
    metadata_df <- metadata_df %>% dplyr::select(-all_of(hide_annotations))
  }
  metadata_df <- metadata_df %>%
    mutate_if(is.factor, as.character) %>%
    replace(is.na(.), "NA")
  for (column in colnames(metadata_df)) {
    if (missing(numericMetadataColumns)) {
      remaining <- unique(metadata_df[column]) %>%
        pull()
      colours[[column]] <- (colours[column] %>%
        unname() %>%
        unlist())[remaining]
    } else if (!column %in% numericMetadataColumns) {
      remaining <- unique(metadata_df[column]) %>% pull()
      colours[[column]] <- (colours[column] %>%
        unname() %>%
        unlist())[remaining]
    }
  }
  if (verbose) {
    print("Rows of metadata:")
    print(nrow(metadata_df))
  }
  if (nrow(mat) == 0) {
    stop(paste(
      "There are no genes left.",
      "Try reducing minMutationPercent"
    ))
  }

  if (simplify_annotation) {
    plot_type <- "simplify"
  } else {
    plot_type <- "original"
  }
  if (verbose) {
    message("calling make_prettyoncoplot()")
  }

  returned <- make_prettyoncoplot(
    mat_input,
    metadata_df,
    these_samples_metadata,
    cluster_oncoplot_rows = cluster_rows,
    cluster_oncoplot_cols = cluster_cols,
    clustering_distance_cols = clustering_distance_cols,
    splitColumnName = splitColumnName,
    column_split = column_split,
    row_split = row_split,
    alter_fun = alter_fun,
    top_annotation = top_annotation,
    right_annotation = right_annotation,
    plot_type = plot_type,
    cnv_df = cnv_df,
    col = col,
    annotation_row = annotation_row,
    annotation_col = annotation_col,
    legend_direction = legend_direction,
    legendFontSize = legendFontSize,
    fontSizeGene = fontSizeGene,
    gene_order = gene_order,
    col_order = col_order,
    legend_row = legend_row,
    legend_col = legend_col,
    col_fun = col_fun,
    colours = colours,
    row_dend = row_dend,
    col_dend = col_dend, column_title = column_title,
    show_legend = show_legend,
    legend_position = legend_position,
    metadataBarHeight = metadataBarHeight,
    metadataBarFontsize = metadataBarFontsize,
    showTumorSampleBarcode = showTumorSampleBarcode,
    return_inputs = return_inputs,
    split_rows_kmeans = split_rows_kmeans,
    split_columns_kmeans = split_columns_kmeans,
    verbose = verbose, use_raster = use_raster,
    pw = plot_width,
    ph = plot_height,
    pct_side = pct_side,
    pctFontSize = pctFontSize,
    row_names_side = row_names_side,
    show_pct = show_pct,
    show_any_legend = show_any_legend,
    metadataSide = metadataSide,
    hide_annotation_name = hide_annotation_name,
    summarizeByColumns = summarizeByColumns,
    longest_label = longest_label,
    right_anno_column = right_anno_column,
    hide_left_annotations = hide_left_annotations,
    coloured_genes = coloured_genes,
    annotation_name_side = annotation_name_side,
    axis_font_size = axis_font_size,
    right_anno_width = right_anno_width
  )
  if (return_inputs) {
    if (simplify_annotation) {
      returned[["all_hit"]] <- all_hit
    }
    returned[["args"]] = all_args
    return(returned)
  }
}

make_prettyoncoplot <- function(
    mat_input,
    metadata_df,
    original_metadata,
    cluster_oncoplot_rows,
    cluster_oncoplot_cols,
    clustering_distance_cols,
    splitColumnName, # nolint: object_name_linter.
    column_split, row_split, alter_fun, top_annotation, right_annotation,
    plot_type = "original", cnv_df,
    heat_mat, col, annotation_row, annotation_col, legend_direction,
    legendFontSize, # nolint: object_name_linter.
    fontSizeGene, # nolint: object_name_linter.
    gene_order,
    col_order,
    legend_row,
    legend_col,
    col_fun,
    colours,
    row_dend,
    col_dend,
    column_title,
    show_legend,
    legend_position,
    metadataBarHeight, # nolint: object_name_linter.
    metadataBarFontsize, # nolint: object_name_linter.
    showTumorSampleBarcode, # nolint: object_name_linter.
    return_inputs,
    split_rows_kmeans,
    split_columns_kmeans,
    verbose, use_raster,
    pw,
    ph,
    show_any_legend,
    pct_side, row_names_side,
    pctFontSize, # nolint: object_name_linter.
    show_pct,
    metadataSide, # nolint: object_name_linter.
    hide_annotation_name,
    summarizeByColumns,
    longest_label,
    right_anno_column,
    hide_left_annotations,
    coloured_genes,
    annotation_name_side,
    axis_font_size,
    right_anno_width) {
  if (plot_type == "simplify") {
    ## Speedier, less detailed plot (3 categories of mutations)
    at <- c("Missense", "Truncating", "Splice_Site")
    if (!missing(cnv_df)) {
      at <- c(at, "CNV")
    }
    heatmap_legend_param <- list(
      title = "Alterations", at = at,
      labels = at, nrow = annotation_row,
      ncol = annotation_col,
      legend_direction = legend_direction,
      labels_gp = gpar(fontsize = legendFontSize)
    )
  } else { # end simplify, set up legend for basic plot type instead
    heatmap_legend_param <- list(
      title = "Alterations", at = c(
        "RNA",
        "3'UTR", "Nonsense_Mutation", "Splice_Site", "Splice_Region",
        "Nonstop_Mutation", "Translation_Start_Site", "In_Frame_Ins",
        "In_Frame_Del", "Frame_Shift_Ins", "Frame_Shift_Del",
        "Multi_Hit", "Missense_Mutation", "Silent", "hot_spot"
      ),
      labels = c(
        "RNA", "3'UTR", "Nonsense Mutation", "Splice Site",
        "Splice Region", "Nonstop Mutation", "Translation Start Site",
        "In Frame Insertion", "In Frame Deletion", "Frame Shift Insertion",
        "Frame Shift Deletion", "Multi Hit", "Missense Mutation",
        "HotSpot", "Silent"
      ), nrow = annotation_row,
      ncol = annotation_col, legend_direction = legend_direction,
      labels_gp = gpar(fontsize = legendFontSize)
    )
  }
  modify_na_elements <- function(x) {
    if (is.character(x)) {
      x[is.na(names(x))] <- "white"
      names(x) <- ifelse(is.na(names(x)), "NA", names(x))
    }
    return(x)
  }
  colours <- lapply(colours, modify_na_elements)
  oncoprint_args <- list(
    mat = mat_input,
    alter_fun = alter_fun,
    top_annotation = top_annotation,
    right_annotation = right_annotation,
    col = col,
    row_order = gene_order,
    column_order = col_order,
    column_labels = NULL,
    show_column_names = showTumorSampleBarcode,
    column_title = column_title,
    row_title = NULL,
    show_heatmap_legend = show_any_legend,
    heatmap_legend_param = heatmap_legend_param,
    
    pct_gp = gpar(fontsize = pctFontSize),
    use_raster = use_raster,
    row_names_side = row_names_side,
    pct_side = pct_side
  )
  
  row_labels = rownames(mat_input[[1]])
  gene_label_colour = setNames(rep("black",length(row_labels)),row_labels)
  if(!missing(coloured_genes)){
    gene_label_colour[coloured_genes] = "red"
    oncoprint_args[['row_names_gp']] =  gpar(fontsize = fontSizeGene, col = gene_label_colour)
  }else{
    oncoprint_args[['row_names_gp']] =  gpar(fontsize = fontSizeGene)
  }

  if(verbose){
    print(gene_label_colour)
  }

  
  
  oncoprint_args[["show_pct"]] <- show_pct
  if (!is.null(pw)) {
    oncoprint_args[["width"]] <- unit(pw, "cm")
  }
  if (!is.null(ph)) {
    oncoprint_args[["height"]] <- unit(ph, "cm")
  }
  #tell oncoprint how much space will be needed for gene names
  if(missing(longest_label)){
    oncoprint_args[["all_row_labels"]] = rownames(mat_input)
  }else{
    oncoprint_args[["all_row_labels"]] = longest_label
  }
  colanno <- ComplexHeatmap::HeatmapAnnotation(
    df = metadata_df,
    show_legend = show_legend,
    show_annotation_name = !hide_annotation_name,
    col = colours, simple_anno_size = unit(
      metadataBarHeight,
      "mm"
    ), gap = unit(0.25 * metadataBarHeight, "mm"),
    annotation_name_gp = gpar(fontsize = metadataBarFontsize),
    annotation_legend_param = list(
      nrow = legend_row, col_fun = col_fun,
      ncol = legend_col, direction = legend_direction,
      labels_gp = gpar(fontsize = legendFontSize)
    )
  )
  if (plot_type == "simplify") {
    any_mut <- mat_input$Missense
    any_mut[] <- 0
    any_mut[mat_input$Missense] <- 1
    any_mut[mat_input$Truncating] <- 1
    any_mut[mat_input$Splice_Site] <- 1
    if ("HotSpot" %in% names(mat_input$mut_mat)) {
      if (verbose) {
        print("Counting hot spots!")
      }
      any_mut[mat_input$HotSpot] <- 1
    }
    mut_status <- as.data.frame(t(any_mut))

    if (!missing(summarizeByColumns)) {
      gg <- rownames(mat_input[[1]])

      # for each row, show a barplot summarizing the incidence of mutations
      # in each category in this column (e.g. pathology)


      get_side_anno = function(summarizeBy,proportional=TRUE,right_anno=FALSE){
          mut_status_long <- rownames_to_column(mut_status, "sample_id") %>%
            pivot_longer(
              cols = -sample_id,
              names_to = "Hugo_Symbol",
              values_to = "mutation"
            )
          mut_status_long <- left_join(
            rownames_to_column(metadata_df, "sample_id"),
            mut_status_long
          ) %>% dplyr::filter(!!sym(summarizeBy) != "NA")

          mut_status_sum <- mut_status_long %>%
            group_by(Hugo_Symbol, .data[[summarizeBy]]) %>%
            summarise(total = sum(mutation), .groups = "drop")

          mut_status_sum <- filter(
            mut_status_sum,
            !!sym(summarizeBy) != "NA"
          )
          scaled = TRUE
          proportional = FALSE
          if(proportional){

            totals = metadata_df %>%
              group_by(.data[[summarizeBy]]) %>%
              count() %>% filter(!!sym(summarizeBy) != "NA")
            gene_totals = mut_status_long %>% group_by(Hugo_Symbol) %>%
              summarise(gene_total=sum(mutation))
            mut_status_sum = left_join(mut_status_sum,gene_totals)

            mut_status_sum = left_join(mut_status_sum,totals) %>%
              mutate(proportion=total/gene_total) #incidence of mutations *in that group*

            mut_status_sum <- mut_status_sum %>%
            select(.,-total:-n) %>%
            pivot_wider(.,
              names_from = {{ summarizeBy }},
              values_from = proportion,
              values_fill = 0
            )

          } else if(scaled){

            totals = original_metadata %>%
              group_by(.data[[summarizeBy]]) %>%
              count() %>% filter(!!sym(summarizeBy) != "NA")
            #get the total for scaling
            mut_status_sum = left_join(mut_status_sum,totals) %>%
              mutate(proportion = total/n)
            
            mut_status_sum_raw = mut_status_sum %>%
            mutate(percent = 100*proportion)
            #stop("===================")
            mut_status_sum <- group_by(mut_status_sum,Hugo_Symbol) %>%
              mutate(total_proportions = sum(proportion)) %>%
              ungroup() %>%
              mutate(adj_proportion = proportion/total_proportions)
            #dplyr::filter(mut_status_sum,Hugo_Symbol %in% c("SLC30A4","MAP2K1","MEF2C")) %>% print()
            mut_status_sum <- mut_status_sum %>%
            select(.,-total,-n,-total_proportions,-proportion) %>%
            rename("proportion"="adj_proportion")
            if(right_anno){
              if(verbose){
                print(head(mut_status_sum_raw))
              }
              
            }
            
            mut_status_sum_raw = mut_status_sum_raw %>%
              select(-total:-proportion) %>%
              pivot_wider(.,
                names_from = {{ summarizeBy }},
                values_from = percent,
                values_fill = 0
              )
            #print(head(mut_status_sum_raw))
            mut_status_sum = mut_status_sum %>%
            pivot_wider(.,
              names_from = {{ summarizeBy }},
              values_from = proportion,
              values_fill = 0
            )
          } else {
            mut_status_sum <- pivot_wider(mut_status_sum,
              names_from = {{ summarizeBy }},
              values_from = total
            )
          }
          mut_status_sum_raw = mut_status_sum_raw %>%
            column_to_rownames("Hugo_Symbol")
          cn = colnames(mut_status_sum)
                mut_status_sum[cn] <- lapply(mut_status_sum[cn], function(col) {
            col[is.nan(col)] <- 0
            col
          })



          mut_status_sum = mut_status_sum %>%
            column_to_rownames("Hugo_Symbol")
          cn = colnames(mut_status_sum)

          mut_status_sum[cn] <- lapply(mut_status_sum[cn], function(col) {
            col[is.nan(col)] <- 0
            col
          })

          anno_col <- map_metadata_to_colours(summarizeBy,metadata_df)[[1]]

          anno_col <- anno_col[colnames(mut_status_sum)]
            if(right_anno){
              return(list(df=mut_status_sum_raw,colours=anno_col))
            }else{
              return(list(df=mut_status_sum,colours=anno_col))
            }
            
      }
      #get the right annotation
      
      if(!missing(right_anno_column)){
        anno_args = get_side_anno(right_anno_column,right_anno = TRUE)
        #print(head(anno_args$df))
        percentages = as.matrix(anno_args$df[gg,])
        right_annotation = HeatmapAnnotation(
          Incidence = anno_barplot(percentages,
            beside=TRUE,
            gp = gpar(fill = anno_args$colours,
            
            col = anno_args$colours),
            axis_param = list(
              #direction = "reverse",
              labels_rot = 0,
              side = "top",
              gp = gpar(fontsize = axis_font_size)
            )),
          width = unit(right_anno_width, "cm"),
          #annotation_name_gp = gpar(fontsize = 0),
          annotation_name_rot = 0,
          annotation_name_side = annotation_name_side,
          annotation_name_offset = unit(1,"cm"),
          annotation_label = NULL,
          which="row"
        )
        oncoprint_args[['right_annotation']]=right_annotation
        oncoprint_args[['pct_gp']]=gpar(fontsize=fontSizeGene)
      }
      #This makes the set of left annotations
      anno_args_list <- lapply(summarizeByColumns, function(col) {
        get_side_anno(col)
      })
      names(anno_args_list) <- paste0("anno_args", seq_along(anno_args_list))
      i = 1
      for(new_anno in names(anno_args_list)){
        oncoprint_args[['col']] = c(oncoprint_args[['col']],
                                    anno_args_list[[new_anno]]$colours)
      }
      if(!hide_left_annotations){
        barplot_annos <- setNames(lapply(seq_along(anno_args_list), function(i) {
          anno_args <- anno_args_list[[i]]
          row_anno_barplot(
            anno_args$df[gg, ], # ensure 'gg' is defined appropriately
            gp = gpar(fill = anno_args$colours, col = anno_args$colours),
            bar_width = 0.5,
            border = FALSE,
            #width = unit(2, "cm"),
            width = unit(right_anno_width, "cm"),
            axis_param = list(
              direction = "reverse",
              labels_rot = 0,
              #side = "top",
              gp = gpar(fontsize = 4)
            )
          )
        }), paste0("anno", seq_along(anno_args_list)))
          left_anno <- do.call(HeatmapAnnotation, c(
            barplot_annos,
            list(
              which = "row",
              gap = unit(0.2, "cm"),
              show_annotation_name = FALSE
            )
          ))
        oncoprint_args[["left_annotation"]] <- left_anno
      }
    }
  }


  if (metadataSide == "bottom") {
    oncoprint_args[["bottom_annotation"]] <- colanno
  } else {
    oncoprint_args[["top_annotation"]] <- colanno
  }
  if (cluster_oncoplot_rows) {
    oncoprint_args[["cluster_rows"]] <- row_dend
  } else {
    # I think row_split is incompatible with row clustering. TBD
    oncoprint_args[["row_split"]] <- row_split
  }
  if (cluster_oncoplot_cols) {
    if (missing(splitColumnName)) {
      oncoprint_args[["cluster_columns"]] <- col_dend
    } else {
      stop(paste(
        "splitColumnName is incompatible with cluster_cols = TRUE",
        "Try prettyStackedOncoplot instead"
      ))
    }
  } else {
    oncoprint_args[["column_split"]] <- column_split
  }
  ch <- do.call("oncoPrint", oncoprint_args)
  ch_ord <- ch@column_order
  names(ch_ord) <- rownames(metadata_df)
  if (return_inputs) {
    if (is.null(gene_order)) {
      # get gene order from Heatmap object
      to_return <- list("gene_order" = ch@row_names_param$labels[ch@row_order])
    } else {
      to_return <- list("gene_order" = gene_order)
    }
    if (plot_type == "simplify") {
      metadata_df <- bind_cols(metadata_df, mut_status)
      return(c(
        to_return,
        list(
          Heatmap = ch,
          mut_mat = mat_input,
          mut_status = mut_status,
          metadata = metadata_df,
          sample_order = ch_ord
        )
      ))
    }
    return(c(
      to_return,
      list(
        Heatmap = ch, mut_mat = mat_input, metadata = metadata_df,
        sample_order = ch_ord
      )
    ))
  }
  suppressMessages(draw(ch,
    heatmap_legend_side = legend_position,
    annotation_legend_side = legend_position
  ))
}


