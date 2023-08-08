#' @title Oncooplot
#'
#' @description Create a highly customizable oncoplot.
#'
#' @details Make an oncoplot that is pretty using ComplexHeatmap. The metadata is expected to follow the structure and column naming used in GAMBL.
#' If you provide your own non-GAMBL samples and metadata, you must include at least the following columns with these names.
#' The first one should match the Tumor_Sample_Barcode in the MAF object or onco_matrix you provide.
#' sample_id, pathology
#'
#' @param maftools_obj A maftools object containing the mutations you want to plot.
#' @param onco_matrix_path Provide a path to an onco_matrix file instead of a MAF object if the former is unavailable (this limits functionality a bit).
#' @param genes An optional vector of genes to restrict your plot to.
#' @param include_noncoding List of non-coding regions to be included, default is NULL. Specify like this: include_noncoding=list("NFKBIZ" = c("3'UTR"), "HNRNPH1" = "Splice_Region")
#' @param keepGeneOrder Set to TRUE if you want to preserve the gene order specified.
#' @param keepSampleOrder Set to TRUE if you want to preserve the sample order specified.
#' @param highlightHotspots Set to TRUE to highlight hot spots. Default is FALSE.
#' @param these_samples_metadata Data frame containing metadata for your samples.
#' @param metadataColumns A vector containing the categorical column names you want to plot below.
#' @param numericMetadataColumns A vector containing the numeric columns you want to plot below.
#' @param expressionColumns Optional variable for retreiving expression values for a specific gene(s).
#' @param numericMetadataMax A numeric vector of cutoffs to apply to numeric columns above.
#' @param sortByColumns A vector containing the column names you want to sort columns (patients) on.
#' @param arrange_descending A Boolean parameter. Set to TRUE to sort metadata in descending fashion. Default is FALSE.
#' @param removeNonMutated Set to TRUE to drop unmutated cases.
#' @param minMutationPercent Only genes mutated in more than minMutationPercent \% patients will be included.
#' @param fontSizeGene Font size for gene labels (default 6).
#' @param annoAlpha Optional alpha to apply to annotation colours.
#' @param mutAlpha Optional alpha to apply to mutation colours.
#' @param recycleOncomatrix Set to TRUE most of the time to reuse the oncomatrix saved by maftools.
#' @param box_col Colour of boxes for outlining mutations (can be problematic with larger oncoprints).
#' @param metadataBarHeight Optional argument to adjust the height of bar with annotations. The default is 1.5.
#' @param metadataBarFontsize Optional argument to control for the font size of metadata annotations. The default is 5.
#' @param hideTopBarplot Optional argument for removing top bar plot. Default value is TRUE.
#' @param tally_all_mutations Optional argument. Set to TRUE to tally all mutations. Default is FALSE.
#' @param tally_all_mutations_max Optional argument. Default is 1000.
#' @param hideSideBarplot Optional argument for removing side bar plot. Default value is FALSE.
#' @param splitColumnName Optional argument to indicate which metadata column to split on. Default is set to pathology.
#' @param splitGeneGroups Split genes into groups for better seperation (between different gene-groups) in prettyOncoplot.
#' @param legend_row Fiddle with these to widen or narrow your legend.
#' @param legend_col Fiddle with these to widen or narrow your legend.
#' @param showTumorSampleBarcode Optional argument for showing tumor barcode. Default is FALSE.
#' @param groupNames optional vector of group names to be displayed above heatmap. Should be the same length as the number of groups that will be shown. Default is NULL (no labels).
#' @param verbose Set to TRUE to enable verbose mode (debugging messages.
#' @param hide_annotations Hide annotations for specifc ashms. argument takes a list with annotations.
#' @param annotate_specific_genes Optional argument, specifying whether the features should be labelled according to their significance in one of the pathologies. Default is FALSE (no annotation).
#' @param this_forest_object If annotate_specific_genes is specified, this arguments takes the output of GAMBLR::prettyForestPlot directly to determine the annotations.
#' @param custom_colours Provide named vector (or named list of vectors) containing custom annotation colours if you do not want to use standartized pallette.
#' @param legend_direction Direction of lgend, defualt is "horizontal".
#' @param ylim Limit for y-axis.
#' @param legend_position Position of legend, default is "bottom".
#' @param annotation_row Row for annotations, default is 2.
#' @param annotation_col Column for annotations, default is 1.
#' @param legendFontSize Font size for legend, default is 10.
#'
#' @return Nothing
#'
#' @import tidyr dplyr circlize ComplexHeatmap ggsci ggplot2 maftools tibble
#' @importFrom utils getFromNamespace
#' @export
#'
#' @examples
#' library(grid)
#'
#' #get some data
#' maf_data = get_coding_ssm(seq_type = "genome")
#' maf_metadata = GAMBLR.data::gambl_metadata
#' maf = maftools::read.maf(maf_data, clinicalData = maf_metadata)
#'
#' #define some genes of interest
#' bl_genes = c("NFKBIZ", "ID3", "TP53", "ARID1A", "FBXO11",
#'              "GNA13", "TCF3", "TFAP4", "HNRNPU", "FOXO1",
#'              "CCND3", "SMARCA4", "DDX3X")
#'
#' dlbcl_genes = c("EZH2", "KMT2D", "MEF2B", "CREBBP", "MYD88")
#'
#' genes = c(bl_genes, dlbcl_genes)
#'
#' #define gene groups
#' gene_groups = c(rep("BL", length(bl_genes)), rep("DLBCL", length(dlbcl_genes)))
#' names(gene_groups) = genes
#'
#' #filter metadata
#' maf_metadata = dplyr::filter(maf_metadata,!lymphgen %in% c("COMPOSITE"))
#'
#' #convert metadata column into factor
#' maf_metadata$pathology = as.factor(maf_metadata$pathology)
#'
#' #define order of factors for selected metadata column
#' maf_metadata$pathology = factor(maf_metadata$pathology,
#'                                 levels = c("DLBCL", "BL",
#'                                            "B-ALL", "CLL",
#'                                            "COMFL", "DLBCL-BL-like",
#'                                            "FL", "HGBL",
#'                                            "MCL", "PBL",
#'                                            "SCBC", "UNSPECIFIED"))
#'
#' maf_metadata = with(maf_metadata, maf_metadata[order(pathology),])
#'
#' #create prettyOncoplot
#' prettyOncoplot(maftools_obj = maf,
#'                genes = genes,
#'                these_samples_metadata = maf_metadata,
#'                splitGeneGroups = gene_groups,
#'                keepGeneOrder = TRUE,
#'                splitColumnName = "pathology",
#'                metadataBarHeight = 5,
#'                metadataBarFontsize = 8,
#'                legend_row = 2,
#'                fontSizeGene = 11,
#'                metadataColumns = c("pathology", "lymphgen", "sex", "EBV_status_inf", "cohort"),
#'                sortByColumns = c("pathology", "lymphgen", "sex", "EBV_status_inf", "cohort"))
#'
prettyOncoplot = function(maftools_obj,
                          onco_matrix_path,
                          genes,
                          include_noncoding = NULL,
                          keepGeneOrder = FALSE,
                          keepSampleOrder = TRUE,
                          highlightHotspots = FALSE,
                          these_samples_metadata,
                          metadataColumns,
                          numericMetadataColumns,
                          expressionColumns = c(),
                          numericMetadataMax,
                          sortByColumns,
                          arrange_descending = FALSE,
                          removeNonMutated = FALSE,
                          minMutationPercent,
                          mutAlpha = 1,
                          recycleOncomatrix = FALSE,
                          splitColumnName,
                          splitGeneGroups,
                          showTumorSampleBarcode = FALSE,
                          groupNames,
                          hide_annotations,
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
                          legendFontSize = 10,
                          fontSizeGene = 6,
                          annotation_row = 2,
                          annotation_col = 1,
                          verbose = FALSE){

  # Use non-exported function from maftools
  createOncoMatrix <- getFromNamespace("createOncoMatrix", "maftools")

  patients = pull(these_samples_metadata, sample_id)
  #ensure patients not in metadata get dropped up-front to ensure mutation frequencies are accurate
  if(!recycleOncomatrix & missing(onco_matrix_path)){
    onco_matrix_path = "onco_matrix.txt"
  #order the data frame the way you want the patients shown
    maf_patients = unique(as.character(maftools_obj@data$Tumor_Sample_Barcode))
    if(any(!maf_patients %in% patients)){
      extra = maf_patients[which(!maf_patients %in% patients)]
      patients = maf_patients[which(maf_patients %in% patients)]
      n_drop = length(extra)
      message(paste(n_drop, "patients are not in your metadata, will drop them from the data before displaying"))
      maftools_obj = subsetMaf(maf = maftools_obj, tsb = patients)
    }
    if(missing(genes)){
      #check that our MAFtools object only contains samples in the supplied metadata
      genes = maftools::getGeneSummary(x = maftools_obj)[order(MutatedSamples, decreasing = TRUE)][,.(Hugo_Symbol, MutatedSamples)]
      colnames(genes)[2] = "mutload"
      totSamps = as.numeric(maftools_obj@summary[3, summary])
      genes[,fractMutated := mutload / totSamps]
      genes = genes[fractMutated * 100 >= minMutationPercent, Hugo_Symbol]
      lg = length(genes)
      message(paste("creating oncomatrix with", lg, "genes"))
      om = createOncoMatrix(m = maftools_obj, g = genes, add_missing = TRUE)
      mat_origin = om$oncoMatrix
      tsbs = levels(maftools::getSampleSummary(x = maftools_obj)[,Tumor_Sample_Barcode])
      print(paste("numcases:", length(tsbs)))

      if(!removeNonMutated){
        tsb.include = matrix(data = 0, nrow = nrow(mat_origin), ncol = length(tsbs[!tsbs %in% colnames(mat_origin)]))
        colnames(tsb.include) = tsbs[!tsbs %in% colnames(mat_origin)]
        rownames(tsb.include) = rownames(mat_origin)
        mat_origin = cbind(mat_origin, tsb.include)
      }
      write.table(mat_origin, file = onco_matrix_path, quote = F, sep = "\t")
    }else{
      if(any(duplicated(genes))){
        stop("There are duplicated elements in the provided gene list (@param genes). Please ensure only unique entries are present in this list.")
      }
      om = createOncoMatrix(m = maftools_obj, g = genes, add_missing = TRUE)
      mat_origin = om$oncoMatrix
      tsbs = levels(maftools::getSampleSummary(x = maftools_obj)[,Tumor_Sample_Barcode])
      print(paste("numcases:",length(tsbs)))
      print(paste("numgenes:",length(mat_origin[,1])))
      if(!removeNonMutated){
        tsb.include = matrix(data = 0, nrow = nrow(mat_origin), ncol = length(tsbs[!tsbs %in% colnames(mat_origin)]))
        colnames(tsb.include) = tsbs[!tsbs %in% colnames(mat_origin)]
        rownames(tsb.include) = rownames(mat_origin)
        mat_origin = cbind(mat_origin, tsb.include)
      }
      write.table(mat_origin, file = onco_matrix_path, quote = F, sep = "\t")
    }
  }
  if(missing(onco_matrix_path)){
    onco_matrix_path = "onco_matrix.txt"
  }
  if(!missing(numericMetadataColumns)){
    message(paste0("The column(s) ", numericMetadataColumns, " specified both in metadata and numeric metadata. Plotting as numeric values..."))
    metadataColumns = metadataColumns[!metadataColumns %in% numericMetadataColumns]
  }
  patients = pull(these_samples_metadata, sample_id)
  #because the way MAFtools writes this file out is the absolute worst for compatability
  old_style_mat = read.table(onco_matrix_path, sep = "\t", stringsAsFactors = FALSE)
  mat = read.table(onco_matrix_path, sep = "\t", header = TRUE, check.names = FALSE, row.names = 1, fill = TRUE, stringsAsFactors = F, na.strings = c("NA", ""))
  colnames(old_style_mat) = colnames(mat)
  mat = old_style_mat
  mat[mat==0]=""
  #add the noncoding mutations to this if requested (only for genes and types specified)
  if(length(include_noncoding) > 0){
    all_genes_df = data.frame(Hugo_Symbol = rownames(mat))
    all_samples_df = data.frame(Tumor_Sample_Barcode = colnames(mat))
    for(gene in names(include_noncoding)){
      for(this_vc in unname(include_noncoding[[gene]])){
        message(paste(gene, "and", this_vc))
        these_samples = dplyr::filter(maftools_obj@maf.silent,Hugo_Symbol == gene & Variant_Classification == this_vc) %>%
          dplyr::select(Tumor_Sample_Barcode, Variant_Classification) %>%
          unique() %>%
          pull(Tumor_Sample_Barcode)
        for(samp in these_samples){
          if(samp %in% colnames(mat)){
            if(mat[gene, samp] == ""){
              mat[gene, samp] = this_vc
            }else{
              mat[gene, samp] = paste0(this_vc, ";", mat[gene, samp])
            }
          }
        }
      }
    }
  }
  #annotate hot spots if necessary
  if(missing(metadataColumns)){
    message("you should name at least one metadata column to show as an annotation. Defaulting to pathology")
    metadataColumns = c("pathology")
  }
  if(missing(genes)){
    genes = rownames(mat)
  }
  col = get_gambl_colours("mutation", alpha = mutAlpha)
  mat[mat == 0]=""
  patients_kept = patients[which(patients %in% colnames(mat))]
  patients_dropped = patients[which(!patients %in% colnames(mat))]
  if(verbose){
    message("====DROPPED=====")
    message(patients_dropped)
  }
  genes_kept = genes[which(genes %in% rownames(mat))]
  genes_dropped = genes[which(!genes %in% maftools_obj@gene.summary$Hugo_Symbol)]
  for (g in genes_dropped) {
    maftools_obj@gene.summary = dplyr::add_row(maftools_obj@gene.summary, Hugo_Symbol = g)
  }
  maftools_obj@gene.summary <- maftools_obj@gene.summary %>% replace(is.na(.), 0)
  if(!missing(minMutationPercent)){
    if(! onco_matrix_path == "onco_matrix.txt"){

      warning("mintMutationPercent option is not available when you provide your own oncomatrix. Feel free to implement this if you need it")
      return()
    }
    mutation_counts <- maftools_obj@gene.summary %>%
      select(Hugo_Symbol, MutatedSamples)

    numpat = length(patients)
    mutation_counts = mutate(mutation_counts, percent_mutated = 100 * MutatedSamples / numpat)
    genes_keep = mutation_counts %>%
      dplyr::filter(percent_mutated >= minMutationPercent) %>%
      pull(Hugo_Symbol)

    genes_kept = genes[genes %in% genes_keep]
  }
  mat = mat[,patients_kept]
  mat = mat[which(rownames(mat) %in% genes_kept),]
  spacing = 0
  height_scaling = 1
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
                gp = gpar(fill = "#e6e6e6", col = box_col))
    },
    RNA = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = "#F2ED36", col = box_col))
    },
    `3'UTR` = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = "#F2ED36", col = box_col))
    },
    Intron = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), 0.75* h*height_scaling,
                gp = gpar(fill = col["Nonsense_Mutation"], col = box_col))
    },
    #big blue
    Nonsense_Mutation = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = "#D8A7CA", col = box_col))
    },
    Splice_Site = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Splice_Site"], col = box_col))
    },
    Splice_Region = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Splice_Region"], col = box_col))
    },
    Nonstop_Mutation = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Nonstop_Mutation"], col = box_col))
    },
    Translation_Start_Site = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Translation_Start_Site"], col = box_col))
    },
    In_Frame_Ins = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["In_Frame_Ins"], col = box_col))
    },
    In_Frame_Del = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["In_Frame_Del"], col = box_col))
    },
    #all frame shifts will be the same colour, magenta
    Frame_Shift_Del = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Frame_Shift_Del"], col = box_col))
    },
    Frame_Shift_Ins = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Frame_Shift_Ins"], col = box_col))
    },
    #big red
    Multi_Hit = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Multi_Hit"], col = box_col))
    },
    #small green
    Missense_Mutation = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Missense_Mutation"], col = box_col))
    },
    hot_spot = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), (height_scaling/5)*h,
                gp = gpar(fill = "white", col = box_col))
    },
    Silent = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), (height_scaling/5)*h,
                gp = gpar(fill = col["Silent"], col = box_col))
    }
  )
  #automagically assign colours for other metadata columns.
  #TO DO: convert the loop below into a "map_metadata_to_colours" function HAS THIS BEEN RESOLVED?
  blood_cols = get_gambl_colours("blood", alpha = annoAlpha)
  colours = list()
  clinical_colours = ggsci::get_ash("clinical")
  all_gambl_colours = get_gambl_colours()
  for(column in metadataColumns){
    these_samples_metadata[[column]] = factor(these_samples_metadata[[column]], levels = unique(these_samples_metadata[[column]]))
    options = these_samples_metadata %>%
      arrange(column) %>%
      dplyr::filter(!is.na(column)) %>%
      pull(column) %>%
      unique()

    options = options[!is.na(options)]
    if(verbose){
      print(">>>>>>>")
      print(levels(options))
      print("<<<<<<<")
    }
    if(column == "sex"){
      these = get_gambl_colours("sex", alpha = annoAlpha)
      these = these[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(sum(levels(options) %in% names(clinical_colours)) == length(levels(options))){
      #we have a way to map these all to colours!
      if(verbose){
        message(paste("found colours for", column, "here"))
      }
      these = clinical_colours[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(("positive" %in% options | "POS" %in% options | "yes" %in% options) & length(options) < 4){
      if(verbose){
        print("using pos_neg")
      }
      these = get_gambl_colours("pos_neg", alpha = annoAlpha)
      these = these[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if("GCB" %in% options){
      these = get_gambl_colours("COO", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(column %in% c("pathology")){
      these = get_gambl_colours(column, alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these,"NA" = "white")
      }
      colours[[column]] = these
    }else if(grepl("lymphgen", column, ignore.case = TRUE)){
      these = get_gambl_colours("lymphgen", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(column == "HMRN"){
      these = get_gambl_colours("hmrn", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(sum(levels(options) %in% names(all_gambl_colours)) == length(levels(options))){
      if(verbose){
        message(paste("found colours for", column, "in all_gambl_colours"))
      }
      these = all_gambl_colours[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(length(levels(options)) > 15){
      these = rainbow(length(levels(options)), alpha = annoAlpha)
      names(these) = levels(options)
        colours[[column]] = these
    }else{
      these = blood_cols[sample(c(1:length(blood_cols)), size = length(levels(options)))]
      names(these) = levels(options)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }
  }
  if (! is.null(custom_colours)){
    colours = custom_colours
  }

  if(highlightHotspots){
    hot_samples = dplyr::filter(maftools_obj@data, hot_spot == TRUE & Hugo_Symbol %in% genes) %>%
      dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
      mutate(mutated = "hot_spot") %>%
      unique()

    all_genes_df = data.frame(Hugo_Symbol = rownames(mat))
    all_samples_df = data.frame(Tumor_Sample_Barcode = colnames(mat))
    hs = left_join(all_samples_df, hot_samples)
    hot_mat = hs %>%
      pivot_wider(names_from = "Tumor_Sample_Barcode", values_from = "mutated") %>%
      left_join(all_genes_df,.) %>%
      column_to_rownames("Hugo_Symbol") %>%
      as.matrix()

    #annotate hotspots in matrix
    for (i in colnames(mat)){
      mat[genes, i][!is.na(hot_mat[genes, i])] = paste0(mat[genes, i][!is.na(hot_mat[genes, i])], ";", hot_mat[genes, i][!is.na(hot_mat[genes, i])])
    }
    colours[["hot_spots"]] = c("hot_spot" = "magenta")
  }
  if(verbose){
    print(colours) #eventually get rid of this once the bugs are gone
  }
  if(!missing(numericMetadataColumns)){
    metadata_df = dplyr::filter(these_samples_metadata, sample_id %in% patients_kept) %>%
      column_to_rownames("sample_id") %>%
      dplyr::select(all_of(c(metadataColumns, numericMetadataColumns, expressionColumns)))

    if(!missing(numericMetadataMax)){
      max_list = setNames(numericMetadataMax, numericMetadataColumns)

      metadata_df = metadata_df %>%
        mutate(across(names(max_list), ~ ifelse(.x > max_list[[cur_column()]], max_list[[cur_column()]], .x)))
    }
  }else{
    metadata_df = dplyr::filter(these_samples_metadata, sample_id %in% patients_kept) %>%
      column_to_rownames("sample_id") %>%
      dplyr::select(all_of(c(metadataColumns, expressionColumns)))
  }
  metadata_df = metadata_df %>%
    mutate(across(all_of(expressionColumns), ~ trim_scale_expression(.x)))
  if(!missing(sortByColumns)){
    if (arrange_descending) {
      metadata_df = arrange(metadata_df, across(sortByColumns, desc))
    } else {
      metadata_df = arrange(metadata_df, across(sortByColumns))
    }
    patients_kept = rownames(metadata_df)
  }
  if(verbose){
    print(genes_kept)
  }
  col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  for(exp in expressionColumns){
    colours[[exp]] = col_fun
  }
  if(missing(splitColumnName)){
    column_split = rep("", length(patients_kept))
  }else{
    column_split = factor(metadata_df[patients_kept, splitColumnName])
  }
  if(missing(splitGeneGroups)){
    row_split = rep("", length(genes))
  }else{
    row_split = factor(splitGeneGroups[genes], levels = unique(splitGeneGroups[genes]))
  }
  if(!missing(groupNames)){
    column_title = groupNames
  }else{
    column_title = NULL
  }
  if(keepGeneOrder){
    gene_order = genes
  }else{
    gene_order = NULL
  }
  if(missing(hide_annotations)){
    show_legend = rep(TRUE, length(colnames(metadata_df)))
  }else{
    show_legend = rep(TRUE, length(colnames(metadata_df)))
    names(show_legend) = colnames(metadata_df)
    show_legend[hide_annotations] = FALSE
  }
  if(missing(sortByColumns)){
    column_order = NULL
  }else{
    column_order = patients_kept
  }
  heatmap_legend_param = list(title = "Alterations",
                         at = c("RNA", "3'UTR" , "Nonsense_Mutation", "Splice_Site","Splice_Region", "Nonstop_Mutation", "Translation_Start_Site",
                         "In_Frame_Ins", "In_Frame_Del", "Frame_Shift_Ins", "Frame_Shift_Del", "Multi_Hit", "Missense_Mutation", "Silent", "hot_spot"),
                         labels = c("RNA", "3'UTR", "Nonsense Mutation", "Splice Site","Splice Region", "Nonstop Mutation", "Translation Start Site",
                         "In Frame Insertion", "In Frame Deletion", "Frame Shift Insertion", "Frame Shift Deletion",
                         "Multi Hit", "Missense Mutation", "Silent", "Hotspot"),
                         nrow = annotation_row, ncol = annotation_col,
                         legend_direction = legend_direction,
                         labels_gp = gpar(fontsize = legendFontSize))
  if(hideTopBarplot){
    top_annotation = NULL
  }else{
    tally_mutations = maftools_obj@data %>%
      dplyr::filter(Tumor_Sample_Barcode %in% patients_kept) %>%
      group_by(Tumor_Sample_Barcode) %>%
      summarize(n_mutations = n()) %>%
      ungroup %>%
      arrange(match(Tumor_Sample_Barcode, patients_kept)) %>%
      select(n_mutations) %>%
      mutate(n_mutations = ifelse(n_mutations > tally_all_mutations_max,
                                  tally_all_mutations_max,
                                  n_mutations))

    if(is.null(ylim) & ! tally_all_mutations){
      top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot())

    }else if (!is.null(ylim) & ! tally_all_mutations){
      top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(ylim=ylim))

    } else if (is.null(ylim) & tally_all_mutations) {
      top_annotation = columnAnnotation(" " = anno_barplot(tally_mutations))
    } else if (! is.null(ylim) & tally_all_mutations) {
      top_annotation = columnAnnotation(" " = anno_barplot(tally_mutations, ylim=ylim))
    }
  }

  # Handle right annotation for specific genes
  if (annotate_specific_genes & is.null(this_forest_object)) {
    message("WARNING: You requested right annotation, but forgot to provide output of GAMBLR::prettyForestPlot")
    message("No right annotation will be drawn.")
    right_annotation = NULL
  } else if (annotate_specific_genes) {

    these_comparisons = this_forest_object$mutmat$comparison %>% levels

    enrichment_label =
      mat[intersect(genes, genes_kept),patients_kept] %>%
      rownames_to_column("gene") %>%
      select(gene) %>%
      left_join(this_forest_object$fisher %>% select(gene, estimate, q.value)) %>%
      mutate("Enriched in" = case_when(
        estimate == "Inf" & q.value <= 0.1 ~ these_comparisons[1],
        estimate == "-Inf" & q.value <= 0.1 ~ these_comparisons[2],
        is.na(estimate) ~ "NA",
        estimate<=1 & q.value <= 0.1 ~ these_comparisons[2],
        estimate > 1 & q.value <= 0.1 ~ these_comparisons[1],
        TRUE ~ "Both"
      )) %>%
      pull("Enriched in")

    right_annotation = rowAnnotation(" " = enrichment_label,
                             col = list(" " = c(get_gambl_colours()[these_comparisons], Both = "#ACADAF", "NA" = "#000000")),
                             simple_anno_size = unit(metadataBarHeight, "mm"),
                             annotation_legend_param =
                               list(title = "Enriched in",
                                    nrow=legend_row,
                                    ncol = legend_col,
                                    direction=legend_direction,
                                    labels_gp = gpar(fontsize = legendFontSize)))
  } else {
      if(hideSideBarplot){
        right_annotation = NULL
      }else{
        right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot())
      }
  }

  ch = ComplexHeatmap::oncoPrint(mat[intersect(genes, genes_kept),patients_kept],
                                 alter_fun = alter_fun,
                                 top_annotation = top_annotation,
                                 right_annotation = right_annotation,
                                 col = col,
                                 row_order = gene_order,
                                 column_order = column_order,
                                 column_labels = NULL,
                                 show_column_names = showTumorSampleBarcode,
                                 column_split = column_split,
                                 column_title = column_title,
                                 row_title = NULL,
                                 row_split = row_split[intersect(genes, genes_kept)],
                                 heatmap_legend_param = heatmap_legend_param,
                                 row_names_gp = gpar(fontsize = fontSizeGene),
                                 pct_gp = gpar(fontsize = fontSizeGene),
                                 bottom_annotation = ComplexHeatmap::HeatmapAnnotation(df = metadata_df,
                                                                                       show_legend = show_legend,
                                                                                       col = colours,
                                                                                       simple_anno_size = unit(metadataBarHeight, "mm"),
                                                                                       gap = unit(0.25 * metadataBarHeight, "mm"),
                                                                                       annotation_name_gp = gpar(fontsize = metadataBarFontsize),
                                                                                       annotation_legend_param = list(nrow = legend_row,
                                                                                                                      col_fun = col_fun,
                                                                                                                      ncol = legend_col,
                                                                                                                      direction = legend_direction,
                                                                                                                      labels_gp = gpar(fontsize = legendFontSize))))

    draw(ch, heatmap_legend_side = legend_position, annotation_legend_side = legend_position)
}
