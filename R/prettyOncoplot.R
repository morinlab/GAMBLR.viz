#' @title Oncooplot
#'
#' @description Create a highly customizable oncoplot.
#'
#' @details Make an oncoplot that is pretty using ComplexHeatmap. The metadata is expected to follow the structure and column naming used in GAMBL.
#' If you provide your own non-GAMBL samples and metadata, you must include at least the following columns with these names.
#' The first one should match the Tumor_Sample_Barcode in the MAF object or onco_matrix you provide.
#' sample_id, pathology
#'
#' @param maf_df A maf as data frame containing the mutations you want to plot.
#' @param onco_matrix_path Provide a path to an onco_matrix file instead of a MAF object if the former is unavailable (this limits functionality a bit).
#' @param genes An optional vector of genes to restrict your plot to.
#' @param include_noncoding List of non-coding regions to be included, default is NULL. Specify like this: include_noncoding=list("NFKBIZ" = c("3'UTR"), "HNRNPH1" = "Splice_Region")
#' @param keepGeneOrder Set to TRUE if you want to preserve the gene order specified.
#' @param keepSampleOrder Set to TRUE if you want to preserve the sample order specified. The default value is FALSE and respects all of the specified ordering.
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
#' @param splitColumnName Optional argument to indicate which metadata column to split on.
#' @param splitGeneGroups Split genes into groups for better seperation (between different gene-groups) in prettyOncoplot.
#' @param legend_row Fiddle with these to widen or narrow your legend.
#' @param legend_col Fiddle with these to widen or narrow your legend.
#' @param showTumorSampleBarcode Optional argument for showing tumor barcode. Default is FALSE.
#' @param groupNames optional vector of group names to be displayed above heatmap. Should be the same length as the number of groups that will be shown. Default is NULL (no labels).
#' @param verbose Set to TRUE to enable verbose mode (debugging messages.
#' @param hide_annotations Hide annotations for specifc ashms. argument takes a list with annotations.
#' @param hide_annotations_tracks When hide_annotations is supplied with a list of columns, this parameter can optionally also not display those columns as the annotation track. Accepts TRUE and FALSE (default).
#' @param annotate_specific_genes Optional argument, specifying whether the features should be labelled according to their significance in one of the pathologies. Default is FALSE (no annotation).
#' @param this_forest_object If annotate_specific_genes is specified, this arguments takes the output of GAMBLR::prettyForestPlot directly to determine the annotations.
#' @param custom_colours Provide named vector (or named list of vectors) containing custom annotation colours if you do not want to use standartized pallette.
#' @param legend_direction Direction of legend, default is "horizontal".
#' @param ylim Limit for y-axis.
#' @param legend_position Position of legend, default is "bottom".
#' @param annotation_row Row for annotations, default is 2.
#' @param annotation_col Column for annotations, default is 1.
#' @param legendFontSize Font size for legend, default is 10.
#' @param cluster_rows Enables clustering of genes with correlated mutation patterns 
#' @param cluster_cols Enables clustering of samples with correlated mutation patterns 
#' @param clustering_distance_rows Distance metric used for clustering when cluster_rows = TRUE
#' @param dry_run Set to TRUE to more efficiently view the clustering result while debugging cluster_rows/clustering_distance_rows
#' @param simplify Collapse/group the variant effect categories to only 3 options. This is a much faster option for when many patients/genes are included.
#'
#' @return When using the simplify option, the function returns a logical matrix indicating the mutation status of each gene and patient shown in the heatmap.
#'
#' @import tidyr dplyr circlize ComplexHeatmap ggplot2 GAMBLR.helpers tibble
#' @export
#'
#' @examples
#' #load packages
#' library(grid)
#' library(GAMBLR.data)
#' library(dplyr)
#'
#' maf_metadata <- get_gambl_metadata(seq_type_filter = "genome") %>%
#'     dplyr::filter(pathology %in% c("FL", "DLBCL"))
#'
#' maf_data <- get_ssm_by_samples(
#'     these_samples_metadata = maf_metadata
#' )
#'
#define some genes of interest
#' fl_genes = c("RRAGC", "CREBBP", "VMA21", "ATP6V1B2")
#'
#' dlbcl_genes = c("EZH2", "KMT2D", "MEF2B", "CD79B", "MYD88", "TP53")
#'
#' genes = c(fl_genes, dlbcl_genes)
#'
#define gene groups
#' gene_groups = c(rep("FL", length(fl_genes)), rep("DLBCL", length(dlbcl_genes)))
#' names(gene_groups) = genes
#'
#' prettyOncoplot(
#'     maf_df = maf_data,
#'     genes = genes,
#'     these_samples_metadata = maf_metadata %>%
#'         arrange(patient_id),
#'     splitGeneGroups = gene_groups,
#'     keepGeneOrder = TRUE,
#'     splitColumnName = "pathology",
#'     metadataBarHeight = 5,
#'     metadataBarFontsize = 8,
#'     legend_row = 2,
#'     fontSizeGene = 11,
#'     metadataColumns = c("pathology", "lymphgen", "sex"),
#'     sortByColumns = c("pathology", "lymphgen", "sex")
#' )
#'
prettyOncoplot = function(maf_df,
                          onco_matrix_path,
                          genes,
                          include_noncoding = NULL,
                          keepGeneOrder = FALSE,
                          keepSampleOrder = FALSE,
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
                          legendFontSize = 10,
                          fontSizeGene = 6,
                          annotation_row = 2,
                          annotation_col = 1,
                          verbose = FALSE,
                          cluster_rows = FALSE,
                          cluster_cols = FALSE,
                          clustering_distance_rows = "binary",
                          dry_run = FALSE,
                          simplify= TRUE){
  
  patients = pull(these_samples_metadata, sample_id)
  if(missing(maf_df)){
    stop(
      "You must provide maf data frame."
    )
  }
  onco_matrix_coding <- GAMBLR.helpers::coding_class[
    !GAMBLR.helpers::coding_class %in% c("Silent", "Splice_Region", "Targeted_Region")
  ]
  #ensure patients not in metadata get dropped up-front to ensure mutation frequencies are accurate
  if(missing(onco_matrix_path)){
    onco_matrix_path = "onco_matrix.txt"
    #order the data frame the way you want the patients shown
    maf_patients = unique(as.character(maf_df$Tumor_Sample_Barcode))
    if(any(!maf_patients %in% patients)){
      extra = maf_patients[which(!maf_patients %in% patients)]
      patients = maf_patients[which(maf_patients %in% patients)]
      n_drop = length(extra)
      message(paste(n_drop, "patients are not in your metadata, will drop them from the data before displaying"))
      maf_df = filter(maf_df, Tumor_Sample_Barcode %in% patients)
    }
    if(missing(genes)){
      #check that our MAFtools object only contains samples in the supplied metadata
      gene_summary = maf_df %>%
        distinct(
          Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification,
          Start_Position, End_Position
        ) %>%
        filter(
          Tumor_Sample_Barcode %in% patients,
          Variant_Classification %in% onco_matrix_coding
        ) %>%
        distinct(Tumor_Sample_Barcode, Hugo_Symbol) %>%
        group_by(Hugo_Symbol) %>%
        summarize(MutatedSamples = n(), .groups = "drop") %>%
        arrange(desc(MutatedSamples))
      genes = gene_summary
      colnames(genes)[2] = "mutload"
      totSamps = as.numeric(length(unique(maf_df$Tumor_Sample_Barcode)))
      genes$fractMutated = genes$mutload / totSamps
      genes = genes %>% filter(fractMutated * 100 >= minMutationPercent) %>% pull(Hugo_Symbol)
      
      lg = length(genes)
      if(!recycleOncomatrix){
        message(paste("creating oncomatrix with", lg, "genes"))
        mat_origin = GAMBLR.helpers::create_onco_matrix(maf_df, genes)
        mat_origin = mat_origin[,!colSums(mat_origin=="") == nrow(mat_origin)]
        tsbs = maf_df %>%
          distinct(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Start_Position, End_Position) %>%
          filter(
            Tumor_Sample_Barcode %in% patients,
            Variant_Classification %in% onco_matrix_coding
          ) %>%
          group_by(Tumor_Sample_Barcode) %>%
          summarize(total = n(), .groups = "drop") %>%
          arrange(desc(total)) %>% pull(Tumor_Sample_Barcode)
        if(!removeNonMutated){
          tsb.include = matrix(data = 0, nrow = nrow(mat_origin), ncol = length(tsbs[!tsbs %in% colnames(mat_origin)]))
          colnames(tsb.include) = tsbs[!tsbs %in% colnames(mat_origin)]
          rownames(tsb.include) = rownames(mat_origin)
          mat_origin = cbind(mat_origin, tsb.include)
        }
        write.table(mat_origin, file = onco_matrix_path, quote = F, sep = "\t")
        if(verbose){
          print(paste("numcases:", length(tsbs)))
        }
      }
      
    }else{
      if(any(duplicated(genes))){
        stop("There are duplicated elements in the provided gene list (@param genes). Please ensure only unique entries are present in this list.")
      }
      gene_summary = maf_df %>%
        distinct(
          Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification,
          Start_Position, End_Position
        ) %>%
        filter(
          Hugo_Symbol %in% genes,
          Tumor_Sample_Barcode %in% patients,
          Variant_Classification %in% onco_matrix_coding
        ) %>%
        distinct(Tumor_Sample_Barcode, Hugo_Symbol) %>%
        group_by(Hugo_Symbol) %>%
        summarize(MutatedSamples = n(), .groups = "drop") %>%
        arrange(desc(MutatedSamples))
      if(!recycleOncomatrix){
        mat_origin = GAMBLR.helpers::create_onco_matrix(maf_df, genes)
        mat_origin = mat_origin[,!colSums(mat_origin=="") == nrow(mat_origin)]
        tsbs = maf_df %>%
          distinct(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Start_Position, End_Position) %>%
          filter(
            Tumor_Sample_Barcode %in% patients,
            Variant_Classification %in% onco_matrix_coding
          ) %>%
          group_by(Tumor_Sample_Barcode) %>%
          summarize(total = n(), .groups = "drop") %>%
          arrange(desc(total)) %>% pull(Tumor_Sample_Barcode)
        if(verbose){
          print(paste("numcases:",length(tsbs)))
          print(paste("numgenes:",length(mat_origin[,1])))
        }
      }
      
      if(!removeNonMutated & !recycleOncomatrix){
        tsb.include = matrix(data = 0, nrow = nrow(mat_origin), ncol = length(tsbs[!tsbs %in% colnames(mat_origin)]))
        colnames(tsb.include) = tsbs[!tsbs %in% colnames(mat_origin)]
        rownames(tsb.include) = rownames(mat_origin)
        mat_origin = cbind(mat_origin, tsb.include)
      }else if(length(include_noncoding) > 0){
        print(
          "You requested to include noncoding mutations and remove non-mutated patients ..."
        )
        these_have_noncoding <- maf_df %>%
          filter(
            Tumor_Sample_Barcode %in% patients,
            Hugo_Symbol %in% names(include_noncoding),
            Variant_Classification %in% unlist(unname(include_noncoding))
          ) %>%
          distinct(
            Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Start_Position, End_Position
          ) %>%
          pull(Tumor_Sample_Barcode)
        tsb.include = matrix(data = 0, nrow = nrow(mat_origin), ncol = length(these_have_noncoding[!these_have_noncoding %in% colnames(mat_origin)]))
        colnames(tsb.include) = these_have_noncoding[!these_have_noncoding %in% colnames(mat_origin)]
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
        these_samples = dplyr::filter(
          maf_df,
          Hugo_Symbol == gene & Variant_Classification == this_vc
        ) %>%
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
  col = GAMBLR.helpers::get_gambl_colours("mutation", alpha = mutAlpha)
  mat[mat == 0]=""
  patients_kept = patients[which(patients %in% colnames(mat))]
  patients_dropped = patients[which(!patients %in% colnames(mat))]
  if(verbose){
    message("====DROPPED=====")
    message(patients_dropped)
  }
  genes_kept = genes[which(genes %in% rownames(mat))]
  if(recycleOncomatrix){
    gene_summary_list = apply(mat,1,function(x){sum(!is.na(x))})
    gene_summary = data.frame(Hugo_Symbol = names(gene_summary_list),MutatedSamples=as.numeric(unname(gene_summary_list)))
    gene_summary = arrange(gene_summary,desc(MutatedSamples))
  }else{
    genes_dropped = genes[which(!genes %in% gene_summary$Hugo_Symbol)]
    for (g in genes_dropped) {
      gene_summary = dplyr::add_row(gene_summary, Hugo_Symbol = g)
    }
    gene_summary <- gene_summary %>% replace(is.na(.), 0)
  }
  if(!missing(minMutationPercent)){
    if(recycleOncomatrix){
      
      warning("mintMutationPercent option is not available when you provide your own oncomatrix. Feel free to implement this if you need it")
      return()
    }
    mutation_counts <- gene_summary %>%
      select(Hugo_Symbol, MutatedSamples)
    numpat = length(patients_kept)
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
  if(simplify){
    #make oncomatrix individually from the MAF
    
    summarize_mutation_by_class = function(mutation_set){
      snv_maf = filter(maf_df,Variant_Classification %in% mutation_set) %>%
        filter(Hugo_Symbol %in% rownames(mat)) %>%
        select(Hugo_Symbol,Tumor_Sample_Barcode) %>% unique()  %>% #keep at most one per group/gene combo
        mutate(Mutated = TRUE)
      
      snv_wide = pivot_wider(snv_maf,names_from = "Tumor_Sample_Barcode",values_from = "Mutated",values_fill = FALSE) %>%
        column_to_rownames("Hugo_Symbol")
      missing_g = rownames(mat)[which(!rownames(mat) %in% rownames(snv_wide))]
      missing_s = colnames(mat)[which(!colnames(mat) %in% colnames(snv_wide))]
      #add missing rows
      for(g in missing_g){
        snv_wide[g,] = FALSE
      }
      #add missing cols
      for(s in missing_s){
        snv_wide[,s] = FALSE
      }
      return(snv_wide[rownames(mat),colnames(mat)])
    }
    
    col["Missense"] = col["Missense_Mutation"]
    col["Truncating"] = col["Nonsense_Mutation"]
    
    
    alter_fun = function(x, y, w, h, v) {
      n = sum(v)  # how many alterations for current gene in current sample
      h = h*0.9
      # use `names(which(v))` to correctly map between `v` and `col`
      if(n) grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h, 
                      gp = gpar(fill = col[names(which(v))], col = NA), just = "top")
    }

    snv_df = summarize_mutation_by_class(mutation_set=c("Missense_Mutation","In_Frame_Del", "In_Frame_Ins","Translation_Start_Site"))
  
    trunc_df = summarize_mutation_by_class(mutation_set=c("Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","Nonstop_Mutation"))
    splice_df = summarize_mutation_by_class(mutation_set = "Splice_Site")
    
    snv_df[trunc_df==TRUE | splice_df==TRUE] = FALSE
    splice_df[trunc_df==TRUE] = FALSE
    any_hit = trunc_df
    any_hit[snv_df == TRUE] = TRUE
    any_hit[splice_df == TRUE] = TRUE
    if(cluster_rows | cluster_cols){
      if(!cluster_cols){
        h_obj = pheatmap(any_hit,
                       clustering_distance_rows = clustering_distance_rows,
                       cluster_cols=F,
                       fontsize_row = 6,show_colnames = F)
      }else if(!cluster_rows){
        h_obj = pheatmap(any_hit,
                         clustering_distance_rows = clustering_distance_rows,
                         cluster_cols=F,
                         fontsize_row = 6,show_colnames = F)
      }else{
        h_obj = pheatmap(any_hit,
                         clustering_distance_rows = clustering_distance_rows,
                         fontsize_row = 6,show_colnames = F)
      }
      if(dry_run){
        print(h_obj)
        return(h_obj)
      }
      row_dend = row_dend(h_obj)
      col_dend = column_dend(h_obj)
    }else{
      row_dend = NULL
    }
    
   
  }else{
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
      }
      ,
      hot_spot = function(x, y, w, h) {
        grid.rect(x, y, w-unit(spacing, "pt"), (height_scaling/5)*h,
                  gp = gpar(fill = "white", col = box_col))
      },
      Silent = function(x, y, w, h) {
        grid.rect(x, y, w-unit(spacing, "pt"), (height_scaling/5)*h,
                  gp = gpar(fill = col["Silent"], col = box_col))
      }
    )
  }
  
  #automagically assign colours for other metadata columns.
  #TO DO: convert the loop below into a "map_metadata_to_colours" function HAS THIS BEEN RESOLVED?
  blood_cols = GAMBLR.helpers::get_gambl_colours("blood", alpha = annoAlpha)
  colours = list()
  clinical_colours = GAMBLR.helpers::get_gambl_colours("clinical")
  all_gambl_colours = GAMBLR.helpers::get_gambl_colours()
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
      these = GAMBLR.helpers::get_gambl_colours("sex", alpha = annoAlpha)
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
      these = GAMBLR.helpers::get_gambl_colours("pos_neg", alpha = annoAlpha)
      these = these[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if("GCB" %in% options){
      these = GAMBLR.helpers::get_gambl_colours("COO", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(column %in% c("pathology")){
      these = GAMBLR.helpers::get_gambl_colours(column, alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these,"NA" = "white")
      }
      colours[[column]] = these
    }else if(grepl("lymphgen", column, ignore.case = TRUE)){
      these = GAMBLR.helpers::get_gambl_colours("lymphgen", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(column == "HMRN"){
      these = GAMBLR.helpers::get_gambl_colours("hmrn", alpha = annoAlpha)
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
    hot_samples = dplyr::filter(maf_df, hot_spot == TRUE & Hugo_Symbol %in% genes) %>%
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
  if(simplify){

    heatmap_legend_param = list(title = "Alterations",
                                at = c("Missense","Truncating","Splice_Site"),
                                labels = c("Missense","Truncating","Splice_Site"),
                                nrow = annotation_row, ncol = annotation_col,
                                legend_direction = legend_direction,
                                labels_gp = gpar(fontsize = legendFontSize))
    
  }else{
    heatmap_legend_param = list(title = "Alterations",
                                at = c("RNA", "3'UTR" , "Nonsense_Mutation", "Splice_Site","Splice_Region", "Nonstop_Mutation", "Translation_Start_Site",
                                       "In_Frame_Ins", "In_Frame_Del", "Frame_Shift_Ins", "Frame_Shift_Del", "Multi_Hit", "Missense_Mutation", "Silent", "hot_spot"),
                                labels = c("RNA", "3'UTR", "Nonsense Mutation", "Splice Site","Splice Region", "Nonstop Mutation", "Translation Start Site",
                                           "In Frame Insertion", "In Frame Deletion", "Frame Shift Insertion", "Frame Shift Deletion",
                                           "Multi Hit", "Missense Mutation", "Silent", "Hotspot"),
                                nrow = annotation_row, ncol = annotation_col,
                                legend_direction = legend_direction,
                                labels_gp = gpar(fontsize = legendFontSize))
  }
  
  if(hideTopBarplot){
    top_annotation = NULL
  }else{
    tally_mutations = maf_df %>%
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
        TRUE ~ "Neither"
      )) %>%
      pull("Enriched in")
    
    right_annotation = rowAnnotation(" " = enrichment_label,
                                     col = list(" " = c(GAMBLR.helpers::get_gambl_colours()[these_comparisons], 
                                                        Neither = "#ACADAF", "NA" = "#000000")),
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
  
  if(missing(hide_annotations)){
    metadata_df = metadata_df
  }else if (hide_annotations_tracks){
    metadata_df = metadata_df %>%
      dplyr:: select(-all_of(hide_annotations))
  }
  
  if(keepSampleOrder){
    patients_kept <- patients_kept[order(
      match(
        patients_kept,
        these_samples_metadata %>%
          filter(Tumor_Sample_Barcode %in% patients_kept) %>%
          pull(Tumor_Sample_Barcode)
      )
    )]
    metadata_df <- metadata_df[
      order(match(rownames(metadata_df), patients_kept)),
      ,
      drop = FALSE
    ]
  }
  metadata_df <- metadata_df %>%
    mutate_if(is.factor, as.character) %>%
    replace(is.na(.), "NA")
  # Only keep the annotation colors for the remaining patients
  for(column in colnames(metadata_df)){
    if(missing(numericMetadataColumns)){
      remaining <- unique(metadata_df[column]) %>% pull()
      colours[[column]] <- (colours[column] %>% unname %>% unlist)[remaining]
    }else if(!column %in% numericMetadataColumns){
      remaining <- unique(metadata_df[column]) %>% pull()
      colours[[column]] <- (colours[column] %>% unname %>% unlist)[remaining]
    }
  }
  if(verbose){
    print("Calling ComplexHeatmap::oncoPrint")
  }
  if(simplify){
    mat_list = list(Missense=as.matrix(snv_df[,patients_kept]),
                    Truncating=as.matrix(trunc_df[,patients_kept]),
                    Splice_Site = as.matrix(splice_df[,patients_kept]))
    mat_input = mat_list
  }else{
    mat_input = mat[intersect(genes, genes_kept),patients_kept]
    any_hit = mat_input
    any_hit[] = FALSE
    any_hit[mat_input != ""] = TRUE
    if(cluster_rows | cluster_cols){
      if(!cluster_cols){
        h_obj = pheatmap(any_hit,
                         clustering_distance_rows = clustering_distance_rows,
                         cluster_cols=F,
                         fontsize_row = 6,show_colnames = F)
      }else if(!cluster_rows){
        h_obj = pheatmap(any_hit,
                         clustering_distance_rows = clustering_distance_rows,
                         cluster_cols=F,
                         fontsize_row = 6,show_colnames = F)
      }else{
        h_obj = pheatmap(any_hit,
                         clustering_distance_rows = clustering_distance_rows,
                         fontsize_row = 6,show_colnames = F)
      }
      if(dry_run){
        print(h_obj)
        return(h_obj)
      }
      row_dend = row_dend(h_obj)
      col_dend = column_dend(h_obj)
    }else{
      row_dend = NULL
    }
  }
  if(cluster_rows){
    if(cluster_cols){
      if(!missing(splitColumnName)){
        message("ignoring splitColumnName because it is incompatible with cluster_cols")
        ch = ComplexHeatmap::oncoPrint(mat_input,
                                       alter_fun = alter_fun,
                                       top_annotation = top_annotation,
                                       right_annotation = right_annotation,
                                       col = col,
                                       row_order = gene_order,
                                       column_order = column_order,
                                       column_labels = NULL,
                                       show_column_names = showTumorSampleBarcode,
                                       #column_split = column_split,
                                       column_title = column_title,
                                       row_title = NULL,
                                       heatmap_legend_param = heatmap_legend_param,
                                       row_names_gp = gpar(fontsize = fontSizeGene),
                                       pct_gp = gpar(fontsize = fontSizeGene),
                                       cluster_rows = row_dend,
                                       cluster_columns = col_dend,
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
      }
    }else{
      ch = ComplexHeatmap::oncoPrint(mat_input,
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
                                   heatmap_legend_param = heatmap_legend_param,
                                   row_names_gp = gpar(fontsize = fontSizeGene),
                                   pct_gp = gpar(fontsize = fontSizeGene),
                                   cluster_rows = row_dend,
                                   #cluster_columns = col_dend,
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
    }
  }else{
    ch = ComplexHeatmap::oncoPrint(mat_input,
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
  }
  
  draw(ch, heatmap_legend_side = legend_position, annotation_legend_side = legend_position)
  if(simplify){
    return(any_hit)
  }
}
