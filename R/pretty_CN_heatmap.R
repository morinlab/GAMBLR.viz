


#' Pretty Copy Number Heatmap
#'
#' @details
#' **Correcting for ploidy.** To show copy number relative to each sample's
#' own ploidy, build the matrix with
#' `GAMBLR.utils::segmented_data_to_cn_matrix(..., adjust_for_ploidy = TRUE)`
#' and leave `scale_by_sample = FALSE` (the default). That adjustment
#' recentres each sample's segments on its dominant CN level *before* the
#' bins are made, so bins with no data are filled after the adjustment and
#' stay neutral.
#'
#' Avoid `scale_by_sample = TRUE` on a matrix built with
#' `fill_missing_with = "diploid"` (the default fill): the filled bins are a
#' fixed CN of 2, so after subtracting a high-ploidy sample's mean they
#' appear as deletions in that sample only (commonly on chrY, pericentromeric
#' 9q and the acrocentric p arms). The mean is also pulled around by large
#' gains and losses. Don't combine the two approaches either, as that shifts
#' samples twice.
#'
#' @param cn_state_matrix Copy number matrix with samples as rows and bins as
#' columns, e.g. from [GAMBLR.utils::segmented_data_to_cn_matrix]. Build it
#' with `adjust_for_ploidy = TRUE` to correct for sample ploidy (see Details).
#' @param scale_by_sample Set to TRUE to recentre each sample on 2 by
#' subtracting its mean CN. Not recommended: use
#' `adjust_for_ploidy = TRUE` when building the matrix instead, and don't
#' combine the two (see Details). Default FALSE.
#' @param these_samples_metadata The output of get_gambl_metadata
#' @param keep_sample_order FALSE. Set to TRUE to ensure samples are in the same
#' order as in the metadata
#' @param metadataColumns One or more columns from the metadata you want to
#' display beside the heatmap
#' @param numericMetadataColumns One or more columns from the metadata that
#' should be considered numeric
#' @param expressionColumns Optional: One or more columns from the metadata that include gene expression values you want shown
#' @param geneBoxPlot Optional: Specify the Hugo symbol of a single gene to embed box plots adjacent to the heatmap. Expression data for this gene must be present in the metadata in a column of the same name.
#' @param show_column_names Set to TRUE to display the sample_id of every sample shown in the heatmap
#' @param show_row_names Set to TRUE to display the ID of every bin (region) shown in the heatmap
#' @param keep_these_chromosomes A vector of chromosome names to include (all others will be excluded)
#' @param hide_these_chromosomes A vector of chromosome names to exclude (all others will be included unless keep_these_chromosomes is specified)
#' @param keep_these_bins A vector of bin names to include (all others will be excluded)
#' @param hide_annotations A vector of annotation names to suppress from legends in the plot
#' @param cluster_columns Set to TRUE to enable clustering of genomic regions (columns) based on their CN value across all patients in the heatmap
#' @param cluster_regions More intuitive alias for cluster_columns
#' @param cluster_samples More intuitive alias for cluster_rows, especially when combining with rotate = TRUE
#' @param cluster_rows Set to TRUE to enable clustering of genomic regions (columns) based on their CN value across all regions in the heatmap
#' @param clusterRowsByGenes Optional: A vector of Hugo gene symbols. When provided (and cluster_rows is TRUE), the rows (samples) are clustered using only the bins that contain these genes rather than all bins. Genes are mapped to bins the same way as for labelTheseGenes. Bins outside keep_these_chromosomes can still be used for clustering.
#' Alternatively, a data frame with columns "gene_id" and "cn_thresh" (the same format as genes_and_cn_threshs in get_cnv_and_ssm_status) restricts each gene to its relevant direction of CN change.
#' A cn_thresh below 2 is the maximum CN considered a loss for that gene (tumour suppressor) and above 2 is the minimum CN considered a gain (oncogene).
#' Any value in that gene's bin not passing its threshold is set to 2 (neutral) before clustering. A cn_thresh of 2 uses the bin without any direction filter. The displayed heatmap values are not altered.
#' @param sortByGenes Optional: genes whose bins (mapped as for labelTheseGenes) order the rows by their CN value, highest first, so gains sort first and losses last. Row clustering is turned off. For direction-aware sorting use sortByFeatures.
#' @param sortByFeatures Optional: names of CN features (columns of `cnv_status`) to order the rows by, in priority order. Samples are sorted by the level of each feature (strongest first: 2, then 1, then 0, using the absolute value so gains and losses both sort ahead of unaltered samples), then by the next feature, and so on. Remaining ties keep the clustering order (or the input order if cluster_rows = FALSE). Row clustering is turned off.
#' @param cnv_status Sample-by-feature CN calls used by sortByFeatures, with sample IDs as row names, e.g. from `GAMBLR.utils::segmented_data_to_cnv_status()`. Samples missing from it count as unaltered.
#' @param sortByBins Optional: A vector containing one or more names of genomic bins that will be used to order the heatmap rows.
#' @param sortByPGA Optional: Sort the rows based on percent genome altered (PGA) instead of the other options
#' @param splitByBinState Optional: A single genomic bin that will be used to split the heatmap based on the CN state of that bin
#' @param sortByMetadataColumns A vector containing one or more names of columns from your metadata that will be used to order the rows overall or within slices (if combined with splitByBinState or geneBoxPlot)
#' @param labelTheseGenes A vector of Hugo gene symbols whose location will be indicated on the top of the heatmap
#' @param labelTheseFeatures Optional data frame of CNV features to label, in the format returned by [GAMBLR.utils::get_cnv_features] (subset it to the features you want). Each feature is labelled at the bin containing its midpoint; features sharing a bin are joined with "/". Can be combined with labelTheseGenes.
#' @param feature_label What to show for each feature in labelTheseFeatures: "gene" (default; representative_gene, or region_name where there is none), "name" (e.g. "2p16.1_gain") or "region" (e.g. "2p16.1").
#' @param colour_feature_labels Colour feature labels by direction (gains red, losses blue; black where features of both directions share a bin). Default TRUE.
#' @param heatmap_name Name of the ComplexHeatmap object (default "CN"). Heatmaps combined into one figure need distinct names.
#' @param heatmap_row_title Optional title shown to the left of the heatmap rows (not used with rotate = TRUE).
#' @param cumulative_as_fraction Plot the gain/loss tracks as the fraction of samples shown instead of counts. Default FALSE.
#' @param cumulative_ylim Optional y-axis limits for the gain/loss tracks, e.g. to share one scale across several heatmaps.
#' @param cumulative_position Where to put the gain/loss tracks: "bottom" (default, below the heatmap with the chromosome strip) or "top" (above the heatmap, under any labels, so labels point at the peaks). "top" isn't supported with rotate = TRUE.
#' @param show_mean_cn Show the mean CN strip above the heatmap. Default TRUE.
#' @param cumulative_axis Axis on the gain/loss tracks: "full" (default), "minimal" (only 0 and the top of the scale, legible on short tracks) or "none".
#' @param show_chromosome_strip Show the chromosome colour strip below the heatmap. Default TRUE.
#' @param heatmap_row_title_fontsize,heatmap_row_title_rot Font size (default 13) and rotation (default 90) of heatmap_row_title.
#' @param cumulative_track_height Height in cm of each of the gain and loss tracks (their width with rotate = TRUE). Default 1, ComplexHeatmap's default.
#' @param annotation_colours Optional named list of colour mappings for the sample annotations (named vectors for categorical columns, circlize::colorRamp2 functions for numeric ones such as PGA). These replace the automatically chosen colours for those annotations.
#' @param show_legends Set to FALSE to omit all legends. Default TRUE.
#' @param label_style How gene, feature and bin labels are drawn. "mark" (default) uses a ComplexHeatmap anno_mark annotation: labels are aligned along one line, spread apart where crowded with connector lines to their bins, and given their own space above the heatmap (or to the left with rotate = TRUE). "classic" draws them on the heatmap body as before, positioned with bin_label_nudge.
#' @param bin_label_fontsize Font size for the gene labels (default 5)
#' @param bin_label_nudge Increase or decrease this value to shift the gene labels up/down (default 1.03). Only used with label_style = "classic".
#' @param bin_label_rotation Rotate the direction of the bin label. Default is 45.
#' @param drop_if_PGA_below Lower limit for proportion of genome altered (PGA). Samples below this value will be dropped (default 0)
#' @param drop_if_PGA_above Upper limit for proportion of genome altered (PGA). Samples above this value will be dropped (default 1)
#' @param show_bottom_annotation_name set to TRUE to label the bottom annotation tracks with more details
#' @param bottom_annotation_name_side If using show_bottom_annotation_name, set this to "left" or "right" to relocate the names
#' @param bin_labels Instead of automatically labeling genes, you can instead explicitly provide a list of labels for any bins in the heatmap. The names of each element should match a bin. The value of each element is the label that will be shown. This option can be used to skip gene location look-ups (see examples).
#' @param focus_on_these_bins Mask all regions outside these bins (set CN to 0). Useful for visualizing GISTIC results.
#' @param legend_direction Which orientation to use for the legend
#' @param legend_position Where to put the legend
#' @param legend_row How many rows for the legend layout
#' @param legend_col How many columns for the legend layout
#' @param boxplot_orientation Either "horizontal" or "vertical" (default: horizontal)
#' @param left_annotation_name_side Which side to put the name of the metadata annotations (top or bottom)
#' @param drop_bin_if_sd_below Force bins with standard deviation below this value to be excluded. Set to -Inf to keep all bins.
#' @param return_data Specify TRUE to get some of the internal data back including the heatmap object
#' @param rotate Optionally, flip the rows/columns of resulting heatmap. Default is FALSE.
#' @param width Width of the heatmap in cm when rotate = TRUE. Default 15; NULL fills the available space.
#' @param height Height of the heatmap in cm when rotate = TRUE. Default 6; NULL fills the available space.
#' @param verbose Control verbosity of the console output. Default is FALSE.
#'
#' @return list (when return_data = TRUE)
#' @import GAMBLR.helpers dplyr
#' @export
#'
#' @examples
#' cat("Running example for function: pretty_CN_heatmap\n")
#' suppressMessages(library(dplyr))
#' suppressPackageStartupMessages(library(GAMBLR.open))
#' #get some metadata for subsetting the data to just one pathology (DLBCL)
#' dlbcl_genome_meta = suppressMessages(get_gambl_metadata()) %>%
#'                     filter(pathology=="DLBCL",
#'                     seq_type=="genome")
#' 
#' #remove any duplicate sample_id/seq_type combinations
#' meta_clean = check_and_clean_metadata(dlbcl_genome_meta,
#'                                       duplicate_action = "keep_first")
#' 
#' # Create the copy number matrix using the helper functions.
#' # adjust_for_ploidy = TRUE shows CN relative to each sample's own ploidy;
#' # use it instead of scale_by_sample (see Details).
#' all_segments = get_cn_segments(these = meta_clean)
#' dlbcl_cn_binned = segmented_data_to_cn_matrix(
#'                                   seg_data = all_segments,
#'                                   strategy="auto_split",
#'                                   n_bins_split=1300,
#'                                   adjust_for_ploidy = TRUE,
#'                                   these_samples_metadata = meta_clean)
#'
#' # Generate a basic genome-wide CN heatmap
#' pretty_CN_heatmap(cn_state_matrix=dlbcl_cn_binned,
#'                   these_samples_metadata = meta_clean,
#'                   hide_annotations = "chromosome")
#'
#' # Disable row (sample) clustering and restrict to a few chromosomes
#' # and highlight some genes of interest
#' pretty_CN_heatmap(cn_state_matrix=dlbcl_cn_binned,
#'   these_samples_metadata = meta_clean,
#'   hide_annotations = "chromosomes",
#'   keep_these_chromosomes = c("9","17"),
#'   cluster_rows=FALSE,
#'   labelTheseGenes = c("CDKN2A","TP53"))
#'
#' \dontrun{
#'  # get gene expression data
#'  gene_exp_all = get_gene_expression(all_genes=T,
#'                                     lazy_join=T,
#'                                      arbitrarily_pick = TRUE,
#'                                      HGNC=T,format="wide")
#' 
#'  genome_meta_exp = left_join(get_gambl_metadata() %>%
#'     dplyr::filter(seq_type=="genome") %>%
#'     dplyr::select(sample_id,pathology,lymphgen),
#'     dplyr::select(gene_exp_all,-sample_id),
#'     by=c("sample_id"="genome_sample_id")) %>%
#'      filter(!is.na(MYC))
#' }
#' # Include gene expression data and embed a box plot showing the expression of one gene across different CN states
#'
#' \dontrun{
#'  pretty_CN_heatmap(cn_state_matrix=all_states_binned,
#'    these_samples_metadata = filter(genome_meta_exp,pathology=="DLBCL"),
#'    hide_annotations = "chromosomes",
#'    cluster_rows=F,
#'    geneBoxPlot = "TP53",
#'    boxplot_orientation="horizontal",bin_label_fontsize = 9,bin_label_nudge = 19
#'  )
#' }
#'
#'
pretty_CN_heatmap = function(cn_state_matrix,
                             scale_by_sample=FALSE,
                             these_samples_metadata,
                             keep_sample_order = FALSE,
                             metadataColumns = c("pathology"),
                             numericMetadataColumns,
                             expressionColumns,
                             genome_build = "grch37",
                             cluster_columns=FALSE,
                             cluster_rows=TRUE,
                             clusterRowsByGenes,
                             show_row_names=FALSE,
                             show_column_names=FALSE,
                             keep_these_chromosomes,
                             hide_these_chromosomes,
                             keep_these_bins,
                             hide_annotations,
                             hide_cumulative_plot = FALSE,
                             sortByBins,
                             sortByGenes,
                             sortByFeatures,
                             cnv_status,
                             splitByBinState,
                             sortByPGA=FALSE,
                             sortByMetadataColumns,
                             labelTheseGenes,
                             labelTheseCytobands,
                             labelTheseFeatures,
                             feature_label = "gene",
                             colour_feature_labels = TRUE,
                             label_style = "mark",
                             heatmap_name = "CN",
                             heatmap_row_title = NULL,
                             cumulative_as_fraction = FALSE,
                             cumulative_ylim = NULL,
                             cumulative_track_height = 1,
                             cumulative_position = "bottom",
                             show_mean_cn = TRUE,
                             cumulative_axis = "full",
                             show_chromosome_strip = TRUE,
                             heatmap_row_title_fontsize = 13,
                             heatmap_row_title_rot = 90,
                             annotation_colours = NULL,
                             show_legends = TRUE,
                             highlightTheseRegions,
                             bin_label_fontsize=5,
                             bin_label_nudge = 1.03,
                             bin_label_rotation=45,
                             drop_if_PGA_below=0,
                             drop_if_PGA_above=1,
                             focus_on_these_bins,
                             geneBoxPlot,
                             show_bottom_annotation_name = FALSE,
                             bottom_annotation_name_side = "left",
                             left_annotation_name_side = "top",
                             bin_labels,
                             legend_direction = "horizontal",
                             legend_position = "bottom",
                             legend_row=2,
                             legend_col=3,
                             metadataBarFontsize = 5,
                             metadataBarHeight = 1.5,
                             boxplot_orientation="vertical",
                             return_data = FALSE,
                             rounded = TRUE,
                             drop_bin_if_sd_below=0,
                             flip=FALSE,
                             max_CN_allowed = 6,
                             verbose=FALSE,
                             rotate = FALSE,
                             width = 15,
                             height = 6,
                             cluster_samples,
                             cluster_regions){
  if(!cumulative_position %in% c("bottom", "top")){
    stop("cumulative_position must be \"bottom\" or \"top\"")
  }
  if(!cumulative_axis %in% c("full", "minimal", "none")){
    stop("cumulative_axis must be \"full\", \"minimal\" or \"none\"")
  }
  if(cumulative_position == "top" && rotate){
    stop("cumulative_position = \"top\" isn't supported with rotate = TRUE")
  }
  if(!label_style %in% c("mark", "classic")){
    stop("label_style must be \"mark\" or \"classic\"")
  }
  # distances can't use missing bins (e.g. from fill_missing_with = "nothing"),
  # so cluster with them treated as neutral; they are still drawn as missing
  neutral_for_clustering = function(m){
    m = as.matrix(m)
    m[is.na(m)] = 2
    m
  }
  #aliases
  if(!missing(cluster_samples)){
    cluster_rows = cluster_samples
  }
  if(!missing(cluster_regions)){
    cluster_columns = cluster_regions
  }

  cn_state_matrix[cn_state_matrix>max_CN_allowed] = max_CN_allowed
  #determine variance of each bin across the entire data set
  bin_vars <- apply(cn_state_matrix, 2, function(x) sd(x, na.rm = TRUE))
  
  #bin_vars = apply(mat1,2,sd)
  old_col = ncol(cn_state_matrix)
  # -Inf switches the filter off (e.g. when bins were chosen for a whole
  # cohort before splitting it, or for groups of one sample, whose sd is NA)
  if(is.finite(drop_bin_if_sd_below)){
    cn_state_matrix = cn_state_matrix[,which(bin_vars > drop_bin_if_sd_below)]
  }
  new_col = ncol(cn_state_matrix)
  if(verbose){
    print(paste("original:",old_col,"new:",new_col))
  }
  if(scale_by_sample){
    cn_state_matrix = cn_state_matrix - rowMeans(cn_state_matrix, na.rm = TRUE) + 2
    if(rounded){
      cn_state_matrix = round(cn_state_matrix)
    }
    cn_state_matrix[cn_state_matrix<0]=0
  }
  if(!missing(focus_on_these_bins)){
    message("focus_on_these_bins provided. All other bins will be set to diploid!")
    other_bins = colnames(cn_state_matrix)[!colnames(cn_state_matrix) %in% focus_on_these_bins]
    cn_state_matrix[,other_bins] = 2
  }
  map_bin_to_bin = function(query_region,
                            regions=colnames(cn_state_matrix),
                            first=TRUE){
    these_coords = suppressMessages(region_to_chunks(query_region))
    these_coords$chromosome = gsub("chr", "", these_coords$chromosome)
    all_matches = c()
    for(r in regions){
      region_coords = region_to_chunks(r)
      region_coords$chromosome = gsub("chr", "", region_coords$chromosome)
      if(these_coords$chromosome == region_coords$chromosome){
        if(((as.integer(these_coords$start) > as.integer(region_coords$start)) &
           (as.integer(these_coords$start) < as.integer(region_coords$end))) ||
           (as.integer(these_coords$start) < as.integer(region_coords$start) &
           as.integer(these_coords$end) > as.integer(region_coords$end))){
          if(first){
            #just return the first match
            return(r)
          }else{
            all_matches = c(all_matches,r)
          }
        }
      }
    }
    if(first){
      return(NA)
    }else{
      return(all_matches)
    }

  }
  if(missing(bin_labels)){
    bin_labels = list()
    cytoband_labels = list()
  }

  if(!missing(geneBoxPlot)){
    if(missing(labelTheseGenes)){
      labelTheseGenes = geneBoxPlot
      expressionColumns = geneBoxPlot
      message(paste("setting expressionColumns and labelTheseGenes to",geneBoxPlot))

    }
  }
  regions_highlight = c()
  if(!missing(highlightTheseRegions)){
    for(a in c(1:nrow(highlightTheseRegions))){
      this_region_bins = map_bin_to_bin(highlightTheseRegions[a,"name"],first=FALSE)

      regions_highlight = c(regions_highlight,this_region_bins)

    }
  }
  gene_bins = list()
  cluster_gene_df = NULL
  if(!missing(clusterRowsByGenes)){
    if(is.data.frame(clusterRowsByGenes)){
      stopifnot('clusterRowsByGenes data frame must have columns "gene_id" (characters) and "cn_thresh" (numeric)' = {
        all(c("gene_id","cn_thresh") %in% names(clusterRowsByGenes)) &&
          is.character(clusterRowsByGenes$gene_id) &&
          is.numeric(clusterRowsByGenes$cn_thresh)
      })
      cluster_genes = clusterRowsByGenes$gene_id
      cluster_threshs = clusterRowsByGenes$cn_thresh
    }else{
      cluster_genes = clusterRowsByGenes
      cluster_threshs = rep(2,length(clusterRowsByGenes))
    }
    message("mapping clustering genes to bins")
    for(j in seq_along(cluster_genes)){
      g = cluster_genes[j]
      gene_region = suppressMessages(gene_to_region(g))
      if(length(gene_region)==0){
        message(paste("no region for",g))
        next
      }
      this_gene_region = map_bin_to_bin(gene_region)
      if (is.na(this_gene_region)){
        message(paste("no region for",g))
        next
      }
      cluster_gene_df = rbind(cluster_gene_df,
                              data.frame(gene_id=g,
                                         bin=this_gene_region,
                                         cn_thresh=cluster_threshs[j]))
    }
    #genes sharing a bin and threshold would otherwise double that bin's weight
    cluster_gene_df = cluster_gene_df[!duplicated(cluster_gene_df[,c("bin","cn_thresh")]),,drop=FALSE]
    if(is.null(cluster_gene_df)||nrow(cluster_gene_df)==0){
      stop("None of the genes in clusterRowsByGenes could be mapped to a bin")
    }
  }
  if(!missing(labelTheseGenes) | !missing(labelTheseCytobands) | !missing(sortByGenes)){

    if(!missing(sortByGenes)){

      message("mapping genes to bins")
      for(g in sortByGenes){
      gene_region = suppressMessages(gene_to_region(g))
      this_gene_region = map_bin_to_bin(gene_region)
      if(!missing(geneBoxPlot)){
        if(g == geneBoxPlot){
          splitByBinState = this_gene_region
          if(missing(keep_these_chromosomes)){
            keep_these_chromosomes = gsub(":.+","",this_gene_region)
          }
        }
      }
      #print(paste(g,gene_region))
      if (is.na(this_gene_region)){
        message(paste("no region for",g))
        next
      }
      #if (this_gene_region %in% names(bin_labels)){
      #  message(paste("Gene", g, "and gene", bin_labels[[this_gene_region]],
      #  "share a region. I will only show the last one!"))
      #}
      #bin_labels[[this_gene_region]]=g
      gene_bins[[g]] = this_gene_region
      }
    } 

    if(!missing(labelTheseGenes)){

      message("mapping genes to bins")
      for(g in labelTheseGenes){
      gene_region = suppressMessages(gene_to_region(g))
      this_gene_region = map_bin_to_bin(gene_region)
      if(!missing(geneBoxPlot)){
        if(g == geneBoxPlot){
          splitByBinState = this_gene_region
          if(missing(keep_these_chromosomes)){
            keep_these_chromosomes = gsub(":.+","",this_gene_region)
          }
        }
      }
      #print(paste(g,gene_region))
      if (is.na(this_gene_region)){
        message(paste("no region for",g))
        next
      }
      if (this_gene_region %in% names(bin_labels)){
        message(paste("Gene", g, "and gene", bin_labels[[this_gene_region]],
        "share a region. I will only show the last one!"))
      }
      bin_labels[[this_gene_region]]=g
      #gene_bins[[g]] = this_gene_region
      }
    } 
    if(!missing(labelTheseCytobands)){
      message("mapping cytobands to bins")
      if(genome_build=="grch37"){
        cytobands = cytobands_grch37 %>% 
          mutate(name=paste0(cb.chromosome,cb.name),
          region=paste0(cb.chromosome,":",cb.start,"-",cb.end))
      }else if(genome_build == "hg38"){
        cytobands = cytobands_hg38 %>%
          mutate(name=paste0(cb.chromosome,cb.name),
          region=paste0(cb.chromosome,":",cb.start,"-",cb.end))
      }
      for(g in labelTheseCytobands){
        cytoband_region = dplyr::filter(cytobands,name==g) %>% pull(region)
        this_cyto_region = map_bin_to_bin(cytoband_region)
      
        if (is.na(this_cyto_region)){
          message(paste("no region for",g))
          next
        }
        if (this_cyto_region %in% names(cytoband_labels)){
          message(paste("Cytoband", g, "and region", cytoband_labels[[this_cyto_region]],
          "share a region. I will only show the last one!"))
        }
        cytoband_labels[[this_cyto_region]]=g
      }
    }
  }else {
    if (!missing(geneBoxPlot)){
      gene_region = suppressMessages(gene_to_region(g))
      this_gene_region = map_bin_to_bin(gene_region)

      #print(paste(g,gene_region))
      if(is.na(this_gene_region)){
        message(paste("no region for",g))
      }
      if(this_gene_region %in% names(bin_labels)){
        message(paste("Gene",g, "and gene",bin_labels[[this_gene_region]], "share a region. I will only show the last one!"))
      }
      if(!missing(geneBoxPlot)){
        splitByBinState = this_gene_region
      }
      bin_labels[[this_gene_region]]=g
      if(missing(keep_these_chromosomes)){
        keep_these_chromosomes = gsub(":.+","",this_gene_region)
      }
    }
  }
  # label colours by bin; bins without an entry are drawn in black
  bin_label_colours = list()
  if(!missing(labelTheseFeatures)){
    required_cols = c("name", "chrom", "start", "end", "direction")
    if(!all(required_cols %in% colnames(labelTheseFeatures))){
      stop(paste("labelTheseFeatures must have columns:",
                 paste(required_cols, collapse = ", ")))
    }
    if(!feature_label %in% c("gene", "name", "region")){
      stop("feature_label must be \"gene\", \"name\" or \"region\"")
    }
    features = as.data.frame(labelTheseFeatures)
    blank = function(x){ is.null(x) || all(is.na(x)) }
    fallback = if(!blank(features$region_name)) features$region_name else features$name
    label_text = switch(
      feature_label,
      gene = if(blank(features$representative_gene)) fallback
             else ifelse(is.na(features$representative_gene) |
                           features$representative_gene == "",
                         fallback, features$representative_gene),
      name = features$name,
      region = fallback
    )
    # place each feature at the bin containing its midpoint
    bin_names = colnames(cn_state_matrix)
    bin_chrom = gsub("^chr", "", sub(":.*", "", bin_names))
    bin_start = as.numeric(sub(".*:([0-9]+)-.*", "\\1", bin_names))
    bin_end = as.numeric(sub(".*-([0-9]+)$", "\\1", bin_names))
    midpoint = (features$start + features$end) / 2
    feature_chrom = gsub("^chr", "", as.character(features$chrom))
    for(i in seq_len(nrow(features))){
      hit = which(bin_chrom == feature_chrom[i] & bin_start <= midpoint[i] &
                    bin_end >= midpoint[i])
      if(length(hit) == 0){
        message(paste("no bin for feature", features$name[i]))
        next
      }
      bin = bin_names[hit[1]]
      colour = if(colour_feature_labels){
        ifelse(features$direction[i] == "gain", "red", "blue")
      }else{
        "black"
      }
      if(bin %in% names(bin_labels)){
        # features sharing a bin are joined rather than overwritten, unless
        # the bin already shows this label (e.g. the same gene from
        # labelTheseGenes)
        if(!label_text[i] %in% strsplit(bin_labels[[bin]], "/")[[1]]){
          bin_labels[[bin]] = paste(bin_labels[[bin]], label_text[i], sep = "/")
        }
        if(!identical(bin_label_colours[[bin]], colour)){
          bin_label_colours[[bin]] = "black"
        }
      }else{
        bin_labels[[bin]] = label_text[i]
        bin_label_colours[[bin]] = colour
      }
    }
  }



  if(!missing(these_samples_metadata)){
    keep_samples = pull(these_samples_metadata, sample_id)
    all_samples = rownames(cn_state_matrix)
    overlap_samples = all_samples[all_samples %in% keep_samples]
    if(keep_sample_order){
      cn_state_matrix = cn_state_matrix[keep_samples,]
    }

    

  }
  if(!missing(numericMetadataColumns)){
  	not_numeric_cols = metadataColumns[!metadataColumns %in% numericMetadataColumns]
  }else{
	not_numeric_cols = metadataColumns
  }
  colours = map_metadata_to_colours(metadataColumns = not_numeric_cols,
                                    these_samples_metadata = these_samples_metadata,
                                    verbose=verbose)
  if(!missing(expressionColumns)){
    these_samples_metadata = these_samples_metadata %>%
      mutate(across(all_of(expressionColumns), ~ trim_scale_expression(.x)))
    exp_col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("blue","white", "red"))

    for(exp in expressionColumns){
      colours[[exp]] = exp_col_fun
    }
  }
  # caller-supplied mappings win, e.g. to keep colours identical across the
  # groups drawn by stacked_CN_heatmap
  for(anno_name in names(annotation_colours)){
    colours[[anno_name]] = annotation_colours[[anno_name]]
  }

  heatmap_args = list()
  cn_state_matrix[cn_state_matrix>5]=5
  col_fun = circlize::colorRamp2(c(1, 2, 5), c("blue", "white", "red"))
  chrom_col = c("chr1"="#555FAB",
                "chr2"="#CE3D31",
                "chr3"="#749B58",
                "chr4"="#F0E584",
                "chr5"="#476A85",
                "chr6"="#BA6338",
                "chr7"="#5CB1DD",
                "chr8"="#7F2368",
                "chr9"="#77C269",
                "chr10"="#D595A6",
                "chr11"="#934823",
                "chr12"="#857B8E",
                "chr13"="#C85328",
                "chr14"="#D58F5B",
                "chr15"="#7A65A7",
                "chr16"="#E3AF69",
                "chr17"="#3C1C54",
                "chr18"="#CEDEB7",
                "chr19"="#612B79",
                "chr20"="#AF2064",
                "chr21"="#E6C66F",
                "chr22"="#5B665D",
                "chrX"="#CA992C",
                "chrY"="white")

  #assume format is chrX:XXX-XXX
  column_chromosome = gsub(":.+", "", colnames(cn_state_matrix))


  if(!missing(hide_these_chromosomes)){
    keep_cols = which(!column_chromosome %in% hide_these_chromosomes)
    cn_state_matrix = cn_state_matrix[,keep_cols]
    column_chromosome = column_chromosome[keep_cols]

  }
  if(!grepl("chr",column_chromosome[1])){
    #remove chr prefix from colours
    names(chrom_col) = gsub("chr", "", names(chrom_col))
  }
  
  splits = NULL
  if(!missing(sortByFeatures)){
    if(missing(cnv_status)){
      stop("sortByFeatures needs cnv_status, e.g. from GAMBLR.utils::segmented_data_to_cnv_status()")
    }
    unknown = setdiff(sortByFeatures, colnames(cnv_status))
    if(length(unknown) > 0){
      stop(paste("sortByFeatures not found in cnv_status:",
                 paste(unknown, collapse = ", ")))
    }
    # ties keep the clustering order (or the input order if not clustering)
    tie_order = seq_len(nrow(cn_state_matrix))
    if(isTRUE(cluster_rows)){
      tie_order = hclust(dist(neutral_for_clustering(cn_state_matrix)))$order
    }
    status = as.data.frame(cnv_status)[rownames(cn_state_matrix), sortByFeatures,
                                        drop = FALSE]
    # rank by the strength of the event in its feature's direction; samples
    # with no status count as unaltered
    status[] = lapply(status, function(x){
      x = abs(as.numeric(x))
      x[is.na(x)] = 0
      x
    })
    keys = c(lapply(status, function(x) -x[tie_order]), list(seq_along(tie_order)))
    new_order = tie_order[do.call(order, keys)]
    cn_state_matrix = cn_state_matrix[new_order, , drop = FALSE]
    cluster_rows = FALSE
  }else if(!missing(sortByGenes)){
    # the sort only survives if rows aren't clustered afterwards
    cluster_rows = FALSE
    sortByBins = unlist(gene_bins)
    cn_state_matrix = arrange(cn_state_matrix,desc(across(sortByBins)))

  } else if(!missing(sortByBins)){
    cluster_rows = FALSE
    cn_state_matrix = arrange(cn_state_matrix,desc(across(sortByBins)))
  }else if(!missing(sortByMetadataColumns)){
    cluster_rows = FALSE
    new_order = arrange(these_samples_metadata,across(sortByMetadataColumns)) %>% pull(sample_id)
    cn_state_matrix = cn_state_matrix[new_order,]

  }

    gain_state = cn_state_matrix
    gain_state[gain_state<=2]=0
    gain_state[gain_state>0]=1
    total_gain = colSums(gain_state, na.rm = TRUE)
    loss_state = cn_state_matrix
    loss_state[loss_state>=2]=0
    loss_state[loss_state!=0]=1
    total_loss = colSums(loss_state, na.rm = TRUE)
    if(!missing(expressionColumns)){
      metadataColumns = c(metadataColumns,expressionColumns)
    }
    row_df = select(these_samples_metadata,sample_id,all_of(unique(c(metadataColumns)))) %>%
      column_to_rownames("sample_id")
    either_state = loss_state
    either_state[gain_state==1]=1
    sample_average = rowMeans(either_state, na.rm = TRUE)
    samples_keep_PGA = which(sample_average >= drop_if_PGA_below & sample_average <= drop_if_PGA_above)

    cn_state_matrix = cn_state_matrix[samples_keep_PGA,]

    sample_average = sample_average[samples_keep_PGA]

    #REDO
    gain_state = cn_state_matrix
    gain_state[gain_state<=2]=0
    gain_state[gain_state>0]=1
    total_gain = colSums(gain_state, na.rm = TRUE)
    loss_state = cn_state_matrix
    loss_state[loss_state>=2]=0
    loss_state[loss_state!=0]=1
    total_loss = colSums(loss_state, na.rm = TRUE)
    #sample_CN_anno = HeatmapAnnotation(PGA=sample_average,which="row")


    if(!is.null(cluster_gene_df)){
      #snapshot before chromosome/bin restriction so genes outside the
      #displayed chromosomes can still drive the sample clustering
      bin_present = cluster_gene_df$bin %in% colnames(cn_state_matrix)
      if(any(!bin_present)){
        message(paste("these clustering bins were removed by other filters:",
                      paste(unique(cluster_gene_df$bin[!bin_present]),collapse=", ")))
      }
      cluster_gene_df = cluster_gene_df[bin_present,,drop=FALSE]
      if(nrow(cluster_gene_df)==0){
        stop("None of the bins for clusterRowsByGenes remain after filtering")
      }
      #one column per gene so a bin shared by genes with opposing
      #thresholds contributes each relevant direction separately
      cluster_gene_mat = as.matrix(cn_state_matrix[,cluster_gene_df$bin,drop=FALSE])
      colnames(cluster_gene_mat) = cluster_gene_df$gene_id
      for(j in seq_len(nrow(cluster_gene_df))){
        thresh = cluster_gene_df$cn_thresh[j]
        vals = cluster_gene_mat[,j]
        if(thresh < 2){
          vals[!is.na(vals) & vals > thresh] = 2
        }else if(thresh > 2){
          vals[!is.na(vals) & vals < thresh] = 2
        }
        cluster_gene_mat[,j] = vals
      }
      if(cluster_rows){
        #A function (rather than a precomputed hclust) keeps this compatible
        #with row_split, since ComplexHeatmap calls it on each slice
        cluster_rows = function(m){
          these_ids = rownames(m)
          if(!all(these_ids %in% rownames(cluster_gene_mat))){
            these_ids = colnames(m)
          }
          hclust(dist(neutral_for_clustering(cluster_gene_mat[these_ids,,drop=FALSE])))
        }
      }else{
        message("clusterRowsByGenes is ignored because row clustering is disabled")
      }
    }

    if(!missing(keep_these_chromosomes)){
      keep_cols = which(column_chromosome %in% keep_these_chromosomes)
      cn_state_matrix = cn_state_matrix[,keep_cols]
      column_chromosome = column_chromosome[keep_cols]
      total_gain = total_gain[keep_cols]
      total_loss = total_loss[keep_cols]
    }

    #erase (set to diploid) any value opposing the most common event in that bin
    if(flip){
    flipped_data = cn_state_matrix
    for(i in c(1:length(total_gain))){
      bin_name = names(total_gain)[i]
      if(total_gain[i] > total_loss[i]){
        flipped_data[cn_state_matrix[,i]<2,i]=2
      }else{
        flipped_data[cn_state_matrix[,i]>2,i]=2
      }
    }
      cn_state_matrix = flipped_data
    }



    if(!missing(keep_these_bins)){
      available_bins = keep_these_bins[which(keep_these_bins %in% colnames(cn_state_matrix))]
      cn_state_matrix = cn_state_matrix[,available_bins]
      regions_highlight = regions_highlight[which(regions_highlight %in% keep_these_bins)]
      column_chromosome = gsub(":.+", "", colnames(cn_state_matrix))
      total_gain = total_gain[available_bins]
      total_loss = total_loss[available_bins]
    }



    anno_rows = row_df[rownames(cn_state_matrix),,drop=FALSE]

    anno_rows$PGA = sample_average

    if(sortByPGA){
      anno_rows = arrange(anno_rows,PGA)
      cn_state_matrix = cn_state_matrix[rownames(anno_rows),]
    }

    if(!missing(splitByBinState)){
      splits = round(cn_state_matrix[,splitByBinState])
      anno_rows$CN_state = splits


      if(!missing(geneBoxPlot)){
        rg = c(0,1)
        panel_fun = function(index, nm) {
          if(boxplot_orientation == "vertical"){
            pushViewport(viewport(yscale = rg, xscale = c(0, 2)))
            grid.rect()
            grid.yaxis(gp = gpar(fontsize = 5))
            grid.boxplot(anno_rows[index,geneBoxPlot],
                         gp = gpar(fill = colours[["CN"]][as.character(splits[index][1])]),
                         pos = 1,
                         direction = "vertical"
                         ,outline = T)
          }else{
            pushViewport(viewport(xscale = rg, yscale = c(0, 2)))
            grid.rect()
            grid.xaxis(gp = gpar(fontsize = 5))
            grid.boxplot(anno_rows[index,geneBoxPlot],
                         gp = gpar(fill = colours[["CN"]][as.character(splits[index][1])]),
                         pos = 1,
                         direction = "horizontal"
                         ,outline = T)
          }

          popViewport()
        }
        anno = anno_link(align_to = splits,
                         which = "row",
                         panel_fun = panel_fun,
                         size = unit(2, "cm"),
                         gap = unit(1, "cm"),
                         width = unit(3, "cm"))
      }
      coldf = get_gambl_colours("copy_number",as_dataframe = T) %>%
        filter(name %in% splits)
      colvec = pull(coldf,colour)
      names(colvec)=pull(coldf,name)
      colours[["CN"]] = colvec
      colours[["CN_state"]] = colvec


      if(!missing(hide_annotations)){
        show_legend = rep(TRUE, length(colnames(anno_rows)))
        names(show_legend) = colnames(anno_rows)
        show_legend[hide_annotations] = FALSE
      }else{
        show_legend = rep(TRUE, length(colnames(anno_rows)))
      }

      left_anno = HeatmapAnnotation(df=anno_rows,
                                    col=colours,
                                    simple_anno_size = unit(metadataBarHeight, "mm"),
                                    which="row",
                                    show_legend=show_legend,
                                    annotation_name_side = left_annotation_name_side,
                                    annotation_name_gp = gpar(fontsize = metadataBarFontsize),
                                    annotation_legend_param =
                                      list(
                                           by_row=F,
                                           nrow=legend_row,
                                           ncol = legend_col,
                                           direction=legend_direction))

    }else{ #no splitByBinState

      if(!missing(hide_annotations)){
        show_legend = rep(TRUE, length(colnames(anno_rows)))
        names(show_legend) = colnames(anno_rows)
        show_legend[hide_annotations] = FALSE
      }else{
        show_legend = rep(TRUE, length(colnames(anno_rows)))
      }
      if(rotate){
        bottom_annotation = HeatmapAnnotation(df=anno_rows,
                                      col=colours,
                                      simple_anno_size = unit(metadataBarHeight, "mm"),
                                      which="column",
                                      show_legend=show_legend,
                                      annotation_name_gp = gpar(fontsize = metadataBarFontsize),
                                      annotation_name_side = "left",
                                      annotation_legend_param =
                                        list(
                                             by_row=F,
                                             nrow=legend_row,
                                             ncol = legend_col,
                                             direction=legend_direction))
      }else{
        left_anno = HeatmapAnnotation(df=anno_rows,
                                      col=colours,
                                      simple_anno_size = unit(metadataBarHeight, "mm"),
                                      which="row",
                                      show_legend=show_legend,
                                      annotation_name_side = left_annotation_name_side,
                                      annotation_name_gp = gpar(fontsize = metadataBarFontsize),
                                      annotation_legend_param =
                                        list(
                                             by_row=F,
                                             nrow=legend_row,
                                             ncol = legend_col,
                                             direction=legend_direction))
      }
    }

    if(!missing(hide_annotations)){
      if("chromosome" %in% hide_annotations){
        show_legend = FALSE
      }else{
        show_legend = TRUE
      }
    }else{
      show_legend=TRUE
    }
    bin_average = colMeans(cn_state_matrix, na.rm = TRUE)
    cn_av_df = tibble(x=bin_average,coordinate=names(bin_average))
    cn_av_df = mutate(cn_av_df,
                      local.minima=ifelse(lag(x,n=15)>x & lead(x,n=15)  > x & x < 2,TRUE,FALSE),
                      local.maxima=ifelse(lag(x,n=15)<x & lead(x,n=15)< x & x > 2,TRUE,FALSE),
                      extreme=ifelse(local.maxima | local.minima,TRUE,FALSE))
    if(rotate){
      average_anno = HeatmapAnnotation(Mean_CN = bin_average,
                                     which = "row",
                                     show_annotation_name = show_bottom_annotation_name,
                                     #annotation_name_side = "left",
                                     show_legend=show_legend,
                                     col=list(Mean_CN=col_fun))
    }else{
      average_anno = HeatmapAnnotation(Mean_CN = bin_average,
                                     which = "column",
                                     show_annotation_name = show_bottom_annotation_name,
                                     #annotation_name_side = "left",
                                     show_legend=show_legend,
                                     col=list(Mean_CN=col_fun))
    }

    # "mark" labels are an annotation laid out with the heatmap, so they are
    # aligned, kept apart and given their own space (rows when rotated)
    label_mark_anno = NULL
    if(label_style == "mark" && length(bin_labels) > 0){
      label_bins = names(bin_labels)
      label_at = match(label_bins, colnames(cn_state_matrix))
      present = !is.na(label_at)
      if(any(present)){
        label_cols = vapply(label_bins[present], function(b){
          if(is.null(bin_label_colours[[b]])) "black" else bin_label_colours[[b]]
        }, character(1))
        label_text = unname(unlist(bin_labels[present]))
        label_at = label_at[present]
        # anno_mark orders labels by position but applies gpar values in the
        # order given, so put everything in position order first
        o = order(label_at)
        label_at = label_at[o]
        label_text = label_text[o]
        label_cols = unname(label_cols)[o]
        label_gp = gpar(fontsize = bin_label_fontsize, col = label_cols)
        if(rotate){
          label_mark_anno = rowAnnotation(
            labels = anno_mark(at = label_at, labels = label_text,
                               which = "row", side = "left",
                               labels_gp = label_gp)
          )
        }else{
          # column anno_mark applies per-label gpar values in reverse order
          # (seen in ComplexHeatmap 2.23.1; the row version is unaffected)
          label_mark_anno = HeatmapAnnotation(
            labels = anno_mark(at = label_at, labels = label_text,
                               which = "column", side = "top",
                               labels_rot = bin_label_rotation,
                               labels_gp = gpar(fontsize = bin_label_fontsize,
                                                col = rev(label_cols))),
            which = "column",
            show_annotation_name = FALSE
          )
        }
      }
    }
    

    # bar heights: counts, or the fraction of samples shown, so groups of
    # different sizes are comparable
    gain_bars = total_gain
    loss_bars = total_loss
    if(cumulative_as_fraction){
      # per bin, out of the samples with data there (NA = no segments)
      n_with_data = pmax(colSums(!is.na(cn_state_matrix)), 1)
      gain_bars = total_gain / n_with_data
      loss_bars = total_loss / n_with_data
    }
    # "minimal" axes label only 0 and the top of the scale, which stays
    # legible when the tracks are short
    show_bar_axis = cumulative_axis != "none"
    bar_axis_param = function(base, bars){
      if(cumulative_axis != "minimal"){
        return(base)
      }
      top = if(!is.null(cumulative_ylim)) cumulative_ylim[2] else max(bars, na.rm = TRUE)
      if(!is.finite(top) || top <= 0){
        return(base)
      }
      c(base, list(at = c(0, top), labels = c("0", format(signif(top, 2)))))
    }
    #need to create a data frame with the bin names as rownames and a different value depending on whether
    # the bin is in regions_highlight
    if(rotate){
      cumulative_anno = HeatmapAnnotation(chromosome=column_chromosome,
                                        Gain = anno_barplot(gain_bars,
                                                            ylim = cumulative_ylim,
                                                            width = unit(cumulative_track_height, "cm"),
                                                            gp=gpar(fill="red",col="red"),
                                                            axis = show_bar_axis,
                                                            axis_param = bar_axis_param(list(direction = "reverse",side="top"), gain_bars)),
                                        loss=anno_barplot(loss_bars,
                                                          ylim = cumulative_ylim,
                                                          width = unit(cumulative_track_height, "cm"),
                                                          gp=gpar(fill="blue",col="blue"),
                                                          axis = show_bar_axis,
                                                          axis_param = bar_axis_param(list(side = "top"), loss_bars)),
                                        
                                                          #),
                                        annotation_legend_param =
                                          list(title = "Chromosome",
                                               nrow=legend_row,
                                               ncol = legend_col,
                                               direction=legend_direction,
                                               by_row=F),
                                 which = "row",
                                 show_annotation_name = show_bottom_annotation_name,
                                 #annotation_name_side = "right",
                                 annotation_name_gp = gpar(fontsize = metadataBarFontsize),
                                 show_legend=show_legend,
                                 col=list(chromosome=chrom_col))

    }else{
      cumulative_anno = HeatmapAnnotation(chromosome=column_chromosome,
                                        #annotation_name_gp = gpar(fontsize = metadataBarFontsize),
                                        Gain = anno_barplot(gain_bars,
                                                            ylim = cumulative_ylim,
                                                            height = unit(cumulative_track_height, "cm"),
                                                            gp=gpar(fill="red",col="red"),
                                                            axis = show_bar_axis,
                                                            axis_param = bar_axis_param(list(), gain_bars)),
                                        loss=anno_barplot(loss_bars,
                                                          ylim = cumulative_ylim,
                                                          height = unit(cumulative_track_height, "cm"),
                                                          gp=gpar(fill="blue",col="blue"),
                                                          axis = show_bar_axis,
                                                          axis_param = bar_axis_param(list(direction = "reverse"), loss_bars)),
                                        annotation_legend_param =
                                          list(title = "Chromosome",
                                               nrow=legend_row,
                                               ncol = legend_col,
                                               direction=legend_direction,
                                               by_row=F),
                                 which = "column",
                                 show_annotation_name = show_bottom_annotation_name,
                                 #annotation_name_side = "right",
                                 annotation_name_gp = gpar(fontsize = metadataBarFontsize),
                                 show_legend=show_legend,
                                 col=list(chromosome=chrom_col))

    }
    
    bars_on_top = !rotate && cumulative_position == "top" && !hide_cumulative_plot
    if(hide_cumulative_plot || bars_on_top){
      cumulative_anno = HeatmapAnnotation(chromosome=column_chromosome,
        annotation_legend_param = list(title = "Chromosome",
                                               nrow=legend_row,
                                               ncol = legend_col,
                                               direction=legend_direction,
                                               by_row=F),
                                 which = ifelse(rotate,"row","column"),
                                 show_annotation_name = show_bottom_annotation_name,
                                 #annotation_name_side = "right",
                                 annotation_name_gp = gpar(fontsize = metadataBarFontsize),
                                 show_legend=show_legend,
                                 col=list(chromosome=chrom_col))

    }
    if(!show_chromosome_strip){
      remaining = setdiff(names(cumulative_anno@anno_list), "chromosome")
      cumulative_anno = if(length(remaining) > 0) cumulative_anno[, remaining] else NULL
    }
    # top annotation, top to bottom: labels, gain/loss tracks (when placed on
    # top, so labels point at the peaks), then the mean CN strip
    top_parts = list()
    if(!rotate){
      if(!is.null(label_mark_anno)){
        top_parts = c(top_parts, list(label_mark_anno))
      }
      if(bars_on_top){
        top_parts = c(top_parts, list(HeatmapAnnotation(
          Gain = anno_barplot(gain_bars,
                              ylim = cumulative_ylim,
                              height = unit(cumulative_track_height, "cm"),
                              gp = gpar(fill = "red", col = "red"),
                              axis = show_bar_axis,
                              # right side: the left is taken by the sample
                              # annotation names drawn above the heatmap
                              axis_param = bar_axis_param(list(side = "right"), gain_bars)),
          loss = anno_barplot(loss_bars,
                              ylim = cumulative_ylim,
                              height = unit(cumulative_track_height, "cm"),
                              gp = gpar(fill = "blue", col = "blue"),
                              axis = show_bar_axis,
                              axis_param = bar_axis_param(list(direction = "reverse", side = "right"), loss_bars)),
          which = "column",
          show_annotation_name = show_bottom_annotation_name,
          annotation_name_gp = gpar(fontsize = metadataBarFontsize)
        )))
      }
      if(show_mean_cn){
        top_parts = c(top_parts, list(average_anno))
      }
    }
    top_anno = if(length(top_parts) > 0) Reduce(c, top_parts) else NULL
    if(anyNA(cn_state_matrix)){
      na_safe_clustering = function(m) hclust(dist(neutral_for_clustering(m)))
      if(isTRUE(cluster_rows)){
        cluster_rows = na_safe_clustering
      }
      if(isTRUE(cluster_columns)){
        cluster_columns = na_safe_clustering
      }
    }
    heatmap_legend_param = list(title = "Copy Number",

                                by_row=F,
                                legend_direction = legend_direction)
    #bottom_anno = HeatmapAnnotation()
    if(!missing(geneBoxPlot)){
      ho = Heatmap(cn_state_matrix,
                   name=heatmap_name,
                   column_title=" ",
                   row_title = NULL,
                   cluster_columns=cluster_columns,
                   cluster_rows=cluster_rows,
                   show_row_names = show_row_names,
                   show_column_names = show_column_names,
                   col = col_fun,
                   bottom_annotation = cumulative_anno,
                   top_annotation = top_anno,
                   left_annotation=left_anno,
                   row_split=splits,
                   right_annotation = rowAnnotation(Expression=anno),
                   #column_split  = column_chromosome,
                   heatmap_legend_param = heatmap_legend_param
      )
    }else{
      if(rotate){
        if(!is.null(width)){
          w = unit(width, "cm")
        }else{
          w = NULL
        }

        if(!is.null(height)){
          h = unit(height, "cm")
        }else{
          h = NULL
        }
        ho = Heatmap(t(cn_state_matrix),
                   name=heatmap_name,
                   column_title=" ",
                   cluster_columns=cluster_rows,
                   cluster_rows=cluster_columns,
                   show_row_names = show_row_names,
                   show_column_names = show_column_names,
                   col = col_fun,
                   bottom_annotation = bottom_annotation,
                   #top_annotation = average_anno,
                   left_annotation = label_mark_anno,
                   #row_split=splits,
                   width = w,
                   height = h,
                   right_annotation = cumulative_anno,
                   #column_split  = column_chromosome,
                   heatmap_legend_param = heatmap_legend_param)

      }else{

        ho = Heatmap(cn_state_matrix,name=heatmap_name,
                   column_title=" ",
                   cluster_columns=cluster_columns,
                   cluster_rows=cluster_rows,
                   show_row_names = show_row_names,
                   show_column_names = show_column_names,
                   col = col_fun,
                   bottom_annotation = cumulative_anno,
                   top_annotation = top_anno,
                   left_annotation=left_anno,
                   row_split=splits,
                   row_title = heatmap_row_title,
                   row_title_gp = gpar(fontsize = heatmap_row_title_fontsize),
                   row_title_rot = heatmap_row_title_rot,
                   #right_annotation = sample_CN_anno,
                   #column_split  = column_chromosome,
                   heatmap_legend_param = heatmap_legend_param)
      }
  }
  if(!return_data){
    draw(ho,heatmap_legend_side=legend_position,annotation_legend_side=legend_position,
         show_heatmap_legend = show_legends,
         show_annotation_legend = show_legends)
  }
  # labels decorate the drawn heatmap, so skip them when it isn't drawn
  if(!return_data &&
     (!missing(labelTheseGenes)|!missing(labelTheseCytobands)|length(bin_labels)>0)){
    if(length(bin_labels) > 0 && label_style == "classic"){
      if(verbose){
        print("gene labeling")
      }
      tmat = t(cn_state_matrix)
      for(i in c(1:length(bin_labels))){
          gene_region = names(bin_labels)[i]
          gene_name = unname(bin_labels[[i]])
          label_colour = bin_label_colours[[gene_region]]
          if(is.null(label_colour)){
            label_colour = "black"
          }
          # the bin may have been removed by filtering after labels were mapped
          if(!gene_region %in% colnames(cn_state_matrix)){
            next
          }
          region_chunks = region_to_chunks(gene_region)
          if(region_chunks$chromosome %in% column_chromosome){
            if(rotate){
              decorate_heatmap_body(heatmap_name,
              {
                i=which(rownames(tmat)==gene_region)
                y=1-i/nrow(tmat)
                grid.text(gene_name,y=y,gp=gpar(fontsize=bin_label_fontsize,
                                                col=label_colour),
                  unit(bin_label_nudge,"npc"),
                      rot=bin_label_rotation,just="top")
                grid.lines(c(1.01, 1), c(y,y), gp = gpar(lwd = 1))
              })
              }else{
              decorate_heatmap_body(heatmap_name,
              {
                i=which(colnames(cn_state_matrix)==gene_region)
                x=i/ncol(cn_state_matrix)
                grid.text(gene_name,x,gp=gpar(fontsize=bin_label_fontsize,
                                              col=label_colour),
                      unit(bin_label_nudge,"npc"),
                      rot=bin_label_rotation,just="top")
                grid.lines(c(x, x), c(1.01, 1), gp = gpar(lwd = 1))
              })
              }
          }
      }

    }
    if(verbose){
      print(cytoband_labels)
    }
    if(!missing(labelTheseCytobands)){
      for(i in c(1:length(cytoband_labels))){
        gene_region = names(cytoband_labels)[i]
        gene_name = unname(cytoband_labels[[i]])
        region_chunks = region_to_chunks(gene_region)
        if(region_chunks$chromosome %in% column_chromosome){          
            if(verbose){
              print("cytoband labeling")
            }
            decorate_annotation("Gain",
                              {i=which(colnames(cn_state_matrix)==gene_region)
                               x=i/ncol(cn_state_matrix)
                               x2 = (i+1)/ncol(cn_state_matrix)
                               xsize = x2 - x
                               xmid = x + xsize/2
                               grid.text(gene_name,x,1.04,gp=gpar(fontsize=bin_label_fontsize),
                                   unit(bin_label_nudge,"npc"),
                                   rot=bin_label_rotation,just="top")
                               grid.lines(c(xmid, xmid), c(-1, 1), gp = gpar(lwd = 1,
                                                                        col = "#00FF0040"))
                              })
          }
        }
    }

  }

  if(return_data){
    if(length(bin_labels)==0){
      bin_labels=NULL
    }
    return(list(heatmap_object=ho,
                data=cn_state_matrix,
                cumulative_gain = total_gain,
                cumulative_loss = total_loss,
                labels=bin_labels,
                chromosome_columns=column_chromosome,
                bin_means=bin_average,
                local_optima=cn_av_df,
                row_anno = anno_rows))
  }

}
