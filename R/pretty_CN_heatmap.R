


#' Pretty Copy Number Heatmap
#'
#' @param cn_state_matrix The output of get_cn_states
#' @param scale_by_sample Set to TRUE to scale CN values within each sample_id
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
#' @param sortByBins Optional: A vector containing one or more names of genomic bins that will be used to order the heatmap rows.
#' @param sortByPGA Optional: Sort the rows based on percent genome altered (PGA) instead of the other options
#' @param splitByBinState Optional: A single genomic bin that will be used to split the heatmap based on the CN state of that bin
#' @param sortByMetadataColumns A vector containing one or more names of columns from your metadata that will be used to order the rows overall or within slices (if combined with splitByBinState or geneBoxPlot)
#' @param labelTheseGenes A vector of Hugo gene symbols whose location will be indicated on the top of the heatmap
#' @param bin_label_fontsize Font size for the gene labels (default 5)
#' @param bin_label_nudge Increase or decrease this value to shift the gene labels up/down (default 1.03)
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
#' @param drop_bin_if_sd_below Force bins with standard deviation below this value to be excluded
#' @param return_data Specify TRUE to get some of the internal data back including the heatmap object
#' @param rotate Optionally, flip the rows/columns of resulting heatmap. Default is FALSE.
#' @param width Set the width of the heatmap. Default is 10.
#' @param verbose Control verbosity of the console output. Default is FALSE.
#'
#' @return list (when return_data = TRUE)
#' @import GAMBLR.helpers dplyr
#' @export
#'
#' @examples
#' cat("Running example for function: pretty_CN_heatmap\n")
#' suppressMessages(library(dplyr))
#' suppressMessages(library(GAMBLR.open))
#' #get some metadata for subsetting the data to just one pathology (DLBCL)
#' dlbcl_genome_meta = suppressMessages(get_gambl_metadata()) %>%
#'                     filter(pathology=="DLBCL",
#'                     seq_type=="genome")
#' 
#' #remove any duplicate sample_id/seq_type combinations
#' meta_clean = check_and_clean_metadata(dlbcl_genome_meta,
#'                                       duplicate_action = "keep_first")
#' 
#' # Create the copy number matrix using the helper functions
#' all_segments = get_cn_segments(these = meta_clean)
#' dlbcl_cn_binned = segmented_data_to_cn_matrix(
#'                                   seg_data = all_segments,
#'                                   strategy="auto_split",
#'                                   n_bins_split=1300,
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
                             show_row_names=FALSE,
                             show_column_names=FALSE,
                             keep_these_chromosomes,
                             hide_these_chromosomes,
                             keep_these_bins,
                             hide_annotations,
                             sortByBins,
                             sortByGenes,
                             splitByBinState,
                             sortByPGA=FALSE,
                             sortByMetadataColumns,
                             labelTheseGenes,
                             labelTheseCytobands,
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
                             drop_bin_if_sd_below=0,
                             flip=FALSE,
                             max_CN_allowed = 6,
                             verbose=FALSE,
                             rotate = FALSE,
                             width = 15,
                             height = 6,
                             cluster_samples,
                             cluster_regions){
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
  cn_state_matrix = cn_state_matrix[,which(bin_vars > drop_bin_if_sd_below)]
  new_col = ncol(cn_state_matrix)
  if(verbose){
    print(paste("original:",old_col,"new:",new_col))
  }
  if(scale_by_sample){
    cn_state_matrix = cn_state_matrix - rowMeans(cn_state_matrix, na.rm = TRUE) + 2
    cn_state_matrix = round(cn_state_matrix)
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
  if(!missing(sortByGenes)){
    cluster_columns = FALSE
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
    

    #need to create a data frame with the bin names as rownames and a different value depending on whether
    # the bin is in regions_highlight
    if(rotate){
      cumulative_anno = HeatmapAnnotation(chromosome=column_chromosome,
                                        Gain = anno_barplot(total_gain,
                                                            gp=gpar(fill="red",col="red"),
                                                            axis_param = list(direction = "reverse",side="top")),
                                        loss=anno_barplot(total_loss,
                                                          gp=gpar(fill="blue",col="blue"),
                                                          axis_param = list(side = "top")),
                                        
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
                                        Gain = anno_barplot(total_gain,
                                                            gp=gpar(fill="red",col="red")),
                                        loss=anno_barplot(total_loss,
                                                          gp=gpar(fill="blue",col="blue"),
                                                          axis_param = list(direction = "reverse")),
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
        heatmap_legend_param = list(title = "Copy Number",

                                by_row=F,
                                legend_direction = legend_direction)
    #bottom_anno = HeatmapAnnotation()
    if(!missing(geneBoxPlot)){
      ho = Heatmap(cn_state_matrix,
                   name="CN",
                   column_title=" ",
                   row_title = NULL,
                   cluster_columns=cluster_columns,
                   cluster_rows=cluster_rows,
                   show_row_names = show_row_names,
                   show_column_names = show_column_names,
                   col = col_fun,
                   bottom_annotation = cumulative_anno,
                   top_annotation = average_anno,
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
                   name="CN",
                   column_title=" ",
                   cluster_columns=cluster_rows,
                   cluster_rows=cluster_columns,
                   show_row_names = show_row_names,
                   show_column_names = show_column_names,
                   col = col_fun,
                   bottom_annotation = bottom_annotation,
                   #top_annotation = average_anno,
                   #left_annotation=left_anno,
                   #row_split=splits,
                   width = w,
                   height = h,
                   right_annotation = cumulative_anno,
                   #column_split  = column_chromosome,
                   heatmap_legend_param = heatmap_legend_param)

      }else{

        ho = Heatmap(cn_state_matrix,name="CN",
                   column_title=" ",
                   cluster_columns=cluster_columns,
                   cluster_rows=cluster_rows,
                   show_row_names = show_row_names,
                   show_column_names = show_column_names,
                   col = col_fun,
                   bottom_annotation = cumulative_anno,
                   top_annotation = average_anno,
                   left_annotation=left_anno,
                   row_split=splits,
                   #right_annotation = sample_CN_anno,
                   #column_split  = column_chromosome,
                   heatmap_legend_param = heatmap_legend_param)
      }
  }
  if(!return_data){
    draw(ho,heatmap_legend_side=legend_position,annotation_legend_side=legend_position)
  }
  if(!missing(labelTheseGenes)|!missing(labelTheseCytobands)|length(bin_labels)>0){
    if(!missing(labelTheseGenes)){
      if(verbose){
        print("gene labeling")
      }
      tmat = t(cn_state_matrix)
      for(i in c(1:length(bin_labels))){
          gene_region = names(bin_labels)[i]
          gene_name = unname(bin_labels[[i]])
          region_chunks = region_to_chunks(gene_region)
          if(region_chunks$chromosome %in% column_chromosome){           
            if(rotate){
              decorate_heatmap_body("CN",
              {
                i=which(rownames(tmat)==gene_region)
                y=1-i/nrow(tmat)
                grid.text(gene_name,y=y,gp=gpar(fontsize=bin_label_fontsize),
                  unit(bin_label_nudge,"npc"),
                      rot=bin_label_rotation,just="top")
                grid.lines(c(1.01, 1), c(y,y), gp = gpar(lwd = 1))
              })
              }else{
              decorate_heatmap_body("CN",
              {
                i=which(colnames(cn_state_matrix)==gene_region)
                x=i/ncol(cn_state_matrix)
                grid.text(gene_name,x,gp=gpar(fontsize=bin_label_fontsize),
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
