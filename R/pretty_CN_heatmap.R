

#' Pretty Copy Number Heatmap
#'
#' @param cn_state_matrix The output of get_cn_states
#' @param scale_by_sample Set to TRUE to scale CN values within each sample_id
#' @param these_samples_metadata The output of get_gambl_metadata
#' @param metadataColumns One or more columns from the metadata you want to display beside the heatmap
#' @param expressionColumns Optional: One or more columns from the metadata that include gene expression values you want shown
#' @param geneBoxPlot Optional: Specify the Hugo symbol of a single gene to embed box plots adjacent to the heatmap. Expression data for this gene must be present in the metadata in a column of the same name. 
#' @param show_column_names Set to TRUE to display the sample_id of every sample shown in the heatmap
#' @param show_row_names Set to TRUE to display the ID of every bin (region) shown in the heatmap
#' @param keep_these_chromosomes A vector of chromosome names to include (all others will be excluded)
#' @param hide_these_chromosomes A vector of chromosome names to exclude (all others will be included unless keep_these_chromosomes is specified)
#' @param keep_these_bins A vector of bin names to include (all others will be excluded)
#' @param hide_annotations A vector of annotation names to suppress from legends in the plot
#' @param cluster_columns Set to TRUE to enable clustering of genomic regions (columns) based on their CN value across all patients in the heatmap
#' @param cluster_rows Set to TRUE to enable clustering of genomic regions (columns) based on their CN value across all regions in the heatmap
#' @param sortByBins Optional: A vector containing one or more names of genomic bins that will be used to order the heatmap rows.
#' @param splitByBinState Optional: A single genomic bin that will be used to split the heatmap based on the CN state of that bin
#' @param sortByMetadataColumns A vector containing one or more names of columns from your metadata that will be used to order the rows overall or within slices (if combined with splitByBinState or geneBoxPlot)
#' @param labelTheseGenes A vector of Hugo gene symbols whose location will be indicated on the top of the heatmap
#' @param bin_label_fontsize Font size for the gene labels (default 5)
#' @param bin_label_nudge Increase or decrease this value to shift the gene labels up/down (default 1.03)
#' @param drop_if_PGA_below Lower limit for proportion of genome altered (PGA). Samples below this value will be dropped (default 0)
#' @param drop_if_PGA_above Upper limit for proportion of genome altered (PGA). Samples above this value will be dropped (default 1)
#' @param show_bottom_annotation_name set to TRUE to label the bottom annotation tracks with more details
#' @param bottom_annotation_name_side If using show_bottom_annotation_name, set this to "left" or "right" to relocate the names
#' @param bin_labels Instead of automatically labeling genes, you can instead explicitly provide a list of labels for any bins in the heatmap. The names of each element should match a bin. The value of each element is the label that will be shown. This option can be used to skip gene location look-ups (see examples). 
#' @param legend_direction Which orientation to use for the legend
#' @param legend_position Where to put the legend
#' @param legend_row How many rows for the legend layout
#' @param legend_col How many columns for the legend layout
#' @param boxplot_orientation Either "horizontal" or "vertical" (default: horizontal)
#' @param left_annotation_name_side Which side to put the name of the metadata annotations (top or bottom)
#' @param drop_bin_if_sd_below Force bins with standard deviation below this value to be excluded
#' @param return_data Specify TRUE to get some of the internal data back including the heatmap object
#'
#' @return list (when return_data = TRUE)
#' @export
#'
#' @examples
#' 
#' # Create the copy number matrix using the helper functions
#' all_segments = get_cn_segments()
#' all_states_binned = get_cn_states(n_bins_split=2500,
#'                                  missing_data_as_diploid = T,
#'                                  seg_data = all_segments)
#'
#' #get some metadata for subsetting the data to just one pathology (DLBCL)
#' dlbcl_genome_meta = get_gambl_metadata() %>% 
#'                     filter(pathology=="DLBCL",
#'                     seq_type=="genome")
#'
#' # Generate a basic genome-wide CN heatmap
#' pretty_CN_heatmap(cn_state_matrix=all_states_binned,
#'                   these_samples_metadata = dlbcl_genome_meta,
#'                   hide_annotations = "chromosome")
#' 
#' # Disable row (sample) clustering and restrict to a few chromosomes 
#' # and highlight some genes of interest
#' pretty_CN_heatmap(cn_state_matrix=all_states_binned,
#'   these_samples_metadata = dlbcl_genome_meta,
#'   hide_annotations = "chromosomes",
#'   keep_these_chromosomes = c("chr9","chr17"),
#'   cluster_rows=F,
#'   labelTheseGenes = c("CDKN2A","TP53"))
#' 
#' 
#' # get gene expression data
#' gene_exp_all = get_gene_expression(all_genes=T,lazy_join=T,arbitrarily_pick = T,HGNC=T,format="wide")
#' 
#' genome_meta_exp = left_join(get_gambl_metadata() %>% 
#'     filter(seq_type=="genome") %>% 
#'     select(sample_id,pathology,lymphgen),
#'     select(gene_exp_all,-sample_id),
#'     by=c("sample_id"="genome_sample_id")) %>% 
#'      filter(!is.na(MYC))
#' 
#' # Include gene expression data and embed a box plot showing the expression of one gene across different CN states
#' 
#' 
#' pretty_CN_heatmap(cn_state_matrix=all_states_binned,
#'   these_samples_metadata = filter(genome_meta_exp,pathology=="DLBCL"),
#'   hide_annotations = "chromosomes",
#'   cluster_rows=F,
#'   geneBoxPlot = "TP53",
#'   boxplot_orientation="horizontal",bin_label_fontsize = 9,bin_label_nudge = 19
#' )
#' 
#' 
pretty_CN_heatmap = function(cn_state_matrix,
                             scale_by_sample=FALSE,
                             these_samples_metadata,
                             metadataColumns = c("pathology"),
                             expressionColumns,
                             cluster_columns=FALSE,
                             cluster_rows=TRUE,
                             show_row_names=FALSE,
                             show_column_names=FALSE,
                             keep_these_chromosomes,
                             hide_these_chromosomes,
                             keep_these_bins,
                             hide_annotations,
                             sortByBins,
                             splitByBinState,
                             sortByMetadataColumns,
                             labelTheseGenes,
                             bin_label_fontsize=5,
                             bin_label_nudge = 1.03,
                             drop_if_PGA_below=0,
                             drop_if_PGA_above=1,
                             geneBoxPlot,
                             show_bottom_annotation_name = FALSE,
                             bottom_annotation_name_side = "left",
                             left_annotation_name_side = "top",
                             bin_labels,
                             legend_direction = "horizontal",
                             legend_position = "bottom",
                             legend_row=2,
                             legend_col=3,
                             boxplot_orientation="vertical",
                             return_data = FALSE,
                             drop_bin_if_sd_below=0,
                             flip=FALSE){
  cn_state_matrix[cn_state_matrix>5] = 5
  #determine variance of each bin across the entire data set
  bin_vars = apply(cn_state_matrix,2,sd)
  old_col = ncol(cn_state_matrix)
  cn_state_matrix = cn_state_matrix[,which(bin_vars > drop_bin_if_sd_below)]
  new_col = ncol(cn_state_matrix)
  print(paste("original:",old_col,"new:",new_col))
  if(scale_by_sample){
    cn_state_matrix = cn_state_matrix - rowMeans(cn_state_matrix) + 2
    cn_state_matrix = round(cn_state_matrix)
    cn_state_matrix[cn_state_matrix<0]=0
  }

  map_bin_to_bin = function(query_region,regions=colnames(cn_state_matrix),first=TRUE){
    these_coords = suppressMessages(region_to_chunks(query_region))
    these_coords$chromosome = str_remove(these_coords$chromosome,"chr")
    all_matches = c()
    for(r in regions){
      region_coords = region_to_chunks(r)
      region_coords$chromosome = str_remove(region_coords$chromosome,"chr")
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
  bin_labels = list()
  if(!missing(geneBoxPlot)){
    if(missing(labelTheseGenes)){
      labelTheseGenes = geneBoxPlot
      expressionColumns = geneBoxPlot
      message(paste("setting expressionColumns and labelTheseGenes to",geneBoxPlot))
      
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
            keep_these_chromosomes = str_remove(this_gene_region,":.+")
          }
        }
      }
      #print(paste(g,gene_region))
      if(is.na(this_gene_region)){
        message(paste("no region for",g))
        next
      }
      if(this_gene_region %in% names(bin_labels)){
        message(paste("Gene",g, "and gene",bin_labels[[this_gene_region]], "share a region. I will only show the last one!"))
      }
      bin_labels[[this_gene_region]]=g
    }
  }else{
    if(!missing(geneBoxPlot)){
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
        keep_these_chromosomes = str_remove(this_gene_region,":.+")
      }
    }
  }
  

  
  if(!missing(these_samples_metadata)){
    keep_samples = pull(these_samples_metadata,sample_id)
    all_samples = rownames(cn_state_matrix)
    cn_state_matrix = cn_state_matrix[all_samples[all_samples %in% keep_samples],]
  }
  
  colours = map_metadata_to_colours(metadataColumns = metadataColumns, 
                                    these_samples_metadata = these_samples_metadata)
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
                "chrX"="#CA992C")
  
  #assume format is chrX:XXX-XXX
  column_chromosome = str_remove(colnames(cn_state_matrix),":.+")
  

  if(!missing(hide_these_chromosomes)){
    keep_cols = which(!column_chromosome %in% hide_these_chromosomes)
    cn_state_matrix = cn_state_matrix[,keep_cols]
    column_chromosome = column_chromosome[keep_cols]
    
  }
  splits = NULL
  if(!missing(sortByBins)){
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
    total_gain = colSums(gain_state)
    loss_state = cn_state_matrix
    loss_state[loss_state>=2]=0
    loss_state[loss_state!=0]=1
    total_loss = colSums(loss_state)
    if(!missing(expressionColumns)){
      metadataColumns = c(metadataColumns,expressionColumns)
    }
    row_df = select(these_samples_metadata,sample_id,all_of(unique(c(metadataColumns)))) %>% 
      column_to_rownames("sample_id")
    either_state = loss_state
    either_state[gain_state==1]=1
    sample_average = rowMeans(either_state)
    samples_keep_PGA = which(sample_average > drop_if_PGA_below & sample_average < drop_if_PGA_above)
    cn_state_matrix = cn_state_matrix[samples_keep_PGA,]
    sample_average = sample_average[samples_keep_PGA]
    
    #REDO
    gain_state = cn_state_matrix
    gain_state[gain_state<=2]=0
    gain_state[gain_state>0]=1
    total_gain = colSums(gain_state)
    loss_state = cn_state_matrix
    loss_state[loss_state>=2]=0
    loss_state[loss_state!=0]=1
    total_loss = colSums(loss_state)
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
      if(bin_name=="chr8:59271381-60481000"){
        print(paste(total_gain[i],total_loss[i],i))
        print(flipped_data[,i])
      }
      if(total_gain[i] > total_loss[i]){
        flipped_data[cn_state_matrix[,i]<2,i]=2
      }else{
        flipped_data[cn_state_matrix[,i]>2,i]=2
      }
      if(bin_name=="chr8:59271381-60481000"){
        print(paste(total_gain[i],total_loss[i],i))
        print(flipped_data[,i])
        print(table(flipped_data[,i],cn_state_matrix[,i]))
      }
    }

    
      cn_state_matrix = flipped_data
    }
    
    keep_rows = rownames(cn_state_matrix)
    colours[["pathology"]]=get_gambl_colours("pathology")
    colours[["lymphgen"]]=get_gambl_colours("lymphgen")
    
    #print(colours)
    if(!missing(keep_these_bins)){
      available_bins = keep_these_bins[which(keep_these_bins %in% colnames(cn_state_matrix))]
      cn_state_matrix = cn_state_matrix[,available_bins]
      print(colnames(cn_state_matrix))
      column_chromosome = str_remove(colnames(cn_state_matrix),":.+")
      total_gain = total_gain[available_bins]
      total_loss = total_loss[available_bins]
    }
    if(!missing(splitByBinState)){
      splits = round(cn_state_matrix[,splitByBinState])
      
      #anno_rows = row_df[rownames(cn_state_matrix),,drop=FALSE]
      anno_rows = row_df[rownames(cn_state_matrix),,drop=FALSE]
      anno_rows$CN_state = splits
      anno_rows$PGA = sample_average
      
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
      #print(show_legend)
      left_anno = HeatmapAnnotation(df=anno_rows,
                                    col=colours,
                                    which="row",
                                    show_legend=show_legend,
                                    annotation_name_side = left_annotation_name_side,
                                    annotation_legend_param =
                                      list(
                                           by_row=F,
                                           nrow=legend_row,
                                           ncol = legend_col,
                                           direction=legend_direction))
      
    }else{
      anno_rows = row_df[rownames(cn_state_matrix),,drop=FALSE]
      anno_rows$PGA = sample_average
      
      if(!missing(hide_annotations)){
        show_legend = rep(TRUE, length(colnames(anno_rows)))
        names(show_legend) = colnames(anno_rows)
        show_legend[hide_annotations] = FALSE
      }else{
        show_legend = rep(TRUE, length(colnames(anno_rows)))
      }
      left_anno = HeatmapAnnotation(df=anno_rows,
                                    col=colours,
                                    which="row",
                                    show_legend = show_legend,
                                    annotation_name_side = left_annotation_name_side,
                                    annotation_legend_param =
                                      list(
                                           by_row=F,
                                           nrow=legend_row,
                                           ncol = legend_col,
                                           direction=legend_direction))
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
    bin_average = colMeans(cn_state_matrix)
    cn_av_df = tibble(x=bin_average,coordinate=names(bin_average))
    cn_av_df = mutate(cn_av_df,
                      local.minima=ifelse(lag(x,n=15)>x & lead(x,n=15)  > x & x < 2,TRUE,FALSE),
                      local.maxima=ifelse(lag(x,n=15)<x & lead(x,n=15)< x & x > 2,TRUE,FALSE),
                      extreme=ifelse(local.maxima | local.minima,TRUE,FALSE))
    
    average_anno = HeatmapAnnotation(Mean_CN = bin_average,
                                     which = "column",
                                     show_annotation_name = show_bottom_annotation_name,
                                     annotation_name_side = bottom_annotation_name_side,
                                     show_legend=show_legend,
                                     col=list(Mean_CN=col_fun))
    cumulative_anno = HeatmapAnnotation(chromosome=column_chromosome,
                                        
                                        
                                        
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
                                 annotation_name_side = bottom_annotation_name_side,
                                 show_legend=show_legend,
                                 col=list(chromosome=chrom_col))
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
                   heatmap_legend_param = heatmap_legend_param
      )
  }
    
  draw(ho,heatmap_legend_side=legend_position,annotation_legend_side=legend_position)
  if(!missing(labelTheseGenes)){
    for(i in c(1:length(bin_labels))){
        gene_region = names(bin_labels)[i]
        gene_name = unname(bin_labels[[i]])
        gene_region_chunks = region_to_chunks(gene_region)
      if(gene_region_chunks$chromosome %in% column_chromosome){
          decorate_heatmap_body("CN",
                              {i=which(colnames(cn_state_matrix)==gene_region)
                              x=i/ncol(cn_state_matrix)
                              grid.text(gene_name,x,gp=gpar(fontsize=bin_label_fontsize),
                                        unit(bin_label_nudge,"npc"),rot=45,just="top")
                              grid.lines(c(x, x), c(1.01, 1), gp = gpar(lwd = 1))})
      }
        #decorate_column_title("CN",{grid.rect(gp = gpar(fill = "#00FF0040"))})
    }
      
  }

  if(return_data){
    return(list(heatmap_object=ho,
                data=cn_state_matrix,
                cumulative_gain = total_gain,
                cumulative_loss = total_loss,
                labels=bin_labels,
                chromosome_columns=column_chromosome,
                bin_means=bin_average,
                local_optima=cn_av_df))
  }
  
}

