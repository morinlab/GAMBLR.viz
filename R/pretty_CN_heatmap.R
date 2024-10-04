
#' Pretty CN heatmap
#'
#' @param cn_state_matrix 
#' @param these_samples_metadata 
#' @param top_barplot 
#' @param cluster_columns 
#' @param show_row_names 
#' @param show_column_names 
#'
#' @return
#' @export
#'
#' @examples
pretty_CN_heatmap = function(cn_state_matrix,
                             these_samples_metadata,
                             top_barplot = TRUE,
                             cluster_columns=FALSE,
                             cluster_rows=TRUE,
                             show_row_names=FALSE,
                             show_column_names=FALSE,
                             infer_chromosomes_from_regions=TRUE,
                             keep_these_chromosomes,
                             hide_these_chromosomes,
                             split_chromosome=FALSE,
                             column_chromosome,
                             sort_by_bins,
                             sort_by_columns,
                             label_bins,
                             bin_label_fontsize=5,
                             bin_label_nudge = 1.03){
  if(!missing(these_samples_metadata)){
    keep_samples = pull(these_samples_metadata,sample_id)
    all_samples = rownames(cn_state_matrix)
    cn_state_matrix = cn_state_matrix[all_samples[all_samples %in% keep_samples],]
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
  if(infer_chromosomes_from_regions){
    #assume format is chrX:XXX-XXX
    column_chromosome = str_remove(colnames(cn_state_matrix),":.+")
    print(head(column_chromosome))
  }else{
    if(missing(column_chromosome)){
      stop("you must provide column_chromosome if you disable chromosome inference")
    }
    
  }
  if(!missing(keep_these_chromosomes)){
    keep_cols = which(column_chromosome %in% keep_these_chromosomes)
    cn_state_matrix = cn_state_matrix[,keep_cols]
    column_chromosome = column_chromosome[keep_cols]

  }
  if(!missing(hide_these_chromosomes)){
    keep_cols = which(!column_chromosome %in% hide_these_chromosomes)
    cn_state_matrix = cn_state_matrix[,keep_cols]
    column_chromosome = column_chromosome[keep_cols]
    
  }
  
  if(!missing(sort_by_bins)){
    cluster_rows = FALSE
    cn_state_matrix = arrange(cn_state_matrix,desc(across(sort_by_bins)))
  }else if(!missing(sort_by_columns)){
    cluster_rows = FALSE
    new_order = arrange(these_samples_metadata,lymphgen) %>% pull(sample_id)
    cn_state_matrix = cn_state_matrix[new_order,]
  }
  
  if(top_barplot){
    gain_state = cn_state_matrix
    gain_state[gain_state<=2]=0
    gain_state[gain_state>0]=1
    total_gain = colSums(gain_state)
    loss_state = cn_state_matrix
    loss_state[loss_state>=2]=0
    loss_state[loss_state!=0]=1
    total_loss = colSums(loss_state)
    print(length(total_gain))
    print(length(total_loss))
    row_df = select(these_samples_metadata,sample_id,pathology,lymphgen,COO_consensus,DHITsig_consensus) %>% 
      column_to_rownames("sample_id")
    left_anno = HeatmapAnnotation(df=row_df[rownames(cn_state_matrix),],
                                  col=list(pathology=get_gambl_colours("pathology"),
                                           lymphgen=get_gambl_colours("lymphgen"),
                                           COO_consensus=get_gambl_colours("COO")),
                                  which="row")
    top_anno = HeatmapAnnotation(chromosome=column_chromosome,
                                 which = "column",col=list(chromosome=chrom_col),
                                 Gain = anno_barplot(total_gain,
                                                     gp=gpar(fill="red",col="red")),
                                 loss=anno_barplot(total_loss,
                                                   gp=gpar(fill="blue",col="blue"),
                                                   axis_param = list(direction = "reverse"))
                                                   )
    #bottom_anno = HeatmapAnnotation()
    ho = Heatmap(cn_state_matrix,name="CN",column_title=" ",
            cluster_columns=cluster_columns,
            cluster_rows=cluster_rows,
            show_row_names = show_row_names,show_column_names = show_column_names,
            col = col_fun,
            #top_annotation = bottom_anno,
            bottom_annotation = top_anno,
            left_annotation=left_anno,
            #column_split  = column_chromosome,
            heatmap_legend_param = list(
              title = "CN"))
    draw(ho)
    if(!missing(label_bins)){
      for(i in c(1:length(label_bins))){
        gene_region = names(label_bins)[i]
        gene_name = unname(label_bins[[i]])

        decorate_heatmap_body("CN",
                              {i=which(colnames(cn_state_matrix)==gene_region)
                              x=i/ncol(cn_state_matrix)
                              grid.text(gene_name,x,gp=gpar(fontsize=bin_label_fontsize),
                                        unit(bin_label_nudge,"npc"),rot=45,just="top")
                              grid.lines(c(x, x), c(1.01, 1), gp = gpar(lwd = 1))})
        #decorate_column_title("CN",{grid.rect(gp = gpar(fill = "#00FF0040"))})
      }
      
    }
  }
  else{
    Heatmap(cn_state_matrix,cluster_columns=cluster_columns,show_row_names = show_row_names,show_column_names = show_column_names,
               col = col_fun)
  }
  
  return(list(heatmap_object=ho,data=cn_state_matrix))
}