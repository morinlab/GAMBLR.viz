#' pretty_circular_mutation_frequency_heatmap
#'
#' @param prettyOncoplot_output The output of the prettyOncoplot function
#' @param cn_status_matrix The output of get_cn_states
#' @param collated_results A list of data frames with sample_id as rownames and features as column names
#' @param these_samples_metadata A data frame with metadata. Usually the output of [GAMBLR.results::get_gambl_metadata].
#' @param genes A vector of genes to label
#' @param cluster Whether to perform clustering. Default is TRUE (clustering is performed).
#' @param keep_these_pathologies A vector of pathology values to show in the plot. All the remaining rows will be ignored.
#' @param min_sample_num Minimum number of samples in a pathology to be considered for the plot. Pathologies with less than this number will be excluded. (20)
#' @param col_fun Color function to modify the default color pallette of the heatmap.
#' @param col_theme Alternatively, provide the color theme instead of  `col_fun` to change the default colors of the heatmap.
#' @param return_data Conditionally return the formatted data used for the plotting. Default is FALSE (only image is plotted and no data is returned).
#' @param dend_location Location of the dendrogram. Default is "inside".
#' @param clustering_distance_method Clustering method. Default is "euclidean".
#' @param border Whether to draw border around heatmap. Default is TRUE (with border).
#' @param split_by_type Whether to split the mutations by type. Default is FALSE (no splitting).
#' @param rotate_degrees Rotate labels. Default is 0 (no rotation).
#' @param gap.degree Gap degree. Default is 15.
#' @param show.sector.labels Show labels for each sector of the heatmap. Default is FALSE (no labels).
#' @param label_cex Number indicating the amount by which plotting text and symbols should be scaled relative to the default when displaying the labels. Default is 0.5.
#' @param rownames_cex Number indicating the amount by which plotting text and symbols should be scaled relative to the default when displaying the rownames. Default is 0.5.
#' @param include_legend Whether to include the legend. Default is FALSE (no legend).
#' @param colour_labels Optionally color labels. Default is FALSE (no coloring).
#' @param label_group How to group the labels. Default is "text".
#' @param label_alpha Value from 0 to 1 to control alpha of the label.
#'
#' @return Nothing or a list of data frames (when return_data = TRUE)
#' @export
#'
#' @examples
#'
#' library(GAMBLR.data)
#' 
#' metadata <- get_gambl_metadata() %>%
#'    filter(!seq_type == "mrna") %>%
#'    filter(pathology %in% names(get_gambl_colours("pathology"))) %>%
#'    distinct(sample_id, .keep_all = TRUE)
#' 
#' all_coding <- get_coding_ssm(these_samples_metadata = metadata)
#' 
#' genes <- lymphoma_genes %>%
#'     filter(DLBCL|FL|BL) %>%
#'     pull(Gene) %>%
#'     unique %>%
#'     sort
#' 
#' oncoplot_output <- prettyOncoplot(
#'     all_coding,
#'     genes = genes,
#'     minMutationPercent = 2,
#'     these_samples_metadata = metadata,
#'     simplify = TRUE,
#'     return_inputs = TRUE
#' )
#' 
#' # Basic plot
#' pretty_circular_mutation_frequency_heatmap(
#'     prettyOncoplot_output = oncoplot_output,
#'     keep_these_pathologies = c(
#'         "FL", "DLBCL", "PMBCL", "BL", "HGBL"
#'     )
#' )
#' 
#' # Add sv layer
#' all_sv <- get_manta_sv(these_samples_metadata = metadata)
#' annotated_sv <- annotate_sv(all_sv) %>%
#'     filter(gene %in% genes, !is.na(partner)) %>%
#'     select(sample_id = tumour_sample_id, gene)
#' 
#' # This is to replicate the output format of collate_sv
#' sv_collated <- annotated_sv %>%
#'     mutate(
#'         gene = paste("manta", gene, "sv", sep = "_"),
#'         mutated = "POS"
#'     ) %>%
#'     distinct %>%
#'     pivot_wider(
#'         names_from = gene,
#'         values_from = mutated
#'     ) %>%
#'     replace(is.na(.), "NEG")
#' 
#' # Plot SSM + SVs
#' pretty_circular_mutation_frequency_heatmap(
#'     collated_results = list(sv_collated),
#'     prettyOncoplot_output = oncoplot_output,
#'     these_samples_metadata = metadata
#' )
#' 
#' # Add aSHM data
#' ashm_freq <- get_ashm_count_matrix(
#'     regions_bed = mutate(
#'         grch37_ashm_regions,
#'         name = paste(gene, region, sep = "_")
#'     ),
#'     this_seq_type = "genome"
#' )
#' 
#' ashm_freq_collated <- mutate(ashm_freq,across(,~ifelse(.x>0,1,0)))
#' 
#' ashm_freq_collated <- ashm_freq_collated[,colSums(ashm_freq_collated) >130]
#' ashm_freq_collated <- rownames_to_column(ashm_freq_collated,"sample_id")
#' 
#' # Comprehensive plot with SSM + SV + aSHM and some non-default arguments
#' pretty_circular_mutation_frequency_heatmap(
#'     collated_results = list(sv_collated, ashm_freq_collated),
#'     prettyOncoplot_output = oncoplot_output,
#'     these_samples_metadata = metadata,
#'     keep_these_pathologies = c("DLBCL", "FL", "BL"),
#'     split_by_type = TRUE,
#'     colour_labels = TRUE,
#'     label_cex = 0.4,
#'     rownames_cex = 0.4,
#'     include_legend = TRUE
#' )
#'
pretty_circular_mutation_frequency_heatmap = function(prettyOncoplot_output,
                                                      cn_status_matrix,
                                                      collated_results,
                                                      these_samples_metadata,
                                                      genes,
                                                      cluster=T,
                                                      keep_these_pathologies,
                                                      min_sample_num=20,
                                                      col_fun,
                                                      col_theme,
                                                      return_data=FALSE,
                                                      dend_location = "inside",
                                                      clustering_distance_method="euclidean",
                                                      border=T,
                                                      split_by_type=FALSE,
                                                      rotate_degrees=0,
                                                      gap.degree=15,
                                                      show.sector.labels=FALSE,
                                                      label_cex=0.5,
                                                      rownames_cex=0.5,
                                                      include_legend=F,
                                                      colour_labels=F,
                                                      label_group="text",
                                                      label_alpha){
  if(!missing(col_theme)){
    if(col_theme == "sunset"){
      col_fun=colorRamp2(c(0,2.5,5,10,20,30,40,50,80),hcl_palette = "ag_Sunset",reverse=T)
    }else if(col_theme=="lajolla"){
      col_fun=colorRamp2(c(0,2.5,5,10,20,30,40,50,80),hcl_palette = "Lajolla",reverse=F)
    }else if(col_theme=="inferno"){
      col_fun=colorRamp2(c(0,2.5,5,10,20,30,40,50,80),hcl_palette = "Inferno",reverse=T)
    }else if(col_theme=="red"){
      col_fun=colorRamp2(c(0,2.5,5,10,20,30,40,50,80),hcl_palette = "Reds",reverse=T)
    }else if(col_theme=="blue"){
      col_fun=colorRamp2(c(0,2.5,5,10,20,30,40,50,80),hcl_palette = "PuBu",reverse=T)
    }
  }
  if(missing(col_fun)){
    col_fun = colorRamp2(c(0,2.5,5,10,20,30,40,50,80),c("white",
                                                        "#E8DACD80",
                                                    "#E1816588",
                                                     "#f0817f",

                                                     "#B0174A99",
                                                     "#91282699",
                                                     "#790821",
                                                     "#5D2E46",
                                                     "black"))
  }
  if(missing(label_alpha)){
    label_alpha=1
  }
  #if(!missing(prettyOncoplot_output) & !missing(cn_status_matrix)){
  #
  #  cn_status_rn = rownames_to_column(cn_status_matrix,"sample_id")
  #  mut_status_rn = left_join(mut_status_rn,cn_status_rn,by="sample_id")
  #  mut_status_rn[is.na(mut_status_rn)]=0
  splits = c()
  if(missing(these_samples_metadata)){
    if(!missing(prettyOncoplot_output)){
      message("defaulting to metadata in prettyOncoplot output")
      these_samples_metadata = rownames_to_column(prettyOncoplot_output$metadata,"sample_id")
      mut_status = rownames_to_column(prettyOncoplot_output$mut_status,"sample_id")
      mut_status_meta = left_join(select(these_samples_metadata,sample_id,pathology),
                                  mut_status)
    }else{
      message("getting metadata for samples in the matrix provided. For more control over this, supply these_samples_metadata")
      these_samples_metadata = get_gambl_metadata() %>%
        filter(sample_id %in% mut_status_rn$sample_id)
      mut_status_meta = select(these_samples_metadata,sample_id,pathology)
    }
  }else{
    mut_status_meta = select(these_samples_metadata,sample_id,pathology)
    if(!missing(prettyOncoplot_output)){
      mut_status = rownames_to_column(prettyOncoplot_output$mut_status,"sample_id")
      mut_status_meta = left_join(select(these_samples_metadata,sample_id,pathology),
                                  mut_status)
    }
  }
  if(ncol(mut_status_meta)>2){
    splits = rep("Mutation",ncol(mut_status_meta)-2)

  }


  if(!missing(collated_results)){
    #the main (possibly only) source of data to include in the plot
    fix_pos_neg = function(x){
      x=case_when(is.numeric(x) & x > 0 ~ 1,
                  is.numeric(x) ~ 0,
                  x=="POS" ~ 1,
                  TRUE ~ 0)
    }
    for(i in c(1:length(collated_results))){
      result_df = collated_results[[i]]
      result_name = names(collated_results)[i]
      num_col = ncol(result_df)-1

      splits = c(splits,rep(result_name,num_col))
      result_df = mutate(result_df,across(-sample_id,~fix_pos_neg(.x)))
      mut_status_meta = left_join(mut_status_meta,result_df)
    }
    #mut_status_rn = column_to_rownames(mut_status,"sample_id")



  }
  genes = colnames(mut_status_meta)[c(3:ncol(mut_status_meta))]

  mut_status_meta = mutate(mut_status_meta,across(all_of(genes),~as.numeric(.)))

  calc_percent = function(x){
    percent = 100*sum(x,na.rm=T)/sum(!is.na(x))
    if(is.nan(percent)){
      percent = 0
    }
    return(percent)
  }
  mut_sums = group_by(mut_status_meta,pathology) %>%
    mutate(total=n()) %>%
    group_by(total,add=T) %>%
    summarise(across(all_of(genes),~calc_percent(.x))) %>%
    as.data.frame()



  if(!missing(keep_these_pathologies)){
    mut_sums = filter(mut_sums,pathology %in% keep_these_pathologies) %>%
      select(-total)
    keep_these_pathologies= keep_these_pathologies[keep_these_pathologies %in% mut_sums$pathology]
  }else{
    mut_sums = filter(mut_sums,total >= min_sample_num) %>%
      select(-total)
  }
  circos_mat = t(column_to_rownames(mut_sums,"pathology"))
  #make rows genes, columns pathologies
  circos_mat = circos_mat[,keep_these_pathologies]



  cn = rev(colnames(circos_mat))
  n = length(cn)

  cn_col = get_gambl_colours("pathology",alpha = label_alpha)
  cn_col = cn_col[rev(colnames(circos_mat))]



  if(colour_labels){
    ncat = length(unique(splits))
    label_pal = brewer.pal(ncat,"Dark2")
    names(label_pal) = unique(splits)
    rownames.col = label_pal[splits]
  }else{
    rownames.col = "black"
  }

  if(border){
    cell.border=0.2
  }else{
    cell.border=NA
  }
  if(dend_location == "outside"){
    circos.clear()
    circos.par(start.degree=75)
    circos.par(gap.degree = c(15))
    if(!split_by_type){
      splits= NULL
    }
    circos.heatmap(circos_mat,
                   col = col_fun,
                   cluster = cluster,
                   dend.side = "outside",
                   split=splits,
                   rownames.cex = rownames_cex,
                   rownames.side = "inside",
                   track.height=0.4,distance.method = clustering_distance_method)


    circos.track(track.index = get.current.track.index(),
                 panel.fun = function(x, y) {
                   if(CELL_META$sector.numeric.index == 1) { # the last sector
                     cn = rev(colnames(circos_mat))
                     n = length(cn)

                     circos.text(rep(CELL_META$cell.xlim[2], n) +
                                   convert_x(1, "mm"),
                                 1:n/3 + 1,
                                 cn,
                                 cex = 0.5,
                                 adj = c(0, 1), facing = "inside")
                   }
                 }, bg.border = NA)


  }else{
    circos.clear()

    circos.par(start.degree=75+ rotate_degrees)


    circos.par(gap.degree = gap.degree)
    if(!split_by_type){
      splits= NULL
    }


    circos.heatmap(circos_mat,
                   col = col_fun,
                   cluster = cluster,
                   dend.side = "inside",
                   split=splits,
                   rownames.cex = rownames_cex,
                   rownames.side = "outside",
                   rownames.col = rownames.col,
                   track.height=0.4,
                   cell.border=cell.border,
                   show.sector.labels=show.sector.labels)


    circos.track(track.index = 2,
                 panel.fun = function(x, y) {
                   #if(CELL_META$sector.numeric.index == 1) { # the last sector

                     if(label_group == "colour"){
                       circos.rect(rep(CELL_META$cell.xlim[2], n) +
                                     convert_x(1, "mm")+6.5,
                                   1:n-1,
                                   rep(CELL_META$cell.xlim[2], n) +
                                     convert_x(1, "mm")-0.5,
                                   1:n,col=cn_col,border=FALSE)

                     }else if(label_group == "text"){
                       circos.text(rep(CELL_META$cell.xlim[2], n) +
                                     convert_x(1, "mm"),
                                   1:n-0.25,
                                   cn,
                                   cex = label_cex,
                                   adj = c(0, 1), facing = "inside",niceFacing=T)
                     }

                   #}
                 }, bg.border = NA)


  }
  if(include_legend){
    lgd = Legend(at = c(0,2.5,5,10,20,30,40,50,80),
                 col_fun = col_fun,
                 title_position = "topcenter",
                 title = "Frequency")

    if(colour_labels){

      lgd2 = Legend(at = names(label_pal),
                                 legend_gp = gpar(col = label_pal),
                                 labels_gp = gpar(col = label_pal),
                                 title_position = "topcenter",
                                 title = "Type")
      if(label_group=="colour"){

        lgd3 = Legend(at = names(cn_col), type="points",
                      legend_gp = gpar(col = cn_col),
                      #labels_gp = gpar(col = cn_col),
                      title_position = "topcenter",
                      background="white",
                      title = "Group")
        lgd_list = packLegend(lgd,lgd2,lgd3,direction = "horizontal")
      }else{
        lgd_list = packLegend(lgd,lgd2,direction = "horizontal")
      }

    }else{
      if(label_group=="colour"){

        lgd3 = Legend(at = names(cn_col), type="points",
                      legend_gp = gpar(col = cn_col),
                      #labels_gp = gpar(col = cn_col),
                      title_position = "topcenter",
                      background="white",
                      title = "Group")
        lgd_list = packLegend(lgd,lgd3,direction = "horizontal")
      }else{
        lgd_list = packLegend(lgd)
      }

    }
    draw(lgd_list, x = unit(12, "mm"), y = unit(12, "mm"), just = c("left", "bottom"))
  }
  if(return_data){
    return(list(mut_status=mut_status_meta,mat=circos_mat,mut_sums=mut_sums))
  }
}
