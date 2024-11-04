#' pretty_circular_mutation_frequency_heatmap
#'
#' @param prettyOncoplot_output The output of the prettyOncoplot function
#' @param cn_status_matrix The output of get_cn_states
#' @param collated_results A list of data frames with sample_id as rownames and features as column names
#' @param genes A vector of genes to label
#' @param keep_these_pathologies A vector of pathology values to show in the plot. All the remaining rows will be ignored.
#' @param min_sample_num Minimum number of samples in a pathology to be considered for the plot. Pathologies with less than this number will be excluded. (20)
#'
#' @return Nothing or a list of data frames (when return_data = TRUE)
#' @export
#'
#' @examples
#'
#' all_gambl_meta = get_gambl_metadata() %>%
#'     filter(!seq_type == "mrna") %>%
#'     filter(pathology %in% names(get_gambl_colours("pathology")))
#'
#' all_coding = get_all_coding_ssm(these_samples_metadata = all_gambl_meta)
#'
#' genes = filter(GAMBLR.data::lymphoma_genes_dlbcl_v_latest,curated==TRUE) %>%
#'     pull(Gene)
#' genes = unique(c(genes,filter(GAMBLR.data::lymphoma_genes_mcl_v_latest,,curated==TRUE) %>% pull(Gene)))
#' genes = unique(c(genes,filter(GAMBLR.data::lymphoma_genes_bl_v_latest,,curated==TRUE) %>% pull(Gene)))
#' oncoplot_output = prettyOncoplot(all_coding,
#'                                 genes=genes,
#'                                 minMutationPercent = 2,
#'                                 these_samples_metadata = all_gambl_meta,
#'                                 simplify = T,return_inputs = T)
#'
#' pretty_circular_mutation_frequency_heatmap(prettyOncoplot_output = oncoplot_output,
#'                                            keep_these_pathologies = c("FL",
#'                                                "DLBCL",
#'                                                "PMBCL",
#'                                                "BL",
#'                                                "HGBL"))
#'
#' genes_and_cn_threshs =
#'     data.frame(gene_id=c("MYC",
#'     "MIR17HG", "TNFAIP3","TCF4",
#'     "TNFRSF14","REL","CD274",
#'     "FAS","HNRNPD","B2M","TP53",
#'     "CDKN2A","RB1","PTEN","CD70",
#'     "TOX","RHOA","CD58","MCL1",
#'     "ETV6","TMEM30A"),
#'     cn_thresh=c(3,3,1,3,1,3,3,1,1,1,1,1,1,1,1,1,1,1,3,1,1)) %>%
#'     mutate(name=ifelse(cn_thresh>2,paste0(gene_id,"_gain"),paste0(gene_id,"_loss")))
#'
#' seg_data = get_cn_segments()
#'
#' cn_status = get_cnv_and_ssm_status(only_cnv="all",
#'                         these_samples_metadata = all_gambl_meta,
#'                         genes_and_cn_threshs = genes_and_cn_threshs,
#'                         adjust_for_ploidy=T)
#'
#'
#' pretty_circular_mutation_frequency_heatmap(cn_status_matrix = cn_status,
#'                                            prettyOncoplot_output = oncoplot_output)
#'
#' sv_collated = GAMBLR.results:::collate_sv_results() %>%
#'     select(sample_id,ends_with("sv"))
#'
#' NFKBIZ_genome = GAMBLR.results:::collate_nfkbiz_results() %>%
#'     select(sample_id,NFKBIZ_UTR)
#' NFKBIZ_capture = GAMBLR.results:::collate_nfkbiz_results(seq_type_filter="capture") %>%
#'     select(sample_id,NFKBIZ_UTR)
#'
#' HNRNPH1_genome = GAMBLR.results:::collate_hnrnph1_mutations() %>%
#'     select(sample_id,HNRNPH1_splice)
#' HNRNPH1_capture = GAMBLR.results:::collate_hnrnph1_mutations(seq_type_filter="capture") %>%
#'     select(sample_id,HNRNPH1_splice)
#'
#' NFKBIZ = bind_rows(NFKBIZ_genome,NFKBIZ_capture)
#' HNRNPH1 = bind_rows(HNRNPH1_genome,HNRNPH1_capture)
#'
#' pretty_circular_mutation_frequency_heatmap(cn_status_matrix = cn_status,
#'                                            collated_results = list(sv_collated,
#'                                                NFKBIZ,
#'                                                HNRNPH1),
#'                                            prettyOncoplot_output = oncoplot_output,
#'                                            these_samples_metadata = all_gambl_meta)
#'
#'
#' all_states_binned = get_cn_states(n_bins_split=2500,
#'     missing_data_as_diploid = T,
#'     seg_data = seg_data)
#'
#' CN_out = pretty_CN_heatmap(cn_state_matrix=all_states_binned,
#'     these_samples_metadata = all_gambl_meta,
#'     hide_annotations = "chromosome",
#'     scale_by_sample=T,
#'     return_data = T)
#'
#' arm_level_events = categorize_CN_events(CN_out)
#'
#'
#' arm_level_annotated = rownames_to_column(arm_level_events$arm_events_simplified,"sample_id")
#'
#' pretty_circular_mutation_frequency_heatmap(cn_status_matrix = cn_status,
#'                                            collated_results = list(sv_collated,
#'                                                NFKBIZ,
#'                                                HNRNPH1,
#'                                                arm_level_annotated),
#'                                            prettyOncoplot_output = oncoplot_output,
#'                                            these_samples_metadata = all_gambl_meta)
#'
#'
#' ashm_freq = get_ashm_count_matrix(
#'                  regions_bed = dplyr::mutate(GAMBLR.data::grch37_ashm_regions,
#'                                name = paste(gene, region, sep = "_")),
#'                  this_seq_type = "genome"
#'                  )
#' ashm_freq_collated = mutate(ashm_freq,across(,~ifelse(.x>0,1,0)))
#'
#'  ashm_freq_collated = ashm_freq_collated[,colSums(ashm_freq_collated) >130]
#'  ashm_freq_collated = rownames_to_column(ashm_freq_collated,"sample_id")
#'
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
