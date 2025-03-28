
prettyJackedOncoplot = function(these_samples_metadata,
                                all_exp,
                                cn_mat,
                                genes,
                                existing_heatmaps,
                                height,
                                width,
                                heatmap_location = "top",
                                returnEverything = FALSE){
  existing_1 = existing_heatmaps[[1]]
  #ashm = pretty_stacked_oncoplot$second_heatmap
  #onco = pretty_stacked_oncoplot$oncoplot_heatmap
  
  samp_order = existing_1@column_names_param$labels
  if(!missing(all_exp)){
    #check for missing samples
    if(any(!samp_order %in% rownames(all_exp))){
      lmis = length(samp_order[!samp_order %in% rownames(all_exp)])
      print(lmis)
      stop("there are samples in your heatmap that are missing from your expression matrix")
      
    }
    some_mat = select(all_exp,any_of(genes)) %>% t()
  
  
    new_heatmap = Heatmap(some_mat[,samp_order]-rowMeans(some_mat[,samp_order]),
                         cluster_columns = T,
                         clustering_distance_columns = "euclidean",
                         show_column_names = F,
                         height = unit(height,"cm"),
                         width = unit(width, "cm"),
                         show_row_names = F,
                         row_dend_side = "right")
  }else if(!missing(cn_mat)){
    #check for missing samples
    if(any(!samp_order %in% rownames(cn_mat))){
      lmis = length(samp_order[!samp_order %in% rownames(cn_mat)])
      print(lmis)
      stop("there are samples in your heatmap that are missing from your copy number matrix")
      
    }
    #ensure metadata is arranged by samp_order
    these_samples_metadata = filter(these_samples_metadata,sample_id %in% samp_order)
    metadata_df <- these_samples_metadata[order(match(rownames(these_samples_metadata), 
                                                      samp_order)), , drop = FALSE]
    print(dim(metadata_df))
    cout = pretty_CN_heatmap(cn_mat[metadata_df$sample_id,],
                                    scale_by_sample = T,
                                    these_samples_metadata = metadata_df,
                                    legend_position = "top",
                                    rotate = T,
                                    keep_sample_order = T,
                             return_data = T,
                             hide_annotations = "chromosome",
                             width=width,height=height
                                    )
    new_heatmap = cout$heatmap_object
    #draw(new_heatmap)
    
  }
  if(heatmap_location == "top"){
    ht_list = NULL %v% new_heatmap
    for(i in c(1:length(existing_heatmaps))){
      ht_list = ht_list %v% existing_heatmaps[[i]]
    }
    #ht_list = 
    #  new_heatmap %v%
    #  ashm %v%
    #  onco
  }else{
    print("bottom")
    ht_list = NULL 
    for(i in c(1:length(existing_heatmaps))){
      ht_list = ht_list %v% existing_heatmaps[[i]]
    }
    ht_list = ht_list %v% new_heatmap
    #ht_list = 
      
    #  ashm %v%
    #  onco %v%
    #  new_heatmap
  }
  
  if(returnEverything){
    if(heatmap_location == "bottom"){
        return(c(existing_heatmaps,list(new=new_heatmap)))
    }else{
      return(c(list(new=new_heatmap),existing_heatmaps))
    }
  }else{
    draw(ht_list,
         heatmap_legend_side = "bottom",
         merge_legend = TRUE)
  }


}
