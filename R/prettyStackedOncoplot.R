#' Pretty stacked oncoplot
#'
#' @param these_samples_metadata A metadata data frame with sample_id as
#' a column. The order of sample IDs in the rows of this data frame will
#' dictate the order of samples in the oncoplot.
#' @param maf_data A data frame with maf data that will be used to populate
#' the oncoplot. Required parameter.
#' @param metadataColumns A character vector of column names that specifies
#' which columns from these_samples_metadata will be displayed below the
#' oncoplot.
#' @param sortByMetadataColumns
#' @param seg_data
#' @param sortByPGA
#' @param cn_state_matrix
#' @param genes
#' @param sortByGenes
#' @param genesCNVinOncoplot
#' @param secondPlotType Defaults to pretty_CN_heatmap, which is currently
#' the only option tested with this function.
#' @param oncoplot_location
#' @param cluster_samples
#' @param secondPlotArgs
#' @param oncoplotArgs
#' @param verbose Set to TRUE for ultra chatty mode
#'
#' @return
#' @export
#'
#' @examples
#' # Prepare some metadata
#'  dlbcl_genome_meta = get_gambl_metadata() %>% 
#'   dplyr::filter(pathology=="DLBCL",seq_type=="genome")
#' # Get CN segments for these samples
#' dlbcl_seg = get_cn_segments(dlbcl_genome_meta)
#' # Prepare CN matrix
#' cn_mat = segmented_data_to_cn_matrix(dlbcl_seg,
#'                                      these = dlbcl_genome_meta,
#'                                       adjust_for_ploidy = TRUE)
#' dlbcl_maf = get_all_coding_ssm(dlbcl_genome_meta)
#' 
#' genes=c("KMT2D","BCL2","CREBBP","EZH2","MYD88","CD79B","TP53",
#'         "PIM1","CARD11","SGK1","SOCS1",'TET2',"SPEN",
#'         "ETV6","CD83","B2M","S1PR2","GNA13","BTG1",
#'         "BTG2","DDX3X","KLHL6","HIST1H1E","TBL1XR1","SMARCA4")
#' 
#' sortGenes = c("TP53","KMT2D","BCL2","EZH2","MYD88","CD79B")
#' 
#' # oncoplot on top with column order dicted by the mutation status
#' # of sortGenes (in the order of the genes appearance in that vector)
#' 
#' prettyStackedOncoplot(these_samples_metadata = dlbcl_genome_meta,  
#'                         maf_data = dlbcl_maf,
#'                         metadataColumns = c("pathology","lymphgen"),
#'                         sortByMetadataColumns = c("pathology","lymphgen"),
#'                         cn_state_matrix = cn_mat,
#'                         genes = genes)
#' 
#' 
#' # oncoplot on top. Clustering of mutations is used to order the columns.
#' 
#' prettyStackedOncoplot(these_samples_metadata = dlbcl_genome_meta,  
#'                       maf_data = dlbcl_maf,
#'                       metadataColumns = c("pathology","lymphgen"),
#'                       cluster_samples = TRUE,
#'                       cn_state_matrix = cn_mat,
#'                       genes = genes)
#' # make a list of arguments for the second (here, the upper) plot
#' CN_args = list("keep_these_chromosomes"=c("2"),
#'                "sortByGenes" = "REL",
#'                "scale_by_sample" = TRUE,
#'                "height"=3)
#' prettyStackedOncoplot(these_samples_metadata = dlbcl_genome_meta,
#'                       maf_data = dlbcl_maf,
#'                       metadataColumns = c("pathology","lymphgen"),
#'                       oncoplot_location = "bottom",
#'                       secondPlotArgs = CN_args,
#'                       cn_state_matrix = cn_mat,
#'                       genes = genes)
#' 
#' # make a list of arguments for the second (here, the upper) plot
#' CN_args = list(
#'               "sortByGenes" = "TP53",
#'               "scale_by_sample" = TRUE,
#'               "height"=8,
#'               "hide_these_chromosomes" = "X")
#' # Specifying sortByGenes for both plots will arrange samples with relevant
#' # CN events or mutations towards one end where possible. 
#' The order in the upper plot restricts the order of the lower plot.
#' 
#' prettyStackedOncoplot(these_samples_metadata = dlbcl_genome_meta,
#'                      maf_data = dlbcl_maf,
#'                      sortByGenes = "TP53",
#'                      metadataColumns = c("pathology","lymphgen"),
#'                      oncoplot_location = "bottom",
#'                      secondPlotArgs = CN_args,
#'                      cn_state_matrix = cn_mat,
#'                      genes = genes)
#'
#' 
prettyStackedOncoplot <- function(these_samples_metadata,
                                  maf_data,
                                  metadataColumns = "pathology",
                                  sortByMetadataColumns,
                                  seg_data,
                                  sortByPGA = FALSE,
                                  cn_state_matrix,
                                  ashm_matrix,
                                  regions_bed,
                                  genes,
                                  sortByGenes,
                                  genes_CN_thresh,
                                  #genes_cn_thresh, #1 for deleted gene, 3 for amplified gene
                                  secondPlotType = "pretty_CN_heatmap", #also allowed: prettyMutationDensity
                                  oncoplot_location = "top",
                                  cluster_samples = FALSE,
                                  secondPlotArgs,
                                  oncoplotArgs,
                                  returnEverything = FALSE,
                                  plot_width = 15,
                                  oncoplotHeight = 6,
                                  secondPlotHeight = 6,
                                  verbose = FALSE) {
  plot_flavour = NULL
  if(grepl("CN",secondPlotType)){
    plot_flavour = "CN"
    # get sample_id that actually have CN data
    if (!missing(cn_state_matrix)) {
      other_samples <- rownames(cn_state_matrix)
    } else {
      print("seg_data will be used. To avoid this, provide cn_state_matrix next time")
      stop()
    }
  }else if(secondPlotType=="prettyMutationDensity"){
    plot_flavour = "mutation"
    if (!missing(ashm_matrix)) {
      other_samples <- rownames(ashm_matrix)
    } else {
      print("ashm_matrix will be generated for you To avoid this, provide ashm_matrix next time")
      stop()
    }
  }
  if (missing(maf_data)) {
    stop("maf_data must be provided")
  } else {
    maf_samples <- unique(maf_data$Tumor_Sample_Barcode)
  }

  other_samples <- dplyr::filter(
      these_samples_metadata,
      sample_id %in% other_samples) %>% 
      pull(sample_id)
    if(verbose){
        print(paste(length(other_samples), "samples with", plot_flavour, "data"))
    }
  
  maf_samples <- dplyr::filter(these_samples_metadata,
                               sample_id %in% maf_samples) %>% pull(sample_id)
  if(verbose){
    print(paste(length(maf_samples), "samples with SSM data"))
  }


  # intersect
  core_samples <- maf_samples[maf_samples %in% other_samples]
  if(verbose){
     print(paste("proceeding with", length(core_samples), "core samples"))
  }
 
  these_samples_metadata <- dplyr::filter(these_samples_metadata, 
                                          sample_id %in% core_samples)


  if (oncoplot_location == "top") {
    if(!missing(sortByGenes) | !missing(genes_CN_thresh)){
      if(!missing(genes_CN_thresh)){
        if(verbose){
          print("getting CN status for sorting")
        }

        if(!missing(sortByGenes)){
          #assume diploid if not specified
          missing_g = sortByGenes[!sortByGenes %in% names(genes_CN_thresh)]
          if(length(missing_g)){
            genes_CN_thresh[missing_g] = 2
          }
        }
        cn_thresh = data.frame(gene_id=names(genes_CN_thresh),cn_thresh = genes_CN_thresh)
        cn_state_mat = get_cnv_and_ssm_status(genes_and_cn_threshs = cn_thresh,
                              #seg_data=dlbcl_seg,
                              cn_matrix = cn_state_matrix,
                              these_samples_metadata = these_samples_metadata,
                              maf_df = maf_data,
                              only_cnv="all") 
        
        old_names= colnames(cn_state_mat)
        cn_name = sapply(old_names,function(x){paste0(x,"_cn")})
        if(!missing(sortByGenes)){
          sort_cn_name = sapply(sortByGenes,function(x){paste0(x,"_cn")})
        }
        
        colnames(cn_state_mat) = c(cn_name)
        cn_state_mat = rownames_to_column(cn_state_mat,"sample_id")
        if(verbose){
          print(head(cn_state_mat))
        }
        
      }else{
        if(!missing(sortByGenes)){
          if(verbose){
            print("making fake CN status for completeness")
          }
          if(plot_flavour=="CN"){
              cn_state_mat = data.frame(sample_id = these_samples_metadata$sample_id,
                             matrix(0,nrow=length(these_samples_metadata$sample_id),
                                    ncol=length(sortByGenes)))
            cn_name = sapply(sortByGenes,function(x){paste0(x,"_cn")})
            colnames(cn_state_mat) = c("sample_id",cn_name)
            sort_cn_name = sapply(sortByGenes,function(x){paste0(x,"_cn")})
          }
          
        
        }

      }
      if(plot_flavour=="CN"){
        these_samples_metadata = left_join(these_samples_metadata,
                                            cn_state_mat) 
      }
      if(!missing(sortByGenes)){
        if(verbose){
          print("getting mutation status for sorting")
        }
          
          ssm_state_mat = get_coding_ssm_status(gene_symbols = sortByGenes,
                                              maf_data = maf_data,
                                              these = these_samples_metadata,
                                              include_hotspots = FALSE)

          these_samples_metadata = these_samples_metadata %>%
                                  left_join(.,ssm_state_mat)
          if(plot_flavour == "CN"){
              these_samples_metadata = these_samples_metadata %>% 
                  arrange(across(any_of(sortByGenes),desc),across(any_of(sort_cn_name),desc))
          }else{
              these_samples_metadata = these_samples_metadata %>% 
                  arrange(across(any_of(sortByGenes),desc))
          }


      }
      if(plot_flavour=="CN"){
        cnv_df = select(these_samples_metadata,sample_id,all_of(cn_name)) 
        cnv_df = column_to_rownames(cnv_df,"sample_id")

        colnames(cnv_df) = gsub("_cn","",colnames(cnv_df))
      }
      if(verbose){
        print(head(cnv_df))
      }
      #select(these_samples_metadata,sample_id,all_of(c(cn_name,sortByGenes))) %>% head() %>% print()
    }
    if(verbose){
      print("making oncoplot")  
    }
    
    oncoplot_args <- list(
      maf_df = maf_data,
      these_samples_metadata = these_samples_metadata,
      keepSampleOrder = FALSE,
      cluster_cols = cluster_samples,
      plot_width = plot_width,
      simplify_annotation = TRUE,
      return_inputs = TRUE
    )

    if (!missing(genes)) {
      oncoplot_args[["genes"]] <- genes
    }
    oncoplot_args[["plot_height"]] = oncoplotHeight
    oncoplot_args[["plot_width"]] = plot_width
    if (!missing(genes_CN_thresh)) {
      oncoplot_args[['cnv_df']] = cnv_df
    }
    if(!missing(sortByGenes)){
      if(verbose){
        print("sortByGenes provided, setting cluster_samples to FALSE")
      }

      oncoplot_args[["cluster_cols"]] <- FALSE
      oncoplot_args[["keepSampleOrder"]] <- TRUE  #test
    } else if (!missing(sortByMetadataColumns)) {
      if(verbose){
        print("sortByMetadataColumns provided, setting cluster_samples to FALSE, sortByColumns")
      }
        
        oncoplot_args[["cluster_cols"]] <- FALSE
        oncoplot_args[["sortByColumns"]] <- sortByMetadataColumns
        oncoplot_args[["metadataColumns"]] <- metadataColumns
    } else {
      oncoplot_args[["cluster_cols"]] <- cluster_samples
    }
    if (!missing(oncoplotArgs)) {
      # heatmap_args <- lapply(split(c(pretty_CN_heatmap_args, heatmap_args), names(c(pretty_CN_heatmap_args, heatmap_args))), function(x) unlist(x, use.names = FALSE))
      oncoplot_args <- lapply(
        split(
          c(oncoplot_args, oncoplotArgs),
          names(c(oncoplot_args, oncoplotArgs))
        ),
        function(x) {
          if (any(sapply(x, is.data.frame))) {
            # If any element is a data frame, keep only the first data frame.
            df_indices <- which(sapply(x, is.data.frame))
            return(x[[df_indices[1]]])
          } else {
            # Otherwise, merge the atomic vectors.
            return(unlist(x, use.names = FALSE))
          }
        }
      )
      
    }
    # stop()
    if(verbose){
      print("ARGS:")
      print(names(oncoplot_args))
    }
    
    oncoplot_out <- do.call("prettyOncoplot", oncoplot_args)
    onco_heatmap <- oncoplot_out$Heatmap

    samples_order <- onco_heatmap@column_names_param$labels
    if(verbose){
        print(head(samples_order))
        print("reordering samples based on prettyOncoplot Heatmap")
    }
    
    metadata_df <- column_to_rownames(these_samples_metadata, "sample_id")
    metadata_df <- metadata_df[order(match(rownames(metadata_df), samples_order)), , drop = FALSE]
    if(verbose){
      print(dim(metadata_df))
      print(length(samples_order))
    }

    # order of oncoplot will dictate sample order
    these_samples_metadata <- rownames_to_column(metadata_df, "sample_id")
    if(verbose){
      print("making other plot")
      print("NOT clustering columns")
    }
    if(plot_flavour == "CN"){
      second_plot_args <- list(
        cn_state_matrix = cn_state_matrix[core_samples, ],
        these_samples_metadata = these_samples_metadata,
        hide_annotations = "chromosome",
        metadataColumns = metadataColumns,
        rotate = TRUE,
        sortByPGA = F,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        return_data = T,
        width = plot_width,
        height = secondPlotHeight
      )
    } else if (plot_flavour == "mutation") {
      second_plot_args <- list(
        mut_freq_matrix = ashm_matrix[core_samples, ],
        regions_bed = regions_bed,
        these_samples_metadata = these_samples_metadata,
        metadataColumns = metadataColumns,
        orientation = "sample_cols",
        cluster_samples = FALSE,
        cluster_regions = FALSE,
        width = plot_width,
        height = secondPlotHeight
      )
    }
    second_plot_args <- list(
      cn_state_matrix = cn_state_matrix[samples_order, ],
      these_samples_metadata = these_samples_metadata,
      hide_annotations = "chromosome",
      metadataColumns = metadataColumns,
      rotate = TRUE,
      sortByPGA = F,
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      return_data = T,
      width = plot_width,
      height = secondPlotHeight
    )
    if (!missing(secondPlotArgs)) {
      second_plot_args <- lapply(
        split(
          c(second_plot_args, secondPlotArgs),
          names(c(second_plot_args, secondPlotArgs))
        ),
        function(x) {
          if (any(sapply(x, is.data.frame))) {
            # If any element is a data frame, keep only the first data frame.
            df_indices <- which(sapply(x, is.data.frame))
            return(x[[df_indices[1]]])
          } else {
            # Otherwise, merge the atomic vectors.
            return(unlist(x, use.names = FALSE))
          }
        }
      )
      if(verbose){
        print(names(second_plot_args))
      }
      
    }
    # disable sortByGenes for the second plot because it will break the linkage
    if(plot_flavour=="CN"){
      second_plot_args[['sortByGenes']] = NULL
    }
    
    pcn <- do.call(secondPlotType, second_plot_args)

    second_heatmap <- pcn$heatmap_object

    ht_list <- onco_heatmap %v% second_heatmap
    if(!returnEverything){
      draw(ht_list, heatmap_legend_side = "bottom", merge_legend = TRUE)
    }

  } else { 
    # oncoplot on the bottom, other image will specify the sample order
    # get a basic mutation status matrix to allow ordering of samples to be 
    # driven by information in both plot types
        if(!missing(sortByGenes) | !missing(genes_CN_thresh)){

      if(!missing(genes_CN_thresh)){
        if(verbose){
          print("getting CN status for sorting")
        }
        
        if(!missing(sortByGenes)){
          #assume diploid if not specified
          missing_g = sortByGenes[!sortByGenes %in% names(genes_CN_thresh)]
          if(length(missing_g)){
            genes_CN_thresh[missing_g] = 2
          }
        }
        

        cn_thresh = data.frame(gene_id=names(genes_CN_thresh),cn_thresh = genes_CN_thresh)
        cn_state_mat = get_cnv_and_ssm_status(genes_and_cn_threshs = cn_thresh,
                              #seg_data=dlbcl_seg,
                              cn_matrix = cn_state_matrix,
                              these_samples_metadata = these_samples_metadata,
                              maf_df = maf_data,
                              only_cnv="all") 
        
        old_names= colnames(cn_state_mat)
        cn_name = sapply(old_names,function(x){paste0(x,"_cn")})
        if(!missing(sortByGenes)){
          sort_cn_name = sapply(sortByGenes,function(x){paste0(x,"_cn")})
        }
        
        colnames(cn_state_mat) = c(cn_name)
        cn_state_mat = rownames_to_column(cn_state_mat,"sample_id")

        
      }else{
        if(!missing(sortByGenes)){
          if(verbose){
            print("making fake CN status for completeness")
          }
          
          cn_state_mat = data.frame(sample_id = these_samples_metadata$sample_id,
                             matrix(0,nrow=length(these_samples_metadata$sample_id),
                                    ncol=length(sortByGenes)))
          cn_name = sapply(sortByGenes,function(x){paste0(x,"_cn")})
          colnames(cn_state_mat) = c("sample_id",cn_name)
          sort_cn_name = sapply(sortByGenes,function(x){paste0(x,"_cn")})
        
        }

      }
      these_samples_metadata = left_join(these_samples_metadata,
                                            cn_state_mat) 
      if(!missing(sortByGenes)){
          if(verbose){
            print("getting mutation status for sorting")
          }

          ssm_state_mat = get_coding_ssm_status(gene_symbols = sortByGenes,
                                              maf_data = maf_data,
                                              these = these_samples_metadata,
                                              include_hotspots = FALSE)

          these_samples_metadata = these_samples_metadata %>%
                                  left_join(.,ssm_state_mat) 
            these_samples_metadata = these_samples_metadata %>% 
                                  arrange(across(any_of(sort_cn_name),desc),across(any_of(sortByGenes),desc))

      }
      cnv_df = select(these_samples_metadata,sample_id,all_of(cn_name)) 
      cnv_df = column_to_rownames(cnv_df,"sample_id")
      colnames(cnv_df) = gsub("_cn","",colnames(cnv_df))
      if(verbose){
        print(head(cnv_df))
      }


      #select(these_samples_metadata,sample_id,all_of(c(cn_name,sortByGenes))) %>% head() %>% print()
    }
    if(verbose){
      print(paste("making CN heatmap, clustering:",cluster_samples))
    }
    
    second_plot_args <- list(
      hide_annotations = "chromosome",
      metadataColumns = metadataColumns,
      rotate = TRUE,
      # scale_by_sample = TRUE,
      "sortByPGA" = sortByPGA,
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      return_data = T,
      width = plot_width,
      height = secondPlotHeight
    )
    if(!missing(sortByGenes)){
      second_plot_args[["keep_sample_order"]] <- TRUE
    } else if (cluster_samples) {
      second_plot_args[["cluster_rows"]] <- TRUE
    } else if (!missing(sortByMetadataColumns)) {
      second_plot_args[["sortByMetadataColumns"]] <- sortByMetadataColumns
    }
    if(verbose){
      print(second_plot_args)
    }
    
    second_plot_args[['cn_state_matrix']] = cn_state_matrix[core_samples, ]
    second_plot_args[['these_samples_metadata']] = these_samples_metadata

    if (!missing(secondPlotArgs)) {
      # heatmap_args <- lapply(split(c(pretty_CN_heatmap_args, heatmap_args), names(c(pretty_CN_heatmap_args, heatmap_args))), function(x) unlist(x, use.names = FALSE))
      second_plot_args <- lapply(
        split(
          c(second_plot_args, secondPlotArgs),
          names(c(second_plot_args, secondPlotArgs))
        ),
        function(x) {
          if (any(sapply(x, is.data.frame))) {
            # If any element is a data frame, keep only the first data frame.
            df_indices <- which(sapply(x, is.data.frame))
            return(x[[df_indices[1]]])
          } else {
            # Otherwise, merge the atomic vectors.
            return(unlist(x, use.names = FALSE))
          }
        }
      )

    }
    pcn <- do.call(secondPlotType, second_plot_args)
    second_heatmap <- pcn$heatmap_object
    samples_order <- second_heatmap@column_names_param$labels
    if(verbose){
      print("reordering samples")
    }
    


    metadata_df <- column_to_rownames(these_samples_metadata, "sample_id")
    metadata_df <- metadata_df[order(match(rownames(metadata_df), samples_order)), , drop = FALSE]
    if(verbose){
      print("making oncoplot")
    }
    
    these_samples_metadata <- rownames_to_column(metadata_df, "sample_id")
    
    oncoplot_args = list("maf_df" = maf_data,
                          "these_samples_metadata" = these_samples_metadata,
                          "genes" = genes,
                          "keepSampleOrder" = TRUE,
                          "cluster_cols" = cluster_samples,
                          "plot_width" = plot_width,
                          "plot_height" = oncoplotHeight,
                          "simplify_annotation" = TRUE,
                          "return_inputs" = TRUE
    )

    if(!missing(genes_CN_thresh)){
      oncoplot_args[["cnv_df"]] = cnv_df
    }
    if (!missing(oncoplotArgs)) {
      oncoplot_args <- lapply(
        split(
          c(oncoplot_args, oncoplotArgs),
          names(c(oncoplot_args, oncoplotArgs))
        ),
        function(x) {
          if (any(sapply(x, is.data.frame))) {
            # If any element is a data frame, keep only the first data frame.
            df_indices <- which(sapply(x, is.data.frame))
            return(x[[df_indices[1]]])
          } else {
            # Otherwise, merge the atomic vectors.
            return(unlist(x, use.names = FALSE))
          }
        }
      )

    }
    if(verbose){
      print(names(oncoplot_args))
    }
    oncoplot_out <- do.call("prettyOncoplot", oncoplot_args)

    onco_heatmap <- oncoplot_out$Heatmap

    ht_list <- second_heatmap %v% onco_heatmap
    if(!returnEverything){
      draw(ht_list, heatmap_legend_side = "bottom", merge_legend = TRUE)
    }
    
  }
  cn_val <- pcn$data
  cn_val[cn_val != 2] <- 1
  cn_val[cn_val == 2] <- 0
  if(returnEverything){


    to_return = list(
      second_heatmap = second_heatmap,
      oncoplot_heatmap = onco_heatmap,
      ordered_metadata = these_samples_metadata,
      oncoplot_mut_status = oncoplot_out$mut_status,
      CN_values = pcn$data,
      CN_status = cn_val)
    if(!missing(sortByGenes) | !missing(genes_CN_thresh)){
      to_return[["cnv_df"]] = cnv_df
      # make a combined status data frame for convenience
      combined_df = oncoplot_out$mut_status
      for(gene in colnames(cnv_df)){
        combined_df[cnv_df[,gene]==1,gene] = 1
      }
      to_return[["combined_mut_status"]] = combined_df
    }
    return(to_return)
  }
}
