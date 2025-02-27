#' Pretty mutual exclusivity plot
#'
#' @param maf_data 
#' @param mutmat 
#' @param genes 
#' @param these_samples_metadata 
#' @param q_threshold 
#' @param exclude_insignificant_genes 
#' @param engine 
#' @param use_MCC 
#' @param font_size 
#' @param use_alpha 
#' @param clustering_distance 
#' @param gene_anno_df 
#' @param size_factor 
#' @param split 
#' @param return_data Default False
#' @param include_silent Default False
#' @param include_hotspots Default False
#' @param review_hotspots Default False
#' @param bonferroni Default False
#'
#' @returns
#' @export
#'
#' @examples
#' dlbcl_meta = get_gambl_metadata() %>% dplyr::filter(pathology=="DLBCL",seq_type!="mrna")
#' all_coding = GAMBLR.open::get_all_coding_ssm(these = dlbcl_meta)
#' 
#' pretty_mutual_exclusivity(maf_data = all_coding, 
#'                                genes = lg_genes,
#'                                these=dlbcl_meta,
#'                                engine = "ComplexHeatmap",
#'                                font_size = 6,
#'                                use_alpha = T,
#'                                clustering_distance="binary",
#'                                include_hotspots = T)
pretty_mutual_exclusivity = function(maf_data,
                                          mutmat,
                                          genes,
                                          these_samples_metadata,
                                          q_threshold = 0.1,
                                          exclude_insignificant_genes=TRUE,
                                          engine = "ggcorrplot",
                                          use_MCC=FALSE,
                                          font_size=7,use_alpha=FALSE,
                                          clustering_distance="binary",
                                          gene_anno_df,
                                          size_factor=0.01,
                                          split,
                                          return_data=F,
                                          include_silent = F,
                                          include_hotspots = F,
                                          review_hotspots = F,
                                          bonferroni = FALSE){
  if(missing(mutmat)){
    if(missing(genes)){
      stop("genes is a required argument")
    }
    print("generating mutation matrix and will return it for convenience")
    
    mutmat = get_coding_ssm_status(gene_symbols = genes,
                                   include_hotspots = include_hotspots,
                                   these_samples_metadata = these_samples_metadata,
                                   include_silent = include_silent,
                                   maf_data = maf_data,review_hotspots = review_hotspots) 
    mutmat = column_to_rownames(mutmat,"sample_id")
    print(dim(mutmat))
  }else{
    if(missing(genes)){
      genes = colnames(mutmat)
    }
    
  }
  mcor = cor(mutmat)
  cpmat = cor_pmat(mutmat)
  #bonferroni correction
  if(bonferroni){
    cpmat = cpmat *nrow(cpmat)^2
  }
  
  M = mcor
  M[cpmat>q_threshold] = 0
  diag(M) <- 0
  M[is.na(M)]=0
  
  
  
  if(engine == "ggcorrplot"){
    
    
    pp = ggcorrplot(M,
                    method = "circle",
                    #size=1,
                    #hc.order = TRUE,
                    #type="upper",
                    lab = FALSE
    ) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                           #breaks = c(-lim, round(-lim/2), 0, round(lim/2), lim),
                           name = "Odds Ratio") +
      ggtitle("Mutation Correlation and Exclusivity (pAdj < 0.1)") +
      theme_minimal() +
      theme(legend.position = "none") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    return(list(plot=pp))
  }
  
  else if(engine=="ComplexHeatmap"){
    # Determine the maximum absolute value in M for scaling:
    M[is.na(M)]=0
    #min_obs = abs(min(M))
    #M[ M> min_obs ]=min_obs
    #M[ M< -min_obs]= -min_obs
    lim <- max(abs(M), na.rm = TRUE)  
    
    
    
    
    #print(gene_anno)
    if(!missing(gene_anno_df)){
      gene_anno_df = gene_anno_df %>% 
        distinct(Feature,.keep_all=TRUE)  %>% 
        #select(Feature,Class) %>%
        filter(Feature %in% rownames(M)) 
      
      missing_anno = rownames(M)[!rownames(M) %in% gene_anno_df$Feature]
      
      
      if(length(missing_anno)>0){
        missing_df = data.frame(Feature=missing_anno,Class="NA")
        gene_anno = bind_rows(gene_anno_df,missing_df) %>% 
          distinct(Feature,.keep_all=TRUE)  %>% 
          column_to_rownames("Feature")
      }else{
        gene_anno = gene_anno_df %>%
          distinct(Feature,.keep_all=TRUE)  %>% 
          column_to_rownames("Feature")
      }
      
      cols =c(get_gambl_colours("lymphgen"),
              get_gambl_colours("genetic_subgroup"),
              "NA"="grey")
      gene_anno = gene_anno[rownames(M),]
      rann = HeatmapAnnotation(df=gene_anno,which="row",
                               col=list(Class=cols,
                                        DLBCL=c("1"=unname(get_gambl_colours("pathology")["DLBCL"]),"0"="white"),
                                        FL=c("1"=unname(get_gambl_colours("pathology")["FL"]),"0"="white"),
                                        BL=c("1"=unname(get_gambl_colours("pathology")["BL"]),"0"="white")))
      #cann = HeatmapAnnotation(df=gene_anno,which="column",
      #                         col=list(Class=cols,
      #                                  DLBCL=get_gambl_colours("pos_neg"),
      #                                  FL = get_gambl_colours("pos_neg"),
      #                                  BL = get_gambl_colours("pos_neg")))
    }
    
    
    
    if(use_alpha){
      # Define your cell function with transparency added to the fill color
      cell_fun <- function(j, i, x, y, width, height, fill) {
        # Get the current value for the cell.
        value <- M[i, j]
        if (!is.na(value)) {
          # Compute a radius proportional to the magnitude of the value.
          # Here we set the radius in "npc" units; you might adjust the factor as needed.
          r <- unit(size_factor * abs(value) / lim, "npc")
          # Adjust the fill color to include transparency (alpha = 0.5 means 50% opaque)
          grid.circle(x, y, r = r, gp = gpar(fill = adjustcolor(fill, alpha.f = 0.5), col = fill))
        }
      }
      
      # Create the heatmap using the modified cell_fun
      if(missing(gene_anno_df)){
        rann = NULL
      }
      heatmap_args = list(matrix=M,
                          name = "Overlap/\nExclusivity",
                          col = colorRamp2(c(-lim, 0, lim), c("blue", "white", "red")),
                          cell_fun = cell_fun,
                          rect_gp = gpar(fill = NA, col = NA),
                          na_col = "white",  # Cells that are NA are shown in white.
                          cluster_rows = TRUE,
                          cluster_columns = TRUE,
                          clustering_distance_rows = clustering_distance,
                          clustering_distance_columns = clustering_distance,
                          show_row_names = TRUE,
                          show_column_names = TRUE,
                          #bottom_annotation = cann,
                          right_annotation = rann,
                          row_names_gp = gpar(fontsize = font_size),
                          column_names_gp = gpar(fontsize = font_size))
      if(!missing(split)){
        heatmap_args[["row_split"]] = split
        heatmap_args[["column_split"]] = split
      }
      ht <- do.call("Heatmap",heatmap_args)
      draw(ht)
      
      return(list(corr_mat = M,mut_mat = mutmat,plot=ht))
      
    }else{
      # Define a cell function that draws a circle in each cell based on the value
      cell_fun <- function(j, i, x, y, width, height, fill) {
        # Get the current value for the cell.
        value <- M[i, j]
        if (!is.na(value)) {
          # Compute a radius proportional to the magnitude of the value.
          # Here we set the radius in "npc" units; you might adjust the 0.4 factor.
          r <- unit(size_factor * abs(value) / lim, "npc")
          grid.circle(x, y, r = r, gp = gpar(fill = fill, col = fill))
        }
      }
      
      heatmap_args = list(matrix=M,
                          name = "Overlap/\nExclusivity",
                          col = colorRamp2(c(-lim, 0, lim), c("blue", "white", "red")),
                          cell_fun = cell_fun,
                          rect_gp = gpar(fill = NA, col = NA), 
                          na_col = "white",  # Cells that are NA are shown in grey.
                          cluster_rows = TRUE,
                          cluster_columns = TRUE,
                          show_row_names = TRUE,
                          show_column_names = TRUE,
                          #bottom_annotation = cann,
                          right_annotation = rann,
                          row_names_gp = gpar(fontsize=font_size),
                          column_names_gp = gpar(fontsize=font_size))
      if(!missing(split)){
        heatmap_args[["row_split"]] = split
        heatmap_args[["column_split"]] = split
      }
      ht= do.call("Heatmap",heatmap_args)
      # Draw the heatmap.
      draw(ht)
      print(range(M))
      if(return_data){
        return(list(corr_mat = M,mut_mat = mutmat,plot=ht))
      }
      
    }
  }
  else{
    stop("engine must be ComplexHeatmap or ggcorrplot")
  } 
}

