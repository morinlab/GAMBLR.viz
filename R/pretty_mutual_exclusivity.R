#' Pretty mutual exclusivity plot
#'
#' @param maf_data
#' @param mutmat
#' @param genes
#' @param these_samples_metadata
#' @param q_threshold
#' @param exclude_insignificant_genes
#' @param engine
#' @param font_size
#' @param use_alpha
#' @param clustering_distance
#' @param gene_anno_df
#' @param metadataBarHeight
#' @param metadataBarFontsize
#' @param legend_direction
#' @param size_factor
#' @param split
#' @param return_data Default False
#' @param include_silent Default False
#' @param include_hotspots Default False
#' @param review_hotspots Default False
#' @param bonferroni Default False
#'
#' @import ggcorrplot
#'
#' @returns
#' @export
#'
#' @examples
#' library(GAMBLR)
#' 
#' bl_fl_dlbcl_meta = get_gambl_metadata() %>% 
#'   dplyr::filter(pathology %in% c("DLBCL","FL","BL"), seq_type != "mrna")
#' dlbcl_meta = dplyr::filter(bl_fl_dlbcl_meta,pathology=="DLBCL")
#' 
#' all_coding <- get_all_coding_ssm(bl_fl_dlbcl_meta)
#' 
#' \dontrun{
#' lymphgens = get_lymphgen(flavour = "no_cnvs.no_sv.with_A53")
#' lg_feats = lymphgens$feature_annotation
#' lg_genes = unique(lg_feats$Feature)
#' 
#' pretty_mutual_exclusivity(
#'    maf_data = all_coding,
#'    genes = lg_genes,
#'    these = dlbcl_meta,
#'    size_factor =  0.007,
#'    engine = "ComplexHeatmap",
#'    font_size = 6,
#'    use_alpha = F,
#'    clustering_distance = "binary",
#'    include_hotspots = T
#' ) 
#' }
#' fl_bl_dlbcl_genes = dplyr::filter(GAMBLR.data::lymphoma_genes,
#'   FL_Tier == 1 | BL_Tier == 1 | DLBCL_Tier ==1) %>%
#'   pull(Gene)
#' 
#' # because the first steps of this are slow we can
#' # store the output matrix as a shortcut for subsequent runs
#' Sys.Date()
#' 
#' outs = pretty_mutual_exclusivity(
#'   maf_data = all_coding,
#'   genes = fl_bl_dlbcl_genes,
#'   these = bl_fl_dlbcl_meta,
#'   engine = "ComplexHeatmap",
#'   font_size = 5,
#'   use_alpha = T,
#'   clustering_distance = "binary",
#'   include_hotspots = F,
#'   return_data = TRUE
#' )
#' Sys.Date()
#' 
#' outs = pretty_mutual_exclusivity(mut_mat=outs$mut_mat,
#'   corr_mat = outs$corr_mat,
#'   p_mat = outs$p_mat,
#'   maf_data = all_coding,
#'   genes = fl_bl_dlbcl_genes,
#'   these = bl_fl_dlbcl_meta,
#'   engine = "ComplexHeatmap",
#'   font_size = 5,
#'   use_alpha = T,
#'   size_factor = 0.004,
#'   clustering_distance = "euclidean",
#'   include_hotspots = F
#' )
#' Sys.Date()
#' 
#' outs = pretty_mutual_exclusivity(
#'   p_mat = outs$p_mat,
#'   maf_data = all_coding,
#'   genes = fl_bl_dlbcl_genes,
#'   these = dlbcl_meta,
#'   engine = "ComplexHeatmap",
#'   font_size = 5,
#'   use_alpha = T,
#'   size_factor = 0.004,
#'   clustering_distance = "euclidean",
#'   legend_direction = "vertical",
#'   width = 15,
#'   include_hotspots = F)
#'
pretty_mutual_exclusivity <- function(maf_data,
                                      mut_mat,
                                      cn_mat,
                                      corr_mat,
                                      p_mat,
                                      min_mutation_percent=2,
                                      genes,
                                      these_samples_metadata,
                                      q_threshold = 0.05,
                                      drop_positive_correlations = FALSE,
                                      exclude_insignificant_genes = TRUE,
                                      engine = "ggcorrplot",
                                      font_size = 7,
                                      use_alpha = FALSE,
                                      clustering_distance = "binary",
                                      gene_anno_df,
                                      size_factor = 0.01,
                                      split,
                                      return_data = FALSE,
                                      include_silent = F,
                                      include_hotspots = F,
                                      review_hotspots = F,
                                      bonferroni = FALSE, 
                                      verbose = FALSE,
                                      metadataBarHeight = 3,
                                      metadataBarFontsize = 4,
                                      legend_direction = "horizontal",
                                      annotate_by_pathology = TRUE,
                                      show_heatmap_legend = TRUE,
                                      cut_k,
                                      width = 10) {
  print("STARTING")
  if (missing(mut_mat)) {
    if (missing(genes)) {
      stop("genes is a required argument")
    }
    if(verbose){
      print("generating mutation matrix and will return it for convenience")
    }
    mutmat <- get_coding_ssm_status(
      gene_symbols = genes,
      include_hotspots = include_hotspots,
      these = these_samples_metadata,
      include_silent = include_silent,
      maf_data = maf_data,
      review_hotspots = review_hotspots
    )
    mutmat <- column_to_rownames(mutmat, "sample_id")
    if(verbose){
      print(dim(mutmat))
    }
    if(!missing(cn_mat)){
      print(dim(mutmat))
      print(dim(cn_mat))
      mutmat = bind_cols(cn_mat[rownames(cn_mat) %in% rownames(mutmat), ],mutmat[rownames(mutmat) %in% rownames(cn_mat),])
      print(dim(mutmat))
      #return(mutmat)
    }
  } else {
    if(missing(cn_mat)){
      mutmat = mut_mat
    }else{
      mutmat = bind_cols(cn_mat,mutmat)
      print(dim(mutmat))
    }
    
    if (missing(genes)) {
      genes <- colnames(mutmat)
    }
  }
  mut_percents = 100*colSums(mutmat)/nrow(mutmat)
  keep_g = names(which(mut_percents > min_mutation_percent))
  print(keep_g)
  drop_g = names(which(!mut_percents > min_mutation_percent))
  print(drop_g)
  mutmat = mutmat[,keep_g]
  if(missing(corr_mat) | missing(p_mat)){
    print("calculating correlation")
    
    mcor <- cor(mutmat)
    cpmat <- cor_pmat(mutmat)
    print("done")
  }else{
    mcor = corr_mat[keep_g,keep_g]
    cpmat = p_mat[keep_g,keep_g]
  }
  
  

  M <- mcor
  M[cpmat > q_threshold] <- 0
  diag(M) <- 0
  diag(cpmat) <- 1
  M[is.na(M)] <- 0
  best_p = apply(cpmat,1,function(x){min(x)})
  best_p = best_p * nrow(cpmat)^2
  print(best_p)
  print(table(best_p < q_threshold))
  # kick out insignificant genes in both dimensions
  if(exclude_insignificant_genes){
    good_g = names(which(best_p< q_threshold))
  #print(good_g)
    M = M[good_g,good_g]
    cpmat = cpmat[good_g,good_g]
  }
  drop_pos = TRUE
  if(drop_positive_correlations){
    M[M>0]=0
    print(rowSums(M))
    not_empty = rowSums(M)<0
    M = M[not_empty,not_empty]
    
  }
  print(dim(M))
  gene_lowest  = apply(M,1,function(x){min(x)})
  gene_highest = apply(M,1,function(x){max(x)})
  
    # bonferroni correction
  if (bonferroni) {
    cpmat <- cpmat * nrow(cpmat)^2
  }
  if (engine == "ggcorrplot") {
    print("GGCORRPLOT")
    pp <- ggcorrplot(M,
      method = "circle",
      # size=1,
      # hc.order = TRUE,
      # type="upper",
      lab = FALSE
    ) +
      scale_fill_gradient2(
        low = "blue", mid = "white", high = "red", midpoint = 0,
        # breaks = c(-lim, round(-lim/2), 0, round(lim/2), lim),
        name = "Odds Ratio"
      ) +
      ggtitle("Mutation Correlation and Exclusivity (pAdj < 0.1)") +
      theme_minimal() +
      theme(legend.position = "none") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    if(return_data){
      return(list(plot = pp,mut_mat=mutmat,corr_mat = mcor, p_mat = cpmat))
    }
    
  } else if (engine == "ComplexHeatmap") {
    ######################################################
    ################## ComplexHeatmap ####################
    ######################################################
    # Determine the maximum absolute value in M for scaling:
    
    lim <- max(abs(M), na.rm = TRUE)
    print("gene_anno_df check")
    if (missing(gene_anno_df)) {

      if(annotate_by_pathology){
        FL_genes = filter(GAMBLR.data::lymphoma_genes,FL_Tier==1) %>% pull(Gene)
        BL_genes = filter(GAMBLR.data::lymphoma_genes,BL_Tier==1) %>% pull(Gene)
        DLBCL_genes = filter(GAMBLR.data::lymphoma_genes,DLBCL_Tier==1) %>% pull(Gene)

        
        gene_anno_df = data.frame(Feature=unique(c(FL_genes,BL_genes,DLBCL_genes))) %>%
          dplyr::filter(Feature %in% colnames(M)) %>%
          mutate(DLBCL=ifelse(Feature %in% DLBCL_genes,1,0)) %>%
          mutate(FL=ifelse(Feature %in% FL_genes,1,0)) %>%
          mutate(BL=ifelse(Feature %in% BL_genes,1,0))
        missing_anno <- rownames(M)[!rownames(M) %in% gene_anno_df$Feature]
        if (length(missing_anno) > 0) {
          missing_df <- data.frame(Feature = missing_anno, DLBCL = 0,FL=0,BL=0)
          gene_anno <- bind_rows(gene_anno_df, missing_df) %>%
            distinct(Feature, .keep_all = TRUE) %>%
            column_to_rownames("Feature")
        } else {
        gene_anno <- gene_anno_df %>%
          distinct(Feature, .keep_all = TRUE) %>%
          column_to_rownames("Feature")
        }
        gene_anno <- gene_anno[rownames(M), ]
        rann <- HeatmapAnnotation(
          simple_anno_size = unit(metadataBarHeight, "mm"),
          annotation_name_gp = gpar(fontsize = metadataBarFontsize),
          show_legend = FALSE,
          annotation_name_side = "top",
          df = gene_anno, which = "row",
          col = list(
            
            DLBCL = c("1" = unname(get_gambl_colours("pathology")["DLBCL"]), "0" = "white"),
            FL = c("1" = unname(get_gambl_colours("pathology")["FL"]), "0" = "white"),
            BL = c("1" = unname(get_gambl_colours("pathology")["BL"]), "0" = "white")
          )
        )
        cann <- HeatmapAnnotation(
          simple_anno_size = unit(metadataBarHeight, "mm"),
          annotation_name_gp = gpar(fontsize = metadataBarFontsize),
          df = gene_anno,
          show_legend = FALSE,
          which = "column",
          col = list(
            
            DLBCL = c("1" = unname(get_gambl_colours("pathology")["DLBCL"]), "0" = "white"),
            FL = c("1" = unname(get_gambl_colours("pathology")["FL"]), "0" = "white"),
            BL = c("1" = unname(get_gambl_colours("pathology")["BL"]), "0" = "white")
          )
        )
      }
    }else{
      gene_anno_df <- gene_anno_df %>%
        distinct(Feature, .keep_all = TRUE) %>%
        # select(Feature,Class) %>%
        filter(Feature %in% rownames(M))

      missing_anno <- rownames(M)[!rownames(M) %in% gene_anno_df$Feature]


      if (length(missing_anno) > 0) {
        missing_df <- data.frame(Feature = missing_anno, Class = "NA")
        gene_anno <- bind_rows(gene_anno_df, missing_df) %>%
          distinct(Feature, .keep_all = TRUE) %>%
          column_to_rownames("Feature")
      } else {
        gene_anno <- gene_anno_df %>%
          distinct(Feature, .keep_all = TRUE) %>%
          column_to_rownames("Feature")
      }

      cols <- c(get_gambl_colours("lymphgen"),
        get_gambl_colours("genetic_subgroup"),
        "NA" = "grey"
      )
      gene_anno <- gene_anno[rownames(M), ]
      rann <- HeatmapAnnotation(
        simple_anno_size = unit(metadataBarHeight, "mm"),
        annotation_name_gp = gpar(fontsize = metadataBarFontsize),
        annotation_name_side = "top",
        df = gene_anno, which = "row",
        show_legend=FALSE,
        col = list(
          Class = cols,
          DLBCL = c("1" = unname(get_gambl_colours("pathology")["DLBCL"]), "0" = "white"),
          FL = c("1" = unname(get_gambl_colours("pathology")["FL"]), "0" = "white"),
          BL = c("1" = unname(get_gambl_colours("pathology")["BL"]), "0" = "white")
        )
      )
       cann = HeatmapAnnotation(df=gene_anno,
         simple_anno_size = unit(metadataBarHeight, "mm"),
         annotation_name_gp = gpar(fontsize = metadataBarFontsize),
         which="column",
         show_legend = FALSE,
         col = list(
          Class = cols,
          DLBCL = c("1" = unname(get_gambl_colours("pathology")["DLBCL"]), "0" = "white"),
          FL = c("1" = unname(get_gambl_colours("pathology")["FL"]), "0" = "white"),
          BL = c("1" = unname(get_gambl_colours("pathology")["BL"]), "0" = "white")
        ))
    }


    if (use_alpha) {
      # Define your cell function with transparency added to the fill color
      cell_fun <- function(j, i, x, y, width, height, fill) {
        # Get the current value for the cell.
        value <- M[i, j]
        if (!is.na(value)) {
          # Compute a radius proportional to the magnitude of the value.
          # Here we set the radius in "npc" units; you might adjust the factor as needed.
          r <- unit(size_factor * abs(value) / lim, "npc")
          #rint(paste(value,r))
          # Adjust the fill color to include transparency (alpha = 0.5 means 50% opaque)
          grid.circle(x, y, r = r, gp = gpar(fill = adjustcolor(fill, alpha.f = 0.5), col = fill))
          #grid.text(sprintf("%.1f", M[i, j]), x, y, gp = gpar(fontsize = 10))
        }
      }

      # Create the heatmap using the modified cell_fun
      if (missing(gene_anno_df)) {
        rann <- NULL
      }

      heatmap_args <- list(
        matrix = M,
        name = "Overlap/\nExclusivity",
        col = colorRamp2(c(-lim, 0, lim), c("blue", "white", "red")),
        cell_fun = cell_fun,
        rect_gp = gpar(fill = NA, col = NA),
        na_col = "white", # Cells that are NA are shown in white.
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        clustering_distance_rows = clustering_distance,
        clustering_distance_columns = clustering_distance,
        show_row_names = TRUE,
        show_column_names = TRUE,
        bottom_annotation = cann,
        right_annotation = rann,
        row_names_gp = gpar(fontsize = font_size),
        column_names_gp = gpar(fontsize = font_size)
      )
      if(!missing(cut_k)){
        heatmap_args["row_km"] = cut_k
      }
      heatmap_args[["width"]] = unit(width,"cm")
      heatmap_args[['heatmap_legend_param']] = list(direction=legend_direction)
      if (!missing(split)) {
        heatmap_args[["row_split"]] <- split
        heatmap_args[["column_split"]] <- split
      }
      ht <- do.call("Heatmap", heatmap_args)
      draw(ht)
      if(return_data){
        #get the order of rows in the heatmap
        rord = row_order(ht)
        ordered_gene_ind = c()
        if("list" %in% class(rord)){
          for(clust in names(rord)){
            ordered_gene_ind = c(ordered_gene_ind,rord[[clust]])
          }
        }else{
          ordered_gene_ind = rord
        }

        ordered_genes = rownames(M)[ordered_gene_ind]
        return(list(plot = ht,
                    M_filt = M,
                    mut_mat=mutmat[,rownames(M)],
                    corr_mat = mcor,
                    p_mat = cpmat,
                    final_order = ordered_genes,
                    original_order = rownames(M),
                    gene_lowest_corr = sort(gene_lowest),
                    gene_highest_corr = sort(gene_highest)) 
                    )
                 
      }
    } else {
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
      if (missing(gene_anno_df)) {
        rann <- NULL
      }
      heatmap_args <- list(
        matrix = M,
        name = "Overlap/\nExclusivity",
        col = colorRamp2(c(-lim, 0, lim), c("blue", "white", "red")),
        cell_fun = cell_fun,
        rect_gp = gpar(fill = NA, col = NA),
        na_col = "white", # Cells that are NA are shown in grey.
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        bottom_annotation = cann,
        right_annotation = rann,
        row_names_gp = gpar(fontsize = font_size),
        column_names_gp = gpar(fontsize = font_size),
        show_heatmap_legend = show_heatmap_legend

      )
      if(!missing(cut_k)){
        heatmap_args["row_km"] = cut_k
      }
      heatmap_args[["width"]] = unit(width,"cm")
      if (!missing(split)) {
        heatmap_args[["row_split"]] <- split
        heatmap_args[["column_split"]] <- split
      }
      heatmap_args[['heatmap_legend_param']] = list(direction=legend_direction)
      ht <- do.call("Heatmap", heatmap_args)
      # Draw the heatmap.
      draw(ht)
      print(range(M))
      if (return_data) {
        #get the order of rows in the heatmap
        rord = row_order(ht)
        ordered_gene_ind = c()
        if("list" %in% class(rord)){
          for(clust in names(rord)){
            ordered_gene_ind = c(ordered_gene_ind,rord[[clust]])
          }
        }else{
          ordered_gene_ind = rord
        }
        
        ordered_genes = rownames(M)[ordered_gene_ind]
        return(list(plot = ht,
                    M_filt = M,
                    mut_mat=mutmat[,rownames(M) ],
                    corr_mat = mcor,
                    p_mat = cpmat,
                    original_order = rownames(M) ,
                    final_order = ordered_genes,
                    gene_lowest_corr = sort(gene_lowest)))
      }
    }
  } else {
    stop("engine must be ComplexHeatmap or ggcorrplot")
  }
}
