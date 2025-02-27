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
#' dlbcl_meta <- get_gambl_metadata() %>% dplyr::filter(pathology == "DLBCL", seq_type != "mrna")
#' all_coding <- get_all_coding_ssm()
#' bl_fl_dlbcl_meta = get_gambl_metadata() %>% 
#'   dplyr::filter(pathology %in% c("DLBCL","FL","BL"), seq_type != "mrna")
#' all_coding <- get_all_coding_ssm(get_gambl_metadata)
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
#'   font_size = 6,
#'   use_alpha = T,
#'   clustering_distance = "binary",
#'   include_hotspots = F,
#'   return_data = TRUE
#' )
#' Sys.Date()
#' 
#' outs = pretty_mutual_exclusivity(mutmat=outs$mut_mat,
#'   corr_mat = outs$corr_mat,
#'   p_mat = outs$p_mat,
#'   maf_data = all_coding,
#'   genes = fl_bl_dlbcl_genes,
#'   these = bl_fl_dlbcl_meta,
#'   engine = "ComplexHeatmap",
#'   font_size = 6,
#'   use_alpha = T,
#'   size_factor = 0.004,
#'   clustering_distance = "euclidean",
#'   include_hotspots = F
#' )
#' Sys.Date()
pretty_mutual_exclusivity <- function(maf_data,
                                      mut_mat,
                                      corr_mat,
                                      p_mat,
                                      genes,
                                      these_samples_metadata,
                                      q_threshold = 0.1,
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
                                      annotate_by_pathology = TRUE,
                                      show_heatmap_legend = TRUE,
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
    
  } else {
    mutmat = mut_mat
    if (missing(genes)) {
      genes <- colnames(mutmat)
    }
  }
  if(missing(corr_mat) | missing(p_mat)){
    print("calculating correlation")
    mcor <- cor(mutmat)
    cpmat <- cor_pmat(mutmat)
    print("done")
  }else{
    mcor = corr_mat
    cpmat = p_mat
  }
  
  # bonferroni correction
  if (bonferroni) {
    cpmat <- cpmat * nrow(cpmat)^2
  }

  M <- mcor
  M[cpmat > q_threshold] <- 0
  diag(M) <- 0
  M[is.na(M)] <- 0



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
       which="column",
       show_legend = FALSE,
                               col=list(Class=cols,
                                        DLBCL=get_gambl_colours("pos_neg"),
                                        FL = get_gambl_colours("pos_neg"),
                                        BL = get_gambl_colours("pos_neg")))
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
      heatmap_args[["width"]] = unit(width,"cm")
      if (!missing(split)) {
        heatmap_args[["row_split"]] <- split
        heatmap_args[["column_split"]] <- split
      }
      ht <- do.call("Heatmap", heatmap_args)
      draw(ht)
      if(return_data){
        #get the order of rows in the heatmap
        ordered_genes = rownames(M)[row_order(ht)]
        return(list(plot = ht,
                    mut_mat=mutmat,
                    corr_mat = mcor,
                    p_mat = cpmat,
                    original_order = rownames(M) , order = ordered_genes))
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
      heatmap_args[["width"]] = unit(width,"cm")
      if (!missing(split)) {
        heatmap_args[["row_split"]] <- split
        heatmap_args[["column_split"]] <- split
      }
      ht <- do.call("Heatmap", heatmap_args)
      # Draw the heatmap.
      draw(ht)
      print(range(M))
      if (return_data) {
        #get the order of rows in the heatmap
        ordered_genes = rownames(M)[row_order(ht)]
        return(list(plot = ht,
                    mut_mat=mutmat,
                    corr_mat = mcor,
                    p_mat = cpmat,
                    original_order = rownames(M) , order = ordered_genes))
      }
    }
  } else {
    stop("engine must be ComplexHeatmap or ggcorrplot")
  }
}
