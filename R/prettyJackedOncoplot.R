#' Add a top/bottom heatmap to an existing ComplexHeatmap stack
#'
#' @description
#' Builds and inserts a new heatmap (expression- or copy-number–based) into an
#' existing list of `ComplexHeatmap` objects, aligning columns to the sample
#' order of the first heatmap in `existing_heatmaps`. The new heatmap can be
#' placed at the top (default) or bottom of the stack, and either drawn
#' immediately or returned for further composition.
#'
#' @details
#' Exactly one of `all_exp` **or** `cn_mat` must be supplied:
#'
#' - **Expression mode** (`all_exp` provided): selects `genes` from the expression
#'   table (samples in rows, genes in columns), transposes to genes × samples,
#'   centers each gene by its row mean, and constructs a `ComplexHeatmap::Heatmap`
#'   with Euclidean column clustering and no row/column labels.
#'
#' - **Copy-number mode** (`cn_mat` provided): validates the same sample order
#'   against `cn_mat` (samples in rows), reorders `these_samples_metadata` to
#'   match the heatmap’s column order, and builds a CNV heatmap via
#'   `pretty_CN_heatmap()` (external helper), returning its `heatmap_object`.
#'
#' In both modes, the column order is taken from
#' `existing_heatmaps[[1]]@column_names_param$labels`. If samples present in the
#' existing heatmap are missing from the provided data, the function stops with
#' an error.
#'
#' @param these_samples_metadata A data frame of sample-level metadata. Must
#'   contain a column `sample_id` matching the sample IDs/column labels of the
#'   existing heatmap stack. Used for CNV mode and for preserving sample order.
#' @param all_exp Optional data frame or matrix of expression values with
#'   **samples in rows** and **genes in columns**. Provide this for expression
#'   mode; do not provide `cn_mat` simultaneously.
#' @param cn_mat Optional matrix or data frame of copy-number features with
#'   **samples in rows** (row names as sample IDs). Provide this for CNV mode;
#'   do not provide `all_exp` simultaneously.
#' @param genes Character vector of gene symbols/feature names to extract from
#'   `all_exp` (ignored in CNV mode).
#' @param existing_heatmaps A list of `ComplexHeatmap` objects that will be
#'   concatenated with the new heatmap using `%v%`. The first element's column
#'   labels define the sample order.
#' @param height Numeric, height of the new heatmap in centimeters.
#' @param width Numeric, width of the new heatmap in centimeters.
#' @param heatmap_location One of `"top"` (default) or `"bottom"`, indicating
#'   where to insert the new heatmap relative to the existing stack.
#' @param returnEverything Logical; if `FALSE` (default) the function calls
#'   `ComplexHeatmap::draw()` on the composed stack and returns `NULL`. If
#'   `TRUE`, no drawing is performed and a list is returned containing the new
#'   heatmap plus the original `existing_heatmaps` (order depends on
#'   `heatmap_location`).
#'
#' @returns
#' - If `returnEverything = FALSE`: invisibly returns `NULL` and draws the
#'   combined heatmap stack to the active graphics device.
#' - If `returnEverything = TRUE`: a list of heatmap objects. When
#'   `heatmap_location = "top"`, the list is `list(new = <Heatmap>, ...)`
#'   followed by `existing_heatmaps`; when `"bottom"`, it is `c(existing_heatmaps, list(new = <Heatmap>))`.
#'
#' @section Requirements & assumptions:
#' - Column order is taken from `existing_heatmaps[[1]]`; all provided data must
#'   contain exactly those samples (extra samples are ignored; missing samples
#'   trigger an error).
#' - Expression mode centers genes by row mean (no scaling to variance).
#' - CNV mode relies on an external helper `pretty_CN_heatmap()` that should
#'   return a list with a `heatmap_object` element.
#'
#' @examples
#' \dontrun{
#' # Suppose 'stack' is a list of ComplexHeatmap objects already composed for the same samples
#' # Expression mode:
#' prettyJackedOncoplot(
#'   these_samples_metadata = md,                 # must contain sample_id
#'   all_exp = expr_df,                           # samples x genes
#'   cn_mat = NULL,
#'   genes = c("MYC","BCL2","BCL6"),
#'   existing_heatmaps = stack,
#'   height = 6, width = 12,
#'   heatmap_location = "top",
#'   returnEverything = FALSE
#' )
#'
#' # CNV mode:
#' out <- prettyJackedOncoplot(
#'   these_samples_metadata = md,
#'   all_exp = NULL,
#'   cn_mat = cn_matrix,                          # samples x CN features
#'   genes = NULL,
#'   existing_heatmaps = stack,
#'   height = 6, width = 12,
#'   heatmap_location = "bottom",
#'   returnEverything = TRUE
#' )
#' out$new  # the new Heatmap object
#' }
#'
#' @seealso
#' [ComplexHeatmap::Heatmap()], [ComplexHeatmap::draw()], the concatenation
#' operator `%v%` from **ComplexHeatmap**, and `pretty_CN_heatmap()` (package helper).
#'
#' @importFrom dplyr select any_of filter
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom grid unit
#' @keywords visualization heatmap genomics
#' @export

prettyJackedOncoplot = function(these_samples_metadata,
                                all_exp,
                                cn_mat,
                                genes,
                                existing_heatmaps,
                                height,
                                width,
                                heatmap_location = "top",
                                returnEverything = FALSE,
                                cn_args = list()) {

  existing_1 = existing_heatmaps[[1]]
  samp_order = existing_1@column_names_param$labels

  if (!missing(all_exp)) {
    # check for missing samples
    if (any(!samp_order %in% rownames(all_exp))) {
      lmis = length(samp_order[!samp_order %in% rownames(all_exp)])
      print(lmis)
      stop("there are samples in your heatmap that are missing from your expression matrix")
    }

    some_mat = dplyr::select(all_exp, dplyr::any_of(genes)) %>% t()

    new_heatmap = ComplexHeatmap::Heatmap(
      some_mat[, samp_order] - rowMeans(some_mat[, samp_order]),
      cluster_columns = TRUE,
      clustering_distance_columns = "euclidean",
      show_column_names = FALSE,
      height = grid::unit(height, "cm"),
      width  = grid::unit(width,  "cm"),
      show_row_names = FALSE,
      row_dend_side = "right"
    )

  } else if (!missing(cn_mat)) {
    # check for missing samples
    if (any(!samp_order %in% rownames(cn_mat))) {
      lmis = length(samp_order[!samp_order %in% rownames(cn_mat)])
      print(lmis)
      stop("there are samples in your heatmap that are missing from your copy number matrix")
    }

    # ensure metadata is arranged by samp_order (match on sample_id)
    these_samples_metadata = dplyr::filter(these_samples_metadata, sample_id %in% samp_order)
    metadata_df <- these_samples_metadata[match(samp_order, these_samples_metadata$sample_id), , drop = FALSE]
    print(dim(metadata_df))

    # defaults for pretty_CN_heatmap; user can override via cn_args
    cn_defaults <- list(
      cn_state_matrix = cn_mat[metadata_df$sample_id, ],
      scale_by_sample = TRUE,
      these_samples_metadata = metadata_df,
      legend_position = "top",
      rotate = TRUE,
      keep_sample_order = TRUE,
      return_data = TRUE,
      hide_annotations = "chromosome",
      width = width,
      height = height
    )

    # merge defaults with user overrides (cn_args wins on conflicts)
    cn_call <- utils::modifyList(cn_defaults, cn_args, keep.null = TRUE)

    cout <- do.call(pretty_CN_heatmap, cn_call)
    new_heatmap = cout$heatmap_object
  }

  # compose the stack
  if (identical(heatmap_location, "top")) {
    ht_list = NULL %v% new_heatmap
    for (i in seq_along(existing_heatmaps)) {
      ht_list = ht_list %v% existing_heatmaps[[i]]
    }
  } else {
    message("bottom")
    ht_list = NULL
    for (i in seq_along(existing_heatmaps)) {
      ht_list = ht_list %v% existing_heatmaps[[i]]
    }
    ht_list = ht_list %v% new_heatmap
  }

  if (isTRUE(returnEverything)) {
    if (identical(heatmap_location, "bottom")) {
      return(c(existing_heatmaps, list(new = new_heatmap)))
    } else {
      return(c(list(new = new_heatmap), existing_heatmaps))
    }
  } else {
    ComplexHeatmap::draw(
      ht_list,
      heatmap_legend_side = "bottom",
      merge_legend = TRUE
    )
    invisible(NULL)
  }
}



og_prettyJackedOncoplot = function(these_samples_metadata,
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
      print(paste(samp_order[!samp_order %in% rownames(all_exp)],collapse=","))
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
    message("arranging metadata order to match layout of original heatmap")
    these_samples_metadata = filter(these_samples_metadata,sample_id %in% samp_order)
    metadata_df <- these_samples_metadata[order(match(rownames(these_samples_metadata), 
                                                      samp_order)), , drop = FALSE]
    
    cout = pretty_CN_heatmap(cn_mat[metadata_df$sample_id,],
                            scale_by_sample = T,
                            these_samples_metadata = metadata_df,
                            legend_position = "top",
                            rotate = T,
                            keep_sample_order = T,
                            return_data = T,
                            hide_annotations = "chromosome",
                            width=width,
                            height=height
                            )
    new_heatmap = cout$heatmap_object
    #draw(new_heatmap)
    
  }
  if(heatmap_location == "top"){
    ht_list = NULL %v% new_heatmap
    for(i in c(1:length(existing_heatmaps))){
      ht_list = ht_list %v% existing_heatmaps[[i]]
    }

  }else{
    print("bottom")
    ht_list = NULL 
    for(i in c(1:length(existing_heatmaps))){
      ht_list = ht_list %v% existing_heatmaps[[i]]
    }
    ht_list = ht_list %v% new_heatmap

  }
  
  if(returnEverything){
    if(heatmap_location == "bottom"){
        #return(c(existing_heatmaps,list(new=new_heatmap)))
        #return(list(c(existing=existing_heatmaps,new=list(new=new_heatmap))))
        return(list(new=new_heatmap,heatmap_list=ht_list,existing=existing_heatmaps))
    }else{
      #return(c(list(new=new_heatmap),existing_heatmaps))
      return(list(new=new_heatmap,heatmap_list=ht_list,existing=existing_heatmaps))
    }
  }else{
    draw(ht_list,
         heatmap_legend_side = "bottom",
         merge_legend = TRUE)
  }


}
