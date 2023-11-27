#' @title Structural Variants Size Plot.
#'
#' @description Generate plot visualizing SV sizes. Subset on variant type, filter on VAF, size etc.
#'
#' @details Plot sample-level SV sizes across selected chromosomes. This function also has a variety of filtering parameters available.
#' For example, it is possible to subset the included variants to a specific VAF threshold with `VAF_cutoff`. The `size_cutoff` is another parameter
#' for filtering the variants on set variant sizes, the default for this parameter is to only include variants of at least 50bp.
#' This function takes either a sample ID (`this_sample_id`) or an already loaded data frame (`maf_data` or a path to a maf-like file with `maf_path`).
#' If `this_sample_id` is called, the function will run [GAMBLR::get_combined_sv] to retrieve SV calls.
#' If either of the `maf` parameters are used, note that it's possible to specify the columns of interest;
#' (`chrom_a_col`, `start_a_col`, `end_a_col` and `variant_type_col`), allowing this function to work with any maf-like data frames.
#' This function also allows the user to customize the returned plot. For more info on how to do this, please refer to the aesthetic
#' parameters; `hide_legend`, `plot_title`, `plot_subtitle`, `adjust_value` and `trim`.
#'
#' @param this_sample_id Sample to be plotted.
#' @param maf_data Optional parameter with copy number df already loaded into R.
#' @param maf_path Optional parameter with path to external cn file.
#' @param chrom_a_col Index of column holding chromosome (to be used with either maf_data or maf_path).
#' @param start_a_col Index of column holding start coordinates (to be used with either maf_data or maf_path).
#' @param end_a_col Index of column holding end coordinates (to be used with either maf_data or maf_path).
#' @param variant_type_col Index of column holding variant type information (to be used with either maf_data or maf_path).
#' @param vaf_cutoff Threshold for filtering variants on VAF (events with a VAF > cutoff will be retained).
#' @param size_cutoff Threshold for filtering variants on size, default is 50bp.
#' @param adjust_value A multiplicate bandwidth adjustment. This makes it possible to adjust the bandwidth while still using the bandwidth estimator. For example, adjust = 1/2 means use half of the default bandwidth.
#' @param trim If FALSE, the default, each density is computed on the full range of the data.
#' @param chr_select Optional argument for subsetting on selected chromosomes, default is all autosomes.
#' @param hide_legend Set to True to remove legend from plot, default is FALSE.
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Subtitle for created plot.
#' @param projection Genomic projection for SVs and circos plot. Accepted values are grch37 and hg38.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import dplyr ggplot2 cowplot stringr
#' @export
#'
#' @examples
#' \dontrun{
#' library(GAMBLR.data)
#' 
#' #build plot sith default parameters
#' fancy_sv_sizedens(this_sample_id = "SP116715")
#'
#' #restrict plot to only chromosome 1 and 2
#' fancy_sv_sizedens(this_sample_id = "SP116715",
#'                   size_cutoff = 0,
#'                   chr_select = c("chr1", "chr2"))
#' }
#'
fancy_sv_sizedens = function(this_sample_id,
                             maf_data,
                             maf_path = NULL,
                             chrom_a_col = 3,
                             start_a_col = 4,
                             end_a_col = 5,
                             variant_type_col = 9,
                             vaf_cutoff = 0,
                             size_cutoff = 50,
                             adjust_value = 1,
                             trim = FALSE,
                             hide_legend = FALSE,
                             chr_select = paste0("chr", c(1:22)),
                             plot_title = paste0(this_sample_id),
                             plot_subtitle = paste0("SV sizes for Manta calls. 
                                                    Dashed line annotates mean 
                                                    variant size.\nVAF cut off: 
                                                    ", vaf_cutoff,", SV size 
                                                    cut off: ", size_cutoff),
                             projection = "grch37"){
  if(!missing(maf_data)){
    svs = maf_data
    svs = as.data.frame(svs)
    colnames(svs)[chrom_a_col] = "CHROM_A"
    colnames(svs)[start_a_col] = "START_A"
    colnames(svs)[end_a_col] = "END_A"
    colnames(svs)[variant_type_col] = "manta_name"

  }else if (!is.null(maf_path)){
    svs = maf_data
    svs = as.data.frame(svs)
    colnames(svs)[chrom_a_col] = "CHROM_A"
    colnames(svs)[start_a_col] = "START_A"
    colnames(svs)[end_a_col] = "END_A"
    colnames(svs)[variant_type_col] = "manta_name"
  }

  #get variants, filter and subset
  if(missing(maf_data) && is.null(maf_path)){
    svs = get_manta_sv(these_sample_ids = this_sample_id, projection = projection, min_vaf = vaf_cutoff) %>%
      dplyr::filter(VAF_tumour > vaf_cutoff) %>%
      dplyr::select(CHROM_A, START_A, END_A, manta_name)
  }

  #split manta_name variable
  svs_df = data.frame( svs$CHROM_A,   svs$START_A, svs$END_A,
                       sub("^(.+?):.*", "\\1", svs$manta_name) )
  
  #rename variables
  names(svs_df)[1] = "chrom"
  names(svs_df)[2] = "start"
  names(svs_df)[3] = "end"
  names(svs_df)[4] = "type"

  #subset df on SV type
  manta_sv = dplyr::filter(svs_df, type %in% c("MantaDEL", "MantaDUP")) %>%
    dplyr::select(chrom, start, end, type)
  
  # check whether enough variants
  manta_sv = mutate( manta_sv, type = factor(type, levels = c("MantaDEL", "MantaDUP")) )
  type_table = table(manta_sv$type)
  if( missing(maf_data) && is.null(maf_path) ){
    k <- gettextf("%i MantaDEL and %i MantaDUP variants were found after filtering by vaf_cutoff.", 
                  type_table["MantaDEL"], type_table["MantaDUP"])
  }else{
    k <- gettextf("%i MantaDEL and %i MantaDUP variants were found.", 
                  type_table["MantaDEL"], type_table["MantaDUP"])
  }
  message(k)
  stopifnot("Plot couldn't be made. At least 2 variants of either type is needed." =
              any(type_table > 1))
  
  #calculate sizes
  manta_sv$size = manta_sv$end - manta_sv$start
  
  #add chr prefix, if missing
  if(!str_detect(manta_sv$chrom, "chr")[1]){
    manta_sv = mutate(manta_sv, chrom = paste0("chr", chrom))
  }

  #subset on selected chromosomes
  manta_sv = manta_sv[manta_sv$chrom %in% chr_select, ]

  #filter out variants < 50 bp
  manta_sv = dplyr::filter(manta_sv, size >= size_cutoff)
  
  # check whether enough variants
  type_table = table(manta_sv$type)
  k <- gettextf("%i MantaDEL and %i MantaDUP variants were left after filtering by size_cutoff.", 
                type_table["MantaDEL"], type_table["MantaDUP"])
  message(k)
  stopifnot("Plot couldn't be made. At least 2 variants of either type is needed." =
              any(type_table > 1))
  
  # groups (MantaDEL or MantaDUP) with only 1 variant are dropped.
  if(type_table["MantaDEL"] == 1){
    manta_sv = filter(manta_sv, type != "MantaDEL")
    message("Warning: At least 2 data points are needed to calculate density estimates. MantaDEL group was dropped because it contains only 1 variant.")
  }
  if(type_table["MantaDUP"] == 1){
    manta_sv = filter(manta_sv, type != "MantaDUP")
    message("Warning: At least 2 data points are needed to calculate density estimates. MantaDUP group was dropped because it contains only 1 variant.")
  }
  manta_sv = mutate(manta_sv, type = droplevels(type))
  
  del_col = GAMBLR.helpers::get_gambl_colours("indels")[[1]]
  dup_col = GAMBLR.helpers::get_gambl_colours("indels")[[2]]

  #plotting
  p = ggplot(manta_sv, aes(x = size, fill = type)) +
        geom_density(alpha = 0.7, color = NA, adjust = adjust_value, trim = trim) +
        labs(title = plot_title, subtitle = plot_subtitle, x = "Size (bp)", y = "") +
        scale_fill_manual(values = c(del_col, dup_col)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        {if(hide_legend)theme(legend.position = "none")} +
        theme_cowplot()

  return(p)
}
