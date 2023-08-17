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
#' @param this_bedpe Parameter with maf like df already loaded into R.
#' @param this_bedpe_path Parameter with path to external maf like file.
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
#' @import dplyr ggplot2 cowplot
#' @export
#'
#' @examples
#' #get data
#' dohh2_bedpe = GAMBLR.data::sample_data$grch37$bedpe %>% dplyr::filter(tumour_sample_id == "DOHH-2")
#' 
#' #build plot
#' fancy_sv_sizedens(this_bedpe = dohh2_bedpe)
#'
fancy_sv_sizedens = function(this_bedpe,
                             this_bedpe_path = NULL,
                             vaf_cutoff = 0,
                             size_cutoff = 50,
                             adjust_value = 1,
                             trim = FALSE,
                             hide_legend = FALSE,
                             chr_select = paste0("chr", c(1:22)),
                             plot_title = "SV Sizes",
                             plot_subtitle = paste0("SV sizes for Manta calls. Dashed line annotates mean variant size.\nVAF cut off: ", vaf_cutoff,", SV size cut off: ", size_cutoff)){
  
  if(!missing(this_bedpe)){
    svs = this_bedpe
  }else if(!is.null(this_bedpe_path)){
    svs = fread_maf(this_bedpe_path)
  }else{
    stop("Please provide either a bedpe file (this_bedpe) or a path to a bedpe file (this_bedpe_path)...")
  }

  #filter and subset
  svs = svs %>%
      dplyr::filter(VAF_tumour > vaf_cutoff) %>%
      dplyr::select(CHROM_A, START_A, END_A, manta_name)

  #split manta_name variable
  svs_df = data.frame(svs$CHROM_A, svs$START_A, svs$END_A, do.call(rbind, strsplit(svs$manta_name, split = ":", fixed = TRUE)))

  #rename variables
  names(svs_df)[1] = "chrom"
  names(svs_df)[2] = "start"
  names(svs_df)[3] = "end"
  names(svs_df)[4] = "type"

  #subset df on SV type
  manta_sv = dplyr::filter(svs_df, type %in% c("MantaDEL", "MantaDUP")) %>%
    dplyr::select(chrom, start, end, type)

  #calculate sizes
  manta_sv$size = manta_sv$end - manta_sv$start
  manta_sv$type = as.factor(manta_sv$type)

  #add chr prefix, if missing
  if(!str_detect(manta_sv$chrom, "chr")[1]){
    manta_sv = mutate(manta_sv, chrom = paste0("chr", chrom))
  }

  #subset on selected chromosomes
  manta_sv = manta_sv[manta_sv$chrom %in% chr_select, ]

  #filter out varaints < 50 bp
  manta_sv = dplyr::filter(manta_sv, size >= size_cutoff)
  manta_sv$row_num = seq.int(nrow(manta_sv))

  del_col = get_gambl_colours("indels")[[1]]
  dup_col = get_gambl_colours("indels")[[2]]

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
