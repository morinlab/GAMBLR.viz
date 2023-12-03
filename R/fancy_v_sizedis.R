#' @title Variant size distribution plot
#'
#' @description Generate a violin plot showing variant (SSM or SVs) size distributions for selected contigs.
#'
#' @details Function for plotting variant size distributions. This function takes either a sample ID given to the `this_sample` parameter.
#' In addition, the function can also accept an already loaded MAF or MAF-like object given to the `maf_data` parameter.
#' As a third option, the function can also read a maf from disk (provide path to maf with `maf_path`).
#' A collection of convenient filtering and data subsetting parameters are also available for this function.
#' This plot can also deal with SVs as well as SSM data. To control this, please use the `ssm` parameter. If set to TRUE and if `this_sample` is called,
#' the function gets data with get_ssm_by_sample and if set to FALSE, the function calls `get_manta_sv` to get SV calls for plotting.
#' If the user calls either `maf_data` or `maf_path`, there are a collection of parameters available for specifying
#' the relevant columns in the given data frame (`variant_type_col`, `chhromosome_col`, `start_col`, `end_col`). It is also possible to
#' restrict the returned plot to any given chromosome. This is done with the `chr_select` parameter (default is all autosomes).
#' In addition, plot aesthetics can also be controlled with `plot_title`, `plot_subtitle`, `scale_value`, `log10`, and `trim`.
#' For more info on how to run with these parameters, refer to the parameter descriptions.
#'
#'
#' @param this_sample_id Sample to be plotted.
#' @param maf_data Optional parameter with maf like df already loaded into R.
#' @param maf_path Optional parameter with path to external maf like file.
#' @param ssm Set to FALSE to get plotting data from get_combined_sv (SVs). Default value is TRUE (plots SSM retrieved from annotate_cn_by_ssm$maf).
#' @param projection Genome build for returned variants (only applicable for ssm = FALSE).
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Recommended: 0 (only applicable for ssm = FALSE).
#' @param variant_type_col Index of column holding Variant Type (to be used with either maf_data or maf_path).
#' @param chromosome_col Index of column holding Chromosome (to be used with either maf_data or maf_path).
#' @param start_col Index of column with variant start coordinates (to be used with either maf_data or maf_path).
#' @param end_col Index of column with variant end coordinates (to be used with either maf_data or maf_path).
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Subtitle for created plot.
#' @param scale_value Scale type for violin plot, accepted values are "area", "width", and "count", default is "count.
#' @param log_10 Boolean statement for y-axis, default is TRUE.
#' @param plot_trim If TRUE, trim the tails of the violins to the range of the data. If FALSE (default), don't trim the tails.
#' @param chr_select vector of chromosomes to be included in plot, defaults to autosomes.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import ggplot2 dplyr cowplot
#' @export
#'
#' @examples
#' #plot SSM size distributions:
#' fancy_v_sizedis(this_sample_id = "DOHH-2")
#'
fancy_v_sizedis = function(this_sample_id,
                           maf_data,
                           maf_path = NULL,
                           ssm = TRUE,
                           projection = "grch37",
                           this_seq_type = "genome",
                           min_vaf = 0,
                           variant_type_col = 10,
                           chromosome_col = 5,
                           start_col = 6,
                           end_col = 7,
                           plot_title = paste0(this_sample_id),
                           plot_subtitle = "Variant Size Distribution",
                           scale_value = "width",
                           log_10 = TRUE,
                           plot_trim = FALSE,
                           chr_select = paste0("chr", c(1:22))){

  if(!missing(maf_data)){
    plot_data = maf_data
    plot_data = as.data.frame(plot_data)
    colnames(plot_data)[variant_type_col] = "Variant_Type"
    colnames(plot_data)[chromosome_col] = "Chromosome"
    colnames(plot_data)[start_col] = "Start_Position"
    colnames(plot_data)[end_col] = "End_Position"

  }else if (!is.null(maf_path)){
    plot_data = fread_maf(maf_path)
    plot_data = as.data.frame(plot_data)
    colnames(plot_data)[variant_type_col] = "Variant_Type"
    colnames(plot_data)[chromosome_col] = "Chromosome"
    colnames(plot_data)[start_col] = "Start_Position"
    colnames(plot_data)[end_col] = "End_Position"
  }

  #get maf data for a specific sample.
  if(missing(maf_data) && is.null(maf_path)){
    if(ssm){
      plot_data = get_ssm_by_sample(this_sample_id = this_sample_id, 
                                    this_seq_type = this_seq_type, 
                                    projection = projection)
    }else{
      plot_data = get_manta_sv(these_sample_ids  = this_sample_id, projection = projection, min_vaf = min_vaf) %>%
        dplyr::select(CHROM_A, START_A, END_A, manta_name)

      #get manta results in required format
      plot_data = data.frame(plot_data$CHROM_A, plot_data$START_A, plot_data$END_A, do.call(rbind, strsplit(plot_data$manta_name, split = ":", fixed = TRUE)))

      #rename variables
      names(plot_data)[1] = "Chromosome"
      names(plot_data)[2] = "Start_Position"
      names(plot_data)[3] = "End_Position"
      names(plot_data)[4] = "Variant_Type"

      #filter out translocations and set order of variables
      plot_data = dplyr::filter(plot_data, Variant_Type %in% c("MantaDEL", "MantaDUP")) %>%
        dplyr::select(Chromosome, Start_Position, End_Position, Variant_Type)

      #remove "Manta" from Variant_Type string
      plot_data$Variant_Type = gsub("^.{0,5}", "", plot_data$Variant_Type)
    }
    
    if(nrow(plot_data) == 0){
      stop("No variants found for the selected sample")
    }
  }

  #add chr prefix if missing
  if(!str_detect(plot_data$Chromosome, "chr")[1]){
    plot_data = mutate(plot_data, Chromosome = paste0("chr", Chromosome))
  }

  #read maf into R and select relevant variables and transform to factor.
  plot_data_df = dplyr::select(plot_data, Chromosome, Start_Position, End_Position, Variant_Type) %>%
    mutate_at(vars(Chromosome, Variant_Type), list(factor))

  #calculate variant size
  plot_data_df$Size = plot_data_df$End_Position - plot_data_df$Start_Position

  if(ssm){
    plot_data_df = plot_data_df[plot_data_df$Variant_Type %in% c("DEL", "INS"), ]
    levels(plot_data_df$Size)[levels(plot_data_df$Size) == "0"] = "1"
    plot_data_df[,5][plot_data_df[,5] == 0] <- 1
  }else{
    plot_data_df = plot_data_df[plot_data_df$Variant_Type %in% c("DEL", "DUP"), ]
  }

  plot_data_df$Size = as.integer(plot_data_df$Size)

  #sub-setting maf based on user-defined parameters
  plot_data_df = plot_data_df[plot_data_df$Chromosome %in% chr_select, ]

  p = ggplot(plot_data_df, aes(x = Variant_Type, y = Size, fill = Variant_Type)) +
        labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "Variant Size (bp)") +
        geom_violin(trim = plot_trim, scale = scale_value, color = NA) +
        stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +
        {if(ssm)scale_fill_manual(values = GAMBLR.helpers::get_gambl_colours("indels"))} +
        {if(!ssm)scale_fill_manual(values = GAMBLR.helpers::get_gambl_colours("svs"))} +
        {if(log_10)scale_y_log10()} +
        theme_cowplot() +
        theme(legend.position = "none")

  return(p)
}
