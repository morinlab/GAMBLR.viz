#' @title Variant size distribution plot
#'
#' @description Generate a violin plot showing variant (SSM or SVs) size distributions for selected contigs.
#'
#' @details Function for plotting variant size distributions. This function takes either a sample ID given to the `this_sample` parameter.
#' In addition, the function can also accept an already loaded MAF or MAF-like object given to the `maf_data` parameter.
#' As a third option, the function can also read a maf from disk (provide path to maf with `maf_path`).
#' A collection of convenient filtering and data subsetting parameters are also available for this function.
#' For restricting your data (if plotting data retrieved with `this_sample_id`), the user can choose to only plot coding mutations with setting `coding_only` to TRUE.
#' This plot can also deal with SVs as well as SSM data. To control this, please use the `ssm` parameter. If set to TRUE and if `this_sample` is called,
#' the function gets data with annotate_cn_by_ssm and if set to FALSE, the function calls `get_combined_sv` to get SV calls for plotting.
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
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param from_flatfile If set to true the function will use flat files instead of the database.
#' @param use_augmented_maf Boolean statement if to use augmented maf, default is FALSE.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import ggplot2 dplyr cowplot
#' @export
#'
#' @examples
#' #plot SSM size distributions:
#' fancy_v_sizedis(this_sample_id = "HTMCP-01-06-00422-01A-01D")
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
                           chr_select = paste0("chr", c(1:22)),
                           coding_only = FALSE,
                           from_flatfile = TRUE,
                           use_augmented_maf = TRUE){

  if(!missing(maf_data)){
    maf = maf_data
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col] = "Variant_Type"
    colnames(maf)[chromosome_col] = "Chromosome"
    colnames(maf)[start_col] = "Start_Position"
    colnames(maf)[end_col] = "End_Position"

  }else if (!is.null(maf_path)){
    maf = fread_maf(maf_path)
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col] = "Variant_Type"
    colnames(maf)[chromosome_col] = "Chromosome"
    colnames(maf)[start_col] = "Start_Position"
    colnames(maf)[end_col] = "End_Position"
  }

  #get maf data for a specific sample.
  if(missing(maf_data) && is.null(maf_path)){
    if(ssm){
      maf = assign_cn_to_ssm(
        this_sample_id = this_sample_id,
        coding_only = coding_only,
        from_flatfile = from_flatfile,
        use_augmented_maf = use_augmented_maf,
        this_seq_type = this_seq_type)$maf
    }else{
      maf = get_combined_sv(these_sample_ids  = this_sample_id, projection = projection, min_vaf = min_vaf) %>%
        dplyr::select(CHROM_A, START_A, END_A, manta_name)

      #get manta results in required format
      maf = data.frame(maf$CHROM_A, maf$START_A, maf$END_A, do.call(rbind, strsplit(maf$manta_name, split = ":", fixed = TRUE)))

      #rename variables
      names(maf)[1] = "Chromosome"
      names(maf)[2] = "Start_Position"
      names(maf)[3] = "End_Position"
      names(maf)[4] = "Variant_Type"

      #filter out translocations and set order of variables
      maf = dplyr::filter(maf, Variant_Type %in% c("MantaDEL", "MantaDUP")) %>%
        dplyr::select(Chromosome, Start_Position, End_Position, Variant_Type)

      #remove "Manta" from Variant_Type string
      maf$Variant_Type = gsub("^.{0,5}", "", maf$Variant_Type)
    }
  }

  #add chr prefix if missing
  if(!str_detect(maf$Chromosome, "chr")[1]){
    maf = mutate(maf, Chromosome = paste0("chr", Chromosome))
  }

  #read maf into R and select relevant variables and transform to factor.
  maf_df = dplyr::select(maf, Chromosome, Start_Position, End_Position, Variant_Type) %>%
    mutate_at(vars(Chromosome, Variant_Type), list(factor))

  #calculate variant size
  maf_df$Size = maf_df$End_Position - maf_df$Start_Position

  if(ssm){
    maf_df = maf_df[maf_df$Variant_Type %in% c("DEL", "INS"), ]
    levels(maf_df$Size)[levels(maf_df$Size) == "0"] = "1"
    maf_df[,5][maf_df[,5] == 0] <- 1
  }else{
    maf_df = maf_df[maf_df$Variant_Type %in% c("DEL", "DUP"), ]
  }

  maf_df$Size = as.integer(maf_df$Size)

  #sub-setting maf based on user-defined parameters
  maf_df = maf_df[maf_df$Chromosome %in% chr_select, ]

  p = ggplot(maf_df, aes(x = Variant_Type, y = Size, fill = Variant_Type)) +
        labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "Variant Size (bp)") +
        geom_violin(trim = plot_trim, scale = scale_value, color = NA) +
        stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +
        {if(ssm)scale_fill_manual(values = get_gambl_colours("indels"))} +
        {if(!ssm)scale_fill_manual(values = get_gambl_colours("svs"))} +
        {if(log_10)scale_y_log10()} +
        theme_cowplot() +
        theme(legend.position = "none")

  return(p)
}
