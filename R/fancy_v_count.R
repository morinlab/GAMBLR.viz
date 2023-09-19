#' @title Total n variants count plot.
#'
#' @description Generate a bar plot visualizing total variant (SSM or SVs) count for selected contigs.
#'
#' @details This function creates a barplot showing the total number of variants for a selected sample.
#' Convenience parameters for restricting the returned plot are available. For example, with `ssm` (Boolean)
#' you can toggle if the plot will be in respect to SSM (`ssm = TRUE`) or if you wish to count SVs (`ssm = FALSE`).
#' In addition, this plot can also accept a variety of incoming data types. Either, you supply the function with a sample ID
#' (`this_sample_id`) and the function will retrieve data using [GAMBLR::assign_cn_to_ssm] or [GAMBLR::get_combined_sv] (depending on how the `ssm` parameter is used).
#' This function also supports a maf or maf-like data frame directly, this is done with `maf_data` or `maf_path`. If data is supplied with either of these parameters,
#' the user can specify what column holds the variant type information as well as chromosome information (`variant_type_col` and `chromosome_col`).
#' Restricting the plot to coding mutations is done with `coding_only = TRUE`. Flat-file and augmented maf options can be toggled with `from_flatfile`
#' and `use_augmented_maf`. Both are TRUE by default and should rarely be set to FALSE. Lastly, this plotting function also have convenient parameters for
#' customizing the returned plot, e.g `plot_title`, `y_interval`, `hide_legend`, and`plot_subtitle` and `snp_colours`. lastly, it is also possible
#' to control what variants are to be counted with `variant_select`. Default is deletions, insertions and duplications, c("DEL", "DUP", "INS"). Not that
#' the variant types specified in this parameter must match with whatever is present in the corresponding `variant_type_col`.
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
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Subtitle for created plot.
#' @param chr_select vector of chromosomes to be included in plot, defaults to autosomes.
#' @param variant_select Subtypes of SVs to be included in plot, default is DEL, INS and DUP.
#' @param snp_colours Optional vector with colours for SNPs (DNP and TNP).
#' @param hide_legend Set to True to remove legend from plot, default is FALSE.
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param log10_y Set to TRUE to force y axis to be in log10.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import ggplot2 dplyr cowplot
#' @export
#'
#' @examples
#' #count all variants for one sample (default parameters)
#' fancy_v_count(this_sample_id = "DOHH-2")
#'
fancy_v_count = function(this_sample_id,
                         maf_data,
                         maf_path = NULL,
                         ssm = TRUE,
                         projection = "grch37",
                         this_seq_type = "genome",
                         min_vaf = 0,
                         variant_type_col = 10,
                         chromosome_col = 5,
                         plot_title = paste0(this_sample_id),
                         plot_subtitle = "Variant Count For Selected Contigs",
                         chr_select = paste0("chr", c(1:22)),
                         variant_select = c("DEL", "INS", "DUP"),
                         snp_colours = c("SNP" = "#2B9971", "DNP" = "#993F2B", "TNP" = "#A62656"),
                         hide_legend = FALSE,
                         coding_only = FALSE,
                         log10_y = FALSE){

  if(!missing(maf_data)){
    maf = maf_data
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col] = "Variant_Type"
    colnames(maf)[chromosome_col] = "Chromosome"

  }else if (!is.null(maf_path)){
    maf = fread_maf(maf_path)
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col] = "Variant_Type"
    colnames(maf)[chromosome_col] = "Chromosome"
  }

  #get maf data for a specific sample.
  if(missing(maf_data) && is.null(maf_path)){
    if(ssm){
      maf = assign_cn_to_ssm(
        this_sample_id = this_sample_id,
        coding_only = coding_only,
        this_seq_type = this_seq_type)$maf
    }else{
      maf = get_manta_sv(these_sample_ids = this_sample_id, projection = projection, min_vaf = min_vaf) %>%
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

  #sub-setting maf based on user-defined parameters
  maf_df = maf_df[maf_df$Chromosome %in% chr_select, ]
  maf_df = maf_df[maf_df$Variant_Type %in% variant_select, ]

  #subset variant data
  sv_count = maf_df %>%
    group_by(Variant_Type) %>%
    summarize(count = n())

  #get colours
  indels_cols = get_gambl_colours("indels")
  colours = append(indels_cols, snp_colours)

  #plot
  p = ggplot(sv_count, aes(x = Variant_Type, y = count, fill = Variant_Type, label = count)) +
    geom_bar(position = "stack", stat = "identity") +
    {if(log10_y)labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "Variants n (log10)", fill = "")} +
    {if(!log10_y)labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "Variants (n)", fill = "")} +
    {if(ssm)scale_fill_manual(values = colours)} +
    {if(!ssm)scale_fill_manual(values = get_gambl_colours("svs"))} +
    geom_text(size = 5, position = position_stack(vjust = 0.5)) +
    {if(log10_y)scale_y_log10(expand = c(0, 0))} +
    {if(!log10_y)scale_y_continuous(expand = c(0, 0))} +
    {if(hide_legend)theme(legend.position = "none")} +
    theme_cowplot()

  return(p)
}
