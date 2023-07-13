#' @title n SNVs per chromosome plot
#'
#' @description Visualizing the number of SNVs per chromosome.
#'
#' @details This function takes on an already loaded maf-like data frame, or a path to the maf file of interest.
#' In addition, the user can also give this function a sample ID and the function will run [GAMBLR::assign_cn_to_ssm]
#' to get data for plotting. If a maf file or data frame is used, the user has the chance to specify what column
#' that holds the Variant Type information (`variant_type_col`), in addition the user can also specify what column
#' in the incoming maf that is corresponding to the chromosome annotations. This function also includes useful subsetting
#' options. For example, `chr_select` allows the user to restrict the plot to specific chromosomes. `include_dnp` is an optional
#' argument (Boolean) for if variants of this subtype should be included or not. The plot can also be restricted to only
#' counting coding mutations (`coding_only`). Flat-file and augmented maf options can be toggled with `from_flatfile`
#' and `use_augmented_maf`. Both are TRUE by default and should rarely be set to FALSE. Lastly, this plotting function
#' also have convenient parameters for customizing the returned plot, e.g `plot_title`, `y_interval`, `hide_legend`, and`plot_subtitle`.
#'
#' @param this_sample_id Sample to be plotted.
#' @param maf_data Optional parameter with maf like df already loaded into R.
#' @param maf_path Optional parameter with path to external maf like file.
#' @param variant_type_col Index of column holding Variant Type (to be used with either maf_data or maf_path).
#' @param chromosome_col Index of column holding Chromosome (to be used with either maf_data or maf_path).
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Subtitle for created plot.
#' @param chr_select vector of chromosomes to be included in plot, defaults to autosomes.
#' @param include_dnp Optional argument for including DNPs. Default is FALSE.
#' @param hide_legend Set to True to remove legend from plot, default is FALSE.
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param from_flatfile If set to true the function will use flat files instead of the database.
#' @param use_augmented_maf Boolean statement if to use augmented maf, default is FALSE.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import ggplot2 dplyr cowplot
#' @export
#'
#' @examples
#' #plot SNVs
#' fancy_snv_chrdistplot(this_sample_id = "HTMCP-01-06-00422-01A-01D")
#'
#' #plot SNVs and DNPs
#' fancy_snv_chrdistplot(this_sample_id = "HTMCP-01-06-00422-01A-01D",
#'                       include_dnp = TRUE,
#'                       plot_subtitle = "SNV + DNP Distribution Per Chromosome")
#'
fancy_snv_chrdistplot = function(this_sample_id,
                                 maf_data,
                                 maf_path = NULL,
                                 variant_type_col = 10,
                                 chromosome_col = 5,
                                 plot_title = paste0(this_sample_id),
                                 plot_subtitle = "SNV Distribution Per Chromosome",
                                 chr_select = paste0("chr", c(1:22)),
                                 include_dnp = FALSE,
                                 hide_legend = FALSE,
                                 coding_only = FALSE,
                                 from_flatfile = TRUE,
                                 use_augmented_maf = TRUE,
                                 this_seq_type = "genome"){

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
    maf = assign_cn_to_ssm(
      this_sample_id = this_sample_id,
      coding_only = coding_only,
      from_flatfile = from_flatfile,
      use_augmented_maf = use_augmented_maf,
      this_seq_type = this_seq_type)$maf
  }

  #add chr prefix if missing
  if(!str_detect(maf$Chromosome[1], "chr")){
    maf = mutate(maf, Chromosome = paste0("chr", Chromosome))
  }

  #subset data frame on snp sub type
  maf_snp = dplyr::filter(maf, Variant_Type == "SNP") %>%
    add_count(Chromosome) %>%
    distinct(Chromosome, .keep_all = TRUE) %>%
    dplyr::select(Chromosome, Variant_Type, n)

  if(!include_dnp){
    #get max number of SNP for the chromosome harboring most variants (for setting y-axis value).
    ymax = max(maf_snp$n)

    #plot
    ggplot(maf_snp, aes(x = Chromosome, y = n)) +
      labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "Count (n)", fill = "") +
      scale_x_discrete(expand = c(0, 0.7), limits = chr_select) +
      geom_bar(position = "stack", stat = "identity", fill = "#2B9971", width = 0.75) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_cowplot() +
      coord_flip()
  }

  else{
    #subset data frame on snp sub type
    maf_dnp = dplyr::filter(maf, Variant_Type == "DNP") %>%
      add_count(Chromosome) %>%
      distinct(Chromosome, .keep_all = TRUE) %>%
      dplyr::select(Chromosome, Variant_Type, n)

    #combine data frames
    maf.count = rbind(maf_snp, maf_dnp)

    #get max number of mutations for the chromosome harboring most variants (for setting y-axis value).
    ymax = max(maf_snp$n) + max(maf_dnp$n)

    #plot
    ggplot(maf.count, aes(x = Chromosome, y = n, fill = Variant_Type)) +
      labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "SNV Count (n)", fill = "") +
      scale_x_discrete(expand = c(0, 0.7), limits = chr_select) +
      geom_bar(position = "stack", stat = "identity", width = 0.75) +
      scale_fill_manual("", values = c("SNP" = "#2B9971", "DNP" = "#993F2B")) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      {if(hide_legend)theme(legend.position = "none")} +
      coord_flip()
  }
}
