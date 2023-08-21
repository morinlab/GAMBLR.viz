#' @title n SNVs per chromosome plot
#'
#' @description Visualizing the number of SNVs per chromosome.
#'
#' @details This function takes on an already loaded maf-like data frame, or a path to the maf file of interest.
#' This function also includes useful subsetting options. For example, `chr_select` allows the user to restrict the plot to specific chromosomes. `include_dnp` is an optional
#' argument (Boolean) for if variants of this subtype should be included or not. Lastly, this plotting function
#' also have convenient parameters for customizing the returned plot, e.g `plot_title`, `hide_legend`, and`plot_subtitle`.
#'
#' @param this_maf Parameter with maf like df already loaded into R.
#' @param this_maf_path Parameter with path to external maf like file.
#' @param variant_type_col Index of column holding Variant Type (to be used with either maf_data or maf_path).
#' @param chromosome_col Index of column holding Chromosome (to be used with either maf_data or maf_path).
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Subtitle for created plot.
#' @param chr_select vector of chromosomes to be included in plot, defaults to autosomes.
#' @param include_dnp Optional argument for including DNPs. Default is FALSE.
#' @param hide_legend Set to True to remove legend from plot, default is FALSE.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import ggplot2 dplyr cowplot
#' @export
#'
#' @examples
#' #get data
#' dohh2_maf = GAMBLR.data::sample_data$grch37$maf %>% dplyr::filter(Tumor_Sample_Barcode == "DOHH-2")
#' 
#' #build plot
#' fancy_snv_chrdistplot(this_maf = dohh2_maf)
#'
fancy_snv_chrdistplot = function(this_maf,
                                 this_maf_path = NULL,
                                 variant_type_col = 10,
                                 chromosome_col = 5,
                                 plot_title = "",
                                 plot_subtitle = "SNV Distribution Per Chromosome",
                                 chr_select = paste0("chr", c(1:22)),
                                 include_dnp = FALSE,
                                 hide_legend = FALSE){

  #get maf
  if(!missing(this_maf)){
    maf = this_maf
  }else if(!is.null(this_maf_path)){
    maf = fread_maf(this_maf_path)
  }else{
    stop("Please provide either a maf file (this_maf) or a path to a maf file (this_maf_path)...")
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
