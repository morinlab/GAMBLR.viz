#' @title Variant size distribution plot
#'
#' @description Generate a violin plot showing variant (SSM or SVs) size distributions for selected contigs.
#'
#' @details Function for plotting variant size distributions. This function takes either a dataf frame with a maf like objecct (`this_maf`),
#' or an absolute pathh to such a file. In addition, if `ssm =. TRUE` the user can also provide a bedpe (`this_bedpe` or `this_bedpe_path`) with SVs to be inlcuded in the returned plot.
#' For more information on how to run this function, see parameter descriptions and examples.
#'
#'
#' @param this_maf Parameter with maf like df already loaded into R.
#' @param this_maf_path Parameter with path to external maf like file.
#' @param this_bedpe Parameter with bedpe like df already loaded into R.
#' @param this_bedpe_path Parameter with path to external bedpe like file.
#' @param ssm Set to FALSE to get plotting data from get_combined_sv (SVs). Default value is TRUE (plots SSM retrieved from annotate_cn_by_ssm$maf).
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
#' #get data
#' dohh2_bedpe = GAMBLR.data::sample_data$grch37$bedpe %>% dplyr::filter(tumour_sample_id == "DOHH-2")
#' dohh2_maf = GAMBLR.data::sample_data$grch37$maf %>% dplyr::filter(Tumor_Sample_Barcode == "DOHH-2")
#' 
#' #build plot
#' fancy_v_sizedis(this_maf = dohh2_maf, this_bedpe = dohh2_bedpe)
#'
fancy_v_sizedis = function(this_maf,
                           this_maf_path = NULL,
                           this_bedpe,
                           this_bedpe_path = NULL,
                           ssm = TRUE,
                           plot_title = "",
                           plot_subtitle = "Variant Size Distribution",
                           scale_value = "width",
                           log_10 = TRUE,
                           plot_trim = FALSE,
                           chr_select = paste0("chr", c(1:22))){

  if(ssm){
    if(!missing(this_maf)){
      maf = this_maf
    }else if(!is.null(this_maf_path)){
      maf = fread_maf(this_maf_path)
    }else{
      stop("Please provide either a maf file (this_maf) or a path to a maf file (this_maf_path)...")
    }  
  }else{
    if(!missing(this_bedpe)){
      maf = this_bedpe
    }else if(!is.null(this_bedpe_path)){
      maf = fread_maf(this_bedpe_path)
    }else{
      stop("Please provide either a bedpe file (this_bedpe) or a path to a bedpe file (this_bedpe_path)...")
    }
    #get manta results in required format
    maf = data.frame(maf$CHROM_A, maf$START_A, maf$END_A, do.call(rbind, strsplit(maf$manta_name, split = ":", fixed = TRUE)))
    
    #rename variables
    names(maf)[1:4] = c("Chromosome", "Start_Position", "End_Position","Variant_Type")
    
    #filter out translocations and set order of variables
    maf = dplyr::filter(maf, Variant_Type %in% c("MantaDEL", "MantaDUP")) %>%
      dplyr::select(Chromosome, Start_Position, End_Position, Variant_Type)
    
    #remove "Manta" from Variant_Type string
    maf$Variant_Type = gsub("^.{0,5}", "", maf$Variant_Type)
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
