#' @title Total n variants count plot.
#'
#' @description Generate a bar plot visualizing total variant (SSM or SVs) count for selected contigs.
#'
#' @details This function creates a bar plot showing the total number of variants for a selected sample.
#' Convenience parameters for restricting the returned plot are available. For example, with `ssm` (Boolean)
#' you can toggle if the plot will be in respect to SSM (`ssm = TRUE`) or if you wish to count SVs (`ssm = FALSE`).
#' Lastly, this plotting function also have convenient parameters for customizing the returned plot, 
#' e.g `plot_title`, `hide_legend`, and`plot_subtitle` and `snp_colours`. 
#' It is also possible to control what variants are to be counted with `variant_select`. 
#' Default is deletions, insertions and duplications, c("DEL", "DUP", "INS").
#'
#' @param this_maf Parameter with maf like df already loaded into R.
#' @param this_maf_path Parameter with path to external maf like file.
#' @param this_bedpe Parameter with bedpe like df already loaded into R.
#' @param this_bedpe_path Parameter with path to external bedpe like file.
#' @param ssm Set to FALSE to get plotting data from get_combined_sv (SVs). Default value is TRUE (plots SSM retrieved from annotate_cn_by_ssm$maf).
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Subtitle for created plot.
#' @param chr_select vector of chromosomes to be included in plot, defaults to autosomes.
#' @param variant_select Subtypes of SVs to be included in plot, default is DEL, INS and DUP.
#' @param snp_colours Optional vector with colours for SNPs (DNP and TNP).
#' @param hide_legend Set to True to remove legend from plot, default is FALSE.
#' @param log10_y Set to TRUE to force y axis to be in log10.
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
#' fancy_v_count(this_maf = dohh2_maf, this_bedpe = dohh2_bedpe)
#'
fancy_v_count = function(this_maf,
                         this_maf_path = NULL,
                         this_bedpe,
                         this_bedpe_path = NULL,
                         ssm = TRUE,
                         plot_title = "",
                         plot_subtitle = "Variant Count For Selected Contigs",
                         chr_select = paste0("chr", c(1:22)),
                         variant_select = c("DEL", "INS", "DUP"),
                         snp_colours = c("SNP" = "#2B9971", "DNP" = "#993F2B", "TNP" = "#A62656"),
                         hide_legend = FALSE,
                         log10_y = FALSE){
  
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
