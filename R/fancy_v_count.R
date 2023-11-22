#' @title Total n variants count plot.
#'
#' @description Generate a bar plot visualizing total variant (SSM or SVs) count for selected contigs.
#'
#' @details This function creates a barplot showing the total number of variants for a selected sample.
#' Convenience parameters for restricting the returned plot are available. For example, with `ssm` (Boolean)
#' you can toggle if the plot will be in respect to SSM (`ssm = TRUE`) or if you wish to count SVs (`ssm = FALSE`).
#' In addition, this plot can also accept a variety of incoming data types. Either, you supply the function with a sample ID
#' (`this_sample_id`) and the function will retrieve data using get_ssm_by_sample or get_manta_sv (depending on how the `ssm` parameter is used).
#' This function also supports a maf or maf-like data frame directly, this is done with `maf_data` or `maf_path`. If data is supplied with either of these parameters,
#' the user can specify what column holds the variant type information as well as chromosome information (`variant_type_col` and `chromosome_col`).
#' Lastly, this plotting function also have convenient parameters for
#' customizing the returned plot, e.g `plot_title`, `y_interval`, `hide_legend`, and`plot_subtitle` and `snp_colours`. lastly, it is also possible
#' to control what variants are to be counted with `variant_select`. Default is deletions, insertions and duplications, c("DEL", "DUP", "INS"). Not that
#' the variant types specified in this parameter must match with whatever is present in the corresponding `variant_type_col`.
#'
#' @param this_sample_id Sample to be plotted.
#' @param maf_data Optional parameter with maf like df already loaded into R.
#' @param maf_path Optional parameter with path to external maf like file.
#' @param ssm Set to FALSE to get plotting data from get_manta_sv (SVs). Default value is TRUE (plots SSM retrieved from get_ssm_by_sample).
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
                         log10_y = FALSE){

  if(!missing(maf_data)){
    plot_data = maf_data
    plot_data = as.data.frame(plot_data)
    colnames(plot_data)[variant_type_col] = "Variant_Type"
    colnames(plot_data)[chromosome_col] = "Chromosome"

  }else if (!is.null(maf_path)){
    plot_data = fread_maf(maf_path)
    plot_data = as.data.frame(plot_data)
    colnames(plot_data)[variant_type_col] = "Variant_Type"
    colnames(plot_data)[chromosome_col] = "Chromosome"
  }

  #get maf data for a specific sample.
  if(missing(maf_data) && is.null(maf_path)){
    if(ssm){
      plot_data = get_ssm_by_sample(these_sample_ids = this_sample_id, 
                                    this_seq_type = this_seq_type, 
                                    projection = projection)
    }else{
      plot_data = get_manta_sv(these_sample_ids = this_sample_id, projection = projection, min_vaf = min_vaf) %>%
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

  #sub-setting maf based on user-defined parameters
  plot_data_df = plot_data_df[plot_data_df$Chromosome %in% chr_select, ]
  plot_data_df = plot_data_df[plot_data_df$Variant_Type %in% variant_select, ]

  #subset variant data
  mutation_count = plot_data_df %>%
    group_by(Variant_Type) %>%
    summarize(count = n())

  #get colours
  indels_cols = GAMBLR.helpers::get_gambl_colours("indels")
  colours = append(indels_cols, snp_colours)

  #plot
  p = ggplot(mutation_count, aes(x = Variant_Type, y = count, fill = Variant_Type, label = count)) +
    geom_bar(position = "stack", stat = "identity") +
    {if(log10_y)labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "Variants n (log10)", fill = "")} +
    {if(!log10_y)labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "Variants (n)", fill = "")} +
    {if(ssm)scale_fill_manual(values = colours)} +
    {if(!ssm)scale_fill_manual(values = GAMBLR.helpers::get_gambl_colours("svs"))} +
    geom_text(size = 5, position = position_stack(vjust = 0.5)) +
    {if(log10_y)scale_y_log10(expand = c(0, 0))} +
    {if(!log10_y)scale_y_continuous(expand = c(0, 0))} +
    {if(hide_legend)theme(legend.position = "none")} +
    theme_cowplot()

  return(p)
}
