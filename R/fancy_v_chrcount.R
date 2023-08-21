#' @title n variants per chromosome plot.
#'
#' @description Visualizing variant (SSM or SVs) counts per chromosome.
#'
#' @details Takes a maf data frame (or path to a maf), counts the number of variants per chromosome.
#' Selected chromosomes (`chr_select`) are plotted along the x-axis and the variant counts are represented on the y-axis.
#' This function can plot both Structural Variants (SV) and Simple Shared Motifs (SSM).
#' It plots SVs per default and SSM can be added with setting `ssm = TRUE`.
#' In addition, the returned plot can also be superimposed with a sample-specific mean coverage (from [GAMBLR::collate_results]).
#' To do so, set `add_qc_metric` to TRUE. A collection of parameters for customizing the returned plot are also available.
#' e.g `plot_title`, `y_interval`, `hide_legend`, and `plot_subtitle`.
#'
#' @param this_maf Parameter with maf like df already loaded into R.
#' @param this_maf_path Parameter with path to external maf like file.
#' @param this_bedpe Parameter with bedpe like df already loaded into R.
#' @param this_bedpe_path Parameter with path to external bedpe like file.
#' @param collated_results If supplied, the function will superimpose QC metrics to thee returned plot. Preferably the return from [GAMBLR.results::collate_results].
#' @param add_qc_metric Boolean statement, if set to TRUE specified QC metric will be added (second y-axis).
#' @param ssm Seet to FALSE to plot SVs instead of SSM. If set to FALSE, this function expects a bedpe with SVs with `this_bedpe` or an absolute path to such a file with `this_bedpe_path`.
#' @param variant_type_col Index of column holding Variant Type (to be used with either maf_data or maf_path).
#' @param chromosome_col Index of column holding Chromosome (to be used with either maf_data or maf_path).
#' @param plot_title Title of plot (default to sample ID).
#' @param y_interval Optional parameter for specifying intervals on y-axis.
#' @param hide_legend Set to True to remove legend from plot, default is FALSE.
#' @param plot_subtitle Subtitle for created plot.
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
#' fancy_v_chrcount(this_maf = dohh2_maf, this_bedpe = dohh2_bedpe, y_interval = 10)
#'
fancy_v_chrcount = function(this_maf,
                            this_maf_path = NULL,
                            this_bedpe,
                            this_bedpe_path = NULL,
                            collated_results,
                            add_qc_metric = FALSE,
                            ssm = TRUE,
                            plot_title = "",
                            y_interval = 1,
                            hide_legend = FALSE,
                            plot_subtitle = "Variant Count Distribution Per Chromosome",
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

  #convert variables to factors
  maf$Variant_Type = as.factor(maf$Variant_Type)
  maf$Chromosome = as.factor(maf$Chromosome)

  #add chr prefix if missing
  if(!str_detect(maf$Chromosome[1], "chr")){
    maf = mutate(maf, Chromosome = paste0("chr", Chromosome))
  }

  #subset data frame on sv sub type
  maf_del = dplyr::filter(maf, Variant_Type == "DEL") %>%
    add_count(Chromosome) %>%
    distinct(Chromosome, .keep_all = TRUE) %>%
    dplyr::select(Chromosome, Variant_Type, n)

  if(ssm){
    maf_ins = dplyr::filter(maf, Variant_Type == "INS") %>%
      add_count(Chromosome) %>%
      distinct(Chromosome, .keep_all = TRUE) %>%
      dplyr::select(Chromosome, Variant_Type, n)

    #combine data frames
    maf.count = rbind(maf_del, maf_ins)

    #get max number of mutations for the chromosome harboring most variants (for setting y-axis value).
    ymax = max(maf_del$n) + max(maf_ins$n)

  }else{
    maf_dup = dplyr::filter(maf, Variant_Type == "DUP") %>%
      add_count(Chromosome) %>%
      distinct(Chromosome, .keep_all = TRUE) %>%
      dplyr::select(Chromosome, Variant_Type, n)

    #combine data frames
    maf.count = rbind(maf_del, maf_dup)

    #get max number of mutations for the chromosome harboring most variants (for setting y-axis value).
    ymax = max(maf_del$n) + max(maf_dup$n)
  }

  if(!missing(collated_results)){
    #get qc data for selected samples
    qc_metrics = collated_results %>%
      dplyr::select(MeanCorrectedCoverage)
  }

  #plot
  p = ggplot(maf.count, aes(x = Chromosome, y = n, fill = Variant_Type, label = n)) +
        labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "Variants (n)", fill = "") +
        scale_x_discrete(expand = c(0, 0.58), limits = chr_select) +
        geom_bar(position = "stack", stat = "identity") +
        {if(ssm)scale_fill_manual(values = get_gambl_colours("indels"))} +
        {if(!ssm)scale_fill_manual(values = get_gambl_colours("svs"))} +
        {if(add_qc_metric)geom_hline(qc_metrics, mapping = aes(yintercept = MeanCorrectedCoverage / 10), linetype = "dashed", group = 2)} +
        {if(!add_qc_metric)scale_y_continuous(expand = c(0, 0), breaks = seq(0, ymax + 2, by = y_interval))} +
        {if(add_qc_metric)scale_y_continuous(expand = c(0, 0), breaks = seq(0, ymax + 2, by = y_interval), sec.axis = sec_axis(~.*10, name = "Mean Corrected Coverage (X)", breaks = seq(0, 100, by = 10)))} +
        theme_cowplot() +
        {if(hide_legend)theme(legend.position = "none")} +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}
