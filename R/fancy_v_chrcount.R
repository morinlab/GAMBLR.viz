#' @title n variants per chromosome plot.
#'
#' @description Visualizing variant (SSM or SVs) counts per chromosome.
#'
#' @details Takes a maf data frame (or path to a maf), counts the number of variants per chromosome.
#' Selected chromosomes (`chr_select`) are plotted along the x-axis and the variant counts are represented on the y-axis.
#' This function can plot both Structural Variants (SV) and Simple Shared Motifs (SSM).
#' It plots SVs per default and SSM can be added with setting `ssm = TRUE`.
#' This plot can also be restricted to only show coding mutations. To do so, set `coding_only` to TRUE.
#' In addition, the returned plot can also be superimposed with a sample-specific mean coverage (from [GAMBLR::collate_results]).
#' To do so, set `add_qc_metric` to TRUE. A collection of parameters for customizing the returned plot are also available.
#' e.g `plot_title`, `y_interval`, `hide_legend`, and `plot_subtitle`.
#'
#' @param this_sample_id Sample to be plotted.
#' @param maf_data Optional parameter with maf like df already loaded into R.
#' @param maf_path Optional parameter with path to external maf like file.
#' @param ssm Set to FALSE to get plotting data from [GAMBLR::get_combined_sv] (SVs). Default value is TRUE (plots SSM retrieved from annotate_cn_by_ssm$maf)
#' @param projection Genome build for returned variants (only applicable for ssm = FALSE)
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Recommended: 0 (only applicable for ssm = FALSE).
#' @param variant_type_col Index of column holding Variant Type (to be used with either maf_data or maf_path).
#' @param chromosome_col Index of column holding Chromosome (to be used with either maf_data or maf_path).
#' @param plot_title Title of plot (default to sample ID).
#' @param y_interval Optional parameter for specifying intervals on y-axis.
#' @param hide_legend Set to True to remove legend from plot, default is FALSE.
#' @param plot_subtitle Subtitle for created plot.
#' @param chr_select vector of chromosomes to be included in plot, defaults to autosomes.
#' @param add_qc_metric Boolean statement, if set to TRUE specified QC metric will be added (second y-axis).
#' @param this_seq_type Default is "genome".
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import ggplot2 dplyr cowplot
#' @export
#'
#' @examples
#' #plot ssm
#' fancy_v_chrcount(this_sample_id = "DOHH-2",
#'                  y_interval = 10)
#'
fancy_v_chrcount = function(this_sample_id,
                            maf_data,
                            maf_path = NULL,
                            ssm = TRUE,
                            projection = "grch37",
                            min_vaf = 0,
                            variant_type_col = 10,
                            chromosome_col = 5,
                            plot_title = paste0(this_sample_id),
                            y_interval = 1,
                            hide_legend = FALSE,
                            plot_subtitle = "Variant Count Distribution Per Chromosome",
                            chr_select = paste0("chr", c(1:22)),
                            add_qc_metric = FALSE,
                            this_seq_type = "genome"){

  if(!missing(maf_data)){
    plot_data = maf_data
    plot_data = as.data.frame(plot_data)
    colnames(plot_data)[variant_type_col] = "Variant_Type"
    colnames(plot_data)[chromosome_col] = "Chromosome"

  }else if(!is.null(maf_path)){
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
      plot_data = get_manta_sv(these_sample_ids = this_sample_id, 
                               projection = projection, 
                               min_vaf = min_vaf) %>%
        dplyr::select(CHROM_A, START_A, END_A, manta_name)

      #get manta results in required format
      plot_data = data.frame(plot_data$CHROM_A, plot_data$START_A, plot_data$END_A, do.call(rbind, strsplit(plot_data$manta_name, split = ":", fixed = TRUE)))

      #rename variables
      names(plot_data)[1:4] = c("Chromosome", "Start_Position", "End_Position","Variant_Type")

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

  #convert variables to factors
  plot_data$Variant_Type = as.factor(plot_data$Variant_Type)
  plot_data$Chromosome = as.factor(plot_data$Chromosome)

  #add chr prefix if missing
  if(!str_detect(plot_data$Chromosome[1], "chr")){
    plot_data = mutate(plot_data, Chromosome = paste0("chr", Chromosome))
  }

  #subset data frame on sv sub type
  plot_del = dplyr::filter(plot_data, Variant_Type == "DEL") %>%
    add_count(Chromosome) %>%
    distinct(Chromosome, .keep_all = TRUE) %>%
    dplyr::select(Chromosome, Variant_Type, n)

  if(ssm){
    plot_ins = dplyr::filter(plot_data, Variant_Type == "INS") %>%
      add_count(Chromosome) %>%
      distinct(Chromosome, .keep_all = TRUE) %>%
      dplyr::select(Chromosome, Variant_Type, n)

    #combine data frames
    mut.count = rbind(plot_del, plot_ins)

    #get max number of mutations for the chromosome harboring most variants (for setting y-axis value).
    ymax = max(plot_del$n) + max(plot_ins$n)

  }else{
    plot_dup = dplyr::filter(plot_data, Variant_Type == "DUP") %>%
      add_count(Chromosome) %>%
      distinct(Chromosome, .keep_all = TRUE) %>%
      dplyr::select(Chromosome, Variant_Type, n)

    #combine data frames
    mut.count = rbind(plot_del, plot_dup)

    #get max number of mutations for the chromosome harboring most variants (for setting y-axis value).
    ymax = max(plot_del$n) + max(plot_dup$n)
  }

  if(add_qc_metric){
    #get qc data for selected samples
    sample_df = data.frame(sample_id = this_sample_id)
    qc_metrics = collate_results(sample_table = sample_df, seq_type_filter = seq_type) %>%
      dplyr::select(MeanCorrectedCoverage)
    if(nrow(qc_metrics) < 1){
      message("No QC metrics available for selected sample...")
    }
  }

  #plot
  p = ggplot(mut.count, aes(x = Chromosome, y = n, fill = Variant_Type, label = n)) +
        labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "Variants (n)", fill = "") +
        scale_x_discrete(expand = c(0, 0.58), limits = chr_select) +
        geom_bar(position = "stack", stat = "identity") +
        {if(ssm)scale_fill_manual(values = GAMBLR.helpers::get_gambl_colours("indels"))} +
        {if(!ssm)scale_fill_manual(values = GAMBLR.helpers::get_gambl_colours("svs"))} +
        {if(add_qc_metric)geom_hline(qc_metrics, mapping = aes(yintercept = MeanCorrectedCoverage / 10), linetype = "dashed", group = 2)} +
        {if(!add_qc_metric)scale_y_continuous(expand = c(0, 0), breaks = seq(0, ymax + 2, by = y_interval))} +
        {if(add_qc_metric)scale_y_continuous(expand = c(0, 0), breaks = seq(0, ymax + 2, by = y_interval), sec.axis = sec_axis(~.*10, name = "Mean Corrected Coverage (X)", breaks = seq(0, 100, by = 10)))} +
        theme_cowplot() +
        {if(hide_legend)theme(legend.position = "none")} +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}
