#' @title Plot Alignment Metrics
#'
#' @description Visualize (stacked barplot) genomic read-subsets (metrics) across a selection of samples.
#'
#' @details This function is available for plotting relevant alignment metrics (read-subsets) for a selection of samples. Per default, this plot returns the following read-metrics;
#' total n reads, total n uniquely mapped reads, total n duplicated reads. This plot can also be superimposed with read metrics from additional samples,
#' allowing for easy comparisons between different sample populations. To run this function, simply specify the sample IDs you are interested in with `these_sample_ids`.
#' This parameter expects a data frame with sample IDs in the first column. Optionally, the user can also provide an already subset (with the sample IDS of interest)
#' metadata table with `these_samples_metadata`. For adding a comparison group to the returned plot, simply give another cohort/set of samples to the `comparison_group` parameter.
#' Similarly to `these_sample_ids`, this parameter also expects a data frame with sample IDs in the first column. In addition, this plot can also add additional read-metrics such as
#' mean values for all plotted metrics and corrected coverage. To enable these features, simply set `add_mean` and `add_corrected_coverage` to TRUE (default).
#'
#' @param these_sample_ids Data frame with sample IDs (to be plotted) in the first column.
#' @param metadata Optional argument, used to derive sample IDs if sample_table is Null.
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process.
#' @param comparison_group Optional argument for plotting mean alignment metrics. Default is plotting the mean for samples provided. This parameter takes a vector of sample IDs.
#' @param seq_type Subset qc metrics to a specific seq_type, default is genome.
#' @param add_mean Set to TRUE to superimpose mean values of plotted variables. Default is TRUE.
#' @param add_corrected_coverage Set to TRUE to add corrected coverage for selected samples.
#' @param keep_cohort If no df with sample IDs is supplied (these_sample_ids = NULL) the function returns metadata and subsets on selected cohorts.
#' @param keep_pathology If no df with sample IDs is supplied (these_sample_ids = NULL) the function returns metadata and subsets on selected pathology.
#' @param this_color_palette Optional parameter that holds the selected colours for the plotted bars.
#' @param plot_sub Optional parameter, add a subtitle to the alignment metric plot.
#'
#' @return A plot as a ggplot object (grob).
#' 
#' @import ggplot2 GAMBLR.helpers dplyr tidyr
#' @export
#'
#' @examples
#' #load packages
#' \dontrun{
#' suppressPackageStartupMessages(library(dplyr))
#' suppressPackageStartupMessages(library(GAMBLR.data))
#'
#' #get sample IDs for available genome samples
#' genome_collated = collate_results(seq_type_filter = "genome") %>%
#'   pull(sample_id)
#'
#' #subset the collated samples on BL samples
#' my_samples = get_gambl_metadata() %>%
#'   dplyr::filter(sample_id %in% genome_collated) %>%
#'   dplyr::filter(pathology == "BL") %>% pull(sample_id)
#'
#' fancy_alignment_plot(these_sample_ids = my_samples)
#' }
fancy_alignment_plot = function(these_sample_ids,
                                metadata,
                                these_samples_metadata,
                                comparison_group,
                                seq_type = "genome",
                                add_mean = TRUE,
                                add_corrected_coverage = TRUE,
                                keep_cohort,
                                keep_pathology,
                                this_color_palette = c("TotalReads" = "#3D405B",
                                                       "TotalUniquelyMapped" = "#81B29A",
                                                       "TotalDuplicatedreads" = "#E07A5F"),
                                plot_sub = ""){
  .Deprecated(msg = "fancy_alignment_plot() is deprecated and will be removed in a future version of GAMBLR.viz.")

  #get gambl metadata (if not supplied)
  if(missing(metadata)){
    this_meta = GAMBLR.helpers::handle_metadata(this_seq_type = seq_type)
  }else{
    this_meta = metadata
  }

  if(!missing(these_samples_metadata)){
    these_sample_ids = dplyr::select(these_samples_metadata, sample_id) %>%
      as.data.frame(strings.as.factors = FALSE)
  }


  #filter metadata on selected cohort/pathology
  if(missing(these_sample_ids)){
    if(!missing(keep_cohort) && missing(keep_pathology)){
      these_sample_ids = dplyr::filter(this_meta, cohort == keep_cohort)
    }

    if(!missing(keep_pathology) && missing(keep_cohort)){
      these_sample_ids = dplyr::filter(this_meta, pathology == keep_pathology)
    }

    if(!missing(keep_cohort) && !missing(keep_pathology)){
      these_sample_ids = dplyr::filter(this_meta, pathology == keep_pathology, cohort == keep_cohort)
    }

    if(missing(keep_cohort) && missing(keep_pathology)){
      these_sample_ids = dplyr::select(this_meta, sample_id)
    }
  }

  #get qc data for selected samples
  qc_metrics = collate_results(sample_table = these_sample_ids, seq_type_filter = seq_type)

  message(paste0("QC Metric successfully retreived for ", nrow(qc_metrics),
                 " samples out of a total of ", nrow(these_sample_ids), " samples in input sample table."))

  #subset alignment metrics
  melt_align = dplyr::select(qc_metrics, c(sample_id, TotalReads, TotalUniquelyMapped, TotalDuplicatedreads)) %>%
    as.data.frame() %>%
    tidyr::pivot_longer(
        !sample_id,
        names_to = "variable",
        values_to = "value"
    ) %>%
    arrange(sample_id)

  mean_cov_df = data.frame(Metric = c("TotalReads", "TotalUniquelyMapped", "TotalDuplicatedreads"),
                           Value = c (mean(melt_align$value[melt_align$variable == "TotalReads"]),
                                      mean(melt_align$value[melt_align$variable == "TotalUniquelyMapped"]),
                                      mean(melt_align$value[melt_align$variable == "TotalDuplicatedreads"])))

  if(!missing(comparison_group)){
    comp_data = collate_results(sample_table = comparison_group, seq_type_filter = seq_type) %>%
      dplyr::select(sample_id, TotalReads, TotalUniquelyMapped, TotalDuplicatedreads) %>%
      as.data.frame() %>%
      tidyr::pivot_longer(
        !sample_id,
        names_to = "variable",
        values_to = "value"
      ) %>%
      arrange(sample_id)

    mean_cov_df_comp = data.frame(Metric = c("TotalReads", "TotalUniquelyMapped", "TotalDuplicatedreads"),
                                  Value = c (mean(comp_data$value[comp_data$variable == "TotalReads"]),
                                             mean(comp_data$value[comp_data$variable == "TotalUniquelyMapped"]),
                                             mean(comp_data$value[comp_data$variable == "TotalDuplicatedreads"])))
  }

  #corrected mean coverage
  melt_cov = dplyr::select(qc_metrics, c(sample_id, MeanCorrectedCoverage)) %>%
    as.data.frame() %>%
    tidyr::pivot_longer(
        !sample_id,
        names_to = "variable",
        values_to = "value"
    ) %>%
    arrange(sample_id)

  #plot alignment data
  p = ggplot() +
    geom_bar(melt_align, mapping = aes(x = sample_id, y = value, fill = variable), position = "dodge", stat = "identity") +
    {if(add_mean)geom_hline(mean_cov_df, mapping = aes(yintercept = Value, color = Metric))} +
    {if(!missing(comparison_group)) geom_hline(mean_cov_df_comp, mapping = aes(yintercept = Value, linetype = Metric), color = "#54555E")} +
    {if(add_corrected_coverage)geom_point(melt_cov, mapping = aes(x = sample_id, y = value * 25000000, shape = variable), fill = "#A892B3", color = "#5C4966", size = 3, position = "dodge", stat = "identity")} +
    {if(add_corrected_coverage)scale_y_continuous(expand = c(0, 0), breaks = seq(0, max(melt_align$value), by = 3e+08), sec.axis = sec_axis(~ . / 2500000000, name = "", labels = function(b){paste0(round(b * 100, 0), "X")}))} +
    {if(!add_corrected_coverage)scale_y_continuous(expand = c(0, 0), breaks = seq(0, max(melt_align$value), by = 3e+08))} +
    labs(title = "Alignment Summary", subtitle = plot_sub, x = "", y = "Reads (n)") +
    scale_fill_manual(values = this_color_palette) +
    scale_linetype_manual(values = c("TotalReads" = "solid", "TotalUniquelyMapped" = "dashed", "TotalDuplicatedreads" = "dotted")) +
    scale_shape_manual(values = c("MeanCorrectedCoverage" = 21)) +
    scale_color_manual(values = c(this_color_palette)) +
    theme_Morons() +
    labs(linetype = "Comparison Group", shape = "Corrected Coverage (right y-axis)", fill = "Alignment Metrics", color = "Alignment Metrics (Mean)")

  return(p)
}
#' @title Copy Number states barplot
#'
#' @description Plot sample-specific CN states and affected bases for each segment
#'
#' @details `fancy_cnbar` visualizes copy number (CN) states on sample-level. Similarly to other fancy_x_plots this function
#' accepts either a sample ID, for which the function will get copy number states with [GAMBLR::get_sample_cn_segments]. The function
#' can also accept an already loaded seg file given to the `seg_data` parameter. It can also load a seg file with the `seg_path`
#' parameter. If the user calls either `seg_data` or `seg_path`, there are a collection of parameters available for specifying
#' the relevant columns in the given data frame (`chrom_col`, `starat_col`, `end_col`, `cn_col`). It is also possible to
#' restrict the returned plot to any given chromosome. This is done with the `chr_select` parameter (default is all autosomes).
#' For further control of the returned plot, it is also possible to set the threshold for maximum CN states to be returned (default is 15).
#' With `include_cn2` (Boolean) the user can control if CN segments = 2 should be added to the plot, default is TRUE.
#' The user can also control the annotations of the returned plot with `plot_title` and `plot_subtitle`. Lastly,
#' This function also computes the number of affected bases for each copy number state and plots these values on a secondary y-axis (right),
#' useful for overviewing the extent of each copy number state, in the context of the full genome.
#'
#' @param this_sample_id Sample to be plotted.
#' @param seg_data Optional parameter with copy number df already loaded into R.
#' @param seg_path Optional parameter with path to external cn file.
#' @param chrom_col Index of column with chromosome annotations (to be used with either maf_data or maf_path).
#' @param start_col Index of column with copy number start coordinates (to be used with either maf_data or maf_path).
#' @param end_col Index of column with copy number end coordinates (to be used with either maf_data or maf_path).
#' @param cn_col Index of column holding copy number information (to be used with either maf_data or maf_path).
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Subtitle for created plot.
#' @param chr_select Vector of chromosomes to be included in plot, defaults to autosomes.
#' @param cutoff Set threshold for maximum CN state to be retrieved.
#' @param include_cn2 Optional boolean statement for including CN = 2 states in plot.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import ggplot2 dplyr GAMBLR.helpers stringr
#' @export
#'
#' @examples
#' \dontrun{
#' suppressPackageStartupMessages(library(GAMBLR.data))
#'
#' #Return a plot for one sample, with default parameters.
#' fancy_cnbar(this_sample_id = "DOHH-2")
#' }
fancy_cnbar = function(this_sample_id,
                       seg_data,
                       seg_path = NULL,
                       chrom_col = 2,
                       start_col = 3,
                       end_col = 4,
                       cn_col = 7,
                       plot_title = paste0(this_sample_id),
                       plot_subtitle = "n CNV Segments (barplots, left y-axis), n Affected bases for each CN state",
                       chr_select = paste0("chr", c(1:22)),
                       cutoff = 15,
                       include_cn2 = TRUE,
                       this_seq_type = "genome") {
  .Deprecated(msg = "fancy_cnbar() is deprecated and will be removed in a future version of GAMBLR.viz.")

  if(!missing(seg_data)){
    seg = seg_data
    seg = as.data.frame(seg)
    colnames(seg)[chrom_col] = "chrom"
    colnames(seg)[start_col] = "start"
    colnames(seg)[end_col] = "end"
    colnames(seg)[cn_col] = "CN"

  }else if (!is.null(seg_path)){
    seg = read.table(seg_path, sep = "\t", header = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
    seg = as.data.frame(seg)
    colnames(seg)[chrom_col] = "chrom"
    colnames(seg)[start_col] = "start"
    colnames(seg)[end_col] = "end"
    colnames(seg)[cn_col] = "CN"
  }

  #get seg data for a specific sample.
  if(missing(seg_data) && is.null(seg_path)){
    seg = get_sample_cn_segments(these_sample_ids = this_sample_id,
                                 streamlined = FALSE,
                                 this_seq_type = this_seq_type
    )
  }

  #add chr prefix if missing
  if(!str_detect(seg$chrom, "chr")[2]){
    seg = mutate(seg, chrom = paste0("chr", chrom))
  }

  #read maf into R and select relevant variables and transformt to factor.
  seg_df = dplyr::select(seg, chrom, start, end, CN) %>%
    mutate_at(vars(chrom), list(factor))

  #subsetting maf based on user-defined parameters
  seg_df = seg_df[seg_df$chrom %in% chr_select, ]

  #transform data type
  seg_df$CN = as.factor(seg_df$CN)

  #count levels of factor
  if(include_cn2){
    cn_states = c(0:cutoff)
    cns_count = dplyr::filter(seg_df, CN %in% cn_states) %>%
      group_by(CN) %>%
      summarize(count = n())
  }else{
    cn_states = c(0:1, 3:cutoff)
    cns_count = dplyr::filter(seg_df, CN %in% cn_states) %>%
      group_by(CN) %>%
      summarize(count = n())
  }

  cns_count$Type = paste0(cns_count$CN)
  cns_count = dplyr::select(cns_count, count, Type)

  #compute lenght of cn segments and transform for plotting
  l_cn_seg = seg
  l_cn_seg$lenght = l_cn_seg$end - l_cn_seg$start
  l_cn_seg$CN = as.factor(l_cn_seg$CN)
  l_cn_seg = dplyr::filter(l_cn_seg, CN %in% cn_states)

  if(!include_cn2){
    l_cn_seg = dplyr::filter(l_cn_seg, CN != 2)
  }

  cn_seg_lenghts = aggregate(l_cn_seg$lenght, list(l_cn_seg$CN), sum)
  colnames(cn_seg_lenghts) = c("CN", "lenght")

  joined_cn = cbind(cns_count, cn_seg_lenghts) %>%
    as.data.frame() %>%
    dplyr::select(CN, count, lenght)

  #get levels of cn states for plotting
  cn_levels = cns_count$Type

  #plot
  p = ggplot(joined_cn, aes(x = CN)) +
    geom_segment(aes(y = 1, yend = lenght/500000, x = CN, xend = CN)) +
    geom_point(aes(y = lenght/500000), colour = "#E6B315", size = 3, group = 2) +
    geom_bar(aes(y = count, fill = CN), position = "stack", stat = "identity") +
    scale_y_log10(limits = c(1, max(joined_cn$count) + 5000), sec.axis = sec_axis(~.*500000, name = "Nucleotides (n)")) +
    labs(title = plot_title, subtitle = plot_subtitle, x = "CN States", y = "CN Segments (n)", fill = "Legend") +
    scale_fill_manual(values = GAMBLR.helpers::get_gambl_colours("copy_number")) +
    scale_x_discrete(limits = cn_levels) +
    geom_text(aes(x = CN, y = count, label = count), colour = "#000000", size = 5, position = position_stack(vjust = 0.5)) +
    theme_Morons() +
    theme(legend.position = "none")

  return(p)
}
#' @title genome-wide ideogram annotated with SSM and CN information
#'
#' @description Generate sample-level ideogram with copy number information, ssm and gene annotations, etc.
#'
#' @details This function generates genome-wide ideograms, visualizing SSM data as well as CN segments.
#' It is also possible to superimpose the plot with gene annotations. Offering a comprehensive overview of all SSM and CN segments of different aneuploidy.
#' The plotting of SSM can be toggled with setting `include_ssm` to TRUE. If so, it is also possible to count the number of SSMs per chromosome with `ssm_count = TRUE`.
#' To get data for plotting, there are a few different options available; like all `fanncy_x_plots` a sample ID can be provided to the `this_sample_id`
#' parameter. If done so, the function will retrieve data (SSm and CN segments) by wrapping the appropriate functions.
#' This data can also be provided with `seg_data`, `seg_path`, `maf_data` and `maf_path`.
#' For more info on how to run with these parameters, refer to the parameter descriptions.
#' In order to annotate the ideogram with genes, simply give the `gene_annotations` parameter a set of genes as a vector of characters or a data frame with gene names in the first column.
#' Another useful parameter for restricting the plotted regions is to call the function with `intersect_regions`.
#' This parameter takes a vector of characters or a data frame with regions that the plotted calls are restricted to.
#'
#' @param this_sample_id Sample to be plotted (for multiple samples, see fancy_multisample_ideogram.
#' @param gene_annotation Annotate ideogram with a set of genes. These genes can either be specified as a vector of characters or a data frame.
#' @param seg_data Optional parameter with copy number df already loaded into R.
#' @param seg_path Optional parameter with path to external seg like file.
#' @param maf_data Optional parameter with maf like df already loaded into R.
#' @param maf_path Optional parameter with path to external maf like file.
#' @param variant_type_col_maf Index of column holding Variant Type (to be used with either maf_data or maf_path).
#' @param chromosome_col_maf Index of column holding Chromosome (to be used with either maf_data or maf_path).
#' @param start_col_maf Index of column with variant start coordinates (to be used with either maf_data or maf_path).
#' @param end_col_maf Index of column with variant end coordinates (to be used with either maf_data or maf_path).
#' @param chrom_col_seg Index of column with chromosome annotations (to be used with either maf_data or maf_path).
#' @param start_col_seg Index of column with copy number start coordinates (to be used with either maf_data or maf_path).
#' @param end_col_seg Index of column with copy number end coordinates (to be used with either maf_data or maf_path).
#' @param cn_col_seg Index of column holding copy number information (to be used with either maf_data or maf_path).
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Optional argument for plot subtitle.
#' @param intersect_regions Optional parameter for subset variant calls to specific regions. Should be either a vector of characters (chr:start-end) or data frame with regions.
#' @param include_ssm Set to TRUE to plot SSMs (dels and ins).
#' @param ssm_count Optional parameter to summarize n variants per chromosome, inlcude_ssm must be set to TRUE.
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import ggplot2 dplyr GAMBLR.helpers GAMBLR.utils
#' @export
#'
#' @examples
#' #build plot
#' \dontrun{
#' fancy_ideogram(this_sample_id = "DOHH-2",
#'                gene_annotation = "MYC",
#'                plot_title = "Sample-level Ideogram Example",
#'                plot_subtitle = "grch37")
#' }
#'
fancy_ideogram = function(this_sample_id,
                          gene_annotation,
                          seg_data,
                          seg_path = NULL,
                          maf_data,
                          maf_path = NULL,
                          variant_type_col_maf = 10,
                          chromosome_col_maf = 5,
                          start_col_maf = 6,
                          end_col_maf = 7,
                          chrom_col_seg = 2,
                          start_col_seg = 3,
                          end_col_seg = 4,
                          cn_col_seg = 7,
                          plot_title = paste0(this_sample_id),
                          plot_subtitle = "Genome-wide Ideogram (grch37).",
                          intersect_regions,
                          include_ssm = TRUE,
                          ssm_count = TRUE,
                          coding_only = FALSE,
                          this_seq_type = "genome"){
  .Deprecated(msg = "fancy_ideogram() is deprecated and will be removed in a future version of GAMBLR.viz.")

  #plot theme
  ideogram_theme = function(){
    theme(legend.position = "bottom", axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank())
  }

  #grch37 coordinates
  grch37_end = dplyr::filter(GAMBLR.data::chromosome_arms_grch37, arm == "q") %>%
    dplyr::filter(chromosome != "Y") %>%
    pull(end)
  grch37_cent_start = dplyr::filter(GAMBLR.data::chromosome_arms_grch37, arm == "p") %>%
    dplyr::filter(chromosome != "Y") %>%
    pull(end)
  grch37_cent_end = dplyr::filter(GAMBLR.data::chromosome_arms_grch37, arm == "q") %>%
    dplyr::filter(chromosome != "Y") %>%
    pull(start)

  #additional regions to plot
  if(!missing(gene_annotation)){
    gene = GAMBLR.utils::gene_to_region(gene_symbol = gene_annotation, projection = "grch37", return_as = "df")
    gene.annotate = gene[gene$chr %in% paste0(c(1:22)), ]
    cols.int = c("chromosome", "start", "end")
    gene.annotate[cols.int] = sapply(gene.annotate[cols.int], as.integer)
  }

  #build chr table for segment plotting
  chr = c( paste0("chr", c(1:22)), "chrX" )
  chr_start = c(0)
  chr_end = grch37_end
  cent_start = grch37_cent_start
  cent_end =  grch37_cent_end
  y = c(1:23)
  yend = c(1:23)

  #transform to data frame
  segment_data = data.frame(chr, chr_start, chr_end, cent_start, cent_end, y, yend)

  #load CN data
  if(!missing(seg_data)){
    cn_states = seg_data
    cn_states = as.data.frame(cn_states)
    colnames(cn_states)[chrom_col_seg] = "chrom"
    colnames(cn_states)[start_col_seg] = "start"
    colnames(cn_states)[end_col_seg] = "end"
    colnames(cn_states)[cn_col_seg] = "CN"

  }else if(!is.null(seg_path)){
    cn_states = read.table(seg_path, sep = "\t", header = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
    cn_states = as.data.frame(cn_states)
    colnames(cn_states)[chrom_col_seg] = "chrom"
    colnames(cn_states)[start_col_seg] = "start"
    colnames(cn_states)[end_col_seg] = "end"
    colnames(cn_states)[cn_col_seg] = "CN"
  }

  #get maf data for a specific sample.
  if(missing(seg_data) && is.null(seg_path)){
    cn_states = get_sample_cn_segments(
      these_sample_ids = this_sample_id,
      with_chr_prefix = FALSE,
      streamlined = FALSE,
      this_seq_type = this_seq_type
    )
  }

  # ignore y chromosome
  cn_states = dplyr::filter(cn_states, chrom != "Y")

  #convert chr into y coordinates
  cn_states$ycoord = cn_states$chrom

  #paste chr in chromosomecolumn, if not there
  if(!str_detect(cn_states$chrom[1], "chr")){
    cn_states = mutate(cn_states, chrom = paste0("chr", chrom))
  }

  if(!missing(intersect_regions)){
    #filter CN states on intersecting regions
    #transform regions to data tables
    #convenience function for converting intersect regions to a df (if it's a string) and renaming columns to match required format.
    if(is.list(intersect_regions)){
      colnames(intersect_regions)[1] = "chrom"
      colnames(intersect_regions)[2] = "start"
      colnames(intersect_regions)[3] = "end"
      if(!str_detect(intersect_regions$chrom, "chr")){
        intersect_regions = mutate(intersect_regions, chrom = paste0("chr", chrom))
      }
      intersect_regions$start = as.numeric(intersect_regions$start)
      intersect_regions$end = as.numeric(intersect_regions$end)
    }

    if(is.character(intersect_regions)){
      if(length(intersect_regions) > 1){
        message("Please only enter one region, only first region will be regarded. For mutiple regions, kindly provide a data frame with regions of interest")
      }

      split_chunks = unlist(strsplit(intersect_regions, ":"))
      split_chunks = unlist(strsplit(split_chunks, "-"))
      chrom = split_chunks[1]
      start = split_chunks[2]
      end = split_chunks[3]
      intersect_regions = cbind(chrom, start, end) %>%
        as.data.frame()

      intersect_regions$start = as.numeric(intersect_regions$start)
      intersect_regions$end = as.numeric(intersect_regions$end)

      if(!str_detect(intersect_regions$chrom[1], "chr")){
        intersect_regions = mutate(intersect_regions, chrom = paste0("chr", chrom))
      }
    }

    #intersect regions
    intersect = cool_overlaps(
        regions_sub, incoming_cn,
        columns1 = c("chrom", "start", "end"),
        columns2 = c("chrom", "start", "end")
    )

    #transform object to data frame
    inter_df = as.data.frame(intersect)

    #organize columns to match the expected format
    cn_states = select(inter_df, ID, chrom, start, end, LOH_flag, log.ratio, CN, ycoord)
  }

  #convert data types
  cn_states$chrom = as.factor(cn_states$chrom)

  # correct y coordinate for X chromosome
  cn_states$ycoord = dplyr::recode(cn_states$ycoord, X="23") %>%
    as.integer

  #subset on CN state
  cn_states$CN[cn_states$CN > 6] = 6
  cn_states = subset(cn_states, CN != 2)
  cn_states$CN = paste0("cn_", cn_states$CN)
  cn_states$CN = as.factor(cn_states$CN)
  cn_states = droplevels(cn_states)
  l = split(cn_states, cn_states$CN)
  list2env(l, envir = .GlobalEnv)

  #load maf data
  if(include_ssm){
    if(!missing(maf_data)){
      maf = maf_data
      maf = as.data.frame(maf)
      colnames(maf)[variant_type_col_maf] = "Variant_Type"
      colnames(maf)[chromosome_col_maf] = "Chromosome"
      colnames(maf)[start_col_maf] = "Start_Position"
      colnames(maf)[end_col_maf] = "End_Position"
    }else if (!is.null(maf_path)){
      maf = GAMBLR.utils::fread_maf(maf_path)
      maf = as.data.frame(maf)
      colnames(maf)[variant_type_col_maf] = "Variant_Type"
      colnames(maf)[chromosome_col_maf] = "Chromosome"
      colnames(maf)[start_col_maf] = "Start_Position"
      colnames(maf)[end_col_maf] = "End_Position"
    }

    if(missing(maf_data) && is.null(maf_path)){
      maf = get_ssm_by_sample(this_sample_id = this_sample_id,
                              projection = "grch37", #currently only supported reference build.
                              this_seq_type = this_seq_type)
    }

    #transform maf data
    maf_trans = dplyr::select(maf, Chromosome, Start_Position, End_Position, Variant_Type)

    #calculate mid points of variants
    maf_trans$mid = ((maf_trans$End_Position - maf_trans$Start_Position) / 2) + maf_trans$Start_Position

    #convert chr into y coordinates
    maf_trans$ystart = maf_trans$yend = dplyr::recode(maf_trans$Chromosome, X="23")

    #paste chr in maf, if not there
    if(!str_detect(maf_trans$Chromosome[1], "chr")){
      maf_trans = mutate(maf_trans, Chromosome = paste0("chr", Chromosome))
    }

    #convert data types
    maf_trans$Start_Position = as.double(maf_trans$Start_Position)
    maf_trans$End_Position = as.double(maf_trans$End_Position)
    maf_trans$ystart = as.integer(maf_trans$ystart)
    maf_trans$yend = as.integer(maf_trans$yend)

    if(!missing(intersect_regions)){
      #filter CN states on intersecting regions
      maf_tmp = maf_trans
      colnames(maf_tmp)[1] = "chrom"
      colnames(maf_tmp)[2] = "start"
      colnames(maf_tmp)[3] = "end"
      maf_tmp = dplyr::select(maf_tmp, chrom, start, end)
      maf.table = as.data.frame(maf_tmp)

      #intersect regions
      intersect_maf = cool_overlaps(
        regions_sub, maf.table,
        columns1 = c("chrom", "start", "end"),
        columns2 = c("chrom", "start", "end")
        )

      #transform object to data frame
      inter_maf_df = as.data.frame(intersect_maf)

      #rename columns
      colnames(inter_maf_df)[1] = "Chromosome"
      colnames(inter_maf_df)[2] = "Start_Position"
      colnames(inter_maf_df)[3] = "End_Position"

      #subset
      inter_maf_df = dplyr::select(inter_maf_df, Chromosome, Start_Position, End_Position)

      #perform a semi join with all cn states (to retain necessary columns)
      maf_trans = dplyr::semi_join(maf_trans, inter_maf_df)
    }

    #subset on variant type
    if(nrow(maf_trans) > 0){
      maf_del = dplyr::filter(maf_trans, Variant_Type == "DEL")
      maf_ins = dplyr::filter(maf_trans, Variant_Type == "INS")

      if(ssm_count){
        del_count = dplyr::filter(maf_trans, Variant_Type == "DEL") %>%
          add_count(Chromosome) %>%
          distinct(Chromosome, .keep_all = TRUE) %>%
          dplyr::select(Chromosome, Variant_Type, yend, n)

        ins_count = dplyr::filter(maf_trans, Variant_Type == "INS") %>%
          add_count(Chromosome) %>%
          distinct(Chromosome, .keep_all = TRUE) %>%
          dplyr::select(Chromosome, Variant_Type, yend, n)
      }
    }
  }

  #get colours and combine palette for indels and cn states
  ideogram_palette = c(GAMBLR.helpers::get_gambl_colours("indels"), GAMBLR.helpers::get_gambl_colours("copy_number"))
  selected_colours = ideogram_palette[c(1,2,19,18,16:13)]
  names(selected_colours)[c(3:8)] = c("CN0", "CN1", "CN3", "CN4", "CN5", "CN6+")

  #plot
  p = ggplot() +
    {if(include_ssm && nrow(maf_del) > 0) geom_segment(data = maf_del, aes(x = mid - 100000, xend = mid + 100000, y = ystart - 0.27, yend = yend - 0.27), color = "#53B1FC", size = 5, stat = "identity", position = position_dodge(width = 0))} + #del
    {if(include_ssm && nrow(maf_del) > 0) geom_segment(data = maf_del, aes(x = mid, xend = mid, y = ystart - 0.35, yend = yend - 0.35), color = "black", lineend = "round", size = 3.5, stat = "identity", position = position_dodge(width = 0))} + #del
    {if(include_ssm && nrow(maf_del) > 0) geom_segment(data = maf_del, aes(x = mid, xend = mid, y = ystart - 0.35, yend = yend - 0.35, color = "DEL"), lineend = "round", size = 3, stat = "identity", position = position_dodge(width = 0))} + #del
    {if(include_ssm && nrow(maf_ins) > 0) geom_segment(data = maf_ins, aes(x = mid - 100000, xend = mid + 100000, y = ystart - 0.27, yend = yend - 0.27), color = "#FC9C6D", size = 5, stat = "identity", position = position_dodge(width = 0))} + #ins
    {if(include_ssm && nrow(maf_ins) > 0) geom_segment(data = maf_ins, aes(x = mid, xend = mid, y = ystart - 0.35, yend = yend - 0.35), color = "black", lineend = "round", size = 3.5, stat = "identity", position = position_dodge(width = 0))} + #ins
    {if(include_ssm && nrow(maf_ins) > 0) geom_segment(data = maf_ins, aes(x = mid, xend = mid, y = ystart - 0.35, yend = yend - 0.35, color = "INS"), lineend = "round", size = 3, stat = "identity", position = position_dodge(width = 0))} + #ins
    {if(ssm_count && nrow(maf_del) > 0) annotate(geom = "text", x = -4000000, y = del_count$yend, label = del_count$n, color = "#3A8799", size = 3)} + #count del
    {if(ssm_count && nrow(maf_trans) > 0) annotate(geom = "text", x = -2300000, y = segment_data$y, label = " | ", color = "black", size = 3)} + #count sep
    {if(ssm_count && nrow(maf_ins) > 0) annotate(geom = "text", x = -1000000, y = ins_count$yend, label = ins_count$n, color = "#E6856F", size = 3)} + #count ins
    geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y, yend = yend), color = "#99A1A6", lineend = "butt", size = 5, stat = "identity", position = position_dodge(width = 0)) + #chr contigs
    geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y, yend = yend), color = "white", size = 6, stat = "identity", position = position_dodge(width = 0)) + #centromeres
    {if("cn_0" %in% levels(cn_states$CN)) geom_segment(data = cn_0, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN0"), size = 4.7, stat = "identity", position = position_dodge(width = 0))} + #cn3
    {if("cn_1" %in% levels(cn_states$CN)) geom_segment(data = cn_1, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN1"), size = 4.7, stat = "identity", position = position_dodge(width = 0))} + #cn3
    {if("cn_3" %in% levels(cn_states$CN)) geom_segment(data = cn_3, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN3"), size = 4.7, stat = "identity", position = position_dodge(width = 0))} + #cn3
    {if("cn_4" %in% levels(cn_states$CN)) geom_segment(data = cn_4, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN4"), size = 4.7, stat = "identity", position = position_dodge(width = 0))} + #cn4
    {if("cn_5" %in% levels(cn_states$CN)) geom_segment(data = cn_5, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN5"), size = 4.7, stat = "identity", position = position_dodge(width = 0))} + #cn5
    {if("cn_6" %in% levels(cn_states$CN)) geom_segment(data = cn_6, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN6+"), size = 4.7, stat = "identity", position = position_dodge(width = 0))} + #cn6 and more
    {if(!missing(gene_annotation)) geom_point(data = gene.annotate, aes(x = ((end - start) / 2) + start, y = chromosome - 0.28), shape = 25, color = "#A63932", fill = "#A63932", stat = "identity", position = position_dodge(width = 0))} + #gene annotation
    {if(!missing(gene_annotation)) geom_label(data = gene.annotate, aes((x = end - start) / 2 + start, y = chromosome - 0.52, label = hugo_symbol), fontface = "bold", color = "white", fill = "#A63932", size = 3)} + #gene annotation text
    geom_text(aes(x = -10000000 , y = yend, label = segment_data$chr), color = "black", size = 5) + #chr labels
    labs(title = plot_title, subtitle = plot_subtitle) + #plot titles
    scale_colour_manual(name = "", values = selected_colours) + #legend/colours
    scale_x_continuous(breaks = seq(0, max(segment_data$chr_end), by = 30000000)) + #set x-axis boundaries
    scale_y_reverse() + #reverse y axis
    theme_Morons()

  return(p)
}
#' @title Plot Quality Control Metrics.
#'
#' @description Plot for visualizing QC metrics and allowing for grouping by different metadata columns.
#'
#' @details This function is readily available for visualizing a variety of quality control metrics. To get started, the user can easily overview all the available metrics with `return_plotdata = TRUE`.
#' When this parameter is set to TRUE, a vector of characters will be returned detailing all the, for this plot, available metrics. After deciding what metric to plot, simply give the metric of choice to the `plot_data` parameter.
#' This function also lets the user provide a data frame with sample IDs to be included in the plot. Optionally, the user can also provide an already filtered metadata table with sample IDs of interest to the `these_samples_metadata`.
#' If none of the two parameters are supplied, the user can easily restrict the plot to any cohort and/or pathology of their liking. This is done by calling `keep_cohort` and `keep_pathology`.
#' If these parameters are used, the function will retrieve metadata for all available GAMBL sample IDs and then subset to the specified cohort or pathology.
#' The layout of the returned plot can also be further customized with `sort_by`. This parameter controls the order in which samples would appear. Similarly, `fill_by` allows the user to control on what factor the plot will be filled by.
#' Sometimes it can also be useful to see how a subset of samples compares to another group; to do this one could call the function with a vector of additional sample IDs given to the `comparison_samples` parameter (see examples for more information).
#' lastly, the plot can also be configured with custom plot title and axis labels (`plot_title` and `y_axis_lab`). For more information, see examples and parameter descriptions.
#'
#' @param these_sample_ids Data frame with sample IDs (to be plotted) in the first column (has to be named sample_id).
#' @param keep_cohort Optional parameter to be used when these_sample is NULL. Returns metadata and filters on the cohort supplied in this parameter.
#' @param keep_pathology Optional parameter to be used when these_sample is NULL. Returns metadata and filters on the pathology supplied in this parameter.
#' @param seq_type Selected seq type for incoming QC metrics.
#' @param metadata Optional, user can provide a metadata df to subset sample IDs from.
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process.
#' @param plot_data Plotting parameter, define the data type to be plotted.
#' @param fill_by Parameter for specifying fill variable for grouped bar plot. Can be any factor from incoming metadata, e.g pathology, cohort, etc.
#' @param comparison_samples Optional parameter, give the function a vector of sample IDs to be compared against the main plotting group. Pathology is default.
#' @param plot_title Plotting parameter, plot title.
#' @param y_axis_lab Plotting parameter, label of y-axis.
#' @param return_plotdata Optional parameter, if set to TRUE a vector of acceptable data types for plotting will be returned, and nothing else.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import dplyr ggplot2 GAMBLR.helpers ggbeeswarm
#' @export
#'
#' @examples
#' \dontrun{
#' #load packages
#' suppressPackageStartupMessages(library(dplyr))
#' suppressPackageStartupMessages(library(GAMBLR.open))
#'
#' #get sample IDs for available genome samples
#' genome_collated = collate_results(seq_type_filter = "genome") %>%
#'   pull(sample_id)
#'
#' #subset the collated samples on BL samples
#' my_samples = get_gambl_metadata() %>%
#'   dplyr::filter(sample_id %in% genome_collated) %>%
#'   dplyr::filter(pathology == "BL") %>% pull(sample_id)
#'
#' fancy_qc_plot(these_sample_ids = my_samples, plot_data = "AverageBaseQuality")
#'}
fancy_qc_plot = function(these_sample_ids,
                         keep_cohort,
                         keep_pathology,
                         seq_type = "genome",
                         metadata,
                         these_samples_metadata,
                         plot_data,
                         fill_by = "pathology",
                         comparison_samples,
                         plot_title = "",
                         y_axis_lab = "",
                         return_plotdata = FALSE){
  .Deprecated(msg = "fancy_qc_plot() is deprecated and will be removed in a future version of GAMBLR.viz.")

  #return a list of acceptable data types for plotting
  if(return_plotdata){
    plotting_variables = c("AverageBaseQuality", "AverageInsertSize", "AverageReadLength",
                           "PairsOnDiffCHR", "TotalReads", "TotalUniquelyMapped",
                           "TotalUnmappedreads", "TotalDuplicatedreads", "ProportionReadsDuplicated",
                           "ProportionReadsMapped", "MeanCorrectedCoverage", "ProportionCoverage10x", "ProportionCoverage30x")

    return(plotting_variables)
  }

  #get gambl metadata (if not supplied)
  if(missing(metadata)){
    this_meta = GAMBLR.helpers::handle_metadata(this_seq_type = seq_type)
  }else{
    this_meta = metadata
  }

  if(!missing(these_samples_metadata)){
    these_sample_ids = dplyr::select(these_samples_metadata, sample_id) %>%
      as.data.frame(strings.as.factors = FALSE)
  }

  #filter metadata on selected cohort/pathology
  if(missing(these_sample_ids) && missing(these_samples_metadata)){
    if(!missing(keep_cohort) && missing(keep_pathology)){
      these_sample_ids = dplyr::filter(this_meta, cohort == keep_cohort)
    }

    if(!missing(keep_pathology) && missing(keep_cohort)){
      these_sample_ids = dplyr::filter(this_meta, pathology == keep_pathology)
    }

    if(!missing(keep_cohort) && !missing(keep_pathology)){
      these_sample_ids = dplyr::filter(this_meta, pathology == keep_pathology, cohort == keep_cohort)
    }

    if(missing(keep_cohort) && missing(keep_pathology)){
      these_sample_ids = dplyr::select(this_meta, sample_id)
    }
  }

  #get QC data for selected samples
  qc_metrics = collate_results(sample_table = these_sample_ids, seq_type_filter = seq_type)
  message(paste0("QC Metric successfully retreived for ", nrow(qc_metrics), " samples out of a total of ", nrow(these_sample_ids), " samples in input sample table."))

  #aggregate sample list with metadata columns
  qc_meta = qc_metrics %>%
    inner_join(this_meta) %>%
    mutate_if(is.integer, as.factor) %>%
    mutate_if(is.character, as.factor)

  qc_meta$group = "main_sample"

  #Retrieve QC metrics for comparison samples, if provided.
  if(!missing(comparison_samples)){
    comp_data = collate_results(sample_table = comparison_samples, seq_type_filter = seq_type)

    #aggregate sample list with metadata columns
    comp_meta = comp_data %>% inner_join(this_meta)
    comp_meta = mutate_if(comp_meta, is.integer, as.factor)
    comp_meta = mutate_if(comp_meta, is.character, as.factor)
    comp_meta$group = "comparison_sample"
    qc_meta = rbind(qc_meta, comp_meta)
  }

  #get gambl colours for selected fill and subset to levels in selected factor
  col_gambl = GAMBLR.helpers::get_gambl_colours(fill_by) %>%
    as.data.frame()

  col_gambl$factors = rownames(col_gambl)
  colnames(col_gambl)[1] = "hex"
  row.names(col_gambl) <- NULL

  levels_fill = levels(qc_meta[[fill_by]]) %>%
    as.data.frame()

  colnames(levels_fill)[1] = "factors"
  sub_cols = dplyr::left_join(levels_fill, col_gambl, by = "factors")
  list_col = as.list(sub_cols$hex)

  #plotting
  p = ggplot(qc_meta) +
    aes_string(x = paste0("group"), y = plot_data, fill = fill_by) +
    geom_boxplot(mapping = aes(x = group), outlier.shape = NA) +
    geom_quasirandom() +
    labs(title = plot_title, x = "", y = y_axis_lab) +
    theme_Morons() +
    scale_fill_manual(values = c(list_col))

  return(p)
}
#' @title n SNVs per chromosome plot
#'
#' @description Visualizing the number of SNVs per chromosome.
#'
#' @details This function takes on an already loaded maf-like data frame, or a path to the maf file of interest.
#' In addition, the user can also give this function a sample ID and the function will run get_ssm_by_sample
#' to get data for plotting. If a maf file or data frame is used, the user has the chance to specify what column
#' that holds the Variant Type information (`variant_type_col`), in addition the user can also specify what column
#' in the incoming maf that is corresponding to the chromosome annotations. This function also includes useful subsetting
#' options. For example, `chr_select` allows the user to restrict the plot to specific chromosomes. `include_dnp` is an optional
#' argument (Boolean) for if variants of this subtype should be included or not. Lastly, this plotting function
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
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param projection Genome build for returned variants. Default is grch37.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import ggplot2 dplyr GAMBLR.helpers GAMBLR.utils
#' @export
#'
#' @examples
#' \dontrun{
#' suppressPackageStartupMessages(library(GAMBLR.open))
#' 
#' #plot SNVs
#' fancy_snv_chrdistplot(this_sample_id = "DOHH-2")
#'
#' #plot SNVs and DNPs
#' fancy_snv_chrdistplot(this_sample_id = "DOHH-2",
#'                       include_dnp = TRUE,
#'                       plot_subtitle = "SNV + DNP Distribution Per Chromosome")
#'}
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
                                 this_seq_type = "genome",
                                 projection = "grch37"){
  .Deprecated(msg = "fancy_snv_chrdistplot() is deprecated and will be removed in a future version of GAMBLR.viz.")

  if(!missing(maf_data)){
    maf = maf_data
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col] = "Variant_Type"
    colnames(maf)[chromosome_col] = "Chromosome"

  }else if (!is.null(maf_path)){
    maf = GAMBLR.utils::fread_maf(maf_path)
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col] = "Variant_Type"
    colnames(maf)[chromosome_col] = "Chromosome"
  }

  #get maf data for a specific sample.
  if(missing(maf_data) && is.null(maf_path)){
    maf = get_ssm_by_sample(this_sample_id = this_sample_id,
                            this_seq_type = this_seq_type,
                            projection = projection)
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
      theme_Morons() +
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
      theme_Morons() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      {if(hide_legend)theme(legend.position = "none")} +
      coord_flip()
  }
}
#' @title Structural Variants Size Plot.
#'
#' @description Generate plot visualizing SV sizes. Subset on variant type, filter on VAF, size etc.
#'
#' @details Plot sample-level SV sizes across selected chromosomes. This function also has a variety of filtering parameters available.
#' For example, it is possible to subset the included variants to a specific VAF threshold with `VAF_cutoff`. The `size_cutoff` is another parameter
#' for filtering the variants on set variant sizes, the default for this parameter is to only include variants of at least 50bp.
#' This function takes either a sample ID (`this_sample_id`) or an already loaded data frame (`maf_data` or a path to a maf-like file with `maf_path`).
#' If `this_sample_id` is called, the function will run [GAMBLR::get_combined_sv] to retrieve SV calls.
#' If either of the `maf` parameters are used, note that it's possible to specify the columns of interest;
#' (`chrom_a_col`, `start_a_col`, `end_a_col` and `variant_type_col`), allowing this function to work with any maf-like data frames.
#' This function also allows the user to customize the returned plot. For more info on how to do this, please refer to the aesthetic
#' parameters; `hide_legend`, `plot_title`, `plot_subtitle`, `adjust_value` and `trim`.
#'
#' @param this_sample_id Sample to be plotted.
#' @param maf_data Optional parameter with copy number df already loaded into R.
#' @param maf_path Optional parameter with path to external cn file.
#' @param chrom_a_col Index of column holding chromosome (to be used with either maf_data or maf_path).
#' @param start_a_col Index of column holding start coordinates (to be used with either maf_data or maf_path).
#' @param end_a_col Index of column holding end coordinates (to be used with either maf_data or maf_path).
#' @param variant_type_col Index of column holding variant type information (to be used with either maf_data or maf_path).
#' @param vaf_cutoff Threshold for filtering variants on VAF (events with a VAF > cutoff will be retained).
#' @param size_cutoff Threshold for filtering variants on size, default is 50bp.
#' @param adjust_value A multiplicate bandwidth adjustment. This makes it possible to adjust the bandwidth while still using the bandwidth estimator. For example, adjust = 1/2 means use half of the default bandwidth.
#' @param trim If FALSE, the default, each density is computed on the full range of the data.
#' @param chr_select Optional argument for subsetting on selected chromosomes, default is all autosomes.
#' @param hide_legend Set to True to remove legend from plot, default is FALSE.
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Subtitle for created plot.
#' @param projection Genomic projection for SVs and circos plot. Accepted values are grch37 and hg38.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import dplyr ggplot2 GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#' suppressPackageStartupMessages(library(GAMBLR.data))
#'
#' #build plot sith default parameters
#' fancy_sv_sizedens(this_sample_id = "SP116715")
#'
#' #restrict plot to only chromosome 1 and 2
#' fancy_sv_sizedens(this_sample_id = "SP116715",
#'                   size_cutoff = 0,
#'                   chr_select = c("chr1", "chr2"))
#' }
#'
fancy_sv_sizedens = function(this_sample_id,
                             maf_data,
                             maf_path = NULL,
                             chrom_a_col = 3,
                             start_a_col = 4,
                             end_a_col = 5,
                             variant_type_col = 9,
                             vaf_cutoff = 0,
                             size_cutoff = 50,
                             adjust_value = 1,
                             trim = FALSE,
                             hide_legend = FALSE,
                             chr_select = paste0("chr", c(1:22)),
                             plot_title = paste0(this_sample_id),
                             plot_subtitle = paste0("SV sizes for Manta calls.
                                                    Dashed line annotates mean
                                                    variant size.\nVAF cut off:
                                                    ", vaf_cutoff,", SV size
                                                    cut off: ", size_cutoff),
                             projection = "grch37"){
  .Deprecated(msg = "fancy_sv_sizedens() is deprecated and will be removed in a future version of GAMBLR.viz.")
  if(!missing(maf_data)){
    svs = maf_data
    svs = as.data.frame(svs)
    colnames(svs)[chrom_a_col] = "CHROM_A"
    colnames(svs)[start_a_col] = "START_A"
    colnames(svs)[end_a_col] = "END_A"
    colnames(svs)[variant_type_col] = "manta_name"

  }else if (!is.null(maf_path)){
    svs = maf_data
    svs = as.data.frame(svs)
    colnames(svs)[chrom_a_col] = "CHROM_A"
    colnames(svs)[start_a_col] = "START_A"
    colnames(svs)[end_a_col] = "END_A"
    colnames(svs)[variant_type_col] = "manta_name"
  }

  #get variants, filter and subset
  if(missing(maf_data) && is.null(maf_path)){
    svs = get_manta_sv(these_sample_ids = this_sample_id, projection = projection, min_vaf = vaf_cutoff) %>%
      dplyr::filter(VAF_tumour > vaf_cutoff) %>%
      dplyr::select(CHROM_A, START_A, END_A, manta_name)
  }

  #split manta_name variable
  svs_df = data.frame( svs$CHROM_A,   svs$START_A, svs$END_A,
                       sub("^(.+?):.*", "\\1", svs$manta_name) )

  #rename variables
  names(svs_df)[1] = "chrom"
  names(svs_df)[2] = "start"
  names(svs_df)[3] = "end"
  names(svs_df)[4] = "type"

  #subset df on SV type
  manta_sv = dplyr::filter(svs_df, type %in% c("MantaDEL", "MantaDUP")) %>%
    dplyr::select(chrom, start, end, type)

  # check whether enough variants
  check_whether_enough_vars = function(manta_type, string1, string2){
    type_table = table(manta_type)
    var_num_message = "%i MantaDEL and %i MantaDUP variants were"
    if( is.null(string2) ){
      var_num_message = gettextf("%s %s.", var_num_message, string1)
    }else{
      var_num_message = gettextf("%s %s after filtering by %s.",
                                 var_num_message, string1, string2)
    }
    type_table = table(manta_type)
    k <- gettextf(var_num_message, type_table["MantaDEL"], type_table["MantaDUP"])
    message(k)
    stopifnot("Plot couldn't be made. At least 2 variants of either type are needed." =
                any(type_table > 1))
  }

  manta_sv = mutate( manta_sv, type = factor(type, levels = c("MantaDEL", "MantaDUP")) )

  if( missing(maf_data) && is.null(maf_path) ){
    check_whether_enough_vars(manta_type = manta_sv$type, string1 = "found", string2 = "vaf_cutoff")
  }else{
    check_whether_enough_vars(manta_type = manta_sv$type, string1 = "found", string2 = NULL)
  }

  #calculate sizes
  manta_sv$size = manta_sv$end - manta_sv$start

  #add chr prefix, if missing
  if(!grepl("chr", manta_sv$chrom[1])){
    manta_sv = mutate(manta_sv, chrom = paste0("chr", chrom))
  }

  #subset on selected chromosomes
  manta_sv = manta_sv[manta_sv$chrom %in% chr_select, ]

  #filter out variants < 50 bp
  manta_sv = dplyr::filter(manta_sv, size >= size_cutoff)

  # check whether enough variants
  check_whether_enough_vars(manta_type = manta_sv$type, string1 = "left", string2 = "size_cutoff")

  # groups (MantaDEL or MantaDUP) with only 1 variant are dropped.
  if(type_table["MantaDEL"] == 1){
    manta_sv = filter(manta_sv, type != "MantaDEL")
    message("Warning: At least 2 data points are needed to calculate density estimates. MantaDEL group was dropped because it contains only 1 variant.")
  }
  if(type_table["MantaDUP"] == 1){
    manta_sv = filter(manta_sv, type != "MantaDUP")
    message("Warning: At least 2 data points are needed to calculate density estimates. MantaDUP group was dropped because it contains only 1 variant.")
  }
  manta_sv = mutate(manta_sv, type = droplevels(type))

  del_col = GAMBLR.helpers::get_gambl_colours("indels")[[1]]
  dup_col = GAMBLR.helpers::get_gambl_colours("indels")[[2]]

  #plotting
  p = ggplot(manta_sv, aes(x = size, fill = type)) +
        geom_density(alpha = 0.7, color = NA, adjust = adjust_value, trim = trim) +
        labs(title = plot_title, subtitle = plot_subtitle, x = "Size (bp)", y = "") +
        scale_fill_manual(values = c(del_col, dup_col)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        {if(hide_legend)theme(legend.position = "none")} +
        theme_Morons()

  return(p)
}
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
#' @import ggplot2 dplyr GAMBLR.helpers GAMBLR.utils
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' suppressPackageStartupMessages(library(GAMBLR.open))
#'
#' #plot ssm
#' fancy_v_chrcount(this_sample_id = "DOHH-2",
#'                  y_interval = 10)
#'}
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
  .Deprecated(msg = "fancy_v_chrcount() is deprecated and will be removed in a future version of GAMBLR.viz.")

  if(!missing(maf_data)){
    plot_data = maf_data
    plot_data = as.data.frame(plot_data)
    colnames(plot_data)[variant_type_col] = "Variant_Type"
    colnames(plot_data)[chromosome_col] = "Chromosome"

  }else if(!is.null(maf_path)){
    plot_data = GAMBLR.utils::fread_maf(maf_path)
    plot_data = as.data.frame(plot_data)
    colnames(plot_data)[variant_type_col] = "Variant_Type"
    colnames(plot_data)[chromosome_col] = "Chromosome"
  }

  #get maf data for a specific sample.
  if(missing(maf_data) && is.null(maf_path)){
    if(ssm){
      plot_data = get_ssm_by_sample(this_sample_id = this_sample_id,
                                    this_seq_type = this_seq_type,
                                    projection = projection)
    }else{
      plot_data = get_manta_sv(these_sample_ids = this_sample_id,
                               projection = projection,
                               min_vaf = min_vaf) %>%
        dplyr::select(CHROM_A, START_A, END_A, manta_name)

      #get manta results in required format
      plot_data = data.frame( plot_data$CHROM_A, plot_data$START_A, plot_data$END_A,
                              sub("^(.+?):.*", "\\1", plot_data$manta_name) )

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
        theme_Morons() +
        {if(hide_legend)theme(legend.position = "none")} +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}
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
#' @import ggplot2 dplyr GAMBLR.utils GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#' suppressPackageStartupMessages(library(GAMBLR.data))
#'
#' #count all variants for one sample (default parameters)
#' fancy_v_count(this_sample_id = "DOHH-2")
#'}
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
  .Deprecated(msg = "fancy_v_count() is deprecated and will be removed in a future version of GAMBLR.viz.")

  if(!missing(maf_data)){
    plot_data = maf_data
    plot_data = as.data.frame(plot_data)
    colnames(plot_data)[variant_type_col] = "Variant_Type"
    colnames(plot_data)[chromosome_col] = "Chromosome"

  }else if (!is.null(maf_path)){
    plot_data = GAMBLR.utils::fread_maf(maf_path)
    plot_data = as.data.frame(plot_data)
    colnames(plot_data)[variant_type_col] = "Variant_Type"
    colnames(plot_data)[chromosome_col] = "Chromosome"
  }

  #get maf data for a specific sample.
  if(missing(maf_data) && is.null(maf_path)){
    if(ssm){
      plot_data = get_ssm_by_sample(this_sample_id = this_sample_id,
                                    this_seq_type = this_seq_type,
                                    projection = projection)
    }else{
      plot_data = get_manta_sv(these_sample_ids = this_sample_id, projection = projection, min_vaf = min_vaf) %>%
        dplyr::select(CHROM_A, START_A, END_A, manta_name)

      #get manta results in required format
      plot_data = data.frame( plot_data$CHROM_A, plot_data$START_A, plot_data$END_A,
                              sub("^(.+?):.*", "\\1", plot_data$manta_name) )

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
    theme_Morons()

  return(p)
}
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
#' @import ggplot2 dplyr GAMBLR.helpers GAMBLR.utils
#' @export
#'
#' @examples
#' \dontrun{
#' suppressPackageStartupMessages(library(GAMBLR.data))
#'
#' #plot SSM size distributions:
#' fancy_v_sizedis(this_sample_id = "DOHH-2")
#'}
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
  .Deprecated(msg = "fancy_v_sizedis() is deprecated and will be removed in a future version of GAMBLR.viz.")

  if(!missing(maf_data)){
    plot_data = maf_data
    plot_data = as.data.frame(plot_data)
    colnames(plot_data)[variant_type_col] = "Variant_Type"
    colnames(plot_data)[chromosome_col] = "Chromosome"
    colnames(plot_data)[start_col] = "Start_Position"
    colnames(plot_data)[end_col] = "End_Position"

  }else if (!is.null(maf_path)){
    plot_data = GAMBLR.utils::fread_maf(maf_path)
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
      plot_data = data.frame( plot_data$CHROM_A, plot_data$START_A, plot_data$END_A,
                              sub("^(.+?):.*", "\\1", plot_data$manta_name) )

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
        theme_Morons() +
        theme(legend.position = "none")

  return(p)
}
#' @title sample-level SV/SSM/CN reports in PDF
#'
#' @description Construct pdf with sample-level plots, using minimum of arguments
#'
#' @details This function runs the complete collection of `fancy_x_plots` for a specific sample ID (`this_sample`), with default parameters.
#' The generated plots are put together into a two-page PDF. In addition, it is also possible to export all individual plots.
#' This can be done by setting `export_individual_plots` to TRUE. It is also possible to use an already loaded seg file instead of using the
#' `this_sample_id` parameter, this is done with the `seg_data` and `maf_data` parameters. Similarly, you can also point this function to a local
#' file on disk with the `seg_path` and `maf_path` parameters.
#'
#' @param this_sample_id Sample ID to be plotted in report.
#' @param export_individual_plots Boolean parameter, set to TRUE to export individual plots.
#' @param out Path to output folder. The default is the working directory.
#' @param seg_data Optional parameter with copy number df already loaded into R.
#' @param seg_path Optional parameter with path to external cn file.
#' @param maf_data Optional parameter with maf like df already loaded into R.
#' @param maf_path Optional parameter with path to external maf like file.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param projection Specify the projection you want the returned plots to be in reference to.
#' Possible values are "grch37" and "hg38". Default is grch37.
#'
#' @return Nothing.
#'
#' @rawNamespace import(ggpubr, except = "get_legend")
#' @import ggplot2 dplyr GAMBLR.utils
#' @export
#'
#' @examples
#' \dontrun{
#' # create a PDF report for one sample, as well as exporting all 
#' # individual plots.
#' comp_report(this_sample_id = "HTMCP-01-06-00422-01A-01D",
#'             out = "./",
#'             export_individual_plots = TRUE)
#' }
#'
#' @keywords internal
comp_report = function(this_sample_id,
                       export_individual_plots = FALSE,
                       out = "./",
                       seg_data,
                       seg_path = NULL,
                       maf_data,
                       maf_path = NULL,
                       this_seq_type = "genome",
                       projection = "grch37"){
  .Deprecated(msg = "comp_report() is deprecated and will be removed in a future version of GAMBLR.viz.")

  if(!missing(maf_data)){
    maf = maf_data
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col_maf] = "Variant_Type"
    colnames(maf)[chromosome_col_maf] = "Chromosome"
    colnames(maf)[start_col_maf] = "Start_Position"
    colnames(maf)[end_col_maf] = "End_Position"

  }else if (!is.null(maf_path)){
    maf = GAMBLR.utils::fread_maf(maf_path)
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col_maf] = "Variant_Type"
    colnames(maf)[chromosome_col_maf] = "Chromosome"
    colnames(maf)[start_col_maf] = "Start_Position"
    colnames(maf)[end_col_maf] = "End_Position"
  }

  if(!missing(seg_data)){
    seg = seg_data
    seg = as.data.frame(seg)
    colnames(seg)[chrom_col_seg] = "chrom"
    colnames(seg)[start_col_seg] = "start"
    colnames(seg)[end_col_seg] = "end"
    colnames(seg)[cn_col_seg] = "CN"

  }else if (!is.null(seg_path)){
    seg = read.table(seg_path, sep = "\t", header = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
    seg = as.data.frame(seg)
    colnames(seg)[chrom_col_seg] = "chrom"
    colnames(seg)[start_col_seg] = "start"
    colnames(seg)[end_col_seg] = "end"
    colnames(seg)[cn_col_seg] = "CN"
  }

  #read maf and seg data into r (avoid calling assign_cn_to_ssm and get_cn_segments for every plotting function)
  if(missing(maf_data) && is.null(maf_path)){
    maf = get_ssm_by_sample(this_sample_id = this_sample_id,
                            this_seq_type = this_seq_type,
                            projection = projection)
  }

  if(missing(seg_data) && is.null(seg_path)){
    seg = get_sample_cn_segments(
      these_sample_ids = this_sample_id,
      streamlined = FALSE,
      this_seq_type = this_seq_type,
      projection = projection
    )
  }

  # check whether maf and seg are empty
  is_maf_seq_empty <- c( nrow(maf), nrow(seg) ) %>%
    {. == 0}
  if( all(is_maf_seq_empty) ){
    stop("The report could not be generated because neither SSMs nor CN segments were found for the specified sample.")
  }
  if( sum(is_maf_seq_empty) == 1 ){
    k <- c("SSMs", "CN segments") [is_maf_seq_empty] %>%
      gettextf("The report could not be generated because %s were not found for the specified sample.", .)
    stop(k)
  }

  #execute a collection of sample-level plots with default parameters
  #page 1
  ssm_chr = fancy_v_chrcount(this_sample_id = this_sample_id, maf_data = maf, plot_title = "", plot_subtitle = "A. SSM Distribution Per Chromosome.", hide_legend = TRUE)
  sv_chr = fancy_v_chrcount(this_sample_id = this_sample_id, plot_title = "", plot_subtitle = "B. SV Distribution Per Chromosome.", ssm = FALSE, hide_legend = TRUE)
  ssm_count = fancy_v_count(this_sample_id = this_sample_id,  maf_data = maf, plot_title = "", plot_subtitle = "C. SSM Counts.", hide_legend = TRUE)
  violine_plot = fancy_v_sizedis(this_sample_id = this_sample_id,  maf_data = maf, plot_title = "", plot_subtitle = "D. SSM Size Distributions.")
  sv_count = fancy_v_count(this_sample_id = this_sample_id, plot_title = "", plot_subtitle = "E. SV Counts.", ssm = FALSE, variant_select = c("DEL", "DUP"), hide_legend = TRUE)
  sv_size = fancy_sv_sizedens(this_sample_id = this_sample_id, plot_title = "", plot_subtitle = "F. SV Size Density.", hide_legend = TRUE)
  snv_plot = fancy_snv_chrdistplot(this_sample_id = this_sample_id,  maf_data = maf, plot_title = "", plot_subtitle = "G. SNV Distribution Per Chromosome.")
  cns = fancy_cnbar(this_sample_id = this_sample_id, seg_data = seg, plot_title = "", plot_subtitle = "H. CN states.")

  #page 2 ideogram
  cnv_ideogram = fancy_ideogram(this_sample_id = this_sample_id, seg_data = seg, maf_data = maf, plot_title = "", plot_subtitle = "F. Ideogram.")

  #build pdf report
  pdf(
    paste0(out, this_sample_id, "_report.pdf"),
    width = 17,
    height = 12,
    onefile = TRUE
)
  # Create a layour for 8 plots where the top row is 2 plots and the rest are
  # in 3-column setup
  ggarrange(
    ggarrange(
        ssm_chr, sv_chr,
        nrow = 1,
        ncol = 2
    ),
    ggarrange(
        ssm_count, violine_plot, sv_count,
        sv_size, snv_plot, cns,
        nrow = 2,
        ncol = 3
    ),
    heights = c(1, 2)
  )
  cnv_ideogram
  dev.off()

  #export individual plots
  if(export_individual_plots){
    ggsave(ssm_chr, filename = paste0(out, this_sample_id, "_ssm_dist_chr.pdf"), limitsize = FALSE, width = 17, height = 12, units = c("in"), dpi = 300)
    ggsave(sv_chr, filename = paste0(out, this_sample_id, "_sv_dist_chr.pdf"), limitsize = FALSE, width = 17, height = 12, units = c("in"), dpi = 300)
    ggsave(snv_plot, filename = paste0(out, this_sample_id, "_snv_dist_chr.pdf"), limitsize = FALSE, width = 17, height = 12, units = c("in"), dpi = 300)
    ggsave(ssm_count, filename = paste0(out, this_sample_id, "_ssm_counts.pdf"), limitsize = FALSE, width = 12, height = 12, units = c("in"), dpi = 300)
    ggsave(sv_count, filename = paste0(out, this_sample_id, "_sv_counts.pdf"), limitsize = FALSE, width = 12, height = 12, units = c("in"), dpi = 300)
    ggsave(sv_size, filename = paste0(out, this_sample_id, "_sv_size_dens.pdf"), limitsize = FALSE, width = 12, height = 12, units = c("in"), dpi = 300)
    ggsave(cns, filename = paste0(out, this_sample_id, "_cn_states.pdf"), limitsize = FALSE, width = 12, height = 12, units = c("in"), dpi = 300)
    ggsave(violine_plot, filename = paste0(out, this_sample_id, "_sv_size_dist.pdf"), limitsize = FALSE, width = 12, height = 12, units = c("in"), dpi = 300)
    ggsave(cnv_ideogram, filename = paste0(out, this_sample_id, "_cnv_ideo.pdf"), limitsize = FALSE, width = 17, height = 12, units = c("in"), dpi = 300)
  }
  return()
}
