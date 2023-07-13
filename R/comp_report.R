#' @title sample-level SV/SSM/CN reports in PDF
#'
#' @description Construct pdf with sample-level plots, using minimum of arguments
#'
#' @details This function runs the complete collection of `fancy_x_plots` for a specific sample ID (`this_sample`), with default parameters.
#' The generated plots are put together into a two-page PDF. In addition, it is also possible to export all individual plots.
#' This can be done by setting `export_individual_plots` to TRUE. It is also possible to use an already loaded seq file instead of using the
#' `this_sample_id` parameter, this is done with the `seq_data` and `maf_data` parameters. Similarly, you can also point this function to a local
#' file on disk with the `seq_path` and `maf_path` parameters.
#'
#' @param this_sample_id Sample ID to be plotted in report.
#' @param export_individual_plots Boolean parameter, set to TRUE to export individual plots.
#' @param out Path to output folder.
#' @param seq_data Optional parameter with copy number df already loaded into R.
#' @param seq_path Optional parameter with path to external cn file.
#' @param maf_data Optional parameter with maf like df already loaded into R.
#' @param maf_path Optional parameter with path to external maf like file.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#'
#' @return Nothing.
#'
#' @rawNamespace import(gridExtra, except = "combine")
#' @import ggplot2 dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' #create a PDF report for one sample, as well as exporting all individual plots.
#' comp_report(this_sample = "HTMCP-01-06-00422-01A-01D",
#'             out = "reports/",
#'             export_individual_plots = TRUE)
#' }
#'
comp_report = function(this_sample_id,
                       export_individual_plots = FALSE,
                       out,
                       seq_data,
                       seq_path = NULL,
                       maf_data,
                       maf_path = NULL,
                       this_seq_type = "genome"){

  if(!missing(maf_data)){
    maf = maf_data
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col_maf] = "Variant_Type"
    colnames(maf)[chromosome_col_maf] = "Chromosome"
    colnames(maf)[start_col_maf] = "Start_Position"
    colnames(maf)[end_col_maf] = "End_Position"

  }else if (!is.null(maf_path)){
    maf = fread_maf(maf_path)
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col_maf] = "Variant_Type"
    colnames(maf)[chromosome_col_maf] = "Chromosome"
    colnames(maf)[start_col_maf] = "Start_Position"
    colnames(maf)[end_col_maf] = "End_Position"
  }

  if(!missing(seq_data)){
    seq = seq_data
    seq = as.data.frame(seq)
    colnames(seq)[chrom_col_seq] = "chrom"
    colnames(seq)[start_col_seq] = "start"
    colnames(seq)[end_col_seq] = "end"
    colnames(seq)[cn_col_seq] = "CN"

  }else if (!is.null(seq_path)){
    seq = read.table(seq_path, sep = "\t", header = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
    seq = as.data.frame(seq)
    colnames(seq)[chrom_col_seq] = "chrom"
    colnames(seq)[start_col_seq] = "start"
    colnames(seq)[end_col_seq] = "end"
    colnames(seq)[cn_col_seq] = "CN"
  }

  #read maf and seq data into r (avoid calling assign_cn_to_ssm and get_cn_segments for every plotting function)
  if(missing(maf_data) && is.null(maf_path)){
    maf = assign_cn_to_ssm(
      this_sample_id = this_sample_id,
      coding_only = FALSE,
      from_flatfile = TRUE,
      use_augmented_maf = TRUE,
      this_seq_type = this_seq_type)$maf
  }

  if(missing(seq_data) && is.null(seq_path)){
    seq = get_sample_cn_segments(
      this_sample_id = this_sample_id,
      multiple_samples = FALSE,
      streamlined = FALSE,
      from_flatfile = TRUE,
      this_seq_type = this_seq_type
    )
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
  cns = fancy_cnbar(this_sample_id = this_sample_id, seq_data = seq, plot_title = "", plot_subtitle = "H. CN states.")

  #page 2 ideogram
  cnv_ideogram = fancy_ideogram(this_sample_id = this_sample_id, seq_data = seq, maf_data = maf, plot_title = "", plot_subtitle = "F. Ideogram.")

  #build pdf report
  pdf(paste0(out, this_sample_id, "_report.pdf"), width = 17, height = 12)
  page1 = grid.arrange(ssm_chr, sv_chr, ssm_count, violine_plot, sv_count, sv_size, snv_plot, cns, nrow = 3, ncol = 6, name = "Report", top = textGrob(paste0(this_sample, " - Report"), gp = gpar(fontsize = 15, fontface = "bold")), bottom = "Page 1", layout_matrix = rbind(c(1,1,1,2,2,2), c(3,3,4,4,5,5), c(6,6,7,7,8,8)))
  page2 = grid.arrange(cnv_ideogram,  nrow = 4, ncol = 4, name = "Report", bottom = "Page 2", layout_matrix = rbind(c(1,1,1,1), c(1,1,1,1), c(1,1,1,1), c(1,1,1,1)))
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
