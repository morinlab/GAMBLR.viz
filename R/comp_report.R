#' @title sample-level SV/SSM/CN reports in PDF
#'
#' @description Construct pdf with sample-level plots, using minimum of arguments
#'
#' @details This function runs the complete collection of `fancy_x_plots` on an incoming MAF (`this_maf`) and SEG file (`this_seg`), with default parameters.
#' The generated plots are put together into a two-page PDF. In addition, it is also possible to export all individual plots.
#' This can be done by setting `export_individual_plots` to TRUE.
#'
#' @param this_seg Required parameter. A SEG file loaded into R.
#' @param this_seg_path Required, if `this_seg` is not provided. An absolute path to the seg file of interest.
#' @param this_maf Required parameter. A MAF file loaded into R.
#' @param this_maf_path Required, if `this_maf` is not provided. An absolute path to the maf file of interest.
#' @param export_individual_plots Boolean parameter, set to TRUE to export individual plots.
#' @param out Path to output folder.
#'
#' @return Nothing.
#'
#' @rawNamespace import(gridExtra, except = "combine")
#' @import ggplot2 dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' #get data for plotting
#' dohh2_maf = GAMBLR.data::sample_data$grch37$maf %>% dplyr::filter(Tumor_Sample_Barcode == "DOHH-2")
#' dohh2_seg = GAMBLR.data::sample_data$grch37$seg %>% dplyr::filter(ID == "DOHH-2")
#' 
#' #create a PDF report for one sample, as well as exporting all individual plots.
#' comp_report(this_seg = dohh2_seg,
#'             this_maf = dohh2_maf,
#'             out = "",
#'             export_individual_plots = TRUE)
#' }
#'
comp_report = function(this_seg,
                       this_seg_path,
                       this_maf,
                       this_maf_path,
                       export_individual_plots = FALSE,
                       out){

  if(!missing(this_maf)){
    maf = this_maf
  }else if(!is.null(this_maf_path)){
    maf = fread_maf(this_maf_path)
  }else{
    stop("Please provide either a maf file (this_maf) or a path to a maf file (this_maf_path)...")
  }

  if(!missing(this_seg)){
    seg = this_seg
  }else if(!is.null(this_seg_path)){
    seg = fread_maf(this_seg_path)
  }else{
    stop("Please provide either a seg file (this_seg) or a path to a seg file (this_seg_path)...")
  }
  
  #ensure only one sample is available in each of the supplied data frames
  maf_sample = unique(maf$Tumor_Sample_Barcode)
  if(length(maf_sample >1)){
    stop("Looks like the provided MAF has more than one smaple ID in it, ensure only one sample ID is present...")
  }
  
  seg_sample = unique(seg$Tumor_Sample_Barcode)
  if(length(seg_sample >1)){
    stop("Looks like the provided SEG has more than one smaple ID in it, ensure only one sample ID is present...")
  }
  
  if(maf_sample != seg_sample){
    stop("The sample ID in thiis_maf does not match the sample ID in this_seg!")
  }

  #execute a collection of sample-level plots with default parameters
  #page 1
  ssm_chr = fancy_v_chrcount(this_maf = maf, plot_title = "", plot_subtitle = "A. SSM Distribution Per Chromosome.", hide_legend = TRUE)
  sv_chr = fancy_v_chrcount(this_maf = maf, plot_title = "", plot_subtitle = "B. SV Distribution Per Chromosome.", ssm = FALSE, hide_legend = TRUE)
  ssm_count = fancy_v_count(this_maf = maf, plot_title = "", plot_subtitle = "C. SSM Counts.", hide_legend = TRUE)
  violine_plot = fancy_v_sizedis(this_maf = maf, plot_title = "", plot_subtitle = "D. SSM Size Distributions.")
  sv_count = fancy_v_count(this_maf = maf, plot_title = "", plot_subtitle = "E. SV Counts.", ssm = FALSE, variant_select = c("DEL", "DUP"), hide_legend = TRUE)
  sv_size = fancy_sv_sizedens(this_maf = maf, plot_title = "", plot_subtitle = "F. SV Size Density.", hide_legend = TRUE)
  snv_plot = fancy_snv_chrdistplot(this_maf = maf, plot_title = "", plot_subtitle = "G. SNV Distribution Per Chromosome.")
  cns = fancy_cnbar(this_maf = maf, this_seg = seg, plot_title = "", plot_subtitle = "H. CN states.")

  #page 2 ideogram
  cnv_ideogram = fancy_ideogram(this_maf = maf, this_seg = seg, plot_title = "", plot_subtitle = "F. Ideogram.")

  #build pdf report
  pdf(paste0(out, maf_sample, "_report.pdf"), width = 17, height = 12)
  page1 = grid.arrange(ssm_chr, sv_chr, ssm_count, violine_plot, sv_count, sv_size, snv_plot, cns, nrow = 3, ncol = 6, name = "Report", top = textGrob(paste0(maf_sample, " - Report"), gp = gpar(fontsize = 15, fontface = "bold")), bottom = "Page 1", layout_matrix = rbind(c(1,1,1,2,2,2), c(3,3,4,4,5,5), c(6,6,7,7,8,8)))
  page2 = grid.arrange(cnv_ideogram,  nrow = 4, ncol = 4, name = "Report", bottom = "Page 2", layout_matrix = rbind(c(1,1,1,1), c(1,1,1,1), c(1,1,1,1), c(1,1,1,1)))
  dev.off()

  #export individual plots
  if(export_individual_plots){
    ggsave(ssm_chr, filename = paste0(out, maf_sample, "_ssm_dist_chr.pdf"), limitsize = FALSE, width = 17, height = 12, units = c("in"), dpi = 300)
    ggsave(sv_chr, filename = paste0(out, maf_sample, "_sv_dist_chr.pdf"), limitsize = FALSE, width = 17, height = 12, units = c("in"), dpi = 300)
    ggsave(snv_plot, filename = paste0(out, maf_sample, "_snv_dist_chr.pdf"), limitsize = FALSE, width = 17, height = 12, units = c("in"), dpi = 300)
    ggsave(ssm_count, filename = paste0(out, maf_sample, "_ssm_counts.pdf"), limitsize = FALSE, width = 12, height = 12, units = c("in"), dpi = 300)
    ggsave(sv_count, filename = paste0(out, maf_sample, "_sv_counts.pdf"), limitsize = FALSE, width = 12, height = 12, units = c("in"), dpi = 300)
    ggsave(sv_size, filename = paste0(out, maf_sample, "_sv_size_dens.pdf"), limitsize = FALSE, width = 12, height = 12, units = c("in"), dpi = 300)
    ggsave(cns, filename = paste0(out, maf_sample, "_cn_states.pdf"), limitsize = FALSE, width = 12, height = 12, units = c("in"), dpi = 300)
    ggsave(violine_plot, filename = paste0(out, maf_sample, "_sv_size_dist.pdf"), limitsize = FALSE, width = 12, height = 12, units = c("in"), dpi = 300)
    ggsave(cnv_ideogram, filename = paste0(out, maf_sample, "_cnv_ideo.pdf"), limitsize = FALSE, width = 17, height = 12, units = c("in"), dpi = 300)
  }
  return()
}
