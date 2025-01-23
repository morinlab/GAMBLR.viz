#' @title Sample-level Circos Plot
#'
#' @description Plot a sample-centric circos overview.
#'
#' @details This function takes a sample ID in the `this_sample_id` parameter.
#' Optionally, the user can supply already loaded data frames (SV, CNV, SSM) with the `sv_df`, `cnv_df` and `ssm_df` parameters.
#' Convenient Boolean parameteers are also avaialble for restricting the plot to specific mutation types (`include_sv`, `include_cnv`, and `include_ssm`).
#'
#' @param this_sample_id Sample ID for the sample to plot.
#' @param sv_df Optional data frame of SVs. If not provided this function will run `get_manta_sv` to retrieve SVs.
#' @param cnv_df Optional data frame of CNVs. If not provided, this function will run `get_sample_cn_segments` to retrieve CNVs.
#' @param ssm_df This parameter does not do anything yet. Maybe it was meant to be implemented.
#' @param include_sv Default TRUE. (does not do anything yet).
#' @param include_cnv Default TRUE. (does not do anything yet).
#' @param this_projection The selected projection, default is grch37 and it's the only supported peojection.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param include_ssm Defaul FALSE. (does not do anything yet).
#' @param legend_metadata_columns Column names from metadata
#' @param legend_metadata_names List of metadata names to be plotted.
#' @param chrom_list List of chromosomes to be plotted. If not stated, chr1-22+X will be used.
#' @param label_genes Gene labels (df, list or what type?)
#' @param auto_label_sv Default is FALSE
#' @param auto_colour_links Whether to apply authomatic coloring of the links. Default is FALSE.
#' @param hide_legend Set to TRUE if you want to suppress the legend. Particularly useful if you are not using GAMBL data/metadata
#'
#' @return Nothing
#'
#' @import dplyr circlize ComplexHeatmap ggplot2
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#'
#' plot_sample_circos(this_sample_id = "02-13135T",
#'                    legend_metadata_columns = c("pathology",
#'                                                "lymphgen",
#'                                                "COO_consensus",
#'                                                "DHITsig_consensus"),
#'                    legend_metadata_names = c("pathology",
#'                                              "LymphGen",
#'                                              "COO",
#'                                              "DHITsig"),
#'                    chrom_list = c("chr3",
#'                                   "chr8",
#'                                   "chr14",
#'                                   "chr18"))
#'
plot_sample_circos = function(this_sample_id,
                              sv_df,
                              cnv_df,
                              ssm_df,
                              include_sv = TRUE,
                              include_ssm = FALSE,
                              legend_metadata_columns,
                              legend_metadata_names = c(),
                              include_cnv = TRUE,
                              this_projection = "grch37",
                              this_seq_type = "genome",
                              chrom_list,
                              label_genes,
                              auto_label_sv = FALSE,
                              auto_colour_links = FALSE,
                              hide_legend = FALSE){



  add_cnv = function(cnv_df){
    bed = data.frame(cnv_df[,c("chrom", "start", "end", "log.ratio")])
    colnames(bed) = c("chr", "start", "end", "value1")
    col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    circos.genomicTrackPlotRegion(bed, stack = TRUE, panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = col_fun(value), border = NA, posTransform = NULL, ...)
      i = getI(...)
      cell.xlim = get.cell.meta.data("cell.xlim")
    }, bg.border = NA)
  }
  if(include_cnv & missing(cnv_df)){
    cnv_df = get_sample_cn_segments(
      these_sample_ids = this_sample_id,
      with_chr_prefix = TRUE,
      this_seq_type = this_seq_type
    )
  }
  if(missing(chrom_list)){

  #should we add chr list for males as well? Sex could then also be added as a paramter for the function were the appropiate chr list is called if not stated?
   chrom_list = paste0("chr", c(1:22,"X"))
  }
  if(!missing(label_genes)){
    if(this_projection=="grch37"){
      gene_bed = GAMBLR.data::grch37_gene_coordinates %>%
        dplyr::filter(gene_name %in% label_genes) %>%
        dplyr::select(chromosome, start, end, gene_name) %>%
        dplyr::mutate(chromosome = paste0("chr", chromosome))
    }else if(this_projection=="hg38"){
      gene_bed = GAMBLR.data::hg38_gene_coordinates %>%
        dplyr::filter(gene_name %in% label_genes) %>%
        dplyr::select(chromosome, start, end, gene_name)
    }

  }
  if(include_sv){
    if(missing(sv_df)){
      sv_df = get_manta_sv(verbose = FALSE) %>%
        dplyr::filter(tumour_sample_id == this_sample_id)
    }else{
      sv_df = sv_df %>%
        dplyr::filter(tumour_sample_id == this_sample_id)
    }
  }


  #add chr prefixes if grch37 is selcted (expected by circlize)
  if(this_projection == "grch37"){
    sv_df = sv_df %>%
      dplyr::mutate(CHROM_A = paste0("chr", CHROM_A)) %>%
      dplyr::mutate(CHROM_B = paste0("chr", CHROM_B))
  }

  sv_df = sv_df %>%
    dplyr::filter(CHROM_A %in% chrom_list) %>%
    dplyr::filter(CHROM_B %in% chrom_list)

  if(auto_label_sv){

    #annotate oncogene SVs and label them
    annotated_sv  = annotate_sv(sv_df, with_chr_prefix = TRUE,genome_build = this_projection) %>%
      dplyr::filter(!is.na(partner)) %>%
      dplyr::filter(tumour_sample_id == this_sample_id)

    these_oncogenes = unique(pull(annotated_sv, gene))
    these_partners = unique(pull(annotated_sv, partner))
    if("IGH" %in% these_partners){
      these_partners = c(these_partners, "IGHV3-62")
    }
    anno_bed1 = annotated_sv %>%
      dplyr::select(chrom1, start1, end1, tumour_sample_id)

    anno_bed2 = annotated_sv %>%
      dplyr::select(chrom2, start2, end2, tumour_sample_id)

    colnames(anno_bed1) = c("chrom", "start", "end", "sample_id")
    colnames(anno_bed2) = c("chrom", "start", "end", "sample_id")

    if(this_projection=="grch37"){
      bed_mut_partner = GAMBLR.data::grch37_partners %>%
        dplyr::filter(gene %in% these_partners) %>%
        mutate(chrom = paste0("chr", chrom))

      bed_mut_onco = GAMBLR.data::grch37_oncogene %>%
        dplyr::filter(gene %in% these_oncogenes) %>%
        mutate(chrom = paste0("chr", chrom))
    }else if(this_projection=="hg38"){
      bed_mut_partner = GAMBLR.data::hg38_partners %>%
        dplyr::filter(gene %in% these_partners) %>%
        mutate(chrom = paste0("chr", chrom))

      bed_mut_onco = GAMBLR.data::hg38_oncogene %>%
        dplyr::filter(gene %in% these_oncogenes) %>%
        mutate(chrom = paste0("chr", chrom))
    }


    bed_mut = bind_rows(bed_mut_partner, bed_mut_onco)
    print(bed_mut)
  }
  bed1 = sv_df %>%
    dplyr::select(CHROM_A, START_A, END_A, tumour_sample_id)

  bed2 = sv_df %>%
    dplyr::select(CHROM_B, START_B, END_B, tumour_sample_id)

  colnames(bed1) = c("chrom", "start", "end", "sample_id")
  colnames(bed2) = c("chrom", "start", "end", "sample_id")
  circos.clear()
  circos.initializeWithIdeogram(chromosome.index = chrom_list)
  if(include_cnv){
    add_cnv(cnv_df)
  }

  if(auto_colour_links){
    bed2 = decorate_bed(bed2,colour_by="chrom",colour_mapping = get_gambl_colours("chromosome"))
    circos.genomicLink(bed1, bed2, col = bed2$color)
  }else{
    circos.genomicLink(bed1, bed2, col = "#bdbdc1")
  }

  if(!missing(label_genes)){
    circos.genomicLabels(gene_bed, labels.column = "gene_name")
  }
  if(auto_label_sv){
    circos.genomicLink(anno_bed1, anno_bed2,col = 'red')
    circos.genomicLabels(bed_mut, labels.column = "gene")
  }
  text(0.75, this_sample_id, cex = 0.8)
  if(hide_legend){
    return()
  }
  if(!missing(legend_metadata_columns)){
    samp_meta = GAMBLR.helpers::handle_metadata(this_seq_type = this_seq_type) %>%
      dplyr::filter(sample_id == this_sample_id)

    these_meta = samp_meta[legend_metadata_columns]
    these_cols = GAMBLR.helpers::get_gambl_colours()

    these_meta <- mutate_if(these_meta, is.factor, as.character)
    vals = as.character(these_meta)

    all_cols = map_metadata_to_colours(legend_metadata_columns, these_meta, verbose = T)

    cols = all_cols[vals]
    print(cols)
    if(length(legend_metadata_names) == length(legend_metadata_columns)){
      for(i in c(1:length(vals))){
        if(!legend_metadata_names[i] == ""){
          vals[i] = paste(legend_metadata_names[i], vals[i])
        }
      }
    }

    lgd_discrete = Legend(labels = vals, title_position = "topleft", legend_gp = gpar(fill = cols))
    draw(lgd_discrete, x = unit(3, "mm"), y = unit(1, "npc") - unit(5, "mm"), just = c("left", "top"))
  }

  #continuous
  lgd_cnv = Legend(at = c(-2, -1, 0, 1, 2),
                   col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                   title_position = "topleft",
                   title = "log\nratio")

  draw(lgd_cnv, x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))
}
