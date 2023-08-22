#' @title Sample-level Circos Plot
#'
#' @description Plot a sample-centric circos overview.
#'
#' @details This function lets the user supply already loaded data frames (SV, CNV, SSM) with the `this_bedpe` and`this_seg` parameters.
#' Convenient Boolean parameters are also available for restricting the plot to specific mutation types (`include_sv` and `include_cnv`).
#'
#' @param this_bedpe Required data frame of SVs. If not provided this function will run `get_manta_sv` to retrieve SVs.
#' @param this_seg Required data frame of CNVs. If not provided, this function will run `get_sample_cn_segments` to retrieve CNVs.
#' @param this_maf This parameter does not do anything yet. Maybe it was meant to be implemented.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param legend_metadata_columns Column names from metadata
#' @param legend_metadata_names List of metadata names to be plotted.
#' @param chrom_list List of chromosomes to be plotted. If not stated, chr1-22+X will bes used.
#' @param label_genes Gene labels.
#' @param auto_label_sv Default is FALSE
#'
#' @return Nothing
#'
#' @import dplyr circlize ComplexHeatmap ggplot2
#' @export
#'
#' @examples
#' #get data
#' dohh2_seg = GAMBLR.data::sample_data$grch37$seg %>% dplyr::filter(ID == "DOHH-2")
#' dohh2_bedpe = GAMBLR.data::sample_data$grch37$bedpe %>% dplyr::filter(tumour_sample_id == "DOHH-2")
#' 
#' #build plot
#' plot_sample_circos(this_bedpe = dohh2_bedpe,
#'                    this_seg = dohh2_seg,
#'                    chrom_list = c("chr2",
#'                                   "chr3",
#'                                   "chr8",
#'                                   "chr14",
#'                                   "chr18"))
#'
plot_sample_circos = function(this_bedpe,
                              this_seg,
                              legend_metadata_columns,
                              legend_metadata_names = c(),
                              this_seq_type = "genome",
                              chrom_list,
                              label_genes,
                              auto_label_sv = FALSE){

  add_cnv = function(this_seg){
    bed = data.frame(this_seg[,c("chrom", "start", "end", "log.ratio")])
    colnames(bed) = c("chr", "start", "end", "value1")
    col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    circos.genomicTrackPlotRegion(bed, stack = TRUE, panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = col_fun(value), border = NA, posTransform = NULL, ...)
      i = getI(...)
      cell.xlim = get.cell.meta.data("cell.xlim")
    }, bg.border = NA)
  }
  if(missing(this_seg)){
    stop("Please provide a data frame with CNVs...")
  }
  if(missing(chrom_list)){

  #should we add chr list for males as well? Sex could then also be added as a paramter for the function were the appropiate chr list is called if not stated?
   chrom_list = paste0("chr", c(1:22,"X"))
  }
  if(!missing(label_genes)){
    gene_bed = GAMBLR.data::grch37_gene_coordinates %>%
      dplyr::filter(gene_name %in% label_genes) %>%
      dplyr::select(chromosome, start, end, gene_name) %>%
      dplyr::mutate(chromosome = paste0("chr", chromosome))
  }
  if(missing(this_bedpe)){
    stop("Please provide a data frame with SVs...")
  }
  
  #add chr prefixes if grch37 is selcted (expected by circlize)
  if(all(!str_detect(this_bedpe$CHROM_A, "chr"))){
    sv_df = this_bedpe %>%
      dplyr::mutate(CHROM_A = paste0("chr", CHROM_A)) %>%
      dplyr::mutate(CHROM_B = paste0("chr", CHROM_B))
  }

  sv_df = sv_df %>%
    dplyr::filter(CHROM_A %in% chrom_list) %>%
    dplyr::filter(CHROM_B %in% chrom_list)

  if(auto_label_sv){
    #annotate oncogene SVs and label them
    annotated_sv  = GAMBLR.utils::annotate_sv(sv_df, with_chr_prefix = TRUE) %>%
      dplyr::filter(!is.na(partner))

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

    bed_mut_partner = GAMBLR.data::grch37_partners %>%
      dplyr::filter(gene %in% these_partners) %>%
      mutate(chrom = paste0("chr", chrom))

    bed_mut_onco = GAMBLR.data::grch37_oncogene %>%
      dplyr::filter(gene %in% these_oncogenes) %>%
      mutate(chrom = paste0("chr", chrom))

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
  add_cnv(this_seg)
  circos.genomicLink(bed1, bed2, col = "#bdbdc1")
  if(!missing(label_genes)){
    circos.genomicLabels(gene_bed, labels.column = "gene_name")
  }
  if(auto_label_sv){
    circos.genomicLink(anno_bed1, anno_bed2,col = 'red')
    circos.genomicLabels(bed_mut, labels.column = "gene")
  }
  
  this_sample = unique(this_bedpe$tumour_sample_id)
  text(c(0.75, 0.75), this_sample, cex = 0.8)
  
  if(!missing(legend_metadata_columns)){
      samp_meta = GAMBLR.helpers::handle_metadata(this_seq_type = this_seq_type) %>%
        dplyr::filter(sample_id == this_sample)
    
    these_meta = samp_meta[legend_metadata_columns]
    these_cols = get_gambl_colours()
    vals = as.character(these_meta)
    names = colnames(these_meta)

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
