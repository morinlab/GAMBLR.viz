#' @title CN VAF Plot
#'
#' @description Create a genome-wide copy number plot for one sample and (optionally) display mutation VAF.
#'
#' @details This function plots copy number segments for one sample and, by default, the VAF of its
#' mutations coloured by the copy number state they fall in. The data are not retrieved internally:
#' supply `maf_data`, `seg_data` and `these_samples_metadata` (e.g. from the accessors in GAMBLR.results or
#' GAMBLR.open). The copy number state of each mutation is assigned with `assign_cn_to_ssm()`.
#' Set `just_segments` to TRUE to plot only the segments. To restrict to coding mutations set
#' `coding_only = TRUE`, optionally with genes to label in `genes_to_label`.
#'
#' @param this_sample_id The sample_id for the sample to plot.
#' @param maf_data Required. A maf_data object (e.g. from `get_ssm_by_samples()`) containing this sample's mutations.
#' @param seg_data Required. A seg_data object (e.g. from `get_cn_segments()`) containing this sample's copy number segments,
#'   in the same genome build as `maf_data`.
#' @param these_samples_metadata Required. Metadata with a row for `this_sample_id`.
#' @param just_segments Specify whether only the segments will be plotted (instead of mutation VAF). Default is FALSE.
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param one_chrom Subset plot to one chromosome.
#' @param genes_to_label Optional. Provide a vector of genes to label (if mutated). Can only be used with coding_only (see above).
#' @param add_chr_prefix If TRUE, "chr" prefix will be added to chr column. Default is FALSE.
#'
#' @return Nothing
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' meta = get_gambl_metadata() %>% dplyr::filter(sample_id == "DOHH-2")
#' maf = get_ssm_by_samples(these_samples_metadata = meta)
#' seg = get_cn_segments(these_samples_metadata = meta)
#' #build plot
#' copy_number_vaf_plot(this_sample_id = "DOHH-2",
#'                      maf_data = maf,
#'                      seg_data = seg,
#'                      these_samples_metadata = meta)
#'
#' #coding only
#' copy_number_vaf_plot(this_sample_id = "DOHH-2",
#'                      maf_data = maf,
#'                      seg_data = seg,
#'                      these_samples_metadata = meta,
#'                      coding_only = TRUE)
#' }
#'
copy_number_vaf_plot = function(this_sample_id,
                                maf_data,
                                seg_data,
                                these_samples_metadata,
                                just_segments = FALSE,
                                coding_only = FALSE,
                                one_chrom,
                                genes_to_label,
                                add_chr_prefix = FALSE){

  if(missing(this_sample_id) || length(this_sample_id) != 1){
    stop("Provide a single this_sample_id.")
  }
  if(missing(maf_data) || missing(seg_data) || missing(these_samples_metadata)){
    stop("maf_data, seg_data and these_samples_metadata are all required ",
         "(e.g. from get_ssm_by_samples(), get_cn_segments() and get_gambl_metadata() in GAMBLR.results or GAMBLR.open).")
  }
  sample_meta = dplyr::filter(these_samples_metadata, sample_id == this_sample_id)
  if(nrow(sample_meta) == 0){
    stop(this_sample_id, " is not in these_samples_metadata.")
  }
  maf_data = dplyr::filter(maf_data, Tumor_Sample_Barcode == this_sample_id)
  seg_data = dplyr::filter(seg_data, ID == this_sample_id)
  if(nrow(seg_data) == 0){
    stop("seg_data has no segments for ", this_sample_id, ".")
  }

  chrom_order = factor(c(1:22, "X"))
  if(add_chr_prefix){
    chrom_order = c(1:22, "X")
    chrom_order = factor(unlist(lapply(chrom_order, function(x){paste0("chr", x)})))
  }
  cn_colours = GAMBLR.helpers::get_gambl_colours(classification = "copy_number")
  maf_and_seg = assign_cn_to_ssm(
    these_samples_metadata = sample_meta,
    maf_data = maf_data,
    seg_data = seg_data,
    coding_only = coding_only
  )
  vaf_cn_maf = maf_and_seg[["maf"]]
  if(!just_segments && nrow(vaf_cn_maf) == 0){
    stop("None of the mutations in maf_data for ", this_sample_id, " fall within a segment in seg_data.")
  }
  vaf_cn_maf = mutate(vaf_cn_maf, CN = case_when(LOH == "1" & CN == 2 ~ "nLOH", TRUE ~ as.character(CN)))
  if(!missing(one_chrom)){
    vaf_cn_maf = dplyr::filter(vaf_cn_maf, Chromosome == one_chrom)
  }
  if(just_segments){
    cn_seg = maf_and_seg[["seg"]]
    cn_seg = mutate(cn_seg, CN_segment = as.numeric(CN), CN = as.character(CN))
    if(!missing(one_chrom)){
      cn_seg = dplyr::filter(cn_seg, Chromosome %in% one_chrom)
    }
    ggplot(cn_seg) +
      geom_segment(data = cn_seg, aes(x = Start_Position, xend = End_Position, y = CN_segment, yend = CN_segment, colour = CN)) +
      facet_wrap(~Chromosome, scales = "free_x") +
      scale_colour_manual(values = cn_colours) +
      theme_minimal() +
      guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3)))
  }else{
    if(coding_only){
      if(missing(genes_to_label)){
        p = mutate(vaf_cn_maf, vaf = t_alt_count / (t_ref_count + t_alt_count)) %>%
        ggplot() +
          geom_point(aes(x = Start_Position, y = vaf, colour = CN), alpha = 0.6, size = 2) +
          scale_colour_manual(values = cn_colours) +
          facet_wrap(~factor(Chromosome, levels = chrom_order), scales = "free_x") +
          theme_minimal() +
          guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3)))
        p + ggtitle(this_sample_id)
      }else{
        #label any mutations that intersect with our gene list
        plot_genes = vaf_cn_maf %>%
          dplyr::filter(Hugo_Symbol %in% genes_to_label)

        p = mutate(vaf_cn_maf, vaf = t_alt_count / (t_ref_count + t_alt_count)) %>%
        ggplot() +
          geom_point(aes(x = Start_Position, y = vaf, colour = CN), size = 2) +
          geom_text(data = plot_genes, aes(x = Start_Position, y = 0.8, label = Hugo_Symbol), size = 3, angle = 90) +
          scale_colour_manual(values = cn_colours) +
          facet_wrap(~factor(Chromosome, levels = chrom_order), scales = "free_x") +
          ylim(c(0,1)) +
          theme_minimal() +
          guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3)))
        p + ggtitle(this_sample_id)
      }
    }else{
      p = mutate(vaf_cn_maf, vaf = t_alt_count / (t_ref_count + t_alt_count)) %>%
      ggplot() +
        geom_point(aes(x = Start_Position, y = vaf, colour = CN), alpha = 0.6, size = 0.2) +
        scale_colour_manual(values = cn_colours) +
        facet_wrap(~factor(Chromosome, levels = chrom_order), scales = "free_x") +
        theme_minimal() +
        guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3)))
      p + ggtitle(this_sample_id)
    }
  }
}
