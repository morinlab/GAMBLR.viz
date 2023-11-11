#' @title CN VAF Plot
#'
#' @description Create a genome-wide copy number plot for one sample and (optionally) display mutation VAF.
#'
#' @details This function takes a sample ID and internally calls [GAMBLR::assign_cn_to_ssm] to get copy number segments for plotting.
#' This plot is visualizing mutation VAFs per default, this can be turned off with setting `just_segments` to TRUE.
#' This only plots the segments. The user can also restrict the plotted segments to coding regions. To do so, set `coding_only= TRUE`,
#' and then specify the genes of interest (coding regions) with the `genes_to_label` (vector of genes).
#'
#' @param this_sample_id The sample_id for the sample to plot.
#' @param just_segments Specify whether only the segments will be plotted (instead of mutation VAF). Default is FALSE.
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param one_chrom Subset plot to one chromosome.
#' @param genes_to_label Optional. Provide a vector of genes to label (if mutated). Can only be used with coding_only (see above).
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param from_flatfile If set to true the function will use flatfiles instead of the database.
#' @param use_augmented_maf Boolean statement if to use augmented maf, default is TRUE.
#' @param add_chr_prefix If TRUE, "chr" prefix will be added to chr column. Default is FALSE.
#'
#' @return Nothing
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' #build plot
#' copy_number_vaf_plot(this_sample_id = "DOHH-2")
#'
#' #coding only
#' copy_number_vaf_plot(this_sample_id = "DOHH-2",
#'                      coding_only = TRUE)
#' }
#'
copy_number_vaf_plot = function(this_sample_id,
                                just_segments = FALSE,
                                coding_only = FALSE,
                                one_chrom,
                                genes_to_label,
                                this_seq_type = "genome",
                                from_flatfile = TRUE,
                                use_augmented_maf = TRUE,
                                add_chr_prefix = FALSE){

  chrom_order = factor(c(1:22, "X"))
  if(add_chr_prefix){
    chrom_order = c(1:22, "X")
    chrom_order = factor(unlist(lapply(chrom_order, function(x){paste0("chr", x)})))
  }
  cn_colours = GAMBLR.helpers::get_gambl_colours(classification = "copy_number")
  maf_and_seg = assign_cn_to_ssm(
    this_sample_id = this_sample_id,
    coding_only = coding_only,
    this_seq_type = this_seq_type
  )
  vaf_cn_maf = maf_and_seg[["maf"]]
  vaf_cn_maf = mutate(vaf_cn_maf, CN = case_when(LOH == "1" & CN == 2 ~ "nLOH", TRUE ~ as.character(CN)))
  if(!missing(one_chrom)){
    vaf_cn_maf = dplyr::filter(vaf_cn_maf, Chromosome == one_chrom)
  }
  if(just_segments){
    cn_seg = maf_and_seg[["seg"]]
    cn_seg = mutate(cn_seg, CN_segment = as.numeric(CN), CN = as.character(CN))
    print(head(cn_seg))
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
