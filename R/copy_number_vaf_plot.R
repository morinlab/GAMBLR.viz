#' @title CN VAF Plot
#'
#' @description Create a genome-wide copy number plot for one sample and (optionally) display mutation VAF.
#'
#' @details This function requires a maf (`this_maf`) and a seg file (`this_seg`) to construct the returned plot.
#' This plot is visualizing mutation VAFs per default, this can be turned off with setting `just_segments` to TRUE.
#' This only plots the segments.
#'
#' @param this_maf A maf loaded into R. This is a required paraemter.
#' @param this_seg A seg loaded into R. This is a required paraemter.
#' @param just_segments Specify whether only the segments will be plotted (instead of mutation VAF). Default is FALSE.
#' @param one_chrom Subset plot to one chromosome.
#' @param genes_to_label Optional. Provide a vector of genes to label (if mutated).
#' @param add_chr_prefix If TRUE, "chr" prefix will be added to chr column. Default is FALSE.
#' @param plot_title The title for the created plot.
#'
#' @return Nothing
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr ggplot2 
#' @export
#'
#' @examples
#' #get data
#' dohh2_maf = GAMBLR.data::sample_data$grch37$maf %>% dplyr::filter(Tumor_Sample_Barcode == "DOHH-2")
#' dohh2_seg = GAMBLR.data::sample_data$grch37$seg %>% dplyr::filter(ID == "DOHH-2")
#' 
#' #build plot
#' copy_number_vaf_plot(this_maf = dohh2_maf,this_seg = dohh2_seg, plot_title = "DOHH-2")
#' 
#'
copy_number_vaf_plot = function(this_maf,
                                this_seg,
                                just_segments = FALSE,
                                coding_only = TRUE,
                                one_chrom,
                                genes_to_label,
                                add_chr_prefix = FALSE,
                                plot_title = "My Plot"){

  chrom_order = factor(c(1:22, "X"))
  
  if(add_chr_prefix){
    chrom_order = c(1:22, "X")
    chrom_order = factor(unlist(lapply(chrom_order, function(x){paste0("chr", x)})))
  }
  
  cn_colours = get_gambl_colours(classification = "copy_number")
  
  #get the maf and seg
  if(!missing(this_maf) && !missing(this_seg)){
    maf_sample = as.data.table(this_maf) %>%
      data.table::setkey(Chromosome, Start_Position, End_Position)
    
    seg_sample = as.data.table(this_seg) %>%
      data.table::setkey(chrom, start, end)
  }else{
    stop("Please provide a MAF and SEG file used for plotting...")
  }
  
  maf_with_segs = data.table::foverlaps(maf_sample, seg_sample, type = "any")
  
  #add CN columns to MAF file
  maf_with_segs$log.ratio = tmp_maf$log.ratio
  maf_with_segs$LOH = tmp_maf$LOH_flag
  maf_with_segs$CN = tmp_maf$CN
  
  maf_and_seg = list(maf = maf_with_segs, seg = seg_sample)
  
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
      geom_segment(data = cn_seg, aes(x = start, xend = end, y = CN_segment, yend = CN_segment, colour = CN)) +
      facet_wrap(~chrom, scales = "free_x") +
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
        p + ggtitle(plot_title)
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
        p + ggtitle(plot_title)
      }
    }else{
      p = mutate(vaf_cn_maf, vaf = t_alt_count / (t_ref_count + t_alt_count)) %>%
      ggplot() +
        geom_point(aes(x = Start_Position, y = vaf, colour = CN), alpha = 0.6, size = 0.2) +
        scale_colour_manual(values = cn_colours) +
        facet_wrap(~factor(Chromosome, levels = chrom_order), scales = "free_x") +
        theme_minimal() +
        guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3)))
      p + ggtitle(plot_title)
    }
  }
}

