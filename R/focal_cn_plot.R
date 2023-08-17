#' @title Copy Number Segments Plot
#'
#' @description Generates a plot of all CN segments for a specified region.
#'
#' @details This function visualizes all CN segments for a defined region, colours the returned segments based on lymphgen information.
#' In addition, this function takes either a specified region (chr:start-end format). If no region is supplied, the user can give the function a gene symbol
#' with `gene`. If so, the function will internally retrieve the region for the specified gene.
#' Sample IDs are specified along the y-axis and the genomic position is visualized along the x-axis.
#'
#' @param this_seg
#' @param this_seg_path
#' @param region Genomic region for plotting in bed format.
#' @param gene Optional variable, converts gene to region if region not supplied.
#' @param these_samples_metadata Required parameter. GAMBL metadata subset to the cases you want to process (or full metadata).
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param type Type of CN segment to be plotted. Default is gain (CN > 2).
#' @param segment_size This parameter controls the size of the segment plotted with ggplot2, default is 1.
#' @param crop_segments Boolean statement that crops segment by first checking if crop segment is smaller than lef/right distance, then adds or subtracts  crop distance to end/start coordinates. Default is TRUE.
#' @param sort_by_annotation Sort CN by annotation, default is "pathology".
#' @param crop_distance Crop distance for cropping segments. Default value is 10000000 bp.
#'
#' @return Nothing
#'
#' @import dplyr cowplot tidyr ggplot2
#' @export
#'
#' @examples
#' #get data
#' 
#' #build plot
#' 
#'
focal_cn_plot = function(this_seg,
                         this_seg_path = NULL,
                         region,
                         gene,
                         these_samples_metadata,
                         this_seq_type = "genome",
                         type = "gain",
                         segment_size = 1,
                         crop_segments = TRUE,
                         sort_by_annotation = c('pathology'),
                         crop_distance = 100000000){
  
  if(missing(these_samples_metadata)){
    stop("Please provide metadata with `these_samples_metadata`...")
  }
  
  if(!missing(gene)){
    region = gene_to_region(gene)
    chunks = region_to_chunks(region)
  }else{
    chunks = region_to_chunks(region)
  }
  
  if(!is.null(this_seg_path)){
    this_seg = read.table(this_seg_path, sep = "\t", header = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  }

  if(!missing(this_seg)){
    if(type == "gain"){
      all_not_dip = this_seg %>%
        mutate(size = end - start) %>%
        dplyr::filter(CN>2)
    }else{
      all_not_dip = this_seg %>%
        mutate(size = end - start) %>%
        dplyr::filter(CN<2)
    }
  }else{
    stop("Please provide a seg file with either `this_seg` or aan absolute path to such a file with `this_seg_path`...")
  }
    

  #crop start and end if they're further than crop_distance from your region
  all_not_dip = mutate(all_not_dip, left_distance = as.numeric(chunks$start) - start)
  all_not_dip = mutate(all_not_dip, right_distance = end - as.numeric(chunks$end))
  
  if(crop_segments){
    all_not_dip = mutate(all_not_dip, end = ifelse(right_distance > crop_distance, as.numeric(chunks$end) + crop_distance, end))
    all_not_dip = mutate(all_not_dip, start = ifelse(left_distance > crop_distance, as.numeric(chunks$start) - crop_distance, start))
  }
  
  all_not_dip = left_join(all_not_dip, these_samples_metadata, by = c("ID" = "sample_id")) %>%
    dplyr::filter(!is.na(pathology))

  all_not_dip = all_not_dip %>%
    arrange(across(all_of(c(sort_by_annotation, "size"))))

  all_not_dip$ID = factor(all_not_dip$ID, levels = unique(all_not_dip$ID))

  ggplot(all_not_dip, aes(x = start, xend = end, y = ID, yend = ID, colour = lymphgen)) +
    geom_vline(aes(xintercept = as.numeric(chunks$start)), alpha = 0.5, colour = get_gambl_colours()[type]) +
    geom_segment(size = segment_size) + theme_cowplot() +
    theme(axis.text.y = element_blank())
}
