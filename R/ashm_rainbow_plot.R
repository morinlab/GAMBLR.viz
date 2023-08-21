#' @title ASHM Rainbow Plot
#'
#' @description Make a rainbow plot of all mutations in a region, ordered and coloured by metadata.
#'
#' @details This function creates a rainbow plot for all mutations in a region. Region can either be specified with the `region` parameter,
#' or the user can provide a maf that has already been subset to the region(s) of interest with `this_maf`.
#' Regardless, the user has to provide a maf or maf-like object with `this_maf` as an required parameter.
#' As a third alternative, the regions can also be specified as a bed file with `bed`.
#' Lastly, this function has a variety of parameters that can be used to further customize the returned plot in many different ways.
#' Refer to the parameter descriptions, examples as well as the vignettes for more demonstrations how this function can be called.
#'
#' @param this_maf Required parameter. A data frame containing mutations (MAF format) within a region of interest.
#' @param metadata A data frame with sample_id as a column.
#' @param exclude_classifications Optional argument for excluding specific classifications from a metadata file.
#' @param drop_unmutated Boolean argument for removing unmutated sample ids in mutated cases.
#' @param classification_column The name of the metadata column to use for ordering and colouring samples.
#' @param bed Optional data frame specifying the regions to annotate (required columns: start, end, name).
#' @param region Genomic region for plotting in bed format.
#' @param custom_colours Provide named vector (or named list of vectors) containing custom annotation colours if you do not want to use standartized pallette.
#' @param hide_ids Boolean argument, if TRUE, ids will be removed.
#'
#' @return ggplot2 object
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' #get data
#' dohh2_maf = GAMBLR.data::sample_data$grch37$maf %>% dplyr::filter(Tumor_Sample_Barcode == "DOHH-2")
#' dohh2_meta = GAMBLR.data::gambl_metadata %>% dplyr::filter(sample_id == "DOHH-2")
#'
#' #build plot
#' ashm_rainbow_plot(this_maf = dohh2_maf,
#'                   metadata = dohh2_meta, 
#'                   region = "chr6:90975034-91066134")
#'
ashm_rainbow_plot = function(this_maf,
                             metadata,
                             exclude_classifications,
                             drop_unmutated = FALSE,
                             classification_column,
                             bed,
                             region,
                             custom_colours,
                             hide_ids = TRUE){
  
  if(!region == ""){
    region = gsub(",", "", region)
    split_chunks = unlist(strsplit(region, ":"))
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = as.numeric(startend[1])
    qend = as.numeric(startend[2])
  }
  
  if(!missing(this_maf)){
    regions_maf = GAMBLR.helpers::handle_ssm_by_region(this_maf = this_maf, region = region)
  }else{
    stop("Please provide a MAF, or MAF-like object with the `this_maf` parameter...")
  }

  if(!missing(classification_column)){
    meta_arranged = arrange(metadata, pathology_rank, lymphgen)
    if(!missing(exclude_classifications)){
      meta_arranged = dplyr::filter(meta_arranged,!get(classification_column) %in% exclude_classifications)
    }
  }else{
    classification_column = "lymphgen"
    meta_arranged = metadata
  }
  
  mutation_positions = regions_maf %>%
    dplyr::select(Tumor_Sample_Barcode, Start_Position) %>%
    as.data.frame()

  mutated_cases = pull(mutation_positions, Tumor_Sample_Barcode) %>%
    unique()

  if(drop_unmutated){
    meta_arranged = meta_arranged %>%
      dplyr::filter(sample_id %in% mutated_cases)
  }
  
  #add a fake mutation at the start position for each sample to ensure every sample shows up
  fake_mutations = data.frame(Tumor_Sample_Barcode = pull(metadata, sample_id), Start_Position = qstart - 1000)
  mutation_positions = rbind(mutation_positions, fake_mutations)

  meta_arranged$classification = meta_arranged[[classification_column]] %>%
    as.factor()

  muts_anno = dplyr::left_join(mutation_positions, meta_arranged, by = c("Tumor_Sample_Barcode" = "sample_id")) %>%
    subset(!is.na(classification))

  muts_anno$sample_id = factor(muts_anno$Tumor_Sample_Barcode, levels = unique(meta_arranged$sample_id))

  if(missing(custom_colours)){
    p = ggplot(muts_anno) +
      geom_point(aes(x = Start_Position, y = sample_id, colour = classification), alpha = 0.4)
  }else{
    p = ggplot(muts_anno) +
      geom_point(aes(x = Start_Position, y = sample_id, colour = classification), alpha = 0.4) +
      scale_colour_manual(values = custom_colours)
  }
  
  if(missing(bed)){
    p + guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3)))
  }else{
    bed = bed %>%
      mutate(size = end - start) %>%
      mutate(midpoint = start + size / 2)
    height = length(unique(meta_arranged$sample_id)) + 10
    p = p + geom_rect(data = bed, aes(xmin = start, xmax = end, ymin = 0, ymax = height + 20), alpha = 0.1) +
      geom_text(data = bed, aes(x = midpoint, y = height, label = name), size = 2.5, angle = 90) +
      guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3)))
  }

  p = p +
    labs(y = "Sample") +
    theme_Morons() +
    theme(plot.margin = ggplot2::margin(1,1,1,1, "cm"), title = element_blank(), plot.subtitle = element_blank(), axis.title.x = element_blank())

  if(hide_ids){
    p = p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }else{
    p = p + theme(axis.text.y = element_text(size = 5))
  }
  return(p)
}
