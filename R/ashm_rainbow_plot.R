#' @title ASHM Rainbow Plot
#'
#' @description Make a rainbow plot of all mutations in a region, ordered and coloured by metadata.
#'
#' @details This function creates a rainbow plot for all mutations in a region specified with 
#'  the `region` parameter for samples specified in `these_samples_metadata`. If 
#'  `these_samples_metadata` is not specified, it will use the samples returned by 
#'  `get_gambl_metadata()` of the appropriate seq types. To highlight and annotate specific 
#'  regions on the plot, an optional bed file of the regions can be specified with `bed`.
#'  An optional MAF-data dataframe with mutations in the `region` can be supplied with
#'  `maf_data`.
#'
#' @param region Genomic region for plotting as a string in the format "chr:start-end".
#' @param maf_data A optional MAF-data data frame containing mutations within the region of 
#'  interest for the samples in `these_samples_metadata`. Ensure it's genome build matches the `projection` parameter.If not supplied, mutations will be retrieved 
#'  from `get_ssm_by_regions` using `these_samples_metadata` and `region` as inputs. 
#' @param these_samples_metadata A data frame with at least columns `sample_id` and the corresponding 
#'  columns given with `classification_column` and `sortByColumns`.
#' @param classification_column The name of the metadata column to use for colouring samples. 
#'  Default: lymphgen.
#' @param exclude_classifications Optional argument for excluding specific classifications within
#'  the `classification_column` values.
#' @param sortByColumns A vector containing the column names you want to sort samples on. The default
#'  is to sort by pathology then `classification_column`.
#' @param projection Plot variants projected to this reference, one of grch37 (default) or hg38. 
#'  Must match the projection of `maf_data` MAF data object if provided.
#' @param drop_unmutated Boolean argument for removing unmutated samples. Default: FALSE.
#' @param bed Optional data frame specifying the regions to annotate. Minimum required column are
#'  start, end, and name with name being used to label the region.
#' @param custom_colours Provide named vector (or named list of vectors) containing custom annotation 
#'  colours if you do not want to use standartized pallette.
#' @param hide_ids Boolean argument, if TRUE (default), sample ids will not appear on the plot.
#' 
#' @return ggplot2 object
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' suppressMessages(library(GAMBLR.open))
#' 
#' # Get lymphgen colours
#' lymphgen_colours = GAMBLR.helpers::get_gambl_colours("lymphgen")
#' 
#' # Basic usage
#' this_region = "chr6:90975034-91066134"
#' this_metadata = get_gambl_metadata()
#' ashm_rainbow_plot(region = this_region,
#'                   these_samples_metadata = metadata, 
#'                   custom_colours = lymphgen_colours)
#' }
#'
ashm_rainbow_plot = function(region,
                             maf_data = NULL,
                             these_samples_metadata,
                             classification_column = "lymphgen",
                             exclude_classifications = NULL,
                             sortByColumns = "pathology",
                             projection = "grch37",
                             drop_unmutated = FALSE,
                             bed,
                             custom_colours,
                             hide_ids = TRUE){

  if(!projection %in% c("grch37", "hg38")){
    stop("projection must be either grch37 or hg38")
  }

  if(missing(these_samples_metadata)){
    # kept for legacy, assumes user provided these_sample_ids
    message("CAUTION! these_samples_metadata was not provided. Using all of get_gambl_metadata().")
    these_samples_metadata = get_gambl_metadata() %>%
      dplyr::filter(seq_type!="mrna")
  }else{
    #drop unsupported seq_type and samples to exclude
    these_samples_metadata = dplyr::filter(these_samples_metadata, seq_type!="mrna")
  }

  # check that classification_column and sortByColumn are in the metadata
  if(!classification_column %in% colnames(these_samples_metadata)){
    stop("classification_column is not present in these_samples_metadata")
  }
  if(!all(sortByColumns %in% colnames(these_samples_metadata))){
    stop("sortByColumns are not all present in these_samples_metadata")
  }

  if(!missing(region)){
    split_chunks = unlist(strsplit(region, ":"))
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = as.numeric(startend[1])
    qend = as.numeric(startend[2])
    if(is.null(maf_data)){
      maf_data = get_ssm_by_regions(regions_list = region,
                                         these_samples_metadata=these_samples_metadata,
                                         projection = projection,
                                         streamlined = TRUE)
    }else{
      # Check that projection and maf_data have the same projection
      if(is.null(get_genome_build(maf_data))){
        stop("No genome_build found for maf_data! Is it a genomic_data or maf_data Object?")
      }
      if(!get_genome_build(maf_data)==projection){
        stop(paste("requested projection:",projection,"and genome_build of maf_data:", 
              get_genome_build(maf_data), "don't match"))
      }
      # Subset MAF data to only the region of interest
      maf_data = get_ssm_by_regions(regions_list = region,
                                        these_samples_metadata=these_samples_metadata,
                                        projection = projection,
                                        streamlined = TRUE, 
                                        maf_data = maf_data)
    }
  }else if(class(region) == "list"| (sum((c("bed_data", "genomic_data", "data.frame") %in% class(region))) >= 1) | length(region)>1 | missing(region)){
    stop("You must supply a single region as a character string")
  }

  sorting_columns = unique(c(sortByColumns, classification_column))
  if(!missing(exclude_classifications)){ 
    meta_arranged = these_samples_metadata %>%
      dplyr::filter(!get(classification_column) %in% exclude_classifications) %>%
      dplyr::arrange(pick(sorting_columns))
  }else{
    meta_arranged = these_samples_metadata %>%
       dplyr::arrange(pick(sorting_columns))
  }

  # Streamlined = TRUE returns col `sample_id` not `Tumor_Sample_Barcode`
  # Renaming it here so that the samples are plotted in the correct order later
  # which is more easily done below when mutations data has a different ID col name 
  # than the arranged metadata
  maf_data <- maf_data %>%
    dplyr::rename(Tumor_Sample_Barcode = sample_id, Start_Position = start)
  mutation_positions = maf_data %>%
    dplyr::select(Tumor_Sample_Barcode, Start_Position) %>%
    as.data.frame()

  mutated_cases = pull(mutation_positions, Tumor_Sample_Barcode) %>%
    unique()

  if(drop_unmutated){
    meta_arranged = meta_arranged %>%
      dplyr::filter(sample_id %in% mutated_cases)
  }
  # Add a fake mutation at the start position to ensure every sample shows
  fake_mutations = data.frame(Tumor_Sample_Barcode = pull(meta_arranged, sample_id), 
                              Start_Position = qstart - 1000)
  mutation_positions = rbind(mutation_positions, fake_mutations)

  meta_arranged$classification = meta_arranged[[classification_column]] %>%
    as.factor()

  muts_anno = dplyr::left_join(
      mutation_positions, 
      meta_arranged, 
      by = c("Tumor_Sample_Barcode" = "sample_id")) %>%
    subset(!is.na(classification))

  muts_anno$sample_id = factor(muts_anno$Tumor_Sample_Barcode, 
                               levels = unique(meta_arranged$sample_id))

  if(missing(custom_colours)){
    p = ggplot(muts_anno) +
      geom_point(aes(x = Start_Position, y = sample_id, colour = classification), 
                alpha = 0.4)
  }else{
    p = ggplot(muts_anno) +
      geom_point(aes(x = Start_Position, y = sample_id, colour = classification), 
                alpha = 0.4) +
      scale_colour_manual(values = custom_colours)
  }
  if(missing(bed)){
    p + guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3)))
  }else{
    bed = bed %>%
      mutate(size = end - start) %>%
      mutate(midpoint = start + size / 2)
    height = length(unique(meta_arranged$sample_id)) + 10
    p = p + geom_rect(data = bed, 
                      aes(xmin = start, xmax = end, ymin = 0, ymax = height + 20), 
                      alpha = 0.1) +
      geom_text(data = bed, 
                aes(x = midpoint, y = height, label = name), 
                size = 2.5, angle = 90) +
      guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3)))
  }

  p = p +
    labs(y = "Sample") +
    GAMBLR.helpers::theme_Morons() +
    theme(plot.margin = ggplot2::margin(1,1,1,1, "cm"), 
          title = element_blank(), 
          plot.subtitle = element_blank(), 
          axis.title.x = element_blank())

  if(hide_ids){
    p = p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }else{
    p = p + theme(axis.text.y = element_text(size = 5))
  }
  return(p)
}
