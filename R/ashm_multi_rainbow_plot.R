#' @title ASHM Multi-panel Rainbow Plot
#'
#' @description Generates a colourful multi-panel overview of hypermutation in regions of interest across many samples.
#'
#' @details Creates a rainbow plot for mutations within regions given in avdata frame in BED format with the following columns; chr, start, end, name.
#' Note that for this function to work, the column names must be exactly this.
#' The user can also to specify a vector of names (`regions_to_display`) to further control what regions are to be displayed on the returned plot.
#' Mutations are coloured by `classification_column`, and certain classifications can be excluded with the parameter `exclude_classifications`.
#' Samples will be plotted in the order they are arranged in `these_samples_metadata`. Otherwise they will be plotted by `classification_column` order.
#' This function will try to obtain mutations internally if `maf_data` is not given.
#' For more info, refer to the parameter descriptions of this function.
#'
#' @param regions_bed A data frame in BED format with chromosome coordinates, must contain columns chr, start, end, name (with these exact names) with maximum 10 regions to plot. Not required if specifying `regions_to_display`. with values from the `GAMBLR.data` ashm regions for the given projection, I.e values in `GAMBLR.data::grch37_ashm_regions` or `GAMBLR.data::hg38_ashm_regions`.
#' @param these_samples_metadata A optional metadata file already subsetted and arranged on the order you want the samples vertically displayed.
#' @param maf_data An optional MAF-data data frame. If not provided, this function will call `get_ssm_by_regions`, using the regions supplied into `regions_bed` and `these_samples_metadata` if provided. Ensure your maf genome build matches the `projection` parameter.
#' @param regions_to_display Optional vector of regions in the format "gene-region" that match values of the gene and region columns in the `GAMBLR.data` ashm regions for the given projection, I.e gene and region columns of `GAMBLR.data::grch37_ashm_regions` or `GAMBLR.data::hg38_ashm_regions`.
#' @param this_seq_type The seq type you want results back for if `maf_data` is not provided. Default: genome.
#' @param classification_column The name of the metadata column to use for colouring samples. 
#'  Default: lymphgen.
#' @param exclude_classifications Optional argument for excluding specific classifications within
#'  the `classification_column` values.
#' @param sortByColumns A vector containing the column names you want to sort samples on. If `these_samples_metadata` is not provided, the default is to sort by pathology then `classification_column`.
#' @param projection Plot variants projected to this reference, one of grch37 (default) or hg38. Bonus regions are only available in grch37.
#' @param custom_colours Provide named vector (or named list of vectors) containing custom annotation colours if you do not want to use standardized pallette.
#' @param verbose  Set to TRUE to maximize the output to console. Default is FALSE.
#'
#' @return Nothing
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' cat("Running example for function: ashm_multi_rainbow_plot\n") 
#' suppressMessages(library(GAMBLR.open))
#' # Get lymphgen colours
#' lymphgen_colours = GAMBLR.helpers::get_gambl_colours("lymphgen")
#' 
#' metadata = suppressMessages(GAMBLR.open::get_gambl_metadata()) %>% 
#'           dplyr::filter(pathology=="DLBCL",
#'                  seq_type=="genome") %>% 
#'           check_and_clean_metadata(.,duplicate_action="keep_first") %>%
#'           dplyr::arrange(lymphgen)
#' regions_bed = GAMBLR.utils::create_bed_data(GAMBLR.data::grch37_ashm_regions,
#'                               fix_names = "concat",
#'                               concat_cols = c("gene","region"),
#'                               sep="-")
#' regions_bed = dplyr::filter(regions_bed,grepl("BCL6",name))
#' ashm_multi_rainbow_plot(regions_bed = regions_bed,
#'                         these_samples_metadata = metadata,
#'                         custom_colours = lymphgen_colours,
#'                         verbose = TRUE)
#' 
#' #build plot
#' \dontrun{
#' ashm_multi_rainbow_plot(regions_to_display = c("BCL2-TSS",
#'                                                "MYC-TSS",
#'                                                "SGK1-TSS",
#'                                                "IGL"),
#'                         these_samples_metadata = metadata,
#'                         custom_colours = lymphgen_colours,
#'                         this_seq_type = "genome")
#' }
#'
ashm_multi_rainbow_plot = function(regions_bed,
                                   these_samples_metadata,
                                   maf_data = NULL,
                                   regions_to_display = NULL,
                                   this_seq_type="genome",
                                   classification_column = "lymphgen",
                                   exclude_classifications,
                                   sortByColumns = "pathology",
                                   projection = "grch37",
                                   custom_colours,
                                   verbose = FALSE) {
  
  if(verbose){
    print("ashm_multi_rainbow_plot")
  }
  if(missing(these_samples_metadata)){
    if(verbose) {
      print("these_samples_metadata was not provided. Using all of get_gambl_metadata().")
    }
    metadata = get_gambl_metadata() %>% 
      dplyr::filter(seq_type %in% this_seq_type)
    
    sorting_columns = unique(c(sortByColumns, classification_column))
    if(!all(sorting_columns %in% colnames(metadata))){
      stop("classification_column or sortByColumns do not exist in the metadata")
    }else{
      meta_arranged = metadata %>% 
        dplyr::arrange(pick(sorting_columns))
    }
  }else{
    metadata = these_samples_metadata  %>% 
      dplyr::filter(seq_type %in% this_seq_type)
    meta_arranged = metadata #assume the user already arranged it the way they wanted
  }

  if(!missing(exclude_classifications)){ 
    meta_arranged = meta_arranged %>%
      dplyr::filter(!get(classification_column) %in% exclude_classifications)
  }

  if(!missing(regions_bed)) {
    print("regions_bed provided")
    if(nrow(regions_bed)>10){
      stop("reduce the regions to display to at most 10")
    }
  }else{
    if(projection == "grch37") {
      regions_bed = GAMBLR.utils::create_bed_data(GAMBLR.data::grch37_ashm_regions,
       fix_names="concat",
       concat_cols=c("gene","region"),
       sep="-"
      )
      if(!is.null(regions_to_display)){
        regions_bed <- regions_bed %>%
          dplyr::filter(name %in% regions_to_display)
        if(nrow(regions_bed) == 0){
        stop("regions_to_display do not appear in GAMBLR.data::grch37_ashm_regions")
        }
      }else{
        # subset to the first 10 rows
        regions_bed <- head(regions_bed, 10)
      }
    }else if(projection == "hg38") {
      regions_bed = GAMBLR.utils::create_bed_data(GAMBLR.data::grch37_ashm_regions,
       fix_names="concat",
       concat_cols=c("gene","region"),
       sep="-"
      )
      if(!is.null(regions_to_display)){
        regions_bed <- regions_bed %>%
          dplyr::filter(name %in% regions_to_display)
        if(nrow(regions_bed) == 0){
        stop("regions_to_display do not appear in GAMBLR.data::grch37_ashm_regions")
        }
      }else{
        # subset to the first 10 rows
        regions_bed <- head(regions_bed, 10)
      }
    }else{
      stop("Please specify one of grch37 or hg38 projections")
    }
  }
  if(verbose){
    print(regions_bed)
  }

  regions_bed = mutate(regions_bed, 
                       region_name=paste0(chrom,":",start,"-",end))
  regions = pull(regions_bed, region_name)
  names = pull(regions_bed, name)
  
  if(missing(maf_data)){
    print("No maf_data provided, using get_ssm_by_regions")
    region_mafs = get_ssm_by_regions(
      regions_bed=regions_bed,
      these_samples_metadata = metadata,
      projection = projection,
      use_name_column = T,
      streamlined = T
    )
  }else{
    region_mafs = get_ssm_by_regions(maf_data=maf_data,
      regions_bed=regions_bed,
      these_samples_metadata = metadata,
      projection = projection,
      use_name_column = T,
      streamlined = T
    )
  }

  if(verbose){
    print(head(region_mafs))
  }
  if(nrow(region_mafs) == 0){
    stop(
      "There are no mutations in the region/sample combination provided."
    )
  }
  
  meta_arranged = meta_arranged %>% mutate_if(is.factor, as.character)
  meta_arranged = meta_arranged %>% mutate(classification = factor(!!sym(classification_column)))
  if(verbose){
    print(head(meta_arranged[,c(1:8)]))
  }
  
  muts_anno = dplyr::left_join(
      region_mafs, 
      meta_arranged)
  muts_first = muts_anno %>%
    dplyr::select(start, region_name) %>%
    dplyr::group_by(region_name) %>%
    dplyr::arrange(start) %>%
    dplyr::filter(row_number() == 1)

  eg = expand_grid(start = pull(muts_first, start),
                   sample_id = pull(meta_arranged, sample_id))
  eg = left_join(eg, muts_first)
  
  #concatenate expanded frame of points with original mutation data
  real_and_fake = bind_rows(region_mafs, eg)
  muts_anno = left_join(real_and_fake, meta_arranged)
  if(any(duplicated(meta_arranged$sample_id))){
    group_by(meta_arranged, sample_id) %>%
      filter(n() > 1) %>%
      print()
    stop("There are duplicated sample ids in the metadata file. Please ensure that each sample id is unique.")
  }
  if(verbose){
    print(head(muts_anno[,c(1:10)]))
  }
  
  muts_anno$sample_id = factor(muts_anno$sample_id, levels = meta_arranged$sample_id)
  if(verbose){
    print("plotting")
  }
  #make the plot
  p = muts_anno %>%
    ggplot() +
    geom_point(aes(x = start,
                   y = sample_id, 
                   colour = classification), alpha = 0.4, size = 0.6) +
    labs(title = "", subtitle = "", x = "", y = "Sample") +
    GAMBLR.helpers::theme_Morons() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = ggplot2::margin(1,1,1,1, "cm"),
          title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title.x = element_blank()) +
    facet_wrap(~region_name, scales = "free_x") +
    guides(color = guide_legend(reverse = TRUE,
                                override.aes = list(size = 3),
                                title={{ classification_column }}))
  
  if(! missing(custom_colours)){
    p = p +
      scale_colour_manual(values = custom_colours)
  }
  if(verbose){
    print("done")
  }
  return(p)
  
}
