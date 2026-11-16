#' @title ASHM Multi-panel Rainbow Plot
#'
#' @description Generates a colourful multi-panel overview of hypermutation in regions of interest across many samples.
#'
#' @details Creates a rainbow plot for mutations within regions given in avdata frame in BED format with the following columns; chr, start, end, name.
#' Note that for this function to work, the column names must be exactly this.
#' The user can also to specify a vector of names (`regions_to_display`) to further control what regions are to be displayed on the returned plot.
#' Mutations are coloured by `classification_column`, and certain classifications can be excluded with the parameter `exclude_classifications`.
#' Samples will be plotted in the order they are arranged in `these_samples_metadata`. Otherwise they will be plotted by `classification_column` order.
#' `these_samples_metadata` and `maf_data` are required; the function does not retrieve data itself.
#' For more info, refer to the parameter descriptions of this function.
#'
#' @param regions_bed A data frame in BED format with chromosome coordinates, must contain columns chr, start, end, name (with these exact names) with maximum 10 regions to plot. Not required if specifying `regions_to_display`. with values from the `GAMBLR.data` ashm regions for the given projection, I.e values in `GAMBLR.data::grch37_ashm_regions` or `GAMBLR.data::hg38_ashm_regions`.
#' @param these_samples_metadata Required. A metadata data frame already subsetted and arranged in the order you want the samples vertically displayed.
#' @param maf_data Required. A MAF data frame (e.g. from GAMBLR.results or GAMBLR.open `get_ssm_by_regions()`) covering the regions to plot. Its genome build must match the `projection` parameter. Calls with fewer than 3 alt reads (`t_alt_count`) are not plotted.
#' @param regions_to_display Optional vector of regions in the format "gene-region" that match values of the gene and region columns in the `GAMBLR.data` ashm regions for the given projection, I.e gene and region columns of `GAMBLR.data::grch37_ashm_regions` or `GAMBLR.data::hg38_ashm_regions`.
#' @param this_seq_type Samples in `these_samples_metadata` with this seq_type are plotted. Default: genome.
#' @param classification_column The name of the metadata column to use for colouring samples. 
#'  Default: lymphgen.
#' @param exclude_classifications Optional argument for excluding specific classifications within
#'  the `classification_column` values.
#' @param sortByColumns Ignored; kept for backward compatibility. Samples are plotted in the order of `these_samples_metadata`.
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
#' suppressPackageStartupMessages(library(GAMBLR.open))
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
#' region_maf = GAMBLR.open::get_ssm_by_regions(regions_bed = regions_bed,
#'                                              these_samples_metadata = metadata)
#' ashm_multi_rainbow_plot(regions_bed = regions_bed,
#'                         these_samples_metadata = metadata,
#'                         maf_data = region_maf,
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
#'                         maf_data = genome_maf,
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
  if(missing(these_samples_metadata) || is.null(these_samples_metadata)){
    stop("these_samples_metadata is required (e.g. from get_gambl_metadata() in GAMBLR.results or GAMBLR.open).")
  }
  if(is.null(maf_data)){
    stop("maf_data is required (e.g. from get_ssm_by_regions() in GAMBLR.results or GAMBLR.open).")
  }
  if(!classification_column %in% colnames(these_samples_metadata)){
    stop("classification_column is not present in these_samples_metadata")
  }
  metadata = these_samples_metadata  %>% 
    dplyr::filter(seq_type %in% this_seq_type)
  meta_arranged = metadata #assume the user already arranged it the way they wanted

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
      regions_bed = GAMBLR.utils::create_bed_data(GAMBLR.data::hg38_ashm_regions,
       fix_names="concat",
       concat_cols=c("gene","region"),
       sep="-"
      )
      if(!is.null(regions_to_display)){
        regions_bed <- regions_bed %>%
          dplyr::filter(name %in% regions_to_display)
        if(nrow(regions_bed) == 0){
        stop("regions_to_display do not appear in GAMBLR.data::hg38_ashm_regions")
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
  
  region_mafs = .subset_maf_to_regions(maf_data,
                                       regions_df = as.data.frame(regions_bed),
                                       sample_ids = metadata$sample_id)

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
