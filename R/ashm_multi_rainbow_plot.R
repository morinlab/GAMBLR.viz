#' @title ASHM Multi-panel Rainbow Plot
#'
#' @description Generates a colourful multi-panel overview of hypermutation in regions of interest across many samples.
#'
#' @details The input for this function is a bed-file with the following columns; chr, start, end, name.
#' Note that for this function to work, the column names must be exactly this.
#' The user also needs to specify a vector of names (`regions_to_display`) to further control what regions are to be displayed on the returned plot.
#' It is also possible to exclude specific classifications from the metadata file. This is achieved with `exclude_classifications`.
#' In addition the user can also use the `metadata` parameter to use an already subset and arranged metadata table.
#' This function will try to obtain mutations internally if `maf_data` is not given.
#' For more info, refer to the parameter descriptions of this function.
#'
#' @param regions_bed Bed file with chromosome coordinates, should contain columns chr, start, end, name (with these exact names). Not required if selecting from many common regions; bonus regions also exist in grch37.
#' @param these_samples_metadata A metadata file already subsetted and arranged on the order you want the samples vertically displayed.
#' @param this_seq_type the seqtype you want results back for if `maf_data` is not provided.
#' @param custom_colours Provide named vector (or named list of vectors) containing custom annotation colours if you do not want to use standardized pallette.
#' @param classification_column Optional. Override default column for assigning the labels used for colouring in the figure.
#' @param maf_data An already loaded maf. If not provided, this function will call `get_ssm_by_region`, using the regions supplied into `regions_bed`. Ensure your maf matches the genome projection.
#' @param projection Provide genome build; default is grch37. Bonus regions are only available in grch37.
#' @param verbose Set to FALSE to prevent printing the full regions bed file to the console. Default is TRUE.
#'
#' @return Nothing
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' cat("Running example for function: ashm_multi_rainbow_plot\n") 
#' suppressMessages(library(GAMBLR.open))
#' #get lymphgen colours
#' lymphgen_colours = GAMBLR.helpers::get_gambl_colours("lymphgen")
#' 
#' metadata = suppressMessages(GAMBLR.open::get_gambl_metadata()) %>% 
#'           dplyr::filter(pathology=="DLBCL",
#'                  seq_type=="genome") %>% 
#'           check_and_clean_metadata(.,duplicate_action="keep_first") %>%
#'           dplyr::arrange(lymphgen)
#' regions_bed = GAMBLR.utils::create_bed_data(grch37_ashm_regions,
#'                               fix_names = "concat",
#'                               concat_cols = c("gene","region"))
#' regions_bed = dplyr::filter(regions_bed,grepl("BCL6",name))
#' ashm_multi_rainbow_plot(regions_bed,
#'                         metadata,
#'                         custom_colours = lymphgen_colours,
#'                         verbose = TRUE)
#' 
#' #build plot
#' \dontrun{
#' ashm_multi_rainbow_plot(regions_to_display = c("BCL2-TSS",
#'                                                "MYC-TSS",
#'                                                "SGK1-TSS",
#'                                                "IGL"),
#'                         custom_colours = lymphgen_colours,
#'                         this_seq_type = "genome")
#' }
#'
ashm_multi_rainbow_plot = function(regions_bed,
                                   these_samples_metadata,
                                   this_seq_type="genome",
                                   custom_colours,
                                   classification_column = "lymphgen",
                                   maf_data,
                                   projection = "grch37",
                                   verbose = FALSE) {
  
  if(verbose){
    print("ashm_multi_rainbow_plot")
  }
  if(missing(these_samples_metadata)){
    if(verbose) {
      print("finding metadata")
    }
    metadata = get_gambl_metadata() %>% 
      dplyr::filter(seq_type %in% this_seq_type)
    meta_arranged = arrange(metadata, pathology, lymphgen)
  }else{
    metadata = these_samples_metadata  %>% 
      dplyr::filter(seq_type %in% this_seq_type)
    meta_arranged = metadata #assume the user already arranged it the way they wanted
  }

  
  if(!missing(regions_bed)) {
    print("regions_bed provided")
    if(nrow(regions_bed)>10){
      stop("reduce the regions to display to at most 10")
    }
  } else {
    if (projection == "grch37") {
      regions_bed = GAMBLR.utils::create_bed_data()(GAMBLR.data::grch37_ashm_regions,
       fix_names="concat",
       concat_cols=c("gene","region"),
       sep="-"
      )
    } else if (projection == "hg38") {
      regions_bed = GAMBLR.utils::create_bed_data(GAMBLR.data::grch37_ashm_regions,
       fix_names="concat",
       concat_cols=c("gene","region"),
       sep="-"
      )
    } else {
      stop(
        "Please specify one of grch37 or hg38 projections"
      )
    }
  }
  if(verbose){
    print(regions_bed)
  }

  
  
  if(nrow(regions_bed) == 0){
    stop("Region to display doesn't have coordinates in supplied bed or GAMBLR.data::somatic_hypermutation_locations_{genome_build}_v_latest."
    )
  }
  regions_bed = mutate(regions_bed,
    region_name=paste0(chrom,":",start,"-",end))
  regions = pull(regions_bed, region_name)
  names = pull(regions_bed, name)
  
  
  if(missing(maf_data)){
    print("get_ssm_by_regions")
    region_mafs = get_ssm_by_regions(
      regions_bed=regions_bed,
      these_samples_metadata = metadata,
      projection = projection,
      use_name_column = T
    ) 
    #%>% strip_genomic_classes() 
  }else{
    region_mafs = get_ssm_by_regions(maf_data=maf_data,
      regions_bed=regions_bed,
      these_samples_metadata = metadata,
      projection = projection,
      use_name_column = T,
      this_seq_type = this_seq_type
    ) %>% strip_genomic_classes()
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
  
  
  muts_anno = dplyr::left_join(region_mafs, meta_arranged)
  muts_first = dplyr::select(muts_anno, start, region_name) %>%
    group_by(region_name) %>%
    arrange(start) %>%
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
