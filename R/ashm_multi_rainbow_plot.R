#' @title ASHM Multi-panel Rainbow Plot
#'
#' @description Generates a colourful multi-panel overview of hypermutation in regions of interest across many samples.
#'
#' @details The input for this function is a bed-file with the following columns; chr, start, end, name.
#' Note that for this function to work, the column names must be exactly this.
#' The user also needs to specify a vector of names (`regions_to_display`) to further control what regions are to be displayed on the returned plot.
#' It is also possible to exclude specific classifications from the metadata file. This is achieved with `exclude_classifications`.
#' In addition the user can also use the `metadata` parameter to use an already subset and arranged metadata table.
#' This function will call [GAMBLR::get_ssm_by_region] if `maf_data` is not called. For more info, refer to the parameter descriptions of this function.
#'
#' @param regions_bed Bed file with chromosome coordinates, should contain columns chr, start, end, name (with these exact names).
#' @param regions_to_display Optional vector of names from default regions_bed to use.
#' @param exclude_classifications Optional argument for excluding specific classifications from a metadeta file.
#' @param metadata A metadata file already subsetted and arranged on the order you want the samples vertically displayed.
#' @param seq_type the seqtype you want results back for if `maf_data` is not provided.
#' @param custom_colours Provide named vector (or named list of vectors) containing custom annotation colours if you do not want to use standartized pallette.
#' @param classification_column Optional. Override default column for assigning the labels used for colouring in the figure.
#' @param maf_data An already loaded maf, if no provided, this function will call `get_ssm_by_region`, using the regions supplied into `regions_bed`.
#' @param verbose Set to FALSE to rpevent printing the full regions bed file to the console. Default is TRUE.
#'
#' @return Nothing
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' #get lymphgen colours
#' lymphgen_colours = get_gambl_colours(classification = "lymphgen")
#'
#' #build plot
#' ashm_multi_rainbow_plot(regions_to_display = c("BCL2-TSS",
#'                                                "MYC-TSS",
#'                                                "SGK1-TSS",
#'                                                "IGL"),
#'                         custom_colours = lymphgen_colours,
#'                         seq_type = "genome")
#'
ashm_multi_rainbow_plot = function(regions_bed,
                                   regions_to_display,
                                   exclude_classifications,
                                   metadata,
                                   seq_type,
                                   custom_colours,
                                   classification_column = "lymphgen",
                                   maf_data,
                                   verbose = TRUE){

  table_name = check_config_value(config::get("results_tables")$ssm)
  db = check_config_value(config::get("database_name"))
  #get the mutations for each region and combine
  #regions_bed should contain chr, start, end, name (with these exact names)
  if(missing(metadata)){
    metadata = get_gambl_metadata()
    meta_arranged = arrange(metadata, pathology_rank, lymphgen)
  }else{
    meta_arranged = metadata #assume the user already arranged it the way they wanted
  }
  if(!missing(exclude_classifications)){
    meta_arranged = dplyr::filter(meta_arranged, !get(classification_column) %in% exclude_classifications)
  }
  if(missing(regions_bed)){
    regions_bed = GAMBLR.data::somatic_hypermutation_locations_GRCh37_v_latest
    regions_bed = mutate(regions_bed, regions = paste0(chr_name, ":", hg19_start, "-", hg19_end))
    regions_bed = mutate(regions_bed, name = paste0(gene, "-", region))
  }else{
    regions_bed = mutate(regions_bed,regions=paste0(chr,":",start,"-",end))
    #if name column is missing, add it
    if(!"name" %in% colnames(regions_bed))
    {
      regions_bed$name = regions_bed$regions
    }
  }

  if(verbose){
    print(regions_bed)
  }

  names = pull(regions_bed, name)
  names = c(names, "NFKBIZ-UTR", "MAF", "PAX5", "WHSC1", "CCND1",
                   "FOXP1-TSS1", "FOXP1-TSS2", "FOXP1-TSS3", "FOXP1-TSS4",
                   "FOXP1-TSS5", "BCL6", "IGH", "IGL", "IGK", "PVT1", "BCL2") #add some additional regions of interest
  regions = pull(regions_bed, regions)
  regions = c(regions,"chr3:101578214-101578365", "chr16:79627745-79634622", "chr9:36898851-37448583", "chr4:1867076-1977887", "chr11:69451233-69460334", "chr3:71623481-71641671",
                      "chr3:71532613-71559445", "chr3:71343345-71363145", "chr3:71167050-71193679", "chr3:71105715-71118362", "chr3:187406804-188522799","chr14:106144562-106344765",
                      "chr22:23217074-23250428","chr2:89073691-89320640", "chr8:128774985-128876311","chr18:60982124-60990180")
  regions_bed = data.frame(regions = regions, names = names)
  regions_bed = dplyr::filter(regions_bed, names %in% regions_to_display)
  regions = pull(regions_bed, regions)
  names = pull(regions_bed, names)
  if(missing(maf_data)){
    region_mafs = lapply(regions, function(x){get_ssm_by_region(region = x, streamlined = TRUE,seq_type = seq_type)})
  }else{
    region_mafs = lapply(regions, function(x){get_ssm_by_region(region = x, streamlined = TRUE, maf_data = maf_data)})
  }
  tibbled_data = tibble(region_mafs, region_name = names)
  unnested_df = tibbled_data %>%
    unnest_longer(region_mafs)

  unlisted_df = mutate(unnested_df, start = region_mafs$Start_Position, sample_id = region_mafs$Tumor_Sample_Barcode) %>%
    dplyr::select(start,sample_id, region_name)

  meta_arranged = meta_arranged %>% mutate_if(is.factor, as.character)
  meta_arranged = meta_arranged %>% mutate(classification = factor(!!sym(classification_column)))


  muts_anno = left_join(unlisted_df, meta_arranged)
  muts_first = dplyr::select(muts_anno, start, region_name) %>%
    group_by(region_name) %>%
    arrange(start) %>%
    dplyr::filter(row_number() == 1)

  eg = expand_grid(start = pull(muts_first, start), sample_id = pull(meta_arranged, sample_id))
  eg = left_join(eg, muts_first)

  #concatenate expanded frame of points with original mutation data
  real_and_fake = bind_rows(unlisted_df, eg)
  muts_anno = left_join(real_and_fake, meta_arranged)

  muts_anno$sample_id = factor(muts_anno$sample_id, levels = meta_arranged$sample_id)

  if(!missing(regions_to_display)){
    muts_anno = dplyr::filter(muts_anno, region_name %in% regions_to_display)
  }
  #make the plot
  p = muts_anno %>%
        ggplot() +
        geom_point(aes(x = start, y = sample_id, colour = classification), alpha = 0.4, size = 0.6) +
        labs(title = "", subtitle = "", x = "", y = "Sample") +
        theme_Morons() +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = ggplot2::margin(1,1,1,1, "cm"), title = element_blank(), plot.subtitle = element_blank(), axis.title.x = element_blank()) +
        facet_wrap(~region_name, scales = "free_x") +
        guides(color = guide_legend(reverse = TRUE,
                                    override.aes = list(size = 3),
                                    title={{ classification_column }}))

  if(! missing(custom_colours)){
    # ensure only relevant color keys are present
    custom_colours = custom_colours[intersect(names(custom_colours), pull(meta_arranged[,classification_column]))]
    p = p +
        scale_colour_manual(values = custom_colours)
  }

  print(p)

}
