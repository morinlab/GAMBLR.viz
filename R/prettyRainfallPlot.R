#' @title Rainfall Plot
#'
#' @description Plot a rainfall plot for one sample. This function takes in the MAF data frame, or path to a custom MAF file.
#'
#' @details Create a sample-level rainfall plot visualizing single nucleotide substitutions mutations for selected chromosomes.
#'
#' @param this_sample_id Sample id for the sample to display. This is argument is not required if you want a multi-sample plot but is otherwise needed.
#' @param label_ashm_genes Boolean argument indicating whether the aSHM regions will be labeled or not.
#' @param projection Specify projection (grch37 or hg38) of mutations. Default is grch37.
#' @param chromosome Provide one or more chromosomes to plot. The chr prefix can be inconsistent with projection and will be handled.
#' @param this_maf Specify custom MAF data frame of mutations.
#' @param maf_path Specify path to MAF file if it is not already loaded into data frame.
#' @param zoom_in_region Provide a specific region in the format "chromosome:start-end" to zoom in to a specific region.
#' @param label_sv Boolean argument to specify whether label SVs or not. Only supported if a specific chromosome or zoom in region are specified.
#' @param seq_type Specify one of "genome" or "capture" when relying on the function to obtain mutations from a region (i.e. if you haven't provided a MAF or single sample_id)
#'
#' @return a ggplot2 plot. Print it using print() or save it using ggsave()
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import ggplot2 dplyr readr stringr tidyr
#' @export
#'
#' @examples
#' prettyRainfallPlot(this_sample_id = "DOHH-2",
#'                    seq_type = "genome",
#'                    zoom_in_region = "8:125252796-135253201",
#'                    label_sv = TRUE)
#'
prettyRainfallPlot = function(this_sample_id,
                              label_ashm_genes = TRUE,
                              projection = "grch37",
                              chromosome,
                              this_maf,
                              maf_path,
                              zoom_in_region,
                              seq_type,
                              label_sv = FALSE) {
  if (missing(this_sample_id)) {
    warning("No sample_id was provided. Using all mutations in the MAF within your region!")
    if(missing(zoom_in_region)){
      stop("Must provide a zoom_in_region to plot when showing data from more than one patient")
    }
  }

  # allow user to specify chromosome prefix inconsistent with chromosome names
  if (!missing(chromosome)) {
    chromosome = GAMBLR.helpers::standardize_chr_prefix(incoming_vector = chromosome, projection = projection)
  }

  # allow to zoom in to a specific region
  if (!missing(zoom_in_region)) {
    region = zoom_in_region
    zoom_in_region = GAMBLR.data::region_to_chunks(zoom_in_region)
    zoom_in_region$chromosome = GAMBLR.helpers::standardize_chr_prefix(incoming_vector = zoom_in_region$chromosome,
                                                                       projection = projection)
    zoom_in_region$start = as.numeric(zoom_in_region$start)
    zoom_in_region$end = as.numeric(zoom_in_region$end)
  }

  if (label_ashm_genes) {
    if (projection == "grch37") {
      ashm_regions = GAMBLR.data::grch37_ashm_regions %>%
        dplyr::rename("start" = "hg19_start",
                      "end" = "hg19_end",
                      "Chromosome" = "chr_name") %>%
        dplyr::mutate(Chromosome = str_remove(Chromosome, pattern = "chr"))
    } else if (projection == "hg38") {
      ashm_regions = GAMBLR.data::hg38_ashm_regions %>%
        rename("start" = "hg38_start",
               "end" = "hg38_end",
               "Chromosome" = "chr_name")
    } else {
      stop("Please specify one of grch37 or hg38 projections")
    }
    if (!missing(chromosome)) {
      ashm_regions = dplyr::filter(ashm_regions, Chromosome %in% chromosome)
    }
    if (!missing(zoom_in_region)) {
      ashm_regions = dplyr::filter(
        ashm_regions,
        (
          Chromosome %in% zoom_in_region$chromosome &
            start >= zoom_in_region$start &
            end <= zoom_in_region$end
        )
      )
    }
    ashm_regions = ashm_regions %>%
      group_by(gene) %>%
      slice_head() %>%
      ungroup()

    # this will be needed for consistent labeling with rainfall plots
    ashm_regions = ashm_regions %>%
      arrange(match(
        Chromosome,
        str_sort(ashm_regions$Chromosome, numeric = TRUE)
      ))
    ashm_regions = ashm_regions %>%
      mutate(Chromosome_f = factor(Chromosome, levels = unique(ashm_regions$Chromosome)))
  }

  # if user is subsetting by chromosome or zooming in to a specific region, it is possible there are no aSHM features to show
  # handle this case separately
  if (nrow(ashm_regions) == 0) {
    message(
      "Warning: after subsetting to a regions you requested to plot, there are no aSHM features to overlap on the final graph."
    )
    label_ashm_genes = FALSE
  }

  # get ssm for the requested sample
  if (!missing(this_maf)) {
    if(missing(this_sample_id)){
      these_ssm=this_maf
      this_sample_id = "all samples"
    }else{
      message ("Using the suppplied MAF df to obrain ser of SSM for the specified sample ...")
      these_ssm = this_maf %>%
        dplyr::filter(Tumor_Sample_Barcode %in% this_sample_id)
    }
  } else if (!missing (maf_path)) {
    message ("Path to custom MAF file was provided, reading SSM using the custom path ...")

    this_maf = suppressMessages(read_tsv(maf_path))
    if(!missing(this_sample_id)){
      this_maf = this_maf %>% dplyr::filter(Tumor_Sample_Barcode %in% this_sample_id)
    }else{
      this_sample_id = "all samples"
    }
  } else if(!missing(this_sample_id)) {
    message ("MAF df or path to custom MAF file was not provided, getting SSM using GAMBLR ...")
    these_ssm = get_ssm_by_sample(this_sample_id,
                                  projection = projection, 
                                  this_seq_type = seq_type)
  }else if(!missing(seq_type)){
    if(missing(this_sample_id)){
      this_sample_id = "all samples"
    }
    message(paste("Will use all mutations for",seq_type, "in this region:",zoom_in_region))
    these_ssm = get_ssm_by_region(region = region,seq_type = seq_type,projection=projection)
  }

  # do rainfall calculation using lag
  rainfall_points = dplyr::select(
    these_ssm,
    Hugo_Symbol,
    Chromosome,
    Start_Position,
    End_Position,
    Reference_Allele,
    Tumor_Seq_Allele2
  ) %>%
    arrange(Chromosome, Start_Position) %>%
    group_by(Chromosome) %>%  # group by chromosome to calculate lag per chromosome
    dplyr::mutate(
      IMD = Start_Position - dplyr::lag(Start_Position),
      # used for coloring
      # all indels are squished to the same color
      Substitution = ifelse((
        Reference_Allele %in% c("A", "T",  "C", "G") &
          Tumor_Seq_Allele2 %in% c("A", "T",  "C", "G")
      ),
      paste(Reference_Allele, Tumor_Seq_Allele2, sep = '>'),
      "InDel"
      )
    ) %>%
    dplyr::mutate(IMD = log(IMD)) %>%
    ungroup() %>%
    drop_na(IMD) # for the first point of each chromosome, NAs are produced generating a warning message

  # collapse substitutions into classes
  rainfall_points$Substitution = rainfall_conv[as.character(rainfall_points$Substitution)]

  # ensure order of grids in the plot is sorted
  rainfall_points = rainfall_points %>%
    arrange(match(
      Chromosome,
      str_sort(rainfall_points$Chromosome, numeric = TRUE)
    ))
  rainfall_points = rainfall_points %>%
    mutate(Chromosome_f = factor(Chromosome, levels = unique(rainfall_points$Chromosome)))
  if (!missing(chromosome)) {
    rainfall_points = dplyr::filter(rainfall_points, Chromosome %in% chromosome)
  }
  if (!missing(zoom_in_region)) {
    rainfall_points = dplyr::filter(
      rainfall_points,
      (
        Chromosome %in% zoom_in_region$chromosome &
          Start_Position >= zoom_in_region$start &
          End_Position <= zoom_in_region$end
      )
    )
  }

  # if user is subsetting by chromosome or zooming in to a specific region, are there any SSM left to plot?
  if (nrow(rainfall_points) == 0) {
    stop("After subsetting to a regions you requested to plot, there are no SSM to display.")
  }

  # label SVs if user wants to overlap this data
  if (!missing(chromosome) & label_sv) {
    sv_chromosome = chromosome
  } else if (!missing(zoom_in_region) & label_sv) {
    sv_chromosome = zoom_in_region$chromosome
  } else if (label_sv) {
    stop(
      "Labeling SV is only supported when a particular chromosome or zoomed region is plotted."
    )
  }

  if (label_sv) {
    message("Getting combined manta + GRIDSS SVs using GAMBLR ...")
    these_sv = get_manta_sv(these_sample_ids  = this_sample_id)
    if ("SCORE" %in% colnames(these_sv)) {
      these_sv = these_sv %>%
        rename("SOMATIC_SCORE" = "SCORE")
    }
    # annotate SV
    these_sv = annotate_sv(these_sv)

    # make SVs a long df with 1 record per SV corresponding to the strand
    sv_to_label =
      melt(
        these_sv %>% select(
          chrom1,
          start1,
          end1,
          chrom2,
          start2,
          end2,
          tumour_sample_id,
          gene,
          partner,
          fusion
        ),
        id.vars = c(
          "tumour_sample_id",
          "gene",
          "partner",
          "fusion",
          "start1",
          "end1",
          "start2",
          "end2"
        ),
        variable.name = "chromosomeN",
        value.name = "Chromosome"
      ) %>%
      dplyr::filter(Chromosome %in% sv_chromosome)

    # are there any SVs on this chromosome/region?
    if (nrow(sv_to_label) > 0) {
      sv_to_label =
        sv_to_label %>%
        melt(
          .,
          id.vars = c(
            "tumour_sample_id",
            "gene",
            "partner",
            "fusion",
            "chromosomeN",
            "Chromosome"
          )
        ) %>%
        group_by(fusion, chromosomeN) %>%
        dplyr::filter(if (grepl("1", chromosomeN))
          variable %in% c("start1", "end1")
          else
            variable %in% c("start2", "end2")) %>%
        dplyr::mutate(variable = gsub("1|2", "", variable)) %>%
        distinct(fusion, Chromosome, variable, .keep_all = TRUE) %>%
        spread(., variable, value) %>%
        dplyr::rename("End_Position" = "end",
                      "Start_Position" = "start") %>%
        ungroup
    } else {
      message(
        "Warning: after subsetting to a regions you requested to plot, there are no SV features to overlap on the final graph."
      )
      label_sv = FALSE
    }

    # when we are plotting region and not whole chromosome, ensure SV is within that region
    if (!missing(zoom_in_region) & label_sv) {
      sv_to_label = dplyr::filter(
        sv_to_label,
        (
          Start_Position >= zoom_in_region$start &
            End_Position <= zoom_in_region$end
        )
      )
      # When we did filtering to start/end for a region, are there any SV to plot?
      if (nrow(sv_to_label) == 0) {
        message(
          "Warning: after subsetting to a regions you requested to plot, there are no SV features to overlap on the final graph."
        )
        label_sv = FALSE
      }
    }

    sv_to_label = sv_to_label %>%
      mutate(Chromosome_f = factor(Chromosome))
  }

  p = ggplot(rainfall_points) +
    geom_point(aes(x = Start_Position, y = IMD, color = Substitution)) +
    scale_color_manual(values = GAMBLR.helpers::get_gambl_colours("rainfall")) +
    ylab("log(IMD)") +
    GAMBLR.helpers::theme_Morons() +
    facet_wrap( ~ Chromosome_f, scales = "free_x") +
    ggtitle(this_sample_id) +
    theme(plot.title = element_text(hjust = 0)) # left-align title plot

  if (label_ashm_genes) {
    p = p +
      ggrepel::geom_text_repel(
        data = ashm_regions,
        aes(start, 1, label = gene),
        size = 4,
        segment.size = 0.5,
        segment.color = "#000000",
        force_pull = 0,
        arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
        max.overlaps = 15,
        segment.curvature = 0.25,
        segment.ncp = 4,
        segment.angle = 25
      )
  }

  if (label_sv) {
    p = p +
      geom_vline(
        data = sv_to_label,
        aes(xintercept = Start_Position),
        color = "lightgreen",
        alpha = .7
      ) +
      geom_text(data = sv_to_label,
                aes(End_Position, 15, label = fusion, color = "lightgreen"))
  }

  # show x-axis coordinates if zooming in to a specific region, but not if looking chromosome/genome-wide
  if (missing(zoom_in_region)) {
    p = p + guides(x = "none")
  }

  return(p)
}
