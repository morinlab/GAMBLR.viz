#' @title Rainfall Plot
#'
#' @description Plot a rainfall plot for one sample. This function takes in the
#'      MAF data frame, or path to a custom MAF file.
#'
#' @details Create a sample-level rainfall plot visualizing single nucleotide
#'      substitutions mutations for selected chromosomes.
#'
#' @param this_sample_id Sample id for the sample to display. This is argument
#'      is not required if you want a multi-sample plot but is otherwise needed.
#' @param label_ashm_genes Boolean argument indicating whether the aSHM regions
#'      will be labeled or not.
#' @param projection Specify projection (grch37 or hg38) of mutations. Default
#'      is grch37.
#' @param chromosome Provide one or more chromosomes to plot. The chr prefix can
#'      be inconsistent with projection and will be handled.
#' @param this_maf Specify custom MAF data frame of mutations.
#' @param maf_path Specify path to MAF file if it is not already loaded into
#'      data frame.
#' @param zoom_in_region Provide a specific region in the format
#'      "chromosome:start-end" to zoom in to a specific region.
#' @param label_sv Boolean argument to specify whether label SVs or not with
#'      green line on rainfall plot.
#' @param this_seq_type Specify one of "genome" or "capture" when relying on the
#'      function to obtain mutations from a region (i.e. if you haven't provided
#'      a MAF or single sample_id)
#' @param plot_title Specify the title for the returned plot, default is
#'      "my_plot".
#' @param annotate_sv Optionally annotate intrachromosomal SVs to label the gene
#'      and partner information on the plot. Default is TRUE (perform
#'      annotation).
#'
#' @return a ggplot2 object
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import ggplot2 dplyr readr stringr tidyr GAMBLR.data GAMBLR.helpers
#' @export
#'
#' @examples
#' # Will annotate and label SVs
#' prettyRainfallPlot(
#'      this_sample_id = "DOHH-2",
#'      this_seq_type = "genome",
#'      zoom_in_region = "8:125252796-135253201",
#'      label_sv = TRUE
#' )
#'
#' # Will not annotate SVs (use raw bedpe) but still label them
#' prettyRainfallPlot(
#'      this_sample_id = "DOHH-2",
#'      this_seq_type = "genome",
#'      zoom_in_region = "8:125252796-135253201",
#'      label_sv = TRUE,
#'      annotate_sv = FALSE
#' )
#'
#' # Will use user-specified SV data
#' sv <- get_manta_sv(
#'     these_sample_ids = "DOHH-2"
#' )
#'
#' prettyRainfallPlot(
#'      this_sample_id = "DOHH-2",
#'      this_seq_type = "genome",
#'      zoom_in_region = "8:125252796-135253201",
#'      sv_data = sv,
#'      label_sv = TRUE
#' )
#'
prettyRainfallPlot = function(
    this_sample_id = NULL,
    label_ashm_genes = TRUE,
    projection = "grch37",
    chromosome,
    sv_data = NULL,
    this_maf = NULL,
    maf_path,
    zoom_in_region,
    this_seq_type,
    label_sv = FALSE,
    plot_title = "my_plot",
    annotate_sv = TRUE
){

    if(is.null(this_sample_id)){
        warning(
            "No sample_id was provided. Using all mutations in the MAF within your region!"
        )
        if(
            !is.null(this_maf) &&
            length(unique(this_maf$Tumor_Sample_Barcode)) > 1 &&
            missing(zoom_in_region)
        ){
            stop(
                "Must provide a zoom_in_region to plot when showing data from more than one sample ID"
            )
        }
    }

    # allow user to specify chromosome prefix inconsistent with chromosome names
    if (!missing(chromosome)) {
        chromosome <- GAMBLR.helpers::standardize_chr_prefix(
            incoming_vector = chromosome,
            projection = projection
        )
    }

    # allow to zoom in to a specific region
    if (!missing(zoom_in_region)) {
        region <- zoom_in_region
        zoom_in_region <- GAMBLR.data::region_to_chunks(zoom_in_region)
        zoom_in_region$chromosome <- GAMBLR.helpers::standardize_chr_prefix(
            incoming_vector = zoom_in_region$chromosome,
            projection = projection
        )
        zoom_in_region$start <- as.numeric(zoom_in_region$start)
        zoom_in_region$end <- as.numeric(zoom_in_region$end)
    }

    if (label_ashm_genes) {
        if (projection == "grch37") {
        ashm_regions <- GAMBLR.data::grch37_ashm_regions %>%
            dplyr::rename(
                "start" = "hg19_start",
                "end" = "hg19_end",
                "Chromosome" = "chr_name"
            ) %>%
            dplyr::mutate(
                Chromosome = str_remove(
                    Chromosome,
                    pattern = "chr"
                )
            )
        } else if (projection == "hg38") {
        ashm_regions <- GAMBLR.data::hg38_ashm_regions %>%
            dplyr::rename(
                "start" = "hg38_start",
                "end" = "hg38_end",
                "Chromosome" = "chr_name"
            )
        } else {
            stop(
                "Please specify one of grch37 or hg38 projections"
            )
        }

        if (!missing(chromosome)) {
            ashm_regions <- dplyr::filter(
                ashm_regions,
                Chromosome %in% chromosome
            )
        }

        if (!missing(zoom_in_region)) {
            ashm_regions <- dplyr::filter(
                ashm_regions,
                (
                    Chromosome %in% zoom_in_region$chromosome &
                        start >= zoom_in_region$start &
                        end <= zoom_in_region$end
                )
            )
        }
        ashm_regions <- ashm_regions %>%
            group_by(gene) %>%
            slice_head() %>%
            ungroup()

        # this will be needed for consistent labeling with rainfall plots
        ashm_regions = ashm_regions %>%
            dplyr::arrange(
                match(
                    Chromosome,
                    str_sort(
                        ashm_regions$Chromosome,
                        numeric = TRUE
                    )
                )
            )
        ashm_regions <- ashm_regions %>%
            dplyr::mutate(
                Chromosome_f = factor(
                    Chromosome,
                    levels = unique(
                        ashm_regions$Chromosome
                    )
                )
            )
        # if user is subsetting by chromosome or zooming in to a specific region,
        # it is possible there are no aSHM features to show
        # handle this case separately
        if (nrow(ashm_regions) == 0) {
            message(
                "Warning: after subsetting to a regions you requested to plot, there are no aSHM features to overlap on the final graph."
            )
            label_ashm_genes = FALSE
        }
    }

    # get ssm for the requested sample
    if (!missing(this_maf)) {
        if(missing(this_sample_id)){
            these_ssm <- this_maf
            this_sample_id <- "all samples"
        }else{
            message(
                "Using the suppplied MAF df to obrain ser of SSM for the specified sample ..."
            )

            these_ssm <- this_maf %>%
                dplyr::filter(
                    Tumor_Sample_Barcode %in% this_sample_id
                )
        }
    } else if (!missing (maf_path)) {
        message (
            "Path to custom MAF file was provided, reading SSM using the custom path ..."
        )

        this_maf <- suppressMessages(read_tsv(maf_path))
        if(!missing(this_sample_id)){
            this_maf <- this_maf %>%
                dplyr::filter(
                    Tumor_Sample_Barcode %in% this_sample_id
                )
        }
    } else if(!missing(this_sample_id)) {
        message (
            "MAF df or path to custom MAF file was not provided, getting SSM using GAMBLR ..."
        )

        these_ssm <- get_ssm_by_sample(
            this_sample_id = this_sample_id,
            projection = projection,
            this_seq_type = this_seq_type
        )
    }else if(!missing(this_seq_type)){
        message(
            paste(
                "Will use all mutations for ",
                this_seq_type,
                "in this region: ",
                zoom_in_region
            )
        )

        these_ssm <- get_ssm_by_region(
            region = region,
            this_seq_type = this_seq_type,
            projection = projection
        )
    }

    # do rainfall calculation using lag
    rainfall_points <- dplyr::select(
        these_ssm,
        Hugo_Symbol,
        Chromosome,
        Start_Position,
        End_Position,
        Reference_Allele,
        Tumor_Seq_Allele2
    ) %>%
        dplyr::arrange(Chromosome, Start_Position) %>%
        group_by(Chromosome) %>%  # group by chromosome to calculate lag per chromosome
        dplyr::mutate(
            IMD = Start_Position - dplyr::lag(Start_Position),
            # used for coloring
            # all indels are squished to the same color
            Substitution = ifelse(
                (
                    Reference_Allele %in% c("A", "T",  "C", "G") &
                    Tumor_Seq_Allele2 %in% c("A", "T",  "C", "G")
                ),
                paste(
                    Reference_Allele,
                    Tumor_Seq_Allele2,
                    sep = '>'
                ),
                "InDel"
            )
        ) %>%
        dplyr::mutate(
            IMD = log10(IMD)
        ) %>%
        ungroup() %>%
        # for the first point of each chromosome, NAs are produced
        # generating a warning message
        drop_na(IMD)

    # collapse substitutions into classes
    rainfall_points$Substitution <- rainfall_conv[
        as.character(rainfall_points$Substitution)
    ]

    # ensure order of grids in the plot is sorted
    rainfall_points <- rainfall_points %>%
        dplyr::arrange(
            match(
                Chromosome,
                str_sort(rainfall_points$Chromosome, numeric = TRUE)
            )
        )
    rainfall_points <- rainfall_points %>%
        dplyr::mutate(
            Chromosome_f = factor(
                Chromosome,
                levels = unique(
                    rainfall_points$Chromosome
                )
            )
        )
    if (!missing(chromosome)) {
        rainfall_points <- dplyr::filter(
            rainfall_points,
            Chromosome %in% chromosome
        )
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

    # if user is subsetting by chromosome or zooming in to a specific region,
    # are there any SSM left to plot?
    if (nrow(rainfall_points) == 0) {
        stop(
            "After subsetting to a regions you requested to plot, there are no SSM to display."
        )
    }

    # label SVs if user wants to overlap this data
    if (!missing(chromosome) & label_sv) {
        sv_chromosome <- chromosome
    } else if (!missing(zoom_in_region) & label_sv) {
        sv_chromosome <- zoom_in_region$chromosome
    } else if (label_sv) {
        sv_chromosome <- standardize_chr_prefix(
            incoming_vector = c(1:22, "X"),
            projection = projection
        )
    }

    if (label_sv) {

        if(missing(sv_data)){
            message(
                "Getting combined manta + GRIDSS SVs using GAMBLR ..."
            )
            these_sv <- get_manta_sv(
                these_sample_ids = this_sample_id,
                projection = projection
            )
        } else {
            these_sv <- sv_data
        }

        if ("SCORE" %in% colnames(these_sv)) {
            these_sv <- these_sv %>%
                rename("SOMATIC_SCORE" = "SCORE")
        }
        # annotate SV
        if (annotate_sv){
            these_sv <- annotate_sv(
                these_sv,
                genome_build = projection
            )
        } else {
            these_sv <- these_sv %>%
                dplyr::rename(
                    chrom1 = "CHROM_A",
                    start1 = "START_A",
                    end1 = "END_A",
                    chrom2 = "CHROM_B",
                    start2 = "START_B",
                    end2 = "END_B"
                ) %>%
                dplyr::mutate(
                    gene = 1:nrow(these_sv),
                    partner = letters[gene],
                    fusion = paste(
                        gene,
                        partner,
                        sep = "-"
                    )
                )
        }

        # make SVs a long df with 1 record per SV corresponding to the strand
        sv_to_label <- reshape2::melt(
            these_sv %>%
                select(
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
        dplyr::filter(
            Chromosome %in% sv_chromosome
        )

        # are there any SVs on this chromosome/region?
        if (nrow(sv_to_label) > 0) {
        sv_to_label <- sv_to_label %>%
            reshape2::melt(
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
            dplyr::filter(
                if (grepl("1", chromosomeN))
                    variable %in% c("start1", "end1")
                else
                    variable %in% c("start2", "end2")
            ) %>%
            dplyr::mutate(
                variable = gsub("1|2", "", variable)
            ) %>%
            distinct(fusion, Chromosome, variable, .keep_all = TRUE) %>%
            spread(., variable, value) %>%
            dplyr::rename(
                "End_Position" = "end",
                "Start_Position" = "start"
            ) %>%
            ungroup
        } else {
            message(
                "Warning: after subsetting to a regions you requested to plot, there are no SV features to overlap on the final graph."
            )
            label_sv <- FALSE
        }

        # when we are plotting region and not whole chromosome,
        # ensure SV is within that region
        if (!missing(zoom_in_region) & label_sv) {
            sv_to_label <- dplyr::filter(
                sv_to_label,
                (
                    Start_Position >= zoom_in_region$start &
                    End_Position <= zoom_in_region$end
                )
            )
            # When we did filtering to start/end for a region,
            # are there any SV to plot?
            if (nrow(sv_to_label) == 0) {
                message(
                    "Warning: after subsetting to a regions you requested to plot, there are no SV features to overlap on the final graph."
                )
                label_sv <- FALSE
            }
        }

        sv_to_label <- sv_to_label %>%
            dplyr::mutate(
                Chromosome_f = factor(Chromosome)
            )
    }


    p <- ggplot(
            rainfall_points,
            aes(
                x = Start_Position,
                y = IMD
            )
        ) +
        scale_color_manual(
            values = GAMBLR.helpers::get_gambl_colours("rainfall")
        ) +
        ylab(expression(log[10](IMD))) +
        GAMBLR.helpers::theme_Morons() +
        facet_grid(
            . ~ Chromosome_f,
            scales = "free_x",
            space = "free_x",
            switch = "x"
        ) +
        ggtitle(plot_title) +
        theme(
            plot.title = element_text(hjust = 0),  # left-align title plot
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 16, colour = "black"),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(colour = "black"),
            panel.spacing.x = unit(0.1, "lines"),
            panel.border = element_blank(),
            text = element_text(size = 16, colour = "black", family="sans"),
            strip.background = element_blank(),
            strip.placement = "outside",
            panel.grid = element_blank()
        )

    if (label_ashm_genes) {
        p <- p +
        ggrepel::geom_text_repel(
            data = ashm_regions,
            aes(
                start,
                1,
                label = gene
            ),
            size = 4,
            segment.size = 0.5,
            segment.color = "#000000",
            force_pull = 0,
            arrow = arrow(
                length = unit(0.05, "inches"),
                type = "closed"
            ),
            max.overlaps = 15,
            segment.curvature = 0.25,
            segment.ncp = 4,
            segment.angle = 25
        )
    }

    if (label_sv) {
        p <- p +
        geom_vline(
            data = sv_to_label,
            aes(
                xintercept = Start_Position
            ),
            color = "lightgreen",
            alpha = .7
        )
    }

    if(annotate_sv) {
        max_val <- max(rainfall_points$IMD)
        p <- p +
        geom_text(
            data = sv_to_label,
            aes(
                End_Position,
                max_val+1,
                label = fusion,
                color = "lightgreen"
            ),
            show.legend = FALSE
        )
    }

    p <- p +
    geom_point(
        inherit.aes = TRUE,
        aes(color = Substitution)
    )

    return(p)
}