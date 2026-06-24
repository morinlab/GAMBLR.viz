#' @title Pretty Lollipop Plot.
#'
#' @description Generates a visually appealing lollipop plot.
#'
#' @details Retrieve maf data of a specific sample or a set of samples. A gene of interest
#' can then be visualized with the given maf data. Silent mutations can be visualized setting
#' include_silent to TRUE.
#'
#' @param maf_df A data frame containing the mutation data. If the `RefSeq` and `Protein_position` columns are absent, a warning is issued and the plot is produced without protein domain annotation.
#' @param these_samples_metadata A data.frame with metadata. Only required when `show_rate = TRUE` (used to derive the mutation rate denominator). Default `NULL`.
#' @param gene The gene symbol to plot.
#' @param plot_title Optional, the title of the plot. Default is gene.
#' @param refseq_id Insert a specific NM_xxx value of interest
#' @param by_allele Set to FALSE to consider all mutations at the same codon as equivalent.
#' When FALSE, and combined with labelPos, the labels will only indicate the amino acid number. Default is TRUE.
#' @param labelPos A vector of amino acid positions that should be labeled on the plot (i.e. c("28", "315")). If unset, label_threshold is used to determine which positions to label.
#' @param labelHGVSp A character vector of values from HGVSp_Short (i.e. c("p.G28D", "p.E315K")) that will be labeled on the plot.
#' @param label_threshold Minimum mutation recurrence for labels to appear on plot. Superseded by either labelPos or labelHGVSp. Default set to 5.
#' @param show_rate Boolean parameter indicating whether to show the somatic mutation rate in the plot title. Default is FALSE.,
#' @param max_count Sets the mutation count for the largest lollipop size. Default is 10.
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is FALSE.
#' @param font_size Customizable font size for the CoLollipop plot somatic mutation rate and comparison labels. Default is 11pt font.
#' @param font Customizable font for all plotting elements. Default is "sans" which uses the system default font. Use `systemfonts::match_fonts(font)` to confirm if your requested font is available. 
#' @param title_size Font size for the plot title. Default is 14.
#' @param x_axis_size Font size for the x-axis labels. Default is 10.
#' @param domain_label_size Font size for the protein domain labels. Default is 4.
#' @param domain_label_orientation Manually override the orientation of the protein domain labels with "horizontal" or "vertical". Default is "auto", inferred from the size of the smallest domain relative to the protein length. 
#' @param aa_label_size Font size for the amino acid labels. Default is 3.
#' @param point_alpha Alpha to apply to the lollipop points. Default is 1 (no transparency).
#' @param point_size_range Vector of length 2 specifiying the size range of points on the lollipops. Default is c(2,8).
#' @param highlight_driver_regions Logical. If `TRUE` and the input `maf_df`
#'   contains a `mutation_alias` column (added by
#'   [GAMBLR.utils::annotate_curated_drivers()]), draws semi-transparent
#'   coloured rectangles spanning the AA coordinate range of each distinct
#'   driver region. Excludes `_other` aliases. Default `FALSE`.
#' @param driver_region_alpha Transparency of the driver region highlight
#'   boxes (0 = fully transparent, 1 = opaque). Default `0.15`.
#' @param driver_region_padding Number of AA positions added to each side of
#'   a region's observed range. Prevents single-point hotspots from collapsing
#'   to a line. Default `5`.
#' @param limit_driver_regions Optional character vector of region names
#'   (gene-prefix-stripped, e.g. `c("KAT", "HAT")`) used to restrict which
#'   driver region boxes are drawn. Regions not listed are silently dropped.
#'   Default `NULL` (show all non-`_other`/non-`_trunc` regions).
#' @param region_df Optional data frame of pre-computed regions (columns
#'   `mutation_alias`, `xmin`, `xmax`, `label`, `x_mid`). When supplied,
#'   overrides internal computation. Intended for use by
#'   [pretty_colollipop_plot()] to enforce identical region boundaries across
#'   both panels. Default `NULL`.
#'
#' @return A list of plot and data objects.
#'
#' @import dplyr ggplot2 ggrepel ggpubr GAMBLR.data
#' @export
#'
#' @examples
#' suppressPackageStartupMessages(library(GAMBLR.open))
#' suppressMessages(
#'   suppressWarnings({
#' 
#' # Get meta data (BL)
#' # Using GAMBLR.open:  
#' metadata <- GAMBLR.open::get_gambl_metadata() %>% 
#'     filter(
#'         seq_type == "genome", 
#'         cohort %in% c("BL_Pediatric", "BL_Adult"),
#'     ) 
#' 
#' maf_df <- GAMBLR.open::get_coding_ssm(
#'         these_samples_metadata = metadata 
#' )
#'
#' # Construct pretty_lollipop_plot.
#' lollipop_result <- pretty_lollipop_plot(
#'     maf_df, 
#'     these_samples_metadata = metadata, 
#'     gene = "DDX3X", 
#'     point_size_range = c(5, 15), 
#'     point_alpha = 0.8,
#'     label_threshold = 3
#'     )
#' lollipop_result$plot
#' }))
#' 
#' \dontrun{
#'  # Or, with GAMBLR.results:
#'  suppressPackageStartupMessages(library(GAMBLR.results))
#'  metadata <- GAMBLR.results::get_gambl_metadata() %>%
#'     dplyr::filter(cohort %in% c("BL_Pediatric", "BL_Adult"))
#' 
#' maf_df <- GAMBLR.results::get_all_coding_ssm(
#'     these_samples_metadata = metadata, 
#'     basic_columns = FALSE
#' )
#' lollipop_result <- pretty_lollipop_plot(
#'     maf_df, 
#'     these_samples_metadata = metadata, 
#'     gene = "DDX3X", 
#'     point_size_range = c(5, 15), 
#'     point_alpha = 0.8,
#'     label_threshold = 3
#'     )
#' lolipop_result$plot
#'}
pretty_lollipop_plot <- function(
    maf_df,
    these_samples_metadata = NULL,
    gene,
    plot_title,
    refseq_id = NULL,
    by_allele = TRUE,
    max_count = 10,
    include_silent = FALSE,
    label_threshold = 5,
    labelPos = NULL, 
    labelHGVSp = NULL,
    mirror = FALSE, 
    font_size = 11,
    font = "sans",
    show_rate = FALSE,
    title_size = 14,
    x_axis_size = 10,
    domain_label_size = 4,
    domain_label_orientation = "auto", 
    aa_label_size = 3, 
    point_alpha = 1,
    point_size_range = c(2, 8),
    highlight_driver_regions = FALSE,
    driver_region_alpha      = 0.15,
    driver_region_padding    = 5,
    limit_driver_regions     = NULL,
    region_df                = NULL,
    coef_limits              = NULL
){
    ##### Input checks #####
    if(missing(gene)){
        stop("Please provide a gene...")
    }

    if(missing(plot_title)){
        plot_title = gene
    }
    
    no_domain_mode <- !"RefSeq" %in% colnames(maf_df) || !"Protein_position" %in% colnames(maf_df)
    if (no_domain_mode) {
        warning("maf_df is missing 'RefSeq' and/or 'Protein_position'; plotting without protein domain information.")
    }
    
    if(!is.null(labelPos) & !is.null(labelHGVSp)){
        stop("Please provide either labelPos or labelHGVSp, not both.")
    }
    
    if(!is.null(labelHGVSp) & !by_allele){
        stop("Please set by_allele = TRUE when using labelHGVSp.
        When by_allele = FALSE, all alleles at the same position are counted together and HGVSp information is lost. ")
    }

    
    ##### Subset the maf to gene and variants of interest #####
    
    if(include_silent){
        variants <- coding_class
    } else {
        variants <- coding_class[!coding_class %in% c(
            "Silent",
            "Splice_Region" # Is this useful? Aren't these all intronic? Won't affect cds
        )]
    }
    
    maf_df <- as.data.frame(maf_df) %>%
        filter(Hugo_Symbol == gene, Variant_Classification %in% variants)

    if (!no_domain_mode) {
        maf_df <- maf_df %>%
            # Separate RefSeq column if it is comma-delimited
            separate_longer_delim(RefSeq, delim = ",") %>%
            mutate(
                RefSeq = str_remove(RefSeq, "[.].*"),
                aa.length = str_remove(Protein_position, ".*/")
            )

        # Verify that all variants belong to a single gene model
        if (length(unique(maf_df$aa.length)) > 1) {
            maf_df %>%
                count(Hugo_Symbol, Transcript_ID, RefSeq, aa.length) %>%
                print()
            warning("The provided maf_df has more than one gene model for the specified gene.
            Selecting the RefSeqID that appears more frequently. Some variants may be dropped ")
            maf_df <- maf_df %>%
                group_by(RefSeq) %>%
                filter(n() == max(n())) %>%
                ungroup()
        }

        ##### Load the protein domain table for the current gene #####

        if (!is.null(refseq_id)) {
            refseq_id <- str_remove(refseq_id, "[.].*")
            protein_domain_subset <- protein_domains %>%
                filter(HGNC == gene, refseq.ID == refseq_id)
            if (length(unique(protein_domain_subset$refseq.ID)) == 0) {
                warning(paste("The provided refseq_id", refseq_id, "does not match any RefSeq IDs for the specified gene."))
                protein_domain_subset <- data.frame()
            }
        } else {
            protein_domain_subset <- protein_domains %>%
                filter(HGNC == gene) %>%
                filter(refseq.ID %in% unique(maf_df$RefSeq)) %>%
                filter(aa.length == unique(maf_df$aa.length))
        }

        if (length(unique(protein_domain_subset$refseq.ID)) > 1) {
            print(protein_domain_subset)
            protein_domain_subset <- protein_domain_subset %>%
                filter(refseq.ID == unique(protein_domain_subset$refseq.ID[1]))
            message(paste("There is more than one RefSeq model matching the maf_df for the specified gene.
        Arbitrarily selecting", unique(protein_domain_subset$refseq.ID), "to plot. Some variants may be dropped. "))
        } else if (length(unique(protein_domain_subset$refseq.ID)) == 0) {
            warning("None of the protein models matches the provided maf_df. Check the Protein_position and RefSeq columns. ")
            protein_domain_subset <- data.frame(
                HGNC = gene,
                Start = -1,
                End = -1,
                Label = "",
                refseq.ID = NA_character_,
                aa.length = as.numeric(unique(maf_df$aa.length))
            )
        }
    }

    ##### Count mutations according to user specified options #####
    gene_df <- maf_df

    if (!no_domain_mode && !is.na(protein_domain_subset$refseq.ID[1])) {
        gene_df <- gene_df %>% filter(RefSeq == unique(protein_domain_subset$refseq.ID))
    }

    gene_df <- gene_df %>%
        mutate(
            AA = as.numeric(
                gsub(
                    "[^0-9]+",
                    "",
                    gsub("([0-9]+).*", "\\1", HGVSp_Short)
                )
            )
        ) %>%
        mutate(AA_position = gsub("^p\\.", "", HGVSp_Short)) %>%
        arrange(AA)
        
    if(by_allele){
        # Keep different HGVSp values even if AA position is the same
        gene_counts <- gene_df %>%
            group_by(
                HGVSp_Short,
                AA,
                AA_position,
                Variant_Classification
            ) %>%
            arrange(AA) %>%
            summarise(mutation_count = n()) %>%
            ungroup() %>% 
            mutate(size = ifelse(
                mutation_count > max_count, 
                max_count, 
                mutation_count
                ))     
    }else{
        # Collapse HGVSp values if AA is the same
        gene_counts <- gene_df %>%
            group_by(
                AA, 
                Variant_Classification
            ) %>%
            arrange(AA) %>%
            summarise(mutation_count = n()) %>% 
            mutate(size = ifelse(
                mutation_count > max_count, 
                max_count, 
                mutation_count
                )) 
    }
    
    if(!is.null(labelHGVSp)){
        # Only label user-specified HGVSp values
        gene_counts <- gene_counts %>% 
            mutate(label = ifelse(
                HGVSp_Short %in% labelHGVSp, 
                AA_position, 
                ""
            ), 
            size = ifelse(mutation_count > max_count, max_count, mutation_count)) 

    } else if(!is.null(labelPos)){
        # Only label user-specified AA positions
        gene_counts <- gene_counts %>% 
            mutate(label = ifelse(
                AA %in% labelPos, 
                as.character(AA), 
                ""
            ), 
            size = ifelse(mutation_count > max_count, max_count, mutation_count)) 

    } else {
        # Only label variants that exceed label_threshold recurrence
        gene_counts <- gene_counts %>% 
            mutate(label = ifelse(
                mutation_count >= label_threshold, 
                str_remove(HGVSp_Short, "p[.]"), 
                ""
            ), 
            size = ifelse(mutation_count > max_count, max_count, mutation_count)) 

    }

    ##### Generate the final domain_data object and make the domain_plot #####
    if (!no_domain_mode) {
        domain_data <- protein_domain_subset %>%
            data.frame(
                start.points = protein_domain_subset$Start,
                end.points = protein_domain_subset$End,
                text.label = protein_domain_subset$Label,
                color = protein_domain_subset$Label
            )
        domain_data$text.position <- (domain_data$start.points + domain_data$end.points) / 2

        domain_plot <- draw_domain_plot(
            domain_data,
            font = font,
            domain_label_size = domain_label_size,
            domain_label_orientation = domain_label_orientation,
            x_axis_size = x_axis_size
        )
        protein_length <- domain_data$aa.length[1]
    } else {
        domain_data <- NULL
        domain_plot <- NULL
        protein_length <- max(gene_counts$AA, na.rm = TRUE)
    }

    ##### Resolve driver region highlight data #####
    final_region_df <- NULL
    if (!is.null(region_df)) {
        final_region_df <- region_df
    } else if (highlight_driver_regions) {
        if ("mutation_alias" %in% colnames(gene_df)) {
            final_region_df <- compute_driver_regions(
                gene_df, gene, driver_region_padding
            )
            if (nrow(final_region_df) == 0) final_region_df <- NULL
        } else {
            warning(
                "highlight_driver_regions = TRUE but 'mutation_alias' column ",
                "not found in maf_df. Run annotate_curated_drivers() first."
            )
        }
    }

    if (!is.null(final_region_df) && !is.null(limit_driver_regions)) {
        final_region_df <- dplyr::filter(
            final_region_df, label %in% limit_driver_regions
        )
        if (nrow(final_region_df) == 0) final_region_df <- NULL
    }

    ##### Generate the final mutation_plot object #####
    mutation_plot <- draw_mutation_plot(
        gene_counts          = gene_counts,
        plot_title           = plot_title,
        protein_length       = protein_length,
        colours_manual       = get_gambl_colours("mutation"),
        font                 = font,
        mirror               = mirror,
        aa_label_size        = aa_label_size,
        title_size           = title_size,
        max_count            = max_count,
        point_size_range     = point_size_range,
        point_alpha          = point_alpha,
        region_df            = final_region_df,
        region_alpha         = driver_region_alpha,
        coef_limits          = coef_limits
    )

    if (show_rate) {
        if (is.null(these_samples_metadata)) {
            stop("these_samples_metadata is required when show_rate = TRUE.")
        }
        denominator <- length(unique(these_samples_metadata$sample_id))
        if (!no_domain_mode) {
            numerator <- length(unique(gene_df[gene_df$RefSeq == protein_domain_subset$refseq.ID[1], ]$Tumor_Sample_Barcode))
        } else {
            numerator <- length(unique(gene_df$Tumor_Sample_Barcode))
        }
        somatic_mutation_rate <- round((numerator / denominator) * 100, 1)

        mutation_plot <- mutation_plot +
            ggtitle(
                plot_title,
                subtitle = paste0("Somatic mutation rate: ", somatic_mutation_rate, "%")
            ) +
            theme(
                text = element_text(family = font),
                plot.subtitle = element_text(hjust = 0.5)
            )
    }

    ##### Combine domain and mutation plots #####
    if (!no_domain_mode) {
        combined_plot <- ggpubr::ggarrange(
            mutation_plot,
            domain_plot,
            ncol = 1,
            align = "v",
            heights = c(3, 1),
            common.legend = TRUE,
            legend = "right"
        )
    } else {
        combined_plot <- mutation_plot
    }

    to_return <- list(
        plot = combined_plot,
        domain_plot = domain_plot,
        mutation_plot = mutation_plot,
        gene_counts = gene_counts,
        domain_data = domain_data
    )
    return(to_return)
}


draw_domain_plot <- function(
    domain_data, 
    font = "sans", 
    x_axis_size = 10, 
    domain_label_size = 4, 
    domain_label_orientation = "auto"
    ){
    domain_palette <- unname(get_gambl_colours("domains"))[1:length(unique(domain_data$text.label))]
    
    # Determine angle of domain labels based on width of domains relative to protein
    if(domain_label_orientation == "auto"){
        angle = ifelse(
            min((domain_data$End - domain_data$Start)/domain_data$aa.length) >= 0.1,
            0, 90
        )
    } else if(domain_label_orientation == "horizontal"){
        angle = 0
    } else if(domain_label_orientation == "vertical"){
        angle = 90
    } else {
        stop("Invalid domain_label_orientation. Use 'auto', 'horizontal', or 'vertical'.")
    }
    
    
    # Add the p-value to the text label if it exists in the provided domain_data
    if("p_value" %in% colnames(domain_data)){
        domain_data$text.label <- ifelse(
            domain_data$p_value != "", 
            str_c(
                domain_data$text.label,
                "\n",
                domain_data$p_value
            ), 
            domain_data$text.label
        )
    } 
    
    # Fix labels that are too close together
    # TODO: Fix it to retain significant p-values if they exist
    domain_data <- domain_data %>% 
        arrange(text.position) %>% 
        mutate(
            text.label = ifelse(
                text.position - lag(text.position) >= 10 | is.na(lag(text.position)), 
                text.label, 
                ""
            )
        )
    # Specify plot start and end points
    x_min <- 1 # Set 1 as the start for any protein (obviously)
    x_max = unique(domain_data$aa.length)[1] # Assuming aa.length is a single value for the gene
    domain_plot <- ggplot() + 
            geom_rect(
                aes(
                    xmin = x_min, 
                    xmax = x_max, 
                    ymin = -1, 
                    ymax = 1
                    ), 
                fill = "lightgrey", 
                color = "lightgrey"
            ) + 
            geom_rect(
                data = domain_data, 
                aes(
                    xmin = start.points, 
                    xmax = end.points, 
                    ymin = -1.2, 
                    ymax = 1.2, 
                    fill = color
                    ), 
                color = "black", 
                show.legend = FALSE
            ) +
            geom_text(
                data = domain_data, 
                aes(
                    x = text.position, 
                    y = 0,
                    angle = angle, 
                    label = text.label
                ), 
                hjust = 0.5, 
                vjust = 0.5,
                size = domain_label_size, 
                family = font
            ) +
            scale_fill_manual(values = domain_palette) + 
            scale_x_continuous(
                breaks = seq(
                    0,
                    x_max,
                    by = 5*round(x_max/25) # Breaks in multiples of 5
                ), 
                limits = c(0, x_max)
            ) +
            theme_void() +
            theme(
                axis.text.x = element_text(
                    size = x_axis_size, 
                    vjust = rel(0.5)
                    ), 
                text = element_text(family = font)
            ) 
    return(domain_plot)
}

draw_mutation_plot <- function(
    gene_counts,
    plot_title,
    protein_length,
    font = "sans",
    colours_manual,
    mirror = FALSE,
    max_count,
    point_size_range,
    title_size,
    aa_label_size,
    point_alpha,
    region_df    = NULL,
    region_alpha = 0.15,
    coef_limits  = NULL
    ){

    nudge_factor <- ifelse(mirror, -0.5, 0.5)

    mutation_plot <- ggplot()

    if (!is.null(region_df) && nrow(region_df) > 0) {
        has_coef <- "coef" %in% colnames(region_df)
        if (has_coef) {
            coef_limits_vec <- if (!is.null(coef_limits)) {
                if (length(coef_limits) == 2L) coef_limits else c(-coef_limits, coef_limits)
            } else {
                lim <- max(abs(region_df$coef), na.rm = TRUE)
                if (lim == 0) lim <- 1
                c(-lim, lim)
            }
            mutation_plot <- mutation_plot +
                geom_rect(
                    data        = region_df,
                    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                        fill = coef),
                    alpha       = region_alpha,
                    inherit.aes = FALSE
                ) +
                scale_fill_gradient2(
                    low      = "#2166AC",
                    mid      = "white",
                    high     = "#D6604D",
                    midpoint = 0,
                    limits   = coef_limits_vec,
                    oob      = scales::squish,
                    name     = "Coef"
                ) +
                geom_text(
                    data        = if ("show_label" %in% colnames(region_df))
                                      dplyr::filter(region_df, show_label)
                                  else region_df,
                    aes(x = x_mid,
                        y     = if (mirror) -Inf else Inf,
                        label = label),
                    vjust       = if (mirror) -0.5 else 1.5,
                    size        = aa_label_size,
                    family      = font,
                    inherit.aes = FALSE
                )
        } else {
            n_reg   <- nrow(region_df)
            reg_pal <- stats::setNames(
                grDevices::hcl.colors(n_reg, palette = "Dark 3"),
                region_df$label
            )
            mutation_plot <- mutation_plot +
                geom_rect(
                    data        = region_df,
                    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                        fill = label),
                    alpha       = region_alpha,
                    inherit.aes = FALSE
                ) +
                scale_fill_manual(name = "Driver\nRegion", values = reg_pal)
        }
    }

    mutation_plot <- mutation_plot +
        geom_segment(
            data = gene_counts, 
            aes(
                x = AA, 
                xend = AA, 
                y = 0, 
                yend = mutation_count
                )
        ) +
        geom_point(
            data = gene_counts, 
            aes(
                x = AA, 
                y = mutation_count, 
                color = Variant_Classification, 
                size = size
                ), 
            alpha = point_alpha
        ) +
        ggrepel::geom_text_repel(data = gene_counts, aes(
                x = AA, 
                y = mutation_count, 
                label = label, 
                point.size = size
            ), 
            min.segment.length = 0, 
            nudge_y = nudge_factor,
            nudge_x = 0,
            direction = "y", 
            size = aa_label_size, 
            family = font
        ) +
        labs(
            x = "AA Position", 
            y = "Mutation Count", 
            title = paste0(
                plot_title
            )
        ) + 
        scale_y_continuous(
            breaks = seq(
                0,
                max(gene_counts$mutation_count) * 1.2,
                by = ifelse(max(gene_counts$mutation_count) * 1.2 <= 5, 1, round(max(gene_counts$mutation_count) * 1.2 / 5))
            ), 
            limits = c(0, max(gene_counts$mutation_count) * 1.2)
        )  +
        xlim(1, protein_length) +
        scale_size_continuous(
            limits = c(1, max_count),
            breaks = seq(
                1, 
                max_count,
                by = round(max_count / 5)
            ), 
            range = point_size_range,
            guide = "none"  
        )  +
        ggpubr::theme_pubr() +
        theme(
            plot.title = element_text(
                hjust = 0.5, 
                size = title_size
                ),
            # axis.text.x = element_blank(),
            # axis.line.x = element_blank(),
            # axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            text = element_text(family = font), 
            legend.position = "right"
        ) +
        scale_color_manual(
            name = "Variant\nClassification", 
            values = colours_manual
        )+
    guides(color = guide_legend(override.aes = list(size = 3)))
    return(mutation_plot)
}

compute_driver_regions <- function(maf_with_aa, gene, padding = 5) {
    maf_with_aa %>%
        dplyr::filter(
            !grepl("_other$|_trunc$", mutation_alias),
            !is.na(mutation_alias),
            !is.na(AA)
        ) %>%
        dplyr::group_by(mutation_alias) %>%
        dplyr::summarise(
            xmin = min(AA, na.rm = TRUE) - padding,
            xmax = max(AA, na.rm = TRUE) + padding,
            .groups = "drop"
        ) %>%
        dplyr::mutate(
            xmin  = pmax(xmin, 1L),
            label = sub(paste0("^", gene, "_"), "", mutation_alias),
            x_mid = (xmin + xmax) / 2
        )
}
