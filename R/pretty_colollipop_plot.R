#' @title Pretty CoLollipop Plot.
#'
#' @description Generates a ggplot-compatible figure of 2 [GAMBLR.viz::pretty_lollipop_plot] mirrored.
#'
#' @details Retrieve maf data of a specific sample or a set of samples for comparison. A gene of interest 
#' can then be visualized with the given maf data files, and comparison commands. Silent mutations can be 
#' visualized setting include_silent to TRUE. 
#'
#' @param maf_df A data frame containing the mutation data. Required when using metadata-based stratification; omit when supplying `maf_top`/`maf_bottom`.
#' @param these_samples_metadata A data.frame with metadata for the CoLollipop Plot. Required when using metadata-based stratification; omit when supplying `maf_top`/`maf_bottom`.
#' @param comparison_column The name of the metadata column containing the comparison values. Required when using metadata-based stratification; omit when supplying `maf_top`/`maf_bottom`.
#' @param comparison_values Optional: If the comparison column contains more than two values or is not a factor, specify a character vector of length two in the order you would like the factor levels to be set, reference group first. When using `maf_top`/`maf_bottom`, serves only as plot labels and defaults to `c("Group 1", "Group 2")`.
#' @param maf_top A data frame of mutations to use for the upper lollipop panel. When provided together with `maf_bottom`, bypasses metadata-based stratification entirely.
#' @param maf_bottom A data frame of mutations to use for the lower lollipop panel. When provided together with `maf_top`, bypasses metadata-based stratification entirely.
#' @param name_top Display label for the upper panel group when using `maf_top`/`maf_bottom`. Used in the plot title and y-axis label. Ignored when `comparison_values` is also supplied. Default `"Group 1"`.
#' @param name_bottom Display label for the lower panel group when using `maf_top`/`maf_bottom`. Used in the plot title and y-axis label. Ignored when `comparison_values` is also supplied. Default `"Group 2"`.
#' @param compare_distributions Logical parameter indicating whether to compare the distribution of variants along the length of the gene and its domains between the groups. Default is FALSE.
#' @param gene The gene symbol to plot.
#' @param label_threshold Threshold for labels to appear on plot. Default 5.
#' @param plot_title Optional, the title of the plot. Default is `{gene}: {comparison_values[1]} vs. {comparison_values[2]}`.
#' @param refseq_id Insert a specific NM_xxx value of interest.
#' @param forestarg Logical parameter indicating whether to plot the colollipopplot with or without the forest plot. Default is TRUE. Ignored (with a warning) when using `maf_top`/`maf_bottom`.
#' @param highlight_driver_regions Logical. When `TRUE` and the input MAF
#'   contains a `mutation_alias` column (from
#'   [GAMBLR.utils::annotate_curated_drivers()]), draws semi-transparent
#'   region boxes on both panels. Region boundaries are derived from the
#'   *combined* mutations of both groups so that both panels share identical
#'   boxes. Default `FALSE`.
#' @param driver_region_alpha Transparency of driver region boxes. Default `0.15`.
#' @param driver_region_padding AA positions of padding on each side of a
#'   region's observed range. Default `5`.
#' @param limit_driver_regions Optional character vector of stripped region
#'   names (e.g. `c("KAT", "HAT")`) to restrict which boxes are drawn in
#'   both panels. Applied to the combined region set before passing to each
#'   panel. Default `NULL` (show all regions).
#' @param glm_results Optional data frame returned by
#'   [GAMBLR.utils::select_informative_features], with at least columns
#'   `feature`, `coef`, and `comparison`. When supplied, driver region boxes
#'   are coloured by their GLM coefficient using the same diverging blue–white–
#'   red scale as [gene_feat_heatmap()], replacing the default qualitative
#'   palette. Feature names must match the `mutation_alias` values in the MAF
#'   (e.g. `CREBBP_KAT`). Default `NULL`.
#' @param coef_scale_limits Optional override for the GLM coefficient colour
#'   scale. Supply a single positive number (symmetric, e.g. `0.3` gives
#'   `c(-0.3, 0.3)`) or a length-2 vector (`c(-0.1, 0.5)`) for asymmetric
#'   limits. Values outside the range are squished to the nearest endpoint
#'   colour. When `NULL` (default), limits are `±max(|coef|)` across both
#'   panels.
#' @param glm_comparison Override for the comparison name used to colour the
#'   **top** panel (after stripping `" vs rest"`, e.g. `"FL"`). When `NULL`
#'   (default), it is inferred from `comparison_values[1]` / `name_top`. The
#'   bottom panel always uses `comparison_values[2]` / `name_bottom` if that
#'   label exists in `glm_results`; otherwise it falls back to qualitative
#'   colouring. Both panels share the same colour-scale limits so their
#'   coefficients are directly comparable.
#' @param ... `pretty_lollipop_plot` arguments, see that function for more info on available parameters.
#'
#' @return A mirrored lollipop plot.
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' suppressPackageStartupMessages(library(GAMBLR.open))
#' 
#' #get meta data (DLBCL_Hilton)
#' meta = GAMBLR.open::get_gambl_metadata()
#' metadata = dplyr::filter(meta, cohort %in% "DLBCL_GenomeCanada")
#' maf_df = GAMBLR.open::get_coding_ssm(
#'   these_samples_metadata = metadata
#' )
#' pretty_colollipop_plot_result <- pretty_colollipop_plot(
#'  maf_df = maf_df,
#'  these_samples_metadata = metadata,
#'  comparison_column = "COO_consensus",
#'  comparison_values = c("ABC", "GCB"),
#'  gene = "PIM1",
#'  label_threshold = 2, 
#'  forestarg = TRUE, 
#'  show_rate = TRUE, 
#'  compare_distributions = TRUE
#' )
#' 
#' pretty_colollipop_plot_result$plot
#' 
#' # Load all the metadata available and restrict to FL and DLBCL
#' meta = GAMBLR.open::get_gambl_metadata()
#' metadata = dplyr::filter(meta,pathology %in% c("FL","DLBCL"))
#'
#' #Load mutations for all genome and capture samples in the metadata
#' maf_df = GAMBLR.open::get_coding_ssm(
#'    these_samples_metadata = metadata
#' )
#' # Compare the mutation profile and pattern for CREBBP in FL and DLBCL
#' pretty_colollipop_plot_result <- pretty_colollipop_plot(
#'   maf_df = maf_df,
#'    these_samples_metadata = metadata,
#'    comparison_column = "pathology",
#'    comparison_values = c("FL", "DLBCL"),
#'    gene = "CREBBP",
#'    label_threshold = 2, 
#'    forestarg = TRUE, 
#'    show_rate = TRUE, 
#'    compare_distributions = TRUE
#')
#'
#'pretty_colollipop_plot_result$plot

pretty_colollipop_plot <- function(
    maf_df = NULL,
    these_samples_metadata = NULL,
    comparison_column = NULL,
    comparison_values,
    gene,
    maf_top = NULL,
    maf_bottom = NULL,
    name_top = "Group 1",
    name_bottom = "Group 2",
    include_silent = FALSE,
    label_threshold = 5,
    plot_title,
    refseq_id,
    forestarg = FALSE,
    compare_distributions = FALSE,
    show_rate = FALSE,
    font = "sans",
    highlight_driver_regions = FALSE,
    driver_region_alpha      = 0.15,
    driver_region_padding    = 5,
    limit_driver_regions     = NULL,
    glm_results              = NULL,
    glm_comparison           = NULL,
    coef_scale_limits        = NULL,
    ...
) {
    dual_maf_mode <- !is.null(maf_top) && !is.null(maf_bottom)

    if (dual_maf_mode) {
        if (!missing(comparison_values) && length(comparison_values) != 2) {
            stop("comparison_values must be a character vector of length 2.")
        }
        if (missing(comparison_values)) {
            comparison_values <- c(name_top, name_bottom)
        }
    } else {
        # metadata-based mode: maf_df, these_samples_metadata, comparison_column required
        defined <- names(as.list(match.call())[-1])
        required <- c("maf_df", "these_samples_metadata", "comparison_column")
        if (any(!required %in% defined)) {
            stop("Please provide maf_df, these_samples_metadata, and comparison_column to the function (or provide maf_top and maf_bottom instead).")
        }

        # If no comparison_values are specified, derive them from comparison_column
        if (missing(comparison_values)) {
            if (class(these_samples_metadata[[comparison_column]]) == "factor") {
                comparison_values <- levels(these_samples_metadata[[comparison_column]])
            } else {
                comparison_values <- unique(these_samples_metadata[[comparison_column]])
            }
        }

        if (length(comparison_values) != 2) {
            stop("Your comparison must have two values. \nEither specify comparison_values as a vector of length 2 or subset your metadata so your comparison_column has only two unique values or factor levels.")
        }
    }

    ##### Subset the maf to gene and variants of interest and split by comparison values #####

    if (include_silent) {
        variants <- coding_class
    } else {
        variants <- coding_class[!coding_class %in% c(
            "Silent",
            "Splice_Region"
        )]
    }

    if (dual_maf_mode) {
        ssm1 <- maf_top %>%
            dplyr::filter(Hugo_Symbol == gene, Variant_Classification %in% variants)
        ssm2 <- maf_bottom %>%
            dplyr::filter(Hugo_Symbol == gene, Variant_Classification %in% variants)
    } else {
        # Subset the metadata to the specified comparison_values
        meta1 <- these_samples_metadata[these_samples_metadata[[comparison_column]] %in% comparison_values[1], ]
        meta2 <- these_samples_metadata[these_samples_metadata[[comparison_column]] %in% comparison_values[2], ]

        # Subset maf to only samples in the comparison values
        ssm1 <- maf_df %>%
            dplyr::filter(Tumor_Sample_Barcode %in% meta1$Tumor_Sample_Barcode) %>%
            dplyr::filter(Hugo_Symbol == gene, Variant_Classification %in% variants)

        ssm2 <- maf_df %>%
            dplyr::filter(Tumor_Sample_Barcode %in% meta2$Tumor_Sample_Barcode) %>%
            dplyr::filter(Hugo_Symbol == gene, Variant_Classification %in% variants)
    }

    ##### Shared geometry from combined mutations (both groups) #####
    # used both for the protein-length axis floor (below) and, when
    # highlight_driver_regions = TRUE, for the driver region boxes — computed
    # once here so the top/bottom panels are sized/shaded identically rather
    # than each pretty_lollipop_plot() call deriving it from its own subset
    combined_ssm <- dplyr::bind_rows(ssm1, ssm2) %>%
        dplyr::mutate(
            AA = as.numeric(gsub(
                "[^0-9]+", "",
                gsub("([0-9]+).*", "\\1", HGVSp_Short)
            ))
        )
    protein_length_min <- suppressWarnings(max(combined_ssm$AA, na.rm = TRUE))
    if (!is.finite(protein_length_min)) protein_length_min <- NULL

    region_df <- NULL
    if (highlight_driver_regions) {
        if ("mutation_alias" %in% colnames(combined_ssm)) {
            region_df <- compute_driver_regions(
                combined_ssm, gene, driver_region_padding
            )
            if (nrow(region_df) == 0) region_df <- NULL
        } else {
            warning(
                "highlight_driver_regions = TRUE but 'mutation_alias' column ",
                "not found. Run annotate_curated_drivers() first."
            )
        }
    }
    if (!is.null(region_df) && !is.null(limit_driver_regions)) {
        region_df <- dplyr::filter(region_df, label %in% limit_driver_regions)
        if (nrow(region_df) == 0) region_df <- NULL
    }

    region_df_top <- region_df
    region_df_bot <- region_df
    shared_coef_lim <- NULL

    if (!is.null(region_df) && !is.null(glm_results)) {
        if (!all(c("feature", "coef", "comparison") %in% colnames(glm_results))) {
            warning(
                "glm_results must contain columns 'feature', 'coef', and ",
                "'comparison'; ignoring glm_results."
            )
        } else {
            glm_results$comparison <- sub(" vs rest$", "", glm_results$comparison)
            comps <- unique(glm_results$comparison)

            # Resolve top comparison
            if (!is.null(glm_comparison)) {
                comp_top <- glm_comparison
            } else {
                candidate <- comparison_values[1]
                if (candidate %in% comps) {
                    comp_top <- candidate
                } else if (length(comps) == 1) {
                    comp_top <- comps[1]
                } else {
                    stop(
                        "glm_results contains multiple comparisons (",
                        paste(comps, collapse = ", "),
                        "). Could not auto-detect from name_top ('", candidate,
                        "'). Specify one via glm_comparison."
                    )
                }
            }

            # Resolve bottom comparison
            comp_bot <- {
                candidate <- comparison_values[2]
                if (candidate %in% comps) candidate else NULL
            }
            if (is.null(comp_bot)) {
                message(
                    "Comparison '", comparison_values[2],
                    "' not found in glm_results; bottom panel will use ",
                    "qualitative region colouring."
                )
            }

            attach_glm_coefs <- function(base_df, comp) {
                glm_sub <- dplyr::filter(glm_results, comparison == comp)
                result <- dplyr::left_join(
                    base_df,
                    dplyr::select(glm_sub, feature, coef),
                    by = c("mutation_alias" = "feature")
                )
                result$show_label <- !is.na(result$coef)
                result$coef[is.na(result$coef)] <- 0
                result
            }

            region_df_top <- attach_glm_coefs(region_df, comp_top)
            if (!is.null(comp_bot)) {
                region_df_bot <- attach_glm_coefs(region_df, comp_bot)
            }

            if (!is.null(coef_scale_limits)) {
                shared_coef_lim <- if (length(coef_scale_limits) == 2L)
                    coef_scale_limits
                else
                    c(-coef_scale_limits, coef_scale_limits)
            } else {
                all_coefs <- c(
                    region_df_top$coef,
                    if (!is.null(comp_bot)) region_df_bot$coef else NULL
                )
                lim <- max(abs(all_coefs), na.rm = TRUE)
                if (lim == 0) lim <- 1
                shared_coef_lim <- c(-lim, lim)
            }
        }
    }

    # Ensure dimensions are greater than zero otherwise no variants are returned.
    dim1 <- nrow(ssm1)
    dim2 <- nrow(ssm2)

    if (dim1 == 0 | dim2 == 0) {
        stop("There are no variants returned for one or both of the subsets defined by the specified comparison values.")
    }
    
    # Generate a unified plot title
    if(missing(plot_title)){
        plot_title <- paste0(
            gene, 
            ": ", 
            comparison_values[1], 
            " vs. ", 
            comparison_values[2]
        )
    } 

    # Arguments to pass into pretty_lollipop_plot
    lollipopplot_args <- list(...)
    # Arguments to pass into domain_plot
    domain_matched <-  names(lollipopplot_args)[names(lollipopplot_args) %in% c(
            "domain_label_size", 
            "domain_label_orientation",
            "x_axis_size"
        )]
    if(length(domain_matched) > 0){
        domainplot_args <- lollipopplot_args[domain_matched]
    } else{
        domainplot_args <- list()
    }
    

    # Get gene_counts data for the first plot (without generating the plot)
    lp1 <- do.call(pretty_lollipop_plot, c(
        list(
            maf_df               = ssm1,
            these_samples_metadata = these_samples_metadata,
            gene                 = gene,
            plot_title           = paste0(gene, " in ", comparison_values[1]),
            region_df            = region_df_top,
            driver_region_alpha  = driver_region_alpha,
            coef_limits          = shared_coef_lim,
            protein_length_min   = protein_length_min
        ),
        lollipopplot_args
    ))
    lp1_gene_counts_data <- as.data.frame(lp1$gene_counts) %>%
        mutate(source = "Sample 1")

    # Get gene_counts data for the second plot
    lp2 <- do.call(pretty_lollipop_plot, c(
        list(
            maf_df               = ssm2,
            these_samples_metadata = these_samples_metadata,
            gene                 = gene,
            plot_title           = paste0(gene, " in ", comparison_values[2]),
            mirror               = TRUE,
            region_df            = region_df_bot,
            driver_region_alpha  = driver_region_alpha,
            coef_limits          = shared_coef_lim,
            protein_length_min   = protein_length_min
        ),
        lollipopplot_args
    ))
    lp2_gene_counts_data <- as.data.frame(lp2$gene_counts)%>% 
        mutate(source = "Sample 2")

    # Combine data for both plots
    combined_gene_counts <- rbind(
        lp1_gene_counts_data, 
        lp2_gene_counts_data 
    )

    # Harmonize the y-axis scales for plotting
    y_limit <- max(combined_gene_counts$mutation_count) * 1.2
    y_breaks <- seq(
            0,
            y_limit,
            by = ifelse(y_limit <= 5, 1, round(y_limit * 1.2 / 5))
        )
    
    if(compare_distributions == TRUE) {

        # Gene KS-Test
        # Compare the distribution of variants along the length of the genes between the groups
        ks_test_result <- suppressWarnings(ks.test(
            rep(lp1$gene_counts$AA, lp1$gene_counts$mutation_count),
            rep(lp2$gene_counts$AA, lp2$gene_counts$mutation_count)
        ))

        gene_p_value <- ks_test_result$p.value

        if (!is.null(lp1$domain_data)) {
            # Domain(s) KS-Test
            domain_data <- lp1$domain_data
            domain_data$p.value <- NA
            # Within each domain, use the KS test to determine if the distribution of variants is different
            for (i in seq_len(nrow(domain_data))) {
                domain_name <- domain_data$text.label[i]
                min_val <- domain_data$start.points[i]
                max_val <- domain_data$end.points[i]

                domain_subset1 <- lp1$gene_counts %>%
                    filter(AA >= min_val & AA <= max_val)
                domain_subset2 <- lp2$gene_counts %>%
                    filter(AA >= min_val & AA <= max_val)

                if (length(domain_subset1$AA) == 0 | length(domain_subset2$AA) == 0) {
                    message("Skipping domain ", domain_name, " as it does not have data for both subgroups")
                } else {
                    domain_ks_test <- suppressWarnings(ks.test(
                        rep(domain_subset1$AA, domain_subset1$mutation_count),
                        rep(domain_subset2$AA, domain_subset2$mutation_count)
                    ))
                    domain_data$p.value[i] <- domain_ks_test$p.value
                }
            }

            domain_data <- domain_data %>%
                mutate(
                    p_value = case_when(
                        p.value < 0.001 ~ "***",
                        p.value < 0.01 ~ "**",
                        p.value < 0.05 ~ "*",
                        TRUE ~ ""
                    )
                )

            domain_plot <- do.call(draw_domain_plot, c(
                list(domain_data = domain_data, font = font),
                domainplot_args
            ))
        } else {
            warning("compare_distributions requires protein domain information; domain KS-test skipped.")
            domain_data <- NULL
            domain_plot <- NULL
        }

        plot_subtitle <- paste0("KS Test P=", signif(gene_p_value, digits = 2))
    } else {
        domain_plot <- lp1$domain_plot
        domain_data <- lp1$domain_data
        plot_subtitle = ""
    }
    
    if (show_rate == TRUE) {
        Somatic_Mutation_Numerator <- c(
            length(unique(ssm1$Tumor_Sample_Barcode)),
            length(unique(ssm2$Tumor_Sample_Barcode))
        )
        if (dual_maf_mode) {
            Somatic_Mutation_Denominator <- c(
                length(unique(maf_top$Tumor_Sample_Barcode)),
                length(unique(maf_bottom$Tumor_Sample_Barcode))
            )
        } else {
            Somatic_Mutation_Denominator <- c(
                length(unique(meta1$sample_id)),
                length(unique(meta2$sample_id))
            )
        }
        Somatic_Mutation_Rate <- round(Somatic_Mutation_Numerator / Somatic_Mutation_Denominator * 100, 1)

        yaxis_labels <- paste0(
            comparison_values,
            " Mutation Count\nRate ",
            Somatic_Mutation_Rate,
            "% (",
            Somatic_Mutation_Numerator,
            "/",
            Somatic_Mutation_Denominator,
            ")"
        )
    } else {
        yaxis_labels <- paste(comparison_values, "\nMutation Count")
    }
    
    # Arrange lollipop and domain plots together
    top_plot <- lp1$mutation_plot +
        scale_y_continuous(breaks = y_breaks, limits = c(0, y_limit)) +
        ylab(yaxis_labels[1]) +
        ggtitle(plot_title, subtitle = plot_subtitle) +
        theme(plot.subtitle = element_text(hjust = 0.5))
    bottom_plot <- lp2$mutation_plot +
        scale_y_continuous(breaks = y_breaks, limits = c(y_limit, 0), transform = "reverse") +
        theme(plot.title = element_blank(), plot.subtitle = element_blank()) +
        ylab(yaxis_labels[2])

    if (!is.null(domain_plot)) {
        colollipop_plot <- ggpubr::ggarrange(
            top_plot, domain_plot, bottom_plot,
            ncol = 1, align = "v", heights = c(3, 1, 3),
            common.legend = TRUE, legend = "right"
        )
    } else {
        colollipop_plot <- ggpubr::ggarrange(
            top_plot, bottom_plot,
            ncol = 1, align = "v", heights = c(1, 1),
            common.legend = TRUE, legend = "right"
        )
    }
    
    
    outputs <- list(
        plot = colollipop_plot,
        combined_gene_counts = combined_gene_counts, 
        lollipop1 = lp1$mutation_plot, 
        lollipop2 = lp2$mutation_plot, 
        domain_data = domain_data, 
        domain_plot= domain_plot
    )

    if (forestarg == TRUE) {
        if (dual_maf_mode) {
            warning("forestarg is not supported when using maf_top/maf_bottom; skipping forest plot.")
        } else {
        forest_plot <- prettyForestPlot(
            maf = maf_df,
            metadata = these_samples_metadata,
            genes = gene,
            comparison_column = comparison_column,
            comparison_values = comparison_values,
            separate_hotspots = FALSE,
            comparison_name = paste0(
                comparison_values[1],
                " vs. ",
                comparison_values[2]
            )
        )

    # arrange colollipop plot and forest plot together one over the other
        combine_lollipop_forest <- ggpubr::ggarrange(
            colollipop_plot, 
            ggpubr::ggarrange(
                forest_plot$forest + 
                    ggtitle(
                        paste0("Fisher's Exact Test P=", signif(forest_plot$fisher$p.value, 2), ", OR=", signif(forest_plot$fisher$estimate, 2) )
                        ) + 
                    theme(
                        text = element_text(family = font, face = "plain"), 
                        axis.text = element_text(family = font, face = "plain"), 
                        axis.title.x = element_text(family = font, face = "plain"), 
                        plot.title = element_text(family = font, face = "plain"),
                        axis.title.y = element_blank()
                        ) ,        
                forest_plot$bar + 
                    theme(
                        text = element_text(family = font, face = "plain"), 
                        axis.text = element_text(family = font, face = "plain"), 
                        axis.title.x = element_text(family = font, face = "plain")
                        ), 
                widths = c(1, 0.6),
                common.legend = TRUE,
                legend = "bottom",
                align = "hv"),          
            ncol = 1, 
            nrow = 2,
            heights = c(3, 1), 
            align = "v"     
        )

    
        outputs$plot <- combine_lollipop_forest
        outputs$forest_plot <- forest_plot
        } # end else (not dual_maf_mode)
    }
    return(outputs)
}
