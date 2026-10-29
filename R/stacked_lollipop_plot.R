#' @title Stacked Lollipop Plot
#'
#' @description Generates a vertically stacked set of lollipop plots, one per
#'   group, sharing a common x-axis (protein length) and y-axis scale. Unlike
#'   [pretty_colollipop_plot()], panels are never mirrored and any number of
#'   groups can be compared. The protein domain bar is placed at the bottom of
#'   the stack.
#'
#' @param maf_list A named list of MAF data frames, one per group. The names
#'   are used as panel y-axis labels and are matched against
#'   `glm_results$comparison` when colouring driver regions.
#'   Example: `list(FL = fl_maf, BL = bl_maf, DLBCL = dlbcl_maf)`.
#' @param gene The gene symbol to plot.
#' @param glm_results Optional data frame with at least columns `feature`,
#'   `coef`, and `comparison` (as returned by
#'   [GAMBLR.utils::select_informative_features()]). When supplied together
#'   with `highlight_driver_regions = TRUE`, region boxes are coloured per
#'   panel by the matching GLM coefficient using a diverging blue–white–red
#'   scale with limits shared across all panels. Each panel is matched by its
#'   name in `maf_list` after stripping `" vs rest"` from the comparison
#'   column. Default `NULL`.
#' @param highlight_driver_regions Logical. When `TRUE` and the MAFs contain a
#'   `mutation_alias` column (added by
#'   [GAMBLR.utils::annotate_curated_drivers()]), draws semi-transparent
#'   driver region boxes. Region boundaries are derived from the union of all
#'   groups' mutations so that all panels share identical box extents. Default
#'   `FALSE`.
#' @param driver_region_alpha Transparency of region boxes (0–1). Default
#'   `0.15`.
#' @param driver_region_padding AA positions of padding added to each side of
#'   a region's observed range. Default `5`.
#' @param limit_driver_regions Optional character vector of stripped region
#'   names (e.g. `c("KAT", "HAT")`) used to restrict which boxes are drawn
#'   across all panels. Default `NULL` (show all regions).
#' @param label_threshold Minimum mutation recurrence for position labels.
#'   Default `5`.
#' @param plot_title Overall plot title placed on the top panel. Defaults to
#'   `"{gene}: {group1} vs. {group2} vs. ..."`.
#' @param coef_scale_limits Optional override for the GLM coefficient colour
#'   scale. Supply either a single positive number (used symmetrically, e.g.
#'   `0.3` gives limits `c(-0.3, 0.3)`) or a length-2 vector
#'   (`c(-0.1, 0.5)`) for asymmetric limits. Values outside the range are
#'   squished to the nearest endpoint colour. When `NULL` (default), limits
#'   are set to `±max(|coef|)` across all panels automatically.
#' @param ... Additional arguments forwarded to [pretty_lollipop_plot()] for
#'   every panel (e.g. `by_allele`, `point_size_range`, `font`,
#'   `refseq_id`). Do **not** pass `mirror`, `region_df`, `coef_limits`,
#'   `protein_length_min`, `highlight_driver_regions`,
#'   `driver_region_padding`, `limit_driver_regions`, or
#'   `driver_region_alpha` here; use the dedicated
#'   parameters above instead.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{`plot`}{The assembled `ggpubr` figure.}
#'     \item{`lollipop_list`}{Named list of [pretty_lollipop_plot()] outputs,
#'       one per group.}
#'     \item{`domain_plot`}{The protein domain bar plot (taken from the first
#'       group's result), or `NULL` if domain data are unavailable.}
#'   }
#'
#' @import dplyr ggplot2 ggpubr
#' @export
stacked_lollipop_plot <- function(
    maf_list,
    gene,
    glm_results              = NULL,
    highlight_driver_regions = FALSE,
    driver_region_alpha      = 0.15,
    driver_region_padding    = 5,
    limit_driver_regions     = NULL,
    label_threshold          = 5,
    plot_title               = NULL,
    coef_scale_limits        = NULL,
    ...
) {
    # ── input checks ─────────────────────────────────────────────────────────
    if (!is.list(maf_list) || is.null(names(maf_list))) {
        stop("'maf_list' must be a named list of MAF data frames.")
    }
    if (missing(gene)) stop("Please provide a gene.")
    group_names <- names(maf_list)
    n_groups    <- length(group_names)
    if (n_groups < 2) stop("'maf_list' must contain at least 2 groups.")

    if (is.null(plot_title)) {
        plot_title <- paste0(gene, ": ", paste(group_names, collapse = " vs. "))
    }

    # Strip args we control so they can't be double-passed through ...
    lollipop_args <- list(...)
    for (arg in c("maf_df", "gene", "plot_title", "label_threshold",
                  "mirror", "region_df", "coef_limits", "protein_length_min",
                  "highlight_driver_regions", "driver_region_alpha",
                  "driver_region_padding", "limit_driver_regions")) {
        lollipop_args[[arg]] <- NULL
    }

    include_silent <- if ("include_silent" %in% names(lollipop_args))
        isTRUE(lollipop_args[["include_silent"]]) else FALSE
    variants <- if (include_silent) {
        coding_class
    } else {
        coding_class[!coding_class %in% c("Silent", "Splice_Region")]
    }

    # ── shared geometry from the union of all groups' mutations ────────────
    # used both for the protein-length axis floor (below) and, when
    # highlight_driver_regions = TRUE, for the driver region boxes — computed
    # once here so every panel is sized/shaded identically rather than each
    # pretty_lollipop_plot() call deriving it from its own group's subset
    combined_ssm <- dplyr::bind_rows(lapply(maf_list, function(maf) {
        dplyr::filter(maf,
                      Hugo_Symbol == gene,
                      Variant_Classification %in% variants)
    }))
    combined_ssm <- dplyr::mutate(
        combined_ssm,
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

    # ── per-group GLM coefs ───────────────────────────────────────────────
    region_df_list  <- stats::setNames(
        replicate(n_groups, region_df, simplify = FALSE),
        group_names
    )
    shared_coef_lim <- NULL

    if (!is.null(region_df) && !is.null(glm_results)) {
        if (!all(c("feature", "coef", "comparison") %in% colnames(glm_results))) {
            warning(
                "glm_results must contain columns 'feature', 'coef', and ",
                "'comparison'; ignoring glm_results."
            )
        } else {
            glm_results$comparison <- sub(" vs rest$", "", glm_results$comparison)
            comps     <- unique(glm_results$comparison)
            all_coefs <- numeric(0)

            for (grp in group_names) {
                if (grp %in% comps) {
                    glm_sub <- dplyr::filter(glm_results, comparison == grp)
                    rdf <- dplyr::left_join(
                        region_df,
                        dplyr::select(glm_sub, feature, coef),
                        by = c("mutation_alias" = "feature")
                    )
                    rdf$show_label <- !is.na(rdf$coef)
                    rdf$coef[is.na(rdf$coef)] <- 0
                    region_df_list[[grp]] <- rdf
                    all_coefs <- c(all_coefs, rdf$coef)
                } else {
                    message(
                        "Group '", grp, "' not found in glm_results; ",
                        "that panel will use qualitative region colouring."
                    )
                }
            }
            if (!is.null(coef_scale_limits)) {
                shared_coef_lim <- if (length(coef_scale_limits) == 2L)
                    coef_scale_limits
                else
                    c(-coef_scale_limits, coef_scale_limits)
            } else if (length(all_coefs) > 0) {
                lim <- max(abs(all_coefs), na.rm = TRUE)
                if (lim == 0) lim <- 1
                shared_coef_lim <- c(-lim, lim)
            }
        }
    }

    # ── call pretty_lollipop_plot for each group ──────────────────────────
    results <- vector("list", n_groups)
    names(results) <- group_names

    for (grp in group_names) {
        results[[grp]] <- do.call(pretty_lollipop_plot, c(
            list(
                maf_df              = maf_list[[grp]],
                gene                = gene,
                plot_title          = paste0(gene, " in ", grp),
                label_threshold     = label_threshold,
                region_df           = region_df_list[[grp]],
                driver_region_alpha = driver_region_alpha,
                coef_limits         = shared_coef_lim,
                protein_length_min  = protein_length_min
            ),
            lollipop_args
        ))
    }

    # ── unified y-axis ────────────────────────────────────────────────────
    all_counts <- unlist(lapply(results, function(r) r$gene_counts$mutation_count))
    y_limit  <- max(all_counts) * 1.2
    y_breaks <- seq(0, y_limit,
                    by = ifelse(y_limit <= 5, 1, round(y_limit * 1.2 / 5)))

    # ── assemble stacked panels ───────────────────────────────────────────
    mutation_plots <- vector("list", n_groups)
    for (i in seq_along(group_names)) {
        grp <- group_names[i]
        p <- suppressMessages(
            results[[grp]]$mutation_plot +
                scale_y_continuous(breaks = y_breaks, limits = c(0, y_limit)) +
                ylab(paste0(grp, "\nMutation Count"))
        )
        if (i == 1) {
            p <- p + ggtitle(plot_title)
        } else {
            p <- p + theme(plot.title = element_blank())
        }
        mutation_plots[[i]] <- p
    }

    domain_plot <- results[[group_names[1]]]$domain_plot

    if (!is.null(domain_plot)) {
        all_plots <- c(mutation_plots, list(domain_plot))
        heights   <- c(rep(3, n_groups), 1)
    } else {
        all_plots <- mutation_plots
        heights   <- rep(1, n_groups)
    }

    stacked_plot <- ggpubr::ggarrange(
        plotlist      = all_plots,
        ncol          = 1,
        align         = "v",
        heights       = heights,
        common.legend = TRUE,
        legend        = "right"
    )

    return(list(
        plot          = stacked_plot,
        lollipop_list = results,
        domain_plot   = domain_plot
    ))
}
