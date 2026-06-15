#' @title Plot feature association results as a forest plot
#'
#' @description Visualise the output of
#' [GAMBLR.utils::test_feature_associations] as a paired forest plot
#' (log odds ratio with 95\% CI) and bar plot (mutation frequency per group).
#'
#' @details This function is the plotting companion to
#' [GAMBLR.utils::test_feature_associations]. It accepts the tidy results
#' tibble directly and produces the same paired forest/bar layout as
#' [prettyForestPlot], but with the statistical computation already done
#' upstream. This separation allows the test results to be inspected,
#' filtered, or modified before plotting.
#'
#' Features are ordered by odds ratio by default. Set
#' \code{keepFeatureOrder = TRUE} to preserve the row order in \code{results}
#' (e.g. to honour an ordering applied by
#' [GAMBLR.utils::order_features_by_gene]).
#'
#' For the bar plot, per-group mutation frequencies are computed from the
#' \code{n_mutated_<group>} columns in \code{results} and the total sample
#' counts supplied via \code{group_sizes} or derived from \code{metadata}.
#' If neither is provided the bar plot is omitted and only the forest plot
#' is returned.
#'
#' @param results A tibble returned by
#' [GAMBLR.utils::test_feature_associations]. Must contain at least
#' \code{feature}, \code{OR}, \code{conf_low}, \code{conf_high},
#' \code{p_value}, and \code{q_value}.
#' @param group_sizes Optional named numeric vector giving the total number
#' of samples per group (e.g. \code{c(FL = 219, DLBCL = 2317)}). Used to
#' compute percentage mutated for the bar plot. When \code{results} contains
#' \code{n_total_<group>} columns (produced by
#' [GAMBLR.utils::test_feature_associations]), those values are used
#' automatically and \code{group_sizes} is only needed to override them.
#' @param metadata Optional metadata data frame. If provided together with
#' \code{comparison_column}, group sizes are derived from it. Superseded by
#' \code{n_total_*} columns in \code{results} when present.
#' @param comparison_column Name of the column in \code{metadata} used to
#' define groups. Required when \code{metadata} is supplied.
#' @param max_q Numeric. Features with \code{q_value} above this threshold
#' are excluded from the plot. Default is 1 (show all).
#' @param keepFeatureOrder Logical. If TRUE, features are plotted in the
#' order they appear in \code{results} rather than sorted by odds ratio.
#' Default is FALSE.
#' @param comparison_name Optional string for the legend title. Defaults to
#' the names of the two groups joined by " vs ".
#' @param custom_colours Optional named vector of colours matching the group
#' labels.
#' @param custom_labels Optional character vector of length 2 providing
#' display labels for the two groups in legend order.
#' @param base_size Numeric font size for axis text. Default is 10.
#' @param bar_label Y-axis label for the bar plot. Default is
#' \code{"percent Mutated"}.
#' @param show_legend Logical. Set to FALSE to suppress the legend.
#' Default is TRUE.
#'
#' @return A named list with elements:
#' \code{forest} (ggplot object), \code{bar} (ggplot object or NULL),
#' \code{arranged} (combined plot or just the forest plot if no bar data),
#' and \code{results} (the filtered results tibble used for plotting).
#'
#' @import dplyr ggplot2 tidyr GAMBLR.helpers
#' @rawNamespace import(ggpubr, except = "get_legend")
#' @export
#'
#' @examples
#' \dontrun{
#' suppressPackageStartupMessages(library(GAMBLR.open))
#' suppressPackageStartupMessages(library(GAMBLR.utils))
#'
#' meta = get_gambl_metadata()
#' fl_dlbcl_meta = dplyr::filter(meta, pathology %in% c("FL", "DLBCL","BL")) %>%
#'   GAMBLR.helpers::check_and_clean_metadata(duplicate_action = "keep_first")
#'
#' maf = get_all_coding_ssm(fl_dlbcl_meta)
#' maf = annotate_curated_drivers(
#'   maf,
#'   genes_of_interest = c("EZH2", "CREBBP", "TP53", "MYD88", "NOTCH1", "NOTCH2","DDX3X","MYC")
#' )
#'
#' results = test_feature_associations(
#'   maf = maf,
#'   metadata = fl_dlbcl_meta,
#'   comparison_column = "pathology",
#'   comparison_values = c("FL", "DLBCL","BL"),
#'   genes = c("EZH2", "CREBBP", "TP53", "MYD88", "NOTCH1", "NOTCH2","DDX3X","MYC"),
#'   maf_column = "mutation_alias"
#' )
#'
#' # group_sizes derived automatically from metadata
#' plots = plot_feature_associations(
#'   results = results,
#'   metadata = fl_dlbcl_meta,
#'   comparison_column = "pathology"
#' )
#' plots$arranged
#'
#' # or supply group sizes directly
#' plots = plot_feature_associations(
#'   results = results,
#'   group_sizes = c(FL = 219, DLBCL = 2317),
#'   comparison_name = "FL vs DLBCL",
#'   max_q = 0.1
#' )
#' plots$arranged
#' }
#'
plot_feature_associations = function(results,
                                     group_sizes = NULL,
                                     metadata = NULL,
                                     comparison_column = NULL,
                                     max_q = 1,
                                     keepFeatureOrder = FALSE,
                                     comparison_name = NULL,
                                     custom_colours = NULL,
                                     custom_labels = NULL,
                                     base_size = 10,
                                     bar_label = "% Mutated",
                                     show_legend = TRUE) {

  required_cols = c("feature", "OR", "conf_low", "conf_high", "p_value", "q_value")
  missing_cols = setdiff(required_cols, colnames(results))
  if (length(missing_cols) > 0) {
    stop("results is missing required columns: ", paste(missing_cols, collapse = ", "),
         "\nDid you use test_feature_associations() with test = 'fisher'?")
  }

  pairwise_mode = "comparison" %in% colnames(results) &&
                  dplyr::n_distinct(results$comparison) > 1

  # detect group names from n_mutated_* columns
  group_count_cols = grep("^n_mutated_[^$]", colnames(results), value = TRUE)
  group_names = sub("^n_mutated_", "", group_count_cols)

  # resolve group_sizes: n_total_* from results takes priority,
  # then metadata, then the caller-supplied group_sizes
  total_cols = paste0("n_total_", group_names)
  if (all(total_cols %in% colnames(results)) && is.null(group_sizes)) {
    group_sizes = vapply(group_names, function(g) {
      as.integer(results[[paste0("n_total_", g)]][[1]])
    }, integer(1))
    names(group_sizes) = group_names
  } else if (!is.null(metadata) && !is.null(comparison_column) &&
             is.null(group_sizes)) {
    gs = table(metadata[[comparison_column]])
    group_sizes = as.integer(gs[group_names])
    names(group_sizes) = group_names
  }

  # filter by max_q and drop infinite OR (complete separation)
  plot_data = dplyr::filter(results, q_value <= max_q, is.finite(OR))

  if (nrow(plot_data) == 0) {
    stop("No features remain after applying max_q = ", max_q,
         ". Try a less stringent threshold.")
  }

  # order features
  if (keepFeatureOrder) {
    feature_levels = unique(as.character(plot_data$feature))
  } else if (pairwise_mode) {
    feature_levels = plot_data %>%
      dplyr::group_by(feature) %>%
      dplyr::summarise(mean_logOR = mean(log(OR), na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(mean_logOR) %>%
      dplyr::pull(feature)
  } else {
    feature_levels = plot_data$feature[order(plot_data$OR)]
  }
  plot_data = dplyr::mutate(plot_data,
                            feature = factor(feature, levels = feature_levels))

  # sizing (base on unique features, not total rows)
  n_features = dplyr::n_distinct(plot_data$feature)
  point_size = min(5, max(1, 50 / n_features))
  error_width = ifelse(point_size >= 5, 0.1, 0.2)
  col_width   = ifelse(point_size >= 5, 0.3, 0.5)

  # colours for bar chart (group fills) — same for both modes
  if (is.null(custom_colours)) {
    gambl_cols = GAMBLR.helpers::get_gambl_colours()
    if (length(intersect(group_names, names(gambl_cols))) == length(group_names)) {
      colours = gambl_cols[group_names]
    } else {
      blood_cols = GAMBLR.helpers::get_gambl_colours(classification = "blood")
      colours = blood_cols[seq_len(length(group_names))]
      names(colours) = group_names
    }
  } else {
    colours = custom_colours
  }

  # labels
  if (is.null(custom_labels)) {
    labels = stats::setNames(group_names, group_names)
  } else if (length(custom_labels) == length(group_names)) {
    labels = stats::setNames(custom_labels, group_names)
  } else {
    warning("custom_labels length does not match number of groups; using group names.")
    labels = stats::setNames(group_names, group_names)
  }

  legend_pos = if (show_legend) "bottom" else "none"

  if (pairwise_mode) {
    comp_levels = unique(as.character(plot_data$comparison))
    plot_data = plot_data %>%
      dplyr::mutate(
        comparison  = factor(comparison, levels = comp_levels),
        comp_group1 = sub(" vs .*", "", as.character(comparison)),
        comp_group2 = sub(".* vs ", "", as.character(comparison)),
        enriched_group = dplyr::case_when(
          conf_low  > 1 ~ comp_group1,
          conf_high < 1 ~ comp_group2,
          TRUE          ~ "none"
        ),
        depleted_group = dplyr::case_when(
          conf_low  > 1 ~ comp_group2,
          conf_high < 1 ~ comp_group1,
          TRUE          ~ "none"
        ),
        enriched_group = factor(enriched_group, levels = c(group_names, "none")),
        depleted_group = factor(depleted_group, levels = c(group_names, "none"))
      )

    comp_shapes = stats::setNames(
      c(15L, 16L, 17L, 18L)[seq_along(comp_levels)],
      comp_levels
    )

    if (is.null(comparison_name)) comparison_name = "Comparison"

    forest = plot_data %>%
      ggplot2::ggplot(ggplot2::aes(x = feature, y = log(OR),
                                   group = comparison)) +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = log(conf_low), ymax = log(conf_high),
                     colour = enriched_group),
        position = ggplot2::position_dodge(width = 0.6),
        width = error_width
      ) +
      ggplot2::geom_point(
        ggplot2::aes(shape = comparison, colour = depleted_group),
        position = ggplot2::position_dodge(width = 0.6),
        size = point_size
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_colour_manual(
        name   = "Enriched in",
        values = c(colours, none = "grey60"),
        breaks = group_names,
        guide  = ggplot2::guide_legend(
          override.aes = list(shape = NA, linewidth = 1)
        )
      ) +
      ggplot2::scale_shape_manual(
        name   = comparison_name,
        values = comp_shapes,
        guide  = ggplot2::guide_legend(
          override.aes = list(colour = "black", linetype = 0)
        )
      ) +
      ggplot2::ylab("ln(Odds Ratio)") +
      ggplot2::xlab("Feature\n") +
      theme_Morons(base_size = base_size) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = base_size))
  } else {
    if (is.null(comparison_name)) {
      comparison_name = paste(group_names, collapse = " vs ")
    }

    # forest plot (single comparison)
    forest = plot_data %>%
      ggplot2::ggplot(ggplot2::aes(x = feature, y = log(OR))) +
      ggplot2::geom_point(size = point_size, shape = "square") +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
      ggplot2::coord_flip() +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = log(conf_low), ymax = log(conf_high)),
        width = error_width
      ) +
      ggplot2::ylab("ln(Odds Ratio)") +
      ggplot2::xlab("Feature\n") +
      theme_Morons(base_size = base_size) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = base_size))
  }

  # bar plot — only if group sizes are available
  if (!is.null(group_sizes) && length(group_count_cols) > 0) {
    bar_data = plot_data %>%
      dplyr::select(feature, dplyr::all_of(group_count_cols)) %>%
      dplyr::distinct() %>%
      tidyr::pivot_longer(
        dplyr::all_of(group_count_cols),
        names_to = "group",
        values_to = "n_mutated"
      ) %>%
      dplyr::mutate(
        group = sub("^n_mutated_", "", group),
        n_total = group_sizes[group],
        percent_mutated = n_mutated / n_total * 100,
        group = factor(group, levels = group_names)
      )

    bar_legend_title = if (pairwise_mode) "Group" else comparison_name

    bar = bar_data %>%
      ggplot2::ggplot(ggplot2::aes(x = feature, y = percent_mutated, fill = group)) +
      ggplot2::geom_col(position = "dodge", width = col_width) +
      ggplot2::xlab("") +
      ggplot2::ylab(bar_label) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(
        name   = bar_legend_title,
        values = colours,
        labels = labels[group_names]
      ) +
      theme_Morons(base_size = base_size) +
      ggplot2::theme(
        axis.text.y  = ggplot2::element_blank(),
        legend.position = legend_pos
      )

    if (pairwise_mode) {
      # Bar Group colours are identical to forest Enriched-in colours so the
      # bar legend is redundant. Use the forest legend as the common legend
      # placed at the bottom of the arranged figure.
      arranged = ggpubr::ggarrange(
        forest,
        bar + ggplot2::theme(legend.position = "none"),
        widths        = c(1, 0.6),
        common.legend = show_legend,
        legend        = "bottom",
        align         = "h"
      )
    } else {
      arranged = ggpubr::ggarrange(
        forest, bar,
        widths = c(1, 0.6),
        common.legend = show_legend,
        align = "h"
      )
    }
  } else {
    if (is.null(group_sizes)) {
      message("No group_sizes or metadata supplied — returning forest plot only.")
    }
    bar = NULL
    arranged = forest
  }

  list(forest = forest, bar = bar, arranged = arranged, results = plot_data)
}
