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
#' Features are ordered by odds ratio by default (mean log OR across
#' comparisons in pairwise mode). Set \code{keepFeatureOrder = TRUE} to
#' preserve the row order in \code{results} (e.g. to honour an ordering
#' applied by [GAMBLR.utils::order_features_by_gene]). In pairwise mode,
#' use \code{sort_by} to order by a specific comparison or group instead.
#'
#' For the bar plot, per-group mutation frequencies are computed from the
#' \code{n_mutated_<group>} columns in \code{results} and the total sample
#' counts supplied via \code{group_sizes} or derived from \code{metadata}.
#' If neither is provided the bar plot is omitted and only the forest plot
#' is returned.
#'
#' @param results A tibble returned by
#' [GAMBLR.utils::test_feature_associations] or
#' [GAMBLR.utils::select_informative_features]. Must contain at least
#' \code{feature} and \code{OR} (or \code{OR_glmnet}). \code{conf_low},
#' \code{conf_high}, and \code{q_value} are optional; when absent, error bars
#' and the \code{max_q} filter are skipped automatically.
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
#' @param max_q Numeric. Rows with \code{q_value} strictly above this
#' threshold are excluded from the plot entirely; rows at or below it are
#' always shown and coloured by enriched/comparator group. Default is 1
#' (show all).
#' @param positive_only Logical. If \code{TRUE}, restricts the plotted
#' features to those with a positive \code{coef} (glmnet coefficient) in at
#' least one comparison, mirroring the
#' \code{select_informative_features(positive_only = TRUE)} filter but
#' applied at plot time instead of at training time — useful for dropping
#' all-negative (depletion-only) features from an already-trained model
#' without re-running \code{cv.glmnet}. Requires a \code{coef} column in
#' \code{results} (present when \code{select_informative_features} was
#' called with \code{return_models = TRUE}, or more generally whenever
#' glmnet was used). Default is \code{FALSE}.
#' @param keepFeatureOrder Logical. If TRUE, features are plotted in the
#' order they appear in \code{results} rather than sorted by odds ratio.
#' Default is FALSE.
#' @param sort_by Optional character scalar or vector of group names
#' controlling feature ordering in pairwise mode. A single group name (e.g.
#' \code{"FL"}) places features most enriched in that group at the top. An
#' ordered vector (e.g. \code{c("FL", "DLBCL", "BL")}) produces a
#' hierarchical sort: features enriched in the first group sort to the top,
#' then features enriched in the second group, and so on; features not
#' enriched in any listed group sort to the bottom. In all cases the sort
#' key is the maximum \code{log(reciproc_OR)} across the relevant
#' comparisons. \code{NULL} (default) sorts by the maximum
#' \code{log(reciproc_OR)} across all comparisons. Ignored when
#' \code{keepFeatureOrder = TRUE} or in non-pairwise mode. In
#' \code{mode = "split"}, the same vector also controls left-to-right panel
#' order: each comparison (e.g. \code{"FL vs rest"}) is ranked by the
#' earliest position of either side in \code{sort_by}; comparisons
#' involving groups not listed in \code{sort_by} are placed after the
#' ordered ones, in their original relative order.
#' @param box_colour String controlling what the filled square (point) colour
#' represents in pairwise mode. \code{"enriched"} (default): square fill shows
#' the numerator (enriched group) and CI line shows the depleted group.
#' \code{"depleted"}: square fill shows the comparator group and CI line shows
#' the enriched group. The legends swap titles accordingly.
#' @param comparison_name Optional string for the bar chart legend title in
#' non-pairwise mode. Defaults to the two group names joined by " vs ".
#' @param custom_colours Optional named vector of colours matching the group
#' labels.
#' @param custom_labels Optional character vector of length 2 providing
#' display labels for the two groups in legend order.
#' @param base_size Numeric font size for axis text. Default is 10.
#' @param bar_label Y-axis label for the bar plot. Default is
#' \code{"percent Mutated"}.
#' @param show_legend Logical. Set to FALSE to suppress the legend.
#' Default is TRUE.
#' @param return_data Logical. Set to TRUE to include \code{plot_data} (the
#' processed data frame used for plotting, with reciprocated OR display columns)
#' in the returned list. Default is FALSE.
#' @param sort_comparison Optional string giving a single comparison label
#' (e.g. \code{"FL vs DLBCL"}) to use as the reference when computing the
#' feature sort order in pairwise mode. When supplied, the sort key is derived
#' only from rows belonging to that comparison rather than taking the maximum
#' across all comparisons. Particularly useful in \code{mode = "split"} to
#' anchor the shared feature order to one specific pair. Ignored when
#' \code{keepFeatureOrder = TRUE} or in non-pairwise mode.
#' @param use_glmnet Logical. Set to \code{TRUE} to use the \code{OR_glmnet}
#' column (exponentiated regularized logistic regression coefficient) from
#' [GAMBLR.utils::select_informative_features] output as the plotted value
#' instead of a Fisher OR. Error bars are suppressed (no CI available for
#' regularized estimates), and the x-axis is labelled \code{"Regularized
#' ln(OR)"}. Auto-detected when \code{OR} is absent but \code{OR_glmnet}
#' is present, so passing \code{select_informative_features()} output directly
#' works without setting this flag. Default is \code{FALSE}.
#' @param mode One of \code{"combined"} (default) or \code{"split"}.
#' In \code{"combined"} mode all pairwise comparisons are overlaid on a single
#' dodged forest plot using the two-colour encoding. In \code{"split"} mode a
#' separate paired forest+bar plot is produced for each comparison; the feature
#' order is derived from the full combined data (respecting \code{sort_by} and
#' \code{keepFeatureOrder}) and is shared across all plots so rows are aligned
#' when panels are placed side by side. The left-to-right panel order also
#' follows \code{sort_by} (see above). Ignored in non-pairwise mode.
#' @param coef_limits Numeric vector of length 2 giving the lower and upper
#' bound used to truncate the plotted \code{log(OR)} (or, in pairwise mode,
#' \code{log(reciproc_OR)}) value before plotting. This prevents a single
#' extreme value (e.g. a near-complete-separation feature in regularized
#' glmnet output, which can produce a regularized ln(OR) far outside the
#' range of all other features) from compressing the visual scale for
#' every other point in the panel. Truncation is applied only to the
#' plotted point position; the underlying \code{results} values (and the
#' feature sort order) are unaffected. Set to \code{NULL} to disable
#' truncation. Default is \code{c(-5, 5)}.
#'
#' @return In \code{"combined"} mode (or non-pairwise mode): a named list with
#' elements \code{forest} (ggplot object), \code{bar} (ggplot object or NULL),
#' \code{arranged} (combined plot or just the forest plot if no bar data),
#' and \code{results} (the filtered results tibble used for plotting).
#' If \code{return_data = TRUE}, also includes \code{plot_data}.
#' In \code{"split"} mode: a named list keyed by comparison label
#' (e.g. \code{"FL vs DLBCL"}), where each element is the same list structure
#' described above, plus a \code{combined} element containing all forest plots
#' arranged side-by-side with a single shared bar plot on the right.
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
#' genes_driver = dplyr::filter(GAMBLR.data::lymphoma_genes,
#'                              DLBCL_Tier == 1 | FL_Tier == 1 | BL_Tier == 1) %>%
#'                              pull(Gene)
#' genes_hotspot = c("EZH2", "CREBBP", "TP53", "MYD88", "NOTCH1", "NOTCH2","DDX3X","MYC")
#' annotated_maf= annotate_curated_drivers(maf_data = maf,
#'                                           genes_of_interest = unique(c(genes_hotspot,
#'                                                                        genes_driver)))
#'
#'  
#' results = test_feature_associations(
#'   maf = annotated_maf,
#'   metadata = fl_dlbcl_meta,
#'   comparison_column = "pathology",
#'   comparison_values = c("FL", "DLBCL","BL"),
#'   maf_column = "mutation_alias"
#' )
#'
#' # group_sizes derived automatically from metadata
#' plots = plot_feature_associations(
#'   results = results,
#'   comparison_column = "pathology",
#'   max_q = 0.1,
#'   sort_by = c("FL", "DLBCL","BL"),
#'   base_size = 6
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
#' 
#' #combined mode with separate panels per comparison
#' plot_feature_associations(results,max_q = 0.00001,
#'  sort_by = "FL",
#'  mode="split",
#'   base_size = 6)
#' }
#' 
#'
plot_feature_associations = function(results,
                                     group_sizes = NULL,
                                     metadata = NULL,
                                     comparison_column = NULL,
                                     max_q = 1,
                                     positive_only = FALSE,
                                     keepFeatureOrder = FALSE,
                                     sort_by = NULL,
                                     box_colour = "enriched",
                                     comparison_name = NULL,
                                     custom_colours = NULL,
                                     custom_labels = NULL,
                                     base_size = 10,
                                     bar_label = "% Mutated",
                                     show_legend = TRUE,
                                     return_data = FALSE,
                                     sort_comparison = NULL,
                                     use_glmnet = FALSE,
                                     mode = "combined",
                                     coef_limits = c(-5, 5)) {

  mode = match.arg(mode, c("combined", "split"))

  # truncate a log(OR)-scale value to coef_limits; NULL disables truncation
  clip_log_or = function(x) {
    if (is.null(coef_limits)) return(x)
    pmin(pmax(x, coef_limits[1]), coef_limits[2])
  }

  # auto-detect glmnet output: OR absent but OR_glmnet present
  glmnet_mode = use_glmnet ||
    (!"OR" %in% colnames(results) && "OR_glmnet" %in% colnames(results))

  if (glmnet_mode) {
    if (!"OR_glmnet" %in% colnames(results)) {
      stop("use_glmnet = TRUE but 'OR_glmnet' column not found. ",
           "Pass output from select_informative_features() with run_fisher = FALSE.")
    }
    results = dplyr::mutate(results, OR = OR_glmnet)
  }

  if (!"feature" %in% colnames(results) || !"OR" %in% colnames(results)) {
    stop("results must contain at least 'feature' and 'OR' (or 'OR_glmnet') columns.",
         if (!glmnet_mode) "\nDid you use test_feature_associations() with test = 'fisher'?" else "")
  }

  has_ci     = !glmnet_mode && all(c("conf_low", "conf_high") %in% colnames(results))
  has_qvalue = "q_value" %in% colnames(results)

  pairwise_mode = "comparison" %in% colnames(results) &&
                  dplyr::n_distinct(results$comparison) > 1

  # detect group names from n_mutated_* columns
  group_count_cols = grep("^n_mutated_[^$]", colnames(results), value = TRUE)
  group_names = sub("^n_mutated_", "", group_count_cols)

  # When run_fisher = FALSE there are no n_mutated_* columns; derive group names
  # from the comparison labels so colours and factor levels are populated correctly.
  if (length(group_names) == 0 && "comparison" %in% colnames(results)) {
    comp_strings = unique(as.character(results$comparison))
    all_sides = unique(c(
      sub(" vs .*", "", comp_strings),
      sub(".* vs ", "", comp_strings)
    ))
    group_names = setdiff(all_sides, "rest")
  }
  # In split-mode recursive calls the comparison column is stripped from
  # sub_results, so the block above never fires.  Fall back to comparison_name.
  if (length(group_names) == 0 && !is.null(comparison_name)) {
    group_names = setdiff(
      unique(c(sub(" vs .*", "", comparison_name),
               sub(".* vs ", "", comparison_name))),
      "rest"
    )
  }

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

  # filter by max_q when q_value is available; in glmnet mode it is absent so
  # all selected features pass through (selection already serves as the filter).
  # Drop Inf OR (complete separation) but keep NA OR — ggplot2 silently omits
  # them, leaving blank space at that feature position in split mode.
  plot_data = if (has_qvalue) {
    dplyr::filter(results, q_value <= max_q, is.finite(OR) | is.na(OR))
  } else {
    dplyr::filter(results, is.finite(OR) | is.na(OR))
  }

  # drop features that are never positive (i.e. depletion-only across every
  # comparison) — same rule as select_informative_features(positive_only),
  # applied here so it doesn't require re-fitting the (slow) glmnet models
  if (positive_only) {
    if (!"coef" %in% colnames(plot_data)) {
      stop("positive_only = TRUE requires a 'coef' column (glmnet coefficient) ",
           "in results. Did you train with select_informative_features()?")
    }
    keep_features = plot_data %>%
      dplyr::group_by(feature) %>%
      dplyr::summarise(any_positive = any(coef > 0, na.rm = TRUE), .groups = "drop") %>%
      dplyr::filter(any_positive) %>%
      dplyr::pull(feature)
    plot_data = dplyr::filter(plot_data, feature %in% keep_features)
  }

  if (nrow(plot_data) == 0) {
    stop("No features remain after applying max_q = ", max_q,
         if (positive_only) " and positive_only = TRUE" else "",
         ". Try a less stringent threshold.")
  }

  # pairwise: add display columns where OR < 1 is reciprocated so all plotted
  # values exceed 1.  Original OR/CI are preserved for sorting (signed, so
  # sort_by="FL" puts FL-enriched at top) and for enriched/depleted colour logic.
  if (pairwise_mode) {
    plot_data = plot_data %>%
      dplyr::mutate(
        comp_group1 = sub(" vs .*", "", as.character(comparison)),
        comp_group2 = sub(".* vs ", "", as.character(comparison)),
        reciproc_OR = dplyr::if_else(OR < 1, 1 / OR, OR)
      )
    if (has_ci) {
      plot_data = plot_data %>%
        dplyr::mutate(
          .cl                = conf_low,
          reciproc_conf_low  = dplyr::if_else(OR < 1, 1 / conf_high, conf_low),
          reciproc_conf_high = dplyr::if_else(OR < 1, 1 / .cl,       conf_high),
          .cl                = NULL
        )
    }
  }

  # order features
  if (keepFeatureOrder) {
    feature_levels = unique(as.character(plot_data$feature))
  } else if (pairwise_mode) {
    # optionally restrict sort-key computation to a single reference comparison
    if (!is.null(sort_comparison)) {
      if (!sort_comparison %in% as.character(plot_data$comparison)) {
        warning("sort_comparison '", sort_comparison,
                "' not found in results; using all comparisons for sort.")
        sort_data = plot_data
      } else {
        sort_data = dplyr::filter(plot_data,
                                  as.character(comparison) == sort_comparison)
      }
    } else {
      sort_data = plot_data
    }
    if (is.null(sort_by)) {
      # sort by max log(reciproc_OR) across sort_data comparisons
      feature_levels = sort_data %>%
        dplyr::group_by(feature) %>%
        dplyr::summarise(sort_key = max(log(reciproc_OR), na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(sort_key) %>%
        dplyr::pull(feature) %>%
        as.character()
    } else {
      invalid = setdiff(sort_by, group_names)
      if (length(invalid) > 0) {
        warning("sort_by value(s) '", paste(invalid, collapse = "', '"),
                "' do not match any group name; ignoring them.")
        sort_by = intersect(sort_by, group_names)
      }
      # hierarchical sort: assign each feature to the first group in sort_by
      # where it has an enriched comparison.
      # log(reciproc_OR) is always >= 0 in the enriched subset, so OR = 0
      # (complete separation) correctly maps to Inf rather than -Inf.
      tier_keys = lapply(seq_along(sort_by), function(i) {
        grp = sort_by[i]
        grp_data = dplyr::filter(
          sort_data,
          grepl(grp, as.character(comparison), fixed = TRUE) &
            ((comp_group1 == grp & OR > 1) | (comp_group2 == grp & OR < 1))
        )
        if (nrow(grp_data) == 0) return(NULL)
        grp_data %>%
          dplyr::group_by(feature) %>%
          dplyr::summarise(sort_key = max(log(reciproc_OR), na.rm = TRUE), .groups = "drop") %>%
          dplyr::mutate(tier = i)
      })
      all_tiers = dplyr::bind_rows(tier_keys) %>%
        dplyr::group_by(feature) %>%
        dplyr::slice_min(tier, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(dplyr::desc(tier), sort_key)
      feature_levels = as.character(all_tiers$feature)
      # features not enriched in any listed group go to the bottom,
      # sorted by min log(OR) ascending (most depleted at bottom,
      # least depleted just below the enriched block)
      remaining = setdiff(unique(as.character(plot_data$feature)), feature_levels)
      if (length(remaining) > 0) {
        remainder_order = sort_data %>%
          dplyr::filter(as.character(feature) %in% remaining) %>%
          dplyr::group_by(feature) %>%
          dplyr::summarise(sort_key = min(log(OR), na.rm = TRUE),
                           .groups = "drop") %>%
          dplyr::arrange(sort_key) %>%
          dplyr::pull(feature) %>%
          as.character()
        # features absent from sort_data (e.g. only in other comparisons) go lowest
        remainder_order = c(
          setdiff(remaining, remainder_order),
          remainder_order
        )
        feature_levels = c(remainder_order, feature_levels)
      }
    }
  } else {
    feature_levels = plot_data$feature[order(plot_data$OR)]
  }
  plot_data = dplyr::mutate(plot_data,
                            feature = factor(feature, levels = feature_levels))

  # sizing, colours, and labels — computed here so they are available inside
  # the split block for building the shared bar from the full plot_data
  n_features = dplyr::n_distinct(plot_data$feature)
  point_size = min(5, max(1, 50 / n_features))
  error_width = ifelse(point_size >= 5, 0.1, 0.2)
  col_width   = ifelse(point_size >= 5, 0.3, 0.5)

  if (is.null(custom_colours)) {
    if (!is.null(comparison_column)) {
      # Build a minimal stub so map_metadata_to_colours can pick the right
      # palette by column name (e.g. "lymphgen" → lymphgen colours)
      stub_meta          = data.frame(
        v = factor(group_names, levels = group_names),
        stringsAsFactors = FALSE
      )
      names(stub_meta)   = comparison_column
      col_map  = map_metadata_to_colours(
        metadataColumns        = comparison_column,
        these_samples_metadata = stub_meta
      )
      colours = col_map[[comparison_column]][group_names]
    } else {
      gambl_cols = GAMBLR.helpers::get_gambl_colours()
      if (length(intersect(group_names, names(gambl_cols))) == length(group_names)) {
        colours = gambl_cols[group_names]
      } else {
        blood_cols     = GAMBLR.helpers::get_gambl_colours(classification = "blood")
        colours        = blood_cols[seq_len(length(group_names))]
        names(colours) = group_names
      }
    }
  } else {
    colours = custom_colours
  }

  if (is.null(custom_labels)) {
    labels = stats::setNames(group_names, group_names)
  } else if (length(custom_labels) == length(group_names)) {
    labels = stats::setNames(custom_labels, group_names)
  } else {
    warning("custom_labels length does not match number of groups; using group names.")
    labels = stats::setNames(group_names, group_names)
  }

  legend_pos = if (show_legend) "bottom" else "none"

  # split mode: recurse once per comparison using the shared feature order
  if (pairwise_mode && mode == "split") {
    drop_cols = c("comparison", "comp_group1", "comp_group2",
                  "reciproc_OR", "reciproc_conf_low", "reciproc_conf_high")
    # union of all features across comparisons (in shared sort order)
    all_features = levels(plot_data$feature)
    comparisons = unique(as.character(plot_data$comparison))
    # order panels by sort_by: each comparison is ranked by the earliest
    # position (of either side, e.g. "FL vs rest") in sort_by; comparisons
    # involving groups absent from sort_by sort last, in their original order
    if (!is.null(sort_by)) {
      comp_side1 = sub(" vs .*", "", comparisons)
      comp_side2 = sub(".* vs ", "", comparisons)
      rank1 = match(comp_side1, sort_by)
      rank2 = match(comp_side2, sort_by)
      panel_rank = pmin(ifelse(is.na(rank1), Inf, rank1),
                        ifelse(is.na(rank2), Inf, rank2))
      comparisons = comparisons[order(panel_rank, seq_along(comparisons))]
    }
    split_plots = lapply(setNames(comparisons, comparisons), function(comp) {
      pair = strsplit(comp, " vs ")[[1]]
      sub_results = plot_data %>%
        dplyr::filter(as.character(comparison) == comp) %>%
        dplyr::select(-dplyr::any_of(drop_cols)) %>%
        # expand to full union: NA rows for features absent from this comparison
        dplyr::mutate(feature = factor(feature, levels = all_features)) %>%
        tidyr::complete(feature) %>%
        # give NA rows q_value = 0 so they survive the max_q filter in the
        # recursive call; OR stays NA so ggplot2 leaves the row blank
        { if (has_qvalue) dplyr::mutate(., q_value = dplyr::coalesce(q_value, 0)) else . } %>%
        dplyr::arrange(feature)
      sub_meta = if (!is.null(metadata) && !is.null(comparison_column))
        dplyr::filter(metadata, .data[[comparison_column]] %in% pair)
      else
        metadata
      plot_feature_associations(
        results           = sub_results,
        # pass already-resolved group_sizes so the recursive call doesn't try
        # to read from NA-filled n_total_* rows
        group_sizes       = group_sizes,
        metadata          = sub_meta,
        comparison_column = comparison_column,
        max_q             = max_q,
        keepFeatureOrder  = TRUE,
        sort_by           = NULL,
        box_colour        = box_colour,
        comparison_name   = comp,
        custom_colours    = custom_colours,
        custom_labels     = custom_labels,
        base_size         = base_size,
        bar_label         = bar_label,
        show_legend       = FALSE,
        return_data       = return_data,
        use_glmnet        = glmnet_mode,
        mode              = "combined",
        coef_limits       = coef_limits
      )
    })

    # combined panel: all forests side-by-side, one shared bar on the right
    legend_pos_split = if (show_legend) "bottom" else "none"
    # When a shared bar chart is available, feature labels are shown there
    # instead (see bar_shared below) so every forest panel can be blanked
    # identically — no panel carries extra margin, so all render the same
    # width. Without a bar, panel 1 keeps the labels since nothing else
    # can carry them.
    has_shared_bar = !is.null(group_sizes) && length(group_count_cols) > 0
    forests = lapply(seq_along(split_plots), function(i) {
      f = split_plots[[i]]$forest +
        ggplot2::ggtitle(names(split_plots)[i]) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                          size = base_size))
      if (has_shared_bar || i > 1) {
        f = f + ggplot2::theme(axis.text.y  = ggplot2::element_blank(),
                               axis.title.y = ggplot2::element_blank(),
                               axis.ticks.y = ggplot2::element_blank())
      }
      f
    })
    # build the shared bar from the full union plot_data so features only
    # significant in one comparison still get bars from all groups; it also
    # carries the feature labels when present (see has_shared_bar above)
    bar_shared = NULL
    if (has_shared_bar) {
      bar_shared = plot_data %>%
        dplyr::select(feature, dplyr::all_of(group_count_cols)) %>%
        dplyr::distinct() %>%
        tidyr::pivot_longer(
          dplyr::all_of(group_count_cols),
          names_to  = "group",
          values_to = "n_mutated"
        ) %>%
        dplyr::mutate(
          group           = sub("^n_mutated_", "", group),
          n_total         = group_sizes[group],
          percent_mutated = n_mutated / n_total * 100,
          group           = factor(group, levels = group_names)
        ) %>%
        ggplot2::ggplot(ggplot2::aes(x = feature, y = percent_mutated,
                                     fill = group)) +
        ggplot2::geom_col(position = "dodge", width = col_width) +
        ggplot2::xlab("Feature\n") +
        ggplot2::ylab(bar_label) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(
          name   = "Group",
          values = colours,
          labels = labels[group_names]
        ) +
        theme_Morons(base_size = base_size) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = base_size),
                       legend.position = legend_pos)
    }
    forests_no_legend = lapply(forests,
      function(f) f + ggplot2::theme(legend.position = "none"))
    if (!is.null(bar_shared)) {
      top_row = ggpubr::ggarrange(
        plotlist = c(forests_no_legend,
                     list(bar_shared + ggplot2::theme(legend.position = "none"))),
        nrow     = 1,
        widths   = c(rep(1, length(forests)), 0.6),
        align    = "h"
      )
      if (show_legend) {
        bar_legend = ggpubr::get_legend(
          bar_shared + ggplot2::theme(legend.position = "bottom")
        )
        split_plots$combined = ggpubr::ggarrange(
          top_row,
          ggpubr::as_ggplot(bar_legend),
          nrow    = 2,
          heights = c(20, 2)
        )
      } else {
        split_plots$combined = top_row
      }
    } else {
      split_plots$combined = ggpubr::ggarrange(
        plotlist = forests_no_legend,
        nrow     = 1,
        align    = "h"
      )
    }
    return(split_plots)
  }

  if (pairwise_mode) {
    comp_levels = unique(as.character(plot_data$comparison))
    # scale cap size so it stays visible regardless of feature count
    error_width = min(0.6, max(0.3, 10 / n_features))
    plot_data = plot_data %>%
      dplyr::mutate(
        comparison     = factor(comparison, levels = comp_levels),
        # enriched/depleted derived from OR direction; q_value gate only when available
        enriched_group = dplyr::case_when(
          OR > 1 & (!has_qvalue | q_value < max_q) ~ comp_group1,
          OR < 1 & (!has_qvalue | q_value < max_q) ~ comp_group2,
          TRUE                                      ~ "none"
        ),
        depleted_group = dplyr::case_when(
          OR > 1 & (!has_qvalue | q_value < max_q) ~ comp_group2,
          OR < 1 & (!has_qvalue | q_value < max_q) ~ comp_group1,
          TRUE                                      ~ "none"
        ),
        enriched_group = factor(enriched_group, levels = c(group_names, "none")),
        depleted_group = factor(depleted_group, levels = c(group_names, "none"))
      ) %>%
      # fill missing feature × comparison combos with NA so position_dodge
      # always allocates equal slot width regardless of which ORs were finite
      tidyr::complete(feature, comparison) %>%
      dplyr::mutate(
        feature                  = factor(feature, levels = feature_levels),
        # truncate the plotted value so one extreme estimate (e.g. a
        # near-complete-separation glmnet coefficient) doesn't compress the
        # visual scale for every other feature; underlying OR is untouched
        y_reciproc_OR            = clip_log_or(log(reciproc_OR)),
        y_reciproc_conf_low      = if (has_ci) clip_log_or(log(reciproc_conf_low))  else NA_real_,
        y_reciproc_conf_high     = if (has_ci) clip_log_or(log(reciproc_conf_high)) else NA_real_
      )

    if (box_colour == "enriched") {
      line_var  = "depleted_group";  line_title = "Comparator"
      point_var = "enriched_group";  point_title = "Enriched in"
    } else {
      line_var  = "enriched_group";  line_title = "Enriched in"
      point_var = "depleted_group";  point_title = "Comparator"
    }

    forest = plot_data %>%
      ggplot2::ggplot(ggplot2::aes(x = feature, y = y_reciproc_OR,
                                   group = comparison)) +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
      { if (has_ci)
          ggplot2::geom_errorbar(
            ggplot2::aes(ymin = y_reciproc_conf_low, ymax = y_reciproc_conf_high,
                         colour = .data[[line_var]]),
            position = ggplot2::position_dodge(width = 0.6),
            width = error_width
          )
        else
          ggplot2::geom_blank()
      } +
      ggplot2::geom_point(
        ggplot2::aes(fill = .data[[point_var]],
                     colour = ggplot2::after_scale(fill)),
        position = ggplot2::position_dodge(width = 0.6),
        size = point_size,
        shape = 22
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_colour_manual(
        name   = line_title,
        values = c(colours, none = "grey60"),
        breaks = group_names,
        guide  = if (has_ci)
          ggplot2::guide_legend(override.aes = list(shape = NA, linewidth = 1))
        else
          "none"
      ) +
      ggplot2::scale_fill_manual(
        name   = point_title,
        values = c(colours, none = "grey60"),
        breaks = group_names,
        guide  = ggplot2::guide_legend(
          override.aes = list(shape = 22, size = 3)
        )
      ) +
      ggplot2::ylab(if (glmnet_mode) "Regularized ln(OR)" else "ln(Odds Ratio)") +
      ggplot2::xlab("Feature\n") +
      theme_Morons(base_size = base_size) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = base_size))
  } else {
    if (is.null(comparison_name)) {
      comparison_name = paste(group_names, collapse = " vs ")
    }

    # extract the two groups from comparison_name rather than group_names[1/2]
    # because in split mode sub-results still carry all groups' n_mutated_* cols
    pair = strsplit(comparison_name, " vs ")[[1]]
    g1 = pair[1]; g2 = pair[2]

    plot_data = plot_data %>%
      dplyr::mutate(
        enriched_group = dplyr::case_when(
          OR > 1 ~ g1,
          OR < 1 ~ g2,
          TRUE   ~ "none"
        ),
        depleted_group = dplyr::case_when(
          OR > 1 ~ g2,
          OR < 1 ~ g1,
          TRUE   ~ "none"
        ),
        enriched_group = factor(enriched_group, levels = c(group_names, "none")),
        depleted_group = factor(depleted_group, levels = c(group_names, "none")),
        # truncate the plotted value so one extreme estimate (e.g. a
        # near-complete-separation glmnet coefficient) doesn't compress the
        # visual scale for every other feature; underlying OR is untouched
        y_OR            = clip_log_or(log(OR)),
        y_conf_low      = if (has_ci) clip_log_or(log(conf_low))  else NA_real_,
        y_conf_high     = if (has_ci) clip_log_or(log(conf_high)) else NA_real_
      )

    if (box_colour == "enriched") {
      line_var   = "depleted_group"; line_title  = "Comparator"
      point_var  = "enriched_group"; point_title = "Enriched in"
    } else {
      line_var   = "enriched_group"; line_title  = "Enriched in"
      point_var  = "depleted_group"; point_title = "Comparator"
    }

    forest = plot_data %>%
      ggplot2::ggplot(ggplot2::aes(x = feature, y = y_OR)) +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
      { if (has_ci)
          ggplot2::geom_errorbar(
            ggplot2::aes(ymin = y_conf_low, ymax = y_conf_high,
                         colour = .data[[line_var]]),
            width = error_width
          )
        else
          ggplot2::geom_blank()
      } +
      ggplot2::geom_point(
        ggplot2::aes(fill = .data[[point_var]],
                     colour = ggplot2::after_scale(fill)),
        size = point_size, shape = 22
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_colour_manual(
        name     = line_title,
        values   = c(colours, none = "grey60"),
        breaks   = group_names,
        na.value = "grey60",
        drop     = FALSE,
        guide    = if (has_ci)
          ggplot2::guide_legend(override.aes = list(shape = NA, linewidth = 1))
        else
          "none"
      ) +
      ggplot2::scale_fill_manual(
        name     = point_title,
        values   = c(colours, none = "grey60"),
        breaks   = group_names,
        na.value = "grey60",
        drop     = FALSE,
        guide    = ggplot2::guide_legend(
          override.aes = list(shape = 22, size = 3)
        )
      ) +
      ggplot2::ylab(if (glmnet_mode) "Regularized ln(OR)" else "ln(Odds Ratio)") +
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

  out = list(forest = forest, bar = bar, arranged = arranged, results = plot_data)
  if (return_data) out$plot_data = plot_data
  out
}
