#' @title Compare classifications of the same samples across runs
#'
#' @description Alluvial plot of how each sample's class (e.g. its DLBClass
#' cluster or LymphGen subgroup) changes from a reference run to one or more
#' other runs, e.g. classifications from different copy number inputs. Each
#' ribbon is a sample, coloured by its class in the reference, and flows
#' through the other runs in the order given, so consistent samples run
#' straight across and the runs that move samples stand out.
#'
#' @details Classes are stacked in the same order on every axis, with
#' "Other" (and "Missing") at the bottom. Within each class, ribbons are
#' ordered by their reference class. Samples whose class is the same in every
#' run are faded (`fade_consistent`). Each run's axis label gives the number
#' and percentage of samples whose class differs from the reference.
#'
#' @param calls The classifications: a named list of data frames, one per
#' run, each with sample_id and the column named in `class_column` (e.g.
#' [GAMBLR.predict::run_dlbclass] output); or one data frame with sample_id
#' and one column per run.
#' @param reference The run that colours and orders the ribbons. Default the
#' first.
#' @param order The other runs, left to right. Runs not named are left out.
#' Default: all, in the order given.
#' @param reference_position Where the reference's axis goes among the
#' others: 1 (default, leftmost), any position up to the number of runs
#' shown, or "middle" (e.g. two runs on either side of the reference).
#' Ribbons are still coloured by the reference class and changes counted
#' against it.
#' @param pairwise Compare each run with the reference separately: its
#' ribbons join only its axis and the reference's, ordered independently of
#' the other runs, and are faded where its class matches the reference's
#' (with `fade_consistent`), so strong ribbons are exactly that run's
#' differences. Otherwise each sample is one ribbon through every axis, faded
#' when its class is the same in every run. Default NULL: TRUE when the
#' reference has runs on both sides (e.g. `reference_position = "middle"`).
#' @param class_column With a list of data frames, the column holding each
#' sample's class. Default "DLBClass", else "PredictedCluster".
#' @param missing Samples absent from some runs: "drop" (default; only
#' samples in every run are shown, with a message) or "show" (as "Missing").
#' @param fade_consistent Fade the ribbons of samples with the same class in
#' every run (or, with `pairwise`, the same class as the reference in that
#' run). Default TRUE.
#' @param class_order Optional order of the classes (top to bottom). Default:
#' sorted, with "Other" and "Missing" last.
#' @param colours Optional named vector of class colours. Default
#' `GAMBLR.helpers::get_gambl_colours()`, grey for classes it lacks.
#' @param label_strata Label each class box with its name and sample count.
#' Default TRUE.
#' @param legend_position Position of the legend, as in
#' [ggplot2::theme()]: "right" (default), "left", "top", "bottom", or "none"
#' to hide it.
#' @param base_size Base font size. Default 11.
#'
#' @return A ggplot object. The per-run summary (n_samples, n_changed,
#' fraction_changed, relative to the reference) is in `attr(p, "runs")`, and
#' each sample's classes and the number of runs in which it differs from the
#' reference in `attr(p, "samples")`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' calls = GAMBLR.predict::run_dlbclass(list(published = gsm_published,
#'                                           cnv_master = gsm_default,
#'                                           qc_score = gsm_qc))
#' plot_classification_alluvial(calls)
#' # qc_score next to the reference, then cnv_master
#' plot_classification_alluvial(calls, order = c("qc_score", "cnv_master"))
#' # the published calls in the middle, one run on each side, no legend
#' plot_classification_alluvial(calls, reference_position = "middle", legend_position = "none")
#' attr(last_plot(), "runs")
#'
#' # one data frame with a column per run, e.g. LymphGen calls
#' plot_classification_alluvial(lymphgen_calls, class_column = NULL)
#' }
plot_classification_alluvial = function(calls,
                                        reference = NULL,
                                        order = NULL,
                                        reference_position = 1,
                                        pairwise = NULL,
                                        class_column = "DLBClass",
                                        missing = "drop",
                                        fade_consistent = TRUE,
                                        class_order = NULL,
                                        colours = NULL,
                                        label_strata = TRUE,
                                        legend_position = "right",
                                        base_size = 11){
  if(!requireNamespace("ggalluvial", quietly = TRUE)){
    stop("plot_classification_alluvial needs the ggalluvial package: install.packages(\"ggalluvial\")")
  }
  missing = match.arg(missing, c("drop", "show"))

  # one column per run, one row per sample
  if(is.data.frame(calls)){
    if(!"sample_id" %in% colnames(calls)) stop("calls needs a sample_id column")
    wide = as.data.frame(calls, stringsAsFactors = FALSE)
    runs = setdiff(colnames(wide), "sample_id")
  }else{
    if(is.null(names(calls)) || any(names(calls) == "")) stop("give calls as a named list")
    runs = names(calls)
    columns = lapply(runs, function(r){
      x = as.data.frame(calls[[r]])
      col = if(!is.null(class_column) && class_column %in% colnames(x)) class_column
            else if("PredictedCluster" %in% colnames(x)) "PredictedCluster"
            else stop("run ", r, " has no ", class_column, " column")
      if(!"sample_id" %in% colnames(x)) stop("run ", r, " has no sample_id column")
      stats::setNames(data.frame(x$sample_id, as.character(x[[col]]), stringsAsFactors = FALSE),
                      c("sample_id", r))
    })
    wide = Reduce(function(a, b) merge(a, b, by = "sample_id", all = TRUE), columns)
  }
  if(length(runs) < 2) stop("calls needs at least two runs")
  dup = unique(wide$sample_id[duplicated(wide$sample_id)])
  if(length(dup) > 0){
    stop(length(dup), " sample_id value(s) appear more than once (e.g. \"", dup[1], "\"); ",
         "give each sample one row per run (empty or NA sample_ids count as one value)")
  }
  if(is.null(reference)) reference = runs[1]
  if(!reference %in% runs) stop("reference must be one of: ", paste(runs, collapse = ", "))
  others = if(is.null(order)) setdiff(runs, reference) else order
  unknown = setdiff(others, runs)
  if(length(unknown) > 0) stop("order names unknown runs: ", paste(unknown, collapse = ", "))
  others = setdiff(others, reference)
  if(length(others) == 0) stop("order leaves no runs to compare with the reference")
  # the reference's axis among the others
  if(identical(reference_position, "middle")) reference_position = ceiling((length(others) + 1) / 2)
  if(!is.numeric(reference_position) || length(reference_position) != 1 ||
     reference_position < 1 || reference_position > length(others) + 1 || reference_position %% 1 != 0){
    stop("reference_position must be \"middle\" or a whole number from 1 to ", length(others) + 1)
  }
  shown = append(others, reference, after = reference_position - 1)
  wide = wide[, c("sample_id", shown), drop = FALSE]
  wide[shown] = lapply(wide[shown], as.character)

  incomplete = !stats::complete.cases(wide[shown])
  if(any(incomplete)){
    if(missing == "drop"){
      message(sum(incomplete), " sample(s) not in every run are left out (missing = \"show\" keeps them)")
      wide = wide[!incomplete, , drop = FALSE]
    }else{
      wide[shown] = lapply(wide[shown], function(v) ifelse(is.na(v), "Missing", v))
    }
  }
  if(nrow(wide) == 0) stop("no samples are in every run")

  # how each run differs from the reference
  changed = sapply(others, function(r) wide[[r]] != wide[[reference]])
  changed = matrix(changed, nrow = nrow(wide), dimnames = list(NULL, others))
  wide$n_runs_changed = rowSums(changed)
  run_summary = data.frame(run = others, n_samples = nrow(wide),
                           n_changed = colSums(changed),
                           fraction_changed = colMeans(changed),
                           stringsAsFactors = FALSE)
  rownames(run_summary) = NULL

  # class order: sorted, then Other and Missing
  classes = unique(unlist(wide[shown]))
  if(is.null(class_order)){
    tail_classes = intersect(c("Other", "Missing"), classes)
    class_order = c(sort(setdiff(classes, tail_classes)), tail_classes)
  }else{
    class_order = c(class_order, setdiff(classes, class_order))
  }
  if(is.null(colours)) colours = GAMBLR.helpers::get_gambl_colours()
  colours = colours[!duplicated(names(colours))]
  palette = stats::setNames(rep("grey70", length(class_order)), class_order)
  known = intersect(class_order, names(colours))
  palette[known] = colours[known]
  if("Missing" %in% class_order && !"Missing" %in% names(colours)) palette["Missing"] = "white"

  axis_labels = ifelse(shown == reference, paste0(reference, "\n(reference)"),
                        sprintf("%s\n%d changed (%.0f%%)", shown,
                                run_summary$n_changed[match(shown, run_summary$run)],
                                100 * run_summary$fraction_changed[match(shown, run_summary$run)]))
  # pairwise: each other run is compared with the reference on its own (its
  # ribbons join only its axis and the reference's, stacked independently,
  # and fade when the sample's class there is the reference's); otherwise
  # each sample is one ribbon through every axis, faded when its class is
  # the same in every run
  if(is.null(pairwise)) pairwise = reference_position > 1 && reference_position < length(shown)
  fade_levels = if(pairwise) c("differs from reference", "same as reference") else c("changes", "same in every run")
  axis_rows = function(axes, fade){
    do.call(rbind, lapply(axes, function(a){
      data.frame(sample_id = wide$sample_id, run = axis_labels[match(a, shown)], class = wide[[a]],
                 reference_class = wide[[reference]], fade = fade, stringsAsFactors = FALSE)
    }))
  }
  as_factors = function(d){
    d$run = factor(d$run, levels = axis_labels)
    d$class = factor(d$class, levels = class_order)
    d$reference_class = factor(d$reference_class, levels = class_order)
    d$fade = factor(d$fade, levels = fade_levels)
    d
  }
  long = as_factors(axis_rows(shown, ifelse(wide$n_runs_changed == 0, fade_levels[2], fade_levels[1])))
  ribbons = if(pairwise){
    lapply(others, function(r) as_factors(axis_rows(c(reference, r),
                                                    ifelse(wide[[r]] == wide[[reference]], fade_levels[2], fade_levels[1]))))
  } else list(long)

  p = ggplot2::ggplot(long, ggplot2::aes(x = run, stratum = class, alluvium = sample_id))
  # one ribbon per sample per layer (geom_flow draws nothing with ggplot2 4);
  # each layer is stacked on its own, so pairwise layers don't interact
  for(d in ribbons){
    p = p + ggalluvial::geom_alluvium(data = d, ggplot2::aes(fill = reference_class,
                                                             alpha = if(fade_consistent) fade else NULL),
                                      lode.guidance = "frontback", aes.bind = "alluvia",
                                      colour = NA, width = 1/4)
  }
  p = p +
    ggalluvial::geom_stratum(ggplot2::aes(fill = class), width = 1/4, colour = "grey20",
                             linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = palette, breaks = class_order, name = paste0(reference, " class"),
                               drop = FALSE) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.3)) +
    ggplot2::labs(x = NULL, y = "samples",
                  subtitle = sprintf("%d samples; ribbons coloured by the %s class%s", nrow(wide), reference,
                                     if(pairwise) ", each run compared with it separately" else "")) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = legend_position)
  if(fade_consistent){
    p = p + ggplot2::scale_alpha_manual(values = stats::setNames(c(0.85, 0.2), fade_levels),
                                        name = NULL, drop = FALSE)
  }
  if(label_strata){
    p = p + ggplot2::geom_text(stat = ggalluvial::StatStratum,
                               ggplot2::aes(label = ggplot2::after_stat(
                                 ifelse(count >= 0.02 * nrow(wide), paste0(stratum, " (", count, ")"), ""))),
                               size = base_size / 3.6, colour = "white", fontface = "bold")
  }
  attr(p, "runs") = run_summary
  attr(p, "samples") = wide
  p
}
