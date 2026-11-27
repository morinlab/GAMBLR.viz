#' @title Kaplan-Meier plot.
#'
#' @description Plot Kaplan-Meier survival curves, overall or by group, with
#' censoring marks, a numbers-at-risk table and the log-rank p-value (and
#' optionally Cox hazard ratios).
#'
#' @details The function takes a data frame with one row per patient (e.g.
#' `get_gambl_outcomes()` joined to sample metadata by patient_id); it
#' doesn't fetch any data. Rows without a time, event or group are dropped,
#' with a message.
#'
#' Colours: `custom_colours` if given; otherwise the GAMBL standard colours
#' (`GAMBLR.helpers::get_gambl_colours()`) when every group has one (e.g.
#' cell of origin, LymphGen); otherwise a fixed-order colourblind-safe
#' palette (up to 6 groups). Curves are also labelled at their ends (up to 4
#' groups) and in the risk table, so groups never rely on colour alone.
#'
#' With `strata_col` (e.g. study), the log-rank test and hazard ratios are
#' stratified, so groups are compared within strata; the curves pool the
#' strata.
#'
#' @param data Data frame with one row per patient.
#' @param time_col Column with the follow-up time. Default "OS_YEARS".
#' @param event_col Column with the event code (1 = event, 0 = censored).
#'   Default "CODE_OS".
#' @param group_col Optional column to draw one curve per value of.
#' @param strata_col Optional column to stratify the log-rank test and hazard
#'   ratios by (e.g. "study").
#' @param group_order Optional order of the groups (legend, colours, risk
#'   table); the first is the hazard ratio reference.
#' @param custom_colours Optional named vector of colours for the groups.
#' @param show_risk_table Draw the numbers-at-risk table. Default TRUE.
#' @param risk_times Times for the risk table and x-axis breaks (default:
#'   about 5 even steps over the follow-up).
#' @param show_censor Mark censored patients with "+". Default TRUE.
#' @param show_ci Draw 95% confidence bands. Default FALSE.
#' @param show_pvalue Print the log-rank p-value. Default TRUE.
#' @param show_hr Print Cox hazard ratios versus the first group. Default
#'   FALSE.
#' @param time_label X-axis label. Default "Years".
#' @param survival_label Y-axis label. Default "Overall survival".
#' @param legend_title Legend title. Default: `group_col`.
#' @param title Optional plot title.
#' @param base_size Base font size. Default 12.
#'
#' @return A ggplot object (the curves and the risk table combined with
#'   cowplot when `show_risk_table = TRUE`).
#'
#' @import ggplot2 dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' meta <- GAMBLR.open::get_gambl_metadata(with_outcomes = TRUE) %>%
#'   dplyr::filter(pathology == "DLBCL", study %in% c("dlbcl_chapuy", "dlbcl_schmitz")) %>%
#'   dplyr::distinct(patient_id, .keep_all = TRUE)
#' prettyKaplanMeier(meta, group_col = "COO_consensus", strata_col = "study",
#'                   group_order = c("GCB", "ABC", "UNCLASS"), show_hr = TRUE)
#' }
prettyKaplanMeier = function(data,
                             time_col = "OS_YEARS",
                             event_col = "CODE_OS",
                             group_col = NULL,
                             strata_col = NULL,
                             group_order = NULL,
                             custom_colours = NULL,
                             show_risk_table = TRUE,
                             risk_times = NULL,
                             show_censor = TRUE,
                             show_ci = FALSE,
                             show_pvalue = TRUE,
                             show_hr = FALSE,
                             time_label = "Years",
                             survival_label = "Overall survival",
                             legend_title = group_col,
                             title = NULL,
                             base_size = 12){
  needed = c(time_col, event_col, group_col, strata_col)
  missing_cols = setdiff(needed, colnames(data))
  if(length(missing_cols) > 0) stop("column(s) not in data: ", paste(missing_cols, collapse = ", "))
  d = as.data.frame(data)
  keep = stats::complete.cases(d[, needed, drop = FALSE])
  if(any(!keep)) message(sum(!keep), " row(s) without a time, event", if(!is.null(group_col)) ", group",
                         if(!is.null(strata_col)) " or stratum", " dropped")
  d = d[keep, , drop = FALSE]
  if(nrow(d) == 0) stop("no rows with both a time and an event")
  d$.time = as.numeric(d[[time_col]])
  d$.event = as.numeric(d[[event_col]])
  grouped = !is.null(group_col)
  d$.group = if(grouped) as.character(d[[group_col]]) else "All patients"
  groups = if(!is.null(group_order)) {
    c(intersect(group_order, unique(d$.group)), setdiff(sort(unique(d$.group)), group_order))
  } else sort(unique(d$.group))
  d$.group = factor(d$.group, levels = groups)

  # colours: custom, else GAMBL standard colours if every group has one, else the fixed-order fallback
  fallback = c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
  colours = if(!is.null(custom_colours)) {
    if(!all(groups %in% names(custom_colours))) stop("custom_colours has no colour for: ",
                                                     paste(setdiff(groups, names(custom_colours)), collapse = ", "))
    custom_colours[groups]
  } else {
    gambl = tryCatch(GAMBLR.helpers::get_gambl_colours("all"), error = function(e) NULL)
    if(grouped && !is.null(gambl) && all(groups %in% names(gambl)) && !anyDuplicated(toupper(gambl[groups]))) {
      gambl[groups]
    } else {
      if(length(groups) > length(fallback)) {
        stop(length(groups), " groups is more than the ", length(fallback), " distinguishable colours: combine ",
             "small groups (e.g. into 'Other') or give custom_colours")
      }
      stats::setNames(fallback[seq_along(groups)], groups)
    }
  }
  if(!grouped) colours = stats::setNames("#3A3A3A", groups)

  # Kaplan-Meier step curves, starting at 100% at time 0
  fit = survival::survfit(survival::Surv(.time, .event) ~ .group, data = d)
  s = summary(fit, censored = TRUE)
  curve = data.frame(time = s$time, surv = s$surv, lower = s$lower, upper = s$upper, n.censor = s$n.censor,
                     group = if(is.null(s$strata)) groups[1] else sub("^\\.group=", "", as.character(s$strata)),
                     stringsAsFactors = FALSE)
  starts = data.frame(time = 0, surv = 1, lower = 1, upper = 1, n.censor = 0, group = groups)
  curve = rbind(starts, curve)
  curve = curve[order(match(curve$group, groups), curve$time), , drop = FALSE]
  curve$group = factor(curve$group, levels = groups)

  max_time = max(d$.time)
  # right margin wide enough for the direct labels (about 2% of the x range per character)
  right_pad = if(grouped && length(groups) <= 4) 0.06 + 0.022 * max(nchar(groups)) else 0.03
  if(is.null(risk_times)) risk_times = pretty(c(0, max_time), n = 5)
  risk_times = risk_times[risk_times <= max_time]

  # statistics (GAMBLR.utils): log-rank p, optionally stratified; Cox HRs vs the first group
  annotation = character(0)
  if(grouped && length(groups) > 1){
    lr = GAMBLR.utils::summarize_survival(d, ".time", ".event", group_col = ".group", strata_col = strata_col)$logrank
    if(show_pvalue) annotation = c(annotation, paste0("Log-rank", if(!is.null(strata_col)) " (stratified)", " p ",
                                                      format_p(lr$p)))
    if(show_hr){
      d$.group_chr = as.character(d$.group)
      cox = GAMBLR.utils::survival_cox(d, ".time", ".event", covariates = ".group_chr", strata_col = strata_col,
                                       reference = list(.group_chr = groups[1]))$coefficients
      cox$term = sub("^\\.group_chr", "", cox$term)
      annotation = c(annotation, sprintf("%s vs %s: HR %.2f (95%% CI %.2f-%.2f)", cox$term, groups[1], cox$hr,
                                         cox$hr_lower, cox$hr_upper))
    }
  }

  ink = "#2B2B2B"
  muted = "#6B6B6B"
  p = ggplot(curve, aes(x = time, y = surv, colour = group))
  if(show_ci){
    p = p + geom_ribbon(data = step_ribbon(curve), aes(x = time, ymin = lower, ymax = upper, fill = group),
                        inherit.aes = FALSE, colour = NA, alpha = 0.12, show.legend = FALSE) +
      scale_fill_manual(values = colours)
  }
  p = p + geom_step(linewidth = 0.8, direction = "hv")
  if(show_censor){
    p = p + geom_point(data = curve[curve$n.censor > 0, , drop = FALSE], shape = 3, size = 1.8, stroke = 0.6,
                       show.legend = FALSE)
  }
  # direct labels at each curve's end (up to 4 groups), in text ink next to the coloured line end
  if(grouped && length(groups) <= 4){
    ends = curve %>% dplyr::group_by(group) %>% dplyr::slice_max(time, n = 1, with_ties = FALSE) %>% dplyr::ungroup()
    # labels sit in the right margin, beyond every curve, joined to their curve's end by a short leader
    p = p + ggrepel::geom_text_repel(data = ends, aes(label = group), colour = ink, hjust = 0,
                                     xlim = c(max_time * 1.05, NA), direction = "y", size = base_size / 3.4,
                                     segment.colour = muted, segment.size = 0.3, min.segment.length = 0,
                                     box.padding = 0.15, show.legend = FALSE)
  }
  p = p + scale_colour_manual(values = colours, name = if(grouped) legend_title else NULL) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1),
                       breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02))) +
    scale_x_continuous(breaks = risk_times, expand = expansion(mult = c(0.045, right_pad))) +
    coord_cartesian(xlim = c(0, max_time), clip = "off") +
    labs(x = time_label, y = survival_label, title = title) +
    theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.3),
          axis.line = element_line(colour = muted, linewidth = 0.3), axis.ticks = element_line(colour = muted),
          axis.text = element_text(colour = muted), axis.title = element_text(colour = ink),
          legend.position = if(grouped) "top" else "none", legend.justification = "left",
          legend.text = element_text(colour = ink), legend.title = element_text(colour = ink),
          plot.title = element_text(colour = ink))
  if(length(annotation) > 0){
    p = p + annotate("text", x = 0.02 * max_time, y = 0.04, label = paste(annotation, collapse = "\n"),
                     hjust = 0, vjust = 0, size = base_size / 3.6, colour = ink, lineheight = 1.1)
  }
  if(!show_risk_table) return(p)

  # numbers at risk under the curves, one row per group (labels in text ink; the legend carries colour)
  risk = summary(fit, times = risk_times, extend = TRUE)
  risk_df = data.frame(time = risk$time, n = risk$n.risk,
                       group = if(is.null(risk$strata)) groups[1] else sub("^\\.group=", "", as.character(risk$strata)))
  risk_df$group = factor(risk_df$group, levels = rev(groups))
  r = ggplot(risk_df, aes(x = time, y = group, label = n)) +
    geom_text(size = base_size / 3.6, colour = ink) +
    scale_x_continuous(breaks = risk_times, limits = c(0, max_time),
                       expand = expansion(mult = c(0.045, right_pad))) +
    coord_cartesian(xlim = c(0, max_time), clip = "off") +
    labs(x = NULL, y = NULL, subtitle = "Number at risk") +
    theme_minimal(base_size = base_size) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(colour = ink, margin = margin(r = base_size * 0.6)), plot.subtitle = element_text(colour = muted, size = base_size * 0.85))
  cowplot::plot_grid(p, r, ncol = 1, align = "v", axis = "lr",
                     rel_heights = c(1, 0.12 + 0.07 * length(groups)))
}

# Confidence limits as steps (each value holds until the next event time), for geom_ribbon
step_ribbon = function(curve){
  do.call(rbind, lapply(split(curve, curve$group, drop = TRUE), function(x){
    x = x[order(x$time), , drop = FALSE]
    n = nrow(x)
    if(n < 2) return(x[, c("time", "lower", "upper", "group")])
    data.frame(time = c(rbind(x$time[-n], x$time[-1]), x$time[n]),
               lower = c(rbind(x$lower[-n], x$lower[-n]), x$lower[n]),
               upper = c(rbind(x$upper[-n], x$upper[-n]), x$upper[n]),
               group = x$group[1])
  }))
}

# p-value for display: "< 0.001" or 2 significant digits
format_p = function(p){
  if(is.na(p)) return("NA")
  if(p < 0.001) "< 0.001" else paste("=", signif(p, 2))
}

utils::globalVariables(c("time", "surv", "lower", "upper", "group", "n", ".time", ".event", ".group"))
