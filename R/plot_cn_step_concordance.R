#' @title Plot copy number steps at two positions against each other
#'
#' @description Scatter plot of each sample's copy number step at one
#' position against its step at another, from
#' `GAMBLR.utils::cn_step_concordance`, one panel per group, e.g. to show
#' the opposite-sign steps an unbalanced t(14;18) leaves at BCL2 and IGH in
#' translocated cases.
#'
#' @details With `axes = "step"`, points are coloured by their `pattern`
#' (both steps in the expected directions, both the other way, a step at one
#' position only, or none); dashed lines mark `min_step` and the expected
#' quadrants are shaded. With `axes = "position"` (loci given as ranges),
#' each sample with a transition of at least `min_step` at both loci is a
#' point at the two transitions' coordinates, coloured by their directions
#' (reds: the expected combinations) and sized by the copies changed;
#' samples with a transition at one locus only are ticks along that axis.
#' The coordinates are segment boundaries, not breakpoints. With several
#' partners (`b` a list in `cn_step_concordance`), each sample is shown with
#' its best-supported partner: by point shape on step axes (the partner's
#' step oriented), and in one row of panels per partner on position axes.
#' Each panel's title gives the group, the number of samples measured at
#' both positions, how many of those with steps at both go the expected way
#' (binomial p), and the Spearman correlation of the steps.
#'
#' @param concordance Output of `GAMBLR.utils::cn_step_concordance`.
#' @param a_label,b_label Axis labels for the two loci. Default "step at a"
#' and "step at b".
#' @param axes "step" (each sample's step at one locus against the other)
#' or "position" (where in each range the transitions are, for loci given
#' to `cn_step_concordance` as ranges). Default: "position" for ranges,
#' otherwise "step".
#' @param expect "opposite" or "same". Default NULL: the partners' steps are
#' oriented (`b_step_oriented`) so the expected combination is always
#' opposite, whatever the arms involved.
#' @param min_step By default that given to `cn_step_concordance`.
#' @param limit Steps beyond this are drawn at it. Default 4.
#' @param ncol Number of panel columns. Default: all groups side by side.
#' @param base_size Base font size. Default 9.
#'
#' @return A ggplot.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' conc = GAMBLR.utils::cn_step_concordance(seg,
#'          a = list(chrom = "18", pos = 60793500), b = list(chrom = "14", pos = 106329000),
#'          these_samples_metadata = meta, group_by = "bcl2_ba")
#' plot_cn_step_concordance(conc, a_label = "BCL2 step (18q tel - cen)",
#'                          b_label = "IGH step (14q tel - cen)")
#' }
plot_cn_step_concordance = function(concordance, a_label = "step at a", b_label = "step at b",
                                    axes = NULL, expect = NULL, min_step = NULL, limit = 4,
                                    ncol = NULL, base_size = 9){
  # partners' steps are oriented (b_step_oriented) so the expected
  # combination is always opposite in sign; a given expect applies as is
  oriented = is.null(expect) && "b_step_oriented" %in% colnames(concordance$steps)
  expect = if(oriented || is.null(expect)) "opposite" else match.arg(expect, c("opposite", "same"))
  if(is.null(min_step)) min_step = if(is.null(concordance$min_step)) 0.5 else concordance$min_step
  loci = concordance$loci
  partners = if(!is.null(loci$b[["chrom"]])) list(loci$b) else loci$b
  ranged = !is.null(loci) && !is.null(loci$a[["from"]]) &&
    all(vapply(partners, function(p) !is.null(p[["from"]]), logical(1)))
  if(is.null(axes)) axes = if(ranged) "position" else "step"
  axes = match.arg(axes, c("step", "position"))
  if(axes == "position" && !ranged) stop("axes = \"position\" needs loci given as ranges (from, to)")
  d = concordance$steps
  d = d[!is.na(d$pattern), , drop = FALSE]
  s = concordance$summary
  # short lines, so the titles fit over square panels
  s$label = sprintf("%s (n = %d)\n%d of %d with both steps as expected%s\nrho = %s",
                    s$group, s$n_measured, s$n_expected, s$n_expected + s$n_other,
                    ifelse(is.na(s$binomial_p), "", sprintf(" (p = %.2g)", s$binomial_p)),
                    ifelse(is.na(s$spearman_rho), "NA",
                           sprintf("%.2f (p = %.2g)", s$spearman_rho, s$spearman_p)))
  d$panel = factor(s$label[match(d$group, s$group)], levels = s$label)
  theme_steps = theme_bw(base_size = base_size) +
    theme(strip.background = element_blank(), strip.text = element_text(hjust = 0),
          panel.grid.minor = element_blank())
  facets = facet_wrap(~ panel, nrow = if(is.null(ncol)) 1 else NULL, ncol = ncol, drop = FALSE)
  if(axes == "position") return(.plot_step_positions(d, loci$a, partners, a_label, b_label, min_step,
                                                     theme_steps, ncol))
  several = "b_partner" %in% colnames(d) && length(unique(stats::na.omit(d$b_partner))) > 1
  d$x = pmax(-limit, pmin(limit, d$a_step))
  d$y = pmax(-limit, pmin(limit, if(oriented) d$b_step_oriented else d$b_step))
  if(oriented && any(unlist(concordance$expect) == "same")) b_label = paste0(b_label, " (oriented)")
  quadrants = if(expect == "opposite"){
    data.frame(xmin = c(min_step, -Inf), xmax = c(Inf, -min_step),
               ymin = c(-Inf, min_step), ymax = c(-min_step, Inf))
  }else{
    data.frame(xmin = c(min_step, -Inf), xmax = c(Inf, -min_step),
               ymin = c(min_step, -Inf), ymax = c(Inf, -min_step))
  }
  colours = c(expected = "#B2182B", other = "#2166AC", one = "grey55", none = "grey80")
  d$pattern = factor(d$pattern, levels = names(colours))
  d = d[order(d$pattern, decreasing = TRUE), , drop = FALSE]   # expected on top
  ggplot(d, aes(x = x, y = y)) +
    geom_rect(data = quadrants, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = "#B2182B", alpha = 0.06) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = c(-min_step, min_step), linetype = "dashed", colour = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = c(-min_step, min_step), linetype = "dashed", colour = "grey70", linewidth = 0.3) +
    (if(several) geom_point(aes(colour = pattern, shape = b_partner), size = 1.4, alpha = 0.8)
     else geom_point(aes(colour = pattern), size = 1.2, alpha = 0.8)) +
    (if(several) labs(shape = "partner")) +
    scale_colour_manual(values = colours, drop = FALSE, name = NULL,
                        labels = c(expected = "both, expected", other = "both, other",
                                   one = "one position", none = "neither")) +
    facets +
    coord_equal(xlim = c(-limit, limit), ylim = c(-limit, limit)) +
    labs(x = a_label, y = b_label) +
    theme_steps
}

# Transition positions plotted against each other: samples with transitions
# at both loci as points, coloured by the directions of the two steps;
# transitions at one locus only as ticks along that axis
.plot_step_positions = function(d, a_locus, partners, a_label, b_label, min_step, theme_steps,
                                ncol = NULL){
  partner_names = vapply(seq_along(partners), function(k){
    n = partners[[k]][["name"]]
    if(is.null(n)) { if(length(partners) == 1) "b" else paste0("b", k) } else as.character(n)
  }, "")
  if(!"b_partner" %in% colnames(d)) d$b_partner = partner_names[1]
  big_a = !is.na(d$a_step) & abs(d$a_step) >= min_step & !is.na(d$a_pos)
  big_b = !is.na(d$b_step) & abs(d$b_step) >= min_step & !is.na(d$b_pos)
  dir = function(v) ifelse(v > 0, "up", "down")
  both = d[big_a & big_b, , drop = FALSE]
  # colours by whether the directions are the expected ones (per partner)
  both$directions = paste0("a ", dir(both$a_step), ", b ", dir(both$b_step))
  both$fit = ifelse(both$pattern == "expected",
                    ifelse(both$a_step < 0, "expected (a down)", "expected (a up)"),
                    ifelse(both$a_step < 0, "other (a down)", "other (a up)"))
  colours = c("expected (a down)" = "#B2182B", "expected (a up)" = "#F4A582",
              "other (a down)" = "#2166AC", "other (a up)" = "#92C5DE")
  both$fit = factor(both$fit, levels = names(colours))
  a_only = d[big_a & !big_b, , drop = FALSE]
  b_only = d[big_b & !big_a, , drop = FALSE]
  mb = function(v) v / 1e6
  xlim = mb(as.numeric(c(a_locus[["from"]], a_locus[["to"]])))
  # one row of panels per partner, each on its own chromosome's range
  limits = do.call(rbind, lapply(seq_along(partners), function(k)
    data.frame(b_partner = partner_names[k],
               y = mb(as.numeric(c(partners[[k]][["from"]], partners[[k]][["to"]]))))))
  limits = merge(limits, unique(d[, "panel", drop = FALSE]))
  for(x in c("both", "a_only", "b_only", "limits")){
    v = get(x); v$b_partner = factor(v$b_partner, levels = partner_names); assign(x, v)
  }
  # a-only transitions belong to no partner: shown in every partner's row
  a_only = do.call(rbind, lapply(partner_names, function(n) transform(a_only, b_partner = n)))
  a_only$b_partner = factor(a_only$b_partner, levels = partner_names)
  p = ggplot() +
    geom_blank(data = limits, aes(y = y)) +
    geom_rug(data = a_only, aes(x = mb(a_pos)), sides = "b", colour = "grey50", alpha = 0.6,
             length = grid::unit(0.03, "npc")) +
    geom_rug(data = b_only, aes(y = mb(b_pos)), sides = "l", colour = "grey50", alpha = 0.6,
             length = grid::unit(0.03, "npc")) +
    geom_point(data = both, aes(x = mb(a_pos), y = mb(b_pos), colour = fit,
                                size = pmin(abs(a_step) + abs(b_step), 6)), alpha = 0.85) +
    scale_colour_manual(values = colours, drop = FALSE, name = "steps") +
    scale_size_continuous(range = c(1, 3), name = "copies changed\n(both loci)") +
    scale_x_continuous(limits = xlim) +
    labs(x = paste0(a_label, " (Mb)"), y = paste0(if(length(partners) > 1) "partner" else b_label,
                                                  " transition (Mb)")) +
    theme_steps
  if(length(partners) > 1){
    p + facet_grid(b_partner ~ panel, scales = "free_y", switch = "y") +
      theme(strip.placement = "outside")
  }else{
    p + facet_wrap(~ panel, nrow = if(is.null(ncol)) 1 else NULL, ncol = ncol, drop = FALSE)
  }
}
