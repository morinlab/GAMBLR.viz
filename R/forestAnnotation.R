#' @title Add Forest Plot Annotations to an Oncoplot
#'
#' @description Add one or more forest-plot style row annotations to a
#' `prettyOncoplot` or `prettyCoOncoplot` result.
#'
#' @details This helper is designed to combine the output of
#' [prettyOncoplot()] or [prettyCoOncoplot()] with one or more
#' [prettyForestPlot()] result objects. The forest statistics are redrawn
#' directly inside ComplexHeatmap row annotations so the odds ratios align with
#' the oncoplot gene rows.
#'
#' The `prettyoncoplot_output` argument must be the returned object from
#' `prettyOncoplot(return_inputs = TRUE)` or
#' `prettyCoOncoplot(returnEverything = TRUE)`, rather than only the drawn
#' heatmap. The `forestplot_output` argument must be a named list of
#' `prettyForestPlot()` outputs.
#'
#' @param prettyoncoplot_output A list returned by [prettyOncoplot()] with
#' `return_inputs = TRUE`, or by [prettyCoOncoplot()] with
#' `returnEverything = TRUE`.
#' @param forestplot_output A named list of [prettyForestPlot()] outputs. Each
#' element should contain a `fisher` data frame and a `bar` ggplot with a
#' `data` slot containing the comparison groups.
#' @param width Width, in cm, of each forest annotation panel.
#' @param side Side on which to place annotation names. If `prettyoncoplot_output`
#' contains saved arguments, `annotation_name_side` is reused from that object.
#' @param axis_font_size Font size for the forest annotation x-axis. If
#' `prettyoncoplot_output` contains saved arguments, `axis_font_size` is reused
#' from that object.
#' @param name_font_size Font size for annotation names.
#' @param point_size Point size, in mm, for the odds-ratio points.
#' @param asterisks If TRUE, add p-value significance symbols above points.
#' @param return_everything If TRUE, return a list containing the drawn heatmap,
#' the assembled heatmap list, row annotations, and processed forest data.
#' Otherwise, the drawn heatmap object is returned invisibly.
#' @param gap Gap, in mm, between heatmap components.
#'
#' @return Invisibly returns the drawn ComplexHeatmap object. If
#' `return_everything = TRUE`, returns a list with `drawn`, `heatmap_list`,
#' `annotations`, and `forest_data`.
#'
#' @import ComplexHeatmap grid
#' @export
#'
#' @examples
#' cat("Running example for function: forestAnnotation\n")
#' library(GAMBLR.open)
#' library(dplyr)
#'
#' metadata <- get_gambl_metadata(seq_type_filter = "genome") %>%
#'     dplyr::filter(
#'         pathology %in% c("FL", "DLBCL"),
#'         study == "FL_Dreval"
#'     )
#' maf <- get_coding_ssm(these_samples_metadata = metadata)
#' genes <- c("CREBBP","KMT2D","BCL2","SGK1","SOCS1",
#'             "EZH2", "TNFRSF14", "RRAGC")
#'
#' oncoplot <- prettyOncoplot(
#'     maf_df = maf,
#'     these_samples_metadata = metadata,
#'     genes = genes,
#'     keepGeneOrder = TRUE,
#'     minMutationPercent = 0,
#'     simplify_annotation = TRUE,
#'     return_inputs = TRUE
#' )
#' forest <- prettyForestPlot(
#'     maf = maf,
#'     metadata = metadata,
#'     genes = genes,
#'     comparison_column = "pathology",
#'     comparison_values = c("DLBCL", "FL"),
#'     keepGeneOrder = TRUE,
#'     max_q = 1
#' )
#' forestAnnotation(oncoplot, list(pathology = forest), point_size=3)
#'
#' # Side-by-side oncoplot workflow.
#' # co <- prettyCoOncoplot(..., returnEverything = TRUE)
#' # fl_forest <- prettyForestPlot(...)
#' # forestAnnotation(co, list(FL = fl_forest))
forestAnnotation <- function(prettyoncoplot_output,
                             forestplot_output,
                             width = 4,
                             side = "bottom",
                             axis_font_size = 9,
                             name_font_size = 12,
                             point_size = 1,
                             asterisks = FALSE,
                             return_everything = FALSE,
                             gap = 5) {
  validate_forest_annotation_inputs(prettyoncoplot_output, forestplot_output)

  if ("args" %in% names(prettyoncoplot_output)) {
    args <- prettyoncoplot_output$args
    if (!is.null(args$annotation_name_side)) {
      side <- args$annotation_name_side
    }
    if (!is.null(args$axis_font_size)) {
      axis_font_size <- args$axis_font_size
    }
  }

  heatmap_obj <- prettyoncoplot_output$Heatmap
  heatmap_annotations <- list()

  for (res_name in names(forestplot_output)) {
    empty_anno <- ComplexHeatmap::anno_empty(
      border = FALSE,
      width = grid::unit(width, "cm"),
      show_name = TRUE,
      which = "row"
    )

    anno_list <- stats::setNames(list(empty_anno), res_name)
    heatmap_annotations[[res_name]] <- do.call(
      ComplexHeatmap::rowAnnotation,
      c(
        anno_list,
        list(
          annotation_name_rot = 0,
          annotation_name_gp = grid::gpar(fontsize = name_font_size),
          annotation_name_side = side,
          annotation_name_offset = grid::unit(1, "cm")
        )
      )
    )
  }

  heatmap_list <- Reduce(`+`, heatmap_annotations, init = heatmap_obj)
  drawn_heatmap <- ComplexHeatmap::draw(
    heatmap_list,
    merge_legends = TRUE,
    annotation_legend_side = "bottom",
    heatmap_legend_side = "bottom",
    gap = grid::unit(gap, "mm")
  )

  is_main_axis <- side != "top"
  genes_ordered <- rev(prettyoncoplot_output$gene_order)
  forest_data <- list()

  for (anno_name in names(heatmap_annotations)) {
    plot_data <- forestplot_output[[anno_name]]$bar$data
    fisher <- prepare_forest_annotation_data(
      fisher = forestplot_output[[anno_name]]$fisher,
      plot_data = plot_data,
      genes_ordered = genes_ordered
    )
    forest_data[[anno_name]] <- fisher

    draw_forest_annotation(
      anno_name = anno_name,
      fisher = fisher,
      axis_font_size = axis_font_size,
      point_size = point_size,
      asterisks = asterisks,
      is_main_axis = is_main_axis
    )
  }

  if (return_everything) {
    return(
      list(
        drawn = drawn_heatmap,
        heatmap_list = heatmap_list,
        annotations = heatmap_annotations,
        forest_data = forest_data
      )
    )
  }

  invisible(drawn_heatmap)
}

validate_forest_annotation_inputs <- function(prettyoncoplot_output,
                                              forestplot_output) {
  if (!is.list(prettyoncoplot_output)) {
    stop("prettyoncoplot_output must be a list returned by prettyOncoplot() or prettyCoOncoplot().")
  }
  if (!"Heatmap" %in% names(prettyoncoplot_output)) {
    stop("prettyoncoplot_output must contain a Heatmap element.")
  }
  if (!"gene_order" %in% names(prettyoncoplot_output)) {
    stop("prettyoncoplot_output must contain a gene_order element.")
  }
  if (!is.list(forestplot_output) || length(forestplot_output) == 0) {
    stop("forestplot_output must be a non-empty named list of prettyForestPlot() outputs.")
  }
  if (is.null(names(forestplot_output)) || any(names(forestplot_output) == "")) {
    stop("forestplot_output must be a named list.")
  }

  required_fisher_cols <- c("gene", "estimate", "conf.low", "conf.high", "p.value")
  for (res_name in names(forestplot_output)) {
    forest_result <- forestplot_output[[res_name]]
    if (!is.list(forest_result) || !"fisher" %in% names(forest_result)) {
      stop("Each forestplot_output element must contain a fisher data frame.")
    }
    missing_cols <- setdiff(required_fisher_cols, colnames(forest_result$fisher))
    if (length(missing_cols) > 0) {
      stop(
        paste0(
          "forestplot_output[['", res_name, "']]$fisher is missing columns: ",
          paste(missing_cols, collapse = ", ")
        )
      )
    }
    if (!"bar" %in% names(forest_result) || is.null(forest_result$bar$data)) {
      stop("Each forestplot_output element must contain a bar plot with a data slot.")
    }
    if (!"comparison" %in% colnames(forest_result$bar$data)) {
      stop("Each forestplot_output bar data frame must contain a comparison column.")
    }
  }
}

prepare_forest_annotation_data <- function(fisher, plot_data, genes_ordered) {
  fisher <- fisher[order(match(fisher$gene, genes_ordered)), , drop = FALSE]

  fisher$conf.low <- log(fisher$conf.low)
  fisher$conf.high <- log(fisher$conf.high)
  fisher$estimate <- log(fisher$estimate)

  finite_estimates <- fisher$estimate[is.finite(fisher$estimate)]
  finite_conf_low <- fisher$conf.low[is.finite(fisher$conf.low)]
  finite_conf_high <- fisher$conf.high[is.finite(fisher$conf.high)]

  if (length(finite_conf_high) > 0) {
    fisher$conf.high[is.infinite(fisher$conf.high)] <- max(finite_conf_high)
  }
  if (length(finite_conf_low) > 0) {
    fisher$conf.low[is.infinite(fisher$conf.low)] <- min(finite_conf_low)
  }
  if (length(finite_estimates) > 0) {
    fisher$estimate[fisher$estimate == -Inf] <- min(finite_estimates)
    fisher$estimate[fisher$estimate == Inf] <- max(finite_estimates)
  }

  fisher$symbol <- dplyr::case_when(
    fisher$p.value < 0.001 ~ "***",
    fisher$p.value < 0.01 ~ "**",
    fisher$p.value < 0.05 ~ "*",
    TRUE ~ ""
  )

  if (is.factor(plot_data$comparison)) {
    group_levels <- levels(plot_data$comparison)
  } else {
    group_levels <- unique(as.character(plot_data$comparison))
  }
  group_levels <- group_levels[!is.na(group_levels)]
  if (length(group_levels) < 2) {
    stop("forestplot_output bar data must contain two comparison groups.")
  }
  group_levels <- group_levels[1:2]

  fisher$colour <- ifelse(fisher$estimate > 0, group_levels[1], group_levels[2])
  cols <- map_metadata_to_colours(
    metadataColumns = "colour",
    these_samples_metadata = fisher
  )
  colours <- cols$colour[group_levels]

  attr(fisher, "group_levels") <- group_levels
  attr(fisher, "colours") <- colours
  fisher
}

draw_forest_annotation <- function(anno_name,
                                   fisher,
                                   axis_font_size,
                                   point_size,
                                   asterisks,
                                   is_main_axis) {
  group_levels <- attr(fisher, "group_levels")
  colours <- attr(fisher, "colours")

  value <- fisher$estimate
  line_start <- fisher$conf.low
  line_end <- fisher$conf.high
  min_value <- min(c(value, line_start), na.rm = TRUE)
  max_value <- max(c(value, line_end), na.rm = TRUE)
  if (!is.finite(min_value) || !is.finite(max_value)) {
    stop("Forest annotation data must contain at least one finite estimate or confidence limit.")
  }
  if (min_value == max_value) {
    min_value <- min_value - 1
    max_value <- max_value + 1
  }
  max_axis_value <- max(c(abs(round(min_value)), abs(round(max_value))))
  if (max_axis_value == 0) {
    max_axis_value <- 1
  }
  rownum <- nrow(fisher)

  ComplexHeatmap::decorate_annotation(anno_name, {
    y <- seq_len(rownum)

    grid::pushViewport(grid::viewport(
      yscale = c(0.5, rownum + 0.5),
      xscale = c(min_value, max_value)
    ))

    grid::grid.segments(
      x0 = 0, x1 = 0,
      y0 = 0.5, y1 = rownum + 0.5,
      default.units = "native"
    )
    grid::grid.segments(
      x0 = line_start, x1 = line_end,
      y0 = y, y1 = y,
      default.units = "native",
      gp = grid::gpar(col = "gray")
    )
    grid::grid.segments(
      x0 = line_start, x1 = line_start,
      y0 = y - 0.05, y1 = y + 0.05,
      default.units = "native",
      gp = grid::gpar(col = "gray")
    )
    grid::grid.segments(
      x0 = line_end, x1 = line_end,
      y0 = y - 0.05, y1 = y + 0.05,
      default.units = "native",
      gp = grid::gpar(col = "gray")
    )

    if (asterisks) {
      grid::grid.text(
        label = fisher$symbol,
        x = value,
        y = y + 0.25,
        gp = grid::gpar(fontsize = 7, col = "black"),
        default.units = "native"
      )
    }

    grid::grid.points(
      value,
      y,
      pch = 16,
      size = grid::unit(point_size, "mm"),
      gp = grid::gpar(
        col = ifelse(
          value > 0,
          colours[group_levels[1]],
          colours[group_levels[2]]
        )
      ),
      default.units = "native"
    )

    grid::grid.xaxis(
      at = c(
        -1 * max_axis_value,
        -1 * round(max_axis_value / 2),
        0,
        round(max_axis_value / 2),
        max_axis_value
      ),
      label = TRUE,
      name = "log(odds ratio)",
      main = is_main_axis,
      gp = grid::gpar(fontsize = axis_font_size, col = "black")
    )
    grid::popViewport()
  })
}
