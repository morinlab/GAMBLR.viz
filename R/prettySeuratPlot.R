#' @title Plot Seurat-style umap plot
#'
#' @description This function uses the output of the `make_and_annotate_umap` to
#'      generate a facet of plots in the single-cell data style of Seurat
#'      package.
#'
#' @param umap_data The output list of the `make_and_annotate_umap` function.
#'      Must contain `df` (a data frame with sample_id, V1, and V2 columns of
#'      the umap coordinates) and `features` (a data frame with feature matrix,
#'      where the rownames are sample identifiers and each column corresponds to
#'      an individual feature).
#' @param features A vector of characters specifying the features to plot. The
#'      order in which the features are specified in this argument will be
#'      respected during the facet plotting. Required argument.
#' @param ncol Integer specifying how many columns should be in the final facet
#'      plot. Default is 3 (the plots will be arranged in 3 columns).
#' @param colours For a custom color palette, provide a named list of colours,
#'      where the item name corresponds to the value from feature matrix and the
#'      item itself is the hex code for the desired colour.
#' @param round_matrix Whether to round values in the provided feature matrix.
#'      Default is FALSE.
#' @import ggplot2 dplyr tibble GAMBLR.helpers
#' @export
#' @return plot
#' @examples
#'
#' \dontrun{
#' gambl_mu_core_all <- make_and_annotate_umap(...)
#'
#' prettySeuratPlot(
#'     gambl_mu_core_all,
#'     features = c(
#'         "EZH2", "CD79B"
#'     ),
#'     ncol = 2
#' )
#'
#' prettySeuratPlot(
#'     gambl_mu_core_all,
#'     features = c(
#'         "EZH2", "CD79B"
#'     ),
#'     ncol = 2,
#'     colours = list(
#'         "0" = "red",
#'         "2" = "gold"
#'     )
#' )
#'
#' prettySeuratPlot(
#'     gambl_mu_core_all,
#'     features = c(
#'         "EZB_feats", "MCD_feats",
#'         "ST2_feats", "BN2_feats"
#'     ),
#'     ncol = 4,
#'     round_matrix = TRUE
#' )
#'
#' }
#'
prettySeuratPlot <- function(
    umap_data,
    features,
    ncol = 3,
    colours = NULL,
    round_matrix = FALSE
){

    if(round_matrix){
        umap_data$features <- round(umap_data$features)
    }

    plot_data <- umap_data$df %>%
    select(sample_id, V1, V2) %>%
    left_join(
        .,
        umap_data$features %>%
            as.data.frame %>%
            rownames_to_column("sample_id") %>%
            select(
                sample_id,  all_of(features)
            )
    ) %>%
    select(sample_id, V1, V2, all_of(features)) %>%
    pivot_longer(
        !c(sample_id, V1, V2)
    ) %>%
    mutate(
        name = factor(name, levels = features),
        value = factor(
            value, levels = sort(unique(value))
        )
    )

    default_colours <- list(
            "0" = "#d2d2d0",
            "1" = "#9870f2",
            "2" = "#4569fd",
            "3" = "#00a6fb"
        )

    allowed <- as.numeric(
        names(default_colours)
    )

    data_vals <- plot_data %>%
        pull(value) %>%
        as.character() %>%
        as.numeric()

    if (!all(data_vals %in% allowed)) {
        bad <- unique(data_vals[!data_vals %in% allowed])
        stop(
            paste0(
                "Invalid values found in the feature matrix:",
                paste(bad, collapse = ", "),
                ".\n",
                "Please provide a named list of colours with the `colours` argument."
            )
        )
    }

    if(is.null(colours)){
        colours_list <- default_colours
    }else{
        colours_list <- colours
    }


    plot <- plot_data %>%
    ggplot(
        aes(
            x = V1, y = V2,
            color = value
        )
    ) +
    geom_point(
        alpha = 0.5,
        stroke = 0.3
    ) +
    scale_colour_manual(
        values = colours_list
    ) +
    facet_wrap(vars(name), ncol = ncol) +
    theme_Morons() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(
            face = "bold", family = "sans"
        )
    )

    return(plot)
}
