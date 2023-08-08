#' @title Morons Theme
#'
#' @description Define function for consistent plot theme.
#'
#' @description This function was set up to have a standardized ggplot theme for creating consistent-styled plots.
#' The parameters in this function lets the user control theme parameters such as; font size, font style and legend positions.
# For more info, refer to the parameter descriptions.
#'
#' @param base_size Size of the font on the plot. Defaults to 14.
#' @param base_family Font family to be used on the plot. Defaults to Arial. Always use cairo device when saving the resulting plot!
#' @param my_legend_position Where to draw the legend? Defaults to the bottom of the plot.
#' @param my_legend_direction Which direction to draw the legend? Defaults to horizontal.
#'
#' @return Nothing.
#'
#' @import ggplot2
#' @importFrom ggthemes theme_foundation
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(mpg, aes(displ, hwy, colour = class)) +
#' geom_point() +
#' theme_Morons()
#'
theme_Morons = function(base_size = 14,
                        base_family = "Arial",
                        my_legend_position = "bottom",
                        my_legend_direction = "horizontal"){

  (ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
   theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
         text = element_text(colour = "black"),
         panel.background = element_rect(colour = NA),
         plot.background = element_rect(colour = NA),
         panel.border = element_rect(colour = NA),
         axis.title = element_text(face = "bold", size = rel(1.2)),
         axis.title.y = element_text(angle = 90, vjust = 2),
         axis.title.x = element_text(vjust = -0.2),
         axis.text = element_text(size = base_size, family = base_family),
         axis.line = element_line(colour = "black", size = rel(0.8)),
         axis.ticks = element_line(),
         panel.grid.major = element_line(colour = "#f0f0f0"),
         panel.grid.minor = element_blank(),
         legend.key = element_rect(colour = NA),
         legend.position = my_legend_position,
         legend.direction = my_legend_direction,
         legend.title = element_text(face = "italic"),
         strip.background = element_rect(color = "black", fill = "white", size = 1, linetype = "solid"),
         strip.text = element_text(face = "bold")
    ))
}
