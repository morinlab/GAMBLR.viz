#' @title Lollipop Plot
#'
#' @description Generates a visually appealing and interactive lollipop plot.
#'
#' @details This function is depending on a modified version of `readMAF` from the `g3viz` package.
#' Returned plot is interactive, meaning the user can hover over individual points in the plot to reveal more information.
#' The plot can also be exported in a variety of file formats inside the interactive view of the lollipop plot.
#'
#' @param maf_df A data frame containing the mutation data (from a MAF).
#' @param gene The gene symbol to plot.
#' @param plot_title Optional (defaults to gene name).
#' @param plot_theme Options: cbioportal(default), blue, simple, nature, nature2, ggplot2, and dark.
#' @param out_name Optional, set the file name of the plot, if you export it to disk. Default name is my_lollipop_plot_{gene}.
#'
#' @return Nothing.
#'
#' @import g3viz dplyr
#' @export
#'
#' @examples
#' #get metadata (Fl and DLBCL)
#' metadata = GAMBLR.data::gambl_metadata
#' this_metadata = dplyr::filter(metadata, consensus_pathology %in% c("FL", "DLBCL"))
#'
#' #get maf data for returned samples
#' maf = get_coding_ssm(limit_samples = this_metadata$sample_id,
#'                      seq_type = "genome")
#'
#' #construct pretty_lollipop_plot.
#' pretty_lollipop_plot(maf_df = maf,
#'                      gene = "MYC",
#'                      plot_title = "Mutation data for MYC",
#'                      plot_theme = "nature2")
#'
pretty_lollipop_plot = function(maf_df,
                                gene,
                                plot_title,
                                plot_theme = "cbioportal",
                                out_name = paste0("my_lollipop_plot_", gene)){
  if(missing(gene)){
    stop("Plese provide a gene...")
  }

  if(missing(plot_title)){
    plot_title = gene
  }

  maf_df = maf_df %>%
    dplyr::filter(Hugo_Symbol == gene)

  #use the readMAF function (modified by Ryan) to parse/convert
  maf_df = g3viz::readMAF(maf.df = maf_df)
  chart.options = g3Lollipop.theme(theme.name = plot_theme, title.text = plot_title)
  g3Lollipop(maf_df,
             gene.symbol = gene,
             plot.options = chart.options,
             output.filename = out_name)
}
