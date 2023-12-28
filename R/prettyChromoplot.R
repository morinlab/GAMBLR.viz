#' @title Chromosome Plot
#'
#' @description Use GISTIC2.0 scores output to reproduce maftools::chromoplot with more flexibility.
#'
#' @details This function uses GISTIC2.0 scores to create a chromosome plot, based on a similar plotting function from `maftools`.
#' The only required parameter for this function is `scores`, which is the path to a file with GISTIC2.0 scores.
#' Other parameters are all optional. For a detailed explanation of how to use these, refer to the parameter descriptions.
#'
#' @param scores Output file scores.gistic from the run of GISTIC2.0
#' @param genes_to_label Optional. Provide a data frame of genes to label (if mutated). The first 3 columns must contain chromosome, start, and end coordinates. Another required column must contain gene names and be named `gene`. All other columns are ignored. If no data frame provided, oncogenes from GAMBLR packages are used by default to annotate on the plot.
#' @param cutoff Optional. Used to determine which regions to color as aberrant. Must be float in the range between 0 and 1. The higher the number, the less regions will be considered as aberrant. The default is 0.5.
#' @param adjust_amps Optional. The value of G-score for highest amplification peak will be multiplied by this value to determine how far up the gene label will be displayed. Default 0.5.
#' @param adjust_dels Optional. The value of G-score for highest deletion peak will be multiplied by this value to determine how far down the gene label will be displayed. Default 2.75.
#' @param label_size Optional. The font size for the gene label to be displayed. Default 3.
#' @param force_pull Optional. How strong the gene name label will be pulled towards a data point. Default 0 (no pulling).
#' @param segment.curvature Optional. Indicates whether arrow to the data point should be curved. Accepts numeric value, where negative is for left-hand and positive for right-hand curves, and 0 for straight lines. Default 0.25.
#' @param segment.ncp Optional. Indicates number of control points to make a smoother curve. Higher value allows for more flexibility for the curve. Default 4.
#' @param segment.angle Optional. Numeric value in the range 0-180, where less than 90 skews control points of the arrow from label to data point toward the start point. Default 25.
#'
#' @return plot
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose", "melt", "dcast"))
#' @import dplyr ggplot2 ggrepel
#' @export
#'
#' @examples
#' \dontrun{
#' #basic usage
#' prettyChromoplot("path_to_gistic_results/scores.gistic")
#'
#' #advanced usages
#' prettyChromoplot(scorees = "path_to_gistic_results/scores.gistic",
#'                  genes_to_label = "path_to_gene_coordinates_table.tsv",
#'                  cutoff = 0.75) +
#'                   ... #any ggplot options to customize plot appearance
#' }
#'
prettyChromoplot = function(scores,
                            genes_to_label,
                            cutoff = 0.5,
                            adjust_amps = 0.5,
                            adjust_dels = 2.75,
                            label_size = 3,
                            force_pull = 0,
                            segment.curvature = 0.25,
                            segment.ncp = 4,
                            segment.angle = 25){

  #read GISTIC scores file, convert G-score to be negative for deletions, and relocate chromosome, start, and end columns to be the first three
  scores = data.table::fread(scores) %>%
    dplyr::mutate(`G-score` = ifelse(Type == "Amp", `G-score`, - 1 * `G-score`)) %>%
    dplyr::relocate(Type, .after = frequency)

  #annotate each region with direction of changes - used for coloring
  scores$fill = ifelse(scores$Type == "Amp" & scores$`-log10(q-value)` > cutoff, "up",
                ifelse(scores$Type == "Del" & scores$`-log10(q-value)` > cutoff, "down", "neutral"))

  #colors to plot
  cnv_palette = c("up" = "#bd0000", "down" = "#2e5096", "neutral" = "#D2D2D3")

  #if no file is provided, annotate with oncogenes in GAMBLR package
  if(missing(genes_to_label)){
    genes_to_label = GAMBLR.data::grch37_oncogene %>%
      dplyr::mutate(across(c(chrom, start, end), as.integer)) %>%
      data.table::as.data.table()
  }else{
    genes_to_label = data.table::fread(genes_to_label)
    colnames(genes_to_label)[1:3] = c("chrom", "start", "end")
    genes_to_label = genes_to_label %>%

      #for now, drop the X chromosome since GISTIC runs without sex chromosmes
      dplyr::filter(!grepl("X", chrom)) %>%
      dplyr::mutate(across(c(chrom, start, end), as.integer)) %>%
      data.table::as.data.table()
  }
  #overlap scores with genes to annotate
  scores = data.table::foverlaps(scores %>%
    data.table::setkey(., Chromosome, Start, End), genes_to_label %>%
    data.table::setkey(., chrom, start, end), by.x = c("Chromosome", "Start", "End"), by.y = c("chrom", "start", "end"), type = "within") %>%

    #if gene to annotate is provided, but it is in region with no CNV, do not label it
    dplyr::mutate(gene=ifelse(!is.na(gene) & fill=="neutral", NA, gene)) %>%
    #if gene is covering multiple adjacent regions, label only once
    dplyr::group_by(gene) %>%
    dplyr::mutate(newcol = ifelse(!is.na(gene) & !duplicated(gene), gene, NA), gene = newcol) %>%
    dplyr::select(-newcol)

  #get coordinates to label chromosome numbers
  xses = scores %>%
    dplyr::group_by(Chromosome) %>%
    dplyr::mutate(End = max(End) / 2) %>%
    dplyr::pull(End)

  #main plotting
  ggplot(data = scores, aes(x = Start, y = `G-score`, color = fill, label = gene)) +
    geom_bar(size = 0.2, stat = 'identity', position = "dodge") +
    ylab('G-score') +
    ggrepel::geom_text_repel(data = subset(scores, !is.na(gene) & Type == "Amp"),
                             nudge_y = max(subset(scores, !is.na(gene) & Type == "Amp")$`G-score`) * adjust_amps,
                             size = label_size,
                             segment.size = 0.5,
                             segment.color = "#000000",
                             force_pull = force_pull,
                             arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
                             segment.curvature = segment.curvature,
                             segment.ncp = segment.ncp,
                             segment.angle = segment.angle) +
    ggrepel::geom_text_repel(data = subset(scores, !is.na(gene) & Type == "Del"),
                             nudge_y = min(subset(scores, !is.na(gene) & Type == "Del")$`G-score`) * adjust_dels,
                             nudge_x = subset(scores, !is.na(gene) & Type == "Del")$Start,
                             size = label_size,
                             segment.size = 0.5,
                             segment.color = "#000000",
                             force_pull = force_pull,
                             arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
                             segment.curvature = segment.curvature,
                             segment.ncp = segment.ncp,
                             segment.angle = segment.angle) +
    facet_grid(. ~ Chromosome, scales = "free") +
    scale_color_manual(values = cnv_palette) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 16, colour = "black"),
          axis.ticks.x = element_blank(), axis.ticks.y = element_line(colour = "black"), legend.position = "none",
          panel.spacing.x = unit(0.1, "lines"), panel.border = element_blank(), text = element_text(size = 16, colour = "black"),
          strip.background = element_blank(), strip.text.x = element_blank(), panel.grid = element_blank()) +
    geom_hline(yintercept = 0, size = 7) +
    geom_text(aes(label = Chromosome, x = xses, y = 0), size = 4, color = "white")
}
