#' @title Chromosome Plot
#'
#' @description Use GISTIC2.0 scores output to reproduce maftools::chromoplot with more flexibility.
#'
#' @details This function uses GISTIC2.0 scores to create a chromosome plot, based on a similar plotting function from `maftools`.
#' The only required parameter for this function is `scores`, which is the path to a file with GISTIC2.0 scores.
#' Other parameters are all optional. For a detailed explanation of how to use these, refer to the parameter descriptions.
#'
#' @param scores_path Output file scores.gistic from the run of GISTIC2.0
#' @param scores_df Optional. Instead of specifying scores_path pass a pre-loaded scores file as a data frame using scores_df
#' @param labels_bed Optional. A bed_data object specifying the regions to apply labels.
#' @param genome_build Defines the chr prefix and the coordinates of the default genes to label if `genes_to_label` is not provided.
#' Automatically set if labels_bed is provided
#' @param cutoff Optional. Used to determine which regions to color as aberrant. Must be float in the range between 0 and 1. The higher the number, the less regions will be considered as aberrant. The default is 0.5.
#' @param adjust_amps Optional. The value of G-score for highest amplification peak will be multiplied by this value to determine how far up the gene label will be displayed. Default 0.5.
#' @param adjust_dels Optional. The value of G-score for highest deletion peak will be multiplied by this value to determine how far down the gene label will be displayed. Default 2.75.
#' @param label_size Optional. The font size for the gene label to be displayed. Default 3.
#' @param force_pull Optional. How strong the gene name label will be pulled towards a data point. Default 0 (no pulling).
#' @param segment.curvature Optional. Indicates whether arrow to the data point should be curved. Accepts numeric value, where negative is for left-hand and positive for right-hand curves, and 0 for straight lines. Default 0.25.
#' @param segment.ncp Optional. Indicates number of control points to make a smoother curve. Higher value allows for more flexibility for the curve. Default 4.
#' @param segment.angle Optional. Numeric value in the range 0-180, where less than 90 skews control points of the arrow from label to data point toward the start point. Default 25.
#' @param hide_neutral Optional. Set to TRUE to hide all neutral (insignificant) regions instead of plotting them in grey
#' @param verbose
#' @return plot
#'
#' @import dplyr ggplot2 ggrepel readr GAMBLR.helpers
#' @export
#'
#' @examples
#' \dontrun{
#' genes = c("MYC","FCGR2B","TNFRSF14","FAS","PTEN","B2M","RB1","TCL1A","CD70",
#'   "BCL2","KLHL14","TCF4","REL","BCL6","HIST1H1C","SMARCA4","CDKN2A","RHOA",
#'   "TNFAIP3","TP53","CDK14","RELN","ETS1","MDM1","MIR17HG","CD58","HNRNPD",
#'   "TOX","PRAME","CD38")
#' gene_bed = select(grch37_gene_coordinates,-1) %>% #remove ensembl ID column
#'   dplyr::filter(hugo_symbol %in% genes) %>% #keep genes of interest
#'   mutate(length = end - start,mid = start + length/2) %>%
#'   mutate(start = mid,end=start+1) %>%
#'   unique() %>%
#'   create_bed_data(genome_build = "grch37") #convert to bed_data format
#' # GISTIC run using grch37
#' prettyChromoplot(scores_path = "scores.gistic",
#'                  labels_bed = gene_bed)
#' #NOTE: genome build is inferred from gene_bed
#' 
#'  # GISTIC run using hg38 data
#' prettyChromoplot(scores_path="scores.gistic",
#'                    cutoff = 0.9,
#'                    label_size=2,
#'                    adjust_amps = 0.5,
#'                    adjust_dels = 0.8,
#'                    genome_build="hg38",
#'                    hide_neutral = T)
#' }
prettyChromoplot <- function(scores_path,
                             scores_df,
                             labels_bed,
                             genome_build,
                             cutoff = 0.5,
                             adjust_amps = 0.5,
                             adjust_dels = 2.75,
                             label_size = 3,
                             force_pull = 0,
                             segment.curvature = 0.25,
                             segment.ncp = 4,
                             segment.angle = 25,
                             hide_neutral = FALSE,
                             verbose = FALSE) {
  if (!missing(scores_path)) {
    # read GISTIC scores file, convert G-score to be negative for deletions, and relocate chromosome, start, and end columns to be the first three
    scores <- read_tsv(scores_path) %>%
      dplyr::mutate(`G-score` = ifelse(Type == "Amp", `G-score`, -1 * `G-score`)) %>%
      dplyr::relocate(Type, .after = frequency) %>%
      mutate(Chromosome = as.character(Chromosome))
  } else {
    if (!missing(scores_df)) {
      if (!"G-score" %in% colnames(scores_df)) {
        stop("scores_df doesn't contain the expected columns e.g. G-score")
      }
      scores <- scores_df %>%
        dplyr::mutate(`G-score` = ifelse(Type == "Amp", `G-score`, -1 * `G-score`)) %>%
        dplyr::relocate(Type, .after = frequency) %>%
        mutate(Chromosome = as.character(Chromosome))
    }
  }

  # annotate each region with direction of changes - used for coloring
  scores$fill <- ifelse(scores$Type == "Amp" & scores$`-log10(q-value)` > cutoff, "up",
    ifelse(scores$Type == "Del" & scores$`-log10(q-value)` > cutoff, "down", "neutral")
  )
  # drop all neutral to avoid genes being dropped when a del and amp peak both overlap
  if (hide_neutral) {
    scores_labeling <- dplyr::filter(scores, fill != "neutral")
  } else {
    scores_labeling <- scores
  }

  # colors to plot
  cnv_palette <- c("up" = "#bd0000", "down" = "#2e5096", "neutral" = "#D2D2D389")


  if (missing(labels_bed)) {
    if (missing(genome_build)) {
      stop("genome_build is required when labels_bed is not used")
    }
    if (genome_build == "grch37") {
      genes_to_label <- create_bed_data(GAMBLR.data::grch37_oncogene)
    } else {
      genes_to_label <- create_bed_data(GAMBLR.data::hg38_oncogene)
    }
  } else {
    if (!"bed_data" %in% class(labels_bed)) {
      stop("labels_bed must be a bed_data object. ?create_bed_data for more details")
    }
    genome_build <- get_genome_build(labels_bed)
    genes_to_label <- labels_bed %>%
      # drop the X chromosome since GISTIC runs without sex chromosmes
      dplyr::filter(!grepl("X", chrom)) %>%
      dplyr::mutate(across(c(start, end), as.integer)) %>%
      strip_genomic_classes()
  }
  # if no file is provided, annotate with oncogenes in GAMBLR package
  if (genome_build == "grch37") {
    chromord <- c(1:22) %>% as.character()
    cytobands <- cytobands_grch37
  } else {
    chromord <- paste0("chr", c(1:22))
    cytobands <- cytobands_hg38
    if (!any(grepl("chr", scores$Chromosome))) {
      if (verbose) {
        print("adding chr prefix")
      }

      scores <- mutate(scores, Chromosome = paste0("chr", Chromosome))
      scores_labeling <- mutate(scores_labeling, Chromosome = paste0("chr", Chromosome))
    }
  }

  cytobands <- dplyr::filter(cytobands, !grepl("X", cb.chromosome)) %>%
    dplyr::filter(!grepl("Y", cb.chromosome))
  # overlap scores with genes to annotate
  scores <- cool_overlaps(
    scores_labeling,
    genes_to_label,
    columns1 = c("Chromosome", "Start", "End"),
    columns2 = c("chrom", "start", "end"),
    type = "any",
    nomatch = TRUE
  )
  print(scores)
  print(genes_to_label)
  scores <- scores %>%
    # if gene to annotate is provided, but it is in region with no CNV, do not label it
    dplyr::mutate(name = ifelse(!is.na(name) & fill == "neutral", NA, name)) %>%
    # if gene is covering multiple adjacent regions, label only once
    dplyr::group_by(name) %>%
    dplyr::mutate(newcol = ifelse(!is.na(name) & !duplicated(name), name, NA), name = newcol) %>%
    dplyr::select(-newcol)
  scores <- mutate(scores,
    Chromosome = factor(Chromosome, levels = chromord)
  )
  # get coordinates to label chromosome numbers

  cytobands <- mutate(cytobands, cb.chromosome = factor(cb.chromosome, levels = chromord))
  xses <- cytobands %>%
    arrange(cb.chromosome) %>%
    dplyr::group_by(cb.chromosome) %>%
    dplyr::summarise(End = max(cb.end), Mid = End / 2) %>%
    dplyr::rename("Chromosome" = "cb.chromosome")

  # main plotting
  ggplot(data = dplyr::filter(scores, fill == "neutral"), aes(x = Start, y = `G-score`, color = fill, label = name)) +
    geom_bar(size = 0.2, stat = "identity", position = "dodge") +
    geom_bar(data = dplyr::filter(scores, fill != "neutral"), size = 0.2, stat = "identity", position = "dodge") +
    ylab("G-score") +
    ggrepel::geom_text_repel(
      data = subset(scores, !is.na(name) & Type == "Amp"),
      nudge_y = max(subset(scores, !is.na(name) & Type == "Amp")$`G-score`) * adjust_amps,
      size = label_size,
      angle = 90,
      segment.size = 0.5,
      segment.color = "#000000",
      force_pull = force_pull,
      arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
      segment.curvature = segment.curvature,
      segment.ncp = segment.ncp,
      segment.angle = segment.angle
    ) +
    ggrepel::geom_text_repel(
      data = subset(scores, !is.na(name) & Type == "Del"),
      nudge_y = min(subset(scores, !is.na(name) & Type == "Del")$`G-score`) * adjust_dels,
      nudge_x = subset(scores, !is.na(name) & Type == "Del")$Start,
      size = label_size,
      angle = 90,
      segment.size = 0.5,
      segment.color = "#000000",
      force_pull = force_pull,
      arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
      segment.curvature = segment.curvature,
      segment.ncp = segment.ncp,
      segment.angle = segment.angle
    ) +
    facet_grid(. ~ Chromosome, scales = "free") +
    scale_color_manual(values = cnv_palette) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 16, colour = "black"),
      axis.ticks.x = element_blank(), axis.ticks.y = element_line(colour = "black"), legend.position = "none",
      panel.spacing.x = unit(0.1, "lines"), panel.border = element_blank(), text = element_text(size = 16, colour = "black"),
      strip.background = element_blank(), strip.text.x = element_blank(), panel.grid = element_blank()
    ) +
    # geom_hline(data=xses,aes(x=),yintercept = 0, size = 7) +
    geom_rect(data = xses, aes(ymin = -0.05, ymax = 0.05, xmin = 1, xmax = End), inherit.aes = FALSE, colour = "black") +
    geom_text(data = xses, aes(label = Chromosome, x = Mid, y = 0), size = 4, color = "white")
}