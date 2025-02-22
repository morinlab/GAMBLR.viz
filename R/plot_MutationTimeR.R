#' @title Plot MutationTimeR Results
#'
#' @description This function plots the timing information from MutationTimeR for SSMs and CNAs.
#'
#' @param this_sample_metadata Metadata with one row containing your sample details
#' @param projection Genome build projection (e.g., "hg38" or "grch37").
#' @param verbose Set to TRUE for a chattier experience. Default is FALSE.
#' @import ggplot2 dplyr ggpubr
#' @export
#' @examples
#'
#' my_meta = suppressMessages(get_gambl_metadata()) %>%
#'   dplyr::filter(sample_id=="01-20985T",
#'                seq_type=="genome")
#' timed = GAMBLR.results::get_timed_mutations(my_meta,"hg38")
#' print(timed$SSM)
#'
#' print(timed$CNA)
#'
#' all_plots = plot_MutationTimeR(my_meta,timed$CNA,timed$SSM)
#'
#' all_plots$full
#'
#' all_plots$minimal
plot_MutationTimeR <- function(this_sample_metadata,
                              timed_cna,
                              timed_ssm,
                              genome_build,
                              verbose = FALSE) {
  if (nrow(this_sample_metadata) > 1) {
    stop("this_sample_metadata must contain exactly one row")
  }
  if(missing(genome_build)){
    genome_build = get_genome_build(timed_cna)
  }
  normal_sample_id <- pull(this_sample_metadata, normal_sample_id)
  sample_id <- pull(this_sample_metadata, sample_id)
  tumour_sample_id <- sample_id

  if (verbose) {
    cat("Checking for CNAs with time info...\n")
    # Check for CNAs with time info
    if (dim(timed_cna %>% filter(!is.na(time)))[1] == 0) {
      cat("No timed CNAs found. Min plot will be empty.\n")
    } else {
      cat("Time info for CNAs found.\n")
    }
  }
  #get chr lengths
  if(missing(genome_build)){
    stop("You must provide a genome_build")
  }else if(genome_build == "grch37"){
    chr_levels = unique(GAMBLR.data::chromosome_arms_grch37$chromosome)
    length_df = GAMBLR.data::chromosome_arms_grch37 %>%
      group_by(chromosome) %>%
      dplyr::filter(chromosome !="Y") %>%
      summarize(length = max(end))
    length_df$chromosome = factor(length_df$chromosome, levels = chr_levels)
    length_df = arrange(length_df, chromosome)
  }else if(genome_build == "hg38"){
    chr_levels = unique(GAMBLR.data::chromosome_arms_hg38$chromosome)
    length_df = GAMBLR.data::chromosome_arms_hg38 %>%
      group_by(chromosome) %>%
      dplyr::filter(chromosome !="chrY") %>%
      summarize(length = max(end))
    length_df$chromosome = factor(length_df$chromosome, levels = chr_levels)
    length_df = arrange(length_df, chromosome)
  }else{
    stop("Unsupported genome_build")
  }
  chrom_lengths = pull(length_df,length)
  names(chrom_lengths) = pull(length_df,chromosome)
  # Plotting Functions
  plot_timed_SSM <- function(
      timed_ssm,
      timed_cna,
      all_ssm = FALSE,
      genome_build = "hg38",
      base_size = 12,
      point_size = 0.5,
      chrom_lengths) {
    s <- names(chrom_lengths)
    timed_ssm <- timed_ssm %>%
      mutate(
        VAF = t_alt_count / t_depth,
        start = Start_Position,
        end = End_Position,
        chr = Chromosome
      )
    timed_ssm$Chromosome <- factor(timed_ssm$Chromosome, levels = s)
    timed_ssm$cls <- factor(timed_ssm$CLS_time_label,
                           levels = c("clonal [early]",
                                      "clonal [late]",
                                      "clonal [NA]",
                                      "subclonal"))

    if (all_ssm) {
      p <- ggplot(data = timed_ssm, aes(x = Start_Position, y = VAF, color = cls)) +
        geom_point(alpha = 0.7, size = point_size, show.legend = TRUE) +
        facet_wrap(~Chromosome, scales = "free_x", nrow = 1) +
        ylim(c(0, 1)) +
        theme_Morons(base_size = base_size) +
        theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
        ylab("VAF") +
        scale_colour_manual(
          values = c("clonal [early]" = "#C77CFF",
                     "clonal [late]" = "#7CAE00",
                     "clonal [NA]" = "#00BFC4",
                     "subclonal" = "#F8766D"),
          drop = FALSE
        )
    } else {
      timed_cna$chr <- factor(timed_cna$chr, levels = s)
      just_timed <- timed_cna %>%
        filter(!is.na(time)) %>%
        mutate(start = startpos, end = endpos)
      #x <- data.table(timed_ssm)
      #y <- data.table(just_timed)
      #setkey(y, chr, start, end)
      print("RUNNING COOL_OVERLAPS")
      print(head(timed_ssm %>% arrange(Chromosome)))
      print(head(just_timed %>% arrange(chr)))
      just_timed = mutate(just_timed,start = as.numeric(start),
      end = as.numeric(end))
      timed_ssm = mutate(timed_ssm,Start_Position = as.numeric(Start_Position),
      End_Position = as.numeric(End_Position))

      kept_ssm <- cool_overlaps(timed_ssm,
                    just_timed,
                    columns1 = c("chr", "start", "end"),
                    columns2 = c("chr","start","end"),type="any") %>%
                    dplyr::filter(!is.na(startpos))


      p <- ggplot(data = kept_ssm,
                  aes(x = Start_Position,
                      y = VAF, color = cls)) +
        geom_point(alpha = 0.7, size = point_size, show.legend = TRUE) +
        facet_wrap(~Chromosome, scales = "free_x", nrow = 1) +
        ylim(c(0, 1)) +
        theme_Morons(base_size = base_size) +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank()) +
        ylab("VAF") +
        scale_colour_manual(
          values = c("clonal [early]" = "#C77CFF",
                     "clonal [late]" = "#7CAE00",
                     "clonal [NA]" = "#00BFC4",
                     "subclonal" = "#F8766D"),
          drop = FALSE
        )
    }

    return(p)
  }

  plot_timed_cna <- function(timed_cna,
                             genome_build = "hg38",
                             base_size = 12,
                             all_chrom = FALSE,
                             chrom_lengths) {
    s <- names(chrom_lengths)
    l <- as.numeric(chrom_lengths)
    cn_bounds <- data.frame(chr = s, start = 1, end = l)

    timed_cna$chr <- factor(timed_cna$chr, levels = s)
    cn_bounds$chr <- factor(cn_bounds$chr, levels = s)

    if (all_chrom) {
      if (dim(timed_cna %>% filter(!is.na(time)))[1] == 0) {
        p <- ggplot() +
          ggtitle("No CNAs were timed")
      } else {
        p <- ggplot(timed_cna) +
          geom_segment(aes(y = time, yend = time, x = startpos, xend = endpos), colour = "red") +
          geom_rect(aes(xmin = endpos, xmax = startpos, ymin = time.lo, ymax = time.up), alpha = 0.2) +
          geom_point(data = cn_bounds, aes(x = start, y = 1), colour = "white") +
          geom_point(data = cn_bounds, aes(x = end, y = 1), colour = "white") +
          ylim(c(0, 1)) +
          facet_wrap(~chr, scales = "free_x", nrow = 1) +
          theme_Morons(base_size = base_size) +
          theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
          ylab("Time")
      }
    } else {
      just_timed <- timed_cna %>% filter(!is.na(time))
      cn_bounds <- cn_bounds %>% filter(chr %in% just_timed$chr)

      p <- ggplot(just_timed) +
        geom_segment(aes(y = time, yend = time, x = startpos, xend = endpos), colour = "red") +
        geom_rect(aes(xmin = endpos, xmax = startpos, ymin = time.lo, ymax = time.up), alpha = 0.2) +
        geom_point(data = cn_bounds, aes(x = start, y = 1), colour = "white") +
        geom_point(data = cn_bounds, aes(x = end, y = 1), colour = "white") +
        ylim(c(0, 1)) +
        facet_wrap(~chr, scales = "free_x", nrow = 1) +
        theme_Morons(base_size = base_size) +
        theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
        ylab("Time")
    }

    return(p)
  }

  plot_total_CN <- function(timed_cna, genome_build = "hg38", base_size = 12, all_chrom = FALSE, chrom_lengths) {
    s <- names(chrom_lengths)
    l <- as.numeric(chrom_lengths)

    timed_cna$chr <- factor(timed_cna$chr, levels = s)

    timed_cna_total <- timed_cna %>%
      mutate(size = endpos - startpos + 1) %>%
      mutate(centre = startpos + size / 2) %>%
      mutate(CN = major_cn + 0.1, CN2 = minor_cn - 0.1)

    just_timed_chrom <- timed_cna_total %>%
      filter(!is.na(time)) %>%
      pull(chr)
    just_timed <- timed_cna_total %>%
      filter(chr %in% just_timed_chrom)

    if (all_chrom) {
      p <- ggplot(timed_cna_total) +
        geom_tile(aes(x = centre, y = CN, width = size, height = 0.2, alpha = clonal_frequency), fill = "lightblue") +
        geom_tile(aes(x = centre, y = CN2, width = size, height = 0.2, alpha = clonal_frequency), fill = "orange") +
        facet_wrap(~chr, scales = "free_x", nrow = 1) +
        ylim(c(-1, 4)) +
        theme_Morons(base_size = base_size) +
        theme(
          axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), axis.title.y = element_text(margin = margin(l = 16))
        )
    } else {
      p <- ggplot(just_timed) +
        geom_tile(aes(x = centre, y = CN, width = size, height = 0.2, alpha = clonal_frequency), fill = "lightblue") +
        geom_tile(aes(x = centre, y = CN2, width = size, height = 0.2, alpha = clonal_frequency), fill = "orange") +
        facet_wrap(~chr, scales = "free_x", nrow = 1) +
        ylim(c(-1, 4)) +
        theme_Morons(base_size = base_size) +
        theme(
          axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), axis.title.y = element_text(margin = margin(l = 16))
        )
    }

    return(p)
  }

  plot_all <- function(timed_ssm, timed_cna, genome_build = "hg38", sample_ids = "", title = "", base_size = 12, point_size = 0.8, all_ssm = FALSE, all_chrom = FALSE) {
    if (all_ssm & all_chrom) {
      CN_total <- plot_total_CN(timed_cna,
                               genome_build = genome_build,
                               base_size = base_size,
                               all_chrom = TRUE,
                               chrom_lengths = chrom_lengths)
      CN_timing <- plot_timed_cna(timed_cna,
                                 genome_build = genome_build,
                                 base_size = base_size,
                                 all_chrom = TRUE,
                                 chrom_lengths = chrom_lengths)
      SSM_timing <- plot_timed_SSM(timed_ssm,
                                  timed_cna,
                                  all_ssm = TRUE,
                                  point_size = point_size,
                                  genome_build = genome_build,
                                  base_size = base_size,
                                  chrom_lengths = chrom_lengths)
      p <- ggarrange(SSM_timing, CN_timing, CN_total, ncol = 1) %>%
        annotate_figure(., top = text_grob(title, face = "bold", size = 12), fig.lab = sample_ids, fig.lab.pos = "top.left")
    } else if (all_ssm != all_chrom) {
      stop("all_ssm and all_chrom are not the same. Set both to either FALSE or TRUE.")
    } else if (dim(timed_cna %>% filter(!is.na(time)))[1] == 0) {
      title <- "No CNAs were timed"
      p <- ggarrange(ggplot()) %>%
        annotate_figure(., top = text_grob(title, face = "bold", size = 12), fig.lab = sample_ids, fig.lab.pos = "top.left")
    } else {
      CN_total <- plot_total_CN(timed_cna, genome_build = genome_build, base_size = base_size, all_chrom = FALSE, chrom_lengths = chrom_lengths)
      CN_timing <- plot_timed_cna(timed_cna, genome_build = genome_build, base_size = base_size, all_chrom = FALSE, chrom_lengths = chrom_lengths)
      SSM_timing <- plot_timed_SSM(timed_ssm, timed_cna, all_ssm = FALSE, point_size = point_size, genome_build = genome_build, base_size = base_size, chrom_lengths = chrom_lengths)
      p <- ggarrange(SSM_timing, CN_timing, CN_total, ncol = 1) %>%
        annotate_figure(., top = text_grob(title, face = "bold", size = 12), fig.lab = sample_ids, fig.lab.pos = "top.left")
    }

    return(p)
  }

  # Creating and saving the plotting functions
  cat("Making full plot...\n")
  plot_full <- suppressWarnings(plot_all(timed_ssm, timed_cna,
    all_ssm = TRUE, all_chrom = TRUE, genome_build = genome_build,
    title = "All SSM and Copy Number States", sample_ids = paste0(tumour_sample_id, "--", normal_sample_id)
  ))

  cat("Making minimum plot...\n")
  plot_min <- suppressWarnings(plot_all(timed_ssm, timed_cna,
    all_ssm = FALSE, all_chrom = FALSE, genome_build = genome_build,
    title = "SSM and Copy Number in Timed CNA Regions", sample_ids = paste0(tumour_sample_id, "--", normal_sample_id)
  ))


  return(list(full = plot_full, minimal = plot_min))
}
