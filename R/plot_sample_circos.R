#' @title Sample-level Circos Plot
#'
#' @description Plot a sample-centric circos overview.
#'
#' @details This function takes a sample ID in the `this_sample_id` parameter
#' and draws one ring per genomic data type you provide: CNV segments
#' (`cnv_df`), SV breakpoint links (`sv_df`) and/or SSM point mutations
#' (`ssm_df`). This function does not fetch any of these on your behalf --
#' retrieve them yourself (e.g. with `get_sample_cn_segments()`,
#' `get_manta_sv()`, `get_ssm_by_sample()` or equivalent) and pass in
#' whichever data types you want plotted. A ring is included automatically
#' when its data frame is supplied; set the corresponding `include_*`
#' argument to `FALSE` to suppress a ring even when data was provided, or
#' to `TRUE` to get an explicit error if you forgot to pass it.
#'
#' @param this_sample_id Sample ID for the sample to plot.
#' @param sv_df Optional data frame of SVs for this sample (e.g. from
#' `get_manta_sv()`). Expected to have `tumour_sample_id`, `CHROM_A`,
#' `START_A`, `END_A`, `CHROM_B`, `START_B`, `END_B`. Rows not matching
#' `this_sample_id` are filtered out automatically.
#' @param cnv_df Optional data frame of CN segments for this sample (e.g.
#' from `get_sample_cn_segments()`). Expected to have `chrom`, `start`,
#' `end`, `log.ratio`.
#' @param ssm_df Optional data frame of SSMs for this sample in MAF-like
#' format (e.g. from `get_ssm_by_sample()`). Expected to have
#' `Tumor_Sample_Barcode`, `Chromosome`, `Start_Position`,
#' `Variant_Classification`. Rows not matching `this_sample_id` are
#' filtered out automatically.
#' @param include_sv Whether to draw the SV link ring. Defaults to `TRUE`
#' when `sv_df` is supplied and `FALSE` otherwise. Set explicitly to `TRUE`
#' to get a clear error if `sv_df` was left out by mistake.
#' @param include_cnv Whether to draw the CNV ring. Defaults to `TRUE` when
#' `cnv_df` is supplied and `FALSE` otherwise. Set explicitly to `TRUE` to
#' get a clear error if `cnv_df` was left out by mistake.
#' @param include_ssm Whether to draw the SSM point-mutation ring. Defaults
#' to `TRUE` when `ssm_df` is supplied and `FALSE` otherwise. Set explicitly
#' to `TRUE` to get a clear error if `ssm_df` was left out by mistake.
#' @param genome_build genome_build to use (this is usually inferred
#' automatically from your data).
#' @param this_seq_type Seq type to use when looking up this sample's
#' metadata for the legend (only relevant when `legend_metadata_columns` is
#' set). Default "genome".
#' @param legend_metadata_columns Column names from metadata
#' @param legend_metadata_names List of metadata names to be plotted.
#' @param chrom_list List of chromosomes to be plotted. If not stated, chr1-22+X will be used.
#' @param label_genes Gene labels (df, list or what type?)
#' @param colour_labels_by_cn Set to TRUE to colour each `label_genes`
#' entry's connector line (not the text itself, which stays black and
#' readable) by the copy number state at that gene's location: red if the
#' overlapping `cnv_df` segment's `log.ratio` is above 0 (gain), blue if
#' below 0 (loss). A flat two-colour rule rather than a continuous scale,
#' since intermediate colours along a gradient are hard to make out against
#' the plot's white background. Requires `include_cnv` (i.e. `cnv_df`). Genes with no
#' overlapping segment in `cnv_df` get a black connector line. Default FALSE.
#' @param auto_label_sv Default is FALSE. Requires `include_sv`.
#' @param link_colour_column Which column from sv_df to use for link colouring (default: NULL)
#' @param cn_colour_range Numeric vector of length 3, `c(low, mid, high)`
#' log ratio breakpoints for the blue-white-red CNV colour scale, shared by
#' the plot and its legend. Defaults to a symmetric range covering the
#' largest absolute log ratio observed in `cnv_df`.
#' @param hide_legend Set to TRUE if you want to suppress the legend. Particularly useful if you are not using GAMBL data/metadata
#' @param verbose Set to TRUE for a chattier experience
#'
#' @return Nothing
#'
#' @import dplyr circlize ComplexHeatmap ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' suppressPackageStartupMessages(library(GAMBLR.open))
#'
#' this_sample_cnv = get_sample_cn_segments(these_sample_ids = "02-13135T",
#'                                          with_chr_prefix = TRUE)
#' this_sample_sv = get_manta_sv() %>%
#'   dplyr::filter(tumour_sample_id == "02-13135T")
#'
#' plot_sample_circos(this_sample_id = "02-13135T",
#'                    cnv_df = this_sample_cnv,
#'                    sv_df = this_sample_sv,
#'                    legend_metadata_columns = c("pathology",
#'                                                "lymphgen",
#'                                                "COO_consensus",
#'                                                "DHITsig_consensus"),
#'                    legend_metadata_names = c("pathology",
#'                                              "LymphGen",
#'                                              "COO",
#'                                              "DHITsig"),
#'                    chrom_list = c("chr3",
#'                                   "chr8",
#'                                   "chr14",
#'                                   "chr18"))
#'}
plot_sample_circos = function(this_sample_id,
                              sv_df,
                              cnv_df,
                              ssm_df,
                              include_sv,
                              include_cnv,
                              include_ssm,
                              legend_metadata_columns,
                              legend_metadata_names = c(),
                              genome_build,
                              this_seq_type = "genome",
                              chrom_list,
                              label_genes,
                              colour_labels_by_cn = FALSE,
                              auto_label_sv = FALSE,
                              link_colour_column = NULL,
                              cn_colour_range,
                              hide_legend = FALSE,
                              verbose = FALSE){

  # a ring is included whenever its data was supplied; include_* lets the
  # caller force inclusion (-> a clear error if the data was forgotten) or
  # force exclusion even when data was passed
  if(missing(include_cnv)){
    include_cnv = !missing(cnv_df)
  } else if(include_cnv && missing(cnv_df)){
    stop("include_cnv = TRUE but cnv_df was not provided. Supply cnv_df yourself (e.g. from get_sample_cn_segments()).")
  }
  if(missing(include_sv)){
    include_sv = !missing(sv_df)
  } else if(include_sv && missing(sv_df)){
    stop("include_sv = TRUE but sv_df was not provided. Supply sv_df yourself (e.g. from get_manta_sv()).")
  }
  if(missing(include_ssm)){
    include_ssm = !missing(ssm_df)
  } else if(include_ssm && missing(ssm_df)){
    stop("include_ssm = TRUE but ssm_df was not provided. Supply ssm_df yourself (e.g. from get_ssm_by_sample()).")
  }
  if(auto_label_sv && !include_sv){
    stop("auto_label_sv = TRUE requires include_sv = TRUE (and sv_df to be supplied).")
  }
  if(colour_labels_by_cn && !include_cnv){
    stop("colour_labels_by_cn = TRUE requires include_cnv = TRUE (and cnv_df to be supplied).")
  }
  if(colour_labels_by_cn && missing(label_genes)){
    stop("colour_labels_by_cn = TRUE has nothing to colour without label_genes.")
  }

  add_cnv = function(cnv_df, col_fun){
    bed = data.frame(cnv_df[,c("chrom", "start", "end", "log.ratio")])
    colnames(bed) = c("chr", "start", "end", "value1")
    if(!any(grepl("chr",bed$chr))){
      bed = mutate(bed,chr = paste0("chr",chr))
    }
    circos.genomicTrackPlotRegion(bed, stack = TRUE, panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = col_fun(value), border = NA, posTransform = NULL, ...)
    }, bg.border = NA)
  }

  add_ssm = function(ssm_df, mutation_colours){
    bed = ssm_df %>%
      dplyr::transmute(
        chr = Chromosome,
        start = Start_Position,
        end = Start_Position,
        y = 0.5,
        colour = unname(mutation_colours[Variant_Classification])
      )
    bed$colour[is.na(bed$colour)] = "grey50"
    if(!any(grepl("chr", bed$chr))){
      bed = mutate(bed, chr = paste0("chr", chr))
    }
    circos.genomicTrackPlotRegion(bed, ylim = c(0, 1), panel.fun = function(region, value, ...) {
      circos.genomicPoints(region, value, numeric.column = "y", col = value$colour, pch = 16, cex = 0.4, ...)
    }, bg.border = NA, track.height = 0.08)
  }

  # one colour per row of gene_bed, from the log.ratio of whichever cnv_df
  # segment overlaps that gene's midpoint (first match if more than one);
  # genes with no overlapping segment stay black, same as the default label
  # colour when colour_labels_by_cn isn't used at all. col_fun is expected
  # to be a label-specific scale (see label_col_fun below), not the ring's
  # own cn_col_fun -- that one passes through white at 0, which is
  # unreadable as text on this plot's white background
  gene_cn_colours = function(gene_bed, cnv_df, col_fun, verbose = FALSE){
    # strip any "chr" prefix from both sides rather than assuming which one
    # has it (gene_bed is chr-prefixed for grch37 but not hg38; cnv_df's
    # prefix depends on how the caller built it) -- comparing bare
    # chromosome names avoids a silent, total match failure either way.
    # start/end are also coerced to numeric: if either came in as character
    # (e.g. read from a file without type conversion), "<=" would do a
    # lexicographic string comparison instead of a numeric one and silently
    # miss most/all overlaps.
    gene_chr = gsub("^chr", "", as.character(gene_bed$chromosome))
    cnv_chr = gsub("^chr", "", as.character(cnv_df$chrom))
    cnv_start = as.numeric(cnv_df$start)
    cnv_end = as.numeric(cnv_df$end)
    gene_mid = (as.numeric(gene_bed$start) + as.numeric(gene_bed$end)) / 2
    vapply(seq_len(nrow(gene_bed)), function(i){
      hits = which(cnv_chr == gene_chr[i] &
                     cnv_start <= gene_mid[i] &
                     cnv_end >= gene_mid[i])
      if(length(hits) == 0){
        if(verbose){
          message(paste0(
            "colour_labels_by_cn: no cnv_df segment overlaps ",
            gene_bed$gene_name[i], " (", gene_bed$chromosome[i], ":",
            round(gene_mid[i]), "); leaving its connector line black."
          ))
        }
        return("black")
      }
      col_fun(cnv_df$log.ratio[hits[1]])
    }, character(1))
  }

  genomic_data = list()
  if(include_cnv){
    genomic_data[["CN"]] = cnv_df
  }
  if(include_sv){
    genomic_data[["SV"]] = sv_df
  }
  if(include_ssm){
    genomic_data[["SSM"]] = ssm_df
  }
  if(length(genomic_data) == 0){
    stop("Nothing to plot: supply at least one of cnv_df, sv_df or ssm_df.")
  }

  if(missing(chrom_list)){
    #should we add chr list for males as well? Sex could then also be added as a paramter for the function were the appropiate chr list is called if not stated?
    chrom_list = paste0("chr", c(1:22,"X"))
  }

  # trust an explicit genome_build outright: check_get_projection()'s
  # cross-validation against each data frame's own genome_build attribute
  # only works for data fetched via GAMBLR wrappers that set it (e.g.
  # get_sample_cn_segments()); plain user-supplied data frames won't have
  # it and would otherwise always report a spurious mismatch
  if(!missing(genome_build)){
    this_projection = genome_build
  }else{
    this_projection = check_get_projection(genomic_data)
  }
  if(!is.character(this_projection) || length(this_projection) != 1 ||
       !this_projection %in% c("grch37", "hg38")){
    stop("Could not determine a genome build. Pass genome_build = \"grch37\" or \"hg38\" explicitly, or supply data frames that carry a genome_build attribute (e.g. from get_sample_cn_segments()/get_manta_sv()).")
  }
  if(verbose){
    message(paste("USING:",this_projection))
  }

  if(!missing(label_genes)){
    if(this_projection=="grch37"){
      gene_bed = GAMBLR.data::grch37_gene_coordinates %>%
        dplyr::filter(gene_name %in% label_genes) %>%
        dplyr::select(chromosome, start, end, gene_name)
    }else if(this_projection=="hg38"){
      gene_bed = GAMBLR.data::hg38_gene_coordinates %>%
        dplyr::filter(gene_name %in% label_genes) %>%
        dplyr::select(chromosome, start, end, gene_name)
    }
    # normalize chr-prefixing regardless of build (grch37_gene_coordinates
    # has bare chromosome names; hg38_gene_coordinates may or may not,
    # depending on source) before comparing against chrom_list below
    gene_bed$chromosome = ifelse(
      grepl("^chr", gene_bed$chromosome),
      gene_bed$chromosome,
      paste0("chr", gene_bed$chromosome)
    )
    # drop genes whose chromosome isn't part of this plot: circos.genomicLabels()
    # errors deep inside circlize ("'names' attribute ... must be the same
    # length as the vector") if asked to place a label on an uninitialized
    # sector, so catch it here with a clear, actionable warning instead
    outside_chrom_list = !gene_bed$chromosome %in% chrom_list
    if(any(outside_chrom_list)){
      warning(paste0(
        "label_genes: dropping gene(s) not on a chromosome in chrom_list (",
        paste(gene_bed$gene_name[outside_chrom_list], collapse = ", "),
        "). Add their chromosome to chrom_list to include them."
      ))
      gene_bed = gene_bed[!outside_chrom_list, , drop = FALSE]
    }
  }

  if(include_cnv){
    if(missing(cn_colour_range)){
      max_abs = suppressWarnings(max(abs(cnv_df$log.ratio), na.rm = TRUE))
      if(!is.finite(max_abs) || max_abs == 0){
        max_abs = 1
      }
      cn_colour_range = c(-max_abs, 0, max_abs)
    }
    cn_col_fun = circlize::colorRamp2(cn_colour_range, c("blue", "white", "red"))
  }

  if(colour_labels_by_cn && nrow(gene_bed) > 0){
    # a two-colour (gain/loss) rule rather than reusing the ring's
    # continuous scale: a connector line sits directly on the plot's white
    # background, where intermediate colours along a diverging gradient
    # (including near-white or near-black ones close to 0) are hard to
    # read, whereas a flat red/blue split stays unambiguous at a glance
    label_col_fun = function(x) ifelse(x > 0, "red", ifelse(x < 0, "blue", "black"))
    label_colours = gene_cn_colours(gene_bed, cnv_df, label_col_fun, verbose = verbose)
  }

  if(include_ssm){
    ssm_df = ssm_df %>%
      dplyr::filter(Tumor_Sample_Barcode == this_sample_id)
    mutation_colours = GAMBLR.helpers::get_gambl_colours("mutation")
  }

  if(include_sv){
    sv_df = sv_df %>%
      dplyr::filter(tumour_sample_id == this_sample_id)

    #add chr prefixes if grch37 is selcted (expected by circlize)
    if(this_projection == "grch37"){
      sv_df = sv_df %>%
        dplyr::mutate(CHROM_A = paste0("chr", CHROM_A)) %>%
        dplyr::mutate(CHROM_B = paste0("chr", CHROM_B))
    }
    sv_df = sv_df %>%
      dplyr::filter(CHROM_A %in% chrom_list) %>%
      dplyr::filter(CHROM_B %in% chrom_list)

    if(auto_label_sv){
      #annotate oncogene SVs and label them
      annotated_sv  = annotate_sv(sv_df, with_chr_prefix = TRUE,genome_build = this_projection) %>%
        dplyr::filter(!is.na(partner)) %>%
        dplyr::filter(tumour_sample_id == this_sample_id)

      these_oncogenes = unique(pull(annotated_sv, gene))
      these_partners = unique(pull(annotated_sv, partner))
      if("IGH" %in% these_partners){
        these_partners = c(these_partners, "IGHV3-62")
      }
      anno_bed1 = annotated_sv %>%
        dplyr::select(chrom1, start1, end1, tumour_sample_id)

      anno_bed2 = annotated_sv %>%
        dplyr::select(chrom2, start2, end2, tumour_sample_id)

      colnames(anno_bed1) = c("chrom", "start", "end", "sample_id")
      colnames(anno_bed2) = c("chrom", "start", "end", "sample_id")

      if(this_projection=="grch37"){
        bed_mut_partner = GAMBLR.data::grch37_partners %>%
          dplyr::filter(gene %in% these_partners) %>%
          mutate(chrom = paste0("chr", chrom))

        bed_mut_onco = GAMBLR.data::grch37_oncogene %>%
          dplyr::filter(gene %in% these_oncogenes) %>%
          mutate(chrom = paste0("chr", chrom))
      }else if(this_projection=="hg38"){
        bed_mut_partner = GAMBLR.data::hg38_partners %>%
          dplyr::filter(gene %in% these_partners) %>%
          mutate(chrom = paste0("chr", chrom))

        bed_mut_onco = GAMBLR.data::hg38_oncogene %>%
          dplyr::filter(gene %in% these_oncogenes) %>%
          mutate(chrom = paste0("chr", chrom))
      }

      bed_mut = bind_rows(bed_mut_partner, bed_mut_onco)
      if(verbose){
        print(bed_mut)
      }
    }
    bed1 = sv_df %>%
      dplyr::select(CHROM_A, START_A, END_A, tumour_sample_id)

    bed2 = sv_df %>%
      dplyr::select(CHROM_B, START_B, END_B, tumour_sample_id)

    colnames(bed1) = c("chrom", "start", "end", "sample_id")
    colnames(bed2) = c("chrom", "start", "end", "sample_id")
  }

  circos.clear()
  circos.initializeWithIdeogram(chromosome.index = chrom_list)
  if(include_cnv){
    add_cnv(cnv_df, cn_col_fun)
  }
  if(include_ssm){
    add_ssm(ssm_df, mutation_colours)
  }

  if(include_sv){
    if(is.null(link_colour_column)){
      circos.genomicLink(bed1, bed2, col = "#bdbdc1")
    } else if(link_colour_column == "CHROM_B"){
      bed2 = decorate_bed(bed2,colour_by="chrom",colour_mapping = get_gambl_colours("chromosome"))
      circos.genomicLink(bed1, bed2, col = bed2$color)
    }else if(link_colour_column == "CHROM_A"){
      bed1 = decorate_bed(bed1,colour_by="chrom",colour_mapping = get_gambl_colours("chromosome"))
      circos.genomicLink(bed1, bed2, col = bed1$color)
    }else{
      #assume values in link_colour_column are already colours
      if(!link_colour_column %in% colnames(sv_df)){
        stop(paste(link_colour_column,"is not an existing column in your sv data frame"))
      }
      if(verbose){
        message(paste("colour values:"))
        message(paste(unique(sv_df[,link_colour_column]), collapse = ", "))
      }
      circos.genomicLink(bed1, bed2, col = sv_df[,link_colour_column])
    }
    if(auto_label_sv){
      circos.genomicLink(anno_bed1, anno_bed2,col = 'red')
      circos.genomicLabels(bed_mut, labels.column = "gene")
    }
  }

  if(!missing(label_genes) && nrow(gene_bed) > 0){
    if(colour_labels_by_cn){
      circos.genomicLabels(gene_bed, labels.column = "gene_name", line_col = label_colours)
    }else{
      circos.genomicLabels(gene_bed, labels.column = "gene_name")
    }
  }
  text(0, 0, this_sample_id, cex = 0.8)
  if(hide_legend){
    return(invisible(NULL))
  }
  if(!missing(legend_metadata_columns)){
    samp_meta = GAMBLR.helpers::handle_metadata(this_seq_type = this_seq_type) %>%
      dplyr::filter(sample_id == this_sample_id)

    these_meta = samp_meta[legend_metadata_columns]
    these_meta <- mutate_if(these_meta, is.factor, as.character)
    vals = as.character(these_meta)

    all_cols = map_metadata_to_colours(legend_metadata_columns, these_meta, verbose = verbose)

    cols = all_cols[vals]
    if(verbose){
      print(cols)
    }
    if(length(legend_metadata_names) == length(legend_metadata_columns)){
      for(i in c(1:length(vals))){
        if(!legend_metadata_names[i] == ""){
          vals[i] = paste(legend_metadata_names[i], vals[i])
        }
      }
    }

    lgd_discrete = Legend(labels = vals, title_position = "topleft", legend_gp = gpar(fill = cols))
    draw(lgd_discrete, x = unit(3, "mm"), y = unit(1, "npc") - unit(5, "mm"), just = c("left", "top"))
  }

  if(include_cnv){
    legend_breaks = round(seq(cn_colour_range[1], cn_colour_range[3], length.out = 5), 2)
    lgd_cnv = Legend(at = legend_breaks,
                     col_fun = cn_col_fun,
                     title_position = "topleft",
                     title = "log\nratio")

    draw(lgd_cnv, x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))
  }
}
