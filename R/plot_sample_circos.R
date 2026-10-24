#' @title Sample-level Circos Plot
#'
#' @description Plot a sample-centric circos overview.
#'
#' @details This function takes a sample ID in the `this_sample_id` parameter
#' and draws one ring per genomic data type you provide: CNV segments
#' (`cnv_df`), SV breakpoint links (`sv_df`) and/or SSM point mutations
#' (`ssm_df`). This function does not fetch any of these on your behalf --
#' retrieve them yourself (e.g. with `get_cn_segments()`, `get_combined_sv()`
#' or `get_manta_sv()`, `get_ssm_by_samples()` or equivalent) and pass in
#' whichever data types you want plotted. A ring is included automatically
#' when its data frame is supplied; set the corresponding `include_*`
#' argument to `FALSE` to suppress a ring even when data was provided, or
#' to `TRUE` to get an explicit error if you forgot to pass it.
#'
#' Only this sample's data, on the chromosomes in `chrom_list`, are used:
#' rows of `cnv_df` (by ID), `sv_df` (by tumour_sample_id) and `ssm_df` (by
#' Tumor_Sample_Barcode) for other samples or chromosomes are dropped, with a
#' message. Chromosome names may have a "chr" prefix or not.
#'
#' \strong{Oncogene SVs.} With `sv_df`, SVs involving known lymphoma
#' oncogenes are found with `GAMBLR.utils::annotate_sv()`. By default each
#' oncogene's SV links are drawn in their own colour on top of the other SVs
#' (grey), with a legend (`highlight_oncogene_sv`), and labelled at their
#' breakpoints: the oncogene at its end and the partner (e.g. IGH, or a
#' partner gene; when known) at the other (`label_oncogene_sv`). `oncogenes`
#' limits this to some genes.
#'
#' @param this_sample_id Sample ID for the sample to plot.
#' @param sv_df Optional data frame of SVs (e.g. from `get_combined_sv()`
#' or `get_manta_sv()`), with `tumour_sample_id`, `CHROM_A`, `START_A`,
#' `END_A`, `CHROM_B`, `START_B`, `END_B` (and `STRAND_A`, `STRAND_B` for
#' oncogene annotation).
#' @param cnv_df Optional data frame of CN segments (e.g. from
#' `get_cn_segments()`), with `chrom`, `start`, `end`, `CN` and/or
#' `log.ratio` (see `cn_value`), and `ID` if it holds more than one sample.
#' @param ssm_df Optional data frame of SSMs in MAF-like format (e.g. from
#' `get_ssm_by_samples()`), with `Tumor_Sample_Barcode`, `Chromosome`,
#' `Start_Position`, `Variant_Classification`.
#' @param include_sv Whether to draw the SV link ring. Defaults to `TRUE`
#' when `sv_df` is supplied and `FALSE` otherwise. Set explicitly to `TRUE`
#' to get a clear error if `sv_df` was left out by mistake.
#' @param include_cnv Whether to draw the CNV ring. Defaults to `TRUE` when
#' `cnv_df` is supplied and `FALSE` otherwise. Set explicitly to `TRUE` to
#' get a clear error if `cnv_df` was left out by mistake.
#' @param include_ssm Whether to draw the SSM point-mutation ring. Defaults
#' to `TRUE` when `ssm_df` is supplied and `FALSE` otherwise. Set explicitly
#' to `TRUE` to get a clear error if `ssm_df` was left out by mistake.
#' @param these_samples_metadata Optional metadata containing this sample,
#' for the legend (`legend_metadata_columns`).
#' @param legend_metadata_columns Columns of `these_samples_metadata` to
#' show in the legend.
#' @param legend_metadata_names Optional names for those columns in the
#' legend.
#' @param genome_build genome_build to use (this is usually inferred
#' automatically from your data).
#' @param this_seq_type When `these_samples_metadata` has more than one row
#' for this sample, the seq_type of the row to use. Default "genome".
#' @param chrom_list Chromosomes to plot. Default chr1-22 and X.
#' @param label_genes Optional gene symbols to label at their locations.
#' @param colour_labels_by_cn Set to TRUE to colour each `label_genes`
#' entry's connector line (not the text itself, which stays black and
#' readable) by the copy number state at that gene's location: red if the
#' overlapping `cnv_df` segment's value (`cn_value`) is above neutral (2
#' for CN, 0 for log.ratio), blue if below. Requires `cnv_df`. Genes with no overlapping segment get
#' a black connector line. Default FALSE.
#' @param highlight_oncogene_sv Colour SVs involving oncogenes, one colour
#' per oncogene, on top of the other SVs, with a legend. Default TRUE (with
#' `sv_df`).
#' @param label_oncogene_sv Label oncogene SVs at their breakpoints: the
#' oncogene and its partner. Default TRUE (with `sv_df`).
#' @param oncogenes Optional oncogene symbols to limit the highlighting and
#' labels to (e.g. c("MYC", "BCL2", "BCL6")). Default: all found.
#' @param auto_label_sv Deprecated: the same as `highlight_oncogene_sv` and
#' `label_oncogene_sv`.
#' @param cnv_features Optional CNV features to call in this sample and
#' label: a feature table (e.g. from `GAMBLR.utils::get_cnv_features()`) or
#' TRUE for the default curated set. The sample's segments are called with
#' `GAMBLR.utils::segmented_data_to_cnv_status()` on the values plotted
#' (`cn_value`'s scale, with `centering = "none"`, so a feature is called
#' against the ring's neutral level: centre the segments first, e.g. with
#' `GAMBLR.utils::center_segments()`, for a sample-specific baseline) and
#' each feature called in its direction is labelled with its name, its
#' connector red for gains and blue for losses (darker for high-level
#' calls). Requires `cnv_df`.
#' @param cnv_feature_args Optional list of further arguments for
#' `segmented_data_to_cnv_status()` (e.g. thresholds, centering or a
#' blacklist).
#' @param link_colour_column Which column from sv_df to use for colouring
#' the other SV links: "CHROM_A" or "CHROM_B" (by chromosome) or a column of
#' colours. Default: grey.
#' @param cn_value The `cnv_df` column the CNV ring (and `colour_labels_by_cn`)
#' shows: "CN" (default; absolute copy number, e.g. centred with
#' `GAMBLR.utils::center_segments` so a sample's neutral level is 2) or
#' "log.ratio".
#' @param show_loh Draw regions with LOH (`LOH_flag` of 1 in `cnv_df`, as
#' reported for Battenberg and controlfreec segments) as a thin ring just
#' inside the CN ring, coloured by the copy number change with the LOH:
#' copy-neutral LOH, deletion with LOH or gain with LOH (see `loh_types`;
#' copy number against
#' 2, so centre high-ploidy samples first; see `cn_value`). Default: when
#' `cnv_df` has LOH_flag values for this sample. Samples without (e.g.
#' capture data) are noted as "not assessed".
#' @param loh_types The kinds of LOH drawn: any of "copy-neutral"
#' (default; grey), "deletion" (deletion with LOH) and "gain" (gain with
#' LOH).
#' @param label_cnloh_features Label recurrently deleted regions (the loss
#' features of `cnv_features`, or of the core set from
#' `GAMBLR.utils::get_cnv_features()` without it) covered by copy-neutral
#' LOH in this sample, by at least the feature's `min_fraction`: e.g. "TP73
#' cnLOH" (the representative gene, or the region without one), with grey
#' connectors. Default FALSE.
#' @param loh_track_height Height of the LOH ring, as a fraction of the
#' circle's radius. Default 0.03.
#' @param cn_colour_range Numeric vector of length 3, `c(low, mid, high)`
#' breakpoints for the blue-white-red CNV colour scale (values beyond them
#' take the end colours), shared by the plot and its legend. Default:
#' `c(0, 2, 6)` for CN; for log.ratio a symmetric range covering the
#' largest absolute log ratio in this sample's `cnv_df`.
#' @param sample_id_position Where to write the sample ID: "topleft"
#' (default; above the metadata legend), "centre" (inside the ring, where
#' SV links can make it hard to read) or "none".
#' @param hide_legend Set to TRUE to suppress the legends.
#' @param return_data Set to TRUE to return (invisibly) the data plotted.
#' @param verbose Set to TRUE for a chattier experience
#'
#' @return Nothing, or with `return_data = TRUE` (invisibly) a list:
#' `cnv`, `sv` and `ssm` (this sample's rows plotted), `oncogene_sv`
#' (the annotated oncogene SVs, with which end is the oncogene),
#' `cnv_features` (the CNV features called and labelled), `loh` (the LOH
#' segments drawn, with `loh_type`) and `cnloh_features` (the loss features
#' labelled as under copy-neutral LOH, with the fraction covered).
#'
#' @import dplyr circlize ComplexHeatmap ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' meta = get_gambl_metadata() %>% dplyr::filter(sample_id == "00-14595_tumorC",
#'                                               seq_type == "genome")
#' plot_sample_circos(this_sample_id = "00-14595_tumorC",
#'                    cnv_df = get_cn_segments(meta),
#'                    sv_df = get_combined_sv(these_samples_metadata = meta),
#'                    these_samples_metadata = meta,
#'                    legend_metadata_columns = c("pathology", "lymphgen",
#'                                                "COO_consensus"),
#'                    legend_metadata_names = c("pathology", "LymphGen", "COO"))
#'}
plot_sample_circos = function(this_sample_id,
                              sv_df,
                              cnv_df,
                              ssm_df,
                              include_sv,
                              include_cnv,
                              include_ssm,
                              these_samples_metadata = NULL,
                              legend_metadata_columns,
                              legend_metadata_names = c(),
                              genome_build,
                              this_seq_type = "genome",
                              chrom_list,
                              label_genes,
                              colour_labels_by_cn = FALSE,
                              highlight_oncogene_sv = TRUE,
                              label_oncogene_sv = TRUE,
                              oncogenes = NULL,
                              auto_label_sv = NULL,
                              cnv_features = NULL,
                              cnv_feature_args = list(),
                              link_colour_column = NULL,
                              cn_value = "CN",
                              show_loh = NULL,
                              loh_types = "copy-neutral",
                              loh_track_height = 0.03,
                              label_cnloh_features = FALSE,
                              cn_colour_range,
                              sample_id_position = "topleft",
                              hide_legend = FALSE,
                              return_data = FALSE,
                              verbose = FALSE){

  # a ring is included whenever its data was supplied; include_* lets the
  # caller force inclusion (-> a clear error if the data was forgotten) or
  # force exclusion even when data was passed
  if(missing(include_cnv)){
    include_cnv = !missing(cnv_df)
  } else if(include_cnv && missing(cnv_df)){
    stop("include_cnv = TRUE but cnv_df was not provided. Supply cnv_df yourself (e.g. from get_cn_segments()).")
  }
  if(missing(include_sv)){
    include_sv = !missing(sv_df)
  } else if(include_sv && missing(sv_df)){
    stop("include_sv = TRUE but sv_df was not provided. Supply sv_df yourself (e.g. from get_combined_sv()).")
  }
  if(missing(include_ssm)){
    include_ssm = !missing(ssm_df)
  } else if(include_ssm && missing(ssm_df)){
    stop("include_ssm = TRUE but ssm_df was not provided. Supply ssm_df yourself (e.g. from get_ssm_by_samples()).")
  }
  if(!is.null(auto_label_sv)){
    highlight_oncogene_sv = auto_label_sv
    label_oncogene_sv = auto_label_sv
  }
  if(colour_labels_by_cn && !include_cnv){
    stop("colour_labels_by_cn = TRUE requires include_cnv = TRUE (and cnv_df to be supplied).")
  }
  if(colour_labels_by_cn && missing(label_genes)){
    stop("colour_labels_by_cn = TRUE has nothing to colour without label_genes.")
  }
  if(!missing(legend_metadata_columns) && !hide_legend){
    if(is.null(these_samples_metadata)){
      stop("legend_metadata_columns needs these_samples_metadata (with this sample's row); ",
           "this function doesn't fetch metadata")
    }
    if(!this_sample_id %in% these_samples_metadata$sample_id){
      stop(this_sample_id, " is not in these_samples_metadata")
    }
  }

  if(missing(chrom_list)){
    chrom_list = paste0("chr", c(1:22, "X"))
  }
  chrom_list = paste0("chr", gsub("^chr", "", chrom_list))
  # (paste0 would turn an empty vector into "chr")
  with_chr = function(x) if(length(x) == 0) character(0) else paste0("chr", gsub("^chr", "", as.character(x)))

  # this sample's rows on the chromosomes plotted
  only_this = function(df, id_col, chrom_cols, what){
    df = as.data.frame(df)
    n = nrow(df)
    if(id_col %in% colnames(df)){
      df = df[as.character(df[[id_col]]) %in% this_sample_id, , drop = FALSE]
    }
    for(col in chrom_cols) df = df[with_chr(df[[col]]) %in% chrom_list, , drop = FALSE]
    if(nrow(df) < n && verbose){
      message(what, ": using ", nrow(df), " of ", n, " rows (this sample, chromosomes plotted)")
    }
    if(id_col %in% colnames(df) && nrow(df) == 0){
      warning("no ", what, " rows for ", this_sample_id, " on the chromosomes plotted")
    }
    df
  }
  if(include_cnv) cnv_df = only_this(cnv_df, "ID", "chrom", "cnv_df")
  if(include_sv) sv_df = only_this(sv_df, "tumour_sample_id", c("CHROM_A", "CHROM_B"), "sv_df")
  if(include_ssm) ssm_df = only_this(ssm_df, "Tumor_Sample_Barcode", "Chromosome", "ssm_df")

  genomic_data = list()
  if(include_cnv) genomic_data[["CN"]] = cnv_df
  if(include_sv) genomic_data[["SV"]] = sv_df
  if(include_ssm) genomic_data[["SSM"]] = ssm_df
  if(length(genomic_data) == 0){
    stop("Nothing to plot: supply at least one of cnv_df, sv_df or ssm_df.")
  }

  # trust an explicit genome_build outright: check_get_projection()'s
  # cross-validation against each data frame's own genome_build attribute
  # only works for data fetched via GAMBLR wrappers that set it
  if(!missing(genome_build)){
    this_projection = genome_build
  }else{
    this_projection = check_get_projection(genomic_data)
  }
  if(!is.character(this_projection) || length(this_projection) != 1 ||
       !this_projection %in% c("grch37", "hg38")){
    stop("Could not determine a genome build. Pass genome_build = \"grch37\" or \"hg38\" explicitly, or supply data frames that carry a genome_build attribute (e.g. from get_cn_segments()/get_manta_sv()).")
  }
  if(verbose){
    message(paste("USING:", this_projection))
  }

  if(!missing(label_genes)){
    gene_coords = if(this_projection == "grch37") GAMBLR.data::grch37_gene_coordinates
                  else GAMBLR.data::hg38_gene_coordinates
    gene_bed = gene_coords %>%
      dplyr::filter(gene_name %in% label_genes) %>%
      dplyr::select(chromosome, start, end, gene_name) %>%
      as.data.frame()
    gene_bed$chromosome = with_chr(gene_bed$chromosome)
    # circos.genomicLabels() errors deep inside circlize when asked to label
    # a chromosome that isn't plotted, so leave those out with a warning
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
    cn_value = match.arg(cn_value, c("CN", "log.ratio"))
    if(!cn_value %in% colnames(cnv_df)){
      stop("cnv_df has no ", cn_value, " column; use cn_value = \"",
           setdiff(c("CN", "log.ratio"), cn_value), "\"")
    }
    cnv_df$.value = as.numeric(cnv_df[[cn_value]])
    # neutral: 2 copies, or a log ratio of 0
    neutral = if(cn_value == "CN") 2 else 0
    if(missing(cn_colour_range)){
      if(cn_value == "CN"){
        cn_colour_range = c(0, 2, 6)
      }else{
        max_abs = suppressWarnings(max(abs(cnv_df$.value), na.rm = TRUE))
        if(!is.finite(max_abs) || max_abs == 0){
          max_abs = 1
        }
        cn_colour_range = c(-max_abs, 0, max_abs)
      }
    }
    cn_col_fun = circlize::colorRamp2(cn_colour_range, c("blue", "white", "red"))
  }
  # LOH segments (LOH_flag 1) typed by the copy number change with the LOH:
  # copy number against 2 (rounded), or without CN the log ratio against
  # +/-0.2
  loh_colours = c("copy-neutral LOH" = "grey40", "deletion with LOH" = "#92C5DE", "gain with LOH" = "#B2182B")
  loh_types = match.arg(loh_types, c("copy-neutral", "deletion", "gain"), several.ok = TRUE)
  loh_shown = paste0(c("copy-neutral" = "copy-neutral LOH", "deletion" = "deletion with LOH",
                       "gain" = "gain with LOH")[loh_types])
  loh_assessed = include_cnv && "LOH_flag" %in% colnames(cnv_df) && any(!is.na(cnv_df$LOH_flag))
  if(is.null(show_loh)) show_loh = loh_assessed
  loh_typed = NULL
  if(loh_assessed){
    loh = cnv_df[cnv_df$LOH_flag %in% 1, , drop = FALSE]
    if(nrow(loh) > 0){
      change = if("CN" %in% colnames(loh) && !all(is.na(loh$CN))) sign(round(as.numeric(loh$CN)) - 2)
               else ifelse(as.numeric(loh$log.ratio) > 0.2, 1, ifelse(as.numeric(loh$log.ratio) < -0.2, -1, 0))
      loh$loh_type = c("deletion with LOH", "copy-neutral LOH", "gain with LOH")[change + 2]
      loh_typed = loh
    }
  }

  if(colour_labels_by_cn && nrow(gene_bed) > 0){
    # a two-colour (gain/loss) rule: intermediate colours of the ring's
    # diverging scale are unreadable as lines on the white background
    label_col_fun = function(x) ifelse(x > neutral, "red", ifelse(x < neutral, "blue", "black"))
    label_colours = .circos_gene_cn_colours(gene_bed, cnv_df, label_col_fun, verbose = verbose)
  }

  # oncogene SVs: annotated on the unprefixed coordinates, with the end in
  # (or nearest) the oncogene identified
  oncogene_sv = NULL
  if(include_sv && (highlight_oncogene_sv || label_oncogene_sv) && nrow(sv_df) > 0){
    oncogene_sv = .circos_oncogene_sv(sv_df, this_projection, oncogenes)
    if(!is.null(oncogene_sv)){
      oncogene_sv = oncogene_sv[with_chr(oncogene_sv$onco_chrom) %in% chrom_list &
                                  with_chr(oncogene_sv$partner_chrom) %in% chrom_list, , drop = FALSE]
      if(nrow(oncogene_sv) == 0) oncogene_sv = NULL
    }
    if(verbose){
      message(if(is.null(oncogene_sv)) 0 else nrow(oncogene_sv), " oncogene SV(s)",
              if(!is.null(oncogene_sv)) paste0(": ", paste(unique(oncogene_sv$gene), collapse = ", ")) else "")
    }
  }

  if(include_ssm){
    mutation_colours = GAMBLR.helpers::get_gambl_colours("mutation")
  }

  # one colour per oncogene with SVs
  onco_colours = NULL
  if(highlight_oncogene_sv && !is.null(oncogene_sv)){
    genes = unique(oncogene_sv$gene)
    palette = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF",
                "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E")
    onco_colours = stats::setNames(rep(palette, length.out = length(genes)), genes)
  }
  # labels: oncogene SV breakpoints (oncogene and partner) and label_genes,
  # in one call so they don't overlap
  labels = NULL
  if(label_oncogene_sv && !is.null(oncogene_sv)){
    # one colour per label (also when there are none, e.g. no named partners)
    label_colour = function(genes){
      if(is.null(onco_colours)) rep("black", length(genes)) else unname(onco_colours[genes])
    }
    onco_lab = unique(data.frame(chr = with_chr(oncogene_sv$onco_chrom), start = oncogene_sv$onco_pos,
                                 end = oncogene_sv$onco_pos + 1, label = oncogene_sv$gene,
                                 colour = label_colour(oncogene_sv$gene)))
    partner_lab = oncogene_sv[!is.na(oncogene_sv$partner), , drop = FALSE]
    partner_lab = unique(data.frame(chr = with_chr(partner_lab$partner_chrom), start = partner_lab$partner_pos,
                                    end = partner_lab$partner_pos + 1, label = partner_lab$partner,
                                    colour = label_colour(partner_lab$gene)))
    # one label per gene and chromosome (several breakpoints of one event)
    lab = rbind(onco_lab, partner_lab)
    labels = lab[!duplicated(paste(lab$label, lab$chr)), , drop = FALSE]
  }
  if(!missing(label_genes) && nrow(gene_bed) > 0){
    gl = data.frame(chr = gene_bed$chromosome, start = gene_bed$start, end = gene_bed$end,
                    label = gene_bed$gene_name,
                    colour = if(colour_labels_by_cn) label_colours else "black")
    labels = rbind(labels, gl[!gl$label %in% labels$label, , drop = FALSE])
  }
  # CNV features called in this sample (GAMBLR.utils::segmented_data_to_cnv_status),
  # labelled in their direction: red for gains, blue for losses, darker for
  # high-level calls
  called_features = NULL
  if(!is.null(cnv_features) && !isFALSE(cnv_features)){
    if(!include_cnv) stop("cnv_features needs cnv_df")
    seg_one = cnv_df
    if(!"ID" %in% colnames(seg_one)) seg_one$ID = this_sample_id
    # called on the values plotted, without re-centring (centering "none"),
    # so the labels describe what the ring shows; cnv_feature_args can
    # override any of these
    status_args = utils::modifyList(
      c(list(seg_data = seg_one, genome_build = this_projection,
             value_scale = if(cn_value == "CN") "CN" else "log_ratio",
             centering = "none", output_format = "long"),
        if(is.data.frame(cnv_features)) list(cnv_features = cnv_features)),
      cnv_feature_args)
    called_features = tryCatch(do.call(GAMBLR.utils::segmented_data_to_cnv_status, status_args),
                               error = function(e){ warning("CNV features couldn't be called: ",
                                                            conditionMessage(e)); NULL })
    if(!is.null(called_features)){
      called_features = as.data.frame(called_features)
      called_features = called_features[called_features$level > 0 &
                                          with_chr(called_features$chrom) %in% chrom_list, , drop = FALSE]
      if(nrow(called_features) == 0) called_features = NULL
    }
    if(verbose) message(if(is.null(called_features)) 0 else nrow(called_features), " CNV feature(s) called")
    if(!is.null(called_features)){
      loss = grepl("loss|del", called_features$name, ignore.case = TRUE)
      high = called_features$level >= 2
      fl = data.frame(chr = with_chr(called_features$chrom), start = called_features$start,
                      end = called_features$end, label = called_features$name,
                      colour = ifelse(loss, ifelse(high, "#08306B", "#4393C3"),
                                      ifelse(high, "#67001F", "#D6604D")))
      labels = rbind(labels, fl[!fl$label %in% labels$label, , drop = FALSE])
    }
  }
  # recurrently deleted regions (loss features) under copy-neutral LOH,
  # labelled by their representative gene (or region)
  cnloh_features = NULL
  if(label_cnloh_features){
    if(!loh_assessed){
      message("label_cnloh_features: no LOH_flag values for this sample, so nothing to label")
    }else{
      # the features given, or the core set
      feats = if(is.data.frame(cnv_features)) as.data.frame(cnv_features)
              else as.data.frame(GAMBLR.utils::get_cnv_features(this_projection))
      feats = feats[feats$direction %in% "loss", , drop = FALSE]
      feats$chrom = gsub("^chr", "", as.character(feats$chrom))
      cn = if(is.null(loh_typed)) NULL else loh_typed[loh_typed$loh_type == "copy-neutral LOH", , drop = FALSE]
      if(!is.null(cn) && nrow(cn) > 0 && nrow(feats) > 0){
        cn_chrom = gsub("^chr", "", as.character(cn$chrom))
        covered = vapply(seq_len(nrow(feats)), function(i){
          ov = pmin(cn$end, feats$end[i]) - pmax(cn$start, feats$start[i])
          sum(ov[cn_chrom == feats$chrom[i] & ov > 0]) / (feats$end[i] - feats$start[i])
        }, numeric(1))
        need = if("min_fraction" %in% colnames(feats)) ifelse(is.na(feats$min_fraction), 0.2, feats$min_fraction) else 0.2
        hit = feats[covered >= need & with_chr(feats$chrom) %in% chrom_list, , drop = FALSE]
        if(nrow(hit) > 0){
          what = if("representative_gene" %in% colnames(hit)) hit$representative_gene else rep(NA, nrow(hit))
          what[is.na(what)] = if("region_name" %in% colnames(hit)) hit$region_name[is.na(what)] else hit$name[is.na(what)]
          cnloh_features = data.frame(name = hit$name, label = paste(what, "cnLOH"), chrom = hit$chrom,
                                      start = hit$start, end = hit$end,
                                      fraction_cnloh = round(covered[covered >= need & with_chr(feats$chrom) %in% chrom_list], 3))
          fl = data.frame(chr = with_chr(cnloh_features$chrom), start = cnloh_features$start,
                          end = cnloh_features$end, label = cnloh_features$label, colour = "grey40")
          labels = rbind(labels, fl[!fl$label %in% labels$label, , drop = FALSE])
        }
      }
      if(verbose) message(if(is.null(cnloh_features)) 0 else nrow(cnloh_features),
                          " recurrently deleted region(s) under copy-neutral LOH")
    }
  }
  # labels go outside the ideogram: the plot is set up without its default
  # ideogram, the labels drawn on the outermost ring, then the chromosome
  # names and the ideogram, then the data rings inside
  species = if(this_projection == "hg38") "hg38" else "hg19"
  circos.clear()
  if(!is.null(labels) && nrow(labels) > 0){
    circos.initializeWithIdeogram(chromosome.index = chrom_list, species = species, plotType = NULL)
    # ring heights relative to the circle (not mm, which on large devices is
    # smaller than circlize's cell padding)
    circos.genomicLabels(labels[, c("chr", "start", "end", "label")], labels.column = "label",
                         side = "outside", line_col = labels$colour, cex = 0.7,
                         connection_height = 0.04)
    circos.track(ylim = c(0, 1), track.height = 0.04, bg.border = NA, cell.padding = c(0, 0, 0, 0),
                 panel.fun = function(x, y){
                   circos.text(CELL_META$xcenter, 0.4, sub("^chr", "", CELL_META$sector.index),
                               cex = 0.7, facing = "inside", niceFacing = TRUE)
                 })
    circos.genomicIdeogram(species = species, track.height = 0.03)
  }else{
    circos.initializeWithIdeogram(chromosome.index = chrom_list, species = species)
  }
  if(include_cnv && nrow(cnv_df) > 0){
    bed = data.frame(chr = with_chr(cnv_df$chrom), start = cnv_df$start, end = cnv_df$end,
                     value1 = cnv_df$.value)
    circos.genomicTrackPlotRegion(bed, stack = TRUE, panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = cn_col_fun(value), border = NA, posTransform = NULL, ...)
    }, bg.border = NA)
  }
  # LOH: a thin ring just inside the CN ring, filled where LOH_flag is 1 with
  # the types in loh_types, coloured by the copy number change with it
  loh_segments = NULL
  if(show_loh && include_cnv){
    if(loh_assessed){
      loh = if(is.null(loh_typed)) NULL else loh_typed[loh_typed$loh_type %in% loh_shown, , drop = FALSE]
      if(!is.null(loh) && nrow(loh) > 0){
        loh_segments = loh
        bed = data.frame(chr = with_chr(loh$chrom), start = loh$start, end = loh$end,
                         colour = unname(loh_colours[loh$loh_type]), stringsAsFactors = FALSE)
        circos.genomicTrackPlotRegion(bed, ylim = c(0, 1), track.height = loh_track_height,
                                      panel.fun = function(region, value, ...){
                                        circos.genomicRect(region, value, ytop = 1, ybottom = 0,
                                                           col = value$colour, border = NA, ...)
                                      }, bg.border = "grey85", bg.lwd = 0.3, cell.padding = c(0, 0, 0, 0))
      }else{
        # assessed, none of these types: an empty ring
        circos.track(ylim = c(0, 1), track.height = loh_track_height, bg.border = "grey85", bg.lwd = 0.3,
                     cell.padding = c(0, 0, 0, 0))
      }
    }else{
      message("cnv_df has no LOH_flag values for this sample (e.g. capture data), so LOH isn't shown")
    }
  }
  if(include_ssm && nrow(ssm_df) > 0){
    bed = data.frame(chr = with_chr(ssm_df$Chromosome), start = ssm_df$Start_Position,
                     end = ssm_df$Start_Position, y = 0.5,
                     colour = unname(mutation_colours[as.character(ssm_df$Variant_Classification)]))
    bed$colour[is.na(bed$colour)] = "grey50"
    circos.genomicTrackPlotRegion(bed, ylim = c(0, 1), panel.fun = function(region, value, ...) {
      circos.genomicPoints(region, value, numeric.column = "y", col = value$colour, pch = 16, cex = 0.4, ...)
    }, bg.border = NA, track.height = 0.08)
  }

  if(include_sv && nrow(sv_df) > 0){
    bed1 = data.frame(chrom = with_chr(sv_df$CHROM_A), start = sv_df$START_A, end = sv_df$END_A)
    bed2 = data.frame(chrom = with_chr(sv_df$CHROM_B), start = sv_df$START_B, end = sv_df$END_B)
    if(is.null(link_colour_column)){
      # lighter when the oncogene SVs are drawn on top
      link_col = if(highlight_oncogene_sv && !is.null(oncogene_sv)) "#bdbdc180" else "#bdbdc1"
    }else if(link_colour_column %in% c("CHROM_A", "CHROM_B")){
      chrom_cols = GAMBLR.helpers::get_gambl_colours("chromosome")
      link_col = unname(chrom_cols[if(link_colour_column == "CHROM_A") bed1$chrom else bed2$chrom])
      link_col[is.na(link_col)] = "#bdbdc1"
    }else{
      if(!link_colour_column %in% colnames(sv_df)){
        stop(paste(link_colour_column, "is not an existing column in your sv data frame"))
      }
      link_col = sv_df[[link_colour_column]]
    }
    circos.genomicLink(bed1, bed2, col = link_col)
    if(highlight_oncogene_sv && !is.null(oncogene_sv)){
      circos.genomicLink(
        data.frame(chrom = with_chr(oncogene_sv$onco_chrom), start = oncogene_sv$onco_pos,
                   end = oncogene_sv$onco_pos + 1),
        data.frame(chrom = with_chr(oncogene_sv$partner_chrom), start = oncogene_sv$partner_pos,
                   end = oncogene_sv$partner_pos + 1),
        col = onco_colours[oncogene_sv$gene], lwd = 2)
    }
  }

  sample_id_position = match.arg(sample_id_position, c("topleft", "centre", "none"))
  if(sample_id_position == "centre") text(0, 0, this_sample_id, cex = 0.8)
  # the sample ID in the top left corner (the centre is busy with links),
  # above the metadata legend
  id_height = grid::unit(0, "mm")
  if(sample_id_position == "topleft"){
    # base graphics like the rest of the plot, at the plot area's top-left
    # corner without clipping (grid text there can be clipped on some
    # devices, e.g. RStudio's plot pane)
    usr = graphics::par("usr")
    graphics::text(usr[1], usr[4], this_sample_id, adj = c(0, 1), font = 2, cex = 0.9, xpd = NA)
    id_height = grid::unit(7, "mm")
  }

  out = if(return_data) list(cnv = if(include_cnv) cnv_df, sv = if(include_sv) sv_df,
                             ssm = if(include_ssm) ssm_df, oncogene_sv = oncogene_sv,
                             cnv_features = called_features, loh = loh_segments,
                             cnloh_features = cnloh_features)
  if(hide_legend){
    return(if(return_data) invisible(out) else invisible(NULL))
  }
  if(!missing(legend_metadata_columns)){
    samp_meta = as.data.frame(these_samples_metadata)
    samp_meta = samp_meta[samp_meta$sample_id == this_sample_id, , drop = FALSE]
    if(nrow(samp_meta) > 1 && "seq_type" %in% colnames(samp_meta) && this_seq_type %in% samp_meta$seq_type){
      samp_meta = samp_meta[samp_meta$seq_type == this_seq_type, , drop = FALSE]
    }
    these_meta = samp_meta[1, legend_metadata_columns, drop = FALSE]
    these_meta = mutate_if(these_meta, is.factor, as.character)
    vals = as.character(these_meta)
    # missing values (e.g. an NA factor level) as "NA", named by their column
    # unless names are given, since Legend() can't take NA labels
    missing_value = is.na(vals)
    vals[missing_value] = "NA"
    # each value's colour within its column's mapping (grey when there's none)
    all_cols = map_metadata_to_colours(legend_metadata_columns, these_meta)
    cols = vapply(seq_along(vals), function(i){
      if(missing_value[i]) return("grey80")   # (GAMBL's NA colour is white)
      v = all_cols[[legend_metadata_columns[i]]][vals[i]]
      if(length(v) == 0 || is.na(v)) "grey70" else unname(v)
    }, character(1))
    if(length(legend_metadata_names) == length(legend_metadata_columns)){
      for(i in seq_along(vals)){
        if(!legend_metadata_names[i] == ""){
          vals[i] = paste(legend_metadata_names[i], vals[i])
        }else if(missing_value[i]){
          vals[i] = paste(legend_metadata_columns[i], vals[i])
        }
      }
    }else{
      vals[missing_value] = paste(legend_metadata_columns[missing_value], "NA")
    }
    lgd_discrete = Legend(labels = vals, title_position = "topleft", legend_gp = gpar(fill = cols))
    draw(lgd_discrete, x = unit(3, "mm"), y = unit(1, "npc") - unit(5, "mm") - id_height,
         just = c("left", "top"))
  }
  if(!is.null(onco_colours)){
    lgd_onco = Legend(labels = names(onco_colours), title = "oncogene SV", type = "lines",
                      legend_gp = gpar(col = onco_colours, lwd = 2), title_position = "topleft")
    # right-aligned by its widest part: the title can be wider than the keys
    title_width = grid::convertWidth(grid::grobWidth(grid::textGrob(
      "oncogene SV", gp = grid::gpar(fontsize = 10, fontface = "bold"))), "mm", valueOnly = TRUE)
    legend_width = grid::convertWidth(grid::grobWidth(lgd_onco@grob), "mm", valueOnly = TRUE)
    draw(lgd_onco, x = unit(1, "npc") - unit(3 + max(0, title_width - legend_width), "mm"),
         y = unit(1, "npc") - unit(5, "mm"), just = c("right", "top"))
  }
  if(show_loh && include_cnv){
    lgd_loh = if(loh_assessed){
      shown = intersect(names(loh_colours), if(is.null(loh_segments)) character(0) else loh_segments$loh_type)
      if(length(shown) == 0) Legend(labels = paste(paste(loh_shown, collapse = ", "), "\nnot detected"),
                                    title = "LOH", legend_gp = gpar(fill = "white"))
      else Legend(labels = shown, title = "LOH (inner ring)", legend_gp = gpar(fill = loh_colours[shown]))
    }else{
      Legend(labels = "not assessed", title = "LOH", legend_gp = gpar(fill = "white"))
    }
    draw(lgd_loh, x = unit(3, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
  }
  if(include_cnv){
    legend_breaks = round(seq(cn_colour_range[1], cn_colour_range[3], length.out = 5), 2)
    lgd_cnv = Legend(at = legend_breaks,
                     col_fun = cn_col_fun,
                     title_position = "topleft",
                     title = if(cn_value == "CN") "CN" else "log\nratio")
    draw(lgd_cnv, x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))
  }
  if(return_data) invisible(out) else invisible(NULL)
}

# One colour per row of gene_bed, from the value (.value) of the cnv_df segment
# overlapping the gene's midpoint (first match); black without one
.circos_gene_cn_colours = function(gene_bed, cnv_df, col_fun, verbose = FALSE){
  gene_chr = gsub("^chr", "", as.character(gene_bed$chromosome))
  cnv_chr = gsub("^chr", "", as.character(cnv_df$chrom))
  cnv_start = as.numeric(cnv_df$start)
  cnv_end = as.numeric(cnv_df$end)
  gene_mid = (as.numeric(gene_bed$start) + as.numeric(gene_bed$end)) / 2
  vapply(seq_len(nrow(gene_bed)), function(i){
    hits = which(cnv_chr == gene_chr[i] & cnv_start <= gene_mid[i] & cnv_end >= gene_mid[i])
    if(length(hits) == 0){
      if(verbose){
        message("colour_labels_by_cn: no cnv_df segment overlaps ", gene_bed$gene_name[i],
                "; leaving its connector line black.")
      }
      return("black")
    }
    col_fun(cnv_df$.value[hits[1]])
  }, character(1))
}

# SVs involving oncogenes (GAMBLR.utils::annotate_sv), with the oncogene's
# end and the partner's end: annotate_sv doesn't keep the oncogene on a
# fixed end, so the end in (or nearest) the oncogene is found from the
# oncogene coordinates
.circos_oncogene_sv = function(sv_df, genome_build, oncogenes = NULL){
  sv = as.data.frame(sv_df)
  needed = c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "tumour_sample_id",
             "STRAND_A", "STRAND_B")
  if(!all(needed %in% colnames(sv))){
    warning("sv_df needs ", paste(setdiff(needed, colnames(sv)), collapse = ", "),
            " to find oncogene SVs; none are highlighted")
    return(NULL)
  }
  sv$CHROM_A = gsub("^chr", "", as.character(sv$CHROM_A))
  sv$CHROM_B = gsub("^chr", "", as.character(sv$CHROM_B))
  a = tryCatch(suppressMessages(GAMBLR.utils::annotate_sv(sv, genome_build = genome_build)),
               error = function(e){ warning("annotate_sv failed: ", conditionMessage(e)); NULL })
  if(is.null(a) || nrow(a) == 0) return(NULL)
  a = as.data.frame(a)
  a = a[!is.na(a$gene), , drop = FALSE]
  if(!is.null(oncogenes)) a = a[a$gene %in% oncogenes, , drop = FALSE]
  if(nrow(a) == 0) return(NULL)
  onco = as.data.frame(if(genome_build == "hg38") GAMBLR.data::hg38_oncogene else GAMBLR.data::grch37_oncogene)
  onco$chrom = gsub("^chr", "", as.character(onco$chrom))
  gap = function(chrom, pos, gene){
    vapply(seq_along(chrom), function(i){
      g = onco[onco$gene == gene[i] & onco$chrom == chrom[i], , drop = FALSE]
      if(nrow(g) == 0) return(Inf)
      min(pmax(0, g$start - pos[i], pos[i] - g$end))
    }, numeric(1))
  }
  c1 = gsub("^chr", "", as.character(a$chrom1)); c2 = gsub("^chr", "", as.character(a$chrom2))
  p1 = as.numeric(a$start1); p2 = as.numeric(a$start2)
  first = gap(c1, p1, a$gene) <= gap(c2, p2, a$gene)
  data.frame(tumour_sample_id = a$tumour_sample_id, gene = a$gene, partner = a$partner,
             onco_chrom = ifelse(first, c1, c2), onco_pos = ifelse(first, p1, p2),
             partner_chrom = ifelse(first, c2, c1), partner_pos = ifelse(first, p2, p1),
             fusion = a$fusion, stringsAsFactors = FALSE)
}
