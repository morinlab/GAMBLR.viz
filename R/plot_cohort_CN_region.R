#' @title Plot the overlap of copy number events in a region
#'
#' @description Show where copy number events in a set of samples overlap
#' within one genomic region, to find the minimal common region: a frequency
#' track (the fraction of samples with an event at each position, stacked by
#' event type) above a tiling plot (one row per sample with an event, its
#' events drawn as bars and sorted so that shared regions line up). Samples
#' can be split into groups (e.g. subtypes), each with its own frequency
#' track, tiling panel and minimal common region.
#'
#' @details \strong{Inputs.} `seg_data` can be generic segmented data (e.g.
#' with columns ID, chrom, start, end and CN (total
#' copy number) or log.ratio, and optionally LOH_flag (1 when only one
#' parental allele remains). It can also be the output of
#' [GAMBLR.utils::summarize_battenberg], whose allele-specific segment
#' properties are used: segments whose state is unresolved (e.g. disagreeing
#' subclonal states) are left out, as gaps.
#'
#' \strong{Filler segments.} Some segmented data contain placeholder rows
#' for uncalled regions (e.g. chromosome ends, acrocentric p arms), reported
#' at the diploid level; in a high-ploidy sample they look like deletions.
#' With `ignore_filler = TRUE` (default), diploid rows (CN 2, or log.ratio
#' 0) whose exact coordinates are shared by at least `filler_min_samples`
#' samples are treated as placeholders: they never count as events and are
#' left out of each sample's baseline, but the sample still counts in the
#' denominator (as with uncalled regions elsewhere in GAMBLR). Rows with
#' `dummy_segment == 1` (as flagged by `GAMBLR.results::get_cn_segments`)
#' are always treated this way, whatever `ignore_filler` is.
#'
#' \strong{Event types} (legend labels in brackets). With CN, each segment is compared with the sample's
#' baseline (`baseline` and `baseline_method`, or 2).
#' With LOH_flag too, losses and copy-neutral segments are split by LOH:
#' \itemize{
#'   \item `homdel` (HOMDEL): homozygous deletion (CN 0).
#'   \item `deletion_loh` (DEL-LOH): deletion with LOH (e.g. 1+0 in a diploid genome).
#'   \item `deletion` (DEL): deletion; with LOH_flag, one that keeps both alleles
#'     (e.g. 2+1 or 1+1 in a tetraploid genome).
#'   \item `deletion_loh_unknown` (DEL-unk): a deletion whose LOH_flag is NA (e.g.
#'     samples from a caller that doesn't report LOH, mixed with ones that
#'     do).
#'   \item `cnloh` (CN-LOH): copy-neutral LOH.
#'   \item `amp` (AMP): amplification, CN at least `amp_threshold` copies
#'     above the baseline (whatever the LOH status).
#'   \item `gain_loh` (GAIN-LOH), `gain` (GAIN), `gain_loh_unknown`
#'     (GAIN-unk): other gains, with, without or with unknown LOH.
#' }
#' Without CN, log.ratio is used: below `-log_ratio_threshold` is a deletion
#' (at or below `homdel_log_ratio` a homozygous deletion), above
#' `log_ratio_threshold` a gain and at or above `amp_log_ratio` an
#' amplification; LOH can't be assessed this way. `direction`
#' chooses the losses (homdel, deletion_loh, deletion, cnloh) or the gains
#' (amp, gain_loh, gain);
#' `event_types` can restrict them further.
#'
#' \strong{Frequencies and the minimal common region.} The frequency at a
#' position is the fraction of samples with an event of a shown type there;
#' each sample counts once. With `missing_data = "unknown"` (default), the
#' denominator at each position is the group's samples with data there (a
#' non-placeholder segment), so a stretch a sample has no call for (e.g. a
#' gap in capture data) is unknown for it rather than counted as no event,
#' and doesn't create a dip. With "no_event", every sample with data
#' anywhere in the window is in the denominator everywhere. The bar chart's
#' denominator follows the same rule (samples with data in the minimal
#' common region). The minimal common region is the stretch
#' (or stretches) of at least `min_mcr_size` over which the frequency stays
#' highest: the highest frequency held across `min_mcr_size` or more, so a
#' sliver created by a few samples' breakpoints can't define it. By default it comes from all the samples given, pooled
#' across groups (`mcr_by = "all"`), and the same region is outlined for
#' every group; `mcr_by = "group"` finds each group's own. It is calculated
#' from the data, not taken from a curated feature table, and depends on the
#' event types shown (e.g. long copy-neutral LOH raises the frequency across
#' a whole arm). It is outlined on top of both panels, at least 0.5\% of the
#' region wide so that it stays visible. With two groups and `mirror = TRUE`, the second group's
#' frequencies are drawn downward.
#'
#' \strong{Tiling order.} Within each group, samples whose events cover the
#' minimal common region (or, with `tiling_focus`, the genes or a given
#' region) come first, followed by samples with
#' events elsewhere in the region. With `tiling_order = "type"` (default),
#' the covering samples are grouped by their most severe event type in the
#' minimal common region (e.g. AMP, then GAIN), and within each type ordered
#' by where the event starts and then by its end, so each type forms its own
#' staircase around the shared region. With "position", they form a single
#' staircase.
#' Starts and ends are compared at `tiling_resolution` (2\% of the window
#' by default), and ties (e.g. events spanning the whole window) are broken
#' by the most severe event type in the minimal common region (or in the
#' event), so samples with similar extents are grouped by type, and then by
#' how much of the event is of that type (e.g. amplified across the whole
#' window before amplified only around the minimal common region).
#' Event segments of a sample less than `min_mcr_size` apart count as one
#' event for this ordering (each segment is still coloured by its type).
#'
#' @param seg_data Segmented data or [GAMBLR.utils::summarize_battenberg]
#' output (see Details). Use whole-genome data so each sample's baseline can
#' be estimated.
#' @param region Region as "chrom:start-end" (e.g. "9:20000000-23000000" or
#' "chr9:20,000,000-23,000,000"). Alternatively give `chrom`, `start` and
#' `end`, or only `genes` (the region then spans them plus `flank`).
#' @param chrom,start,end Region coordinates, instead of `region`.
#' @param these_samples_metadata Optional metadata with sample_id (and the
#' `split_by` column). When given, only its samples are shown. When it has
#' a `cn_available` column (as from `GAMBLR.results::get_data_by_region`),
#' only samples with TRUE are plotted, with a message.
#' @param split_by Optional metadata column defining groups.
#' @param split_order Optional order of the groups (top to bottom).
#' @param direction "loss" (default) or "gain".
#' @param event_types Optional subset of event types to show (see Details).
#' @param baseline "sample" (default; each sample's own baseline copy
#' number, see `baseline_method`) or "diploid" (CN 2; use this for
#' segments that are already centred, which is chosen automatically for
#' output of `GAMBLR.utils::center_segments`). Ignored for
#' log.ratio input. With "sample", a `baseline_cn` column in `seg_data`
#' (e.g. from `GAMBLR.results::get_cn_segments(region = ...)`, which
#' computes it genome-wide) is used; otherwise the baseline is estimated
#' from the segments, which must then cover each sample's genome (see
#' `GAMBLR.utils::check_segment_coverage`): segments filtered to a region,
#' or very sparse samples, stop with an error naming them.
#' @param baseline_method How each sample's baseline is defined (see
#' `GAMBLR.utils::add_segment_baseline`): "mode" (default; the mean copy
#' number within the sample's dominant state, so fractional offsets such as
#' a neutral level at 1.9 are removed), "mode_rounded" (that state as a
#' whole number) or "mean". Events are called on the copy number shifted by
#' the baseline (`round(CN - baseline + 2)`); with "mode_rounded" this is
#' the same as comparing the rounded copy number with the baseline. A
#' `baseline_cn` column of another method stops with an error, except that
#' "mode" baselines also serve "mode_rounded".
#' @param show_loh Set to FALSE to ignore LOH status: deletions (DEL-LOH,
#' DEL, DEL-unk) are shown as one type (DEL), as are gains (GAIN-LOH, GAIN,
#' GAIN-unk as GAIN), and copy-neutral LOH, which changes no copy number, is
#' not shown. Homozygous deletions and amplifications are kept. Default TRUE.
#' @param amp_threshold Copies above the baseline for an amplification.
#' Default 2 (e.g. CN 4 or more in a diploid genome, 6 or more in a
#' tetraploid one).
#' @param log_ratio_threshold,homdel_log_ratio,amp_log_ratio Thresholds for
#' log.ratio input. Defaults 0.2, -2 and 1.
#' @param genes Optional genes to mark: gene symbols (located with
#' [GAMBLR.utils::gene_to_region]) or a data frame with gene, chrom, start
#' and end.
#' @param label_genes Optional other genes to mark and label when they fall
#' in the window (gene symbols, or a data frame like `genes`), e.g. a list of
#' lymphoma genes. They are drawn a shade lighter than `genes` and don't
#' affect the region or `tiling_focus`. Genes not found or outside the
#' window are skipped.
#' @param gene_track Draw the genes (`genes`, `label_genes`) as a track at
#' the top: each from its start to its end, with strand when given as a
#' data frame with a strand column ("+"/"-"), and its name. Default TRUE;
#' with FALSE, genes are labelled in the frequency track instead. Either way
#' they are shaded in the data panels.
#' @param gene_models Gene models for the gene track (exons, with UTRs
#' thinner than coding parts, introns and strand). By default (NULL) they are retrieved for the window from
#' mygene.info with [GAMBLR.utils::get_gene_models] (needs internet; without
#' it genes are drawn as plain bars, with a message). Give a data frame in
#' that format to use your own (e.g. edited or cached) models, or FALSE to
#' skip them. Genes in the window that aren't in `genes` or `label_genes`
#' are added to the track (lighter; when retrieved automatically, only if
#' there are 30 or fewer) but not shaded in the plot.
#' @param genes_show_isoforms Optional gene symbols whose transcripts are all drawn
#' in the gene track, one row each (labelled with the transcript), instead
#' of the gene's merged model. With automatic gene models they are retrieved
#' too; with your own `gene_models`, give per-transcript rows for them (e.g.
#' from `GAMBLR.utils::get_gene_models(..., transcripts = "all")`). Can make
#' the track tall for genes with many isoforms.
#' @param isoforms Which transcripts of the `genes_show_isoforms` genes are
#' drawn (it has no effect without `genes_show_isoforms`):
#' "distinct_ends" (default): one per distinct transcription start site or
#' 3' end (strand-aware, ends within 100 bp treated as the same), i.e.
#' versions differing only internally (e.g. skipped exons) are combined,
#' drawn as the one with the most exons and labelled with how many it stands
#' for, e.g. "FOXP1 NM_032682 (+4)"; "all": every transcript.
#' @param other_genes Genes in the window besides `genes` and `label_genes`
#' (from the gene models) in the gene track: "auto" (default; added when
#' there are 30 or fewer, or always with your own `gene_models`), "none"
#' (only the genes you named) or "all".
#' @param shade_genes Shade the genes you named (`genes`, `label_genes`) in
#' the data panels, from start to end (narrow ones as a dashed line).
#' Default TRUE; FALSE leaves them only in the gene track.
#' @param flank Distance added on each side of `genes` when they define the
#' region. Default 5e6.
#' @param genome_build Genome build, for locating `genes`. Default "grch37".
#' @param mirror With two groups, draw the second group's frequencies
#' downward on one track (default TRUE); otherwise one track per group.
#' @param min_mcr_size Smallest stretch (bp) that can be a minimal common
#' region: it is where the highest frequency is held over at least this
#' length. Default 1e4.
#' @param missing_data How stretches a sample has no data for count:
#' "unknown" (default; left out of the denominator there) or "no_event".
#' See Details.
#' @param ignore_filler,filler_min_samples Treat placeholder rows as no
#' event (see Details). Defaults TRUE and 3.
#' @param show_mcr Outline the minimal common region. Default TRUE.
#' @param mcr_by Find the minimal common region from "all" samples (default)
#' or within each "group". See Details.
#' @param mcr_colour Colour of the minimal common region outline. Default
#' "grey10".
#' @param colours Optional named vector of colours for the event types,
#' replacing the defaults (GAMBLR's copy number colours: homdel "#4393C3",
#' deletion_loh "#92C5DE", deletion "#D1E5F0", deletion_loh_unknown
#' "#BFC8CF", cnloh "#AFBBF2", amp "#67001F", gain_loh "#B2182B", gain
#' "#D6604D", gain_loh_unknown "#E3C4BC").
#' @param base_size Base font size. Default 9.
#' @param rel_heights Relative heights of the frequency track and the tiling
#' plot. Default c(1, 2).
#' @param show_tiling Show the lower section, one row per sample's events.
#' Default TRUE; FALSE leaves only the upper section (frequency track, MCR
#' summary and any SV/SSM density tracks), e.g. when a region has too many
#' events for the rows to be readable. The tiling order and plot are still
#' computed and returned with `return_data = TRUE`.
#' @param sv_data,ssm_data Optional SV calls and mutations (MAF) to show as
#' density tracks between the frequency track and the tiling plot, built as
#' by [plot_sv_density] and [plot_ssm_density] (see there for the columns
#' needed). Each bin shows the fraction of each group's samples with data
#' that have a breakpoint or mutation there.
#' @param sv_samples,ssm_samples The samples with SV or mutation data (each
#' track's denominator). Default: the plotted samples with `sv_available` /
#' `ssm_available` TRUE when `these_samples_metadata` has those columns (as
#' from `GAMBLR.results::get_data_by_region`), otherwise all the samples
#' plotted, i.e. the data are assumed to cover every sample. Give e.g. only
#' the genomes when that isn't so: capture data have no SV calls and cover
#' only their targets.
#' @param sv_height,ssm_height Height of each track relative to the
#' frequency track. Default 0.35.
#' @param sv_args,ssm_args Optional lists of other arguments for the tracks,
#' as in [plot_sv_density] and [plot_ssm_density]: `measure` (default
#' "samples" here, matching the copy number tracks; "patients" counts each
#' patient once per group; for SSMs also "mutations"), `bin_size`,
#' `bin_fraction`, `max_bins`, `colour` (default: each group in its GAMBLR
#' palette colour), `group_colours`, `class_colours`, `y_scale`, the outlier
#' settings (`outliers`, `outlier_ratio`, `outlier_k`, `outlier_min_events`,
#' `outlier_scope`; isolated outliers are excluded by default), and
#' `colour_by`
#' ("type" or "partner" for SVs; see [plot_ssm_density] for SSMs) and
#' `max_classes` (SSM).
#' @param tiling_label Y-axis label of the tiling plot. Default NULL: chosen
#' from the event types present, e.g. "deletion and LOH events" (losses with
#' LOH information), "deletion events" (losses without it, or with only
#' deletion types shown) or "gain events".
#' @param tiling_min_rows Minimum height of each group's tiling panel, in
#' sample rows, so that small groups keep room for their label. Default
#' NULL: 15\% of the largest group's rows, and at least 12.
#' @param tiling_order How samples whose events cover the minimal common
#' region are ordered in the tiling plot: "type" (default) groups them by
#' their most severe event type there (e.g. all AMP, then all GAIN), each
#' type as its own staircase; "position" makes one staircase, using the type
#' only to break ties. See Details.
#' @param tiling_focus What the tiling plot is ordered around: "mcr"
#' (default; the minimal common region), "gene" (the span of `genes`) or
#' `c(start, end)` coordinates. Samples whose events cover it come first,
#' grouped and staircased by their event there (see `tiling_order`). The
#' minimal common region is still outlined and summarised.
#' @param tiling_resolution Event starts and ends within this fraction of
#' the window are treated as the same when ordering the tiling plot, so
#' that similar events are grouped by type. Default 0.02.
#' @param frequency_label Y-axis label of the frequency track. Default
#' "samples".
#' @param show_mcr_summary Add a stacked bar per group (on one 0-100\%
#' scale) of the share of its samples with each event type in the minimal
#' common region, beside the frequency track. Default TRUE.
#' @param mcr_summary_include_none With TRUE (default), each bar is all of
#' the group's samples, including those with no event in the minimal common
#' region (white), so bar heights compare how often each group is affected.
#' With FALSE, each bar is only the samples with an event there (their
#' number in the bar's label), to compare the relative contribution of the event
#' types between groups.
#' @param mcr_summary_width Width of that bar chart relative to the main
#' plot. Default 0.08. The bars have no axis: each is labelled inside with
#' its group and number of samples.
#' @param legend_position "right" (default): the legend is stacked
#' vertically in the right margin, below the bar chart. "top": a horizontal
#' legend above the frequency track. Without the bar chart
#' (`show_mcr_summary = FALSE`) the legend is always on top.
#' @param legend_width Width of the right margin (bar chart and legend)
#' relative to the main plot when `legend_position = "right"`. Default 0.12.
#' @param return_data Set to TRUE to also return the frequencies, the
#' minimal common regions and the plotted events.
#'
#' @return The combined plot (cowplot). With `return_data = TRUE`, a list:
#' \itemize{
#'   \item `plot`: the combined plot.
#'   \item `region`: the window shown, as "chrN:start-end" (see
#'     [region_window]).
#'   \item `events`: the rows of `seg_data` (or of its `segments` table)
#'     behind the plot: segments in the window with an event of a shown
#'     type, as given (not clipped), with their `group`, `baseline_cn`,
#'     `type` and `event` (legend label).
#'   \item `sv_events`, `ssm_events`: with `sv_data` or `ssm_data`, their rows
#'     behind the tracks (as from [plot_sv_density] and [plot_ssm_density]).
#'   \item `data`: every segment in the region (clipped to it), one row per
#'     segment: ID, group, chrom, start, end, size, CN and/or log.ratio,
#'     `baseline_cn` (the sample's baseline), LOH_flag, `type` (event type,
#'     NA if none), `event` (its legend label, for shown types), `shown`,
#'     `in_mcr` (overlaps the minimal common region used for its group) and
#'     `row` (its row in the tiling plot).
#'   \item `mcr`: the minimal common region(s), with the frequency held
#'     across it (`fraction`), the highest frequency within it
#'     (`peak_fraction`) and the group ("all" when pooled).
#'   \item `mcr_summary`: per group, the number and fraction of samples by
#'     event type in the minimal common region (each sample counted once, by
#'     its most severe type: homozygous deletion, then deletion with LOH,
#'     copy-neutral LOH, other deletions; gains with LOH before other gains),
#'     and "none"; `fraction` is of all the group's samples and
#'     `fraction_of_events` of those with an event there
#'     (`n_with_event`).
#'   \item `tiling_order`: the tiling plot's samples, top to bottom per
#'     group, with what they are sorted by: whether their event covers the
#'     focus (`covers_focus`: the minimal common region, or as set by
#'     `tiling_focus`), that event's start and end (`block_start`,
#'     `block_end`), its most severe type there (`block_type`) and the share
#'     of the event at that type (`block_type_share`).
#'   \item `sv_density`, `ssm_density`: with `sv_data` or `ssm_data`, per bin,
#'     group and class, the number (`n`) and share (`value`) of samples with
#'     data (`n_units`) that have a breakpoint or mutation there.
#'   \item `frequencies`: per interval, group and event type, the number of
#'     samples with the event (`n`), the number with data there
#'     (`n_with_data`, the denominator with `missing_data = "unknown"`) and
#'     the fraction.
#'   \item `frequency_plot`, `tiling_plot`, `mcr_summary_plot`: the panels
#'     as ggplot objects.
#' }
#'
#' @import ggplot2 cowplot dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' seg = get_cn_segments(these = meta, projection = "grch37")
#' # 9p21 losses by mechanism, by COO
#' plot_cohort_CN_region(seg, genes = "CDKN2A", flank = 4e6,
#'                    these_samples_metadata = meta, split_by = "COO_consensus")
#'
#' # from Battenberg results (unresolved subclonal segments left out)
#' bb = GAMBLR.utils::summarize_battenberg(battenberg_dir = "battenberg_results")
#' plot_cohort_CN_region(bb, region = "1:0-15000000", genes = "TNFRSF14",
#'                    these_samples_metadata = meta, split_by = "group")
#' }
plot_cohort_CN_region = function(seg_data,
                              region = NULL,
                              chrom = NULL,
                              start = NULL,
                              end = NULL,
                              these_samples_metadata = NULL,
                              split_by = NULL,
                              split_order = NULL,
                              direction = "loss",
                              event_types = NULL,
                              baseline = "sample",
                              baseline_method = "mode",
                              log_ratio_threshold = 0.2,
                              homdel_log_ratio = -2,
                              amp_threshold = 2,
                              amp_log_ratio = 1,
                              show_loh = TRUE,
                              genes = NULL,
                              flank = 5e6,
                              label_genes = NULL,
                              gene_track = TRUE,
                              gene_models = NULL,
                              genes_show_isoforms = NULL,
                              isoforms = "distinct_ends",
                              shade_genes = TRUE,
                              other_genes = "auto",
                              genome_build = "grch37",
                              mirror = TRUE,
                              show_mcr = TRUE,
                              mcr_by = "all",
                              min_mcr_size = 1e4,
                              mcr_colour = "grey10",
                              ignore_filler = TRUE,
                              missing_data = "unknown",
                              filler_min_samples = 3,
                              colours = NULL,
                              base_size = 9,
                              rel_heights = c(1, 2),
                              show_tiling = TRUE,
                              sv_data = NULL,
                              sv_samples = NULL,
                              sv_height = 0.35,
                              sv_args = list(),
                              ssm_data = NULL,
                              ssm_samples = NULL,
                              ssm_height = 0.35,
                              ssm_args = list(),
                              tiling_label = NULL,
                              tiling_min_rows = NULL,
                              tiling_resolution = 0.02,
                              tiling_order = "type",
                              tiling_focus = "mcr",
                              show_mcr_summary = TRUE,
                              mcr_summary_include_none = TRUE,
                              mcr_summary_width = 0.08,
                              legend_position = "right",
                              legend_width = 0.12,
                              frequency_label = "samples",
                              return_data = FALSE){
  direction = match.arg(direction, c("loss", "gain"))
  baseline = match.arg(baseline, c("sample", "diploid"))
  # segments already centred (by GAMBLR.utils::center_segments) are on a
  # diploid scale; their baseline_cn is the uncentred one
  if(baseline == "sample" && "CN_uncentered" %in% colnames(seg_data)){
    message("seg_data is already centred (it has a CN_uncentered column), so events are ",
            "called relative to CN 2 (baseline = \"diploid\")")
    baseline = "diploid"
  }
  baseline_method = match.arg(baseline_method, GAMBLR.utils::baseline_methods())
  mcr_by = match.arg(mcr_by, c("all", "group"))
  missing_data = match.arg(missing_data, c("unknown", "no_event"))
  legend_position = match.arg(legend_position, c("right", "top"))
  tiling_order = match.arg(tiling_order, c("type", "position"))
  isoforms_given = !missing(isoforms)
  isoforms = match.arg(isoforms, c("distinct_ends", "all"))
  if(isoforms_given && length(genes_show_isoforms) == 0){
    warning("isoforms only applies to genes_show_isoforms, which is not set")
  }
  if(is.character(tiling_focus)){
    tiling_focus = match.arg(tiling_focus, c("mcr", "gene"))
  }else if(!is.numeric(tiling_focus) || length(tiling_focus) != 2 ||
           tiling_focus[1] >= tiling_focus[2]){
    stop("tiling_focus must be \"mcr\", \"gene\" or c(start, end)")
  }
  seg = .common_region_segments(seg_data)
  seg$filler = if(ignore_filler) .filler_segments(seg, filler_min_samples)
               else rep(FALSE, nrow(seg))
  # rows flagged as dummy (placeholder) segments are always treated this way
  if("dummy_segment" %in% colnames(seg)){
    seg$filler = seg$filler | seg$dummy_segment %in% 1
  }

  # the window and the genes to mark (index genes and label_genes)
  w = .resolve_window(region, chrom, start, end, genes, flank, label_genes, genome_build)
  w = .add_gene_models(w, .resolve_gene_models(w, gene_models, genome_build, gene_track,
                                               genes_show_isoforms),
                       max_extra_genes = .other_gene_cap(other_genes, gene_models),
                       genes_show_isoforms = genes_show_isoforms, isoforms = isoforms)
  if(!shade_genes && !is.null(w$marks)) w$marks$shade = FALSE
  chrom = w$chrom; start = w$start; end = w$end
  gene_df = w$genes
  marks = w$marks
  if(identical(tiling_focus, "gene") && (is.null(gene_df) || nrow(gene_df) == 0)){
    stop("tiling_focus = \"gene\" needs genes in the region")
  }

  # samples and groups
  if(!is.null(these_samples_metadata)){
    meta_df = as.data.frame(these_samples_metadata)
    # copy number availability annotated on the metadata (e.g. by
    # GAMBLR.results::get_data_by_region): only those samples are plotted
    if("cn_available" %in% colnames(meta_df)){
      unusable = meta_df$sample_id[!meta_df$cn_available %in% TRUE]
      if(length(unusable) > 0){
        message(sum(meta_df$cn_available %in% TRUE), " of ", nrow(meta_df),
                " samples have usable copy number data (cn_available) and are plotted")
      }
      meta_df = meta_df[meta_df$cn_available %in% TRUE, , drop = FALSE]
      these_samples_metadata = meta_df
    }
    seg = seg[seg$ID %in% meta_df$sample_id, , drop = FALSE]
  }
  if(!is.null(split_by)){
    if(is.null(these_samples_metadata) || !split_by %in% colnames(these_samples_metadata)){
      stop("split_by must name a column of these_samples_metadata")
    }
    group_of = stats::setNames(as.character(these_samples_metadata[[split_by]]),
                               these_samples_metadata$sample_id)
    seg$group = unname(group_of[seg$ID])
    seg = seg[!is.na(seg$group), , drop = FALSE]
  }else{
    seg$group = "all samples"
  }
  if(nrow(seg) == 0){
    stop("no segments for the samples given")
  }

  # classify every segment (baseline from the whole genome)
  typed = .common_region_types(seg, baseline, log_ratio_threshold,
                               homdel_log_ratio, amp_threshold, amp_log_ratio,
                               genome_build = genome_build,
                               baseline_method = baseline_method)
  seg$type = typed$type
  seg$baseline_cn = typed$baseline
  has_loh = "LOH_flag" %in% colnames(seg) && any(!is.na(seg$LOH_flag))
  if(!show_loh){
    # LOH status ignored: deletions and gains are not split by it, and
    # copy-neutral LOH (no copy number change) is not an event
    seg$type[seg$type %in% c("deletion_loh", "deletion_loh_unknown")] = "deletion"
    seg$type[seg$type %in% c("gain_loh", "gain_loh_unknown")] = "gain"
    seg$type[seg$type %in% "cnloh"] = NA_character_
    has_loh = FALSE
  }
  # most to least severe, for ordering and summarising
  severity = c("homdel", "deletion_loh", "cnloh", "deletion",
               "deletion_loh_unknown", "amp", "gain_loh", "gain",
               "gain_loh_unknown")
  shown = if(direction == "loss") c("homdel", "deletion_loh", "deletion",
                                    "deletion_loh_unknown", "cnloh")
          else c("amp", "gain_loh", "gain", "gain_loh_unknown")
  # unknown-LOH types only when some segments have them, and none split by
  # LOH when it is ignored
  shown = shown[!grepl("_unknown$", shown) | shown %in% seg$type]
  if(!show_loh){
    shown = setdiff(shown, c("deletion_loh", "cnloh", "gain_loh"))
  }
  if(!is.null(event_types)){
    bad = setdiff(event_types, shown)
    if(length(bad) > 0){
      stop(paste0("event_types for direction \"", direction, "\" are: ",
                  paste(shown, collapse = ", ")))
    }
    shown = intersect(shown, event_types)
  }
  # legend labels: DEL and GAIN keep both alleles when LOH is known
  labels = c(homdel = "HOMDEL", deletion_loh = "DEL-LOH", deletion = "DEL",
             deletion_loh_unknown = "DEL-unk", cnloh = "CN-LOH",
             amp = "AMP", gain_loh = "GAIN-LOH", gain = "GAIN",
             gain_loh_unknown = "GAIN-unk")
  # GAMBLR copy number colours (CN 0, CN 1, CN 3/4), with a lighter
  # step for deletions that keep both alleles
  default_colours = c(homdel = "#4393C3", deletion_loh = "#92C5DE",
                      deletion = "#D1E5F0", deletion_loh_unknown = "#BFC8CF",
                      cnloh = "#AFBBF2", amp = "#67001F", gain_loh = "#B2182B",
                      gain = "#D6604D", gain_loh_unknown = "#E3C4BC")
  if(!is.null(colours)){
    default_colours[names(colours)] = colours
  }

  # clip to the region; samples with data there form the denominators
  seg = seg[seg$chrom == chrom & seg$end > start & seg$start < end, , drop = FALSE]
  if(nrow(seg) == 0){
    stop("no segments in the region")
  }
  seg$start = pmax(seg$start, start); seg$end = pmin(seg$end, end)
  groups = unique(seg$group)
  if(!is.null(split_order)){
    groups = c(intersect(split_order, groups), setdiff(groups, split_order))
  }else if(!is.null(split_by) && is.factor(these_samples_metadata[[split_by]])){
    groups = intersect(levels(these_samples_metadata[[split_by]]), groups)
  }else{
    groups = sort(groups)
  }
  n_samples = vapply(groups, function(g) length(unique(seg$ID[seg$group == g])), numeric(1))
  events = seg[!is.na(seg$type) & seg$type %in% shown, , drop = FALSE]
  events$type = factor(events$type, levels = shown)

  # frequency per interval between breakpoints, per group and type
  bp = sort(unique(c(start, end, seg$start, seg$end)))
  n_int = length(bp) - 1
  # samples with data (a non-placeholder segment) at each interval
  with_data = function(g){
    d = seg[seg$group == g & !seg$filler, , drop = FALSE]
    n = numeric(n_int)
    for(id in unique(d$ID)){
      x = d[d$ID == id, , drop = FALSE]
      hit = rep(FALSE, n_int)
      for(i in seq_len(nrow(x))){
        hit[match(x$start[i], bp):(match(x$end[i], bp) - 1)] = TRUE
      }
      n = n + hit
    }
    n
  }
  frequencies = do.call(rbind, lapply(groups, function(g){
    counts = matrix(0, n_int, length(shown), dimnames = list(NULL, shown))
    e = events[events$group == g, , drop = FALSE]
    for(i in seq_len(nrow(e))){
      ix = match(e$start[i], bp):(match(e$end[i], bp) - 1)
      counts[ix, as.character(e$type[i])] = counts[ix, as.character(e$type[i])] + 1
    }
    # denominator: all of the group's samples, or only those with data there
    denom = if(missing_data == "unknown") with_data(g) else rep(n_samples[[g]], n_int)
    data.frame(group = g, start = rep(bp[-length(bp)], length(shown)),
               end = rep(bp[-1], length(shown)),
               type = rep(shown, each = n_int), n = as.vector(counts),
               n_with_data = rep(denom, length(shown)),
               fraction = ifelse(rep(denom, length(shown)) > 0,
                                 as.vector(counts) / rep(denom, length(shown)), 0),
               stringsAsFactors = FALSE)
  }))

  # minimal common region(s): where the total frequency peaks, over all
  # samples (mcr_by = "all") or within each group
  peak = function(total, label, n){
    size = bp[-1] - bp[-length(bp)]
    if(max(total) == 0) return(NULL)
    min_size = min(min_mcr_size, sum(size))
    # contiguous stretches of intervals at or above a frequency level
    runs_at = function(level){
      above = total >= level
      run = cumsum(c(1, diff(above) != 0))
      list(run = run, lens = tapply(size[above], run[above], sum))
    }
    # the highest level held over a stretch of at least min_size (lowering
    # the level only lengthens stretches, so binary search applies)
    levels = sort(unique(total[total > 0]), decreasing = TRUE)
    lo = 1; hi = length(levels)
    while(lo < hi){
      mid = (lo + hi) %/% 2
      if(any(runs_at(levels[mid])$lens >= min_size)) hi = mid else lo = mid + 1
    }
    level = levels[lo]
    r = runs_at(level)
    keep = as.numeric(names(r$lens)[r$lens >= min_size])
    do.call(rbind, lapply(keep, function(k){
      ix = which(r$run == k)
      data.frame(group = label, chrom = chrom, start = bp[min(ix)],
                 end = bp[max(ix) + 1], fraction = level,
                 peak_fraction = max(total[ix]), n_samples = n,
                 stringsAsFactors = FALSE)
    }))
  }
  # frequency of any shown event per interval: events over samples with data
  interval_fraction = function(f){
    key = factor(f$start, levels = bp[-length(bp)])
    n = as.numeric(tapply(f$n, key, sum))
    first_type = f$type == f$type[1]
    denom = as.numeric(tapply(f$n_with_data[first_type], key[first_type], sum))
    ifelse(denom > 0, n / denom, 0)
  }
  if(mcr_by == "all"){
    mcr = peak(interval_fraction(frequencies), "all", sum(n_samples))
  }else{
    mcr = do.call(rbind, lapply(groups, function(g){
      peak(interval_fraction(frequencies[frequencies$group == g, ]), g, n_samples[[g]])
    }))
  }
  # the region(s) used for each group's sorting and plotting
  mcr_of = function(g){
    if(is.null(mcr)) return(NULL)
    if(mcr_by == "all") mcr else mcr[mcr$group == g, , drop = FALSE]
  }

  # the region the tiling plot is ordered around: the minimal common
  # region, the genes, or given coordinates
  focus_of = function(g){
    if(identical(tiling_focus, "mcr")) return(mcr_of(g))
    if(identical(tiling_focus, "gene")){
      return(data.frame(start = min(gene_df$start), end = max(gene_df$end)))
    }
    data.frame(start = tiling_focus[1], end = tiling_focus[2])
  }
  # minimum panel height (rows): given, or scaled with the largest group
  rows_per_group = if(nrow(events) > 0) tapply(events$ID, events$group,
                                               function(v) length(unique(v))) else 0
  min_rows = if(is.null(tiling_min_rows)) max(12, ceiling(0.15 * max(rows_per_group)))
             else tiling_min_rows
  # tiling order: samples covering the group's minimal common region first
  events = events[order(events$ID, events$start), , drop = FALSE]
  sample_order = do.call(rbind, lapply(groups, function(g){
    e = events[events$group == g, , drop = FALSE]
    if(nrow(e) == 0) return(NULL)
    m = focus_of(g)
    keys = do.call(rbind, lapply(split(e, e$ID), function(s){
      # contiguous blocks of event segments: segments separated by less
      # than min_mcr_size (e.g. the 1 bp between consecutive 1-based
      # segments) belong to one block
      block = cumsum(c(1, s$start[-1] - s$end[-nrow(s)] > min_mcr_size))
      b_start = tapply(s$start, block, min); b_end = tapply(s$end, block, max)
      hit = if(!is.null(m) && nrow(m) > 0)
        which(b_start < m$end[1] & b_end > m$start[1]) else integer(0)
      k = if(length(hit) > 0) hit[1] else 1
      # tie-breaker: the most severe type in the region (or in the block)
      in_block = block == as.numeric(names(b_start)[k])
      if(length(hit) > 0){
        in_block = in_block & s$start < m$end[1] & s$end > m$start[1]
      }
      key_type = min(match(as.character(s$type[in_block]), severity))
      # share of the whole event at that type (solid rows before mixed ones)
      whole = block == as.numeric(names(b_start)[k])
      len = s$end - s$start
      data.frame(ID = s$ID[1], covers = length(hit) > 0,
                 key_start = b_start[k], key_end = b_end[k], key_type = key_type,
                 key_share = sum(len[whole & match(as.character(s$type), severity) == key_type]) /
                   sum(len[whole]))
    }))
    # positions compared at the plot's resolution (tiling_resolution of the
    # window), so events that look the same are grouped by type rather than
    # ordered by small breakpoint differences
    bin = max(1, tiling_resolution * (end - start))
    start_bin = floor((keys$key_start - start) / bin)
    end_bin = -floor((keys$key_end - start) / bin)
    keys = if(tiling_order == "type"){
      # grouped by type at the minimal common region, a staircase within each
      keys[order(!keys$covers, keys$key_type, start_bin, end_bin,
                 -round(keys$key_share, 2), keys$key_start, -keys$key_end), ]
    }else{
      # one staircase; type only breaks ties
      keys[order(!keys$covers, start_bin, end_bin, keys$key_type,
                 -round(keys$key_share, 2), keys$key_start, -keys$key_end), ]
    }
    # first sample at the top; small groups padded below to min_rows
    data.frame(ID = keys$ID, group = g,
               row = rev(seq_len(nrow(keys))) + max(0, min_rows - nrow(keys)),
               covers_focus = keys$covers, block_start = keys$key_start,
               block_end = keys$key_end,
               block_type = severity[keys$key_type],
               block_type_share = keys$key_share)
  }))
  events$row = rep(NA_real_, nrow(events))
  if(!is.null(sample_order)){
    events$row = sample_order$row[match(paste(events$group, events$ID),
                                        paste(sample_order$group, sample_order$ID))]
  }

  # plots
  if(is.null(tiling_label)){
    present = unique(as.character(events$type))
    loh_shown = any(c("cnloh", "deletion_loh", "gain_loh") %in% present)
    tiling_label = if(direction == "loss"){
      if(loh_shown) "deletion and LOH events" else "deletion events"
    }else{
      if(loh_shown) "gain and LOH events" else "gain events"
    }
  }
  group_label = stats::setNames(paste0(groups, "\n(n=", n_samples, ")"), groups)
  x_lim = c(start, end) / 1e6
  x_scale = scale_x_continuous(limits = x_lim, expand = c(0, 0))
  fill_scale = scale_fill_manual(values = default_colours[shown],
                                 labels = labels[shown], drop = FALSE,
                                 name = NULL)
  # genes in the data panels (shaded spans); labelled there only without a
  # gene track
  gene_layers = function(with_labels, label_facet = NULL){
    .gene_layers(marks, start, end, base_size, with_labels = with_labels && !gene_track,
                 label_facet = label_facet, facet_levels = group_label,
                 avoid = if(show_mcr) mcr else NULL)
  }
  # minimal common region drawn as an outline on top of the bars
  mcr_outline = function(d, ymin, ymax){
    d$ymin = ymin; d$ymax = ymax
    geom_rect(data = d, aes(xmin = start / 1e6, xmax = end / 1e6,
                            ymin = ymin, ymax = ymax),
              fill = NA, colour = mcr_colour, linewidth = 0.5,
              inherit.aes = FALSE)
  }
  two_mirrored = mirror && length(groups) == 2
  freq_data = frequencies[frequencies$n > 0, ]
  freq_data$type = factor(freq_data$type, levels = shown)
  freq_data = freq_data[order(freq_data$group, freq_data$start, freq_data$type), ]
  freq_data = dplyr::mutate(dplyr::group_by(freq_data, group, start),
                            ymax = cumsum(fraction), ymin = ymax - fraction)
  freq_data = as.data.frame(freq_data)
  mcr_plot = NULL
  if(show_mcr && !is.null(mcr)){
    mcr_plot = do.call(rbind, lapply(groups, function(g){
      m = mcr_of(g)
      if(is.null(m) || nrow(m) == 0) return(NULL)
      m$group = g
      m
    }))
  }
  if(!is.null(mcr_plot)){
    # keep narrow regions visible
    min_width = 0.005 * (end - start)
    centre = (mcr_plot$start + mcr_plot$end) / 2
    narrow = mcr_plot$end - mcr_plot$start < min_width
    mcr_plot$start[narrow] = centre[narrow] - min_width / 2
    mcr_plot$end[narrow] = centre[narrow] + min_width / 2
  }
  if(two_mirrored){
    flip = freq_data$group == groups[2]
    freq_data$ymin[flip] = -freq_data$ymin[flip]
    freq_data$ymax[flip] = -freq_data$ymax[flip]
    y_max = max(abs(c(freq_data$ymin, freq_data$ymax)), 0.05)
    p_freq = ggplot(freq_data) +
      geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = ymin,
                    ymax = ymax, fill = type))
    if(!is.null(mcr_plot)){
      p_freq = p_freq + mcr_outline(mcr_plot,
                                    ifelse(mcr_plot$group == groups[2], -y_max, 0),
                                    ifelse(mcr_plot$group == groups[2], 0, y_max))
    }
    p_freq = p_freq +
      geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
      annotate("text", x = x_lim[1], y = c(y_max, -y_max),
               label = group_label[groups], hjust = -0.02,
               vjust = c(1.2, -0.4), size = base_size / 3) +
      scale_y_continuous(limits = c(-y_max, y_max),
                         labels = function(v) paste0(round(abs(v) * 100), "%"))
  }else{
    freq_data$group = factor(group_label[freq_data$group], levels = group_label)
    p_freq = ggplot(freq_data) +
      geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = ymin,
                    ymax = ymax, fill = type))
    if(!is.null(mcr_plot)){
      mcr_facet = mcr_plot
      mcr_facet$group = factor(group_label[mcr_facet$group], levels = group_label)
      p_freq = p_freq + mcr_outline(mcr_facet, 0, Inf)
    }
    p_freq = p_freq +
      facet_grid(group ~ ., drop = FALSE, switch = "y") +
      scale_y_continuous(labels = function(v) paste0(round(v * 100), "%"),
                         expand = expansion(mult = c(0, 0.1)))
  }
  p_freq = .behind(p_freq, gene_layers(TRUE, if(two_mirrored) NULL else group_label[[1]])) +
    x_scale + fill_scale +
    labs(x = NULL, y = frequency_label) +
    theme_classic(base_size = base_size) +
    theme(legend.position = "top", axis.text.x = element_blank(),
          strip.background = element_blank(), strip.placement = "outside",
          strip.text.y.left = element_text(angle = 90))

  tile_data = events[!is.na(events$row), , drop = FALSE]
  tile_data$group = factor(group_label[tile_data$group], levels = group_label)
  # every group's panel at least min_rows tall, so its label fits
  min_rows = data.frame(group = factor(rep(group_label, each = 2), levels = group_label),
                        row = rep(c(1, min_rows), length(group_label)))
  p_tile = ggplot(tile_data) +
    geom_blank(data = min_rows, aes(y = row), inherit.aes = FALSE) +
    geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = row - 0.42,
                  ymax = row + 0.42, fill = type))
  p_tile = .behind(p_tile, gene_layers(FALSE))
  if(!is.null(mcr_plot)){
    mcr_tile = mcr_plot
    mcr_tile$group = factor(group_label[mcr_tile$group], levels = group_label)
    p_tile = p_tile + mcr_outline(mcr_tile, -Inf, Inf)
  }
  p_tile = p_tile +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", drop = FALSE,
               switch = "y") +
    x_scale + fill_scale +
    scale_y_continuous(breaks = NULL, expand = expansion(add = 0.6)) +
    labs(x = paste0("chr", chrom, " (Mb)"), y = tiling_label) +
    theme_classic(base_size = base_size) +
    theme(legend.position = "none", strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 90),
          panel.spacing = unit(2, "mm"))

  # each group's samples by their (most severe) event in the minimal common
  # region, as a share of the group, on one scale
  mcr_summary = NULL
  p_bar = NULL
  if(!is.null(mcr)){
    mcr_summary = do.call(rbind, lapply(groups, function(g){
      m = mcr_of(g)
      e = events[events$group == g, , drop = FALSE]
      hit = rep(FALSE, nrow(e))
      for(k in seq_len(nrow(m))){
        hit = hit | (e$start < m$end[k] & e$end > m$start[k])
      }
      e = e[hit, , drop = FALSE]
      worst = if(nrow(e) == 0) character(0) else
        tapply(match(as.character(e$type), severity), e$ID,
               function(v) severity[min(v)])
      counts = table(factor(worst, levels = shown))
      # denominator: the group, or its samples with data in the region
      n_group = n_samples[[g]]
      if(missing_data == "unknown"){
        d = seg[seg$group == g & !seg$filler, , drop = FALSE]
        has = rep(FALSE, nrow(d))
        for(k in seq_len(nrow(m))){
          has = has | (d$start < m$end[k] & d$end > m$start[k])
        }
        n_group = length(unique(d$ID[has]))
      }
      data.frame(group = g, event = c(shown, "none"),
                 n = c(as.numeric(counts), n_group - sum(counts)),
                 n_samples = n_group, stringsAsFactors = FALSE)
    }))
    mcr_summary$fraction = mcr_summary$n / mcr_summary$n_samples
    # share among the group's samples with an event there
    with_event = tapply(mcr_summary$n * (mcr_summary$event != "none"),
                        mcr_summary$group, sum)
    mcr_summary$n_with_event = as.numeric(with_event[mcr_summary$group])
    mcr_summary$fraction_of_events = ifelse(
      mcr_summary$event == "none" | mcr_summary$n_with_event == 0, NA_real_,
      mcr_summary$n / mcr_summary$n_with_event)
    if(show_mcr_summary){
      bar = mcr_summary
      bar$group = factor(bar$group, levels = groups)
      if(mcr_summary_include_none){
        n_shown = tapply(mcr_summary$n_samples, mcr_summary$group, `[`, 1)[groups]
        bar_label = data.frame(group = groups,
                               label = paste0(groups, " (n=", n_shown, ")"))
      }else{
        # composition among samples with an event, labelled with their number
        bar = bar[bar$event != "none", , drop = FALSE]
        bar$fraction = ifelse(is.na(bar$fraction_of_events), 0, bar$fraction_of_events)
        n_with = tapply(bar$n_with_event, as.character(bar$group), `[`, 1)[groups]
        bar_label = data.frame(group = groups,
                               label = paste0(groups, " (n=", n_with, ")"))
      }
      bar_label$group = factor(bar_label$group, levels = groups)
      bar$event = factor(bar$event, levels = c("none", rev(shown)))
      # no axes: group names (and numbers) run up inside each bar
      p_bar = ggplot(bar, aes(x = group, y = fraction, fill = event)) +
        geom_col(width = 0.85, colour = "grey60", linewidth = 0.2) +
        # on a translucent white background, readable over dark segments
        geom_label(data = bar_label, aes(x = group, y = 0.02, label = label),
                   inherit.aes = FALSE, angle = 90, hjust = 0,
                   size = base_size / 3.4, colour = "grey10",
                   fill = grDevices::adjustcolor("white", alpha.f = 0.75),
                   label.size = 0, label.padding = unit(0.6, "mm")) +
        scale_fill_manual(values = c(default_colours[shown], none = "white"),
                          guide = "none") +
        # zoom rather than limit the scale: a stack summing to 1 + rounding
        # error would otherwise lose its top segment
        coord_cartesian(ylim = c(0, 1), expand = FALSE) +
        theme_void(base_size = base_size) +
        theme(plot.margin = margin(0, 2, 0, 4))
    }
  }

  # samples with each track's data: given, the metadata's availability
  # column, or every plotted sample
  track_samples = function(given, column){
    if(!is.null(given)) return(given)
    if(!is.null(these_samples_metadata) && column %in% colnames(these_samples_metadata)){
      return(these_samples_metadata$sample_id[these_samples_metadata[[column]] %in% TRUE])
    }
    unique(seg$ID)
  }
  # breakpoint (SV) and mutation (SSM) density tracks, built as by
  # plot_sv_density() and plot_ssm_density()
  density_track = function(points, profiled, args, y_title, all_ids = points$ID){
    a = utils::modifyList(list(measure = "samples", bin_size = NULL, bin_fraction = 0.01,
                               max_bins = 2000, colour = NULL, group_colours = NULL,
                               class_colours = NULL, y_scale = "shared",
                               outliers = "gap", outlier_ratio = 3, outlier_k = 1.5,
                               outlier_min_events = 20, outlier_scope = "window"),
                          args)
    profiled = intersect(as.character(profiled), unique(seg$ID))
    if(length(profiled) == 0){
      warning("no plotted samples have ", y_title, " data; that track is left out")
      return(NULL)
    }
    units = .region_units(profiled, these_samples_metadata, split_by, groups, a$measure)
    units$groups = groups
    outlier_table = .find_outliers(points$ID[points$ID %in% units$map$ID],
                                   all_ids[all_ids %in% units$map$ID], a$outliers,
                                   a$outlier_k, a$outlier_min_events, a$outlier_scope,
                                   a$outlier_ratio)
    dropped = .drop_outliers(units, points, outlier_table, y_title)
    units = dropped$units; points = dropped$points
    breaks = .region_bins(start, end, a$bin_size, a$bin_fraction, a$max_bins)
    density = .binned_density(points, units, breaks, a$measure)
    n_units = tapply(density$n_units, density$group, `[`, 1)[groups]
    n_units[is.na(n_units)] = 0
    axis_title = c("SV breakpoints" = "SV frequency", "SSMs" = "SSM frequency")[[y_title]]
    y = paste0(axis_title, "\n(n=", paste(n_units, collapse = "/"), ")")
    p = .density_track(density, groups, group_label, start, end, marks, y_label = y,
                       measure = a$measure, fill = if(is.null(a$colour)) "grey35" else a$colour,
                       group_colours = if(is.null(a$colour)) .group_colours(groups, these_samples_metadata,
                                                                            split_by, a$group_colours),
                       class_colours = a$class_colours, y_scale = a$y_scale,
                       mirror = mirror, base_size = base_size, gene_labels = FALSE,
                       strips = FALSE) +
      theme(axis.text.x = element_blank(), axis.title.y = element_text(size = base_size * 0.8),
            legend.position = "inside", legend.position.inside = c(0.995, 0.98),
            legend.justification = c(1, 1), legend.direction = "horizontal",
            legend.key.size = unit(2.5, "mm"), legend.text = element_text(size = base_size * 0.7),
            legend.background = element_blank())
    list(plot = p, density = density, outliers = outlier_table)
  }
  sv_track = NULL
  if(!is.null(sv_data)){
    points = .sv_breakpoints(sv_data, chrom, start, end)
    classes = .sv_classes(points, sv_args$colour_by, chrom)
    if(!is.null(classes$class)) points$class = classes$class
    if(!is.null(classes$priority)) points$priority = classes$priority
    if(is.null(sv_args$class_colours) && !is.null(classes$colours)){
      sv_args$class_colours = classes$colours
    }
    sv_points = points
    points$type = NULL; points$partner = NULL
    sv_track = density_track(points,
                             track_samples(sv_samples, "sv_available"),
                             sv_args[setdiff(names(sv_args), "colour_by")], "SV breakpoints",
                             all_ids = .sv_sample_ids(sv_data))
  }
  ssm_track = NULL
  if(!is.null(ssm_data)){
    points = .ssm_points(ssm_data, chrom, start, end, ssm_args$colour_by,
                         if(is.null(ssm_args$max_classes)) 8 else ssm_args$max_classes)
    ssm_points = points
    ssm_track = density_track(points,
                              track_samples(ssm_samples, "ssm_available"),
                              ssm_args[setdiff(names(ssm_args), c("colour_by", "max_classes"))],
                              "SSMs", all_ids = .ssm_sample_ids(ssm_data))
  }
  sv_density = if(is.null(sv_track)) NULL else sv_track$density
  track_units = function(profiled, args){
    a = utils::modifyList(list(measure = "samples"), args)
    .region_units(intersect(as.character(profiled), unique(seg$ID)),
                  these_samples_metadata, split_by, groups,
                  if(a$measure == "mutations") "samples" else a$measure)
  }
  sv_events = if(is.null(sv_track)) NULL else
    .contributing_rows(sv_data, sv_points[!sv_points$ID %in% sv_track$outliers$ID, , drop = FALSE],
                       track_units(track_samples(sv_samples, "sv_available"), sv_args),
                       extra = list(breakpoint_in_window = "pos"))
  ssm_events = if(is.null(ssm_track)) NULL else
    .contributing_rows(ssm_data, ssm_points[!ssm_points$ID %in% ssm_track$outliers$ID, , drop = FALSE],
                       track_units(track_samples(ssm_samples, "ssm_available"), ssm_args))
  ssm_density = if(is.null(ssm_track)) NULL else ssm_track$density
  # without the tiling plot, the lowest track gets the x axis
  if(!show_tiling){
    with_axis = function(p) p + labs(x = paste0("chr", chrom, " (Mb)")) +
      theme(axis.text.x = element_text(), axis.title.x = element_text())
    if(!is.null(ssm_track)){
      ssm_track$plot = with_axis(ssm_track$plot)
    }else if(!is.null(sv_track)){
      sv_track$plot = with_axis(sv_track$plot)
    }else{
      p_freq = with_axis(p_freq)
    }
  }
  # left column: frequency track, SV and SSM tracks (optional), tiling plot
  # (optional)
  left_plots = c(list(p_freq), if(!is.null(sv_track)) list(sv_track$plot),
                 if(!is.null(ssm_track)) list(ssm_track$plot), if(show_tiling) list(p_tile))
  left_heights = c(rel_heights[1], if(!is.null(sv_track)) sv_height,
                   if(!is.null(ssm_track)) ssm_height, if(show_tiling) rel_heights[2])
  genes_plot = if(gene_track) .gene_track(marks, start, end, base_size, w$exons) else NULL
  track_height = if(is.null(genes_plot)) 0 else .gene_track_height(genes_plot$rows) * rel_heights[1]

  if(!is.null(genes_plot)){
    left_plots = c(list(genes_plot$plot), left_plots)
    left_heights = c(track_height, left_heights)
  }
  if(is.null(p_bar)){
    combined = cowplot::plot_grid(plotlist = left_plots, ncol = 1, align = "v",
                                  axis = "lr", rel_heights = left_heights)
  }else{
    legend_right = legend_position == "right"
    legend = NULL
    if(legend_right){
      # the legend, stacked vertically, goes below the bar chart
      legend = cowplot::get_plot_component(
        p_freq + theme(legend.position = "right", legend.direction = "vertical",
                        legend.justification = c(0, 0.5),
                        legend.key.spacing.y = unit(1, "mm")),
        "guide-box-right", return_all = FALSE)
      p_freq = p_freq + theme(legend.position = "none")
    }
    # the bar chart beside the frequency track: padded to match its panel's
    # top and bottom, in a column with the same row heights
    left_plots[[if(is.null(genes_plot)) 1 else 2]] = p_freq
    left = cowplot::plot_grid(plotlist = left_plots, ncol = 1, align = "v",
                              axis = "lr", rel_heights = left_heights)
    bar_aligned = cowplot::align_plots(p_freq, p_bar, align = "h", axis = "tb")[[2]]
    k = if(is.null(genes_plot)) 0 else 1
    legend_panel = if(legend_right) cowplot::ggdraw() +
      cowplot::draw_grob(legend, x = 0.02, y = 0.5, width = 0.96, height = 0.96,
                         hjust = 0, vjust = 0.5) else NULL
    # nothing below the frequency track (no tiling plot or density tracks):
    # the legend goes beside the bar chart rather than below it
    legend_beside = legend_right && length(left_heights) == k + 1
    side_width = if(legend_beside) mcr_summary_width + legend_width
                 else if(legend_right) max(legend_width, mcr_summary_width) else mcr_summary_width
    top_right = cowplot::plot_grid(bar_aligned, if(legend_beside) legend_panel, nrow = 1,
                                   rel_widths = c(mcr_summary_width,
                                                  side_width - mcr_summary_width))
    # the legend centred vertically in the space below the bar chart
    bottom_right = if(legend_right && !legend_beside) legend_panel else NULL
    # rows: (gene track), frequency track, the rest
    right = if(legend_beside){
      if(k == 1) cowplot::plot_grid(NULL, top_right, ncol = 1, rel_heights = left_heights)
      else top_right
    }else if(k == 1){
      cowplot::plot_grid(NULL, top_right, bottom_right, ncol = 1,
                         rel_heights = c(left_heights[1], left_heights[2],
                                         sum(left_heights[-(1:2)])))
    }else{
      cowplot::plot_grid(top_right, bottom_right, ncol = 1,
                         rel_heights = c(left_heights[1], sum(left_heights[-1])))
    }
    combined = cowplot::plot_grid(left, right, nrow = 1,
                                  rel_widths = c(1, side_width))
  }
  if(!return_data){
    return(combined)
  }
  # every segment in the region with its classification
  data = seg
  data$size = data$end - data$start
  data$shown = !is.na(data$type) & data$type %in% shown
  data$event = ifelse(data$shown, unname(labels[data$type]), NA_character_)
  data$in_mcr = FALSE
  for(g in groups){
    m = mcr_of(g)
    if(is.null(m)) next
    ix = which(data$group == g)
    for(k in seq_len(nrow(m))){
      data$in_mcr[ix] = data$in_mcr[ix] | (data$start[ix] < m$end[k] & data$end[ix] > m$start[k])
    }
  }
  data$row = events$row[match(paste(data$ID, data$start), paste(events$ID, events$start))]
  data = data[, intersect(c("ID", "group", "chrom", "start", "end", "size",
                            "CN", "log.ratio", "baseline_cn", "LOH_flag",
                            "type", "event", "shown", "in_mcr", "row"),
                          colnames(data))]
  rownames(data) = NULL
  tiling_order = sample_order
  if(!is.null(tiling_order)){
    tiling_order = tiling_order[order(tiling_order$group, -tiling_order$row), ]
    rownames(tiling_order) = NULL
  }
  # the input segments behind the plot: events of a shown type in the window
  # (original, unclipped rows), with their group, type and legend label
  input_df = if(is.list(seg_data) && !is.data.frame(seg_data)) as.data.frame(seg_data$segments)
             else as.data.frame(seg_data)
  ev = seg[!is.na(seg$type) & seg$type %in% shown, , drop = FALSE]
  ev = ev[!duplicated(ev$input_row), , drop = FALSE]
  events_out = input_df[ev$input_row, , drop = FALSE]
  events_out$group = ev$group
  events_out$baseline_cn = ev$baseline_cn
  events_out$type = ev$type
  events_out$event = unname(labels[ev$type])
  rownames(events_out) = NULL
  list(plot = combined, region = w$region, data = data, events = events_out,
       sv_events = sv_events, ssm_events = ssm_events,
       sv_outliers = if(is.null(sv_track)) NULL else sv_track$outliers,
       ssm_outliers = if(is.null(ssm_track)) NULL else ssm_track$outliers,
       mcr = mcr, mcr_summary = mcr_summary,
       sv_density = sv_density, ssm_density = ssm_density,
       tiling_order = tiling_order,
       frequencies = frequencies, frequency_plot = p_freq,
       tiling_plot = p_tile, mcr_summary_plot = p_bar)
}

# Segments with ID, chrom, start, end and CN and/or log.ratio (and LOH_flag
# if known), from generic segmented data or summarize_battenberg() output.
.common_region_segments = function(seg_data){
  if(is.list(seg_data) && !is.data.frame(seg_data)){
    if(is.null(seg_data$segments)){
      stop("seg_data must be segmented data or summarize_battenberg() output")
    }
    s = seg_data$segments
    # unresolved (unassessable) segments are left out as gaps
    s$input_row = seq_len(nrow(s))
    s = s[!is.na(s$total_cn) & !is.na(s$homologues), , drop = FALSE]
    return(data.frame(input_row = s$input_row, ID = s$sample, chrom = as.character(s$chrom),
                      start = s$start, end = s$end, CN = s$total_cn,
                      LOH_flag = as.integer(s$homologues != "both"),
                      stringsAsFactors = FALSE))
  }
  s = as.data.frame(seg_data)
  needed = c("ID", "chrom", "start", "end")
  if(!all(needed %in% colnames(s))){
    stop("seg_data needs columns ID, chrom, start and end (and CN or log.ratio)")
  }
  if(!any(c("CN", "log.ratio") %in% colnames(s))){
    stop("seg_data needs a CN or log.ratio column")
  }
  s$input_row = seq_len(nrow(s))
  s$ID = as.character(s$ID)
  s$chrom = gsub("^chr", "", as.character(s$chrom))
  s$start = as.numeric(s$start); s$end = as.numeric(s$end)
  s = s[order(s$ID, s$chrom, s$start), , drop = FALSE]
  # overlapping segments within a sample (e.g. two profiles under one ID,
  # such as genome and capture) would be double counted and mis-ordered
  n = nrow(s)
  overlap = n > 1 & c(FALSE, s$ID[-1] == s$ID[-n] & s$chrom[-1] == s$chrom[-n] &
                        s$start[-1] < s$end[-n])
  if(any(overlap)){
    ids = unique(s$ID[overlap])
    warning(length(ids), " sample(s) have overlapping segments, e.g. more than ",
            "one copy number profile per ID (check seq_type); they are double ",
            "counted and their tiling rows overlap: ",
            paste(utils::head(ids, 10), collapse = ", "),
            if(length(ids) > 10) ", ..." else "", call. = FALSE)
  }
  s
}

# Each segment's sample baseline: the segments' own baseline_cn of
# baseline_method (e.g. computed genome-wide for region segments by
# GAMBLR.results::get_cn_segments), or estimated from whole-genome segments
# with GAMBLR.utils::add_segment_baseline; filler rows never count
.sample_baselines = function(seg, baseline_method, genome_build){
  have = if("baseline_method" %in% colnames(seg)) unique(stats::na.omit(seg$baseline_method))
  if(length(have) == 1 && have != baseline_method &&
     !(have == "mode" && baseline_method == "mode_rounded")){
    stop("the segments' baselines (baseline_cn) are \"", have, "\", not \"", baseline_method,
         "\": use baseline_method = \"", have, "\", or get segments with that baseline, e.g. ",
         "GAMBLR.results::get_cn_segments(region = ..., baseline_method = \"", baseline_method,
         "\")", call. = FALSE)
  }
  s = seg
  s$dummy_segment = as.integer(seg$filler)
  s = GAMBLR.utils::add_segment_baseline(s, method = baseline_method, genome_build = genome_build)
  as.numeric(s$baseline_cn)
}

# Event type of every segment (NA if none), relative to each sample's
# baseline (see .sample_baselines) or CN 2.
.common_region_types = function(seg, baseline, log_ratio_threshold,
                                homdel_log_ratio, amp_threshold = 2,
                                amp_log_ratio = 1, genome_build = "grch37",
                                baseline_method = "mode"){
  has_cn = "CN" %in% colnames(seg) && any(!is.na(seg$CN))
  has_flag = "LOH_flag" %in% colnames(seg) && any(!is.na(seg$LOH_flag))
  # NA: LOH unknown (no LOH_flag for this segment)
  loh = if(has_flag) seg$LOH_flag == 1 else rep(FALSE, nrow(seg))
  if(has_cn){
    cn_raw = as.numeric(seg$CN)
    if(baseline == "sample"){
      base = .sample_baselines(seg, baseline_method, genome_build)
    }else{
      base = rep(2, nrow(seg))
    }
    # copy number relative to the baseline, on a diploid scale: with a
    # whole-number baseline this is round(CN) compared with the baseline
    cn = round(cn_raw - base + 2)
    type = ifelse(is.na(cn), NA_character_,
            ifelse(round(cn_raw) == 0, "homdel",
            ifelse(cn < 2, ifelse(is.na(loh), "deletion_loh_unknown",
                                  ifelse(loh, "deletion_loh", "deletion")),
            ifelse(cn >= 2 + amp_threshold, "amp",
            ifelse(cn > 2, ifelse(is.na(loh), "gain_loh_unknown",
                                  ifelse(loh, "gain_loh", "gain")),
            ifelse(!is.na(loh) & loh, "cnloh", NA_character_))))))
  }else{
    base = rep(NA_real_, nrow(seg))
    lr = as.numeric(seg$log.ratio)
    type = ifelse(is.na(lr), NA_character_,
            ifelse(lr <= homdel_log_ratio, "homdel",
            ifelse(lr < -log_ratio_threshold, "deletion",
            ifelse(lr >= amp_log_ratio, "amp",
            ifelse(lr > log_ratio_threshold, "gain", NA_character_)))))
  }
  # placeholder rows are never events
  type[seg$filler] = NA_character_
  list(type = type, baseline = base)
}

# Placeholder rows: diploid rows (CN 2 or log.ratio 0) at exact
# coordinates shared by min_samples or more samples.
.filler_segments = function(seg, min_samples){
  diploid = if("CN" %in% colnames(seg)) !is.na(seg$CN) & seg$CN == 2
            else !is.na(seg$log.ratio) & seg$log.ratio == 0
  key = paste(seg$chrom, seg$start, seg$end)
  shared = key %in% names(which(tapply(seg$ID, key, function(v) length(unique(v))) >= min_samples))
  diploid & shared
}

