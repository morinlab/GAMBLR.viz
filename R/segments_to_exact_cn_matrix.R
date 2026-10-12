#' @title Segmented data to an exact (breakpoint-resolution) CN matrix
#'
#' @description Convert segmented CN data into a sample x region matrix the
#' same way [segmented_data_to_cn_matrix()] does, except the column grid is
#' derived from the cohort's own segment breakpoints instead of a fixed
#' number of equal-width bins. This is the slow, run-once-or-a-few-times
#' step; its output is meant to be reused across many calls to
#' [pretty_CN_cohort_heatmap()] while you tweak the plot.
#'
#' @details The algorithm (all vectorized, no per-region loop over the full
#' segment table -- that nested-scan pattern is what makes the existing
#' bin-count-driven approach slow at this scale):
#' \enumerate{
#'   \item Build a `GRanges` from every segment in `seg_data`. If
#'     `merge_tolerance` > 0, round each breakpoint to the nearest multiple
#'     of it first -- this is what keeps the column count tractable, at the
#'     cost of collapsing breakpoints that fall within `merge_tolerance` bp
#'     of each other. `merge_tolerance = 0` gives true single-bp resolution
#'     (likely impractical at cohort scale; see Examples).
#'   \item `IRanges::disjoin()` the (optionally rounded) ranges once, cohort
#'     wide. This is the breakpoint-exact column grid: the minimal set of
#'     non-overlapping intervals such that no input segment boundary falls
#'     strictly inside one of them (the same idea `CGHregions` uses).
#'   \item `GenomicRanges::findOverlaps()` the segments against those bins
#'     once, cohort-wide, and resolve each (bin, sample) pair to a CN value
#'     (overlap-width-weighted average if more than one segment from the
#'     same sample touches a bin, which can happen after rounding).
#'   \item Pivot to a wide matrix with `data.table::dcast()` (much faster
#'     than `tidyr::pivot_wider()` at this width) and fill gaps per
#'     `fill_missing_with`, mirroring [segmented_data_to_cn_matrix()].
#' }
#'
#' @param seg_data A data frame of segments with columns `ID`, `chrom`,
#' `start`, `end`, `CN` (as returned by `get_cn_segments()`/
#' `get_sample_cn_segments()`).
#' @param these_samples_metadata Optional metadata table (with `sample_id`)
#' to subset to and to guarantee a row for every listed sample, even ones
#' with zero segments. If missing, every `ID` in `seg_data` gets a row.
#' @param merge_tolerance Breakpoints within this many bp of each other are
#' collapsed before disjoining, capping the number of resulting columns.
#' Default 10000 (10kb). Set to 0 for true breakpoint-exact resolution --
#' see the Examples for why that's usually impractical at >1000 samples.
#' @param max_segments_per_sample Samples with more total segments than this
#' are flagged with a `warning()` naming them, as likely noisy/low-quality
#' CN profiles -- a high segment count often means the caller fit lots of
#' spurious small breakpoints rather than real CNAs, and such samples also
#' disproportionately widen the cohort-wide bin grid for no biological
#' benefit. This is a heads up only: nothing is excluded automatically, and
#' the choice of whether to drop them (from `seg_data`/
#' `these_samples_metadata`) before calling this function again is yours.
#' Default 500. Set to `Inf` to disable the check.
#' @param fill_missing_with How to fill (sample, bin) combinations with no
#' overlapping segment: "diploid" (CN = 2), "avg_ploidy" (that sample's own
#' width-weighted mean CN), or "nothing" (leave `NA`). Default "diploid".
#' @param adjust_for_ploidy Set to TRUE to recentre each sample's CN values
#' on 2 using its own baseline ploidy. The dominant CN level -- whichever
#' (rounded) CN level covers the most total genomic length -- is found the
#' same way a mode would be, rather than using the length-weighted mean like
#' [segmented_data_to_cn_matrix()]'s argument of the same name (the mean is
#' skewed by any sample with a substantial but non-dominant burden of gains
#' or losses, which over-corrects an otherwise-diploid sample's baseline).
#' The shift applied is then the continuous, length-weighted mean CN of the
#' segments within that dominant level, not the rounded level itself -- a
#' sample whose dominant level rounds to 2 but whose actual values there
#' cluster at, say, 1.9 (common for low-signal/"boring" samples, where the
#' caller's own neutral reference may be off by a similar amount) still gets
#' shifted back to a true 2.0, instead of silently receiving no correction.
#' See `ploidy_baseline_autosomes_only`.
#' @param ploidy_baseline_autosomes_only When `adjust_for_ploidy = TRUE`,
#' exclude chrX/chrY segments from the baseline (mode) calculation. Default
#' TRUE: sex chromosome copy number reflects patient sex, not tumour
#' biology, and would bias the baseline for roughly half of a mixed-sex
#' cohort. The resulting (autosome-derived) shift is still applied to X/Y
#' segments in the output, same as everywhere else.
#' @param max_CN_allowed CN values above this are capped to it. Default 6.
#' @param rounded Set to TRUE to round the final matrix to integers. Default
#' FALSE -- unlike [segmented_data_to_cn_matrix()], this function defaults
#' to keeping the continuous values, since preserving them is the point.
#' @param verbose Default TRUE. Prints a timing breakdown for each step,
#' since the point of this function is to let you judge feasibility on your
#' own hardware before committing to it.
#'
#' @return A numeric matrix: `sample_id` as rows, `"chrom:start-end"` bin
#' names as columns (in genomic order), continuous CN values.
#'
#' @import dplyr data.table
#' @export
#'
#' @examples
#' \dontrun{
#' dlbcl_genome_meta = get_gambl_metadata() %>%
#'   dplyr::filter(pathology == "DLBCL", seq_type == "genome")
#' all_segments = get_cn_segments(these_samples_metadata = dlbcl_genome_meta)
#'
#' # a coarse first pass to see the shape/timing before going finer
#' cn_mat_coarse = segments_to_exact_cn_matrix(
#'   seg_data = all_segments,
#'   these_samples_metadata = dlbcl_genome_meta,
#'   merge_tolerance = 100000
#' )
#'
#' # true breakpoint resolution -- only try this after merge_tolerance = 1e5
#' # or 1e4 has told you roughly how the column count/runtime scales
#' cn_mat_exact = segments_to_exact_cn_matrix(
#'   seg_data = all_segments,
#'   these_samples_metadata = dlbcl_genome_meta,
#'   merge_tolerance = 0
#' )
#' }
segments_to_exact_cn_matrix = function(seg_data,
                                       these_samples_metadata,
                                       merge_tolerance = 10000,
                                       max_segments_per_sample = 500,
                                       fill_missing_with = "diploid",
                                       adjust_for_ploidy = FALSE,
                                       ploidy_baseline_autosomes_only = TRUE,
                                       max_CN_allowed = 6,
                                       rounded = FALSE,
                                       verbose = TRUE){
  required_cols = c("ID", "chrom", "start", "end", "CN")
  if(!all(required_cols %in% colnames(seg_data))){
    stop(paste(
      "seg_data must have columns:", paste(required_cols, collapse = ", ")
    ))
  }
  if(!fill_missing_with %in% c("diploid", "avg_ploidy", "nothing")){
    stop("fill_missing_with must be 'diploid', 'avg_ploidy' or 'nothing'")
  }

  t_start = Sys.time()
  tick = function(label, from){
    if(verbose){
      message(sprintf(
        "[segments_to_exact_cn_matrix] %s: %.2fs",
        label, as.numeric(difftime(Sys.time(), from, units = "secs"))
      ))
    }
  }

  if(!missing(these_samples_metadata)){
    all_samples = unique(these_samples_metadata$sample_id)
    seg_data = dplyr::filter(seg_data, ID %in% all_samples)
  }else{
    all_samples = unique(seg_data$ID)
  }
  seg_data = as.data.frame(seg_data)

  if(is.finite(max_segments_per_sample)){
    segment_counts = table(seg_data$ID)
    noisy_samples = names(segment_counts)[segment_counts > max_segments_per_sample]
    if(length(noisy_samples) > 0){
      warning(sprintf(
        paste(
          "%d sample(s) have more than %d total segments and likely",
          "represent noisy/low-quality CN profiles: %s.",
          "Consider excluding them before building the matrix",
          "(adjust with max_segments_per_sample)."
        ),
        length(noisy_samples), max_segments_per_sample,
        paste(noisy_samples, collapse = ", ")
      ))
    }
  }

  if(adjust_for_ploidy){
    t0 = Sys.time()
    seg_data = seg_data %>%
      dplyr::mutate(
        CN = ifelse(CN > max_CN_allowed, max_CN_allowed, CN),
        .size = end - start
      )

    # Baseline ploidy per sample: the length-weighted MODE of (rounded) CN,
    # not the length-weighted mean. The mean isn't robust to a skewed CN
    # distribution -- a sample that's genuinely diploid at baseline but
    # carries a substantial burden of gains (very common) gets its mean
    # pulled above 2 by those gains alone, even though most of its genome
    # (by length) is still truly CN=2; shifting by that inflated mean then
    # pushes the actually-neutral majority of the genome into "loss"
    # territory. The mode -- whichever (rounded) CN level covers the most
    # total genomic length -- is far more robust to that skew.
    #
    # Sex chromosomes are excluded from the baseline calculation (but not
    # from the matrix/plot -- the resulting autosome-derived shift is still
    # applied to X/Y) by default: X/Y copy number reflects patient sex, not
    # tumour biology, and would bias the baseline for roughly half of any
    # mixed-sex cohort.
    #
    # The dominant level is *identified* by rounding (robust to skew -- see
    # above), but the shift itself is computed from the continuous,
    # length-weighted mean CN of the segments inside that dominant bucket,
    # not from the rounded integer label. Using the rounded label directly
    # would silently skip correction for exactly the samples that most need
    # it: a low-signal/"boring" sample whose caller-assigned neutral CN is
    # consistently ~1.9 (not 2.0) still has its dominant bucket round to 2,
    # which would give a baseline of 2 and a shift of 0 -- leaving the whole
    # sample under-corrected even though adjust_for_ploidy ran.
    baseline_segs = seg_data
    if(ploidy_baseline_autosomes_only){
      baseline_segs = dplyr::filter(
        baseline_segs, !grepl("^(chr)?[XY]$", chrom, ignore.case = TRUE)
      )
    }
    baseline_df = baseline_segs %>%
      dplyr::mutate(CN_level = round(CN)) %>%
      dplyr::group_by(ID, CN_level) %>%
      dplyr::summarise(
        total_length = sum(.size),
        level_CN = stats::weighted.mean(CN, .size),
        .groups = "drop"
      ) %>%
      dplyr::arrange(CN_level) %>%
      dplyr::group_by(ID) %>%
      dplyr::slice_max(total_length, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(ID, .baseline = level_CN)

    seg_data = seg_data %>%
      dplyr::left_join(baseline_df, by = "ID") %>%
      dplyr::mutate(
        # a sample with no segments left after the autosome filter (shouldn't
        # normally happen) falls back to no adjustment rather than erroring
        .baseline = ifelse(is.na(.baseline), 2, .baseline),
        CN = 2 + CN - .baseline
      ) %>%
      dplyr::select(-.size, -.baseline)
    tick("ploidy adjustment (length-weighted mode baseline)", t0)
  }

  # 1. disjoint, breakpoint-exact (or merge_tolerance-exact) bins, cohort-wide
  t0 = Sys.time()
  chroms_present = unique(seg_data$chrom)
  has_prefix = any(grepl("^chr", chroms_present))
  chr_order = c(1:22, "X", "Y", "M")
  if(has_prefix) chr_order = paste0("chr", chr_order)
  chr_order = c(chr_order[chr_order %in% chroms_present],
               setdiff(chroms_present, chr_order))

  seg_start = seg_data$start
  seg_end = seg_data$end
  if(merge_tolerance > 0){
    seg_start = round(seg_start / merge_tolerance) * merge_tolerance
    seg_end = round(seg_end / merge_tolerance) * merge_tolerance
    seg_end = pmax(seg_end, seg_start + merge_tolerance)
  }
  seg_gr = GenomicRanges::GRanges(
    seqnames = factor(seg_data$chrom, levels = chr_order),
    ranges = IRanges::IRanges(start = pmax(seg_start, 1), end = seg_end)
  )
  bins_gr = GenomicRanges::sort(IRanges::disjoin(seg_gr))
  tick(sprintf("disjoin (%d bins from %d segments)", length(bins_gr), nrow(seg_data)), t0)

  # the output matrix is dense (every sample x every bin), so its cell
  # count -- not the input segment count -- is what actually drives memory
  # use in the steps below. Warn before attempting the expensive overlap
  # step rather than after it exhausts memory minutes in.
  projected_cells = length(bins_gr) * length(all_samples)
  if(verbose){
    message(sprintf(
      "[segments_to_exact_cn_matrix] projected matrix: %d samples x %d bins = %s cells",
      length(all_samples), length(bins_gr), format(projected_cells, big.mark = ",", scientific = FALSE)
    ))
  }
  if(projected_cells > 2e8){
    warning(sprintf(
      paste(
        "This will produce a dense matrix of ~%s cells (%d samples x %d bins).",
        "That's likely to be very slow and/or exhaust memory.",
        "Consider a larger merge_tolerance."
      ),
      format(projected_cells, big.mark = ","), length(all_samples), length(bins_gr)
    ))
  }

  # 2. resolve each (bin, sample) pair to a CN value in one vectorized pass
  t0 = Sys.time()
  hits = GenomicRanges::findOverlaps(bins_gr, seg_gr)
  overlap_width = IRanges::width(IRanges::pintersect(
    bins_gr[S4Vectors::queryHits(hits)], seg_gr[S4Vectors::subjectHits(hits)]
  ))
  hit_dt = data.table::data.table(
    bin_idx = S4Vectors::queryHits(hits),
    sample_id = seg_data$ID[S4Vectors::subjectHits(hits)],
    CN = seg_data$CN[S4Vectors::subjectHits(hits)],
    w = overlap_width
  )
  # the vast majority of (bin, sample) pairs have exactly one overlapping
  # segment -- those need no aggregation at all. Only a bin that spans more
  # than one of the *same* sample's segments (possible after merge_tolerance
  # rounding) needs resolving, and that subset is normally tiny. Doing a
  # by-group weighted.mean() over all pairs is enormously slower, since
  # weighted.mean() isn't one of data.table's GForce-optimized aggregates
  # and so gets called as a plain R function once per group -- in the
  # hundreds-of-millions-of-pairs regime that dominates runtime.
  hit_dt[, n := .N, by = list(bin_idx, sample_id)]
  if(all(hit_dt$n == 1L)){
    resolved = hit_dt[, list(bin_idx, sample_id, CN)]
  }else{
    singles = hit_dt[n == 1L, list(bin_idx, sample_id, CN)]
    multis = hit_dt[n > 1L]
    multis_resolved = multis[, list(CN = stats::weighted.mean(CN, w)),
                             by = list(bin_idx, sample_id)]
    resolved = data.table::rbindlist(list(singles, multis_resolved))
  }
  tick(sprintf("overlap + resolve (%d bin x sample pairs)", nrow(resolved)), t0)

  # 3. pivot to a wide sample x bin matrix
  t0 = Sys.time()
  bin_names = paste0(
    as.character(GenomicRanges::seqnames(bins_gr)), ":",
    IRanges::start(bins_gr), "-", IRanges::end(bins_gr)
  )
  resolved[, bin_name := bin_names[bin_idx]]
  wide = data.table::dcast(resolved, sample_id ~ bin_name, value.var = "CN")
  cn_matrix = as.matrix(wide[, -1, with = FALSE])
  rownames(cn_matrix) = wide$sample_id
  cn_matrix = cn_matrix[, bin_names, drop = FALSE]
  tick("pivot to wide matrix", t0)

  # 4. ensure every requested sample is a row (even ones with 0 segments),
  #    fill gaps, cap and (optionally) round
  t0 = Sys.time()
  missing_samples = setdiff(all_samples, rownames(cn_matrix))
  if(length(missing_samples) > 0){
    pad = matrix(NA_real_, nrow = length(missing_samples), ncol = ncol(cn_matrix),
                dimnames = list(missing_samples, colnames(cn_matrix)))
    cn_matrix = rbind(cn_matrix, pad)
  }
  cn_matrix = cn_matrix[all_samples, , drop = FALSE]

  if(fill_missing_with == "diploid"){
    cn_matrix[is.na(cn_matrix)] = 2
  }else if(fill_missing_with == "avg_ploidy"){
    row_means = rowMeans(cn_matrix, na.rm = TRUE)
    row_means[is.na(row_means)] = 2
    na_idx = which(is.na(cn_matrix), arr.ind = TRUE)
    if(nrow(na_idx) > 0){
      cn_matrix[na_idx] = row_means[na_idx[, 1]]
    }
  }
  # "nothing" -> leave NA as-is

  cn_matrix[cn_matrix > max_CN_allowed] = max_CN_allowed
  if(rounded){
    cn_matrix = round(cn_matrix)
  }
  tick("fill/finalize matrix", t0)
  tick(sprintf("TOTAL (%d samples x %d bins)", nrow(cn_matrix), ncol(cn_matrix)), t_start)

  return(cn_matrix)
}
