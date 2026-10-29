#' @title Copy number of several regions side by side
#'
#' @description Show the copy number of two or more (non-contiguous) regions
#' as a heatmap, one block of columns per region and one row per sample, in
#' the same sample order in every block, so that how a sample's copy number
#' in one region relates to another is visible. Samples can be ordered by
#' metadata, by their copy number in one of the regions, by clustering, or in
#' a given order.
#'
#' @details Each region is cut into bins of `bin_size` (the last one may be
#' shorter), so a region's block is as wide as it is long and every bin
#' covers the same length in every block; regions up to about 1 Mb at the
#' default 50 kb make 20 bins each. A bin's value is the length-weighted mean
#' copy number of the segments overlapping it (`bin_value = "mean"`), or the
#' value furthest from 2 (`"extreme"`), which keeps events smaller than a bin
#' (e.g. a focal homozygous deletion) visible. Bins with no segment are NA
#' (grey), except in controlfreec samples (a `module` column of
#' "controlfreec"), whose gaps are copy neutral and are filled with 2.
#'
#' Copy number is shown as given, or centred first (`center`, with
#' `GAMBLR.utils::center_segments`); segments already centred by it are
#' recognised. Colours are centred on 2.
#'
#' \strong{Sample order.} `sort_by` takes one or more keys, applied in turn
#' (later keys break ties), within the groups of `split_by`:
#' \itemize{
#'   \item a metadata column (in its factor or alphabetical order);
#'   \item the name of a region: by the sample's mean copy number in that
#'     region, highest first;
#'   \item "cluster": hierarchical clustering on all the regions' bins
#'     (as the last key; NA counts as 2).
#' }
#' Prefix a key with "-" to reverse it (e.g. "-MYC": lowest first). For any
#' other order, give `sample_order` (sample IDs, top to bottom); the
#' order used is returned with `return_data = TRUE`, so it can be reused.
#'
#' @param seg_data Segments with ID, chrom, start, end and CN (e.g. from
#' `GAMBLR.results::get_cn_segments`).
#' @param regions The regions, as a named list (names are the block titles).
#' Each element is a region string ("chr8:127000000-128000000"), one or more
#' gene symbols (the window spans them +/- `flank`), or a list with `genes`
#' and `flank` (or `region`). Unnamed elements are named after their genes
#' or coordinates.
#' @param these_samples_metadata Metadata (sample_id and any columns used by
#' `split_by`, `sort_by` or `annotation_columns`). Only its samples are
#' shown (those with segments in at least one of the regions; see
#' `keep_empty`), in its row order unless sorted. When it has a `cn_available`
#' column, only samples with TRUE are shown.
#' @param bin_size Bin size in bases. Default 5e4.
#' @param bin_value "mean" (default) or "extreme". See Details.
#' @param flank Bases added on each side of genes given as a region. Default
#' 2.5e5.
#' @param center NULL (default) shows CN as given; a baseline method
#' ("mode", "mode_rounded" or "mean") centres it first with
#' `GAMBLR.utils::center_segments` (not needed for segments it already
#' centred).
#' @param sort_by Sort keys (see Details). Default "cluster".
#' @param sample_order Optional sample IDs in the order to show them (top to
#' bottom); samples not listed follow. Overrides `sort_by`.
#' @param split_by Optional metadata column: samples are split into groups
#' (row blocks), sorted within each.
#' @param annotation_columns Metadata columns shown as annotation bars beside
#' the heatmap, replacing the default: `split_by` and the metadata columns
#' in `sort_by`.
#' @param additional_columns Further metadata columns to show as annotation
#' bars, after those of `annotation_columns` (or the default), e.g.
#' c("lymphgen", "pairing_status", "purity"). Categorical columns get GAMBL
#' colours where they exist; numeric ones a continuous scale.
#' @param marks Optional positions to mark with a dashed vertical line
#' through the heatmap (e.g. breakpoints): a data frame (or list of lists)
#' with chrom and pos, and optionally label (shown above). Marks outside
#' every region are ignored.
#' @param mark_genes Optional gene symbols whose start and end are marked
#' with dotted vertical lines, with the gene's name above its bins
#' (coordinates for `genome_build`).
#' @param sv_data Optional structural variants (bedpe-like, as from
#' `GAMBLR.results::get_combined_sv` or `get_manta_sv`: CHROM_A, START_A,
#' CHROM_B, START_B and a sample column). Each breakpoint in a region is
#' drawn as a point in its sample's row, where it is in the region: black
#' when the other end is in another of the regions (e.g. an IGH-BCL2
#' translocation), grey when both ends are in the same region, white when
#' the other end is elsewhere. An "SV" annotation column gives each
#' sample's most notable class.
#' @param sv_samples Optional IDs of the samples with SV data; by default
#' those with `sv_available` TRUE when the metadata has that column. Samples
#' without SV data are marked in an "SV_data" annotation column, since they
#' can't show any.
#' @param sv_point_size Size of the breakpoint points in mm. Default 1.5.
#' @param sv_show The breakpoint classes to show (points and the "SV"
#' annotation): any of "linked" (other end in another of the regions),
#' "within" (other end in the same region, e.g. class switch or V(D)J
#' deletions in IGH) and "other" (other end elsewhere). Default all three.
#' Either one set for every region, or a list by region name, e.g.
#' `list(BCL2 = c("linked", "other"), IGH = "linked")` (regions not in the
#' list show all three).
#' Linked breakpoints are drawn by the orientation of their junction (from
#' STRAND_A and STRAND_B, "+" keeping the side before the breakpoint, the
#' strands of the ends in the order of `regions`): an up triangle for "+-",
#' down for "-+", a diamond for "++" and a square for "--", the same for
#' both of a junction's points (other breakpoints are circles). A balanced
#' reciprocal translocation gives both complementary orientations (e.g. for
#' t(14;18) with BCL2 listed first, "+-" is der(18) and "-+" der(14)); a
#' `linked_junctions` annotation (also a metadata column) says whether each
#' sample has both ("reciprocal") or only one ("one orientation"). A message
#' gives what each shape means for the regions plotted.
#' The number of breakpoints shown in each region per sample is added to the
#' metadata as `<region>_SVs` columns, which can be used in `sort_by`, e.g.
#' `sort_by = c("-BCL2_SVs", "cluster")` for samples with BCL2 breakpoints
#' (of the classes shown) first.
#' To filter by SV type or size instead, subset `sv_data` first.
#' @param ssm_data Optional mutations (MAF: Tumor_Sample_Barcode,
#' Chromosome, Start_Position), as one data frame or a list of them (e.g.
#' one `GAMBLR.results::get_ssm_by_region()` call per region). The number of
#' mutations each sample has in each region is shown as a bar chart beside
#' the heatmap, one per region, and added to the metadata as
#' `<region>_SSMs` columns (which can be used in `sort_by`). Duplicated
#' mutations (same sample, position and alleles) count once.
#' @param ssm_samples Optional IDs of the samples whose mutation data cover
#' the regions; by default those with `ssm_available` TRUE when the
#' metadata has that column. Samples without are marked in an "SSM_data"
#' annotation column, since their count of 0 means no data (e.g. capture
#' samples outside their targets).
#' @param ssm_bar_width Width of each SSM bar chart in mm. Default 12.
#' @param concordance Optional output of `GAMBLR.utils::cn_step_concordance`
#' for these samples: its per-sample columns (step_type, pattern, a_pos,
#' b_pos, a_step, b_step and b_partner) are added to the metadata, only
#' samples with a step type in `step_types` are shown, and unless given,
#' `split_by` is "step_type" and `sort_by` c("a_pos", "b_pos", "cluster").
#' @param step_types With `concordance`, the step types to show (e.g.
#' c("BCL2 down, IGH up", "BCL2 up, IGH down")). Default: all but "none".
#' @param keep_empty Show samples with no segments in any of the regions
#' (as grey rows). Default FALSE: they are left out, with a message.
#' @param cn_max Copy number at which the colour scale saturates. Default 6.
#' @param max_bins Stop if the regions would make more bins than this in
#' total. Default 2000.
#' @param show_row_names Show sample IDs. Default: when there are at most
#' 100 samples.
#' @param genome_build "grch37" (default) or "hg38", for genes.
#' @param return_data Set to TRUE to also return the matrix, the sample
#' order and the regions.
#'
#' @return A ComplexHeatmap heatmap (drawn when printed), or with
#' `return_data = TRUE` a list: `plot` (the heatmap; also as `heatmap`),
#' the input rows behind it: `events` (the `seg_data` rows overlapping each
#' region for the samples shown, with a `region` column), `sv_events` (the
#' `sv_data` rows of the breakpoints shown, with `sv_class`, `sv_pair` and
#' `sv_orientation`) and `ssm_events` (the `ssm_data` rows counted, with a
#' `region` column); and `matrix` (samples by bins, in the order shown),
#' `sample_order`, `metadata` (the samples shown, in that order,
#' with the columns added: `<region>_SVs`, `<region>_SSMs`,
#' `linked_junctions`, step columns), `regions` (name, region, chrom, start,
#' end, bins) and `bins` (region, bin, chrom, start, end) and, with `ssm_data`, `ssm_counts`
#' (one row per sample shown, one column per region).
#'
#' @import ComplexHeatmap
#' @export
#'
#' @examples
#' \dontrun{
#' meta = get_gambl_metadata() %>% dplyr::filter(pathology == "DLBCL", seq_type == "genome")
#' meta = GAMBLR.results::annotate_cn_availability(meta)
#' seg = GAMBLR.results::get_cn_segments(dplyr::filter(meta, cn_available))
#' plot_multi_region_CN(seg,
#'                      regions = list(MYC = "MYC", BCL6 = "BCL6", REL = "REL",
#'                                     CDKN2A = list(genes = "CDKN2A", flank = 5e5)),
#'                      these_samples_metadata = meta, center = "mode",
#'                      split_by = "COO_consensus", sort_by = c("REL", "cluster"))
#' }
plot_multi_region_CN = function(seg_data,
                                regions,
                                these_samples_metadata,
                                bin_size = 5e4,
                                bin_value = "mean",
                                flank = 2.5e5,
                                center = NULL,
                                sort_by = "cluster",
                                sample_order = NULL,
                                split_by = NULL,
                                annotation_columns = NULL,
                                additional_columns = NULL,
                                marks = NULL,
                                mark_genes = NULL,
                                sv_data = NULL,
                                sv_samples = NULL,
                                sv_point_size = 1.5,
                                sv_show = c("linked", "within", "other"),
                                ssm_data = NULL,
                                ssm_samples = NULL,
                                ssm_bar_width = 12,
                                concordance = NULL,
                                step_types = NULL,
                                keep_empty = FALSE,
                                cn_max = 6,
                                max_bins = 2000,
                                show_row_names = NULL,
                                genome_build = "grch37",
                                return_data = FALSE){
  bin_value = match.arg(bin_value, c("mean", "extreme"))
  if(missing(these_samples_metadata)) stop("these_samples_metadata is required")
  meta = as.data.frame(these_samples_metadata)
  meta = meta[!duplicated(meta$sample_id), , drop = FALSE]
  # copy number steps (GAMBLR.utils::cn_step_concordance): their columns are
  # added to the metadata, samples are split by step type and sorted by the
  # transitions' positions unless told otherwise
  if(!is.null(concordance)){
    steps = as.data.frame(concordance$steps)
    step_cols = intersect(c("step_type", "pattern", "b_partner", "a_pos", "b_pos", "a_step", "b_step"),
                          colnames(steps))
    replaced = intersect(step_cols, colnames(meta))
    if(length(replaced) > 0) message("replacing metadata column(s) with the concordance's: ",
                                     paste(replaced, collapse = ", "))
    meta = meta[, setdiff(colnames(meta), step_cols), drop = FALSE]
    meta = merge(meta, steps[, c("sample_id", step_cols), drop = FALSE], by = "sample_id", sort = FALSE)
    if(is.null(step_types)){
      keep = !is.na(meta$step_type) & meta$step_type != "none"
    }else{
      keep = as.character(meta$step_type) %in% step_types
    }
    message(sum(keep), " of ", nrow(meta), " samples have a shown step type",
            if(is.null(step_types)) " (any but \"none\"; see step_types)" else "")
    meta = meta[keep, , drop = FALSE]
    meta$step_type = droplevels(factor(meta$step_type, levels = levels(steps$step_type)))
    if(is.null(split_by)) split_by = "step_type"
    if(missing(sort_by)) sort_by = c("a_pos", "b_pos", "cluster")
  }
  if("cn_available" %in% colnames(meta)){
    n_all = nrow(meta)
    meta = meta[meta$cn_available %in% TRUE, , drop = FALSE]
    if(nrow(meta) < n_all){
      message(nrow(meta), " of ", n_all, " samples have usable copy number data (cn_available) and are shown")
    }
  }
  seg = as.data.frame(seg_data)
  if(!all(c("ID", "chrom", "start", "end", "CN") %in% colnames(seg))){
    stop("seg_data needs columns ID, chrom, start, end and CN")
  }
  dropped = setdiff(unique(seg$ID), meta$sample_id)
  if(length(dropped) > 0){
    message("dropping ", length(dropped), " sample(s) not in these_samples_metadata")
  }
  seg = seg[seg$ID %in% meta$sample_id, , drop = FALSE]
  if(!is.null(center) && !"CN_uncentered" %in% colnames(seg)){
    center = match.arg(center, GAMBLR.utils::baseline_methods())
    seg = as.data.frame(GAMBLR.utils::center_segments(seg, method = center))
  }

  # the regions, and their bins
  wins = .multi_region_windows(regions, flank, genome_build)
  n_bins = ceiling((wins$end - wins$start) / bin_size)
  if(sum(n_bins) > max_bins){
    stop("these regions make ", sum(n_bins), " bins of ", bin_size, " bp (more than max_bins = ",
         max_bins, "): use a larger bin_size or smaller regions")
  }
  bins = do.call(rbind, lapply(seq_len(nrow(wins)), function(i){
    st = seq(wins$start[i], by = bin_size, length.out = n_bins[i])
    data.frame(region = wins$name[i], bin = seq_along(st), chrom = wins$chrom[i],
               start = st, end = pmin(st + bin_size, wins$end[i]), stringsAsFactors = FALSE)
  }))
  bins$column = paste0(bins$region, ":", bins$bin)
  wins$bins = n_bins

  # samples by bins
  ids = meta$sample_id
  mat = matrix(NA_real_, nrow = length(ids), ncol = nrow(bins), dimnames = list(ids, bins$column))
  seg$chrom = gsub("^chr", "", as.character(seg$chrom))
  for(i in seq_len(nrow(wins))){
    b = bins[bins$region == wins$name[i], , drop = FALSE]
    s = seg[seg$chrom == wins$chrom[i] & seg$end > wins$start[i] & seg$start < wins$end[i] &
              !is.na(seg$CN), , drop = FALSE]
    if(nrow(s) == 0) next
    # each segment against each bin it overlaps
    nb = nrow(b)
    first = pmax(1, pmin(nb, findInterval(pmax(s$start, wins$start[i]), b$start)))
    last = pmax(1, pmin(nb, findInterval(pmin(s$end, wins$end[i]) - 1, b$start)))
    k = last - first + 1
    row = rep(seq_len(nrow(s)), k)
    col = first[row] + sequence(k) - 1
    overlap = pmin(s$end[row], b$end[col]) - pmax(s$start[row], b$start[col])
    keep = overlap > 0
    row = row[keep]; col = col[keep]; overlap = overlap[keep]
    cn = as.numeric(s$CN[row])
    key = paste(s$ID[row], col, sep = "\r")
    value = if(bin_value == "mean"){
      tapply(cn * overlap, key, sum) / tapply(overlap, key, sum)
    }else{
      tapply(cn, key, function(v) v[which.max(abs(v - 2))])
    }
    parts = do.call(rbind, strsplit(names(value), "\r", fixed = TRUE))
    mat[cbind(match(parts[, 1], ids), match(b$column[as.integer(parts[, 2])], colnames(mat)))] =
      as.numeric(value)
  }
  # controlfreec doesn't report copy-neutral segments: its gaps are 2
  if("module" %in% colnames(seg)){
    modules = tapply(as.character(seg$module), seg$ID, function(m) paste(unique(m), collapse = "+"))
    neutral_gaps = names(modules)[modules == "controlfreec"]
    rows = ids %in% neutral_gaps
    filled = mat[rows, , drop = FALSE]
    filled[is.na(filled)] = 2
    mat[rows, ] = filled
  }
  no_data = ids[rowSums(!is.na(mat)) == 0]
  if(length(no_data) > 0){
    if(keep_empty){
      message(length(no_data), " sample(s) have no segments in any of the regions (shown grey)")
    }else{
      message(length(no_data), " sample(s) have no segments in any of the regions and are left out ",
              "(keep_empty = TRUE shows them)")
      keep = !ids %in% no_data
      mat = mat[keep, , drop = FALSE]
      meta = meta[keep, , drop = FALSE]
      ids = ids[keep]
    }
  }
  if(length(ids) == 0) stop("no samples have segments in any of the regions")

  # SSMs per sample in each region, as metadata columns <region>_SSMs (so
  # they can also be sort keys)
  ssm_cols = character(0)
  if(!is.null(ssm_data)){
    counts = .multi_region_ssm_counts(ssm_data, wins, ids)
    ssm_cols = paste0(wins$name, "_SSMs")
    for(k in seq_along(ssm_cols)) meta[[ssm_cols[k]]] = counts[, k]
    message(sum(counts), " SSM(s) in the regions, in ", sum(rowSums(counts) > 0), " of the samples shown")
  }

  # SV breakpoints in the regions (before ordering, so their counts can be
  # sort keys)
  sv_fill = c(linked = "black", within = "grey45", other = "white")
  svs = if(is.null(sv_data)) NULL else .multi_region_svs(sv_data, wins, bins, ids)
  # the classes shown in each region: one set for all, or a list by region
  # name (regions not in it show all classes)
  show_by_region = stats::setNames(lapply(wins$name, function(n){
    v = if(is.list(sv_show)) { if(n %in% names(sv_show)) sv_show[[n]] else names(sv_fill) } else sv_show
    match.arg(v, names(sv_fill), several.ok = TRUE)
  }), wins$name)
  if(is.list(sv_show)){
    unknown = setdiff(names(sv_show), wins$name)
    if(length(unknown) > 0) stop("sv_show names no region(s): ", paste(unknown, collapse = ", "))
  }
  sv_shown = intersect(names(sv_fill), unique(unlist(show_by_region)))
  if(!is.null(svs)){
    keep = mapply(function(region, link) link %in% show_by_region[[region]], svs$region, svs$link)
    if(any(!keep)) message(sum(!keep), " breakpoint(s) of the classes not in sv_show are hidden")
    svs = svs[keep, , drop = FALSE]
    if(nrow(svs) == 0) svs = NULL
  }
  # linked junctions per sample: both complementary orientations of a
  # region pair (reciprocal), or only one
  if(!is.null(sv_data)){
    junctions = rep("none", length(ids))
    lk = if(is.null(svs)) NULL else unique(svs[svs$link == "linked" & !is.na(svs$orientation),
                                               c("ID", "pair", "orientation"), drop = FALSE])
    if(!is.null(lk) && nrow(lk) > 0){
      flip = function(o) chartr("+-", "-+", o)
      by_sample = split(lk, lk$ID)
      cls = vapply(by_sample, function(x){
        recip = any(vapply(split(x$orientation, x$pair), function(o) any(flip(o) %in% o), logical(1)))
        if(recip) "reciprocal" else "one orientation"
      }, "")
      junctions[match(names(cls), ids)] = cls
      # what the shapes mean, for each region pair, with the numbers of SVs
      # (calls) and of samples of each orientation, and of samples with both
      # of a complementary pair
      sv_lk = unique(svs[svs$link == "linked" & !is.na(svs$orientation), c("ID", "sv_id", "pair", "orientation")])
      shape_name = c("+-" = "up triangle", "-+" = "down triangle", "++" = "diamond", "--" = "square")
      for(pr in unique(lk$pair)){
        rg = strsplit(pr, "-", fixed = TRUE)[[1]]
        x = sv_lk[sv_lk$pair == pr, , drop = FALSE]
        o = sort(unique(x$orientation))
        n_sv = vapply(o, function(k) length(unique(x$sv_id[x$orientation == k])), integer(1))
        n_smp = vapply(o, function(k) length(unique(x$ID[x$orientation == k])), integer(1))
        both = vapply(split(x$orientation, x$ID), function(v) any(chartr("+-", "-+", v) %in% v), logical(1))
        message("linked ", pr, " junctions (", rg[1], " strand/", rg[2], " strand): ",
                paste0(shape_name[o], " = ", substr(o, 1, 1), "/", substr(o, 2, 2), " (", n_sv,
                       " SVs in ", n_smp, " samples)", collapse = "; "),
                ". Samples with both orientations of a pair (reciprocal): ", sum(both), " of ",
                length(both))
      }
    }
    meta$linked_junctions = factor(junctions, levels = c("reciprocal", "one orientation", "none"))
    if(any(junctions != "none")) annotation_columns = unique(c(annotation_columns, "linked_junctions"))
  }
  # shown breakpoints per sample in each region, as metadata columns
  # <region>_SVs (e.g. sort_by = "-BCL2_SVs" puts SV-positive samples first)
  if(!is.null(sv_data)){
    for(n in wins$name){
      hits = if(is.null(svs)) character(0) else svs$ID[svs$region == n]
      meta[[paste0(n, "_SVs")]] = as.integer(table(factor(hits, levels = ids)))
    }
  }
  # sample order
  if(!is.null(split_by) && !split_by %in% colnames(meta)){
    stop("these_samples_metadata has no column ", split_by,
         if(split_by %in% c("step_type", "pattern", "b_partner"))
           " (it comes from GAMBLR.utils::cn_step_concordance(): pass its result as concordance)" else "",
         call. = FALSE)
  }
  groups = if(is.null(split_by)) factor(rep("all", length(ids))) else {
    g = as.character(meta[[split_by]])
    g[is.na(g)] = "NA"
    factor(g, levels = unique(c(.level_order(meta[[split_by]]), "NA")[c(.level_order(meta[[split_by]]), "NA") %in% g]))
  }
  ord = if(!is.null(sample_order)){
    given = intersect(as.character(sample_order), ids)
    c(match(given, ids), setdiff(seq_along(ids), match(given, ids)))
  }else{
    .multi_region_order(mat, meta, bins, sort_by, groups)
  }
  if(!is.null(split_by)) ord = ord[order(as.integer(groups[ord]))]
  mat = mat[ord, , drop = FALSE]
  meta = meta[ord, , drop = FALSE]
  groups = groups[ord]

  # annotations
  sort_keys = sub("^-", "", sort_by)
  if(is.null(annotation_columns)){
    annotation_columns = unique(c(split_by, intersect(sort_keys, colnames(meta))))
  }
  annotation_columns = unique(c(annotation_columns, additional_columns))
  missing_columns = setdiff(annotation_columns, colnames(meta))
  if(length(missing_columns) > 0){
    stop("these_samples_metadata has no column(s): ", paste(missing_columns, collapse = ", "))
  }
  # the samples' rows in the order shown
  if(!is.null(svs)) svs$row = match(svs$ID, rownames(mat))
  if(!is.null(sv_data)){
    message(if(is.null(svs)) 0 else nrow(svs), " SV breakpoint(s) in the regions, in ",
            if(is.null(svs)) 0 else length(unique(svs$row)), " of the samples shown")
    per_sample = rep("none", nrow(mat))
    if(!is.null(svs) && nrow(svs) > 0){
      rank = match(svs$link, names(sv_fill))
      best = tapply(rank, svs$row, min)
      per_sample[as.integer(names(best))] = names(sv_fill)[best]
    }
    meta$SV = factor(per_sample, levels = c(sv_shown, "none"))
    annotation_columns = unique(c(annotation_columns, "SV"))
  }
  # which samples have SV data (when known): those without can't show any
  if(!is.null(sv_data)){
    profiled = if(!is.null(sv_samples)) meta$sample_id %in% sv_samples
               else if("sv_available" %in% colnames(meta)) meta$sv_available %in% TRUE
               else NULL
    if(!is.null(profiled) && any(!profiled)){
      meta$SV_data = ifelse(profiled, "yes", "no")
      annotation_columns = unique(c(annotation_columns, "SV_data"))
    }
  }
  # which samples have SSM data covering the regions (when known)
  if(!is.null(ssm_data)){
    profiled = if(!is.null(ssm_samples)) meta$sample_id %in% ssm_samples
               else if("ssm_available" %in% colnames(meta)) meta$ssm_available %in% TRUE
               else NULL
    if(!is.null(profiled) && any(!profiled)){
      meta$SSM_data = ifelse(profiled, "yes", "no")
      annotation_columns = unique(c(annotation_columns, "SSM_data"))
    }
  }
  left = NULL
  if(length(annotation_columns) > 0){
    anno = meta[, annotation_columns, drop = FALSE]
    # factors keep their level order (in the legend too)
    anno[] = lapply(anno, function(v) if(is.numeric(v) || is.factor(v)) v else as.character(v))
    categorical = names(anno)[!vapply(anno, is.numeric, logical(1))]
    colours = if(length(categorical) > 0)
      map_metadata_to_colours(categorical, meta[, categorical, drop = FALSE]) else list()
    colours = colours[intersect(names(colours), categorical)]
    if("SV_data" %in% names(anno)) colours$SV_data = c(yes = "grey20", no = "grey90")
    if("SSM_data" %in% names(anno)) colours$SSM_data = c(yes = "grey20", no = "grey90")
    if("linked_junctions" %in% names(anno)){
      colours$linked_junctions = c(reciprocal = "black", "one orientation" = "#E08214", none = "white")
    }
    if("SV" %in% names(anno)) colours$SV = c(linked = "black", within = "grey45", other = "grey80",
                                             none = "white")[levels(meta$SV)]
    left = ComplexHeatmap::rowAnnotation(df = anno, col = colours,
                                         annotation_name_side = "top",
                                         annotation_name_gp = grid::gpar(fontsize = 8))
  }
  # genes marked above their bins
  gene_bins = rep("", nrow(bins))
  for(i in seq_len(nrow(wins))){
    g = wins$genes[[i]]
    if(is.null(g) || nrow(g) == 0) next
    for(j in seq_len(nrow(g))){
      hit = bins$region == wins$name[i] & bins$end > g$start[j] & bins$start < g$end[j]
      gene_bins[hit] = g$gene[j]
    }
  }
  # marks: positions (marks; dashed lines) and gene starts and ends
  # (mark_genes; dotted lines, the gene's name above its bins)
  mark_tab = data.frame(chrom = character(0), pos = numeric(0), label = character(0), lty = numeric(0))
  if(!is.null(marks)){
    m = if(is.data.frame(marks)) marks else do.call(rbind, lapply(marks, as.data.frame))
    mark_tab = rbind(mark_tab, data.frame(chrom = as.character(m$chrom), pos = as.numeric(m$pos),
                                          label = if("label" %in% colnames(m)) as.character(m$label)
                                                  else NA_character_, lty = 2))
  }
  label_at = integer(0); label_text = character(0)
  if(!is.null(mark_genes)){
    g = as.data.frame(suppressMessages(GAMBLR.utils::gene_to_region(
      gene_symbol = mark_genes, projection = genome_build, return_as = "df")))
    colnames(g)[seq_len(3)] = c("chrom", "start", "end")
    gene_col = grep("hugo|gene", colnames(g), ignore.case = TRUE)[1]
    g$gene = g[[gene_col]]
    g = g[!duplicated(g$gene), , drop = FALSE]
    not_found = setdiff(mark_genes, g$gene)
    if(length(not_found) > 0) message("not found, so not marked: ", paste(not_found, collapse = ", "))
    g$chrom = gsub("^chr", "", as.character(g$chrom))
    mark_tab = rbind(mark_tab, data.frame(chrom = rep(g$chrom, 2), pos = c(g$start, g$end),
                                          label = NA_character_, lty = 3))
    # each gene's name over the middle of the bins it overlaps
    for(j in seq_len(nrow(g))){
      k = which(bins$chrom == g$chrom[j] & bins$end > g$start[j] & bins$start < g$end[j])
      if(length(k) == 0) next
      label_at = c(label_at, k[ceiling(length(k) / 2)]); label_text = c(label_text, g$gene[j])
    }
  }
  mark_tab$chrom = gsub("^chr", "", mark_tab$chrom)
  hits = do.call(rbind, lapply(seq_len(nrow(mark_tab)), function(m){
    k = which(bins$chrom == mark_tab$chrom[m] & bins$start <= mark_tab$pos[m] & bins$end > mark_tab$pos[m])
    if(length(k) == 0) return(NULL)
    data.frame(col = k, frac = (mark_tab$pos[m] - bins$start[k]) / (bins$end[k] - bins$start[k]),
               lty = mark_tab$lty[m], label = mark_tab$label[m])
  }))
  if(!is.null(hits) && nrow(hits) > 0){
    labelled = hits[!is.na(hits$label), , drop = FALSE]
    label_at = c(label_at, labelled$col); label_text = c(label_text, labelled$label)
  }
  layer = NULL
  if((!is.null(hits) && nrow(hits) > 0) || (!is.null(svs) && nrow(svs) > 0)){
    layer = function(j, i, x, y, w, h, fill){
      if(!is.null(hits) && nrow(hits) > 0){
        for(col in intersect(unique(j), hits$col)){
          k = which(j == col)[1]
          for(r in which(hits$col == col)){
            xx = x[k] - w[k] * 0.5 + w[k] * hits$frac[r]
            grid::grid.lines(grid::unit.c(xx, xx), grid::unit(c(0, 1), "npc"),
                             gp = grid::gpar(col = if(hits$lty[r] == 3) "grey25" else "black",
                                             lwd = 0.8, lty = hits$lty[r]))
          }
        }
      }
      if(!is.null(svs) && nrow(svs) > 0){
        # the cell of each breakpoint in this slice, and where in the bin
        k = match(paste(svs$row, svs$col), paste(i, j))
        here = which(!is.na(k))
        if(length(here) > 0){
          kk = k[here]
          xx = x[kk] - w[kk] * 0.5 + w[kk] * svs$frac[here]
          grid::grid.points(xx, y[kk], pch = svs$pch[here], size = grid::unit(sv_point_size, "mm"),
                            gp = grid::gpar(col = "black", lwd = 0.6,
                                            fill = sv_fill[svs$link[here]]))
        }
      }
    }
  }
  top_parts = list(gene = ComplexHeatmap::anno_simple(ifelse(gene_bins == "", "none", "gene"),
                                                      col = c(none = "white", gene = "grey20"),
                                                      height = grid::unit(1.5, "mm")))
  if(length(label_at) > 0){
    # labels in the same bin are shown together
    label_text = tapply(label_text, label_at, function(v) paste(unique(v), collapse = "/"))
    label_at = as.integer(names(label_text)); label_text = unname(as.character(label_text))
    # each label centred over its bin (so it sits on its line), in two rows
    # where neighbours are close
    o = order(label_at); label_at = label_at[o]; label_text = label_text[o]
    row = rep(1L, length(label_at))
    for(k in seq_along(label_at)[-1]){
      close = label_at[k] - label_at[k - 1] < max(nchar(label_text[c(k, k - 1)])) * 0.6 * (ncol(mat) / 100)
      if(close && row[k - 1] == 1L) row[k] = 2L
    }
    label_rows = lapply(1:2, function(r){
      v = rep("", ncol(mat)); v[label_at[row == r]] = label_text[row == r]; v
    })
    text_row = function(v) ComplexHeatmap::anno_text(v, which = "column", rot = 0, just = "center",
                                                     location = grid::unit(0.5, "npc"),
                                                     gp = grid::gpar(fontsize = 8),
                                                     height = grid::unit(3.5, "mm"))
    top_parts = c(if(any(row == 2L)) list(labels2 = text_row(label_rows[[2]])),
                  list(labels = text_row(label_rows[[1]])), top_parts)
  }
  top = do.call(ComplexHeatmap::HeatmapAnnotation, c(top_parts, list(show_annotation_name = FALSE)))
  # SSM counts per region: a bar chart beside the heatmap for each
  right = NULL
  if(length(ssm_cols) > 0){
    bars = lapply(ssm_cols, function(col){
      ComplexHeatmap::anno_barplot(meta[[col]], which = "row", border = FALSE, bar_width = 1,
                                   gp = grid::gpar(fill = "grey30", col = NA),
                                   axis_param = list(gp = grid::gpar(fontsize = 6)),
                                   width = grid::unit(ssm_bar_width, "mm"))
    })
    names(bars) = sub("_SSMs$", " SSMs", ssm_cols)
    right = do.call(ComplexHeatmap::rowAnnotation,
                    c(bars, list(annotation_name_side = "top", annotation_name_rot = 90,
                                 annotation_name_gp = grid::gpar(fontsize = 8),
                                 gap = grid::unit(1.5, "mm"))))
  }
  col_fun = circlize::colorRamp2(c(0, 1, 2, 3, cn_max),
                                 c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B"))
  titles = paste0(wins$name, "\n", wins$region)
  if(is.null(show_row_names)) show_row_names = length(ids) <= 100
  ht = ComplexHeatmap::Heatmap(
    mat, name = if(!is.null(center) || "CN_uncentered" %in% colnames(seg)) "CN\n(centred)" else "CN",
    col = col_fun, na_col = "grey85",
    cluster_rows = FALSE, cluster_columns = FALSE,
    column_split = factor(bins$region, levels = wins$name), column_title = titles,
    column_title_gp = grid::gpar(fontsize = 9), column_gap = grid::unit(2, "mm"),
    row_split = if(is.null(split_by)) NULL else groups,
    row_title_gp = grid::gpar(fontsize = 9), row_gap = grid::unit(1, "mm"), row_title_rot = 0,
    show_column_names = FALSE, show_row_names = show_row_names,
    row_names_gp = grid::gpar(fontsize = 6),
    left_annotation = left, top_annotation = top, right_annotation = right, layer_fun = layer,
    use_raster = length(ids) > 500)
  if(!return_data) return(ht)
  # the input rows behind the plot
  shown = rownames(mat)
  segs = as.data.frame(seg_data)
  seg_chrom = gsub("^chr", "", as.character(segs$chrom))
  events = do.call(rbind, lapply(seq_len(nrow(wins)), function(w){
    hit = as.character(segs$ID) %in% shown & seg_chrom == wins$chrom[w] &
      segs$end > wins$start[w] & segs$start < wins$end[w]
    if(!any(hit)) return(NULL)
    data.frame(region = wins$name[w], segs[hit, , drop = FALSE], row.names = NULL, check.names = FALSE)
  }))
  sv_events = NULL
  if(!is.null(svs) && nrow(svs) > 0){
    shown_sv = unique(svs[, c("sv_id", "link", "pair", "orientation")])
    shown_sv = shown_sv[!duplicated(shown_sv$sv_id), , drop = FALSE]
    sv_events = data.frame(as.data.frame(sv_data)[shown_sv$sv_id, , drop = FALSE],
                           sv_class = shown_sv$link, sv_pair = shown_sv$pair,
                           sv_orientation = shown_sv$orientation, row.names = NULL, check.names = FALSE)
  }
  ssm_events = if(is.null(ssm_data)) NULL else .multi_region_ssm_rows(ssm_data, wins, shown)
  list(plot = ht, heatmap = ht, events = events, sv_events = sv_events, ssm_events = ssm_events,
       matrix = mat, sample_order = shown, metadata = meta,
       ssm_counts = if(length(ssm_cols) > 0) data.frame(sample_id = meta$sample_id, meta[, ssm_cols, drop = FALSE],
                                                        check.names = FALSE, row.names = NULL) else NULL,
       regions = wins[, c("name", "region", "chrom", "start", "end", "bins")], bins = bins)
}

# The regions as a table: name, region string, chrom, start, end and the
# genes in each (for marking), from region strings, genes or lists
.multi_region_windows = function(regions, flank, genome_build){
  if(!is.list(regions)) regions = as.list(regions)
  if(length(regions) < 1) stop("give at least one region")
  rows = lapply(seq_along(regions), function(i){
    r = regions[[i]]
    if(is.list(r)){
      w = .resolve_window(region = r$region, genes = r$genes,
                          flank = if(is.null(r$flank)) flank else r$flank,
                          genome_build = genome_build, quiet = TRUE)
      label = if(!is.null(r$genes)) paste(r$genes, collapse = "+") else r$region
    }else if(length(r) == 1 && grepl(":", r)){
      w = .resolve_window(region = r, genome_build = genome_build, quiet = TRUE)
      label = r
    }else{
      w = .resolve_window(genes = r, flank = flank, genome_build = genome_build, quiet = TRUE)
      label = paste(r, collapse = "+")
    }
    name = names(regions)[i]
    if(is.null(name) || is.na(name) || name == "") name = label
    list(name = name, region = w$region, chrom = w$chrom, start = w$start, end = w$end,
         genes = w$genes)
  })
  out = data.frame(name = vapply(rows, `[[`, "", "name"), region = vapply(rows, `[[`, "", "region"),
                   chrom = vapply(rows, `[[`, "", "chrom"),
                   start = vapply(rows, function(x) as.numeric(x$start), 0),
                   end = vapply(rows, function(x) as.numeric(x$end), 0), stringsAsFactors = FALSE)
  if(any(duplicated(out$name))) stop("region names must be unique")
  out$genes = lapply(rows, `[[`, "genes")
  out
}

# Row order from sort keys (metadata columns, region names, "cluster"),
# applied within groups; "-" reverses a key
.multi_region_order = function(mat, meta, bins, sort_by, groups){
  n = nrow(mat)
  if(is.null(sort_by) || length(sort_by) == 0) return(seq_len(n))
  keys = list()
  cluster = FALSE
  for(k in sort_by){
    reverse = startsWith(k, "-")
    name = sub("^-", "", k)
    if(name == "cluster"){
      cluster = TRUE
      next
    }
    if(name %in% bins$region){
      v = rowMeans(mat[, bins$region == name, drop = FALSE], na.rm = TRUE)
      v[is.nan(v)] = NA
      keys[[k]] = if(reverse) v else -v                   # highest first by default
    }else if(name %in% colnames(meta)){
      v = meta[[name]]
      v = if(is.numeric(v)) v else as.integer(factor(v, levels = .level_order(v)))
      keys[[k]] = if(reverse) -v else v
    }else{
      stop("sort_by key \"", name, "\" is neither a metadata column, a region name nor \"cluster\"",
           if(name %in% c("a_pos", "b_pos", "a_step", "b_step", "step_type"))
             paste0(" (", name, " comes from GAMBLR.utils::cn_step_concordance(): pass its result ",
                    "as concordance, or join its steps to the metadata)") else "", call. = FALSE)
    }
  }
  base = if(length(keys) > 0) do.call(order, c(keys, list(na.last = TRUE))) else seq_len(n)
  if(!cluster) return(base)
  # clustering breaks the remaining ties: within groups of equal keys
  tie = if(length(keys) > 0) do.call(paste, c(lapply(keys, function(v) v[base]), sep = "\r"))
        else rep("", n)
  tie = paste(groups[base], tie, sep = "\r")
  out = integer(0)
  for(t in unique(tie)){
    idx = base[tie == t]
    if(length(idx) > 2){
      m = mat[idx, , drop = FALSE]
      m[is.na(m)] = 2
      idx = idx[stats::hclust(stats::dist(m))$order]
    }
    out = c(out, idx)
  }
  out
}

# Levels in a sensible order: a factor's own, otherwise sorted (NA last)
.level_order = function(v){
  if(is.factor(v)) levels(v) else sort(unique(as.character(v[!is.na(v)])))
}


# SV breakpoints falling in the regions' bins: one row per breakpoint with
# its sample's row in the matrix, its bin and where in it, and whether the
# other end is in another region ("linked"), the same region ("within") or
# elsewhere ("other")
.multi_region_svs = function(sv_data, wins, bins, sample_ids){
  sv = as.data.frame(sv_data)
  needed = c("CHROM_A", "START_A", "CHROM_B", "START_B")
  if(!all(needed %in% colnames(sv))) stop("sv_data needs columns CHROM_A, START_A, CHROM_B and START_B")
  # only SVs of the samples shown with an end on a plotted chromosome
  ids = .sv_sample_ids(sv)
  keep = ids %in% sample_ids & (gsub("^chr", "", as.character(sv$CHROM_A)) %in% wins$chrom |
                                  gsub("^chr", "", as.character(sv$CHROM_B)) %in% wins$chrom)
  input_row = which(keep)
  sv = sv[keep, , drop = FALSE]
  ids = ids[keep]
  if(nrow(sv) == 0) return(NULL)
  ends = function(side, other){
    pos = function(sd){
      st = as.numeric(sv[[paste0("START_", sd)]])
      en = if(paste0("END_", sd) %in% colnames(sv)) as.numeric(sv[[paste0("END_", sd)]]) else st
      (st + en) / 2
    }
    strand = function(sd){
      col = paste0("STRAND_", sd)
      if(col %in% colnames(sv)) as.character(sv[[col]]) else rep(NA_character_, nrow(sv))
    }
    data.frame(ID = ids, sv_id = input_row,
               chrom = gsub("^chr", "", as.character(sv[[paste0("CHROM_", side)]])),
               pos = pos(side), strand = strand(side),
               other_chrom = gsub("^chr", "", as.character(sv[[paste0("CHROM_", other)]])),
               other_pos = pos(other), other_strand = strand(other), stringsAsFactors = FALSE)
  }
  bp = rbind(ends("A", "B"), ends("B", "A"))
  bp = bp[bp$ID %in% sample_ids, , drop = FALSE]
  region_of = function(chrom, pos){
    r = rep(NA_character_, length(chrom))
    for(w in seq_len(nrow(wins))){
      hit = chrom == wins$chrom[w] & pos >= wins$start[w] & pos < wins$end[w]
      r[hit & is.na(r)] = wins$name[w]
    }
    r
  }
  bp$region = region_of(bp$chrom, bp$pos)
  bp = bp[!is.na(bp$region), , drop = FALSE]
  if(nrow(bp) == 0) return(NULL)
  other_region = region_of(bp$other_chrom, bp$other_pos)
  bp$link = ifelse(is.na(other_region), "other", ifelse(other_region == bp$region, "within", "linked"))
  # linked junctions: the region pair (in the order of the regions) and the
  # strands of the ends in that order (BEDPE: "+" keeps the side before the
  # breakpoint), so both points of a junction share its orientation; the
  # reciprocal junctions of a balanced translocation have complementary
  # orientations (e.g. "+-" and "-+")
  order_of = match(bp$region, wins$name); other_order = match(other_region, wins$name)
  first_is_this = is.na(other_order) | order_of <= other_order
  bp$pair = ifelse(bp$link == "linked",
                   ifelse(first_is_this, paste(bp$region, other_region, sep = "-"),
                          paste(other_region, bp$region, sep = "-")), NA_character_)
  bp$orientation = ifelse(bp$link == "linked",
                          ifelse(first_is_this, paste0(bp$strand, bp$other_strand),
                                 paste0(bp$other_strand, bp$strand)), NA_character_)
  bp$orientation[grepl("NA", bp$orientation)] = NA_character_
  shapes = c("+-" = 24, "-+" = 25, "++" = 23, "--" = 22)
  bp$pch = ifelse(bp$link == "linked" & !is.na(bp$orientation), shapes[bp$orientation], 21)
  bp$col = vapply(seq_len(nrow(bp)), function(k){
    m = which(bins$region == bp$region[k] & bins$start <= bp$pos[k] & bins$end > bp$pos[k])
    if(length(m) == 0) NA_integer_ else m[1]
  }, integer(1))
  bp = bp[!is.na(bp$col), , drop = FALSE]
  bp$frac = (bp$pos - bins$start[bp$col]) / (bins$end[bp$col] - bins$start[bp$col])
  bp$row = match(bp$ID, sample_ids)
  bp
}

# Mutations per sample in each region (samples by regions), from one MAF or
# a list of them; a mutation repeated (e.g. in two region queries) counts once
.multi_region_ssm_counts = function(ssm_data, wins, sample_ids){
  mafs = if(is.data.frame(ssm_data)) list(ssm_data) else ssm_data
  cols = c("Tumor_Sample_Barcode", "Chromosome", "Start_Position")
  maf = do.call(rbind, lapply(mafs, function(m){
    m = as.data.frame(m)
    if(!all(cols %in% colnames(m))) stop("ssm_data needs columns ", paste(cols, collapse = ", "))
    alleles = intersect(c("Reference_Allele", "Tumor_Seq_Allele2"), colnames(m))
    m[, c(cols, alleles), drop = FALSE]
  }))
  maf = maf[!duplicated(maf), , drop = FALSE]
  maf = maf[as.character(maf$Tumor_Sample_Barcode) %in% sample_ids, , drop = FALSE]
  chrom = gsub("^chr", "", as.character(maf$Chromosome))
  pos = as.numeric(maf$Start_Position)
  counts = matrix(0L, nrow = length(sample_ids), ncol = nrow(wins),
                  dimnames = list(sample_ids, wins$name))
  for(w in seq_len(nrow(wins))){
    hit = chrom == wins$chrom[w] & pos >= wins$start[w] & pos < wins$end[w]
    if(!any(hit)) next
    n = table(factor(as.character(maf$Tumor_Sample_Barcode[hit]), levels = sample_ids))
    counts[, w] = as.integer(n)
  }
  counts
}

# The mutations counted by .multi_region_ssm_counts: MAF rows of the samples
# shown in each region (duplicates once), with the region
.multi_region_ssm_rows = function(ssm_data, wins, sample_ids){
  mafs = if(is.data.frame(ssm_data)) list(ssm_data) else ssm_data
  maf = do.call(rbind, lapply(mafs, function(m){
    m = as.data.frame(m)
    m[, intersect(colnames(m), colnames(as.data.frame(mafs[[1]]))), drop = FALSE]
  }))
  key_cols = intersect(c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", "Reference_Allele",
                         "Tumor_Seq_Allele2"), colnames(maf))
  maf = maf[!duplicated(maf[, key_cols]), , drop = FALSE]
  maf = maf[as.character(maf$Tumor_Sample_Barcode) %in% sample_ids, , drop = FALSE]
  chrom = gsub("^chr", "", as.character(maf$Chromosome))
  pos = as.numeric(maf$Start_Position)
  do.call(rbind, lapply(seq_len(nrow(wins)), function(w){
    hit = chrom == wins$chrom[w] & pos >= wins$start[w] & pos < wins$end[w]
    if(!any(hit)) return(NULL)
    data.frame(region = wins$name[w], maf[hit, , drop = FALSE], row.names = NULL, check.names = FALSE)
  }))
}
