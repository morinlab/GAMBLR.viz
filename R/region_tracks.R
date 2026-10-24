# Shared building blocks for region plots: plot_sv_density(),
# plot_ssm_density() and plot_cohort_CN_region() use these so their windows,
# gene marks, groups, bins and tracks match and can be stacked.

# The window (chrom, start, end) from region, chrom/start/end or genes +
# flank, and the gene marks in it: `genes` (index = TRUE) and `label_genes`.
.resolve_window = function(region = NULL, chrom = NULL, start = NULL, end = NULL,
                           genes = NULL, flank = 5e6, label_genes = NULL,
                           genome_build = "grch37", quiet = FALSE){
  locate = function(x){
    if(is.null(x) || length(x) == 0) return(NULL)
    if(is.character(x)){
      g = as.data.frame(suppressMessages(GAMBLR.utils::gene_to_region(
        gene_symbol = x, projection = genome_build, return_as = "df")))
      if(nrow(g) == 0) return(NULL)
      colnames(g)[seq_len(3)] = c("chrom", "start", "end")
      gene_col = grep("hugo|gene", colnames(g), ignore.case = TRUE)[1]
      g = data.frame(gene = g[[gene_col]], chrom = g$chrom, start = g$start, end = g$end,
                     strand = NA_character_)
    }else{
      x = as.data.frame(x)
      g = x[, c("gene", "chrom", "start", "end")]
      g$strand = if("strand" %in% colnames(x)) as.character(x$strand) else NA_character_
    }
    g$chrom = gsub("^chr", "", as.character(g$chrom))
    g
  }
  gene_df = locate(genes)
  derived = FALSE
  if(!is.null(region)){
    clean = gsub(",", "", region)
    parts = regmatches(clean, regexec("^(chr)?([0-9XY]+):([0-9.e+]+)-([0-9.e+]+)$", clean))[[1]]
    if(length(parts) == 0){
      stop("region must look like \"9:20000000-23000000\"")
    }
    chrom = parts[3]; start = as.numeric(parts[4]); end = as.numeric(parts[5])
  }else if(is.null(chrom) && !is.null(gene_df)){
    if(length(unique(gene_df$chrom)) != 1){
      stop("genes are on different chromosomes; give region")
    }
    chrom = gene_df$chrom[1]
    start = max(0, min(gene_df$start) - flank)
    end = max(gene_df$end) + flank
    derived = TRUE
  }
  if(is.null(chrom) || is.null(start) || is.null(end)){
    stop("give region, chrom/start/end, or genes")
  }
  chrom = gsub("^chr", "", as.character(chrom))
  in_window = function(g) if(is.null(g)) NULL else
    g[g$chrom == chrom & g$end > start & g$start < end, , drop = FALSE]
  gene_df = in_window(gene_df)
  extra = in_window(locate(label_genes))
  if(!is.null(extra)){
    extra = extra[!extra$gene %in% gene_df$gene & !duplicated(extra$gene), , drop = FALSE]
  }
  marks = rbind(if(!is.null(gene_df) && nrow(gene_df) > 0)
                  data.frame(gene_df[, c("gene", "start", "end", "strand")], index = TRUE),
                if(!is.null(extra) && nrow(extra) > 0)
                  data.frame(extra[, c("gene", "start", "end", "strand")], index = FALSE))
  region = .region_string(chrom, start, end)
  if(derived && !quiet){
    message("region: ", region, " (", paste(unique(gene_df$gene), collapse = ", "), " +/- ",
            format(flank, big.mark = ",", scientific = FALSE), " bp)")
  }
  list(chrom = chrom, start = start, end = end, genes = gene_df, marks = marks, region = region)
}

# "chrN:start-end" for a window (GAMBLR getters accept it for either build)
.region_string = function(chrom, start, end){
  paste0("chr", gsub("^chr", "", chrom), ":", format(floor(start), scientific = FALSE, trim = TRUE),
         "-", format(ceiling(end), scientific = FALSE, trim = TRUE))
}

#' @title The region a region plot shows
#'
#' @description Get the window (as "chrN:start-end") that the region plots
#' ([plot_cohort_CN_region], [plot_sv_density], [plot_ssm_density],
#' [plot_region_genes]) show for the same `genes` and `flank` (or `region`,
#' or `chrom`/`start`/`end`), e.g. to get the data for that window first:
#' the plots compute it the same way.
#'
#' @param genes Gene symbols (or a data frame with gene, chrom, start, end).
#' The window spans them, +/- `flank`.
#' @param flank Bases added on each side of the genes. Default 5e6, as in
#' [plot_cohort_CN_region]; [plot_sv_density] defaults to 1e6 and
#' [plot_ssm_density] to 5e4, so give the same value as in the plot call.
#' @param region,chrom,start,end Alternatively the window itself.
#' @param genome_build "grch37" (default) or "hg38".
#'
#' @return A region string, "chrN:start-end".
#'
#' @export
#'
#' @examples
#' \dontrun{
#' region = region_window("MALT1", flank = 5e6)
#' ssm = GAMBLR.results::get_ssm_by_region(region = region, these_samples_metadata = meta)
#' plot_cohort_CN_region(seg, genes = "MALT1", flank = 5e6, ssm_data = ssm,
#'                       these_samples_metadata = meta)
#' }
region_window = function(genes = NULL, flank = 5e6, region = NULL, chrom = NULL, start = NULL,
                         end = NULL, genome_build = "grch37"){
  .resolve_window(region, chrom, start, end, genes, flank, NULL, genome_build, quiet = TRUE)$region
}

# Add gene models (as from GAMBLR.utils::get_gene_models: gene, chrom,
# start, end, strand, exon_start, exon_end) to a resolved window: fills in
# strand for marked genes, adds the models' other genes in the window as
# lighter marks, and keeps their exons (w$exons) for the gene track.
# The gene models to use for a window: given ones (a data frame), none
# (FALSE), or by default fetched for the window from mygene.info via
# GAMBLR.utils::get_gene_models() (NULL if that fails, e.g. offline)
.resolve_gene_models = function(w, gene_models, genome_build, needed = TRUE,
                                genes_show_isoforms = NULL){
  if(isFALSE(gene_models) || !needed) return(NULL)
  if(!is.null(gene_models)) return(as.data.frame(gene_models))
  region = paste0(w$chrom, ":", format(floor(w$start), scientific = FALSE), "-",
                  format(ceiling(w$end), scientific = FALSE))
  fail = function(e){
    message("gene models not retrieved (", conditionMessage(e),
            "); genes are drawn without exons. Use gene_models = FALSE to skip ",
            "this, or give your own")
    NULL
  }
  models = tryCatch(suppressMessages(GAMBLR.utils::get_gene_models(region = region,
                                                                   genome_build = genome_build)),
                    error = fail)
  if(is.null(models) || length(genes_show_isoforms) == 0) return(models)
  # every transcript of the genes asked for, in place of their merged models
  iso = tryCatch(suppressMessages(GAMBLR.utils::get_gene_models(genes = genes_show_isoforms,
                                                                genome_build = genome_build,
                                                                transcripts = "all")),
                 error = function(e) NULL)
  if(is.null(iso) || nrow(iso) == 0) return(models)
  rbind(models[!models$gene %in% iso$gene, , drop = FALSE], iso[, colnames(models)])
}

# One transcript per distinct pair of ends: transcripts whose
# transcription start (strand-aware) and 3' end are both within
# end_tolerance of a kept transcript's are dropped (e.g. skipped-exon
# versions); kept ones are those with the most exons, then the longest, and
# their label notes how many they stand for.
.distinct_end_isoforms = function(iso, end_tolerance){
  out = lapply(split(iso, iso$parent), function(x){
    tx = x[!duplicated(x$gene), c("gene", "start", "end", "strand"), drop = FALSE]
    minus = !is.na(tx$strand) & tx$strand == "-"
    tx$tss = ifelse(minus, tx$end, tx$start)
    tx$tes = ifelse(minus, tx$start, tx$end)
    tx$n_exons = as.numeric(table(x$gene)[tx$gene])
    tx = tx[order(-tx$n_exons, -(tx$end - tx$start)), , drop = FALSE]
    kept = integer(0); covers = integer(0)
    for(i in seq_len(nrow(tx))){
      same = kept[abs(tx$tss[kept] - tx$tss[i]) <= end_tolerance &
                    abs(tx$tes[kept] - tx$tes[i]) <= end_tolerance]
      if(length(same) == 0){
        kept = c(kept, i); covers = c(covers, 0)
      }else{
        covers[match(same[1], kept)] = covers[match(same[1], kept)] + 1
      }
    }
    label = stats::setNames(ifelse(covers > 0, paste0(tx$gene[kept], " (+", covers, ")"),
                                   tx$gene[kept]), tx$gene[kept])
    x = x[x$gene %in% names(label), , drop = FALSE]
    x$gene = unname(label[x$gene])
    x
  })
  do.call(rbind, out)
}

# For the genes in genes_show_isoforms, one entry per transcript (gene label
# "GENE transcript", extent from its exons) in place of the gene, in the
# models and the marks; an isoform keeps the gene's index status, and only
# the first is shaded in the plot.
.expand_isoforms = function(m, marks, genes_show_isoforms, isoforms = "distinct_ends",
                            end_tolerance = 100){
  if(length(genes_show_isoforms) == 0 || !"transcript" %in% colnames(m)) return(list(m = m, marks = marks))
  sel = m$gene %in% genes_show_isoforms & !is.na(m$transcript) & m$transcript != "union"
  missing = setdiff(intersect(genes_show_isoforms, m$gene), m$gene[sel])
  if(length(missing) > 0){
    message("no per-transcript models for ", paste(missing, collapse = ", "),
            " (genes_show_isoforms needs gene_models with transcripts = \"all\")")
  }
  if(!any(sel)) return(list(m = m, marks = marks))
  parents = unique(m$gene[sel])
  iso = m[sel, , drop = FALSE]
  iso$parent = iso$gene
  iso$gene = paste(iso$gene, iso$transcript)
  iso$start = as.numeric(tapply(iso$exon_start, iso$gene, min)[iso$gene])
  iso$end = as.numeric(tapply(iso$exon_end, iso$gene, max)[iso$gene])
  if(isoforms == "distinct_ends"){
    iso = .distinct_end_isoforms(iso, end_tolerance)
  }
  m = rbind(m[!m$gene %in% parents, , drop = FALSE], iso[, colnames(m)])
  for(parent in parents){
    tx = iso[iso$parent == parent & !duplicated(iso$gene), , drop = FALSE]
    tx = tx[order(tx$start), , drop = FALSE]
    at = match(parent, marks$gene)
    index = if(is.na(at)) FALSE else marks$index[at]
    shade = if(is.na(at)) FALSE else marks$shade[at]
    if(is.null(marks$track)) marks$track = rep(TRUE, nrow(marks))
    # one row per isoform in the gene track, and the gene's whole extent
    # (all isoforms shown) for the shading only
    entries = data.frame(gene = tx$gene, start = tx$start, end = tx$end, strand = tx$strand,
                         index = index, shade = FALSE, track = TRUE)
    span = data.frame(gene = parent, start = min(tx$start), end = max(tx$end),
                      strand = tx$strand[1], index = index, shade = shade, track = FALSE)
    marks = rbind(if(is.na(at)) marks else marks[-at, , drop = FALSE], entries, span)
  }
  list(m = m, marks = marks)
}

.add_gene_models = function(w, gene_models, max_extra_genes = Inf, genes_show_isoforms = NULL,
                            isoforms = "distinct_ends"){
  if(is.null(gene_models) || nrow(gene_models) == 0) return(w)
  m = as.data.frame(gene_models)
  needed = c("gene", "chrom", "start", "end")
  if(!all(needed %in% colnames(m))){
    stop("gene_models needs columns gene, chrom, start and end (see GAMBLR.utils::get_gene_models)")
  }
  m$chrom = gsub("^chr", "", as.character(m$chrom))
  m = m[m$chrom == w$chrom & m$end > w$start & m$start < w$end, , drop = FALSE]
  if(nrow(m) == 0) return(w)
  if(is.null(m$strand)) m$strand = NA_character_
  genes = m[!duplicated(m$gene), c("gene", "start", "end", "strand"), drop = FALSE]
  marks = w$marks
  if(is.null(marks)) marks = data.frame(gene = character(0), start = numeric(0),
                                        end = numeric(0), strand = character(0),
                                        index = logical(0))
  if(is.null(marks$shade)) marks$shade = rep(TRUE, nrow(marks))
  expanded = .expand_isoforms(m, marks, genes_show_isoforms, isoforms)
  m = expanded$m; marks = expanded$marks
  genes = m[!duplicated(m$gene), c("gene", "start", "end", "strand"), drop = FALSE]
  known = match(marks$gene, genes$gene)
  fill = is.na(marks$strand) & !is.na(known)
  marks$strand[fill] = genes$strand[known[fill]]
  # genes with models take their extent from them, so the drawn model, the
  # shading and the coordinates agree
  has = !is.na(known)
  marks$start[has] = genes$start[known[has]]
  marks$end[has] = genes$end[known[has]]
  new = genes[!genes$gene %in% marks$gene, , drop = FALSE]
  if(nrow(new) > max_extra_genes){
    if(max_extra_genes > 0){
      message(nrow(new), " other genes in the window are not drawn (more than ",
              max_extra_genes, "); name them in genes/label_genes, or use other_genes = \"all\"")
    }
    new = new[0, , drop = FALSE]
    m = m[m$gene %in% marks$gene, , drop = FALSE]
  }
  # genes only from the models: in the gene track, not shaded in the panels
  if(nrow(new) > 0){
    extra = data.frame(new, index = FALSE, shade = FALSE)
    if(!is.null(marks$track)) extra$track = TRUE
    marks = rbind(marks, extra)
  }
  w$marks = marks
  if(all(c("exon_start", "exon_end") %in% colnames(m))){
    if(is.null(m$feature)) m$feature = "coding"
    w$exons = m[!is.na(m$exon_start), c("gene", "exon_start", "exon_end", "feature"), drop = FALSE]
  }
  w
}

# Bin edges for a window: bin_size (bp) if given, otherwise bin_fraction of
# the window; refuses more than max_bins bins (costly and unreadable).
.region_bins = function(start, end, bin_size = NULL, bin_fraction = 0.01,
                        max_bins = 2000){
  width = end - start
  if(!is.null(bin_size)){
    if(!is.numeric(bin_size) || length(bin_size) != 1 || bin_size <= 0){
      stop("bin_size must be a positive number of bp")
    }
    bin = bin_size
  }else{
    if(!is.numeric(bin_fraction) || length(bin_fraction) != 1 ||
       bin_fraction <= 0 || bin_fraction > 1){
      stop("bin_fraction must be between 0 and 1")
    }
    bin = width * bin_fraction
  }
  n = ceiling(width / bin)
  if(n > max_bins){
    stop(sprintf(paste0("%s bins of %s bp across this %s Mb window is more than ",
                        "max_bins (%s); use a larger bin_size or bin_fraction ",
                        "(at least %s bp here)"),
                 format(n, big.mark = ","), format(round(bin), big.mark = ","),
                 round(width / 1e6, 2), max_bins,
                 format(ceiling(width / max_bins), big.mark = ",")))
  }
  unique(c(seq(start, end, by = bin), end))
}

# Sample -> group (and counting unit) for the samples with data. Units are
# samples, or with measure "patients", patients within a group (a patient's
# samples in different groups count once in each).
.region_units = function(ids, these_samples_metadata = NULL, split_by = NULL,
                         split_order = NULL, measure = "patients"){
  ids = unique(as.character(ids))
  meta = if(is.null(these_samples_metadata)) NULL else as.data.frame(these_samples_metadata)
  if(!is.null(meta)){
    ids = intersect(ids, as.character(meta$sample_id))
  }
  group = rep("all samples", length(ids))
  if(!is.null(split_by)){
    if(is.null(meta) || !split_by %in% colnames(meta)){
      stop("split_by must name a column of these_samples_metadata")
    }
    group = as.character(meta[[split_by]][match(ids, meta$sample_id)])
  }
  unit = ids
  counted = if(measure == "patients") "patients" else "samples"
  if(measure == "patients"){
    if(is.null(meta) || !"patient_id" %in% colnames(meta)){
      message("no patient_id in these_samples_metadata; counting samples")
      counted = "samples"
    }else{
      patient = as.character(meta$patient_id[match(ids, meta$sample_id)])
      unit = ifelse(is.na(patient), ids, patient)
    }
  }
  map = data.frame(ID = ids, group = group, unit = paste(group, unit, sep = "\r"),
                   stringsAsFactors = FALSE)
  map = map[!is.na(map$group), , drop = FALSE]
  groups = unique(map$group)
  if(!is.null(split_order)){
    groups = c(intersect(split_order, groups), setdiff(groups, split_order))
  }else if(!is.null(split_by) && is.factor(meta[[split_by]])){
    groups = intersect(levels(meta[[split_by]]), groups)
  }else{
    groups = sort(groups)
  }
  list(map = map, groups = groups, counted = counted)
}

# Per bin, group and class: with measure "patients"/"samples" the number and
# fraction of units with at least one event (a unit counts once per bin, in
# its first class by level order); with "mutations", events per sample per
# kb. points: ID, pos and optional class (factor).
.binned_density = function(points, units, breaks, measure = "patients"){
  map = units$map
  groups = units$groups
  points = points[points$ID %in% map$ID, , drop = FALSE]
  if(is.null(points$class)) points$class = factor(rep("all", nrow(points)), levels = "all")
  classes = levels(points$class)
  points$bin = findInterval(points$pos, breaks, rightmost.closed = TRUE, all.inside = TRUE)
  points$group = map$group[match(points$ID, map$ID)]
  points$unit = map$unit[match(points$ID, map$ID)]
  n_bins = length(breaks) - 1
  out = do.call(rbind, lapply(groups, function(g){
    p = points[points$group == g, , drop = FALSE]
    n_units = length(unique(map$unit[map$group == g]))
    n_samples = sum(map$group == g)
    counts = matrix(0, n_bins, length(classes))
    if(nrow(p) > 0){
      if(measure == "mutations"){
        counts = unclass(table(factor(p$bin, levels = seq_len(n_bins)), p$class))
      }else{
        # one class per unit and bin: the first by priority (default: level
        # order)
        prio = if(!is.null(p$priority)) p$priority else as.integer(p$class)
        p = p[order(p$unit, p$bin, prio), , drop = FALSE]
        p = p[!duplicated(paste(p$unit, p$bin)), , drop = FALSE]
        counts = unclass(table(factor(p$bin, levels = seq_len(n_bins)), p$class))
      }
    }
    denom = if(measure == "mutations") n_samples * diff(breaks) / 1e3 else rep(n_units, n_bins)
    data.frame(group = g, start = rep(breaks[-length(breaks)], length(classes)),
               end = rep(breaks[-1], length(classes)),
               class = factor(rep(classes, each = n_bins), levels = classes),
               n = as.vector(counts), n_units = if(measure == "mutations") n_samples else n_units,
               value = ifelse(rep(denom, length(classes)) > 0,
                              as.vector(counts) / rep(denom, length(classes)), 0),
               stringsAsFactors = FALSE)
  }))
  out
}

# Genes in a data panel: a faint span from start to end (a dashed line for
# genes too narrow to see), optionally labelled. Labels flip left when less
# of them would then sit on the `avoid` boxes (start, end in bp). With a gene
# track, panels are drawn without labels.
.gene_layers = function(marks, start, end, base_size, with_labels = FALSE,
                        label_facet = NULL, facet_levels = NULL, avoid = NULL){
  if(!is.null(marks) && !is.null(marks$shade)) marks = marks[marks$shade, , drop = FALSE]
  if(is.null(marks) || nrow(marks) == 0) return(list())
  lo = pmax(marks$start, start) / 1e6
  hi = pmin(marks$end, end) / 1e6
  mid = (marks$start + marks$end) / 2e6
  narrow = (hi - lo) < 0.003 * (end - start) / 1e6
  out = list()
  if(any(!narrow)){
    span = data.frame(lo = lo, hi = hi, alpha = ifelse(marks$index, 0.12, 0.06))[!narrow, ]
    out = c(out, list(geom_rect(data = span, aes(xmin = lo, xmax = hi, ymin = -Inf,
                                                 ymax = Inf, alpha = I(alpha)),
                                fill = "grey20", inherit.aes = FALSE)))
  }
  if(any(narrow)){
    line = data.frame(x = mid, colour = ifelse(marks$index, "grey30", "grey60"))[narrow, ]
    out = c(out, list(geom_vline(data = line, aes(xintercept = x, colour = I(colour)),
                                 linetype = "22", linewidth = 0.3, inherit.aes = FALSE)))
  }
  if(!with_labels) return(out)
  o = order(mid)
  vjust = numeric(length(mid))
  vjust[o] = c(1.2, 2.6, 4)[(seq_along(o) - 1) %% 3 + 1]
  text_colour = ifelse(marks$index, "black", "grey40")
  hjust = rep(-0.05, length(mid))
  if(!is.null(avoid) && nrow(avoid) > 0){
    width = nchar(marks$gene) * 0.012 * (end - start) / 1e6
    covered = function(lo, hi) vapply(seq_along(lo), function(i)
      sum(pmax(0, pmin(avoid$end / 1e6, hi[i]) - pmax(avoid$start / 1e6, lo[i]))),
      numeric(1))
    hjust[covered(mid - width, mid) < covered(mid, mid + width)] = 1.05
  }
  if(is.null(label_facet)){
    c(out, list(annotate("text", x = mid, y = Inf, label = marks$gene,
                         vjust = vjust, hjust = hjust, size = base_size / 3.2,
                         fontface = "italic", colour = text_colour)))
  }else{
    lab = data.frame(x = mid, gene = marks$gene, vjust = vjust, hjust = hjust,
                     colour = text_colour,
                     group = factor(label_facet, levels = facet_levels))
    c(out, list(geom_text(data = lab, aes(x = x, y = Inf, label = gene,
                                          vjust = vjust, hjust = hjust,
                                          colour = I(colour)),
                          size = base_size / 3.2, fontface = "italic",
                          inherit.aes = FALSE)))
  }
}

# A gene track: each gene as a bar over its extent (at least a sliver wide),
# with strand chevrons when strand is known and its name above; genes whose
# bars or labels would overlap go on separate rows. Returns list(plot, rows).
.gene_track = function(marks, start, end, base_size = 9, exons = NULL){
  if(!is.null(marks) && !is.null(marks$track)) marks = marks[marks$track, , drop = FALSE]
  if(is.null(marks) || nrow(marks) == 0) return(NULL)
  w = end - start
  g = marks
  if(is.null(g$strand)) g$strand = NA_character_
  g$lo = pmax(g$start, start); g$hi = pmin(g$end, end)
  min_w = 0.004 * w
  short = g$hi - g$lo < min_w
  centre = (g$lo + g$hi) / 2
  g$lo[short] = centre[short] - min_w / 2; g$hi[short] = centre[short] + min_w / 2
  g$mid = (g$lo + g$hi) / 2
  label_w = nchar(g$gene) * 0.012 * w
  ext_lo = pmin(g$lo, g$mid - label_w / 2); ext_hi = pmax(g$hi, g$mid + label_w / 2)
  # greedy row packing on the extent of bar and label
  o = order(ext_lo)
  g = g[o, , drop = FALSE]; ext_lo = ext_lo[o]; ext_hi = ext_hi[o]
  row_end = numeric(0); g$row = NA_integer_
  for(i in seq_len(nrow(g))){
    r = which(row_end < ext_lo[i])[1]
    if(is.na(r)) r = length(row_end) + 1
    g$row[i] = r
    row_end[r] = ext_hi[i] + 0.005 * w
  }
  n = max(g$row)
  g$y = n - g$row + 1
  g$fill = ifelse(g$index, "grey25", "grey60")
  g$text = ifelse(g$index, "black", "grey35")
  # genes with an exon structure: intron line and exon blocks; others: a bar
  ex = NULL
  if(!is.null(exons) && nrow(exons) > 0){
    ex = exons[exons$gene %in% g$gene & exons$exon_end > start & exons$exon_start < end, , drop = FALSE]
    ex$lo = pmax(ex$exon_start, start); ex$hi = pmin(ex$exon_end, end)
    # keep exons visible at this scale
    tiny = ex$hi - ex$lo < 0.0015 * w
    centre = (ex$lo + ex$hi) / 2
    ex$lo[tiny] = centre[tiny] - 0.00075 * w; ex$hi[tiny] = centre[tiny] + 0.00075 * w
    ex$y = g$y[match(ex$gene, g$gene)]
    ex$fill = g$fill[match(ex$gene, g$gene)]
    # coding parts full height; UTRs and non-coding exons thinner
    if(is.null(ex$feature)) ex$feature = "coding"
    ex$half = ifelse(is.na(ex$feature) | ex$feature == "coding", 0.16, 0.09)
  }
  structured = g$gene %in% ex$gene
  step = w / 60
  chev = do.call(rbind, lapply(seq_len(nrow(g)), function(i){
    st = g$strand[i]
    if(is.na(st) || !st %in% c("+", "-", "1", "-1")) return(NULL)
    len = g$hi[i] - g$lo[i]
    x = if(len > 2 * step) seq(g$lo[i] + step / 2, g$hi[i] - step / 2, by = step) else g$mid[i]
    if(structured[i]){
      # on the intron line only, away from exon blocks
      e = ex[ex$gene == g$gene[i], , drop = FALSE]
      inside = vapply(x, function(v) any(v > e$lo - step / 3 & v < e$hi + step / 3), logical(1))
      x = x[!inside]
      if(length(x) == 0) return(NULL)
    }
    data.frame(x = x, y = g$y[i], label = if(st %in% c("+", "1")) ">" else "<",
               colour = if(structured[i]) g$fill[i] else "white")
  }))
  bars = g[!structured, , drop = FALSE]
  lines = g[structured, , drop = FALSE]
  p = ggplot(g)
  if(nrow(bars) > 0){
    p = p + geom_rect(data = bars, aes(xmin = lo / 1e6, xmax = hi / 1e6, ymin = y - 0.14,
                                       ymax = y + 0.14, fill = I(fill)))
  }
  if(nrow(lines) > 0){
    p = p + geom_segment(data = lines, aes(x = lo / 1e6, xend = hi / 1e6, y = y, yend = y,
                                           colour = I(fill)), linewidth = 0.4) +
      geom_rect(data = ex, aes(xmin = lo / 1e6, xmax = hi / 1e6, ymin = y - half,
                               ymax = y + half, fill = I(fill)), inherit.aes = FALSE)
  }
  p = p + geom_text(aes(x = mid / 1e6, y = y + 0.22, label = gene, colour = I(text)),
                    vjust = 0, size = base_size / 3.2, fontface = "italic")
  if(!is.null(chev) && nrow(chev) > 0){
    p = p + geom_text(data = chev, aes(x = x / 1e6, y = y, label = label, colour = I(colour)),
                      size = base_size / 3.8, fontface = "bold", vjust = 0.45,
                      inherit.aes = FALSE)
  }
  p = p + scale_x_continuous(limits = c(start, end) / 1e6, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.75, n + 0.8), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = base_size) +
    theme(plot.margin = margin(2, 5.5, 0, 5.5))
  list(plot = p, rows = n)
}

# Relative height for a gene track of `rows` rows (the data panels are ~1)
.gene_track_height = function(rows) 0.06 + 0.09 * rows

# A density track: two groups mirrored, otherwise one panel per group,
# stacked by class. density from .binned_density(); group_label names the
# groups (facet strips, or text inside a mirrored track).
.density_track = function(density, groups, group_label, start, end, marks,
                          y_label, measure, fill = "grey35", class_colours = NULL,
                          group_colours = NULL,
                          mirror = TRUE, base_size = 9, gene_labels = FALSE,
                          strips = TRUE, y_scale = "shared"){
  two_mirrored = mirror && length(groups) == 2
  free = y_scale == "free"
  # only classes that appear (a class can lose every bin to an earlier one,
  # since a sample counts once per bin; an empty legend key would be drawn)
  shown_classes = levels(density$class)[levels(density$class) %in%
                                          as.character(density$class[density$n > 0])]
  if(length(shown_classes) == 0) shown_classes = levels(density$class)[1]
  density$class = factor(as.character(density$class), levels = shown_classes)
  density = density[!is.na(density$class), , drop = FALSE]
  d = density[density$n > 0, , drop = FALSE]
  d = d[order(d$group, d$start, d$class), , drop = FALSE]
  if(nrow(d) > 0){
    d$ymax = stats::ave(d$value, d$group, d$start, FUN = cumsum)
  }else{
    d$ymax = numeric(0)
  }
  d$ymin = d$ymax - d$value
  label_fun = if(measure == "mutations") function(v) format(signif(abs(v), 2))
              else function(v) paste0(round(abs(v) * 100), "%")
  y_max = max(c(d$ymax, if(measure == "mutations") 1e-3 else 0.01))
  stacked = length(levels(density$class)) > 1
  # one colour, or (not stacked by class) each group in its own colour
  d$bar_fill = if(!stacked && !is.null(group_colours)){
    col = unname(group_colours[d$group]); col[is.na(col)] = fill; col
  }else rep(fill, nrow(d))
  rect = if(stacked) geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = ymin,
                                   ymax = ymax, fill = class))
         else geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = ymin, ymax = ymax,
                            fill = I(bar_fill)))
  floor_max = if(measure == "mutations") 1e-3 else 0.01
  if(two_mirrored){
    flip = d$group == groups[2]
    top_max = y_max; bottom_max = y_max
    if(free){
      # each half scaled to its own maximum (drawn to the same height),
      # each labelled with its own value
      top_max = max(c(d$ymax[!flip], floor_max)); bottom_max = max(c(d$ymax[flip], floor_max))
      d$ymin[flip] = d$ymin[flip] * y_max / bottom_max
      d$ymax[flip] = d$ymax[flip] * y_max / bottom_max
      d$ymin[!flip] = d$ymin[!flip] * y_max / top_max
      d$ymax[!flip] = d$ymax[!flip] * y_max / top_max
    }
    d$ymin[flip] = -d$ymin[flip]; d$ymax[flip] = -d$ymax[flip]
    p = ggplot(d) + rect +
      geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
      scale_y_continuous(limits = c(-y_max, y_max), breaks = c(-y_max, 0, y_max),
                         labels = c(label_fun(bottom_max), label_fun(0), label_fun(top_max)))
    if(strips){
      p = p + annotate("text", x = start / 1e6, y = c(y_max, -y_max),
                       label = group_label[groups], hjust = -0.02,
                       vjust = c(1.1, -0.1), size = base_size / 3.2)
    }
  }else{
    d$group = factor(group_label[d$group], levels = group_label[groups])
    p = ggplot(d) + rect
    if(free){
      # each panel scaled to its own maximum, which is its axis label
      p = p + geom_blank(data = data.frame(group = factor(group_label[groups],
                                                          levels = group_label[groups]),
                                           ymax = floor_max, start = start, end = start),
                         aes(x = start / 1e6, y = ymax), inherit.aes = FALSE) +
        facet_grid(group ~ ., drop = FALSE, switch = "y", scales = "free_y") +
        # the limits seen here include the 10% headroom above the maximum
        scale_y_continuous(breaks = function(l) max(l) / 1.1, labels = label_fun,
                           expand = expansion(mult = c(0, 0.1)))
    }else{
      p = p + facet_grid(group ~ ., drop = FALSE, switch = "y") +
        scale_y_continuous(limits = c(0, y_max),
                           breaks = if(length(groups) > 1) y_max else c(0, y_max),
                           labels = label_fun, expand = expansion(mult = c(0, 0.1)))
    }
  }
  if(stacked){
    values = class_colours
    if(is.null(values)){
      values = stats::setNames(grDevices::hcl.colors(length(levels(density$class)), "Dark 3"),
                               levels(density$class))
    }
    values = values[names(values) %in% levels(density$class)]
    p = p + scale_fill_manual(values = values, drop = FALSE, name = NULL)
  }
  facet_first = if(two_mirrored) NULL else group_label[[groups[1]]]
  .behind(p, .gene_layers(marks, start, end, base_size, with_labels = gene_labels,
                          label_facet = facet_first, facet_levels = group_label[groups])) +
    scale_x_continuous(limits = c(start, end) / 1e6, expand = c(0, 0)) +
    labs(x = NULL, y = y_label) +
    theme_classic(base_size = base_size) +
    theme(strip.background = element_blank(), strip.placement = "outside",
          strip.text.y.left = if(strips) element_text(angle = 90) else element_blank(),
          legend.position = "right")
}

# Breakpoints (both ends of each SV) in a window: ID, pos and SV type
.sv_breakpoints = function(sv_data, chrom, start, end){
  sv_data = as.data.frame(sv_data)
  needed = c("CHROM_A", "START_A", "CHROM_B", "START_B")
  if(!all(needed %in% colnames(sv_data))){
    stop("sv_data needs columns CHROM_A, START_A, CHROM_B and START_B")
  }
  ids = .sv_sample_ids(sv_data)
  chrom_a = gsub("^chr", "", as.character(sv_data$CHROM_A))
  chrom_b = gsub("^chr", "", as.character(sv_data$CHROM_B))
  # SV type: a type column, the Manta name, or translocation/other
  type_col = intersect(c("type", "SVTYPE", "sv_type"), colnames(sv_data))[1]
  type = if(!is.na(type_col)) toupper(as.character(sv_data[[type_col]]))
         else if("manta_name" %in% colnames(sv_data))
           sub("^Manta([A-Z]+).*$", "\\1", as.character(sv_data$manta_name))
         else rep("OTHER", nrow(sv_data))
  type[chrom_a != chrom_b] = "TRA"
  end_of = function(side, chr, partner){
    s = as.numeric(sv_data[[paste0("START_", side)]])
    e = if(paste0("END_", side) %in% colnames(sv_data))
      as.numeric(sv_data[[paste0("END_", side)]]) else s
    data.frame(ID = ids, chrom = chr, pos = (s + e) / 2, type = type,
               partner = partner, input_row = seq_along(ids), stringsAsFactors = FALSE)
  }
  bp = rbind(end_of("A", chrom_a, chrom_b), end_of("B", chrom_b, chrom_a))
  bp[bp$chrom == chrom & bp$pos >= start & bp$pos <= end,
     c("ID", "pos", "type", "partner", "input_row"), drop = FALSE]
}

# Class (and colours) for SV breakpoints: "type" (SV type) or "partner"
# (chromosome of the other end, in GAMBLR's chromosome colours)
.sv_classes = function(points, colour_by, chrom = NULL){
  if(is.null(colour_by)) return(list(class = NULL, colours = NULL))
  if(colour_by == "type"){
    types = intersect(c("DEL", "DUP", "INV", "INS", "TRA", "BND", "OTHER"),
                      c(unique(points$type), "OTHER"))
    return(list(class = factor(ifelse(points$type %in% types, points$type, "OTHER"),
                               levels = types), colours = NULL))
  }
  if(colour_by != "partner"){
    stop("colour_by must be \"type\" or \"partner\"")
  }
  chr = paste0("chr", points$partner)
  order_key = suppressWarnings(as.numeric(sub("^chr", "", chr)))
  order_key[is.na(order_key)] = c(X = 23, Y = 24)[sub("^chr", "", chr[is.na(order_key)])]
  # legend in chromosome order (chr1, chr2, ... chr10, ... chrX, chrY)
  levels = unique(chr[order(order_key, chr, na.last = TRUE)])
  palette = GAMBLR.helpers::get_gambl_colours("chromosome")
  colours = stats::setNames(unname(palette[levels]), levels)
  # chromosomes without a usable colour (missing or NA in the palette) in grey
  colours[is.na(colours) | colours == ""] = "grey60"
  # counting priority: the window's own chromosome (intrachromosomal SVs)
  # last, so a sample with both counts under its translocation partner
  own = paste0("chr", chrom)
  rank = match(chr, c(setdiff(levels, own), intersect(own, levels)))

  list(class = factor(chr, levels = levels), colours = colours, priority = rank)
}

.sv_sample_ids = function(sv_data){
  sv_data = as.data.frame(sv_data)
  col = intersect(c("tumour_sample_id", "Tumor_Sample_Barcode", "sample_id", "ID"),
                  colnames(sv_data))[1]
  if(is.na(col)){
    stop("sv_data needs a sample column (tumour_sample_id, sample_id or ID)")
  }
  as.character(sv_data[[col]])
}

# Mutation positions in a window, with a class for colour_by
.ssm_points = function(maf_data, chrom, start, end, colour_by = NULL, max_classes = 8){
  maf = as.data.frame(maf_data)
  # the reduced format of get_ssm_by_regions(): start, sample_id, region_name
  if(!"Start_Position" %in% colnames(maf) &&
     all(c("start", "sample_id", "region_name") %in% colnames(maf))){
    maf$Tumor_Sample_Barcode = maf$sample_id
    maf$Start_Position = maf$start
    maf$Chromosome = sub("^(chr)?([^:]+):.*$", "\\2", as.character(maf$region_name))
  }
  needed = c("Tumor_Sample_Barcode", "Chromosome", "Start_Position")
  if(!all(needed %in% colnames(maf))){
    stop("maf_data needs columns Tumor_Sample_Barcode, Chromosome and Start_Position ",
         "(a MAF), or start, sample_id and region_name (as from get_ssm_by_regions())")
  }
  keep = gsub("^chr", "", as.character(maf$Chromosome)) == chrom &
    maf$Start_Position >= start & maf$Start_Position <= end
  input_row = which(keep)
  maf = maf[keep, , drop = FALSE]
  class = NULL
  if(!is.null(colour_by)){
    class = .ssm_classes(maf, colour_by)
    # keep the most common classes; the rest are "other"
    top = names(sort(table(class), decreasing = TRUE))
    if(length(top) > max_classes){
      class[!class %in% top[seq_len(max_classes - 1)]] = "other"
    }
    levels = c(setdiff(names(sort(table(class), decreasing = TRUE)), c("other", "NA")),
               intersect(c("other", "NA"), unique(class)))
    class = factor(class, levels = levels)
  }
  out = data.frame(ID = as.character(maf$Tumor_Sample_Barcode),
                   pos = as.numeric(maf$Start_Position), input_row = input_row,
                   stringsAsFactors = FALSE)
  if(!is.null(class)) out$class = class
  out
}

# Classes for colour_by: "coding", "substitution" or a MAF column (e.g. a
# motif column added by GAMBLR.utils::annotate_ssm_context)
.ssm_classes = function(maf, colour_by){
  if(identical(colour_by, "coding")){
    noncoding = c("Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank", "IGR", "RNA",
                  "lincRNA", "Targeted_Region")
    if(!"Variant_Classification" %in% colnames(maf)){
      stop("colour_by = \"coding\" needs Variant_Classification")
    }
    vc = as.character(maf$Variant_Classification)
    return(ifelse(vc %in% noncoding, "non-coding", "coding"))
  }
  if(identical(colour_by, "substitution")){
    if(!all(c("Reference_Allele", "Tumor_Seq_Allele2") %in% colnames(maf))){
      stop("colour_by = \"substitution\" needs Reference_Allele and Tumor_Seq_Allele2")
    }
    ref = toupper(as.character(maf$Reference_Allele))
    alt = toupper(as.character(maf$Tumor_Seq_Allele2))
    snv = nchar(ref) == 1 & nchar(alt) == 1 & ref %in% c("A", "C", "G", "T") &
      alt %in% c("A", "C", "G", "T")
    comp = c(A = "T", C = "G", G = "C", T = "A")
    purine = ref %in% c("A", "G")
    ref[snv & purine] = comp[ref[snv & purine]]
    alt[snv & purine] = comp[alt[snv & purine]]
    return(ifelse(snv, paste0(ref, ">", alt), "indel"))
  }
  if(!colour_by %in% colnames(maf)){
    stop("colour_by must be \"coding\", \"substitution\" or a column of maf_data")
  }
  v = as.character(maf[[colour_by]])
  v[is.na(v)] = "NA"
  v
}

#' @title Plot structural variant breakpoint density in a region
#'
#' @description Show, for one genomic window, how often samples (or
#' patients) have a structural variant breakpoint in each bin, optionally
#' split into groups. A standalone version of the breakpoint track of
#' [plot_cohort_CN_region]; it can be stacked with [plot_ssm_density] or other
#' region tracks using [stack_region_tracks].
#'
#' @details Either end of an SV in the window counts, so translocations into
#' the window are included. With `measure = "patients"` (default), each
#' patient counts once per group (their samples in the same group are
#' combined; samples in different groups count in each), using the
#' patient_id column of `these_samples_metadata`; "samples" counts samples.
#' The denominator is the samples (patients) with SV data in each group. By
#' default all samples in `these_samples_metadata` are assumed to have SV
#' data, including those with no SVs; give `sv_samples` to say otherwise
#' (e.g. only the genomes, since capture data have no SV calls). Without
#' metadata, the samples appearing in `sv_data` are used. SV types come from a type column, the
#' Manta name, or are inferred (different chromosomes: TRA).
#'
#' @param sv_data SV calls: one row per SV with the sample (tumour_sample_id,
#' sample_id or ID) and CHROM_A, START_A, (END_A), CHROM_B, START_B,
#' (END_B), as from GAMBLR's SV functions (e.g. `get_manta_sv()`). This
#' function doesn't retrieve data itself.
#' @param region,chrom,start,end,genes,flank,label_genes The window and the
#' genes to mark, as in [plot_cohort_CN_region]; data frames of genes may have a
#' strand column.
#' @param gene_models Gene models for the gene track (exons, introns and
#' strand). By default (NULL) they are retrieved for the window from
#' mygene.info with [GAMBLR.utils::get_gene_models] (needs internet; without
#' it genes are drawn as plain bars, with a message). Give a data frame in
#' that format to use your own (e.g. edited or cached) models, or FALSE to
#' skip them. Genes in the window that aren't in `genes` or `label_genes`
#' are added to the track (lighter; when retrieved automatically, only if
#' there are 30 or fewer) but not shaded in the plot.
#' @param these_samples_metadata Optional metadata (sample_id, patient_id and
#' the `split_by` column); when given, only its samples are used.
#' @param split_by,split_order Optional metadata column defining groups, and
#' their order.
#' @param sv_samples Samples with SV data (the denominator). Default: the
#' samples in `these_samples_metadata` with `sv_available` TRUE when it has
#' that column (as from `GAMBLR.results::get_data_by_region`), otherwise all
#' of them. See Details.
#' @param measure "patients" (default) or "samples".
#' @param bin_size Bin width in bp. Overrides `bin_fraction`.
#' @param bin_fraction Bin width as a fraction of the window. Default 0.01.
#' @param max_bins Refuse to use more bins than this (a tiny bin size over a
#' large window is slow and unreadable). Default 2000.
#' @param colour_by Optionally stack the bars by "type" (SV type) or
#' "partner" (the chromosome the other end of each SV lands on, in GAMBLR's
#' IGV-style chromosome colours; intrachromosomal SVs get the window's own
#' chromosome). The legend is in chromosome order. A sample (patient) counts
#' once per bin, in one class: for "partner", a translocation partner takes
#' precedence over the window's own chromosome, then chromosome order; for
#' "type", legend order. Default NULL (one colour).
#' @param outliers How samples with outlying numbers of events (e.g. a
#' hypermutated sample, which would dominate its group) are handled. They
#' are excluded from the plot and its denominator, with a message, by:
#' "gap" (default): samples (within the top 5\\%) whose count is more than
#' `outlier_ratio` times that of the next sample below them, i.e. isolated
#' far above everyone else. A heavy but continuous tail, as at aSHM loci
#' (e.g. BCL2, MYC), is kept.
#' "fence": counts above Q3 + `outlier_k` x IQR (on a log scale) among
#' samples with at least one event; stricter, and can remove real tails.
#' "none": keep every sample.
#' Counts are per sample (see `outlier_scope`); samples with
#' `outlier_min_events` or fewer are never outliers.
#' @param outlier_ratio For "gap": how many times the next count a sample
#' must exceed. Default 3.
#' @param outlier_k For "fence": the multiplier. Default 1.5.
#' @param outlier_min_events Samples with this many events or fewer are never
#' outliers. Default 20.
#' @param outlier_scope Count each sample's events in the "window" (default)
#' or in all the "data" given (e.g. genome-wide, when a whole MAF is given).
#' @param colour One colour for all the bars (when not stacked by class).
#' Default NULL: each group in its colour from GAMBLR's palette for the
#' `split_by` column ([map_metadata_to_colours]), or `group_colours`, and
#' grey for groups without one (and without groups).
#' @param group_colours Optional named colours for the groups, overriding the
#' palette.
#' @param class_colours Optional named colours for the classes (SV types, or
#' partner chromosomes named like "chr14").
#' @param mirror With two groups, draw the second downward. Default TRUE.
#' @param y_scale "shared" (default): all groups on one y scale, to compare
#' frequencies directly; "free": each group's panel (or, mirrored, each
#' half) scaled to its own maximum, which its axis label shows.
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
#' the plot, from start to end (narrow ones as a dashed line). Default TRUE;
#' FALSE leaves them only in the gene track.
#' @param genome_build For locating genes. Default "grch37".
#' @param base_size Base font size. Default 9.
#' @param gene_track Draw the genes (`genes`, `label_genes`) as a track above
#' the plot: each from its start to its end, with strand when `genes` is a
#' data frame with a strand column ("+"/"-"). Default TRUE. Set to FALSE
#' when stacking with [stack_region_tracks], which can add one gene track
#' for all the plots. Genes are also shaded in the plot.
#' @param return_data Set to TRUE to return a list with `plot`,
#' `density_plot` (without the gene track), `density` (per bin, group and
#' class: n, n_units and value), `outliers` (samples excluded: ID, events,
#' threshold) and `events`: the rows of `sv_data` behind
#' the plot (SVs with a breakpoint in the window, of counted samples), with
#' their `group`, `class` (with `colour_by`) and `breakpoint_in_window`, and
#' with `measure = "patients"` the `patient_id` they count under (each
#' sample's SVs are kept, since breakpoints of the same event can differ
#' between samples).
#'
#' @return A ggplot, or with `return_data = TRUE` a list with `plot`, `region`
#' (the window shown, as "chrN:start-end"; see [region_window]), `density_plot`,
#' `density`, `events` and `outliers`.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' sv = get_manta_sv(these_samples_metadata = genome_meta)
#' plot_sv_density(sv, genes = "BCL6", flank = 2e6,
#'                 these_samples_metadata = meta, split_by = "COO_consensus",
#'                 sv_samples = genome_meta$sample_id)
#' }
plot_sv_density = function(sv_data, region = NULL, chrom = NULL, start = NULL,
                           end = NULL, genes = NULL, flank = 1e6, label_genes = NULL,
                           gene_models = NULL, genes_show_isoforms = NULL,
                           isoforms = "distinct_ends", shade_genes = TRUE, other_genes = "auto",
                           these_samples_metadata = NULL, split_by = NULL,
                           split_order = NULL, sv_samples = NULL,
                           measure = "patients", bin_size = NULL,
                           bin_fraction = 0.01, max_bins = 2000,
                           colour_by = NULL, colour = NULL, group_colours = NULL,
                           y_scale = "shared",
                           outliers = "gap", outlier_ratio = 3, outlier_k = 1.5,
                           outlier_min_events = 20, outlier_scope = "window",
                           class_colours = NULL, mirror = TRUE,
                           genome_build = "grch37", base_size = 9,
                           gene_track = TRUE, return_data = FALSE){
  measure = match.arg(measure, c("patients", "samples"))
  y_scale = match.arg(y_scale, c("shared", "free"))
  isoforms_given = !missing(isoforms)
  isoforms = match.arg(isoforms, c("distinct_ends", "all"))
  if(isoforms_given && length(genes_show_isoforms) == 0){
    warning("isoforms only applies to genes_show_isoforms, which is not set")
  }
  outliers = match.arg(outliers, c("gap", "fence", "none"))
  outlier_scope = match.arg(outlier_scope, c("window", "data"))
  w = .resolve_window(region, chrom, start, end, genes, flank, label_genes, genome_build)
  w = .add_gene_models(w, .resolve_gene_models(w, gene_models, genome_build, gene_track,
                                               genes_show_isoforms),
                       max_extra_genes = .other_gene_cap(other_genes, gene_models),
                       genes_show_isoforms = genes_show_isoforms, isoforms = isoforms)
  if(!shade_genes && !is.null(w$marks)) w$marks$shade = FALSE
  breaks = .region_bins(w$start, w$end, bin_size, bin_fraction, max_bins)
  profiled = .profiled_samples(sv_samples, these_samples_metadata, .sv_sample_ids(sv_data),
                               column = "sv_available")
  units = .region_units(profiled, these_samples_metadata, split_by, split_order, measure)
  if(nrow(units$map) == 0){
    stop("no samples with SV data (among these_samples_metadata)")
  }
  points = .sv_breakpoints(sv_data, w$chrom, w$start, w$end)
  classes = .sv_classes(points, colour_by, w$chrom)
  if(!is.null(classes$class)) points$class = classes$class
  if(!is.null(classes$priority)) points$priority = classes$priority
  if(is.null(class_colours)) class_colours = classes$colours
  outlier_table = .find_outliers(points$ID[points$ID %in% units$map$ID],
                                 .sv_sample_ids(sv_data)[.sv_sample_ids(sv_data) %in% units$map$ID],
                                 outliers, outlier_k, outlier_min_events, outlier_scope,
                                 outlier_ratio)
  dropped = .drop_outliers(units, points, outlier_table, "SV breakpoints")
  units = dropped$units; points = dropped$points
  events = .contributing_rows(sv_data, points, units,
                              extra = list(breakpoint_in_window = "pos"))
  points$type = NULL; points$partner = NULL
  density = .binned_density(points, units, breaks, measure)
  n_units = tapply(density$n_units, density$group, `[`, 1)[units$groups]
  group_label = stats::setNames(paste0(units$groups, "\n(n=", n_units, ")"), units$groups)
  p = .density_track(density, units$groups, group_label, w$start, w$end, w$marks,
                     y_label = "SV frequency", measure = measure,
                     fill = if(is.null(colour)) "grey35" else colour,
                     group_colours = if(is.null(colour)) .group_colours(units$groups, these_samples_metadata,
                                                                        split_by, group_colours),
                     class_colours = class_colours, mirror = mirror, y_scale = y_scale,
                     base_size = base_size) +
    labs(x = paste0("chr", w$chrom, " (Mb)"))
  combined = .with_gene_track(p, w, gene_track, base_size)
  if(return_data) list(plot = combined, region = w$region, density_plot = p, density = density,
                       events = events, outliers = outlier_table) else combined
}

#' @title Plot simple somatic mutation density in a region
#'
#' @description Show, for one genomic window, how often samples (or
#' patients) have a mutation in each bin, or how many mutations they carry,
#' optionally split into groups and coloured by a mutation class. It can be
#' stacked with [plot_sv_density] or other region tracks using
#' [stack_region_tracks], and is also used by [plot_cohort_CN_region].
#'
#' @details Measures:
#' \itemize{
#'   \item "patients" (default): the fraction of patients in each group with
#'     at least one mutation in the bin. A patient's samples in the same group
#'     are combined (samples in different groups count in each), using the
#'     patient_id column of `these_samples_metadata`.
#'   \item "samples": the same, counting samples.
#'   \item "mutations": mutations per sample per kb, e.g. for aSHM loci,
#'     where one sample can carry many mutations.
#' }
#' The denominator is the samples with mutation data covering the window. By
#' default all samples in `these_samples_metadata` are assumed to have it,
#' including those with no mutations there; give `ssm_samples` to say
#' otherwise. Exome and capture data only cover their targets, so for
#' non-coding sequence give the genomes as `ssm_samples`. Without metadata,
#' the samples appearing in `maf_data` are used (with a message).
#'
#' `colour_by` stacks the bars by: "coding" (coding vs non-coding, from
#' Variant_Classification), "substitution" (the six pyrimidine-centred
#' substitution classes, and indels) or any column of `maf_data`, such as the
#' motif column added by [GAMBLR.utils::annotate_ssm_context] (e.g. WRCY:
#' SITE / MOTIF / FALSE) or a trinucleotide context column. A unit counts
#' once per bin, in the first class by legend order (most common first);
#' classes beyond `max_classes` are combined as "other".
#'
#' @param maf_data Mutations: a MAF (Tumor_Sample_Barcode, Chromosome,
#' Start_Position, plus the columns `colour_by` needs), e.g. from
#' `get_ssm_by_regions(..., streamlined = FALSE)`, or the reduced format of
#' `get_ssm_by_regions()` (start, sample_id, region_name; no `colour_by`).
#' This function doesn't retrieve data itself.
#' @param region,chrom,start,end,genes,flank,label_genes The window and the
#' genes to mark, as in [plot_cohort_CN_region]; data frames of genes may have a
#' strand column.
#' @param gene_models Gene models for the gene track (exons, introns and
#' strand). By default (NULL) they are retrieved for the window from
#' mygene.info with [GAMBLR.utils::get_gene_models] (needs internet; without
#' it genes are drawn as plain bars, with a message). Give a data frame in
#' that format to use your own (e.g. edited or cached) models, or FALSE to
#' skip them. Genes in the window that aren't in `genes` or `label_genes`
#' are added to the track (lighter; when retrieved automatically, only if
#' there are 30 or fewer) but not shaded in the plot.
#' @param these_samples_metadata Optional metadata (sample_id, patient_id and
#' the `split_by` column); when given, only its samples are used.
#' @param split_by,split_order Optional metadata column defining groups, and
#' their order.
#' @param ssm_samples Samples with mutation data covering the window (the
#' denominator). Default: the samples in `these_samples_metadata` with
#' `ssm_available` TRUE when it has that column (as from
#' `GAMBLR.results::get_data_by_region`), otherwise all of them. See
#' Details.
#' @param measure "patients" (default), "samples" or "mutations".
#' @param colour_by Optional class to stack by. See Details.
#' @param max_classes Most classes shown; the rest are "other". Default 8.
#' @param outliers How samples with outlying numbers of events (e.g. a
#' hypermutated sample, which would dominate its group) are handled. They
#' are excluded from the plot and its denominator, with a message, by:
#' "gap" (default): samples (within the top 5\\%) whose count is more than
#' `outlier_ratio` times that of the next sample below them, i.e. isolated
#' far above everyone else. A heavy but continuous tail, as at aSHM loci
#' (e.g. BCL2, MYC), is kept.
#' "fence": counts above Q3 + `outlier_k` x IQR (on a log scale) among
#' samples with at least one event; stricter, and can remove real tails.
#' "none": keep every sample.
#' Counts are per sample (see `outlier_scope`); samples with
#' `outlier_min_events` or fewer are never outliers.
#' @param outlier_ratio For "gap": how many times the next count a sample
#' must exceed. Default 3.
#' @param outlier_k For "fence": the multiplier. Default 1.5.
#' @param outlier_min_events Samples with this many events or fewer are never
#' outliers. Default 20.
#' @param outlier_scope Count each sample's events in the "window" (default)
#' or in all the "data" given (e.g. genome-wide, when a whole MAF is given).
#' @param bin_size Bin width in bp. Overrides `bin_fraction`.
#' @param bin_fraction Bin width as a fraction of the window. Default 0.01.
#' @param max_bins Refuse to use more bins than this. Default 2000.
#' @param colour One colour for all the bars (when not stacked by class).
#' Default NULL: each group in its colour from GAMBLR's palette for the
#' `split_by` column ([map_metadata_to_colours]), or `group_colours`, and
#' grey for groups without one (and without groups).
#' @param group_colours Optional named colours for the groups, overriding the
#' palette.
#' @param class_colours Optional named colours for the classes.
#' @param mirror With two groups, draw the second downward. Default TRUE.
#' @param y_scale "shared" (default): all groups on one y scale, to compare
#' frequencies directly; "free": each group's panel (or, mirrored, each
#' half) scaled to its own maximum, which its axis label shows.
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
#' the plot, from start to end (narrow ones as a dashed line). Default TRUE;
#' FALSE leaves them only in the gene track.
#' @param genome_build For locating genes. Default "grch37".
#' @param base_size Base font size. Default 9.
#' @param gene_track Draw the genes as a track above the plot (see
#' [plot_sv_density]). Default TRUE.
#' @param return_data Set to TRUE to return a list with `plot`,
#' `density_plot` (without the gene track), `density`, `outliers` (samples
#' excluded: ID, events, threshold) and `events`: the rows
#' of `maf_data` behind the plot (mutations in the window, of counted
#' samples), with their `group` and `class` (with `colour_by`). With
#' `measure = "patients"`, also `patient_id`, and a variant found in more
#' than one of a patient's samples in the same group (same position and
#' alleles) is listed once (the first sample's row), so the rows are the
#' patients' distinct mutations.
#'
#' @return A ggplot, or with `return_data = TRUE` a list with `plot`, `region`
#' (the window shown, as "chrN:start-end"; see [region_window]), `density_plot`,
#' `density`, `events` and `outliers`.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' plot_ssm_density(maf, genes = "BCL6", flank = 5e4, measure = "mutations",
#'                  these_samples_metadata = meta, split_by = "COO_consensus",
#'                  ssm_samples = genome_ids, colour_by = "coding")
#' }
plot_ssm_density = function(maf_data, region = NULL, chrom = NULL, start = NULL,
                            end = NULL, genes = NULL, flank = 5e4, label_genes = NULL,
                            gene_models = NULL, genes_show_isoforms = NULL,
                           isoforms = "distinct_ends", shade_genes = TRUE, other_genes = "auto",
                            these_samples_metadata = NULL, split_by = NULL,
                            split_order = NULL, ssm_samples = NULL,
                            measure = "patients", colour_by = NULL, max_classes = 8,
                            outliers = "gap", outlier_ratio = 3, outlier_k = 1.5,
                            outlier_min_events = 20, outlier_scope = "window",
                            bin_size = NULL, bin_fraction = 0.01, max_bins = 2000,
                            colour = NULL, group_colours = NULL, class_colours = NULL,
                            mirror = TRUE, y_scale = "shared",
                            genome_build = "grch37", base_size = 9,
                            gene_track = TRUE, return_data = FALSE){
  measure = match.arg(measure, c("patients", "samples", "mutations"))
  y_scale = match.arg(y_scale, c("shared", "free"))
  isoforms_given = !missing(isoforms)
  isoforms = match.arg(isoforms, c("distinct_ends", "all"))
  if(isoforms_given && length(genes_show_isoforms) == 0){
    warning("isoforms only applies to genes_show_isoforms, which is not set")
  }
  outliers = match.arg(outliers, c("gap", "fence", "none"))
  outlier_scope = match.arg(outlier_scope, c("window", "data"))
  w = .resolve_window(region, chrom, start, end, genes, flank, label_genes, genome_build)
  w = .add_gene_models(w, .resolve_gene_models(w, gene_models, genome_build, gene_track,
                                               genes_show_isoforms),
                       max_extra_genes = .other_gene_cap(other_genes, gene_models),
                       genes_show_isoforms = genes_show_isoforms, isoforms = isoforms)
  if(!shade_genes && !is.null(w$marks)) w$marks$shade = FALSE
  breaks = .region_bins(w$start, w$end, bin_size, bin_fraction, max_bins)
  profiled = .profiled_samples(ssm_samples, these_samples_metadata,
                               .ssm_sample_ids(maf_data), column = "ssm_available")
  units = .region_units(profiled, these_samples_metadata, split_by, split_order,
                        if(measure == "mutations") "samples" else measure)
  if(nrow(units$map) == 0){
    stop("no samples with mutation data (among these_samples_metadata)")
  }
  points = .ssm_points(maf_data, w$chrom, w$start, w$end, colour_by, max_classes)
  all_ids = .ssm_sample_ids(maf_data)
  outlier_table = .find_outliers(points$ID[points$ID %in% units$map$ID],
                                 all_ids[all_ids %in% units$map$ID],
                                 outliers, outlier_k, outlier_min_events, outlier_scope,
                                 outlier_ratio)
  dropped = .drop_outliers(units, points, outlier_table, "mutations")
  units = dropped$units; points = dropped$points
  events = .contributing_rows(maf_data, points, units,
                              dedupe = intersect(c("Chromosome", "Start_Position", "End_Position",
                                                   "Reference_Allele", "Tumor_Seq_Allele2"),
                                                 colnames(as.data.frame(maf_data))))
  density = .binned_density(points, units, breaks, measure)
  n_units = tapply(density$n_units, density$group, `[`, 1)[units$groups]
  group_label = stats::setNames(paste0(units$groups, "\n(n=", n_units, ")"), units$groups)
  y_label = if(measure == "mutations") "mutations per sample per kb"
            else "SSM frequency"
  p = .density_track(density, units$groups, group_label, w$start, w$end, w$marks,
                     y_label = y_label, measure = measure,
                     fill = if(is.null(colour)) "grey35" else colour,
                     group_colours = if(is.null(colour)) .group_colours(units$groups, these_samples_metadata,
                                                                        split_by, group_colours),
                     class_colours = class_colours, mirror = mirror, y_scale = y_scale,
                     base_size = base_size) +
    labs(x = paste0("chr", w$chrom, " (Mb)"))
  combined = .with_gene_track(p, w, gene_track, base_size)
  if(return_data) list(plot = combined, region = w$region, density_plot = p, density = density,
                       events = events, outliers = outlier_table) else combined
}

#' @title Stack region tracks
#'
#' @description Stack plots of the same window, e.g. from [plot_sv_density]
#' and [plot_ssm_density], vertically with their x axes aligned. The x axis
#' is shown on the bottom plot only.
#'
#' @param ... ggplots of the same window (or a single list of them), e.g.
#' made with `gene_track = FALSE`.
#' @param gene_track Optional gene track for the top, from
#' [plot_region_genes] (with `return_data = TRUE`, or its plot).
#' @param rel_heights Relative heights of the plots in `...`. Default equal.
#'
#' @return A combined plot (cowplot).
#'
#' @import cowplot ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' stack_region_tracks(plot_sv_density(sv, genes = "BCL2", flank = 1e5, ...),
#'                     plot_ssm_density(maf, genes = "BCL2", flank = 1e5, ...))
#' }
stack_region_tracks = function(..., gene_track = NULL, rel_heights = NULL){
  plots = list(...)
  if(length(plots) == 1 && is.list(plots[[1]]) && !inherits(plots[[1]], "gg")){
    plots = plots[[1]]
  }
  n = length(plots)
  for(i in seq_len(n - 1)){
    plots[[i]] = plots[[i]] + theme(axis.text.x = element_blank(),
                                    axis.title.x = element_blank())
  }
  heights = if(is.null(rel_heights)) rep(1, n) else rel_heights
  if(!is.null(gene_track)){
    rows = if(is.list(gene_track) && !inherits(gene_track, "gg")) gene_track$rows else 1
    track = if(is.list(gene_track) && !inherits(gene_track, "gg")) gene_track$plot else gene_track
    plots = c(list(track), plots)
    heights = c(.gene_track_height(rows) * stats::median(heights), heights)
  }
  cowplot::plot_grid(plotlist = plots, ncol = 1, align = "v", axis = "lr",
                     rel_heights = heights)
}

#' @title Plot the genes in a region
#'
#' @description A gene track for a genomic window: each gene drawn from its
#' start to its end, with strand when known and its name above, overlapping
#' genes on separate rows. Used at the top of [plot_cohort_CN_region],
#' [plot_sv_density] and [plot_ssm_density], and with [stack_region_tracks].
#'
#' @param region,chrom,start,end,genes,flank,label_genes The window and the
#' genes, as in [plot_cohort_CN_region]. `genes` (drawn darker) and
#' `label_genes` can be gene symbols or data frames with gene, chrom, start,
#' end and optionally strand ("+"/"-").
#' @param gene_models Gene models for the gene track (exons, introns and
#' strand). By default (NULL) they are retrieved for the window from
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
#' the plot, from start to end (narrow ones as a dashed line). Default TRUE;
#' FALSE leaves them only in the gene track.
#' @param genome_build For locating genes given as symbols. Default "grch37".
#' @param base_size Base font size. Default 9.
#' @param return_data Set to TRUE to return a list with `plot` and `rows`
#' (the number of rows, for sizing it in a stack).
#'
#' @return A ggplot (NULL if no genes are in the window), or a list.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' genes = data.frame(gene = c("BCL6", "LPP"), chrom = "3",
#'                    start = c(187439165, 187871072), end = c(187463515, 188608460),
#'                    strand = c("-", "+"))
#' plot_region_genes(region = "3:187000000-189000000", genes = genes)
#' }
plot_region_genes = function(region = NULL, chrom = NULL, start = NULL, end = NULL,
                             genes = NULL, flank = 1e6, label_genes = NULL,
                             gene_models = NULL, genes_show_isoforms = NULL,
                           isoforms = "distinct_ends", shade_genes = TRUE, other_genes = "auto",
                             genome_build = "grch37", base_size = 9,
                             return_data = FALSE){
  isoforms_given = !missing(isoforms)
  isoforms = match.arg(isoforms, c("distinct_ends", "all"))
  if(isoforms_given && length(genes_show_isoforms) == 0){
    warning("isoforms only applies to genes_show_isoforms, which is not set")
  }
  w = .resolve_window(region, chrom, start, end, genes, flank, label_genes, genome_build)
  w = .add_gene_models(w, .resolve_gene_models(w, gene_models, genome_build,
                                               genes_show_isoforms = genes_show_isoforms),
                       max_extra_genes = .other_gene_cap(other_genes, gene_models),
                       genes_show_isoforms = genes_show_isoforms, isoforms = isoforms)
  if(!shade_genes && !is.null(w$marks)) w$marks$shade = FALSE
  track = .gene_track(w$marks, w$start, w$end, base_size, w$exons)
  if(return_data) track else if(is.null(track)) NULL else track$plot
}

# A plot with a gene track above it (when asked for and there are genes)
.with_gene_track = function(p, w, gene_track, base_size){
  if(!gene_track) return(p)
  track = .gene_track(w$marks, w$start, w$end, base_size, w$exons)
  if(is.null(track)) return(p)
  cowplot::plot_grid(track$plot, p, ncol = 1, align = "v", axis = "lr",
                     rel_heights = c(.gene_track_height(track$rows), 1))
}

# The input rows behind a density track: those of samples counted (with
# data, in a group) that fall in the window, annotated with their group and
# class (and for SVs the breakpoint position in the window)
.contributing_rows = function(input, points, units, extra = NULL, dedupe = NULL){
  input = as.data.frame(input)
  p = points[points$ID %in% units$map$ID, , drop = FALSE]
  if(nrow(p) == 0) return(input[0, , drop = FALSE])
  rows = unique(p$input_row)
  out = input[rows, , drop = FALSE]
  first = p[match(rows, p$input_row), , drop = FALSE]
  out$group = units$map$group[match(first$ID, units$map$ID)]
  # the unit each row was counted under: its patient (or sample)
  unit = units$map$unit[match(first$ID, units$map$ID)]
  if(units$counted == "patients"){
    out$patient_id = sub("^[^\r]*\r", "", unit)
    # a variant found in several samples of a patient (same group) once
    if(!is.null(dedupe) && all(dedupe %in% colnames(out))){
      key = do.call(paste, c(list(unit), out[dedupe], sep = "\r"))
      keep = !duplicated(key)
      out = out[keep, , drop = FALSE]; first = first[keep, , drop = FALSE]
    }
  }
  if(!is.null(p$class)) out$class = as.character(first$class)
  for(col in names(extra)) out[[col]] = first[[extra[[col]]]]
  rownames(out) = NULL
  out
}

# Add layers underneath a plot's existing layers (e.g. gene spans behind data)
.behind = function(p, layers){
  if(length(layers) == 0) return(p)
  n = length(p$layers)
  p = p + layers
  k = length(p$layers)
  if(k > n && n > 0) p$layers = c(p$layers[(n + 1):k], p$layers[seq_len(n)])
  p
}

# Samples with data (a track's denominator): those given, otherwise all
# samples in the metadata (the data are assumed to cover them, including
# samples without events), otherwise the samples appearing in the data
.profiled_samples = function(given, these_samples_metadata, data_ids, column = NULL){
  if(!is.null(given)) return(unique(as.character(given)))
  if(!is.null(these_samples_metadata)){
    meta = as.data.frame(these_samples_metadata)
    # availability annotated on the metadata (e.g. by
    # GAMBLR.results::get_data_by_region): only those samples
    if(!is.null(column) && column %in% colnames(meta)){
      use = meta[[column]] %in% TRUE
      if(any(!use)){
        message(sum(use), " of ", nrow(meta), " samples have this data type (", column,
                "); the others are left out of the denominator")
      }
      return(unique(as.character(meta$sample_id[use])))
    }
    return(unique(as.character(meta$sample_id)))
  }
  message("no these_samples_metadata: the denominator is the samples appearing in ",
          "the data (samples without events are missed)")
  unique(as.character(data_ids))
}

# Sample IDs in mutation data (MAF, or the reduced get_ssm_by_regions() format)
.ssm_sample_ids = function(maf_data){
  maf = as.data.frame(maf_data)
  col = intersect(c("Tumor_Sample_Barcode", "sample_id"), colnames(maf))[1]
  if(is.na(col)) return(character(0))
  as.character(maf[[col]])
}

# Samples with outlying numbers of events (e.g. a hypermutated sample),
# from their counts (in the window, or all of the data given) among samples
# with at least one event. "gap" (default): the leading samples (within the
# top 5%) whose count is more than `ratio` times that of the next sample
# below them, i.e. isolated far above the rest, as opposed to a heavy but
# continuous tail (as at aSHM loci). "fence": above Q3 + k * IQR on a log
# scale. Either way only counts above min_events. Returns ID, events and
# threshold (the count above which samples were excluded).
.find_outliers = function(ids_in_window, ids_in_data, method = "gap", k = 1.5,
                          min_events = 20, scope = "window", ratio = 3){
  empty = data.frame(ID = character(0), events = numeric(0), threshold = numeric(0))
  if(method == "none") return(empty)
  ids = if(scope == "data") ids_in_data else ids_in_window
  counts = sort(table(as.character(ids)), decreasing = TRUE)
  n = length(counts)
  if(n < 4) return(empty)
  x = as.numeric(counts)
  if(method == "gap"){
    top = min(n - 1, max(1, ceiling(0.05 * n)))
    cut = 0
    for(i in seq_len(top)){
      if(x[i] > min_events && x[i] > ratio * x[i + 1]) cut = i
    }
    if(cut == 0) return(empty)
    threshold = ratio * x[cut + 1]
    return(data.frame(ID = names(counts)[seq_len(cut)], events = x[seq_len(cut)],
                      threshold = round(max(threshold, min_events), 1),
                      stringsAsFactors = FALSE))
  }
  lx = log1p(x)
  q = stats::quantile(lx, c(0.25, 0.75), names = FALSE)
  threshold = max(expm1(q[2] + k * (q[2] - q[1])), min_events)
  out = x > threshold
  if(!any(out)) return(empty)
  data.frame(ID = names(counts)[out], events = x[out], threshold = round(threshold, 1),
             stringsAsFactors = FALSE)
}

# Drop outlier samples from the units (denominator) and points, with a message
.drop_outliers = function(units, points, outliers, what){
  if(nrow(outliers) == 0) return(list(units = units, points = points))
  message(nrow(outliers), " outlier sample(s) excluded (more than ", outliers$threshold[1],
          " ", what, "): ", paste(utils::head(paste0(outliers$ID, " (", outliers$events, ")"), 10),
                                  collapse = ", "),
          "; see outliers = \"none\"")
  units$map = units$map[!units$map$ID %in% outliers$ID, , drop = FALSE]
  list(units = units, points = points[!points$ID %in% outliers$ID, , drop = FALSE])
}

# Colours for the groups of split_by: given ones, or GAMBLR's palette for
# that column (map_metadata_to_colours); NULL without groups or a palette
.group_colours = function(groups, these_samples_metadata, split_by, given = NULL){
  if(!is.null(given)) return(given)
  if(is.null(split_by) || is.null(these_samples_metadata)) return(NULL)
  col = tryCatch(suppressMessages(map_metadata_to_colours(
    metadataColumns = split_by, these_samples_metadata = as.data.frame(these_samples_metadata)))[[split_by]],
    error = function(e) NULL)
  if(is.null(col) || !any(groups %in% names(col))) return(NULL)
  col[names(col) %in% groups]
}

# How many genes besides the named ones the gene track may add
.other_gene_cap = function(other_genes, gene_models){
  other_genes = match.arg(other_genes, c("auto", "none", "all"))
  switch(other_genes, none = 0, all = Inf,
         auto = if(is.null(gene_models)) 30 else Inf)
}
