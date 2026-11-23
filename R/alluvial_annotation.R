# Helpers for prettyOncoplot(alluvial_columns = ...): a stacked alluvial drawn as a column annotation of the
# oncoplot, so its axes and ribbons share the oncoplot's columns, and the sample order that goes with it.

# The comparable key of each class label: `equivalents` maps labels that mean the same class onto one key
# (e.g. c(C1 = "BN2", C3 = "EZB")); unmapped labels are their own key
.alluvial_key = function(v, equivalents = NULL){
  v = as.character(v)
  v[is.na(v)] = "NA"
  if(!is.null(equivalents)){
    hit = v %in% names(equivalents)
    v[hit] = unname(equivalents[v[hit]])
  }
  v
}

# Rank of each key in the class order: the order given, then any other keys alphabetically, then "Other" and
# "NA" last
.alluvial_rank = function(keys, class_order = NULL){
  default = c("BN2", "EZB", "ST2", "MCD", "N1", "A53", "C2", "Composite")
  if(is.null(class_order)) class_order = default
  rest = sort(setdiff(unique(keys), c(class_order, "Other", "NA")))
  full = c(class_order, rest, "Other", "NA")
  match(keys, full)
}

# Sample order for an alluvial over `columns` (listed top to bottom; the last one sits next to the oncoplot):
# by the class of the last axis; within it, for each axis further away (bottom to top), by `sort` for that axis:
#   "class"      samples whose class there matches the last axis's first, then the other classes in order
#   "frequency"  the classes shared by most samples within the block (all axes nearer the oncoplot) first,
#                ties in class order
# `sort` is one value for every axis above the last, or one per axis (top to bottom, the last excluded)
# `tiebreak`: optional numeric column of `metadata`; the last sort key, highest first (missing values last)
.alluvial_sample_order = function(metadata, columns, class_order = NULL, equivalents = NULL, sort = "class",
                                  tiebreak = NULL){
  k = length(columns)
  if(k > 1){
    if(length(sort) == 1) sort = rep(sort, k - 1)
    if(length(sort) != k - 1 || !all(sort %in% c("class", "frequency"))){
      stop("alluvial sort must be \"class\" or \"frequency\", once or once per axis above the last")
    }
  }
  keys = lapply(columns, function(cl) .alluvial_key(metadata[[cl]], equivalents))
  base = keys[[k]]
  sort_keys = list(.alluvial_rank(base, class_order))
  block = base
  if(k > 1){
    for(j in (k - 1):1){
      cls_rank = .alluvial_rank(keys[[j]], class_order)
      sort_keys[[length(sort_keys) + 1]] = if(sort[j] == "class") ifelse(keys[[j]] == base, 0L, cls_rank) else {
        n_in_block = stats::ave(seq_along(block), block, keys[[j]], FUN = length)
        -n_in_block
      }
      sort_keys[[length(sort_keys) + 1]] = cls_rank   # ties
      block = paste(block, keys[[j]], sep = "\r")
    }
  }
  if(!is.null(tiebreak)){
    tb = suppressWarnings(as.numeric(metadata[[tiebreak]]))
    tb[is.na(tb)] = -Inf
    sort_keys[[length(sort_keys) + 1]] = -tb
  }
  do.call(order, sort_keys)
}

# Consistency of each sample across the axes, relative to the axis next to the oncoplot (the last): "consistent"
# when every other axis has the same class, "inconsistent" when none does, "semi-consistent" otherwise
.alluvial_consistency = function(metadata, columns, equivalents = NULL){
  keys = lapply(columns, function(cl) .alluvial_key(metadata[[cl]], equivalents))
  k = length(keys)
  base = keys[[k]]
  if(k == 1) return(rep("consistent", length(base)))
  agree = Reduce(`+`, lapply(keys[-k], function(x) as.integer(x == base)))
  out = ifelse(agree == k - 1, "consistent", ifelse(agree == 0, "inconsistent", "semi-consistent"))
  factor(out, levels = c("consistent", "semi-consistent", "inconsistent"))
}

# The alluvial as a ComplexHeatmap column annotation. `values`: data frame with one column per axis (top to
# bottom), rows in the order of the oncoplot matrix's columns. Each axis is a bar of sample-wide cells coloured
# by class; the last axis follows the oncoplot's column order, and each axis above it orders its samples by its
# own class, then by the class and position on the axis below (so ribbons from one class bundle together).
# Ribbons join each sample's positions on neighbouring axes, coloured by the class on the upper axis
# (ribbon_colour = "upper") or the lower one ("lower"), or, given one of the axis names, every ribbon by that axis's
# class (so one classification's classes can be followed through all the axes).
.alluvial_annotation = function(values, colours, class_order = NULL, equivalents = NULL, height = 4,
                                bar_fraction = 0.1, ribbon_alpha = 0.45, ribbon_colour = "upper",
                                axis_label_fontsize = 7, show_axis_names = TRUE){
  values = as.data.frame(values, stringsAsFactors = FALSE)
  axes = colnames(values)
  if(length(ribbon_colour) != 1 || !ribbon_colour %in% c("upper", "lower", axes)){
    stop("ribbon_colour must be \"upper\", \"lower\" or one of the alluvial columns: ", paste(axes, collapse = ", "))
  }
  k = length(axes)
  n_all = nrow(values)
  colour_of = function(j, v){
    pal = colours[[axes[j]]]
    v = as.character(v); v[is.na(v)] = "NA"
    out = if(is.null(pal)) rep(NA_character_, length(v)) else unname(pal[v])
    out[is.na(out)] = "grey70"
    out
  }
  draw_fun = function(index, k_slice, n_slice){
    n = length(index)
    v = values[index, , drop = FALSE]
    keys = lapply(axes, function(a) .alluvial_key(v[[a]], equivalents))
    # positions on each axis: the last follows the column order; each one above by its own class, then the
    # class and position on the axis below
    pos = vector("list", k)
    pos[[k]] = seq_len(n)
    if(k > 1){
      for(j in (k - 1):1){
        o = order(.alluvial_rank(keys[[j]], class_order), .alluvial_rank(keys[[j + 1]], class_order), pos[[j + 1]])
        pos[[j]] = integer(n); pos[[j]][o] = seq_len(n)
      }
    }
    bh = bar_fraction
    yc = if(k == 1) 1 - bh / 2 else 1 - bh / 2 - (seq_len(k) - 1) * (1 - bh) / (k - 1)
    grid::pushViewport(grid::viewport(xscale = c(0.5, n + 0.5), yscale = c(0, 1)))
    # ribbons between neighbouring axes: one polygon per sample, sigmoid edges
    t = seq(0, 1, length.out = 16); s = t * t * (3 - 2 * t)
    if(k > 1){
      for(j in seq_len(k - 1)){
        y_top = yc[j] - bh / 2; y_bot = yc[j + 1] + bh / 2
        yy = y_top + (y_bot - y_top) * t
        ja = if(ribbon_colour %in% axes) match(ribbon_colour, axes) else if(ribbon_colour == "upper") j else j + 1
        fill = colour_of(ja, v[[axes[ja]]])
        # samples adjacent on both axes with the same colour form one ribbon (no seams between them)
        o = order(pos[[j]]); xt = pos[[j]][o]; xb = pos[[j + 1]][o]; fl = fill[o]
        brk = c(TRUE, diff(xt) != 1 | diff(xb) != 1 | fl[-1] != fl[-length(fl)])
        run = cumsum(brk)
        px = c(); py = c(); id = c(); rf = c()
        for(r in unique(run)){
          w = which(run == r)
          left = (xt[w[1]] - 0.5) + ((xb[w[1]] - 0.5) - (xt[w[1]] - 0.5)) * s
          right = (xt[w[length(w)]] + 0.5) + ((xb[w[length(w)]] + 0.5) - (xt[w[length(w)]] + 0.5)) * s
          px = c(px, left, rev(right)); py = c(py, yy, rev(yy)); id = c(id, rep(r, 2 * length(t))); rf = c(rf, fl[w[1]])
        }
        grid::grid.polygon(x = px, y = py, id = id, default.units = "native",
                           gp = grid::gpar(fill = grDevices::adjustcolor(rf, alpha.f = ribbon_alpha), col = NA))
      }
    }
    # axis bars: one cell per sample, at its position on that axis
    for(j in seq_len(k)){
      # consecutive samples of one class as one block
      o = order(pos[[j]]); xs = pos[[j]][o]; fl = colour_of(j, v[[axes[j]]])[o]
      run = cumsum(c(TRUE, fl[-1] != fl[-length(fl)]))
      starts = tapply(xs, run, min); ends = tapply(xs, run, max); rf = fl[match(unique(run), run)]
      grid::grid.rect(x = (starts + ends) / 2, y = yc[j], width = ends - starts + 1, height = bh, default.units = "native",
                      gp = grid::gpar(fill = rf, col = NA))
      if(show_axis_names){
        grid::grid.text(axes[j], x = grid::unit(1, "npc") + grid::unit(1, "mm"), y = grid::unit(yc[j], "native"),
                        just = "left", gp = grid::gpar(fontsize = axis_label_fontsize))
      }
    }
    grid::popViewport()
  }
  fun = ComplexHeatmap::AnnotationFunction(fun = draw_fun, n = n_all, which = "column", subsettable = FALSE,
                                           height = grid::unit(height, "cm"),
                                           # AnnotationFunction gives the drawing function its own environment:
                                           # everything it uses is passed here
                                           var_import = list(values = values, axes = axes, k = k, class_order = class_order,
                                                             equivalents = equivalents, bar_fraction = bar_fraction,
                                                             ribbon_alpha = ribbon_alpha, ribbon_colour = ribbon_colour,
                                                             axis_label_fontsize = axis_label_fontsize,
                                                             show_axis_names = show_axis_names, colour_of = colour_of))
  ComplexHeatmap::HeatmapAnnotation(alluvial = fun, show_annotation_name = FALSE)
}
