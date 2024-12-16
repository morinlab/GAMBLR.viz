#' @title Heatmap
#'
#' @description Create a highly customizable heatmap using the ComplexHeatmap package.
#'
#' @details Make an heatmap that is looking cute using ComplexHeatmap. The metadata is expected to follow the structure and column naming used in GAMBL.
#' If you provide your own non-GAMBL samples and metadata, you must include at least the columns with names corresponding to annotation tracks and column "Tumor_Sample_Barcode".
#' showing sample ids. The metadata can contain numeric columns, which will be plotted as numeric variables in the annotation. The feature matrix is supplied in `this_matrix` argument.
#' and is expected to have samples in rows, and features in columns. The argument importance_values is similar to the widths of NMF object or importance values for feature/group from RF models.
#' It is also expected to have column names (having names of the groups that will be shown on heatmap) and rownames (corresponding to feature ids).
#'
#' @param this_matrix A data frame with column Tumor_Sample_Barcode and a column for each feature. Can be binary. Expected to not contain negative values.
#' @param importance_values Provide a data frame of feature (in rows) by group (in columns) with numeric values representative of feature importance. Can be obtained from rf$inportance or basis(NMF).
#' @param these_samples_metadata Data frame containing metadata for your samples.
#' @param max_number_of_features_per_group Optional argument to indicate how many features from each group to be considered for display. Default is 10.
#' @param splitColumnName Optional argument to indicate which metadata column to split on. Default is set to pathology.
#' @param metadataColumns A vector containing the categorical column names you want to plot below.
#' @param numericMetadataColumns A vector containing the numeric columns you want to plot below.
#' @param numericMetadataMax A numeric vector of cutoffs to apply to numeric columns above.
#' @param custom_colours Provide named vector (or named list of vectors) containing custom annotation colours if you do not want to use standartized pallette.
#' @param prioritize_ordering_on_numeric Logical argument specifying whether to sort on numeric metadata first or other metadata columns. Default is TRUE (sort on numeric metadata, then on other columns).
#' @param legend_direction Optional argument to indicate whether legend should be in horizontal (default) or vertical position.
#' @param legend_position Optional argument to indicate where the legend should be drawn. The default is set to bottom, but can also accept top, right, and left.
#' @param legend_row Fiddle with these to widen or narrow your legend (default 3).
#' @param legend_col Fiddle with these to widen or narrow your legend (default 3).
#' @param fontSizeGene Font size for gene labels (default 6).
#' @param metadataBarHeight Optional argument to adjust the height of bar with annotations. The default is 1.5
#' @param leftStackedWidth Optional argument to control how wide should the stacked plot on the left be. The default is 4.
#' @param metadataBarFontsize Optional argument to control for the font size of metadata annotations. The default is 5.
#' @param groupNames optional vector of group names to be displayed above heatmap. Should be the same length as the number of groups that will be shown. Default is NULL (no labels).
#'
#' @return Nothing
#'
#' @import dplyr circlize ComplexHeatmap ggplot2 tibble
#' @export
#'
#' @examples
#' \dontrun{
#' splendidHeatmap(this_matrix = data,
#'                 importance_values = rf$importance[,c(1:3)],
#'                 these_samples_metadata = MASTER.METADATA,
#'                 splitColumnName = "pathology",
#'                 metadataColumns = c("cohort",
#'                                     "pathology",
#'                                     "sex",
#'                                     ".",
#'                                     "COO_consensus",
#'                                     "DHITsig_consensus",
#'                                     "seq_type"),
#'                 numericMetadataColumns = ".",
#'                 numericMetadataMax = 0.7,
#'                 custom_colours = custom_colours)
#' }
#'
splendidHeatmap = function(this_matrix,
                           importance_values,
                           these_samples_metadata,
                           max_number_of_features_per_group = 10,
                           splitColumnName = "pathology",
                           metadataColumns = c("pathology"),
                           numericMetadataColumns = NULL,
                           numericMetadataMax = NULL,
                           prioritize_ordering_on_numeric = TRUE,
                           custom_colours = NULL,
                           legend_direction = "horizontal",
                           legend_position = "bottom",
                           legend_row = 3,
                           legend_col = 3,
                           fontSizeGene = 6,
                           metadataBarHeight = 1.5,
                           leftStackedWidth = 4,
                           metadataBarFontsize = 5,
                           groupNames = NULL){
  comparison_groups = colnames(importance_values)

  if(!is.null(splitColumnName) & (splitColumnName %in% metadataColumns)){
    metadataColumns = c(splitColumnName, metadataColumns[!metadataColumns == splitColumnName])
  }

  if(!is.null(numericMetadataColumns) & length(intersect(numericMetadataColumns, metadataColumns))>0){
    message(paste0("The column(s) ", numericMetadataColumns, " specified both in metadata and numeric metadata. Plotting as numeric values..."))
    metadataColumns = metadataColumns[!metadataColumns %in% numericMetadataColumns]
  }

  #get which group samples belong to
  metadata_df = these_samples_metadata[,c("Tumor_Sample_Barcode", metadataColumns, numericMetadataColumns)] %>%
    as.data.frame() %>%
    column_to_rownames(., "Tumor_Sample_Barcode")

  if(!is.null(numericMetadataMax)){
      max_list = setNames(numericMetadataMax, numericMetadataColumns)
      metadata_df = metadata_df %>%
        dplyr::mutate(across(names(max_list), ~ ifelse(.x > max_list[[cur_column()]], max_list[[cur_column()]], .x)))
  }

  # count N of features for every dsample and add it to metadata
  metadata_df =
  this_matrix %>%
    as.data.frame %>%
    column_to_rownames("Tumor_Sample_Barcode") %>%
    rowSums %>%
    as.data.frame %>%
    `names<-`("N_features") %>%
    rownames_to_column ("Tumor_Sample_Barcode") %>%
    base::merge(metadata_df %>%
                  rownames_to_column ("Tumor_Sample_Barcode"),
                .) %>%
    column_to_rownames("Tumor_Sample_Barcode")

  my_colours = NULL
  these_names = NULL
  for (i in 1:length(metadataColumns)){
    this_metadata_column = GAMBLR.helpers::get_gambl_colours(metadataColumns[i])
    if (sum(is.na(names(this_metadata_column[unlist(c(unique(these_samples_metadata[,metadataColumns[i]])))]))) <= 1 &
        nrow(unique(these_samples_metadata[,metadataColumns[i]])) > 1){
      these_names = c(these_names, metadataColumns[i])
      my_colours = append(my_colours, list(c(this_metadata_column, "NA" = "#BDBDC1FF")))
      names(my_colours) = these_names
    }
  }

  my_colours = c(custom_colours, my_colours)

  col_fun=circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  for(exp in numericMetadataColumns){
    my_colours[[exp]] = col_fun
  }

  #get all features
  w = importance_values[,comparison_groups]
  w = as.data.frame(w) %>%
    dplyr::mutate_if(is.character,as.numeric)

  # extract most important features, while taking the feature with highest weight for a particular cluster if it was seen before for other cluster with lower weight
  FEATURES <- w[,1] %>%
    as.data.frame() %>%
    `rownames<-`(rownames(w)) %>%
    `names<-`("importance") %>%
    dplyr::arrange(desc(importance)) %>%
    head(., max_number_of_features_per_group) %>%
    rownames_to_column(., var = "Feature") %>%
    dplyr::mutate(group = comparison_groups[1])
  for (i in 2:length(comparison_groups)){
    FEATURES = rbind(as.data.frame(FEATURES), w[,i] %>%
      as.data.frame() %>%
      `rownames<-`(rownames(w)) %>%
      `names<-`("importance") %>%
       arrange(desc(importance)) %>%
       head(., max_number_of_features_per_group + 3) %>%
       rownames_to_column(., var = "Feature") %>%
       dplyr::mutate(group = comparison_groups[i])) %>%
       dplyr::mutate(importance=as.numeric(importance)) %>%
       dplyr::group_by(Feature) %>%
       dplyr::filter(importance == max(importance)) %>%
       dplyr::arrange(group)
  }
  FEATURES = as.data.frame(FEATURES)

  mat = this_matrix %>%
    base::merge(., metadata_df %>%
    rownames_to_column(., "Tumor_Sample_Barcode") %>%
    dplyr::select(Tumor_Sample_Barcode, splitColumnName)) %>%
    as.data.frame()
  mat[,splitColumnName] = factor(mat[,splitColumnName])

  #breaks used to display groups with different colors on heatmap
  bk = c(0, seq(0.5, length(comparison_groups) + 0.5, 1))

  #colors used to show on the heatmap body. Starts with white - the color of feature absence
  my_palette = c("white", rev(unlist(my_colours[splitColumnName])))
  my_palette = unname(my_palette)

  #get each group and label the events for each feature with group number
  mat_2 = mat[,-ncol(mat)]
  #subset samples of each group
  MY.LIST = list()
  for (i in 1:length(comparison_groups)){
    MY.LIST[[i]] = assign(comparison_groups[i], mat_2 %>%
      as.data.frame(.) %>%
      column_to_rownames(., var = "Tumor_Sample_Barcode") %>%
      t(.) %>%
      as.data.frame(.) %>%
      dplyr::select(metadata_df %>%
      dplyr::filter(base::get(splitColumnName) == comparison_groups[i]) %>%
      rownames))
  }

  #assign numbers - used for coloring of heatmap body
  for(i in 1:length(comparison_groups)){
    MY.LIST[[i]][MY.LIST[[i]] > 0] = i
  }

  #bind them all together for plotting
  mat_2 = do.call(cbind, MY.LIST) %>%
    as.data.frame(.) %>%
    t(.) %>%
    as.data.frame(.) %>%
    rownames_to_column(., var = "Tumor_Sample_Barcode") %>%
    base::merge(., metadata_df %>%
    rownames_to_column(., "Tumor_Sample_Barcode") %>%
    dplyr::select(Tumor_Sample_Barcode, splitColumnName), by = "Tumor_Sample_Barcode")

  #specify where row breaks should be on heatmap
  sorted_groups <- sort(
    FEATURES$group,
    method = "radix",
    decreasing = FALSE,
    na.last = TRUE
  )
  sorted_groups <- sorted_groups[
    order(as.numeric(sorted_groups))
  ]
  FEATURES = FEATURES %>%
    arrange(match(
      group,
      sorted_groups
    ))
  breaks = 0
  for (this_group in comparison_groups){
    N = (nrow(FEATURES %>%
      dplyr::filter(group == this_group)))
    breaks = c(breaks, N)
  }

  #second, make a vector that will be supplied to ComplexHeatmap
  my_vector = NULL
  for (i in 1:(length(breaks))){
    my_vector = c(my_vector, rep(i - 1, breaks[i]))
  }

  #prepare matrix for stacked barplots on the left
  STACKED = data.frame(matrix(NA, ncol = 1, nrow = nrow(FEATURES)))[-1]
  rownames(STACKED) = FEATURES$Feature
  for (i in 1:length(comparison_groups)) {
  STACKED = cbind(STACKED, mat_2[,c("Tumor_Sample_Barcode", FEATURES$Feature)] %>%
    base::merge(., metadata_df %>%
    rownames_to_column(., "Tumor_Sample_Barcode") %>%
    dplyr::select(Tumor_Sample_Barcode, splitColumnName), by = "Tumor_Sample_Barcode") %>%
    dplyr::arrange(!!sym(splitColumnName)) %>%
    dplyr::filter(base::get(splitColumnName) == comparison_groups[i]) %>%
    dplyr::select(-Tumor_Sample_Barcode, -splitColumnName) %>%
    dplyr::summarise_all(funs(sum)) %>%
    t(.) %>%
    `colnames<-`(comparison_groups[i]) %>%
    as.data.frame(.) %>%
    dplyr::mutate_all(~(./i) / nrow(metadata_df)))
  }

  m = t(apply(STACKED, 1, function(x) x/sum(x)))
  m[is.na(m)] <- 0

  if(prioritize_ordering_on_numeric & ! is.null(numericMetadataColumns)){ # numeric metadata is provided and is prioritized for column sorting
    used_for_ordering_df = t(base::merge(mat_2 %>%
    dplyr::select(-splitColumnName), metadata_df %>%
    rownames_to_column(., "Tumor_Sample_Barcode"), by = "Tumor_Sample_Barcode") %>%
    column_to_rownames(., var = "Tumor_Sample_Barcode") %>%
      dplyr::arrange(desc(!!!syms(numericMetadataColumns)),
        !!!syms(metadataColumns)) %>%
      dplyr::select(FEATURES$Feature))

    this_is_ordered_df = metadata_df[ (order(match(rownames(metadata_df), colnames(used_for_ordering_df)))), ] %>%
      dplyr::arrange(desc(!!!syms(numericMetadataColumns)),
        !!!syms(metadataColumns))
  }else if(! is.null(numericMetadataColumns)){ # numeric metadata is provided, but is not prioritized for column sorting
    used_for_ordering_df = t(base::merge(mat_2 %>%
    dplyr::select(-splitColumnName), metadata_df %>%
    rownames_to_column(., "Tumor_Sample_Barcode"), by = "Tumor_Sample_Barcode") %>%
    column_to_rownames(., var = "Tumor_Sample_Barcode") %>%
      dplyr::arrange(!!!syms(metadataColumns),
        desc(!!!syms(numericMetadataColumns))) %>%
      dplyr::select(FEATURES$Feature))

    this_is_ordered_df = metadata_df[ (order(match(rownames(metadata_df), colnames(used_for_ordering_df)))), ] %>%
      dplyr::arrange(!!!syms(metadataColumns),
        desc(!!!syms(numericMetadataColumns)))
  }else{ # no numeric metadata is proveded to plot
    used_for_ordering_df = t(base::merge(mat_2 %>%
    dplyr::select(-splitColumnName), metadata_df %>%
    rownames_to_column(., "Tumor_Sample_Barcode"), by = "Tumor_Sample_Barcode") %>%
    column_to_rownames(., var = "Tumor_Sample_Barcode") %>%
      dplyr::arrange(!!!syms(metadataColumns)) %>%
      dplyr::select(FEATURES$Feature))

    this_is_ordered_df = metadata_df[ (order(match(rownames(metadata_df), colnames(used_for_ordering_df)))), ] %>%
      dplyr::arrange(!!!syms(metadataColumns))
  }

  # left annotation: stacked feature weights
  ha = rowAnnotation(`feature abundance` = anno_barplot(m, gp = gpar(fill = my_palette[1:length(comparison_groups)+1]),
                                                      bar_width = 1, width = unit(leftStackedWidth, "cm"),
                                                      axis_param = list(side = legend_position, labels_rot = 0)))

  #bottom annotation: tracks indicating metadata
  ha_bottom = HeatmapAnnotation(
    df = this_is_ordered_df %>% dplyr::select(-c(splitColumnName, N_features)),
    col = my_colours,
    simple_anno_size = unit(metadataBarHeight, "mm"),
    gap = unit(0.25 * metadataBarHeight, "mm"),
    annotation_name_gp = gpar(fontsize = metadataBarFontsize),
    annotation_legend_param = list(nrow = legend_row, ncol = legend_col, direction = legend_direction)
  )

  #top annotation: groups of interest to split on
  top_bar_colors = list(my_colours[[splitColumnName]] %>% rev)
  names(top_bar_colors) = splitColumnName
  names(top_bar_colors[[splitColumnName]]) = names(top_bar_colors[[splitColumnName]]) %>% rev()

  ha_top = HeatmapAnnotation(
    group = anno_block(gp = gpar(fill = top_bar_colors[[1]], fontsize = fontSizeGene * 1.5), labels = groupNames),
    "N of features" = anno_barplot(this_is_ordered_df$N_features)
  )

  splendidHM = ComplexHeatmap::Heatmap(used_for_ordering_df,
                                       col = my_palette,
                                       show_column_names = FALSE,
                                       cluster_columns = FALSE,
                                       cluster_rows = FALSE,
                                       row_names_gp = gpar(fontsize = fontSizeGene),
                                       show_heatmap_legend = FALSE,
                                       row_split = my_vector,
                                       row_title = NULL,
                                       left_annotation = ha,
                                       bottom_annotation = ha_bottom,
                                       top_annotation = ha_top,
                                       column_split = dplyr::pull(metadata_df[(order(match(rownames(metadata_df), colnames(used_for_ordering_df)))), ], splitColumnName),
                                       column_title = NULL)

  draw(splendidHM, heatmap_legend_side = legend_position, annotation_legend_side = legend_position)
}
