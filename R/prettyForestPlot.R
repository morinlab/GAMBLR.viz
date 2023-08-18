#' @title Forest Plot.
#'
#' @description Create a forest plot comparing mutation frequencies for a set of genes between two groups.
#'
#' @details This function returns two types of plot (box plot and forest plot), the user can either view them separately or arranged on the same grob.
#' In addition this function also lets the user control the mutation frequencies of the plotted genes.
#' If no `genes` is provided, this function auto-defaults to all genes in the incoming maf.
#' The user can then control the minimum number of mutations requirement for a gene to be included in the plot.
#' This is done with the `min_mutations` parameter.
#' For extended examples on how to use this function, refer to the example inside the function documentation or the vignettes.
#'
#' @param this_maf A maf data frame. Minimum required columns are Hugo_Symbol and Tumor_Sample_Barcode.
#' @param mutmat Optional argument for binary mutation matrix. If not supplied, function will generate this matrix from the file used in argument "maf".
#' @param these_samples_metadata Metadata for the comparisons. Minimum required columns are Tumor_Sample_Barcode and the column assigning each case to one of two groups.
#' @param genes An optional vector of genes to restrict your plot to. If no gene-list is supplied, the function will extract all mutated genes from the incoming maf. See min_mutations parameter for more info.
#' @param min_mutations Optional parameter for when no `genes` is provided. This parameter ensures only genes with n mutations are kept in `genes`. Default value is 1, this means all genes in the incoming maf will be plotted.
#' @param comparison_column Mandatory: the name of the metadata column containing the comparison values.
#' @param rm_na_samples Set to TRUE to remove 0 mutation samples. Default is FALSE.
#' @param comparison_values Optional: If the comparison column contains more than two values or is not a factor, specify a character vector of length two in the order you would like the factor levels to be set, reference group first.
#' @param separate_hotspots Optional: If you would like to treat hotspots separately from other mutations in any gene. Requires that the maf file is annotated with [GAMBLR::annotate_hotspots].
#' @param comparison_name Optional: Specify the legend title if different from the comparison column name.
#' @param custom_colours Optional: Specify a named vector of colours that match the values in the comparison column.
#' @param custom_labels Optional: Specify custom labels for the legend categories. Must be in the same order as comparison_values.
#' @param max_q cut off for q values to be filtered in fish test
#'
#' @return A convenient list containing all the data frames that were created in making the plot, including the mutation matrix. It also produces (and returns) ggplot object with a side-by-side forest plot and bar plot showing mutation incidences across two groups.
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @rawNamespace import(ggpubr, except = "get_legend")
#' @import dplyr cowplot forcats ggplot2 purrr tidyr broom
#' @export
#'
#' @examples
#' #get data
#' my_meta = GAMBLR.data::gambl_metadata
#' my_maf = GAMBLR.data::sample_data$grch37$maf
#'
#' prettyForestPlot(this_maf = my_maf,
#'                  these_samples_metadata = my_meta,
#'                  min_mutations = 50,
#'                  genes = c("ATP6V1B2",
#'                            "EZH2",
#'                            "TNFRSF14",
#'                            "RRAGC"),
#'                  comparison_column = "pathology",
#'                  comparison_values = c("DLBCL",
#'                                        "FL"),
#'                  separate_hotspots = FALSE,
#'                  comparison_name = "FL vs DLBCL")
#'
prettyForestPlot = function(this_maf,
                            mutmat,
                            these_samples_metadata,
                            genes,
                            min_mutations = 1,
                            comparison_column,
                            rm_na_samples = FALSE,
                            comparison_values = FALSE,
                            separate_hotspots = FALSE,
                            comparison_name = FALSE,
                            custom_colours = FALSE,
                            custom_labels = FALSE,
                            max_q = 1){

  #If no comparison_values are specified, derive the comparison_values from the specified comparison_column
  if(comparison_values[1] == FALSE){
    if(class(these_samples_metadata[[comparison_column]]) == "factor"){
      comparison_values = levels(these_samples_metadata[[comparison_column]])
    } else {
      comparison_values = unique(these_samples_metadata[[comparison_column]])
    }
  }
  #Ensure there are only two comparison_values
  {
    if(length(comparison_values) != 2)
      stop(paste0("Your comparison must have two values. \nEither specify comparison_values as a vector of length 2 or subset your metadata so your comparison_column has only two unique values or factor levels."))
  }
  #Subset the metadata to the specified comparison_values and the maf to the remaining sample_ids
  these_samples_metadata = these_samples_metadata[these_samples_metadata[[comparison_column]] %in% comparison_values, ]

  #Ensure the metadata comparison column is a factor with levels matching the input
  these_samples_metadata$comparison = factor(these_samples_metadata[[comparison_column]], levels = comparison_values)

  #read maf into r
  if(!missing(this_maf)){
    #extract gene symbols from maf with minimum N mutations (if no genes list is provided)
    if(missing(genes)){
      genes = this_maf %>%
        dplyr::select(Hugo_Symbol) %>%
        add_count(Hugo_Symbol) %>%
        distinct(Hugo_Symbol, .keep_all = TRUE) %>%
        dplyr::filter(n >= min_mutations) %>%
        pull(Hugo_Symbol)
    }
    maf = this_maf[this_maf$Hugo_Symbol %in% genes, ]
    maf = this_maf[this_maf$Tumor_Sample_Barcode %in% these_samples_metadata$Tumor_Sample_Barcode, ]
  }else{
    stop("Please provide a maf with `this_maf` paraemter...")
  }

  #If separate_hotspots = true, confirm the input maf is hotspot annotated
  if(!missing(mutmat)){
    #add the required columns from the metadata and make the names consistent
    mutmat = left_join(dplyr::select(these_samples_metadata, sample_id, comparison),mutmat) %>%
      dplyr::rename("Tumor_Sample_Barcode"="sample_id")
  }else if(!missing(this_maf)){
    if(separate_hotspots){
      if(!"hot_spot" %in% colnames(this_maf))
        stop("No \"hot_spot\" column in maf file. Annotate your maf file with GAMBLR::annotate_hot_spots() first. ")
      maf$Hugo_Symbol = ifelse(!is.na(maf$hot_spot), paste0(maf$Hugo_Symbol, "_hotspot"), maf$Hugo_Symbol)
    }
    #Convert the maf file to a binary matrix
    mutmat = this_maf %>%
      dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
      full_join(dplyr::select(these_samples_metadata, Tumor_Sample_Barcode, comparison), by = "Tumor_Sample_Barcode") %>%
      distinct() %>%
      dplyr::mutate(is_mutated = 1) %>%
      pivot_wider(names_from = Hugo_Symbol, values_from = is_mutated, values_fill = 0) %>%
      dplyr::mutate(across(where(is.numeric), ~replace_na(., 0)))
    if(rm_na_samples && "NA" %in% colnames(mutmat)){
      mutmat = mutmat %>%
        dplyr::filter(`NA` == 0) #filter out all samples that show no mutations in the selected genes (i.e na_samples = 1).
    }
    if("NA" %in% colnames(mutmat)){
      mutmat = dplyr::select(mutmat, -`NA`) #remove NA column (if there).
    }
  }else{
    message("provide a MAF or mutation matrix")
    return()
  }

  fish_test = mutmat %>%
    pivot_longer(-c(Tumor_Sample_Barcode, comparison), names_to = "gene", values_to = "is_mutated") %>%
    dplyr::mutate(is_mutated = factor(is_mutated, levels = c("1", "0"))) %>%
    group_by(gene) %>%
    dplyr::summarise(table = list(table(is_mutated, comparison))) %>%
    dplyr::mutate(test = map(table, fisher.test), tidy = map(test, broom::tidy)) %>%
    unnest(tidy) %>%
    dplyr::mutate(q.value = p.adjust(p.value, "BH")) %>%
    dplyr::select(-c(test, method, alternative)) %>%
    dplyr::filter(q.value <= max_q) %>%
    dplyr::mutate(gene = fct_reorder(gene, estimate))

  flatten_table <- function(Row){

    mut_n <- Row[2] %>%
    as.data.frame %>%
    `colnames<-`(
      gsub(
        "table.",
        "",
        colnames(.)
      )
    ) %>%
    mutate(
      is_mutated = ifelse(
        is_mutated == 1,
        "mutated",
        "non-mutated"
      ),
      group = paste(
        is_mutated,
        comparison,
        sep = "_"
      )
    ) %>%
    select(group, Freq) %>%
    pivot_wider(
        names_from = group,
        values_from = Freq
    ) %>%
    t %>%
    as.data.frame

    Row <- Row[-2] %>%
      do.call(cbind, .) %>%
      as.data.frame %>%
      t %>% as.data.frame

    rbind(Row, mut_n)

  }


  fish_test <- apply(
    fish_test,
    1,
    flatten_table
  ) %>%
  do.call(cbind, .) %>%
  t %>%
  as.data.frame %>%
  `rownames<-`(NULL) %>%
  mutate_at(c(2:10), as.numeric) %>%
  arrange(estimate)


  point_size = 50 / round(length(fish_test$gene))
  if(point_size < 1){
    point_size = 1
  }
  font_size = 360 / round(length(fish_test$gene))
  if(font_size < 4){
    font_size = 4
  }else if(font_size > 20){
    font_size = 20
  }
  message(paste("FONT:", font_size, "POINT:", point_size, length(fish_test$gene)))
  forest = fish_test %>%
    dplyr::mutate(gene = factor(gene, levels = fish_test$gene)) %>%
    ggplot(aes(x = gene, y = log(estimate))) +
    geom_point(size = point_size, shape = "square") +
    geom_hline(yintercept = 0, lty = 2) +
    coord_flip() +
    geom_errorbar(aes(ymin = log(conf.low), ymax = log(conf.high), width = 0.2)) +
    ylab("ln(Odds Ratio)") +
    xlab("Mutated Genes") +
    cowplot::theme_cowplot() +
    theme(axis.text.y = element_text(size = font_size))

  if(comparison_name == FALSE){
    comparison_name = comparison_column
  }

  if(custom_colours[1] == FALSE){
    if(length(levels(these_samples_metadata$comparison)[levels(these_samples_metadata$comparison) %in% names(get_gambl_colours())]) == 2){
      colours = get_gambl_colours()[levels(these_samples_metadata$comparison)]
    } else {
      colours = get_gambl_colours(classification = "blood")[c("Red", "Blue")]
      names(colours) = levels(these_samples_metadata$comparison)
    }
  } else {
    colours = custom_colours
  }

  if(custom_labels[1] == FALSE){
    labels = levels(these_samples_metadata$comparison)
    names(labels) = levels(these_samples_metadata$comparison)
  } else if(length(custom_labels) != 2) {
    labels = levels(these_samples_metadata$comparison)
    names(labels) = levels(these_samples_metadata$comparison)
    print("Provided custom labels is not a character vector of length 2. Defaulting to comparison factor levels as labels. ")
  } else {
    labels = custom_labels
    names(labels) = comparison_values
  }

  bar = mutmat %>%
    dplyr::select(-Tumor_Sample_Barcode) %>%
    pivot_longer(
		  !comparison,
		  names_to = "gene",
		  values_to = "is_mutated"
	  ) %>%
    group_by(gene, comparison) %>%
    drop_na() %>%
    summarise(percent_mutated = sum(is_mutated) / n() * 100) %>%
    dplyr::filter(gene %in% fish_test$gene) %>%
    dplyr::mutate(gene = factor(gene, levels = fish_test$gene)) %>%
    ggplot(aes(x = gene, y = percent_mutated, fill = comparison)) +
    geom_col(position = "dodge", width = 0.5) +
    xlab("") + ylab("% Mutated") +
    coord_flip() +
    scale_fill_manual(name = comparison_name, values = colours, labels = labels[levels(these_samples_metadata$comparison)]) +
    cowplot::theme_cowplot() +
    theme(axis.text.y = element_blank(), legend.position = "bottom")

  legend = cowplot::get_legend(bar)

  plots = plot_grid(forest, bar +
                      theme(legend.position = "none"), rel_widths = c(1, 0.6), nrow = 1)

  arranged_plot = cowplot::plot_grid(plot_grid(NULL, legend, NULL, nrow = 1), plots, nrow = 2, rel_heights = c(0.1, 1))

  return(list(fisher = fish_test, forest = forest, bar = bar, legend = legend, arranged = arranged_plot, mutmat = mutmat))
}
