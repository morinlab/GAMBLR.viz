#' @title Forest Plot.
#'
#' @description Create a forest plot comparing mutation frequencies
#' for a set of genes between two groups.
#'
#' @details This function returns two types of plot (box plot and
#' forest plot), the user can either view them separately or arranged
#' on the same grob.
#' In addition this function also lets the user control the mutation
#' frequencies of the plotted genes.
#' If no `genes` is provided, this function auto-defaults to all genes
#' in the incoming maf.
#' The user can then control the minimum number of mutations requirement
#' for a gene to be included in the plot.
#' This is done with the `min_mutations` parameter.
#' For extended examples on how to use this function, refer to the example
#' inside the function documentation or the vignettes.
#'
#' @param maf A maf data frame. Minimum required columns are Hugo_Symbol
#' and Tumor_Sample_Barcode.
#' @param mutmat Optional argument for binary mutation matrix. If not
#' supplied, function will generate this matrix from the file used in
#' argument "maf".
#' @param metadata Metadata for the comparisons. Minimum required
#' columns are Tumor_Sample_Barcode and the column assigning each case
#' to one of two groups.
#' @param genes An optional vector of genes to restrict your plot to.
#' If no gene-list is supplied, the function will extract all mutated
#' genes from the incoming maf. See min_mutations parameter for more info.
#' @param keepGeneOrder Set to TRUE if you want to preserve the gene
#' order specified.
#' @param min_mutations Optional parameter for when no `genes` is provided.
#' This parameter ensures only genes with n mutations are kept in `genes`.
#' Default value is 1, this means all genes in the incoming maf will be plotted.
#' @param comparison_column Mandatory: the name of the metadata column
#' containing the comparison values.
#' @param rm_na_samples Set to TRUE to remove 0 mutation samples.
#' Default is FALSE.
#' @param comparison_values Optional: If the comparison column contains
#' more than two values or is not a factor, specify a character vector
#' of length two in the order you would like the factor levels to be set,
#' reference group first.
#' @param separate_hotspots Deprecated. Use `maf_column = "hot_spot"` instead.
#' @param maf_column The MAF column to use as the feature identifier for each
#' row in the mutation matrix. Default is `"Hugo_Symbol"` (one test per gene).
#' Set to `"mutation_alias"` or `"driver_alias"` to use sub-gene categories
#' produced by [GAMBLR.utils::annotate_curated_drivers] (e.g. `TP53_trunc`,
#' `TP53_R175`, `MYD88_L265P`). Set to `"hot_spot"` for backward-compatible
#' hotspot separation: TRUE rows become `<Hugo_Symbol>HOTSPOT`, FALSE rows
#' keep `Hugo_Symbol`. The `genes` parameter always filters by `Hugo_Symbol`
#' first, regardless of `maf_column`.
#' @param comparison_name Optional: Specify the legend title if different
#' from the comparison column name.
#' @param custom_colours Optional: Specify a named vector of colours that
#' match the values in the comparison column.
#' @param custom_labels Optional: Specify custom labels for the legend
#' categories. Must be in the same order as comparison_values.
#' @param max_q cut off for q values to be filtered in fish test
#' @param base_size Numeric value to specify font size of the gene labels on the plot. Default is 10.
#' @param bar_label Optional: Specify label for forest plot. Default label is "\% Mutated".
#' @param show_legend Set to FALSE to hide legend from final figure. Default is TRUE.
#'
#' @return A convenient list containing all the data frames that were
#' created in making the plot, including the mutation matrix. It also
#' produces (and returns) ggplot object with a side-by-side forest plot
#' and bar plot showing mutation incidences across two groups.
#'
#' @rawNamespace import(ggpubr, except = "get_legend")
#' @import dplyr ggplot2 purrr tidyr GAMBLR.helpers
#' @export
#'
#' @examples
#' cat("Running example for function: prettyForestPlot\n")
#'
#' # Compare FL vs DLBCL using fine-grained driver/hotspot categories from
#' # annotate_curated_drivers() instead of plain gene-level mutation status.
#' # This reveals whether specific mutation subtypes (e.g. EZH2 hotspot vs
#' # truncating vs other missense) drive the association rather than the
#' # gene overall.
#'
#' suppressPackageStartupMessages(library(GAMBLR.open))
#' suppressPackageStartupMessages(library(GAMBLR.utils))
#' suppressWarnings(
#'   suppressMessages({
#'
#' meta = get_gambl_metadata()
#' fl_dlbcl_meta = dplyr::filter(meta, pathology %in% c("FL", "DLBCL")) %>%
#'   check_and_clean_metadata(duplicate_action = "keep_first")
#'
#' # Get all coding SSMs (genome + capture) for these samples
#' fl_dlbcl_maf = get_all_coding_ssm(fl_dlbcl_meta)
#'
#' # Annotate driver/hotspot categories for selected genes
#' fl_dlbcl_maf = annotate_curated_drivers(
#'   maf_data = fl_dlbcl_maf,
#'   genes_of_interest = c("EZH2", "CREBBP", "TP53",
#'                         "MYD88", "NOTCH1", "NOTCH2")
#' )
#'
#' # mutation_alias gives the finest resolution: each named region or hotspot
#' # codon becomes its own feature (e.g. EZH2_Y641, CREBBP_KAT, TP53_trunc,
#' # TP53_R175). One Fisher test is run per category.
#' alias_plots = prettyForestPlot(
#'   maf = fl_dlbcl_maf,
#'   metadata = fl_dlbcl_meta,
#'   genes = c("EZH2", "CREBBP", "TP53", "MYD88", "NOTCH1", "NOTCH2"),
#'   maf_column = "mutation_alias",
#'   comparison_column = "pathology",
#'   comparison_values = c("DLBCL", "FL"),
#'   comparison_name = "FL vs DLBCL"
#' )
#' alias_plots$arranged
#'
#' # driver_alias collapses to just two categories per gene: GENE_driver (any
#' # mutation matching a curated driver region, regardless of specific position)
#' # vs GENE_other (everything else). This is coarser than mutation_alias but
#' # useful when you only care whether a driver-class mutation is enriched,
#' # not which specific hotspot or region it falls in.
#' driver_plots = prettyForestPlot(
#'   maf = fl_dlbcl_maf,
#'   metadata = fl_dlbcl_meta,
#'   genes = c("EZH2", "CREBBP", "TP53", "MYD88", "NOTCH1", "NOTCH2"),
#'   maf_column = "driver_alias",
#'   comparison_column = "pathology",
#'   comparison_values = c("DLBCL", "FL"),
#'   comparison_name = "FL vs DLBCL"
#' )
#' driver_plots$arranged
#'
#' # hotspot_alias captures only the strictest hotspot class: named recurrent
#' # codons (e.g. EZH2_hotspot, MYD88_hotspot). Truncating mutations and
#' # general LOF regions do not produce a hotspot_alias entry, so genes like
#' # TP53 whose driver biology is dominated by truncations will have few or no
#' # rows here. Use mutation_alias if you want those included.
#' hotspot_plots = prettyForestPlot(
#'   maf = fl_dlbcl_maf,
#'   metadata = fl_dlbcl_meta,
#'   genes = c("EZH2", "CREBBP", "TP53", "MYD88", "NOTCH1", "NOTCH2"),
#'   maf_column = "hotspot_alias",
#'   comparison_column = "pathology",
#'   comparison_values = c("DLBCL", "FL"),
#'   comparison_name = "FL vs DLBCL"
#' )
#' hotspot_plots$arranged
#'
#' }))
#'
#' suppressPackageStartupMessages(library(GAMBLR.open))
#' suppressWarnings(
#'   suppressMessages({
#'
#' metadata = get_gambl_metadata()
#' this_meta = dplyr::filter(metadata, pairing_status == "matched")
#' this_meta = dplyr::filter(this_meta, pathology %in% c("FL", "DLBCL")) %>%
#'             check_and_clean_metadata(.,duplicate_action="keep_first")
#'
#' maf = get_coding_ssm(these_samples_metadata = this_meta)
#'
#' plots = prettyForestPlot(maf = maf,
#'                  metadata = this_meta,
#'                  genes = c("ATP6V1B2",
#'                            "EZH2",
#'                            "TNFRSF14",
#'                            "RRAGC"),
#'                  comparison_column = "pathology",
#'                  comparison_values = c("DLBCL",
#'                                        "FL"),
#'                  separate_hotspots = FALSE,
#'                  comparison_name = "FL vs DLBCL")
#' plots$arranged
#'
#' }))
prettyForestPlot = function(maf,
                            mutmat,
                            metadata,
                            genes,
                            keepGeneOrder = FALSE,
                            min_mutations = 1,
                            comparison_column,
                            rm_na_samples = FALSE,
                            comparison_values = FALSE,
                            separate_hotspots = FALSE,
                            maf_column = "Hugo_Symbol",
                            comparison_name = FALSE,
                            custom_colours = FALSE,
                            custom_labels = FALSE,
                            max_q = 1,
                            base_size = 10,
                            bar_label = "% Mutated",
                            show_legend = TRUE
                            ){

  #If no comparison_values are specified, derive the comparison_values
  # from the specified comparison_column
  if(comparison_values[1] == FALSE){
    if("factor" %in% class(metadata[[comparison_column]])){
      comparison_values = levels(metadata[[comparison_column]])
    } else {
      comparison_values = unique(metadata[[comparison_column]])
    }
  }
  #Ensure there are only two comparison_values
  {
    if(length(comparison_values) != 2)
      stop(paste0("Your comparison must have two values. \nEither specify comparison_values as a vector of length 2 or subset your metadata so your comparison_column has only two unique values or factor levels."))
  }

  # Need to properly sort the genes later
  if(missing(genes)){
    user_provided_gene_order = FALSE
  }else{
    user_provided_gene_order = TRUE
  }

  # Subset the metadata to the specified comparison_values
  # and the maf to the remaining sample_ids
  metadata = metadata[metadata[[comparison_column]] %in% comparison_values, ]

  #Ensure the metadata comparison column is a factor with levels matching the input
  metadata$comparison = factor(metadata[[comparison_column]],
                               levels = comparison_values)

  #read maf into r
  if(!missing(maf)){
    maf <- GAMBLR.utils::strip_genomic_classes(maf)
    #extract gene symbols from maf with minimum N mutations (if no genes list is provided)
    if(missing(genes)){
      genes = maf %>%
        dplyr::select(Hugo_Symbol) %>%
        add_count(Hugo_Symbol) %>%
        distinct(Hugo_Symbol, .keep_all = TRUE) %>%
        dplyr::filter(n >= min_mutations) %>%
        pull(Hugo_Symbol)
    }
    maf = maf[maf$Hugo_Symbol %in% genes, ]
    maf = maf[maf$Tumor_Sample_Barcode %in% metadata$Tumor_Sample_Barcode, ]
  }

  if(!missing(mutmat)){
    #add the required columns from the metadata and make the names consistent
    mutmat = left_join(dplyr::select(metadata, sample_id, comparison),mutmat) %>%
      dplyr::rename("Tumor_Sample_Barcode"="sample_id")
  }else if(!missing(maf)){
    # backward compat: separate_hotspots = TRUE is equivalent to maf_column = "hot_spot"
    if(separate_hotspots && maf_column == "Hugo_Symbol"){
      .Deprecated(msg = "separate_hotspots is deprecated; use maf_column = \"hot_spot\" instead.")
      maf_column = "hot_spot"
    }

    # derive the feature column used to name columns in the mutation matrix
    if(maf_column == "hot_spot"){
      if(!"hot_spot" %in% colnames(maf))
        stop("No \"hot_spot\" column in maf. Run annotate_curated_drivers() first.")
      maf$.feature = ifelse(
        maf$hot_spot %in% c(TRUE, "TRUE"),
        paste0(maf$Hugo_Symbol, "HOTSPOT"),
        maf$Hugo_Symbol
      )
    } else if(maf_column == "Hugo_Symbol"){
      maf$.feature = maf$Hugo_Symbol
    } else {
      if(!maf_column %in% colnames(maf))
        stop(paste0("Column \"", maf_column, "\" not found in maf. Run annotate_curated_drivers() first."))
      maf$.feature = maf[[maf_column]]
    }

    #Convert the maf file to a binary matrix
    mutmat = maf %>%
      dplyr::select(.feature, Tumor_Sample_Barcode) %>%
      full_join(dplyr::select(metadata, Tumor_Sample_Barcode, comparison), by = "Tumor_Sample_Barcode") %>%
      distinct() %>%
      dplyr::mutate(is_mutated = 1) %>%
      pivot_wider(names_from = .feature, values_from = is_mutated, values_fill = 0) %>%
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
    dplyr::mutate(
        test = map(table, fisher.test),
        estimate = map_dbl(test, ~ .x$estimate),
        p.value = map_dbl(test, ~ .x$p.value),
        conf.low = map_dbl(test, ~ .x$conf.int[1]),
        conf.high = map_dbl(test, ~ .x$conf.int[2]),
        method = map_chr(test, ~ .x$method),
        alternative = map_chr(test, ~ .x$alternative)
    ) %>%
    dplyr::mutate(q.value = p.adjust(p.value, "BH")) %>%
    dplyr::select(-c(test, method, alternative)) %>%
    dplyr::filter(q.value <= max_q)

  if(nrow(fish_test) < 1){
    stop(
        "There are no differences between the comparison groups that pass the max_q value!"
    )
  }

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
  mutate_at(c(2:10), as.numeric)


  if(keepGeneOrder & user_provided_gene_order){
    # account for possible gene loss due to the max_q cutoff
    message(
        "Supporting the custom-provided gene order instead of ordering on OR"
    )
    genes <- genes[genes %in% fish_test$gene]
    # now order
    # need to actually reverse the order here because there is coord_flip in ggplot
    fish_test <- fish_test %>%
        arrange(match(gene, rev(genes)))
  }else{
    fish_test <- fish_test %>%
        arrange(estimate)
  }


  point_size = 50 / round(length(fish_test$gene))
  error_width = 0.2
  geom_col_width = 0.5
  if(point_size < 1){
    point_size = 1
  }
  if(point_size >= 5){
    point_size = 5
    error_width = 0.1
    geom_col_width = 0.3
  }

  font_size <- base_size
  message(paste("FONT:", font_size, "POINT:", point_size, length(fish_test$gene)))
  forest = fish_test %>%
    dplyr::mutate(gene = factor(gene, levels = fish_test$gene)) %>%
    ggplot(aes(x = gene, y = log(estimate))) +
    geom_point(size = point_size, shape = "square") +
    geom_hline(yintercept = 0, lty = 2) +
    coord_flip() +
    geom_errorbar(aes(ymin = log(conf.low), ymax = log(conf.high), width = error_width)) +
    ylab("ln(Odds Ratio)") +
    xlab("Mutated Genes\n") +
    theme_Morons(base_size = base_size) +
    theme(axis.text.y = element_text(size = font_size))

  if(comparison_name == FALSE){
    comparison_name = comparison_column
  }

  if(custom_colours[1] == FALSE){
    if(length(levels(metadata$comparison)[levels(metadata$comparison) %in% names(GAMBLR.helpers::get_gambl_colours())]) == 2){
      colours = GAMBLR.helpers::get_gambl_colours()[levels(metadata$comparison)]
    } else {
      colours = GAMBLR.helpers::get_gambl_colours(classification = "blood")[c("Red", "Blue")]
      names(colours) = levels(metadata$comparison)
    }
  } else {
    colours = custom_colours
  }

  if(custom_labels[1] == FALSE){
    labels = levels(metadata$comparison)
    names(labels) = levels(metadata$comparison)
  } else if(length(custom_labels) != 2) {
    labels = levels(metadata$comparison)
    names(labels) = levels(metadata$comparison)
    print("Provided custom labels is not a character vector of length 2. Defaulting to comparison factor levels as labels. ")
  } else {
    labels = custom_labels
    names(labels) = comparison_values
  }

  if (show_legend == FALSE) {
    legend_pos = "none"
  } else {
    legend_pos = "bottom"
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
    geom_col(position = "dodge", width = geom_col_width) +
    xlab("") + ylab(bar_label) +
    coord_flip() +
    scale_fill_manual(name = comparison_name, values = colours, labels = labels[levels(metadata$comparison)]) +
    theme_Morons(base_size = base_size) +
    theme(axis.text.y = element_blank(), legend.position = legend_pos)

  arranged_plot = ggarrange(
    forest,
    bar,
    widths = c(1, 0.6),
    common.legend = show_legend,
    align = "h"
  )

  return(list(fisher = fish_test, forest = forest, bar = bar, arranged = arranged_plot, mutmat = mutmat))
}
