#' @title Side-by-side Oncoplots
#'
#' @description `prettyCoOncoplot` returns ggplot-compatible figure of 2 [GAMBLR::prettyOncoplot] side-by-side.
#'
#' @details This function will generate a graphic displaying 2 oncoplots side-by-side. Optionally user can
#' annotate each oncoplot with it's own title that will be displayed at the top. All the arguments
#' recognized by [GAMBLR::prettyOncoplot] are supported and can be specified when calling this function.
#' For both oncoplots the same specified parameters will be applied (e.g. genes to display, split columns,
#' font size, top annotation etc). If the provided argument is not recognized by prettyOncoplot,
#' it will be discarded. If you want a specific order of oncoplots on the left and right, please
#' ensure the argument `comparison_column` is a factor with first level being the group
#' you want to be plotted on the left side. For developers: new arguments added to prettyOncoplot in the future
#' are expected to be out-of-the-box compatible with this function nd would not need code modifications.
#'
#' @param maf Required argument. A maftools object containing the mutations you want to plot on both oncoplots.
#' @param metadata Required argument. A data.frame with metadata for both oncoplots.
#' @param comparison_values Optional: If the comparison column contains more than two values or is not a factor, specify a character vector of length two in the order you would like the factor levels to be set, reference group first.
#' @param comparison_column Required: the name of the metadata column containing the comparison values.
#' @param label1 Optional argument. Label to be shown as a title for the oncoplot #1.
#' @param label2 Optional argument. Label to be shown as a title for the oncoplot #2.
#' @param ... `prettyOncoplot` arguments, see that function for more info on avaialble parameters.
#'
#' @return A ggplot object with 2 oncoplots side-by-side.
#'
#' @rawNamespace import(ggpubr, except = "get_legend")
#' @import ComplexHeatmap dplyr maftools
#' @export
#'
#' @examples
#' #get data for plotting
#' ssm = get_coding_ssm(limit_cohort = c("BL_Adult", "BL_Pediatric"), seq_type = "genome")
#' ssm = maftools::read.maf(ssm)
#' meta = get_gambl_metadata()
#' meta = dplyr::filter(meta, cohort %in% c("BL_Adult", "BL_Pediatric"))
#'
#' #build plot
#' prettyCoOncoplot(maf = ssm,
#'                  metadata = meta,
#'                  comparison_column = "cohort",
#'                  include_noncoding = NULL,
#'                  minMutationPercent = 0,
#'                  genes = c("MYC",
#'                            "TET2",
#'                            "TP53",
#'                            "DDX3X",
#'                            "ID3"),
#'                  metadataColumns = c("pathology",
#'                                      "EBV_status_inf",
#'                                      "pairing_status",
#'                                      "cohort"),
#'                  splitColumnName = "EBV_status_inf",
#'                  metadataBarHeight = 10,
#'                  fontSizeGene = 12,
#'                  metadataBarFontsize = 10,
#'                  legend_row = 2,
#'                  label1 = "Adult",
#'                  label2 = "Pediatric")
#'
prettyCoOncoplot = function(maf,
                            metadata,
                            comparison_column,
                            comparison_values,
                            label1,
                            label2,
                            ...){

    # check for required arguments
    required = c("maf", "metadata", "comparison_column")

    defined = names(as.list(match.call())[-1])

    if (any(!required %in% defined)) {
      stop("Please provide mutation data and metadata for 2 pretty Oncoplots with specified comparison_column.")
    }

    #If no comparison_values are specified, derive the comparison_values from the specified comparison_column
    if(missing(comparison_values)){
      if(class(metadata[[comparison_column]]) == "factor"){
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

    #Subset the metadata to the specified comparison_values and the maf to the remaining sample_ids
    meta1 = metadata[metadata[[comparison_column]] %in% comparison_values[1], ]
    meta2 = metadata[metadata[[comparison_column]] %in% comparison_values[2], ]

    # Subset maf to only samples in the comparison values
    ssm1 = maftools::subsetMaf(maf,tsb=pull(meta1, Tumor_Sample_Barcode))
    ssm2 = maftools::subsetMaf(maf,tsb=pull(meta2, Tumor_Sample_Barcode))

    # Arguments to pass into prettyOncoplot
    oncoplot_args = list(...)
    # Discard any arguments not supported by prettyOncoplot
    oncoplot_args = oncoplot_args[names(oncoplot_args) %in% intersect(names(oncoplot_args),
                                                                      formalArgs(prettyOncoplot))]
    # Build oncoplot No1
    op1 = do.call(prettyOncoplot, c(
      list(
        maftools_obj = ssm1,
        these_samples_metadata = meta1
      ),
      oncoplot_args
    ))
    # convert it to ggplot object
    op1 = grid.grabExpr(draw(op1), width = 10, height = 17)
    # if user provided annotation label, place it as a name for oncoplot No1
    if (!missing(label1)) {
      op1 = annotate_figure(op1,
                            top = text_grob(label1,
                                            face = "bold",
                                            size = 20))
    }
    # Build oncoplot No2
    op2 = do.call(prettyOncoplot, c(
      list(
        maftools_obj = ssm2,
        these_samples_metadata = meta2
      ),
      oncoplot_args
    ))
    # convert it to ggplot object
    op2 = grid.grabExpr(draw(op2), width = 10, height = 17)
    # if user provided annotation label, place it as a name for oncoplot No2
    if (!missing(label2)) {
      op2 = annotate_figure(op2,
                            top = text_grob(label2,
                                            face = "bold",
                                            size = 20))
    }
    # arrange 2 oncoplots together side by side
    p = ggarrange(op1,
                  op2,
                  ncol = 2,
                  nrow = 1)

    return(p)
  }
