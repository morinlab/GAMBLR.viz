#' @title Side-by-side Oncoplots
#'
#' @description `prettyCoOncoplot` returns figure of an arbitrary number of
#' [GAMBLR.viz::prettyOncoplot] images placed side-by-side with the same gene order.
#'
#' @details This function will generate a graphic displaying 2 or more oncoplots
#' side-by-side. Optionally user can annotate each oncoplot with it's own title
#' that will be displayed at the top via the labels argument.
#' All the arguments recognized by
#' [GAMBLR.viz::prettyOncoplot] are supported and can be specified when calling
#' this function. For both oncoplots the same specified parameters will be
#' applied (e.g. genes to display, split columns, font size, top annotation
#' etc). When the provided argument is not recognized by
#' [GAMBLR.viz::prettyOncoplot], it will be ignored. Some features of the oncoplots
#' will be dropped from all but the right-most oncoplot (e.g. gene names).
#' The order of the oncoplots is defined (left-to-right) by either the order of
#' the factors in the `comparison_column` when it's a factor, the order of
#' appearance of the values in the `comparison_column`, or by the order of
#' elements in the `comparison_values` when the comparison column contains more
#' than 2 levels. The left oncoplot will dictate the order of genes to be
#' displayed on the right oncoplot. If you want a specific order of oncoplots on
#' the left and right, please ensure the argument `comparison_column` is a
#' factor with first level being the group you want to be plotted on the left
#' side.
#' For developers: new arguments added to [GAMBLR.viz::prettyOncoplot] in the
#' future are expected to be out-of-the-box compatible with this function and
#' would not need code modifications in this function.
#'
#' @param maf Required argument. A data frame containing the mutations you want
#'      to plot on both oncoplots.
#' @param metadata Required argument. A data.frame with metadata for both
#'      oncoplots.
#' @param comparison_values Optional: If the comparison column contains more
#'      than two values or is not a factor, specify a character vector of length
#'      two in the order you would like the factor levels to be set, reference
#'      group first.
#' @param comparison_column Required: the name of the metadata column containing
#'      the comparison values.
#' @param labels Optional: A named vector with the value representing the name
#' you want displayed on the top of each oncoplot. The names of the vector
#'      should match the values in the `comparison_column` of the metadata.
#' @param keepGeneOrder Optional: If TRUE, the genes will be displayed
#' in the order they are found in the `genes` argument.
#' @param individualOncoplotArgs Optional: A named list of arguments to be
#' passed to [GAMBLR.viz::prettyOncoplot] for customizing each oncoplot.
#' The names of the list must match at least one of the values in the 
#' `comparison_column` of the metadata.
#' @param returnEverything Optional: If TRUE, the function will return a list of
#'     the inputs used to generate the oncoplots.
#'     Useful for subsequent functions that require the same inputs.
#' @param ... `prettyOncoplot` arguments, see that function for more info on
#'      avaialble parameters.
#'
#' @return A complexheatmap object with 2 oncoplots side-by-side.
#'
#' @import ComplexHeatmap dplyr
#' @export
#'
#' @examples
#' cat("Running example for function: prettyCoOncoplot\n")
#' library(GAMBLR.open)
#' #get data for plotting
#' meta <- get_gambl_metadata()
#' meta <- meta %>%
#'     dplyr::filter(
#'         study == "FL_Dreval",
#'         pathology %in% c("DLBCL", "FL")
#'     )
#' ssm <- get_coding_ssm(
#'     these_samples_metadata = meta
#' )
#'
#' #build plot
#' prettyCoOncoplot(
#'     maf = ssm,
#'     metadata = meta,
#'     comparison_column = "pathology",
#'     genes=dplyr::filter(lymphoma_genes,
#'                         FL_Tier==1 | DLBCL_Tier==1) %>% 
#'                         dplyr::pull(Gene),
#'     metadataColumns = c(
#'         "pathology",
#'         "lymphgen",
#'         "pairing_status"
#'     ),
#'     metadataBarHeight = 10,
#'     fontSizeGene = 12,
#'     metadataBarFontsize = 10,
#'     legend_row = 2,
#'     label1 = "DLBCL",
#'     label2 = "FL",
#'     simplify_annotation =TRUE,
#'     minMutationPercent = 5
#' )
#'
prettyCoOncoplot = function(maf,
                            metadata,
                            comparison_column,
                            comparison_values,
                            labels,
                            genes,
                            keepGeneOrder = FALSE,
                            returnEverything = FALSE,
                            individualOncoplotArgs,
                            ...){

    # check for required arguments
    required = c("maf", "metadata", "comparison_column")

    defined = names(as.list(match.call())[-1])
    all_args <- as.list(match.call(expand.dots = TRUE)[-1])
    #all_args <- mget(names(formals()),             # the formal names
    #               sys.frame(sys.nframe()),      # the current call frame
    #               inherits = FALSE)             # don’t look up the stack
    if (any(!required %in% defined)) {
      stop("Please provide mutation data and metadata for 2 pretty Oncoplots with specified comparison_column.")
    }

    #If no comparison_values are specified, derive the comparison_values from the specified comparison_column
    if(missing(comparison_values)){
      if("factor" %in% class(metadata[[comparison_column]])){
        comparison_values = levels(metadata[[comparison_column]])
      } else {
        comparison_values = unique(metadata[[comparison_column]])
      }
    }
    #warn the user about how many comparison values are being used
    if(length(comparison_values) > 1){
      message(paste("Number of comparison values is:",
                    length(comparison_values),
                    "Proceeding with all of them."))
    }else{
        stop("comparison_values has only one value.")
    }
    #Subset the metadata to the specified comparison_values and the maf to the remaining sample_ids
 
    #list of metadata split by comparison_column
    meta_split = split(metadata, metadata[[comparison_column]])
    # Subset maf to only samples in the comparison values
    maf = left_join(maf,select(metadata, sample_id, all_of(comparison_column)),
                    by=c("Tumor_Sample_Barcode"="sample_id"))
    maf_split = split(maf, maf[[comparison_column]])


    # Arguments to pass into prettyOncoplot
    oncoplot_args = list(...)
    # Discard any arguments not supported by prettyOncoplot
    oncoplot_args = oncoplot_args[names(oncoplot_args) %in% intersect(names(oncoplot_args),
                                                                      formalArgs(prettyOncoplot))]
    # Ignore genes argument for the right oncoplot so it matches the left side
    message(
        "Setting up the order of genes to match between the two oncoplots..."
    )
    oncoplot_args_rerun <- oncoplot_args[
        !names(oncoplot_args) %in% c(
            "genes", "keepGeneOrder", "minMutationPercent"
        )
    ]
    if(!keepGeneOrder && !missing(genes)){
        op1 <- do.call(
            prettyOncoplot,
            c(
                list(
                    maf_split[[comparison_values[1]]],
                    these_samples_metadata = meta_split[[comparison_values[1]]],
                    return_inputs = TRUE,
                    genes = genes
                ),
            oncoplot_args
            )
        )
        gene_order = op1$gene_order
    }else{
        gene_order = genes
    }
    ht_list = list()
    i = 1
    suppressMessages(
        suppressWarnings({
            for(group_name in comparison_values){
                if(i == length(comparison_values)){
                    hide_gene_names = FALSE
                    hide_annotation_name = FALSE
                    show_any_legend = TRUE
                } else {
                    hide_gene_names = TRUE
                    hide_annotation_name = TRUE
                    show_any_legend = FALSE
                }
                if(!missing(labels)){
                    label = labels[group_name]
                }else{
                    label = group_name
                }
                oncoArgs = c(list(
                                maf_df = maf_split[[group_name]],
                                these_samples_metadata = meta_split[[group_name]],
                                genes = gene_order,
                                keepGeneOrder = TRUE,
                                
                                minMutationPercent = 0,
                                return_inputs = TRUE,
                                groupNames = label,
                                hide_gene_names = hide_gene_names,
                                hide_annotation_name = hide_annotation_name,
                                show_any_legend = show_any_legend
                                ),
                            oncoplot_args_rerun)
                            #override arguments specified in individualOncoplotArgs if applicable
                if(!missing(individualOncoplotArgs)){
                    if(group_name %in% names(individualOncoplotArgs)){
                        #print(group_name)
                        for(argname in names(individualOncoplotArgs[[group_name]])){
                            #if the argument is not in the oncoplot_args, add it
                            if(!argname %in% names(oncoArgs)){
                                message(paste("Adding", argname, "for", group_name))
                                oncoArgs[[argname]] = individualOncoplotArgs[[group_name]][[argname]]
                            }else{
                                #if the argument is in the oncoplot_args, override it
                                message(paste("Overriding", argname, "for", group_name))
                                oncoArgs[[argname]] = individualOncoplotArgs[[group_name]][[argname]]
                            }
                        }
                    }
                }
                if(!any(grepl("sort", names(oncoArgs)))){
                    oncoArgs[["sortByGenes"]] = gene_order
                }
                i = i + 1
                op <- do.call(
                    prettyOncoplot,
                    oncoArgs
                )
                ht_list[[group_name]] <- op$Heatmap
            }
        
            
            heatmap_list <- Reduce(`+`, ht_list) 
            
                
            p <- draw(
                heatmap_list,
                merge_legends = TRUE,
                heatmap_legend_side = "bottom",
                annotation_legend_side = "bottom"
            )
        })
    )
    if(returnEverything){
        return(
            list(Combined = p,
                args = all_args,
                Heatmap = heatmap_list,
                gene_order = gene_order)
            )
    }
    return(p)
  }
