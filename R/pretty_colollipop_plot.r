#' @title Pretty CoLollipop Plot.
#'
#' @description Generates a ggplot-compatible figure of 2 [GAMBLR.viz::pretty_lollipop_plot] mirrored, and a 
#' [GAMBLR.viz::prettyForestPlot] displayed below.
#'
#' @details Retrieve maf data of a specific sample or a set of samples for comparison. A gene of interest 
#' can then be visualized with the given maf data files, and comparison commands. Silent mutations can be 
#' visualized setting include_silent to TRUE. 
#'
#' @param maf_df A data frame containing the mutation data.
#' @param metadata Required argument. A data.frame with metadata for the CoLollipop Plot.
#' @param comparison_column Required: the name of the metadata column containing the comparison values.
#' @param comparison_values Optional: If the comparison column contains more than two values or is not a factor, specify a character vector of length two in the order you would like the factor levels to be set, reference group first.
#' @param gene The gene symbol to plot.
#' @param plot_title Optional, the title of the plot. Default is gene.
#' @param ... `pretty_lollipop_plot` arguments, see that function for more info on avaialble parameters.
#' 
#' @return A mirrored lollipop plot.
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#' 
#' #get meta data (DLBCL_Hilton)
#' meta = GAMBLR.data::get_gambl_metadata()
#' metadata = dplyr::filter(meta, cohort %in% "DLBCL_Hilton")
#' maf_df = GAMBLR.data::get_coding_ssm(
#'   these_samples_metadata = metadata
#' )
#' pretty_colollipop_plot_result <- pretty_colollipop_plot(maf_df = maf_df,
#'                                                         metadata = metadata,
#'                                                         comparison_column = "sex",
#'                                                         comparison_values = c("M", "F"),
#'                                                         gene = "IGLL5")
pretty_colollipop_plot <- function(
    maf_df,
    metadata,
    comparison_column,
    comparison_values,
    gene = NULL,
    plot_title,
    ...
) {

    # check for required arguments
    required <- c(
        "maf_df", 
        "metadata", 
        "comparison_column"
        )

    defined <- names(as.list(match.call())[-1])

    if(any(!required %in% defined)) {
        stop("Please provide mutation data and metadata for 2 pretty Oncoplots with specified comparison_column.")
    }

    # If no comparison_values are specified, derive the comparison_values from the specified comparison_column
    if(missing(comparison_values)){
        if(class(metadata[[comparison_column]]) == "factor"){
            comparison_values = levels(metadata[[comparison_column]])
        } else {
            comparison_values = unique(metadata[[comparison_column]])
        }
    }

    # Ensure there are only two comparison_values
    {
        if(length(comparison_values) != 2)
        stop(paste0("Your comparison must have two values. \nEither specify comparison_values as a vector of length 2 or subset your metadata so your comparison_column has only two unique values or factor levels."))
    }

    # Subset the metadata to the specified comparison_values and the maf to the remaining sample_ids
    meta1 <- metadata[metadata[[comparison_column]] %in% comparison_values[1], ]
    meta2 <- metadata[metadata[[comparison_column]] %in% comparison_values[2], ]

    # Subset maf to only samples in the comparison values
    ssm1 <- maf_df %>%
        dplyr::filter(
            Tumor_Sample_Barcode %in% meta1$Tumor_Sample_Barcode
        )

    ssm2 <- maf_df %>%
        dplyr::filter(
            Tumor_Sample_Barcode %in% meta2$Tumor_Sample_Barcode
        )

    # Ensure dimensions are greater than zero otherwise no variants are returned.
    dim1 <- dim(ssm1)[1]  
    dim2 <- dim(ssm2)[1]

    if(dim1 == 0 | dim2 == 0) {
        stop(paste0("Ensure all variants in metadata are accounted for. \nEither ensure all variants are loaded into metadata as one dataframe, or subset your metadata so that all variants are included."))
    } 

    # Arguments to pass into pretty_lollipop_plot
    lollipopplot_args <- list(...)

    # Get gene_counts data for the first plot (without generating the plot)
    lp1_gene_counts <- do.call(pretty_lollipop_plot, c(
        list(
            maf_df = ssm1,
            gene = gene,
            plotarg = FALSE
        ), 
        lollipopplot_args
    ))
    lp1_gene_counts_data <- as.data.frame(lp1_gene_counts)

    # Get gene_counts data for the second plot
    lp2_gene_counts <- do.call(pretty_lollipop_plot, c(
        list(
            maf_df = ssm2,
            gene = gene,
            plotarg = FALSE
        ), 
        lollipopplot_args
    ))
    lp2_gene_counts_data <- as.data.frame(lp2_gene_counts)
    lp2_gene_counts_data <- lp2_gene_counts_data %>%
        mutate(mutation_count = mutation_count * -1)

    # Combine data for both plots
    lp1_gene_counts_data <- lp1_gene_counts_data %>% 
        mutate(source = "Sample 1")
    lp2_gene_counts_data <- lp2_gene_counts_data %>% 
        mutate(source = "Sample 2")

    combined_gene_counts <- rbind(
        lp1_gene_counts_data, 
        lp2_gene_counts_data 
    )

    # Setting Somatic Mutation Statistic
    meta1_counter <- length(unique(meta1$Tumor_Sample_Barcode))
    meta2_counter <- length(unique(meta2$Tumor_Sample_Barcode))

    # Pass combined_gene_counts to pretty_lollipop_plot for plotting
    colollipop_plot <- pretty_lollipop_plot(
        maf_df = maf_df, 
        gene = gene, 
        plot_title = plot_title,
        include_silent = FALSE,
        plotarg = TRUE,
        mirrorarg = TRUE,
        combined_gene_counts = combined_gene_counts,
        meta1_counter = meta1_counter,
        meta2_counter = meta2_counter,
        Sample1 = comparison_values[1],
        Sample2 = comparison_values[2]
    )

    forest_plot <- prettyForestPlot(
        maf = maf_df,
        metadata = metadata,
        genes = gene,
        comparison_column = comparison_column,
        comparison_values = comparison_values,
        separate_hotspots = FALSE,
        comparison_name = paste0(
            comparison_values[1],
            " vs. ",
            comparison_values[2]
        ),
        mirrorarg = TRUE
    )

    # arrange colollipop plot and forest plot together one over the other
    plot <- ggarrange(
        colollipop_plot,              
        ggarrange(
            forest_plot$forest,        
            forest_plot$bar, 
            widths = c(1, 0.6),
            common.legend = TRUE,
            align = "h"),          
        ncol = 1, 
        nrow = 2,
        heights = c(2, 1)            
    )
    return(plot)
}