#' @title Pretty CoLollipop Plot.
#'
#' @description Generates a ggplot-compatible figure of 2 [GAMBLR.viz::pretty_lollipop_plot] mirrored.
#'
#' @details Retrieve maf data of a specific sample or a set of samples for comparison. A gene of interest 
#' can then be visualized with the given maf data files, and comparison commands. Silent mutations can be 
#' visualized setting include_silent to TRUE. 
#'
#' @param maf_df A data frame containing the mutation data.
#' @param these_samples_metadata Required argument. A data.frame with metadata for the CoLollipop Plot.
#' @param comparison_column Required: the name of the metadata column containing the comparison values.
#' @param comparison_values Optional: If the comparison column contains more than two values or is not a factor, specify a character vector of length two in the order you would like the factor levels to be set, reference group first.
#' @param compare_distributions Logical parameter indicating whether to compare the distribution of variants along the length of the gene and its domains between the groups. Default is FALSE. 
#' @param gene The gene symbol to plot.
#' @param label_threshold Threshold for labels to appear on plot. Default 5. 
#' @param plot_title Optional, the title of the plot. Default is `{gene}: {comparison_values[1]} vs. {comparison_values[2]}`.
#' @param refseq_id Insert a specific NM_xxx value of interest.
#' @param forestarg Logical parameter indicating whether to plot the colollipopplot with or without the forest plot. Default is TRUE.
#' @param ... `pretty_lollipop_plot` arguments, see that function for more info on available parameters.
#' 
#' @return A mirrored lollipop plot.
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' library(GAMBLR.open)
#' 
#' #get meta data (DLBCL_Hilton)
#' meta = GAMBLR.open::get_gambl_metadata()
#' metadata = dplyr::filter(meta, cohort %in% "DLBCL_GenomeCanada")
#' maf_df = GAMBLR.open::get_coding_ssm(
#'   these_samples_metadata = metadata
#' )
#' pretty_colollipop_plot_result <- pretty_colollipop_plot(
#'  maf_df = maf_df,
#'  these_samples_metadata = metadata,
#'  comparison_column = "COO_consensus",
#'  comparison_values = c("ABC", "GCB"),
#'  gene = "PIM1",
#'  label_threshold = 2, 
#'  forestarg = TRUE, 
#'  show_rate = TRUE, 
#'  compare_distributions = TRUE
#' )
#' 
#' pretty_colollipop_plot_result$plot

pretty_colollipop_plot <- function(
    maf_df,
    these_samples_metadata,
    comparison_column,
    comparison_values,
    gene,
    include_silent = FALSE,
    label_threshold = 5,
    plot_title,
    refseq_id,
    forestarg = FALSE,
    compare_distributions = FALSE,
    show_rate = FALSE,
    font = "sans",
    domain_label_size = 3,
    x_axis_size = 10,
    ...
) {

    # check for required arguments
    required <- c(
        "maf_df", 
        "these_samples_metadata", 
        "comparison_column"
        )

    defined <- names(as.list(match.call())[-1])

    if(any(!required %in% defined)) {
        stop("Please provide maf_df, these_samples_metadata, and comparison_column to the function.")
    }

    # If no comparison_values are specified, derive the comparison_values from the specified comparison_column
    if(missing(comparison_values)){
        if(class(these_samples_metadata[[comparison_column]]) == "factor"){
            comparison_values = levels(these_samples_metadata[[comparison_column]])
        } else {
            comparison_values = unique(these_samples_metadata[[comparison_column]])
        }
    }

    # Ensure there are only two comparison_values
    if(length(comparison_values) != 2) stop(paste0("Your comparison must have two values. \nEither specify comparison_values as a vector of length 2 or subset your metadata so your comparison_column has only two unique values or factor levels."))

    # Subset the metadata to the specified comparison_values and the maf to the remaining sample_ids
    meta1 <- these_samples_metadata[these_samples_metadata[[comparison_column]] %in% comparison_values[1], ]
    meta2 <- these_samples_metadata[these_samples_metadata[[comparison_column]] %in% comparison_values[2], ]

    ##### Subset the maf to gene and variants of interest and split by comparison values #####
    
    if(include_silent){
        variants <- coding_class
    } else {
        variants <- coding_class[!coding_class %in% c(
            "Silent",
            "Splice_Region" # Is this useful? Aren't these all intronic? Won't affect cds
        )]
    }

    # Subset maf to only samples in the comparison values
    ssm1 <- maf_df %>%
        dplyr::filter(
            Tumor_Sample_Barcode %in% meta1$Tumor_Sample_Barcode
        ) %>% 
        filter(
            Hugo_Symbol == gene, 
            Variant_Classification %in% variants
        )

    ssm2 <- maf_df %>%
        dplyr::filter(
            Tumor_Sample_Barcode %in% meta2$Tumor_Sample_Barcode
        )%>% 
        filter(
            Hugo_Symbol == gene, 
            Variant_Classification %in% variants
        )

    # Ensure dimensions are greater than zero otherwise no variants are returned.
    dim1 <- dim(ssm1)[1]  
    dim2 <- dim(ssm2)[1]

    if(dim1 == 0 | dim2 == 0) {
        stop(paste0("There are no variants returned for one or both of the subsets defined by the specified comparison values. "))
    }
    
    # Generate a unified plot title
    if(missing(plot_title)){
        plot_title <- paste0(
            gene, 
            ": ", 
            comparison_values[1], 
            " vs. ", 
            comparison_values[2]
        )
    } 

    # Arguments to pass into pretty_lollipop_plot
    lollipopplot_args <- list(...)

    # Get gene_counts data for the first plot (without generating the plot)
    lp1 <- do.call(pretty_lollipop_plot, c(
        list(
            maf_df = ssm1,
            these_samples_metadata = these_samples_metadata,
            gene = gene, 
            plot_title = paste0(gene, " in ", comparison_values[1])
        ), 
        lollipopplot_args
    ))
    lp1_gene_counts_data <- as.data.frame(lp1$gene_counts)%>% 
        mutate(source = "Sample 1")

    # Get gene_counts data for the second plot
    lp2 <- do.call(pretty_lollipop_plot, c(
        list(
            maf_df = ssm2,
            these_samples_metadata = these_samples_metadata,
            gene = gene, 
            plot_title = paste0(gene, " in ", comparison_values[2]), 
            mirror = TRUE
        ), 
        lollipopplot_args
    ))
    lp2_gene_counts_data <- as.data.frame(lp2$gene_counts)%>% 
        mutate(source = "Sample 2")

    # Combine data for both plots
    combined_gene_counts <- rbind(
        lp1_gene_counts_data, 
        lp2_gene_counts_data 
    )

    # Harmonize the y-axis scales for plotting
    y_limit <- max(combined_gene_counts$mutation_count) * 1.2
    y_breaks <- seq(
            0,
            y_limit,
            by = ifelse(y_limit <= 5, 1, round(y_limit * 1.2 / 5))
        )
    
    if(compare_distributions == TRUE) {
        
        # Gene KS-Test
        # Compare the distribution of variants along the length of the genes between the groups
        ks_test_result <- suppressWarnings(ks.test(
            rep(lp1$gene_counts$AA, lp1$gene_counts$mutation_count),
            rep(lp2$gene_counts$AA, lp2$gene_counts$mutation_count) 
        ))

        gene_p_value <- ks_test_result$p.value
        

        # Domain(s) KS-Test
        domain_names <- c()
        domain_list <- list()
        domain_data <- lp1$domain_data
        domain_data$p.value <- NA
        # Within each domain, use the KS test to determine if the distribution of variants is different
        for (i in 1:nrow(domain_data)) {
            domain_name <- domain_data$text.label[i]
            min_val <- domain_data$start.points[i]
            max_val <- domain_data$end.points[i]

            # Subset data for the current domain
            domain_subset1 <- lp1$gene_counts %>% 
                filter(AA >= min_val & AA <= max_val)
                
            domain_subset2 <- lp2$gene_counts %>% 
                filter(AA >= min_val & AA <= max_val)
                
            if(length(domain_subset1$AA) == 0 | length(domain_subset2$AA) == 0) {
                message("Skipping domain ", domain_name, " as it does not have data for both subgroups")
            }else{
                # Compare the distribution of variants along the length of the domains between the groups
                domain_ks_test <- suppressWarnings(ks.test(
                        rep(domain_subset1$AA, domain_subset1$mutation_count),
                        rep(domain_subset2$AA, domain_subset2$mutation_count) 
                    ))
                domain_data$p.value[i] <- domain_ks_test$p.value
            }
        }
        
        domain_data <- domain_data %>%         
            mutate(
                p_value = case_when(
                    p.value < 0.001 ~ "***",
                    p.value < 0.01 ~ "**",
                    p.value < 0.05 ~ "*", 
                    TRUE ~ ""
                )
            )

        domain_plot <- draw_domain_plot(
            domain_data, 
            domain_label_size = domain_label_size,
            x_axis_size = x_axis_size,
            font = font
            )
        plot_subtitle <- paste0("KS Test P=", signif(gene_p_value, digits = 2))
    } else{
        domain_plot <- lp1$domain_plot
        domain_data <- lp1$domain_data
        plot_subtitle = ""
    }
    
    if(show_rate == TRUE){
        Somatic_Mutation_Numerator <- c(length(unique(ssm1$Tumor_Sample_Barcode)), length(unique(ssm2$Tumor_Sample_Barcode)))
        Somatic_Mutation_Denominator <- c(length(unique(meta1$sample_id)), length(unique(meta2$sample_id)))
        Somatic_Mutation_Rate <- round(Somatic_Mutation_Numerator / Somatic_Mutation_Denominator * 100, 1)
        
        yaxis_labels <- paste0(
            comparison_values,
            " Mutation Count\nRate ", 
            Somatic_Mutation_Rate, 
            "% (", 
            Somatic_Mutation_Numerator, 
            "/", 
            Somatic_Mutation_Denominator, 
            ")"
        )
    } else {
        yaxis_labels <- paste(comparison_values, "\nMutation Count")
    }
    
    # Arrange lollipop and domain plots together
    colollipop_plot <- ggpubr::ggarrange(
        lp1$mutation_plot + scale_y_continuous(breaks = y_breaks, limits = c(0, y_limit))+ 
            ylab(yaxis_labels[1]) + 
            ggtitle(plot_title, subtitle = plot_subtitle) + 
            theme(plot.subtitle = element_text(hjust = 0.5)), 
        domain_plot, 
        lp2$mutation_plot + 
            scale_y_continuous(breaks = y_breaks, limits = c(y_limit, 0), transform = "reverse") + 
            theme(
                plot.title = element_blank(), 
                plot.subtitle = element_blank()
                ) + 
            ylab(yaxis_labels[2]), 
        ncol = 1, 
        align = "v", 
        heights = c(3, 1, 3),
        common.legend = TRUE, 
        legend = "right"
    ) 
    
    
    outputs <- list(
        plot = colollipop_plot,
        combined_gene_counts = combined_gene_counts, 
        lollipop1 = lp1, 
        lollipop2 = lp2, 
        domain_data = domain_data
    )

    if (forestarg == TRUE) {
        forest_plot <- prettyForestPlot(
            maf = maf_df,
            metadata = these_samples_metadata,
            genes = gene,
            comparison_column = comparison_column,
            comparison_values = comparison_values,
            separate_hotspots = FALSE,
            comparison_name = paste0(
                comparison_values[1],
                " vs. ",
                comparison_values[2]
            )
        )

    # arrange colollipop plot and forest plot together one over the other
        combine_lollipop_forest <- ggpubr::ggarrange(
            colollipop_plot, 
            ggpubr::ggarrange(
                forest_plot$forest + 
                    ggtitle(
                        paste0("Fisher's Exact Test P=", signif(forest_plot$fisher$p.value, 2))
                        ) + 
                    theme(
                        text = element_text(family = font, face = "plain"), 
                        axis.text = element_text(family = font, face = "plain"), 
                        axis.title.x = element_text(family = font, face = "plain"), 
                        plot.title = element_text(family = font, face = "plain"),
                        axis.title.y = element_blank()
                        ) ,        
                forest_plot$bar + 
                    theme(
                        text = element_text(family = font, face = "plain"), 
                        axis.text = element_text(family = font, face = "plain"), 
                        axis.title.x = element_text(family = font, face = "plain")
                        ), 
                widths = c(1, 0.6),
                common.legend = TRUE,
                legend = "bottom",
                align = "hv"),          
            ncol = 1, 
            nrow = 2,
            heights = c(3, 1), 
            align = "v"     
        )

    
        outputs$plot <- combine_lollipop_forest
        outputs$forest_plot <- forest_plot
    } 
    return(outputs)
}
