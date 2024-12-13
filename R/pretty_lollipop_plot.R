#' @title Pretty Lollipop Plot.
#'
#' @description Generates a visually appealing lollipop plot.
#'
#' @details Retrieve maf data of a specific sample or a set of samples. A gene of interest
#' can then be visualized with the given maf data. Silent mutations can be visualized setting
#' include_silent to TRUE.
#'
#' @param maf_df A data frame containing the mutation data.
#' @param gene The gene symbol to plot.
#' @param plot_title Optional, the title of the plot. Default is gene.
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is FALSE.
#'
#' @return A lollipop plot.
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#'
#' # get meta data (BL_Thomas)
#' metadata <- get_gambl_metadata() %>%
#'     filter(cohort == "BL_Thomas")
#'
#' maf_df <- get_ssm_by_samples(
#'     these_samples_metadata = metadata
#' )
#'
#' # construct pretty_lollipop_plot.
#' lolipop_result <- pretty_lollipop_plot(maf_df, "MYC")
#'
pretty_lollipop_plot <- function(
    maf_df,
    gene = NULL,
    refseq = NULL,
    metadata = NULL,
    split_by_column = NULL,
    split_column_values = NULL,
    plot_title,
    include_silent = FALSE,
    label_threshold = 5) {
    if (missing(gene)) {
        stop("Please provide a gene...")
    }
    if (missing(plot_title)) {
        plot_title <- gene
    }

    # Load the maf_df as a data.frame
    maf_df <- as.data.frame(maf_df)
    if (include_silent) {
        variants <- coding_class
    } else {
        variants <- coding_class[!coding_class %in% c(
            "Silent",
            "Splice_Region"
        )]
    }
    # Subset metadata to include stratifications of interest
    if (!is.null(metadata)) {
        if (length(split_by_column) != 1) {
            stop("split_by_column should be a character vector of length 1. ")
        } else if (length(split_column_values) != 2) {
            stop("split_column_values should be a character vector of length 2.")
        }
        metadata <- metadata[metadata[[split_by_column]] %in% split_column_values, ]
        maf_df <- maf_df %>%
            filter(Tumor_Sample_Barcode %in% metadata$sample_id) %>%
            left_join(select(metadata, Tumor_Sample_Barcode = sample_id, all_of(split_by_column)), by = "Tumor_Sample_Barcode")
    } else {
        # make a dummy grouping column
        split_by_column <- "dummy"
        split_column_values <- "dummy"
        maf_df$dummy <- "dummy"
    }
    # Subset the maf to gene and variants of interest
    nc_maf_df <- maf_df %>%
        filter(Hugo_Symbol == gene) %>%
        filter(Variant_Classification %in% variants)
    gene_df <- nc_maf_df %>%
        mutate(position = as.numeric(gsub(
            "[^0-9]+",
            "", gsub("([0-9]+).*", "\\1", HGVSp_Short)
        ))) %>%
        mutate(AA = gsub("^p\\.", "", HGVSp_Short)) %>%
        arrange(position)
    # Count mutations at each position
    gene_counts <- gene_df %>%
        group_by(
            AA, position, Start_Position, End_Position,
            Variant_Classification, Reference_Allele, !!sym(split_by_column)
        ) %>%
        summarize(mutation_count = n()) %>%
        arrange(position) %>%
        # Specify the amino acid change as a label based on the label_threshold value
        mutate(label = ifelse(mutation_count >= label_threshold, AA, ""))
    # Subset the protein domains table
    protein_domain_subset <- subset(protein_domains, HGNC ==
        gene)
    unique_refseq <- unique(protein_domain_subset$refseq.ID)
    if (length(unique_refseq) > 1 & is.null(refseq)) {
        refseq <- unique_refseq[1]
        warning(paste(
            "Multiple RefSeq IDs found for gene ", gene,
            ". Using first RefSeq ID: ", refseq
        ))
        warning(paste(
            "To specify a different RefSeq ID, use the refseq argument. Available RefSeq IDs: ", unique_refseq
        ))
    }
    # Subset to a single refseq ID
    protein_domain_subset <- subset(protein_domain_subset, refseq.ID == refseq)
    # Process protein domains
    domain_data <- protein_domain_subset %>% data.frame(
        start.points = protein_domain_subset$Start,
        end.points = protein_domain_subset$End, text.label = protein_domain_subset$Label,
        color = protein_domain_subset$Label
    )
    domain_data$text.position <- (domain_data$start.points +
        domain_data$end.points) / 2

    # Set x_max to length of the protein model
    x_max <- unique(domain_data$aa.length)
    x_min <- 0
    # Set y_max to the nearest multiple of 5 above the maximum mutation count
    y_max <- ceiling(max(gene_counts$mutation_count) / 5) * 5

    # Load gambl mutation colours
    colours_manual <- get_gambl_colours("mutation")

    # Create the domain plot
    domain_plot <- ggplot() +
        geom_rect(
            aes(
                xmin = x_min,
                xmax = x_max, ymin = -0.2, ymax = 0.2
            ),
            fill = "lightgrey",
            color = "lightgrey"
        ) +
        geom_rect(
            data = domain_data, aes(
                xmin = start.points,
                xmax = end.points, ymin = -0.4, ymax = 0.4, fill = color
            ),
            color = "black", show.legend = FALSE
        ) +
        geom_text(
            data = domain_data,
            aes(x = text.position, y = 0, label = text.label)
        ) +
        scale_fill_brewer("set3") +
        theme_void() +
        theme(
            axis.text.x = element_text(vjust = rel(0.5))
        )

    # Function for mutation lollipop plot
    mutation_plot_fun <- function(gene_counts) {
        ggplot() +
            geom_segment(data = gene_counts, aes(
                x = position,
                xend = position, y = 0, yend = mutation_count
            )) +
            geom_point(
                data = gene_counts,
                aes(
                    x = position, y = mutation_count, color = Variant_Classification,
                    size = mutation_count
                )
            ) +
            geom_text(data = gene_counts, aes(
                x = position, y = mutation_count, label = label, angle = 45, hjust = -0.25
            )) +
            ggpubr::theme_pubr() +
            scale_color_manual(name = "Mutation Type", values = colours_manual) +
            ylab("Mutation Count") +
            xlab("") +
            xlim(0, x_max) +
            theme(
                axis.line.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                legend.position = "right"
            ) +
            guides(size = "none")
    }

    if (length(split_column_values) == 2) {
        # Co-lollipop plot when split_by_column is provided
        mutation_plot1 <- mutation_plot_fun(gene_counts %>% filter(!!sym(split_by_column) == split_column_values[1])) +
            ylim(0, y_max)
        mutation_plot2 <- mutation_plot_fun(gene_counts %>% filter(!!sym(split_by_column) == split_column_values[2])) +
            scale_y_reverse(limits = c(ymax, 0))
        plot <- ggpubr::ggarrange(
            mutation_plot1, domain_plot, mutation_plot2,
            ncol = 1, align = "hv",
            heights = c(3, 1, 3),
            labels = c(split_column_values[1], "", split_column_values[2]),
            label.y = c(1, 0, 0.2),
            label.x = c(0.8, 0, 0.8),
            common.legend = TRUE,
            legend = "right"
        )
    } else {
        # Single lollipop plot when split_by_column is not provided
        mutation_plot <- mutation_plot_fun(gene_counts) +
            ggtitle(plot_title) +
            ylim(0, y_max)
        plot <- ggpubr::ggarrange(
            mutation_plot, domain_plot,
            ncol = 1, align = "hv",
            heights = c(5, 1)
        )
    }
    return(plot)
}
