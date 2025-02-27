#' @title Pretty Lollipop Plot.
#'
#' @description Generates a visually appealing lollipop plot.
#'
#' @details Retrieve maf data of a specific sample or a set of samples. A gene of interest
#' can then be visualized with the given maf data. Silent mutations can be visualized setting
#' include_silent to TRUE. 
#'
#' @param maf_df A data frame containing the mutation data.
#' @param these_samples_metadata Required argument. A data.frame with metadata.
#' @param gene The gene symbol to plot.
#' @param plot_title Optional, the title of the plot. Default is gene.
#' @param refseq_id Insert a specific NM_xxx value of interest
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is FALSE. 
#' @param label_threshold Threshold for labels to appear on plot. Default set to 5.
#' @param plotarg Logical parameter indicating whether to plot the lollipopplot or return the data in data frame format. Default is TRUE.
#' @param mirrorarg Logical paramter for when mirroring lollipop data in prety_co_lollipop plot. Default is FALSE.
#' @param combined_gene_counts A dataframe containing data for a mirrored lollipop analysis.
#' @param meta1_counter A dataframe for calculating Somatic Mutation Rate in `pretty_lollipop_plot`.
#' @param meta2_counter A dataframe for calculating Somatic Mutation Rate in `pretty_lollipop_plot`.
#' @param Sample1 A label for displaying Somatic Mutation Rate in `pretty_lollipop_plot`.
#' @param Sample2 A label for displaying Somatic Mutation Rate in `pretty_lollipop_plot`.
#' @param font Customizable font size for the CoLollipop plot somatic mutation rate and comparison labels. Default is 11pt font.
#' 
#' @return A lollipop plot.
#'
#' @import dplyr ggplot2 GAMBLR.data
#' @export
#'
#' @examples
#' library(GAMBLR.open)

#' #get meta data (BL_Thomas)
#' metadata <- get_gambl_metadata() %>%
#'     dplyr::filter(cohort == "BL_Thomas")
#' 
#' maf_df <- get_coding_ssm(
#'     these_samples_metadata = metadata
#' )
#'
#' #construct pretty_lollipop_plot.
#' lolipop_result <- pretty_lollipop_plot(maf_df, these_samples_metadata = metadata, "DDX3X")
#' print(lolipop_result)
#' 
#' \dontrun{
#'  # Or, with GAMBLR.results:
#'  library(GAMBLR.results)
#'  metadata <- get_gambl_metadata() %>%
#'     dplyr::filter(pathology == "BL")
#' 
#' maf_df <- get_all_coding_ssm(
#'     these_samples_metadata = metadata
#' )
#' lolipop_result <- pretty_lollipop_plot(maf_df, these_samples_metadata = metadata, "DDX3X")
#' lolipop_result
#'}
pretty_lollipop_plot <- function(
    maf_df,
    these_samples_metadata,
    gene,
    plot_title,
    refseq_id,
    include_silent = FALSE,
    label_threshold = 5,
    plotarg = TRUE,
    mirrorarg = FALSE,
    combined_gene_counts,
    meta1_counter,
    meta2_counter,
    Sample1 = Sample1,
    Sample2 = Sample2,
    font
){
    if(missing(gene)){
        stop("Please provide a gene...")
    }

    if(missing(plot_title)){
        plot_title = gene
    }

    maf_df <- as.data.frame(maf_df)

    # Specifying noncoding regions with coding_class as a bundled object with GAMBLR.helpers
    if(include_silent){
        variants <- coding_class
    } else {
        variants <- coding_class[!coding_class %in% c(
            "Silent",
            "Splice_Region"
        )]
    }

    if (!"RefSeq" %in% colnames(maf_df)) {
        stop("Error: The provided maf_df is missing the 'RefSeq' column. 
              Ensure that get_ssm_by_samples() is run with basic_columns = FALSE.")
    }

    max_splits <- maf_df %>%
        pull(RefSeq) %>%
        strsplit(",") %>%
        sapply(length) %>%
        max()

    new_columns <- paste0("RefSeq_", seq_len(max_splits))

    maf_df <- separate(
        maf_df,
        col = RefSeq,
        into = new_columns,
        sep = ",",
        remove = TRUE,
        convert = FALSE
    )

    maf_df <- maf_df %>%
        pivot_longer(
        cols = new_columns,       
        names_to = "RefSeq_col", # Name for the column containing old column names
        values_to = "RefSeq", # Name for the new single column with values
        values_drop_na = TRUE      
        ) %>%
        select(-RefSeq_col) # Remove the column names   

    maf_df <- maf_df %>% 
        mutate(RefSeq = sub("\\.\\d+$", "", RefSeq)) # Remove decimal places to match with protein_domains

    # Subset the maf to gene and variants of interest
    nc_maf_df <- maf_df %>%
        filter(
            Hugo_Symbol == gene
        ) %>% 
        filter(
            Variant_Classification %in% variants 
        )

    # Filter for the specific gene 
    gene_df <- nc_maf_df %>% 
        mutate(
            AA = as.numeric(
                gsub(
                    "[^0-9]+", 
                    "", 
                    gsub(
                        "([0-9]+).*", 
                        "\\1", 
                        HGVSp_Short
                    )
                )
            )
        ) %>%
        mutate(AA_position = gsub("^p\\.", "", HGVSp_Short)) %>%
        arrange(AA)

    repetition_factor <- length(unique(gene_df$RefSeq)) 

    # Count mutations at each position
    gene_counts <- gene_df %>% 
        group_by(
            AA, 
            Start_Position, 
            End_Position, 
            Variant_Classification, 
            Reference_Allele, 
            Tumor_Seq_Allele2
        ) %>%
        arrange(AA) %>%
        mutate(mutation_count = n()/repetition_factor) 

    if (mirrorarg == TRUE){
        gene_counts <- combined_gene_counts %>%
        # Specify the amino acid change as a label based on the label_threshold value
        mutate(label = as.character(ifelse(
            mutation_count >= label_threshold, 
            AA_position, 
            ""
        )))
    } else {
        gene_counts <- gene_counts %>%
        # Specify the amino acid change as a label based on the label_threshold value
        mutate(label = as.character(ifelse(
            mutation_count >= label_threshold, 
            AA_position, 
            ""
        )))
    }

    # protein_domains a bundled object with GAMBLR.data
    protein_domain_subset <- subset(
        protein_domains, 
        HGNC == gene
    )

    # Initial check: Ensure NM_XXX in protein_domain_subset matches gene_df$RefSeq
    matching_values <- protein_domain_subset$refseq.ID[protein_domain_subset$refseq.ID %in% unique(gene_df$RefSeq)]
    # Trigger error only if there are no matches
    if (length(matching_values) == 0) {
        stop(paste(
            "Error: None of the RefSeq IDs in gene_df match any refseq.ID in protein_domain_subset",
            "\nThe following refseq IDs are available in protein_domain_subset:",
            paste(unique(protein_domain_subset$refseq.ID), collapse = ", ")
        ))
    }

    if (missing(refseq_id)){
        if (length(unique(protein_domain_subset$refseq.ID)) == 1){ # Case 1: No refseq_id provided. Gene has only one NM_XXX
            protein_domain_subset <- protein_domain_subset 
        } else if (length(unique(protein_domains$refseq.ID)) > 1){ # Case 2: No refseq_id provided. Gene has multiple NM_XXX; choose the first one
            protein_domain_subset <- protein_domain_subset %>%
                filter(refseq.ID == protein_domain_subset$refseq.ID[1])  
        }
    } else if (refseq_id %in% protein_domain_subset$refseq.ID){ # Case 3: User-provided refseq_id matches one in protein_domains
        protein_domain_subset <- protein_domain_subset %>%
            filter(refseq.ID == refseq_id)
    } else { # Case 4: refseq_id not found in protein_domains
        stop(paste("Error: refseq_id", refseq_id, "is not found in protein_domains.Rd", 
            "\nThe following refseq IDs are available in protein_domains.Rd:",
            paste(unique(protein_domain_subset$refseq.ID), collapse = ", ")))
    }

    domain_data <- protein_domain_subset %>% 
        data.frame(
            start.points = protein_domain_subset$Start,
            end.points = protein_domain_subset$End,
            text.label = protein_domain_subset$Label,
            color = protein_domain_subset$Label
        )

    domain_data$text.position <- (domain_data$start.points + domain_data$end.points) / 2

    # Determine the x-axis range
    x_max <- max(max(domain_data$end.points), 
        max(
            gene_counts$AA, 
            na.rm = TRUE
        )
    )
    x_min <- 0

    # get_gambl_colours() from GAMBLR.helpers
    colours_manual <- get_gambl_colours("mutation")
    domain_palette <- unname(get_gambl_colours("domains"))[1:length(unique(domain_data$text.label))]

    # Somatic mutation statistic
    if (mirrorarg == TRUE){
        
        # Initialize vectors
        Somatic_Mutation_Numerator <- numeric(2)
        Somatic_Mutation_Denominator <- numeric(2)
        Somatic_Mutation_Rate <- numeric(2)

        # Somatic Mutation rate for lp1
        Somatic_Mutation_Numerator[1] <- combined_gene_counts %>%
            filter(
                Hugo_Symbol == gene,
                source == "Sample 1"
            ) %>% 
            filter(
                Variant_Classification %in% variants 
            ) %>% 
            distinct(Tumor_Sample_Barcode) %>% 
            nrow()

        Somatic_Mutation_Denominator[1] <- meta1_counter

        Somatic_Mutation_Rate[1] <- Somatic_Mutation_Numerator[1]/Somatic_Mutation_Denominator[1] *100
        Somatic_Mutation_Rate[1] <- round(Somatic_Mutation_Rate[1], 2)

        # Somatic Mutation rate for lp2
        Somatic_Mutation_Numerator[2] <- combined_gene_counts %>%
            filter(
                Hugo_Symbol == gene,
                source == "Sample 2"
            ) %>% 
            filter(
                Variant_Classification %in% variants 
            ) %>% 
            distinct(Tumor_Sample_Barcode) %>% 
            nrow()

        Somatic_Mutation_Denominator[2] <- meta2_counter

        Somatic_Mutation_Rate[2] <- Somatic_Mutation_Numerator[2]/Somatic_Mutation_Denominator[2] *100
        Somatic_Mutation_Rate[2] <- round(Somatic_Mutation_Rate[2], 2)
    } else {
        Somatic_Mutation_Numerator <- maf_df %>%
            filter(
                Hugo_Symbol == gene
            ) %>% 
            filter(
                Variant_Classification %in% variants 
            ) %>% 
            distinct(Tumor_Sample_Barcode) %>% 
            nrow()

        Somatic_Mutation_Denominator <- length(unique(these_samples_metadata$Tumor_Sample_Barcode))

        Somatic_Mutation_Rate <- Somatic_Mutation_Numerator/Somatic_Mutation_Denominator *100
        Somatic_Mutation_Rate <- round(Somatic_Mutation_Rate, 2)
    }

    # KS-Test 
    if(mirrorarg == TRUE){
        # Gene KS-Test
        aa_frequency_data <- combined_gene_counts %>%
            group_by(AA, source) %>%
            summarise(
                freq = n()/repetition_factor,
                .groups = "drop"
            )

        ks_test_result <- ks.test(
            aa_frequency_data$freq[aa_frequency_data$source == "Sample 1"], 
            aa_frequency_data$freq[aa_frequency_data$source == "Sample 2"]
        )

        gene_p_value <- ks_test_result$p.value
        
        if(gene_p_value <= 0.05){
            gene_p_value <- gene_p_value
        } else {
            gene_p_value <- NA
        }

        # Domain(s) KS-Test
        domain_names <- c()
        domain_list <- list()
        domain_data$p_value <- NA
    
        for (i in 1:nrow(domain_data)) {
            domain_name <- domain_data$text.label[i]
            min_val <- domain_data$start.points[i]
            max_val <- domain_data$end.points[i]

            # Subset data for the current domain
            domain_subset <- combined_gene_counts %>%
                filter(AA >= min_val & AA <= max_val)

            if (nrow(domain_subset) > 0) {

                domain_frequency_data <- domain_subset %>%
                    group_by(
                        AA, 
                        source
                    ) %>%
                    summarise(
                        freq = n()/repetition_factor,
                        .groups = "drop"
                    )

                domain_ks_test <- tryCatch({
                    ks.test(
                        domain_frequency_data$freq[domain_frequency_data$source == "Sample 1"],
                        domain_frequency_data$freq[domain_frequency_data$source == "Sample 2"]
                    )
                }, error = function(e) {
                    list(p.value = NA)
                })

                # Store the p-value in domain_data
                domain_data$p_value[i] <- domain_ks_test$p.value

                if (!is.na(domain_data$p_value[i])) {
                    if (domain_data$p_value[i] <= 0.001) {
                        domain_data$p_value[i] <- "***"
                    } else if (domain_data$p_value[i] <= 0.01) {
                        domain_data$p_value[i] <- "**"
                    } else if (domain_data$p_value[i] <= 0.05) {
                        domain_data$p_value[i] <- "*"
                    } else {
                        domain_data$p_value[i] <- NA
                    }
                }            

            } else {
                print(paste(
                    "No mutation data for domain", 
                    domain_name
                ))
            }
        }
    }

    if (mirrorarg == TRUE) {

        # Processing p-values for domain plot
        for (i in 1:nrow(domain_data)) {
            if (!is.na(domain_data$p_value[i])) {
                domain_data$p_value[i] <- as.character(domain_data$p_value[i])
            } else {
                domain_data$p_value[i] <- " "
            }
        }

        # Create the domain plot
        domain_plot <- ggplot() + 
            geom_rect(
                aes(
                    xmin = x_min, 
                    xmax = x_max, 
                    ymin = -0.2, 
                    ymax = 0.2
                    ), 
                fill = "lightgrey", 
                color = "lightgrey"
            ) + 
            geom_rect(
                data = domain_data, 
                aes(
                    xmin = start.points, 
                    xmax = end.points, 
                    ymin = -0.4, 
                    ymax = 0.4, 
                    fill = color
                    ), 
                color = "black", 
                show.legend = FALSE
            ) +
            geom_text(
                data = domain_data, 
                aes(
                    x = text.position, 
                    y = 0,
                    angle = 90, 
                    label = paste0(
                        text.label,
                        "\n",
                        domain_data$p_value
                        )
                )
            ) +
            scale_fill_manual(values = domain_palette) +
            theme_void() +
            theme(
                axis.text.x = element_text(vjust = rel(0.5)),
                text = element_text(family = "arial")
            )

        # Function for mutation lollipop plot
        mutation_plot <- function(gene_counts){
            ggplot() + 
            geom_segment(
                data = gene_counts, 
                aes(
                    x = AA, 
                    xend = AA, 
                    y = 0, 
                    yend = mutation_count
                )
            ) +
            geom_point(
                data = gene_counts, 
                aes(
                    x = AA, 
                    y = mutation_count, 
                    color = Variant_Classification, 
                    size = mutation_count
                )
            ) + 
            ggpubr::theme_pubr() 
        }

        # common legend
        co_legend_data <- mutation_plot(gene_counts) +
            guides(size = "none") +  
            scale_color_manual(name = "Mutation Type", values = colours_manual) +
            theme(legend.position = "right") +
            theme(text = element_text(family = "arial"))
        gg_co_legend <- as_ggplot(ggpubr::get_legend(co_legend_data))
      
        mutation_plot1 <- mutation_plot(gene_counts %>% filter(source == "Sample 1")) +
            scale_y_continuous(
                breaks = seq(
                    0,
                    max(gene_counts$mutation_count),
                    by = 1
                )
            ) +
            scale_color_manual(name = "Mutation Type", values = colours_manual) +
            ylab("Mutation Count") +
            xlab("") +
            xlim(0, x_max) +
            theme(
                axis.line.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                text = element_text(family = "arial")
            ) 

        mutation_plot2 <- mutation_plot(gene_counts %>% filter(source == "Sample 2")) +
            scale_y_reverse(
                breaks = seq(
                    0,
                    max(gene_counts$mutation_count),
                    by = 1
                )
            ) +
            scale_color_manual(name = "Mutation Type", values = colours_manual) +
            ylab("Mutation Count") +
            xlab("") +
            xlim(0, x_max) +
            theme(
                axis.line.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                text = element_text(family = "arial")
            ) 
        
        domain_plot <- domain_plot + 
            theme(
                axis.title.x = element_blank(),
                axis.text.x = element_text(),
                axis.ticks.x = element_line(color = "black"),
                plot.margin = margin(t = -20), # Reduce space between plots
                text = element_text(family = "arial")
            )

        # Combine mutation plots and domain plot
        plot <- ggpubr::ggarrange(ggpubr::ggarrange(
                mutation_plot1, 
                domain_plot, 
                mutation_plot2,
                ncol = 1, 
                heights = c(3, 1, 3),
                labels = c(
                    paste0(
                            '               "', 
                            Sample1, 
                            '"',
                            "\n",
                            "Somatic Mutation Rate ", 
                            Somatic_Mutation_Rate[1], 
                            "%",
                            "\n",
                            "             N = ", 
                            Somatic_Mutation_Denominator[1]
                    ), 
                    "", 
                    paste0(
                            '               "', 
                            Sample2, 
                            '"',
                            "\n",
                            "Somatic Mutation Rate ", 
                            Somatic_Mutation_Rate[2],
                            "%",
                            "\n", 
                            "              N = ",
                            Somatic_Mutation_Denominator[2]
                    )   
                ),
                font.label = list(size = font),  # Adjust font size 
                label.y = c(0.85, 0, 0.35),
                label.x = c(0.7, 0, 0.7),
                legend = "none"
            ),
        gg_co_legend, 
        ncol = 2,
        widths = c(5,1)
        ) %>% 
        ggpubr::annotate_figure(
            plot,
            top = text_grob(
                paste0(
                    plot_title,
                    "\n",
                    ifelse(
                        !is.na(gene_p_value), 
                        paste0("p = ", round(gene_p_value, 3)), 
                        ""
                    )
                ),
                color = "black", 
                face = "bold", 
                size = 14
            ) 
            
        )
    } else {
        plot <- ggplot() +
            geom_segment(
                data = gene_counts, 
                aes(
                    x = AA, 
                    xend = AA, 
                    y = 0, 
                    yend = mutation_count
                    )
            ) +
            geom_point(
                data = gene_counts, 
                aes(
                    x = AA, 
                    y = mutation_count, 
                    color = Variant_Classification, 
                    size = mutation_count
                    )
            ) +
            geom_text(data = gene_counts, aes(
                x = AA, 
                y = mutation_count, 
                label = label, 
                angle = 45, 
                hjust = -0.25
            )) +
            # Background rectangle for regions without domain data
            geom_rect(
                aes(
                    xmin = x_min, 
                    xmax = x_max, 
                    ymin = -0.2, 
                    ymax = 0.2
                    ), 
                fill = "darkgrey", 
                color = "darkgrey"
            ) + 
            # Domain rectangles
            geom_rect(
                data = domain_data, 
                aes(
                    xmin = start.points, 
                    xmax = end.points, 
                    ymin = -0.4, 
                    ymax = 0.4, 
                    fill = color
                    ), 
                color = "black", 
                show.legend = FALSE
            ) +
            scale_fill_manual(values = domain_palette) + 
            geom_text(
                data = domain_data, 
                aes(
                    x = text.position, 
                    y = 0, 
                    angle = 90, 
                    label = text.label
                    )
            ) +
            annotate(
                "text",
                x = 0,
                y = max(gene_counts$mutation_count) * 1.1,
                label = paste0(
                    "Somatic Mutation Rate ", 
                    Somatic_Mutation_Rate, 
                    "%"
                ),
                hjust = 0
            ) +
            labs(
                x = "AA Position", 
                y = "Mutation Count", 
                title = paste0(
                    plot_title
                )
            ) + 
            scale_y_continuous(
                breaks = seq(
                    0,
                    max(gene_counts$mutation_count),
                    by = 1
                )
            )  +
            scale_size_continuous(
                breaks = seq(
                    floor(min(gene_counts$mutation_count)), 
                    ceiling(max(gene_counts$mutation_count)), 
                    by = 1
                ),
                labels = function(x) format(x, scientific = FALSE)  
            )  +
            theme_bw() +
            theme(
                plot.title = element_text(
                    hjust = 0.5
                    ),
                axis.text.x = element_text(
                  angle = 0, 
                  hjust = 1
                  ),
                text = element_text(family = "arial")
            ) +
            scale_color_manual(
              name = "Legend", 
              values = colours_manual
            )
    }


    if (plotarg == TRUE) {
        return(plot)
    } else {
        return(gene_counts)
    }
}
