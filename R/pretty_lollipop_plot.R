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
#' @param refseq_id Insert a specific NM_xxx value of interest
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is FALSE. 
#' @param plotarg Logical parameter indicating whether to plot the lollipopplot or return the data in data frame format. Default is TRUE.
#' @param mirrorarg Logical paramter for when mirroring lollipop data in prety_co_lollipop plot. Default is FALSE.
#' @param combined_gene_counts A dataframe containing data for a mirrored lollipop analysis.
#' @param meta1_counter A dataframe for calculating Somatic Mutation Rate in `pretty_lollipop_plot`.
#' @param meta2_counter A dataframe for calculating Somatic Mutation Rate in `pretty_lollipop_plot`.
#' @param Sample1 A label for displaying Somatic Mutation Rate in `pretty_lollipop_plot`.
#' @param Sample2 A label for displaying Somatic Mutation Rate in `pretty_lollipop_plot`.
#' 
#' @return A lollipop plot.
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#' 
#' #get meta data (BL_Thomas)
#' metadata <- get_gambl_metadata() %>%
#'     filter(cohort == "BL_Thomas")
#' 
#' maf_df <- get_ssm_by_samples(
#'     these_samples_metadata = metadata
#' )
#'
#' #construct pretty_lollipop_plot.
#' lolipop_result <- pretty_lollipop_plot(maf_df, "MYC")
#'
pretty_lollipop_plot <- function(
    maf_df, 
    gene,
    plot_title,
    refseq_id,
    include_silent = FALSE,
    plotarg = TRUE,
    mirrorarg = FALSE,
    combined_gene_counts,
    meta1_counter,
    meta2_counter,
    Sample1 = Sample1,
    Sample2 = Sample2
) {
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
        arrange(AA)

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
        mutate(mutation_count = n()) 

    if (mirrorarg == TRUE){
        gene_counts <- combined_gene_counts
    } else {
        gene_counts <- gene_counts
    }

    # protein_domains a bundled object with GAMBLR.data
    protein_domain_subset <- subset(
        protein_domains, 
        HGNC == gene
    )
   
    if (missing(refseq_id)){
        if (length(unique(protein_domain_subset$refseq.ID)) == 1){ # Gene only has 1 NM_XXX on protein_domains$refseq.ID
            protein_domain_subset <- protein_domain_subset
        } else if (length(unique(protein_domains$refseq.ID)) > 1){ # Gene has >1 NM_XXX on protein_domains$refseq.ID, the first one is selected
            protein_domain_subset <- protein_domain_subset %>%
                filter(refseq.ID == protein_domain_subset$refseq.ID[1])
        }
    } else if (refseq_id %in% protein_domain_subset$refseq.ID){ # refseq_id matches to an NM_XXX on protein_domains$refseq.ID
        protein_domain_subset <- protein_domain_subset %>%
            filter(refseq.ID == refseq_id)
    } else { # refseq_id has no matches to protein_domains.Rd
        stop(paste("Error: refseq_id", refseq_id, "is not found in protein_domains", 
            "\nThe following refseq IDs are available in protein_domains.RD:",
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

        Somatic_Mutation_Denominator <- length(unique(maf_df$Tumor_Sample_Barcode))

        Somatic_Mutation_Rate <- Somatic_Mutation_Numerator/Somatic_Mutation_Denominator *100
        Somatic_Mutation_Rate <- round(Somatic_Mutation_Rate, 2)
    }

    # KS-Test 
    if(mirrorarg == TRUE){
        # Gene KS-Test
        aa_frequency_data <- combined_gene_counts %>%
            group_by(AA, source) %>%
            summarise(freq = n())

        ks_test_result <- ks.test(
            aa_frequency_data$freq[aa_frequency_data$source == "Sample 1"], 
            aa_frequency_data$freq[aa_frequency_data$source == "Sample 2"]
        )

        gene_p_value <- ks_test_result$p.value

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

                # Mutation counts for Sample 1 and Sample 2
                domain_mutation1 <- domain_subset %>%
                    filter(
                        Variant_Classification %in% variants, 
                        source == "Sample 1"
                    ) %>%
                    nrow()

                domain_mutation2 <- domain_subset %>%
                    filter(
                        Variant_Classification %in% variants, 
                        source == "Sample 2"
                    ) %>%
                    nrow()

                domain_frequency_data <- domain_subset %>%
                    group_by(
                        AA, 
                        source
                    ) %>%
                    summarise(freq = n())

                domain_ks_test <- tryCatch({
                    ks.test(
                        domain_frequency_data$freq[domain_frequency_data$source == "Sample 1"],
                        domain_frequency_data$freq[domain_frequency_data$source == "Sample 2"]
                    )
                }, error = function(e) {
                    list(p.value = NA)
                })

                # Save domain-specific results in a list
                unique_domain_name <- paste0(
                    "Domain_", 
                    domain_name
                )
                domain_list[[unique_domain_name]] <- list(
                    mutation_count_sample1 = domain_mutation1,
                    mutation_count_sample2 = domain_mutation2,
                    ks_p_value = domain_ks_test$p.value
                )

                # Store the p-value in domain_data
                domain_data$p_value[i] <- domain_ks_test$p.value

                # Store the unique domain name
                domain_names <- c(domain_names, unique_domain_name)
            } else {
                print(paste(
                    "No mutation data for domain", 
                    domain_name
                ))
            }
        }
    }

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
                size = abs(mutation_count)
                )
        ) +
        # Background rectangle for regions without domain data
        geom_rect(
            aes(
                xmin = x_min, 
                xmax = x_max, 
                ymin = -0.2, 
                ymax = 0.2
                ), 
            fill = "black", 
            color = "black"
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
        geom_text(
            data = domain_data, 
            aes(
                x = text.position, 
                y = 0.1, 
                label = text.label
                )
        ) 

        if (mirrorarg == TRUE) {
            plot <- plot +
                annotate(
                    "text",
                    x = 0,
                    y = max(gene_counts$mutation_count) * 1.1,
                    label = paste0(
                        "Somatic Mutation Rate ", 
                        Somatic_Mutation_Rate[1], 
                        "% N = ", 
                        Somatic_Mutation_Denominator[1],
                        "\n",
                        "Comparison Value ", 
                        '"', 
                        Sample1, 
                        '"'
                    ),
                    hjust = 0
                ) +
                annotate(
                    "text",
                    x = 0,
                    y = min(gene_counts$mutation_count) * 1.1,
                    label = paste0(
                        "Somatic Mutation Rate ", 
                        Somatic_Mutation_Rate[2], 
                        "% N = ",
                        Somatic_Mutation_Denominator[2], 
                        "\n", 
                        "Comparison Value ",
                        '"', 
                        Sample2, 
                        '"'
                    ),
                    hjust = 0
                ) +
                geom_text(
                    data = domain_data, 
                    aes(
                        x = text.position, 
                        y = -0.2, 
                        label = paste0(
                            "p = ", 
                            round(p_value, 3)
                        )  # Display domain-specific p-value
                    )
                ) +
                labs(
                    x = "AA Position", 
                    y = "Mutation Count", 
                    title = paste0(
                        plot_title,
                        "\n",
                        "p = ",
                        round(gene_p_value, 3)
                    )
                )
        } else {
            plot <- plot +
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
                )
        }

        plot <- plot +
        scale_size_continuous(
            name = "Mutation Count", 
            labels = function(x) abs(x)  
        ) + 
        theme_bw() +
        theme(
            plot.title = element_text(
                hjust = 0.5
                ),
            axis.text.x = element_text(
              angle = 90, 
              hjust = 1
              )
        ) +
        scale_color_manual(
          name = "Legend", 
          values = colours_manual
        ) 

    if (plotarg == TRUE) {
        return(plot)
    } else {
        return(gene_counts)
    }

}
