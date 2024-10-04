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
    maf_df = NULL, 
    gene = NULL,
    plot_title,
    include_silent = FALSE,
    plotarg = TRUE,
    mirrorarg = FALSE,
    combined_gene_counts = NULL,
    meta1_counter = NULL,
    meta2_counter = NULL,
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
                y = 0, 
                label = text.label
                )
        ) +
        labs(
            x = "AA Position", 
            y = "Mutation Count", 
            title = paste0(
                plot_title
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
                        "%", 
                        "\n", 
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
                        "%", 
                        "\n", 
                        '"', 
                        Sample2, 
                        '"'
                    ),
                    hjust = 0
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
