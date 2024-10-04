#' @title Mirrored Lollipop Plot.
#'
#' @description Generates a visually appealing mirrored lollipop plot.
#'
#' @details Retrieve two maf data files of a specific sample or a set of samples for comparison. 
#' A gene of interest can then be visualized with the given maf data files. Silent mutations can 
#' be visualized setting include_silent to TRUE. 
#'
#' @param maf_df1 A data frame containing the mutation data from a given cohort or pathology, etc.
#' @param maf_df2 A data frame containing the mutation data from an other cohort or pathology, etc.
#' @param gene The gene symbol to plot.
#' @param plot_title Optional, the title of the plot. Default is gene.
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is FALSE. 
#' 
#' @return A mirrored lollipop plot.
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#' 
#' metadata <- get_gambl_metadata()
#' metadata1 <- metadata %>%
#'      filter(cohort == "DLBCL_Hilton")
#' metadata2 <- metadata %>%
#'      filter(cohort == "DLBCL_Thomas")
#' 
#' maf_df1 <- get_ssm_by_samples(
#'      these_samples_metadata = metadata1
#' )
#'
#' maf_df2 <- get_ssm_by_samples(
#'      these_samples_metadata = metadata2
#' )
#' 
#' #construct mirrored_lollipop_plot
#' mirrored_lollipop_result <- mirrored_lollipop_plot(maf_df1, maf_df2, "IGLL5")
#' 

mirrored_lollipop_plot <- function(
    maf_df1,
    maf_df2, 
    gene = NULL,
    plot_title,
    include_silent = FALSE
) {
    if(missing(gene)){
        stop("Please provide a gene...")
    }
  
    if(missing(plot_title)){
        plot_title=gene
    }

    maf_df2 <- as.data.frame(maf_df2)
    maf_df1 <- as.data.frame(maf_df1)

    # Specifying noncoding regions with coding_class as a bundled object with GAMBLR.helpers
    if(include_silent){
        variants <- coding_class
    } else {
        variants <- coding_class[!coding_class %in% c(
            "Silent",
            "Splice_Region"
        )]
    }

    nc_maf_df2 <- maf_df2 %>%
        filter(
            Hugo_Symbol == gene
        ) %>% 
        filter(
            Variant_Classification %in% variants 
        )

    nc_maf_df1 <- maf_df1 %>%
        filter(
            Hugo_Symbol == gene
        ) %>% 
        filter(
            Variant_Classification %in% variants 
        ) 

    # Filter for the specific gene 
    gene_df2 <- nc_maf_df2 %>% 
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

    gene_counts2 <- gene_df2 %>% 
        group_by(
            AA, 
            Start_Position, 
            End_Position, 
            Variant_Classification, 
            Reference_Allele, 
            Tumor_Seq_Allele2
        ) %>%
        arrange(AA) %>%
        summarise(mutation_count = n())

    gene_df1 <- nc_maf_df1 %>% 
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

    gene_counts1 <- gene_df1 %>% 
        group_by(
            AA, 
            Start_Position, 
            End_Position, 
            Variant_Classification, 
            Reference_Allele, 
            Tumor_Seq_Allele2
        ) %>%
        arrange(AA) %>%
        summarise(mutation_count = n())     

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
    x_max <- max(
        max(
            domain_data$end.points
        ), 
        max(
            gene_counts1$AA, 
            na.rm = TRUE
        ),
        max(
            gene_counts2$AA, 
            na.rm = TRUE
        )
    )
    x_min <- 0

    # get_gambl_colours() from GAMBLR.helpers
    colours_manual <- get_gambl_colours("mutation")
 
    # Somatic mutation statistic
    Somatic_Mutation_Numerator2 <- maf_df2 %>%
        filter(
            Hugo_Symbol == gene
        ) %>% 
        filter(
            Variant_Classification %in% variants 
        ) %>% 
        distinct(Tumor_Sample_Barcode) %>% 
        nrow()

    Somatic_Mutation_Denominator2 <- length(unique(maf_df2$Tumor_Sample_Barcode))

    Somatic_Mutation_Rate2 <- Somatic_Mutation_Numerator2/Somatic_Mutation_Denominator2 *100
    Somatic_Mutation_Rate2 <- round(Somatic_Mutation_Rate2, 2)

    Somatic_Mutation_Numerator1 <- maf_df1 %>%
        filter(
            Hugo_Symbol == gene
        ) %>% 
        filter(
            Variant_Classification %in% variants 
        ) %>% 
        distinct(Tumor_Sample_Barcode) %>% 
        nrow()

    Somatic_Mutation_Denominator1 <- length(unique(maf_df1$Tumor_Sample_Barcode))

    Somatic_Mutation_Rate1 <- Somatic_Mutation_Numerator1/Somatic_Mutation_Denominator1 *100
    Somatic_Mutation_Rate1 <- round(Somatic_Mutation_Rate1, 2)

    plot <- ggplot() +
        geom_segment(
            data = gene_counts2, 
            aes(
                x = AA, 
                xend = AA, 
                y = 0, 
                yend = -mutation_count
                )
        ) +
        geom_segment(
            data = gene_counts1, 
            aes(
                x = AA, 
                xend = AA, 
                y = 0, 
                yend = mutation_count
                )
        ) +
        geom_point(
            data = gene_counts2, 
            aes(
                x = AA, 
                y = -mutation_count, 
                color = Variant_Classification, 
                size = mutation_count
                )
        ) + 
        annotate(
                 "text",
                 x = 0,
                 y = max(gene_counts2$mutation_count) * -1.1,
                 label = paste0("Sample 2: ", Somatic_Mutation_Rate2, "%"),
                 hjust = 0
                ) +
        geom_point(
            data = gene_counts1, 
            aes(
                x = AA, 
                y = mutation_count, 
                color = Variant_Classification, 
                size = mutation_count
                )
        ) + 
        annotate(
                 "text",
                 x = 0,
                 y = max(gene_counts1$mutation_count) * 1.1,
                 label = paste0("Sample 1: ", Somatic_Mutation_Rate1, "%"),
                 hjust = 0
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
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(
                hjust = 0.5
                ),
            axis.text.x = element_text(
              angle = 45, 
              hjust = 1
              )
        ) +
        scale_color_manual(
          name = "Legend", 
          values = colours_manual
        )

  return(plot)
}
