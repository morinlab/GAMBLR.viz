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
#' @param by_allele Set to FALSE to consider all mutations at the same codon as equivalent.
#' When FALSE, and combined with labelPos, the labels will only indicate the amino acid number. Default is TRUE.
#' @param labelPos Specify which AA positions to label in the plot (default is no labels).
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is FALSE.
#'
#' @return A lollipop plot.
#'
#' @import dplyr ggplot2 GAMBLR.data
#' @export
#'
#' @examples
#' library(GAMBLR.open)

#' #get meta data (BL_Thomas)
#' metadata <- suppressMessages(get_gambl_metadata()) %>%
#'     filter(cohort == "BL_Thomas")
#'
#' maf_df <- get_coding_ssm(
#'     these_samples_metadata = metadata
#' )
#'
#' #construct pretty_lollipop_plot.
#' lolipop_result <- pretty_lollipop_plot(maf_df, "DDX3X")
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
#' lolipop_result <- pretty_lollipop_plot(maf_df, "DDX3X")
#' lolipop_result
#'}
pretty_lollipop_plot <- function(
    maf_df,
    gene = NULL,
    plot_title,
    by_allele = TRUE,
    max_count = 10,
    include_silent = FALSE,
    labelPos = NULL,
    show_rate = FALSE,
    title_size = 8,
    x_axis_size = 4,
    domain_label_size = 0,
    aa_label_size = 3
) {
    if(missing(gene)){
        stop("Please provide a gene...")
    }

    if(missing(plot_title)){
        plot_title=gene
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
    if(by_allele){
      gene_counts <- gene_df %>%
        group_by(
          HGVSp_Short,
          AA,
          Variant_Classification
        ) %>%
        arrange(AA) %>%
        summarise(mutation_count = n()) %>%
        ungroup()
      gene_counts = mutate(gene_counts,label=HGVSp_Short,size=ifelse(mutation_count>max_count,max_count,mutation_count))
    }else{
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
        summarise(mutation_count = n()) %>%
        ungroup()
      gene_counts = mutate(gene_counts,label=AA,size=ifelse(mutation_count>max_count,max_count,mutation_count))
    }
    max_mutation_count = max(gene_counts$mutation_count)
    #select(gene_counts,AA,HGVSp_Short,mutation_count) %>% arrange(desc(mutation_count)) %>% print()
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
        # Background rectangle for regions without domain data
        geom_rect(
            aes(
                xmin = x_min ,
                xmax = x_max
                
                ),           
                ymin = -0.01, #* max_mutation_count,
                ymax = 0.01, #* max_mutation_count,     

            fill = "black",
            color = "black"
        ) +
        # Domain rectangles
        geom_rect(
            data = domain_data,
            aes(
                xmin = start.points,
                xmax = end.points,
                
                fill = color
                ),
            ymin = -0.02 * max_mutation_count,
            ymax = 0.02 * max_mutation_count,
            color = "black",
            show.legend = FALSE
        ) +
        geom_point(
            data = gene_counts,
            colour = "black",
            shape = 21,   
            aes(
                x = AA,
                y = mutation_count,
                fill = Variant_Classification, 
                size = size,    
                )
        ) 
        if(domain_label_size > 0){
        plot = plot +  
          geom_text(
            data = domain_data,
            aes(
                x = text.position,
                y = 0,
                label = text.label,
                size = domain_label_size
                )
        )
        }
        

        if(!is.null(labelPos)){
          dplyr::filter(gene_counts, AA %in% labelPos) %>%
            unique() %>% print()
          maxAA = dplyr::filter(gene_counts, AA %in% labelPos) %>%
            select(AA,label,mutation_count) %>% pull(AA) %>% max()
          gene_counts = mutate(gene_counts,label=gsub("p.","",label))
          plot = plot + geom_text_repel(data = dplyr::filter(gene_counts, AA %in% labelPos) %>%
                                          select(AA,label,mutation_count) %>%
                                          unique(),
                                   aes(
                                     x = AA,
                                     y = mutation_count,
                                     label = label
                                     
                                   ),
                                   size = aa_label_size,
                                   nudge_y = max_mutation_count * 0.2,
                                   nudge_x = maxAA/3,
                                   max.overlaps=2)
        }
        if(show_rate){
            plot_title = paste0(
                plot_title,
                "\n[Somatic Mutation Rate: ",
                Somatic_Mutation_Rate,
                "%]"
                )
        }
        max_y = max_mutation_count + max_mutation_count * 0.2
        #if(max_mutation_count > max_count){
        #    max_y = max_y + max_count/2
        #}else {
        #    max_y = max_y + max_mutation_count/2
        #}
        plot = plot +
        labs(
            x = "AA Position",
            y = "Mutation Count",
            title = plot_title
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(
                hjust = 0.5,
                size = title_size,
                ),
            axis.text.x = element_text(
              angle = 45,
              hjust = 1
              ),
            axis.title = element_text(
              size = x_axis_size
            )
        ) +
        scale_fill_manual(
          name = "Legend",
          values = colours_manual
        ) + 
        ylim(c(0,max_y))
  #print(colours_manual)
  return(plot)
}
