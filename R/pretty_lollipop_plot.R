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
#' @param by_allele Set to FALSE to consider all mutations at the same codon as equivalent.
#' When FALSE, and combined with labelPos, the labels will only indicate the amino acid number. Default is TRUE.
#' @param labelPos A vector of amino acid positions that should be labeled on the plot (i.e. c("28", "315")). If unset, label_threshold is used to determine which positions to label.
#' @param labelHGVSp A character vector of values from HGVSp_Short (i.e. c("p.G28D", "p.E315K")) that will be labeled on the plot.
#' @param label_threshold Minimum mutation recurrence for labels to appear on plot. Superseded by either labelPos or labelHGVSp. Default set to 5.
#' @param show_rate Boolean parameter indicating whether to show the somatic mutation rate in the plot title. Default is FALSE.,
#' @param max_count Sets the mutation count for the largest lollipop size. Default is 10.
#' @param include_silent Logical parameter indicating whether to include silent mutations into coding mutations. Default is FALSE.
#' @param font_size Customizable font size for the CoLollipop plot somatic mutation rate and comparison labels. Default is 11pt font.
#' @param font Customizable font for all plotting elements. Default is "sans" which uses the system default font. Use `systemfonts::match_fonts(font)` to confirm if your requested font is available. 
#' @param title_size Font size for the plot title. Default is 14.
#' @param x_axis_size Font size for the x-axis labels. Default is 10.
#' @param domain_label_size Font size for the protein domain labels. Default is 4.
#' @param aa_label_size Font size for the amino acid labels. Default is 3.
#' @param point_alpha Alpha to apply to the lollipop points. Default is 1 (no transparency).
#' @param point_size_range Vector of length 2 specifiying the size range of points on the lollipops. Default is c(2,8). 
#'
#' @return A list of plot and data objects.
#'
#' @import dplyr ggplot2 ggrepel ggpubr GAMBLR.data
#' @export
#'
#' @examples
#' library(GAMBLR.open)
#' suppressMessages(
#'   suppressWarnings({
#' 
#' #get meta data (BL_Thomas)
#' metadata <- get_gambl_metadata() %>%
#'     filter(seq_type == "genome") %>%
#'     check_and_clean_metadata(.,duplicate_action="keep_first")
#'
#' maf_df <- get_coding_ssm(
#'     these_samples_metadata = metadata
#' )
#'
#' #construct pretty_lollipop_plot.
#' lolipop_result <- pretty_lollipop_plot(maf_df, these_samples_metadata = metadata, "DDX3X")
#' print(lolipop_result)
#' }))
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
    refseq_id = NULL,
    by_allele = TRUE,
    max_count = 10,
    include_silent = FALSE,
    label_threshold = 5,
    labelPos = NULL, 
    labelHGVSp = NULL,
    mirror = FALSE, 
    font_size = 11,
    font = "sans",
    show_rate = FALSE,
    title_size = 14,
    x_axis_size = 10,
    domain_label_size = 4,
    aa_label_size = 3, 
    point_alpha = 1, 
    point_size_range = c(2, 8)
){
    ##### Input checks #####
    if(missing(gene)){
        stop("Please provide a gene...")
    }

    if(missing(plot_title)){
        plot_title = gene
    }
    
    if (!"RefSeq" %in% colnames(maf_df)) {
        stop("Error: The provided maf_df is missing the 'RefSeq' column. 
              Ensure that get_ssm_by_samples() is run with basic_columns = FALSE, 
              or that you have updated GAMBLR.data to the latest version. .")
    }
    
    if(!is.null(labelPos) & !is.null(labelHGVSp)){
        stop("Please provide either labelPos or labelHGVSp, not both.")
    }
    
    if(!is.null(labelHGVSp) & !by_allele){
        stop("Please set by_allele = TRUE when using labelHGVSp.
        When by_allele = FALSE, all alleles at the same position are counted together and HGVSp information is lost. ")
    }

    
    ##### Subset the maf to gene and variants of interest #####
    
    if(include_silent){
        variants <- coding_class
    } else {
        variants <- coding_class[!coding_class %in% c(
            "Silent",
            "Splice_Region" # Is this useful? Aren't these all intronic? Won't affect cds
        )]
    }
    
    maf_df <- as.data.frame(maf_df) %>%
        filter(
            Hugo_Symbol == gene
        ) %>%
        filter(
            Variant_Classification %in% variants 
        ) %>% 
        # Separate RefSeq column if it is comma-delimited
        separate_longer_delim(
            RefSeq, 
            delim = ","
            ) %>% 
        mutate(
            RefSeq = str_remove(RefSeq, "[.].*"), 
            aa.length = str_remove(Protein_position, ".*/")
        )
        
    # Verify that all variants belong to a single gene model
    if(length(unique(maf_df$aa.length)) > 1){
        maf_df %>% 
            select(
                Hugo_Symbol, 
                Transcript_ID,
                RefSeq, 
                aa.length
            ) %>% 
            distinct() %>% 
            print()
        stop("The provided maf_df has more than one gene model for the specified gene. 
            Check the Transcript_ID and Protein_position columns to ensure one consistent gene model was used. ")
    }

    ##### Load the protein domain table for the current gene #####
    
    if(!is.null(refseq_id)){
        # User-specified refseq_id
        refseq_id = str_remove(refseq_id, "[.].*") # Remove version number if present
        protein_domain_subset <- protein_domains %>% 
            filter(HGNC == gene) %>% 
            filter(refseq.ID == refseq_id)
        if(length(unique(protein_domain_subset$refseq.ID)) == 0){
            stop(paste("The provided refseq_id", refseq_id, "does not match any RefSeq IDs for the specified gene."))
        }
    } else {
        # Infer refseq_id to use from the maf_df
        # Use the length of the gene model from the maf to identify the correct RefSeq ID
        protein_domain_subset <- protein_domains %>% 
            filter(HGNC == gene) %>% 
            filter(refseq.ID %in% unique(maf_df$RefSeq)) %>%
            filter(aa.length == unique(maf_df$aa.length))
    }
    
    # Verify that this reduces down to a single gene model
    if(length(unique(protein_domain_subset$refseq.ID)) > 1){
        print(protein_domain_subset)
        protein_domain_subset <- protein_domain_subset %>%
            filter(refseq.ID == unique(protein_domain_subset$refseq.ID[1]))
        message(paste("There is more than one RefSeq model matching the maf_df for the specified gene.
        Arbitrarily selecting", unique(protein_domain_subset$refseq.ID), "to plot."))
    } else if (length(unique(protein_domain_subset$refseq.ID)) == 0) {
       print(protein_domain_subset)
       stop("None of the protein models matches the provided maf_df. Check the Protein_position and RefSeq columns. ")
    }

    ##### Count mutations according to user specified options #####
    gene_df <- maf_df %>%
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
        
    if(by_allele){
        # Keep different HGVSp values even if AA position is the same
        gene_counts <- gene_df %>%
            group_by(
                HGVSp_Short,
                AA,
                AA_position,
                Variant_Classification
            ) %>%
            arrange(AA) %>%
            summarise(mutation_count = n()) %>%
            ungroup() %>% 
            mutate(size = ifelse(
                mutation_count > max_count, 
                max_count, 
                mutation_count
                ))     
    }else{
        # Collapse HGVSp values if AA is the same
        gene_counts <- gene_df %>%
            group_by(
                AA, 
                Variant_Classification
            ) %>%
            arrange(AA) %>%
            summarise(mutation_count = n()) %>% 
            mutate(size = ifelse(
                mutation_count > max_count, 
                max_count, 
                mutation_count
                )) 
    }
    
    if(!is.null(labelHGVSp)){
        # Only label user-specified HGVSp values
        gene_counts <- gene_counts %>% 
            mutate(label = ifelse(
                HGVSp_Short %in% labelHGVSp, 
                AA_position, 
                ""
            ), 
            size = ifelse(mutation_count > max_count, max_count, mutation_count)) 

    } else if(!is.null(labelPos)){
        # Only label user-specified AA positions
        gene_counts <- gene_counts %>% 
            mutate(label = ifelse(
                AA %in% labelPos, 
                as.character(AA), 
                ""
            ), 
            size = ifelse(mutation_count > max_count, max_count, mutation_count)) 

    } else {
        # Only label variants that exceed label_threshold recurrence
        gene_counts <- gene_counts %>% 
            mutate(label = ifelse(
                mutation_count >= label_threshold, 
                str_remove(HGVSp_Short, "p[.]"), 
                ""
            ), 
            size = ifelse(mutation_count > max_count, max_count, mutation_count)) 

    }

    ##### Generate the final domain_data object and make the domain_plot #####
    domain_data <- protein_domain_subset %>% 
        data.frame(
            start.points = protein_domain_subset$Start,
            end.points = protein_domain_subset$End,
            text.label = protein_domain_subset$Label,
            color = protein_domain_subset$Label
        )

    domain_data$text.position <- (domain_data$start.points + domain_data$end.points) / 2

    domain_plot <- draw_domain_plot(
        domain_data, 
        font = font, 
        domain_label_size = domain_label_size,
        x_axis_size = x_axis_size
        )

    ##### Generate the final mutation_plot object #####
    mutation_plot <- draw_mutation_plot(
        gene_counts = gene_counts, 
        plot_title = plot_title,
        protein_length = domain_data$aa.length[1],
        colours_manual = get_gambl_colours("mutation"), 
        font = font, 
        mirror = mirror, 
        aa_label_size = aa_label_size,
        title_size = title_size, 
        max_count = max_count, 
        point_size_range = point_size_range,
        point_alpha = point_alpha
    )
    
    if(show_rate){
        denominator <- length(unique(these_samples_metadata$sample_id))
        numerator <- length(unique(gene_df[gene_df$RefSeq == protein_domain_subset$refseq.ID[1],]$Tumor_Sample_Barcode))
        somatic_mutation_rate <- round((numerator / denominator) * 100, 1)
        
        mutation_plot <- mutation_plot + 
            ggtitle(
                plot_title, 
                subtitle = paste0(
                    "Somatic mutation rate: ", 
                    somatic_mutation_rate, 
                    "%"
                )
            ) + 
            theme(
                text = element_text(family = font), 
                plot.subtitle = element_text(hjust = 0.5)
            )
    }
    ##### Combine domain and mutation plots #####
    combined_plot <- ggpubr::ggarrange(
        mutation_plot,
        domain_plot, 
        ncol = 1, 
        align = "v", 
        heights = c(3, 1), 
        common.legend = TRUE, 
        legend = "right"
    )

    to_return <- list(
        plot = combined_plot,
        domain_plot = domain_plot, 
        mutation_plot = mutation_plot, 
        gene_counts = gene_counts, 
        domain_data = domain_data
    )
    return(to_return)
}


draw_domain_plot <- function(
    domain_data, 
    font = "sans", 
    x_axis_size = 10, 
    domain_label_size = 4
    ){
    domain_palette <- unname(get_gambl_colours("domains"))[1:length(unique(domain_data$text.label))]
    
    # Determine angle of domain labels based on width of domains relative to protein
    angle = ifelse(
        min((domain_data$End - domain_data$Start)/domain_data$aa.length) >= 0.2,
        0, 90
    )
    
    # Add the p-value to the text label if it exists in the provided domain_data
    if("p_value" %in% colnames(domain_data)){
        domain_data$text.label <- ifelse(
            domain_data$p_value != "", 
            str_c(
                domain_data$text.label,
                "\n",
                domain_data$p_value
            ), 
            domain_data$text.label
        )
    } 
    
    # Fix labels that are too close together
    # TODO: Fix it to retain significant p-values if they exist
    domain_data <- domain_data %>% 
        arrange(text.position) %>% 
        mutate(
            text.label = ifelse(
                text.position - lag(text.position) >= 10 | is.na(lag(text.position)), 
                text.label, 
                ""
            )
        )
    # Specify plot start and end points
    x_min <- 1 # Set 1 as the start for any protein (obviously)
    x_max = unique(domain_data$aa.length)[1] # Assuming aa.length is a single value for the gene
    domain_plot <- ggplot() + 
            geom_rect(
                aes(
                    xmin = x_min, 
                    xmax = x_max, 
                    ymin = -1, 
                    ymax = 1
                    ), 
                fill = "lightgrey", 
                color = "lightgrey"
            ) + 
            geom_rect(
                data = domain_data, 
                aes(
                    xmin = start.points, 
                    xmax = end.points, 
                    ymin = -1.2, 
                    ymax = 1.2, 
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
                    angle = angle, 
                    label = text.label
                ), 
                hjust = 0.5, 
                vjust = 0.5,
                size = domain_label_size, 
                family = font
            ) +
            scale_fill_manual(values = domain_palette) +
            theme_void() +
            theme(
                axis.text.x = element_text(
                    size = x_axis_size, 
                    vjust = rel(0.5)
                    ), 
                text = element_text(family = font)
            ) + 
            xlim(x_min, x_max)
    return(domain_plot)
}

draw_mutation_plot <- function(
    gene_counts, 
    plot_title, 
    protein_length, 
    font = "sans", 
    colours_manual, 
    mirror = FALSE,
    max_count, 
    point_size_range,
    title_size,
    aa_label_size, 
    point_alpha
    ){
    
    nudge_factor <- ifelse(mirror, -0.5, 0.5)
    mutation_plot <- ggplot() +
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
                size = size
                ), 
            alpha = point_alpha
        ) +
        ggrepel::geom_text_repel(data = gene_counts, aes(
                x = AA, 
                y = mutation_count, 
                label = label, 
                point.size = size
            ), 
            min.segment.length = 0, 
            nudge_y = nudge_factor,
            nudge_x = 0,
            direction = "y", 
            size = aa_label_size, 
            family = font
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
                max(gene_counts$mutation_count) * 1.2,
                by = ifelse(max(gene_counts$mutation_count) * 1.2 <= 5, 1, round(max(gene_counts$mutation_count) * 1.2 / 5))
            ), 
            limits = c(0, max(gene_counts$mutation_count) * 1.2)
        )  +
        xlim(1, protein_length) +
        scale_size_continuous(
            limits = c(1, max_count),
            breaks = seq(
                1, 
                max_count,
                by = round(max_count / 5)
            ), 
            range = point_size_range,
            guide = "none"  
        )  +
        ggpubr::theme_pubr() +
        theme(
            plot.title = element_text(
                hjust = 0.5, 
                size = title_size
                ),
            axis.text.x = element_blank(),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            text = element_text(family = font), 
            legend.position = "right"
        ) +
        scale_color_manual(
            name = "Variant\nClassification", 
            values = colours_manual
        )+
    guides(color = guide_legend(override.aes = list(size = 3))) 
    return(mutation_plot)
}
