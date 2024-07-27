#' @title Lollipop Plot
#'
#' @description Generates a visually appealing and interactive lollipop plot.
#'
#' @details Retrieve a maf of a specific sample or a set of samples. Returned plot is interactive, meaning the 
#' user can hover over individual points in the plot to reveal more information. The plot can also be exported 
#' in a variety of file formats inside the interactive view of the lollipop plot.
#'
#' @param maf_data A data frame containing the mutation data (from a MAF).
#' @param gene_symbol The gene symbol to plot.
#' 
#' @return A lolipop plot.
#'
#' @import tidyverse BiocManager maftools
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#' 
#' #get metadata (BL_Thomas)
#' metadata <- get_gambl_metadata() %>%
#'     filter(cohort == "BL_Thomas")
#' 
#' maf_data <- get_ssm_by_samples(
#'     these_samples_metadata = metadata
#' )
#'
#' #construct pretty_lollipop_plot.
#' lolipop_result <- get_lolipop_by_gene_symbol(maf_data, "MYC")
#' print(lolipop_result)
#'
get_lolipop_by_gene_symbol <- function(maf_data, 
                                       gene_symbol = NULL) {
  # Read the MAF object if it's not already in MAF format
  if (!inherits(maf_data, "MAF")) {
    maf_data <- read.maf(maf_data)
  }
  
  maf_data <- maf_data@data 
  # Filter for the specific gene 
  gene_df <- maf_data %>% 
    filter(Hugo_Symbol == gene_symbol) %>%
    mutate(AA = as.numeric(gsub("[^0-9]+", "", gsub("([0-9]+).*", "\\1", HGVSp_Short)))) %>%
    arrange(AA)

  gene_counts <- gene_df %>% 
    group_by(AA, Start_Position, 
             End_Position, 
             Variant_Classification, 
             Reference_Allele, 
             Tumor_Seq_Allele2) %>%
    arrange(AA) %>%
    summarise(mutation_count = n()) 

  protein_domain_subset <- subset(protein_domains, HGNC == gene_symbol)
  if (!is.data.frame(protein_domain_subset)) {
    protein_domain_subset <- as.data.frame(protein_domain_subset)
  }
   
  domain_data <- protein_domain_subset %>% data.frame(start.points = protein_domain_subset$Start,
                                                      end.points = protein_domain_subset$End,
                                                      text.label = protein_domain_subset$Label,
                                                      color = protein_domain_subset$Label)

  domain_data$text.position <- (domain_data$start.points + domain_data$end.points) / 2

  # Determine the x-axis range
  x_max <- max(max(domain_data$end.points), 
               max(gene_counts$AA, na.rm = TRUE))
  x_min <- 0

  colors_manual <- c(setNames(protein_domain_subset$Label, protein_domain_subset$Label), 
                     "Frame_Shift_Del" = "#e50505", 
                     "Frame_Shift_Ins" = "#df9307" ,
                     "In_Frame_Del" = "#f6f60c",
                     "In_Frame_Ins" = "#04e504",
                     "Missense_Mutation" = "#9e12f6", 
                     "Nonsense_Mutation" = "#06ebeb",
                     "Nonstop_Mutation" = "#0505df", 
                     "Splice_Site" = "#035703",
                     "Translation_Start_Site" = "#f3899b",
                     "RNA" = "#146deb")
  
  Somatic_Mutation_Statistic_Numerator <- maf_data %>%
    filter(Hugo_Symbol == gene_symbol) %>% 
    filter(!(Variant_Classification %in% c("Intron",
                                          "Silent",
                                          "3'UTR",
                                          "5'UTR",
                                          "IGR",
                                          "3'Flank",
                                          "5'Flank",
                                          "Splice_Region"))) %>% 
    distinct(Tumor_Sample_Barcode) %>% 
    nrow()
  
  Somatic_Mutation_Statistic_Denominator <- length(unique(maf_data$Tumor_Sample_Barcode))

  Somatic_Mutation_Rate <- Somatic_Mutation_Statistic_Numerator/Somatic_Mutation_Statistic_Denominator *100
  Somatic_Mutation_Rate <- round(Somatic_Mutation_Rate, 2)

  plot <- ggplot() +
          geom_segment(data = gene_counts, aes(x = AA, xend = AA, y = 0, yend = mutation_count)) +
          geom_point(data = gene_counts, aes(x = AA, y = mutation_count, color = Variant_Classification, size = mutation_count)) +
          # Background rectangle for regions without domain data
          geom_rect(aes(xmin = x_min, xmax = x_max, ymin = -0.2, ymax = 0.2), fill = "black", color = "black") + 
          # Domain rectangles
          geom_rect(data = domain_data, aes(xmin = start.points, xmax = end.points, ymin = -0.4, ymax = 0.4, fill = color), color = "black", show.legend = FALSE) +
          geom_text(data = domain_data, aes(x = text.position, y = 0, label = text.label)) +
          labs(x = "AA Position", y = "Mutation Count", title = paste0(gene_symbol, " :[Somatic Mutation Rate: ", Somatic_Mutation_Rate, "%]")) +
          theme_bw() +
          theme(
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1)
          ) +
          scale_color_manual(name = "Legend", values = colors_manual)

  return(plot)
}
