#' @title genome-wide ideogram annotated with SSM and CN information
#'
#' @description Generate sample-level ideogram with copy number information, ssm and gene annotations, etc.
#'
#' @details This function generates genome-wide ideograms, visualizing SSM data as well as CN segments.
#' It is also possible to superimpose the plot with gene annotations. Offering a comprehensive overview of all SSM and CN segments of different aneuploidy.
#' The plotting of SSM can be toggled with setting `include_ssm` to TRUE. If so, it is also possible to count the number of SSMs per chromosome with `ssm_count = TRUE`.
#' To get data for plotting, there are a few different options available; like all `fanncy_x_plots` a sample ID can be provided to the `this_sample_id`
#' parameter. If done so, the function will retrieve data (SSm and CN segments) by wrapping the appropriate functions.
#' This data can also be provided with `seg_data`, `seg_path`, `maf_data` and `maf_path`.
#' For more info on how to run with these parameters, refer to the parameter descriptions.
#' In order to annotate the ideogram with genes, simply give the `gene_annotations` parameter a set of genes as a vector of characters or a data frame with gene names in the first column.
#' Another useful parameter for restricting the plotted regions is to call the function with `intersect_regions`.
#' This parameter takes a vector of characters or a data frame with regions that the plotted calls are restricted to.
#'
#' @param this_sample_id Sample to be plotted (for multiple samples, see fancy_multisample_ideogram.
#' @param gene_annotation Annotate ideogram with a set of genes. These genes can either be specified as a vector of characters or a data frame.
#' @param seg_data Optional parameter with copy number df already loaded into R.
#' @param seg_path Optional parameter with path to external seg like file.
#' @param maf_data Optional parameter with maf like df already loaded into R.
#' @param maf_path Optional parameter with path to external maf like file.
#' @param variant_type_col_maf Index of column holding Variant Type (to be used with either maf_data or maf_path).
#' @param chromosome_col_maf Index of column holding Chromosome (to be used with either maf_data or maf_path).
#' @param start_col_maf Index of column with variant start coordinates (to be used with either maf_data or maf_path).
#' @param end_col_maf Index of column with variant end coordinates (to be used with either maf_data or maf_path).
#' @param chrom_col_seg Index of column with chromosome annotations (to be used with either maf_data or maf_path).
#' @param start_col_seg Index of column with copy number start coordinates (to be used with either maf_data or maf_path).
#' @param end_col_seg Index of column with copy number end coordinates (to be used with either maf_data or maf_path).
#' @param cn_col_seg Index of column holding copy number information (to be used with either maf_data or maf_path).
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Optional argument for plot subtitle.
#' @param intersect_regions Optional parameter for subset variant calls to specific regions. Should be either a vector of characters (chr:start-end) or data frame with regions.
#' @param include_ssm Set to TRUE to plot SSMs (dels and ins).
#' @param ssm_count Optional parameter to summarize n variants per chromosome, inlcude_ssm must be set to TRUE.
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import ggplot2 dplyr cowplot
#' @export
#'
#' @examples
#' #build plot
#' \dontrun{
#' fancy_ideogram(this_sample_id = "DOHH-2",
#'                gene_annotation = "MYC",
#'                plot_title = "Sample-level Ideogram Example",
#'                plot_subtitle = "grch37")
#' }
#'
fancy_ideogram = function(this_sample_id,
                          gene_annotation,
                          seg_data,
                          seg_path = NULL,
                          maf_data,
                          maf_path = NULL,
                          variant_type_col_maf = 10,
                          chromosome_col_maf = 5,
                          start_col_maf = 6,
                          end_col_maf = 7,
                          chrom_col_seg = 2,
                          start_col_seg = 3,
                          end_col_seg = 4,
                          cn_col_seg = 7,
                          plot_title = paste0(this_sample_id),
                          plot_subtitle = "Genome-wide Ideogram (grch37).",
                          intersect_regions,
                          include_ssm = TRUE,
                          ssm_count = TRUE,
                          coding_only = FALSE,
                          this_seq_type = "genome"){

  #plot theme
  ideogram_theme = function(){
    theme(legend.position = "bottom", axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank())
  }

  #grch37 coordinates
  grch37_end = dplyr::filter(GAMBLR.data::chromosome_arms_grch37, arm == "q") %>% 
    dplyr::filter(chromosome != "Y") %>% 
    pull(end)
  grch37_cent_start = dplyr::filter(GAMBLR.data::chromosome_arms_grch37, arm == "p") %>% 
    dplyr::filter(chromosome != "Y") %>% 
    pull(end)
  grch37_cent_end = dplyr::filter(GAMBLR.data::chromosome_arms_grch37, arm == "q") %>% 
    dplyr::filter(chromosome != "Y") %>% 
    pull(start)
  
  #additional regions to plot
  if(!missing(gene_annotation)){
    gene = GAMBLR.utils::gene_to_region(gene_symbol = gene_annotation, genome_build = "grch37", return_as = "df")
    gene.annotate = gene[gene$chr %in% paste0(c(1:22)), ]
    cols.int = c("chromosome", "start", "end")
    gene.annotate[cols.int] = sapply(gene.annotate[cols.int], as.integer)
  }

  #build chr table for segment plotting
  chr = c( paste0("chr", c(1:22)), "chrX" )
  chr_start = c(0)
  chr_end = grch37_end
  cent_start = grch37_cent_start
  cent_end =  grch37_cent_end
  y = c(1:23)
  yend = c(1:23)

  #transform to data frame
  segment_data = data.frame(chr, chr_start, chr_end, cent_start, cent_end, y, yend)

  #load CN data
  if(!missing(seg_data)){
    cn_states = seg_data
    cn_states = as.data.frame(cn_states)
    colnames(cn_states)[chrom_col_seg] = "chrom"
    colnames(cn_states)[start_col_seg] = "start"
    colnames(cn_states)[end_col_seg] = "end"
    colnames(cn_states)[cn_col_seg] = "CN"

  }else if(!is.null(seg_path)){
    cn_states = read.table(seg_path, sep = "\t", header = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
    cn_states = as.data.frame(cn_states)
    colnames(cn_states)[chrom_col_seg] = "chrom"
    colnames(cn_states)[start_col_seg] = "start"
    colnames(cn_states)[end_col_seg] = "end"
    colnames(cn_states)[cn_col_seg] = "CN"
  }

  #get maf data for a specific sample.
  if(missing(seg_data) && is.null(seg_path)){
    cn_states = get_sample_cn_segments(
      these_sample_ids = this_sample_id,
      with_chr_prefix = FALSE,
      streamlined = FALSE,
      this_seq_type = this_seq_type
    )
  }

  # ignore y chromosome
  cn_states = dplyr::filter(cn_states, chrom != "Y")
  
  #convert chr into y coordinates
  cn_states$ycoord = cn_states$chrom

  #paste chr in chromosomecolumn, if not there
  if(!str_detect(cn_states$chrom[1], "chr")){
    cn_states = mutate(cn_states, chrom = paste0("chr", chrom))
  }

  if(!missing(intersect_regions)){
    #filter CN states on intersecting regions
    #transform regions to data tables
    #convenience function for converting intersect regions to a df (if it's a string) and renaming columns to match required format.
    if(is.list(intersect_regions)){
      colnames(intersect_regions)[1] = "chrom"
      colnames(intersect_regions)[2] = "start"
      colnames(intersect_regions)[3] = "end"
      if(!str_detect(intersect_regions$chrom, "chr")){
        intersect_regions = mutate(intersect_regions, chrom = paste0("chr", chrom))
      }
      intersect_regions = as.data.table(intersect_regions)
      intersect_regions$start = as.numeric(intersect_regions$start)
      intersect_regions$end = as.numeric(intersect_regions$end)
    }

    if(is.character(intersect_regions)){
      if(length(intersect_regions) > 1){
        message("Please only enter one region, only first region will be regarded. For mutiple regions, kindly provide a data frame with regions of interest")
      }

      split_chunks = unlist(strsplit(intersect_regions, ":"))
      split_chunks = unlist(strsplit(split_chunks, "-"))
      chrom = split_chunks[1]
      start = split_chunks[2]
      end = split_chunks[3]
      intersect_regions = cbind(chrom, start, end) %>%
        as.data.frame()

      intersect_regions$start = as.numeric(intersect_regions$start)
      intersect_regions$end = as.numeric(intersect_regions$end)

      if(!str_detect(intersect_regions$chrom[1], "chr")){
        intersect_regions = mutate(intersect_regions, chrom = paste0("chr", chrom))
      }
    }

    incoming_cn = as.data.table(cn_states)
    regions_sub = as.data.table(intersect_regions)

    #set keys
    data.table::setkey(incoming_cn, chrom, start, end)
    data.table::setkey(regions_sub, chrom, start, end)

    #intersect regions
    intersect = data.table::foverlaps(regions_sub, incoming_cn, nomatch = 0)

    #transform object to data frame
    inter_df = as.data.frame(intersect)

    #organize columns to match the expected format
    cn_states = select(inter_df, ID, chrom, start, end, LOH_flag, log.ratio, CN, ycoord)
  }

  #convert data types
  cn_states$chrom = as.factor(cn_states$chrom)
  
  # correct y coordinate for X chromosome
  cn_states$ycoord = dplyr::recode(cn_states$ycoord, X="23") %>% 
    as.integer

  #subset on CN state
  cn_states$CN[cn_states$CN > 6] = 6
  cn_states = subset(cn_states, CN != 2)
  cn_states$CN = paste0("cn_", cn_states$CN)
  cn_states$CN = as.factor(cn_states$CN)
  cn_states = droplevels(cn_states)
  l = split(cn_states, cn_states$CN)
  list2env(l, envir = .GlobalEnv)

  #load maf data
  if(include_ssm){
    if(!missing(maf_data)){
      maf = maf_data
      maf = as.data.frame(maf)
      colnames(maf)[variant_type_col_maf] = "Variant_Type"
      colnames(maf)[chromosome_col_maf] = "Chromosome"
      colnames(maf)[start_col_maf] = "Start_Position"
      colnames(maf)[end_col_maf] = "End_Position"
    }else if (!is.null(maf_path)){
      maf = fread_maf(maf_path)
      maf = as.data.frame(maf)
      colnames(maf)[variant_type_col_maf] = "Variant_Type"
      colnames(maf)[chromosome_col_maf] = "Chromosome"
      colnames(maf)[start_col_maf] = "Start_Position"
      colnames(maf)[end_col_maf] = "End_Position"
    }

    if(missing(maf_data) && is.null(maf_path)){
      maf = get_ssm_by_sample(this_sample_id = this_sample_id, 
                              projection = "grch37", #currently only supported reference build.
                              this_seq_type = this_seq_type)
    }

    #transform maf data
    maf_trans = dplyr::select(maf, Chromosome, Start_Position, End_Position, Variant_Type)

    #calculate mid points of variants
    maf_trans$mid = ((maf_trans$End_Position - maf_trans$Start_Position) / 2) + maf_trans$Start_Position

    #convert chr into y coordinates
    maf_trans$ystart = maf_trans$yend = dplyr::recode(maf_trans$Chromosome, X="23")
    
    #paste chr in maf, if not there
    if(!str_detect(maf_trans$Chromosome[1], "chr")){
      maf_trans = mutate(maf_trans, Chromosome = paste0("chr", Chromosome))
    }

    #convert data types
    maf_trans$Start_Position = as.double(maf_trans$Start_Position)
    maf_trans$End_Position = as.double(maf_trans$End_Position)
    maf_trans$ystart = as.integer(maf_trans$ystart)
    maf_trans$yend = as.integer(maf_trans$yend)

    if(!missing(intersect_regions)){
      #filter CN states on intersecting regions
      maf_tmp = maf_trans
      colnames(maf_tmp)[1] = "chrom"
      colnames(maf_tmp)[2] = "start"
      colnames(maf_tmp)[3] = "end"
      maf_tmp = dplyr::select(maf_tmp, chrom, start, end)
      maf.table = as.data.table(maf_tmp)
      data.table::setkey(maf.table, chrom, start, end)

      #intersect regions
      intersect_maf = data.table::foverlaps(regions_sub, maf.table, nomatch = 0)

      #transform object to data frame
      inter_maf_df = as.data.frame(intersect_maf)

      #rename columns
      colnames(inter_maf_df)[1] = "Chromosome"
      colnames(inter_maf_df)[2] = "Start_Position"
      colnames(inter_maf_df)[3] = "End_Position"

      #subset
      inter_maf_df = dplyr::select(inter_maf_df, Chromosome, Start_Position, End_Position)

      #perform a semi join with all cn states (to retain necessary columns)
      maf_trans = dplyr::semi_join(maf_trans, inter_maf_df)
    }

    #subset on variant type
    if(nrow(maf_trans) > 0){
      maf_del = dplyr::filter(maf_trans, Variant_Type == "DEL")
      maf_ins = dplyr::filter(maf_trans, Variant_Type == "INS")

      if(ssm_count){
        del_count = dplyr::filter(maf_trans, Variant_Type == "DEL") %>%
          add_count(Chromosome) %>%
          distinct(Chromosome, .keep_all = TRUE) %>%
          dplyr::select(Chromosome, Variant_Type, yend, n)

        ins_count = dplyr::filter(maf_trans, Variant_Type == "INS") %>%
          add_count(Chromosome) %>%
          distinct(Chromosome, .keep_all = TRUE) %>%
          dplyr::select(Chromosome, Variant_Type, yend, n)
      }
    }
  }

  #get colours and combine palette for indels and cn states
  ideogram_palette = c(GAMBLR.helpers::get_gambl_colours("indels"), GAMBLR.helpers::get_gambl_colours("copy_number"))
  selected_colours = ideogram_palette[c(1,2,19,18,16:13)]
  names(selected_colours)[c(3:8)] = c("CN0", "CN1", "CN3", "CN4", "CN5", "CN6+")

  #plot
  p = ggplot() +
    {if(include_ssm && nrow(maf_del) > 0) geom_segment(data = maf_del, aes(x = mid - 100000, xend = mid + 100000, y = ystart - 0.27, yend = yend - 0.27), color = "#53B1FC", size = 5, stat = "identity", position = position_dodge(width = 0))} + #del
    {if(include_ssm && nrow(maf_del) > 0) geom_segment(data = maf_del, aes(x = mid, xend = mid, y = ystart - 0.35, yend = yend - 0.35), color = "black", lineend = "round", size = 3.5, stat = "identity", position = position_dodge(width = 0))} + #del
    {if(include_ssm && nrow(maf_del) > 0) geom_segment(data = maf_del, aes(x = mid, xend = mid, y = ystart - 0.35, yend = yend - 0.35, color = "DEL"), lineend = "round", size = 3, stat = "identity", position = position_dodge(width = 0))} + #del
    {if(include_ssm && nrow(maf_ins) > 0) geom_segment(data = maf_ins, aes(x = mid - 100000, xend = mid + 100000, y = ystart - 0.27, yend = yend - 0.27), color = "#FC9C6D", size = 5, stat = "identity", position = position_dodge(width = 0))} + #ins
    {if(include_ssm && nrow(maf_ins) > 0) geom_segment(data = maf_ins, aes(x = mid, xend = mid, y = ystart - 0.35, yend = yend - 0.35), color = "black", lineend = "round", size = 3.5, stat = "identity", position = position_dodge(width = 0))} + #ins
    {if(include_ssm && nrow(maf_ins) > 0) geom_segment(data = maf_ins, aes(x = mid, xend = mid, y = ystart - 0.35, yend = yend - 0.35, color = "INS"), lineend = "round", size = 3, stat = "identity", position = position_dodge(width = 0))} + #ins
    {if(ssm_count && nrow(maf_del) > 0) annotate(geom = "text", x = -4000000, y = del_count$yend, label = del_count$n, color = "#3A8799", size = 3)} + #count del
    {if(ssm_count && nrow(maf_trans) > 0) annotate(geom = "text", x = -2300000, y = segment_data$y, label = " | ", color = "black", size = 3)} + #count sep
    {if(ssm_count && nrow(maf_ins) > 0) annotate(geom = "text", x = -1000000, y = ins_count$yend, label = ins_count$n, color = "#E6856F", size = 3)} + #count ins
    geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y, yend = yend), color = "#99A1A6", lineend = "butt", size = 5, stat = "identity", position = position_dodge(width = 0)) + #chr contigs
    geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y, yend = yend), color = "white", size = 6, stat = "identity", position = position_dodge(width = 0)) + #centromeres
    {if("cn_0" %in% levels(cn_states$CN)) geom_segment(data = cn_0, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN0"), size = 4.7, stat = "identity", position = position_dodge(width = 0))} + #cn3
    {if("cn_1" %in% levels(cn_states$CN)) geom_segment(data = cn_1, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN1"), size = 4.7, stat = "identity", position = position_dodge(width = 0))} + #cn3
    {if("cn_3" %in% levels(cn_states$CN)) geom_segment(data = cn_3, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN3"), size = 4.7, stat = "identity", position = position_dodge(width = 0))} + #cn3
    {if("cn_4" %in% levels(cn_states$CN)) geom_segment(data = cn_4, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN4"), size = 4.7, stat = "identity", position = position_dodge(width = 0))} + #cn4
    {if("cn_5" %in% levels(cn_states$CN)) geom_segment(data = cn_5, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN5"), size = 4.7, stat = "identity", position = position_dodge(width = 0))} + #cn5
    {if("cn_6" %in% levels(cn_states$CN)) geom_segment(data = cn_6, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN6+"), size = 4.7, stat = "identity", position = position_dodge(width = 0))} + #cn6 and more
    {if(!missing(gene_annotation)) geom_point(data = gene.annotate, aes(x = ((end - start) / 2) + start, y = chromosome - 0.28), shape = 25, color = "#A63932", fill = "#A63932", stat = "identity", position = position_dodge(width = 0))} + #gene annotation
    {if(!missing(gene_annotation)) geom_label(data = gene.annotate, aes((x = end - start) / 2 + start, y = chromosome - 0.52, label = hugo_symbol), fontface = "bold", color = "white", fill = "#A63932", size = 3)} + #gene annotation text
    geom_text(aes(x = -10000000 , y = yend, label = segment_data$chr), color = "black", size = 5) + #chr labels
    labs(title = plot_title, subtitle = plot_subtitle) + #plot titles
    scale_colour_manual(name = "", values = selected_colours) + #legend/colours
    scale_x_continuous(breaks = seq(0, max(segment_data$chr_end), by = 30000000)) + #set x-axis boundaries
    scale_y_reverse() + #reverse y axis
    theme_cowplot() + #themes
    ideogram_theme() #themes

  return(p)
}
