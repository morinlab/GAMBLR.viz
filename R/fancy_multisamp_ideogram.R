#' @title Genome-wide ideogram (CN segments) for multiple samples.
#'
#' @description Generate ideograms for selected samples, visualizing copy number variation segments.
#'
#' @details To create multi-sample ideograms, i.e showing CN segments across multiple samples, this function was created.
#' This can be used to infer inheritance patterns, hotspots, etc. across multiple samples or the sample ID but for different timepoints.
#' In addition, this plot can also allow to only plot concordant (or discordant) cn segments between two samples. i.e how two samples differ, or are alike.
#' The function automatically detects the number of samples provided and sets the plotting parameters accordingly.
#' The maximum number of samples this plot can deal with is 4 and sample IDs of interest are given to the `these_sample_ids`.
#' In order to only plot segments that are concordant between the selected samples, set `komapre` to TRUE and `concordance` to TRUE.
#' To instead plot discordant CN segments, set this parameter to FALSE.
#'
#' @param these_sample_ids Sample to be plotted (accepts 2, 3 or 4 samples).
#' @param plot_title Main title of plot.
#' @param plot_sub Subtitle of plot.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param chr_anno_dist Optional parameter to adjust chromosome annotations, default value is 3, increase to adjust annotations to left.
#' @param chr_select Optional parameter to subset plot to specific chromosomes. Default value is chr1-22.
#' @param include_cn2 Set to TRUE for plotting CN states == 2.
#' @param kompare Boolean statement, set to TRUE to call cnvKompare on the selected samples for plotting concordant (or discordant) cn segments across selected chromosomes.
#' @param concordance Boolean parameter to be used when kompare = TRUE. Default is TRUE, to plot discordant segments, set the parameter to FALSE.
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#'
#' @return A plot as a ggplot object (grob).
#'
#' @import ggplot2 dplyr stringr GAMBLR.helpers
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#'
#' #two samples ideogram
#' fancy_multisamp_ideogram(these_sample_ids = c("05-18426T", "05-18426T"),
#'                          include_cn2 = TRUE,
#'                          plot_title = "Multi-sample Ideograms Example",
#'                          plot_sub = "grch37",
#'                          chr_anno_dist = 2.5,
#'                          chr_select = paste0("chr", c(1:22)))
#'
fancy_multisamp_ideogram = function(these_sample_ids,
                                    plot_title = "CN Segments Ideogram",
                                    plot_sub = "grch37",
                                    this_seq_type = "genome",
                                    chr_anno_dist = 3,
                                    chr_select = paste0("chr", c(1:22)),
                                    include_cn2 = FALSE,
                                    kompare = FALSE,
                                    concordance = TRUE,
                                    coding_only = FALSE){

  #plot theme
  ideogram_theme = function(){
    theme(legend.position = "bottom", axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank())}

  #transform chr_anno_dist and set some other variables
  anno_dist = chr_anno_dist * -10000000
  if(length(these_sample_ids) == 2){
    seg_dist = 0.18
    seg_size = 4
    seg_size_cent = 5
  }else if(length(these_sample_ids) == 3){
    seg_dist = 0.28
    seg_size = 3.5
    seg_size_cent = 4
  }else if(length(these_sample_ids) == 4){
    seg_dist = 0.12}

  #chr segment coordinates
  grch37_end = GAMBLR.data::chromosome_arms_grch37[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44),3]
  grch37_cent_start = GAMBLR.data::chromosome_arms_grch37[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43),3]
  grch37_cent_end = GAMBLR.data::chromosome_arms_grch37[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44),2]

  #build chr table for segment plotting
  chr = paste0("chr", c(1:22))
  chr_start = c(0)
  chr_end = grch37_end
  cent_start = grch37_cent_start
  cent_end =  grch37_cent_end
  y = c(1:22)
  yend = c(1:22)

  #transform to data frame
  segment_data = data.frame(chr, chr_start, chr_end, cent_start, cent_end, y, yend)
  segment_data$chr = as.factor(segment_data$chr)

  #sub-setting maf based on user-defined parameters
  segment_data = segment_data[segment_data$chr %in% chr_select, ]
  segment_data = droplevels(segment_data)

  if(kompare){
    #call cnvKompare to retreive CN segments shared (or not) shared between selected samples.
    cnv_komp = GAMBLR.utils::cnvKompare(these_sample_ids = these_sample_ids)

    #select concordant or discordant CN segments for plotting.
    if(concordance){
      cnv_cord = cnv_komp$concordant_cytobands
    }else{
      cnv_cord = cnv_komp$discordant_cytobands
    }

    #transform log.ratio to CN states
    cnv_cord$CN_tmp = 2*2^cnv_cord$log.ratio
    cnv_cord$CN = round(cnv_cord$CN_tmp) %>%
      as.factor()

    colnames(cnv_cord)[2] = "chrom"
    colnames(cnv_cord)[3] = "start"
    colnames(cnv_cord)[4] = "end"

    cn_states = dplyr::select(cnv_cord, ID, chrom, start, end, CN)
  }else{
    #load CN data
    cn_states = get_sample_cn_segments(
      these_sample_ids = these_sample_ids,
      streamlined = FALSE,
      this_seq_type = this_seq_type
    )
  }

  #convert chr into y coordinates
  cn_states$ycoord = cn_states$chrom

  #paste chr in chromosome column, if not there
  if(!str_detect(cn_states$chrom[1], "chr")){
    cn_states = mutate(cn_states, chrom = paste0("chr", chrom))
    }

  #transform data types
  cols.int = c("start", "end", "ycoord")
  cn_states[cols.int] = sapply(cn_states[cols.int], as.integer)
  cn_states$chrom = as.factor(cn_states$chrom)
  cn_states$ID = as.factor(cn_states$ID)

  #sub-setting maf based on user-defined parameters
  cn_states = cn_states[cn_states$chrom %in% chr_select, ]
  cn_states = droplevels(cn_states)

  #retrieve sample names
  samples = levels(cn_states$ID)

  #first sample
  sample1 = samples[1]
  sample1_cn = dplyr::filter(cn_states, ID == sample1)
  GAMBLR.helpers::subset_cnstates(cn_segments = sample1_cn, samplen = 1, include_2 = include_cn2)
  sample1_cn$CN = as.factor(sample1_cn$CN)
  sample1_cn = droplevels(sample1_cn)

  #second sample
  sample2 = samples[2]
  sample2_cn = dplyr::filter(cn_states, ID == sample2)
  GAMBLR.helpers::subset_cnstates(cn_segments = sample2_cn, samplen = 2, include_2 = include_cn2)
  sample2_cn$CN = as.factor(sample2_cn$CN)
  sample2_cn = droplevels(sample2_cn)

  #third sample (if provided...)
  if(length(these_sample_ids) > 2){
    sample3 = samples[3]
    sample3_cn = dplyr::filter(cn_states, ID == sample3)
    GAMBLR.helpers::subset_cnstates(cn_segments = sample3_cn, samplen = 3, include_2 = include_cn2)
    sample3_cn$CN = as.factor(sample3_cn$CN)
    sample3_cn = droplevels(sample3_cn)}

  #fourth sample (if provided...)
  if(length(these_sample_ids) > 3){
    sample4 = samples[4]
    sample4_cn = dplyr::filter(cn_states, ID == sample4)
    GAMBLR.helpers::subset_cnstates(cn_segments = sample4_cn, samplen = 4, include_2 = include_cn2)
    sample4_cn$CN = as.factor(sample4_cn$CN)
    sample4_cn = droplevels(sample4_cn)}

  #get colours and combine palette for indels and cn states
  ideogram_palette = GAMBLR.helpers::get_gambl_colours("copy_number")
  if(include_cn2){
    selected_colours = ideogram_palette[c(17:11)]
    names(selected_colours)[c(1:6)] = c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6+")
  }else{
    selected_colours = ideogram_palette[c(17,16,14:11)]
    names(selected_colours)[c(1:6)] = c("CN0", "CN1", "CN3", "CN4", "CN5", "CN6+")
  }

  #plot
  if(length(these_sample_ids) >= 2 & length(these_sample_ids) < 4){
    p = ggplot() + geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y + seg_dist, yend = yend + seg_dist, label = chr), color = "#99A1A6", lineend = "butt", size = seg_size, stat = "identity", position = position_dodge()) + #chr contigs
      geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y + seg_dist, yend = yend + seg_dist), color = "white", size = seg_size_cent, stat = "identity", position = position_dodge()) + #centromeres
      annotate(geom = "text", x = -2000000, y = segment_data$yend + seg_dist, label = sample1, color = "black", size = 3, hjust = 1) + #sample name annotations
      {if("0" %in% levels(sample1_cn$CN)) geom_segment(data = cn_0_sample1, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN0"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn0
      {if("1" %in% levels(sample1_cn$CN)) geom_segment(data = cn_1_sample1, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN1"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn1
      {if("2" %in% levels(sample1_cn$CN)) geom_segment(data = cn_2_sample1, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN2"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn2
      {if("3" %in% levels(sample1_cn$CN)) geom_segment(data = cn_3_sample1, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN3"), size = seg_size, stat = "identity", position = position_dodge())} + #cn3
      {if("4" %in% levels(sample1_cn$CN)) geom_segment(data = cn_4_sample1, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN4"), size = seg_size, stat = "identity", position = position_dodge())} + #cn4
      {if("5" %in% levels(sample1_cn$CN)) geom_segment(data = cn_5_sample1, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN5"), size = seg_size, stat = "identity", position = position_dodge())} + #cn5
      {if("6+" %in% levels(sample1_cn$CN)) geom_segment(data = cn_6_sample1, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN6+"), size = seg_size, stat = "identity", position = position_dodge())} + #second sample
      geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y - seg_dist, yend = yend - seg_dist, label = chr), color = "#99A1A6", lineend = "butt", size = seg_size, stat = "identity", position = position_dodge()) + #chr contigs
      geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y - seg_dist, yend = yend - seg_dist), color = "white", size = seg_size_cent, stat = "identity", position = position_dodge()) + #centromeres
      annotate(geom = "text", x = -2000000, y = segment_data$yend - seg_dist, label = sample2, color = "black", size = 3, hjust = 1) + #sample name annotations
      {if("0" %in% levels(sample2_cn$CN)) geom_segment(data = cn_0_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN0"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn0
      {if("1" %in% levels(sample2_cn$CN)) geom_segment(data = cn_1_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN1"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn1
      {if("2" %in% levels(sample2_cn$CN)) geom_segment(data = cn_2_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN2"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn2
      {if("3" %in% levels(sample2_cn$CN)) geom_segment(data = cn_3_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN3"), size = seg_size, stat = "identity", position = position_dodge())} + #cn3
      {if("4" %in% levels(sample2_cn$CN)) geom_segment(data = cn_4_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN4"), size = seg_size, stat = "identity", position = position_dodge())} + #cn4
      {if("5" %in% levels(sample2_cn$CN)) geom_segment(data = cn_5_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN5"), size = seg_size, stat = "identity", position = position_dodge())} + #cn5
      {if("6+" %in% levels(sample2_cn$CN)) geom_segment(data = cn_6_sample2, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN6+"), size = seg_size, stat = "identity", position = position_dodge())} #third sample
    if(length(these_sample_ids) > 2){
      p = p + geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y + seg_dist - 0.28, yend = yend + seg_dist - 0.28, label = chr), color = "#99A1A6", lineend = "butt", size = seg_size, stat = "identity", position = position_dodge()) + #chr contigs
        geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y + seg_dist - 0.28, yend = yend + seg_dist - 0.28), color = "white", size = seg_size_cent, stat = "identity", position = position_dodge()) + #centromeres
        annotate(geom = "text", x = -2000000, y = segment_data$yend + seg_dist - 0.28, label = sample3, color = "black", size = 3, hjust = 1) + #sample name annotations
        {if("0" %in% levels(sample3_cn$CN)) geom_segment(data = cn_0_sample3, aes(x = start, xend = end, y = ycoord + seg_dist - 0.28, yend = ycoord + seg_dist - 0.28, color = "CN0"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn0
        {if("1" %in% levels(sample3_cn$CN)) geom_segment(data = cn_1_sample3, aes(x = start, xend = end, y = ycoord + seg_dist - 0.28, yend = ycoord + seg_dist - 0.28, color = "CN1"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn1
        {if("2" %in% levels(sample3_cn$CN)) geom_segment(data = cn_2_sample3, aes(x = start, xend = end, y = ycoord + seg_dist - 0.28, yend = ycoord + seg_dist - 0.28, color = "CN2"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn2
        {if("3" %in% levels(sample3_cn$CN)) geom_segment(data = cn_3_sample3, aes(x = start, xend = end, y = ycoord + seg_dist - 0.28, yend = ycoord + seg_dist - 0.28, color = "CN3"), size = seg_size, stat = "identity", position = position_dodge())} + #cn3
        {if("4" %in% levels(sample3_cn$CN)) geom_segment(data = cn_4_sample3, aes(x = start, xend = end, y = ycoord + seg_dist - 0.28, yend = ycoord + seg_dist - 0.28, color = "CN4"), size = seg_size, stat = "identity", position = position_dodge())} + #cn4
        {if("5" %in% levels(sample3_cn$CN)) geom_segment(data = cn_5_sample3, aes(x = start, xend = end, y = ycoord + seg_dist - 0.28, yend = ycoord + seg_dist - 0.28, color = "CN5"), size = seg_size, stat = "identity", position = position_dodge())} + #cn5
        {if("6+" %in% levels(sample3_cn$CN)) geom_segment(data = cn_6_sample3, aes(x = start, xend = end, y = ycoord + seg_dist - 0.28, yend = ycoord + seg_dist - 0.28, color = "CN6+"), size = seg_size, stat = "identity", position = position_dodge())}} #cn6+
    p = p + geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y - 0.5, yend = yend - 0.5, label = chr), color = "white", lineend = "butt", size = 2, stat = "identity", position = position_dodge()) + #white space between chromosome groups (upper)
      geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y + 0.5, yend = yend + 0.5, label = chr), color = "white", lineend = "butt", size = 2, stat = "identity", position = position_dodge()) + #white space between chromosome groups (bottom)
      annotate(geom = "text", x = anno_dist, y = segment_data$yend, label = segment_data$chr, color = "black", size = 5, hjust = 1) + #chr labels
      labs(title = plot_title, subtitle = plot_sub) + #plot titles
      scale_colour_manual(name = "", values = selected_colours) + #colours and legend
      scale_x_continuous(breaks = seq(0, max(segment_data$chr_end), by = 30000000)) + #x-axis boundaries
      scale_y_reverse() + #flip ideogram
      theme_Morons()

    #plotting has its own plotting chunk, due to re-scaling of geom_segments and segment widths etc.
  }else if(length(these_sample_ids) == 4){
    p = ggplot() + geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y - (seg_dist - 0.48), yend = yend - (seg_dist - 0.48), label = chr), color = "#99A1A6", lineend = "butt", size = 3, stat = "identity", position = position_dodge()) + #chr contigs
      geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y - (seg_dist - 0.48), yend = yend - (seg_dist - 0.48)), color = "white", size = 3, stat = "identity", position = position_dodge()) + #centromeres
      annotate(geom = "text", x = -2000000, y = segment_data$yend - (seg_dist - 0.48), label = sample1, color = "black", size = 3, hjust = 1) + #sample name annotations
      {if("0" %in% levels(sample1_cn$CN)) geom_segment(data = cn_0_sample1, aes(x = start, xend = end, y = ycoord - (seg_dist - 0.48), yend = ycoord - (seg_dist - 0.48), color = "CN0"), size = 3, stat = "identity", position = position_dodge())} +  #cn0
      {if("1" %in% levels(sample1_cn$CN)) geom_segment(data = cn_1_sample1, aes(x = start, xend = end, y = ycoord - (seg_dist - 0.48), yend = ycoord - (seg_dist - 0.48), color = "CN1"), size = 3, stat = "identity", position = position_dodge())} +  #cn1
      {if("2" %in% levels(sample1_cn$CN)) geom_segment(data = cn_2_sample1, aes(x = start, xend = end, y = ycoord - (seg_dist - 0.48), yend = ycoord - (seg_dist - 0.48), color = "CN2"), size = 3, stat = "identity", position = position_dodge())} +  #cn2
      {if("3" %in% levels(sample1_cn$CN)) geom_segment(data = cn_3_sample1, aes(x = start, xend = end, y = ycoord - (seg_dist - 0.48), yend = ycoord - (seg_dist - 0.48), color = "CN3"), size = 3, stat = "identity", position = position_dodge())} + #cn3
      {if("4" %in% levels(sample1_cn$CN)) geom_segment(data = cn_4_sample1, aes(x = start, xend = end, y = ycoord - (seg_dist - 0.48), yend = ycoord - (seg_dist - 0.48), color = "CN4"), size = 3, stat = "identity", position = position_dodge())} + #cn4
      {if("5" %in% levels(sample1_cn$CN)) geom_segment(data = cn_5_sample1, aes(x = start, xend = end, y = ycoord - (seg_dist - 0.48), yend = ycoord - (seg_dist - 0.48), color = "CN5"), size = 3, stat = "identity", position = position_dodge())} + #cn5
      {if("6+" %in% levels(sample1_cn$CN)) geom_segment(data = cn_6_sample1, aes(x = start, xend = end, y = ycoord - (seg_dist - 0.48), yend = ycoord - (seg_dist - 0.48), color = "CN6+"), size = 3, stat = "identity", position = position_dodge())} + #cn6+
      geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y - seg_dist, yend = yend - seg_dist, label = chr), color = "#99A1A6", lineend = "butt", size = 3, stat = "identity", position = position_dodge()) + #chr contigs
      geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y - seg_dist, yend = yend - seg_dist), color = "white", size = 3, stat = "identity", position = position_dodge()) + #centromeres
      annotate(geom = "text", x = -2000000, y = segment_data$yend - seg_dist, label = sample2, color = "black", size = 3, hjust = 1) + #sample name annotations
      {if("0" %in% levels(sample2_cn$CN)) geom_segment(data = cn_0_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN0"), size = 3, stat = "identity", position = position_dodge())} +  #cn0
      {if("1" %in% levels(sample2_cn$CN)) geom_segment(data = cn_1_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN1"), size = 3, stat = "identity", position = position_dodge())} +  #cn1
      {if("2" %in% levels(sample2_cn$CN)) geom_segment(data = cn_2_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN2"), size = 3, stat = "identity", position = position_dodge())} +  #cn2
      {if("3" %in% levels(sample2_cn$CN)) geom_segment(data = cn_3_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN3"), size = 3, stat = "identity", position = position_dodge())} + #cn3
      {if("4" %in% levels(sample2_cn$CN)) geom_segment(data = cn_4_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN4"), size = 3, stat = "identity", position = position_dodge())} + #cn4
      {if("5" %in% levels(sample2_cn$CN)) geom_segment(data = cn_5_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN5"), size = 3, stat = "identity", position = position_dodge())} + #cn5
      {if("6+" %in% levels(sample2_cn$CN)) geom_segment(data = cn_6_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN6+"), size = 3, stat = "identity", position = position_dodge())} + #cn6+
      geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y + seg_dist, yend = yend + seg_dist, label = chr), color = "#99A1A6", lineend = "butt", size = 3, stat = "identity", position = position_dodge()) + #chr contigs
      geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y + seg_dist, yend = yend + seg_dist), color = "white", size = 3, stat = "identity", position = position_dodge()) + #centromeres
      annotate(geom = "text", x = -2000000, y = segment_data$yend + seg_dist, label = sample3, color = "black", size = 3, hjust = 1) + #sample name annotations
      {if("0" %in% levels(sample3_cn$CN)) geom_segment(data = cn_0_sample3, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN0"), size = 3, stat = "identity", position = position_dodge())} +  #cn0
      {if("1" %in% levels(sample3_cn$CN)) geom_segment(data = cn_1_sample3, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN1"), size = 3, stat = "identity", position = position_dodge())} +  #cn1
      {if("2" %in% levels(sample3_cn$CN)) geom_segment(data = cn_2_sample3, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN2"), size = 3, stat = "identity", position = position_dodge())} +  #cn2
      {if("3" %in% levels(sample3_cn$CN)) geom_segment(data = cn_3_sample3, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN3"), size = 3, stat = "identity", position = position_dodge())} + #cn3
      {if("4" %in% levels(sample3_cn$CN)) geom_segment(data = cn_4_sample3, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN4"), size = 3, stat = "identity", position = position_dodge())} + #cn4
      {if("5" %in% levels(sample3_cn$CN)) geom_segment(data = cn_5_sample3, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN5"), size = 3, stat = "identity", position = position_dodge())} + #cn5
      {if("6+" %in% levels(sample3_cn$CN)) geom_segment(data = cn_6_sample3, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN6+"), size = 3, stat = "identity", position = position_dodge())} + #cn6+
      geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y + (seg_dist - 0.48), yend = yend + (seg_dist - 0.48), label = chr), color = "#99A1A6", lineend = "butt", size = 3, stat = "identity", position = position_dodge()) + #chr contigs
      geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y + (seg_dist - 0.48), yend = yend  + (seg_dist - 0.48)), color = "white", size = 3, stat = "identity", position = position_dodge()) + #centromeres
      annotate(geom = "text", x = -2000000, y = segment_data$yend  + (seg_dist - 0.48), label = sample4, color = "black", size = 3, hjust = 1) + #sample name annotations
      {if("0" %in% levels(sample4_cn$CN)) geom_segment(data = cn_0_sample4, aes(x = start, xend = end, y = ycoord + (seg_dist - 0.48), yend = ycoord + (seg_dist - 0.48), color = "CN0"), size = 3, stat = "identity", position = position_dodge())} +  #cn0
      {if("1" %in% levels(sample4_cn$CN)) geom_segment(data = cn_1_sample4, aes(x = start, xend = end, y = ycoord + (seg_dist - 0.48), yend = ycoord + (seg_dist - 0.48), color = "CN1"), size = 3, stat = "identity", position = position_dodge())} +  #cn1
      {if("2" %in% levels(sample4_cn$CN)) geom_segment(data = cn_2_sample4, aes(x = start, xend = end, y = ycoord + (seg_dist - 0.48), yend = ycoord + (seg_dist - 0.48), color = "CN2"), size = 3, stat = "identity", position = position_dodge())} +  #cn2
      {if("3" %in% levels(sample4_cn$CN)) geom_segment(data = cn_3_sample4, aes(x = start, xend = end, y = ycoord + (seg_dist - 0.48), yend = ycoord + (seg_dist - 0.48), color = "CN3"), size = 3, stat = "identity", position = position_dodge())} + #cn3
      {if("4" %in% levels(sample4_cn$CN)) geom_segment(data = cn_4_sample4, aes(x = start, xend = end, y = ycoord + (seg_dist - 0.48), yend = ycoord + (seg_dist - 0.48), color = "CN4"), size = 3, stat = "identity", position = position_dodge())} + #cn4
      {if("5" %in% levels(sample4_cn$CN)) geom_segment(data = cn_5_sample4, aes(x = start, xend = end, y = ycoord + (seg_dist - 0.48), yend = ycoord + (seg_dist - 0.48), color = "CN5"), size = 3, stat = "identity", position = position_dodge())} + #cn5
      {if("6+" %in% levels(sample4_cn$CN)) geom_segment(data = cn_6_sample4, aes(x = start, xend = end, y = ycoord + (seg_dist - 0.48), yend = ycoord + (seg_dist - 0.48), color = "CN6+"), size = 3, stat = "identity", position = position_dodge())} + #cn6+
      geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y - 0.5, yend = yend - 0.5, label = chr), color = "white", lineend = "butt", size = 2, stat = "identity", position = position_dodge()) + #white space between chromosome groups (upper)
      geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y + 0.5, yend = yend + 0.5, label = chr), color = "white", lineend = "butt", size = 2, stat = "identity", position = position_dodge()) + #white space between chromosome groups (bottom)
      annotate(geom = "text", x = anno_dist, y = segment_data$yend, label = segment_data$chr, color = "black", size = 5, hjust = 1) + #chr labels
      labs(title = plot_title, subtitle = plot_sub) + #plot titles
      scale_colour_manual(name = "", values = selected_colours) + #colours and legend
      scale_x_continuous(breaks = seq(0, max(segment_data$chr_end), by = 30000000)) + #x-axis boundaries
      scale_y_reverse() + #flip ideogram
      theme_Morons()
  }
  return(p)
}
