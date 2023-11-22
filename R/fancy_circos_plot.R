#' @title SSM and SV Circos Plot.
#'
#' @description Create a circos plot visualizing SVS and SSM with optional gene annotations.
#'
#' @details This function is using RCircos to create sample-level cirocs plots, annotating SVs and SSM with the potential of adding gene annotations.
#' To control what variants are to be plotted, simply use the two Boolean parameters; `ssm_calls` and `sv_calls` (both TRUE by default).
#' Provide the sample ID of interest in with the `this_sample_id` parameter. This function calls get_ssm_by_sample and get_manta_sv to retrieve data for plotting.
#' Since this function does not create a grob, but rather outputs a rendered PDF/PNG, the user has to provide an output path with the `out` parameter.
#' In addition, the user can control the output format. For PDF, set `pdf` to TRUE (default) and to export the created plot as PNG, set the same parameter to FALSE.
#' This function also has convenient filtering parameters available, see parameter descriptions for more information and how to properly use the filtering parameters.
#' Lastly, this plot can also highlight genes of interest. To do so, provide a data frame (comparable to the return from `gene_to_region(return_as = "bed")`) to the `gene_list` parameter.
#'
#' @param this_sample_id Sample to be plotted.
#' @param gene_list Optional parameter to annotate genes on the circos plot from a data frame of genes. Is compatible with [GAMBLR::gene_to_region] (return_as = "bed") output format. See examples.
#' @param ssm_calls Boolean parameter for plotting ssm. Default is TRUE.
#' @param sv_calls Boolean parameter for plotting SVs, default is TRUE.
#' @param chr_select Optional argument for subset on selected chromosomes, default is all autosomes.
#' @param vaf_cutoff Threshold for filtering variants on VAF (events with a VAF > cutoff will be retained).
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param projection Genomic projection for variants and circos plot. Accepted values are grch37 and hg38, default is grch37.
#' @param this_seq_type Seq type for returned CN segments. One of "genome" (default) or "capture".
#' @param out Path to output folder, to where the plot will be exported.
#' @param plot_title Optional parameter for naming your plot, default is this_sample.
#' @param pdf Set to FALSE for png, default is TRUE (pdf).
#' @param file_name Optional parameter for specifying the file name of generated circos plot, default is "{this_sample}_circos.pdf". If pdf is set to FALSE, a png will be generated, thus the .png extension needs to be attached to the file_name.
#'
#' @return Nothing.
#'
#' @import dplyr RCircos
#' @export
#'
#' @examples
#' \dontrun{
#' #retrieve gene names for FL genes
#' fl_genes = dplyr::filter(GAMBLR.data::lymphoma_genes_lymphoma_genes_v0.0, FL == TRUE) %>%
#'   pull(Gene)
#'
#' # get regions for selected genes
#' fl_genes_list = GAMBLR.utils::gene_to_region(gene_symbol = fl_genes,
#'                                return_as = "bed")
#'
#' fancy_circos_plot(this_sample_id = "DOHH-2",
#'                   ssm_calls = FALSE,
#'                   gene_list = fl_genes_list,
#'                   chr_select = c("chr8",
#'                                  "chr14",
#'                                  "chr18"),
#'                   out = "../../plots/",
#'                   plot_title = "DOHH-2 (SVs) Example Plot",
#'                   pdf = FALSE,
#'                   file_name = "dohh2_example.png")
#' }
#'
fancy_circos_plot = function(this_sample_id,
                             gene_list,
                             ssm_calls = TRUE,
                             sv_calls = TRUE,
                             chr_select = paste0("chr", c(1:22)),
                             vaf_cutoff = 0,
                             coding_only = FALSE,
                             projection = "grch37",
                             this_seq_type = "genome",
                             plot_title = paste0(this_sample_id),
                             out,
                             pdf = TRUE,
                             file_name = paste0(this_sample_id, "_circos.pdf")){

  #set track properties based on selected plotting data
  if(ssm_calls && sv_calls && !missing(gene_list)){
    gene_con_track = 1
    gene_name_track = 2
    sv_del_track = 4
    sv_dup_track = 5
    ssm_snp_track = 6
    ssm_dnp_track = 7
    ssm_del_track = 8
    ssm_ins_track = 9
    trans_track = 10
  }


  if(ssm_calls && sv_calls && missing(gene_list)){
    sv_del_track = 1
    sv_dup_track = 2
    ssm_snp_track = 3
    ssm_dnp_track = 4
    ssm_del_track = 5
    ssm_ins_track = 6
    trans_track = 7
  }

  if(ssm_calls && !sv_calls && !missing(gene_list)){
    gene_con_track = 1
    gene_name_track = 2
    ssm_snp_track = 4
    ssm_dnp_track = 5
    ssm_del_track = 6
    ssm_ins_track = 7
  }

  if(ssm_calls && !sv_calls && missing(gene_list)){
    ssm_snp_track = 1
    ssm_dnp_track = 2
    ssm_del_track = 3
    ssm_ins_track = 4
  }

  if(!ssm_calls && sv_calls && !missing(gene_list)){
    gene_con_track = 1
    gene_name_track = 2
    sv_del_track = 4
    sv_dup_track = 5
    trans_track = 6
  }

  if(!ssm_calls && sv_calls && missing(gene_list)){
    sv_del_track = 1
    sv_dup_track = 2
    trans_track = 3
  }

  #sanity checking incoming gene list and renaming columns if needed
  if(!missing(gene_list)){
    #check type of incoming gene_list
    if(!is.list(gene_list)){
      message("Please ensure that incoming gene list is in fact a data frame with the following columns: chr:start:end:gene")
    }

    #rename columns, if needed (Rcircos is expecting the column names to always be the same...)
    if(!"Chromosome" %in% colnames(gene_list)[1]){
      colnames(gene_list)[1] = "Chromosome"
    }
    if(!"chromStart" %in% colnames(gene_list)[2]){
      colnames(gene_list)[2] = "chromStart"
    }
    if(!"chromEnd" %in% colnames(gene_list)[3]){
      colnames(gene_list)[3] = "chromEnd"
    }
    if(!"Gene" %in% colnames(gene_list)[4]){
      colnames(gene_list)[4] = "Gene"
    }

    #add "chr" prefix, if needed
    if( any( !str_detect(gene_list$Chromosome, "chr") ) ){
      stopifnot( "Inconsistent chromosome names in argument `gene_list`." = all( !str_detect(gene_list$Chromosome, "chr") ) )
      gene_list = mutate(gene_list, Chromosome = paste0("chr", Chromosome))
    }

    #filter gene list on selected chromosomes
    gene_list = gene_list[gene_list$Chromosome %in% chr_select, ]
  }


  #get SSM data
  if(ssm_calls){
    maf = get_ssm_by_sample(these_sample_ids = this_sample_id, projection = projection, this_seq_type = this_seq_type)
    maf_tmp = dplyr::select(maf, Chromosome, Start_Position, End_Position, Variant_Type) #select appropriate columns
    maf_tmp$Variant_Size = maf_tmp$End_Position - maf_tmp$Start_Position # calcualte variant size
    maf_tmp$Variant_Type = as.factor(maf_tmp$Variant_Type) #transform Variant_Type to factor
    maf_tmp[maf_tmp==0] <- 1 #transform all lenght coordinates == 0 to 1

    if( any( !str_detect(maf_tmp$Chromosome, "chr") ) ){ #add chr prefix, if missing...
      maf_tmp = mutate(maf_tmp, Chromosome = paste0("chr", Chromosome))
    }

    maf_tmp = maf_tmp[maf_tmp$Chromosome %in% chr_select, ] #filter incoming maf on selected chromosomes
    ssm_del = dplyr::filter(maf_tmp, Variant_Type == "DEL") #subset on deletions
    ssm_ins = dplyr::filter(maf_tmp, Variant_Type == "INS") #subset on insertions
    ssm_snp = dplyr::filter(maf_tmp, Variant_Type == "SNP") #subset on single nucleotide polymorphism
    ssm_dnp = dplyr::filter(maf_tmp, Variant_Type == "DNP") #subset on dinucleotide polymorphism
    message(paste0(nrow(ssm_del) + nrow(ssm_dnp) + nrow(ssm_ins) + nrow(ssm_snp)), " SSMs found for ", this_sample_id)
  }

  #get SVs
  if(sv_calls){
    svs = get_manta_sv(these_sample_ids = this_sample_id, projection = projection)

    #filter on vaf
    svs = dplyr::filter(svs, VAF_tumour > vaf_cutoff)

    #subset on relevant variables
    svs_df = dplyr::select(svs, CHROM_A, START_A, END_A, CHROM_B, START_B, END_B, manta_name)

    #split manta_name variable
    svs_df = data.frame( svs_df$CHROM_A, svs_df$START_A, svs_df$END_A, 
                         svs_df$CHROM_B, svs_df$START_B, svs_df$END_B, 
                         sub("^(.+?):.*", "\\1", svs_df$manta_name) )

    #rename variables
    colnames(svs_df)[1:7] = c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "TYPE")

    #add chr prefix, if missing
    if(!str_detect(svs_df$CHROM_A, "chr")[1]){
      svs_df = mutate(svs_df, CHROM_A = paste0("chr", CHROM_A))
    }

    if(!str_detect(svs_df$CHROM_B, "chr")[4]){
      svs_df = mutate(svs_df, CHROM_B = paste0("chr", CHROM_B))
    }

    #subset on selected chromosomes
    svs_df = svs_df[svs_df$CHROM_A %in% chr_select, ]
    svs_df = svs_df[svs_df$CHROM_B %in% chr_select, ]

    #subset df on SV type
    sv_trans = dplyr::filter(svs_df, TYPE == "MantaBND") %>%
      dplyr::select(CHROM_A, START_A, END_A, CHROM_B, START_B, END_B)

    sv_del = dplyr::filter(svs_df, TYPE == "MantaDEL") %>%
      dplyr::select(CHROM_A, START_A, END_A)

    sv_dup = dplyr::filter(svs_df, TYPE == "MantaDUP") %>%
      dplyr::select(CHROM_A, START_A, END_A)

    #calculate sizes
    sv_del$SIZE = sv_del$END_A - sv_del$START_A
    sv_dup$SIZE = sv_dup$END_A - sv_dup$START_A

    message(paste0(nrow(sv_trans) + nrow(sv_del) + nrow(sv_dup)), " SVs found for ", this_sample_id)
  }

  #plotting
  #define reference build
  if(projection == "grch37"){
    data(UCSC.HG19.Human.CytoBandIdeogram, package = "RCircos")
    cytobands = UCSC.HG19.Human.CytoBandIdeogram
  }else if(projection == "hg38"){
    data(UCSC.HG38.Human.CytoBandIdeogram, package = "RCircos")
    cytobands = UCSC.HG38.Human.CytoBandIdeogram
  }


  #get chr excluded (reversed of chr included, since RCircos only accept this)
  chr_all = paste0("chr", c(1:22, "X", "Y"))
  chr_exclude = setdiff(chr_all, chr_select)

  #set core components
  suppressMessages(RCircos::RCircos.Set.Core.Components(cyto.info = cytobands, chr.exclude = chr_exclude, tracks.inside = 10))

  #set plot parameters
  RCircos.params = RCircos.Get.Plot.Parameters()

  #define plotting parameters
  out.file = paste0(out, file_name)

  if(pdf){
    pdf(out.file, height = 7, width = 7)
    pdf(out.file, height = 7, width = 7)
  }else{
    png(out.file, height = 7, width = 7, units = "in", res = 300)
    png(out.file, height = 7, width = 7, units = "in", res = 300)
  }

  RCircos.Set.Plot.Area(margins = 0);

  #create empty plot
  RCircos.Chromosome.Ideogram.Plot()

  #add tracks to plot
  #add gene names
  if(!missing(gene_list)){
    RCircos.Gene.Connector.Plot(gene_list, gene_con_track, "in")
    RCircos.Gene.Name.Plot(gene_list, 4, gene_name_track, "in")
  }

  #aadd tracks
  if(ssm_calls){
    #ssm deletions
    RCircos.params$track.background = "steelblue2"
    RCircos.params$max.layers = 1
    RCircos.Reset.Plot.Parameters(RCircos.params)
    RCircos.Tile.Plot(ssm_del, track.num = ssm_del_track, side = "in")

    #ssm insertions
    RCircos.params$track.background = "sienna2"
    RCircos.params$max.layers = 1
    RCircos.Reset.Plot.Parameters(RCircos.params)
    RCircos.Tile.Plot(ssm_ins, track.num = ssm_ins_track, side = "in")

    #ssm snp
    RCircos.params$track.background = "seagreen"
    RCircos.params$max.layers = 1
    RCircos.Reset.Plot.Parameters(RCircos.params)
    RCircos.Tile.Plot(ssm_snp, track.num = ssm_snp_track, side = "in")

    #ssm dnp
    RCircos.params$track.background = "tomato4"
    RCircos.params$max.layers = 1
    RCircos.Reset.Plot.Parameters(RCircos.params)
    RCircos.Tile.Plot(ssm_dnp, track.num = ssm_dnp_track, side = "in")
  }

  if(sv_calls){
    #translocations
    RCircos.Link.Plot(sv_trans, track.num = trans_track, by.chromosome = FALSE)

    #duplications
    RCircos.params$track.background = "sienna2"
    RCircos.params$max.layers = 1
    RCircos.Reset.Plot.Parameters(RCircos.params)
    RCircos.Tile.Plot(sv_dup, track.num = sv_dup_track, side = "in")

    #deletions
    RCircos.params$track.background = "steelblue2"
    RCircos.params$max.layers = 1
    RCircos.Reset.Plot.Parameters(RCircos.params)
    RCircos.Tile.Plot(sv_del, track.num = sv_del_track, side = "in")
  }

  #add plot title and legends
  if(sv_calls && !ssm_calls){
    text(x = 0, y = 2.5, plot_title, font = 2, cex = 1.2)
    legend(x = 2, y = -1.2, legend = c("Del", "Dup"), bty = "n", col = c("steelblue2", "sienna2"), pch = 19, text.col = "black", horiz = FALSE, pt.cex = 2, title = "SVs", inset = c(0.1, 0.1))
    legend(x = 2, y = -1.2, legend = c("Del", "Dup"), bty = "n", col = c("steelblue2", "sienna2"), pch = 19, text.col = "black", horiz = FALSE, pt.cex = 2, title = "SVs", inset = c(0.1, 0.1))
  }

  if(!sv_calls && ssm_calls){
    text(x = 0, y = 2.5, plot_title, font = 2, cex = 1.2)
    legend(x = 2, y = -1.8, legend = c("SNP", "DNP", "Del", "Ins"), bty = "n", col = c("seagreen", "tomato4", "steelblue2", "sienna2"), pch = 19, text.col = "black", horiz = FALSE, pt.cex = 2, title = "SSM", inset = c(0.1, 0.1))
    legend(x = 2, y = -1.8, legend = c("SNP", "DNP", "Del", "Ins"), bty = "n", col = c("seagreen", "tomato4", "steelblue2", "sienna2"), pch = 19, text.col = "black", horiz = FALSE, pt.cex = 2, title = "SSM", inset = c(0.1, 0.1))
  }

  if(sv_calls && ssm_calls){
    text(x = 0, y = 2.5, plot_title, font = 2, cex = 1.2)
    legend(x = 2, y = -1.2, legend = c("Del", "Dup"), bty = "n", col = c("steelblue2", "sienna2"), pch = 19, text.col = "black", horiz = FALSE, pt.cex = 2, title = "SVs", inset = c(0.1, 0.1))
    legend(x = 2, y = -1.8, legend = c("SNP", "DNP", "Del", "Ins"), bty = "n", col = c("seagreen", "tomato4", "steelblue2", "sienna2"), pch = 19, text.col = "black", horiz = FALSE, pt.cex = 2, title = "SSM", inset = c(0.1, 0.1))
    legend(x = 2, y = -1.2, legend = c("Del", "Dup"), bty = "n", col = c("steelblue2", "sienna2"), pch = 19, text.col = "black", horiz = FALSE, pt.cex = 2, title = "SVs", inset = c(0.1, 0.1))
    legend(x = 2, y = -1.8, legend = c("SNP", "DNP", "Del", "Ins"), bty = "n", col = c("seagreen", "tomato4", "steelblue2", "sienna2"), pch = 19, text.col = "black", horiz = FALSE, pt.cex = 2, title = "SSM", inset = c(0.1, 0.1))
  }

  invisible(dev.off())
}
