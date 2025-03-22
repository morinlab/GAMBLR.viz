
#' @title Represent mutation hot spots as a sequence logo
#'
#' @description Generates a summary of the frequency of SSMs in a small
#' region around a mutation hot spot using a "sequence logo"-style stacked plot
#'
#' @param maf_df
#' @param gene_symbol
#' @param hotspot_position
#' @param pad_length
#' @param fasta_path
#' @param include_reference
#' @param include_AA
#' @param group_AA
#' @param annotate_motif
#' @param annotate_only
#' @param return_data
#' @param base_size
#' @param aa_size
#' @rawNamespace import(IRanges, except = c("collapse", "union", "slice", "intersect", "setdiff", "desc", "reduce"))
#' @rawNamespace import(GenomicRanges, except = c("union", "intersect", "setdiff", "reduce"))
#' @import ggseqlogo Rsamtools cowplot BSgenome
#' @returns a named list containing a ggplot object and various processed data
#' @export
#'
#' @examples
#' g_meta = suppressMessages(get_gambl_metadata()) %>% 
#'   dplyr::filter(seq_type=="genome")
#' # CXCR4 non-coding hotspot
#' muts = get_ssm_by_regions(regions_list = "2:136875000-136875097",
#'                                 streamlined = F,
#'                                 projection = "grch37",
#'                                 these_samples_metadata = g_meta
#'                                     )
#'
#'
#' mutation_hotspot_logo(maf_df=muts,
#'                       hotspot_position=136875037,
#'                       include_AA = T,
#'                       pad_length = 51,
#'                       annotate_motif=T,
#'                       base_size = 2,
#'                       prepend_plot_label = "CXCR4")
#'                       
#'                       
#' cap_meta = suppressMessages(get_gambl_metadata()) %>% 
#'   dplyr::filter(seq_type=="capture")
#'   
#' #ID3 hotspot
#' muts = get_ssm_by_regions(regions_list = "1:23885706-23885798",
#'                                 streamlined = F,
#'                                 projection = "grch37",
#'                                 this_seq_type = "capture",
#'                                 these_samples_metadata = cap_meta
#'                                     )
#' 
#' 
#' mutation_hotspot_logo(maf_df=muts,
#'                               hotspot_position=23885750,
#'                               include_AA = T,aa_size = 3,
#'                               pad_length = 42,
#'                               annotate_motif=T,
#'                               base_size = 3,
#'                               prepend_plot_label = "ID3",text_size = 8)
#'
mutation_hotspot_logo = function(maf_df,
                             hotspot_position=128750677,
                             genome_build,
                             pad_length=20,
                             fasta_path,
                             include_AA = FALSE,
                             group_AA = FALSE,
                             annotate_motif=FALSE,
                             annotate_only=FALSE,
                             return_data=FALSE,
                             base_size=5,
                             aa_size=3,
                             text_size=5,
                             prepend_plot_label="",
                             include_title = TRUE,
                             real_coordinates = TRUE,
                             reverse = FALSE,
                             verbose = FALSE){

  bsgenome_loaded = FALSE
  include_reference = TRUE
  # If there is no fastaPath, it will read it from config key

  # Based on the genome_build the fasta file which will be loaded is different
  if (missing(fasta_path)){
    if("maf_data" %in% class(maf_df)){
      genome_build = get_genome_build(maf_df)
    }
    if(missing(genome_build)){
      stop("no genome_build information provided or present in maf")
    }
    base <- check_config_value(config::get("repo_base"))
    fasta_path <- paste0(
      base,
      "ref/lcr-modules-references-STABLE/genomes/",
      genome_build,
      "/genome_fasta/genome.fa"
    )
    if(!file.exists(fasta_path)){
      #try BSgenome
      installed = installed.genomes()
      if(genome_build=="hg38"){
        bsgenome_name = "BSgenome.Hsapiens.UCSC.hg38"
      }else if(genome_build == "grch37"){
        bsgenome_name = "BSgenome.Hsapiens.UCSC.hg19"
      }else{
        stop(paste("unsupported genome:",genome_build))
      }
      if(bsgenome_name %in% installed){
        genome = getBSgenome(bsgenome_name)
        bsgenome_loaded = TRUE
      }else{
        print(installed)
        print(paste("Local Fasta file cannot be found and missing genome_build",bsgenome_name,"Supply a fastaPath for a local fasta file or install the missing BSGenome package and re-run"))
      }
    }
  }
  # It checks for the presence of a local fastaPath
  if(!bsgenome_loaded){
    # Create a reference to an indexed fasta file.
    if (!file.exists(fasta_path)) {
      stop("Failed to find the fasta file and no compatible BSgenome found")
    }
    fasta = Rsamtools::FaFile(file = fasta_path)
  }
  mut_noindel = maf_df %>%
    filter(Variant_Type=="SNP")


  chrom = pull(mut_noindel,Chromosome)[1]
  start = hotspot_position - pad_length
  end = hotspot_position + pad_length
  maf_df_in_region = dplyr::filter(maf_df,Start_Position > start, Start_Position < end)
  if(nrow(maf_df_in_region)==0){
    print(paste0("Padded region: ",chrom,":",start,"-",end))
    print(paste0("mutation coordinate range: ",
                maf_df$Chromosome[1],":",
                min(maf_df$Start_Position),"-",
                max(maf_df$Start_Position)))
    stop("no mutations in maf_df are in the padded region around the specified coordinate")
  }
  
  some_mutations = filter(mut_noindel,
                          Start_Position<=end,
                          Start_Position >= start) %>%
    mutate(rel_start = Start_Position - start+1)
  if(bsgenome_loaded){
    if(genome_build == "grch37"){
      chrom_pre = paste0("chr",chrom)
    }else{
      chrom_pre = chrom
    }
    context = as.character(Rsamtools::getSeq(
                                           genome,
                                           chrom_pre,
                                           start = start,
                                          end = end)
                          )
  }else{
    context = as.character(Rsamtools::getSeq(fasta,
                                             GenomicRanges::GRanges(chrom,
                                                                    IRanges(start = start,
                                                                            end = end)))) %>% unname()


  }

  some_mutations = annotate_ssm_motif_context(some_mutations,
                                              fastaPath = fasta_path,
                                              genome_build = genome_build)

  if(annotate_only){
    return(list(annotated=some_mutations,region=context))
  }

  mod_string <- function(original, index, newbase) {
    # Create a string of N's with the same length as original
    masked_string <- paste(rep("N", nchar(original)), collapse = "")
    # Replace the character at the specified index with newbase
    substr(masked_string, index, index) <- newbase
    return(masked_string)
  }
  count_bases = function(original,index,newbase,base_counts,exclude_unmutated=T){
    all_bases = unlist(strsplit(original,""))
    for (i in c(1:length(all_bases))) {
      if(i == index){
        base = newbase
        base_counts[base,i] = base_counts[base,i] + 1
      }else {
        if (!exclude_unmutated) {
          base = all_bases[i]
          base_counts[base,i] = base_counts[base,i] + 1
        }
      }
    }
    return(base_counts)
  }




  some_mutations <- some_mutations %>%
  mutate(
    unmut = context,
    #mut = mod_string(context, rel_start, Tumor_Seq_Allele2),
    AA = gsub("p\\.", "", HGVSp_Short),
    To_AA = sub(".*(\\w)$", "\\1", AA),
    AA = sub("(\\d+).+", "\\1", AA)
  )

  base_counts = matrix(nrow=4,ncol=nchar(context))
  rownames(base_counts)=c("A","T","G","C")
  base_counts[] =0
  for(i in c(1:nrow(some_mutations))){
    base_counts = count_bases(as.character(some_mutations[i,"unmut"]),
                              some_mutations[i,"rel_start"],
                              as.character(some_mutations[i,"Tumor_Seq_Allele2"]),
                              base_counts)
  }
  num_base = nchar(context)
  if(real_coordinates){
    common_scale <- scale_x_continuous(
      limits = c(0, num_base+1),
      breaks = seq(1, num_base, by = 10),
      expand = c(0.1, 0.1),
      labels = seq(start, end, by = 10)
    )
  }else{
    common_scale <- scale_x_continuous(
      limits = c(0, num_base+1),
      breaks = 1:num_base,
      expand = c(0.1, 0.1)
    )
    
  }
  
  region_name = paste0(chrom,":",
                       hotspot_position - pad_length,
                       "-",
                       hotspot_position + pad_length)

  p1 = ggseqlogo::ggseqlogo(base_counts,
                            method='custom', seq_type='dna') +
    theme(axis.title.y =element_blank(),
          #axis.text.y=element_blank(),
          #axis.ticks.y=element_blank()
          )+
    common_scale +
    theme(
      axis.line.x  = element_line(color = "black"),
      axis.ticks.x = element_line(color = "black"),
      axis.text.x = element_blank()
    )
    
    #if(real_coordinates){
    #  p1 = p1 + theme(axis.text.x = element_text(angle = 45,
    #                                             vjust = 0.5,
    #                                             hjust=1,
    #                                             size = text_size))
    #}else{
    #  p1 = p1 + theme(axis.text.x = element_text(angle = 90,
    #                                             vjust = 0.5,
    #                                             hjust=1,
    #                                             size = text_size))
    #}
  
  if(include_title){
   p1 = p1 + ggtitle(paste(prepend_plot_label,region_name))
  }

  if(include_reference){
    if(annotate_motif){
      aln_df = dplyr::select(some_mutations,
                             rel_start,
                             Reference_Allele,WRCY)
      all_pos = data.frame(rel_start = c(1:nchar(context)))
      all_pos$Reference_Allele = unlist(strsplit(context,""))
      aln_df = left_join(all_pos,aln_df) %>% unique()

      aln_df = mutate(aln_df,WRCY=ifelse(is.na(WRCY),"unmutated",WRCY))
      if(verbose){
        print(aln_df)
      }
      
      p2 = ggplot(aln_df, aes(rel_start,1)) +
        geom_text(aes(label=Reference_Allele, color=WRCY),size=base_size) +
        common_scale +
        xlab('') +
        scale_color_manual(values=c('black', 'orange','red','grey')) +
        theme_logo() +
        #theme(legend.position = "none") +
        theme(axis.title.y =element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()
        ) +
        theme(
          axis.line.x  = element_line(color = "black"),
          axis.ticks.x = element_line(color = "black")
        )
      #if(real_coordinates){
      #  p2 = p2 + theme(axis.text.x = element_text(angle = 45,
      #                                             vjust = 0.5,
      #                                             hjust=1,
      #                                             size = text_size))
      #}else{
      #  p2 = p2 + theme(axis.text.x = element_text(angle = 90,
      #                                             vjust = 0.5,
      ##                                             hjust=1,
      #                                             size = text_size))
      #}
    }else{
      aln = data.frame(
        letter=strsplit(context, "")[[1]],
        x=rep(1:nchar(context))
      )
      pos = some_mutations$rel_start
      aln$mut = 'no'
      aln$mut[ pos ] = 'yes'

      p2 = ggplot(aln, aes(x,1)) +
        geom_text(aes(label=letter, color=mut),size=base_size) +
        common_scale +
                           xlab('') +
        scale_color_manual(values=c('black', 'red')) +
        theme_logo() +
        #theme(legend.position = "none") +
        theme(axis.title.y =element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()
        ) 
      #if(real_coordinates){
      #  p2 = p2 + theme(axis.text.x = element_text(angle = 45,
      #                                             vjust = 0.5,
      ##                                             hjust=1,
      #                                             size = text_size))
      #}else{
      #  p2 = p2 + theme(axis.text.x = element_text(angle = 90,
      #                                             vjust = 0.5,
      #                                             hjust=1,
      #                                             size = text_size))
      #}
    }

    if(include_AA | group_AA){
      if(group_AA){
        aa_df = some_mutations %>% dplyr::select(AA,rel_start) %>% unique()
      }else{
        aa_df = some_mutations %>% dplyr::select(AA,rel_start) %>% unique()
      }

      dummy1 = data.frame(AA="",rel_start=1)
      dummy2 = data.frame(AA="",rel_start=nchar(context))
      aa_df = bind_rows(dummy1,aa_df,dummy2)

      check = dplyr::select(some_mutations,AA,rel_start,HGVSp_Short)

      pAA = ggplot(aa_df, aes(rel_start,1)) +
        geom_text(aes(label=AA),size=aa_size,angle = 90) +
        common_scale +
        scale_color_manual(values=c('black', 'red')) +

        theme_logo() + theme(legend.position = "none") +

        theme(axis.title.y =element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.title.x = element_blank(),
              axis.text.x=element_blank(),
        )
      p1 <- p1 + theme(plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5))
      pAA <- pAA + theme(plot.margin = margin(t = -5, r = 5.5, b = -2, l = 5.5))
      p2 <- p2 + theme(plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5))
      
      #pp = cowplot::plot_grid(pAA,p1,p2,ncol=1,align="v",rel_heights = c(3,14,5))
      pp = cowplot::plot_grid(p1,pAA,p2,ncol=1,align="v",rel_heights = c(14,0.5,3))

      
    }else{
      p1 <- p1 + theme(plot.margin = margin(t = 5.5, r = 5.5, b = -5, l = 5.5))
      p2 <- p2 + theme(plot.margin = margin(t = -4, r = 5.5, b = 5.5, l = 5.5))
      pp = cowplot::plot_grid(p1,p2,ncol=1,align="v",rel_heights = c(7,2))

    }
  }
  if(return_data){
    if(include_AA){
      return(list(df=some_mutations,plot=pp,freqs=base_counts,plot1=p1,plot2=p2,plot3=pAA))
    }else{
      return(list(df=some_mutations,plot=pp,freqs=base_counts,plot1=p1,plot2=p2))
    }

  }else{
    print(pp)
  }

}

plot_signature = function(maf_df,genes="MYC",return_data=FALSE,separate_silent=FALSE){


  all_mutations = expand_grid(ref=c("T","C"),alt=c("A","C","T","G")) %>%
    mutate(mutation=paste0(ref,">",alt))
  all_contexts = expand_grid(up=c("A","C","T","G"),ref=c("T","C"),alt=c("A","C","T","G"),down=c("A","C","T","G"))
  all_contexts = mutate(all_contexts,seq=paste0(up,ref,down),mutation=paste0(ref,">",alt)) %>%
    dplyr::filter(!ref==alt) %>%
    mutate(mutation_context = paste(mutation,seq,sep="_"))
  all_mutations = dplyr::filter(maf_df,Hugo_Symbol %in% genes)
  silent_mutations = dplyr::filter(all_mutations,Variant_Classification %in% c("Silent","5'UTR","3'UTR","5'Flank","3'Flank"))
  all_count = group_by(all_mutations,mutation,seq,Hugo_Symbol) %>% tally()
  silent_count = group_by(silent_mutations,mutation,seq,Hugo_Symbol) %>% tally()
  all_count = left_join(all_contexts,all_count)
  silent_count = left_join(all_contexts,silent_count)
  all_count$Hugo_Symbol = factor(all_count$Hugo_Symbol,levels=genes)
  silent_count$Hugo_Symbol = factor(silent_count$Hugo_Symbol,levels=genes)
  if(return_data){
    all_count = dplyr::select(all_count,mutation_context,Hugo_Symbol,n) %>%
      pivot_wider(id_cols = Hugo_Symbol,names_from = mutation_context,values_from=n,values_fill = 0) %>%
      column_to_rownames("Hugo_Symbol")
    return(all_count)
  }
  if(separate_silent){
    silent_count = dplyr::filter(silent_count,!is.na(Hugo_Symbol))
    silent_count = mutate(silent_count,Hugo_Symbol = paste0(Hugo_Symbol,"-Silent"))

    all_count = bind_rows(silent_count,all_count)
    all_count = filter(all_count,!is.na(Hugo_Symbol))
    all_count = group_by(all_count,Hugo_Symbol) %>%
      mutate(gene_total=sum(n)) %>%
      mutate(gene_norm=n/gene_total) %>%
      ungroup()
  }
  ggplot(all_count,aes(x=mutation_context,y=gene_norm,fill=mutation)) + geom_col() +
    facet_wrap(~Hugo_Symbol,ncol=1) +
    theme_Morons(base_size=8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
  #ggplot(all_count,aes(x=mutation_context,y=n,fill=mutation)) + geom_col() +
  #  facet_wrap(~Hugo_Symbol,scales='free_y',ncol=1) +
  #  theme_Morons(base_size=8) +
  #    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))


}


