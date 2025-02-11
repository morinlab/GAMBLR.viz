
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
#' @import ggseqlogo Rsamtools cowplot
#' @returns a named list containing a ggplot object and various processed data
#' @export
#'
#' @examples
#' meta = get_gambl_metadata() %>% dplyr::filter(seq_type=="genome")
#' muts = get_ssm_by_regions(regions_list = "2:136875000-136875097",
#'                                 streamlined = F,
#'                                 these_samples_metadata = meta
#'                                     )
#'                                     
#'                                     
#' gene_hotspot_logo(maf_df=muts,
#'                       gene_symbol = "CXCR4",
#'                       hotspot_position=136875037,
#'                       include_AA = T,
#'                       pad_length = 51,
#'                       annotate_motif=T,
#'                       base_size = 4)

gene_hotspot_logo = function(maf_df,
                             gene_symbol="MYC",
                             hotspot_position=128750677,
                             genome_build = "grch37",
                             pad_length=20,
                             fasta_path="~/Downloads/genome.fa",
                             include_reference = TRUE,
                             include_AA = FALSE,
                             group_AA = FALSE,
                             annotate_motif=FALSE,
                             annotate_only=FALSE,
                             return_data=TRUE,
                             base_size=5,
                             aa_size=3){
  fasta <- Rsamtools::FaFile(file = fasta_path)
  
  mut_noindel = maf_df %>%
    dplyr::filter(Hugo_Symbol==gene_symbol) %>%
    filter(Variant_Type=="SNP")
  

  chrom = pull(mut_noindel,Chromosome)[1]
  start = hotspot_position - pad_length
  end = hotspot_position + pad_length
  print(paste(chrom,start,end))
  some_mutations = filter(mut_noindel,Start_Position<=end, Start_Position >= start) %>%
    mutate(rel_start = Start_Position - start+1)
  context = as.character(Rsamtools::getSeq(fasta,GenomicRanges::GRanges(chrom,IRanges(start = start,end = end)))) %>% unname()
  
  some_mutations = annotate_ssm_motif_context(some_mutations,fastaPath = fasta_path)
  
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
  
  
  #some_mutations = mutate(some_mutations,
  #                        unmut = context,
  #                        mut = mod_string(context,rel_start,Tumor_Seq_Allele2)) %>%
  # mutate(AA = str_remove(HGVSp_Short,"p\\.")) %>%
  #  mutate(To_AA = str_extract(AA,"\\w$")) %>%
  #  mutate(AA=str_replace(AA,"(\\d+).+","\\1"))
  some_mutations <- some_mutations %>%
  mutate(
    unmut = context,
    mut = mod_string(context, rel_start, Tumor_Seq_Allele2),
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
  

  region_name = paste0(hotspot_position - pad_length,"-",hotspot_position + pad_length)
  
  p1 = ggseqlogo::ggseqlogo(base_counts, method='custom', seq_type='dna') + 
    theme(axis.text.x = element_blank()) + ggtitle(paste(gene_symbol,region_name)) 
  if(include_reference){
    if(annotate_motif){
      aln_df = dplyr::select(some_mutations,rel_start,Reference_Allele,WRCY)
      all_pos = data.frame(rel_start = c(1:nchar(context)))
      all_pos$Reference_Allele = unlist(strsplit(context,""))
      aln_df = left_join(all_pos,aln_df) %>% unique()
      
      aln_df = mutate(aln_df,WRCY=ifelse(is.na(WRCY),"unmutated",WRCY))
      print(aln_df)
      p2 = ggplot(aln_df, aes(rel_start,1)) +
        geom_text(aes(label=Reference_Allele, color=WRCY),size=base_size) + 
        scale_x_continuous(breaks=1:nchar(context), expand = c(0.065, 0)) + xlab('') + 
        scale_color_manual(values=c('black', 'orange','red','grey')) + 
        theme_logo() + 
        #theme(legend.position = "none") + 
        theme(axis.title.y =element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank() 
        ) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
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
        scale_x_continuous(breaks=1:nchar(context), expand = c(0.065, 0)) + xlab('') + 
        scale_color_manual(values=c('black', 'red')) + 
        theme_logo() + 
        #theme(legend.position = "none") + 
        theme(axis.title.y =element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank() 
        ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
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
      print(aa_df)
      check = dplyr::select(some_mutations,AA,rel_start,HGVSp_Short)

      pAA = ggplot(aa_df, aes(rel_start,1)) +
        geom_text(aes(label=AA),size=aa_size,angle = 90) + 
        scale_x_continuous(breaks=1:nchar(context), expand = c(0.065, 0)) + xlab('') + 
        scale_color_manual(values=c('black', 'red')) + 
        
        theme_logo() + theme(legend.position = "none") + 
        
        theme(axis.title.y =element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank() 
        ) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
      #pAA
      pp = cowplot::plot_grid(p1,pAA,p2,ncol=1,align="v",rel_heights = c(7,2,2))
      print(pp)
    }else{
      pp = cowplot::plot_grid(p1,p2,ncol=1,align="v",rel_heights = c(7,2))
      print(pp)
    }
    
  }else{
    print(p1)
  }
  if(return_data){
    return(list(df=some_mutations,plot=pp))
  }
  
}

plot_signature = function(maf_df,genes="MYC",return_data=FALSE,separate_silent=FALSE){

 
  all_mutations = expand_grid(ref=c("T","C"),alt=c("A","C","T","G")) %>% mutate(mutation=paste0(ref,">",alt))
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
                             
  
