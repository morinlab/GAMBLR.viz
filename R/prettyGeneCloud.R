#' @title Gene Cloud Plot
#'
#' @description Make a word cloud of gene names from a MAF file based on mutation frequency.
#'
#' @details Create a wordcloud from an incoming MAF. Required parameter is `maf_df`.
#' Optional parameters are `these_genes`, `other_genes`, `these_genes_colour`, `other_genes_colour` and `colour_index`.
#' If no genes are supplied when calling the function, this function will default to all lymphoma genes.
#'
#' @param maf_df A MAF-format data frame containing the mutations you want to summarize
#' in a gene word cloud
#' @param these_genes An optional vector of gene symbols (defaults to all lymphoma
#' genes)
#' @param other_genes An optional second vector of gene symbols to include in your
#' cloud in a second colour
#' @param these_genes_colour Specify the hex code of a colour to use for the first
#' set of genes
#' @param other_genes_colour Specify another hex code of a colour to use for the
#' second set of genes
#' @param colour_index Optional named character vector with a name for each gene
#' in these_genes and a colour as the value
#' @param return_genes Optionally request a vector of the genes in the plot. 
#' Default is FALSE.
#' @param wordcloud_version Specify whether you want to use the original
#' wordcloud package or wordcloud2 (1 or 2)
#' @param zoomout Increase or decrease to fit more or less of the words
#' into view. Only works when wordcloud_version is 2. Default 1
#'
#' @return data frame with counts for each gene
#'
#' @import wordcloud wordcloud2 RColorBrewer dplyr
#' @export
#'
#' @examples
#' #get all coding SSM directly from GAMBLR.data
#' maf = GAMBLR.data::sample_data$grch37$maf
#'
#'
#' #build wordcloud
#' prettyGeneCloud(maf_df = maf,
#'                 wordcloud_version = 2,
#'                 zoomout=0.3)
#'
prettyGeneCloud = function(maf_df,
                           these_genes,
                           other_genes,
                           these_genes_colour="#B2DF8A",
                           other_genes_colour="#bc42f5",
                           colour_index,
                           return_genes = FALSE,
                           wordcloud_version = 2,
                           zoomout = 1,
                           fontFamily = "Arial"){
  if(missing(these_genes)){
    these_genes = filter(GAMBLR.data::lymphoma_genes,
                         DLBCL_Tier ==1 | FL_Tier == 1 | BL_Tier == 1) %>%
                  pull(Gene)
  }
  #drop genes not in the list then tally only coding variants (by default).
  # TODO: eventually allow an option to collapse samples from the same patient
  if(missing(other_genes)){
    these_genes_maf = dplyr::filter(maf_df,Hugo_Symbol %in% these_genes)
  }else{
    these_genes_maf = dplyr::filter(maf_df,Hugo_Symbol %in% c(these_genes,other_genes))
  }


  #drop non-coding
  these_genes_maf = dplyr::filter(these_genes_maf,Variant_Classification %in% coding_vc)
  these_genes_unique = group_by(these_genes_maf,Hugo_Symbol,Tumor_Sample_Barcode) %>%
    slice_head() %>% ungroup() %>% group_by(Hugo_Symbol) %>% tally()
  these_genes_unique = group_by(these_genes_maf,Hugo_Symbol,Tumor_Sample_Barcode) %>%
    slice_head() %>% ungroup() %>% group_by(Hugo_Symbol) %>% tally()
  make_wordcloud2 = function(genes,
                             frequencies,
                             colours,
                             ordered.colors,
                             scaling,
                             random.order,
                             zoomout,
                             fontFamily){
    gene_df = data.frame(word=genes,freq=frequencies,
                        colour = sample(colours,
                        length(genes),replace=T))
    gene_df = arrange(gene_df,desc(freq))
    wc = wordcloud2::wordcloud2(gene_df,
                          color=gene_df$colour,
                          size=zoomout,
                          fontFamily=fontFamily)
    return(wc)
  }
  if(!missing(other_genes)){
    #assign a colour to each gene list
    these_genes_unique = these_genes_unique %>%
      mutate(this_col=ifelse(Hugo_Symbol %in% these_genes,
                             these_genes_colour,
                             other_genes_colour)) %>% 
      arrange(dplyr::desc(n))
    if(wordcloud_version == 1){
          wc = wordcloud::wordcloud(these_genes_unique$Hugo_Symbol,
                              these_genes_unique$n,colors=these_genes_unique$this_col,
                              ordered.colors = T,
                              scale=c(8,0.3),random.order = F)
    }

  }else{
    if(!missing(colour_index)){
      #use the colours in colour_index to colour the gene names
      if(any(!these_genes %in% names(colour_index))){
        stop("all genes in these_genes must be among the names of colour_index if you specify this variable")
      }
      these_genes_unique$color=colour_index[these_genes_unique$Hugo_Symbol]
      #make cloud with the user-specified colours mapped to the genes
      if(wordcloud_version == 1){
        wc = wordcloud::wordcloud(these_genes_unique$Hugo_Symbol,
                          these_genes_unique$n,
                          random.order=F,
                          ordered.colors=T,
                          colors=these_genes_unique$color)
      }else if(wordcloud_version == 2){
         wc = make_wordcloud2(these_genes_unique$Hugo_Symbol,
                          these_genes_unique$n,
                          random.order=F,
                          ordered.colors=T,
                          colours=these_genes_unique$color,
                          zoomout=zoomout,
                          fontFamily=fontFamily)
      }
    }else{
      if(wordcloud_version == 1){
        wc = wordcloud::wordcloud(these_genes_unique$Hugo_Symbol,
                                  these_genes_unique$n,random.color=TRUE,
                                  colors=RColorBrewer::brewer.pal(12,"Set2"))
      }else if(wordcloud_version == 2){
         wc = make_wordcloud2(these_genes_unique$Hugo_Symbol,
                          these_genes_unique$n,
                          random.order=F,
                          ordered.colors=T,
                          colours=RColorBrewer::brewer.pal(12,"Set2"),
                          zoomout=zoomout,
                          fontFamily=fontFamily)
      }
    }
  }
  these_genes_unique = arrange(these_genes_unique,n)
  these_genes_unique$Hugo_Symbol = factor(these_genes_unique$Hugo_Symbol,
                                          levels=these_genes_unique$Hugo_Symbol)
  if(return_genes){
    return(these_genes_unique)
  }else{
    return(wc)
  }
  
}
