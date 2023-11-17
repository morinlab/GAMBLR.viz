ssm_to_proteinpaint = function(gene_symbols, regions, maf, these_samples_metadata, 
                               print_removed_rows = FALSE, ...){
  
  # check input parameters
  is_missing <- c( gene_symbols = missing(gene_symbols), regions = missing(regions),
                   maf = missing(maf), these_samples_metadata = missing(these_samples_metadata) )
  stopifnot( "None or only one of arguments `gene_symbols`, `regions` or `maf` should be provided." = any( 0:1 %in% sum(!is_missing[c("gene_symbols", "regions", "maf")]) ) )
  
  # get arguments passed via `...`
  arguments <- list(...)
  # basic_columns is always FALSE 
  arguments$basic_columns <- FALSE
  
  # check whether `...` contains the right arguments
  # if `get_coding_ssm` will be used
  all_missing <- all( is_missing[c("gene_symbols", "regions", "maf")] )
  if( !is_missing["gene_symbols"] | all_missing ){
    right_args <- names(arguments) %in% formalArgs(get_coding_ssm)
    if( !all(right_args) ){
      k <- names(arguments)[!right_args] %>% 
        paste(collapse = ", ") %>% 
        gettextf("If argument `gene_symbols` is provided, `get_coding_ssm` is used internally and `...` refers to arguments of this function.\nOffending arguments: %s", .)
      stop(k)
    }
  }
  # if `get_ssm_by_regions` will be used
  if( !is_missing["regions"] ){
    right_args <- names(arguments) %in% formalArgs(get_ssm_by_regions)
    if( !all(right_args) ){
      k <- names(arguments)[!right_args] %>% 
        paste(collapse = ", ") %>% 
        gettextf("If argument `regions` is provided, `get_ssm_by_regions` is used internally and `...` refers to arguments of this function.\nOffending arguments: %s", .)
      stop(k)
    }
  }
  
  if( is_missing["maf"] ){
    # if maf is not provided
    
    # inform that provided these_samples_metadata is ignored if any
    if( !is_missing["these_samples_metadata"] ){
      message("Warning: `these_samples_metadata` argument is used only when `maf` if specified and will be ignorred.")
    }
    
    # to ensure the metadata is compatible with the maf, internally create the metadata
    ggmt_arguments <- names(arguments) %>% 
      recode(seq_type = "seq_type_filter", from_indexed_flatfile = "from_flatfile") %>% 
      `names<-`(arguments, .) %>% 
      "["( ., names(.) %in% formalArgs(get_gambl_metadata) )
    these_samples_metadata <- do.call(get_gambl_metadata, ggmt_arguments)
    
    if( !is_missing["gene_symbols"] ){
      # if gene_symbols is provided
      # get maf with get_coding_ssm and subset by gene names
      maf <- do.call(get_coding_ssm, arguments) %>% 
        filter(Hugo_Symbol == "MYC")
    }else if( !is_missing["regions"] ){
      # if gene_symbols is regions
      # get maf with get_ssm_by_regions
      
      # figure out if `regions` is `regions_list` or `regions_bed`   # <<<<<<<<<<<<<<<<<<<<<<<<<<<< todo
      # and add it to `arguments`
      
      maf <- do.call(get_ssm_by_regions, arguments)
    }else if(all_missing){
      # if no specific gene (`gene_symbols`) neither region (`regions`) is provided
      # get maf with get_coding_ssm from all genes
      maf <- do.call(get_coding_ssm, arguments)
    }
  }
  
  if( !is_missing["maf"] & is_missing["these_samples_metadata"] ){
    # if maf is provided but not these_samples_metadata
    message("Warning: `these_samples_metadata` argument is missing when `maf` was provided. Hence, some optional columns are going to be omitted.")
  }else{
    # add metadata columns to maf
    maf <- left_join(maf, these_samples_metadata, by = join_by(Tumor_Sample_Barcode == sample_id))
  }
  
  # adequate maf to ProteinPaint visualization format
  maf <- rename(maf, gene = Hugo_Symbol, chromosome = Chromosome, 
                start = Start_Position, class = Variant_Classification,
                refseq = RefSeq, aachange = HGVSp_Short, origin = Mutation_Status,
                disease = pathology, mutant_in_tumor = t_alt_count, 
                total_in_tumor = t_depth, mutant_in_normal = n_alt_count, 
                total_in_normal = n_depth, patient = patient_id, 
                sample = Tumor_Sample_Barcode, sampletype = time_point,
                REF = Reference_Allele) %>% 
    mutate(class = recode(class,
                          Missense_Mutation = "missense",
                          In_Frame_Del = "proteinDel",
                          In_Frame_Ins = "proteinIns",
                          Nonsense_Mutation = "nonsense",
                          Nonstop_Mutation = "nonsense",
                          Silent = "silent",
                          Splice_Region = "splice_region",
                          Splice_Site = "splice",
                          Frame_Shift_Del = "frameshift",
                          Frame_Shift_Ins = "frameshift"))
  
  # remove version names in refseq column
  maf$refseq <- maf$refseq %>% 
    gsub("\\.[0-9]+,", ",", .) %>% 
    gsub("\\.[0-9]+$", "", .)
  
  # if alt allele is heterogeneous alternative, show both allele separated 
  # by comma (like in VCF files)
  is_het_alt <- ! maf$REF == maf$Tumor_Seq_Allele1
  maf$ALT <- is_het_alt %>% 
    { paste( maf$Tumor_Seq_Allele1[.], maf$Tumor_Seq_Allele2[.], sep = "," ) } %>% 
    replace(maf$Tumor_Seq_Allele2, is_het_alt, .)
  
  # add VAF (variant allele fraction)
  maf <- mutate(maf, VAF = mutant_in_tumor/total_in_tumor)
  
  # keep only columns that are meaningful for proteinpaint
  maf <- select(maf, gene, refseq, chromosome, start, aachange, class, disease, origin,
                patient, sample, sampletype, mutant_in_tumor, total_in_tumor, 
                mutant_in_normal, total_in_normal, REF, ALT, VAF)
  
  # filter out rows that don't contain required column values and output warning message
  required_cols <- select(maf, gene, refseq, chromosome, start, aachange, class) %>% 
    mutate(i = row_number()) %>% 
    drop_na()
  required_cols <- apply(required_cols != "", 1, all) %>% 
    filter(required_cols, .)
  if(print_removed_rows){  # if only removed rows should be returned (for checking purpose)
    removed_cols <- 1:nrow(maf) %>% 
      "["(! . %in% required_cols$i) %>% 
      slice(maf, .)
    return(removed_cols)
  }
  num_removed_rows <- (nrow(maf) - nrow(required_cols)) %>% 
    format(big.mark="'")
  k <- format( nrow(maf), big.mark="'" ) %>% 
    gettextf("Warning: %s rows out of %s were removed from the output table because there were missing values in required columns. Run `ssm_to_proteinpaint` again with `print_removed_rows = TRUE` to see these rows.",
             num_removed_rows, .)
  message(k)
  maf <- slice(maf, required_cols$i)
  
  # return filnal output table
  maf
}
