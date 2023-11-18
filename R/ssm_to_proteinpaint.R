#' @title MAF to ProteinPaint table
#' 
#' @description This function takes a MAF-like data frame and convert it to the right format 
#' for visualization with ProteinPaint.
#' 
#' @details For visualization with ProteinPaint, the output data frame must be saved to a 
#' file and uploaded to https://proteinpaint.stjude.org/. 
#' 
#' For a valid ProteinPaint table for visualization, the MAF table must contain the required 
#' columns to be converted to the ProteinPaint format. They are: Hugo_Symbol, RefSeq, 
#' Chromosome, Start_Position, HGVSp_Short, Variant_Classification. Optional columns used 
#' from the MAF table are: Mutation_Status, t_alt_count, t_depth, n_alt_count, n_depth, 
#' Tumor_Sample_Barcode, Reference_Allele. In the case of using functions `get_ssm_by_samples`, 
#' `get_coding_ssm`, or `get_ssm_by_patients` to get your MAF table and required/desired 
#' columns are missing, consider using parameter `basic_columns = FALSE` to ensure that your 
#' MAF has all needed columns. Other columns are taken from the metadata table generated 
#' internally in `ssm_to_proteinpaint`.
#' 
#' @param maf A data frame in MAF format.
#' @param print_removed_rows Boolean parameter. Set to TRUE for returning rows from the 
#' incoming MAF that do not contain any values in the required columns. Commonly used for 
#' checking purposes only. Setting this to TRUE, does not produce an output compatible with 
#' Protein paint. The default is FALSE.
#'
#' @return A data frame.
#' 
#' @import tidyr
#' @export
#'
#' @examples
#' library(GAMBLR.data)
#' 
#' # get maf
#' my_maf = get_coding_ssm(basic_columns = FALSE)
#' 
#' # convert maf to ProteinPaint format
#' pp_df = ssm_to_proteinpaint(my_maf)
#' 
ssm_to_proteinpaint = function(maf,
                               print_removed_rows = FALSE){
  
  # check for required columns in maf 
  maf_req_cols = c("Hugo_Symbol", "RefSeq", "Chromosome", "Start_Position", "HGVSp_Short", 
                    "Variant_Classification")
  maf_req_cols_not_present = maf_req_cols %>% 
    "["(! . %in% names(maf))
  if( length(maf_req_cols_not_present) > 0 ){
    k = paste(maf_req_cols_not_present, collapse = ", ") %>% 
      gettextf("Required columns missing in the input MAF data frame: %s.", .)
    stop(k)
  }
  
  # check for optional columns in maf
  maf_opt_cols = c("Mutation_Status", "t_alt_count", "t_depth", "n_alt_count", "n_depth", 
                   "Tumor_Sample_Barcode", "Reference_Allele")
  maf_opt_cols_not_present = maf_opt_cols %>% 
    "["(! . %in% names(maf))
  if( length(maf_opt_cols_not_present) > 0 ){
    k = paste(maf_opt_cols_not_present, collapse = ", ") %>% 
      gettextf("Warning: Otional columns missing in the input MAF data frame: %s.", .)
    message(k)
  }
  
  # add metadata columns to maf
  metadata = get_gambl_metadata(seq_type_filter = c("genome", "capture", "mrna"))
  maf = left_join(maf, metadata, by = join_by(Tumor_Sample_Barcode == sample_id))
  
  # check for optional columns in maf object (including columns from the metadata)
  maf_opt_cols = c("pathology", "patient_id", "time_point")
  maf_opt_cols_not_present = maf_opt_cols %>% 
    "["(! . %in% names(maf))
  if( length(maf_opt_cols_not_present) > 0 ){
    k = paste(maf_opt_cols_not_present, collapse = ", ") %>% 
      gettextf("Warning: Otional columns missing in the metadata (created internally): %s.", .)
    message(k)
  }
  
  # adequate maf to ProteinPaint visualization format
  maf = rename(maf, gene = Hugo_Symbol, chromosome = Chromosome, 
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
  maf$refseq = maf$refseq %>% 
    gsub("\\.[0-9]+,", ",", .) %>% 
    gsub("\\.[0-9]+$", "", .)
  
  # if alt allele is heterogeneous alternative, show both allele separated 
  # by comma (like in VCF files)
  is_het_alt = ! maf$REF == maf$Tumor_Seq_Allele1
  maf$ALT = is_het_alt %>% 
    { paste( maf$Tumor_Seq_Allele1[.], maf$Tumor_Seq_Allele2[.], sep = "," ) } %>% 
    replace(maf$Tumor_Seq_Allele2, is_het_alt, .)
  
  # add VAF (variant allele fraction)
  maf = mutate(maf, VAF = mutant_in_tumor/total_in_tumor)
  
  # keep only columns that are meaningful for proteinpaint
  maf = dplyr::select(maf, gene, refseq, chromosome, start, aachange, class, disease, origin,
                      patient, sample, sampletype, mutant_in_tumor, total_in_tumor, 
                      mutant_in_normal, total_in_normal, REF, ALT, VAF)
  
  # filter out rows that don't contain required column values and output warning message
  required_cols = dplyr::select(maf, gene, refseq, chromosome, start, aachange, class) %>% 
    mutate(i = row_number()) %>% 
    tidyr::drop_na()
  required_cols = apply(required_cols != "", 1, all) %>% 
    dplyr::filter(required_cols, .)
  if(print_removed_rows){  # if only removed rows should be returned (for checking purpose)
    removed_rows = 1:nrow(maf) %>% 
      "["(! . %in% required_cols$i) %>% 
      slice(maf, .)
    return(removed_rows)
  }
  num_removed_rows = (nrow(maf) - nrow(required_cols)) %>% 
    format(big.mark="'")
  k = format( nrow(maf), big.mark="'" ) %>% 
    gettextf("Warning: %s rows out of %s were removed from the output table because there were missing values in required columns. Run `ssm_to_proteinpaint` again with `print_removed_rows = TRUE` to see these rows.",
             num_removed_rows, .)
  message(k)
  maf = slice(maf, required_cols$i)
  
  # return final output table
  maf = as.data.frame(maf)
  maf
}
