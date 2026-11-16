# Internal: subset a user-supplied MAF to regions and samples, returning the
# streamlined columns the rainbow plots use (start, sample_id, region_name).
# Replaces calls to get_ssm_by_regions(maf_data = ...), which GAMBLR.viz does
# not import. Chromosome names are compared without any "chr" prefix.
# regions_df needs columns chrom, start, end, name. Rows with fewer than
# min_read_support alt reads (t_alt_count) are dropped, matching the default
# of get_ssm_by_regions().
.subset_maf_to_regions <- function(maf_data, regions_df, sample_ids, min_read_support = 3) {
    maf <- as.data.frame(maf_data)
    needed <- c("Chromosome", "Start_Position", "Tumor_Sample_Barcode")
    if (!all(needed %in% colnames(maf))) {
        stop("maf_data must contain the columns ", paste(needed, collapse = ", "))
    }
    if (min_read_support > 0 && "t_alt_count" %in% colnames(maf)) {
        maf <- maf[!is.na(maf$t_alt_count) & maf$t_alt_count >= min_read_support, ]
    }
    maf <- maf[maf$Tumor_Sample_Barcode %in% sample_ids, ]
    maf_chrom <- gsub("^chr", "", as.character(maf$Chromosome))
    hits <- lapply(seq_len(nrow(regions_df)), function(i) {
        keep <- maf_chrom == gsub("^chr", "", as.character(regions_df$chrom[i])) &
            maf$Start_Position >= regions_df$start[i] &
            maf$Start_Position <= regions_df$end[i]
        if (!any(keep)) return(NULL)
        data.frame(start = maf$Start_Position[keep],
                   sample_id = maf$Tumor_Sample_Barcode[keep],
                   region_name = regions_df$name[i])
    })
    out <- do.call(rbind, hits)
    if (is.null(out)) out <- data.frame(start = numeric(0), sample_id = character(0), region_name = character(0))
    out
}
