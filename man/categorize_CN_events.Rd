% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/circular_CN_plot.R
\name{categorize_CN_events}
\alias{categorize_CN_events}
\title{Categorize arm-level and chromosomal CNV events}
\usage{
categorize_CN_events(pretty_CN_heatmap_output)
}
\arguments{
\item{pretty_CN_heatmap_output}{The output from running the pretty_CN_heatmap function}
}
\value{
List of data frames
}
\description{
Categorize arm-level and chromosomal CNV events
}
\examples{
\dontrun{
library(GAMBLR.open)
meta = get_gambl_metadata()
meta = check_and_clean_metadata(meta,duplicate_action="keep_first")

all_segments = get_cn_segments(these_samples_metadata = meta)

all_states_binned = segmented_data_to_cn_matrix(
                                 seg_data = all_segments,
                                 strategy="auto_split",
                                 n_bins_split=500,
                                 these_samples_metadata = meta)
labelTheseGenes = c("REL","TP53")
cn_out = pretty_CN_heatmap(cn_state_matrix=all_states_binned,
    scale_by_sample = TRUE,
    these_samples_metadata = meta,
    metadataColumns = c("pathology","seq_type"),
    return_data = TRUE)
    
aneuploidies = categorize_CN_events(cn_out)

select(aneuploidies,1:5) \%>\% head()

}
}
\keyword{internal}
