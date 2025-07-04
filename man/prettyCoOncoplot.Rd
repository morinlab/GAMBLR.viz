% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prettyCoOncoplot.R
\name{prettyCoOncoplot}
\alias{prettyCoOncoplot}
\title{Side-by-side Oncoplots}
\usage{
prettyCoOncoplot(
  maf,
  metadata,
  comparison_column,
  comparison_values,
  labels,
  genes,
  keepGeneOrder = FALSE,
  returnEverything = FALSE,
  individualOncoplotArgs,
  ...
)
}
\arguments{
\item{maf}{Required argument. A data frame containing the mutations you want
to plot on both oncoplots.}

\item{metadata}{Required argument. A data.frame with metadata for both
oncoplots.}

\item{comparison_column}{Required: the name of the metadata column containing
the comparison values.}

\item{comparison_values}{Optional: If the comparison column contains more
than two values or is not a factor, specify a character vector of length
two in the order you would like the factor levels to be set, reference
group first.}

\item{labels}{Optional: A named vector with the value representing the name
you want displayed on the top of each oncoplot. The names of the vector
     should match the values in the `comparison_column` of the metadata.}

\item{keepGeneOrder}{Optional: If TRUE, the genes will be displayed
in the order they are found in the `genes` argument.}

\item{returnEverything}{Optional: If TRUE, the function will return a list of
the inputs used to generate the oncoplots.
Useful for subsequent functions that require the same inputs.}

\item{individualOncoplotArgs}{Optional: A named list of arguments to be
passed to [GAMBLR.viz::prettyOncoplot] for customizing each oncoplot.
The names of the list must match at least one of the values in the 
`comparison_column` of the metadata.}

\item{...}{`prettyOncoplot` arguments, see that function for more info on
avaialble parameters.}
}
\value{
A complexheatmap object with 2 oncoplots side-by-side.
}
\description{
`prettyCoOncoplot` returns figure of an arbitrary number of
[GAMBLR.viz::prettyOncoplot] images placed side-by-side with the same gene order.
}
\details{
This function will generate a graphic displaying 2 or more oncoplots
side-by-side. Optionally user can annotate each oncoplot with it's own title
that will be displayed at the top via the labels argument.
All the arguments recognized by
[GAMBLR.viz::prettyOncoplot] are supported and can be specified when calling
this function. For both oncoplots the same specified parameters will be
applied (e.g. genes to display, split columns, font size, top annotation
etc). When the provided argument is not recognized by
[GAMBLR.viz::prettyOncoplot], it will be ignored. Some features of the oncoplots
will be dropped from all but the right-most oncoplot (e.g. gene names).
The order of the oncoplots is defined (left-to-right) by either the order of
the factors in the `comparison_column` when it's a factor, the order of
appearance of the values in the `comparison_column`, or by the order of
elements in the `comparison_values` when the comparison column contains more
than 2 levels. The left oncoplot will dictate the order of genes to be
displayed on the right oncoplot. If you want a specific order of oncoplots on
the left and right, please ensure the argument `comparison_column` is a
factor with first level being the group you want to be plotted on the left
side.
For developers: new arguments added to [GAMBLR.viz::prettyOncoplot] in the
future are expected to be out-of-the-box compatible with this function and
would not need code modifications in this function.
}
\examples{
cat("Running example for function: prettyCoOncoplot\n")
library(GAMBLR.open)
#get data for plotting
meta <- get_gambl_metadata()
meta <- meta \%>\%
    dplyr::filter(
        study == "FL_Dreval",
        pathology \%in\% c("DLBCL", "FL")
    )
ssm <- get_coding_ssm(
    these_samples_metadata = meta
)

#build plot
prettyCoOncoplot(
    maf = ssm,
    metadata = meta,
    comparison_column = "pathology",
    genes=dplyr::filter(lymphoma_genes,
                        FL_Tier==1 | DLBCL_Tier==1) \%>\% 
                        dplyr::pull(Gene),
    metadataColumns = c(
        "pathology",
        "lymphgen",
        "pairing_status"
    ),
    metadataBarHeight = 10,
    fontSizeGene = 12,
    metadataBarFontsize = 10,
    legend_row = 2,
    label1 = "DLBCL",
    label2 = "FL",
    simplify_annotation =TRUE,
    minMutationPercent = 5
)

}
