#' @title Map Metadata to Colours
#'
#' @description Assign a colour palette to metadata columns automatically and consistently.
#' 
#' @details Internally called by `plot_sample_circos` in GAMBLR.viz
#'
#' @param metadataColumns Names of the metadata columns to assign colours for.
#' @param these_samples_metadata Metadata for just the samples you need colours for.
#' @param column_alias A list of column_names with values indicating what gambl colour scheme they should use (e.g. pos_neg, pathology, lymphgen).
#' @param as_vector Boolean statement that is set to TRUE per default.
#' @param verbose Set to TRUE to enable verbose mode (debugging messages.
#' @param annoAlpha Optional alpha to apply to annotation colours.
#'
#' @return Either a vector or list of colours.
#'
#' @import dplyr ggsci
#' 
#' @export
#'
#' @examples
#' #get colours
#' all_cols = map_metadata_to_colours(metadataColumns = c("lymphgen",
#'                                                        "pathology",
#'                                                        "genetic_subgroup"),
#'                                    these_samples_metadata = GAMBLR.data::gambl_metadata,
#'                                    column_alias = list("nothing" = "FL"),
#'                                    as_vector = F)
#'
map_metadata_to_colours = function(metadataColumns,
                                   these_samples_metadata,
                                   column_alias = list(),
                                   as_vector = TRUE,
                                   verbose = FALSE,
                                   annoAlpha = 1){
  
  clinical_colours = ggsci::get_ash("clinical")
  all_gambl_colours = GAMBLR.viz::get_gambl_colours()
  colours = list()
  colvec = c()
  
  aliases = c(colour_aliases, column_alias)
  for(column in metadataColumns){
    this_value = these_samples_metadata[[column]]
    options = this_value
    if(verbose){
      print(">>>>>>>")
      message("finding colour for", this_value)
      print("<<<<<<<")
    }
    if(column %in% names(aliases)){
      key = aliases[[column]]
      if(verbose){
        print(paste("using alias to look up colours for", column))
        message(paste("using", key, "for", column))
      }
      these = GAMBLR.viz::get_gambl_colours(classification = key)
      colours[[column]] = these
      colvec = c(colvec, these[this_value])
      if(verbose){
        message("adding:", these[this_value])
      }
    }else if(column == "sex"){
      these = GAMBLR.viz::get_gambl_colours("sex", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these= c(these, "NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these[this_value])
      message("adding:", these[this_value])
    }else if(sum(levels(options) %in% names(clinical_colours)) == length(levels(options))){
      
      #we have a way to map these all to colours!
      if(verbose){
        message(paste("found colours for", column, "in clinical"))
      }
      these = clinical_colours[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these,"NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these[this_value])
    }else if(("positive" %in% options | "POS" %in% options) & length(options)<4){
      if(verbose){
        print("using pos_neg")
      }
      
      these = GAMBLR.viz::get_gambl_colours("pos_neg", alpha = annoAlpha)
      these = these[levels(options)]
      
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these[this_value])
    }else if("GCB" %in% options){
      these = GAMBLR.viz::get_gambl_colours("COO", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these,"NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec,these)
    }else if(column %in% c("pathology")){
      these = GAMBLR.viz::get_gambl_colours(column, alpha = annoAlpha)
      
      if(!"NA" %in% names(these)){
        these = c(these,"NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these)
    }else if(grepl("lymphgen", column, ignore.case = TRUE)){
      these = GAMBLR.viz::get_gambl_colours("lymphgen", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these)
    }else if(column == "HMRN"){
      these = GAMBLR.viz::get_gambl_colours("hmrn", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these= c(these,"NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these)
    }else if(sum(levels(options) %in% names(all_gambl_colours)) == length(levels(options))){
      if(verbose){
        message(paste("found colours for", column, "in all_gambl_colours"))
      }
      these = all_gambl_colours[levels(options)]
      if(!"NA" %in% names(these)){
        these= c(these,"NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these)
    }else if(length(levels(options)) > 15){
      
      these = rainbow(length(levels(options)), alpha = annoAlpha)
      names(these) = levels(options)
      
      colours[[column]] = these
      colvec = c(colvec, these)
    }else{
      these = blood_cols[sample(c(1:length(blood_cols)), size = length(levels(options)))]
      names(these) = levels(options)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these)
    }
  }
  if(as_vector){
    return(colvec)
  }
  return(colours)
}
