#' @title Assign a colour palette to metadata columns automatically and consistently.
#'
#' @param metadataColumns Names of the metadata columns to assign colours for.
#' @param these_samples_metadata Metadata for just the samples you need colours for.
#' @param column_alias A list of column_names with values indicating what gambl colour scheme they should use (e.g. pos_neg, pathology, lymphgen).
#' @param verbose Set to TRUE to enable verbose mode (debugging messages.
#' @param annoAlpha Optional alpha to apply to annotation colours.
#'
#' @return Either a vector or list of colours.
#'
#' @import dplyr
#'
#' @examples
#' library(GAMBLR.data)
#' 
#' #get metadata
#' all_meta = suppressMessages(get_gambl_metadata())
#'
#' #get colours
#' all_cols = map_metadata_to_colours(
#'      metadataColumns = c(
#'          "lymphgen",
#'          "pathology",
#'          "genetic_subgroup"),
#'      these_samples_metadata = all_meta,
#'      column_alias = list("nothing" = "FL")
#' )
#'
#' @export
map_metadata_to_colours = function(metadataColumns,
                                   these_samples_metadata,
                                   column_alias = list(),
                                   verbose = FALSE,
                                   annoAlpha = 1){

  #automagically assign colours for other metadata columns.
  blood_cols = GAMBLR.helpers::get_gambl_colours("blood", alpha = annoAlpha)
  colours = list()
  clinical_colours = GAMBLR.helpers::get_gambl_colours("clinical")
  all_gambl_colours = GAMBLR.helpers::get_gambl_colours()
  for(column in metadataColumns){
    these_samples_metadata[[column]] = factor(these_samples_metadata[[column]], levels = unique(these_samples_metadata[[column]]))
    options = these_samples_metadata %>%
      arrange(column) %>%
      dplyr::filter(!is.na(column)) %>%
      pull(column) %>%
      unique()
    
    options = options[!is.na(options)]
    if(verbose){
      print(">>>XX>>>")
      print(levels(options))
      print("<<<XX<<<")
    }
    if(grepl("^chr",column)){
      these = GAMBLR.helpers::get_gambl_colours("aneuploidy", alpha = annoAlpha)
      these = these[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }
    if(column == "sex"){
      these = GAMBLR.helpers::get_gambl_colours("sex", alpha = annoAlpha)
      these = these[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(sum(levels(options) %in% names(clinical_colours)) == length(levels(options))){
      #we have a way to map these all to colours!
      if(verbose){
        message(paste("found colours for", column, "here"))
      }
      these = clinical_colours[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(("positive" %in% options | "POS" %in% options | "yes" %in% options) & length(options) < 4){
      if(verbose){
        print("using pos_neg")
      }
      these = GAMBLR.helpers::get_gambl_colours("pos_neg", alpha = annoAlpha)
      these = these[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if("GCB" %in% options){
      these = GAMBLR.helpers::get_gambl_colours("COO", alpha = annoAlpha)
      these = these[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(column %in% c("pathology")){
      these = GAMBLR.helpers::get_gambl_colours(column, alpha = annoAlpha)
      these = these[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these,"NA" = "white")
      }
      colours[[column]] = these
    }else if(grepl("lymphgen", column, ignore.case = TRUE)){
      these = GAMBLR.helpers::get_gambl_colours("lymphgen", alpha = annoAlpha)
      these = these[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
      
    }else if(column == "HMRN"){
      these = GAMBLR.helpers::get_gambl_colours("hmrn", alpha = annoAlpha)
      these = these[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(sum(levels(options) %in% names(all_gambl_colours)) == length(levels(options))){
      if(verbose){
        message(paste("found colours for", column, "in all_gambl_colours"))
      }
      these = all_gambl_colours[levels(options)]
      these = these[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(length(levels(options)) > 15){
      these = rainbow(length(levels(options)), alpha = annoAlpha)
      names(these) = levels(options)
      colours[[column]] = these
    }else{
      these = blood_cols[sample(c(1:length(blood_cols)), size = length(levels(options)))]
      names(these) = levels(options)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      
      colours[[column]] = these
    }
    #check that our colour set is complete
    
    loptions = levels(options)
    if(any(!loptions %in% names(these))){
      #try the full set of colours  
      if(!any(!loptions %in% names(all_gambl_colours))){
        colours[[column]] = all_gambl_colours[loptions]
      }else{
        missing = options[!loptions %in% names(all_gambl_colours)]
        #try to replace any "/"-containing strings with "-COMP"
        missing_cols = c()
        colours[[column]] = all_gambl_colours[loptions[loptions %in% names(all_gambl_colours)]]
        if(verbose){
          print("XXXXX")
          print(colours[[column]])
        }
        for(colname in missing){
          if(grepl("/",colname)){
            prefix = str_remove(colname,"/.+")
            if(verbose){
              print(paste("trying to map:",colname,"using",prefix))
            }
            composite = paste0(prefix,"-COMP")
            if(composite %in% names(all_gambl_colours)){
              if(verbose){
                print(paste("Success",all_gambl_colours[composite]))
              }
              this_col = all_gambl_colours[composite]
              names(this_col)=colname
              missing_cols = c(missing_cols,this_col)
            }else{
              if(verbose){
                print(paste("FAIL",composite,all_gambl_colours[composite]))
              }
            }
          }
        }
        if(verbose){
          print(missing_cols)
          print(colours[[column]])
        }
        if(length(missing_cols)< length(missing)){
          #last-ditch effort (random)
          fill_in_cols = rainbow(length(missing))
          names(fill_in_cols) = missing
          colours[[column]] = c(colours[[column]],fill_in_cols)
        }else{
          colours[[column]] = c(colours[[column]],missing_cols)
        }
      }
    }
    if(!"NA" %in% names(colours[[column]])){
      colours[[column]] = c(colours[[column]], "NA" = "white")
    }
  }
  return(colours)
}
