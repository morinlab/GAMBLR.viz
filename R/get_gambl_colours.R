#' @title Get GAMBL Colours.
#'
#' @description Get GAMBL colour schemes for annotating figures.
#'
#' @details This function was designed to retrieve specified GAMBL colour palettes.
#' By default, this function returns all the colours currently available.
#' The user can easily specify what classification to return colors for with the `classification` parameter.
#' It is also possible to return any given colour in different formats.
#' To do so, refer to the Boolean arguments; `as_list` and `as_dataframe`.
#' For more information regarding the available colours, refer to the utilities vignette.
#'
#' @param classification Optionally request only colours for pathology, lymphgen, mutation or copy_number.
#' @param alpha Alpha of plotted colours.
#' @param as_list Boolean parameter controlling the format of the return. Default is FALSE.
#' @param as_dataframe Boolean parameter controlling the format of the return. Default is FALSE.
#' @param return_available Set to TRUE for returning all available colours. Default is FALSE.
#' @param verbose Default is FALSE
#'
#' @return A named vector of colour codes for lymphgen classes and pathology.
#'
#' @import dplyr ggsci stringr tidyr
#' @export
#'
#' @examples
#' lymphgen_cols = get_gambl_colours("lymphgen")
#'
#' \dontrun{
#' #be sure to install ggsci from https://github.com/morinlab/ggsci
#' #install_github("morinlab/ggsci")
#' }
#'
get_gambl_colours = function(classification = "all",
                             alpha = 1,
                             as_list = FALSE,
                             as_dataframe = FALSE,
                             return_available = FALSE,
                             verbose = FALSE){
  
  all_colours = list()
  everything = c()
  blood_cols = ggsci::get_ash("blood")
  
  all_colours[["seq_type"]] = c("mrna" = "#E41A1C",
                                "genome" = "#377EB8",
                                "capture" = "#4DAF4A")
  
  all_colours[["type"]] = c("gain" = "blue",
                            "loss" = "red")
  
  all_colours[["hmrn"]] = c("BCL2-MYC" = "#52000F",
                            "BCL2" = "#721F0F",
                            "SOCS1/SGK1" = "#D66B1F",
                            "TET2/SGK1" = "#C41230",
                            "MYD88" = "#3B5FAC",
                            "NOTCH2" = "#7F3293",
                            "NOTCH1" = "#55B55E",
                            "Other" = "#ACADAF")
  
  all_colours[["EBV"]] =  c("EBV-positive" = "#7F055F",
                            "EBV-negative" = "#E5A4CB",
                            "POS" = "#7F055F",
                            "NEG" = "#E5A4CB")
  
  all_colours[["BL"]] = c("Q53-BL" = "#A6CEE3",
                          "M53-BL" = "#A6CEE3", #added because genetic subgroup still refers to it this way
                          "DLBCL-A" = "#721F0F",
                          "IC-BL" = "#45425A",
                          "DGG-BL" = "#E90C8B",
                          "DLBCL-B" = "#FB9A99",
                          "DLBCL-C" = "#C41230")
  
  all_colours[["FL"]] = c(dFL = "#99C1B9", cFL = "#D16666", DLBCL = "#479450")
  
  all_colours[["lymphgenerator"]] = c("MP3"="#5B8565",
                                      "EGB" = "#98622A",
                                      "ETB"="#813F3D",
                                      "aSCI"="#D66B1F",
                                      "aSEL"="#6A0D18",
                                      "MCaP"="#5F8CFF",
                                      "BNZ"="#8870B6",
                                      "EZB"="#721F0F",
                                      "ST2"="#C41230",
                                      "UNCLASS"="#05631E"
  )
  
  all_colours[["chapuy_classifier"]] = c(
    C0 = "#bebebe",
    C1 = "#803D99",
    C2 ="#00A2D2",
    C3 = "#F39123",
    C4 = "#50BFAD",
    C5 = "#DE292A"
  )
  
  all_colours[["lacy_classifier"]] = all_colours[["hmrn"]]
  
  all_colours[["lymphgen"]] = c("EZB-MYC" = "#52000F",
                                "EZB" = "#721F0F",
                                "EZB-COMP" = "#C7371A",
                                "ST2" = "#C41230",
                                "ST2-COMP" = "#EC3251",
                                "MCD" = "#3B5FAC",
                                "MCD-COMP" = "#6787CB",
                                "BN2" =  "#7F3293",
                                "BN2-COMP" = "#A949C1",
                                "N1" = "#55B55E",
                                "N1-COMP" = "#7FC787",
                                "A53" = "#5b6d8a",
                                "Other" = "#ACADAF",
                                "COMPOSITE" = "#ACADAF")
  
  #all_colours[["coding_class"]] = c("Frame_Shift_Del","Frame_Shift_Ins",
  #                 "In_Frame_Del","In_Frame_Ins",
  #                 "Missense_Mutation","Nonsense_Mutation",
  #                 "Nonstop_Mutation","Splice_Region","Splice_Site",
  #                 "Targeted_Region","Translation_Start_Site")
  all_colours[["mutation"]]=
    c(
      "Nonsense_Mutation"=unname(blood_cols["Red"]),
      "Missense_Mutation"=unname(blood_cols["Green"]),
      "Multi_Hit"=unname(blood_cols["Steel Blue"]),
      "Frame_Shift_Ins" = unname(blood_cols["Magenta"]),
      "Frame_Shift_Del" = unname(blood_cols["Magenta"]),
      "In_Frame_Ins" = unname(blood_cols["Brown"]),
      "In_Frame_Del" = unname(blood_cols["Brown"]),
      "Nonstop_Mutation" = unname(blood_cols["Light Blue"]),
      "Translation_Start_Site" = unname(blood_cols["Lavendar"]),
      "Splice_Site" = unname(blood_cols["Orange"]),
      "Splice_Region" = unname(blood_cols["Orange"]),
      "3'UTR" = unname(blood_cols["Yellow"]),
      "Silent" = "#A020F0")
  
  all_colours[["rainfall"]] =
    c(
      "C>A" = "#2196F3FF",
      "C>G" = "#3F51B5FF",
      "C>T" = "#F44336FF",
      "InDel" = "purple",
      "T>A" = "#4CAF50FF",
      "T>C" = "#FFC107FF",
      "T>G" = "#FF9800FF"
    )
  
  all_colours[["pos_neg"]]=c(
    "POS"="#c41230",
    "NEG"="#E88873",
    "PARTIAL"="#E88873",
    "yes"="#c41230",
    "no"="#E88873",
    "YES"="#c41230",
    "NO"="#E88873",
    "FAIL"="#bdbdc1",
    "positive"="#c41230",
    "negative"="#E88873",
    "fail"="#bdbdc1")
  
  all_colours[["copy_number"]]=c(
    "nLOH"="#E026D7",
    "14"="#380015",
    "15"="#380015",
    "13"="#380015",
    "12"="#380015",
    "11"="#380015",
    "10"="#380015",
    "9"="#380015",
    "8"="#380015",
    "7"="#380015",
    "6"="#380015",
    "5"="#67001F",
    "4"="#B2182B",
    "3"="#D6604D",
    "2"="#ede4c7",
    "1"="#92C5DE",
    "0"="#4393C3"
  )
  all_colours[["blood"]] = c(
    "Red" = "#c41230", "Blue"="#115284","Green" = "#39b54b",
    "Purple" = "#5c266c", "Orange"="#fe9003","Green" = "#046852",
    "Lavendar" = "#8781bd", "Steel Blue"= "#455564",
    "Light Blue" = "#2cace3", "Magenta" = "#e90c8b", "Mustard" = "#b76d29",
    "LimeGreen" = "#a4bb87", "Brown" = "#5f3a17", "Gray" = "#bdbdc1",
    "Yellow" = "#f9bd1f"
  )
  all_colours[["sex"]]=c(
    "M"="#118AB2",
    "Male"="#118AB2",
    "male"="#118AB2",
    "F"="#EF476F",
    "Female"="#EF476F",
    "female"="#EF476F")
  all_colours[["clinical"]]=ggsci::get_ash("clinical")
  all_colours[["pathology"]] = c(
    "B-ALL"="#C1C64B",
    "CLL"="#889BE5",
    "MCL"="#F37A20",
    "BL"="#926CAD",
    "mBL"="#34C7F4",
    "tFL"="#FF8595",
    "DLBCL-BL-like"="#34C7F4",
    "pre-HT"="#754F5B",
    "PMBL"= "#227C9D",
    "PMBCL"="#227C9D",
    "FL"="#EA8368",
    "no-HT"="#EA8368",
    "COMFL"="#8BBC98",
    "COM"="#8BBC98",
    "post-HT"="#479450",
    "DLBCL"="#479450",
    "denovo-DLBCL"="#479450",
    "HGBL-NOS"="#294936",
    "HGBL"="#294936",
    "HGBL-DH/TH"="#7A1616",
    "PBL" = "#E058C0",
    "Plasmablastic" = "#E058C0",
    "CNS" = "#E2EF60",
    "THRLBCL" = "#A5F2B3",
    "MM"="#CC9A42",
    "SCBC"="#8c9c90",
    "UNSPECIFIED"="#cfba7c",
    "OTHER"="#cfba7c",
    "MZL"="#065A7F",
    "SMZL"="#065A7F",
    "Prolymphocytic" = "#7842f5"
  )
  all_colours[["coo"]] = c(
    "ABC" = "#05ACEF",
    "UNCLASS" = "#05631E",
    "Unclass" = "#05631E",
    "U" = "#05631E",
    "UNC" = "#05631E",
    "GCB"= "#F58F20",
    "DHITsig-"= "#F58F20",
    "DHITsigNeg"= "#F58F20",
    "DHITsig-IND" = "#003049",
    "DHITsig+" = "#D62828",
    "DHITsigPos" = "#D62828",
    "NA" = "#ACADAF"
  )
  all_colours[["cohort"]] = c("Chapuy"="#8B0000","Chapuy, 2018"="#8B0000",
                              "Arthur"= "#8845A8","Arthur, 2018"= "#8845A8",
                              "Schmitz"= "#2C72B2","Schmitz, 2018"= "#2C72B2",
                              "Reddy" = "#E561C3","Reddy, 2017" = "#E561C3",
                              "Morin"= "#8DB753", "Morin, 2013"= "#8DB753",
                              "Kridel"= "#4686B7", "Kridel, 2016"= "#4686B7",
                              "ICGC"="#E09C3B","ICGC, 2018"="#E09C3B",
                              "Grande"="#e90c8b", "Grande, 2019"="#e90c8b")
  
  all_colours[["indels"]] = c("DEL" = "#53B1FC", "INS" = "#FC9C6D")
  all_colours[["svs"]] = c("DEL" = "#53B1FC", "DUP" = "#FC9C6D")
  all_colours[["genetic_subgroup"]] = c(all_colours[["lymphgen"]],all_colours[["BL"]],all_colours[["FL"]])
  #print(all_colours)
  if(alpha <1){
    for(colslot in names(all_colours)){
      raw_cols = all_colours[[colslot]]
      raw_cols_rgb = col2rgb(raw_cols)
      alpha_cols = rgb(raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ], alpha = alpha * 255L, names = names(raw_cols), maxColorValue = 255L)
      names(alpha_cols) = names(raw_cols)
      all_colours[[colslot]] = alpha_cols
    }
  }
  for(this_group in names(all_colours)){
    everything = c(everything, all_colours[[this_group]])
  }
  #return matching value from lowercase version of the argument if it exists
  lc_class = stringr::str_to_lower(classification)
  if(return_available){
    return(names(all_colours))
  }
  if(classification %in% names(all_colours)){
    if(as_dataframe){
      some_col=all_colours[[classification]]
      df_ugly = data.frame(name=names(some_col),colour=unname(some_col))
      df_tidy = mutate(df_ugly,group=classification)
      return(df_tidy)
    }
    return(all_colours[[classification]])
  }else if(lc_class %in% names(all_colours)){
    return(all_colours[[lc_class]])
  }else if(as_list){
    return(all_colours)
  }else if(as_dataframe){
    df_ugly = data.frame(name = names(unlist(all_colours, use.names = T)), colour = unlist(all_colours, use.names = T))
    df_tidy = separate(df_ugly,name,into=c("group","name"),sep="\\.")
    return(df_tidy)
  }else{
    return(everything)
  }
}
