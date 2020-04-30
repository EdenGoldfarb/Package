#' Add_AAtoCDS - Add aa to the CDS df.
#'
#' This function add aa information to the CDS df. it adds the following data:
#' 1. aa count ; 2.aa positions; 3.nts positions
#' @param  aa amino acid. one capital letter.
#' @param CDS_DB cds df
#' @return  A new dataframe with the add aa columns
#' @export
Add_AAtoCDS <- function(aa,CDS_DB){
  CDS_DB[paste(aa,"count",sep = "")] <- str_count(CDS_DB$AA,aa)
  CDS_DB[paste(aa,"pos",sep = "")]  <- as.data.frame(mapply(function(i) {paste(str_locate_all(i,pattern = aa)[[1]][,1],collapse = ",")}, CDS_DB$AA))[,1]
  CDS_DB[paste(aa,"posnt",sep = "")] <- mapply(function(i) {paste(as.character(as.integer(str_split(i,",")[[1]])*3-1),collapse = ",")}  , CDS_DB[,paste(aa,"pos",sep = "")])

  return(CDS_DB)
}
