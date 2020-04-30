#' Range_creator - creating a ranges of exons/cds of a gene
#'
#' This function creating a ranges of exons/cds of a gene which is needed to to find the aa in find_aa
#' @param  geneid geneid can be knownGene name (UCSC) or geneSymbole
#' @param bed bed file
#' @return  DF with the gene's range
#' @export
Range_creator <- function(geneid,bed,exonsORcds){
  if (exonsORcds=="exons") {
    if (sum(bed$hg19.knownGene.name==geneid)>0) {
      indx <- which(bed$hg19.knownGene.name==geneid)
      start <- as.integer(na.omit(as.integer(str_split(bed$hg19.knownGene.exonStarts[indx],",")[[1]])))
      end <- as.integer(na.omit(as.integer(str_split(bed$hg19.knownGene.exonEnds[indx],",")[[1]])))
      df <- as.data.frame(IRanges(start,end))
      df$cumsum <- cumsum(df[,3])
      return(df)
    }else {if (sum(bed$hg19.kgXref.geneSymbol==geneid)>0) {
      indx <- which(bed$hg19.kgXref.geneSymbol==geneid)
      start <- as.integer(na.omit(as.integer(str_split(bed$hg19.knownGene.exonStarts[indx],",")[[1]])))
      end <- as.integer(na.omit(as.integer(str_split(bed$hg19.knownGene.exonEnds[indx],",")[[1]])))
      df <- as.data.frame(IRanges(start,end))
      df$cumsum <- cumsum(df[,3])
      return(df)

    } else{ return(print("Gene Not Found"))}}


  }else {if (exonsORcds=="cds") {
    if (sum(bed$hg19.knownGene.name==geneid)>0) {
      indx <- which(bed$hg19.knownGene.name==geneid)
      start <- as.integer(na.omit(as.integer(str_split(bed$hg19.knownGene.cdsStart[indx],",")[[1]])))
      end <- as.integer(na.omit(as.integer(str_split(bed$hg19.knownGene.cdsEnd[indx],",")[[1]])))
      df <- as.data.frame(IRanges(start,end))
      df$cumsum <- cumsum(df[,3])
      return(df)
    }else {if (sum(bed$hg19.kgXref.geneSymbol==geneid)>0) {
      indx <- which(bed$hg19.kgXref.geneSymbol==geneid)
      start <- as.integer(na.omit(as.integer(str_split(bed$hg19.knownGene.cdsStart[indx],",")[[1]])))
      end <- as.integer(na.omit(as.integer(str_split(bed$hg19.knownGene.cdsEnd[indx],",")[[1]])))
      df <- as.data.frame(IRanges(start,end))
      df$cumsum <- cumsum(df[,3])
      return(df)

    } else{ return(print("Gene Not Found"))}}

  }}


}



#' cds_range - creating a ranges of cds of a gene
#'
#' This function creating a ranges of cds of a gene which is needed to to find the aa in AA_finder
#' @param  rangeCreator.o range creator output
#' @param CDS
#' @param cdsStart start location of gene
#' @param cdsEnd end location of gene
#' @param  geneid geneid can be knownGene name (UCSC) or geneSymbole
#' @return  DF with the gene's range
#' @export
cds_range <- function(rangeCreator.o,CDS,cdsStart,cdsEnd,geneid){
  if (sum(CDS[CDS$geneSymbol==geneid,"strand"]=="-", CDS[CDS$UCSCid==geneid,"strand"]=="-")>0) {
    rangeCreator.o <- rangeCreator.o[seq(dim(rangeCreator.o)[1],1),]
    bigger <- which(rangeCreator.o$start>cdsEnd)
    if (length(bigger)==0) {
      rangeCreator.o[1,2] <- cdsEnd
      rangeCreator.o[,1] <- rangeCreator.o[,1]+1
      cdsRange <- as.data.frame(IRanges(rangeCreator.o[,1],rangeCreator.o[,2]))
      cdsRange$cumsum <- cumsum(cdsRange[,3])
      return(cdsRange)
    }  else {  rangeCreator.o[1:(length(bigger)+1),2] <- cdsEnd
    rangeCreator.o[,1] <- rangeCreator.o[,1]+1
    cdsRange <- as.data.frame(IRanges(rangeCreator.o[((length(bigger)+1):nrow(rangeCreator.o)),1],rangeCreator.o[((length(bigger)+1):nrow(rangeCreator.o)),2]))
    cdsRange$cumsum <- cumsum(cdsRange[,3])
    return(cdsRange)}

  }else {  smaller <- which(rangeCreator.o$end<cdsStart)
  if (length(smaller)==0) {
    rangeCreator.o[1,1] <- cdsStart
    rangeCreator.o[,1] <- rangeCreator.o[,1]+1
    cdsRange <- as.data.frame(IRanges(rangeCreator.o[,1],rangeCreator.o[,2]))
    cdsRange$cumsum <- cumsum(cdsRange[,3])
    return(cdsRange)
  } else {  rangeCreator.o[1:(length(smaller)+1),1] <- cdsStart
  rangeCreator.o[,1] <- rangeCreator.o[,1]+1
  cdsRange <- as.data.frame(IRanges(rangeCreator.o[((length(smaller)+1):nrow(rangeCreator.o)),1],rangeCreator.o[((length(smaller)+1):nrow(rangeCreator.o)),2]))
  cdsRange$cumsum <- cumsum(cdsRange[,3])
  return(cdsRange)}

  }}
