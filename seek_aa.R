
data("CDS")
data("bed")
data("codonTable")



#' AA_finder - find aa genomic coordinates.
#'
#' The function finds the genomic coordinates of aa in a specified gene.
#' @param  aa amino acid. character of one capital letter.
#' @param i integer - which one is wanted? the second tryptophan(W)? the 4th glutamine (W) ?
#' @param  geneid geneid can be knownGene name (UCSC) or geneSymbole
#' @param bed DF - the bed file
#' @param CDSdb cds df
#' @param  exonsORcds takes "exons" or "cds".
#' @param genesymORucsc "genesym" or "ucsc" used geneid
#' @return  3 strings. start location, end location and chromosome.
#' @export
AA_finder <- function(aa,i,geneid,bed,CDSdb,exonsORcds,genesymORucsc) {
  if (sum(colnames(CDS)==paste(aa,"count",sep = ""))==0) {
    CDS <<- Add_AAtoCDS(aa,CDSdb)
    AA_finder(aa,i,geneid,bed,CDS,exonsORcds,genesymORucsc)

  }
  CDS <- CDS[CDS[paste(aa,"count",sep = "")]>0,]
  if (genesymORucsc == "genesym") {
    chr <- as.character(bed[bed$hg19.kgXref.geneSymbol==geneid,which(colnames(bed)=="hg19.knownGene.chrom")])[1]
    if (geneid %in% as.character(CDS$geneSymbol)) {
      pos <- str_split(as.character(CDS[CDS$geneSymbol==geneid,paste(aa,"pos",sep = "")]),pattern = ",")[[1]] %>% as.integer()
      rangeCreator.o <- Range_creator(geneid,bed,exonsORcds)
      cdsStart <- bed[bed$hg19.kgXref.geneSymbol==geneid,which(colnames(bed)=="hg19.knownGene.cdsStart")][1]
      cdsEnd <- bed[bed$hg19.kgXref.geneSymbol==geneid,which(colnames(bed)=="hg19.knownGene.cdsEnd")][1]
      cdsRange <- cds_range(rangeCreator.o,CDS,cdsStart,cdsEnd,geneid)

      if (CDS[CDS$geneSymbol==geneid,"strand"][1]=="-") {
        if (which(cdsRange$cumsum>(pos[i]*3))[1]==1) {
          firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]-((pos[i]*3)-3)-2 # // 6
          lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]-((pos[i]*3)-3) # // 6
          if (firstnt<cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]) {
            realfirst <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]
            diff <- realfirst-firstnt
            firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1]-1,2]-(diff+1)
            return(c(firstnt,lasttnt,chr))}
          return(c(firstnt,lasttnt,chr))

        }else {
          dis <- pos[i]*3-cdsRange[(which(cdsRange$cumsum>(pos[i]*3))[1]-1),4] # this was done in cases that W is at the beginning of an exome
          if (dis<3) {
            if (dis==0) {
              lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1]-1,1]+2
              firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1]-1,1]
              return(c(firstnt,lasttnt,chr))

            }else {
              lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1]-1,1]+(3-dis-1)
              firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]-((pos[i]*3)-(cdsRange[(which(cdsRange$cumsum>(pos[i]*3))[1]-1),4]+3))-2
              return(c(firstnt,lasttnt,chr))
            }


          }else {
            lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]-((pos[i]*3)-(cdsRange[(which(cdsRange$cumsum>(pos[i]*3))[1]-1),4]+3)) # // 4
            firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]-((pos[i]*3)-(cdsRange[(which(cdsRange$cumsum>(pos[i]*3))[1]-1),4]+3))-2 # // 4
            if (firstnt<cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]) {
              realfirst <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]
              diff <- realfirst-firstnt
              firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1]-1,2]-(diff+1)
              return(c(firstnt,lasttnt,chr))}


          }
          return(c(firstnt,lasttnt,chr))

        }



      }else{  if (which(cdsRange$cumsum>(pos[i]*3))[1]==1) {
        firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]+(pos[i]*3)-3 #this is always get messed up // 2,3,8
        lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]+(pos[i]*3)-1#+2 // 2,3,8
        if (lasttnt>cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]) {
          reallast <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]
          diff <- lasttnt-reallast
          lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1]+1,1]+(diff-1)
          return(c(firstnt,lasttnt,chr))}
        return(c(firstnt,lasttnt,chr))

      }else {
        firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]+((pos[i]*3)-cdsRange[(which(cdsRange$cumsum>(pos[i]*3))[1]-1),4])-3 # also
        lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]+(pos[i]*3-cdsRange[(which(cdsRange$cumsum>(pos[i]*3))[1]-1),4])-1
        if (lasttnt>cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]) {
          reallast <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]
          diff <- lasttnt-reallast
          lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1]+1,1]+(diff-1)
          return(c(firstnt,lasttnt,chr))}
        return(c(firstnt,lasttnt,chr))

      }

      }

    } else {print("Gene not found in CDS DB")}

  } else {
    chr <- as.character(bed[bed$hg19.knownGene.name==geneid,which(colnames(bed)=="hg19.knownGene.chrom")])[1]
    if (geneid %in% as.character(CDS$UCSCid)) {
      pos <- str_split(as.character(CDS[CDS$UCSCid==geneid,paste(aa,"pos",sep = "")]),pattern = ",")[[1]] %>% as.integer()
      rangeCreator.o <- Range_creator(geneid,bed,exonsORcds)
      cdsStart <- bed[bed$hg19.knownGene.name==geneid,which(colnames(bed)=="hg19.knownGene.cdsStart")][1]
      cdsEnd <- bed[bed$hg19.knownGene.name==geneid,which(colnames(bed)=="hg19.knownGene.cdsEnd")][1]
      cdsRange <- cds_range(rangeCreator.o,CDS,cdsStart,cdsEnd,geneid)

      if (CDS[CDS$UCSCid==geneid,"strand"][1]=="-") {
        if (which(cdsRange$cumsum>(pos[i]*3))[1]==1) {
          firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]-((pos[i]*3)-3)-2 # // 6
          lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]-((pos[i]*3)-3) # // 6
          if (firstnt<cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]) {
            realfirst <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]
            diff <- realfirst-W_firstnt
            firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1]-1,2]-(diff+1)
            return(c(firstnt,lasttnt,chr))}
          return(c(firstnt,lasttnt,chr))

        }else {
          dis <- pos[i]*3-cdsRange[(which(cdsRange$cumsum>(pos[i]*3))[1]-1),4] # this was done in cases that W is at the beginning of an exome
          if (dis<3) {
            if (dis==0) {
              lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1]-1,1]+2
              firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1]-1,1]
              return(c(firstnt,lasttnt,chr))

            }else {
              lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1]-1,1]+(3-dis-1)
              firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]-((pos[i]*3)-(cdsRange[(which(cdsRange$cumsum>(pos[i]*3))[1]-1),4]+3))-2
              return(c(firstnt,lasttnt,chr))
            }


          }else {
            lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]-((pos[i]*3)-(cdsRange[(which(cdsRange$cumsum>(pos[i]*3))[1]-1),4]+3)) # // 4
            firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]-((pos[i]*3)-(cdsRange[(which(cdsRange$cumsum>(pos[i]*3))[1]-1),4]+3))-2 # // 4
            if (firstnt<cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]) {
              realfirst <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]
              diff <- realfirst-firstnt
              firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1]-1,2]-(diff+1)
              return(c(firstnt,lasttnt,chr))}


          }
          return(c(firstnt,lasttnt,chr))

        }



      }else{  if (which(cdsRange$cumsum>(pos[i]*3))[1]==1) {
        firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]+(pos[i]*3)-3 #this is always get messed up // 2,3,8
        lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]+(pos[i]*3)-1#+2 // 2,3,8
        if (lasttnt>cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]) {
          reallast <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]
          diff <- lasttnt-reallast
          lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1]+1,1]+(diff-1)
          return(c(firstnt,lasttnt,chr))}
        return(c(firstnt,lasttnt,chr))

      }else {
        firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]+((pos[i]*3)-cdsRange[(which(cdsRange$cumsum>(pos[i]*3))[1]-1),4])-3 # also
        lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]+(pos[i]*3-cdsRange[(which(cdsRange$cumsum>(pos[i]*3))[1]-1),4])-1
        if (lasttnt>cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]) {
          reallast <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],2]
          diff <- lasttnt-reallast
          lasttnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1]+1,1]+(diff-1)
          return(c(firstnt,lasttnt,chr))}
        if (firstnt<cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]) {
          diff <- lasttnt-cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1],1]
          diff <- 3-diff-2
          firstnt <- cdsRange[which(cdsRange$cumsum>(pos[i]*3))[1]-1,2]-diff

        }
        return(c(firstnt,lasttnt,chr))

      }

      }

    } else {print("Gene not found in CDS DB")}}


}

# a <- AA_finder(aa = "W",i = 6,
#                geneid = "CD44",
#                bed = bed,
#                CDSdb = CDS,
#                exonsORcds = "exons",
#                genesymORucsc = "genesym")
# writeClipboard(paste(a[3],a[2],sep = ":"))



#' CODON_finder - find the codon of the specific amino acid wanted.
#'
#' TThe function finds the codon of the specific amino acid wanted.
#' @param  aa amino acid. character of one capital letter.
#' @param i integer - which one is wanted? the second tryptophan(W)? the 4th glutamine (W) ?
#' @param genesymORucsc "genesym" or "ucsc" used geneid
#' @return  list of 3 nts
#' @export
CODON_finder <- function(aa,i,genesymORucsc){
  if (genesymORucsc == "genesym"){
    posnt <- str_split(as.character(CDS[CDS$geneSymbol==geneid,paste(aa,"posnt",sep = "")]),pattern = ",")[[1]][i] %>% as.integer()
    mid <- s2c(as.character(CDS[CDS$geneSymbol==geneid,"seq"]))[posnt]
    fir <- s2c(as.character(CDS[CDS$geneSymbol==geneid,"seq"]))[posnt-1]
    las <- s2c(as.character(CDS[CDS$geneSymbol==geneid,"seq"]))[posnt+1]
    return(list(middle=mid,first=fir,last=las))

  } else {
    posnt <- str_split(as.character(CDS[CDS$UCSCid==geneid,paste(aa,"posnt",sep = "")]),pattern = ",")[[1]][i] %>% as.integer()
    mid <- s2c(as.character(CDS[CDS$geneSymbol==geneid,"seq"]))[posnt]
    fir <- s2c(as.character(CDS[CDS$geneSymbol==geneid,"seq"]))[posnt-1]
    las <- s2c(as.character(CDS[CDS$geneSymbol==geneid,"seq"]))[posnt+1]
    return(list(middle=mid,first=fir,last=las))

  }

}
# unlist(CODON_finder(aa,i,genesymORucsc))













