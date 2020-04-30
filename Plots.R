
#' Plot_Gene - Plot all CDS of entire gene
#'
#' The function plots entire RFP reads of a specified gene.
#' @param  geneid geneid can be knownGene name (UCSC) or geneSymbole
#' @param wig wig file with reads
#' @param bed DF - the bed file
#' @param barsWidth plot's bars width
#' @return  ggplot object
#' @export

Plot_Gene <- function(geneid,wig,bed,barsWidth){
  chr <- as.character(bed[bed$hg19.kgXref.geneSymbol==geneid,which(colnames(bed)=="hg19.knownGene.chrom")])
  start <- bed[bed$hg19.kgXref.geneSymbol==geneid,which(colnames(bed)=="hg19.knownGene.cdsStart")]
  end <- bed[bed$hg19.kgXref.geneSymbol==geneid,which(colnames(bed)=="hg19.knownGene.cdsEnd")]
  subs <- wig[wig$seqnames==chr & wig$start>start & wig$start<end,]

  plotG <- ggplot(subs,mapping = aes(x = subs$start,y = abs(subs$score))) +
    geom_bar(stat="Identity",width = barsWidth) + ylab("Reads")+ xlab("Position")+ ggtitle(geneid) +
    theme(axis.text.x = element_blank())+
    theme_Publication() + scale_fill_Publication()
  return(plotG)
}

#geneid <- "WARS"
#Plot_Gene(geneid,wig,bed,50)


####################################################### Plot aa vacinity #######################################################
# RFPIFNaMD55A.m <- as.data.frame(import.wig("z:/edeng/NSG/runs/MD55/alignments/fp_chx_uninf_7RFPIFNaMD55A_R1/human/wig/unique_align/fp_chx_uninf_7RFPIFNaMD55A_R1_hg19.align.unique.p_offset.5_prime.merged.minus_Copy.wig"))
# RFPIFNaMD55A.p <- as.data.frame(import.wig("z:/edeng/NSG/runs/MD55/alignments/fp_chx_uninf_7RFPIFNaMD55A_R1/human/wig/unique_align/fp_chx_uninf_7RFPIFNaMD55A_R1_hg19.align.unique.p_offset.5_prime.merged.plus_Copy.wig"))
# RFPmIFNaMD55A.m <- as.data.frame(import.wig("z:/edeng/NSG/runs/MD55/alignments/fp_chx_uninf_5RFPmIFNaMD55A_S5_R1_001/human/wig/unique_align/fp_chx_uninf_5RFPmIFNaMD55A_S5_R1_001_hg19.align.unique.p_offset.5_prime.merged.minus_COPY.wig"))
# RFPmIFNaMD55A.p <- as.data.frame(import.wig("z:/edeng/NSG/runs/MD55/alignments/fp_chx_uninf_5RFPmIFNaMD55A_S5_R1_001/human/wig/unique_align/fp_chx_uninf_5RFPmIFNaMD55A_S5_R1_001_hg19.align.unique.p_offset.5_prime.merged.plus_COPY.wig"))



# aa= 'W'
# geneid='WARS'
# wig=RFPIFNaMD55A.m
# bed=bed
# CDSdb=CDS
# numberOfaa = "FULL"
# SINGLE = 1
# across = across
# plot = T
# genesymORucsc = "genesym"
# exonsORcds = "exons"
# data = F



#' Plot aa vacinity
#'
#' plots the vacinity of a specific aa specified.
#' @param aa amino acid. character of one capital letter.
#' @param geneid geneid can be knownGene name (UCSC) or geneSymbole
#' @param wig wig file with reads
#' @param bed DF - the bed file
#' @param CDSdb cds df
#' @param numberOfaa multiplot several aa. takes an integer with max value- number of specified aa apppearences.
#' @param SINGLE if only one aa (one plot) is needed. this iteger specify which aa will be plotted (e.g the third tryptophan)
#' @param across across how many nt down and up stream. integer.
#' @param exonsORcds takes "exons" or "cds".
#' @param genesymORucsc "genesym" or "ucsc" used geneid
#' @param plot boolean. whether to plot or not
#' @param data boolean return the data as DF?
#' @return  ggplot object and/or DF
#' @export
#'
Plot_vacinity <- function(aa,geneid,wig,bed,CDSdb,numberOfaa,SINGLE,across,genesymORucsc,exonsORcds,plot,data){
  if (sum(colnames(CDS)==paste(aa,"count",sep = ""))==0) {
    CDS <<- Add_AAtoCDS(aa,CDSdb)
    Plot_vacinity(aa,geneid,wig,bed,CDS,numberOfaa,SINGLE,across,genesymORucsc,exonsORcds,plot,data)

  }
  CDS <- CDS[CDS[paste(aa,"count",sep = "")]>0,]
  if (genesymORucsc=="genesym") {
    pos <- str_split(as.character(CDS[CDS$geneSymbol==geneid,paste(aa,"pos",sep = "")]),pattern = ",")[[1]] %>% as.integer()
    chr <- as.character(bed[bed$hg19.kgXref.geneSymbol==geneid,which(colnames(bed)=="hg19.knownGene.chrom")])

  } else {
    pos <- str_split(as.character(CDS[CDS$UCSCid==geneid,paste(aa,"pos",sep = "")]),pattern = ",")[[1]] %>% as.integer()
    chr <- as.character(bed[bed$hg19.knownGene.name==geneid,which(colnames(bed)=="hg19.knownGene.chrom")])
  }
  listofplots <- list()
  listofdfs <- list()
  k=1
  if (is.numeric(numberOfaa) & numberOfaa<length(pos)) {
    for (i in 1:numberOfaa) {
      firstnt <- as.integer(AA_finder(aa,i,geneid,bed,CDS,exonsORcds,genesymORucsc)[1])
      lastnt <- as.integer(AA_finder(aa,i,geneid,bed,CDS,exonsORcds,genesymORucsc)[2])
      subs <- wig[wig$seqnames==chr & wig$start>(firstnt-across) & wig$start<(lastnt+across),]
      listofdfs[[i]] <- subs
      ticks <- sort(c(seq((firstnt-across),(lastnt+across),by = 10),(firstnt:lastnt)))
      ticks <- ticks[!duplicated(ticks)]
      ticks <- ticks[ticks>=min(subs$start) & ticks<=max(subs$start)]
      labels <- ticks
      posnt <- which(ticks %in% firstnt:lastnt ,arr.ind = ticks)
      codon <- unlist(CODON_finder(aa,i,genesymORucsc))
      if (CDS[CDS$geneSymbol==geneid,"strand"]=="-") {
        labels[posnt] <- c(as.character(codon["last"]),as.character(codon["middle"]),as.character(codon["first"]))
      }else {labels[posnt] <- c(as.character(codon["first"]),as.character(codon["middle"]),as.character(codon["last"]))}
      colors <- ifelse(labels==as.character(codon["first"]) | labels==as.character(codon["middle"]) | labels==as.character(codon["last"]), "red", "black")
      if (!is.na(sum(subs$start)) & sum(subs$start)>0) {
        listofplots[[k]] <- ggplot(wig[wig$seqnames==chr & wig$start>(firstnt-across) & wig$start<(lastnt+across),]
                                   ,mapping = aes(x = start,y = score)) +
          geom_bar(stat="Identity") + ylab("Reads") + xlab("Position") + ggtitle(paste(geneid,i,sep = "-")) +
          scale_x_continuous(breaks = ticks,labels = labels) + theme(axis.text.x = element_text(colour = colors))+
          theme_Publication() + scale_fill_Publication() + theme(axis.text.x = element_text(angle = 25,colour = colors))
        k=k+1
      }

    }
  } else if (numberOfaa=="FULL") {
    for (i in 1:length(pos)) {
      firstnt <- as.integer(AA_finder(aa,i,geneid,bed,CDS,exonsORcds,genesymORucsc)[1])
      lastnt <- as.integer(AA_finder(aa,i,geneid,bed,CDS,exonsORcds,genesymORucsc)[2])
      subs <- wig[wig$seqnames==chr & wig$start>(firstnt-across) & wig$start<(lastnt+across),]
      listofdfs[[i]] <- subs
      ticks <- sort(c(seq((firstnt-across),(lastnt+across),by = 10),(firstnt:lastnt)))
      ticks <- ticks[!duplicated(ticks)]
      ticks <- ticks[ticks>=min(subs$start) & ticks<=max(subs$start)]
      labels <- ticks
      posnt <- which(ticks %in% firstnt:lastnt ,arr.ind = ticks)
      codon <- unlist(CODON_finder(aa,i,genesymORucsc))
      if (CDS[CDS$geneSymbol==geneid,"strand"]=="-") {
        labels[posnt] <- c(as.character(codon["last"]),as.character(codon["middle"]),as.character(codon["first"]))
      }else {labels[posnt] <- c(as.character(codon["first"]),as.character(codon["middle"]),as.character(codon["last"]))}
      colors <- ifelse(labels==as.character(codon["first"]) | labels==as.character(codon["middle"]) | labels==as.character(codon["last"]), "red", "black")
      if (!is.na(sum(subs$start)) & sum(subs$start)>0) {
        listofplots[[k]] <- ggplot(wig[wig$seqnames==chr & wig$start>(firstnt-across) & wig$start<(lastnt+across),]
                                   ,mapping = aes(x = start,y = score)) +
          geom_bar(stat="Identity") + ylab("Reads") + xlab("Position") + ggtitle(paste(geneid,i,sep = "-")) +
          scale_x_continuous(breaks = ticks,labels = labels) + theme(axis.text.x = element_text(colour = colors))+
          theme_Publication() + scale_fill_Publication() + theme(axis.text.x = element_text(angle = 25,colour = colors))
        k=k+1
      }

    }

  } else if (numberOfaa == "SINGLE") {
    firstnt <- as.integer(AA_finder(aa,SINGLE,geneid,bed,CDS,exonsORcds,genesymORucsc)[1])
    lastnt <- as.integer(AA_finder(aa,SINGLE,geneid,bed,CDS,exonsORcds,genesymORucsc)[2])
    subs <- wig[wig$seqnames==chr & wig$start>(firstnt-across) & wig$start<(lastnt+across),]
    listofdfs[[i]] <- subs
    ticks <- sort(c(seq((firstnt-across),(lastnt+across),by = 10),(firstnt:lastnt)))
    ticks <- ticks[!duplicated(ticks)]
    ticks <- ticks[ticks>=min(subs$start) & ticks<=max(subs$start)]
    labels <- ticks
    posnt <- which(ticks %in% firstnt:lastnt ,arr.ind = ticks)
    codon <- unlist(CODON_finder(aa,i,genesymORucsc))
    if (CDS[CDS$geneSymbol==geneid,"strand"]=="-") {
      labels[posnt] <- c(as.character(codon["last"]),as.character(codon["middle"]),as.character(codon["first"]))
    }else {labels[posnt] <- c(as.character(codon["first"]),as.character(codon["middle"]),as.character(codon["last"]))}
    colors <- ifelse(labels==as.character(codon["first"]) | labels==as.character(codon["middle"]) | labels==as.character(codon["last"]), "red", "black")
    if (!is.na(sum(subs$start)) & sum(subs$start)>0) {
      listofplots[[k]] <- ggplot(wig[wig$seqnames==chr & wig$start>(firstnt-across) & wig$start<(lastnt+across),]
                                 ,mapping = aes(x = start,y = score)) +
        geom_bar(stat="Identity") + ylab("Reads") + xlab("Position") + ggtitle(paste(geneid,SINGLE,sep = "-")) +
        scale_x_continuous(breaks = ticks,labels = labels) + theme(axis.text.x = element_text(colour = colors))+
        theme_Publication() + scale_fill_Publication() + theme(axis.text.x = element_text(angle = 25,colour = colors))
      k=k+1
    }

  } else{print("Number of aa exeeded or wrong input")}
  if (length(listofplots)>0) {
    if (data) {
      return(listofdfs)
    } else {if (plot) {
      gene_plotofplots = ggarrange(plotlist = listofplots,common.legend = T)
      return(gene_plotofplots)
    }
      else{return(listofplots)}}


  }else { print("No Plots")}

}



#' Plot_vacinity_IFNvsControl
#'
#' plots the vacinity of a specific aa specified in the context of two treatments. calculated the ratio of IFN to Control
#' @param IFN_data/Control_data a list of DFs.
#' @param i This iteger specify which aa will be plotted (e.g the third tryptophan)
#' @param aa amino acid. character of one capital letter.
#' @param geneid geneid can be knownGene name (UCSC) or geneSymbole
#' @param bed DF - the bed file
#' @param CDSdb cds df
#' @param across across how many nt down and up stream. integer.
#' @param exonsORcds takes "exons" or "cds".
#' @param genesymORucsc "genesym" or "ucsc" used geneid
#' @return  ggplot object
#' @export
#'
Plot_vacinity_IFNvsControl <- function(IFN_data,Control_data,aa,i,geneid,bed,CDS,across,exonsORcds,genesymORucsc){
  firstnt <- as.integer(AA_finder(aa,i,geneid,bed,CDS,exonsORcds,genesymORucsc)[1])
  lastnt <- as.integer(AA_finder(aa,i,geneid,bed,CDS,exonsORcds,genesymORucsc)[2])
  merged <- merge(as.data.frame(IFN_data[[i]]),as.data.frame(Control_data[[i]][,c("start","score")]),by = "start",all.x = T)
  merged[is.na(merged)] <- 1
  merged[,c("score.x","score.y")] <- abs(merged[,c("score.x","score.y")])
  merged['IFN/Control'] <- merged[,c("score.x")]/merged[,c("score.y")]
  ticks <- sort(c(seq((firstnt-across),(lastnt+across),by = 10),(firstnt:lastnt)))
  ticks <- ticks[!duplicated(ticks)]
  ticks <- ticks[ticks>=min(merged$start) & ticks<=max(merged$start)]
  labels <- ticks
  posnt <- which(ticks %in% firstnt:lastnt ,arr.ind = ticks)
  codon <- unlist(CODON_finder(aa,i,genesymORucsc))
  labels[posnt] <- c("",aa,"")
  aaindex <- which(labels==aa)
  k=1
  if (CDS[CDS$geneSymbol==geneid,"strand"]=="-") {
    for (j in 1:(aaindex-2)) {
      labels[j] <- 10*(aaindex-j-1)}
    for (j in (aaindex+2):length(labels)) {
      labels[j] <- -10*(k)
      k <- k+1

    }
  }else{
    for (j in 1:(aaindex-2)) {
      labels[j] <- -10*(aaindex-j-1)}
    for (j in (aaindex+2):length(labels)) {
      labels[j] <- 10*(k)
      k <- k+1

    }
    if (labels[1]=="0") {
      labels[1]=""

    }

  }
  colors <- ifelse(labels==aa,"#fdb462" , "black")
  plot <- ggplot(merged ,mapping = aes(x = start,y = as.integer(merged$`IFN/Control`))) +
    geom_bar(stat="Identity",fill="#386cb0") + ylab("IFN/Control") + xlab("nt") + ggtitle(paste(geneid,sep = "-")) +
    scale_x_continuous(breaks = ticks,labels = labels) + theme(axis.text.x = element_text(colour = colors))+
    theme_Publication(base_size = 50) + scale_fill_Publication() + theme(axis.text.x = element_text(colour = colors))
  return(plot)

}



# CDS <- CDS_DB
# geneid <- "EMP1"
# across <- 150
# Plot_vacinity(aa = "W",
#               geneid = geneid,
#               wig = RFPIFNaMD55A.p,
#               bed = bed,CDSdb = CDS,
#               numberOfaa = 1,
#               SINGLE = 1,
#               across = across,
#               plot = T,
#               genesymORucsc = "genesym",
#               exonsORcds = "exons"
#               ,data = F)
# Plot_vacinity(aa = "W",
#               geneid = geneid,
#               wig = RFPmIFNaMD55A.p,
#               bed = bed,CDSdb = CDS,
#               numberOfaa = 3,
#               SINGLE = 1,
#               across = across,
#               plot = T,
#               genesymORucsc = "genesym",
#               exonsORcds = "exons",
#               data = F)

# across <- 600
# Plot_vacinity(aa = "W",
#             geneid = 'WARS',
#             wig = RFPIFNaMD55A.m,
#             bed = bed,CDSdb = CDS,
#             numberOfaa = "FULL",
#             SINGLE = 1,
#             across = across,
#             plot = T,
#             genesymORucsc = "genesym",
#             exonsORcds = "exons",
#             data = F)
# #
# IFN <-   Plot_vacinity(aa = "W",
#                             geneid = geneid,
#                             wig = RFPIFNaMD55A.p,
#                             bed = bed,CDSdb = CDS,
#                             numberOfaa = "FULL",
#                             SINGLE = 1,
#                             across = across,
#                             plot = T,
#                             genesymORucsc = "genesym",
#                             exonsORcds = "exons",
#                             data = T)

#                 Plot_vacinity_IFNvsControl(IFN_data = list(RFPIFNaMD55A.m,RFPIFNaMD55A.p),
#                          Control_data = list(RFPmIFNaMD55A.m,RFPmIFNaMD55A.p),
#                          aa = 'W',i = 2,
#                          geneid = geneid,
#                          bed = bed,CDS = CDS,
#                          across = across,
#                          exonsORcds = "exons",
#                          genesymORucsc = "genesym")
# IFN_data = list(RFPIFNaMD55A.m,RFPIFNaMD55A.p)
# Control_data=list(RFPmIFNaMD55A.m,RFPmIFNaMD55A.p)
# aa = 'W'
# i = 2
# geneid = geneid
# bed = bed
# CDS = CDS
# across = across
# exonsORcds = "exons"
# genesymORucsc = "genesym"


