# 
# library("MSnbase")
# library("Peptides")
# library("RGenetics")
# aa <- get.amino.acids()
# 
# aa$Abbrev3 <-  c("Asn","His")
# 
# 
# codonTable <- geneticCodeTable(DNA = FALSE)
# 
# codonTable <- aggregate(codonTable$GeneticCode  ~ codonTable$AminoAcids,codonTable, paste,collapse=",")
# 
# 
# aa <- merge(aa,codonTable,by.x = 'Abbrev3',by.y = 'codonTable$AminoAcids')
# codonTable <- aa
# use_data(codonTable)
# colnames(aa)[ncol(aa)] <- "RNA Codons"
# 
# 
# aa[aa$Abbrev3 %in% c('Asn','His','Lys','Gln'),]
# aa[aa$Abbrev3 %in% c('Leu','Ser'),]
# aa[aa$Abbrev3 %in% c('Asn','Ser'),]
# aa[aa$Abbrev3 %in% c('Asp','Glu'),]
# 
# 
# 
# hydrophobicity('SLREAESLLAK','Roseman')
