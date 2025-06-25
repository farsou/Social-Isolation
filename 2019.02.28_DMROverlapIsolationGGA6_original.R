###############
# Authors: Fabio Pertille;
# Description: MeDIPS Isolation OVerllap
# Date: 07th of july, 2018
# Files: .txt
# Code: GBSMEDIP002_overlap
library(GenomicFeatures)
library(GenomicFeatures)
library("TxDb.Ggallus.UCSC.galGal6.refGene")

library(MEDIPS)
library("MEDIPSData")
library(GenomicRanges)
#setwd("/home/fabio/postdoc/IsolationEarly/RresultsGGA6")
setwd("H:/Medicina Veterinária/2017.POSDOC-ESALQ-USP/Experimento_Isolation/RresultsAAGGA6")
#load(file="mr.edgeR_E_IsolationEarly.rda")
#load(file="mr.edgeR_S_IsolationEarly.rda")
load(file="mr.edgeROI_E_IsolationEarly.rda")
load(file="mr.edgeROI_S_IsolationEarly.rda")

EdgeR_SW_Early_S_ROI<-mr.edgeROI_S; rm (mr.edgeROI_S )
EdgeR_SW_Early_E_ROI<-mr.edgeROI_E; rm (mr.edgeROI_E)

#setwd("/home/fabio/postdoc/IsolationOld/RresultsGGA6")
load(file="mr.edgeROI_E_IsolationOld.rda")
load(file="mr.edgeROI_S_IsolationOld.rda")


#setwd("/home/fabio/postdoc/IsolationOld/RresultsGGA6")
EdgeR_SW_Old_S_ROI<-mr.edgeROI_S; rm (mr.edgeROI_S )
EdgeR_SW_Old_E_ROI<-mr.edgeROI_E; rm (mr.edgeROI_E)

#BR_Early_S$tag <- "BR_E_S"; BR_Early_S$strand <- "*"
#TesteBR_Early_S <- with(BR_Early_S, GRanges(chr, IRanges(start, stop), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, tag))

EdgeR_SW_Early_S_ROI$Test <- "Young_S";  EdgeR_SW_Early_S_ROI$strand <- "*"
EdgeR_SW_Early_E_ROI$Test <- "Young_E";  EdgeR_SW_Early_E_ROI$strand <- "*"
EdgeR_SW_Old_S_ROI$Test <- "Adult_S";  EdgeR_SW_Old_S_ROI$strand <- "*"
EdgeR_SW_Old_E_ROI$Test <- "Adult_E";  EdgeR_SW_Old_E_ROI$strand <- "*"

EdgeR_SW_Early_S_ROI$Location <- paste0(EdgeR_SW_Early_S_ROI$chr, ":", EdgeR_SW_Early_S_ROI$start, "-", EdgeR_SW_Early_S_ROI$stop)
EdgeR_SW_Early_E_ROI$Location <- paste0(EdgeR_SW_Early_E_ROI$chr, ":", EdgeR_SW_Early_E_ROI$start, "-", EdgeR_SW_Early_E_ROI$stop)
EdgeR_SW_Old_S_ROI$Location <- paste0(EdgeR_SW_Old_S_ROI$chr, ":", EdgeR_SW_Old_S_ROI$start, "-", EdgeR_SW_Old_S_ROI$stop)
EdgeR_SW_Old_E_ROI$Location <- paste0(EdgeR_SW_Old_E_ROI$chr, ":", EdgeR_SW_Old_E_ROI$start, "-", EdgeR_SW_Old_E_ROI$stop)

SW_Early_S_ROI <-EdgeR_SW_Early_S_ROI
SW_Early_E_ROI <-EdgeR_SW_Early_E_ROI
SW_Old_S_ROI <-EdgeR_SW_Old_S_ROI
SW_Old_E_ROI <-EdgeR_SW_Old_E_ROI

#GenomicRanges Objetcts
GR_EdgeR_SW_Early_S_ROI<-with(EdgeR_SW_Early_S_ROI, GRanges(chr, IRanges(start, stop), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, Test, Location))
GR_EdgeR_SW_Early_E_ROI<-with(EdgeR_SW_Early_E_ROI, GRanges(chr, IRanges(start, stop), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, Test, Location))
GR_EdgeR_SW_Old_S_ROI<-with(EdgeR_SW_Old_S_ROI, GRanges(chr, IRanges(start, stop), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, Test, Location))
GR_EdgeR_SW_Old_E_ROI<-with(EdgeR_SW_Old_E_ROI, GRanges(chr, IRanges(start, stop), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, Test, Location))

GR_SW_Early_S_ROI <- GR_EdgeR_SW_Early_S_ROI
GR_SW_Early_E_ROI <- GR_EdgeR_SW_Early_E_ROI
GR_SW_Old_S_ROI <- GR_EdgeR_SW_Old_S_ROI
GR_SW_Old_E_ROI <- GR_EdgeR_SW_Old_E_ROI


# setwd("/home/fabio/postdoc/IsolationOverlap/gga6")
output_folder <- "//argos.storage.uu.se/MyFolder$/farso212/others/my papers/Isolation stress/isolation data and R/inta age/try"
library(BSgenome.Ggallus.UCSC.galGal6)
library(seqinr)
library("rGADEM")
library("Biostrings")




#SEQ
Seq_GR_SW_Early_S_ROI = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_SW_Early_S_ROI)
Seq_GR_SW_Early_E_ROI = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_SW_Early_E_ROI)
Seq_GR_SW_Old_S_ROI = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_SW_Old_S_ROI)
Seq_GR_SW_Old_E_ROI = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_SW_Old_E_ROI)

Seq_CG_GR_SW_Early_S_ROI <- vcountPattern("CG", Seq_GR_SW_Early_S_ROI)
Seq_CG_GR_SW_Early_E_ROI <- vcountPattern("CG", Seq_GR_SW_Early_E_ROI)
Seq_CG_GR_SW_Old_S_ROI <- vcountPattern("CG", Seq_GR_SW_Old_S_ROI)
Seq_CG_GR_SW_Old_E_ROI <- vcountPattern("CG", Seq_GR_SW_Old_E_ROI)

Seq_GC_GR_SW_Early_S_ROI <- vcountPattern("GC", Seq_GR_SW_Early_S_ROI)
Seq_GC_GR_SW_Early_E_ROI <- vcountPattern("GC", Seq_GR_SW_Early_E_ROI)
Seq_GC_GR_SW_Old_S_ROI <- vcountPattern("GC", Seq_GR_SW_Old_S_ROI)
Seq_GC_GR_SW_Old_E_ROI <- vcountPattern("GC", Seq_GR_SW_Old_E_ROI)


#add CpG and GpC collumn to the Grange objects
CG<- Seq_CG_GR_SW_Early_S_ROI
GC<- Seq_GC_GR_SW_Early_S_ROI
SW_Early_S_ROI <- cbind(SW_Early_S_ROI, CG,GC)
CG<- Seq_CG_GR_SW_Early_E_ROI
GC<- Seq_GC_GR_SW_Early_E_ROI
SW_Early_E_ROI <- cbind(SW_Early_E_ROI, CG,GC)
CG<- Seq_CG_GR_SW_Old_S_ROI
GC<- Seq_GC_GR_SW_Old_S_ROI
SW_Old_S_ROI <- cbind(SW_Old_S_ROI, CG,GC)
CG<- Seq_CG_GR_SW_Old_E_ROI
GC<- Seq_GC_GR_SW_Old_E_ROI
SW_Old_E_ROI <- cbind(SW_Old_E_ROI, CG,GC)

#Genomic Ranges Object
GR_SW_Early_S_ROI <-with(SW_Early_S_ROI, GRanges(chr, IRanges(start, stop), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, Test, Location, CG, GC))
GR_SW_Early_E_ROI <-with(SW_Early_E_ROI, GRanges(chr, IRanges(start, stop), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, Test, Location, CG, GC))
GR_SW_Old_S_ROI <-with(SW_Old_S_ROI, GRanges(chr, IRanges(start, stop), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, Test, Location, CG, GC))
GR_SW_Old_E_ROI <-with(SW_Old_E_ROI, GRanges(chr, IRanges(start, stop), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, Test, Location, CG, GC))

SW_Young_S_adj0.5 <- subset(GR_SW_Early_S_ROI, edgeR.adj.p.value<.51)
SW_Young_E_adj0.5 <- subset(GR_SW_Early_E_ROI, edgeR.adj.p.value<.51)
SW_Adult_S_adj0.5 <- subset(GR_SW_Old_S_ROI, edgeR.adj.p.value<.51)
SW_Adult_E_adj0.5 <- subset(GR_SW_Old_E_ROI, edgeR.adj.p.value<.51)

SW_All_adj0.5 <- c(SW_Young_S_adj0.5, SW_Young_E_adj0.5, SW_Adult_S_adj0.5, SW_Adult_E_adj0.5)

#Sweden low p-value
#BRSW_EE_S_ROI <- subset(GR_EdgeR_SW_Early_S_ROI, edgeR.p.value<.0005)
#All low p-value 0.0005
GR_SW_Early_S_ROI_p.0005<- subset(GR_SW_Early_S_ROI, edgeR.p.value<.0005)
GR_SW_Early_E_ROI_p.0005<- subset(GR_SW_Early_E_ROI, edgeR.p.value<.0005)
GR_SW_Old_S_ROI_p.0005<- subset(GR_SW_Old_S_ROI, edgeR.p.value<.0005)
GR_SW_Old_E_ROI_p.0005<- subset(GR_SW_Old_E_ROI, edgeR.p.value<.0005)
#p-value0.005
GR_SW_Early_S_ROI_p.005<- subset(GR_SW_Early_S_ROI, edgeR.p.value<.005)
GR_SW_Early_E_ROI_p.005<- subset(GR_SW_Early_E_ROI, edgeR.p.value<.005)
GR_SW_Old_S_ROI_p.005<- subset(GR_SW_Old_S_ROI, edgeR.p.value<.005)
GR_SW_Old_E_ROI_p.005<- subset(GR_SW_Old_E_ROI, edgeR.p.value<.005)

GR_SW_Early_S_ROI_p.05<- subset(GR_SW_Early_S_ROI, edgeR.p.value<.05)
GR_SW_Early_E_ROI_p.05<- subset(GR_SW_Early_E_ROI, edgeR.p.value<.05)
GR_SW_Old_S_ROI_p.05<- subset(GR_SW_Old_S_ROI, edgeR.p.value<.05)
GR_SW_Old_E_ROI_p.05<- subset(GR_SW_Old_E_ROI, edgeR.p.value<.05)




#write ranges
write.table(GR_SW_Early_S_ROI_p.05 , "GR_SW_Early_S_ROI_p.05.txt", sep="\t", col.names=T, row.names=T, quote=F)
write.table(GR_SW_Early_E_ROI_p.05 , "GR_SW_Early_E_ROI_p.05.txt", sep="\t", col.names=T, row.names=T, quote=F)
write.table(GR_SW_Old_S_ROI_p.05 , "GR_SW_Old_S_ROI_p.05.txt", sep="\t", col.names=T, row.names=T, quote=F)
write.table(GR_SW_Old_E_ROI_p.05 , "GR_SW_Old_E_ROI_p.05.txt", sep="\t", col.names=T, row.names=T, quote=F)

write.table(GR_SW_Early_S_ROI_p.005 , "GR_SW_Early_S_ROI_p.005.txt", sep="\t", col.names=T, row.names=T, quote=F)
write.table(GR_SW_Early_E_ROI_p.005 , "GR_SW_Early_E_ROI_p.005.txt", sep="\t", col.names=T, row.names=T, quote=F)
write.table(GR_SW_Old_S_ROI_p.005 , "GR_SW_Old_S_ROI_p.005.txt", sep="\t", col.names=T, row.names=T, quote=F)
write.table(GR_SW_Old_E_ROI_p.005 , "GR_SW_Old_E_ROI_p.005.txt", sep="\t", col.names=T, row.names=T, quote=F)

write.table(GR_SW_Early_S_ROI_p.0005 , "GR_SW_Early_S_ROI_p.0005.txt", sep="\t", col.names=T, row.names=T, quote=F)
write.table(GR_SW_Early_E_ROI_p.0005 , "GR_SW_Early_E_ROI_p.0005.txt", sep="\t", col.names=T, row.names=T, quote=F)
write.table(GR_SW_Old_S_ROI_p.0005 , "GR_SW_Old_S_ROI_p.0005.txt", sep="\t", col.names=T, row.names=T, quote=F)
write.table(GR_SW_Old_E_ROI_p.0005 , "GR_SW_Old_E_ROI_p.0005.txt", sep="\t", col.names=T, row.names=T, quote=F)

#All low p-value 0.005

#Extend range in charge to get a larger interval
#windextend <- BRSW_EE_S
#start(windextend) <- start(windextend) - 300
#end(windextend) <- end(windextend) + 300


#SEQ
setwd("/home/fabio/postdoc/IsolationOverlap/gga6")
seq_windextend = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, windextend)
Biostrings::writeXStringSet(seq_windextend, "seq_windextend.fasta")


Seq_GR_SW_Early_S_ROI_p.0005 = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_SW_Early_S_ROI_p.0005)
Seq_GR_SW_Early_E_ROI_p.0005 = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_SW_Early_E_ROI_p.0005)
Seq_GR_SW_Old_S_ROI_p.0005 = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_SW_Old_S_ROI_p.0005)
Seq_GR_SW_Old_E_ROI_p.0005 = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_SW_Old_E_ROI_p.0005)


names(seq) = paste0("SEQUENCE_", seq_along(seq))
Biostrings::writeXStringSet(seq, "olapsBReSWe.fasta")

names(seq1) = paste0("SEQUENCE_", seq_along(seq1))
Biostrings::writeXStringSet(seq1, "lapsBReSWo.fasta")

names(seq) = paste0("SEQUENCE_", seq_along(seq))
Biostrings::writeXStringSet(seq, "my.fasta")

write.table(Seq_GR_SW_Early_S_ROI_p.0005 , "Seq_GR_SW_Early_S_ROI_p.0005.txt", sep="\t", col.names=T, row.names=T, quote=F)
write.table(Seq_GR_SW_Early_E_ROI_p.0005 , "Seq_GR_SW_Early_E_ROI_p.0005.txt", sep="\t", col.names=T, row.names=T, quote=F)
write.table(Seq_GR_SW_Old_S_ROI_p.0005 , "Seq_GR_SW_Old_S_ROI_p.0005.txt", sep="\t", col.names=T, row.names=T, quote=F)
write.table(Seq_GR_SW_Old_E_ROI_p.0005 , "Seq_GR_SW_Old_E_ROI_p.0005.txt", sep="\t", col.names=T, row.names=T, quote=F)


#to extend the windows
#start(GR_BR_Early_S)<- start(GR_BR_Early_S) - 1000
#stop(GR_BR_Early_S)<- stop(GR_BR_Early_S) + 1000
#start(gr) <- start(gr) - 10
#end(gr) <- end(gr) + 10

#start(GR_SW_Early_S_ROI_p.0005)<- start(GR_SW_Early_S_ROI_p.0005) - 300
#start(GR_SW_Early_E_ROI_p.0005)<- start(GR_SW_Early_E_ROI_p.0005) - 300
#start(GR_SW_Old_S_ROI_p.0005)<- start(GR_SW_Old_S_ROI_p.0005) - 300
#start(GR_SW_Old_E_ROI_p.0005)<- start(GR_SW_Old_E_ROI_p.0005) - 300
#end(GR_SW_Early_S_ROI_p.0005)<- end(GR_SW_Early_S_ROI_p.0005) + 300
#end(GR_SW_Early_E_ROI_p.0005)<- end(GR_SW_Early_E_ROI_p.0005) + 300
#end(GR_SW_Old_S_ROI_p.0005)<- end(GR_SW_Old_S_ROI_p.0005) + 300
#end(GR_SW_Old_E_ROI_p.0005)<- end(GR_SW_Old_E_ROI_p.0005) + 300
#write.table(Seq_GR_SW_Early_S_ROI_p.0005 , "Seq_GR_SW_Early_S_ROI_p.0005_300bp.txt", sep="\t", col.names=T, row.names=T, quote=F)
#write.table(Seq_GR_SW_Early_E_ROI_p.0005 , "Seq_GR_SW_Early_E_ROI_p_300bp.0005.txt", sep="\t", col.names=T, row.names=T, quote=F)
#write.table(Seq_GR_SW_Old_S_ROI_p.0005 , "Seq_GR_SW_Old_S_ROI_p.0005_300bp.txt", sep="\t", col.names=T, row.names=T, quote=F)
#write.table(Seq_GR_SW_Old_E_ROI_p.0005 , "Seq_GR_SW_Old_E_ROI_p.0005_300bp.txt", sep="\t", col.names=T, row.names=T, quote=F)

#extend p.005 windows in 3k
#start(GR_SW_Early_S_ROI_p.005)<- start(GR_SW_Early_S_ROI_p.005) - 20000
#start(GR_SW_Early_E_ROI_p.005)<- start(GR_SW_Early_E_ROI_p.005) - 20000
#start(GR_SW_Old_S_ROI_p.005)<- start(GR_SW_Old_S_ROI_p.005) - 20000
#start(GR_SW_Old_E_ROI_p.005)<- start(GR_SW_Old_E_ROI_p.005) - 20000
#end(GR_BR_Early_S_ROI_p.005)<- end(GR_BR_Early_S_ROI_p.005) + 20000
#end(GR_SW_Early_S_ROI_p.005)<- end(GR_SW_Early_S_ROI_p.005) + 20000
#end(GR_SW_Early_E_ROI_p.005)<- end(GR_SW_Early_E_ROI_p.005) + 20000
#end(GR_SW_Old_S_ROI_p.005)<- end(GR_SW_Old_S_ROI_p.005) + 20000
#end(GR_SW_Old_E_ROI_p.005)<- end(GR_SW_Old_E_ROI_p.005) + 20000


#laps - Brasil Early vs Sweden Early
#SW_E_ES_p.05 <-subsetByOverlaps(GR_SW_Early_S_ROI_p.05, GR_SW_Early_E_ROI_p.05)
#SW_O_ES_p.05 <-subsetByOverlaps(GR_SW_Old_S_ROI_p.05, GR_SW_Old_E_ROI_p.05)

#SW_EO_E_p.05 <-subsetByOverlaps(GR_SW_Early_E_ROI_p.05, GR_SW_Old_E_ROI_p.05)
#SW_EO_S_p.05 <-subsetByOverlaps(GR_SW_Early_S_ROI_p.05, GR_SW_Old_S_ROI_p.05)

SW_E_SE <-subsetByOverlaps(GR_SW_Early_S_ROI, GR_SW_Early_E_ROI)
SW_O_SE <-subsetByOverlaps(GR_SW_Old_S_ROI, GR_SW_Old_E_ROI)
SW_E_ES <-subsetByOverlaps(GR_SW_Early_E_ROI, GR_SW_Early_S_ROI)
SW_O_ES <-subsetByOverlaps(GR_SW_Old_E_ROI, GR_SW_Old_S_ROI)


SW_EO_E <-subsetByOverlaps(GR_SW_Early_E_ROI, GR_SW_Old_E_ROI)
SW_EO_S <-subsetByOverlaps(GR_SW_Early_S_ROI, GR_SW_Old_S_ROI)
SW_OE_E <-subsetByOverlaps(GR_SW_Old_E_ROI, GR_SW_Early_E_ROI)
SW_OE_S <-subsetByOverlaps(GR_SW_Old_S_ROI, GR_SW_Early_S_ROI)

SW_E_SE_p.05 <-subsetByOverlaps(GR_SW_Early_S_ROI_p.05, GR_SW_Early_E_ROI_p.05)
SW_O_SE_p.05 <-subsetByOverlaps(GR_SW_Old_S_ROI_p.05, GR_SW_Old_E_ROI_p.05)
SW_E_ES_p.05 <-subsetByOverlaps(GR_SW_Early_E_ROI_p.05, GR_SW_Early_S_ROI_p.05)
SW_O_ES_p.05 <-subsetByOverlaps(GR_SW_Old_E_ROI_p.05, GR_SW_Old_S_ROI_p.05)


SW_EO_E_p.05 <-subsetByOverlaps(GR_SW_Early_E_ROI_p.05, GR_SW_Old_E_ROI_p.05)
SW_EO_S_p.05 <-subsetByOverlaps(GR_SW_Early_S_ROI_p.05, GR_SW_Old_S_ROI_p.05)
SW_OE_E_p.05 <-subsetByOverlaps(GR_SW_Old_E_ROI_p.05, GR_SW_Early_E_ROI_p.05)
SW_OE_S_p.05 <-subsetByOverlaps(GR_SW_Old_S_ROI_p.05, GR_SW_Early_S_ROI_p.05)

overlap0.5 <- c(SW_E_SE_p.05, SW_E_ES_p.05,  SW_O_SE_p.05, SW_O_ES_p.05, SW_OE_E_p.05, SW_EO_E_p.05, SW_OE_S_p.05, SW_EO_S_p.05)

EOvsSE <- subsetByOverlaps (SW_E_ES, SW_O_ES) #113
OEvsSE <- subsetByOverlaps (SW_O_ES, SW_E_ES) #113

ESvsEO <- subsetByOverlaps (SW_EO_E, SW_EO_S) #121
SEvsEO <- subsetByOverlaps (SW_EO_S, SW_EO_E) #121

EOSE<- subsetByOverlaps(EOvsSE,SEvsEO) #111
SEEO<- subsetByOverlaps(SEvsEO,EOvsSE) #111

#write.table(as.data.frame(SW_DMRROI_EE_S_p.05) , "SW_DMRROI_EE_S_p.05.txt", sep="\t", col.names=T, row.names=T, quote=F)
#write.table(as.data.frame(SW_ROIDMR_EE_S_p.05) , "SW_ROIDMR_EE_S_p.05.txt", sep="\t", col.names=T, row.names=T, quote=F)


########################################
#  http://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html 
########################################
  library(calibrate)
  require(maptools)
setwd("/home/fabio/postdoc/IsolationOverlap/gga6")

  resSW_ES <- GR_SW_Early_S_ROI
  resSW_EE <- GR_SW_Early_E_ROI
  resSW_OS <- GR_SW_Old_S_ROI
  resSW_OE <- GR_SW_Old_E_ROI
  
  resbind <- c(resSW_ES, resSW_EE,resSW_OS, resSW_OE)
  
   # res <- EdgeR_SW_Early_S
#First Question: What are the significant DMR windows in each of the 4 tests? (ES, OS, EE, OS)
  tiff("Vulcano_1-4_E-ESvsO_S.tiff", units="cm", width=30, height=30, res=600, compression = "lzw")
     # Make a basic volcano plot
  with(resbind, plot(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="white", main="Volcano plot", xlim=c(-4,6)))
   # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(resSW_ES, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="yellow"))
  with(subset(resSW_OS, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="yellow3"))
#  with(subset(resSW_EE, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="cyan"))
#  with(subset(resSW_OE, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="cyan3"))
  with(subset(resSW_ES, edgeR.p.value<.0005 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
  with(subset(resSW_OS, edgeR.p.value<.0005 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
#  with(subset(resSW_EE, edgeR.p.value<.0005 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
#  with(subset(resSW_OE, edgeR.p.value<.0005 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
#  with(subset(resSW_ES, edgeR.adj.p.value<.5 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
#  with(subset(resSW_OS, edgeR.adj.p.value<.5 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
#  with(subset(resSW_EE, edgeR.adj.p.value<.5 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
#  with(subset(resSW_OE, edgeR.adj.p.value<.5 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
  
  
  abline(1.30103, 0, col = "red")
  abline(2.30103, 0, col = "yellow")
  abline(3.30103, 0, col = "green")
  dev.off() 
  
#Second Question: Is there DMR overlap between stress and enriched environments when compared to controls? 
  tiff("Vulcano_5_eoSEvsSE.tiff", units="cm", width=30, height=30, res=600, compression = "lzw")
  
  # Make a basic volcano plot
  #with(c(SW_E_ES, SW_O_ES), plot(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="gray", main="Volcano plot", xlim=c(-4,6)))
  with(resbind, plot(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="white", main="Volcano plot", xlim=c(-4,6)))
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(SW_E_ES, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="seagreen1"))
  with(subset(SW_E_SE, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="seagreen1"))
 
   with(subset(SW_O_ES, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="seagreen4"))
  with(subset(SW_O_SE, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="seagreen4"))
  # 
  #with(subset(SW_E_ES, edgeR.p.value<.05 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
#  with(subset(SW_O_ES, edgeR.p.value<.05 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
  abline(1.30103, 0, col = "red")
  abline(2.30103, 0, col = "yellow")
  abline(3.30103, 0, col = "green")
  dev.off()

#Third: Knowing what remains methylated after cell cycle renewal  
  tiff("Vulcano_6_eoSSvsEE.tiff", units="cm", width=30, height=30, res=600, compression = "lzw")
  
  # Make a basic volcano plot
  with(resbind, plot(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="white", main="Volcano plot", xlim=c(-4,6)))
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(SW_EO_E, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="red"))
  with(subset(SW_OE_E, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="red"))
  
   with(subset(SW_EO_S, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="red4"))
   with(subset(SW_OE_S, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="red4"))
  
   # with(subset(SW_EO_E, edgeR.p.value<.05 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
#  with(subset(SW_EO_S, edgeR.p.value<.05 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
  abline(1.30103, 0, col = "red")
  abline(2.30103, 0, col = "yellow")
  abline(3.30103, 0, col = "green")
  dev.off()
  
#Fourth: Place in the same chart what is common across the treatments  (green) and what is common across  the time of living(red) 
  tiff("Vulcano_7_eoSSvsEE.tiff", units="cm", width=30, height=30, res=600, compression = "lzw")
  with(resbind, plot(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="white", main="Volcano plot", xlim=c(-4,6)))
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(EOvsSE, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="green"))
  with(subset(OEvsSE, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="green"))
  
  with(subset(ESvsEO, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="red"))
  with(subset(SEvsEO, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="red"))
 
  #with(subset(EOvsSE, edgeR.p.value<.05 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
  #with(subset(SEvsEO, edgeR.p.value<.05 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
  abline(1.30103, 0, col = "red")
  abline(2.30103, 0, col = "yellow")
  abline(3.30103, 0, col = "green")
  dev.off() 
  
  tiff("Vulcano_8_all.tiff", units="cm", width=30, height=30, res=600, compression = "lzw")
  with(resbind, plot(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="white", main="Volcano plot", xlim=c(-4,6)))
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(EOSE, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="black"))
  with(subset(SEEO, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="black"))
  
  abline(1.30103, 0, col = "red")
  abline(2.30103, 0, col = "yellow")
  abline(3.30103, 0, col = "green")
  dev.off() 
  
  
  
  #write ranges
  write.table(as.data.frame(resSW_ES), "1_resSW_ES.txt", sep="\t", col.names=T, row.names=F, quote=F) #yellow
  write.table(as.data.frame(resSW_OS), "1_resSW_OS.txt", sep="\t", col.names=T, row.names=F, quote=F) #yellow3
  write.table(as.data.frame(resSW_EE), "1_resSW_EE.txt", sep="\t", col.names=T, row.names=F, quote=F) #cyan
  write.table(as.data.frame(resSW_OE), "1_resSW_OE.txt", sep="\t", col.names=T, row.names=F, quote=F) #cyan3
  
  write.table(as.data.frame(subset(SW_E_SE, edgeR.p.value<.05)), "2_SW_E_SE.txt", sep="\t", col.names=T, row.names=F, quote=F) #6
  write.table(as.data.frame(subset(SW_E_ES, edgeR.p.value<.05)), "2_SW_E_ES.txt", sep="\t", col.names=T, row.names=F, quote=F) #11
  write.table(as.data.frame(subset(SW_O_SE, edgeR.p.value<.05)), "3_SW_O_SE.txt", sep="\t", col.names=T, row.names=F, quote=F) #13
  write.table(as.data.frame(subset(SW_O_ES, edgeR.p.value<.05)), "3_SW_O_ES.txt", sep="\t", col.names=T, row.names=F, quote=F) #4
  
  write.table(as.data.frame(subset(SW_EO_E, edgeR.p.value<.05)), "4_SW_EO_E.txt", sep="\t", col.names=T, row.names=F, quote=F) #4
  write.table(as.data.frame(subset(SW_OE_E, edgeR.p.value<.05)), "4_SW_OE_E.txt", sep="\t", col.names=T, row.names=F, quote=F) #5
  write.table(as.data.frame(subset(SW_EO_S, edgeR.p.value<.05)), "5_SW_EO_S.txt", sep="\t", col.names=T, row.names=F, quote=F) #5
  write.table(as.data.frame(subset(SW_OE_S, edgeR.p.value<.05)), "5_SW_OE_S.txt", sep="\t", col.names=T, row.names=F, quote=F) #20
  
  write.table(as.data.frame(subset(EOvsSE, edgeR.p.value<.05)), "6_EOvsSE.txt", sep="\t", col.names=T, row.names=F, quote=F) #2
  write.table(as.data.frame(subset(OEvsSE, edgeR.p.value<.05)), "6_OEvsSE.txt", sep="\t", col.names=T, row.names=F, quote=F) #0
  
  write.table(as.data.frame(subset(ESvsEO, edgeR.p.value<.05)), "7_ESvsEO.txt", sep="\t", col.names=T, row.names=F, quote=F) #3
  write.table(as.data.frame(subset(SEvsEO, edgeR.p.value<.05)), "7_SEvsEO.txt", sep="\t", col.names=T, row.names=F, quote=F) #3
  write.table(as.data.frame(subset(EOSE, edgeR.p.value<.05)), "8_EOSE.txt", sep="\t", col.names=T, row.names=F, quote=F) #15
  write.table(as.data.frame(subset(SEEO, edgeR.p.value<.05)), "8_SEEO.txt", sep="\t", col.names=T, row.names=F, quote=F) #18
  
#
  write.table(as.data.frame(SW_E_SE), "2_SW_E_SE.txt", sep="\t", col.names=T, row.names=F, quote=F) #6
  write.table(as.data.frame(SW_E_ES), "2_SW_E_ES.txt", sep="\t", col.names=T, row.names=F, quote=F) #11
  write.table(as.data.frame(SW_O_SE), "3_SW_O_SE.txt", sep="\t", col.names=T, row.names=F, quote=F) #13
  write.table(as.data.frame(SW_O_ES), "3_SW_O_ES.txt", sep="\t", col.names=T, row.names=F, quote=F) #4
  
  write.table(as.data.frame(SW_EO_E), "4_SW_EO_E.txt", sep="\t", col.names=T, row.names=F, quote=F) #4
  write.table(as.data.frame(SW_OE_E), "4_SW_OE_E.txt", sep="\t", col.names=T, row.names=F, quote=F) #5
  write.table(as.data.frame(SW_EO_S), "5_SW_EO_S.txt", sep="\t", col.names=T, row.names=F, quote=F) #5
  write.table(as.data.frame(SW_OE_S), "5_SW_OE_S.txt", sep="\t", col.names=T, row.names=F, quote=F) #20
  
  write.table(as.data.frame(EOvsSE), "6_EOvsSE.txt", sep="\t", col.names=T, row.names=F, quote=F) #2
  write.table(as.data.frame(OEvsSE), "6_OEvsSE.txt", sep="\t", col.names=T, row.names=F, quote=F) #0
  
  write.table(as.data.frame(ESvsEO), "7_ESvsEO.txt", sep="\t", col.names=T, row.names=F, quote=F) #3
  write.table(as.data.frame(SEvsEO), "7_SEvsEO.txt", sep="\t", col.names=T, row.names=F, quote=F) #3
  write.table(as.data.frame(EOSE), "8_EOSE.txt", sep="\t", col.names=T, row.names=F, quote=F) #15
  write.table(as.data.frame(SEEO), "8_SEEO.txt", sep="\t", col.names=T, row.names=F, quote=F) #18
  

SW_E_ES_2 <- c(SW_E_ES, SW_E_SE) ##seagreen1
SW_E_ES_2 <- sortSeqlevels(SW_E_ES_2)
SW_O_ES_2 <- c(SW_O_ES, SW_O_SE) #seagreen4
SW_O_ES_2 <- sortSeqlevels(W_O_ES_2)
EOvsSE_2 <- c(EOvsSE, OEvsSE) #green
EOvsSE_2 <- sortSeqlevels(EOvsSE_2)
ESvsEO_2 <- c(ESvsEO, SEvsEO) #red
ESvsEO_2 <- sortSeqlevels(ESvsEO_2)
EOSE_2 <- c(EOSE,SEEO)
EOSE_2 <- sortSeqlevels(EOSE_2)

write.table(as.data.frame(SW_E_ES_2), "SW_E_ES_2", sep="\t", col.names=T, row.names=F, quote=F)
write.table(as.data.frame(SW_O_ES_2), "SW_O_ES_2", sep="\t", col.names=T, row.names=F, quote=F)
write.table(as.data.frame(EOvsSE_2), "EOvsSE_2", sep="\t", col.names=T, row.names=F, quote=F)
write.table(as.data.frame(ESvsEO_2), "ESvsEO_2", sep="\t", col.names=T, row.names=F, quote=F)
write.table(as.data.frame(EOSE_2), "EOSE_2", sep="\t", col.names=T, row.names=F, quote=F)





#graph
library(qtl)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)
library(gghighlight)
library(ggbio)
library(tidyverse)

##  coord:genome
resbindS <- c(resSW_ES, resSW_OS)
resbindE <- c(resSW_EE,resSW_OE)

resbindA <- c(resSW_OS, resSW_OE)
resbindY <- c(resSW_ES, resSW_EE)



resbind <- sortSeqlevels(resbind)
resbindS <- sortSeqlevels(resbindS)
resbindE <- sortSeqlevels(resbindE)
resbindA <- sortSeqlevels(resbindA)
resbindY <- sortSeqlevels(resbindY)

SW_E_ES_2<- sortSeqlevels(SW_E_ES_2)
#install.packages("Cairo")
#setHook(packageEvent("grDevices", "onLoad"),
#        function(...) grDevices::X11.options(type='cairo'))
#options(device='x11')


tiff("Chicken_Isolation_manhatahn.tiff", units="cm", width=30.00, height=15.00, res=600, compression = "lzw")

autoplot(resbind, coord = "genome", geom = "point", aes(y = -log10(edgeR.p.value), color=Test), space.skip = 0.01, spaceline = TRUE, legend = TRUE) +

  scale_color_manual(values=c("cyan3", "yellow3", "cyan", "yellow"))+
#                                OE       OS        EE      ES        
#  scale_color_manual(values=c("yellow3", "yellow"))+
# scale_color_manual(values=c("cyan3", "cyan"))+
  
   geom_hline(yintercept= 3.30103, color='seagreen', size=0.2) +
  geom_hline(yintercept= 2.30103, color='yellow4', size=0.2) +
  geom_hline(yintercept= 1.30103, color='red', size=0.2) +

  theme(axis.text.x=element_text(angle=90, hjust=0, size = 6)) 

dev.off()

#resbindS$sig <- "sig"
#resbindS$sig[resbindS$edgeR.p.value>.05] <- "nonsig"

tiff("Chicken_Foldchange_manhatahn_E.tiff", units="cm", width=30.00, height=15.00, res=600, compression = "lzw")

autoplot(resbindE, coord = "genome", geom = "point", aes(y = edgeR.logFC, color=Test), space.skip = 0.01, spaceline = TRUE, legend = TRUE) +
 #  scale_color_manual(values=c("cyan3", "yellow3", "cyan", "yellow"))+
geom_point(color = ifelse(resbindE$edgeR.p.value<.05, "black", NA), shape=23, size = 3) +
#geom_text(aes(label=ifelse(edgeR.adj.p.value<.5, as.character(Location),'')),hjust=0,vjust=0)

  #geom_point(aes(if(resbindS$edgeR.p.value<.05)), color=TRUE, shape=23)+
  
  #                                OE       OS        EE      ES        
  #scale_color_manual(values=c("yellow3", "yellow"))+
                               
  scale_color_manual(values=c("cyan3", "cyan"))+
  
# geom_hline(yintercept= 3.30103, color='seagreen', size=0.2) +
#  geom_hline(yintercept= 2.30103, color='yellow4', size=0.2) +
#  geom_hline(yintercept= 1.30103, color='red', size=0.2) +
  geom_hline(yintercept= 0, color='black', size=0.2) +


 # geom_point(shape=23, aes(color= sig , size=3))+
  
  theme(axis.text.x=element_text(angle=90, hjust=0, size = 6))
  #geom_point(data = ex, mapping = aes(x = Location, y =  -log10(edgeR.p.value)), color="red") +
  

  


dev.off()

tiff("Chicken_Foldchange_manhatahn_S.tiff", units="cm", width=30.00, height=15.00, res=600, compression = "lzw")

autoplot(resbindS, coord = "genome", geom = "point", aes(y = edgeR.logFC, color=Test), space.skip = 0.01, spaceline = TRUE, legend = TRUE) +
  geom_point(color = ifelse(resbindS$edgeR.p.value<.05, "black", NA), shape=23, size = 3) +
  
    #  scale_color_manual(values=c("cyan3", "yellow3", "cyan", "yellow"))+
  #                                OE       OS        EE      ES        
#  scale_color_manual(values=c("yellow3", "yellow"))+     #stress
 scale_color_manual(values=c("cyan3", "cyan"))+           #enrich
  
  # geom_hline(yintercept= 3.30103, color='seagreen', size=0.2) +    
  #  geom_hline(yintercept= 2.30103, color='yellow4', size=0.2) +    
  #  geom_hline(yintercept= 1.30103, color='red', size=0.2) +
  geom_hline(yintercept= 0, color='black', size=0.2) +
  
  theme(axis.text.x=element_text(angle=90, hjust=0, size = 6)) 

dev.off()

tiff("Chicken_Foldchange_manhatahn_A.tiff", units="cm", width=30.00, height=15.00, res=600, compression = "lzw")
autoplot(resbindA, coord = "genome", geom = "point", aes(y = edgeR.logFC, color=Test), space.skip = 0.01, spaceline = TRUE, legend = TRUE) +
  geom_point(color = ifelse(resbindA$edgeR.p.value<.05, "black", NA), shape=23, size = 3) +
  
    #  scale_color_manual(values=c("cyan3", "yellow3", "cyan", "yellow"))+
  #                                OE       OS        EE      ES        
  #  scale_color_manual(values=c("yellow3", "yellow"))+
  #scale_color_manual(values=c("cyan3", "yellow3"))+    #adult
  scale_color_manual(values=c("cyan", "yellow"))+   #young
  
  # geom_hline(yintercept= 3.30103, color='seagreen', size=0.2) +
  #  geom_hline(yintercept= 2.30103, color='yellow4', size=0.2) +
  #  geom_hline(yintercept= 1.30103, color='red', size=0.2) +
  geom_hline(yintercept= 0, color='black', size=0.2) +
  
  theme(axis.text.x=element_text(angle=90, hjust=0, size = 6)) 

dev.off()


#SW_E_ES_2_NA <- na.omit(SW_E_ES_2, cols=c("edgeR.logFC"))
autoplot(SW_E_ES_2, coord = "genome", geom = "point", aes(y = edgeR.logFC, color=Test), space.skip = 0.01, spaceline = TRUE, legend = TRUE) +
  #  scale_color_manual(values=c("cyan3", "yellow3", "cyan", "yellow"))+
  #                                OE       OS        EE      ES        
  #  scale_color_manual(values=c("yellow3", "yellow"))+
  scale_color_manual(values=c("cyan3", "cyan"))+
  
  # geom_hline(yintercept= 3.30103, color='seagreen', size=0.2) +
  #  geom_hline(yintercept= 2.30103, color='yellow4', size=0.2) +
  #  geom_hline(yintercept= 1.30103, color='red', size=0.2) +
  geom_hline(yintercept= 0, color='black', size=0.2) +
  
  theme(axis.text.x=element_text(angle=90, hjust=0, size = 6)) 

dev.off()


################################
###############################
##############################
#some counts
width(resbind)
genome <- BSgenome.Ggallus.UCSC.galGal6
head(seqlengths(genome))
getSeq(genome, as.character=FALSE)

aggregate(width(resbind), by=list(Category=resbind$Test), function(x) c(mean = mean(x), sd = sd(x)))

#Category	x.mean	x.sd
#Adult_E	479.2977	529.2058
#Adult_S	351.3135	381.4918
#Young_E	535.5546	521.3032
#Young_S	513.3682	505.3228
plot(density(as.matrix(resbind$edgeR.p.value)))
plot(density(as.matrix(resbind$edgeR.logFC)))

#size simulation
library('regioneR')
AE <- createRandomRegions(nregions=100000, length.mean=479, length.sd=529, genome=genome)
AS<- createRandomRegions(nregions=100000, length.mean=351, length.sd=381, genome=genome)
YE <- createRandomRegions(nregions=100000, length.mean=535, length.sd=521, genome=genome)
YS<- createRandomRegions(nregions=100000, length.mean=513, length.sd=505, genome=genome)
TT<- createRandomRegions(nregions=100000, length.mean=470, length.sd=484, genome=genome)

TT_df<- as.data.frame (TT)
table(TT_df$seqnames)

resbind_df <- as.data.frame(resbind)
table(resbind_df$seqnames)


#annotating the significant ones
#overlap0.5 - todos os overlaps
#SW_All_adj0.5 - todos os sig
setwd("H:/Medicina Veterinária/2017.POSDOC-ESALQ-USP/Experimento_Isolation/RresultsAAGGA6/")
library (BSgenome.Ggallus.UCSC.galGal6)
library ("ChIPseeker")
library(MEDIPS)
library("MEDIPSData")
library (BSgenome.Ggallus.UCSC.galGal6)
library(readr)
library(RMariaDB)
library(GenomicFeatures)
library(org.Gg.eg.db)
library(ChIPseeker)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rtracklayer)
library("GenomicRanges")
library('seqinr')
library("rGADEM")
library("Biostrings")
library(org.Gg.eg.db)
library(ChIPpeakAnno)
library("Cairo")

supportedUCSCtables(genome = "galGal6")
gg_txdb <- makeTxDbFromUCSC(genome="galGal6", tablename="ensGene")

#DO ONE BY ONE
peak_file<- as.data.frame(overlap0.5)
#peak_file<- as.data.frame (SW_All_adj0.5)


peak_file <- subset(peak_file, select=c(seqnames, start,end))
names(peak_file) <- c("CHR", "BP", "BP2")

peak_anno_function <- function(input_data){
  b2 <- input_data[,c("CHR" , "BP" , "BP2")]
  write.table(b2 , "annotate_file.txt" , sep = "\t" , col.names = TRUE , row.names = FALSE)
  peakAnno <- annotatePeak("H:/Medicina Veterinária/2017.POSDOC-ESALQ-USP/Experimento_Isolation/RresultsAAGGA6/annotate_file.txt",
                           TxDb=gg_txdb, annoDb="org.Gg.eg.db")
  return(peakAnno)
}

peak_chicken <- peak_anno_function(peak_file)
plotAnnoPie(peak_chicken)
info_chicken <- peak_chicken@anno@elementMetadata

overlap0.5  <- cbind(overlap0.5, info_chicken , all=T)
#SW_All_adj0.5 <- cbind(SW_All_adj0.5, info_chicken , all=T)

write.table(overlap0.5  , "overlap0.5 .txt" , sep="\t" , row.names = F , quote = F)
#write.table(SW_All_adj0.5 , "SW_All_adj0.5.txt" , sep="\t" , row.names = F , quote = F)

#Venn Diagramm
resSW_ES05 <- GR_SW_Early_S_ROI_p.05
resSW_EE05 <- GR_SW_Early_E_ROI_p.05
resSW_OS05 <- GR_SW_Old_S_ROI_p.05
resSW_OE05 <- GR_SW_Old_E_ROI_p.05
resbind05 <- c(resSW_ES05, resSW_EE05,resSW_OS05, resSW_OE05)

grl <- splitAsList(resbind05, resbind05$Test)
CairoPDF("vennShort_term05.pdf")
res <- makeVennDiagram(Peaks=grl, NameOfPeaks=names(grl), main="Isolation: Shor-term: Young vs Adults (DMR p<0.05)", col = "black",fill = c("red", "blue", "yellow", "green"), alpha = 0.50, cex= 1.5)
dev.off()

grl <- splitAsList(resbind, resbind$Test)
CairoPDF("venn_Shor_term.pdf")
res <- makeVennDiagram(Peaks=grl, NameOfPeaks=names(grl), main="Isolation: Shor-term: Young vs Adults (peaking called ADJp<0.1)", col = "black",fill = c("red", "blue", "yellow", "green"), alpha = 0.50, cex= 1.5)
dev.off()



###SEQ only significant ones to use in the motif analysis
#setwd("/home/fabio/postdoc/IsolationOverlap/gga6/OldvsYoungControl")
GR_SW_Early_S_ROI <- GR_EdgeR_SW_Early_S_ROI
GR_SW_Early_E_ROI <- GR_EdgeR_SW_Early_E_ROI
GR_SW_Old_S_ROI <- GR_EdgeR_SW_Old_S_ROI
GR_SW_Old_E_ROI <- GR_EdgeR_SW_Old_E_ROI


GR_Y_S_adj <- subset(GR_SW_Early_S_ROI, edgeR.adj.p.value<.5)
GR_Y_E_adj <-  subset(GR_SW_Early_E_ROI, edgeR.adj.p.value<.5)
GR_A_S_adj <-  subset(GR_SW_Old_S_ROI, edgeR.adj.p.value<.5)
GR_A_E_adj <-  subset(GR_SW_Old_E_ROI, edgeR.adj.p.value<.5)

GR_Y_S <- subset(GR_SW_Early_S_ROI, edgeR.p.value<.05)
GR_Y_E <- subset(GR_SW_Early_E_ROI, edgeR.p.value<.05)
GR_A_S <- subset(GR_SW_Old_S_ROI, edgeR.p.value<.05)
GR_A_E <- subset(GR_SW_Old_E_ROI, edgeR.p.value<.05)

library(data.table)
library ("GenomicRanges")
library(BSgenome.Ggallus.UCSC.galGal6)
library(seqinr)
library("rGADEM")
Seq_GR_Y_S_adj = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_Y_S_adj)
Seq_GR_Y_E_adj = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_Y_E_adj)
Seq_GR_A_S_adj = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_A_S_adj)
Seq_GR_A_E_adj = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_A_E_adj)
Seq_GR_Y_S = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_Y_S)
Seq_GR_Y_E = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_Y_E)
Seq_GR_A_S = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_A_S)
Seq_GR_A_E = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_A_E)

names(Seq_GR_Y_S_adj) = GR_Y_S_adj$Location
names(Seq_GR_Y_E_adj) = GR_Y_E_adj$Location
names(Seq_GR_A_S_adj) = GR_A_S_adj$Location
names(Seq_GR_A_E_adj) = GR_A_E_adj$Location

names(Seq_GR_Y_S) = GR_Y_S$Location
names(Seq_GR_Y_E) = GR_Y_E$Location
names(Seq_GR_A_S) = GR_A_S$Location
names(Seq_GR_A_E) = GR_A_E$Location

setwd("/home/fabio/postdoc/IsolationEarly/RresultsGGA6/fasta")
setwd("H:/Medicina Veterinária/2017.POSDOC-ESALQ-USP/Experimento_Isolation/RresultsAAGGA6/fastaIsolationBrasilp005_adj05")
Biostrings::writeXStringSet(Seq_GR_Y_S_adj, "Y_S_adj.fasta")
Biostrings::writeXStringSet(Seq_GR_Y_E_adj, "Y_E_adj.fasta")
Biostrings::writeXStringSet(Seq_GR_A_S_adj, "A_S_adj.fasta")
Biostrings::writeXStringSet(Seq_GR_A_E_adj, "A_E_adj.fasta")

Biostrings::writeXStringSet(Seq_GR_Y_S, "Y_S.fasta")
Biostrings::writeXStringSet(Seq_GR_Y_E, "Y_E.fasta")
Biostrings::writeXStringSet(Seq_GR_A_S, "A_S.fasta")
Biostrings::writeXStringSet(Seq_GR_A_E, "A_E.fasta")

#MOTIFFS
#obtain genomic sequence of flanking regions

#NEW OVERLAP 
overlapALL <- findOverlapsOfPeaks(resSW_ES, resSW_EE, resSW_OS, resSW_OE)
overlapStres <- findOverlapsOfPeaks(resSW_ES, resSW_OS)
overlapEnrich<- findOverlapsOfPeaks(resSW_EE, resSW_OE)

overlapALL <- overlapALL$peaklist[["resSW_ES///resSW_EE///resSW_OS///resSW_OE"]]
overlapStres <- overlapStres$peaklist[["resSW_ES///resSW_OS"]]
overlapEnrich <- overlapEnrich$peaklist[["resSW_EE///resSW_OE"]]

overlapALL$peakNames <- "Overlap"
overlapStres$peakNames <- "OverlapStress"
overlapEnrich$peakNames <- "OverlapEnrich"

#p<0.05

lapGR_Y_S<- subsetByOverlaps(overlapStres, GR_Y_S)
lapGR_Y_E<- subsetByOverlaps(overlapEnrich, GR_Y_E)
lapGR_A_S<- subsetByOverlaps(overlapStres, GR_A_S)
lapGR_A_E<- subsetByOverlaps(overlapEnrich, GR_A_E)
OverlapMerge <- unique(c(lapGR_Y_S, lapGR_Y_E, lapGR_A_S, lapGR_A_E ))
OlapRegions <- OverlapMerge
OlapR <- OverlapMerge
#mcols(OlapR)$peakNames <- NULL
lapGR_Y_S<- subsetByOverlaps(GR_Y_S,overlapStres)
lapGR_Y_E<- subsetByOverlaps(GR_Y_E, overlapEnrich)
lapGR_A_S<- subsetByOverlaps(GR_A_S, overlapStres)
lapGR_A_E<- subsetByOverlaps(GR_A_E, overlapEnrich)

OverlapMerge <- unique(c(lapGR_Y_S, lapGR_Y_E, lapGR_A_S, lapGR_A_E ))
OlapDMR<- OverlapMerge

setwd("H:/Medicina Veterinária/2017.POSDOC-ESALQ-USP/Experimento_Isolation/RresultsAAGGA6/MOTIF")
grl <- splitAsList(OverlapMerge, OverlapMerge$Test)
tiff("Venn_OverlapedAllDMRs.tiff", units="cm", width=30.00, height=30.00, res=600, compression = "lzw")
#CairoPDF("venn.pdf")
res <- makeVennDiagram(Peaks=grl, NameOfPeaks=names(grl), main = "Number of DMRs found within overlapped regions ", col = "black", alpha = 0.50, cex= 2.5, ignore.strand=T, scaled=F, fill = c("red", "blue", "yellow", "green"))
#fill = c("Yellow", "cyan"), 
dev.off()

#all the peaks from all the comparisons (i will run these line to each one of the files using common and unique peaks for each treatment)
peak_file <- as.data.frame(OlapRegions)
peak_file <- as.data.frame(OverlapMerge)

peak_file <- subset(peak_file, select=c(seqnames, start,end))
names(peak_file) <- c("CHR", "BP", "BP2")

peak_anno_function <- function(input_data){
  b2 <- input_data[,c("CHR" , "BP" , "BP2")]
  write.table(b2 , "annotate_file.txt" , sep = "\t" , col.names = TRUE , row.names = FALSE)
  peakAnno <- annotatePeak("H:/Medicina Veterinária/2017.POSDOC-ESALQ-USP/Experimento_Isolation/RresultsAAGGA6/MOTIF/annotate_file.txt",
                           TxDb=gg_txdb, annoDb="org.Gg.eg.db")
  return(peakAnno)
}

peak_chicken <- peak_anno_function(peak_file)
plotAnnoPie(peak_chicken)
info_chicken <- peak_chicken@anno@elementMetadata

OlapRegions <- cbind(OlapRegions, info_chicken , all=T)
OverlapMerge<- cbind(OverlapMerge, info_chicken , all=T)

write.table(OlapRegions , "OlapRegions.txt" , sep="\t" , row.names = F , quote = F)
write.table(OverlapMerge , "OverlapedDMRs.txt" , sep="\t" , row.names = F , quote = F)

OlapR <- as.data.frame(OlapR)
OlapR$Location <- paste0(OlapR$seqnames, ":", OlapR$start, "-", OlapR$end)
OlapR<-with(OlapR, GRanges(seqnames, IRanges(start, end), strand, Location))

library(data.table)
library ("GenomicRanges")
library(BSgenome.Ggallus.UCSC.galGal6)
library(seqinr)
library("rGADEM")
Seq_OlapR = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, OlapR)
Seq_OlapDMR = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, OlapDMR)

names(Seq_OlapR) = OlapR$Location
names(Seq_OlapDMR) = OlapDMR$Location

#setwd("/home/fabio/postdoc/IsolationEarly/RresultsGGA6/fasta")
#setwd("H:/Medicina Veterinária/2017.POSDOC-ESALQ-USP/Experimento_Isolation/RresultsAAGGA6/fastaIsolationBrasilp005_adj05")
Biostrings::writeXStringSet(Seq_OlapR, "Seq_OlapR.fasta")
Biostrings::writeXStringSet(Seq_OlapDMR, "Seq_OlapDMR.fasta")


