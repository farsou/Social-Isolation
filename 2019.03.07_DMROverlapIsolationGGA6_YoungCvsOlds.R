###############
# Authors: Fabio Pertille;
# Description: MeDIPS Isolation OVerllap YEarly C vs Olds
# Date: 07th of march, 2019
# Files: .txt
# Code: GBSMEDIP002_overlap

module load HPC/R-3.5.1
R
library(MEDIPS)
library("MEDIPSData")
library (BSgenome.Ggallus.UCSC.galGal6)
library(readr)
library(RMariaDB)
library(GenomicFeatures)
library(org.Gg.eg.db)
library(ChIPseeker)
library(ChIPpeakAnno)
library(TxDb.Ggallus.UCSC.galGal6.refGene)
library(rtracklayer)
library("GenomicRanges")
library('seqinr')
library("rGADEM")
library("Biostrings")
library(org.Gg.eg.db)

#setwd("/home/fabio/postdoc/IsolationOverlap/gga6/OldvsYoungControl")
setwd("C:/Users/fabpe757/ROCKSCIENCE/2017.POSDOC-ESALQ-USP/Experimento_Isolation/RresultsAAGGA6/EarlyCvsOlds")
load(file="mr.edgeROI_S_IsolationEarlyOLD.rda")
load(file="mr.edgeROI_E_IsolationEarlyOLD.rda")
load(file="mr.edgeROI_C_IsolationEarlyOLD.rda")

EdgeR_EarlyC_OLD_S <-mr.edgeROI_S; rm (mr.edgeROI_S)
EdgeR_EarlyC_OLD_E <-mr.edgeROI_E; rm (mr.edgeROI_E)
EdgeR_EarlyC_OLD_C <-mr.edgeROI_C; rm (mr.edgeROI_C)


EdgeR_EarlyC_OLD_S$Test <- "y_C_a_S";  EdgeR_EarlyC_OLD_S$strand <- "*"
EdgeR_EarlyC_OLD_E$Test <- "y_C_a_E";  EdgeR_EarlyC_OLD_E$strand <- "*"
EdgeR_EarlyC_OLD_C$Test <- "y_C_a_C";  EdgeR_EarlyC_OLD_C$strand <- "*"

EdgeR_EarlyC_OLD_S$Location <- paste0(EdgeR_EarlyC_OLD_S$chr, ":", EdgeR_EarlyC_OLD_S$start, "-", EdgeR_EarlyC_OLD_S$stop)
EdgeR_EarlyC_OLD_E$Location <- paste0(EdgeR_EarlyC_OLD_E$chr, ":", EdgeR_EarlyC_OLD_E$start, "-", EdgeR_EarlyC_OLD_E$stop)
EdgeR_EarlyC_OLD_C$Location <- paste0(EdgeR_EarlyC_OLD_C$chr, ":", EdgeR_EarlyC_OLD_C$start, "-", EdgeR_EarlyC_OLD_C$stop)

EarlyC_OLD_S <-EdgeR_EarlyC_OLD_S
EarlyC_OLD_E <-EdgeR_EarlyC_OLD_E
EarlyC_OLD_C <-EdgeR_EarlyC_OLD_C

GR_EdgeR_EarlyC_OLD_S <-with(EdgeR_EarlyC_OLD_S, GRanges(chr, IRanges(start, stop), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, Test, Location))
GR_EdgeR_EarlyC_OLD_E <-with(EdgeR_EarlyC_OLD_E, GRanges(chr, IRanges(start, stop), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, Test, Location))
GR_EdgeR_EarlyC_OLD_C <-with(EdgeR_EarlyC_OLD_C, GRanges(chr, IRanges(start, stop), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, Test, Location))

EarlyC_OLD_S <-as.data.frame(GR_EdgeR_EarlyC_OLD_S)
EarlyC_OLD_E <-as.data.frame(GR_EdgeR_EarlyC_OLD_E)
EarlyC_OLD_C <-as.data.frame(GR_EdgeR_EarlyC_OLD_C)

GR_EarlyC_OLD_S <- GR_EdgeR_EarlyC_OLD_S
GR_EarlyC_OLD_E <- GR_EdgeR_EarlyC_OLD_E
GR_EarlyC_OLD_C <- GR_EdgeR_EarlyC_OLD_C

#setwd("/home/fabio/postdoc/IsolationOverlap/gga6/OldvsYoungControl")

Seq_GR_EarlyC_OLD_S = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_EarlyC_OLD_S)
Seq_GR_EarlyC_OLD_E = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_EarlyC_OLD_E)
Seq_GR_EarlyC_OLD_C = BSgenome::getSeq(BSgenome.Ggallus.UCSC.galGal6, GR_EarlyC_OLD_C)

Seq_CG_GR_EarlyC_OLD_S <- vcountPattern("CG", Seq_GR_EarlyC_OLD_S)
Seq_CG_GR_EarlyC_OLD_E <- vcountPattern("CG", Seq_GR_EarlyC_OLD_E)
Seq_CG_GR_EarlyC_OLD_C <- vcountPattern("CG", Seq_GR_EarlyC_OLD_C)

Seq_GC_GR_EarlyC_OLD_S <- vcountPattern("GC", Seq_GR_EarlyC_OLD_S)
Seq_GC_GR_EarlyC_OLD_E <- vcountPattern("GC", Seq_GR_EarlyC_OLD_E)
Seq_GC_GR_EarlyC_OLD_C <- vcountPattern("GC", Seq_GR_EarlyC_OLD_C)

#add CpG and GpC collumn to the Grange objects
CG <-Seq_CG_GR_EarlyC_OLD_S
GC <-Seq_GC_GR_EarlyC_OLD_S
EarlyC_OLD_S <- cbind(EarlyC_OLD_S , CG,GC)
CG <-Seq_CG_GR_EarlyC_OLD_E
GC <-Seq_GC_GR_EarlyC_OLD_E
EarlyC_OLD_E <- cbind(EarlyC_OLD_E, CG,GC)
CG <-Seq_CG_GR_EarlyC_OLD_C
GC <-Seq_GC_GR_EarlyC_OLD_C
EarlyC_OLD_C <- cbind(EarlyC_OLD_C, CG,GC)

#
#I am eliminating this DMR that p-value is too much low compared with the others
#seqnames  start    end width strand edgeR.logFC edgeR.p.value edgeR.adj.p.value    Test           Location CG GC
#chr16    520810    521933  1124      *  1.1672522982  8.035507e-50  6.516796e-47 y_C_a_S      chr16:520810-521933 101  131
#EarlyC_OLD_S <- EarlyC_OLD_S[-c(141),] 

#Genomic Ranges Object
GR_EarlyC_OLD_S <-with(EarlyC_OLD_S, GRanges(seqnames, IRanges(start, end), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, Test, Location, CG, GC))
GR_EarlyC_OLD_E <-with(EarlyC_OLD_E, GRanges(seqnames, IRanges(start, end), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, Test, Location, CG, GC))
GR_EarlyC_OLD_C <-with(EarlyC_OLD_C, GRanges(seqnames, IRanges(start, end), strand, edgeR.logFC, edgeR.p.value, edgeR.adj.p.value, Test, Location, CG, GC))

#Sweden low p-value
#BRSW_EE_S_ROI <- subset(GR_EdgeR_SW_Early_S_ROI, edgeR.p.value<.0005)
#All low p-value 0.0005
#write ranges
#write.table(GR_SW_Early_S_ROI_p.05 , "GR_SW_Early_S_ROI_p.05.txt", sep="\t", col.names=T, row.names=T, quote=F)
#laps - Brasil Early vs Sweden Early
#SW_E_ES_p.05 <-subsetByOverlaps(GR_SW_Early_S_ROI_p.05, GR_SW_Early_E_ROI_p.05)


########################################
#  http://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html 
########################################
library(calibrate)
require(maptools)
setwd("/home/fabio/postdoc/IsolationOverlap/gga6/OldvsYoungControl")

resbind <- c(GR_EarlyC_OLD_S, GR_EarlyC_OLD_C, GR_EarlyC_OLD_E)
resbind <- sortSeqlevels(resbind)

tiff("Vulcano_ECvOall.tiff", units="cm", width=30, height=30, res=600, compression = "lzw")
# Make a basic volcano plot
with(resbind, plot(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="white", main="Volcano plot", xlim=c(-4,4)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(GR_EarlyC_OLD_C, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="red"))

with(subset(GR_EarlyC_OLD_S, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="yellow4"))
with(subset(GR_EarlyC_OLD_E, edgeR.p.value<.9), points(edgeR.logFC, -log10(edgeR.p.value), pch=20, col="cyan4"))

#with(subset(GR_EarlyC_OLD_S, edgeR.p.value<.0005 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
#with(subset(GR_EarlyC_OLD_E, edgeR.p.value<.0005 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
#with(subset(GR_EarlyC_OLD_C, edgeR.p.value<.0005 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))
#  with(subset(resSW_EE, edgeR.p.value<.0005 & abs(edgeR.logFC)>1), pointLabel(edgeR.logFC, -log10(edgeR.p.value), labels=Location, cex=.6, offset = 1))

abline(1.30103, 0, col = "red")
abline(2.30103, 0, col = "yellow")
abline(3.30103, 0, col = "green")
dev.off() 

write.table(as.data.frame(GR_EarlyC_OLD_S), "1_EarlyC_OLD_S.txt", sep="\t", col.names=T, row.names=F, quote=F) #yellow
write.table(as.data.frame(GR_EarlyC_OLD_C), "1_EarlyC_OLD_C.txt", sep="\t", col.names=T, row.names=F, quote=F) #cyan
write.table(as.data.frame(GR_EarlyC_OLD_E), "1_EarlyC_OLD_E.txt", sep="\t", col.names=T, row.names=F, quote=F) #blue

#graph
library(qtl)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)
library(gghighlight)
library(tidyverse)
library(ggbio)

tiff("Chicken_Isolation_manhatahnESC.tiff", units="cm", width=30.00, height=15.00, res=600, compression = "lzw")

autoplot(resbind, coord = "genome", geom = "point", aes(y = -log10(edgeR.p.value), color=Test), space.skip = 0.01, spaceline = TRUE, legend = TRUE) +
  
  scale_color_manual(values=c("red3", "cyan3", "yellow3"))+
  #                                OE       OS        EE      ES        
  #  scale_color_manual(values=c("yellow3", "yellow"))+
  # scale_color_manual(values=c("cyan3", "cyan"))+
  
  geom_hline(yintercept= 3.30103, color='seagreen', size=0.2) +
  geom_hline(yintercept= 2.30103, color='yellow4', size=0.2) +
  geom_hline(yintercept= 1.30103, color='red', size=0.2) +
  
  theme(axis.text.x=element_text(angle=90, hjust=0, size = 6)) 

dev.off()

tiff("Chicken_Foldchange_manhatahn_ESC.tiff", units="cm", width=30.00, height=15.00, res=600, compression = "lzw")

autoplot(resbind, coord = "genome", geom = "point", aes(y = edgeR.logFC, color=Test), space.skip = 0.01, spaceline = TRUE, legend = TRUE) +
  #  scale_color_manual(values=c("cyan3", "yellow3", "cyan", "yellow"))+
  geom_point(color = ifelse(resbind$edgeR.p.value<.05, "black", NA), shape=23, size = 3) +
  #geom_text(aes(label=ifelse(edgeR.adj.p.value<.5, as.character(Location),'')),hjust=0,vjust=0)
  
  #geom_point(aes(if(resbindS$edgeR.p.value<.05)), color=TRUE, shape=23)+
  
  #                                OE       OS        EE      ES        
  #scale_color_manual(values=c("yellow3", "yellow"))+
  
  scale_color_manual(values=c("red3", "cyan3", "yellow3"))+
  
  # geom_hline(yintercept= 3.30103, color='seagreen', size=0.2) +
  #  geom_hline(yintercept= 2.30103, color='yellow4', size=0.2) +
  #  geom_hline(yintercept= 1.30103, color='red', size=0.2) +
  geom_hline(yintercept= 0, color='black', size=0.2) +
  
  
  # geom_point(shape=23, aes(color= sig , size=3))+
  
  theme(axis.text.x=element_text(angle=90, hjust=0, size = 6))
#geom_point(data = ex, mapping = aes(x = Location, y =  -log10(edgeR.p.value)), color="red") +


dev.off()

resbind05 <- subset(resbind, edgeR.p.value<.05)
tiff("Chicken_Foldchange_manhatahn_ESC005.tiff", units="cm", width=30.00, height=15.00, res=600, compression = "lzw")

autoplot(resbind, coord = "genome", geom = "point", aes(y = edgeR.logFC, color = ifelse(resbind$edgeR.p.value<.05, Test, NA)), space.skip = 0.01, spaceline = TRUE, legend = FALSE) +
  #  scale_color_manual(values=c("cyan3", "yellow3", "cyan", "yellow"))+
#  geom_point(color = ifelse(resbind$edgeR.p.value<.05, Test, NA), shape=23, size = 3) +
      #geom_text(aes(label=ifelse(edgeR.adj.p.value<.5, as.character(Location),'')),hjust=0,vjust=0)
   #                                OE       OS        EE      ES        
  #scale_color_manual(values=c("yellow3", "yellow"))+
  
  scale_color_manual(values=c("red3", "cyan3", "yellow3"))+
  
  # geom_hline(yintercept= 3.30103, color='seagreen', size=0.2) +
  #  geom_hline(yintercept= 2.30103, color='yellow4', size=0.2) +
  #  geom_hline(yintercept= 1.30103, color='red', size=0.2) +
  geom_hline(yintercept= 0, color='black', size=0.2) +
  
  
  # geom_point(shape=23, aes(color= sig , size=3))+
  
  theme(axis.text.x=element_text(angle=90, hjust=0, size = 6))
#geom_point(data = ex, mapping = aes(x = Location, y =  -log10(edgeR.p.value)), color="red") +


dev.off()


#Venn Diagram
library(ChIPpeakAnno)
library("Cairo")
GR_EarlyC_OLD_S05<-subset(GR_EarlyC_OLD_S, edgeR.p.value<.05)
GR_EarlyC_OLD_C05<-subset(GR_EarlyC_OLD_C, edgeR.p.value<.05)
GR_EarlyC_OLD_E05<-subset(GR_EarlyC_OLD_E, edgeR.p.value<.05)
resbind05 <- c(GR_EarlyC_OLD_S05, GR_EarlyC_OLD_C05, GR_EarlyC_OLD_E05)
#res <- makeVennDiagram(Peaks=list(peaks1, peaks2, peaks3),
#                       NameOfPeaks=c("TF1", "TF2", "TF3"))

#res <- makeVennDiagram(Peaks=list(GR_EarlyC_OLD_S, GR_EarlyC_OLD_C, GR_EarlyC_OLD_E),
#                                              NameOfPeaks=c("y_C_a_S", "y_C_a_E", "y_C_a_C"))

setwd("/home/fabio/postdoc/IsolationOverlap/gga6/OldvsYoungControl")
setwd("H:/Medicina Veterin?ria/2017.POSDOC-ESALQ-USP/Experimento_Isolation/RresultsAAGGA6/EarlyCvsOlds")

grl <- splitAsList(resbind05, resbind05$Test)
CairoPDF("venn.pdf")
res <- makeVennDiagram(Peaks=grl, NameOfPeaks=names(grl), main="Isolation: Young Controls vs Adults", col = "black",fill = c("red", "blue", "yellow"), alpha = 0.50, cex= 1.5)
dev.off()

resbind05<-subset(resbind, edgeR.p.value<.05)
grl <- splitAsList(resbind05, resbind05$Test)

CairoPDF("venn05.pdf")
res <- makeVennDiagram(Peaks=grl, NameOfPeaks=names(grl), main="Isolation: Young Controls vs Adults P<0.05", col = "black",fill = c("red", "blue", "yellow"), alpha = 0.50, cex= 1.5)
#plot(res)
dev.off()


#A pie chart is used to demonstrate the overlap features of the common peaks.
#pie1(table(ol$overlappingPeaks[["GR_EarlyC_OLD_S05///GR_EarlyC_OLD_C05///GR_EarlyC_OLD_E05"]]$overlapFeature))
ol <- findOverlapsOfPeaks(GR_EarlyC_OLD_S05, GR_EarlyC_OLD_C05, GR_EarlyC_OLD_E05)
ol$venn_cnt
ol$peaklist

CairoPDF("OverlapFeatures.pdf")
for (i in 1:1){
  pie1(table(ol$overlappingPeaks[["GR_EarlyC_OLD_S05///GR_EarlyC_OLD_C05"]]$overlapFeature), main="Isolation: YCvsA_Stress vs Control")
  pie1(table(ol$overlappingPeaks[["GR_EarlyC_OLD_S05///GR_EarlyC_OLD_E05"]]$overlapFeature), main="Isolation: YCvsA_Stress vs Enriched")
  pie1(table(ol$overlappingPeaks[["GR_EarlyC_OLD_C05///GR_EarlyC_OLD_E05"]]$overlapFeature), main="Isolation: YCvsA_Enriched vs Control")
}
dev.off()

#After that I will create objects to each part of the Venn I would like to annotate.
y_C_a_SCE <- ol$peaklist[["GR_EarlyC_OLD_S05///GR_EarlyC_OLD_C05///GR_EarlyC_OLD_E05"]]
y_C_a_S <- ol$peaklist[["GR_EarlyC_OLD_S05"]]
y_C_a_C <- ol$peaklist[["GR_EarlyC_OLD_C05"]]
y_C_a_E <- ol$peaklist[["GR_EarlyC_OLD_E05"]]

y_C_a_SC <- ol$peaklist[["GR_EarlyC_OLD_S05///GR_EarlyC_OLD_C05"]]
y_C_a_EC <- ol$peaklist[["GR_EarlyC_OLD_C05///GR_EarlyC_OLD_E05"]]
y_C_a_ES <- ol$peaklist[["GR_EarlyC_OLD_S05///GR_EarlyC_OLD_E05"]]

#MeSH.Gga.eg.db


#https://bioconductor.org/packages/release/bioc/vignettes/ChIPpeakAnno/inst/doc/ChIPpeakAnno.html
setwd("H:/Medicina Veterin?ria/2017.POSDOC-ESALQ-USP/Experimento_Isolation/RresultsAAGGA6/EarlyCvsOlds")
library ("ChIPseeker")
supportedUCSCtables(genome = "galGal6")
gg_txdb <- makeTxDbFromUCSC(genome="galGal6", tablename="ensGene")

#all the peaks from all the comparisons (i will run these line to each one of the files using common and unique peaks for each treatment)
peak_file <- as.data.frame(resbind05)
peak_file <- as.data.frame(y_C_a_SCE)
peak_file <- as.data.frame(y_C_a_S)
peak_file <- as.data.frame(y_C_a_C)
peak_file <- as.data.frame(y_C_a_E)

peak_file <- as.data.frame(y_C_a_SC)
peak_file <- as.data.frame(y_C_a_EC)
peak_file <- as.data.frame(y_C_a_ES)

peak_file <- subset(peak_file, select=c(seqnames, start,end))
                    names(peak_file) <- c("CHR", "BP", "BP2")

peak_anno_function <- function(input_data){
  b2 <- input_data[,c("CHR" , "BP" , "BP2")]
  write.table(b2 , "annotate_file.txt" , sep = "\t" , col.names = TRUE , row.names = FALSE)
  peakAnno <- annotatePeak("H:/Medicina Veterin?ria/2017.POSDOC-ESALQ-USP/Experimento_Isolation/RresultsAAGGA6/EarlyCvsOlds/annotate_file.txt",
                           TxDb=gg_txdb, annoDb="org.Gg.eg.db")
  return(peakAnno)
}

peak_chicken <- peak_anno_function(peak_file)
plotAnnoPie(peak_chicken)
info_chicken <- peak_chicken@anno@elementMetadata


all <- as.data.frame(peak_file)
y_C_a_SCE <- as.data.frame(peak_file)
y_C_a_S <- as.data.frame(peak_file)
y_C_a_C <- as.data.frame(peak_file)
y_C_a_E <- as.data.frame(peak_file)

y_C_a_SC <- as.data.frame(peak_file)
y_C_a_EC <- as.data.frame(peak_file)
y_C_a_ES <- as.data.frame(peak_file)

#info_chicken$loc.id <- rownames(all$loc.id)
all <- cbind(all, info_chicken , all=T)
#all_info <- merge(all, info_chicken , all=T)

y_C_a_SCE  <- cbind(y_C_a_SCE , info_chicken , all=T)
y_C_a_S <- cbind(y_C_a_S, info_chicken , all=T)
y_C_a_C<- cbind(y_C_a_C, info_chicken , all=T)
y_C_a_E<- cbind(y_C_a_E, info_chicken , all=T)

y_C_a_SC<- cbind(y_C_a_SC, info_chicken , all=T)
y_C_a_EC<- cbind(y_C_a_EC, info_chicken , all=T)
y_C_a_ES<- cbind(y_C_a_ES, info_chicken , all=T)


write.table(all , "y_C_a_All.txt" , sep="\t" , row.names = F , quote = F)
write.table(y_C_a_SCE , "y_C_a_SCE.txt" , sep="\t" , row.names = F , quote = F)
write.table(y_C_a_S , "y_C_a_S.txt" , sep="\t" , row.names = F , quote = F)
write.table(y_C_a_C , "y_C_a_C.txt" , sep="\t" , row.names = F , quote = F)
write.table(y_C_a_E , "y_C_a_E.txt" , sep="\t" , row.names = F , quote = F)

write.table(y_C_a_SC , "y_C_a_SC.txt" , sep="\t" , row.names = F , quote = F)
write.table(y_C_a_EC , "y_C_a_EC.txt" , sep="\t" , row.names = F , quote = F)
write.table(y_C_a_ES , "y_C_a_ES.txt" , sep="\t" , row.names = F , quote = F)


plotAnnoPie(all$annotation)


##ONLY DEG


for_path_DEG <- list( "y_C_a_All"= all$ENTREZID, 
                      "y_C_a_SCE"= y_C_a_SCE$ENTREZID, 
                      "y_C_a_S"=  y_C_a_S$ENTREZID, 
                      "y_C_a_C"=  y_C_a_C$ENTREZID, 
                      "y_C_a_E"= y_C_a_E$ENTREZID)
for_path_DEG <- list( "y_C_a_SCE"= y_C_a_SCE$ENTREZID, 
                      "y_C_a_S"=  y_C_a_S$ENTREZID, 
                      "y_C_a_C"=  y_C_a_C$ENTREZID, 
                      "y_C_a_E"= y_C_a_E$ENTREZID)

for_path_DEG <- list( "y_C_a_SCE"= y_C_a_SCE$ENTREZID, 
                      "y_C_a_S"=  y_C_a_S$ENTREZID, 
                      "y_C_a_C"=  y_C_a_C$ENTREZID, 
                      "y_C_a_E"= y_C_a_E$ENTREZID,
                      "y_C_a_SC"= y_C_a_SC$ENTREZID,
                      "y_C_a_EC"= y_C_a_EC$ENTREZID,
                      "y_C_a_ES"= y_C_a_ES$ENTREZID)   


library(clusterProfiler)

GO_analysis <- compareCluster(geneClusters = for_path_DEG , fun = "enrichGO" , ont="BP" ,OrgDb=org.Gg.eg.db, pvalueCutoff = 0.9, pAdjustMethod = "fdr", readable=T)
jpeg("GO_ANALYSIS_genes_ONLY.jpeg" , width = 900 , height = 850)
dotplot(GO_analysis)
dev.off()
kegg_analysis <- compareCluster(geneClusters = for_path_DEG , fun = "enrichKEGG" ,organism="gga", pvalueCutoff = 0.5, pAdjustMethod = "fdr")
jpeg("KEGG_ANALYSIS_genes_ONLY.jpeg" , width = 600 , height = 850)
dotplot(kegg_analysis)
dev.off()

kegg_analysis <- compareCluster(geneClusters = for_path_DEG , fun = "enrichDO" ,OrgDb=org.Gg.eg.db, pvalueCutoff = 0.5, pAdjustMethod = "fdr")

groupGO_analysis <- compareCluster(geneClusters = for_path_DEG , level=15, fun = "groupGO" , ont="BP" ,OrgDb=org.Gg.eg.db)
jpeg("groupGO_ANALYSIS_genes_ONLY_level15.jpeg" , width = 1200 , height = 850)
dotplot(groupGO_analysis , showCategory=30)
dev.off()




re <- makeVennDiagram(Peaks=ol, NameOfPeaks=names(ol), main="Isolation: Young Controls vs Adults P<0.05", col = "black",fill = c("red", "blue", "yellow"), alpha = 0.50)

library(reshape2)

SC<- subsetByOverlaps(GR_EarlyC_OLD_S, GR_EarlyC_OLD_C)                      
#398
SE<- subsetByOverlaps(GR_EarlyC_OLD_S, GR_EarlyC_OLD_E)
#228
CE<- subsetByOverlaps(GR_EarlyC_OLD_C, GR_EarlyC_OLD_E)
#236

GR_EarlyC_OLD_S05<-subset(GR_EarlyC_OLD_S, edgeR.p.value<.05)
GR_EarlyC_OLD_C05<-subset(GR_EarlyC_OLD_C, edgeR.p.value<.05)
GR_EarlyC_OLD_E05<-subset(GR_EarlyC_OLD_E, edgeR.p.value<.05)
resbind05 <-

SC05<- subsetByOverlaps(GR_EarlyC_OLD_S05, GR_EarlyC_OLD_C05)                      
#88
SE05<- subsetByOverlaps(GR_EarlyC_OLD_S05, GR_EarlyC_OLD_E05)
#42
CE05<- subsetByOverlaps(GR_EarlyC_OLD_C05, GR_EarlyC_OLD_E05)
#49



#
#width of each chromossome coverad by the test
library('regioneR')

#some counts
width(resbind)
width(resbind05)

genome <- BSgenome.Ggallus.UCSC.galGal6
head(seqlengths(genome))
getSeq(genome, as.character=FALSE)

aggregate(width(resbind), by=list(Category=resbind$Test), function(x) c(mean = mean(x), sd = sd(x)))

#Category   x.mean     x.sd
#1  y_C_a_C 348.6082 492.6483
#2  y_C_a_E 512.2944 919.9470
#3  y_C_a_S 362.6235 516.8533

plot(density(as.matrix(resbind$edgeR.p.value)))
plot(density(as.matrix(resbind$edgeR.logFC)))

#size simulation
library('regioneR')
TT<- createRandomRegions(nregions=100000, length.mean=407.84, length.sd=643.15, genome=genome)

TT_df<- as.data.frame (TT)
table(TT_df$seqnames)

resbind_df <- as.data.frame(resbind)
table(resbind_df$seqnames)

resbind05_df <- as.data.frame(resbind05)
table(resbind05_df$seqnames)

params <- new("BSParams", X = Celegans, FUN = alphabetFrequency,
              exclude = c("M", "_"))
bsapply(params)


