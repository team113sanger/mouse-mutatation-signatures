library("GenVisR")
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(tidyverse)

#figure4A
main_layer <- theme_grey()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7,color=1))
custom_pallete <- c("grey30","grey90","steelblue1","yellow","pink","green")

genes=read.delim('drivers_lung.txt',header=F,sep='\t')
colnames(genes)=c('Tumor_Sample_Barcode','chr','pos','ref','alt','Hugo_Symbol','Change','Changep','Variant_Classification')
waterfall(genes,mainXlabel=T,main_geneLabSize=12,mainLabelCol = "Changep",mainLabelAngle = 90, mainLabelSize = 2,mainDropMut = TRUE,section_heights = c(0, 1),mainLayer = main_layer,mainPalette=custom_pallete)

#figure4B
genesliv=read.delim('drivers_liver.txt',header=F,sep='\t')
colnames(genesliv)=c('Tumor_Sample_Barcode','chr','pos','ref','alt','Hugo_Symbol','Change','Changep','Variant_Classification')
waterfall(genesliv,mainXlabel=T,main_geneLabSize=12,mainLabelCol = "Changep",mainLabelAngle = 90, mainLabelSize = 2,mainDropMut = TRUE,section_heights = c(0, 1),mainLayer = main_layer,mainPalette=custom_pallete)

#supplementary figure, clustering of CNVs
data=read.delim('segmented_cnvs.txt',header=T)
data=data[,4:137]
samples=colnames(data)
myColor <- colorRampPalette(c("blue", "white", "red"))(60)
myBreaks <- c(seq(0, 4, length.out=ceiling(60)))
pheatmap(data,color=myColor,cluster_rows=F,fontsize_col=7,breaks=myBreaks,clustering_distance_cols='correlation',clustering_method='ward.D',fontsize_row=1)

#Figure4D
cnvs=read.delim('cnvs_anno.txt',header=F)
colnames(cnvs)=c('chromosome','start','end','segmean','sample','sample1','tissue','chemicals')
cnvsk=cnvs %>% filter(tissue=='KIDNEY')
pdf("cnvskidney.pdf",height=4,width=26)
cnFreq(cnvsk,plot_title='KIDNEY TUMOURS',genome='mm10',plotChr=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19'))+scale_y_continuous(limits = c(-1, 1))
dev.off()

cnvsli=cnvs %>% filter(tissue=='LIVER')
pdf("cnvsliver.pdf",height=4,width=26)
cnFreq(cnvsli,plot_title='LIVER TUMOURS',genome='mm10',plotChr=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19'))+scale_y_continuous(limits = c(-1, 1))
dev.off()

cnvslu=cnvs %>% filter(tissue=='LUNG')
pdf("cnvslung.pdf",height=4,width=26)
cnFreq(cnvslu,plot_title='LUNG TUMOURS',genome='mm10',plotChr=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19'))+scale_y_continuous(limits = c(-1, 1))
dev.off() 


