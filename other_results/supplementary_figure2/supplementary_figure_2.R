library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(regioneR)
library(karyoploteR)

uni=readRDS('../starting_data/snvs.rds')
samples=levels(factor(uni$sample))
ref_genome = "BSgenome.Mmusculus.UCSC.mm10"

vargs=GRangesList() 
vargs1=GRangesList() 

#|(alt == 'C' )
for (i in 1:181){
uni1= uni %>% filter(sample == samples[i])  
#uni1= uni1 %>% filter((ref == 'C') & (alt == 'T' ))
uni1= uni1 %>% filter(((ref == 'C') & (alt == 'A' )) | ((ref == 'C') & (alt == 'G' )) | ((ref == 'C') & (alt == 'T' )) | ((ref == 'T') & (alt == 'A' )) | ((ref == 'T') & (alt == 'C' )) | ((ref == 'T') & (alt == 'G' )) )
varg=makeGRangesFromDataFrame(uni1,
                              keep.extra.columns=TRUE,
                              ignore.strand=TRUE,
                              seqinfo=NULL,
                              seqnames.field="chrom",
                              start.field="pos",
                              end.field="pos",
                              starts.in.df.are.0based=FALSE)
vargs[[i]]=varg
}

for (i in 1:181){
uni1= uni %>% filter(sample == samples[i])  
#uni1= uni1 %>% filter((ref == 'G') & (alt == 'A' ))
uni1= uni1 %>% filter(((ref == 'G') & (alt == 'T' )) | ((ref == 'G') & (alt == 'C' )) | ((ref == 'G') & (alt == 'A' )) | ((ref == 'A') & (alt == 'T' )) | ((ref == 'A') & (alt == 'G' )) | ((ref == 'A') & (alt == 'C' )) )
varg1=makeGRangesFromDataFrame(uni1,
                              keep.extra.columns=TRUE,
                              ignore.strand=TRUE,
                              seqinfo=NULL,
                              seqnames.field="chrom",
                              start.field="pos",
                              end.field="pos",
                              starts.in.df.are.0based=FALSE)
vargs1[[i]]=varg1
}


#https://bernatgel.github.io/karyoploter_tutorial//Examples/Rainfall/Rainfall.html

pdf(file="supplementary_figure_2.pdf", onefile=TRUE,height=11,width=17) 

pp <- getDefaultPlotParams(plot.type = 4)
pp$data1inmargin <- 0
pp$bottommargin <- 20

#in furan there is "something for some samples". However it is in the same region for all the samples. This cluster is made up of a group of putative false positive mutations 
#I think they are false positive mutations because this region is characterized by many mutations in common between samples, so they are filtered out. However, few remains. 
for (i in 1:181)
{
sm.gr =vargs1[[i]] 
variant.colors <- getVariantsColors(sm.gr$ref, sm.gr$alt)
kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL,labels.plotter = NULL, plot.params = pp, genome='mm10')
kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt=45)
kpAddMainTitle(kp, main=unique(vargs[[i]]$sample), cex=1.2)

kpPlotRainfall(kp, data = sm.gr, col=variant.colors, r0=0, r1=0.34)
kpAxis(kp, ymax = 7, tick.pos = 1:7, r0=0, r1=0.34)
kpAddLabels(kp, labels = c("G>N and A>N"), srt=90, pos=1, label.margin = 0.03, r0=0, r1=0.34)
kpAddLabels(kp, labels = c("Distance between mutations (log10)"), srt=90, pos=1, label.margin = 0.04, r0=0, r1=0.34)
kpPlotDensity(kp, data = sm.gr, r0=0.35, r1=0.50)
kpAddLabels(kp, labels = c("G>N and A>N"), srt=90, pos=1, label.margin = 0.03, r0=0.35, r1=0.50)
kpAddLabels(kp, labels = c("Density"), srt=90, pos=1, label.margin = 0.04, r0=0.35, r1=0.50)

sm.gr =vargs[[i]]
variant.colors <- getVariantsColors(sm.gr$ref, sm.gr$alt)

kpPlotRainfall(kp, data = sm.gr, col=variant.colors, r0=0.51, r1=0.85)
kpAxis(kp, ymax = 7, tick.pos = 1:7, r0=0.51, r1=0.85)
kpAddLabels(kp, labels = c("C>N and T>N"), srt=90, pos=1, label.margin = 0.03, r0=0.51, r1=0.85)
kpAddLabels(kp, labels = c("Distance between mutations (log10)"), srt=90, pos=1, label.margin = 0.04, r0=0.51, r1=0.85)
kpPlotDensity(kp, data = sm.gr, r0=0.86, r1=1)
kpAddLabels(kp, labels = c("C>N and T>N"), srt=90, pos=1, label.margin = 0.03, r0=0.86, r1=1)
kpAddLabels(kp, labels = c("Density"), srt=90, pos=1, label.margin = 0.04, r0=0.86, r1=1)
legend("bottomleft",inset=.02,c("C>A & G>T","C>G & G>C","C>T & G>A","T>A & A>T","T>C & A>G","T>C & A>C"), fill=c("#4c64ae","#000000","#e40611","#bf4a96","#fbe800","#6eb529"), horiz=F, cex=0.8)
}

dev.off()
















