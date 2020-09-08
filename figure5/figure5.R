#Laura Riva April 2020

library(tidyverse)
library(ggsci)
library(MutationalPatterns)

#identify mutational signatures in human cancers. Data are downloaded from 
#Alexandrov et al., The repertoire of mutational signatures in human cancer, Nature 2020
#accession code syn11801889, available at https://www.synapse.org/#!Synapse:syn11801889
#In particular, I downloaded the following data (on September 13 2019):
#Mutation Catalogs -- Spectra of Individual Tumours: WES_Other.96.csv, WES_TCGA.96.csv, WGS_Other.96.csv, WGS_PCAWG.96.csv 
#SP_Signatures_in_Samples: nonPCAWG_WES_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv, TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv,
#nonPCAWG_WGS_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv, PCAWG_sigProfiler_SBS_signatures_in_samples.csv

############################################################################
#upload the data 
source('../signature_decomposition/signature_decomposition.R')
cancer_signatures60 <- read.table('../starting_data/PCAWG_signatures.txt', sep = "\t", header = TRUE)
cancer_signatures60 <- as.matrix(cancer_signatures60[,2:66])
nm <- colnames(cancer_signatures60)
nhdp <- readRDS('mousesignatures_norm.rds')

names <- c("A[C>A]A", "A[C>A]C",
"A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C",
"G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "A[C>G]A", "A[C>G]C",
"A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A", "G[C>G]C",
"G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C",
"A[C>T]G", "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C",
"G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T", "A[T>A]A", "A[T>A]C",
"A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "G[T>A]A", "G[T>A]C",
"G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "A[T>C]A", "A[T>C]C",
"A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C",
"G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T", "A[T>G]A", "A[T>G]C",
"A[T>G]G", "A[T>G]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C",
"G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")

humanWGS1 <- read.csv('WGS_PCAWG.96.csv',check.names=FALSE)
humanWGS2 <- read.csv('WGS_Other.96.csv',check.names=FALSE)
humanWGS <- cbind(humanWGS1[,3:2782],humanWGS2[,3:1867])
rownames(humanWGS) <- names
#I upload this file because I used a specific order of the samples when I did the analysis 
namessample <- read.delim('nameorder.txt',header=F,check.names=FALSE)
namessample <- as.character(namessample$V1)
humandata <- humanWGS[,namessample]

activeWGS1 <- read.csv('PCAWG_sigProfiler_SBS_signatures_in_samples.csv',check.names=FALSE)
activeWGS2 <- read.csv('nonPCAWG_WGS_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv',check.names=FALSE)
activeWGS <- rbind(activeWGS1,activeWGS2)
colnames(activeWGS) <- c('Cancer.Types','Sample.Names',colnames(activeWGS)[3:68])
activeWGS <- activeWGS %>% mutate(name=paste(Cancer.Types,Sample.Names,sep='::'))
rownames(activeWGS) <- activeWGS$name
activeWGS <- activeWGS[namessample,]
active_signatures <- activeWGS[,1:68]

activeWES1 <- read.csv('nonPCAWG_WES_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv',check.names=FALSE)
activeWES2 <- read.csv('TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv',check.names=FALSE)
humanWES1 <- read.csv('WES_Other.96.csv',check.names=FALSE)
humanWES2 <- read.csv('WES_TCGA.96.csv',check.names=FALSE)

humanWES <- cbind(humanWES1[,3:9693],humanWES2[,3:9495])
rownames(humanWES) <- names
namessample <- sort(colnames(humanWES))
humanWES <- humanWES[,namessample]

activeWES <- rbind(activeWES1,activeWES2)
colnames(activeWES) <- c('Cancer.Types','Sample.Names',colnames(activeWES)[3:68])
activeWES <- activeWES %>% mutate(name=paste(Cancer.Types,Sample.Names,sep='::'))
rownames(activeWES) <- activeWES$name
activeWES <- activeWES[namessample,]
activeWES <- activeWES[,1:68]

sigtocons <- c('SBS1','SBS2','SBS3','SBS4','SBS5','SBS6','SBS7a','SBS7b','SBS7c','SBS7d','SBS8','SBS9','SBS10a','SBS10b','SBS11','SBS12','SBS13','SBS14','SBS15','SBS16','SBS17a',
'SBS17b','SBS18','SBS19','SBS20','SBS21','SBS22','SBS23','SBS24','SBS25','SBS26','SBS28','SBS29','SBS30','SBS31','SBS32','SBS33','SBS34','SBS35','SBS36','SBS37','SBS38','SBS39',
'SBS40','SBS41','SBS42','SBS44')

sigtest <- cbind(cancer_signatures60[,sigtocons],nhdp[,4],nhdp[,7])
colnames(sigtest) <- c(sigtocons,'SBS_N1','SBS_N2')

#use SignatureEstimation to identify SBS19, SBS42, mSBS_N1 and mSBS_N2 in human samples
############################################################################
##this part is computationally intensive and was done on the farm, dividing the original WES and WGS spectra in parts
##WGS
#wgs <- matrix(0,4645,ncol=49)
#colnames(wgs) <- colnames(sigtest)
#for (i in seq(1,4645)){
#	boot <- bootstrapSigExposures(humandata[,i]/sum(humandata[,i]), sigtest, R = 10000, mutation.count = sum(humandata[,i]))
#	signature_fractionhdpT <- apply(as.matrix(boot$exposures),1,function(x) 1-(1+length(which(x>0.05)))/10001)
#	wgs[i,] <- signature_fractionhdpT
#}
#saveRDS(wgs,file='WGS_0.05.rds')
##WES
#wes <- matrix(0,nrow=19184,ncol=49)
#colnames(wes) <- colnames(sigtest)
#for (i in seq(1,19184)){
#	print(i)
#	boot <- bootstrapSigExposures(humanWES[,i]/sum(humanWES[,i]), sigtest, R = 10000, mutation.count = sum(humanWES[,i]))
#	signature_fractionhdpT<-apply(as.matrix(boot$exposures),1,function(x) 1-(1+length(which(x>0.05)))/10001)
#	wes[i,] <- signature_fractionhdpT
#}
#saveRDS(wes,file='WES_0.05.rds')
############################################################################

wgs_0.05 <- readRDS('wgs_0.05.rds')
wes_0.05 <- readRDS('wes_0.05.rds')

tumwgs_19_0.05 <- which(wgs_0.05[,'SBS19']==0 & colSums(humandata)>=200)
tumwgs_42_0.05 <- which(wgs_0.05[,'SBS42']==0 & colSums(humandata)>=200)
tumwgs_N1_0.05 <- which(wgs_0.05[,'SBS_N1']==0 & colSums(humandata)>=200)
tumwgs_N2_0.05 <- which(wgs_0.05[,'SBS_N2']==0 & colSums(humandata)>=200)

tumwes_19_0.05 <- which(wes_0.05[,'SBS19']==0 & colSums(humanWES)>=200)
tumwes_42_0.05 <- which(wes_0.05[,'SBS42']==0 & colSums(humanWES)>=200)
tumwes_N1_0.05 <- which(wes_0.05[,'SBS_N1']==0 & colSums(humanWES)>=200)
tumwes_N2_0.05 <- which(wes_0.05[,'SBS_N2']==0 & colSums(humanWES)>=200)

wgs_0.05_19 <- matrix(1,nrow=length(tumwgs_19_0.05),ncol=1)
rownames(wgs_0.05_19) <- names(tumwgs_19_0.05)
for (i in seq(1,length(tumwgs_19_0.05))){
	wgs_0.05_19[i] <- findSigExposures(humandata[,tumwgs_19_0.05[i]], sigtest)$exposures['SBS19',]
}
wes_0.05_19 <- matrix(1,nrow=length(tumwes_19_0.05),ncol=1)
rownames(wes_0.05_19) <- names(tumwes_19_0.05)
for (i in seq(1,length(tumwes_19_0.05))){
	wes_0.05_19[i] <- findSigExposures(humanWES[,tumwes_19_0.05[i]], sigtest)$exposures['SBS19',]
}
nwgs <- length(tumwgs_19_0.05)
nwes <- length(tumwes_19_0.05)
pp <- rbind(cbind(wgs_0.05_19*colSums(humandata[,tumwgs_19_0.05]),rep('WGS',nwgs),rep('SBS19',nwgs)),cbind(wes_0.05_19*colSums(humanWES[,tumwes_19_0.05]),rep('WES',nwes),rep('SBS19',nwes)))

wgs_0.05_42 <- matrix(1,nrow=length(tumwgs_42_0.05),ncol=1)
rownames(wgs_0.05_42) <- names(tumwgs_42_0.05)
for (i in seq(1,length(tumwgs_42_0.05))){
	wgs_0.05_42[i] <- findSigExposures(humandata[,tumwgs_42_0.05[i]], sigtest)$exposures['SBS42',]
}
wes_0.05_42 <- matrix(1,nrow=length(tumwes_42_0.05),ncol=1)
rownames(wes_0.05_42) <- names(tumwes_42_0.05)
for (i in seq(1,length(tumwes_42_0.05))){
	wes_0.05_42[i] <- findSigExposures(humanWES[,tumwes_42_0.05[i]], sigtest)$exposures['SBS42',]
}
nwgs <- length(tumwgs_42_0.05)
nwes <- length(tumwes_42_0.05)
pp <- rbind(pp,cbind(wgs_0.05_42*sum(humandata[,tumwgs_42_0.05]),rep('WGS',nwgs),rep('SBS42',nwgs)),cbind(wes_0.05_42*colSums(humanWES[,tumwes_42_0.05]),rep('WES',nwes),rep('SBS42',nwes)))

wgs_0.05_N1 <- matrix(1,nrow=length(tumwgs_N1_0.05),ncol=1)
rownames(wgs_0.05_N1) <- names(tumwgs_N1_0.05)
for (i in seq(1,length(tumwgs_N1_0.05))){
	wgs_0.05_N1[i] <- findSigExposures(humandata[,tumwgs_N1_0.05[i]], sigtest)$exposures['SBS_N1',]
}
nwgs <- length(tumwgs_N1_0.05)
pp <- rbind(pp,cbind(wgs_0.05_N1*colSums(humandata[,tumwgs_N1_0.05]),rep('WGS',nwgs),rep('mSBS_N1',nwgs)))

wgs_0.05_N2 <- matrix(1,nrow=length(tumwgs_N2_0.05),ncol=1)
rownames(wgs_0.05_N2) <- names(tumwgs_N2_0.05)
for (i in seq(1,length(tumwgs_N2_0.05))){
	wgs_0.05_N2[i] <- findSigExposures(humandata[,tumwgs_N2_0.05[i]], sigtest)$exposures['SBS_N2',]
}
wes_0.05_N2 <- matrix(1,nrow=length(tumwes_N2_0.05),ncol=1)
rownames(wes_0.05_N2) <- names(tumwes_N2_0.05)
for (i in seq(1,length(tumwes_N2_0.05))){
	wes_0.05_N2[i] <- findSigExposures(humanWES[,tumwes_N2_0.05[i]], sigtest)$exposures['SBS_N2',]
}
nwgs <- length(tumwgs_N2_0.05)
nwes <- length(tumwes_N2_0.05)
pp <- rbind(pp,cbind(wgs_0.05_N2*colSums(humandata[,tumwgs_N2_0.05]),rep('WGS',nwgs),rep('mSBS_N2',nwgs)),cbind(wes_0.05_N2*colSums(humanWES[,tumwes_N2_0.05]),rep('WES',nwes),rep('mSBS_N2',nwes)))

d <- data.frame('mutations'=as.numeric(pp[,1]),'sequencing'=as.factor(pp[,2]),'treatment'=as.factor(pp[,3]),Sample=rownames(pp))
d <- d %>%mutate(tumour=str_extract(Sample,':'))
d <- d %>%mutate(tumour=str_split_fixed(Sample,':',n=2)[,1])
d <- d %>%mutate(tumour=str_replace_all(d$tumour,'CA','Ca'))
d <- d %>%mutate(tumour=str_replace_all(d$tumour,'Sarcoma-bone','Bone-Sarcoma'))
d <- d %>%mutate(tumour=as.factor(tumour))

#Figure 5a
############################################################################
# This script plots the data for Figure 5A as a series of Pie Charts
# Set the run parameters:
params <- list(
    'count.scale' = c('WES'=50, 'WGS'=2500),
    'col.palette' = pal_npg(palette='nrc', alpha=1)(6),
    'shape.palette' = c(15, 16, 17),
    'treatment.order' = c('SBS19', 'SBS42', 'mSBS_N1', 'mSBS_N2'),
    'log.y' = TRUE
)

# Load the cancer types:
cancer.types <- read.table('cancer-types.txt', sep='\t', header=TRUE, row.names=1, colClasses=c('character', rep('factor', 3)))

# Load the raw data:
d$mutations.mb <- d$mutations / params$count.scale[as.character(d$sequencing)]
d$treatment <- factor(d$treatment, levels=params$treatment.order) # Reorder treatments into the correct order

# Reorder the individual samples in each treatment by increasing mutation rate:
d <- do.call(rbind, lapply(levels(d$treatment), function(treatment){
    x <- d[d$treatment == treatment,]
    x <- x[order(x$mutations.mb), ]
    x$index <- 1:nrow(x)
    return(x)
}))

# Build the shape and colour ramps for the different cancer types:
point.types <- expand.grid('colour'=params$col.palette, 'shape'=params$shape.palette, stringsAsFactors=FALSE)
if(nlevels(cancer.types$tissue) > nrow(point.types)) stop('not enough rows in point.types!')
point.types <- point.types[1:nlevels(cancer.types$tissue),]
rownames(point.types) <- levels(cancer.types$tissue)

# Define the colour and shape vectors:
type.colours <- point.types$colour
type.shapes <- point.types$shape
names(type.colours) <- names(type.shapes) <- rownames(point.types)

# Set the cancer tissue and types:
d$tissue <- cancer.types[as.character(d$tumour), 'tissue']

# Build and generate the output plot:
g <- ggplot(d, aes(x=index, y=mutations.mb, colour=tissue, shape=tissue))
g <- g + geom_point()
if(identical(params$log.y, TRUE)){
    g <- g + scale_y_log10(expand=expand_scale(mult=0.1))
} else {
    g <- g + scale_y_continuous(expand=expand_scale(mult=0.1))
}
g <- g + facet_grid(rows=vars(treatment))
g <- g + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
g <- g + theme(panel.background=element_blank(), panel.border=element_rect(fill=NA))
g <- g + scale_colour_manual(values=type.colours)
g <- g + scale_shape_manual(values=type.shapes)
g <- g + labs(y=expression(Mutations/Mb), x='')
g <- g + theme(legend.position='right', legend.direction='vertical', legend.key=element_blank(), legend.key.size=unit(0.1, 'cm'), strip.background=element_blank())
ggsave(plot=g, file='figure-5a.pdf', device=cairo_pdf, width=6, height=3)

############################################################################
#I used SigprofilerPlotting in python, however here I obtain the same results using MutationalPatterns
#for Figure 5b 
#I want to show mSBS19 and some liver samples where I found with SBS19 
signb=cbind(nhdp[,3],humandata[,c('Liver-HCC::SP97771','Liver-HCC::SP112141')])
colnames(signb)=c('mSB19','Liver-HCC:SP97771','Liver-HCC:SP112141')
fig5b<-plot_96_profile(signb,ymax=0.09)
ggsave(plot=fig5b, file='figure-5b.pdf', device=cairo_pdf, width=7, height=5)

#for Figure 5c 
#I want to show mSBS42 and some liver samples where I found with SBS42
signc=cbind(nhdp[,5],humanWES[,c('Biliary-AdenoCa::BD114T','Biliary-AdenoCa::Case4(CHCOSK005)')])
colnames(signc)=c('mSB42','Biliary-AdenoCa:BD114T','Biliary-AdenoCa:Case4(CHCOSK005)')
fig5c<-plot_96_profile(signc,ymax=0.17)
ggsave(plot=fig5c, file='figure-5c.pdf', device=cairo_pdf, width=7, height=5)

############################################################################
#Supplementary Table 11a
tablespectra <- cbind(humandata[,tumwgs_19_0.05],humanWES[,tumwes_19_0.05],
humandata[,tumwgs_42_0.05],humanWES[,tumwes_42_0.05],humandata[,tumwgs_N1_0.05],
humandata[,tumwgs_N2_0.05],humanWES[,tumwes_N2_0.05])
names <- unique(colnames(tablespectra))
tablespectra <- tablespectra[,names]
tablespectra <- tablespectra %>%mutate(Mutation_Type=rownames(tablespectra))
write.table(tablespectra,file='Supplementary_Table_11a.txt',quote=F,sep='\t',row.names=F,col.names=T)

############################################################################
#Supplementary Table 11b
tableact <- rbind(active_signatures[tumwgs_19_0.05,],
activeWES[tumwes_19_0.05,],
active_signatures[tumwgs_42_0.05,],
activeWES[tumwes_42_0.05,],
active_signatures[tumwgs_N1_0.05,],
active_signatures[tumwgs_N2_0.05,],
activeWES[tumwes_N2_0.05,])

signfound <- c(rep('SBS19_wgs',dim(active_signatures[tumwgs_19_0.05,])[1]),
rep('SBS19_wxs',dim(active_signatures[tumwes_19_0.05,])[1]),
rep('SBS42_wgs',dim(active_signatures[tumwgs_42_0.05,])[1]),
rep('SBS42_wxs',dim(active_signatures[tumwes_42_0.05,])[1]),
rep('mSBS_N1_wgs',dim(active_signatures[tumwgs_N1_0.05,])[1]),
rep('mSBS_N2_wgs',dim(active_signatures[tumwgs_N2_0.05,])[1]),
rep('mSBS_N2_wxs',dim(active_signatures[tumwes_N2_0.05,])[1]))

tableact <- tableact%>%mutate(signatures=signfound)
write.table(tableact,file='Supplementary_Table_11b.txt',quote=F,sep='\t',row.names=F,col.names=T)

############################################################################
#Supplementary Table 12
#I had to change the names of the samples manually, because I realized that samples with different names were part of the same tumour type
#for example cns-piloastro and cns-lgg were considered oif the same tumour type, but it was not clear from the starting data
#I considered Classification_of_cancer_types_for_WES_data.xlsx file from 

tumwgs_19_0.05 <- table(tolower(active_signatures$Cancer.Types[which(wgs_0.05[,'SBS19']==0 & colSums(humandata)>=200)]))
tumwgs_42_0.05 <- table(tolower(active_signatures$Cancer.Types[which(wgs_0.05[,'SBS42']==0 & colSums(humandata)>=200)]))
tumwgs_N1_0.05 <- table(tolower(active_signatures$Cancer.Types[which(wgs_0.05[,'SBS_N1']==0 & colSums(humandata)>=200)]))
tumwgs_N2_0.05 <- table(tolower(active_signatures$Cancer.Types[which(wgs_0.05[,'SBS_N2']==0 & colSums(humandata)>=200)]))

tumwes_19_0.05 <- table(tolower(activeWES$Cancer.Types[which(wes_0.05[,'SBS19']==0 & colSums(humanWES)>=200)]))
tumwes_42_0.05 <- table(tolower(activeWES$Cancer.Types[which(wes_0.05[,'SBS42']==0 & colSums(humanWES)>=200)]))
tumwes_N1_0.05 <- table(tolower(activeWES$Cancer.Types[which(wes_0.05[,'SBS_N1']==0 & colSums(humanWES)>=200)]))
tumwes_N2_0.05 <- table(tolower(activeWES$Cancer.Types[which(wes_0.05[,'SBS_N2']==0 & colSums(humanWES)>=200)]))

tum_19_0.05 <- table(c(tolower(active_signatures$Cancer.Types[which(wgs_0.05[,'SBS19']==0 & colSums(humandata)>=200)]),tolower(activeWES$Cancer.Types[which(wes_0.05[,'SBS19']==0 & colSums(humanWES)>=200)])))
tum_42_0.05 <- table(c(tolower(active_signatures$Cancer.Types[which(wgs_0.05[,'SBS42']==0 & colSums(humandata)>=200)]),tolower(activeWES$Cancer.Types[which(wes_0.05[,'SBS42']==0 & colSums(humanWES)>=200)])))
tum_N1_0.05 <- table(c(tolower(active_signatures$Cancer.Types[which(wgs_0.05[,'SBS_N1']==0 & colSums(humandata)>=200)]),tolower(activeWES$Cancer.Types[which(wes_0.05[,'SBS_N1']==0 & colSums(humanWES)>=200)])))
tum_N2_0.05 <- table(c(tolower(active_signatures$Cancer.Types[which(wgs_0.05[,'SBS_N2']==0 & colSums(humandata)>=200)]),tolower(activeWES$Cancer.Types[which(wes_0.05[,'SBS_N2']==0 & colSums(humanWES)>=200)])))

selwes <- activeWES[which(colSums(humanWES)>=200),]
selwgs <- active_signatures[which(colSums(humandata)>=200),]
totalnames <- table(sort(c(tolower(selwgs$Cancer.Types),tolower(selwes$Cancer.Types))))
totalsamples <- sum(totalnames)
#8079 in total

pt1 <- array(0,dim=29)
#SBS19
for (i in seq(1,length(tum_19_0.05))){
test1 <- fisher.test(matrix(c(tum_19_0.05[names(tum_19_0.05)[i]],
(sum(tum_19_0.05)-tum_19_0.05[names(tum_19_0.05)[i]]),
(totalnames[names(tum_19_0.05)[i]]-tum_19_0.05[names(tum_19_0.05)[i]]),
(totalsamples-sum(tum_19_0.05)-totalnames[names(tum_19_0.05)[i]]+tum_19_0.05[names(tum_19_0.05)[i]])),nrow=2,ncol=2),alternative='greater')
pt1[i] <- test1$p.value}
#SBS42
for (i in seq(1,length(tum_42_0.05))){
test1 <- fisher.test(matrix(c(tum_42_0.05[names(tum_42_0.05)[i]],
(sum(tum_42_0.05)-tum_42_0.05[names(tum_42_0.05)[i]]),
(totalnames[names(tum_42_0.05)[i]]-tum_42_0.05[names(tum_42_0.05)[i]]),
(totalsamples-sum(tum_42_0.05)-totalnames[names(tum_42_0.05)[i]]+tum_42_0.05[names(tum_42_0.05)[i]])),nrow=2,ncol=2),alternative='greater')
j=length(tum_19_0.05)+i
pt1[j] <- test1$p.value}
#SBS_N1
for (i in seq(1,length(tum_N1_0.05))){
test1 <- fisher.test(matrix(c(tum_N1_0.05[names(tum_N1_0.05)[i]],
(sum(tum_N1_0.05)-tum_N1_0.05[names(tum_N1_0.05)[i]]),
(totalnames[names(tum_N1_0.05)[i]]-tum_N1_0.05[names(tum_N1_0.05)[i]]),
(totalsamples-sum(tum_N1_0.05)-totalnames[names(tum_N1_0.05)[i]]+tum_N1_0.05[names(tum_N1_0.05)[i]])),nrow=2,ncol=2),alternative='greater')
j=length(tum_19_0.05)+length(tum_42_0.05)+i
pt1[j] <- test1$p.value}
#SBS_N2
for (i in seq(1,length(tum_N2_0.05))){
test1 <- fisher.test(matrix(c(tum_N2_0.05[names(tum_N2_0.05)[i]],
(sum(tum_N2_0.05)-tum_N2_0.05[names(tum_N2_0.05)[i]]),
(totalnames[names(tum_N2_0.05)[i]]-tum_N2_0.05[names(tum_N2_0.05)[i]]),
(totalsamples-sum(tum_N2_0.05)-totalnames[names(tum_N2_0.05)[i]]+tum_N2_0.05[names(tum_N2_0.05)[i]])),nrow=2,ncol=2),alternative='greater')
j=length(tum_19_0.05)+length(tum_42_0.05)+length(tum_N1_0.05)+i
pt1[j] <- test1$p.value}

qt1 <- p.adjust(pt1,method='BH')
names0.05 <- c(names(tum_19_0.05),names(tum_42_0.05),names(tum_N1_0.05),names(tum_N2_0.05))
namesig <- c(rep('SBS19',length(tum_19_0.05)),rep('SBS42',length(tum_42_0.05)),rep('SBS_N1',length(tum_N1_0.05)),rep('SBS_N2',length(tum_N2_0.05)))
t1 <-data.frame(signature=namesig,cancertype=names0.05,nsample=c(tum_19_0.05,tum_42_0.05,tum_N1_0.05,tum_N2_0.05),
ntumours=c(totalnames[names(tum_19_0.05)],totalnames[names(tum_42_0.05)],totalnames[names(tum_N1_0.05)],totalnames[names(tum_N2_0.05)]),pval=pt1,qval=qt1)

write.table(t1,file='Supplementary_Table_12.txt',quote=F,sep='\t',row.names=F,col.names=T)








