library("GenVisR")
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(tidyverse)
library(MutationalPatterns)
library(BSgenome.Mmusculus.UCSC.mm10)

############################################################################
#figure4A
main_layer <- theme_grey()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7,color=1))
custom_pallete <- c("grey30","grey90","steelblue1","yellow","pink","green")
geneslung <- read.delim('drivers_lung.txt',header=F,sep='\t')
colnames(geneslung) <- c('Tumor_Sample_Barcode','chr','pos','ref','alt','Hugo_Symbol','Change','Changep','Variant_Classification')
pdf("Fig-4a.pdf",height=5,width=12)
waterfall(geneslung,mainXlabel=T,main_geneLabSize=12,mainLabelCol = "Changep",mainLabelAngle = 90, mainLabelSize = 2,mainDropMut = TRUE,section_heights = c(0, 1),mainLayer = main_layer,mainPalette=custom_pallete)
dev.off()

############################################################################
#figure4B
genesliv <- read.delim('drivers_liver.txt',header=F,sep='\t')
colnames(genesliv) <- c('Tumor_Sample_Barcode','chr','pos','ref','alt','Hugo_Symbol','Change','Changep','Variant_Classification')
pdf("Fig-4b.pdf",height=5,width=12)
waterfall(genesliv,mainXlabel=T,main_geneLabSize=12,mainLabelCol = "Changep",mainLabelAngle = 90, mainLabelSize = 2,mainDropMut = TRUE,section_heights = c(0, 1),mainLayer = main_layer,mainPalette=custom_pallete)
dev.off()

############################################################################
#figure4C and supplementary Table 9
#test association between specific mutations and mutational signatures
#exposures=NxS exposure matrix with N the number of samples and S the number of signatures 
#withdriver=indices of the samples with the driver mutation
#nodriver=indices of the samples without the driver mutation
#signatures to perform the test. Example: signatures=c('SBS1','SBS5'), only signatures that are 
#present with confidence in the samples are considered
#I tested only mutations with at least 5 samples i a tissue type

DriverAssociationSignatures<-function(exposures,exposures_conf,driversamples,genesel,tissuesel,mutation){
	withdriver <- array(0,dim=length(driversamples))
	for (i in 1:length(driversamples)){withdriver[i] <- str_which(rownames(exposures),driversamples[i])}
	signatures <- which(colSums(exposures_conf[withdriver,])>0)
	ind <- str_which(rownames(exposures),tissuesel) 
	nodriver <- setdiff(ind,withdriver)
	
	pvaltot <- array(0,dim=length(signatures))
	for (i in 1:length(signatures)){
		wtest <- wilcox.test(exposures[withdriver,signatures[i]],exposures[nodriver,signatures[i]],alternative='greater')
		pvaltot[i] <- wtest$p.value
		}
		tablepval <- as.data.frame(pvaltot)
		tablepval <- tablepval %>% mutate(n1=length(withdriver))
		tablepval <- tablepval %>% mutate(n2=length(nodriver))
		tablepval <- tablepval %>% mutate(tissue=tissuesel)
		tablepval <- tablepval %>% mutate(gene=genesel)
		tablepval <- tablepval %>% mutate(change=mutation)
		tablepval <- tablepval %>% mutate(signature=names(signatures))
		return(tablepval)
	}

sbss <- read.delim('../starting_data/SigProfilerMatrixGenerator_matrices/SBS96.all',header=T,sep='\t',check.names=F)
sbss1 <- sbss[,2:182]
rownames(sbss1)<-sbss[,1]

samples1 <- sort(unique(colnames(sbss1)))
samples <- c(samples1[1:36],samples1[42:83],samples1[37:41],samples1[84:181])

ind_basechanges <- c(1:4,25:28,49:52,73:76)
mut_mat_m <- sbss1[rownames(sbss1)[c(ind_basechanges,ind_basechanges+4,ind_basechanges+8,ind_basechanges+12,ind_basechanges+16,ind_basechanges+20)],samples]

rr1 <- readRDS('../figure1/mexposure.rds')
rr1 <- as.matrix(rr1)
rownames(rr1)=samples
qq <- readRDS('../figure1/mexposuresig.rds')
qq <- as.matrix(qq)
rownames(qq)=samples

#test liver
namesp <- genesliv %>% filter(Hugo_Symbol=='Hras' & Changep=='p.Q61K') %>% select(Tumor_Sample_Barcode) %>% mutate(paste('LIVER_',str_replace_all(Tumor_Sample_Barcode,' ','_'),sep=''))
wt <- DriverAssociationSignatures(rr1,qq,namesp[,2],'Hras','LIVER','p.Q61K')
namesp <- genesliv %>% filter(Hugo_Symbol=='Hras' & Changep=='p.Q61L') %>% select(Tumor_Sample_Barcode) %>% mutate(paste('LIVER_',str_replace_all(Tumor_Sample_Barcode,' ','_'),sep=''))
wt1 <- DriverAssociationSignatures(rr1,qq,namesp[,2],'Hras','LIVER','p.Q61L')
totalta <- rbind(wt,wt1)
namesp <- genesliv %>% filter(Hugo_Symbol=='Hras' & Changep=='p.Q61R') %>% select(Tumor_Sample_Barcode) %>% mutate(paste('LIVER_',str_replace_all(Tumor_Sample_Barcode,' ','_'),sep=''))
wt1 <- DriverAssociationSignatures(rr1,qq,namesp[,2],'Hras','LIVER','p.Q61R')
totalta <- rbind(totalta,wt1)
namesp <- genesliv %>% filter(Hugo_Symbol=='Braf' & Changep=='p.V637E') %>% select(Tumor_Sample_Barcode) %>% mutate(paste('LIVER_',str_replace_all(Tumor_Sample_Barcode,' ','_'),sep=''))
wt1 <- DriverAssociationSignatures(rr1,qq,namesp[,2],'Braf','LIVER','p.V637E')
totalta <- rbind(totalta,wt1)
totalta <- totalta %>% mutate(q.value=p.adjust(as.numeric(pvaltot),method='BH'))
totalta <- totalta %>% mutate(tissue = 'LIVER')

#test lung
namesp <- geneslung %>% filter(Hugo_Symbol=='Fgfr2' & Changep=='p.C401R') %>% select(Tumor_Sample_Barcode) %>% mutate(paste('LUNG_',str_replace_all(Tumor_Sample_Barcode,' ','_'),sep=''))
wt <- DriverAssociationSignatures(rr1,qq,namesp[,2],'Fgfr2','LUNG','p.C401R')
namesp <- geneslung %>% filter(Hugo_Symbol=='Kras' & Changep=='p.G12D') %>% select(Tumor_Sample_Barcode) %>% mutate(paste('LUNG_',str_replace_all(Tumor_Sample_Barcode,' ','_'),sep=''))
wt1 <- DriverAssociationSignatures(rr1,qq,namesp[,2],'Kras','LUNG','p.G12D')
totaltb <- rbind(wt,wt1)
namesp <- geneslung %>% filter(Hugo_Symbol=='Braf' & Changep=='p.V637E') %>% select(Tumor_Sample_Barcode) %>% mutate(paste('LUNG_',str_replace_all(Tumor_Sample_Barcode,' ','_'),sep=''))
wt1 <- DriverAssociationSignatures(rr1,qq,namesp[,2],'Braf','LUNG','p.V637E')
totaltb <- rbind(totaltb,wt1)
totaltb <- totaltb %>% mutate(tissue = 'LUNG')
totaltb <- totaltb %>% mutate(q.value=p.adjust(as.numeric(pvaltot),method='BH'))
write.table(rbind(totalta,totaltb),'Supplementary_table_9.txt',quote=F,sep='\t',row.names=F,col.names=T)

############################################################################
#test presence of drivers Supplementary table 8
# I did the test only for chemicals with at least 3 mutations 

samples <- as.character(lapply(samples,function(x) str_replace_all(x,'_',' ')))

DriverAssociationSample<-function(genesel,tissuesel,samples,genes){
	a <- genes %>% filter(Hugo_Symbol==genesel) %>% select(Tumor_Sample_Barcode) %>% mutate(name=paste(tissuesel,Tumor_Sample_Barcode,sep=' ')) %>% select(name)
	a <- as.character(unique(a[,1]))
	b <- as.character(lapply(a,function(x) str_trim(str_sub(x, start = 1L, end =str_locate(x,"\\d$")[1]-2 ))))
	c <- table(b)
	nam=names(which(c>2))
	if (length(nam)>0){
		p=array(0,length(nam))
		n1=array(0,length(nam))
		n2=array(0,length(nam))
		n3=array(0,length(nam))
		n4=array(0,length(nam))
		for (i in seq(1,length(nam))){
			a1=length(str_which(as.character(a),nam[i]))	
			test1 <- fisher.test(matrix(c(a1,(length(str_which(samples,nam[i]))-a1),
			(length(a)-a1),
			(length(str_which(samples,tissuesel))-length(str_which(samples,nam[i]))+a1-length(a))),nrow=2,ncol=2),alternative='greater')
			p[i]=test1$p.value
			n1[i]=a1
			n2[i]=(length(str_which(samples,nam[i]))-a1)
			n3[i]=length(a)-a1
			n4[i]=(length(str_which(samples,tissuesel))-length(str_which(samples,nam[i]))+a1-length(a))
		}	
	}
	tablepval <- as.data.frame(n1)
	tablepval <- tablepval %>% mutate(n2=n2)
	tablepval <- tablepval %>% mutate(n3=n3)
	tablepval <- tablepval %>% mutate(n4=n4)	
	tablepval <- tablepval %>% mutate(p=p)
	tablepval <- tablepval %>% mutate(tissue=tissuesel)
	tablepval <- tablepval %>% mutate(gene=genesel)
	tablepval <- tablepval %>% mutate(chemical=nam)
	return(tablepval)
}

wt1 <- DriverAssociationSample('Kras','LUNG',samples,geneslung)
wt2 <- DriverAssociationSample('Braf','LUNG',samples,geneslung)
wt3 <- DriverAssociationSample('Fgfr2','LUNG',samples,geneslung)

# I decided to test Ctnnb1 for ANTIMONY TRIOXIDE because it is the only chemical with 2 samples out of 3 in lung 
# with mutations
a <- geneslung %>% filter(Hugo_Symbol=='Ctnnb1') %>% select(Tumor_Sample_Barcode) %>% mutate(name=paste('LUNG',Tumor_Sample_Barcode,sep=' ')) %>% select(name)
a <- as.character(unique(a[,1]))
b <- as.character(lapply(a,function(x) str_trim(str_sub(x, start = 1L, end =str_locate(x,"\\d$")[1]-2 ))))
c <- table(b)
nam=names(which(c>1))
wt4=c()
wt4$n1=length(str_which(as.character(a),nam))	
wt4$n2=(length(str_which(samples,nam))-wt4$n1)
wt4$n3=length(a)-wt4$n1
wt4$n4=(length(str_which(samples,'LUNG'))-length(str_which(samples,nam))+wt4$n1-length(a))
wt4$p <- fisher.test(matrix(c(wt4$n1,wt4$n2,wt4$n3,wt4$n4),nrow=2,ncol=2),alternative='greater')$p.value
wt4$tissue <- 'LUNG'
wt4$gene <- 'Ctnnb1'
wt4$chemical <- nam

totalta <- rbind(wt1,wt2,wt3,wt4)
totalta <- totalta %>% mutate(q.value=p.adjust(as.numeric(p),method='BH'))

wt1 <- DriverAssociationSample('Hras','LIVER',samples,genesliv)
wt2 <- DriverAssociationSample('Egfr','LIVER',samples,genesliv)
wt3 <- DriverAssociationSample('Ctnnb1','LIVER',samples,genesliv)
totaltb <- rbind(wt1,wt2,wt3)
totaltb <- totaltb %>% mutate(q.value=p.adjust(as.numeric(p),method='BH'))

write.table(rbind(totalta,totaltb),'Supplementary_table_8.txt',quote=F,sep='\t',row.names=F,col.names=T)

############################################################################
#Figure4D CNVs
cnvs=read.delim('cnvs_anno.txt',header=F)
colnames(cnvs)=c('chromosome','start','end','segmean','sample','sample1','tissue','chemicals')
cnvsk=cnvs %>% filter(tissue=='KIDNEY')
pdf("Fig-4d_kidney.pdf",height=4,width=26)
cnFreq(cnvsk,plot_title='KIDNEY TUMOURS',genome='mm10',plotChr=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19'))+scale_y_continuous(limits = c(-1, 1))
dev.off()

cnvsli=cnvs %>% filter(tissue=='LIVER')
pdf("Fig-4d_liver.pdf",height=4,width=26)
cnFreq(cnvsli,plot_title='LIVER TUMOURS',genome='mm10',plotChr=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19'))+scale_y_continuous(limits = c(-1, 1))
dev.off()

cnvslu=cnvs %>% filter(tissue=='LUNG')
pdf("Fig-4d_lung.pdf",height=4,width=26)
cnFreq(cnvslu,plot_title='LUNG TUMOURS',genome='mm10',plotChr=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19'))+scale_y_continuous(limits = c(-1, 1))
dev.off() 

############################################################################
#supplementary figure 6
#at the end we are considering a total of 132 samples
#for the other we have very few changes so I decided to exclude them
#to perform the clustering I needed to segment the data (I used bedtools and intersectbed).
#I created the segmented_cnvs.txt file using data from Theta2. Theta2 can identify subclonal copy number. 
#In case we had subclonal copy number for a sample, we considered only the biggest clone.

data=read.delim('segmented_cnvs.txt',header=T)
data=data[,4:137]
samples=colnames(data)
myColor <- colorRampPalette(c("blue", "white", "red"))(60)
myBreaks <- c(seq(0, 4, length.out=ceiling(60)))
pdf("Supplementary_figure_6.pdf",height=10,width=26)
pheatmap(data,color=myColor,cluster_rows=F,fontsize_col=7,breaks=myBreaks,clustering_distance_cols='correlation',clustering_method='ward.D',fontsize_row=1)
dev.off() 

############################################################################
#supplementary Table 10a and supplementary Table 10b

DriverAssociationSignatures_ML <- function(namesp,sigs,pps,samples,namesig,trinucleotides,totest){
	ind <- which(samples==namesp)
	j <- which(trinucleotides==totest)
	s1 <- sigs*matrix(rep(pps[ind,],96),96,11,byrow=T)
	ss1 <- s1/rowSums(s1)		
	colnames(ss1) <- namesig
	if (length(which(ss1[j,]>=0.5))>0) {
		l <- names(which(ss1[j,]==max(ss1[j,])))
				}
	else {l <- ''}		
	return(l)
}

DriverAssociationSignatures_ML_total <- function(namesp,sigs,pps,samples,namesig,trinucleotides,totest){
	ind <- which(samples==namesp)
	s1 <- sigs*matrix(rep(pps[ind,],96),96,11,byrow=T)
	ss1 <- s1/rowSums(s1)		
	colnames(ss1) <- namesig
	j <- which(trinucleotides==totest)
	signatureinsample <- names(which(ss1[j,]==max(ss1[j,])))
	return(signatureinsample)
}

complement<-function(base){
	base0 <- vector(mode='character',length(base))
	base0[which(base=='A')] <- 'T'
	base0[which(base=='C')] <- 'G'
	base0[which(base=='G')] <- 'C'
	base0[which(base=='T')] <- 'A'
	return(base0)
}

ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"
trinucleotides <- rownames(mut_mat_m)
summ <- colSums(mut_mat_m)
sigs <- readRDS('../figure1/mSBSs.rds')
sigs <- as.matrix(sigs)
pps <- rr1*matrix(rep(summ,11),nrow=181,ncol=11)
samples <- colnames(mut_mat_m)
namesig <- c("mSBS1","mSBS2","mSBS3" ,"mSBS4","mSBS5","mSBS6","mSBS7","mSBS8","mSBS9","mSBS10","mSBS11")

#upload a file with all the snvs that are in in Tier1 of the cancer gene censusclassified as somatic and having nonsense, missense, splice site or frameshift mutations  (mouse homologous)
genes <- read.delim('snvsindrivers.txt',header=F)
colnames(genes) <- c('sample','MD','chrom','position','ref','alt','gene')
mutation_types <- c('C > A','C > G','C > T','T > A','T > C','T > G')

genes1 <- genes %>% filter(paste(ref,alt,sep=' > ') %in% mutation_types)
genes1 <- genes1 %>% mutate(context = paste(as.character(getSeq(get(ref_genome),paste('chr',as.character(chrom),sep=''),position-1,position-1)),'[',ref,'>',alt,']', as.character(getSeq(get(ref_genome),paste('chr',as.character(chrom),sep=''),position+1,position+1)),sep=''))

genes2 <- genes %>% filter(!(paste(ref,alt,sep=' > ') %in% mutation_types))
genes2 <- genes2 %>% mutate(context = paste(complement(as.character(getSeq(get(ref_genome),paste('chr',as.character(chrom),sep=''),position+1,position+1))),'[',complement(ref),'>',complement(alt),']', complement(as.character(getSeq(get(ref_genome),paste('chr',as.character(chrom),sep=''),position-1,position-1))),sep=''))


genest <- rbind(genes1,genes2)

signature_sig <- vector(mode='character',length=280)
for (i in seq(1:280)){signature_sig[i] <- DriverAssociationSignatures_ML(genest$sample[i],sigs,pps,samples,namesig,trinucleotides,genest$context[i])}
signature <- vector(mode='character',length=280)
for (i in seq(1:280)){signature[i] <- DriverAssociationSignatures_ML_total(genest$sample[i],sigs,pps,samples,namesig,trinucleotides,genest$context[i])}

genest <- cbind(genest,signature,signature_sig)

genest10a <- genest %>% filter(((position==39627783 & alt=='T') | position==130196315 | position==130196323 | position ==141192550 | position==141192551 | (position==145246771 & alt=='T')) & !(str_detect(sample,'STOMACH')))  
write.table(genest10a,'supplementary_table_10a.txt',quote=F,sep='\t',col.names=T,row.names=F)

genest10b <- genest %>% filter(signature_sig!='')  
write.table(genest10b,'supplementary_table_10b.txt',quote=F,sep='\t',col.names=T,row.names=F)

