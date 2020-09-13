library(tidyverse)
options(stringsAsFactors = F)

#upload the "raw" file that I have generated coming from bedtools, intersection with replication and transcription and signature annotations
replication_time_tsb <- readRDS('../starting_data/snvs_sign_replication_transcription.rds')
replication_time_tsb <- replication_time_tsb %>% mutate(mutationtype=str_sub(context,3,5))

replication_time_tsb_transcribed <- replication_time_tsb %>% filter(transcription_strand =='1' | transcription_strand =='2')
replication_time_tsb_transcribed <-replication_time_tsb_transcribed %>% mutate(transcription_strand=as.numeric(transcription_strand))
replication_time_tsb_transcribed1 <- replication_time_tsb_transcribed %>% filter(ref=='C' | ref=='T')
replication_time_tsb_transcribed2 <- replication_time_tsb_transcribed %>% filter(ref=='A' | ref=='G') %>% mutate(transcription_strand =-(transcription_strand-3))
replication_time_tsb_transcribed=rbind(replication_time_tsb_transcribed1,replication_time_tsb_transcribed2)

mutation_types<-c('C>A','C>G','C>T','T>A','T>C','T>G')
tissues <- c('KIDNEY','LIVER','LUNG','STOMACH')
namesig <- c('mSBS5','mSBS40','mSBS19','mSBS_N1','mSBS42','mSBS12','mSBS_N2','mSBS18','mSBS1','mSBS17','mSBS_N3')

############################################################################
#supplementary table 6a
#analyze transcription strand bias per signature per tissue FIGURE2A and supplementary table 6a
replication_time_tsb_transcribed_tis <- replication_time_tsb_transcribed %>% filter(replicationTime!='.')
total_tsb_all <- as.data.frame(table(replication_time_tsb_transcribed_tis[,c("signature_t0.5","transcription_strand","mutationtype","tissue")]))
transcriptionbias <- c()
for (k in 1:4){
	for (i in 1:11){ 
		for (j in 1:6){
			t1 <- total_tsb_all %>% filter(signature_t0.5==namesig[i]) %>% filter(tissue==tissues[k]) %>% filter(mutationtype==mutation_types[j]) %>% filter(transcription_strand==1)
			t2 <- total_tsb_all %>% filter(signature_t0.5==namesig[i]) %>% filter(tissue==tissues[k]) %>% filter(mutationtype==mutation_types[j]) %>% filter(transcription_strand==2)
			if ((t1$Freq>0 | t2$Freq>0) & (t1$Freq+t2$Freq>=50)){
						pval <- binom.test(c(as.numeric(t1$Freq),as.numeric(t2$Freq)))$p.value
						transcriptionbias <- rbind(transcriptionbias,c(tissues[k],as.character(namesig[i]),mutation_types[j],t1$Freq,t2$Freq,pval))
						}
					}
				}
			}
transcriptionbias <- as.data.frame(transcriptionbias)
colnames(transcriptionbias) <- c('tissue','signature_t0.5','mutationtype','transcribed','nontranscribed','p_value')
transcriptionbias$q_value <- p.adjust(transcriptionbias$p_value,method="BH")
transcriptionbias <- transcriptionbias %>% mutate('enrichment'=(as.numeric(transcribed))/as.numeric(nontranscribed))
write.table(transcriptionbias,'Supplementary_table_6a.txt',quote=F,sep='\t',row.names=F,col.names=T)

############################################################################
#supplementary table 6b
samples1 <- sort(unique(replication_time_tsb$sample))
samples <- c(samples1[1:36],samples1[42:83],samples1[37:41],samples1[84:181])
total_tsb_all <- as.data.frame(table(replication_time_tsb_transcribed[,c("sample","signature_t0.5","transcription_strand","mutationtype")]))
transcriptionbias <- c()
for (k in 1:length(samples)){
	for (i in 1:length(namesig)){
		for (j in 1:6){
			t1 <- total_tsb_all %>% filter(signature_t0.5==namesig[i]) %>% filter(sample==samples[k]) %>% filter(mutationtype==mutation_types[j]) %>% filter(transcription_strand==1)
			t2 <- total_tsb_all %>% filter(signature_t0.5==namesig[i]) %>% filter(sample==samples[k]) %>% filter(mutationtype==mutation_types[j]) %>% filter(transcription_strand==2)
			if ((t1$Freq>0 | t2$Freq>0) & (t1$Freq+t2$Freq>=50)){
						pval <- binom.test(c(as.numeric(t1$Freq),as.numeric(t2$Freq)))$p.value
						transcriptionbias <- rbind(transcriptionbias,c(samples[k],as.character(namesig[i]),mutation_types[j],t1$Freq,t2$Freq,pval))
			}
		}
	}
}
transcriptionbias <- as.data.frame(transcriptionbias)
colnames(transcriptionbias) <- c('sample','signature_t0.5','mutationtype','transcribed','nontranscribed','p_value')
transcriptionbias$q_value <- p.adjust(transcriptionbias$p_value,method="BH")
transcriptionbias <- transcriptionbias %>% mutate('enrichment'=(as.numeric(transcribed))/as.numeric(nontranscribed))
write.table(transcriptionbias,'Supplementary_table_6b.txt',quote=F,sep='\t',row.names=F,col.names=T)

############################################################################
#supplementary table 7a
#replication time bias per signature per tissue FIGURE2C and supplementary table 7a
replication_time_tsb <- replication_time_tsb %>%filter(replicationTime!='.')
replication_time_tsb <- replication_time_tsb %>% mutate(replicationTimeplot=-(as.numeric(replicationTime)-11))

total_tsb_all_rep <- replication_time_tsb[,c("signature_t0.5","replicationTimeplot","mutationtype","tissue")]
total_tsb_all_reptest <- total_tsb_all_rep %>% filter(as.numeric(replicationTimeplot) <4 | as.numeric(replicationTimeplot)>7)
total_tsb_all_reptest$replicationTimeplot[which(total_tsb_all_reptest$replicationTimeplot==2)] <- 1
total_tsb_all_reptest$replicationTimeplot[which(total_tsb_all_reptest$replicationTimeplot==3)] <- 1
total_tsb_all_reptest$replicationTimeplot[which(total_tsb_all_reptest$replicationTimeplot==8)] <- 10
total_tsb_all_reptest$replicationTimeplot[which(total_tsb_all_reptest$replicationTimeplot==9)] <- 10
total_tsb_all_reptest1 <- as.data.frame(table(total_tsb_all_reptest))

replicationtimebias <- c()
for (k in 1:4){
	for (i in 1:length(namesig)){
		for (j in 1:6) {
			t1 <- total_tsb_all_reptest1 %>% filter(signature_t0.5==namesig[i]) %>% filter(tissue==tissues[k]) %>% filter(mutationtype==mutation_types[j]) %>% filter(replicationTimeplot=='1')
			t2 <- total_tsb_all_reptest1 %>% filter(signature_t0.5==namesig[i]) %>% filter(tissue==tissues[k]) %>% filter(mutationtype==mutation_types[j]) %>% filter(replicationTimeplot=='10')
			if ((t1$Freq>0 | t2$Freq>0) & (t1$Freq+t2$Freq>=50)){
				pval <- binom.test(c(as.numeric(t1$Freq),as.numeric(t2$Freq)),p=0.46)$p.value
				replicationtimebias <- rbind(replicationtimebias,c(tissues[k],as.character(namesig[i]),mutation_types[j],t1$Freq,t2$Freq,pval))
				}
			}
		}
	}

replicationtimebias <- as.data.frame(replicationtimebias)
colnames(replicationtimebias) <- c('tissue','signature_t0.5','Mutation_type','early','late','p_value')
replicationtimebias$q_value <- p.adjust(replicationtimebias$p_value,method="BH")
replicationtimebias <- replicationtimebias %>% mutate('enrichment'=(as.numeric(early)/as.numeric(late)))
write.table(replicationtimebias,'Supplementary_table_7a.txt',quote=F,sep='\t',row.names=F,col.names=T)

############################################################################
#supplementary table 7b
total_tsb_all_rep <- replication_time_tsb[,c("sample","signature_t0.5","replicationTimeplot","mutationtype")]
total_tsb_all_reptest <- total_tsb_all_rep %>% filter(as.numeric(replicationTimeplot) <4 | as.numeric(replicationTimeplot)>7)
total_tsb_all_reptest$replicationTimeplot[which(total_tsb_all_reptest$replicationTimeplot==2)] <- 1
total_tsb_all_reptest$replicationTimeplot[which(total_tsb_all_reptest$replicationTimeplot==3)] <- 1
total_tsb_all_reptest$replicationTimeplot[which(total_tsb_all_reptest$replicationTimeplot==8)] <- 10
total_tsb_all_reptest$replicationTimeplot[which(total_tsb_all_reptest$replicationTimeplot==9)] <- 10
total_tsb_all_reptest1 <- as.data.frame(table(total_tsb_all_reptest))

replicationtimebias <- c()
for (k in 1:length(samples)){
	for (i in 1:length(namesig)){ 
		for (j in 1:6) {
			t1 <- total_tsb_all_reptest1 %>% filter(sample==samples[k]) %>% filter(signature_t0.5==namesig[i]) %>% filter(mutationtype==mutation_types[j]) %>% filter(replicationTimeplot=='1')
			t2 <- total_tsb_all_reptest1 %>% filter(sample==samples[k]) %>% filter(signature_t0.5==namesig[i]) %>% filter(mutationtype==mutation_types[j]) %>% filter(replicationTimeplot=='10')
			if ((t1$Freq>0 | t2$Freq>0) & ((t1$Freq+t2$Freq)>=50)){
				pval <- binom.test(c(as.numeric(t1$Freq),as.numeric(t2$Freq)),p=0.46)$p.value
				replicationtimebias <- rbind(replicationtimebias,c(as.character(samples[k]),as.character(namesig[i]),mutation_types[j],t1$Freq,t2$Freq,pval))
				}
			}
		}
	}

replicationtimebias <- as.data.frame(replicationtimebias)
colnames(replicationtimebias) <- c('sample','signature_t0.5','Mutation_type','early','late','p_value')
replicationtimebias$q_value <- p.adjust(replicationtimebias$p_value,method="BH")
replicationtimebias <- replicationtimebias %>% mutate('enrichment'=(as.numeric(early)/as.numeric(late)))
write.table(replicationtimebias,'Supplementary_table_7b.txt',quote=F,sep='\t',row.names=F,col.names=T)

############################################################################
#figures 2a and 2C
#read the files that I have saved
limits2 <- function(qval) {
	sig <- vector(mode='numeric',length=length(qval))
	for (i in seq(1,length(qval))){
		if (qval[i]>0.05 ){sig[i] <- 0}
		else if (qval[i]>0.01 & qval[i]<=0.05 ){sig[i] <- 0.25}		
		else if (qval[i]>0.0001 & qval[i]<=0.01 ){sig[i] <- 0.5}	
		else if (qval[i]<=0.0001 ){sig[i]=1}	
	}
	return(sig)
}

consider <- c('LUNG_mSBS1',
'LUNG_mSBS5',
'LUNG_mSBS12',
'LUNG_mSBS17',
'LUNG_mSBS18',
'LUNG_mSBS40',
'LUNG_mSBS_N1',
'LUNG_mSBS_N3',
'LIVER_mSBS1',
'LIVER_mSBS5',
'LIVER_mSBS12',
'LIVER_mSBS17',
'LIVER_mSBS18',
'LIVER_mSBS19',
'LIVER_mSBS40',
'LIVER_mSBS_N1',
'LIVER_mSBS_N3',
'STOMACH_mSBS_N2',
'STOMACH_mSBS19',
'STOMACH_mSBS42',
'KIDNEY_mSBS12',
'KIDNEY_mSBS_N1')

consider <- rev(consider)

rep <- read.delim('Supplementary_table_7a.txt')
tra <- read.delim('Supplementary_table_6a.txt')

tra <- tra %>% mutate(type=paste(tissue,signature_t0.5,sep='_'))
rep <- rep %>% mutate(type=paste(tissue,signature_t0.5,sep='_'))
tra <- tra %>% mutate(significance=limits2(q_value))
tra <- tra %>% filter(type %in% consider)
rep <- rep %>% mutate(significance=limits2(q_value))
rep <- rep %>% filter(type %in% consider)

t1.rect1 <- data.frame (xmin=1.5, xmax=2.5, ymin=0.5, ymax=22.5)
t2.rect1 <- data.frame (xmin=3.5, xmax=4.5, ymin=0.5, ymax=22.5)
t3.rect1 <- data.frame (xmin=5.5, xmax=6.5, ymin=0.5, ymax=22.5)

pdf('Fig-2a.pdf',width=10, height=11)
spot.theme <- list(
    theme_classic(),
    theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 13)),
    theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 12)),
    theme(axis.line=element_blank()),
    theme(text = element_text(size = 22)),
    theme(panel.background = element_rect(fill = 'white')),
    theme(plot.margin = unit(c(10,10,10,10), "mm")),
    theme(legend.box.background = element_rect(color='white')),
    scale_size_continuous(range = c(-1, 8)),
    scale_colour_gradient2(low = "blue",high = "red",	mid = "white",midpoint=0,guide = "colourbar",aesthetics = "colour",limits=c(-1,1)),
    scale_x_discrete(position = "top"))

p <- tra %>% mutate(name=fct_relevel(type, consider)) %>% ggplot(aes(x=mutationtype, y=name)) + geom_point(aes(colour = log2(enrichment), size = significance))+spot.theme+xlab("")+ylab("")

p +
  geom_rect(data=t1.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t2.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t3.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_hline(yintercept=1.5,color = "white",size=1)+  geom_hline(yintercept=2.5,color = "white",size=1)+geom_hline(yintercept=3.5,color = "white",size=1)+
  geom_hline(yintercept=4.5,color = "white",size=1)+  geom_hline(yintercept=4.5,color = "white",size=1)+geom_hline(yintercept=5.5,color = "white",size=1)+
  geom_hline(yintercept=6.5,color = "white",size=1)+  geom_hline(yintercept=7.5,color = "white",size=1)+geom_hline(yintercept=8.5,color = "white",size=1)+
  geom_hline(yintercept=9.5,color = "white",size=1)+  geom_hline(yintercept=10.5,color = "white",size=1)+
  geom_hline(yintercept=11.5,color = "white",size=1)+  geom_hline(yintercept=12.5,color = "white",size=1)+geom_hline(yintercept=13.5,color = "white",size=1)+
  geom_hline(yintercept=14.5,color = "white",size=1)+  geom_hline(yintercept=15.5,color = "white",size=1)+geom_hline(yintercept=16.5,color = "white",size=1)+
  geom_hline(yintercept=17.5,color = "white",size=1)+  geom_hline(yintercept=18.5,color = "white",size=1)+geom_hline(yintercept=19.5,color = "white",size=1)+
  geom_hline(yintercept=20.5,color = "white",size=1)+  geom_hline(yintercept=21.5,color = "white",size=1)+
  geom_vline(xintercept=0.5,color ="gray",size=1)+geom_vline(xintercept=6.5,color = "gray",size=1)
dev.off()

pdf('Fig-2c.pdf',width=10, height=11)
spot.theme <- list(
    theme_classic(),
    theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 13)),
    theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 12)),
    theme(axis.line=element_blank()),
    theme(text = element_text(size = 22)),
    theme(panel.background = element_rect(fill = 'white')),
    theme(plot.margin = unit(c(10,10,10,10), "mm")),
    theme(legend.box.background = element_rect(color='white')),
    scale_size_continuous(range = c(-1, 8)),
    scale_colour_gradient2(low = "blue",high = "red",mid = "white",midpoint=0,guide = "colourbar",aesthetics = "colour",limits=c(-3.1,1)),
    scale_x_discrete(position = "top"))

p <-  rep %>% mutate(name=fct_relevel(type, consider)) %>% ggplot(aes(x=Mutation_type, y=name)) + geom_point(aes(colour = log2(enrichment), size = significance))+spot.theme+xlab("")+ylab("")
p +
  geom_rect(data=t1.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t2.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t3.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_hline(yintercept=1.5,color = "white",size=1)+  geom_hline(yintercept=2.5,color = "white",size=1)+geom_hline(yintercept=3.5,color = "white",size=1)+
  geom_hline(yintercept=4.5,color = "white",size=1)+  geom_hline(yintercept=4.5,color = "white",size=1)+geom_hline(yintercept=5.5,color = "white",size=1)+
  geom_hline(yintercept=6.5,color = "white",size=1)+  geom_hline(yintercept=7.5,color = "white",size=1)+geom_hline(yintercept=8.5,color = "white",size=1)+
  geom_hline(yintercept=9.5,color = "white",size=1)+  geom_hline(yintercept=10.5,color = "white",size=1)+
  geom_hline(yintercept=11.5,color = "white",size=1)+  geom_hline(yintercept=12.5,color = "white",size=1)+geom_hline(yintercept=13.5,color = "white",size=1)+
  geom_hline(yintercept=14.5,color = "white",size=1)+  geom_hline(yintercept=15.5,color = "white",size=1)+geom_hline(yintercept=16.5,color = "white",size=1)+
  geom_hline(yintercept=17.5,color = "white",size=1)+  geom_hline(yintercept=18.5,color = "white",size=1)+geom_hline(yintercept=19.5,color = "white",size=1)+
  geom_hline(yintercept=20.5,color = "white",size=1)+  geom_hline(yintercept=21.5,color = "white",size=1)+
  geom_vline(xintercept=0.5,color ="gray",size=1)+geom_vline(xintercept=6.5,color = "gray",size=1)
dev.off()

############################################################################
#figure 2b examples of transcriptional strand bias, samples are selected from
#../starting_data/plots/SBS_384_plots_ntp_snvs.pdf
#this pdf file contains the transcriptional stand bias spectra for all the samples

#figure 2d samples with mSBS_N3, samples are selected from 
#../starting_data/plots/SBS_96_plots_ntp_snvs.pdf
#this pdf file contains the spectra for all the samples

############################################################################
#figures 2E AID motif 
library(MutationalPatterns)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggseqlogo)
library(gridExtra)
library(grid)

sbss <- read.delim('../starting_data/SigProfilerMatrixGenerator_matrices/SBS96.all',header=T,sep='\t',check.names=F)
samples1 <- sort(unique(colnames(sbss[,2:182])))
samples <- c(samples1[1:36],samples1[42:83],samples1[37:41],samples1[84:181])

uni <- replication_time_tsb %>% select(chrom,position,ref,alt,sample)
uni <- uni %>% mutate(chrom=paste('chr',as.character(chrom),sep=''))

signatures <-readRDS('../figure1/mSBSs.rds')
signatures <- as.matrix(signatures)
activitiessig=readRDS('../figure1/mexposuresig.rds')
activitiessig <- as.matrix(activitiessig)
rownames(activitiessig)=samples

consider <- rownames(activitiessig[which(activitiessig[,11]>0),])
ref_genome <- "BSgenome.Mmusculus.UCSC.mm10"

uni1 <- uni %>% filter(sample %in% consider)
uni1 <- uni1 %>% mutate(vcf_context = as.character(getSeq(get(ref_genome),chrom,position-2,position+ 2)))
uni1 <- uni1 %>% mutate(vcf_context_small = as.character(getSeq(get(ref_genome),chrom,position-1,position+ 1)))
uni2 <- uni1 %>% mutate(vcf_context_smallr=reverse(chartr('ATGC', 'TACG', vcf_context_small))) %>% mutate(vcf_contextr=reverse(chartr('ATGC','TACG',vcf_context)))  

#I will filter in different parts already here
uni_GCC_CG <- uni1 %>% filter(ref=='C' & alt=='G' & vcf_context_small=='GCC')
uni_GCC_CGr <- uni2 %>% filter(ref=='G' & alt=='C' & vcf_context_smallr=='GCC')
context_GCC_CG <- c(uni_GCC_CG$vcf_context,uni_GCC_CGr$vcf_contextr)

uni_GCC_CT <- uni1 %>% filter(ref=='C' & alt=='T' & vcf_context_small=='GCC')
uni_GCC_CTr <- uni2 %>% filter(ref=='G' & alt=='A' & vcf_context_smallr=='GCC')
context_GCC_CT <- c(uni_GCC_CT$vcf_context,uni_GCC_CTr$vcf_contextr)

sel3 <- as.data.frame(context_GCC_CG)%>%mutate(sml=substr(context_GCC_CG, 1, 4))
fig2e1 <- ggseqlogo(sel3$sml,method = 'prob')

sel3 <- as.data.frame(context_GCC_CT)%>%mutate(sml=substr(context_GCC_CT, 1, 4))
fig2e2<- ggseqlogo(sel3$sml,method = 'prob')
ggsave(plot=grid.arrange(fig2e1,fig2e2,nrow=1), file='Fig-2e.pdf', device=cairo_pdf, width=12, height=4)

