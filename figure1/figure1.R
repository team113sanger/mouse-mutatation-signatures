library(MutationalPatterns)
library(tidyverse)
library(pheatmap)
library(hdp)
library(ggsci) 

source('../signature_decomposition/signature_decomposition.R')

cancer_signatures60 <- read.table('../starting_data/PCAWG_signatures.txt', sep = "\t", header = TRUE)
cancer_signatures60 <- as.matrix(cancer_signatures60[,2:66])
nm <- colnames(cancer_signatures60)

sbss <- read.delim('../starting_data/SigProfilerMatrixGenerator_matrices/SBS96.all',header=T,sep='\t',check.names=F)
sbss1 <- sbss[,2:182]
rownames(sbss1)<-sbss[,1]

samples1 <- sort(unique(colnames(sbss1)))
samples <- c(samples1[1:36],samples1[42:83],samples1[37:41],samples1[84:181])

ind_basechanges <- c(1:4,25:28,49:52,73:76)

mut_mat_m <- sbss1[rownames(sbss1)[c(ind_basechanges,ind_basechanges+4,ind_basechanges+8,ind_basechanges+12,ind_basechanges+16,ind_basechanges+20)],samples]
samples <- as.character(lapply(samples,function(x) str_replace_all(x,'_',' ')))

############################################################################
#figure 1A
nsnvs <- colSums(mut_mat_m)
totalmut <- as.data.frame(nsnvs)
totalmut <- totalmut %>% mutate(sample=samples)
colnames(totalmut) <- c('number_of_SNVs','sample')

namem_ord <- c("LUNG SPONTANEOUS",
"LUNG ANTIMONY TRIOXIDE",          
"LUNG COBALT METAL",                           
"LUNG ISOBUTYL NITRITE",           
"LUNG NICKEL OXIDE",                
"LUNG NICKEL SUBSULFIDE",          
"LUNG NICKEL SULFATE HEXAHYDRATE",             
"LUNG SODIUM TUNGSTATE DIHYDRATE",  
"LUNG VANADIUM PENTOXIDE",         
"LUNG VINYLIDENE CHLORIDE",
"LIVER SPONTANEOUS", 
"LIVER ANTHRAQUINONE",
"LIVER ANTIMONY TRIOXIDE",
"LIVER BROMOCHLOROACETIC ACID",   
"LIVER COBALT METAL", "LIVER CUMENE","LIVER DIETHANOLAMINE",            
"LIVER FURAN","LIVER GINKGO BILOBA EXTRACT",      
"LIVER INDIUM PHOSPHIDE","LIVER ISOBUTYL NITRITE","LIVER NICKEL OXIDE",              
"LIVER NICKEL SULFATE HEXAHYDRATE",            
"LIVER NICKEL SUBSULFIDE",        
"LIVER OXAZEPAM",                               
"LIVER DE-71",         
"LIVER PRIMACLONE",                
"LIVER SODIUM TUNGSTATE DIHYDRATE",                   
"LIVER VANADIUM PENTOXIDE",          
"LIVER VINYLIDENE CHLORIDE",           
"LIVER 1,2,3 TRICHLOROPROPANE",
"STOMACH 1,2,3 TRICHLOROPROPANE", 
"KIDNEY VINYLIDENE CHLORIDE")

totalmut <- totalmut %>% mutate(category=as.character(lapply(samples,function(x) str_trim(str_sub(x, start = 1L, end =str_locate(x,"\\d$")[1]-2 )))))
totalmut <- totalmut %>% mutate(tissue=as.character(lapply(samples,function(x) str_sub(x, start = 1L, end =str_locate(x,"\\s")[1] ))))

#with all the points
fig1a <- totalmut %>% mutate(name = fct_relevel(category, namem_ord)) %>% ggplot(aes(x=name, y=log10(number_of_SNVs),fill=tissue)) + geom_boxplot(outlier.shape = NA)+geom_jitter(color="black", size=0.5,alpha=0.95)+scale_fill_npg()+theme(axis.text.x=element_text(size = 9, angle = 45, hjust = 0))+theme(axis.text.y=element_text(size = 12))+scale_x_discrete(position = "top")+theme(text = element_text(size = 18))+xlab("tumour type")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+theme(legend.key = element_rect(fill = "white", colour = "white"))
ggsave(plot=fig1a, file='Fig-1a.pdf', device=cairo_pdf, width=12, height=4.5)

############################################################################
#test the difference in the number of snvs for lung and liver
nsnvs <- colSums(mut_mat_m)

#liver substitutions
spoliver <- which((str_detect(names(nsnvs),'LIVER') & str_detect(names(nsnvs),'SPONTANEOUS'))==TRUE)
chem <- c("ANTHRAQUINONE","ANTIMONY_TRIOXIDE","BROMOCHLOROACETIC_ACID",   
"COBALT", "CUMENE","DIETHANOLAMINE",            
"FURAN","GINKGO_BILOBA_EXTRACT",      
"INDIUM_PHOSPHIDE","ISOBUTYL_NITRITE","NICKEL_OXIDE",              
"NICKEL_SULFATE_HEXAHYDRATE",            
"NICKEL_SUBSULFIDE",        
"OXAZEPAM",                               
"DE-71",         
"PRIMACLONE",                
"SODIUM_TUNGSTATE_DIHYDRATE",                   
"VANADIUM_PENTOXIDE",          
"VINYLIDENE_CHLORIDE",           
"1,2,3_TRICHLOROPROPANE")
pval <- rep(0,length(chem))
for (i in seq(1,length(chem))){
chemtestliver <- which((str_detect(names(nsnvs),'LIVER') & str_detect(names(nsnvs),chem[i]))==TRUE)
pval[i] <- wilcox.test(nsnvs[chemtestliver],nsnvs[spoliver],alternative='greater')$p.value
}
qval <- p.adjust(pval,method='BH')
#0.88461538 0.91499340 0.91499340 0.88461538 0.88461538 0.91499340 0.95768267 0.95768267 0.98201981 0.88461538 0.88461538 0.88461538 0.88461538 0.95768267 0.91499340 0.91499340 0.88461538 0.88461538 0.91499340 0.01465201 --> 1,2,3_TRICHLOROPROPANE
qval_liver <- data.frame(chemical=chem,q.value=qval)
write.table(qval_liver,file='liver_substitutions_qval.txt',quote=F,sep='\t',row.names=F)

#lung substitutions
spolung <- which((str_detect(names(nsnvs),'LUNG') & str_detect(names(nsnvs),'SPONTANEOUS'))==TRUE)
chem <- c("ANTIMONY_TRIOXIDE",          
"COBALT",                           
"ISOBUTYL_NITRITE",           
"NICKEL_OXIDE",                
"NICKEL_SUBSULFIDE",          
"NICKEL_SULFATE_HEXAHYDRATE",             
"SODIUM_TUNGSTATE_DIHYDRATE",  
"VANADIUM_PENTOXIDE",         
"VINYLIDENE_CHLORIDE")
pval <- rep(0,length(chem))
for (i in seq(1,length(chem))){
chemtestlung <- which((str_detect(names(nsnvs),'LUNG') & str_detect(names(nsnvs),chem[i]))==TRUE)
pval[i] <- wilcox.test(nsnvs[chemtestlung],nsnvs[spolung],alternative='greater')$p.value
}
qval <- p.adjust(pval,method='BH')
#0.9247468218 0.0004848093 0.9247468218 0.9247468218 0.7689560440 0.9247468218 0.9247468218 0.0197802198 0.6632191338 --> COBALT, VANADIUM_PENTOXIDE
qval_lung <- data.frame(chemical=chem,q.value=qval)
write.table(qval_lung,file='lung_substitutions_qval.txt',quote=F,sep='\t',row.names=F)

#ndinuc liver and lung difference
pp <- which((str_detect(names(nsnvs),'LUNG'))==TRUE)
pp1 <- which((str_detect(names(nsnvs),'LIVER'))==TRUE)
pval <- wilcox.test(nsnvs[pp1],nsnvs[pp])$p.value
#0.0007545472 (1734,1354 median)

############################################################################
#signature extraction with HDP
#since I am performing the novo analysis I do not normalize mouse data to human data, I will normalize later when I need to compare mouse signatures to human signatures
hdp_mut <- hdp_init(ppindex = c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), # index of parental node
                     cpindex = c(1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2), # index of the CP to use
                     hh = rep(1, 96), # prior is uniform over 96 categories
                     alphaa = rep(1,times=2), # shape hyperparameters for 2 CPs
                     alphab = rep(1,times=2))  # rate hyperparameters for 2 CPs
hdp_mut <- hdp_addconparam(hdp_mut,
                            alphaa = rep(1,33), # shape hyperparams for 17 CPs
                            alphab = rep(1,33)) # rate hyperparams for 17 CPs

pp <- c(rep(2,6),rep(3,4),rep(4,5),rep(5,6),rep(6,5), 
    rep(7,5),rep(8,5),rep(9,5),rep(10,5),rep(11,3),
    rep(12,5),rep(13,6),rep(14,4),rep(15,3),rep(16,6),
    rep(17,5),rep(18,5),rep(19,5),rep(20,6),rep(21,11), 
    rep(22,6),rep(23,7),rep(24,5),rep(25,6),rep(26,6),
    rep(27,5),rep(28,4),rep(29,6),rep(30,6),rep(31,12),
    rep(32,3),rep(33,5),rep(34,5))

cp <- pp+1                            

hdp_mut <- hdp_adddp(hdp_mut,
                     numdp = 181, # add 177 nodes
                     ppindex=pp, 
                     cpindex=cp)
                   
hdp_mut <- hdp_setdata(hdp_mut,
                       dpindex=35:215, # which nodes to add data to (leaves) 9+193 202
                       data=rbind(t(mut_mat_m)))# input data matrix

saveRDS(hdp_mut,file='hdp_sbss.rds')

chlist <- vector("list", 8)
############################################################################
# This part is computationally intensive. 
#for (i in 1:8){
#activated <- dp_activate(hdp_mut,dpindex = 1:numdp(hdp_mut),initcc =20,seed = jobIndex*1000)
#chlist[[i]]<- hdp_posterior(activated,burnin = 200000,n = 100,space = 1000,cpiter = 3,seed = 33333*(i+2)+ 1)
# }

#In practice, I run chains in parallel using a FARM job array
# e.g. with a command line like:
##bsub -q long -oo /path/MutData-out181.%I.txt -J 'MYJOB[1-8]' -R'select[mem>MEMVAL] rusage[mem=MEMVAL]' -M MEMVAL 'R CMD BATCH --vanilla path/denovo_181.R path/${LSB_JOBINDEX}.denovo4.181.Rout'

# where the denovo_181.R is :
#jobIndex <- as.numeric(Sys.getenv("LSB_JOBINDEX"))
#hdp_mut <-readRDS('hdp_sbss.rds')
#activated <- dp_activate(hdp_mut,dpindex = 1:numdp(induced),initcc =20,seed = jobIndex*1000)
#chlist<- hdp_posterior(activated,burnin = 200000,n = 100,space = 1000,cpiter = 3,seed = 33333*(jobIndex+2)+ 1)
#saveRDS(chlist,file=paste('hdp_mut_181_denovo.',jobIndex,'.rds',sep=''))

#and the I uploaded the files that I have obtained
#chlist <- vector("list", 8)
############################################################################
chlist[[1]] <- readRDS('hdp_mut_181_denovo.1.rds')
chlist[[2]] <- readRDS('hdp_mut_181_denovo.2.rds')
chlist[[3]] <- readRDS('hdp_mut_181_denovo.3.rds')
chlist[[4]] <- readRDS('hdp_mut_181_denovo.4.rds')
chlist[[5]] <- readRDS('hdp_mut_181_denovo.5.rds')
chlist[[6]] <- readRDS('hdp_mut_181_denovo.6.rds')
chlist[[7]] <- readRDS('hdp_mut_181_denovo.7.rds')
chlist[[8]] <- readRDS('hdp_mut_181_denovo.8.rds')

#combines all in one hdpMultiChain object, checks diagnostic plots, and then
mut_multi <- hdp_multi_chain(chlist)
#extract_components
mut_multi <- hdp_extract_components(mut_multi,min.sample=3)


############################################################################
#Figure 1B
#compare mouse signatures to human signatures 
#I normalize mouse signatures to compare them to the human signatures
tricountHMnorm <- read.delim('../starting_data/mousetohuman_normalization.txt',sep='\t',header = F)
mSBSs <- t(mut_multi@comp_categ_distn$mean[1:11,])
msig <- mSBSs*tricountHMnorm$V2
shdp_norm <- msig/t(matrix(rep(colSums(msig),96),nrow=11,ncol=96))
namesig <- c("mSBS1","mSBS2","mSBS3" ,"mSBS4","mSBS5","mSBS6","mSBS7","mSBS8","mSBS9","mSBS10","mSBS11")
colnames(shdp_norm) <- namesig
csmap <-cos_sim_compare_multiplesignatures(cancer_signatures60,shdp_norm)
cossim <- array(0,11)
bestsig <- matrix('',nrow=11,ncol=1)
for (i in seq(1,11)){
bestsig[i] <- colnames(cancer_signatures60)[which(csmap[,i]==max(csmap[,i]))];cossim[i]=round(max(csmap[,i]),digit=2)}
d <- data.frame('signature'=namesig,'besthumansignature'=bestsig,'cosine_similarity'=cossim)
#we realized that similarity of mSBS10 with SBS17 is lower than the threshold although visually they are very smilar.
#This is mainly due to the fact that SBS17 has been divided in SBS17a and SBS17b but in our case they always appear together.
#we decided to reconstruct the old signature 17 with SBS17a and SBS17b. We report this cosine similarity.  
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")cancer_signatures <- read.table(sp_url, sep = "\t", header = TRUE)
QP <- findSigExposures(as.numeric(cancer_signatures[,17]),cancer_signatures60[,c('SBS17a','SBS17b')])$exposures
d$cosine_similarity[10] <- round(cos_sim(shdp_norm[,10],((0.41*cancer_signatures60[,'SBS17a'])+(0.59*cancer_signatures60[,'SBS17b']))),digit=2)
write.table(d,file='Fig-2b_table.txt',quote=F,sep='\t',row.names=F)
namesig <- c('mSBS5','mSBS40','mSBS19','mSBS_N1','mSBS42','mSBS12','mSBS_N2','mSBS18','mSBS1','mSBS17','mSBS_N3')

#supplementary table 4, deconstruct mouse signature using the minimum number of human signatures, if more than 3 signatures are needed: it is a new signature
weight <- array(0,11)
cossim <- array(0,11)
bestsig <- matrix('',nrow=11,ncol=1)
for (i in seq(1,11)){
	a <- min_best_cos_sim_signature(shdp_norm[,i],cancer_signatures60)
	if (!(is.na(a[1]))){
		weight[i] <- paste(as.character(a[[1]]),' ')
		bestsig[i] <- as.character(tibble(rownames(as.matrix(a[[1]]))))
		cossim[i] <- as.numeric(a[[2]])
	}
	else {bestsig[i] <- 'NA'
	}
}
d <- data.frame('signature'=namesig,'besthumansignature'=bestsig,'weight'=weight,'cosine_similarity'=cossim)
write.table(d,file='Supplementary_table_4.txt',quote=F,sep='\t',row.names=F)

############################################################################
#save some of the data to use
mSBSs <- t(mut_multi@comp_categ_distn$mean[1:11,])
colnames(mSBSs) <- namesig

mexposure <- mut_multi@comp_dp_distn$mean[35:215,]
mexposure <- as.data.frame(mexposure)
rownames(mexposure) <- samples
colnames(mexposure) <- namesig

mexposuresig <- mut_multi@comp_dp_distn$mean;for (i in 2:dim(mut_multi@comp_dp_distn$mean)[1]) {mexposuresig[i,which(mut_multi@comp_dp_distn$cred.int[[i]][1,]==0)]=0}
mexposuresig <- mexposuresig[35:215,]
mexposuresig <- as.data.frame(mexposuresig)
rownames(mexposuresig) <- samples
colnames(mexposuresig) <- namesig

saveRDS(mSBSs,file='mSBSs.rds')
saveRDS(mexposure,file='mexposure.rds')
saveRDS(mexposuresig,file='mexposuresig.rds')


############################################################################
#Figure 1 D
#plot signatures, I used Sigprofiler plot but I think it is easier here to just plot it in R
#plot with correct names
fig1d <- plot_96_profile(mSBSs,ymax=0.12)
ggsave(plot=fig1d, file='Fig-1d.pdf', device=cairo_pdf, width=7, height=10)

############################################################################
#summplementary figure, clustering
myColor <- colorRampPalette(c("white", "blue"))(60)
myBreaks <- c(seq(0, 1, length.out=ceiling(60)))
oo <- pheatmap(t(mexposure),clustering_distance_cols='correlation',cluster_rows=F,clustering_method='complete',fontsize_col=7,fontsize_row=12,breaks=myBreaks,col=myColor)
SupplementaryFigure <- plot_contribution(t(mexposure[oo$tree_col$order,]),coord_flip=FALSE,palette=c(RColorBrewer::brewer.pal(9, "Set1")[c(9)],RColorBrewer::brewer.pal(12, "Paired")))+
    theme(axis.text.x= element_text(angle = 90, hjust = 1, vjust = 0.5,size=7))
ggsave(plot=SupplementaryFigure, file="Supplementary_Figure_3.pdf", device=cairo_pdf, width=16, height=6)

############################################################################
#Figure 1 C
namem <- sort(unique(as.character(lapply(samples,function(x) str_trim(str_sub(x, start = 1L, end =str_locate(x,"\\d$")[1]-2 ))))))
namem <- c(namem[1:7],namem[9:17],namem[8],namem[18:33])
mexposuretotal <- mexposure*matrix(rep(nsnvs,11),nrow=181,ncol=11)
rate <- matrix(0,nrow=11,ncol=33)
mediana <- matrix(0,nrow=11,ncol=33)

for (i in seq(1,33)){
 	indi <- str_which(samples,namem[i])	
 	for (j in seq(1,11)){
 	rate[j,i] <- length(which(mexposure[indi,j]>=0.1))/length(indi)	
 	if (length(which(mexposure[indi,j]>=0.1))>0){
 	mediana[j,i] <- mean(mexposure[indi[which(mexposure[indi,j]>=0.1)],j])}
 	}
 }

rownames(rate) <- namesig
colnames(rate) <- namem
rownames(mediana) <- namesig
colnames(mediana) <- namem

ratesi <- data.frame(t(rate),stringsAsFactors=F,check.names=F)
ratesi <- ratesi %>% mutate(tumour_type=rownames(ratesi))
ratesi.df <- gather(ratesi, namesig, key='signature',value='rate')

medianasi <- data.frame(t(mediana),stringsAsFactors=F,check.names=F)
medianasi <- medianasi %>% mutate(tumour_type=rownames(medianasi))
medianasi.df <- gather(medianasi, namesig, key='signature',value='median_counts')

ratesi.df <- ratesi.df %>% mutate(mecounts=medianasi.df$median_counts)

t1.rect1 <- data.frame (xmin=1.5, xmax=2.5, ymin=0.5, ymax=11.5)
t2.rect1 <- data.frame (xmin=3.5, xmax=4.5, ymin=0.5, ymax=11.5)
t3.rect1 <- data.frame (xmin=5.5, xmax=6.5, ymin=0.5, ymax=11.5)
t4.rect1 <- data.frame (xmin=7.5, xmax=8.5, ymin=0.5, ymax=11.5)
t5.rect1 <- data.frame (xmin=9.5, xmax=10.5, ymin=0.5, ymax=11.5)
t6.rect1 <- data.frame (xmin=11.5, xmax=12.5, ymin=0.5, ymax=11.5)
t7.rect1 <- data.frame (xmin=13.5, xmax=14.5, ymin=0.5, ymax=11.5)
t8.rect1 <- data.frame (xmin=15.5, xmax=16.5, ymin=0.5, ymax=11.5)
t9.rect1 <- data.frame (xmin=17.5, xmax=18.5, ymin=0.5, ymax=11.5)
t10.rect1 <- data.frame (xmin=19.5, xmax=20.5, ymin=0.5, ymax=11.5)
t11.rect1 <- data.frame (xmin=21.5, xmax=22.5, ymin=0.5, ymax=11.5)
t12.rect1 <- data.frame (xmin=23.5, xmax=24.5, ymin=0.5, ymax=11.5)
t13.rect1 <- data.frame (xmin=25.5, xmax=26.5, ymin=0.5, ymax=11.5)
t14.rect1 <- data.frame (xmin=27.5, xmax=28.5, ymin=0.5, ymax=11.5)
t15.rect1 <- data.frame (xmin=29.5, xmax=30.5, ymin=0.5, ymax=11.5)
t16.rect1 <- data.frame (xmin=31.5, xmax=32.5, ymin=0.5, ymax=11.5)

namesig1 <- c('mSBS_N3','mSBS_N2','mSBS_N1','mSBS42','mSBS40','mSBS19','mSBS18','mSBS17','mSBS12','mSBS5','mSBS1')

spot.theme <- list(
    theme_classic(),
    theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 10, angle = 45, hjust = 0)),
    theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 12)),
    theme(axis.line=element_blank()),
    theme(text = element_text(size = 22)),
    theme(panel.background = element_rect(fill = 'white')),
    theme(plot.margin = unit(c(10,10,10,10), "mm")),
    theme(legend.box.background = element_rect(color='white')),
    scale_size_continuous(range = c(-1, 8)),
    scale_colour_gradient2(low = "#F39B7FB2",high = "black",mid = "#3C5488B2",midpoint=0.4,guide = "colourbar",aesthetics = "colour",limits=c(0.08,0.8)),
    scale_x_discrete(position = "top"))


p <- ratesi.df %>% mutate(name = fct_relevel(tumour_type, namem_ord)) %>% mutate(mSBS = fct_relevel(signature, namesig1)) %>% ggplot( aes(x=name, y=mSBS)) + geom_point(aes(colour = mecounts, size = rate))+xlab("tumour type")+ spot.theme

p <- p +
  geom_rect(data=t1.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t2.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t3.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t4.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t5.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t6.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t7.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t8.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t9.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t10.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t11.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t12.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t13.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t14.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t15.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t16.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_hline(yintercept=1.5,color = "white",size=1)+  geom_hline(yintercept=2.5,color = "white",size=1)+geom_hline(yintercept=3.5,color = "white",size=1)+
  geom_hline(yintercept=4.5,color = "white",size=1)+  geom_hline(yintercept=4.5,color = "white",size=1)+geom_hline(yintercept=5.5,color = "white",size=1)+
  geom_hline(yintercept=6.5,color = "white",size=1)+  geom_hline(yintercept=7.5,color = "white",size=1)+geom_hline(yintercept=8.5,color = "white",size=1)+
  geom_hline(yintercept=9.5,color = "white",size=1)+  geom_hline(yintercept=10.5,color = "white",size=1)+
  geom_vline(xintercept=10.5,color ="black",size=1)+geom_vline(xintercept=31.5,color = "black",size=1)+
  geom_vline(xintercept=0.5,color ="black",size=1)+geom_vline(xintercept=33.5,color = "black",size=1)+
geom_vline(xintercept=1.5,color ="grey",size=1)+geom_vline(xintercept=11.5,color = "grey",size=1)	

ggsave(plot=p, file='Fig-1c.pdf', device=cairo_pdf, width=16, height=7)

############################################################################
#Figure 1E TCP
library(gridExtra)
library(grid)
ind <- which(str_detect(samples,'TRICHLOROPROPANE')=='TRUE')
del <- mexposuretotal[ind,]
p1 <- plot_contribution(t(del),mSBSs,mode = "relative",coord_flip=TRUE,palette=c(pal_npg("nrc",alpha=0.8)(10),'#FFFFFF'))+theme(axis.text.x= element_text(angle = 90, hjust = 1, vjust = 0.5,size=7))
p2 <- plot_contribution(t(del),mSBSs,mode = "absolute",coord_flip=TRUE,palette=c(pal_npg("nrc",alpha=0.8)(10),'#FFFFFF'))+theme(axis.text.x= element_text(angle = 90, hjust = 1, vjust = 0.5,size=7))
ggsave(plot=grid.arrange(p1,p2,nrow=1), file='Fig-1e.pdf', device=cairo_pdf, width=12, height=6)

############################################################################
#Figure 1F VINYLIDENE CHLORIDE
ind1 <- which(str_detect(samples,'VINYLIDENE')=='TRUE')
ind1 <- c(1 ,2 ,3 ,4 , 5 ,  6 ,  114, 116, 112, 113, 115, 117, 118, 173, 176, 172,  174 ,175 )
del <- mexposuretotal[ind1,]
p3 <- plot_contribution(t(del),mSBSs,mode = "relative",coord_flip=TRUE,palette=c(pal_npg("nrc",alpha=0.8)(10),'#FFFFFF'))+theme(axis.text.x= element_text(angle = 90, hjust = 1, vjust = 0.5,size=7))
p4 <- plot_contribution(t(del),mSBSs,mode = "absolute",coord_flip=TRUE,palette=c(pal_npg("nrc",alpha=0.8)(10),'#FFFFFF'))+theme(axis.text.x= element_text(angle = 90,hjust = 1, vjust = 0.5,size=7))
ggsave(plot=grid.arrange(p3,p4,nrow=1), file='Fig-1f.pdf', device=cairo_pdf, width=12, height=6)
