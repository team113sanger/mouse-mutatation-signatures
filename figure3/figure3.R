library(MutationalPatterns)
library(tidyverse)
library(hdp)
library(ggsci) 
options(stringsAsFactors = F)
source('../signature_decomposition/signature_decomposition.R')

############################################################################
#dinucleotides
dinuc <- read.delim('../starting_data/SigProfilerMatrixGenerator_matrices/DBS78.all',header=T,sep='\t',check.names=F)
dinuc1 <- dinuc[,2:182]
rownames(dinuc1) <- dinuc[,1]

samples1 <- sort(unique(colnames(dinuc1)))
samples <- c(samples1[1:36],samples1[42:83],samples1[37:41],samples1[84:181])
dinuc2 <- dinuc1[,samples]

samples <- as.character(lapply(samples,function(x) str_replace_all(x,'_',' ')))

############################################################################
#figure 3A
totaldinuc <- colSums(dinuc2)
totaldinuc <- as.data.frame(totaldinuc)
totaldinuc <- totaldinuc %>% mutate(sample=samples)
colnames(totaldinuc) <- c('number_of_dinucleotides','sample')

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

totaldinuc <- totaldinuc %>% mutate(category=as.character(lapply(samples,function(x) str_trim(str_sub(x, start = 1L, end =str_locate(x,"\\d$")[1]-2 )))))
totaldinuc <- totaldinuc %>% mutate(tissue=as.character(lapply(samples,function(x) str_sub(x, start = 1L, end =str_locate(x,"\\s")[1] ))))

fig3a <- totaldinuc %>% mutate(name = fct_relevel(category, namem_ord)) %>% ggplot(aes(x=name, y=number_of_dinucleotides,fill=tissue)) + geom_boxplot(outlier.shape = NA)+geom_jitter(color="black", size=0.5,alpha=0.95)+xlab("")+scale_fill_npg()+theme_light()+theme(axis.text.x=element_text(size = 9, angle = 45, hjust = 0))+theme(axis.text.y=element_text(size = 12))+scale_x_discrete(position = "top")+theme(text = element_text(size = 18))+xlab("tumour type")+ylab("number of dinucleotides")+theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+theme(legend.key = element_rect(fill = "white", colour = "white"))
ggsave(plot=fig3a, file='Figure-3a.pdf', device=cairo_pdf, width=11, height=5)

############################################################################
#test the difference in the number of snvs for lung and liver
ndinuc <- colSums(dinuc2)

#ndinuc liver
spoliver <- which((str_detect(names(ndinuc),'LIVER') & str_detect(names(ndinuc),'SPONTANEOUS'))==TRUE)
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
chemtestliver <- which((str_detect(names(ndinuc),'LIVER') & str_detect(names(ndinuc),chem[i]))==TRUE)
pval[i] <- wilcox.test(ndinuc[chemtestliver],ndinuc[spoliver],alternative='greater')$p.value
}
qval <- p.adjust(pval,method='BH')
#0.9059340 0.9312023 0.9312023 0.9312023 0.9312023 0.9312023 0.9312023 0.9312023 0.9312023 0.9312023 0.9312023 0.9312023 0.9312023 0.9312023 0.9312023 0.9312023 0.9059340 0.9312023 0.9059340 0.9059340
qval_liver <- data.frame(chemical=chem,q.value=qval)
write.table(qval_liver,file='liver_dinucleotides_qval.txt',quote=F,sep='\t',row.names=F)

#ndinuc lung
spolung <- which((str_detect(names(ndinuc),'LUNG') & str_detect(names(ndinuc),'SPONTANEOUS'))==TRUE)
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
chemtestlung <- which((str_detect(names(ndinuc),'LUNG') & str_detect(names(ndinuc),chem[i]))==TRUE)
pval[i] <- wilcox.test(ndinuc[chemtestlung],ndinuc[spolung],alternative='greater')$p.value
}
qval <- p.adjust(pval,method='BH')
#0.697291551 0.003811582 0.513234391 0.697291551 0.814941591 0.814941591 0.697291551 0.129651900 0.697291551 --> COBALT
qval_lung <- data.frame(chemical=chem,q.value=qval)
write.table(qval_lung,file='lung_dinucleotides_qval.txt',quote=F,sep='\t',row.names=F)

#ndinuc liver and lung difference
pp <- which((str_detect(names(ndinuc),'LUNG'))==TRUE)
pp1 <- which((str_detect(names(ndinuc),'LIVER'))==TRUE)
wilcox.test(ndinuc[pp1],ndinuc[pp])$p.value
#3.216488e-12 (16,6 median)

############################################################################
#figure 3B
dinuc3 <- as.data.frame(t(dinuc2))
colnames(dinuc3) <-c('AC.CA','AC.CG','AC.CT','AC.GA','AC.GG','AC.GT','AC.TA','AC.TG','AC.TT','AT.CA','AT.CC','AT.CG','AT.GA','AT.GC','AT.TA','CC.AA','CC.AG','CC.AT','CC.GA','CC.GG','CC.GT','CC.TA','CC.TG','CC.TT','CG.AT','CG.GC','CG.GT','CG.TA','CG.TC','CG.TT','CT.AA','CT.AC','CT.AG','CT.GA','CT.GC','CT.GG','CT.TA','CT.TC','CT.TG','GC.AA','GC.AG','GC.AT','GC.CA','GC.CG','GC.TA','TA.AT','TA.CG','TA.CT','TA.GC','TA.GG','TA.GT','TC.AA','TC.AG','TC.AT','TC.CA','TC.CG','TC.CT','TC.GA','TC.GG','TC.GT','TG.AA','TG.AC','TG.AT','TG.CA','TG.CC','TG.CT','TG.GA','TG.GC', 'TG.GT','TT.AA','TT.AC','TT.AG','TT.CA','TT.CC','TT.CG','TT.GA','TT.GC','TT.GG')

dinuc3 <- dinuc3%>%mutate(AC=AC.CA+AC.CG+AC.CT+AC.GA+AC.GG+AC.GT+AC.TA+AC.TG+AC.TT)
dinuc3 <- dinuc3%>%mutate(AT=AT.CA+AT.CC+AT.CG+AT.GA+AT.GC+AT.TA)
dinuc3 <- dinuc3%>%mutate(CC=CC.AA+CC.AG+CC.AT+CC.GA+CC.GG+CC.GT+CC.TA+CC.TG+CC.TT)
dinuc3 <- dinuc3%>%mutate(CG=CG.AT+CG.GC+CG.GT+CG.TA+CG.TC+CG.TT)
dinuc3 <- dinuc3%>%mutate(CT=CT.AA+CT.AC+CT.AG+CT.GA+CT.GC+CT.GG+CT.TA+CT.TC+CT.TG)
dinuc3 <- dinuc3%>%mutate(GC=GC.AA+GC.AG+GC.AT+GC.CA+GC.CG+GC.TA)
dinuc3 <- dinuc3%>%mutate(TA=TA.AT+TA.CG+TA.CT+TA.GC+TA.GG+TA.GT )
dinuc3 <- dinuc3%>%mutate(TC=TC.AA+TC.AG+TC.AT+TC.CA+TC.CG+TC.CT+TC.GA+TC.GG+TC.GT)
dinuc3 <- dinuc3%>%mutate(TG=TG.AA+TG.AC+TG.AT+TG.CA+TG.CC+TG.CT+TG.GA+TG.GC+TG.GT )
dinuc3 <- dinuc3%>%mutate(TT=TT.AA+TT.AC+TT.AG+TT.CA+TT.CC+TT.CG+TT.GA+TT.GC+TT.GG)

dinuc4 <-data.frame(Sample=samples,'AC'=dinuc3['AC'],'AT'=dinuc3['AT'],'CC'=dinuc3['CC'],'CG'=dinuc3['CG'],'CT'=dinuc3['CT'],'GC'=dinuc3['GC'],'TA'=dinuc3['TA'],'TC'=dinuc3['TC'],'TG'=dinuc3['TG'],'TT'=dinuc3['TT'])
#rownames(dinuc4)=samples
namem <- sort(unique(as.character(lapply(samples,function(x) str_trim(str_sub(x, start = 1L, end =str_locate(x,"\\d$")[1]-2 ))))))
namem <- c(namem[1:7],namem[9:17],namem[8],namem[18:33])

dinuc4 <- dinuc4 %>% mutate(category=as.character(lapply(samples,function(x) str_trim(str_sub(x, start = 1L, end =str_locate(x,"\\d$")[1]-2 )))))
#str_replace(x,"\\s\\d$",''))))
dinuc4 <- dinuc4 %>% mutate(tissue=as.character(lapply(samples,function(x) str_sub(x, start = 1L, end =str_locate(x,"\\s")[1] ))))
dinuc4 <- dinuc4 %>% mutate(total=AC+AT+CC+CG+CT+GC+TA+TC+TG+TT)
dinucshow <- gather(dinuc4,Dinucleotides,Contribution, AC:TT)

mm <- matrix(0,nrow=33,ncol=10)
dinuctypes <- unique(dinucshow$Dinucleotides)
for (i in seq(1,33)){
	for (j in seq(1,10)){
		cc <- dinucshow %>% filter(category==namem[i]) %>% filter(Dinucleotides==dinuctypes[j]) %>% select(Contribution)
		mm[i,j] <- mean(cc$Contribution)
	}
}

rownames(mm) <- namem
colnames(mm) <- dinuctypes
mm <- as.data.frame(mm)
mm <- mm %>% mutate(sample=namem)
dinucat <- gather(mm,Dinucleotides,Contribution, AC:TT)

fig3b <- dinucat %>% mutate(name = fct_relevel(sample, namem_ord)) %>% ggplot(aes(x =name, y =Contribution, fill =Dinucleotides)) + geom_bar(stat="identity", colour = "black")+xlab("")+scale_fill_npg()+theme_light()+theme(axis.text.x=element_text(size = 9, angle = 45, hjust = 0))+theme(axis.text.y=element_text(size = 12))+scale_x_discrete(position = "top")+theme(text = element_text(size = 18))+xlab("tumour type")+ylab("mean number of dinucleotides")
ggsave(plot=fig3b, file='Figure-3b.pdf', device=cairo_pdf, width=11, height=5)

############################################################################
#signature extraction with HDP, 112 samples
dinuchdp <- dinuc2[,which(colSums(dinuc2)>9)]

#"KIDNEY_VINYLIDENE_CHLORIDE_1" 6 
#"LIVER_1.2.3_TRICHLOROPROPANE_1" 4     
#"LIVER_ANTHRAQUINONE_1" 5  
#"LIVER_ANTIMONY_TRIOXIDE_2"   5               
#"LIVER_BROMOCHLOROACETIC_ACID_2" 3     
#"LIVER_COBALT_1" 4                  
#"LIVER_CUMENE_1" 5                   
#"LIVER_DIETHANOLAMINE_1"  4   
#"LIVER_FURAN_1" 3 
#"LIVER_GINKGO_BILOBA_EXTRACT_2" 2     
#"LIVER_INDIUM_PHOSPHIDE_2" 3          
#"LIVER_ISOBUTYL_NITRITE_1" 5         
#"LIVER_NICKEL_OXIDE_1"  2            
#"LIVER_NICKEL_SUBSULFIDE_1" 3         
#"LIVER_NICKEL_SULFATE_1" 5           
#"LIVER_OXAZEPAM_2" 3                  
#"LIVER_DE-71_1" 4
#"LIVER_PRIMACLONE_1" 3 
#"LIVER_SODIUM_TUNGSTATE_DIHYDRATE_1" 6
#"LIVER_SPONTANEOUS_1" 10                             
#"LIVER_VANDIUM_PENTOXIDE_1" 5 
#"LIVER_VINYLIDENE_CHLORIDE_2" 6    
#"LUNG_COBALT_1"    6  
#"LUNG_NICKEL_SULFATE_5"  1            
#"LUNG_SODIUM_TUNGSTATE_DIHYDRATE_6" 1  
#"LUNG_VANADIUM_PENTOXIDE_2"  2        
#"LUNG_VINYLIDENE_CHLORIDE_4" 1        
#"STOMACH_1.2.3_TRICHLOROPROPANE_1"  5

hdp_mut <- hdp_init(ppindex = c(0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), # index of parental node
                     cpindex = c(1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2), # index of the CP to use
                     hh = rep(1, 78), # prior is uniform over 96 categories
                     alphaa = rep(1,times=2), # shape hyperparameters for 2 CPs
                     alphab = rep(1,times=2))  # rate hyperparameters for 2 CPs
hdp_mut <- hdp_addconparam(hdp_mut,
                            alphaa = rep(1,28), # shape hyperparams for 17 CPs
                            alphab = rep(1,28)) # rate hyperparams for 17 CPs

pp <- c(rep(2,6),rep(3,4),rep(4,5),rep(5,5),rep(6,3), 
    rep(7,4),
    rep(8,5),
    rep(9,4),
    rep(10,3),    
    rep(11,2),
    rep(12,3),
    rep(13,5),rep(14,2),
    rep(15,3), 
    rep(16,5),
    rep(17,3),
    rep(18,4), 
    rep(19,3), 
    rep(20,6),rep(21,10), 
    rep(22,5), 
    rep(23,6),rep(24,6), 
    rep(25,1), rep(26,1), 
    rep(27,2), rep(28,1), 
    rep(29,5))

cp <- pp+1                            

hdp_mut <- hdp_adddp(hdp_mut,
                     numdp = 112, # add 112 nodes
                     ppindex=pp, 
                     cpindex=cp)
                   
hdp_mut <- hdp_setdata(hdp_mut,
                       dpindex=30:141, # which nodes to add data to (leaves) 9+193 202
                       data=rbind(t(dinuchdp)))# input data matrix

saveRDS(hdp_mut,file='hdp_dinuc.rds')

chlist <- vector("list", 8)
############################################################################
# This part is computationally intensive.
#for (i in 1:8){
#activated <- dp_activate(hdp_mut,dpindex = 1:numdp(hdp_mut),initcc =20,seed = jobIndex*1000)
#chlist[[i]]<- hdp_posterior(activated,burnin = 200000,n = 100,space = 1000,cpiter = 3,seed = 33333*(jobIndex)+ 1)
#}

#In practice, I run chains in parallel using a FARM job array
# e.g. with a command line like:
##bsub -q long -oo /path/MutDatadinuc-out181.%I.txt -J 'MYJOB[1-8]' -R'select[mem>MEMVAL] rusage[mem=MEMVAL]' -M MEMVAL 'R CMD BATCH --vanilla path/denovodinuc_181.R path/${LSB_JOBINDEX}.denovodinuc.181.Rout'

# where the denovodinuc_181.R is :
#jobIndex <- as.numeric(Sys.getenv("LSB_JOBINDEX"))
#induced <- readRDS('hdp_dinuc.rds')
#activated <- dp_activate(induced,dpindex = 1:numdp(induced),initcc =20,seed = jobIndex*1000)
#chlist<- hdp_posterior(activated,burnin = 200000,n = 100,space = 1000,cpiter = 3,seed = 33333*(jobIndex)+ 1)
#saveRDS(chlist,file=paste('hdp_mut_dinuc_181_denovo.',jobIndex,'.rds',sep=''))

#and the I uploaded the files that I have obtained
#chlist <- vector("list", 8)
############################################################################
chlist[[1]] <- readRDS('hdp_mut_dinuc_181_denovo.1.rds')
chlist[[2]] <- readRDS('hdp_mut_dinuc_181_denovo.2.rds')
chlist[[3]] <- readRDS('hdp_mut_dinuc_181_denovo.3.rds')
chlist[[4]] <- readRDS('hdp_mut_dinuc_181_denovo.4.rds')
chlist[[5]] <- readRDS('hdp_mut_dinuc_181_denovo.5.rds')
chlist[[6]] <- readRDS('hdp_mut_dinuc_181_denovo.6.rds')
chlist[[7]] <- readRDS('hdp_mut_dinuc_181_denovo.7.rds')
chlist[[8]] <- readRDS('hdp_mut_dinuc_181_denovo.8.rds')
mut_multi <- hdp_multi_chain(chlist)
mut_multi <- hdp_extract_components(mut_multi,min.sample=3)

############################################################################
#compare mouse dinucleotide signatures to human dinucleotide signatures...I do not do any normalization because we are working with indels
sign <- mut_multi@comp_categ_distn$mean
#DBS signatures for mouse (normalized)
dbs <- read.delim('../starting_data/DBS_signatures_genome_builds.txt')
dbs1 <- dbs[1:78,2:12]
rownames(dbs1) <- as.matrix(dbs[1:78,1])
csmap <- cos_sim_compare_multiplesignatures(dbs1,t(sign))
cossim <- array(0,4)
bestsig <- matrix('',nrow=4,ncol=1)
namesig <- c('mDBS1','mDBS2','mDBS3','mDBS4')
for (i in seq(1,4)){
bestsig[i] <- colnames(dbs1)[which(csmap[,i]==max(csmap[,i]))];cossim[i]=round(max(csmap[,i]),digit=2)}
d <- data.frame('signature'=namesig,'besthumansignature'=bestsig,'cosine_similarity'=cossim)
write.table(d,file='compare_dinucleotides.txt',quote=F,sep='\t',row.names=F)
#two signatures are new and we can identify mDBS2
namesig <- c('0','mDBS2','mDBS_N1','mDBS_N2')

############################################################################
#Supplementary_Figure_4b
matt <- mut_multi@comp_dp_distn$mean
matt <- matt[30:141,]
#I used SigProfilerPlotting to plot Supplementary_Figure_4a, sigdinuc contains the dinucleotide signatures 
sigdinuc <- t(mut_multi@comp_categ_distn$mean)

rownames(matt) <- colnames(dinuchdp)
colnames(matt) <- namesig
dinucexposureabsolute <- plot_contribution(t(matt*colSums(dinuchdp)), t(mut_multi@comp_categ_distn$mean),mode = "absolute",coord_flip = F,palette=RColorBrewer::brewer.pal(12, "Paired"))+theme(axis.text.x= element_text(angle = 90, hjust = 1, vjust = 0.5,size=7))
ggsave(plot=dinucexposureabsolute, file='Supplementary_Figure_4b.pdf', device=cairo_pdf, width=18, height=6)


############################################################################
#Supplementary_Figure_4c
#mDBS_N2 is confidentially present only in cobalt, in all the cobalt samples 
#the rest is mainly made up of DBS2
matt1 <- mut_multi@comp_dp_distn$mean
for (i in 2:141) {matt1[i,which(mut_multi@comp_dp_distn$cred.int[[i]][1,]==0)]=0}
matt1 <- matt1[30:141,]
rownames(matt1) <- colnames(dinuchdp)
colnames(matt1) <- namesig
which(matt1[,4]>0)
#LUNG_COBALT_METAL_1 LUNG_COBALT_METAL_2 LUNG_COBALT_METAL_3 LUNG_COBALT_METAL_4 LUNG_COBALT_METAL_5 LUNG_COBALT_METAL_6 

#I used SigProfilerPlotting to plot Supplementary_Figure_4c, data in dinucm
#I do not consider the samples without any dinucleotide mutations
dinucn=t(t(dinuc2)/colSums(dinuc2))
ind=which(str_detect(colnames(dinucn), 'LUNG_COBALT_METAL')==1)
dinucm=rowMeans(dinucn[,ind])




############################################################################
#Now I start analyzing indels
indel <- read.delim('../starting_data/SigProfilerMatrixGenerator_matrices/ind.ID83.all',header=T,sep='\t',check.names=F)
indel1 <- indel[,2:182]
rownames(indel1) <- indel[,1]

samples1 <- sort(unique(colnames(indel1)))
samples <- c(samples1[1:36],samples1[42:83],samples1[37:41],samples1[84:181])
indel2 <- indel1[,samples]

############################################################################
#figure 3C
samples <- as.character(lapply(samples,function(x) str_replace_all(x,'_',' ')))

totalindel <- as.data.frame(colSums(indel2))
colnames(totalindel) <- c('number_of_indels')
totalindel <- totalindel %>% mutate(sample=samples)

totalindel <- totalindel %>% mutate(category=as.character(lapply(samples,function(x) str_trim(str_sub(x, start = 1L, end =str_locate(x,"\\d$")[1]-2 )))))
totalindel <- totalindel %>% mutate(tissue=as.character(lapply(samples,function(x) str_sub(x, start = 1L, end =str_locate(x,"\\s")[1] ))))

#with all the points
fig3c <- totalindel %>% mutate(name = fct_relevel(category, namem_ord)) %>% ggplot(aes(x=name, y=number_of_indels,fill=tissue)) + geom_boxplot(outlier.shape = NA)+geom_jitter(color="black", size=0.5,alpha=0.95)+xlab("")+scale_fill_npg()+theme_light()+theme(axis.text.x=element_text(size = 9, angle = 45, hjust = 0))+theme(axis.text.y=element_text(size = 12))+scale_x_discrete(position = "top")+theme(text = element_text(size = 18))+xlab("tumour type")+ylab("number of indels")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+theme(legend.key = element_rect(fill = "white", colour = "white"))
ggsave(plot=fig3c, file='Figure-3c.pdf', device=cairo_pdf, width=11, height=5)

############################################################################
#test the difference in the number of snvs for lung and liver
nindel <- colSums(indel2)

#nindel liver
spoliver <- which((str_detect(names(nindel),'LIVER') & str_detect(names(nindel),'SPONTANEOUS'))==TRUE)
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
chemtestliver <- which((str_detect(names(nindel),'LIVER') & str_detect(names(nindel),chem[i]))==TRUE)
pval[i] <- wilcox.test(nindel[chemtestliver],nindel[spoliver],alternative='greater')$p.value
}
qval <- p.adjust(pval,method='BH')
#0.9967201 0.9967201 0.9967201 0.7252353 0.4342114 0.9967201 0.5174151 0.9967201 0.9967201 0.7252353 0.9967201 0.9967201 0.5174151 0.9967201 0.9967201 0.9967201 0.2684269 0.4342114 0.9967201 0.4342114
qval_liver <- data.frame(chemical=chem,q.value=qval)
write.table(qval_liver,file='liver_indel_qval.txt',quote=F,sep='\t',row.names=F)

#nindel lung
spolung <- which((str_detect(names(nindel),'LUNG') & str_detect(names(nindel),'SPONTANEOUS'))==TRUE)
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
chemtestlung <- which((str_detect(names(nindel),'LUNG') & str_detect(names(nindel),chem[i]))==TRUE)
pval[i] <- wilcox.test(nindel[chemtestlung],nindel[spolung],alternative='greater')$p.value
}
qval <- p.adjust(pval,method='BH')
#0.9950161 0.9835770 0.9950161 0.9835770 0.9835770 0.9950161 0.9950161 0.3192502 0.9950161
qval_lung <- data.frame(chemical=chem,q.value=qval)
write.table(qval_lung,file='lung_indel_qval.txt',quote=F,sep='\t',row.names=F)

#ndinuc liver and lung difference
pp <- which((str_detect(names(nindel),'LUNG'))==TRUE)
pp1 <- which((str_detect(names(nindel),'LIVER'))==TRUE)
wilcox.test(nindel[pp1],nindel[pp])$p.value
#8.004319e-05 121,178 median

############################################################################
#signature extraction with HDP with prior
#KIDNEY_VINYLIDENE_CHLORIDE 1-6       
#LIVER_1,2,3_TRICHLOROPROPANE 1-4     
#LIVER_ANTHRAQUINONE 1-5                          
#LIVER_ANTIMONY_TRIOXIDE 1-6          
#LIVER_BROMOCHLOROACETIC_ACID 1-5    
#LIVER_COBALT_METAL 1-5
#LIVER_CUMENE 1-5                    
#LIVER_DIETHANOLAMINE 1-5             
#LIVER_FURAN 1-5         
#LIVER_GINKGO_BILOBA_EXTRACT 1-3  
#LIVER_INDIUM_PHOSPHIDE 1-5           
#LIVER_ISOBUTYL_NITRITE 1-6              
#LIVER_NICKEL_OXIDE 1-4         
#LIVER_NICKEL_SUBSULFIDE 1-3          
#LIVER_NICKEL_SULFATE_HEXAHYDRATE 1-6 
#LIVER_OXAZEPAM 1-5                  
#LIVER_DE-71 1-5        
#LIVER_PRIMACLONE 1-5                
#LIVER_SODIUM_TUNGSTATE_DIHYDRATE 1-6
#LIVER_SPONTANEOUS 1-11  
#LIVER_VANDIUM_PENTOXIDE 1-6
#LIVER_VINYLIDENE_CHLORIDE 1-7
#LUNG_ANTIMONY_TRIOXIDE 1-5          
#LUNG_COBALT_METAL 1-6                     
#LUNG_ISOBUTYL_NITRITE 1-6
#LUNG_NICKEL_OXIDE 1-5
#LUNG_NICKEL_SUBSULFIDE 1-4
#LUNG_NICKEL_SULFATE_HEXAHYDRATE 1-6
#LUNG_SODIUM_TUNGSTATE_DIHYDRATE 1-6
#LUNG_SPONTANEOUS 1-12                 
#LUNG_VANADIUM_PENTOXIDE 1-3    
#LUNG_VINYLIDENE_CHLORIDE 1-5
#STOMACH_1,2,3_TRICHLOROPROPANE 1-5

cancer_indel <- read.delim('../starting_data/sigProfiler_ID_signatures.txt')
cancer_indel1 <- cancer_indel[,2:18]
rownames(cancer_indel1) <- cancer_indel[,1]

luad_prior <- hdp_prior_init(prior_distn = as.matrix(cancer_indel1), # matrix of prior sigs
                              prior_pseudoc = rep(1000, 17), # pseudocount weights
                              hh=rep(1, 83), # uniform prior over 83 categories
                              alphaa=c(1, 1), # shape hyperparams for 2 CPs
                              alphab=c(1, 1)) # rate hyperparams for 2 CPs

luad_prior <- hdp_addconparam(luad_prior,
                              alphaa = c(1,1), # shape hyperparams for 2 new CPs
                              alphab = c(1,1)) # rate hyperparams for 2 new CPs
                              
luad_prior <- hdp_adddp(luad_prior,numdp = 34,
					ppindex = c(1,rep(19,33)), # index of parental node
                    cpindex = c(3,rep(4,33)) # index of the CP to use
                    )  

luad_prior <- hdp_addconparam(luad_prior,
                            alphaa = rep(1,33), # shape hyperparams for 17 new CPs
                            alphab = rep(1,33)) # rate hyperparams for 17 new CPs

cp <- c(rep(5,6),rep(6,4),rep(7,5),rep(8,6),rep(9,5), 
    rep(10,5),
    rep(11,5),
    rep(12,5),
    rep(13,5),
    rep(14,3),
    rep(15,5),
    rep(16,6),rep(17,4),
    rep(18,3), 
    rep(19,6),
    rep(20,5),
    rep(21,5), 
    rep(22,5), 
    rep(23,6),rep(24,11), 
    rep(25,6), 
    rep(26,7),rep(27,5), 
    rep(28,6), rep(29,6), 
    rep(30,5), rep(31,4), 
    rep(32,6),rep(33,6),rep(34,12),rep(35,3),rep(36,5),rep(37,5))

pp <- cp+15                            

luad_prior <- hdp_adddp(luad_prior,
                     numdp = 181, 
                     ppindex=pp, 
                     cpindex=cp)
      
luad_prior <- hdp_setdata(luad_prior,
                       dpindex=53:233, # which nodes to add data to (leaves)
                       data=rbind(t(indel2))) # input data matrix

saveRDS(luad_prior,file='hdp_indel_prior.rds')

chlist <- vector("list", 8)
############################################################################
# This part is computationally intensive. 
#for (i in 1:8){
#activated <- dp_activate(indel2,dpindex = 19:numdp(indel2),initcc =20,seed = jobIndex*1000)
#chlist[[i]]<- hdp_posterior(activated,burnin = 200000,n = 100,space =1000,cpiter = 3,seed = 33333*i+ 1)
# }

#In practice, I run chains in parallel using a FARM job array
# e.g. with a command line like:
##bsub -q long -oo /path/MutDataindprior-out.%I.txt -J 'MYJOB[1-8]' -R'select[mem>MEMVAL] rusage[mem=MEMVAL]' -M MEMVAL 'R CMD BATCH --vanilla path/prior181_ind.R ${LSB_JOBINDEX}.priordenovo.181.Rout'

# where the prior181_ind.R is :
#jobIndex <- as.numeric(Sys.getenv("LSB_JOBINDEX"))
#induced_prior <- readRDS('hdp_indel_prior.rds')
#activated <- dp_activate(induced_prior,dpindex = 19:numdp(induced_prior),initcc =20,seed = jobIndex*1000)
#chlist <- hdp_posterior(activated,burnin = 200000,n = 100,space = 1000,cpiter = 3,seed = 33333*(jobIndex)+ 1)
#saveRDS(chlist,file=paste('hdp_ind_181_prior.',jobIndex,'.rds',sep=''))

#and the I uploaded the files that I have obtained
#chlist <- vector("list", 8)
############################################################################
chlist[[1]] <- readRDS('hdp_ind_181_prior.1.rds')
chlist[[2]] <- readRDS('hdp_ind_181_prior.2.rds')
chlist[[3]] <- readRDS('hdp_ind_181_prior.3.rds')
chlist[[4]] <- readRDS('hdp_ind_181_prior.4.rds')
chlist[[5]] <- readRDS('hdp_ind_181_prior.5.rds')
chlist[[6]] <- readRDS('hdp_ind_181_prior.6.rds')
chlist[[7]] <- readRDS('hdp_ind_181_prior.7.rds')
chlist[[8]] <- readRDS('hdp_ind_181_prior.8.rds')
mut_multi <- hdp_multi_chain(chlist)
mut_multi <- hdp_extract_components(mut_multi,min.sample=3)

############################################################################
#I used SigProfilerPlotting to plot Supplementary_Figure_5, sign contains the indel signatures 
sign <- t(mut_multi@comp_categ_distn$mean)
rownames(sign) <- indel[,1]
csmap <- cos_sim_compare_multiplesignatures(cancer_indel1,sign)
cossim <- array(0,6)
bestsig <- matrix('',nrow=4,ncol=1)
namesig<-c('mID1','mID2','mID3','mID4','mID5','mID6')
for (i in seq(1,6)){
bestsig[i] <- colnames(dbs1)[which(csmap[,i]==max(csmap[,i]))];cossim[i]=round(max(csmap[,i]),digit=2)}
d <- data.frame('signature'=namesig,'besthumansignature'=bestsig,'cosine_similarity'=cossim)
write.table(d,file='compare_indels.txt',quote=F,sep='\t',row.names=F)
#considering a thereshold of 0.85, I can match them to human signatures
namesig<-c('0','mID1','mID2','mID9','mID3','mID8')


############################################################################
#Figure3D
matt <- mut_multi@comp_dp_distn$mean
matt <- matt[20:52,]
colnames(matt) <- namesig
matt <- as.data.frame(matt)
matt <- matt %>% mutate(sample=namem)
mattr <- gather(matt,'0','mID1','mID2','mID9','mID3','mID8',key='signature',value='contribution')

fig3d <- mattr %>% mutate(name = fct_relevel(sample, namem_ord)) %>% ggplot(aes(x=factor(name),y=contribution,fill=factor(signature))) +geom_bar(stat="identity", colour = "black")+scale_fill_npg()+theme_light()+theme(axis.text.x=element_text(size = 9, angle = 45, hjust = 0))+theme(axis.text.y=element_text(size = 12))+scale_x_discrete(position = "top")+theme(text = element_text(size = 18))+xlab("tumour type")+ylab("contribution")
ggsave(plot=fig3d, file='Figure-3d.pdf', device=cairo_pdf, width=11, height=5)
