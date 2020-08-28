#Laura Riva April 2020

library(tidyverse)
source('../signature_decomposition/signature_decomposition.R')
cancer_signatures60 = read.table('../starting_data/PCAWG_signatures.txt', sep = "\t", header = TRUE)
cancer_signatures60 = as.matrix(cancer_signatures60[,2:66])
nm=colnames(cancer_signatures60)
nhdp=readRDS('mousesignatures_norm.rds')

names=c("A[C>A]A", "A[C>A]C",
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

humandata1=read.csv('allWGS.96.csv',check.names=FALSE)
humandata=humandata1[,3:4647]
rownames(humandata)=names
active_signatures=readRDS('active_signaturesall.rds')

activeWES1=read.csv('nonPCAWG_WES_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv',check.names=FALSE)
activeWES2=read.csv('TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv',check.names=FALSE)
humanWES1=read.csv('WES_Other.96.csv',check.names=FALSE)
humanWES2=read.csv('WES_TCGA.96.csv',check.names=FALSE)

humanWES=cbind(humanWES1[,3:9693],humanWES2[,3:9495])
rownames(humanWES)=names
namessample=sort(colnames(humanWES))
humanWES=humanWES[,namessample]

activeWES=rbind(activeWES1,activeWES2)
colnames(activeWES)=c('Type','SampleName',colnames(activeWES)[3:68])
activeWES=activeWES %>% mutate(name=paste(Type,SampleName,sep='::'))
rownames(activeWES)=activeWES$name
activeWES=activeWES[namessample,]
activenormwes=activeWES[,4:68]/rowSums(activeWES[,4:68])
rownames(activenormwes)=colnames(humanWES)

sigtocons=c('SBS1','SBS2','SBS3','SBS4','SBS5','SBS6','SBS7a','SBS7b','SBS7c','SBS7d','SBS8','SBS9','SBS10a','SBS10b','SBS11','SBS12','SBS13','SBS14','SBS15','SBS16','SBS17a','SBS17b','SBS18','SBS19','SBS20','SBS21','SBS22','SBS23','SBS24','SBS25','SBS26','SBS28','SBS29','SBS30','SBS31','SBS32','SBS33','SBS34','SBS35','SBS36','SBS37','SBS38','SBS39','SBS40','SBS41','SBS42','SBS44')

sigtest=cbind(cancer_signatures60[,sigtocons],nhdp[,4],nhdp[,7])
colnames(sigtest)=c(sigtocons,'SBS_N1','SBS_N2')

############################################################################
##this part is computationally intensive and was done on the farm, dividing the original WES and WGS spectra in parts
##WGS
#wgs=matrix(0,4645,ncol=49)
#colnames(wgs)=colnames(sigtest)
#for (i in seq(1,4645)){
#	boot = bootstrapSigExposures(humandata[,i]/sum(humandata[,i]), sigtest, R = 10000, mutation.count = sum(humandata[,i]))
#	signature_fractionhdpT=apply(as.matrix(boot$exposures),1,function(x) 1-(1+length(which(x>0.05)))/10001)
#	wgs[i,]=signature_fractionhdpT
#}
#saveRDS(wgs,file='WGS_0.05.rds')
##WES
#wes=matrix(0,nrow=19184,ncol=49)
#colnames(wes)=colnames(sigtest)
#for (i in seq(1,19184)){
#	print(i)
#	boot = bootstrapSigExposures(humanWES[,i]/sum(humanWES[,i]), sigtest, R = 10000, mutation.count = sum(humanWES[,i]))
#	signature_fractionhdpT=apply(as.matrix(boot$exposures),1,function(x) 1-(1+length(which(x>0.05)))/10001)
#	wes[i,]=signature_fractionhdpT
#}
#saveRDS(wes,file='WES_0.05.rds')
############################################################################

wgs_0.05=readRDS('wgs_0.05.rds')
wes_0.05=readRDS('wes_0.05.rds')

tumwgs_19_0.05=which(wgs_0.05[,'SBS19']==0 & colSums(humandata)>=200)
tumwgs_42_0.05=which(wgs_0.05[,'SBS42']==0 & colSums(humandata)>=200)
tumwgs_N1_0.05=which(wgs_0.05[,'SBS_N1']==0 & colSums(humandata)>=200)
tumwgs_N2_0.05=which(wgs_0.05[,'SBS_N2']==0 & colSums(humandata)>=200)

tumwes_19_0.05=which(wes_0.05[,'SBS19']==0 & colSums(humanWES)>=200)
tumwes_42_0.05=which(wes_0.05[,'SBS42']==0 & colSums(humanWES)>=200)
tumwes_N1_0.05=which(wes_0.05[,'SBS_N1']==0 & colSums(humanWES)>=200)
tumwes_N2_0.05=which(wes_0.05[,'SBS_N2']==0 & colSums(humanWES)>=200)

#write table 
#write samples











