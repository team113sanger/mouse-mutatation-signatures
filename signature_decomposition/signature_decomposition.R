# Laura Riva Jan 2019

library(SignatureEstimation)

cos_sim = function(x, y) {
	res = sum(x*y)/(sqrt(sum(x*x))*sqrt(sum(y*y)))
  	res = as.numeric(res)
  	return(res)
  	}


#calculate the cosine similarity between mysignature and cancer_signatures (i.e. COSMIC signatures) without decomposition
cos_sim_comparesignature <- function(signatureinput,cancer_signatures){
	signcompared=matrix(0,1,ncol(cancer_signatures))
	for (i in 1:ncol(cancer_signatures)){
		signcompared[i]=cos_sim(signatureinput,cancer_signatures[,i])
		}
	return(signcompared)
	}


cos_sim_compare_multiplesignatures <- function(signatures1,signatures2){
	signcompared=matrix(nrow=dim(signatures1)[2],ncol=dim(signatures2)[2])
	for (i in 1:dim(signatures1)[2]){
		for (j in 1:dim(signatures2)[2]){	
			signcompared[i,j]=cos_sim(signatures1[,i],signatures2[,j])
			}
		}	
	return(signcompared)
	}


#calculate the cosine similarity between the recostruction of one of my signature (signatureinput) and known input signatures (cancer_signatures)
cos_sim_reconstructed = function(signatureinput,cancer_signatures){
	if (is.null(dim(cancer_signatures))) {
		return(cos_sim(as.numeric(signatureinput),as.numeric(cancer_signatures)))
		}
  	else {QP = findSigExposures(as.numeric(signatureinput),cancer_signatures)$exposures
		reconstructed <- as.matrix(cancer_signatures) %*% as.numeric(QP)
  		return(cos_sim(as.numeric(reconstructed),as.numeric(signatureinput)))
  		}
  	}  	


#calculate the minimum number of known signatures (cancer_signatures) that can reconstruct my signature (signatureinput)
#cancer_signatures should have at least 9 signatures. If a starting signature is made up of more than 3 known signatures, this signature is considered a novel signature. 
#3 is the default value but if can be increased to 4 or 5 but I think that there is no reason to have more than 5 signatures as the majority (84%) of human tumour samples are made up of max of 5 signatures   
min_best_cos_sim_signature=function(signatureinput,cancer_signatures,cs_threshold=0.9,num_sig_threshold=3){
	signatureinput=as.numeric(signatureinput)
	cancer_signatures=as.matrix(cancer_signatures)	
	signames=colnames(cancer_signatures)
	bestonesignature=cos_sim_comparesignature(signatureinput,cancer_signatures)
	compare1=which(bestonesignature>=cs_threshold)
	if (length(compare1)>0) {
		res=data.frame('weight'=1)
		rownames(res)=signames[which(bestonesignature==max(bestonesignature))]
		return(list(res,data.frame('cos_sim'=max(bestonesignature))))
	}
	else {
		allsolutions=matrix(0,nrow=length(colnames(cancer_signatures))-1,ncol=dim(combn(length(colnames(cancer_signatures)),num_sig_threshold))[2])
		for (i in 2:num_sig_threshold){
			regions=combn(c(1:length(colnames(cancer_signatures))),i)
 			for (j in 1:dim(regions)[2]){
 				sel=cancer_signatures[,regions[,j]]
 				allsolutions[i,j]=cos_sim_reconstructed(signatureinput,sel) 
 				}
 			if (length(which(allsolutions[i,]>=cs_threshold))>0){
 				break}
 		}
 		if (length(which(allsolutions[i,]>=cs_threshold))>0){
 			position=which(allsolutions[i,]>=cs_threshold)
			ind=allsolutions[i,which(allsolutions[i,]>=cs_threshold)]
			sigs=signames[combn(c(1:length(colnames(cancer_signatures))),i)[,position[which(ind==max(ind)[1])]]]
			res=findSigExposures(as.numeric(signatureinput),cancer_signatures[,sigs])$exposures
			colnames(res)=c('weight')
			return(list(as.data.frame(res),data.frame('cos_sim'=allsolutions[i,position[which(ind==max(ind)[1])]])))
			}
		else {
			res=NaN
			return(res)
			}
		}
	}