library(tidyverse)
library(BSgenome.Mmusculus.UCSC.mm10)
library(MutationalPatterns)
options(stringsAsFactors = F)

uni=readRDS('../starting_data/snvs.rds')
samples=sort(unique(uni$sample))
ref_genome = "BSgenome.Mmusculus.UCSC.mm10"

vargs=GRangesList() 

for (i in 1:length(samples)){
uni1= uni %>% filter(sample == samples[i])  
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

repli_strand_granges = GRanges(seqnames = repli_strand$chrom,ranges = IRanges(start = repli_strand$star + 1,end = repli_strand$end),strand_info=factor(repli_strand$versus))
mut_mat_s_rep <- mut_matrix_stranded(vargs, ref_genome, repli_strand_granges, mode = "replication")
colnames(mut_mat_s_rep)=samples

strand_counts_rep <- strand_occurrences(mut_mat_s_rep,by=samples)
strand_bias_rep <- strand_bias_test(strand_counts_rep)

#this is the table with all the significant values like the one for transcriptional strand bias, should be corrected for false discovery rate
strand_bias_rep=strand_bias_rep %>% mutate(q_value=p.adjust(p_poisson,method='BH'))
strand_bias_rep=strand_bias_rep %>% select(c(group,type,left,right,total,ratio,p_poisson,q_value))
write.table(strand_bias_rep,'Supplementary_table_5.txt',quote=F,sep='\t',col.names=T,row.names=F)





