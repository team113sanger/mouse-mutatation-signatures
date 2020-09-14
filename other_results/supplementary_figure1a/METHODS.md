We used ENCODE data to look for possible associations with histone markers or open chromatin. We attempted to match sample types as close as possible.

For the open chromatin DNase-seq analysis, only one suitable dataset was found for lung (ENCSR000CNM); containing 3 isogenic replicates. No suitable DNase-seq datasets were found for liver. DNase-seq hotspot data were downloaded for these datasets.

For the histone markers, datasets for 5 markers (H3K27ac, H3K27me3, H3K36me3, H3K79me2, and H3K9ac) were located (ENCSR000CDH, ENCSR000CEN, ENCSR000CEO, ENCSR000CEP, and ENCSR000CEQ respectively) for liver tissue; whilst 2 datasets (H3K4me3 and H3K4me1) were located for lung (ENCSR000CAR and ENCSR000CAQ respectively). 

Once downloaded, all data were processed in R (version 3.6.2) using the MutationalPatterns package (version 1.10.0). Briefly, we first loaded all data into R by converting both our SNV data into a set of genomic ranges, as well as each of the ENCODE datasets. All data used were mapped to the UCSC mm10 reference genome. We used the complete genome as the surveyed region list.

The genomic_distribution test within the MutationalPatterns package was used to determine enrichment or depletion between matching ENCODE data and our SNV list.

After plotting, the differences in Observed/Expected ratios for each marker/tissue pair were assessed using an un-paired Mann–Whitney test. False-discovery rate p-value correction was applied to all p-values.
