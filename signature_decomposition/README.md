<h1>MinSig Minimin Number of Signatures</h1>

<h2>signature_decomposition.R</h2>


signature_decomposition.R can be used to decompose a de novo extracted signature in known signatures.    
A signature is defined as a “new signature” if it cannot be decomposed in a maximum of 3 input signatures.    
It finds the minimum number of known signatures whose reconstructed signature have  
a cosine similarity >=0.9 with the input signature (0.9 is the default threshold but it can be set by the user).    
old and new (PCAWG) COSMIC signatures can be used 

**TO RUN signature_decomposition**
```
source("signature_decomposition.R")
```

Dependencies are **SignatureEstimation** and **tidyverse** (for  sample_decomposition.R)

you can find **SignatureEstimation** here:  
https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/index.cgi#signatureestimation  
https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/software/signatureestimation/SignatureEstimation.tar.gz  
paper https://www.ncbi.nlm.nih.gov/pubmed/29028923)  
