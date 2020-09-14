# This R script runs the SNV correlation analysis on Laura's data
# Alastair Droop, 2020-04-09

# Get the necessary libraries:
library(MutationalPatterns)
library(BSgenome.Mmusculus.UCSC.mm10.masked)
library(parallel)
library(matrixStats) # For the fast "count" function

# Set the run IO:
io <- list(
    'rdata.file' = file.path('..', 'rdata', 'analysis-data.Rdata'),
    'output.dir' = file.path('..', 'output')
)

# Set the run parameters:
params <- list(
    'genome' = 'BSgenome.Mmusculus.UCSC.mm10',
    'n.bootstraps' = 250,
    'cores' = 8,
    'threshold.p' = 0.05
)

# Load the appropriate genome library:
message(sprintf('using genome "%s"...', params$genome))
library(package=params$genome, character.only=TRUE)
ref.genome <- base::get(params$genome)

# Load the input data:
message(sprintf('reading input data from "%s"...', io$rdata.file))
load(io$rdata.file)

# A function to generate a random set of SNVs:
sampleSNVs <- function(s){
    chrs <- seqlevelsInUse(s)
    chr.lens <- seqlengths(seqinfo(s)[chrs])
    chr.n <- as.vector(table(seqnames(s))[chrs])
    r <- c(mapply(sample.int, chr.lens, chr.n), recursive=TRUE)
    chr.d <- GRanges(
        seqnames = Rle(chrs, chr.n),
        ranges = IRanges(start=r, width=1),
        strand = '*',
        seqinfo = seqinfo(s)
    )
    return(chr.d)    
}

# A function to count the number of SNVs that hit a region:
overlapSNVs <- function(s, regions){
    return(count(countOverlaps(s, regions) > 0))
}

# A function to bootstrap random SNV range interactions:
bootstrapSNVs <- function(s, regions, k=1){
    res <- sapply(1:k, function(i){
        return(overlapSNVs(sampleSNVs(s), regions))
    })
    return(res)
}

# Build the possible treatment/region combinations:
rt.combn <- expand.grid(region=names(encode.list), treatment=names(snv.list))
rt.n <- nrow(rt.combn)
message(sprintf('%d (region, treatment) pairs to run', nrow(rt.combn)))

# Bootstrap the possible treatment/region combinations:
treat.test <- do.call(rbind, mclapply(1:rt.n, function(i){
    # Pull out the required region and treatment names:
    region <- rt.combn[i, 'region']
    treatment <- rt.combn[i, 'treatment']
    message(sprintf('[%3d/%3d] running (%s, %s)...', i, rt.n, region, treatment))

    # Pull out the SNV and region data:
    regions <- encode.list[[region]]
    snvs <- snv.list[[treatment]]

    # Save the SNV count:
    n <- length(snvs)

    # Calculate real hits:
    observed.hits <- overlapSNVs(snvs, regions)
    observed.prob <- observed.hits / n

    # Calculate simulated hits:
    expected.hits <- mean(bootstrapSNVs(snvs, regions, k=params$n.bootstraps))
    expected.prob <- expected.hits / n

    # Perform a Bionomial test to see if the observed and expected have the same probability:
    observed.test <- binom.test(x=observed.hits, n=n, p=expected.prob)

    # Build the output data.frame:
    output <- data.frame(
        'treatment' = treatment,
        'region' = region,
        'n.boot' = params$n.bootstraps,
        'n' = n,
        'obs.hits' = observed.hits,
        'exp.hits' = expected.hits,
        'obs.p' = observed.prob,
        'exp.p' = expected.prob,
        'p' = observed.test$p.value
    )

    # Return the output:
    return(output)
}, mc.cores=params$cores))

# Apply p-value adjustment:
treat.test$p.adjust <- p.adjust(treat.test$p, method='fdr')

# Apply significance marker:
treat.test$sig <- c('TRUE'='*', 'FALSE'='')[as.character(treat.test$p.adjust <= params$threshold.p)]

# Sort by p-value:
treat.test <- treat.test[order(treat.test$p), ]

# Save the output data:
output.filename <- file.path(io$output.dir, 'binomial-test.txt')
message(sprintf('saving test data to "%s"...', output.filename))
write.table(treat.test, file=output.filename, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
