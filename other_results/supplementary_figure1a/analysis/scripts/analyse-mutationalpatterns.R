# This R script runs the SNV correlation analysis on Laura's data
# Alastair Droop, 2020-04-09

# Get the necessary libraries:
library(MutationalPatterns)
library(BSgenome.Mmusculus.UCSC.mm10.masked)
library(parallel)

# Set the run IO:
io <- list(
    'rdata.file' = file.path('..', 'rdata', 'analysis-data.Rdata'),
    'output.dir' = file.path('..', 'output')
)

# Set the run parameters:
params <- list(
    'genome' = 'BSgenome.Mmusculus.UCSC.mm10'
)

# Load the appropriate genome library:
message(sprintf('using genome "%s"...', params$genome))
library(package=params$genome, character.only=TRUE)
ref.genome <- base::get(params$genome)

# Load the input data:
message(sprintf('reading input data from "%s"...', io$rdata.file))
load(io$rdata.file)

# Calculate the mutational distribution:
message('calculating genomic distributions...')
snv.dist <- genomic_distribution(snv.list, rep(list(ref.coverage), length(snv.list)), encode.list)
snv.dist$sample.type <- gsub('_\\d+', '', snv.dist$sample, perl=TRUE)

# Test the mutational distributions:
message('testing genomic distributions...')
snv.test <- enrichment_depletion_test(snv.dist, by=snv.dist$sample.type)

# Save the mutation distribution data:
snv.dist.filename <- file.path(io$output.dir, 'mutation-patterns-dist.txt')
message(sprintf('writing mutation distribution table to "%s"...', snv.dist.filename))
write.table(snv.dist, file=snv.dist.filename, col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

# Save the mutation distribution test results:
snv.test.filename <- file.path(io$output.dir, 'mutation-patterns-test.txt')
message(sprintf('writing mutation test table to "%s"...', snv.test.filename))
write.table(snv.test, file=snv.test.filename, col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)
