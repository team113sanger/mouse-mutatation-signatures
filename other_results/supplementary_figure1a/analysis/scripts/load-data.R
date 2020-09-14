# This R script loads the SNV correlation analysis data into R for Laura's data
# Alastair Droop, 2020-04-17

# Get the necessary libraries:
library(MutationalPatterns)
library(BSgenome.Mmusculus.UCSC.mm10.masked)
library(parallel)

# Set the run IO:
io <- list(
    'input.snv.file' = file.path('..', '..', 'input', 'snv', 'snvs_181tumours.txt'),
    'input.encode.dir' = c(
        'histone' = file.path('..', '..', 'input', 'encode-histone', 'bed'),
        'chromatin' = file.path('..', '..', 'input', 'encode-chromatin', 'bed')
    ),
    'output.dir' = file.path('..', 'rdata')
)

# Set the run parameters:
params <- list(
    'snv.datafilename' = 'snvs_181tumours.data',
    'output.filename' = 'analysis-data.RData',
    'genome' = 'BSgenome.Mmusculus.UCSC.mm10',
    'cores' = 8,
    'encode.cols' = list(
        'histone' = c('chrom', 'start', 'end', 'name', 'score', 'strand', 'signal', 'p.value', 'q.value', 'peak'),
        'chromatin' = c('chrom', 'start', 'end', 'name', 'score', 'strand', 'signal', 'p.value', 'q.value')
    )
)

# Load the appropriate genome library:
message(sprintf('using genome "%s"...', params$genome))
library(package=params$genome, character.only=TRUE)
ref.genome <- base::get(params$genome)

# A function to load a specific BED file:
loadBED <- function(filename, cols, type, tissue, marker, rep=1){
    message(sprintf('reading %s BED file from "%s"...', type, filename))
    d <- read.table(filename, sep='\t', header=FALSE, col.names=cols)
    rownames(d) <- d$name
    res <- makeGRangesFromDataFrame(d, seqinfo=seqinfo(ref.genome), keep.extra.columns=FALSE)
    res$type <- type
    res$tissue <- tissue
    res$marker <- marker
    res$replicate <- rep
    return(res)
}

# Build a list of ENCODE input data:
message('identifying ENCODE BED files...')
encode.input <- do.call(rbind, lapply(names(io$input.encode.dir), function(type){
    input.dir <- io$input.encode.dir[[type]]
    input.files <- list.files(input.dir, pattern='*\\.bed$')
    input.data <- data.frame(
        'type' = type,
        'marker' = gsub('^([^-_]+)(:?_r(\\d))?-([^.]+)\\.bed$', '\\4', input.files, perl=TRUE),
        'tissue' = gsub('^([^-_]+)(:?_r(\\d))?-([^.]+)\\.bed$', '\\1', input.files, perl=TRUE),
        'rep' = as.numeric(gsub('^([^-_]+)(:?_r(\\d))?-([^.]+)\\.bed$', '\\3', input.files, perl=TRUE)),
        'filename' = file.path(input.dir, input.files)
    )
    input.data$rep[is.na(input.data$rep)] <- 1
    input.data$type <- factor(input.data$type)
    input.data$marker <- factor(input.data$marker)
    input.data$tissue <- factor(input.data$tissue)
    return(input.data)
}))
rownames(encode.input) <- sprintf('%s-%s-%s-%d', encode.input$type, encode.input$marker, encode.input$tissue, encode.input$rep)
message(sprintf('%d ENCODE BED files located', nrow(encode.input)))

# Load the ENCODE data into a GRangesList object:
encode.list <- lapply(1:nrow(encode.input), function(i){
    input.type <- encode.input[i, 'type']
    return(loadBED(
        filename = encode.input[i, 'filename'],
        cols = params$encode.cols[[input.type]],
        type = input.type,
        tissue = encode.input[i, 'tissue'],
        marker = encode.input[i, 'marker'],
        rep = encode.input[i, 'rep']
    ))
})
names(encode.list) <- rownames(encode.input)
encode.list <- GRangesList(encode.list)

# Load the SNV data:
snv.datafile <- file.path(io$output.dir, params$snv.datafilename)
if(file.exists(snv.datafile)){
    #Read from file if it is present:
    message(sprintf('reading reformatted SNV list from "%s"...', snv.datafile))
    snv.list <- readRDS(snv.datafile)
} else {
    # Load the input SNV data:
    message(sprintf('reading SNVs from "%s"...', io$input.snv.file))
    snvs <- read.table(io$input.snv.file, sep='\t', header=TRUE, colClasses=c(rep('character', 2), 'numeric', rep('character', 2)))
    snvs$sample <- factor(snvs$sample)
    snvs$chrom <- factor(snvs$chrom, levels=sprintf('chr%s', c(1:19, 'X', 'Y')), ordered=TRUE)

    # Convert SNV data into a GRangesList:
    message('reformatting SNVs...')
    snv.list <- GRangesList(mclapply(levels(snvs$sample), function(s){
        message(sprintf('  reformatting %s...', s))
        s.data <- snvs[snvs$sample==s, ]
        ref <- DNAStringSet(s.data$ref)
        alt <- DNAStringSetList(as.list(DNAStringSet(s.data$alt))) # <--- This line needs serious optimisation!
        output <- GRanges(
            seqnames = s.data$chrom,
            ranges = IRanges(start=s.data$pos, width=1),
            strand = '*',
            seqinfo = seqinfo(ref.genome),
            REF = ref,
            ALT = alt
        )
        names(output) <- sprintf('%s:%d_%s/%s', gsub('chr', '', seqnames(output)), start(output), s.data$ref, s.data$alt)
        return(output)
    }, mc.cores=params$cores))
    names(snv.list) <- levels(snvs$sample)

    # Save the SNV list data to file:
    message(sprintf('saving reformatted SNVs to "%s"...', io$input.snv.listfile))
    saveRDS(snv.datafile, io$input.snv.listfile)
}

# Build the coverage object:
message('building reference coverage...')
ref.coverage <- GRanges(seqinfo(ref.genome))
ref.coverage <- ref.coverage[seqlevels(ref.coverage) %in% seqlevelsInUse(snv.list)]

# Save the loaded data:
output.file <- file.path(io$output.dir, params$output.filename)
message(sprintf('saving processed data to "%s"...', output.file))
save(snv.list, encode.list, ref.coverage, encode.input, file=output.file)
