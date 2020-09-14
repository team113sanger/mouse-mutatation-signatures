# This script plots the histone marker data to show if there is a difference between spontaneous & non-spontaneous
# Alastair Droop, 2020-04-16

# Set the run IO:
io <- list(
    'binomial.input' = file.path('..', 'output', 'binomial-test.txt'),
    'mp.input' = file.path('..', 'output', 'mutation-patterns-dist.txt'),
    'output.pdf' = file.path('..', 'output', 'analysis-results.pdf'),
    'output.tests' = file.path('..', 'output', 'analysis-results.csv'),
    'output.counts' = file.path('..', 'output', 'analysis-counts.csv')    
)

# Set the run parameters:
params <- list(
)

# Load the necessary libraries:
library(ggplot2)

# Load the input data:
binom.data <- read.table(io$binomial.input, sep='\t', header=TRUE)
mp.data <- read.table(io$mp.input, sep='\t', header=TRUE)

# Merge test data:
binom.data <- cbind.data.frame('test'='Binomial', binom.data[,c('treatment', 'region', 'n', 'obs.hits', 'exp.hits')])
mp.data <- cbind.data.frame('test'='MutPatterns', mp.data[, c('sample', 'region', 'n_muts', 'observed', 'expected')])
colnames(binom.data) <- c('test', 'treatment', 'region.label', 'n', 'observed', 'expected')
colnames(mp.data) <- c('test', 'treatment', 'region.label', 'n', 'observed', 'expected')
test.data <- rbind.data.frame(binom.data, mp.data)
test.data$treatment.tissue <- tolower(gsub('^([^_]+)_(.*)_(\\d+)$', '\\1', test.data$treatment, perl=TRUE))
test.data$treatment.chemical <- gsub('^([^_]+)_(.*)_(\\d+)$', '\\2', test.data$treatment, perl=TRUE)
test.data$treatment.rep <- as.numeric(gsub('^([^_]+)_(.*)_(\\d+)$', '\\3', test.data$treatment, perl=TRUE))
test.data$type <- gsub('^([^-]+)-([^-]+)-([^-]+)-(\\d+)$', '\\1', test.data$region.label, perl=TRUE)
test.data$marker <- gsub('^([^-]+)-([^-]+)-([^-]+)-(\\d+)$', '\\2', test.data$region.label, perl=TRUE)
test.data$region.tissue <- gsub('^([^-]+)-([^-]+)-([^-]+)-(\\d+)$', '\\3', test.data$region.label, perl=TRUE)
test.data$region.rep <- as.numeric(gsub('^([^-]+)-([^-]+)-([^-]+)-(\\d+)$', '\\4', test.data$region.label, perl=TRUE))
for(i in c('test', 'treatment', 'region.label', 'treatment.tissue', 'treatment.chemical', 'type', 'marker', 'region.tissue')) {test.data[[i]] <- factor(test.data[[i]])}
test.data$induction <- factor(c('Induced', 'Spontaneous')[as.numeric(test.data$treatment.chemical == 'SPONTANEOUS')+1])
test.data <- test.data[, c('test', 'induction', 'treatment.tissue', 'treatment.chemical', 'treatment.rep', 'type', 'marker', 'region.tissue', 'region.rep', 'n', 'observed', 'expected')]

# # Only look at lung & liver treatments:
# test.data <- test.data[test.data$treatment.tissue %in% c('LIVER', 'LUNG'), ]
# test.data$treatment.tissue <- droplevels(test.data$treatment.tissue)

# Calculate the o/e ratios:
test.data$ratio <- test.data$observed / test.data$expected

# Only look at data where the treatment tissue matches the region tissue:
sel.data <- test.data[as.character(test.data$treatment.tissue) == as.character(test.data$region.tissue), ]

# Only look at The MutationPatterns data:
sel.data <- sel.data[sel.data$test == 'MutPatterns', ]

# Plot the data:
# pdf(io$output.pdf, width=18, height=3)
p <- ggplot(sel.data, aes(y=ratio))
p <- p + geom_hline(yintercept=1, colour='grey75', linetype='dotted')
p <- p + geom_boxplot(aes(colour=induction), outlier.shape=1)
p <- p + facet_grid(~marker+treatment.tissue)
p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
p <- p + theme(panel.background=element_blank(), panel.border=element_rect(fill=NA))
p <- p + labs(y=expression(frac(Observed, Expected)))
p <- p + ggtitle('SNV Ratios by Induction Status and Marker')
ggsave(plot=p, file=io$output.pdf, device=cairo_pdf, width=18, height=3)

# Extract the specific comparisons:
comparisons <- data.frame(
    'marker' = c('DNase', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K79me2', 'H3K9ac'),
    'tissue' = c('lung', 'liver', 'liver', 'liver', 'lung', 'lung', 'liver', 'liver')
)

wcTest <- function(i){
    marker <- i[['marker']]
    tissue <- i[['tissue']]
    sel.d <- sel.data[sel.data$marker == marker & sel.data$treatment.tissue == tissue, ]
    x <- sel.d[sel.d$induction == 'Induced', 'ratio']
    y <- sel.d[sel.d$induction == 'Spontaneous', 'ratio']
    test.res <- wilcox.test(x=x, y=y, paired=FALSE, exact=TRUE)
    return(test.res$p.value)
}

comparisons$p.value <- apply(comparisons, 1, wcTest)
comparisons$p.adjust <- p.adjust(comparisons$p.value, method='fdr')
write.csv(comparisons, file=io$output.tests, row.names=FALSE)

# Extract the sample counts from sel.data:
sel.counts <- data.frame(
    'marker' = factor(rep(c('DNase', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K79me2', 'H3K9ac'), each=2), levels=levels(sel.data$marker)),
    'tissue' = factor(rep(c('lung', 'liver', 'liver', 'liver', 'lung', 'lung', 'liver', 'liver'), each=2), levels=levels(sel.data$treatment.tissue)),
    'induction' = factor(rep(c('Induced', 'Spontaneous'), times=8), levels=levels(sel.data$induction))
)

group.data <- sapply(1:nrow(sel.counts), function(i){sel.data[(sel.data$marker==sel.counts[i, 'marker']) & (sel.data$treatment.tissue == sel.counts[i, 'tissue']) & (sel.data$induction == sel.counts[i, 'induction']),]}, simplify=FALSE)
sel.counts$n <- sapply(group.data, nrow)
sel.counts$Q1 <- sapply(group.data, function(i){quantile(i$observed / i$expected, probs=0.25)})
sel.counts$median <- sapply(group.data, function(i){median(i$observed / i$expected)})
sel.counts$Q3 <- sapply(group.data, function(i){quantile(i$observed / i$expected, probs=0.75)})
write.csv(sel.counts, file=io$output.counts, row.names=FALSE)
