# Script to produce figure S1e

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure1Se.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("openxlsx")
  library("GenomicRanges")
  library("annotatr")
  library("ChIPseeker")
  library("BSgenome.Mmusculus.UCSC.mm10")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene") 
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--cre"), type="character", default=NULL, 
              help="cCRE information", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$cre)){
  print_help(opt_parser)
  stop("cCRE information (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load Genome\n")
genome <- GRanges(seqnames=seqnames(BSgenome.Mmusculus.UCSC.mm10), 
  IRanges(start=1, end=seqlengths(BSgenome.Mmusculus.UCSC.mm10)))
genome <- genome[-as.numeric(grep('_',seqnames(genome)))]
genome <- genome[-as.numeric(grep('chrM',seqnames(genome)))]
genome <- unlist(tile(genome, width=1000))

par(mfrow = c(2, 2))

genome$region <- paste(genome)

cat("Annotation to TSS\n")
peakAnno <- annotatePeak(genome, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion=c(-3000, 3000), annoDb="org.Mm.eg.db")

print(paste('Percentage annotated:', round(mean(1:length(genome) %in% subjectHits(findOverlaps(as.GRanges(peakAnno), genome)))*100,2), '%'))

svg(paste0(opt$outdir,"/Fig1e.svg"))
plotAnnoPie(peakAnno, ndigit = 1, main = "Genome")
  
cat("Annotation to ENCODE cCRE\n")
annotated <- readRDS(opt$cre)

print(paste('Percentage annotated:', round(mean(1:length(genome) %in% subjectHits(findOverlaps(annotated, genome)))*100,2), '%'))
genome <- annotated[unique(queryHits(findOverlaps(annotated, genome)))]

slices <- round(table(genome$element) / length(genome)*100,2)
lbls <-  paste0(names(slices),": ", round(slices,1), "%")

pie(slices,labels = lbls, col=rainbow(length(lbls), v=0.75), cex = 0.8, main = "Genome")

dev.off()

sessionInfo()