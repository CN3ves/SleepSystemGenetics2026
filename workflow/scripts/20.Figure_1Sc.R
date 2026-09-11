# Script to produce figure S1c

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure1Sc.log', open = "wt")
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

genome$region <- paste(genome)

cat("Annotation to TSS\n")
peakAnno <- annotatePeak(genome, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion=c(-3000, 3000), annoDb="org.Mm.eg.db")

print(paste('Percentage annotated:', round(mean(1:length(genome) %in% subjectHits(findOverlaps(as.GRanges(peakAnno), genome)))*100,2), '%'))

svg(paste0(opt$outdir,"/FigS1c1.svg"))
plotAnnoPie(peakAnno, ndigit = 1, main = "Genome")
dev.off()

cat("Annotation to ENCODE cCRE\n")
annotated <- readRDS(opt$cre)

print(paste('Percentage annotated:', round(mean(1:length(genome) %in% subjectHits(findOverlaps(annotated, genome)))*100,2), '%'))
genome <- annotated[unique(queryHits(findOverlaps(annotated, genome)))]

slices <- round(table(genome$element) / length(genome)*100,2)
lbls <-  paste0(names(slices),": ", round(slices,1), "%")

svg(paste0(opt$outdir,"/FigS1c2.svg"))
par(mar=c(5.1, 4.1, 4.1, 10.5))
pie(slices,labels = lbls, col=rainbow(length(lbls), v=0.75), cex = 0.8, main = "Genome")

dev.off()

# merge inputs for table
ov <- findOverlaps(as.GRanges(peakAnno),genome)

mcols(genome)$region <- paste(genome)
cre <- as.data.frame(genome[subjectHits(ov)])[,c('region', 'element')]

peakAnno <- as.GRanges(peakAnno)
mcols(peakAnno)$region <- paste(peakAnno)
tss <- as.data.frame(peakAnno[queryHits(ov)])[,c('region', 'annotation')]

names(cre) <- paste0('region:ENCODE_',names(cre))
names(tss) <- paste0('region:TSS_',names(tss))

tab <- cbind(tss, cre)

write.csv(tab, paste0(opt$outdir, "/data/S1c.csv"))

sessionInfo()