# Script to produce figure 3a

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure3a.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("org.Mm.eg.db")
  library("ChIPseeker")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene") 
  library("GenomicRanges")
  library("openxlsx")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--S6"), type="character", default=NULL, 
              help="Table S6", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$S6)){
  print_help(opt_parser)
  stop("Table S6 (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load Table S6\n")
par(mfrow = c(2, 1))
QTL <- read.xlsx(opt$S6, "ATAC_QTL")

cat("TSS annotation\n")
QTL$seqnames <- gsub(":.*","",QTL$region) 
QTL$start <- as.numeric(gsub(".*:(.*)-.*","\\1",QTL$region))
QTL$end <- as.numeric(gsub(".*-","",QTL$region) )

QTL <- GRanges(QTL[,c("seqnames","start", "end", "lod")])
  
#Peak annotation
peakAnno <- annotatePeak(unique(QTL[,-1]), TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion=c(-3000, 3000), annoDb="org.Mm.eg.db")

svg(paste0(opt$outdir,"/Fig3a.svg"))
plotAnnoPie(peakAnno, ndigit = 1, main='up')

cat("ENCODE annotation\n")
cre <- read.delim("encode_files/mm10-cCREs.bed", header= FALSE) 
names(cre) <- c("seqnames", "start", "end", "something", "else","element")
cre <- GRanges(cre$seqnames, IRanges(start=cre$start, end=cre$end), mcols=cre[,"element", drop = FALSE])

original_ranges <- as.GRanges(peakAnno)
annotated = annotate_regions(
  regions = original_ranges,
  annotations = cre,
  ignore.strand = TRUE,
  quiet = FALSE)

mcols(annotated)$element <- annotated$annot$mcols.element
mcols(annotated) <- mcols(annotated)[, -which(names(mcols(annotated))=="annot")]

annotated <- GRanges(annotated[,1:5],mcols=annotated[,-1:-5])
names(mcols(annotated)) <- gsub("mcols.", "", names(mcols(annotated)))

# Pie Chart with Percentages
slices <- round(table(annotated$element) / length(annotated)*100,1)
lbls <-  paste0(names(slices),": ", slices, "%")

par(mar=c(0.5, 10, 0.5, 10))
pie(slices,labels = lbls, col=rainbow(length(lbls), v=0.75), cex = 0.8) 

cat("Save plot\n")
dev.off()

sessionInfo()