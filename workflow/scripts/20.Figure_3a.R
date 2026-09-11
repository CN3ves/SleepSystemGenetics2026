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
  make_option(c("-b", "--cre"), type="character", default=NULL, 
              help="cCRE information", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$S6)){
  print_help(opt_parser)
  stop("Table S6 (-a) is missing", call.=FALSE)
}
if (is.null(opt$cre)){
  print_help(opt_parser)
  stop("cCRE information (-b) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load Table S6\n")

QTL <- read.xlsx(opt$S6, sheet = "ATAC_QTL")
QTL <- QTL[grep("0.05",QTL$test),]

cat("TSS annotation\n")
QTL$seqnames <- gsub(":.*","",QTL$region) 
QTL$start <- as.numeric(gsub(".*:(.*)-.*","\\1",QTL$region))
QTL$end <- as.numeric(gsub(".*-","",QTL$region) )

QTL <- makeGRangesFromDataFrame(QTL, keep.extra.columns=TRUE,seqnames.field="seqnames")
  
#Peak annotation
peakAnno <- annotatePeak(unique(QTL[,-1]), TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion=c(-3000, 3000), annoDb="org.Mm.eg.db")

cat("Save plot 1\n")
svg(paste0(opt$outdir,"/Fig3a1.svg"))
plotAnnoPie(peakAnno, ndigit = 1, main='up')
dev.off()

cat("ENCODE annotation\n")
annotated <- readRDS(opt$cre)

#cCRE annotation up
print(paste('Percentage annotated up:', round(mean(1:length(QTL) %in% subjectHits(findOverlaps(annotated, QTL)))*100,2), '%'))
cre <- annotated[annotated$Peak_ID %in% QTL$regionID]

# Pie Chart with Percentages
slices <- round(table(cre$element) / length(cre)*100,1)
lbls <-  paste0(names(slices),": ", slices, "%")

cat("Save plot 2\n")

svg(paste0(opt$outdir,"/Fig3a2.svg"))
par(mar=c(0.5, 10, 0.5, 10))
pie(slices,labels = lbls, col=rainbow(length(lbls), v=0.75), cex = 0.8) 
dev.off()

# merge inputs for table
ov <- findOverlaps(as.GRanges(peakAnno),cre)

mcols(cre)$region <- paste(cre)
cre <- as.data.frame(cre[subjectHits(ov)])[,c('region', 'element')]
names(cre) <- paste0('region:ENCODE_',names(cre))

peakAnno <- as.GRanges(peakAnno)
mcols(peakAnno)$region <- paste(peakAnno)
tss <- data.frame(peakAnno[queryHits(ov)])[,c('region', 'annotation')]
names(tss) <- paste0('region:TSS_',names(tss))

tab <- cbind(tss, cre)

write.csv(tab, paste0(opt$outdir, "/data/3a.csv"))

sessionInfo()