# Script to produce figure 1b

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure1b.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("openxlsx")
  library("GenomicRanges")
  library("annotatr")
  library("ChIPseeker")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene") 
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--S2"), type="character", default=NULL, 
              help="Table S2", metavar="character"),
  make_option(c("-b", "--cre"), type="character", default=NULL, 
              help="cCRE information", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$S2)){
  print_help(opt_parser)
  stop("Table S2 (-a) is missing", call.=FALSE)
}
if (is.null(opt$cre)){
  print_help(opt_parser)
  stop("cCRE information (-b) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load Table S2\n")
atac <- read.xlsx(opt$S2, sheet="Differential_Accessibility")
par(mfrow = c(2, 2))

atac <- makeGRangesFromDataFrame(atac, keep.extra.columns=TRUE)
atac$region <- paste(atac)

cat("Annotation to TSS\n")
#TSS annotation up
tss_up <- atac[atac$ATAC_FDR <0.05 & atac$ATAC_logFC >0]
peakAnno <- annotatePeak(tss_up, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion=c(-3000, 3000), annoDb="org.Mm.eg.db")
print(paste('Percentage annotated up:', round(mean(1:length(tss_up) %in% subjectHits(findOverlaps(as.GRanges(peakAnno), tss_up)))*100,2), '%'))
print(paste('size of up:', length(tss_up)))

svg(paste0(opt$outdir,"/Fig1b.svg"))
plotAnnoPie(peakAnno, ndigit = 1, main = "up")
  
#TSS annotation down
tss_down <- atac[atac$ATAC_FDR <0.05 & atac$ATAC_logFC < 0]
peakAnno <- annotatePeak(tss_down, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion=c(-3000, 3000), annoDb="org.Mm.eg.db")
print(paste('Percentage annotated down:', round(mean(1:length(tss_down) %in% subjectHits(findOverlaps(as.GRanges(peakAnno), tss_down)))*100,2), '%'))
print(paste('size of down:', length(tss_down)))

plotAnnoPie(peakAnno, ndigit = 1, main = "down")

cat("Annotation to ENCODE cCRE\n")
annotated <- readRDS(opt$cre)

#cCRE annotation up
print(paste('Percentage annotated up:', round(mean(1:length(tss_up) %in% subjectHits(findOverlaps(annotated, tss_up)))*100,2), '%'))
up <- annotated[annotated$Peak_ID %in% tss_up$Region_ID]

slices <- round(table(up$element) / length(up)*100,2)
lbls <-  paste0(names(slices),": ", round(slices,1), "%")

pie(slices,labels = lbls, col=rainbow(length(lbls), v=0.75), cex = 0.8, main = "up")

#cCRE annotation down
down <- annotated[annotated$Peak_ID %in% tss_down$Region_ID]
print(paste('Percentage annotated down:', round(mean(1:length(tss_down) %in% subjectHits(findOverlaps(annotated, tss_down)))*100,2), '%'))

slices <- round(table(down$element) / length(down)*100,2)
lbls <-  paste0(names(slices),": ", round(slices,1), "%")

pie(slices,labels = lbls, col=rainbow(length(lbls), v=0.75), cex = 0.8, main = "down")
dev.off()

sessionInfo()