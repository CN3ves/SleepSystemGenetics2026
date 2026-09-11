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

atac <- makeGRangesFromDataFrame(atac, keep.extra.columns=TRUE)
atac$region <- paste(atac)

cat("Annotation to TSS\n")
#TSS annotation up
tss_up <- atac[atac$ATAC_FDR <0.05 & atac$ATAC_logFC >0]
peakAnno_up <- annotatePeak(tss_up, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion=c(-3000, 3000), annoDb="org.Mm.eg.db")
print(paste('Percentage annotated up:', round(mean(1:length(tss_up) %in% subjectHits(findOverlaps(as.GRanges(peakAnno_up), tss_up)))*100,2), '%'))
print(paste('size of up:', length(tss_up)))

svg(paste0(opt$outdir,"/Fig1b1.svg"))
plotAnnoPie(peakAnno_up, ndigit = 1, main = "up")
dev.off()

#TSS annotation down
tss_down <- atac[atac$ATAC_FDR <0.05 & atac$ATAC_logFC < 0]
peakAnno_down <- annotatePeak(tss_down, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion=c(-3000, 3000), annoDb="org.Mm.eg.db")
print(paste('Percentage annotated down:', round(mean(1:length(tss_down) %in% subjectHits(findOverlaps(as.GRanges(peakAnno_down), tss_down)))*100,2), '%'))
print(paste('size of down:', length(tss_down)))

svg(paste0(opt$outdir,"/Fig1b2.svg"))
plotAnnoPie(peakAnno_down, ndigit = 1, main = "down")
dev.off()

cat("Annotation to ENCODE cCRE\n")
annotated <- readRDS(opt$cre)

#cCRE annotation up
print(paste('Percentage annotated up:', round(mean(1:length(tss_up) %in% subjectHits(findOverlaps(annotated, tss_up)))*100,2), '%'))
cre_up <- annotated[annotated$Peak_ID %in% tss_up$Region_ID]

slices <- round(table(cre_up$element) / length(cre_up)*100,2)
lbls <-  paste0(names(slices),": ", round(slices,1), "%")

svg(paste0(opt$outdir,"/Fig1b3.svg"))
par(mar=c(5.1, 4.1, 4.1, 10.5))
pie(slices,labels = lbls, col=rainbow(length(lbls), v=0.75), cex = 0.8, main = "up")
dev.off()

#cCRE annotation down
cre_down <- annotated[annotated$Peak_ID %in% tss_down$Region_ID]
print(paste('Percentage annotated down:', round(mean(1:length(tss_down) %in% subjectHits(findOverlaps(annotated, tss_down)))*100,2), '%'))

slices <- round(table(cre_down$element) / length(cre_down)*100,2)
lbls <-  paste0(names(slices),": ", round(slices,1), "%")

svg(paste0(opt$outdir,"/Fig1b4.svg"))
par(mar=c(5.1, 4.1, 4.1, 10.5))
pie(slices,labels = lbls, col=rainbow(length(lbls), v=0.75), cex = 0.8, main = "down")
dev.off()

# merge inputs for table
up <- findOverlaps(as.GRanges(peakAnno_up),cre_up)
down <- findOverlaps(as.GRanges(peakAnno_down),cre_down)

mcols(cre_up)$region <- paste(cre_up)
cre_up <- as.data.frame(cre_up[subjectHits(up)])[,c('region', 'element')]
names(cre_up) <- paste0('region:ENCODE_',names(cre_up))

mcols(cre_down)$region <- paste(cre_down)
cre_down <- as.data.frame(cre_down[subjectHits(down)])[,c('region', 'element')]
names(cre_down) <- paste0('region:ENCODE_',names(cre_down))

peakAnno_up <- as.GRanges(peakAnno_up)
mcols(peakAnno_up)$region <- paste(peakAnno_up)
tss_up <- as.data.frame(peakAnno_up[queryHits(up)])[,c('region', 'annotation')]
names(tss_up) <- paste0('region:TSS_',names(tss_up))
tss_up$Direction <- 'More accessible'

peakAnno_down <- as.GRanges(peakAnno_down)
mcols(peakAnno_down)$region <- paste(peakAnno_down)
tss_down <- as.data.frame(peakAnno_down[queryHits(down)])[,c('region', 'annotation')]
names(tss_down) <- paste0('region:TSS_',names(tss_down))
tss_down$Direction <- 'Less accessible'

tab <- cbind(rbind(tss_up,tss_down), rbind(cre_up,cre_down))

write.csv(tab, paste0(opt$outdir, "/data/1b.csv"))

sessionInfo()