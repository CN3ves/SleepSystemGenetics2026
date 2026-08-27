# Script designed to annotate regions to nearest trasncript

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/10-BXD_annotation/gene_10.3.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ggplot2")
  library("ggrepel")
  library("ChIPseeker")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene") 
  library("BSgenome.Mmusculus.UCSC.mm10")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-p", "--peaks"), type="character", default=NULL, 
              help="File with the list of genomic features", metavar="character"), 
  make_option(c("-d", "--dar"), type="character", default=NULL, 
              help="Rdata file with merged region for differential analysis", metavar="character"),
  make_option(c("-q", "--qtl"), type="character", default=NULL, 
              help="Rdata file with merged region for QTL analysis", metavar="character"),
  make_option(c("-i", "--int"), type="character", default=NULL, 
              help="Rdata file with merged region for QTL interactions analysis", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

   
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$peaks)){
  print_help(opt_parser)
  stop("File with the list of genomic features (-p) is missing", call.=FALSE)
}
if (is.null(opt$dar)){
  print_help(opt_parser)
  stop("Rdata file with merged region for differential analysis (-d) is missing", call.=FALSE)
}
if (is.null(opt$qtl)){
  print_help(opt_parser)
  stop("Rdata file with merged region for QTL analysis (-q) is missing", call.=FALSE)
}
if (is.null(opt$int)){
  print_help(opt_parser)
  stop("Rdata file with merged region for QTL interactions analysis (-i) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
peaks <- read.delim(opt$peaks, header = FALSE, col.names = c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame",  "attribute"))
peaks <- makeGRangesFromDataFrame(peaks, keep.extra.columns=TRUE)
peaks <- peaks[-as.numeric(grep('_',seqnames(peaks)))]

all <- GRanges(seqnames=seqnames(BSgenome.Mmusculus.UCSC.mm10), 
  IRanges(start=1, end=seqlengths(BSgenome.Mmusculus.UCSC.mm10)))
all <- all[-as.numeric(grep('_',seqnames(all)))]
all <- unlist(tile(all, width=1000))

peaks <- list("All"= all,
              "Peaks" = reduce(peaks),
              "DAR" = readRDS(opt$dar),
              "QTL" = readRDS(opt$qtl),
              "QTLxSD" = readRDS(opt$int))

#Peak annotation
peakAnno <- lapply(peaks, annotatePeak, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion=c(-3000, 3000), annoDb="org.Mm.eg.db")
  
print(peakAnno)

for(test in names(peakAnno)) {
  cat(paste("Plotting for", test,"\n"))
  png(paste0(opt$outdir,"/Annot_pie_", test,".png"))
  plotAnnoPie(peakAnno[[test]])
  dev.off()
  
  png(paste0(opt$outdir,"/Annot_bar_", test,".png"))
  print(plotAnnoBar(peakAnno[[test]]))
  dev.off()
  
  png(paste0(opt$outdir,"/Annot_vennpie_", test,".png"))
  vennpie(peakAnno[[test]])
  dev.off()
  
  #too long and not much more information 
  #png(paste0(opt$outdir,"/Annot_upset_", test,".png"))
  #print(upsetplot(peakAnno[[test]]))
  #dev.off()
  
  # Visualize distribution of TF-binding loci relative to TSS
  png(paste0(opt$outdir,"/Annot_TSSdist_", test,".png"))
  print(plotDistToTSS(peakAnno[[test]],
              title="Distribution of transcription factor-binding loci\nrelative to TSS"))
  dev.off()
  
  annotated_peaks <- as.GRanges(peakAnno[[test]])
  annotated_peaks <- as.data.frame(annotated_peaks)
  if(!any(names(annotated_peaks) == 'FDR')) {
    cat("No statistics\n")
    annotated_peaks$FDR <- -1
  }
  annotated_peaks <- annotated_peaks[order(annotated_peaks$FDR),]
  
  write.csv(annotated_peaks, paste0(opt$outdir,"/",test,"_annotated.csv"))
}
cat("Saving all annotation\n")
saveRDS(peakAnno,paste0(opt$outdir,"/ranges_annotated.RData"))

# Distant to promoter plot
peaks_TSS <- lapply(peakAnno, function(x) as.GRanges(x)$distanceToTSS)
xlims <- round(max(abs(range(peaks_TSS)))/1000000)*1000000
ylims <- max(sapply(peaks_TSS, function(x) max(density(x)$y)))

png(paste0(opt$outdir,"/distance_TSS.png"))
plot(density(peaks_TSS[['All']]), main="Distribution of distances to TSS", xlab="Distance to TSS", xlim=c(-xlims,xlims), ylim = c(0, ylims), col="black")
lines(density(peaks_TSS[['DAR']]), col = "red")
lines(density(peaks_TSS[['QTL']]), col = "green", lty = "dashed")
lines(density(peaks_TSS[['QTLxSD']]), col = "blue", lty = "dashed")
abline(v=0, col="grey")
legend("topright", legend=c("Tested","DAR", "QTL", "QTLxSD"), lwd =1, col=c("black","red", "green", "blue"))
dev.off()

png(paste0(opt$outdir,"/distance_TSSzoom.png"))
plot(density(peaks_TSS[['All']]), main="Distribution of distances to TSS (withing 50Kb)", xlab="Distance to TSS", xlim=c(-50000,50000), ylim = c(0, ylims), col="black")
lines(density(peaks_TSS[['DAR']]), col = "red")
lines(density(peaks_TSS[['QTL']]), col = "green", lty = "dashed")
lines(density(peaks_TSS[['QTLxSD']]), col = "blue", lty = "dashed")
abline(v=0, col="grey")
legend("topright", legend=c("Tested","DAR", "QTL", "QTLxSD"), lwd =1, col=c("black","red", "green", "blue"))
dev.off()

print("Fraction of significant peaks 50Kb away from TSS:")
print(paste0(names(peaks_TSS) , ":", sapply(peaks_TSS, function(x) round(mean(abs(x) > 50000),4)), collapse = "; "))
print("Fraction of significant within 5Kb from TSS:")
print(paste0(names(peaks_TSS) , ":", sapply(peaks_TSS, function(x) round(mean(abs(x) < 5000),4)), collapse = "; "))
print("Fraction of significant peaks 5Kb upstream of TSS:")
print(paste0(names(peaks_TSS) , ":", sapply(peaks_TSS, function(x) round(mean(x < 5000 &  x >0 ),4)), collapse = "; "))


sessionInfo()