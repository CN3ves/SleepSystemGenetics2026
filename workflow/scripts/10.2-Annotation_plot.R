# Script designed to get description and QC of significant regions

# Redirect all R logs to Snakemake log
log <- file('logs/10-BXD_annotation/plots_10.2.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ChIPseeker")
  library("ggplot2")
  library("dplyr")
  library("BSgenome.Mmusculus.UCSC.mm10")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene") 
  library("GenomicRanges")
  library("optparse")
})

source("workflow/scripts/10.0-Annotation_helper.R")

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
  make_option(c("-b", "--basesnps"), type="character", default=NULL, 
              help="File containing the filtered QTL results", metavar="character"),
  make_option(c("-f", "--fcsnps"), type="character", default=NULL, 
              help="File containing the filtered QTL results", metavar="character"),
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
if (is.null(opt$basesnps)){
  print_help(opt_parser)
  stop("File containing the filtered QTL results (-b) is missing", call.=FALSE)
}
if (is.null(opt$fcsnps)){
  print_help(opt_parser)
  stop("File containing the filtered iteraction QTL results (-f) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
snps_qtl <-  read.csv(opt$basesnps)
snps_int <-  read.csv(opt$fcsnps)

peaks <- read.delim(opt$peaks, header = FALSE, col.names = c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame",  "attribute"))
peaks <- makeGRangesFromDataFrame(peaks, keep.extra.columns=TRUE)
peaks <- peaks[-as.numeric(grep('_',seqnames(peaks)))]

peaks <- list("All" = reduce(peaks),
              "DAR" = readRDS(opt$dar),
              "QTL" = readRDS(opt$qtl),
              "QTLxSD" = readRDS(opt$int))


cat(paste("Number of contiguous regions:\n", paste(names(peaks), sapply(peaks, length), collapse = "\n"),"\n"))

# Percentage of genome changed

cg_text <- function(results, text, direction=NA) {
  if (is.na(direction)) {
    if (!'direction' %in% names(results)) {
        results$direction <- 'all has no have direction'
    } 
    direction = unique(results$direction)    
  }
  direction <- c("", direction) # Make sure it is a vector
  cg <- sum(width(results[results$direction %in% direction]))/sum(seqlengths(BSgenome.Mmusculus.UCSC.mm10)[grep("[_MY]",seqnames(BSgenome.Mmusculus.UCSC.mm10), invert = T)])*100

  print(paste0("Percentage of ", text, " (excluding M/Y): ", round(cg,3), "%"))

}
cg_text(peaks[["All"]], "tested genome")
cg_text(peaks[["DAR"]], "differentially accessible genome")
cg_text(peaks[["DAR"]], "differentially more accessible genome", "up")
cg_text(peaks[["DAR"]], "differentially less accessible genome", "down")
cg_text(peaks[["QTL"]], "regions related to genotype")
cg_text(peaks[["QTLxSD"]], "regions with genotype x treatment interactions")

# CG enrichment
seqs <- list("All" = BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, peaks[["All"]]),
            "DAR" = BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, peaks[["DAR"]]),
            "QTL" =BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, peaks[["QTL"]]),
            "QTLxSD" = BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, peaks[["QTLxSD"]]))
 
cg_content  <- list(
  "All" = as.numeric(Biostrings::letterFrequency(x = seqs[['All']], letters = "GC", as.prob = TRUE)),
  "DAR" = as.numeric(Biostrings::letterFrequency(x = seqs[['DAR']], letters = "GC", as.prob = TRUE)),
  "QTL" = as.numeric(Biostrings::letterFrequency(x = seqs[['QTL']], letters = "GC", as.prob = TRUE)),
  "QTLxSD" = as.numeric(Biostrings::letterFrequency(x = seqs[['QTLxSD']], letters = "GC", as.prob = TRUE)))


cat("Peaks CG content distribution:\n")
dump <- sapply(1:length(cg_content), function(i) {
  cat(names(cg_content)[i], "\n")
  cat(paste0(names(summary(cg_content[[i]])), ": ", round(summary(cg_content[[i]]),3), collapse= "; "), "\n")
})

png(paste0(opt$outdir,"/cg_content.png"))
plot(density(cg_content[['All']]), main="Distribution of GC content for significant reads", xlab="CG content", ylim= c(0,10))
lines(density(cg_content[['DAR']]), col = "red")
lines(density(cg_content[['QTL']]), col = "green")
lines(density(cg_content[['QTLxSD']]), col = "blue", lty = "dashed")
abline(v=0.42, col="grey")
legend("topright", legend=c("Tested","DAR", "QTL", "QTLxSD", "mm10 GC mean"), lwd =1, col=c("black","red", "green", "blue","gray"))
dev.off()


# ChIP peaks coverage plot
cat("Covplot DAR:\n")
png(paste0(opt$outdir,"/covplot_DAR.png"))
plot_covfc(peaks[["DAR"]], title = "Differentially accessible regions")
dev.off()

cat("Covplot QTL:\n")
plot <-  qtl_plot(peaks= peaks[["QTL"]], snps=snps_qtl)
png(paste0(opt$outdir,"/covplot_QTL.png"))
plot_covfc(plot, title = "QTL regions",labels=c("QTL","SNPs","Chr limit"), legend= "LOD score")
dev.off()

cat("Covplot interactions:\n")
plot <-  qtl_plot(peaks= peaks[["QTLxSD"]], snps=snps_int)
png(paste0(opt$outdir,"/covplot_QTLxSD.png"))
plot_covfc(plot, title = "QTLxSD regions", labels=c("QTL","SNPs","Chr limit"), legend= "LOD score")
dev.off()

# Profile plots
cat("Profile:\n")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
tag_plots(3000, peaks[-1], opt$outdir)

sessionInfo()