# Script to produce table S1

# Redirect all R logs to Snakemake log
log <- file('logs/19-Tables/S1_19.1.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("openxlsx")
  library("tidyverse")
  library("org.Mm.eg.db")
  library("ChIPseeker")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
  library("BSgenome.Mmusculus.UCSC.mm10")
  library("foreach")
  library("doParallel")
  library("doSNOW")
  library("annotatr")
  library("optparse")
})

cat("Checking arguments\n")

option_list = list(
  make_option(c("-a", "--pre_align"), type="character", default=NULL, 
              help="Read QC results' table", metavar="character"),
  make_option(c("-b", "--post_align"), type="character", default=NULL, 
              help="Alignment QC results' table", metavar="character"),
  make_option(c("-c", "--bam"), type="character", default=NULL, 
              help="Bam QC results' table", metavar="character"),
  make_option(c("-d", "--peaks"), type="character", default=NULL, 
              help="Peak calling logs' directory", metavar="character"),
  make_option(c("-e", "--call"), type="character", default=NULL, 
              help="Peakcalling QC results' table", metavar="character"),
  make_option(c("-f", "--bed"), type="character", default=NULL, 
              help="Merged peak list", metavar="character"),
  make_option(c("-g", "--gtf"), type="character", default=NULL, 
              help="Genomic regions of interest", metavar="character"),
  make_option(c("-i", "--cre"), type="character", default=NULL, 
              help="ENCODE cCRE list", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$pre_align)){
  print_help(opt_parser)
  stop("Read QC results' table (-a) is missing", call.=FALSE)
}
if (is.null(opt$post_align)){
  print_help(opt_parser)
  stop("Alignment QC results' table (-b) is missing", call.=FALSE)
}
if (is.null(opt$bam)){
  print_help(opt_parser)
  stop("Bam QC results' table (-c) is missing", call.=FALSE)
}
if (is.null(opt$peaks)){
  print_help(opt_parser)
  stop("Peak calling logs' directory (-d) is missing", call.=FALSE)
}
if (is.null(opt$call)){
  print_help(opt_parser)
  stop("Peakcalling QC results' table (-e) is missing", call.=FALSE)
}
if (is.null(opt$bed)){
  print_help(opt_parser)
  stop("Merged peak list (-f) is missing", call.=FALSE)
}
if (is.null(opt$gtf)){
  print_help(opt_parser)
  stop("Genomic regions of interest (-g) is missing", call.=FALSE)
}
if (is.null(opt$cre)){
  print_help(opt_parser)
  stop("ENCODE cCRE list (-i) is missing", call.=FALSE) 
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

# Create a new workbook and add a sheet
wb <- createWorkbook()

cat("Load sequencing QC tab\n")
pre_align <- read.csv(opt$pre_align)
pre_align <- pre_align[-grep("adapter_counts", pre_align$X),]
names(pre_align)[1:2] <- c("Metric", "Summary")
names(pre_align)[-1:-2] <- gsub("^X","",names(pre_align)[-1:-2])
names(pre_align)[-1:-2] <- paste0(sub("_"," (",names(pre_align)[-1:-2]),")")

addWorksheet(wb, "Sequencing_QC")
writeData(wb, "Sequencing_QC", pre_align, rowNames=FALSE)

cat("Load alignmnet QC tab\n")
post_align <- read.csv(opt$post_align)
rownames(post_align) <- post_align$File
bam_qc <- read.csv(opt$bam)
bam_qc <- bam_qc[, -which(names(bam_qc)== 'X')]
rownames(bam_qc) <- bam_qc$sample

bam_qc$effective_depth <- round(as.numeric(post_align[bam_qc$sample,'Reads']) * 
  as.numeric(bam_qc$nonRedundantFraction) * 
  as.numeric(bam_qc$Align_rate_chrM) * 
  (as.numeric(gsub('%','',post_align[bam_qc$sample,'Unaligned']))*100))

# Summarize data
bam_qc <- rbind(c("", apply(bam_qc[,-1], 2, function(col) {
  col <- as.numeric(col)
  col <- summary(col)
  col <- paste0("Min: ", round(col["Min."],3), 
    "; Max: ", round(col["Max."],3), "; Mean: ", round(col["Mean"],3))
  return(col)
})),bam_qc)
  
post_align <- merge(post_align,bam_qc, all.x=TRUE, by.x="File",by.y="sample")
post_align <- t(post_align)
post_align <- cbind(rownames(post_align),post_align)
post_align <- as.data.frame(post_align)
names(post_align) <- post_align['File',]
names(post_align)[1:2] <- c("Metric", "Summary")
names(post_align)[-1:-2] <- paste0(sub("_"," (",names(post_align)[-1:-2]),")")

post_align <- post_align[-1:-2,]

addWorksheet(wb, "Alignment_QC")
writeData(wb, "Alignment_QC", post_align, rowNames=FALSE)

cat("Load peak calling QC tab\n")
peak_qc <- read.csv(opt$call)

names(peak_qc)[1:2] <- c("Metric", "Summary")
names(peak_qc)[-1:-2] <- gsub("^X","",names(peak_qc)[-1:-2])
names(peak_qc)[-1:-2] <- paste0(sub("_"," (",names(peak_qc)[-1:-2]),")")

addWorksheet(wb, "Peak_QC")
writeData(wb, "Peak_QC", peak_qc, rowNames=FALSE)

cat("Load peak calling number\n")
peaks_files <- list.files(opt$peaks, pattern="log",recursive=TRUE, full.name= TRUE)
peaks_call <- lapply(peaks_files[-1:-2], function(file) {
  x <- readLines(file)
  x <- x[grep("Peaks identified",x)]
  x <- gsub(".*: ", "",x)
  x <- data.frame(Sample= gsub(".*/(.*).out", "\\1", file), Peaks = x)#as.numeric(x))
  return(x)
})
peaks_call <- do.call(rbind, peaks_call)
peaks_call$Sample <- basename(peaks_call$Sample)

peaks_call$Set <- gsub("_.*", "", peaks_call$Sample)
n <- strsplit(peaks_call$Sample, '_')
n <- sapply(n, function(x) ifelse(length(x)==4,gsub(".log","",x[4]),''))
peaks_call$Set <- paste0(peaks_call$Set,n)

peaks_call$Strain <- gsub("_*[0-9]*.log", "", peaks_call$Sample)
peaks_call$Treatment <- gsub(".*_", "", peaks_call$Strain)
peaks_call$Strain <- gsub(".*_(.*)_.*", "\\1", peaks_call$Strain)
peaks_call$Size <- gsub(".* \\((.*)bp)", "\\1 bp", peaks_call$Peaks)

genome_size <- sum(seqlengths(BSgenome.Mmusculus.UCSC.mm10))
peaks_call$Perc <- paste0(round(as.numeric(gsub(" bp", "", peaks_call$Size))/genome_size *100,2),"%")
peaks_call$Peaks <- gsub(" .*","",peaks_call$Peaks) 

peaks_call <- peaks_call %>% dplyr::select(!c(Sample)) %>% pivot_wider(names_from=Treatment, values_from=c(Peaks,Size,Perc))

names(peaks_call) <- c('Set', 'Strain', 'CTRL: n of peaks', 'SD: n of peaks', 
  'CTRL: bp covered','SD: bp covered', 'CTRL: % of genome', 'SD: % of genome')
  
addWorksheet(wb, "Peak_Summary")
writeData(wb, "Peak_Summary", peaks_call, rowNames=FALSE)

cat("Load features tab\n") 
merged_peaks <- read.delim(opt$bed, header = FALSE)
names(merged_peaks) <- c("Chromosome", "Start", "End")

ranges <- GRanges(IRanges(start=as.numeric(merged_peaks$Start), 
                                 end=as.numeric(merged_peaks$End)), 
                         seqnames=merged_peaks$Chromosome)
annot <- annotatePeak(ranges, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion=c(-3000, 3000),   annoDb="org.Mm.eg.db")
annot <- as.GRanges(annot)

annot <- c(annot,ranges[which(!ranges %in% annot)])

merged_peaks$'Nearest Gene' <- annot$SYMBOL[match(paste(ranges), paste(annot))]

addWorksheet(wb, "Merged_Peaks")
writeData(wb, "Merged_Peaks", merged_peaks, rowNames=FALSE)

cat("Split regions\n") 
cre <- read.delim(opt$cre, header= FALSE) 
names(cre) <- c("seqnames", "start", "end", "something", "else","element")
cre <- GRanges(cre$seqnames, IRanges(start=cre$start, end=cre$end), mcols=cre[,"element", drop = FALSE])

split_peaks <- read.delim(opt$gtf, header = FALSE)

split_peaks <- split_peaks[,c(1,4,5,9)]
names(split_peaks) <- c("Chromosome", "Start", "End", "Peak_ID")
split_peaks$Peak_ID <- gsub("Peak_ID (.*);","\\1",split_peaks$Peak_ID)
split_peaks$Length <- split_peaks$End - split_peaks$Start

ranges <- GRanges(IRanges(start=as.numeric(split_peaks$Start), 
                                 end=as.numeric(split_peaks$End)), 
                         seqnames=split_peaks$Chromosome,
                  Peak_ID = split_peaks$Peak_ID)

annot <- annotatePeak(ranges, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion=c(-3000, 3000),   annoDb="org.Mm.eg.db")
annot <- as.GRanges(annot)
annot <- c(annot,ranges[which(!ranges %in% annot)])

annot_cre = annotate_regions(
  regions = annot,
  annotations = cre,
  ignore.strand = TRUE,
  quiet = FALSE)

mcols(annot_cre)$element <- annot_cre$annot$mcols.element
mcols(annot_cre) <- mcols(annot_cre)[, -which(names(mcols(annot_cre))=="annot")]

saveRDS(annot_cre, paste0(opt$outdir,"/data/cCRE.RData"))

split_peaks$'Nearest Gene' <- annot$SYMBOL[match(paste(ranges), paste(annot))]
split_peaks$'Annotation' <- annot$annotation[match(paste(ranges), paste(annot))]

cat("Compute cCRE for regions\n") 
cl <- parallel::makeCluster(20)
doParallel::registerDoParallel(cl)
paths <- .libPaths()
clusterExport(cl, "paths")
clusterEvalQ(cl, .libPaths(paths))

cCRE <- foreach(i = 1:nrow(split_peaks), .combine = rbind)  %dopar% {
  idx <- which(annot_cre$Peak_ID %in% split_peaks$Peak_ID[i])
  
  res <- data.frame("Peak_ID"=split_peaks$Peak_ID[i], "cCRE"= paste(sort(unique(annot_cre[idx]$element)), collapse="; "))
    
  return(res)
} 

parallel::stopCluster(cl)

saveRDS(cCRE, paste0(opt$outdir,"/data/cCRE_summary.RData"))

split_peaks <- merge(split_peaks, cCRE)

split_peaks <- split_peaks[c("Peak_ID","Chromosome", "Start", "End", "Length", "Nearest Gene","cCRE","Annotation")]

addWorksheet(wb, "Final_Regions")
writeData(wb, "Final_Regions", split_peaks, rowNames=FALSE)

cat("Save Table\n") 
saveWorkbook(wb, paste0(opt$outdir,"/TableS1-ATAC_QC_summary.xlsx"), overwrite = TRUE)

sessionInfo()