# Script designed to run GRANIE

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/13-BXD_integrate/granie_13.4.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("readr")
  library("GRaNIE")
  library("ggplot2")
  library("ggrepel")
  library("limma")
  library("edgeR")
  library("org.Mm.eg.db")
  library("biomaRt")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-r", "--rna"), type="character", default=NULL, 
              help="RNA raw counts", metavar="character"),
  make_option(c("-a", "--atac"), type="character", default=NULL, 
              help="ATAC raw counts", metavar="character"),
  make_option(c("-p", "--regions"), type="character", default=NULL, 
              help="List of genomice regions", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$rna)){
  print_help(opt_parser)
  stop("RNA raw counts (-r) is missing", call.=FALSE)
}
if (is.null(opt$atac)){
  print_help(opt_parser)
  stop("ATAC raw counts (-a) is missing", call.=FALSE)
}
if (is.null(opt$regions)){
  print_help(opt_parser)
  stop("List of genomice regions (-p) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

# Helper function Granie does not access decimals, so sum of reads is taken
collapse_counts <- function(obj) {
  
  sum_counts <- lapply(unique(obj$samples$group), function(group) {
    samples <- obj$samples$group == group
    
    group <- gsub("BXD0*", "BXD", group)
    group <- gsub("96","48a", group)
    group <- gsub("97","65a", group)
    group <- gsub("103","73b", group)
    
    sum_count <- data.frame(rowSums(edgeR::getCounts(obj)[,samples, drop = FALSE]))
    names(sum_count) <- group
    
    return(sum_count)
  })
  
  as.data.frame(sum_counts) 
}

cat("Reading files\n")
folder_TFBS_6TFs = "https://www.embl.de/download/zaugg/diffTF/TFBS/TFBS_mm10_PWMScan_HOCOMOCOv10.tar.gz"
download.file(folder_TFBS_6TFs, file.path(paste0(opt$outdir,"/HOCOMOCOv10.tar.gz")), quiet = FALSE)
untar(file.path(paste0(opt$outdir,"/HOCOMOCOv10.tar.gz")), exdir = opt$outdir)
    
countsPeaks = readRDS(opt$atac)
countsRNA = readRDS(opt$rna)
countsPeaks <- collapse_counts(countsPeaks)
countsPeaks <- countsPeaks[-grep("chrUn",rownames(countsPeaks)),]
countsPeaks <- countsPeaks[-grep("random",rownames(countsPeaks)),]

#NAMES MUST BE chr:start-end
regions <- read.delim(opt$regions, header =FALSE)
regions$V9 <- gsub("Peak_ID (.*);","\\1",regions$V9)
regions$interval <- paste0(regions[,1],":",regions[,4],"-",regions[,5])
regions <- regions[regions$V9 %in%  rownames(countsPeaks),]
countsPeaks[regions$V9,"PeakID"] <- regions$interval
countsPeaks <- countsPeaks[,c(ncol(countsPeaks),1:(ncol(countsPeaks)-1))]
rownames(countsPeaks) <- NULL

countsRNA <- collapse_counts(countsRNA)
countsRNA$GeneID <- rownames(countsRNA)
countsRNA <- countsRNA[,c(ncol(countsRNA),1:(ncol(countsRNA)-1))]
rownames(countsRNA) <- NULL  

countsRNA$GeneID <- mapIds(org.Mm.eg.db, keys = countsRNA$GeneID, keytype = "SYMBOL", column="ENSEMBL")
countsRNA <- countsRNA[!is.na(countsRNA$GeneID),]

cat("Preparing sample data\n")
sampleMetadata <- c(colnames(countsRNA),colnames(countsPeaks))
sampleMetadata <- unique(grep("ID", sampleMetadata, invert = TRUE, value=T))
sampleMetadata <- data.frame("SampleID" = sampleMetadata,
                             "Strain" = gsub("_.*", "",sampleMetadata),
                             "Treatment" = gsub(".*_", "",sampleMetadata))
 
Peaks <- countsPeaks[,c(TRUE,sapply(colnames(countsPeaks), function(x) x %in% sampleMetadata$SampleID)[-1])]
RNA <-  countsRNA[,c(TRUE,sapply(colnames(countsRNA), function(x) x %in%   sampleMetadata$SampleID)[-1])]
print(dim(Peaks))
print(dim(RNA))

write_tsv(Peaks,paste0(opt$outdir,"/peaks.tsv"))
write_tsv(RNA,paste0(opt$outdir,"/rna.tsv"))
write_tsv(sampleMetadata,paste0(opt$outdir,"/meta.tsv"))

cat("Initialiye GRN\n")
objectMetadata = list(name = opt$outdir,
                      file_peaks = paste0(opt$outdir,"/peaks.tsv"),
                      file_rna = paste0(opt$outdir,"/rna.tsv"),
                      file_sampleMetadata = paste0(opt$outdir,"/meta.tsv"),
                      genomeAssembly = "mm10")

GRN = initializeGRN(objectMetadata = objectMetadata, outputFolder = opt$outdir,
                        genomeAssembly = "mm10")
    
    
print("MART CHECKPOINT") 
GRN = addData(GRN, counts_peaks = Peaks, normalization_peaks = "DESeq2_sizeFactors",
              idColumn_peaks = "PeakID", counts_rna = RNA, normalization_rna = "limma_quantile",
              idColumn_RNA = "GeneID", sampleMetadata = sampleMetadata, forceRerun = TRUE)

GRN = plotPCA_all(GRN, data = c("rna","peaks"), topn = 500, type = "normalized", plotAsPDF = TRUE,   forceRerun   = TRUE)
    
motifFolder = tools::file_path_as_absolute(paste0(opt$outdir,"/PWMScan_HOCOMOCOv10"))

GRN = addTFBS(GRN, motifFolder = motifFolder, TFs = "all", filesTFBSPattern = "_TFBS",
              fileEnding = ".bed.gz", forceRerun = TRUE)

GRN = overlapPeaksAndTFBS(GRN, nCores = 20, forceRerun = TRUE)

chrToKeep_peaks = c(paste0("chr", 1:22), "chrX")
GRN = filterData(GRN, minNormalizedMean_peaks = 5, minNormalizedMeanRNA = 1, chrToKeep_peaks = chrToKeep_peaks, maxSize_peaks = 10000, forceRerun = TRUE)
    
GRN = addConnections_TF_peak(GRN, plotDiagnosticPlots = FALSE, connectionTypes = c("expression"),
                             corMethod = "pearson", forceRerun = TRUE)

GRN = plotDiagnosticPlots_TFPeaks(GRN, dataType = c("real"), plotAsPDF = TRUE)

GRN = AR_classification_wrapper(GRN, significanceThreshold_Wilcoxon = 0.05, outputFolder =     paste0(opt$outdir,"/plots"),plot_minNoTFBS_heatmap = 100, plotDiagnosticPlots = TRUE, forceRerun = TRUE)

print("SAVEPOINT 1")
GRN_file_outputRDS = paste0(opt$outdir,"/GRN.rds")

GRN = addConnections_peak_gene(GRN, corMethod = "pearson", promoterRange = 250000,
                               TADs = NULL, nCores = 20, plotDiagnosticPlots = FALSE, plotGeneTypes = list(c("all")), forceRerun = TRUE) 
    
GRN = plotDiagnosticPlots_peakGene(GRN, gene.types = list(c("protein_coding", "lincRNA")),  plotAsPDF = TRUE)
    
GRN = filterGRNAndConnectGenes(GRN, TF_peak.fdr.threshold = 0.2, peak_gene.fdr.threshold = 0.2, peak_gene.fdr.method = "BH", gene.types = c("protein_coding", "lincRNA"), allowMissingTFs = FALSE, allowMissingGenes = FALSE, forceRerun = TRUE)
    
GRN = add_TF_gene_correlation(GRN, corMethod = "pearson", nCores = 20, forceRerun = TRUE)
    
GRN_connections.all = getGRNConnections(GRN, type = "all.filtered", include_TF_gene_correlations = TRUE, include_geneMetadata = TRUE)
    
print("SAVE NET")
saveRDS(GRN_connections.all, paste0(opt$outdir,"/graph.RData"))

sessionInfo()