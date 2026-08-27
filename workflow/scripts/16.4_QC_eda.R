# Script designed to filter counts from ATAC count files

# Redirect all R logs to Snakemake log
log <- file('logs/16-NRF_bam/eda_16.4.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("XML")
  library("EDASeq")
  library("edgeR")
  library("biomaRt")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-f", "--fastq"), type="character", default=NULL, 
              help="Processed FASTQ read files", metavar="character"),
  make_option(c("-b", "--bam"), type="character", default=NULL, 
              help="Sorted aligmnetd BAM files", metavar="character"),
  make_option(c("-c", "--count"), type="character", default=NULL, 
              help="Gene-level count files", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="Processed metadata file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$fastq)){
  print_help(opt_parser)
  stop("Processed FASTQ read files (-f) is missing", call.=FALSE)
}
if (is.null(opt$bam)){
  print_help(opt_parser)
  stop("Sorted aligmnetd BAM files (-b) is missing", call.=FALSE)
}
if (is.null(opt$count)){
  print_help(opt_parser)
  stop("Gene-level count files (-c) is missing", call.=FALSE)
}
if (is.null(opt$meta)){
  print_help(opt_parser)
  stop("Processed metadata file (-m) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
# Unaligned reads
files <-  strsplit(opt$fastq, " ")[[1]]
names(files) <- gsub("\\.fq.*", "", basename(files))
meta <- read.csv(opt$meta, row.names=1)
rownames(meta) <- meta$SampleName

fastq <- FastqFileList(files[meta$SampleName])
elementMetadata(fastq) <- meta

# Aligned reads
files <- strsplit(opt$bam, " ")[[1]]
names(files) <- gsub(".co.bam", "", basename(files))
bfs <- BamFileList(files[meta$SampleName])
elementMetadata(bfs) <- meta

# Counts
counts_files <- strsplit(opt$count, " ")[[1]]
counts <- lapply(counts_files, function(file) {
  df <- read.delim(file, header=FALSE)[,1:2]
  colnames(df) <-  c('ID', gsub(".*/(.*)_ReadsPerGene.*", "\\1",file))
  df <- df[-grep("N_", df$ID),]
  return(df)
})
counts <- Reduce(function(x, y) merge(x, y, by="ID"), counts)
rownames(counts) <- counts$ID
counts <- counts[,-1]

cat("Read-level EDA\n")
#Numbers of unaligned and aligned reads(((((((())))))))
cols <- as.factor(meta$Group)
names(cols) <- cols
levels(cols) <- rainbow(length(cols))

png(paste0(opt$outdir,"/EDASeq_alined_reads.png"))
barplot(bfs,las=2,col=cols)
legend("top", legend = unique(names(cols)), fill = unique(cols), horiz=TRUE)
dev.off()

png(paste0(opt$outdir,"/EDASeq_total_reads.png"))
barplot(fastq,las=2,col=cols)
legend("top", legend = unique(names(cols)), fill = unique(cols), horiz=TRUE)
dev.off()

png(paste0(opt$outdir,"/EDASeq_percentage_mapped.png"))
plot(x=bfs, y=fastq,las=2,col=cols)
legend("topright", legend = names(cols), fill = cols)
dev.off()

#Read quality scores
png(paste0(opt$outdir,"/EDASeq_read_quality.png"))
plotQuality(bfs,col=cols,lty=1)
legend("topright", legend = unique(names(cols)), fill = unique(cols), horiz=TRUE)
dev.off()

# Individual lane summaries
for (i in 1:length(bfs)) {
  png(paste0(opt$outdir,"/EDASeq_read_quality",i,".png"))
  barplot(bfs[[i]],las=2)
  dev.off()
}

cat("Gene-level EDA\n")
# Get GC content for genes
gc <- getGeneLengthAndGCContent(id=gsub("\\.[0-9]*","",rownames(counts)), org="mm10",mode="org.db")
rownames(gc) <- rownames(counts)

#filter the non-expressed genes
filter <- apply(counts,1,function(x) mean(x)>10)
table(filter)

common <- intersect(rownames(gc),rownames(counts[filter,]))
length(common)

#filter by have CG content information
sum(is.na(gc[common,"gc"]))

rm_cg <- names(gc[common,"gc"][is.na(gc[common,"gc"])])
common <- common[!common%in% rm_cg]

data <- newSeqExpressionSet(counts=as.matrix(counts[common,]),
                            featureData=data.frame(gc=gc[common,"gc"],length=gc[common,"length"]),
                            phenoData=meta)

stopifnot(all(rownames(counts(data)) == rownames(featureData(data))))
stopifnot(all(colnames(counts(data)) == rownames(phenoData(data))))

# Between-lane distribution of gene-level counts
png(paste0(opt$outdir,"/count_dist.png"))
boxplot(data,las=2,col=cols)
dev.off()

cat("Normalize counts\n")
dataWithin <- withinLaneNormalization(data,"gc", which="full",offset=TRUE) # within-lane gene-specific (and possibly lane-specific) effects, e.g., related to gene length or GC-content
dataNorm <- betweenLaneNormalization(dataWithin, which="full",offset=TRUE) # effects related to between-lane distributional differences, e.g., sequencing depth.

png(paste0(opt$outdir,"/bias_gc_norm.png"))
biasPlot(dataNorm, "gc", log=TRUE, col=cols)
dev.off()

png(paste0(opt$outdir,"/bias_length_norm.png"))
biasPlot(dataNorm, "length", log=TRUE, col=cols)
dev.off()

png(paste0(opt$outdir,"/count_dist_norm.png"))
boxplot(dataNorm,las=2,col=cols)
dev.off()

saveRDS(dataNorm, paste0(opt$outdir,"/EDA_normalised.RData"))

png(paste0(opt$outdir,"/MDP.png"))
MDPlot(data, 1: ncol(data))
dev.off()

for (g in unique(phenoData(data)$Group)) {
  png(paste0(opt$outdir,"/MDP_",g,".png"))
  MDPlot(data, which(phenoData(data)$Group == g))
  dev.off()  
}

# Over-dispersion
png(paste0(opt$outdir,"/varPlot.png"))
meanVarPlot(data, log=TRUE, col=cols)
dev.off()

for (g in unique(phenoData(data)$Group)) {
  png(paste0(opt$outdir,"/varPlot_",g,".png"))
  meanVarPlot(data[phenoData(data)$Group == g,], log=TRUE, col=cols)
  dev.off()
}

# Gene-specific effects on read counts
png(paste0(opt$outdir,"/bias_gc.png"))
biasPlot(data, "gc", log=TRUE, col=cols)
dev.off()

png(paste0(opt$outdir,"/bias_length.png"))
biasPlot(data, "length", log=TRUE, col=cols)
dev.off()

sessionInfo()