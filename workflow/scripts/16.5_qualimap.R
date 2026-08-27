# Script designed to filter counts from ATAC count files

# Redirect all R logs to Snakemake log
log <- file('logs/16-NRF_bam/qualimap_16.5.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-b", "--bams"), type="character", default=NULL, 
              help="Alignment BAM files", metavar="character"),
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="Gene-level count files", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="Processed metadata file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bams)){
  print_help(opt_parser)
  stop("Alignment BAM files (-b) is missing", call.=FALSE)
}
if (is.null(opt$counts)){
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
# Counts
counts_files <- strsplit(opt$counts, " ")[[1]]
bam_files <-  strsplit(opt$bams, " ")[[1]] 
metadata <- read.csv(opt$meta, row.names=1)

# Process count matrix
counts_files <- strsplit(opt$count, " ")[[1]]
counts <- lapply(counts_files, function(file) {
  df <- read.delim(file, header=FALSE)[,1:2]
  colnames(df) <-  c('ID', gsub(".*/(.*)_ReadsPerGene.*", "\\1",file))
  df <- df[-grep("N_", df$ID),] # Remove Meta tags
  return(df)
})
counts <- Reduce(function(x, y) merge(x, y, by="ID"), counts)
rownames(counts) <- counts$ID
counts <- counts[,-1]

metadata$Group <- gsub(".*_","",metadata$Group)
metadata$file <- paste0(opt$outdir,"/qualimap_counts.txt")
metadata$index <- match(metadata$SampleName,colnames(counts))+1
metadata <-  metadata[,c("SampleName","Group", "file", "index")]


write.table(counts,  unique(metadata$file), sep="\t", quote = FALSE, col.names=FALSE)

### DO IT FOR CONDITION BECAUSE CORRELATION PLOT IS CRAP
write.table(metadata[grep("FV", metadata$SampleName),], paste0(opt$outdir,"/qualimap_FVsamples.txt"), sep="\t", quote = FALSE,row.names = FALSE, col.names=FALSE)

write.table(metadata[grep("FT", metadata$SampleName),], paste0(opt$outdir,"/qualimap_FTsamples.txt"), sep="\t", quote = FALSE,row.names = FALSE, col.names=FALSE)

write.table(metadata[grep("CT", metadata$SampleName),], paste0(opt$outdir,"/qualimap_CTsamples.txt"), sep="\t", quote = FALSE,row.names = FALSE, col.names=FALSE)

write.table(metadata, paste0(opt$outdir,"/qualimap_samples.txt"), sep="\t", quote = FALSE,row.names = FALSE, col.names=FALSE)

bam_samples <-data.frame("path"= bam_files,
                         "Sample"= gsub(".co.bam","",basename(bam_files)))
  
write.table(bam_samples, paste0(opt$outdir,"/qualimap_multi.txt"), sep="\t", quote = FALSE,row.names = FALSE, col.names=FALSE)

sessionInfo()