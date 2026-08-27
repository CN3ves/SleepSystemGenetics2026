# Script designed to filter counts from ATAC count files

# Redirect all R logs to Snakemake log
log <- file('logs/7-BXD_normalization/object_7.0.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("edgeR")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="Directory with the count files", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="Processed metadata file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("Directory with the count files (-c) is missing", call.=FALSE)
}
if (is.null(opt$meta)){
  print_help(opt_parser)
  stop("Processed metadata file (-m) is missing", call.=FALSE)
}

cat("Reading inputs\n")
counts <- strsplit(opt$counts, " ")[[1]]
metadata <- read.csv(opt$meta, row.names=1)
metadata$Sample <- sprintf("%03d",metadata$Sample)

# Process count matrix
counts <- readDGE(counts)
# Meta tags detected: __no_feature, __ambiguous, __too_low_aQual, __not_aligned, __alignment_not_unique
colnames(counts$counts) <-  gsub(".*/", "",colnames(counts$counts))
rownames(counts$samples) <-  gsub(".*/", "",rownames(counts$samples))
# Remove Meta tags
counts$counts <- counts$counts[-grep("__", rownames(counts$counts)),]

# Add metadata
counts$samples$Treatment <- metadata[match(gsub('_.*','',rownames(counts$samples)),metadata$Sample),"Treatment"]
counts$samples$Strain <- metadata[match(gsub('_.*','',rownames(counts$samples)),metadata$Sample),"Strain"]
counts$samples$Strain[grep("^[0-9]", counts$samples$Strain)] <- paste0("BXD", counts$samples$Strain[grep("^[0-9]", counts$samples$Strain)])

counts$samples$group <- paste(counts$samples$Strain,counts$samples$Treatment, sep = "_")

saveRDS(counts, paste0(opt$outdir,"/atac_counts.RData"))

sessionInfo()