# Script designed to convert the sleep phenotyping results to an R object

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/10-BXD_annotation/reduce_10.1.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("foreach")
  library("doParallel")
  library("doSNOW")
  library("rtracklayer")
  library("csaw")
  library("tidyverse")
  library("GenomicRanges")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="Rdata file with feature counts", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="Processed metadata file", metavar="character"),
  make_option(c("-p", "--peaks"), type="character", default=NULL, 
              help="GTF file with the peaks genomic coordinates", metavar="character"),
  make_option(c("-d", "--dar"), type="character", default=NULL, 
              help="Table with the differential statistics", metavar="character"),
  make_option(c("-q", "--qtl"), type="character", default=NULL, 
              help="Table with the QTL results", metavar="character"),
  make_option(c("-i", "--int"), type="character", default=NULL, 
              help="Table with the QTL interaction results", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
   
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("Rdata file with feature counts (-c) is missing", call.=FALSE)
}
if (is.null(opt$meta)){
  print_help(opt_parser)
  stop("Processed metadata file (-m) is missing", call.=FALSE)
}
if (is.null(opt$peaks)){
  print_help(opt_parser)
  stop("GTF file with the peaks genomic coordinates (-p) is missing", call.=FALSE)
}
if (is.null(opt$dar)){
  print_help(opt_parser)
  stop("Table with the differential statistics (-s) is missing", call.=FALSE)
}
if (is.null(opt$qtl)){
  print_help(opt_parser)
  stop("Table with the QTL results (-q) is missing", call.=FALSE)
}
if (is.null(opt$int)){
  print_help(opt_parser)
  stop("Table with the QTL interaction results (-i) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}


cat("Reading inputs\n")
counts <- readRDS(opt$counts)
metadata <- read.csv(opt$meta, row.names=1)
metadata$Sample <- sprintf("%03d",metadata$Sample)
peaks <- read.delim(opt$peaks, header = FALSE, col.names = c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame",  "attribute"))

annotations <- list(dar = read.csv(opt$dar, row.names=1),
                    qtl = read.csv(opt$qtl, row.names=1),
                    int = read.csv(opt$int, row.names=1))

#Get features
cat("Adding average counts to peaks\n")
peaks <- makeGRangesFromDataFrame(peaks, keep.extra.columns=TRUE)

# Add mean count values
SD <- metadata$Sample[metadata$Treatment=="SD"]
meanSD <- counts$counts[,grep('SD', colnames(counts$counts))]
stopifnot(all(sapply(SD, function(x) sum(grepl(paste0('^',x,'_'), colnames(meanSD)))) == 1))
meanSD <- rowMeans(meanSD)

CTRL <-metadata$Sample[metadata$Treatment=="CTRL"]
meanCTRL <- counts$counts[,grep('CTRL', colnames(counts$counts))]
stopifnot(all(sapply(CTRL, function(x) sum(grepl(paste0('^',x,'_'), colnames(meanCTRL)))) == 1))
meanCTRL <- rowMeans(meanCTRL)

elementMetadata(peaks)["SD"] <- meanSD[gsub("Peak_ID (.*);","\\1",peaks$attribute)]
elementMetadata(peaks)["CTRL"] <- meanCTRL[gsub("Peak_ID (.*);","\\1",peaks$attribute)]

# Reduce the peaks
reduce_metadata <- function(peaks) {
  reduced <- GenomicRanges::reduce(peaks)
  
  overlaps <- findOverlaps(peaks, reduced)
  ids <- subjectHits(overlaps)
  combined <- combineTests(ids, as.data.frame(peaks))
  
  elementMetadata(reduced) <- combined

  cl <- makeCluster(20)
  registerDoParallel(cl)
  
  attributes <- foreach(i = 1:length(reduced), .packages = "GenomicRanges") %dopar% {
    
    idx <- subjectHits(overlaps) == i
    attrs <- queryHits(overlaps)[idx]
    attrs <- peaks[attrs]

    peakID <- unique(gsub(".*(chr.*)\\..*","\\1", attrs$attribute))

    if (length(peakID) > 1) {
      ids <- paste(sort(gsub(".*_","",peakID)), collapse = "_")
      peakID <- paste0(unique(gsub("_.*","_",peakID)),ids)
      
    }
    SD <- mean(attrs$SD)
    CTRL <- mean(attrs$CTRL)
    
    if ("snp" %in% names( elementMetadata(attrs))) {# For QTL
     snp <- paste(unique(attrs$snp), collapse = "; ")
     return(c("peakID"=peakID,"SD"=SD,"CTRL"=CTRL, "SNPs" = snp))
    }
    
    return(c("peakID"=peakID,"SD"=SD,"CTRL"=CTRL))
    
  }
  
  stopCluster(cl)
  
  attributes <- as.data.frame(do.call(rbind, attributes))
  elementMetadata(reduced)$PeakID <- attributes$peakID
  elementMetadata(reduced)$SD <- as.numeric(attributes$SD)
  elementMetadata(reduced)$CTRL <- as.numeric(attributes$CTRL)
  names(elementMetadata(reduced)) <- gsub("rep.logFC","logFC", names(elementMetadata(reduced)))
  elementMetadata(reduced)$NAME <- gsub(":.$","",paste(reduced))
  
  return(reduced)
}

for (tab in names(annotations)) {
  cat(paste("Merging", tab, "\n"))
  regions <- peaks
  stats <- annotations[[tab]]
  # sort regions in correct order
  if(tab == 'dar') {
    cat("Match DAR\n")
    regions <- regions[match(rownames(stats), gsub("Peak_ID (.*);","\\1",regions$attribute)),]
  } else if(tab == 'qtl'){
    cat("Match QTL\n")
    regions <- regions[match(stats$lodcolumn, gsub("Peak_ID (.*);","\\1",regions$attribute)),]
  } else if(tab == 'int'){
    cat("Match QTLxFC\n")
    regions <- regions[match(stats$region, gsub("Peak_ID (.*);","\\1",regions$attribute)),]
  }
  seqlevels(regions) <-  seqlevels(regions)[seqlevels(regions) %in% unique(seqnames(regions))]

  # add stats information
  if(tab == 'dar') {
    cat("Set DAR values\n")
    elementMetadata(regions)["logFC"] <- stats[gsub("Peak_ID (.*);","\\1",regions$attribute), "logFC"]
    elementMetadata(regions)["FDR"] <- stats[gsub("Peak_ID (.*);","\\1",regions$attribute), "FDR"]
    elementMetadata(regions)["PValue"] <- stats[gsub("Peak_ID (.*);","\\1",regions$attribute), "PValue"]

    regions <- regions[regions$FDR<0.05,]
  } else {
    cat("Set QTL values\n")
    stopifnot(all(sort(unique(stats$test)) == c("alpha:0.05|-|alpha:0.1","alpha:0.1")))
    stopifnot(all(gsub("Peak_ID (.*);","\\1",regions$attribute) == stats$lodcolumn))
    elementMetadata(regions)["logFC"] <- stats$lod
    elementMetadata(regions)["PValue"] <- ifelse(grepl('0.05',stats$test),0.05,0.1)
    elementMetadata(regions)["test"] <- stats$test
    elementMetadata(regions)["snp"] <- stats$snp
    elementMetadata(regions)["bp"] <- stats$bp
  } 

  cat("Add metadata \n")
  reduced <- reduce_metadata(regions)

  cat("Save results \n")
  export.bed(reduced, paste0(opt$outdir, '/ranges_',tab,'.bed'))  
  saveRDS(reduced, paste0(opt$outdir, '/ranges_',tab,'.RData'))
}

sessionInfo()