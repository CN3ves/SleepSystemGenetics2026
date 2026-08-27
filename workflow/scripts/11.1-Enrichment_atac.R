# Script designed for enrichment analyses for ATAC regions

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ChIPseeker")
  library("clusterProfiler")
  library("ggplot2")
  library("enrichplot")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene") 
  library("org.Mm.eg.db")
  library("tidyverse")
  library("ReactomePA")
  library("optparse")
})

source("workflow/scripts/11.0-Enrichment_helper.R")

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-a", "--annotation"), type="character", default=NULL, 
              help="Rdata containing annotated peaks information", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="Enrichment database", metavar="character"),
  make_option(c("-s", "--stats"), type="character", default=NULL, 
              help="Statistical results from differential analysis", metavar="character"),
  make_option(c("-p", "--peaks"), type="character", default=NULL, 
              help="Genomic region feature information", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

   
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$annotation)){
  print_help(opt_parser)
  stop("Rdata containing annotated peaks information (-p) is missing", call.=FALSE)
}
if (is.null(opt$type)){
  print_help(opt_parser)
  stop("Enrichment database (-t) is missing (GO, KEEG, Wikipaths, Reactome)", call.=FALSE)
}
if (is.null(opt$stats)){
  print_help(opt_parser)
  stop("Statistical results from differential analysis (-s) is missing", call.=FALSE)
}
if (is.null(opt$peaks)){
  print_help(opt_parser)
  stop("Genomic region feature information (-f) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}


# Redirect all R logs to Snakemake log
type <- opt$type
log <- file(paste0('logs/11-BXD_enrichment/atac_',type,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Reading inputs\n")
enrich_obj <- readRDS(opt$annotation)   
enrich_obj <- enrich_obj[-which(names(enrich_obj) %in% c('All', 'Peaks'))] 

#add FC information for GSEA analyses
stats <- readRDS(opt$stats)$table
peaks <- read.delim(opt$peaks, header = FALSE, col.names = c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame",  "attribute"))
peaks <- makeGRangesFromDataFrame(peaks, keep.extra.columns=TRUE)
peaks <- peaks[-as.numeric(grep('_',seqnames(peaks)))]
peaks$ID <-  gsub(';','',gsub('Peak_ID ','',peaks$attribute))

peaks <- peaks[peaks$ID %in% rownames(stats)]
peaks$logFC <- stats[peaks$ID,'logFC']
peaks <- as.GRanges(annotatePeak(peaks, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion=c(-3000, 3000), annoDb="org.Mm.eg.db"))
  
cat(paste("ATAC seq enrichment analyses", type))
for (run in names(enrich_obj)) {
  print(paste("Running enrichment for", run))
  
  if(grepl("DAR", run, ignore.case = TRUE)) {
    params <- list(ora = c("genes", "up", "down"), gsea = TRUE)
  } else { 
    params <- list(ora = c("genes"), gsea = FALSE)
  }
  
  lsts <- get_lists_from_regions(enrich_obj[[run]], peaks)
  
  run_enrich(lsts, param = params, name=run, type = type, dir=opt$outdir)
}

sessionInfo()
