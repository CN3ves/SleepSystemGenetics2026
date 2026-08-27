# Script designed for enrichment analyses for RNA trasncripts

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
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
  make_option(c("-d", "--deg"), type="character", default=NULL, 
              help="Table with the differential transcriptomic results", metavar="character"),
  make_option(c("-q", "--qtl"), type="character", default=NULL, 
              help="Table with the eQTL results", metavar="character"),
  make_option(c("-i", "--int"), type="character", default=NULL, 
              help="Table with the eQTL interaction results", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="Enchirment database", metavar="character")  ,
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

   
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$deg)){
  print_help(opt_parser)
  stop("Table with the differential transcriptomic results (-p) is missing", call.=FALSE)
}
if (is.null(opt$qtl)){
  print_help(opt_parser)
  stop("Table with the eQTL results (-q) is missing", call.=FALSE)
}
if (is.null(opt$int)){
  print_help(opt_parser)
  stop("Table with the eQTL interaction results (-i) is missing", call.=FALSE)
}
if (is.null(opt$type)){
  print_help(opt_parser)
  stop("Enrichment database (-t) is missing (GO, KEEG, Wikipaths, Reactome)", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}


# Redirect all R logs to Snakemake log
type <- opt$type
log <- file(paste0('logs/11-BXD_enrichment/rna_',type,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")


cat("Reading inputs\n")
deg <- readRDS(opt$deg)$table
qtl <- read.csv(opt$qtl, row.names=1)
int <- read.csv(opt$int, row.names=1)

enrich_obj <- list("DEG" = deg,
                   "QTL" = qtl,
                   "QTLxSD" = int)

for (run in names(enrich_obj)) {
  print(paste("Running enrichment for", run))
  
  if(grepl("DAR", run, ignore.case = TRUE)) {
    params <- list(ora = c("genes", "up", "down"), gsea = TRUE)
  } else { 
    params <- list(ora = c("genes"), gsea = FALSE)
  }
  
  lsts <- get_lists_from_genes(enrich_obj[[run]], universe=keys(TxDb.Mmusculus.UCSC.mm10.knownGene))
  
  run_enrich(lsts, param = params, name=run, type = opt$type, dir=opt$outdir)
}

sessionInfo()