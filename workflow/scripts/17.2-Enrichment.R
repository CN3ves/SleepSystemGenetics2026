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

source("workflow/scripts/17.0-Enrichment_helper.R")

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-d", "--deg"), type="character", default=NULL, 
              help="Table with the differential transcriptomic results", metavar="character"),
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
log <- file(paste0('logs/17-NRF_differential/enrichment_',type,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")


cat("Reading inputs\n")
deg <- readRDS(opt$deg)

for(run in names(deg)[1:3]) {
  stats <- deg[[run]]$table

  params <- list(ora = c("genes", "up", "down"), gsea = TRUE)

  lsts <- get_lists_from_genes(stats, universe=keys(TxDb.Mmusculus.UCSC.mm10.knownGene))

  print(paste("Running enrichment for", run))
  run_enrich(lsts, param = params, name=run, type = type, dir=opt$outdir)
}

sessionInfo()