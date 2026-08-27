# Script designed to load sleep phenotype dtaa into an R object


# Redirect all R logs to Snakemake log
log <- file(paste0('logs/9-BXD_qtl/aspects_9.1_pheno.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("edgeR")
  library("ggplot2")
  library("ComplexHeatmap")
  library("tidyr")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="Rdata file with the normalized feature counts", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("Rdata file with the normalized feature counts (-c) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}


cat("Reading inputs\n")
counts <- read.delim(opt$counts, sep = " ", row.names=1)

counts <- list("counts" = t(counts), 
               "sample" = data.frame("Treatment" = NA, 
                                      "Strain" = rownames(counts),
                                      "group" = rownames(counts)))

rownames(counts$sample) <- counts$sample$Strain
saveRDS(counts, paste0(opt$outdir,"/sleep_counts_disp.RData"))

sessionInfo()