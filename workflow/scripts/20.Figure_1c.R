# Script to produce figure 1c

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure1c.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("clusterProfiler")
  library("ggplot2")
  library("enrichplot")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--enrich"), type="character", default=NULL, 
              help="Enrichment information", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$enrich)){
  print_help(opt_parser)
  stop("Enrichment information (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load Table S3\n")
enrichment <- readRDS(opt$enrich)

cat("Get ATAC enrichment\n")
go <- enrichment[["ATAC_GO"]]
go@result$Description <- paste0(go@result$Description, " (GO-", go@result$ONTOLOGY, ")")

kegg <- enrichment[["ATAC_KEGG"]]
kegg@result$Description <- paste(kegg@result$Description, "(KEGG)")
kegg@result$ONTOLOGY <- NA
kegg@result <- kegg@result[, names(go@result)]

atac <- go
atac@result <- rbind(atac@result,kegg@result)
atac@result <- atac@result[order(atac@result$qvalue),]

cat("Get RNA enrichment\n")
go <- enrichment[["RNA_GO"]]
go@result$Description <- paste0(go@result$Description, " (GO-", go@result$ONTOLOGY, ")")

kegg <- enrichment[["RNA_KEGG"]]
kegg@result$Description <- paste(kegg@result$Description, "(KEGG)")
kegg@result$ONTOLOGY <- NA
kegg@result <- kegg@result[, names(go@result)]

rna <- go
rna@result <- rbind(rna@result,kegg@result)
rna@result <- rna@result[order(rna@result$qvalue),]


cat("Filter overlapping terms\n")

all <- atac
all@result <- all@result[all@result$ID %in% rna@result$ID,]


termsim <- pairwise_termsim(all)
p <- enrichplot::emapplot(termsim, showCategory =  200) + theme(text = element_text(size = 8))


cat("Save plot\n")
svg(paste0(opt$outdir,"/Fig1c.svg"),height=14,width=14)
print(p)
dev.off()

sessionInfo()