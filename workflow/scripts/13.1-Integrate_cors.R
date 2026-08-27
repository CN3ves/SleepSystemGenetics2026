# Script designed to calculate correlated between genomic regions and the transcripts they are annotated to

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/13-BXD_integrate/cors_13.1.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("edgeR")
  library("GenomicRanges")
  library("ChIPseeker")
  library("org.Mm.eg.db")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
  library("foreach")
  library("doParallel")
  library("doSNOW")
  library("ggplot2")
  library("optparse")
})

source("workflow/scripts/13.0-Integrate_helper.R")

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-a", "--atac_counts"), type="character", default=NULL, 
              help="Normalised and filteres ATAC count R object", metavar="character"),
  make_option(c("-r", "--rna_counts"), type="character", default=NULL, 
              help="Normalised and filteres RNA count R object", metavar="character"),
  make_option(c("-c", "--atac_diff"), type="character", default=NULL, 
              help="Differential results for ATAC", metavar="character"),
  make_option(c("-t", "--rna_diff"), type="character", default=NULL, 
              help="Differential results for RNA", metavar="character"),
  make_option(c("-p", "--regions"), type="character", default=NULL, 
              help="Genomic region list", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$atac_counts)){
  print_help(opt_parser)
  stop("Normalised and filteres ATAC count R object (-a) is missing", call.=FALSE)
}
if (is.null(opt$rna_counts)){
  print_help(opt_parser)
  stop("Normalised and filteres RNA count R object (-r) is missing", call.=FALSE)
}
if (is.null(opt$atac_diff)){
  print_help(opt_parser)
  stop("Differential results for ATAC (-c) is missing", call.=FALSE)
}
if (is.null(opt$rna_diff)){
  print_help(opt_parser)
  stop("Differential results for RNA (-t) is missing", call.=FALSE)
}
if (is.null(opt$regions)){
  print_help(opt_parser)
  stop("Genomic region list (-p) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading files\n")
atac_cpm <- collapse_counts(readRDS(opt$atac_counts))
rna_cpm <- collapse_counts(readRDS(opt$rna_counts))

atac_diff <- read.csv(opt$atac_diff)
rna_diff <- read.csv(opt$rna_diff)

cat("Filter data for shared samples\n")
rna_cpm <- as.data.frame(rna_cpm[,colnames(atac_cpm)])
atac_cpm <- as.data.frame(atac_cpm[,colnames(rna_cpm)])

cat("Filter low count number features\n")
atac_sig <- atac_cpm[rownames(atac_cpm) %in% atac_diff$X,]
rna_sig <- rna_cpm[rownames(rna_cpm) %in% rna_diff$X,]
  
cat("Annotate ATAC regions\n")
regions <- read.delim(opt$regions,header= FALSE)
regions <- GRanges(seqnames=regions$V1, 
                   IRanges(start=as.numeric(regions$V4), end=as.numeric(regions$V5)),
                   PeakID = gsub("Peak_ID (.*);","\\1",regions$V9))

regions <- annotatePeak(regions, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion=c(-3000, 3000), annoDb="org.Mm.eg.db")
regions <- as.GRanges(regions)
atac_sig$SYMBOL <-  regions[match(rownames(atac_sig), regions$PeakID)]$SYMBOL
  
cat("Filter data for matching genes\n")
atac_sig <- atac_sig[atac_sig$SYMBOL %in% rownames(rna_sig),]
rna_sig <- rna_sig[unique(atac_sig$SYMBOL),]


cat("Calculate correlations\n")
correlations <- list()

cl <- parallel::makeCluster(40)
doParallel::registerDoParallel(cl)
paths <- .libPaths()
clusterExport(cl, "paths")
clusterEvalQ(cl, .libPaths(paths))

correlations <-foreach(i = 1:nrow(atac_sig), .packages = c("doParallel","GenomicRanges")) %dopar% {

  atac_all <- atac_sig[i,c('SYMBOL',names(atac_sig)[-ncol(atac_sig)])]
  atac_ctrl <- atac_sig[i,c('SYMBOL',grep('CTRL',names(atac_sig),value=TRUE))]
  atac_sd <- atac_sig[i,c('SYMBOL',grep('SD',names(atac_sig),value=TRUE))]
  
  rna_all <- rna_sig[atac_all$SYMBOL,names(atac_all)[-1]]
  rna_ctrl <- rna_sig[atac_ctrl$SYMBOL,names(atac_ctrl)[-1]]
  rna_sd <- rna_sig[atac_sd$SYMBOL,names(atac_sd)[-1]]
  
  cor_all <- cor.test(as.numeric(rna_all),as.numeric(atac_all[colnames(rna_all)]))
  cor_ctrl <- cor.test(as.numeric(rna_ctrl),as.numeric(atac_ctrl[colnames(rna_ctrl)]))
  cor_sd <- cor.test(as.numeric(rna_sd),as.numeric(atac_sd[colnames(rna_sd)]))
  
  data.frame(gene=rownames(rna_all), regionID=rownames(atac_all), 
             region=paste(regions[regions$PeakID %in% rownames(atac_all)]),
             cor_ctrl=cor_ctrl$estimate, p_ctrl=cor_ctrl$p.value,
             cor_sd=cor_sd$estimate, p_sd=cor_sd$p.value,
             cor_all=cor_all$estimate, p_all=cor_all$p.value, row.names = NULL)
}

parallel::stopCluster(cl)
saveRDS(correlations,paste0(opt$outdir,'/cors_annotation.Rdata'))

correlations <- do.call(rbind,correlations)
correlations$padj_all <- p.adjust(correlations$p_all, method = "BH")
correlations$padj_ctrl <- p.adjust(correlations$p_ctrl, method = "BH")
correlations$padj_sd <- p.adjust(correlations$p_sd, method = "BH")
correlations <- correlations[order(correlations$padj_all),]
  
write.csv(correlations, paste0(opt$outdir,"/gene_region_cors_annotation.csv"))

sessionInfo()