# Script designed to calculate correlated between genomic regions and the transcripts with interaction QTLs

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/13-BXD_integrate/qtl_13.2.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("edgeR")
  library("GenomicRanges")
  library("tidyverse")
  library("dplyr")
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
  make_option(c("-c", "--atac_qtl"), type="character", default=NULL, 
              help="QTL t-test interactions screening for ATAC", metavar="character"),
  make_option(c("-t", "--rna_qtl"), type="character", default=NULL, 
              help="QTL t-test interactions screening for RNA", metavar="character"),
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
if (is.null(opt$atac_qtl)){
  print_help(opt_parser)
  stop("QTL t-test interactions creening for ATAC(-c) is missing", call.=FALSE)
}
if (is.null(opt$rna_qtl)){
  print_help(opt_parser)
  stop("QTL t-test interactions screening for RNA (-t) is missing", call.=FALSE)
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

atac_inverted <- classify_qtl(read.csv(opt$atac_qtl, row.names=1))
rna_inverted <- classify_qtl(read.csv(opt$rna_qtl, row.names=1))

# get correlations
regions <- read.delim(opt$regions, header =FALSE)
regions <- GRanges(seqnames=regions$V1, 
                   IRanges(start=as.numeric(regions$V4), end=as.numeric(regions$V5)),
                   region = gsub("Peak_ID (.*);","\\1",regions$V9))

atac_inverted <- match2coords(atac_inverted, regions, filter=TRUE) 
rna_inverted <- rna_inverted[rna_inverted$type != 'rm',]

cat("Find overlapping qtls\n")
ov <- lapply(1:nrow(rna_inverted), function(i) {
  hits <- which(atac_inverted$snp == rna_inverted[i,"snp"])
  if(length(hits) ==0 ) hits=NA
  data.frame(queryHits=i, subjectHits=hits)
  })
ov <- do.call(rbind,ov)
ov <- ov[!is.na(ov$subjectHits),]
ov <- Hits(ov$queryHits, ov$subjectHits, length(atac_inverted), length(atac_inverted), sort.by.query=TRUE)

cat("Calculate correlations\n")
correlations <- list()
for (i in 1:length(ov)) {
  message('\r', paste0(i, " out of ", length(ov),': ', round(i/length(ov)* 100,5),"%          "), appendLF = FALSE)

  rna <- queryHits(ov)[i]
  atac <- subjectHits(ov)[i]
  
  rna <- rna_inverted[rna,]
  atac <- atac_inverted[atac,]
  
  rna_snp <- rna$snp
  atac_snp <- atac$snp
  
  rna <- rna_cpm[rna$region, colnames(rna_cpm) %in% colnames(atac_cpm)]
  atac <- atac_cpm[atac$region, colnames(rna)]
      
  correlation <- cor.test(as.numeric(rna),as.numeric(atac))
  
  region <- regions$region == rownames(atac)
  region <- paste0(regions[region])
  
  correlations[[i]] <- data.frame(gene=rownames(rna),
                    regionID=rownames(atac),
                    region=region,
                    cor=correlation$estimate, 
                    p=correlation$p.value, 
                    rna_snp=rna_snp, 
                    atac_snp = atac_snp, 
                    same_snp = rna_snp==atac_snp,
                    row.names = NULL)
  
}
cat("Save results\n")
saveRDS(correlations,paste0(opt$outdir,"/gene_region_qtl.RData"))

correlations <- do.call(rbind,correlations)
correlations$padj <- p.adjust(correlations$p, method = "BH")
correlations <- correlations[order(correlations$padj),]
  
write.csv(correlations, paste0(opt$outdir,"/gene_region_cors_qtl.csv"))

sessionInfo()