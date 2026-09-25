# Script to produce figure S8a

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure8Sa.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("edgeR")
  library("EDASeq")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--counts"), type="character", default=NULL, 
              help="Nrf1 counts", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("Nrf1 counts (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load normalized counts\n")

counts <-  readRDS(opt$counts)

pData(counts)$NRF <- ifelse(gsub("_.*", "", pData(counts)$Group) == "FT", "negative", "positive")
pData(counts)$NRF <- as.factor(pData(counts)$NRF)
pData(counts)$NRF <- relevel(pData(counts)$NRF, ref="positive")
pData(counts)$Treatment <- as.factor(pData(counts)$Treatment)
pData(counts)$Treatment <- relevel(pData(counts)$Treatment, ref="Vehicle")

y <- DGEList(counts=counts(counts), group=pData(counts)$Group)
y$offset <- -offst(counts)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

logcpm <- cpm(y, log=TRUE)

colors <- c("CT_NSD"="grey", "FT_NSD"="red", "FV_NSD"="green", "CT_SD"="darkgrey", "FT_SD"="darkred", "FV_SD"="darkgreen")
pchs <-  c("CT_NSD"= 16, "FT_NSD"=16, "FV_NSD"=16, "CT_SD"=23, "FT_SD"=23, "FV_SD"=23)
  
svg(paste0(opt$outdir,"/FigS8a.svg"),width=10,height=10)
mds <- plotMDS(logcpm, col = colors[pData(counts)$Group], pch = pchs[pData(counts)$Group],main=paste0("PCA using all (", nrow(logcpm),") features"),top=nrow(logcpm), gene.selection="common", dim.plot = c(1,2))

legend("top", legend=names(colors), col=colors, pch=pchs, ncol=3)
dev.off()

tab <- mds$eigen.vectors[,1:2]
rownames(tab) <- pData(counts)$Group
tab <- as.data.frame(tab)
tab$color <- colors[pData(counts)$Group]
tab$shape <- pchs[pData(counts)$Group]

write.csv(tab, paste0(opt$outdir, "/data/S8a.csv"))

sessionInfo()
