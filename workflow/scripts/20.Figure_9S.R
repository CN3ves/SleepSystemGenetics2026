# Script to produce figure S9

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure9S.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("EDASeq")
  library("edgeR")
  library("tidyverse")
  library("ggplot2")
  library("ggrepel")
  library("ggpubr")
  library("org.Mm.eg.db")
  library("igraph")
  library("openxlsx")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--rna"), type="character", default=NULL, 
              help="RNA-seq counts", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$rna)){
  print_help(opt_parser)
  stop("RNA-seq counts(-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

genes <- c("Ncoa5",  "Dlgap4",  "Arc",  "Grin1", "Gmds",  "Gabrb1", "Tfam",  "Nrf1", "Htr1b")
counts <-  readRDS(opt$rna)

y <- DGEList(counts=counts(counts), group=pData(counts)$Group)
y$offset <- -offst(counts)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

logcpm <- cpm(y, log=TRUE)

cpm <- logcpm %>% as.data.frame %>% mutate(gene = rownames(logcpm)) %>% pivot_longer(!gene, values_to="logcpm", names_to = "Sample") 
cpm$Group <- pData(counts)[cpm$Sample,"Group"]
cpm$SYMBOL <- mapIds(org.Mm.eg.db, gsub("\\..*","",cpm$gene),  "SYMBOL","ENSEMBL")
cpm <- cpm[cpm$SYMBOL %in% genes,]
cpm$Treatment <- gsub(".*_","",cpm$Group)
cpm$Group <- gsub("_.*","",cpm$Group)

cpm$test <- paste(cpm$Group ,cpm$Treatment)

comps <- list(c("CT NSD","CT SD"),
              c("FT NSD","FT SD"),
              c("FV NSD","FV SD"),
              c("CT NSD","FT NSD"),
              c("CT NSD","FV NSD"),
              c("FT NSD","FV NSD"))

g <- ggplot(cpm, aes(x=test, y=logcpm, fill=Treatment,shape=Group)) + 
    geom_boxplot(position=position_dodge(1), outlier.shape = NA) +
    geom_point(aes(fill = Treatment), size=1, position = position_jitterdodge()) + 
    labs(title="",x="Genotype", y = "RPKM") +
    theme_classic(base_size = 10) + stat_compare_means(comparisons =comps, method = "t.test") +
    theme(plot.title=element_text(hjust=0.5),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(~SYMBOL, scales = "free_y")

svg(paste0(opt$outdir, "/FigS9.svg"),width=10,height=10)
print(g)
dev.off()


write.csv(cpm, paste0(opt$outdir, "/data/S9.csv"))

sessionInfo()