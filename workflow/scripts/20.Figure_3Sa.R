# Script to produce figure S3a

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure3Sa.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("EDASeq")
  library("edgeR")
  library("tidyverse")  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--norm"), type="character", default=NULL, 
              help="TMM normalised counts", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$norm)){
  print_help(opt_parser)
  stop("TMM normalised counts (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load counts\n")
counts <- readRDS(opt$norm)

# Get color and shapes
dba <- rgb(226/255,187/255,144/255)
c57 <- rgb(148/255,151/255,152/255)

strains <- as.factor(counts$sample$Strain)
colors <- rainbow(length(levels(strains)))
names(colors) <- levels(strains)
colors["DBA"]<- dba
colors["C57Bl6"]<- c57

points <- as.numeric(gsub("CTRL", 19,  gsub("SD", 17, counts$sample$Treatment)))

cat("Save plot\n")
svg(paste0(opt$outdir,"/FigS3a.svg"))

mds <-plotMDS(counts, col = colors[strains], pch = points, top=100, gene.selection="common", main=paste("ATAC MDS (PCA) using top 100 variable features"))
legend("topleft", legend=c("DBA", "C57Bl6"), pch=19, col=c(dba,c57), ncol=1)
legend("bottomleft", legend=c("SD","CTRL"), pch=c(17,19), col="black", ncol=1)

labs <- data.frame(x=mds$x,y=mds$y,label=strains, col=colors[strains],treatment=points) %>% filter(treatment == 19)  %>% group_by(label,col) %>% summarize(x=mean(x),y=mean(y))

text(labs$x, labs$y, labs$label, cex=1, pos=4, col=labs$col)

dev.off()

sessionInfo()

#DAR PCA  
# library("openxlsx")
# counts <- readRDS('results/7-BXD_normalization/EDA/atac_filtered_normalised_counts.RData')
# diff <-read.xlsx('manuscript/tables/TableS2-Differential_Correlation.xlsx','Differential_Accessibility')
# sigs <- diff[diff$ATAC_FDR < 0.05,'Region_ID']

# strains <- as.factor(counts$sample$Strain)
# colors <- rainbow(length(levels(strains)))
# names(colors) <- levels(strains)
# colors["DBA"]<- dba
# colors["C57Bl6"]<- c57

# points <- as.numeric(gsub("CTRL", 19,  gsub("SD", 17, counts$sample$Treatment)))

# cat("Save plot\n")
# svg("manuscript/figures/FigS3ATAC.svg")

# mds <-plotMDS(counts[sigs,], col = colors[strains], pch = points, top=100, gene.selection="common", main=paste("ATAC MDS (PCA) using top 100 variable features"))
# legend("topleft", legend=c("DBA", "C57Bl6"), pch=19, col=c(dba,c57), ncol=1)
# legend("bottomleft", legend=c("SD","CTRL"), pch=c(17,19), col="black", ncol=1)

# labs <- data.frame(x=mds$x,y=mds$y,label=strains, col=colors[strains],treatment=points) %>% filter(treatment == 19)  %>% group_by(label,col) %>% summarize(x=mean(x),y=mean(y))

# text(labs$x, labs$y, labs$label, cex=1, pos=4, col=labs$col)

# dev.off()