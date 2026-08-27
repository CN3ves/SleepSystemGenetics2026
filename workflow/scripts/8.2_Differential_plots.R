# Script designed for plotting differential analyses results

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("edgeR")
  library("ggplot2")
  library("ComplexHeatmap")
  library("tidyr")
  library("optparse")
})

source("workflow/scripts/8.0_Differential_helper.R")
source("workflow/scripts/7.3_Normalization_helper.R")

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="RDS counts object", metavar="character"),
  make_option(c("-s", "--stats"), type="character", default=NULL, 
              help="Rdata file with the statistical results", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("RDS counts object (-c) is missing", call.=FALSE)
}
if (is.null(opt$stats)){
  print_help(opt_parser)
  stop("Rdata file with the statistical results (-s) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

exp <- gsub('_stat.*','',basename(opt$stats))

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/8-BXD_differential/plots_8.2_', exp,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Reading inputs\n")
original_fit <- readRDS(opt$stats)
counts <- readRDS(opt$counts)
plot_settings <- plot_helper(counts$sample$Strain, counts$sample$Treatment)

# Label significant features
tbl_original <- original_fit$table
coef_original <- as.data.frame(original_fit$coef)
sig_original <-  rownames(tbl_original)[tbl_original$FDR < 0.05]

coef_original$sig <- "Not sig"
coef_original[sig_original,"sig"] <- "Sig"

# Plot model coefficients
coef_box(coef_original, title="Coefficient distribution vs CTRL", file=paste0(opt$outdir, "/", exp, "_coeff_box.png"))

sigs <- rownames(coef_original)[coef_original$sig == "Sig"]
filter <- rownames(tbl_original[sigs,])[order(tbl_original[sigs,]$FDR)][1:min(2000, length(sigs))]
coef_heatmap(coef_original, filter, "BXD line effect (2000 most significant)", paste0(opt$outdir, "/", exp,"_coeff_heatFDR.png"))

filter <- rownames(tbl_original[sigs,])[order(abs(tbl_original[sigs,]$logFC), decreasing = TRUE)][1:min(2000, length(sigs))]
coef_heatmap(coef_original, filter, "BXD line effect (2000 highest effect)", paste0(opt$outdir, "/", exp, "_coeff_heatFC.png"))

#Plot MDS for significant features
# Plot a nice PCA with variance explained 
for (i in c(100,500,1000)) {
  msd <- plotMDS(counts[sig_original,], top=i, gene.selection="common", dim.plot = c(1,2), ndim=4,plot=FALSE)

  pca <- prcomp(as.dist(msd$distance.matrix))
  var <- summary(pca)$importance["Proportion of Variance",1:2]
  var <- paste0(round(var * 100,2), "%")

  png(paste0(opt$outdir, "/", exp, "_sig_PC12.png"))
  plot(pca$rotation[,"PC1"],  pca$rotation[,"PC2"], col = plot_settings[["colors"]], pch =  plot_settings[["shapes"]], 
    xlab = paste0("PC1 (", var[1], ")"), ylab = paste0("PC2 (", var[2], ")"), main=paste("MDS (PCA) using top variable", i, "features"))
  legend("topright", legend=c("DBA", "C57Bl6"), pch=19, col=c(plot_settings[["parents"]][['DBA']], plot_settings[["parents"]][['C57Bl6']]), ncol=1)
  legend("bottomright", legend=c("SD","CTRL"), pch=c(17,19), col="black", ncol=1)
  dev.off()
}

sessionInfo()