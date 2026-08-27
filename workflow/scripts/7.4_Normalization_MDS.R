# Script designed to compute multridimentional scaling 

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("EDASeq")
  library("edgeR")
  library("optparse")
})
source("workflow/scripts/7.3_Normalization_helper.R")

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="RDS counts object", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("RDS counts object (-c) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

exp <- gsub('_filter.*','',basename(opt$counts))

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/7-BXD_normalization/eda_7.4_', exp,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Reading inputs\n")
counts <- readRDS(opt$counts)
plot_settings <- plot_helper(counts$sample$Strain, counts$sample$Treatment)

# TMM Normalization
counts <- calcNormFactors(counts)
norm_counts <- sapply(1:ncol(counts$counts), function(x) counts$counts[,x] /(counts$samples$lib.size[x] * counts$samples$norm.factors[x]) )

png(paste0(opt$outdir, "/", exp, "_RLE_TMM.png"))
plotRLE(norm_counts, outline=FALSE, col=plot_settings[["colors"]], main="RLE plot with TMM normalization")
dev.off()

saveRDS(counts, paste0(opt$outdir, "/",exp,"_filtered_normalised_counts.RData"))

# Screen MDS plots
for (i in c(100,500,1000,5000,10000,15000,20000,50000)) {

  png(paste0(opt$outdir, "/", exp, "_QC_PCA",i,".png"))
  plotMDS(counts, col = plot_settings[["colors"]], pch = plot_settings[["shapes"]], top=i, gene.selection="common", main=paste("MDS (PCA) using top",i, exp, "variable features"))
  legend("topright", legend=c("DBA", "C57Bl6"), pch=19, col=c(plot_settings[["parents"]][['DBA']], plot_settings[["parents"]][['C57Bl6']]), ncol=1)
  legend("bottomright", legend=c("SD","CTRL"), pch=c(17,19), col="black", ncol=1)
  dev.off()
  
  png(paste0(opt$outdir, "/", exp,"_QC_PCoA",i,".png"))
  plotMDS(counts, col = plot_settings[["colors"]], pch = plot_settings[["shapes"]], top=i, gene.selection="pairwise", main=paste("MDS (PCoA) using top",i, exp, "variable features"))
  legend("topright", legend=c("DBA","C57Bl6") , pch=19, col=c(plot_settings[["parents"]][['DBA']], plot_settings[["parents"]][['C57Bl6']]), ncol=1)
  legend("bottomright", legend=c("SD","CTRL"), pch=c(17,19), col="black", ncol=1)
  dev.off()
}

# Plot a nice PCA with variance explained (top 100)
msd <- plotMDS(counts, top=100, gene.selection="common", dim.plot = c(1,2), ndim=4,plot=FALSE)

pca <- prcomp(as.dist(msd$distance.matrix))
var <- summary(pca)$importance["Proportion of Variance",1:4]
var <- paste0(round(var * 100,2), "%")

png(paste0(opt$outdir, "/", exp, "_QC_PCA_100_PC12.png"))
plot(pca$rotation[,"PC1"],  pca$rotation[,"PC2"], col = plot_settings[["colors"]], pch =  plot_settings[["shapes"]], xlab = paste0("PC1 (", var[1], ")"), ylab = paste0("PC2 (", var[2], ")"), main="MDS (PCA) using top variable 100 features")
legend("topright", legend=c("DBA", "C57Bl6"), pch=19, col=c(plot_settings[["parents"]][['DBA']], plot_settings[["parents"]][['C57Bl6']]), ncol=1)
legend("bottomright", legend=c("SD","CTRL"), pch=c(17,19), col="black", ncol=1)
dev.off()

png(paste0(opt$outdir, "/", exp,"_QC_PCA_100_PC34.png"))
plot(pca$rotation[,"PC3"],  pca$rottion[,"PC4"], col = plot_settings[["colors"]], pch =  plot_settings[["shapes"]], xlab = paste0("PC3 (", var[3], ")"), ylab = paste0("PC4 (", var[4], ")"), main="MDS (PCA) using top variable 100 features")
legend("topright", legend=c("DBA", "C57Bl6"), pch=19, col=c(plot_settings[["parents"]][['DBA']], plot_settings[["parents"]][['C57Bl6']]), ncol=1)
legend("bottomright", legend=c("SD","CTRL"), pch=c(17,19), col="black", ncol=1)
dev.off()

# Alternative plots for exploration
# Focus plot on Treatment
for (i in c(5000,10000,20000)) {
  png(paste0(opt$outdir, "/", exp, "_QC_PCA_",i,"treat.png"))
  msd <- plotMDS(counts, col = as.factor( plot_settings[["shapes"]]), pch =  plot_settings[["shapes"]], top=i, gene.selection="common", dim.plot = c(1,2), main=paste("MDS (PCA) using top",i, "variable features"))
  legend("topright", legend=c("SD","CTRL"), pch=c(17,19), col= as.factor( plot_settings[["shapes"]]), ncol=1)
  idx <-  as.factor(counts$samples$Strain)  %in% c(names(plot_settings[["parents"]]),sample(unique(counts$samples$Strain),3))
  cols <- gsub("17", "darkgreen", gsub("19","darkorange",  plot_settings[["shapes"]][idx]))
  labels <- gsub("BXD", "", gsub("C57.*", "C57",counts$samples$Strain[idx]))
  text(msd$x[idx], msd$y[idx], labels= labels, cex= 0.7, pos=2, col= cols)
  dev.off()
}

# Color by strains
sample_lists <- list()
if(exp == "ATAC") {
  sample_lists <- list("0" = c("BXD49", "BXD48", "BXD43", "BXD100"))
}
samples <- sample(unique(counts$samples$Strain))
samples <- samples[!samples %in% sample_lists[["0"]]]

# Group strains for easier coloring
for(i in 1:(length(samples)/5)) {
  s <- (i-1)*5+1
  sample_lists[[as.character(i)]] <- samples[s:(s+4)]
}

# plot for each strain group
for(s in 1:length(sample_lists)) {
  strains <- sample_lists[[s]]
  
  colors <- rep("grey", length(counts$samples$group))
  cols <- c("red", "green", "blue", "orange", "cyan")[1:length(strains)]
  names(cols) <- strains
  idx <- counts$samples$Strain %in% strains
  colors[idx] <- cols[counts$samples$Strain[idx]]

  png(paste0(opt$outdir, "/", exp,"_QC_PCA_strains_",s,".png"))
  msd <- plotMDS(counts, col = colors, pch = 1, top=100, gene.selection="common")
  legend("topright", legend=strains, pch=1, col=cols, ncol=1)
  text(msd$x[idx], msd$y[idx], labels=rownames(counts$samples)[idx], cex= 0.7, pos=2, col = paste0("dark",colors[idx]))
  dev.off()
  
  if (s==1 & exp == "ATAC") { 
    # First group has the mixed samples. Plot with original metadata
    colors[idx] <- cols[original_metadata$Strain[match(original_metadata$Sample,rownames(counts$samples))][idx]]
    
    png(paste0(opt$outdir, "/", exp,"_QC_PCA_strains_",s,"_original.png"))
    msd <- plotMDS(counts, col = colors, pch = 1, top=100, gene.selection="common")
    legend("topright", legend=strains, pch=1, col=cols, ncol=1)
    text(msd$x[idx], msd$y[idx], labels=rownames(counts$samples)[idx], cex= 0.7, pos=2, col = paste0("dark",colors[idx]))
    dev.off()
    }
}

sessionInfo()