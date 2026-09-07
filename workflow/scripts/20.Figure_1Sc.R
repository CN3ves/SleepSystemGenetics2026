# Script to produce figure S1b

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure1Sc.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ggplot2")
  library("openxlsx")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--S2"), type="character", default=NULL, 
              help="Table S2", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$S2)){
  print_help(opt_parser)
  stop("Table S2 (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load Table S2\n")
atac <- read.xlsx(opt$S2, sheet="Differential_Accessibility")

rna <- read.xlsx(opt$S2, sheet="Differential_Expression")

correlations <- atac[atac$ATAC_FDR < 0.05 & atac$Pearson_correlation_all != '-', 
                     c("Region_ID","ATAC_logFC","Nearest.Gene",
                       "Pearson_correlation_all","Correlation_FDR_all")]

cat("Process correlations\n")
correlations$Group <- 'Correlation not significant'
correlations$Group[as.numeric(correlations$Correlation_FDR_all) < 0.05] <- 'Significant correlation'

rna <- rna[rna$RNA_FDR < 0.05,"Gene_ID"]
  
correlations$RNA <- 'No significant transcript response to SD'
correlations$RNA[correlations$Nearest.Gene %in% rna] <- 'Significant transcript response to SD'

correlations$Direction <- 'ERROR'
correlations$Direction[correlations$ATAC_logFC > 0] <- 'More accessible'
correlations$Direction[correlations$ATAC_logFC < 0 ] <- 'Less accessible'

dens <- density(as.numeric(correlations$Pearson_correlation_all))
mode <- which.max(dens$y)
mode <- dens$x[mode]

g <- ggplot(data = correlations,aes(x = as.numeric(Pearson_correlation_all), fill = Group, alpha = RNA)) +
  geom_histogram(colour="grey", lwd=0.1, binwidth=0.01)+
  geom_density(aes(y=0.01*after_stat(count)),fill='grey', alpha = .05) + 
  theme_classic() + 
  ylab("Gene/region pair counts") + xlab("Pearson Correlation") + 
  scale_alpha_manual(values=c(0.5, 1))  +
  geom_vline(xintercept=mode, col='black', lwd=0.1) + 
  scale_x_continuous(breaks = seq(-1, 1, 0.1))
  
cat("Save plot\n")
svg(paste0(opt$outdir,"/FigS1c.svg"), width=14, height=10)
print(g)
dev.off()

sessionInfo()