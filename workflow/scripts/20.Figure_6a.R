# Script to produce figure 6a

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure6a.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ggplot2")
  library("ggrepel")
  library("openxlsx")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--stats"), type="character", default=NULL, 
              help="SD network", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$stats)){
  print_help(opt_parser)
  stop("Nrf1 stats (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load baseline statistics\n")

ct <- read.xlsx(opt$stats, sheet="-Baseline- ciKO vs Genotype")
ct$test <- 'FTvsCT'

fv <- read.xlsx(opt$stats, sheet="-Baseline- ciKO vs Tamoxifen")
fv$test <- 'FTvsFV'

deg <- rbind(ct,fv)
deg$logFC <- -deg$logFC # invert for plotting because the FT is the reference

deg$sig <- "No change"
deg$sig[deg$FDR < 0.05 & deg$logFC > 0] <- "Up"
deg$sig[deg$FDR < 0.05 & deg$logFC < 0] <- "Down"

shared <- ct[ct$FDR <0.05,"ENSEMBL"] %in% fv[fv$FDR <0.05,"ENSEMBL"]
shared <- ct[ct$FDR <0.05,"ENSEMBL"][shared]

deg$label <- deg$SYMBOL
deg$SYMBOL[!deg$ENSEMBL %in% shared] <- ""

lim <-  round(max(abs(deg$logFC)),1)

g <- ggplot(data=deg, aes(x=logFC, y=-log10(FDR), fill=sig, col=test, label = label, shape=test)) +
  geom_line(data= deg[deg$label != "",], aes(group=label), col="grey", linetype="dashed") +
  geom_point(size = 3, alpha=0.75) + 
  geom_text_repel(show.legend = FALSE, size=4,col=ifelse(deg$test=='FTvsCT',"darkgreen","darkolivegreen")) + 
  theme_classic() + xlim(-lim,lim) + 
  guides(alpha= "none", col=guide_legend(title="DEG", override.aes = list(size=5))) + 
  scale_fill_manual(values=c("navy", "grey","#e2725b")) +
  scale_color_manual(values=c("darkgreen","darkolivegreen")) +
  scale_shape_manual(values=c(24,25)) + 
  geom_hline(yintercept=-log10(0.05), col="grey", alpha = 0.75) + 
  ggtitle("Differential expression") +
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 20)) 
  
cat("Save figure\n")
svg(paste0(opt$outdir,"/Fig6a.svg"),width=10, height=10)
print(g)
dev.off()

write.csv(deg, paste0(opt$outdir, "/data/6a.csv"))

sessionInfo()
