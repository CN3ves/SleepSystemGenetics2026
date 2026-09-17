# Script to produce figure 6b

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure6b.log', open = "wt")
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
              help="Nrf1 stats", metavar="character"),
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

cat("Load interactions statistics\n")
ct <- read.xlsx(opt$stats, sheet="-Interaction- SD x Genotype")
ct$test <- 'FTvsCT'

fv <- read.xlsx(opt$stats, sheet="-Interaction- SD x Tamoxifen")
fv$test <- 'FTvsFV'

deg <- rbind(ct,fv)
deg$logFC <- -deg$logFC

deg$sig <- "No change"
deg$sig[deg$PValue < 0.01] <- "Interaction changes"

id <- deg$ENSEMBL[deg$sig=="Interaction changes"]
id <- table(id)
id <- names(which(id==2))
deg$label <- ""
deg$label[deg$ENSEMBL %in% id] <- deg$SYMBOL[deg$ENSEMBL %in% id]

lim <-  round(max(abs(deg$logFC)),1)

g <- ggplot(data=deg, aes(x=logFC, y=-log10(PValue), fill=sig, col=test, label = label, shape=test)) +
  geom_line(data= deg[deg$label != "",], aes(group=label), col="grey", linetype="dashed") +
  geom_point(size = 3, alpha=0.75) + 
  geom_text_repel(show.legend = FALSE, size=4,max.overlaps=100, col=ifelse(deg$test=='FTvsCT',"darkgreen","darkolivegreen")) + 
  theme_classic() + xlim(-lim,lim) + 
  guides(alpha= "none", col=guide_legend(title="Interaction", override.aes = list(size=5))) + 
  scale_fill_manual(values=c("black", "grey","black")) +
  scale_color_manual(values=c("darkgreen","darkolivegreen")) +
  scale_shape_manual(values=c(24,25)) + 
  geom_hline(yintercept=-log10(0.01), col="grey", alpha = 0.75) + 
  ggtitle("Differential interaction") + xlab("Interaction effect") + 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 20)) 
  
cat("Save figure\n")
svg(paste0(opt$outdir,"/Fig6b.svg"),width=10, height=10)
print(g)
dev.off()

write.csv(deg, paste0(opt$outdir, "/data/6b.csv"))

sessionInfo()
