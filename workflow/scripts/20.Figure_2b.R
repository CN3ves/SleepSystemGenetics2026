# Script to produce figure 2b

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure2b.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ggplot2")
  library("ggrepel")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--footprints"), type="character", default=NULL, 
              help="Footprint statistics", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$footprints)){
  print_help(opt_parser)
  stop("Footprint statistics (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load footprints\n")
meta_footprints <- read.csv(opt$footprints)

# Make aggregate plot
df <- meta_footprints[order(meta_footprints$mean_activity),]
df$tf <- rownames(df)
df$Sig <- "Not sig."
df$Sig[df$padj < 0.05] <- "Sig."
df$Sig[df$Sig == "Sig." &  df$tf %in% names(which(tf_sig >= 30))] <- "Sig. 30+ lines"
df$Sig[df$tf %in% names(which(tf_sig == n))] <- "Sig. all lines"
cols <- c("Not sig." = "grey", "Sig." = "darkred", "Sig. 30+ lines" = "darkgreen","Sig. all lines" = "limegreen")

df$labs <- toupper(df$tf)
df$labs[df$Sig == "Not sig."] <- ""
df[names(which(tf_sig[rownames(df)] < 5)),"labs"] <- ""
df$labs <- gsub("MA[0-9]*\\.[0-9]\\.", "", df$labs)
df$labs[df$labs != ""] <- paste0(df$labs[df$labs != ""], " (", tf_sig[df$tf[df$labs != ""]],")")
df$labs[grep("NA",df$labs)] <- "" 

df$tf <- factor(df$tf, levels=sample(df$tf))

y_lims <- round(max(abs(df$mean_activity)),2)+0.01

g <- ggplot(df, aes(x=tf,y=mean_activity, fill=Sig, col = Sig, label = labs)) +
  geom_point(shape = 21, size=1) +
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) + 
  theme_classic() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title=element_text(size=14,face="bold")) +
  geom_text_repel(size=3,show.legend = FALSE,max.overlaps=20) +
  ylab("Mean activity score\nCTRL <-> SD") + xlab("") + 
  ylim(c(-y_lims,y_lims)) + 
  geom_hline(yintercept=0, color = "black", linewidth=0.5) 

cat("Save plot\n")
svg(paste0(opt$outdir,"/Fig2b.svg"), width=14, height=7)
print(g)
dev.off()

sessionInfo()