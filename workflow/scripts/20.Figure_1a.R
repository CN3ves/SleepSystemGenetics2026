# Script to produce figure 1a

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure1a.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("openxlsx")
  library("ggplot2")
  library("ggrepel")
  library("GenomicRanges")
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
# Volcano plot
atac <- read.xlsx(opt$S2, sheet="Differential_Accessibility")

diff <- makeGRangesFromDataFrame(atac, keep.extra.columns=TRUE)
diff$region <- paste(diff)
diff$Correlation <- as.numeric(diff$Correlation_FDR_all) < 0.05
tab <- as.data.frame(diff)[,c("region", "ATAC_logFC", "ATAC_FDR","Nearest.Gene", "Correlation")]

cat("Set plot parameters\n")
tab$color <- "No significant change"
tab$color[tab$ATAC_FDR < 0.05 & tab$ATAC_logFC > 0] <- "More accessible"
tab$color[tab$ATAC_FDR < 0.05 & tab$ATAC_logFC < 0] <- "Less accessible"

tab$Correlation <- ifelse(tab$Correlation, "Correlated","Not correlated")
tab$Correlation[is.na(tab$Correlation)] <- "No expression"
tab$Correlation[tab$Nearest.Gene == "-"] <- "No annotation"

lim <-  round(max(abs(tab$ATAC_logFC)),1)
  
g_lines <- c(1.2,1.5,2)
g_lines <- c(1/g_lines,g_lines)
g_lines <- log(g_lines)

# Plot with gene annotation
idx_up <- which(tab$ATAC_logFC > 0 )
idx_down <- which(tab$ATAC_logFC < 0)
idx <- unique(c(idx_up[order(tab$ATAC_FDR[idx_up])][1:50],idx_up[order(tab$ATAC_logFC[idx_up])][1:50],
    idx_down[order(tab$ATAC_FDR[idx_down])][1:50],idx_down[order(tab$ATAC_logFC[idx_down])][1:50]))

tab$label <- NA
tab$label[idx] <- tab$Nearest.Gene[idx]

#sample non-significant regions to reduce image size
filter_idx <- tab$color != "No significant change" & -log10(tab$ATAC_FDR) < 5

keep <- c(which(tab$color != "No significant change" & -log10(tab$ATAC_FDR) > 5),
         sample(which(filter_idx),sum(filter_idx)*0.4))

cat("Make plot\n")
g <- ggplot(data=tab[keep,], aes(x=ATAC_logFC, y=-log10(ATAC_FDR), col=color, label = label, pch=Correlation)) +
  geom_point(size = 2) + 
  geom_text_repel(show.legend = FALSE, size=10) + 
  theme_classic() + xlim(-lim,lim) + 
  guides(alpha= "none", col=guide_legend(title="DAR", override.aes = list(size=7)),
         shape=guide_legend(override.aes = list(size=7))) + 
  scale_color_manual(values=c("Less accessible" = "navy", "No significant change" = "grey", "More accessible" = "#e2725b"),name="DAR") +
  scale_shape_manual(values=c("Correlated" = 19, "Not correlated" = 21,"No annotation" = 3,"No expression" = 4),name="Annotated transcript") +
  geom_vline(xintercept=g_lines, col="grey", lty="dashed", alpha = 0.5) +
  geom_hline(yintercept=-log10(0.05), col="grey", alpha = 0.75) + 
  annotate("text", x = min(tab$ATAC_logFC), y = max(-log10(tab$ATAC_FDR), na.rm=TRUE),size=8,
           label = paste0("> Acc.: ", table(tab$color)["More accessible"]), color = "navy", hjust = 0) +
  annotate("text", x = min(tab$ATAC_logFC), y = max(-log10(tab$ATAC_FDR)-1.5, na.rm=TRUE),size=8,
           label = paste0("< Acc.: ", table(tab$color)["Less accessible"]), color = "#e2725b", hjust = 0)  +
  ggtitle("Differentially Accessible Regions") +
  theme(plot.title = element_text(hjust = 0.5,size=30), 
        legend.title=element_text(size=25),
        legend.text=element_text(size=25),
        axis.text=element_text(size=20),
        axis.title=element_text(size=30))

cat("Save plot\n")
svg(paste0(opt$outdir,"/Fig1a.svg"),width=20, height=20)
print(g)
dev.off()

write.csv(tab, paste0(opt$outdir, "/data/1a.csv"))

sessionInfo()