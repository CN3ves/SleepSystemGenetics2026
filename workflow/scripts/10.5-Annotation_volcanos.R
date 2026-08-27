# Script designed to get description and QC of significant regions

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/10-BXD_annotation/volcanos_10.5.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ggplot2")
  library("ggrepel")
  library("GenomicRanges")
  library("ChIPseeker")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene") 
  library("org.Mm.eg.db")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-a", "--annot"), type="character", default=NULL, 
              help="ATAC seq results annotation file", metavar="character"),
  make_option(c("-r", "--rna"), type="character", default=NULL, 
              help="RNA seq differential results table", metavar="character"),
  make_option(c("-q", "--qtl"), type="character", default=NULL, 
              help="RNA seq qtl results table", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

   
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$annot)){
  print_help(opt_parser)
  stop("ATAC seq results annotation file (-a) is missing", call.=FALSE)
}
if (is.null(opt$rna)){
  print_help(opt_parser)
  stop("RNA seq differential results table (-d) is missing", call.=FALSE)
}
if (is.null(opt$qtl)){
  print_help(opt_parser)
  stop("RNA seq qtl results table (-q) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

#helper function
volcano <- function(data, name) {

  colors <- c('Down'="navy", "No change"="grey",  "Up"="#e2725b")
  colors <- colors[unique(data$color)]
  g <- ggplot(data=data, aes(x=logFC, y=-log10(FDR), col=color, label = label, alpha=alpha)) +
      geom_point(size = 0.2) + 
    geom_text_repel(show.legend = FALSE, size=3) + 
    theme_classic() + xlim(-lim,lim) + 
      guides(alpha= "none", col=guide_legend(title="DAR", override.aes = list(size=2))) + 
    scale_color_manual(values=colors) +
      geom_vline(xintercept=g_lines, col="grey", lty="dashed", alpha = 0.5) +
    geom_hline(yintercept=-log10(0.05), col="grey", alpha = 0.75) + 
      ggtitle("Differential accessibility") +
    theme(plot.title = element_text(hjust = 0.5))

  png(name)
  print(g)
  dev.off()
}

cat("Reading inputs\n")
for (exp in c('ATAC', 'RNA')) {
  cat(paste("Plotting volcanos for",exp,"\n"))
  if( exp == 'ATAC') {
    atac <- readRDS(opt$annot)
    dar <- atac[['DAR']]
    tab <- as.data.frame(dar)[,c("NAME", "logFC", "FDR","SYMBOL")]
    qtl <- atac[['QTL']]

    ov <- findOverlaps(dar,qtl)
    qtl <- queryHits(ov)
    range <- c(10,5,25)
    }
  if( exp == 'RNA') {
    deg <- readRDS(opt$rna)$table 
    tab <- as.data.frame(deg)[,c('logFC', 'FDR')] 
    tab$SYMBOL <- rownames(tab)
  
    qtl <- read.csv(opt$qtl)
    qtl <- tab$SYMBOL %in% qtl$lodcolumn
    range <- c(50,50,50)
  }

  # Set volcano parameters
  tab$color <- "No change"
  tab$color[tab$FDR < 0.05 & tab$logFC > 0] <- "Up"
  tab$color[tab$FDR < 0.05 & tab$logFC < 0] <- "Down"

  lim <-  round(max(abs(tab$logFC)),1)
  
  g_lines <- c(1.2,1.5,2)
  g_lines <- c(1/g_lines,g_lines)
  g_lines <- log(g_lines)

  # Plot overlap with QTL and gene annotation
  tab$alpha <- 0.25
  tab$alpha[qtl] <- 1

  idx_up <- which(tab$logFC > 0 & tab$alpha==1)
  idx_down <- which(tab$logFC < 0 & tab$alpha==1)
  idx <- c(idx_up[order(tab$FDR[idx_up])[1:range[1]]],
         idx_down[order(tab$FDR[idx_down])[1:range[2]]])

  tab$label <- NA
  tab$label[idx] <- tab$SYMBOL[idx]

  volcano(tab, paste0(opt$outdir,'/',exp,"_volcano_QTLgenes.png"))

  # Plot overlap with QTL and region name
  if( exp == 'ATAC') {  
    tab$label <- NA
    tab$label[idx] <- tab$NAME[idx]
  
    volcano(tab,  paste0(opt$outdir,"/ATAC_volcano_QTLregions.png"))
  }

  # Plot with gene annotation
  idx_up <- which(tab$logFC > 0 )
  idx_down <- which(tab$logFC < 0)
  idx <- c(idx_up[order(tab$FDR[idx_up])[1:range[3]]],
          idx_down[order(tab$FDR[idx_down])[1:range[3]]])
  tab$label <- NA
  tab$label[idx] <- tab$SYMBOL[idx]
  tab$alpha <- 1

  volcano(tab, paste0(opt$outdir,'/',exp,"_volcano_genes.png"))

  # Plot with region name
  if(exp == 'ATAC') {  
    tab$label <- NA
    tab$label[idx] <- tab$NAME[idx]
  
    volcano(tab,  paste0(opt$outdir,"/ATAC_volcano_regions.png"))
  }
}

sessionInfo()