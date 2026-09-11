# Script to produce figure S4a

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure4Sa.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("GenomicRanges")
  library("BSgenome.Mmusculus.UCSC.mm10")
  library("ggplot2")
  library("dplyr")
  library("openxlsx")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--S6"), type="character", default=NULL, 
              help="Table S6", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$S6)){
  print_help(opt_parser)
  stop("Table S6 (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load table S6\n")
QTL <- read.xlsx(opt$S6, "ATAC_QTL")

QTL$seqnames <- gsub(":.*","",QTL$region) 
QTL$start <- as.numeric(gsub(".*:(.*)-.*","\\1",QTL$region))
QTL$end <- as.numeric(gsub(".*-","",QTL$region) )

snps <- QTL[, c("chr","bp", "bp", "lod")]
names(snps) <- c("seqnames","start", "end", "lod")

snps$seqnames <- paste0("chr",snps$seqnames)
snps <- GRanges(snps)
snps$lod <- -snps$lod

QTL <- GRanges(QTL[,c("seqnames","start", "end", "lod")])
  
peaks <- c(QTL,snps)

# Adapted from ChIPseeker
coverage <- coverage(peaks, weight=mcols(peaks)[["lod"]])

pcov <- lapply(coverage, IRanges::slice, lower=0.001)
ncov <- lapply(coverage, IRanges::slice, upper=-0.001)
cov <- lapply(1:length(coverage), function(i) {
  c(pcov[[i]],ncov[[i]])
})
names(cov) <- names(pcov)

# Make data frame
ldf <- lapply(1:length(cov), function(i) {
    
    x <- cov[[i]]
    l <- length(BSgenome.Mmusculus.UCSC.mm10[[names(cov)[i]]])
    
    if (length(x@ranges) == 0) return(NA)
    
    rval <- runValue(x)
    rval <- sapply(rval@listData, mean)
    
    data.frame(chr   = names(cov[i]),
               start = start(x),
               end   = end(x),
               cnt   = rval,
               len = l
              )
  })

ldf <- ldf[!is.na(ldf)]
df <- do.call("rbind", ldf)

#Sort chrs
chr.name <- as.character(unique(df$chr))
noChr <- suppressWarnings(as.numeric(sub("chr", "", chr.name)))
ch.idx <- which(is.na(noChr))
n.idx <- which(!is.na(noChr))
chr.n <- noChr[n.idx]
chr.sorted <- chr.name[n.idx][order(chr.n)]
chr.sorted <- c(chr.sorted, sort(chr.name[ch.idx]))
df$chr <- factor(df$chr, levels=chr.sorted)

# Group data for plotting
tm <- group_by(df, chr, start, end, len) %>% summarise(value=sum(cnt))

#uniform snp values
tm$value[tm$value  <0] <- -mean(tm$value[tm$value  >0])

# Add colors
colors <- rep("darkorange", nrow(tm))
colors[tm$value > 0] <- "darkcyan"
tm$color <- colors
  
# get chr boundaries
limits <- unique(tm[,c("chr","len")])
limits$start <- limits$len # right limit up tick
limits$end <- limits$len
limits$color <- "black"
limits <- rbind(limits, limits)
limits$value <- c(rep(round(mean(abs(tm$value))),nrow(unique(tm[,c("chr","len")]))),
                  rep(-round(mean(abs(tm$value))),nrow(unique(tm[,c("chr","len")]))))
limits <- limits[,names(tm)] #sort limits to fit plot data
limits0 <- limits
limits0$start <- limits0$end <- 0 
tm <- do.call(rbind, list(tm,limits,limits0))

tm$color <- factor(tm$color, levels=c("darkcyan",  "darkorange","black"))

p <- ggplot(tm, aes(start, value, color=color, fill=color))
p <- p + geom_rect(aes(xmin = start, ymin = 0, xmax = end,  ymax = value))
p <- p + facet_grid(chr ~ ., scales = "free")
p <- p + theme_classic()
p <- p + xlab("Chromosome Size (bp)") + ylab("logFC") + ggtitle("aQTLs")
p <- p + scale_y_continuous(breaks=seq(-0.5,0.5,0.5))
p <- p + theme(strip.text.y = element_text(angle = 360),plot.title = element_text(hjust = 0.5))
p <- p  + scale_color_manual(values=c("darkcyan", "darkorange","black"), labels=c("Regions","QTL peak","Chr limit"))
p <- p  + scale_fill_manual(values=c("darkcyan", "darkorange","black"), labels=c("Regions","QTL peak","Chr limit"))
p <- p + guides(color=guide_legend(title="SD vs CTRL"), fill=guide_legend(title="aQTL"))
p <- p  + geom_hline(yintercept=0, linetype="dashed", color = "grey")

cat("Save plot\n")
svg(paste0(opt$outdir,"/FigS3c.svg"),width=16,height=7)
print(p)
dev.off()

write.csv(tm, paste0(opt$outdir, "/data/S4a.csv"))

sessionInfo()
