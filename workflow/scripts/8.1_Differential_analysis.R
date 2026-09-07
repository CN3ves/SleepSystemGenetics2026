# Script designed for differential analyses

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("edgeR")
  library("ggplot2")
  library("tidyr")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="Rdata file with the normalized feature counts", metavar="character"),
  make_option(c("-f", "--frip"), type="character", default=NULL, 
              help="Table with Fraction of reads in peak considering mitochondrial DNA", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("Rdata file with the normalized feature counts (-c) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

exp <- gsub('_filter.*','',basename(opt$counts))
if(!is.null(opt$frip)) exp <- paste0(exp,"_frip")

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/8-BXD_differential/DA_8.1_', exp,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Reading inputs\n")
counts <- readRDS(opt$counts)
if(!is.null(opt$frip)) { # test the impact of low FRIP samples
  frips <- read.csv(opt$frip, row.names=1)
  frips <-frips[grep('full', names(frips))]
  rm <- gsub('_full','',gsub('^X','',names(frips)[frips <= 0.2]))
  print(paste("Excluded sample", colnames(counts)[colnames(counts) %in% rm], "due to low FRiP score"))
  counts <- counts[,!colnames(counts) %in% rm]
}

depth <- colSums(counts$counts)/nrow(counts$counts)
print(paste("Minimum estimated depth of experiment:", min(depth)))
print(paste("Percentage of samples with ok depth of experiment:", mean(depth>3)))
print(paste("Percentage of samples with good depth of experiment:", mean(depth>10)))

#Distribution of depth
png(paste0(opt$outdir, "/", exp, "_depth_experiment.png"))
hist(depth, breaks=30, freq = TRUE, main="Estimated depth of experiments", xlab="Depth")
dev.off()

#  Estimating dispersion
treatment <- as.factor(counts$samples$Treatment)
treatment <- relevel(treatment,"CTRL")
lines <- as.factor(counts$samples$Strain)
lines <- relevel(lines,"C57Bl6")

design <- model.matrix(~0 + lines + treatment)
colnames(design) <- gsub("lines", "", gsub("treatment", "", colnames(design)) )

plot.design <- as.data.frame(design)
plot.design$Sample <- paste0(counts$samples$group, "(", rownames(counts$samples), ")")
plot.design <-  plot.design %>% 
    pivot_longer(1:(ncol(plot.design)-1), names_to = "perc", values_to ="val")

png(paste0(opt$outdir, "/", exp, "_design_main.png"))
ggplot(plot.design, aes(Sample, perc, fill = c("white", "black")[ val + 1 ])) +
  geom_tile() +
  scale_fill_identity() +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(size = 5, angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(size=5)) + 
  xlab("Coefficients") + ylab("Samples") + 
  ggtitle("Model Matrix")
dev.off()

disp <- estimateDisp(counts,design) # May take 3-4 hours
print(paste("Common feature dispersion", round(disp$common.dispersion, 5)))

png(paste0(opt$outdir, "/", exp,"_BCV_treament.png"))
plotBCV(disp)
dev.off()

# Testing for DA
fit <- glmQLFit(disp,design)

qlf <- glmQLFTest(fit,coef="SD")
qlf$table$FDR <- p.adjust(qlf$table$PValue, method='BH')
qlf$table <- qlf$table[order(qlf$table$FDR),]

stats <- summary(decideTests(qlf))
print(paste(rownames(stats), stats), sep = ": ")

changed <- sum(summary(decideTests(qlf))[-2,1])
all <- sum(summary(decideTests(qlf))[,1])
print(paste0("Percentage of significant features ", round(sum(changed)/all*100,2),"%"))

#Distribution of pvalues
png(paste0(opt$outdir, "/", exp, "pvals_treament.png"))
hist(qlf$table$PValue, freq = FALSE, main=" Distribution of p-values (SD vs CTRL)", xlab="p-value")
lines(density(qlf$table$FDR), col = "red")
legend("topright", legend="BH-corrected", lwd =1, col="red")
dev.off()

# MA plot
png(paste0(opt$outdir, "/", exp, "_MA_treatment.png"))
plotMD(qlf)
abline(h=c(-0.5, 0.5), col="green")
abline(h=c(-0.1, 0.1), col="green")
dev.off()

write.csv(qlf$table, paste0(opt$outdir, "/", exp, "_peaks_treament.csv"))
saveRDS(disp, paste0(opt$outdir, "/", exp, "_counts_disp.RData"))
saveRDS(qlf, paste0(opt$outdir, "/", exp, "_stats_treament.RData"))

sessionInfo()