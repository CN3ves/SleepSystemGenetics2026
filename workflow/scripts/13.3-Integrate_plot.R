# Script designed to plot overlaps between RNA and ATAc results

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/13-BXD_integrate/plots_13.3.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ggplot2")
  library("VennDiagram")
  library("GenomicRanges")
  library("RColorBrewer")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-r", "--rna_diff"), type="character", default=NULL, 
              help="Differential results for RNA", metavar="character"),
  make_option(c("-q", "--rna_qtl"), type="character", default=NULL, 
              help="QTL results for RNA", metavar="character"),
  make_option(c("-i", "--rna_int"), type="character", default=NULL, 
              help="Interaction QTL results for RNA", metavar="character"),
  make_option(c("-y", "--atac"), type="character", default=NULL, 
              help="Annotation for ATAC results", metavar="character"),
  make_option(c("-a", "--atac_qtl"), type="character", default=NULL, 
              help="QTL results for ATAC", metavar="character"),
  make_option(c("-c", "--atac_int"), type="character", default=NULL, 
              help="Interaction QTL results for ATAC", metavar="character"),
  make_option(c("-s", "--sleep"), type="character", default=NULL, 
              help="QTL results for sleep phenotypes", metavar="character"),
  make_option(c("-z", "--cors_annot"), type="character", default=NULL, 
              help="Correlations results between transcripts and regions based on annotation", metavar="character"),
  make_option(c("-x", "--cors_int"), type="character", default=NULL, 
              help="Correlations results between transcripts and regions based on interaction QTLs", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$rna_diff)){
  print_help(opt_parser)
  stop("Differential results for RNA (-r) is missing", call.=FALSE)
}
if (is.null(opt$rna_qtl)){
  print_help(opt_parser)
  stop("QTL results for RNA (-q) is missing", call.=FALSE)
}
if (is.null(opt$rna_int)){
  print_help(opt_parser)
  stop("Interaction QTL results for RNA (-i) is missing", call.=FALSE)
}
if (is.null(opt$atac)){
  print_help(opt_parser)
  stop("Annotation for ATAC results (-y) is missing", call.=FALSE)
}
if (is.null(opt$atac_qtl)){
  print_help(opt_parser)
  stop("QTL results for ATAC (-a) is missing", call.=FALSE)
}
if (is.null(opt$atac_int)){
  print_help(opt_parser)
  stop("Interaction QTL results for ATAC (-c) is missing", call.=FALSE)
}
if (is.null(opt$sleep)){
  print_help(opt_parser)
  stop("QTL results for sleep phenotypes (-s) is missing", call.=FALSE)
}
if (is.null(opt$cors_annot)){
  print_help(opt_parser)
  stop("Correlations results between transcripts and regions based on annotation (-z) is missing", call.=FALSE)
}
if (is.null(opt$cors_int)){
  print_help(opt_parser)
  stop("Correlations results between transcripts and regions based on interaction QTLs (-x) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading files\n")
plot <- function(sets,cols,names, file) {

  sets <- lapply(sets, function(x) x[!is.na(x)])
  venn.diagram(
        x = sets,
        category.names = names,
        filename = file,
        output=FALSE,
        disable.logging = TRUE,
        
        # Output features
        imagetype="png" ,
        height = 480 , 
        width = 480 , 
        resolution = 300,
        compression = "lzw",
        
        # Circles
        lwd = 2,
        lty = 'blank',
        fill = cols,
        
        # Numbers
        cex = .6,
        fontface = "bold",
        fontfamily = "sans",

        
        # Set names
        cat.cex = 0.2,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        #cat.pos = c(-10, 180),
        #cat.dist = c(0.02, 0.02),
        cat.fontfamily = "sans",
        cat.col =cols
  )
}

myCol4 <- brewer.pal(4, "Pastel2")
myCol3 <- myCol4[1:3]
myCol2 <- myCol3[-2]

cat("Plot Venn for differential analyses\n")
rna_diff <- read.csv(opt$rna_diff)
atac <- readRDS(opt$atac)

plot(list(unique(rna_diff[rna_diff$FDR<0.05,'X']), unique(atac[['DAR']][atac[['DAR']]$FDR <0.05,]$SYMBOL)), myCol2,
  c("DEG Transcripts" , "DAR Genomic Regions"), paste0(opt$outdir,"/venn_diff.png"))

cat("Plot Venn for QTL analyses\n")
rna_qtl <- read.csv(opt$rna_qtl)
atac_qtl <- read.csv(opt$atac_qtl)
sleep <- read.csv(opt$sleep)

plot(list(unique(rna_qtl$snp), unique(atac_qtl$snp), unique(sleep$snp)), myCol3,
  c("SNPs for Transcripts" , "SNPs for Genomic Regions", "SNPs for sleep phenotypes"), paste0(opt$outdir,"/venn_qtl_snps.png"))

plot(list(unique(rna_qtl[,'lodcolumn']), unique(atac[['QTL']]$SYMBOL)), myCol2,
  c("Transcripts with QTL" , "Genomic Regions with QTL"), paste0(opt$outdir,"/venn_qtl.png"))

cat("Plot Venn for QTL interaction analyses\n")
rna_int <- read.csv(opt$rna_int)
atac_int <- read.csv(opt$atac_int)

plot(list(unique(rna_int$snp), unique(atac_int$snp), unique(sleep$snp)), myCol3,
  c("Interaction SNPs for Transcripts", "Interaction SNPs for Genomic Regions", "SNPs for sleep phenotypes"),paste0(opt$outdir,"/venn_int_snps.png"))

plot(list(unique(rna_int[,'lodcolumn']), unique(atac[['QTLxSD']]$SYMBOL)), myCol2,
  c("Transcripts with interaction QTL", "Genomic Regions with interaction QTL"),paste0(opt$outdir,"/venn_int.png"))

plot(list(unique(rna_qtl$snp), unique(rna_int$snp), unique(atac_qtl$snp), unique(atac_int$snp)), myCol4,
  c("SNPs for Transcripts", "Interaction SNPs for Transcripts" ,  "SNPs for Genomic Regions", "Interaction SNPs for Genomic Regions"),paste0(opt$outdir,"/venn_intqtl_snps.png"))

plot(list(unique(rna_qtl[,'lodcolumn']), unique(rna_int[,'lodcolumn']), unique(atac[['QTL']]$SYMBOL), unique(atac[['QTLxSD']]$SYMBOL)), myCol4,
  c("Transcripts with QTL", "Transcripts with interaction QTL", "Genomic Regions with QTL","Genomic Regions with interaction QTL"), paste0(opt$outdir,"/venn_intqtl.png"))


cat("Plot p-value distributions for correlations\n")

correlations <- read.csv(opt$cors_annot, row.names=1)
correlations$Group <- 'Not sig.'
correlations$Group[correlations$p_all < 0.05] <- 'p < 0.05'
correlations$Group[correlations$padj_all < 0.05] <- 'padj < 0.05'

png(paste0(opt$outdir,"/cors_annotation.png"))
ggplot(data = correlations,aes(x = cor_all, fill = Group)) +
  geom_histogram(colour = "black")+
  theme_classic() + 
  ylab("Number of Gene/region pairs") + xlab("Pearson Correlation")
dev.off()

correlations <- read.csv(opt$cors_int, row.names=1)
correlations$Group <- 'Not sig.'
correlations$Group[correlations$p < 0.05] <- 'p < 0.05'
correlations$Group[correlations$padj < 0.05] <- 'padj < 0.05'

png(paste0(opt$outdir,"/cors_int.png"))
ggplot(data = correlations,aes(x = cor, fill = Group)) +
  geom_histogram(colour = "black")+
  theme_classic() + 
  ylab("Number of Gene/region pairs") + xlab("Pearson Correlation")
dev.off()


sessionInfo()