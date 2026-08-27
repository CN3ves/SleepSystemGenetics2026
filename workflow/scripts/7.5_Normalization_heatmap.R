# Script designed to cluster samples

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ComplexHeatmap")
  library("circlize")
  library("pheatmap")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="RDS counts object", metavar="character"),
  make_option(c("-a", "--annot"), type="character", default=NULL,
              help="Comma-separated metadata columns to annotate heatmap", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL,
              help="Processed metadata file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("RDS counts object (-c) is missing", call.=FALSE)
}
if (is.null(opt$meta)){
  print_help(opt_parser)
  stop("Processed metadata file (-m) is missing", call.=FALSE)
}
if (is.null(opt$annot)){
  print_help(opt_parser)
  stop("Metadata columns to annotate heatmap (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

exp <- gsub('_filter.*','',basename(opt$counts))

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/7-BXD_normalization/eda_7.5_', exp,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Reading inputs\n")

counts <- readRDS(opt$counts)
metadata <- read.csv(opt$meta, row.names=1)
rownames(metadata) <- sprintf("%03d",metadata$Sample)
if (exp == 'rna') metadata <- counts$samples
meta_annot <- strsplit(opt$annot, ",")[[1]]

# Extract the 2000 most variable features
c_mat <- sapply(1:ncol(counts$counts), function(x) counts$counts[,x] /(counts$samples$lib.size[x] * counts$samples$norm.factors[x]) )
colnames(c_mat) <- colnames(counts$counts)
s <- rowMeans((c_mat-rowMeans(c_mat))^2)
o <- order(s,decreasing=TRUE)
c_mat <- c_mat[o[1L:min(2000, nrow(c_mat))],,drop=FALSE]

write.csv(as.data.frame(rownames(c_mat)), paste0(opt$outdir,"/", exp, "2000most_variable.csv"))

# Make heatmap annotation
my_sample_col <- metadata[,meta_annot, drop = FALSE]
my_sample_col <- my_sample_col[gsub('_.*','',colnames(c_mat)),]
rownames(my_sample_col) <- colnames(c_mat)

annot_cols <- list()
f0_color <- list("DBA"=rgb(226/255,187/255,144/255), "C57Bl6"=rgb(148/255,151/255,152/255))

# Column annotations
for (col in meta_annot) {
   my_sample_col[,col] <- as.factor(my_sample_col[,col])
   uniq <- levels(my_sample_col[,col])
   color <-  rainbow(length(uniq))
   names(color) <-uniq
   if (col %in% c("Strain", "genotype")) {
     dba <- unique(grep("^DBA",metadata[,col], value = TRUE))
     dba <- names(sort(sapply(dba, nchar))[1])

     c57 <- unique(grep("^C57",metadata[,col], value = TRUE))
     c57 <- names(sort(sapply(c57, nchar))[1])


     color[dba] <- f0_color[['DBA']]
     color[c57]<-  f0_color[['C57Bl6']]
   }
   annot_cols[[col]] <- sapply(my_sample_col[,col], function(x) color[x])

}

# Row annotation
annotation_row = data.frame(
                  top500 = c(rep("100", 100), rep("500", 400), rep("2000", 1500))
                )
annotation_row <- annotation_row[1:nrow(c_mat),, FALSE]
rownames(annotation_row) <- rownames(c_mat)
annot_cols$top500 <- c("100"="red", "500"="grey", "2000"="white")

if (exp == "atac") {
  annot_cols$chr <- rainbow(21)
  annotation_row$chr = gsub("_.*","",rownames(c_mat))
  names(annot_cols$chr) <- unique(annotation_row$chr)
  annot_cols$chr <- annot_cols$chr[annotation_row$chr]
}


# Plot heatmap
png(paste0(opt$outdir,"/QC_", exp,"_pheat.png"))
hmap <- pheatmap(c_mat, scale = 'row',
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE,show_colnames = FALSE,
         fontsize_col = 2,
         border_color = 'NA',
         annotation_col = my_sample_col,
         annotation_row = annotation_row,
         annotation_colors = annot_cols,
         annotation_legend = FALSE,
         main = 'Top 2000 features counts',
         cellwidth=1.5)
dev.off()

# part 2: check for possible counfounders and metadata associations

# Instantiate association matrix
corr_p <- matrix(NA, nrow = ncol(metadata), ncol = ncol(metadata), dimnames = list(names(metadata), names(metadata)))

# Calculate Chi square p-vals
if("Age" %in% names(metadata)) {
  metadata$Age <- cut(metadata$Age,breaks = seq(11,15, 0.5))
}

for (var in names(metadata)) {
  for (covar in names(metadata)) {
    tab <- table(metadata[,var],metadata[,covar])
    tab <- tab[rowSums(tab)!=0, colSums(tab)!=0]
    corr_p[var,covar] <- chisq.test(tab, simulate.p.value=TRUE)$p.value
    corr_p[covar,var] <- corr_p[var,covar] 
  }
}
col_fun = colorRamp2(c(0, 0.05,0.5,1), c("cyan", "red", "white", "black"))

png(paste0(opt$outdir,"/",exp ,"_QC_covars.png"))
Heatmap(corr_p, border_gp = gpar(col = "black", lty = 1),col = col_fun(seq(0,1,0.05)),
        heatmap_legend_param = list(title = "ChiSq p-value"))
dev.off()

sessionInfo()