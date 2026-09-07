# Script to produce figure 2a

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure2a.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ggplot2")
  library("ggrepel")
  library("ComplexHeatmap")
  library("RColorBrewer")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--footprints"), type="character", default=NULL, 
              help="Footprint statistics", metavar="character"),
  make_option(c("-b", "--heatmap"), type="character", default=NULL, 
              help="Footprint with all samples", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$footprints)){
  print_help(opt_parser)
  stop("Footprint statistics (-a) is missing", call.=FALSE)
}
if (is.null(opt$heatmap)){
  print_help(opt_parser)
  stop("Footprint with all samples (-b) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load footprints\n")
meta <- read.csv(opt$footprints)
df <- read.delim(opt$heatmap)

samples <- unique(gsub("Protection_Score_","",names(df)[grep("Protection_Score_",names(df))]))

for(s in samples) {
  df[,s] <- rowSums(df[,grepl(s, names(df))])
}

df <- df[,-c(grep("Protection_Score_", names(df)), grep("TC_", names(df)))]

rownames(df) <- df$Motif
df$Motif <- NULL
df_plot <- df
df_plot$Num <- NULL

df_plot_scale <-t(apply(df_plot, 1, scale))
colnames(df_plot_scale) <- colnames(df_plot)
rownames(df_plot_scale) <- rownames(df_plot)

rowlabels <- rownames(df_plot_scale)  

sigs <-  meta[meta$padj < 0.05,'X']
idx <- unique(c(
  unlist(sapply(sigs, function(x) grep(paste0('^',x,'$'), rowlabels))),
  seq(1,length(rowlabels),5)))

rowlabels <- gsub("MA[0-9]*\\.[0-9]\\.","",rowlabels)
rowlabels[-idx] <- "" 

ra = rowAnnotation(
  labels = anno_text(rowlabels, which = "row", just="left", gp = gpar(fontsize = 5, col=ifelse(rowlabels %in% sigs, "purple", "black"))),
  show_legend = FALSE
 )

p <- Heatmap(as.matrix(df_plot_scale),
             name = "TF", 
             col = rev(brewer.pal(n = 11, name = "RdBu")),
             left_annotation = ra,
             
             cluster_rows = TRUE,  show_row_names = FALSE,
             show_row_dend = TRUE, row_dend_side = "right", 
             row_names_gp = gpar(fontsize = 4),
             row_split=8, #row_title = NULL, 
             row_gap = unit(2, "mm"),
             column_names_gp = gpar(fontsize = 8, col=ifelse(grepl("_SD", colnames(df_plot_scale)), "purple","darkgreen")),
             show_column_dend = TRUE,
             column_split=8, #column_title = NULL, 
             column_gap = unit(2, "mm"),
             
             #clustering_method_rows = "ward.D2",
             )


cat("Save plot\n")
svg(paste0(opt$outdir,"/Fig2a.svg"),height=10, width=10)
print(p)
dev.off()

sessionInfo()