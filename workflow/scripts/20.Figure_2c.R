# Script to produce figure 2c

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure2c.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ComplexHeatmap")
  library("openxlsx")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--footprints"), type="character", default=NULL, 
              help="Footprint statistics", metavar="character"),
  make_option(c("-b", "--S5"), type="character", default=NULL, 
              help="Table S5", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$footprints)){
  print_help(opt_parser)
  stop("Footprint statistics (-a) is missing", call.=FALSE)
}
if (is.null(opt$S5)){
  print_help(opt_parser)
  stop("Table S5 (-b) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load footprints\n")
footprints <-  read.csv(opt$footprints)
footprints <- footprints[footprints$padj < 0.05,]
tfs <- gsub("MA[0-9]*\\.[0-9]\\.","",footprints$X)
tfs <- unique(toupper(gsub("\\(.*","",unlist(strsplit(tfs,"::")))))

col <- lapply(tfs, function(x) unique(footprints[grep(x,toupper(footprints$X)),"mean_activity"] >0))
col <- sapply(col, function(x) ifelse(length(x) ==1,ifelse(x,"darkgreen","purple"),"black"))
names(col) <- tfs

cat("Load KEA\n")
kea_list <- list(read.xlsx(opt$S5, sheet="Mean_rank"),
                read.xlsx(opt$S5, sheet="Integrated_scaled_rank"))

df <- data.frame(TF=tfs)
df[,unique(unlist(sapply(kea_list, function(x) gsub("\\*","",x$Protein))))] <- 0
rownames(df) <- df$TF
df <- df[,-1]

for(lst in 1:length(kea_list)) {
  kea <- kea_list[[lst]]
  
  for (i in 1:nrow(kea)) {
    targets <- strsplit(kea$Overlapping.Proteins[i],",")[[1]]
    df[targets,gsub("\\*","",kea$Protein[i])] <- 1
  }
}

kinases <- grep("\\*",unlist(lapply(kea_list, function(x) x$Protein)), value = TRUE)
kinases <- c(kinases,unlist(lapply(kea_list, function(x) x$Protein[1:10])))
kinases <- unique(gsub("\\*","",kinases))

p <- Heatmap(as.matrix(df[,colnames(df) %in% kinases]),
             col = c("white","blue"),
             heatmap_legend_param = list(at=0:1, labels = c("Not target", "Known target"),color_bar = "discrete"),
             name = "Target",
             row_names_gp = gpar(fontsize = 10,col=col),
             column_names_gp = gpar(fontsize = 10),
             column_title = "Kinases",
             row_title = "Significant TFs",
             column_title_side = "bottom"#, column_labels = gsub("\\*","",colnames(df[,grep("\\*",colnames(df))]))
              )

cat("Save plot\n")
svg(paste0(opt$outdir,"/Fig2c.svg"), width=15, height=10)
print(p)
dev.off()

write.csv(as.matrix(df[,colnames(df) %in% kinases]), paste0(opt$outdir, "/data/2c.csv"))

sessionInfo()