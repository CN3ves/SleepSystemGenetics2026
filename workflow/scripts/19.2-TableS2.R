# Script to produce table S2

# Redirect all R logs to Snakemake log
log <- file('logs/19-Tables/S2_19.2.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("openxlsx")
  library("edgeR")
  library("optparse")
})

cat("Checking arguments\n")

option_list = list(
  make_option(c("-a", "--dar"), type="character", default=NULL, 
              help="Differential accessibility results' table", metavar="character"),
  make_option(c("-b", "--deg"), type="character", default=NULL, 
              help="Differential expression results", metavar="character"),
  make_option(c("-c", "--cors"), type="character", default=NULL, 
              help="Annotated gene to transcript correlations", metavar="character"),
  make_option(c("-d", "--S1"), type="character", default=NULL, 
              help="Table S1", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dar)){
  print_help(opt_parser)
  stop("Differential accessibility results (-a) is missing", call.=FALSE)
}
if (is.null(opt$deg)){
  print_help(opt_parser)
  stop("Differential expression results (-b) is missing", call.=FALSE)
}
if (is.null(opt$cors)){
  print_help(opt_parser)
  stop("Annotated gene to transcript correlations (-c) is missing", call.=FALSE)
}
if (is.null(opt$S1)){
  print_help(opt_parser)
  stop("Table S1 (-d) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

# Create a new workbook and add a sheet
wb <- createWorkbook()

cat("Load ATAC differential analysis\n")
dar <-  readRDS(opt$dar)
dar <- dar$table
dar$Peak_ID <- rownames(dar)

regions <- read.xlsx(opt$S1, "Final_Regions")
dar <- merge(dar,regions)

cat("Loading correlations\n")
cors <-  read.csv(opt$cors)
cors <- cors[,-1]
names(cors)[1:2] <- c("Nearest.Gene", "Peak_ID")

mean(apply(cors[,c("Nearest.Gene", "Peak_ID")],1,paste, collapse=' ') %in% apply(dar[,c("Nearest.Gene", "Peak_ID")],1,paste, collapse=' '))
#1

cat("Merge results\n")
res <- merge(dar,cors, by=c("Nearest.Gene", "Peak_ID"), all.x = TRUE)
res <- res[,c("Peak_ID", "Chromosome", "Start", "End", "logFC", "PValue", "FDR", "Annotation", "cCRE", "Nearest.Gene","cor_ctrl", "p_ctrl", "padj_ctrl", "cor_sd", "p_sd", "padj_sd", "cor_all", "p_all", "padj_all")]
  
res <- res[order(res$FDR, decreasing= FALSE),]

res$Differential_Correlated_bsl <-  res$FDR < 0.05 & res$padj_ctrl < 0.05
res$Differential_Correlated_sd <-  res$FDR < 0.05 & res$padj_sd < 0.05
res$Differential_Correlated_all <-  res$FDR < 0.05 & res$padj_all < 0.05

names(res) <- c("Region_ID", "Chromosome", "Start", "End", "ATAC_logFC", "ATAC_PValue", "ATAC_FDR","Annotation", "cCRE", "Nearest Gene", 
                "Pearson_correlation_baseline", "Correlation_PValue_baseline","Correlation_FDR_baseline",
                "Pearson_correlation_SD", "Correlation_PValue_SD", "Correlation_FDR_SD",
                "Pearson_correlation_all", "Correlation_PValue_all", "Correlation_FDR_all",
                "Differential and Correlated in baseline", "Differential and Correlated after SD", "Differential and Correlated in all samples")

res[is.na(res)] <- "-"

addWorksheet(wb, "Differential_Accessibility")
writeData(wb, "Differential_Accessibility", res, rowNames=FALSE)

cat("Load RNA differential analysis\n")
deg <-  readRDS(opt$deg)
deg <- deg$table
deg$Gene_ID <- rownames(deg)

deg <- deg[,c("Gene_ID", "logFC", "PValue", "FDR")]
names(deg) <- c("Gene_ID", "RNA_logFC", "RNA_PValue", "RNA_FDR")

addWorksheet(wb, "Differential_Expression")
writeData(wb, "Differential_Expression", deg, rowNames=FALSE)

cat("Save table\n")
saveWorkbook(wb, paste0(opt$outdir,"/TableS2-Differential_Correlation.xlsx"), overwrite = TRUE)

sessionInfo()