# Script to produce table S8

# Redirect all R logs to Snakemake log
log <- file('logs/19-Tables/S8_19.8.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("openxlsx")
  library("optparse")
})

cat("Checking arguments\n")
option_list = list(
  make_option(c("-a", "--meta"), type="character", default=NULL, 
              help="Rna-seq metadata' table", metavar="character"),
  make_option(c("-b", "--pre_align"), type="character", default=NULL, 
              help="Read QC results' table", metavar="character"),
  make_option(c("-c", "--post_align"), type="character", default=NULL, 
              help="Alignment QC results' table", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$meta)){
  print_help(opt_parser)
  stop("RNA-seq metadata (-a) is missing", call.=FALSE)
}
if (is.null(opt$pre_align)){
  print_help(opt_parser)
  stop("Read QC results' table (-b) is missing", call.=FALSE)
}
if (is.null(opt$post_align)){
  print_help(opt_parser)
  stop("Alignment QC results' table (-c) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}


# Create a new workbook and add a sheet
wb <- createWorkbook()

cat("Loading metadata\n")
metadata <- read.csv(opt$meta)
names(metadata)[1] <- "Sample"
rownames(metadata) <- metadata$SampleName

cat("Load sequencing QC tab\n")
pre_align <- read.csv(opt$pre_align)
names(pre_align)[1:2] <- c("Metric", "Summary")

group_name <- metadata[gsub('^X','',names(pre_align)[-1:-2]),c("Sample", "Group")]
group_name <- apply(group_name, 1, paste, collapse=" (")
names(pre_align)[-1:-2] <- paste0(group_name,")")

addWorksheet(wb, "Sequencing_QC")
writeData(wb, "Sequencing_QC", pre_align, rowNames=FALSE)

cat("Load alignmnet QC tab\n")
post_align <- read.csv(opt$post_align)
post_align <- t(post_align)
post_align <- cbind(rownames(post_align),post_align)
post_align <- as.data.frame(post_align)
names(post_align) <- post_align[2,]
names(post_align)[1:2] <- c("Metric", "Summary")
post_align <- post_align[-1:-2,]

group_name <- metadata[names(post_align)[-1:-2],c("Sample", "Group")]
group_name <- apply(group_name, 1, paste, collapse=" (")
names(post_align)[-1:-2] <- paste0(group_name,")")

addWorksheet(wb, "Alignment_QC")
writeData(wb, "Alignment_QC", post_align, rowNames=FALSE)

# Save the workbook
cat("Save table\n")
saveWorkbook(wb, paste0(opt$outdir,"/TableS8-RNAQC_summary.xlsx"), overwrite = TRUE)

sessionInfo()