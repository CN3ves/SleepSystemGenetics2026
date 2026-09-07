# Script designed to aggregate all read QC summaries

# Redirect all R logs to Snakemake log
log <- file('logs/1-BXD_fastq/QCsummary_1.1.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("rjson")
  library("dplyr")
  library("optparse")
})
source("workflow/scripts/1.0-QC_reads_helper.R")

cat("Checking arguments\n")
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="Directory with the Fastp json output", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="Processed metadata file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("Directory containing the Fastp .json output (-d) is missing", call.=FALSE)
}
if (is.null(opt$meta)){
  print_help(opt_parser)
  stop("Processed metadata file (-m) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
reports <- list.files(opt$dir, full.names = TRUE)
metadata <- read.csv(opt$meta, row.names=1)

cat("Collecting QC data summaries\n")
# Collect QC data summaries
qc_table <- list()
for (report in reports) {
  # Read json report
  qc <- fromJSON(file = report)
  
  # Summary QC
  res <- qc_summary(qc)
  # Duplication QC 
  res <- c(res, qc_duplication(qc))
  # Adapter QC   
  res <- c(res,qc_adaptor(qc))
  # Filtering QC   
  for (time in c("before", "after")) {
    res <- c(res, qc_filter(qc,time))
  }
  
  # Save to table
  report <- strsplit(report, "/")[[1]]
  report <-  gsub("\\..*","", report[length(report)])
  qc_table[[report]] <- matrix(as.character(res), ncol = 1, dimnames = list(names(res), report))
  
}

cat("Generating summary table\n")
# Create placeholder table with sorted rownames
rows <- sort(unique(unlist(sapply(qc_table, rownames))))
idx <-  grep("_adapter_counts", rows)
rows <- c(rows[-idx], rows[idx])

summary_tab <- matrix(rows, ncol = 1, dimnames = list(rows, "summary"))

# Add all samples to table
for (tab in qc_table) {
  tab <- tab[match(rownames(summary_tab), rownames(tab)),,drop = FALSE]
  summary_tab <- cbind(summary_tab, tab)
 }

# Summarize results
summary_tab <- as.data.frame(summary_tab, stringsAsFactors=FALSE)
summary_tab$summary <- sapply (1:nrow(summary_tab), summarize, summary_tab)

cat("Checking low read samples\n")
# Samples with reads below 50M
idx <- as.numeric(summary_tab["total_reads (before filtering)",-1]) < 50e6
idx <- summary_tab["total_reads (before filtering)",c(FALSE,idx)]
san_tab <- metadata[metadata$Sample %in% as.numeric(gsub("_.*","",names(idx))), c("Strain", "Mice_ID", "Treatment", "Sample")]
san_tab$Reads <- as.numeric(idx)
print("SAMPLES with low counts (<50M) before filtering")
print(san_tab)

cat("Saving table\n")
# Write read's QC table
write.csv(summary_tab, paste0(opt$outdir,"/QC_fastq.csv"))

sessionInfo()