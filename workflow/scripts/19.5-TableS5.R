# Script to produce table S5

# Redirect all R logs to Snakemake log
log <- file('logs/19-Tables/S5_19.5.log', open = "wt")
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
  make_option(c("-a", "--kea"), type="character", default=NULL, 
              help="Directory with files downloaded from KEA3", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$kea)){
  print_help(opt_parser)
  stop("Directory with files downloaded from KEA3 (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

# Create a new workbook and add a sheet
wb <- createWorkbook()

files <- list.files(opt$kea, full.names=TRUE)

cat("Load footprint analyses\n")
for (table in files){ 
  tab <- read.delim(table)
  sheet <- gsub(" ","_", gsub(".tsv","",basename(table)))
  addWorksheet(wb, sheet)
  writeData(wb, sheet, tab, rowNames=FALSE)

}

cat("Save table\n")
saveWorkbook(wb, paste0(opt$outdir,"/TableS5-KEA_enrichment.xlsx"), overwrite = TRUE)

sessionInfo()