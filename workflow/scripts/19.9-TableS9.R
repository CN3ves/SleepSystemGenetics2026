# Script to produce table S9

# Redirect all R logs to Snakemake log
log <- file('logs/19-Tables/S9_19.9.log', open = "wt")
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
  make_option(c("-a", "--deg"), type="character", default=NULL, 
              help="Differential expression results", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$deg)){
  print_help(opt_parser)
  stop("Differential expression results (-b) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

# Create a new workbook and add a sheet
wb <- createWorkbook()

cat("Load Differential analyses\n")
stats <-  readRDS(opt$deg)

cat("SD vs baseline per group\n")
for (sheet in c("SD effect CT", "SD effect FV","SD effect FT")) {
  i <- which(tolower(gsub(".* ","",sheet)) == tolower(names(stats)))
  print(names(stats)[i])
  mat <- stats[[i]]$table
  mat$ENSEMBLE <- rownames(mat)
  mat <- mat[,c("ENSEMBLE", "SYMBOL", "NAME", "logFC", "logCPM", "F", "PValue", "FDR")]
  
  sheet <- gsub("FT","ciKO",gsub("FV","Tamoxifen Control",gsub("CT","Genotype Control",sheet)))
  sheet <- gsub("SD effect ","SD effect -",sheet)
  addWorksheet(wb, sheet)
  writeData(wb, sheet, mat, rowNames=FALSE)
    
}

cat("Interacition model\n")
for (sheet in c("baselineCT", "baselineFV", "SDFT", "SDCT", "SDFV" )) {
  mat <- stats[[sheet]]$table
  mat$ENSEMBLE <- rownames(mat)
  mat <- mat[,c("ENSEMBLE", "SYMBOL", "NAME", "logFC", "logCPM", "F", "PValue", "FDR")]

  if(sheet == "baselineCT") sheet <- "-Baseline- FT vs CT"
  if(sheet == "baselineFV") sheet <- "-Baseline- FT vs FV"
  if(sheet == "SDFT") sheet <- "-Interaction- SD x FT"
  if(sheet == "SDCT") sheet <- "-Interaction- SD x CT"
  if(sheet == "SDFV") sheet <- "-Interaction- SD x FV"
  
  print(sheet)
  sheet <- gsub("FT","ciKO",gsub("FV","Tamoxifen",gsub("CT","Genotype",sheet)))
  addWorksheet(wb, sheet)
  writeData(wb, sheet, mat, rowNames=FALSE)
    
}

cat("Save table\n")
saveWorkbook(wb, paste0(opt$outdir,"/TableS9-Differential_Expression.xlsx"), overwrite = TRUE)

sessionInfo()