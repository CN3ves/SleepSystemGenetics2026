# Script designed to filter counts from RNA matrix

# Redirect all R logs to Snakemake log
log <- file('logs/7-BXD_normalization/object_7.1.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("edgeR")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="Directory with the count files", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="Processed metadata file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("Directory with the count files (-c) is missing", call.=FALSE)
}
if (is.null(opt$meta)){
  print_help(opt_parser)
  stop("Processed metadata file (-m) is missing", call.=FALSE)
}

cat("Reading inputs\n")
count_mat <- read.delim(opt$counts, row.names=1)
metadata <- read.delim(opt$meta)
metadata$genotype <- gsub("/RwwJ ", "", metadata$genotype)
metadata$genotype <- gsub("/TyJ ", "", metadata$genotype)
metadata$genotype <- gsub(".*F1", "F1", metadata$genotype)
metadata$genotype <- gsub(" \\(.*", "", metadata$genotype)
metadata$genotype <- gsub("C57BL/6J", "C57Bl6",metadata$genotype)
metadata$genotype <- gsub("DBA/2J", "DBA", metadata$genotype)
  
metadata$condition <- gsub("Control", "CTRL", metadata$condition)
metadata$condition <- gsub("Sleep Deprived", "SD", metadata$condition)
  
tissue <- "Cortex" # Liver not used in this work

# Process count matrix
counts <- count_mat[,metadata$tissue == tissue]
meta <- metadata[metadata$tissue == tissue,]
print(table(metadata[metadata$Sample_geo_accession %in% colnames(counts),"tissue"]))
counts <- DGEList(counts)

# Add metadata
counts$samples$Treatment <- meta[match(rownames(counts$samples),meta$Sample_geo_accession),"condition"]
counts$samples$Strain <- metadata[match(rownames(counts$samples),metadata$Sample_geo_accession),"genotype"]
counts$samples$group <- paste(counts$samples$Strain,counts$samples$Treatment, sep = "_")

saveRDS(counts, paste0(opt$outdir,"/rna_counts.RData"))

sessionInfo()