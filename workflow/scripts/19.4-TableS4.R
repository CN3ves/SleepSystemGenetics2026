# Script to produce table S4

# Redirect all R logs to Snakemake log
log <- file('logs/19-Tables/S4_19.4.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("openxlsx")
  library("org.Mm.eg.db")
  library("ChIPseeker")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
  library("tidyverse")
  library("optparse")
})

cat("Checking arguments\n")

option_list = list(
  make_option(c("-a", "--footprints"), type="character", default=NULL, 
              help="Footprint results directory", metavar="character"),
  make_option(c("-b", "--stats"), type="character", default=NULL, 
              help="Footprint meta-analysis", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$footprints)){
  print_help(opt_parser)
  stop("Footprint results directory (-a) is missing", call.=FALSE)
}
if (is.null(opt$stats)){
  print_help(opt_parser)
  stop("Footprint meta-analysis (-b) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

# Create a new workbook and add a sheet
wb <- createWorkbook()

cat("Load footprint analyses\n")
footprints <-  list.files(opt$footprints, pattern="bed", full.name=TRUE)

footprints <- lapply(footprints, function(file) {
  tab <- read.delim(file, header=FALSE)
  names(tab) <- c("seqnames", "start","end","motif","score","strand","something")
  tab$sample <- gsub(".*match//(.*)_mpbs.bed", "\\1", file)
  tab$motif_TF <- gsub(".*\\.[0-9]\\.","",tab$motif)
  tab$motif_ID <- gsub("(.*\\.[0-9])\\..*","\\1",tab$motif)

  tab <- tab[, -which(names(tab) %in% c("motif", "something"))]
  
  return(tab)
  })
  
footprints <- do.call(rbind,footprints)

footprints$sample[grep("BXD",footprints$sample, invert = TRUE)] <- paste0('BXD',footprints$sample[grep("BXD",footprints$sample, invert = TRUE)])

footprints$sample <- gsub("BXD([C|D])","\\1",footprints$sample)

footprints <- makeGRangesFromDataFrame(footprints, keep.extra.columns=TRUE)

cat("Annotate footprints\n")
annot <- annotatePeak(footprints, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion=c(-3000, 3000),   annoDb="org.Mm.eg.db")
annot <- as.GRanges(annot)

annot <- annot[,c("sample", "motif_TF", "motif_ID", "score",   "annotation", "distanceToTSS", "SYMBOL", "GENENAME", "transcriptId",  "ENSEMBL")]

annot <- as.data.frame(annot) %>% pivot_wider(names_from = sample,values_from=score) %>% GRanges

names(mcols(annot)) <- gsub("results/12-BXD_footprints/motifs/","",names(mcols(annot)))
names(mcols(annot)) <- gsub("_mpbs.bed","",names(mcols(annot)))

annot <- annot[, -grep('sub',names(mcols(annot)))]

cat("Filter table\n")
annot <- as.data.frame(annot) %>% select(-c('width',"strand",'ENSEMBL', "GENENAME"))
annot[,grep('_[SC]',names(annot))] <- lapply(grep('_[SC]',names(annot),value=TRUE), function(s) ifelse(is.na(annot[,s]),NA,gsub('^X','',s)))
annot$Strains <- apply(annot[,grep('_[SC]',names(annot))], 1,function(x) paste(sort(unique(x[!is.na(x)])), collapse = '; '))

annot$annotation <- gsub(" \\(.*","",annot$annotation)
annot <- annot[,-grep('_[SC]',names(annot))] 

addWorksheet(wb, "Footprints")
writeData(wb, "Footprints",annot,  rowNames=FALSE)

cat("Load footprint meta-analisys\n")

meta <- read.csv(opt$stats)
meta <- meta[order(meta$padj),]
names(meta)[1] <- 'Motif'

addWorksheet(wb, "Meta-analysis")
writeData(wb, "Meta-analysis", meta, rowNames=FALSE)

cat("Save table\n")
saveWorkbook(wb, paste0(opt$outdir,"/TableS4-Footprints.xlsx"), overwrite = TRUE)
   
sessionInfo()