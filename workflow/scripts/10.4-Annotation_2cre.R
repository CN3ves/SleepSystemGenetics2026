# Script designed to annotate regions to ENCODE

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/10-BXD_annotation/encode_10.4.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("GenomicRanges")
  library("annotatr")
  library("ChIPseeker")
  library("ComplexHeatmap")
  library("circlize")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-e", "--encode"), type="character", default=NULL, 
              help="cCRE regions file", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=NULL, 
              help="Rdata containing peak annotation", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

   
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$annotation)){
  print_help(opt_parser)
  stop("Rdata containing containing peak annotation (-a) is missing", call.=FALSE)
}
if (is.null(opt$encode)){
  print_help(opt_parser)
  stop("ENCODE cSRE annotation (-e) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
peaks_obj <- readRDS(opt$annotation)

cre <- read.delim(opt$encode, header= FALSE) 
names(cre) <- c("seqnames", "start", "end", "something", "else","element")
cre <- GRanges(cre$seqnames, IRanges(start=cre$start, end=cre$end), mcols=cre[,"element", drop = FALSE])

for (test in names(peaks_obj)) {
  cat(paste("Annotate regions",test,"\n"))
  # Intersect the regions we read in with the annotations
  original_ranges <- as.GRanges(peaks_obj[[test]])
  annotated = annotate_regions(
    regions = original_ranges,
    annotations = cre,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  mcols(annotated)$element <- annotated$annot$mcols.element
  mcols(annotated) <- mcols(annotated)[, -which(names(mcols(annotated))=="annot")]
  
  missing <- original_ranges[!original_ranges$PeakID %in% annotated$PeakID]

  if(length(missing) >0 ) {
    mcols(missing)$element <- "NA"
  
    annotated <- unique(rbind(as.data.frame(annotated,row.names=NULL), as.data.frame(missing,row.names=NULL)))
    annotated <- GRanges(annotated[,1:5],mcols=annotated[,-1:-5])
  }

  
  names(mcols(annotated)) <- gsub("mcols.", "", names(mcols(annotated)))

  peaks_obj[[test]] <- annotated
  
  # Pie Charot with Percentages
  slices <- round(table(annotated$element) / length(annotated)*100,2)
  lbls <-  paste0(names(slices),": ", slices, "%")

  png(paste0(opt$outdir,"/screen_pie_", test, ".png"))
  par(mar=c(5.1, 4.1, 4.1, 10.5))
  pie(slices,labels = lbls, col=rainbow(length(lbls), v=0.75),
      main=paste("Distribution of regions for", test), cex = 0.8) 
  dev.off()
  
  tab <- table(gsub(" \\(.*", "",annotated$annotation),gsub(",.*","",annotated$element))
  tab <- matrix(tab, ncol = ncol(tab), dimnames = dimnames(tab))

  
  col_fun = colorRamp2(c(0,0.1, 0.2, 1), c("blue","cyan", "green", "red"))
  png(paste0(opt$outdir,"/ucsc_encode_", test, ".png"))
  print(Heatmap(apply(tab, 2, function(x) x/sum(x)), col =  col_fun(seq(0, 1, 0.1)), name =  paste("% UCSC\n", test), column_title = "ENCODE annotation", row_title = "UCSC TxDb annotation",show_column_dend = FALSE, show_row_dend = FALSE))
  dev.off()
  
  png(paste0(opt$outdir,"/encode_ucsc", test, ".png"))
  print(Heatmap(apply(tab, 1, function(x) x/sum(x)), col =  col_fun(seq(0, 1, 0.1)), name =  paste("% ENCODE\n", test), column_title = "ENCODE annotation", row_title = "UCSC TxDb annotation", show_column_dend = FALSE, show_row_dend = FALSE))
  dev.off() 
}

lapply(peaks_obj, function(x) table(x$element))

cat("Saving annotations\n")
saveRDS(peaks_obj, paste0(opt$outdir,"/cCRE_annotated.RData"))


sessionInfo()