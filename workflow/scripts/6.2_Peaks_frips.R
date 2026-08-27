# Script designed to estimate the Fraction of Reads in Peaks for a given list of features

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("Rsubread")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-f", "--features"), type="character", default=NULL, 
              help="Directory with the bed/gtf feature files", metavar="character"),
  make_option(c("-b", "--bam"), type="character", default=NULL, 
              help="Directory with the bam alignment files", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$features)){
  print_help(opt_parser)
  stop("Directory with the bed/gtf feature files (-f)", call.=FALSE)
}
if (is.null(opt$bam)){
  print_help(opt_parser)
  stop("Directory with the bam aligment files (-b)", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
bams <- list.files(opt$bam, pattern="bam$", full.names = TRUE)
features <- opt$features

# Redirect all R logs to Snakemake log
bname <- gsub('.*features/(.*)/genomic_features([0-9]*).gtf','\\1\\2',features)
log <- file(paste0('logs/6-BXD_features/frip_6.2_', bname,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Use featurecounts for Frip (standard method)
fc <-featureCounts(bams,annot.ext =features, nthreads =30,
    isGTFAnnotationFile=TRUE, GTF.featureType='Peak', GTF.attrType='Peak_ID',
    useMetaFeatures=FALSE, ,allowMultiOverlap =TRUE)

frips <- fc$stat[1,-1]/colSums(fc$stat[,-1])

names(frips) <- gsub('.co.bam', '', names(frips))
rownames(frips) <- bname

write.csv(frips, paste0(opt$outdir, '/frip_', bname, '.csv'))

sessionInfo()
