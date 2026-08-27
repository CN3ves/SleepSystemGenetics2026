# Script designed to aggregate summaries from BAM files

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ATACseqQC")
  library("Rsamtools")
  library("preseqR")
  library("foreach")
  library("doParallel")
  library("ggplot2")
  library("GenomicAlignments" )
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
  library("optparse")
})
source("workflow/scripts/3.0-QC_align_helper.R")

cat("Checking arguments\n")
option_list = list(
  make_option(c("-b", "--bam"), type="character", default=NULL, 
              help="Bam file to QC", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bam)){
  print_help(opt_parser)
  stop("Bam file (-b) missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}


bamfile <- opt$bam
bname <- gsub('\\..*','',basename(bamfile))
print("Reading bam file")

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/3-BXD_bam/bamQC_', bname,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Reading inputs\n")
# Load bam file
possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                 "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                               "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                               "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                               "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                               "U2"))
param <- ScanBamParam(what = c("rname", "strand", "cigar", "pos", "mrnm", "mpos", "isize", "qwidth"), flag = scanBamFlag(isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE, isNotPassingQualityControls = FALSE), tag=unlist(possibleTag))
  
bam <- scanBam(bamfile, index = bamfile, param = param, asMate = FALSE)
  
tags <- names(bam[[1]]$tag)[lengths(bam[[1]]$tag)>0]
gal <- readBamFile(bamfile, tag=tags, bigFile=TRUE)
  
bam[[1]] <- bam[[1]][-which(names(bam[[1]]) == "tag")]
  
# Load transcript models
txs <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
seq_idx <- grep("_", seqlevels(txs), invert = TRUE)
txs <- txs[seqnames(txs) %in% seqlevels(txs)[seq_idx]]
seqlevels(txs) <- seqlevels(txs)[seq_idx]

# General QC 
print("Running General QC")
qc <-  bamQC(bamfile, outPath=NULL)

# Complexity
print("Getting complexity estimates")
histFile <- readsDupFreq(bam)
complex <- estimateLibComplexity(histFile)
rm("histFile")

# Read length
print("Getting Read lengths")
read.len <-table(bam[[1]]$qwidth)
x <- 1:300
frag.len <- read.len[match(x, names(read.len))]
frag.len[is.na(frag.len)] <- 0
y <- frag.len/sum(frag.len) 
y <- as.numeric(y)
read_lens <- data.frame("x"=x, "y"=y* 10^3)

# QC enrichment Scores
print("Getting QC Scores")
pt <- PTscore(gal, txs, seqlev = seqlevels(txs))
pt <- unique(data.frame("log"=pt$log2meanCoverage, "PT"= pt$PT_score, "chr"=seqnames(pt)))

nfr <- NFRscore(gal, txs, seqlev = seqlevels(txs))
nfr <- unique(data.frame("log"=nfr$log2meanCoverage, "NFR"= nfr$NFR_score, "chr"=seqnames(nfr)))

tsse <- TSSEscore(gal, txs, seqlev = grep("chr[0-9YZ]+",seqlevels(txs), value=T))

# Coverage
print("Getting mean values for correlation")
cov <- mean_cov(gal, txs) 

# Save results as RData objects
res <- list("QC" = qc, "Complexity"=complex, "Length"=read_lens, "PT"=pt, "NFR"=nfr, "TSSE"=tsse, "VM"= cov["vms"], "sex"=cov["sex"])
file <- paste0(opt$outdir,"/",bname,".RData")
print(paste("Storing data on", file))
save(res, file=file)

sessionInfo()