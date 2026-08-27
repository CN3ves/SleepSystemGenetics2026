# Script designed to confirm chr10 abnormalisties in mm10 annotation and correct with mm9 

# Redirect all R logs to Snakemake log
log <- file('logs/4-BXD_genotype/validate.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("optparse")
})
source("workflow/scripts/3.0-QC_align_helper.R")

cat("Checking arguments\n")
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="Directory with the BXD genotype files", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="Processed metadata file", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("Directory with the BXD genotype files (-d) is missing", call.=FALSE)
}
if (is.null(opt$meta)){
  print_help(opt_parser)
  stop("Processed metadata file (-m) is missing", call.=FALSE)
}

cat("Reading inputs\n")
metadata <- read.csv(opt$meta, row.names=1)
metadata$Sample <- sprintf("%03d",metadata$Sample)
bxd_gen10 <- read.csv(paste0(opt$dir,"/bxd_geno.csv"), skip=3)
bxd_snps <- read.csv(paste0(opt$dir,"/bxd_pmap.csv"), skip=3)
bxd_gen9 <- read.delim(paste0(opt$dir,"/BXD_mm9.geno.txt"), skip = 14)

# Get list of used BXD strains
strains <- unique(metadata$Strain)
strains <- strains[! strains %in% c("C57Bl6", "DBA")]
#Add BXD label to the "numeric" strains
strains[grep("BXD", strains, invert = TRUE)] <- paste0("BXD", strains[grep("BXD", strains, invert = TRUE)] )

# Filter for dbSNP Reference SNP ans strains
bxd_gen10 <- bxd_gen10[grep("rs",bxd_gen10$marker),]
bxd_gen10 <- bxd_gen10[, c("marker", strains)]

# Correct SNP marker location (from Mb to bp)
bxd <- merge(bxd_gen10, bxd_snps)
bxd$pos <- bxd$pos * 1000000 
rownames(bxd) <- bxd$marker

# Correct issue with BXD100 (10.1371/journal.pcbi.1010552) on mm10!
print(paste("Confirm mm10 issue with BXD100. Unique allele on chr 8, 10-19, X and Y:", unique(bxd[!bxd$chr %in% c(1:7,9), "BXD100"])))

bxd_gen9 <- bxd_gen9[,c("Chr","Locus","BXD100")]
rownames(bxd_gen9) <- bxd_gen9$Locus

# Input allele information from mm9
bad_locus <- rownames(bxd)[!bxd$chr %in% c(1:7,9)]
print(paste("Proportion of problematic mm10 loci in mm9:", round(mean(bad_locus %in% bxd_gen9$Locus) * 100,2), "%"))
locus <-  bad_locus[bad_locus %in% bxd_gen9$Locus]

stopifnot(all(bxd[rownames(bxd) %in% locus,"marker"] == bxd_gen9[locus,"Locus"]))

print("Changes in allelles:")
print(table(bxd[rownames(bxd) %in% locus,"BXD100"], bxd_gen9[locus,"BXD100"]))
bxd[rownames(bxd) %in% locus,"BXD100"] <- bxd_gen9[locus,"BXD100"]

saveRDS(bxd, file = paste0(opt$dir,"/genotypes.Rdata"))

sessionInfo()