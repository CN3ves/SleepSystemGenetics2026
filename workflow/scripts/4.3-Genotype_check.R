# Script designed compare the called variants for the BXD samples to the BXD reference

# Redirect all R logs to Snakemake log
log <- file('logs/4-BXD_genotype/check_4.3.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ComplexHeatmap")
  library("circlize")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="Directory with the sample variants", metavar="character"),
  make_option(c("-g", "--geno"), type="character", default=NULL, 
              help="Rdata file with the genotype object", metavar="character"),
  make_option(c("-b", "--bxd"), type="character", default=NULL, 
              help="File with annotated BXD gentypes", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="Processed metadata file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("Directory with the sample variants (-d) is missing", call.=FALSE)
}
if (is.null(opt$geno)){
  print_help(opt_parser)
  stop("Rdata file with the genotype object (-g) is missing", call.=FALSE)
}
if (is.null(opt$bxd)){
  print_help(opt_parser)
  stop("File with annotated BXD gentypes (-b) is missing", call.=FALSE)
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
bxd <- readRDS(opt$geno)
genos <- list.files(opt$dir, pattern="[0-9].*csv", full.names=T)
metadata <- read.csv(opt$meta, row.names=1)
metadata$Sample <- sprintf("%03d",metadata$Sample)
bxd_snps <- read.csv(paste0(opt$bxd,"/bxd_pmap.csv"), skip=3)

genos <- sapply(genos, read.csv)
genos <- genos[-1,]
colnames(genos) <- gsub(".*/([0-9]*)_.*.csv","\\1", colnames(genos))
genos <- as.data.frame(genos)

strains <- metadata$Strain[metadata$Sample %in% names(genos)]

# Call parent genotypes
call_allele <- function(mat) {
  uniq_mat <- apply(mat, 1, function(x) {
    x <- gsub(",<\\*>", "",x)
    x <- x[!is.na(x) & x != "NA"]
    x <- x[x != "<*>"]
    x <- gsub(",.*", "",x)
    x <- unique(x)
    if(length(x) !=1 ) x <- "X"
    
    return(x)
  })
  names(uniq_mat) <- rownames(mat)
  return(uniq_mat)
}

b_allele <- call_allele(genos[,strains=="C57Bl6"])
d_allele <- call_allele(genos[,strains=="DBA"])

no_info_snps <- which((b_allele == d_allele) | (b_allele == "X") |( d_allele == "X") )
b_allele <- b_allele[-no_info_snps]
d_allele <- d_allele[-no_info_snps]

bxd_sub <- bxd[unique(names(b_allele), names(d_allele)),grep("BXD", names(bxd))]
samples <- genos[rownames(bxd_sub),!strains %in% c("C57Bl6", "DBA")]

# Match genopytes to parent strains
classify <- function(samples, genotypes, subset, metadata) {
 
  # Prepare object for storing call genotypes
  res <- metadata[metadata$Sample %in% names(samples), c("Sample", "Strain")]
  res$assign <- NA
  rownames(res) <- res$Sample
  res$Sample <- NULL
  res$Strain[grep("BXD", res$Strain, invert=TRUE)] <- paste0("BXD", res$Strain[grep("BXD", res$Strain, invert=TRUE)] )

  # Prepare object to store percentage of matching
  mat <- as.data.frame(matrix(NA, nrow=nrow(res), ncol=ncol(subset), dimnames = list(rownames(res),names(subset))))
  
  for (sample in names(samples)) {
    alleles <- call_allele(samples[,sample, drop=FALSE])
    alleles <- alleles[rownames(genotypes)]
    alleles <- alleles[!is.na(alleles)]
    for (a in names(alleles)) {
      if (length(alleles[a]) > 0) {
        b <- alleles[a]==b_allele[a] 
        d <- alleles[a]==d_allele[a] 
        if(length(b) ==0 ) b <- FALSE
        if(length(d) ==0 ) d <- FALSE
        if(b&d) {
          alleles[a] <- "H"
        } else if (b) {
          alleles[a] <- "B"
        } else if (d) {
          alleles[a] <- "D"
        } else {
          alleles[a] <- NA
        }
      }
    }
  
  alleles <- alleles[sapply(alleles,length) > 0]
  strains <- sapply(bxd_sub[names(alleles),], function(x) {
    hs <- x == "H" | as.character(alleles) == "H" | is.na(alleles)
    mean(x[!hs] == as.character(alleles)[!hs])
    })
  res[sample,"assign"] <- names(which.max(strains))
  mat[sample,] <- strains
  }
  
  return(list("assign"=res, "match"=mat))

}

geno_res <- classify(samples, bxd, bxd_sub, metadata)
print("Misclassified samples:")
print(geno_res[["assign"]][geno_res[["assign"]][,1] != geno_res[["assign"]][,2],])

# Try without bad mm10 chromossomes
good_chr <- bxd_snps[bxd_snps$chr %in% c(1:7,9), "marker"]
bxd_sub <- bxd_sub[good_chr, ]
samples <- genos[rownames(bxd_sub),!strains %in% c("C57Bl6", "DBA")]

geno_res_good <- classify(samples, bxd, bxd_sub, metadata)
print("Misclassified samples (no bad mm10 chromosomes):")
print(geno_res_good[["assign"]][geno_res_good[["assign"]][,1] != geno_res_good[["assign"]][,2],])

# Genotype call plots
M <- geno_res[["match"]]
rownames(M) <- paste0(metadata$Strain[match(rownames(M),metadata$Sample)], " (", rownames(M),")")
rownames(M)[grep("BXD", rownames(M), invert = TRUE)] <- paste0("BXD", rownames(M)[grep("BXD", rownames(M), invert = TRUE)] )
  
M <- M[order(rownames(M)),order(colnames(M))]
  
col_fun = colorRamp2(c(0,0.5,1), c("red",  "white", "blue"))

png(paste0(opt$outdir,"/Genotype_heat.png"))
Heatmap(as.matrix(M), border_gp = gpar(col = "black", lty = 1),
        column_title = "BXD Genotyping",
        heatmap_legend_param = list(title = "% SNPs matching"),
        cluster_rows = FALSE, cluster_columns = FALSE,
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 3),
        col = col_fun(seq(0,1,0.5)))
dev.off()

mixed <- apply(M, 1, function(x) sort(x, decreasing=T)[1]-sort(x, decreasing=T)[2]<0.1)
mixed <- unique(c(gsub(" .*", "",names(which(mixed))),"BXD48", "BXD49"))
mixed <- sapply(mixed, function(x) grep(x, rownames(M)))
mixed <- unique(unlist(mixed))
  
png(paste0(opt$outdir,"/Genotype_heat_focus.png"))
Heatmap(as.matrix(M[mixed ,]), border_gp = gpar(col = "black", lty = 1),
        column_title = "BXD Genotyping",
        heatmap_legend_param = list(title = "% SNPs matching (uncertain"),
        cluster_rows = FALSE, cluster_columns = FALSE,
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 7),
        col = col_fun(seq(0,1,0.5)))
dev.off()

# Correct metadata
i39 <- metadata$Sample == "039"
i38 <- metadata$Sample == "038"
metadata[i39, "Sample"] <- "038"
metadata[i38, "Sample"] <- "039"

write.csv(metadata, paste0(opt$outdir,"/metadata_corrected.csv"))

sessionInfo()