# Script to produce figure S1a

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure1Sa.log', open = "wt")
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
  make_option(c("-a", "--dir"), type="character", default=NULL, 
              help="Called variants", metavar="character"),
  make_option(c("-b", "--map"), type="character", default=NULL, 
              help="BXD genetic map directory", metavar="character"),
  make_option(c("-c", "--meta"), type="character", default=NULL, 
              help="ATAC-seq metadata", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("Called variants (-a) is missing", call.=FALSE)
}
if (is.null(opt$map)){
  print_help(opt_parser)
  stop("BXD genetic map directory (-b) is missing", call.=FALSE)
}
if (is.null(opt$meta)){
  print_help(opt_parser)
  stop("ATAC-seq metadata (-c) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load Genotyping data\n")
# Read inputs
bxd <- readRDS(paste0(opt$map,"/genotypes.Rdata"))
geno_files <- list.files(opt$dir, pattern="csv", full.names=T)
metadata <- read.csv(opt$meta, row.names=1)
metadata$Sample <- sprintf("%03d",metadata$Sample)
bxd_snps <- read.csv(paste0(opt$map,"/bxd_pmap.csv"), skip=3)

genos <- lapply(geno_files, read.csv, row.names=1)

names(genos) <- gsub(".*/(.*)_full.csv","\\1", geno_files)

vars <- sort(unique(unlist(lapply(genos,names))))
vars <- as.data.frame(matrix(rep(NA,length(vars)*length(genos)),
                   nrow=length(vars),ncol=length(genos), 
                   dimnames=list(vars,names(genos))))

for(g in names(genos)) {
  vars[,g] <- as.character(genos[[g]][rownames(vars)])
}

strains <- metadata$Strain[metadata$Sample %in% gsub('_.*','',names(genos))]

cat("Get parent strain genotypes\n")
call_allele <- function(mat) {
  uniq_mat <- apply(mat, 1, function(x) {
    x <- gsub(",<\\*>", "",x)
    x <- x[!is.na(x) & x != "NA"]
    x <- x[x != "<*>"]
    x <- gsub(",", "",x)
    x <- unique(x)
    if(length(x) !=1 ) x <- "X"
    if(nchar(x) > 1) x <- "X"
    
    return(x)
  })
  names(uniq_mat) <- rownames(mat)
  return(uniq_mat)
}

b_allele <- call_allele(vars[,strains=="C57Bl6"])
d_allele <- call_allele(vars[,strains=="DBA",])

no_info_snps <- which((b_allele == d_allele) | (b_allele == "X") | ( d_allele == "X") )
b_allele <- b_allele[-no_info_snps]
d_allele <- d_allele[-no_info_snps]
references = list('b6' = b_allele, 'dba'=d_allele)

bxd <- bxd[unique(c(names(b_allele), names(d_allele))),grep("BXD", names(bxd))]
samples <- vars[rownames(bxd),!strains %in% c("C57Bl6", "DBA")]
names(samples) <- gsub("_.*","",names(samples))

cat("Match BXD lines to parents\n")
classify <- function(samples, genotypes, metadata,references) {
 
  # Prepare object for storing call genotypes
  res <- metadata[metadata$Sample %in% gsub('_.*','',names(samples)), c("Sample", "Strain")]
  res$assign <- NA
  rownames(res) <- res$Sample
  res$Sample <- NULL
  res$Strain[grep("BXD", res$Strain, invert=TRUE)] <- paste0("BXD", res$Strain[grep("BXD", res$Strain, invert=TRUE)] )

  # Prepare object to store percentage of matching
  mat <- as.data.frame(matrix(NA, nrow=nrow(res), ncol=ncol(genotypes), dimnames = list(rownames(res),names(genotypes))))
  
  for (sample in names(samples)) {
    alleles <- call_allele(samples[,sample, drop=FALSE])
    alleles <- alleles[rownames(genotypes)]
    alleles <- alleles[!is.na(alleles)]
    for (a in names(alleles)) {
      if (length(alleles[a]) > 0) {
        b <- alleles[a]== references[['b6']][a] 
        d <- alleles[a]== references[['dba']][a] 
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
    
    strains <- sapply(genotypes[names(alleles),], function(x) {
      hs <- x == "H" | as.character(alleles) == "H" | is.na(alleles)
      mean(x[!hs] == as.character(alleles)[!hs])
      })
    
    res[sample,"assign"] <- names(which.max(strains))
    res[sample,"score"] <- max(strains)
    mat[sample,] <- strains
  }
  
  return(list("assign"=res, "match"=mat))

}

geno_res <- classify(samples, bxd, metadata, references)

# Genotype call plots
M <- geno_res[["match"]]
rownames(M) <- paste0(rownames(M), " (", metadata$Strain[match(rownames(M),metadata$Sample)],")")
rownames(M)[grep("BXD", rownames(M), invert = TRUE)] <- gsub(" \\(", " (BXD", rownames(M)[grep("BXD", rownames(M), invert = TRUE)] )

M <- M[order(gsub(".* ","",rownames(M))),order(colnames(M))]
  
col_fun = colorRamp2(c(0,0.5,1), c("red",  "white", "blue"))

named_rows <- c(seq(1,nrow(M),5),grep("03[89]",rownames(M)))
named_rows <- rownames(M)[named_rows] 

cat("Save plot\n")
svg(paste0(opt$outdir,"/FigS1a.svg"), width=10, height=10)  
Heatmap(as.matrix(M), border_gp = gpar(col = "black", lty = 1),
        column_title = "BXD Genotyping",
        heatmap_legend_param = list(title = "% SNPs matching reference"),
        cluster_rows = FALSE, cluster_columns = FALSE,
        column_names_gp = grid::gpar(fontsize = 12),
        row_names_gp = grid::gpar(fontsize = 12),
        row_labels = sapply(rownames(M), function(r) ifelse(r %in% named_rows, r,"")),
        col = col_fun(seq(0,1,0.5)))
dev.off()

write.csv(as.matrix(M), paste0(opt$outdir, "/data/S1a.csv"))

sessionInfo()