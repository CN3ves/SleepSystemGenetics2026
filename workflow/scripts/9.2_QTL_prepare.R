# Script designed to produce the files required to run the QTL analyses

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("rjson")
  library("edgeR")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene") 
  library("org.Mm.eg.db")
  library("optparse")
})
source("workflow/scripts/9.0_QTL_helper.R")

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="Rdata file with the normalized feature counts", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="Directory with the bxd template files", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("Rdata file with the normalized feature counts (-c) is missing", call.=FALSE)
}
if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("Directory with the bxd template files (-d) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

# Redirect all R logs to Snakemake log
exp <- gsub('_count.*','',basename(opt$counts))

log <- file(paste0('logs/9-BXD_qtl/data_9.2_', exp,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Reading inputs\n")
template <- read_template(opt$dir)
counts <- readRDS(opt$counts)

# Make phenotype
metadata <- counts$sample
if(exp=='sleep') {
  counts <- counts$counts
} else {
  counts <- log(cpm(counts))
}

colnames(counts) <- metadata[colnames(counts),"group"]
colnames(counts) <- gsub("BXD0*","BXD", colnames(counts)) # match strain names
colnames(counts) <-  gsub("B6([1-2])", "C57\\1", colnames(counts))
colnames(counts) <-  gsub("DB([1-2])", "DBA\\1", colnames(counts))
colnames(counts) <-  gsub("BXD96", "BXD48a", colnames(counts))
colnames(counts) <-  gsub("BXD97", "BXD65a", colnames(counts))
colnames(counts) <-  gsub("BXD103", "BXD73b", colnames(counts))

# Keep valid strains
strains <- gsub("_.*", "", colnames(counts)) 
strains_idx <- strains %in% names(template[["genotype"]][["data"]])
print(paste("Ignored strains/samples for", exp, ":", paste(sort(unique(strains[!strains_idx])), collapse = ", ")))

# Remove ignored samples
counts <- counts[,strains_idx]
strains <- strains[strains_idx]

# Generate files for QTL analysis
geno <- make_genotype(template[["genotype"]], select = strains, sample_names = colnames(counts), file = paste0(opt$outdir,"/", exp, "_geno.csv"))

cross <- make_cross(genotype = names(geno)[-1], strains = names(geno)[-1], template = template[["cross"]], file=paste0(opt$outdir,"/", exp, "_crossinfo.csv"))

if(exp == "atac") {
  subset = unique(gsub("_.*","",rownames(counts))) # split by chromossome
} else {
  subset = "skip"
}

if (exp != 'sleep') {#Get FC values
    geno_fc <- make_genotype(template[["genotype"]], select = unique(strains), sample_names = unique(strains), file = paste0(opt$outdir,"/", exp, "_genoFC.csv"))
    cross_fc <- make_cross(genotype = names(geno_fc)[-1], strains = names(geno_fc)[-1], template = template[["cross"]], file=paste0(opt$outdir,"/", exp, "_crossinfoFC.csv"))
}

for (chr in subset) {
  
  if(chr != "skip"){
    print(paste("Subset:", chr))
    idx <- gsub("_.*","",rownames(counts)) == chr
    mat <- counts[idx,]
    # average replicates to match other data types
    groups <- unique(colnames(mat)) 
    for (g in groups) {
      dups <- grep(g, colnames(mat))
      tmp <- mat[,dups]
      stopifnot(length(unique(colnames(tmp))) == 1)
      
      tmp <- rowMeans(tmp)[rownames(mat)] # add mean to first rep
      tmp <- matrix(tmp, ncol=1, dimnames=list(names(tmp), g))
      mat <- mat[,-dups] # remove other replicates
      
      mat <- cbind(mat, tmp)

    }
     
  } else {
    mat <- counts
    chr = ""
  }
  pheno <- make_phenotype(counts = mat, template= template[["phenotype"]], file = paste0(opt$outdir,"/", exp, "_pheno",chr, ".csv"))
  pheno$id <- gsub("\\.[0-9]", "", pheno$id)
  
  regions <- get_regions(exp, genes = names(pheno)[-1])                      
  covar <- make_covar(phenotype=names(pheno)[-1], regions = regions, template = template[["covariates"]] , file = paste0(opt$outdir,"/", exp, "_phenocovar",chr, ".csv"))
  
  # Write control file
  json <- template[['control']]
  
  json$geno <- paste0(exp, "_geno.csv")
  json$cross_info$file <- paste0(exp, "_crossinfo.csv")
  json$pheno <-paste0(exp, "_pheno",chr, ".csv")
  json$phenocovar <-  paste0(exp, "_phenocovar",chr, ".csv")
  
  write(toJSON(json),  paste0(opt$outdir,"/", exp, chr, ".json"))

  if (exp != 'sleep') {#Get FC values
    print("Getting Fold change values for interaction model")
    mat <- as.data.frame(mat)
    for (s in unique(strains)) {
      idx <- grep(paste0('^',s,'_'),colnames(mat))
 
      fc <- exp(mat[,idx])
      fc <- rowMeans(exp(fc[,grep('SD',colnames(fc)),drop=FALSE])) / rowMeans(exp(fc[,grep('CTRL',colnames(fc)),drop=FALSE]))
  
      mat <-  mat[,-idx]
      mat[,s] <- log(fc)
    }

    pheno_fc <- make_phenotype(counts = mat, template= template[["phenotype"]], file = paste0(opt$outdir,"/", exp, "_phenoFC",chr, ".csv"))
    pheno_fc$id <- gsub("\\.[0-9]", "", pheno_fc$id)
  
    regions_fc <- get_regions(exp, genes = names(pheno_fc)[-1])                      
    covar_fc <- make_covar(phenotype=names(pheno_fc)[-1], regions = regions_fc, template = template[["covariates"]] , file = paste0(opt$outdir,"/", exp, "_phenocovarFC",chr, ".csv"))
  
    # Write control file
    json_fc <- template[['control']]
  
    json_fc$geno <- paste0(exp, "_genoFC.csv")
    json_fc$cross_info$file <- paste0(exp, "_crossinfoFC.csv")
    json_fc$pheno <-paste0(exp, "_phenoFC",chr, ".csv")
    json_fc$phenocovar <-  paste0(exp, "_phenocovarFC",chr, ".csv")
  
    write(toJSON(json_fc),  paste0(opt$outdir,"/", exp, chr, "FC.json"))
  }
}


sessionInfo()