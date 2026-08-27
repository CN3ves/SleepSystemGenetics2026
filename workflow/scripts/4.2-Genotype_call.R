# Script designed to call variants from the ATAC-seq reads around BXD variants

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("doParallel")
  library("foreach")
  library("optparse")
})

cat("Checking arguments\n")
option_list = list(
  make_option(c("-b", "--bam"), type="character", default=NULL, 
              help="Bam file to genotype", metavar="character"),
  make_option(c("-g", "--geno"), type="character", default=NULL, 
              help="Rdata file with the genotype object", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bam)){
  print_help(opt_parser)
  stop("Bam file to genotype (-b) is missing", call.=FALSE)
}
if (is.null(opt$geno)){
  print_help(opt_parser)
  stop("Rdata file with the genotype object (-g) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
bamfile <- opt$bam
bxd <- readRDS(opt$geno)
bname <- gsub('\\..*','',basename(bamfile))
print("Reading bam file")

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/4-BXD_genotype/call_', bname,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Get variants from experimental data
cat("Start computation\n")
cl <- makeCluster(20)
registerDoParallel(cl)
clusterEvalQ(cl, library(doParallel, lib='Rlibs/'))

var_lst <- foreach(row = 1:nrow(bxd)) %dopar% {
  chr <- bxd[row,"chr"]
  pos <- bxd[row,"pos"]
  filter <- paste0("chr",chr,":", pos-10000, "-", pos+10000)
  
  temp <- gsub(".co.bam", "temp", basename(bamfile))
 
  # Call samtools and bcftools to get variant is it exists in the data
  system(paste0("/dcsrsoft/spack/20260114/spack/opt/spack/linux-zen2/samtools-1.21-wgap6pkyp3icy7jh5yu3ef4m3yycbu3l/bin/samtools view -h ", bamfile, " '", filter, "' | /dcsrsoft/spack/20260114/spack/opt/spack/linux-zen2/samtools-1.21-wgap6pkyp3icy7jh5yu3ef4m3yycbu3l/bin/samtools view -bS - > ", temp,row,".bam"))
  var <-  system(paste0("/dcsrsoft/spack/20260114/spack/opt/spack/linux-zen2/bcftools-1.22-wyegavpf2nyim3vxlyxjt554cm3akbcg/bin/bcftools mpileup --no-reference ", temp,row,".bam | grep ", pos, " | cut -f5"), intern=TRUE) 
  
  #Remove temp file
  system(paste0("rm ",temp, row,".bam"))
  
  if (length(var)>0) {
    return(var)
  } else {
    return(NA)
  }
}

stopCluster(cl)
cat("Computation ended\n")
names(var_lst) <- rownames(bxd)

cat("Save table\n")
file <- gsub(".*/", "", bamfile)
file <- gsub(".co.bam", ".csv", file)
file <- paste0(opt$outdir, "/",file)
write.csv(as.data.frame(var_lst), file)

sessionInfo()