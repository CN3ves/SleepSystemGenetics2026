# Script designed to runnpermutations of the qtl model

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("qtl2")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-q", "--qtl"), type="character", default=NULL, 
              help="Rdata containing the QTL model", metavar="character"),
  make_option(c("-i", "--index"), type="character", default=NULL, 
              help="Rdata file with the normalized feature counts", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$index)){
  print_help(opt_parser)
  stop("Feature index (-i) is missing", call.=FALSE)
}
if (is.null(opt$qtl)){
  print_help(opt_parser)
  stop("Rdata containing the QTL model (-q) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
qtl_file <- opt$qtl
index <- as.numeric(opt$index)
exp <- gsub('_QTL.*','',basename(qtl_file))
perms_n <-100

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/9-BXD_qtl/qtl_9.4_', exp,'_',index,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat(paste("Reading ",qtl_file,"\n"))
qtl_data <- readRDS(qtl_file)
map <- insert_pseudomarkers(map=qtl_data[['qtl_data']]$gmap, step=1) # spaced at 1 cM intervals 



#Performing a permutation test
print(paste("Running permutations on", exp))
seed <- nchar(exp)+index
print(paste("Using index", index," - seed =", seed))
print(paste(perms_n*index, "permutations of 1500:", round(perms_n*index / 1500*100,2), "%"))

set.seed(seed)
n <- ncol(qtl_data[['pheno']])
print(paste("Run first half of features:",1,'-',ceiling(n/2) ,"out of", n))

perm <- scan1perm(genoprobs = qtl_data[['genoprobs']], 
                            pheno = qtl_data[['pheno']][,1:ceiling(n/2)], 
                            addcovar = qtl_data[['addcovar']],  
                            kinship=qtl_data[['kinship']],
                            perm_Xsp=TRUE,
                            chr_lengths=chr_lengths(map),
                            n_perm = perms_n, 
                            cores=10)
print(paste("Run second half of features:",ceiling(n/2),'-',n ,"out of", n))

perm2 <- scan1perm(genoprobs = qtl_data[['genoprobs']], 
                            pheno = qtl_data[['pheno']][,(ceiling(n/2)+1):n], 
                            addcovar = qtl_data[['addcovar']],  
                            kinship=qtl_data[['kinship']],
                            perm_Xsp=TRUE,
                            chr_lengths=chr_lengths(map),
                            n_perm = perms_n, 
                            cores=10)

perm_addictive <- cbind(perm,perm2)

print(paste("Merge both halves:", sapply(perm_addictive,dim)[1,1], 'permutations of', sapply(perm_addictive,dim)[2,1],'features'))

print("Saving results")
file <- paste0(opt$outdir,"/",exp, "_QTL_perm",index,".RData")

saveRDS(perm_addictive, file)

sessionInfo()