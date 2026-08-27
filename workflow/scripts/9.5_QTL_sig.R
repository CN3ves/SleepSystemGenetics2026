# Script designed to gather the main results and permutation value to find statistical significant QTLs

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("qtl2")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-m", "--model"), type="character", default=NULL, 
              help="RData file containignnthe QTL model", metavar="character"),
  make_option(c("-p", "--perms"), type="character", default=NULL, 
              help="RData file containing the permutation results", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$model)){
  print_help(opt_parser)
  stop("RData file containignnthe QTL model (-m) is missing", call.=FALSE)
}
if (is.null(opt$perms)){
  print_help(opt_parser)
  stop("RData file containing the permutation result (-p) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
exp <- gsub('_.*','',basename(opt$model))
qtl_data <- readRDS(opt$model)
map <- insert_pseudomarkers(map=qtl_data[['qtl_data']]$gmap, step=1)

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/9-BXD_qtl/sigs_9.5_', exp,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Reading QTL models and permutations\n")
qtl_fit <- qtl_data[['model']]
perm <- strsplit(opt$perms, " ")[[1]]
perm <- lapply(perm,readRDS)
perm <- do.call(rbind, perm)

step <- 1000
alpha <- c(0.05,0.1)
ci <- 0.95
stopifnot(ci==(1-alpha[1]))

# Finding LOD peaks
add_test <- function(test, name,treatment,data) {#name the test for the results table
  map <- insert_pseudomarkers(map=data$gmap, step=1) # spaced at 1 cM intervals 
  pr <- calc_genoprob(cross=data, map=map, error_prob=0.0001)

    if (nrow(test) !=0) {
      test$test <- name
      
      genos <- apply(test[, c("chr", "pos")], 1, function(x){
        # maxarg has randomness and does not produce the same result
        # given the same input for differen alpha runs in same cases!
        # force seed
        set.seed(ceiling(as.numeric(x["pos"])))
        geno=maxmarg(pr, map, chr=x["chr"], pos=as.numeric(x["pos"]), return_char=TRUE, minprob=0.5)
        
        stopifnot(all(treatment %in% names(geno)))
        geno <- geno[treatment]
        treatment =  gsub(".*_","",treatment)
        tab = table(treatment, geno)
        
        if (any(grepl("SD",names(geno)))) {
          data.frame("BB_Ctrl" = tab["CTRL","BB"], "BB_SD" = tab["SD","BB"], "DD_Ctrl" = tab["CTRL","DD"], "DD_SD" = tab["SD","DD"])
        } else {
          tab = colSums(tab)
          data.frame("BB" = tab["BB"], "DD" = tab["DD"])
        }
        
        })
      test <- cbind(test,do.call(rbind,genos))
    } else {
      test$test <- character()
    }
  rownames(test) <- NULL
  return(test)
}

convert_pos <- function(qtls, rqtl_object) { # convert from morgan (gmap) to bp (pmap) so it's useful
  
  if (nrow(qtls) !=0) {
    qtls$snps <- qtls$bp_hi <-  qtls$bp_lo <- qtls$snp <- qtls$bp<- 0
    for (r in 1:nrow(qtls)) {
      chr <- qtls$chr[r]
    
      #Get closest marker from map (becsuse pseudomakers are not in pmap)
      pos <- which.min(sapply(rqtl_object$gmap[[chr]], function(x) abs(x-qtls[r,"pos"])))
      lo <- which.min(sapply(rqtl_object$gmap[[chr]], function(x) abs(x-qtls[r,"ci_lo"])))
      hi <- which.min(sapply(rqtl_object$gmap[[chr]], function(x) abs(x-qtls[r,"ci_hi"])))
    
      qtls$bp[r] <-  rqtl_object$pmap[[chr]][names(pos)]*1000000
      qtls$bp_lo[r] <- rqtl_object$pmap[[chr]][names(lo)]*1000000
      qtls$bp_hi[r] <- rqtl_object$pmap[[chr]][names(hi)]*1000000
      
      qtls$snp[r] <- names(pos)
      qtls$snps[r] <- paste0(names(rqtl_object$pmap[[chr]][lo:hi]), collapse=";" )
      
      #sanity check
      if (qtls$bp_lo[r] > qtls$bp[r]) print("Lower boundary is after best location!")
      if (qtls$bp_hi[r] < qtls$bp[r]) print("higher boundary is after best location!")
    }
  } else {
     qtls$bp <- qtls$bp_hi <-  qtls$bp_lo <- character()
  }
  
 return(qtls)
}
test <- ifelse(grepl('FC',opt$model), "interaction", "additive")
print(paste("Finding peaks for", test, "model:", exp))

res <- list()
for (a in alpha) {
  cols <- ncol(perm[['A']])
  rows <- nrow(perm[['A']])
  print(paste("Alpha threshold", a, "features:", cols, "permutations:", rows))
  
  peaks <- lapply(seq(1,cols,step), function(i) {
    p <- perm
    #legacy but where splitting the permutation to decrease RAM lead to duplicated columns
    #should be fixed now on line 74 of 9.4_QTL_perms.R
    #with pheno = qtl_data[['pheno']][,(ceiling(n/2)+1):n],
    #but re-running takes a long time 
    p[['A']] <- p[['A']][,!duplicated(colnames(p[['A']]))]
    p[['X']] <- p[['X']][,!duplicated(colnames(p[['X']]))]
    cols <- ncol(p[['A']])

    p[['A']] <- p[['A']][,i:min(i+step-1,cols)]
    p[['X']] <- p[['X']][,i:min(i+step-1,cols)]
    class(p) <- class(perm)

    ths <- summary(p, alpha = a)
    fit <- qtl_fit[,i:min(i+step-1,cols)]

    stopifnot(all(colnames(fit) == colnames(ths[['A']])))

    qtls_A <- find_peaks(scan1_output = fit, map = map, threshold = ths[['A']], 
      prob = ci, peakdrop = 1.8, expand2markers = TRUE) # Autossomes
    qtls_A <- qtls_A[qtls_A$chr!="X",]
    qtls_A <- add_test(qtls_A, paste0(test,'_alpha:', a), rownames(qtl_data[['addcovar']]), data =  qtl_data[['qtl_data']])
    qtls_A <- convert_pos(qtls_A,  qtl_data[['qtl_data']])

    qtls_X <- find_peaks(scan1_output = fit, map = map,threshold = ths[['X']], 
    prob = ci, peakdrop = 1.8, expand2markers = TRUE)  # X chr
    qtls_X <- qtls_X[qtls_X$chr=="X",]
    qtls_X <-  add_test(qtls_X, paste0(test,'_alpha:', a), rownames(qtl_data[['addcovar']]), data =  qtl_data[['qtl_data']])
    qtls_X <- convert_pos(qtls_X,qtl_data[['qtl_data']])

    return(rbind(qtls_A, qtls_X))

  })

  res[[paste('alpha',a)]] <- do.call(rbind, peaks)
}

res <- unique(do.call(rbind,res))

print("Saving results")
write.csv(res, paste0(opt$outdir,"/",exp, "_QTL_sigs.csv"))

sessionInfo()