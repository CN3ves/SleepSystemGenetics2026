# Script designed to run post-hoc t-test to identify QTL baseline and SD response trends

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("foreach")
  library("doParallel")
  library("doSNOW")
  library("tidyverse")
  library("qtl2")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-q", "--qtl"), type="character", default=NULL, 
              help="Folder containing the QTL models", metavar="character"),
  make_option(c("-t", "--table"), type="character", default=NULL, 
              help="Interaction QTL results table", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$qtl)){
  print_help(opt_parser)
  stop("Folder containing the QTL models (-q) is missing", call.=FALSE)
}
if (is.null(opt$table)){
  print_help(opt_parser)
  stop("Interaction QTL results table (-t) is missing",call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

exp <- gsub(".*filter_(.*).csv","\\1",opt$table)
exp <- gsub("FC","",exp)

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/9-BXD_qtl/classify_9.7_',exp,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Reading inputs\n")
test <- read.csv(opt$table,row.names=1)
models <- list.files(opt$qtl, pattern = 'model', full.names=TRUE)

cat("Running t-test for each significant interaction \n")
cl <- makeCluster(20)
registerDoParallel(cl)
paths <- .libPaths()
clusterExport(cl, "paths")

output <- foreach(i = 1:nrow(test), .packages = c("qtl2")) %dopar% {
  region = test$lodcolumn[i]

  # open the data for model with no interatin (contains treatment information)
  if(grepl('chr[0-9XY]{1,2}',region)) {
    chr <- gsub('_.*','',region)
    model <- grep(paste0(chr,'_'), models, value = TRUE)
  } else {
    model <- grep(paste0(exp,'_'), models, value = TRUE)
  }

  geno_info <- readRDS(model)[['qtl_data']]

  map <- insert_pseudomarkers(map=geno_info$gmap, step=1) # spaced at 1 cM intervals 
  pr <- calc_genoprob(cross=geno_info, map=map, error_prob=0.0001)
  
  genotype <-  maxmarg(pr, map, chr=test$chr[i], pos=test$pos[i], return_char=TRUE)
  phenotype <- geno_info$pheno[,region]
  genotype <- genotype[gsub("\\..*","", names(phenotype))]
  
  d <- data.frame(geno=genotype, val = phenotype, treat = gsub("\\..*","", gsub(".*_","",names(phenotype))))
  d$group <- paste(d$geno, d$treat)
  infs <- grep('Inf',d$val) #Fold changes can have infinit values
  if (length(infs) > 0) d <- d[-infs,]

  catch <- function(exp) {
    # function to catch any error in t.test (not enough samples)
    tryCatch({
      exp
    },error = function(e) {
      list(estimate=c(NA,NA), p.val=NA)
    })
  }


  bl <- catch(t.test(d$val[d$group== "BB CTRL"], d$val[d$group=="DD CTRL"]))
  bb <- catch(t.test(d$val[d$group== "BB CTRL"], d$val[d$group=="BB SD"]))
  dd <- catch(t.test(d$val[d$group== "DD CTRL"], d$val[d$group=="DD SD"]))
  res <- data.frame("pval" = paste0(c(bl$p.val,bb$p.val,dd$p.val), ifelse(length(infs)==0, '', ' INF')), 
                    "sig" = c(bl$p.val<0.05,bb$p.val<0.05,dd$p.val<0.05), 
                    "fc"=c(bl$estimate[2]/bl$estimate[1],
                           bb$estimate[2]/bb$estimate[1],
                           dd$estimate[2]/dd$estimate[1]),
                    "sign"=c(ifelse(bl$estimate[2]>bl$estimate[1],"up","down"),
                             ifelse(bb$estimate[2]>bb$estimate[1],"up","down"),
                             ifelse(dd$estimate[2]>dd$estimate[1],"up","down")),
                    "region" = region,
                    "chr" = test$chr[i],
                    "pos" = test$pos[i],
                    "snp" = test$snp[i],
                    row.names= c("baseline-DDvsBB", "BB-SDvsCTRL", "DD-SDvsCTRL"))
  
  res$group <- sapply(1:nrow(res), function(i) ifelse(res$sig[i], paste0(rownames(res)[i],":",res$sign[i]), paste0(rownames(res)[i],":ns")))
  res$INVERTED <- ifelse(all(res$pval[2:3] < 0.05) & (res$sign[2] != res$sign[3]) ,TRUE,FALSE)
  
  return(res)
}

stopCluster(cl)

print(paste("Any value is NULL?",  which(sapply(output, is.null))))
saveRDS(output,paste0(opt$outdir,"/qtl_ttests_",exp,".RData"))

tab <- do.call(rbind,output)

tab <- tab %>% select(-c(fc,sig,INVERTED,pval, sign)) %>% mutate(groups=gsub(":.*","",group)) %>% pivot_wider(names_from=groups, values_from=group)

tab <- lapply(tab, function(x) gsub(".*:", "",x))
tab <- as.data.frame(tab)

test <- test[,-1]
names(test)[1] <- "region"

tab <- merge(tab, test, by=c("region","chr","pos","snp"))

write.csv(tab, paste0(opt$outdir,"/qtl_ttests_",exp,".csv"))

sessionInfo()