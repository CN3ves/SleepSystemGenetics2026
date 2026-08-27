# Script designed to aggregame and filter all QTL results (split due to unwieldy number of permutations)

# Redirect all R logs to Snakemake log
log <- file('logs/9-BXD_qtl/aggregate_9.6.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-t", "--table"), type="character", default=NULL, 
              help="QTL results tables", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$table)){
  print_help(opt_parser)
  stop("QTL results tables (-t) is missing",call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
table_files <- strsplit(opt$table, " ")[[1]]

exps <- c('atac', 'rna', 'sleep')

for (exp in exps) { 

 # Merge QTL table
  exp_files <- grep(exp, table_files, value=TRUE)
  fc <- 1
  if (any(grepl('FC',grep(exp, table_files, value=TRUE)))) fc <- 2 # second iteration looks of FC files
  
  files <- grep('FC', exp_files, value=TRUE, invert=TRUE)
  
  for(tab in 1:fc) {
    print("Merge files")
    print(files)
    tables <- lapply(files, read.csv, row.names=1)
    tables  <- do.call(rbind, tables)

    # Merge alpha thresholds
    qtls <- apply(tables[,c("lodcolumn", "chr","snp")],1, paste, collapse='/')

    dups <- which(duplicated(qtls))
    for (dup in dups) {
      idx <- which(qtls == qtls[dup])       
      stopifnot(length(idx)==2)
      tables[idx[1],] <- sapply(tables[idx,], function(x) paste(unique(x), collapse = '|-|'))
    }
    tables$test <- gsub('additive_','', tables$test)
    tables$test <- gsub('interaction_','', tables$test)
    tables <- tables[-dups,]

    #Save
    f <- paste0(opt$outdir, "/QTL_table_", exp,".csv")
    if(tab ==2) f <- gsub('.csv','FC.csv',f)
    write.csv(tables, f)

    print("Filtering table")
    if (exp != "sleep" & tab == 1) {
      genos <- tables[,c("BB_Ctrl", "BB_SD", "DD_Ctrl", "DD_SD")]
    } else {
      genos <- tables[,c("BB", "DD")]
    }
    genos[] <- lapply(genos, as.numeric)
    print(paste("Check merging tables - NAs:", sum(is.na(genos))))
    samples <- rowSums(genos)

    idx <- rowSums(genos < 5) == 0 

    n <- gsub("QTL_table_","QTL_table_filter_", f)

    write.csv(tables[idx,], n)

    print(paste0(basename(f), " Original size ", nrow(tables), " Filtered size ", nrow(tables[idx,])))
    files <- grep('FC', exp_files, value=TRUE, invert=FALSE) 
  }
  
}

sessionInfo()