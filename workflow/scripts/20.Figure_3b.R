# Script to produce figure 3b

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure3b.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("qtl2")
  library("openxlsx")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--qtl"), type="character", default=NULL, 
              help="Diretory for genetic information", metavar="character"),
  make_option(c("-b", "--models"), type="character", default=NULL, 
              help="Directory for ATAC QTL models", metavar="character"),
  make_option(c("-c", "--S2"), type="character", default=NULL, 
              help="Table S2", metavar="character"),
  make_option(c("-d", "--S6"), type="character", default=NULL, 
              help="Table S6", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$qtl)){
  print_help(opt_parser)
  stop("Diretory for genetic information (-a) is missing", call.=FALSE)
}
if (is.null(opt$models)){
  print_help(opt_parser)
  stop("Directory for ATAC QTL models (-b) is missing", call.=FALSE)
}
if (is.null(opt$S2)){
  print_help(opt_parser)
  stop("Table S2 (-c) is missing", call.=FALSE)
}
if (is.null(opt$S6)){
  print_help(opt_parser)
  stop("Table S6 (-d) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

get_perms <- function(files) {
  perms <- lapply(files, readRDS)
  c <- class(perms[[1]])
  perms <- do.call(rbind, perms)

  class(perms) <- c
  return(perms)
}

cat("Load Table S6\n")
atac_QTL <- read.xlsx(opt$S6, "ATAC_QTL")
rna_QTL <-  read.xlsx(opt$S6, "RNA_QTL")
cat("Load Table S2\n")
atac_diff <- read.xlsx(opt$S2, "Differential_Accessibility")

cat("Select examples\n")
genes <- c('Nrf1', 'Ncald', "Nav1")
genes <- genes[genes %in% rna_QTL$Gene]

regions <- sapply(genes, function(g) {
  region <- atac_QTL[atac_QTL$gene ==g,"regionID"]
  
  if (any(atac_diff$Region_ID %in% region)) {
    region <- atac_diff$Region_ID[atac_diff$Region_ID %in% region][1]
  } else {
    region <- region[1]
  }
  return(region)
  })

cat("Load genomic information\n")
atac_data <- read_cross2(file =  paste0(opt$qtl,"/atacchr6.json"))
rna_data <- read_cross2(file =  paste0(opt$qtl,"/rna.json"))

out_rna <- readRDS(paste0(opt$models,"/rna_QTL_model.RData"))
perms_rna <- list.files(opt$models, pattern="rna_.*_perm", full.names=TRUE)
perms_rna <- get_perms(perms_rna)

atac_files <- list.files(opt$models, pattern="atac.*_model.RData", full.names=TRUE)
atac_files <- grep('FC_QTL', atac_files, invert = TRUE, value=TRUE)
atac_list <- lapply(atac_files, function(file) {
  f <- readRDS(file)
  return(colnames(f[['model']]))
  }  )
names(atac_list) <- atac_files


for (i in 1:length(genes)) {
  gene <- genes[i]
  region <- regions[gene]
  
  cat("Load QTL information for ATAC\n")
  for (file in 1:length(atac_list)){
    if(region %in% atac_list[[file]]) {
      out_atac <- readRDS(names(atac_list)[file])

      perms <- basename(gsub('model.RData','perm.*RData', names(atac_list)[file]))

      perms_atac <- list.files(opt$models, pattern=perms, full.names=TRUE)
      perms_atac <- get_perms(perms_atac)

      break
    }
  }
  
  atac_map <- insert_pseudomarkers(map=atac_data$gmap, step=1)
  rna_map <- insert_pseudomarkers(map=rna_data$gmap, step=1)
  
  atac_idx <- which(colnames(out_atac$model) == region)
  rna_idx <- which(colnames(out_rna$model) == gene)
  
  same <- atac_QTL[atac_QTL$regionID == region,'snp'] %in% rna_QTL[rna_QTL$Gene == gene,'snp'] 
  same <- ifelse(same,'same','')
  
  svg(paste0(opt$outdir, "/Fig3b",i,same,".svg"), height=14, width=10)
  par(mar=c(5.1, 4.1, 1.1, 1.1), mfrow=c(2,1))
  plot(out_atac$model, atac_map, lodcolumn=atac_idx, col="violetred")
  abline(h=summary(perms_atac, alpha = 0.05)[['A']][atac_idx],lwd=2)
  legend("topright", lwd=2, cex=1.5, col="violetred", bg="gray90",unique(atac_QTL[atac_QTL$regionID == colnames(out_atac$model)[atac_idx],'region']))
  
  plot(out_rna$model, rna_map, lodcolumn=rna_idx, col="purple")
  abline(h=summary(perms_rna, alpha = 0.05)[['A']][rna_idx],lwd=2)
  legend("topright", lwd=2, cex=1.5, col="purple", bg="gray90",colnames(out_rna$model)[rna_idx])
  dev.off()
}

sessionInfo()