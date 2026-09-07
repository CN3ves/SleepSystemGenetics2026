# Script to produce table S6

# Redirect all R logs to Snakemake log
log <- file('logs/19-Tables/S6_19.6.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("openxlsx")
  library("GenomicRanges")
  library("tidyverse")
  library("optparse")
})

cat("Checking arguments\n")

option_list = list(
  make_option(c("-a", "--qtl"), type="character", default=NULL, 
              help="Directory with QTL results tables", metavar="character"),
  make_option(c("-b", "--int"), type="character", default=NULL, 
              help="Directory with the interaction trend results", metavar="character"),
  make_option(c("-c", "--sleep"), type="character", default=NULL, 
              help="Directory with phenotype metadata", metavar="character"),
  make_option(c("-d", "--S2"), type="character", default=NULL, 
              help="Table S2", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$qtl)){
  print_help(opt_parser)
  stop("Directory with QTL results tables (-a) is missing", call.=FALSE)
}
if (is.null(opt$int)){
  print_help(opt_parser)
  stop("Directory with the interaction trend results (-b) is missing", call.=FALSE)
}
if (is.null(opt$sleep)){
  print_help(opt_parser)
  stop("Directory with phenotype metadata (-c) is missing", call.=FALSE)
}
if (is.null(opt$S2)){
  print_help(opt_parser)
  stop("Table S2 (-d) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

# Create a new workbook and add a sheet
wb <- createWorkbook()

cat("Load QTL results\n")
files <- list.files(opt$qtl, pattern = ".csv", full.name=TRUE)
files <- grep("filter", files, value=TRUE,invert = TRUE)
qtls <- lapply(files,read.csv)
names(qtls) <- files

inv <- list.files(opt$int, pattern = "RData", recursive = TRUE,full.name=TRUE)
types <- lapply(inv,readRDS)
names(types) <- inv

# Helper functions
classify_invertions <- function(tab) {
  #Load and format data
  tab <- do.call(rbind,tab)
  tab$group <- paste(tab$group, round(tab$fc,2), sep = "_")
  tab <- tab %>% select(-c(fc,sig,INVERTED,pval, sign)) %>% 
    mutate(results= gsub(".*:","",group), groups=gsub(":.*","",group)) %>%
    filter(groups != 'NA_NA') %>%  
    select(-group) %>% 
    pivot_wider(names_from=groups, values_from=results)
  
  tab <- as.data.frame(tab)
  
  # Label test results
  tab$baseline <- !grepl("ns",tab$'baseline-DDvsBB')
  tab$inverted <- (grepl("down",tab$'BB-SDvsCTRL') & grepl("up",tab$'DD-SDvsCTRL')) |
                  (grepl("up",tab$'BB-SDvsCTRL') & grepl("down",tab$'DD-SDvsCTRL'))
  tab$protected <- (grepl("ns",tab$'BB-SDvsCTRL') & !grepl("ns",tab$'DD-SDvsCTRL')) |
                   (!grepl("ns",tab$'BB-SDvsCTRL') & grepl("ns",tab$'DD-SDvsCTRL'))
  
  tab$noSD <- grepl("ns",tab$'BB-SDvsCTRL') & grepl("ns",tab$'DD-SDvsCTRL')
  tab$up_both <- grepl("up",tab$'BB-SDvsCTRL') & grepl("up",tab$'DD-SDvsCTRL')
  tab$down_both <- grepl("down",tab$'BB-SDvsCTRL') & grepl("down",tab$'DD-SDvsCTRL')
  
  tab$BSL_trend <- ifelse(tab$baseline,"Significant in baseline", "No baseline effect")
  
  tab$SD_trend <- ifelse(tab$inverted,"Opposed effect",
                  ifelse(tab$protected,"Protected effect",
                  ifelse(tab$up_both,"SD up-regulated",
                  ifelse(tab$down_both ,"SD down-regulated",
                  ifelse(tab$noSD,"No effect",
                         "ERROR")))))
  return(tab)
}

process_QTL <- function(file, data, invertions, ATACdiff=NA) {
  QTL <- data[[grep(file,names(data))]]

  QTL <- QTL[rowSums(QTL[,grep('BB|DD',names(QTL))]<5)==0,] #remove imbalance
  

  if(!all(is.na(ATACdiff))) {
    QTL$region <- ATACdiff[QTL$lodcolumn,"region"]
    QTL$gene <- ATACdiff[QTL$lodcolumn,"Nearest.Gene"]
    QTL$Annotation <- ATACdiff[QTL$lodcolumn,"Annotation"]
    QTL$cCRE <- ATACdiff[QTL$lodcolumn,"cCRE"]
    
    QTL <- QTL[,c("lodcolumn", "region", "Annotation", "cCRE", "gene", "lod", "chr", "pos", "ci_lo", "ci_hi", "bp","bp_lo", "bp_hi", "test", grep('BB|DD',names(QTL), value= TRUE), "snp", "snps")]
    names(QTL)[1] <- "regionID"
  } else {
    QTL <- QTL[,c("lodcolumn", "lod", "chr", "pos", "ci_lo", "ci_hi", "bp","bp_lo", "bp_hi", "test", grep('BB|DD',names(QTL), value= TRUE), "snp", "snps")]
    names(QTL)[1] <- "Gene"
  }

  if(grepl('FC', file)) {
    inv <- gsub("QTL_table","qtl_ttests",file)
    inv <- gsub("FC","",inv)
    inv <- gsub(".csv",".RData",inv)
  
    inv <- classify_invertions(invertions[[grep(inv,names(invertions))]])
  
    if(!all(is.na(ATACdiff))) {
      names(inv)[1] <- "regionID"
      inv <- inv[,c("regionID", "chr", "pos", "snp", "BSL_trend", "SD_trend")]
    } else {
      names(inv)[1] <- "Gene"
      inv <- inv[,c("Gene", "chr", "pos", "snp", "BSL_trend", "SD_trend")]
    }
  
    QTL <- merge(QTL,inv, all.x=TRUE)
  }

  if(!all(is.na(ATACdiff))) {
    n <- c("regionID", "region",  "Annotation", "cCRE", "gene", "lod", "chr", "pos", "ci_lo", "ci_hi", "bp","bp_lo", "bp_hi", "BSL_trend", "SD_trend", grep('BB|DD',names(QTL), value= TRUE), "snp", "snps")
    if(!grepl('FC', file)) n <- n[!n %in% c("BSL_trend", "SD_trend")]
    QTL <- QTL[,n]
  } else {
    n <- c("Gene", "lod", "chr", "pos", "ci_lo", "ci_hi", "bp","bp_lo", "bp_hi", "BSL_trend", "SD_trend", grep('BB|DD',names(QTL), value= TRUE), "snp", "snps")
    if(!grepl('FC', file)) n <- n[!n %in% c("BSL_trend", "SD_trend")]
    QTL <- QTL[,n]
  }
  
  return(QTL)
}

QLT_overlap <- function(lst) {
  
  ov <- lst[[1]][,c(names(lst[[1]])[1],"snp")]
  
  for(i in 2:length(lst)) {
    n <- which(names(lst[[i]]) %in% c('Gene','region','Motif_TF', 'phenotype_ID'))[1]
    ov <- merge(ov, lst[[i]][,c(names(lst[[i]])[n],"snp")], all=TRUE)
    
  }

  
  return(ov)

}

cat("Process accessibility QTLs\n")
ATACdiff <- read.xlsx(opt$S2, "Differential_Accessibility")
rownames(ATACdiff) <- ATACdiff$Region_ID
ATACdiff$region <- paste0(ATACdiff$Chromosome,":",ATACdiff$Start,"-",ATACdiff$End)

aQTL <- process_QTL("QTL_table_atac.csv", qtls, types, ATACdiff)
aQTL_FC <-  process_QTL("QTL_table_atacFC.csv",qtls,types, ATACdiff)

addWorksheet(wb, "ATAC_QTL")
writeData(wb, "ATAC_QTL", aQTL, rowNames=FALSE)

addWorksheet(wb, "ATACxSD_QTL")
writeData(wb, "ATACxSD_QTL", aQTL_FC, rowNames=FALSE)

cat("Process expressions QTLs\n")
eQTL <- process_QTL("QTL_table_rna.csv", qtls, types)
eQTL_FC <- process_QTL("QTL_table_rnaFC.csv", qtls, types)

addWorksheet(wb, "RNA_QTL")
writeData(wb, "RNA_QTL", eQTL, rowNames=FALSE)

addWorksheet(wb, "RNAxSD_QTL")
writeData(wb, "RNAxSD_QTL", eQTL_FC, rowNames=FALSE)

cat("Process footprints QTLs\n")
fQTL <- process_QTL("QTL_table_footprints.csv", qtls, types)
fQTL_FC <- process_QTL("QTL_table_footprintsFC.csv", qtls, types)
names(fQTL)[1] <- 'Motif_TF'
names(fQTL_FC)[1] <- 'Motif_TF'

addWorksheet(wb, "TF_QTL")
writeData(wb, "TF_QTL", fQTL, rowNames=FALSE)

addWorksheet(wb, "TFxSD_QTL")
writeData(wb, "TFxSD_QTL", fQTL_FC, rowNames=FALSE)

cat("Process sleep QTLs\n")
pQTL <- qtls[[grep("QTL_table_sleep.csv",names(qtls))]]
phenos <- read.xlsx(paste0(opt$sleep,"/General_Information.xlsx"), sheet="Phenotypes")
phenos$PhenotypeID <- trimws(phenos$PhenotypeID)
translation <- read.xlsx(paste0(opt$sleep,"/phenotypes.xlsx"))

idx <- c(which(!is.na(phenos$Description)), nrow(phenos))

for(i in 1:(length(idx)-1)){
  #Fill in descriptions missing from excel format
  phenos$Description[idx[i]:(idx[i+1]-1)] <-  phenos$Description[idx[i]]
  phenos$Comments[idx[i]:(idx[i+1]-1)] <-  phenos$Comments[idx[i]]
  
}

pQTL$lodcolumn <- gsub("quant.out.","", pQTL$lodcolumn )
idx <- match(pQTL$lodcolumn,phenos$PhenotypeID)
pQTL$description <-  phenos[idx,]$Description 
pQTL$comments <-  phenos[idx,]$Comments 
pQTL$heritability<-  phenos[idx,]$Heritability 

idx <- match(pQTL$lodcolumn,translation$Phenotype)
pQTL$phenotype <-  translation[idx,]$translation 

pQTL <- pQTL[,c("lodcolumn","phenotype","description","comments", "heritability", "lod", "chr", "pos", "ci_lo", "ci_hi", "bp","bp_lo", "bp_hi", "test", "BB", "DD", "snp", "snps")]
names(pQTL)[1] <- "phenotype_ID"

addWorksheet(wb, "pheno_QTL")
writeData(wb, "pheno_QTL", pQTL, rowNames=FALSE)

cat("Overlap QTLs\n")
ov <- QLT_overlap(list('eQTL' = eQTL, 'eQTL_FC' = eQTL_FC,
                       'aQTL' = aQTL, 'aQTL_FC' = aQTL_FC,
                       'tfQTL' = fQTL, 'tfQTL_FC' = fQTL_FC, 
                       'pQTL' = pQTL))

ov <- ov %>% 
summarize(Motif_TF=paste(unique(sort(Motif_TF)), collapse='; '),
  Chromatin_region=paste(unique(sort(region)), collapse='; '),
  Gene=paste(unique(sort(Gene)), collapse='; '),
  phenotype=paste(unique(sort(phenotype_ID)), collapse='; '), .by =snp)

addWorksheet(wb, "QTL Peaks Overlaps")
writeData(wb, "QTL Peaks Overlaps", ov, rowNames=FALSE)

cat("Save table\n")
saveWorkbook(wb, paste0(opt$outdir,"/TableS6-QTL.xlsx"), overwrite = TRUE)

sessionInfo()