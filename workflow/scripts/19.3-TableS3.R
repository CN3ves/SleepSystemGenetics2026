# Script to produce table S3

# Redirect all R logs to Snakemake log
log <- file('logs/19-Tables/S3_19.3.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("openxlsx")
  library("org.Mm.eg.db")
  library("ChIPseeker")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
  library("clusterProfiler")
  library("tidyverse")
  library("optparse")
})

cat("Checking arguments\n")

option_list = list(
  make_option(c("-a", "--S2"), type="character", default=NULL, 
              help="Table S2", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$S2)){
  print_help(opt_parser)
  stop("Table S2 (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

# Create a new workbook and add a sheet
wb <- createWorkbook()

cat("Load differential analyses\n")
# Get differential results
atac_stats <- read.xlsx(opt$S2, "Differential_Accessibility")
rna_stats <- read.xlsx(opt$S2, "Differential_Expression")

cat("Prepare RNA results \n")
rna_genes <- rna_stats$RNA_logFC[order(rna_stats$RNA_logFC,decreasing=T)]
names(rna_genes) <- rna_stats$Gene_ID[order(rna_stats$RNA_logFC,decreasing=T)]
  
rna_convert <-  AnnotationDbi::select(org.Mm.eg.db, 
       keys = names(rna_genes),
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

names(rna_genes) <-  rna_convert$ENTREZID
rna_genes <- rna_genes[!is.na(names(rna_genes))]

cat("Prepare ATAC results \n")
stats <- atac_stats %>% filter(Nearest.Gene %in% rna_stats$Gene_ID) %>% group_by(Nearest.Gene) %>% summarise(FC=mean(log2(mean(2**ATAC_logFC))))
  
atac_genes <- stats$FC[order(stats$FC,decreasing=T)]
names(atac_genes) <- stats$Nearest.Gene[order(stats$FC,decreasing=T)]

atac_convert <-   AnnotationDbi::select(org.Mm.eg.db, 
       keys = names(atac_genes),
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")
  
names(atac_genes) <- atac_convert$ENTREZID
atac_genes <- atac_genes[!is.na(names(atac_genes))]

cat("Run GO enrichment\n")
convert2symbol <- function(enrichment,covert_tab) {
  covert_tab <- covert_tab[!is.na(covert_tab$ENTREZID),]
  rownames(covert_tab) <- covert_tab$ENTREZID
  
  symbols <- lapply(strsplit(enrichment@result$core_enrichment, "/"), function(row) {
    paste(sort(covert_tab[row,"SYMBOL"]),collapse = "/")
  })
  
  return(unlist(symbols))
}

rna_GO <- gseGO(rna_genes, ont = "ALL", OrgDb="org.Mm.eg.db")
rna_GO@result$core_enrichment_SYMBOL <- convert2symbol(rna_GO, rna_convert)

atac_GO <- gseGO(atac_genes, ont = "ALL", OrgDb="org.Mm.eg.db")
atac_GO@result$core_enrichment_SYMBOL <- convert2symbol(atac_GO, atac_convert)

cat("Run KEGG enrichment\n")
rna_KEGG <- gseKEGG(rna_genes, organism = "mmu")
rna_KEGG@result$core_enrichment_SYMBOL <- convert2symbol(rna_KEGG, rna_convert)
rna_KEGG@result$Description <- gsub(" - Mus musculus.*","", rna_KEGG@result$Description)

atac_KEGG <- gseKEGG(atac_genes, organism = "mmu")
atac_KEGG@result$core_enrichment_SYMBOL <- convert2symbol(atac_KEGG, atac_convert)
atac_KEGG@result$Description <- gsub(" - Mus musculus.*","", atac_KEGG@result$Description)

cat("Get overlapping enrichment\n")
overlap_enrich <- function(atac, rna) {
  atac_df <- atac@result[atac@result$Description %in% rna@result$Description,]
  rna_df  <- rna@result[rna@result$Description %in% atac@result$Description,]
  
  if(!"ONTOLOGY" %in% names(rna_df)) {
    rna_df$ONTOLOGY <- atac_df$ONTOLOGY <- "KEGG"
  } else {
    rna_df$ONTOLOGY <- paste("GO", rna_df$ONTOLOGY)
    atac_df$ONTOLOGY <- paste("GO", atac_df$ONTOLOGY)
  }
  
  rna_df <- rna_df[rownames(atac_df),c("ID","Description","ONTOLOGY","core_enrichment_SYMBOL")]
  
  rna_genes <- strsplit(rna_df$core_enrichment_SYMBOL, "/")
  atac_genes <- strsplit(atac_df$core_enrichment_SYMBOL, "/")
  
  for(i in 1:length(rna_genes)) {
    rna_df$core_enrichment_SYMBOL[[i]] <-paste(sort(
      rna_genes[[i]][rna_genes[[i]] %in% atac_genes[[i]]]),
      collapse=", ")
    
  }
  
  names(rna_df) <- c("ID","Description","Type","Genes")
  return(rna_df)
}

ov_kegg <- overlap_enrich(atac_KEGG, rna_KEGG)
ov_go <- overlap_enrich(atac_GO, rna_GO)

ov <- rbind(ov_go,ov_kegg)
 
#save for plotting
enrichment <- list("RNA_GO"=rna_GO,"ATAC_GO"=atac_GO,
                   "RNA_KEGG"=rna_KEGG,"ATAC_KEGG"=atac_KEGG, 
                   "Overlap" = ov)

saveRDS(enrichment,paste0(opt$outdir,"/data/enrichment.RData"))

for(enrich in names(enrichment)) {
 
  addWorksheet(wb, enrich)
  writeData(wb, enrich,enrichment[[enrich]], rowNames=FALSE)
}

cat("Save table\n")
saveWorkbook(wb, paste0(opt$outdir,"/TableS3-Enrichment_BXD.xlsx"), overwrite = TRUE)
        
sessionInfo()