# Script to produce figure S4b

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure4Sb.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("qtl2")
  library("edgeR")
  library("ggplot2")
  library("openxlsx")
  library("ggpubr")
  library("tidyverse")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--qtl"), type="character", default=NULL, 
              help="Diretory for genetic information", metavar="character"),
  make_option(c("-b", "--atac"), type="character", default=NULL, 
              help="ATAC counts", metavar="character"),
  make_option(c("-c", "--rna"), type="character", default=NULL, 
              help="RNA counts", metavar="character"),
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
if (is.null(opt$atac)){
  print_help(opt_parser)
  stop("ATAC counts (-b) is missing", call.=FALSE)
}
if (is.null(opt$rna)){
  print_help(opt_parser)
  stop("RNA counts (-c) is missing", call.=FALSE)
}
if (is.null(opt$S6)){
  print_help(opt_parser)
  stop("Table S6 (-d) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load Table S6\n")
atac_QTL <- read.xlsx(opt$S6, "ATACxSD_QTL")
rna_QTL <-  read.xlsx(opt$S6, "RNAxSD_QTL")

collapse_counts <- function(obj) {
  
  mean_counts <- lapply(unique(obj$samples$group), function(group) {
    samples <- obj$samples$group == group
    
    group <- gsub("BXD0*", "BXD", group)
    group <- gsub("96","48a", group)
    group <- gsub("97","65a", group)
    group <- gsub("103","73b", group)
    
    samples_means <- data.frame(rowMeans(cpm(obj)[,samples, drop = FALSE]))
    names(samples_means) <- group
    
    return(samples_means)
  })
  
  as.data.frame(mean_counts)
}

cat("Load counts\n")
atac_cpm <- collapse_counts(readRDS(opt$atac))
rna_cpm <- collapse_counts(readRDS(opt$rna))

rna_cpm <- rna_cpm[,colnames(atac_cpm)]
atac_cpm <- atac_cpm[,colnames(rna_cpm)]

rna_data <- read_cross2(file =  paste0(opt$qtl,"/rna.json"))

select <- list("snp" = c("rs46476185","rs46476185","rs47499225","rs4165503"),
               "region" = c("chr11:113925739-113925989", "chr6:40250165-40250406","chr7:29371318-29371568","chr15:98804780-98805030"),
               "gene" = c("Zfp383", "Zfp383","Tceb3","Dph2"))


points <- list()
for (i in 1:length(select[[1]])) {
  region <- select[['region']][i]
  gene <- select[['gene']][i]
  snp <- select[['snp']][i]

  a_qtl <- atac_QTL[atac_QTL$region %in% region & atac_QTL$snp %in% snp,]
  r_qtl <- rna_QTL[rna_QTL$Gene == gene & rna_QTL$snp %in% snp,]

  c <- cor.test(as.numeric(atac_cpm[a_qtl$regionID,]),as.numeric(rna_cpm[gene,]))
  print(c)
  
  chr <- gsub(':.*','',region)
  atac_data <- read_cross2(file =  paste0(opt$qtl,"/atac",chr,".json"))

  atac_map <- insert_pseudomarkers(map=atac_data$gmap, step=1)
  rna_map <- insert_pseudomarkers(map=rna_data$gmap, step=1)
    
  atac_pr <- calc_genoprob(cross=atac_data, map=atac_map, error_prob=0.0001)
  rna_pr <- calc_genoprob(cross=rna_data, map=rna_map, error_prob=0.0001)

  atac_geno <-  maxmarg(atac_pr, atac_map, chr=a_qtl$chr, pos=a_qtl$pos, return_char=TRUE)
  rna_geno <-  maxmarg(rna_pr, rna_map, chr=r_qtl$chr, pos=r_qtl$pos, return_char=TRUE)[names(atac_geno)]
    
  atac_vals <- as.numeric(atac_cpm[a_qtl$regionID,names(atac_geno)])
  rna_vals <- as.numeric(rna_cpm[r_qtl$Gene,names(rna_geno)])
  
  # plot 
  d <- data.frame(geno=c(atac_geno,rna_geno), val = c(atac_vals,rna_vals), 
                  dataset=c(rep(a_qtl$region,length(atac_geno)),rep(r_qtl$Gene,length(rna_geno))),
                  treat = gsub("\\..*","", gsub(".*_","",c(names(atac_geno),names(rna_geno)))))
  d$group <- paste(d$geno, d$treat)
  d <- d[!is.na(d$geno),]
  d$snp <- snp

  points[[snp]] <- d

  p <- ggplot(d, aes(x=group, y=val, fill=treat)) + 
    geom_boxplot(position=position_dodge(1)) +
    geom_point(aes(fill = treat), shape = 21, position = position_jitterdodge()) + 
    labs(title=unique(a_qtl$snp,r_qtl$snp),x="Genotype", y = "Average counts per million") +
    theme_classic() + stat_compare_means(comparisons = list(c("DD SD", "DD CTRL"), c("BB SD", "BB CTRL"),   c("DD SD", "BB SD"), c("DD CTRL", "BB CTRL")), method = "t.test") +
    theme(plot.title=element_text(hjust=0.5),text=element_text(size=16)) + facet_wrap(~ dataset,scales="free")
    
   svg(paste0(opt$outdir,"/FigS4b_",snp,".svg"), width=14,height=10)
   print(p)
   dev.off()  

}

write.csv(do.call(rbind,points), paste0(opt$outdir, "/data/S4b.csv"), row.names=FALSE)

sessionInfo()
