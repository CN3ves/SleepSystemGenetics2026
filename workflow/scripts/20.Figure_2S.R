# Script to produce figure S2

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure2S.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ggplot2")
  library("edgeR")
  library("openxlsx")
  library("tidyverse")
  library("ggpubr")
  library("GenomicRanges")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--atac"), type="character", default=NULL, 
              help="ATAC counts", metavar="character"),
  make_option(c("-b", "--rna"), type="character", default=NULL, 
              help="RNA counts", metavar="character"),
  make_option(c("-c", "--S2"), type="character", default=NULL, 
              help="Table S2", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$atac)){
  print_help(opt_parser)
  stop("ATAC counts (-a) is missing", call.=FALSE)
}
if (is.null(opt$rna)){
  print_help(opt_parser)
  stop("RNA counts (-b) is missing", call.=FALSE)
}
if (is.null(opt$S2)){
  print_help(opt_parser)
  stop("Table S2 (-c) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load counts\n")
collapse_means <- function(obj) {
  
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

atac_cpm <- collapse_means(readRDS(opt$atac))
rna_cpm <- collapse_means(readRDS(opt$rna))

rna_cpm <- as.data.frame(rna_cpm[,colnames(atac_cpm)])
atac_cpm <- as.data.frame(atac_cpm[,colnames(rna_cpm)])

cat("Load Table S2\n")
atac <- read.xlsx(opt$S2, sheet="Differential_Accessibility")
rna <- read.xlsx(opt$S2, sheet="Differential_Expression")

cors <- atac[atac$ATAC_FDR <0.05 & atac$Pearson_correlation_all != "-",]
cors$Pearson_correlation_all <- as.numeric(cors$Pearson_correlation_all)

cat("Get examples\n")
#genes <- c("Nrf1","Sgsm1", "Fxyd3","Nav1", "Ncald", "Homer1",'Rbfox3','Tshz1','Dcaf7','Klf4','Ano2','Rasl11b','Cadm1')
genes <- c("Hps1", "Nrf1","Tshz1", "Ncald", 'Nav1','Sgsm1')

all(genes %in% rna$Gene_ID[rna$RNA_FDR < 0.05])

region <- cors[cors$Nearest.Gene %in% genes,] %>% group_by(Nearest.Gene) %>% filter(abs(Pearson_correlation_all) == max(abs(Pearson_correlation_all))) %>% pull(Region_ID) #pick best correlating region
region[grep("chr15", region)] <-  cors[1,"Region_ID"] # switch Ncald's region for top significant

cat("Filter dataset\n")
atac_sig <- atac_cpm[rownames(atac_cpm) %in% region,]
rna_sig <- rna_cpm[rownames(rna_cpm) %in% genes,]
 
atac_sig <- atac_sig %>% mutate(id=rownames(atac_sig)) %>% pivot_longer(!id, names_to="sample",values_to="counts") 

sub_cors <- cors[cors$Region_ID %in% atac_sig$id,]
atac_sig$gene <- sub_cors[match(atac_sig$id, sub_cors$Region_ID),]$Nearest.Gene

rna_sig <- rna_sig %>% mutate(id=rownames(rna_sig)) %>% pivot_longer(!id, names_to="sample",values_to="counts") 

df <- merge(atac_sig,rna_sig,by.x=c("gene", "sample"), by.y=c("id","sample"))

df$region <- paste(GRanges(atac[match(df$id, atac$Region_ID),]))

df$label <- paste(df$region, df$gene,sep=" -> ")

df <- df[,c("label", "sample",  "counts.x",  "counts.y")]
names(df) <- c("label", "sample",  "ATAC", "RNA")

df$treatment <- gsub(".*_","", df$sample)
df$strain <- gsub("_.*","", df$sample)
df$ATAC <- as.numeric(df$ATAC )
df$RNA <- as.numeric(df$RNA )

df$label <- factor(df$label, level=c("chr1:135623342-135623558 -> Nav1",
                                     "chr19:42788388-42788638 -> Hps1", 
                                     "chr6:30028536-30028786 -> Nrf1",
                                     "chr5:113297968-113298218 -> Sgsm1", 
                                     "chr18:83895324-83895574 -> Tshz1",
                                     "chr15:37607642-37607892 -> Ncald"))

g <- ggplot(df, aes(x=ATAC,y=RNA, shape=treatment, col = treatment)) +
  geom_point(size=2) +
  geom_line(aes(group=strain), col="grey", linetype="dashed") + 
  theme_classic() + 
  ylab("Average RNA-seq cpm") + xlab("Average ATAC-seq cpm") +
  geom_smooth(aes(group=1), method = 'lm', col="black") +
  stat_cor(aes(group=1), size = 4, show.legend=FALSE) + 
  facet_wrap(~ label, scales="free")


cat("Save plot\n")
svg(paste0(opt$outdir,"/FigS2.svg"), width=17, height= 10)
print(g)
dev.off()

write.csv(df, paste0(opt$outdir, "/data/S2.csv"))

sessionInfo()