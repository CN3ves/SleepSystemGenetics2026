# Script to produce figure S11

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure11S.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("edgeR")
  library("tidyverse")
  library("openxlsx")
  library("ggplot2")
  library("ggpubr")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--atac"), type="character", default=NULL, 
              help="ATAC-seq counts", metavar="character"),
  make_option(c("-b", "--rna"), type="character", default=NULL, 
              help="RNA-seq counts", metavar="character"),
  make_option(c("-c", "--prints"), type="character", default=NULL, 
              help="BXD footprints", metavar="character"),
  make_option(c("-d", "--phenos"), type="character", default=NULL, 
              help="Sleep phenotypes", metavar="character"),
  make_option(c("-e", "--S2"), type="character", default=NULL, 
              help="Table S2", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$atac)){
  print_help(opt_parser)
  stop("ATAC-seq counts (-a) is missing", call.=FALSE)
}
if (is.null(opt$rna)){
  print_help(opt_parser)
  stop("RNA-seq counts(-b) is missing", call.=FALSE)
}
if (is.null(opt$prints)){
  print_help(opt_parser)
  stop("BXD footprints (-c) is missing", call.=FALSE)
}
if (is.null(opt$phenos)){
  print_help(opt_parser)
  stop("Sleep phenotypes (-d) is missing", call.=FALSE)
}
if (is.null(opt$S2)){
  print_help(opt_parser)
  stop("Table S2 (-e) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

# Helper functions
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

fix_names <- function(group) {
  group <- gsub("BXD0*", "BXD", group)
  group <- gsub("96","48a", group)
  group <- gsub("97","65a", group)
  group <- gsub("103","73b", group)

  group <- paste0('BXD', group)

  group <- gsub("BXDC57B[lL]6","C57Bl6", group)
  group <- gsub("BXDDBA","DBA", group)
  group <- gsub("DBA2","DBA", group)
  group <- gsub("BXDBXD","BXD", group)

  return(group)
}

cat("Load phenotypes \n")
sleep <- read.xlsx(opt$phenos)
rownames(sleep) <- fix_names(sleep$strID)
sleep <- sleep[,'delta-2.(2.25-3.25Hz)%', drop = FALSE]
colnames(sleep) <- 'Fast_delta2'

cat("Load footprint meta-analysis \n")
prints <- read.delim(opt$prints)
samples <- unique(gsub("Protection_Score_","",names(prints)[grep("Protection_Score_",names(prints))]))
for(s in samples) {
  prints[,s] <- rowSums(prints[,grepl(s, names(prints))])
}
tfs <- prints[,-c(grep("Protection_Score_", names(prints)), grep("TC_", names(prints)))]
tfs$Motif <- toupper(gsub("MA.*\\.","",gsub("\\(.*","",tfs$Motif)))
tfs <- tfs %>% select(!Num) %>% group_by(Motif) %>% summarise_all(mean) %>% as.data.frame
rownames(tfs) <- tfs$Motif
tfs <- tfs[,-1]
colnames(tfs) <- fix_names(colnames(tfs))

cat("Load mean values \n")
atac_fc <- collapse_means(readRDS(opt$atac))
rna_fc <- collapse_means(readRDS(opt$rna))

samples <-  colnames(atac_fc)[ colnames(atac_fc) %in% colnames(rna_fc)]

rna <- as.data.frame(rna_fc[rownames(rna_fc) %in% c('Nrf1', 'Lgr6','Poln'),samples])
tfs <- as.data.frame(tfs[rownames(tfs) %in% c('FOSL2', 'JUNB','NRF1'),samples])

ATACdiff <- read.xlsx(opt$S2, sheet="Differential_Accessibility" )
ATACdiff <- ATACdiff[ATACdiff$Nearest.Gene %in% 'Nrf1',]
region <- ATACdiff[as.logical(ATACdiff$Differential.and.Correlated.in.all.samples),'Region_ID']

atac <- as.data.frame(atac_fc[rownames(atac_fc) %in% region,samples, drop=FALSE]) 

cat("Merge data \n")
df <- rbind(rna,atac,tfs) %>% t %>%  as.data.frame %>% rownames_to_column('Treatment') %>%  
  mutate(strain = gsub('_.*','',Treatment) ,Treatment = gsub('.*_','',Treatment)) 

sleep <- as.data.frame(sleep[unique(df$strain),,drop=FALSE]) %>%  as.data.frame %>% 
  rownames_to_column('strain') 

df <- merge(df, sleep, by='strain') %>% 
  pivot_longer(!c(strain, Treatment,Fast_delta2),names_to='marker', values_to = 'value')

svg(paste0(opt$outdir, "/FigS11.svg"),width=20,height=20)

ggplot(df, aes(x=value,y=Fast_delta2, shape=Treatment, col = Treatment)) +
  geom_point(size=4) +
  theme_classic() + theme(text = element_text(size = 20)) + 
  xlab("Average marker intensity") + ylab("Average fast delta power after SD") +
  geom_smooth(method = 'lm', se=FALSE) +
  stat_cor(size = 8, show.legend=FALSE) + 
  facet_wrap(marker ~ Treatment, scales="free")

dev.off()

write.csv(as_data_frame(df, what = "edges"), paste0(opt$outdir, "/data/S11.csv"))

sessionInfo()