# Script designed to aggregate all alignment QC summaries

# Redirect all R logs to Snakemake log
log <- file('logs/3-BXD_bam/alignQC_3.2.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ggplot2")
  library("tidyr")
  library("optparse")
})
source("workflow/scripts/3.0-QC_align_helper.R")

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="Directory with the bowtie log output", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="Processed metadata file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("Directory containing the bowtie logs (-d) is missing", call.=FALSE)
}
if (is.null(opt$meta)){
  print_help(opt_parser)
  stop("Processed metadata file (-m) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
logs <- list.files(opt$dir, pattern='align_', full.names = TRUE)
metadata <- read.csv(opt$meta, row.names=1)
metadata$Sample <- sprintf("%03d",metadata$Sample)

# Extract alignment results
align_summary <- alignQC(logs)

# Samples with less than 200M total reads
san_tab <- align_summary[-grep('sub',align_summary$File),c("File","Reads")]
san_tab$File <- gsub("_.*","",san_tab$File)
san_tab <- merge(metadata, san_tab, by.x="Sample", by.y="File")
san_tab$Reads <- as.numeric(san_tab$Reads)
san_tab <- aggregate(san_tab$Reads, by=list(Category=san_tab$Group), FUN=sum)
san_tab <- san_tab[order(san_tab$x, decreasing =TRUE),]
print("GROUPS with low total reads (<200M) - footprinting")
print(paste0(sum(san_tab$x < 200000000),'/',length(san_tab$x)))
print("GROUPS with >200M total reads")
print(san_tab[san_tab$x >= 200000000,])
print("GROUPS with the lowest total reads")
print(tail(san_tab[san_tab$x < 200000000,]))

cat("Generating diagnostic plots\n")
# Plot Summaries
plot_df <- align_summary[-1,]
stopifnot(all(metadata[match(gsub("_.*", "",align_summary$File[-1]),metadata$Sample),"Sample"] == gsub("_.*", "",align_summary$File[-1])))
plot_df$Treatment <- metadata[match(gsub("_.*", "",align_summary$File[-1]),metadata$Sample),"Treatment"] 

plot_df$Reads <- as.numeric(plot_df$Reads)
plot_df$sampling <- gsub('.*_','',plot_df$File)
  
png(paste0(opt$outdir,"/QC_align_reads.png"))
ggplot(plot_df, mapping = aes(x=Treatment, y = Reads, fill=Treatment)) +
  geom_boxplot() + 
  labs(title="Reads considered for alignment") + ylab("Number of reads") +
  theme_classic() + 
  theme( plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  facet_wrap(~ sampling, scales ='free_y')
dev.off()

plot_df <- plot_df[,-1:-2]

plot_df <- plot_df %>% 
  data.frame(stringsAsFactors = FALSE) %>%
  pivot_longer(1:(ncol(plot_df)-2), names_to = "perc", values_to ="val")
plot_df$val <- as.numeric(gsub("%","",plot_df$val))/100

png(paste0(opt$outdir,"/QC_align_mapping.png"))
  ggplot(plot_df , mapping = aes(x= perc, y = val, fill = Treatment)) +
  geom_boxplot() + 
  labs(title="Alignment QC") + ylab("Proportion") +
  theme_classic() + 
  theme( plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  facet_wrap(~ sampling, scales ='free_y')
dev.off()

cat("Generating summary table\n")
# Write alignment QC table
write.csv(align_summary, paste0(opt$outdir,"/QC_align.csv"))

sessionInfo()