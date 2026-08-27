# Script designed to aggregate all alignment QC summaries

# Redirect all R logs to Snakemake log
log <- file('logs/16-NRF_bam/alignQC_16.2.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ggplot2")
  library("tidyr")
  library("optparse")
})
source("workflow/scripts/16.0-QC_align_helper.R")

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
logs <- list.files(opt$dir, pattern='Log.final', full.names = TRUE)
metadata <- read.csv(opt$meta, row.names=1)
rownames(metadata) <- metadata$SampleName

# Extract alignment results
align_summary <- alignQC(logs)

cat("Generating diagnostic plots\n")
# Plot Summaries
plot_df <- align_summary[-1,]
stopifnot(all(rownames(metadata)[match(align_summary$File[-1],rownames(metadata))] == align_summary$File[-1]))
plot_df$Treatment <- metadata[match(align_summary$File[-1],rownames(metadata)),"Group"] 

plot_df$Reads <- as.numeric(plot_df$Reads)
  
png(paste0(opt$outdir,"/QC_align_reads.png"))
ggplot(plot_df, mapping = aes(x=Treatment, y = Reads, fill=Treatment)) +
  geom_boxplot() + 
  labs(title="Reads considered for alignment") + ylab("Number of reads") +
  theme_classic() + 
  geom_hline(yintercept=50000, color = "grey") + 
  theme( plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()

plot_df <- plot_df[,c(grep("rate", names(plot_df),value=TRUE), "Treatment")]

plot_df <- plot_df %>% 
  data.frame(stringsAsFactors = FALSE) %>%
  pivot_longer(1:(ncol(plot_df)-1), names_to = "perc", values_to ="val")
plot_df$val <- as.numeric(gsub("%","",plot_df$val))/100

png(paste0(opt$outdir,"/QC_align_mapping.png"))
ggplot(plot_df , mapping = aes(x= perc, y = val, fill = Treatment)) +
  geom_boxplot() + 
  labs(title="Alignment QC") + ylab("Proportion") +
  theme_classic() + 
  theme( plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()

cat("Generating summary table\n")
# Write alignment QC table
write.csv(align_summary, paste0(opt$outdir,"/QC_align.csv"))

sessionInfo()