# Script designed to plot read QC summaries

# Redirect all R logs to Snakemake log
log <- file('logs/1-BXD_fastq/QCsummary_2.2.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("dplyr")
  library("tidyverse")
  library("ggplot2")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Read's QC file", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="Processed metadata file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Read's QC file (-d) is missing", call.=FALSE)
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
summary_tab <- read.csv(opt$file, row.names=1)
metadata <- read.csv(opt$meta, row.names=1)
metadata$Sample <- sprintf("%03d",as.numeric(metadata$Sample))

# Set parameters for each plot
rows <- c("q.*rate", "content", ".*uplica.*", ".*_reads.*")
labels <- list(c("Read QC summary", "Minimum quality", "Proportion"),
               c("Content summary", "Base", "Percentage"),
               c("Duplication summary", "Duplication", "Proportion"),
               c("Reads summary", "Reads", "M reads"))

cat("Generating boxplot diagnostics\n")
# Box plots
for (i in 1:length(rows)) {
  print(paste("Plot", labels[[i]][1]))
  
  # Subset table for plot
  row <- grep(rows[i], rownames(summary_tab), value=TRUE)
  plot <- summary_tab[row,-1]
  plot["treatment",] <- metadata$Treatment[match(metadata$Sample, gsub("X([0-9]*)_.*","\\1",colnames(plot)))]
  
  # Convert table from wide to long
  plot <- plot %>% 
    t() %>%
    data.frame(stringsAsFactors = FALSE) %>%
    pivot_longer(1:length(row), names_to = "perc", values_to ="val")
  plot <- plot[,c("perc", "val", "treatment")]
  
  # Make colour groups
  plot$group <- paste("Raw", plot$treatment)
  plot$group[grep("after", plot$perc)] <- paste("Filtered", plot$treatment[grep("after", plot$perc)])

  # Correct x categories
  plot$perc <- gsub("q(.*)_rate.*.filtering.", "Bases with quality > \\1", plot$perc) # i=1
  plot$perc <- gsub("content.*", "%", plot$perc) # i=2
  plot$perc <- gsub("\\.with\\.3\\.\\.dupl", ".with.3+.dupl", plot$perc) # i=3
  plot$perc <-gsub("\\.{3}(.*).$", " (%\\1)", plot$perc) # i=3
  plot$perc <- gsub("total_reads.*", "total_reads", plot$perc) # i=4
  plot$perc <- gsub("\\.", " ", plot$perc)
  
  # Correct value types
  plot$perc <- as.factor(plot$perc)
  plot$val <- as.numeric(plot$val)
  plot$group <- factor(plot$group, levels = c("Raw CTRL", "Filtered CTRL", "Raw SD", "Filtered SD"))
  plot <- plot[! is.na(plot$val),]

  # Save plot
  p <- ggplot(plot, aes(x=perc, y=val, fill=group)) + 
    geom_boxplot() + 
    labs(title=labels[[i]][1]) + xlab(labels[[i]][2]) + ylab(labels[[i]][3]) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
          plot.title = element_text(hjust = 0.5))
  
  png(paste0(opt$outdir,"/QC_plot", gsub(" ", "",labels[[i]][1]),".png"))
  print(p)
  dev.off()
}

# Histogram of read lengths
col <- "read1_mean_length (after filtering)"
hist_data <- as.data.frame(t(summary_tab[col,-1]))

hist_data$reads <- as.numeric(hist_data[,col])
hist_data$Treatment <- metadata[match(gsub("X([0-9]*)_.*","\\1",rownames(hist_data)),metadata$Sample),"Treatment"]

cat("Generating histogram diagnostics\n")
png(paste0(opt$outdir,"/QC_plot_density.png"))
ggplot(hist_data, aes(x = reads, col=Treatment, fill=Treatment)) + 
  geom_density(aes(alpha = .25)) +
  labs(title="Distribution of read lengths") + xlab("Mean read length") + ylab("Density")+
  scale_alpha(guide = "none") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

sessionInfo()