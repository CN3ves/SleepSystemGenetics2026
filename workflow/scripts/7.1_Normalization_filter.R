# Script designed to filter counts

# Redirect all R logs to Snakemake log
log <- file('logs/7-BXD_Normalzation/filter_7.1.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("edgeR")
  library("ggplot2")
  library("tidyr")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="RDS counts object", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("RDS counts object (-c) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
counts <- readRDS(opt$counts)
exp <- gsub(".*/(.*?)_original.*", "\\1", opt$counts)

# Remove genomic scaffolds
print(paste(exp, "count dimensions (original):", paste(dim(counts), collapse = ", ")))
counts$counts <- counts$counts[grep("_.*_",rownames(counts$counts), invert = TRUE),]
  
# Filter counts 
print(paste(exp, "count dimensions (only standard chrs):", paste(dim(counts), collapse = ", ")))
keep <- filterByExpr(counts, group=counts$samples$Treatment)

counts_filter <- counts[keep,,keep.lib.sizes=FALSE]
print(paste(exp, "count dimensions (filtered by counts):", paste(dim(counts_filter), collapse = ", "), "| Percentage kept:", round(mean(keep)*100,2), "%"))

# Plot mean reads per feature distribution
mean_counts <- apply(counts,1,mean)
mean_counts <- density(log2(mean_counts+1))
mean_counts_filter <- apply(counts_filter,1,mean)
mean_counts_filter <- density(log2(mean_counts_filter+1))

xlab <- paste0("Unfiltered: N=", mean_counts_filter$n, "; Bandwidth=", signif(mean_counts_filter$bw,3), "\nFiltered: N=", mean_counts$n, "; Bandwidth=", signif(mean_counts$bw,3))
  
png(paste0(opt$outdir, "/QC_", exp, "_count_dist.png"))
plot(mean_counts_filter, col="green", main="Distribution of log2 mean reads per ORC", xlab = xlab)
lines(mean_counts, col="red")
legend("topright", legend=c("pre-filter", "post-filter"), pch=16, col=c("red", "green"))
dev.off()

# check samples with a lot of 0s counts
print(paste(exp, "samples with highest percentages of 0 counts:"))
print(round(sort(apply(counts$counts,2,function(x) mean(x==0)), decreasing = T)[1:15]*100,2))

saveRDS(counts_filter, paste0(opt$outdir, "/", exp, "_counts.RData"))

# Check counts per million
cpms <- cpm(counts_filter)
print(paste(exp, "proportion of values with very low cpm", round(mean(cpms<1), 3)))
print(paste(exp, "proportion of values with low cpm", round(mean(cpms>1 & cpms < 5), 3)))
print(paste(exp, "proportion of values with high cpm", round(mean(cpms>10), 3)))

cpms <- rbind(cpms, counts_filter$samples[colnames(cpms), "Treatment"])
rownames(cpms)[nrow(cpms)] <- "Treatment"
cpms <- rbind(cpms, colnames(cpms))
rownames(cpms)[nrow(cpms)] <- "Sample"

cpms <- cpms %>% 
  t() %>% 
  data.frame(stringsAsFactors = FALSE) %>%
  pivot_longer(1:(nrow(cpms)-2), names_to = "perc", values_to ="val") 

cpms <- as.data.frame(cpms, stringsAsFactors = FALSE)
cpms$val <- as.numeric(cpms$val)
  
png(paste0(opt$outdir, "/QC_", exp, "_count_cpm.png"))
ggplot(cpms, aes(x = Sample, y=val, col=Treatment)) + 
  geom_boxplot()  +
  labs(title="Count per million reads") + xlab("Samples") + ylab("CPM")+
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png(paste0(opt$outdir, "/QC_", exp, "_count_cpm2.png"))
ggplot(cpms, aes(x = Sample, y=val, col=Treatment)) + 
  geom_boxplot()  +
  labs(title="Count per million reads") + xlab("Samples") + ylab("CPM")+
  theme_classic() + ylim(0,10) + 
  geom_hline(yintercept=1, color="red") + 
  geom_hline(yintercept=5, color="yellow") + 
  theme(plot.title = element_text(hjust = 0.5)) 
dev.off()

sessionInfo()