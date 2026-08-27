# Script designed to aggregate all peaks and counts QC summaries

# Redirect all R logs to Snakemake log
log <- file('logs/6-BXD_features/QC_6.3.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("edgeR")
  library("dplyr")
  library("tidyr")
  library("ggplot2")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-l", "--logs"), type="character", default=NULL, 
              help="Directory with the HTseq logs", metavar="character"),
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="Directory with the count files", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="Processed metadata file", metavar="character"),
  make_option(c("-f", "--features"), type="character", default=NULL, 
              help="Directory with the bed/gtf feature files", metavar="character"),
  make_option(c("-r", "--readsinpeaks"), type="character", default=NULL, 
              help="Feature counts Frip summary files", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$logs)){
  print_help(opt_parser)
  stop("Directory with the HTseq logs (-l) is missing", call.=FALSE)
}
if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("Directory with the count files (-c) is missing", call.=FALSE)
}
if (is.null(opt$meta)){
  print_help(opt_parser)
  stop("Processed metadata file (-m) is missing", call.=FALSE)
}
if (is.null(opt$features)){
  print_help(opt_parser)
  stop("Directory with the bed/gtf feature files (-f)", call.=FALSE)
}
if (is.null(opt$readsinpeaks)){
  print_help(opt_parser)
  stop("Feature counts Frip summary files (-r)", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
logs <- list.files(opt$logs, pattern = "peaks", full.names = TRUE)
counts <- strsplit(opt$counts, " ")[[1]]
feature <- opt$features
metadata <- read.csv(opt$meta, row.names=1)
metadata$Sample <- sprintf("%03d",metadata$Sample)
frips <- strsplit(opt$readsinpeaks, " ")[[1]]

# Collect QC data summaries
qc_table <- list()
for (log in logs) {
  qc <- read.delim(file = log, header = FALSE)

  # Make table
  qc <- qc[-c(grep("Warning",qc$V1),grep("Logs",qc$V1)),1]
  qc <- sapply(as.data.frame(do.call(rbind, strsplit(qc, ":"))), trimws)

  # Clean names
  rm_idx <- c()
  for (idx in grep("^\\(", qc[,1])) {
    qc[idx-1,1] <- paste(qc[idx-1,1], qc[idx,1])
    rm_idx <- c(rm_idx, -idx)
  }
  qc <- qc[rm_idx,]
  
  for (idx in grep("^duplicates$", qc[,1])) {
    qc[idx,1] <- paste(qc[idx-1,1], qc[idx,1])
  }
 
  rownames(qc) <- trimws(qc[,1])
  
  # Split individual samples
  idxs <- grep("Processing experimental file", qc[,1])
  samples <- c(idxs, nrow(qc))
  
  for (i in 1:length(idxs)) {
    idx <- samples[i:(i+1)]
    sample_qc <- qc[idx[1]:idx[2],]
    
    # Get sample name
    sample <- grep("Processing experimental file #", sample_qc[,1])[1]
    sample <- gsub("\\..*","", basename(sample_qc[sample,2]))
    
    # Remove uninformative rows
    rm <- grep("Processing experimental file", sample_qc[,1])
    for (idx in 1:nrow(sample_qc)){
      if (sample_qc[idx,1]==sample_qc[idx,2]) {
        rm <- c(rm,idx)
      }
    }

    qc_table[[sample]] <- sample_qc[-rm,-1]
  }
}

# Match names
n <- unique(unlist(lapply(qc_table, names)))
qc_table <- as.data.frame(sapply(qc_table,function(x) {x<-x[n]; names(x)<-n; x}))
# Remove rows with no variations
print(paste("Genome size:", unique(qc_table["Genome length",][!is.na(qc_table["Genome length",])])))
qc_table <- qc_table[apply(qc_table, 1, function(r) length(unique(r[!is.na(r)])) != 1),]

# Estimate FRiP
counts <- readDGE(counts)

reads_in_peak <- colSums(counts$counts[-grep("__", rownames(counts$counts)),])
all_reads <- colSums(counts$counts)
multi_count <- counts$counts["__ambiguous",]

frip <- (reads_in_peak - multi_count) / (all_reads  - 2*multi_count)

print(paste("Frip > 0.3:", mean(frip > 0.3)))
print(paste("Frip > 0.2:", mean(frip > 0.2)))

# load frips from featurecount
frips <-  do.call(rbind, lapply(frips, read.csv, row.names=1))
rownames(frips) <- gsub('list','Final_features',rownames(frips))
rownames(frips) <- gsub('mito','Mitochondrial_features',rownames(frips))
rownames(frips) <- gsub('parent','Parent_line_subset',rownames(frips))
rownames(frips) <- gsub('subsample','Downsampled',rownames(frips))
rownames(frips) <-paste0('FeatureCounts_FRiP:',rownames(frips))
names(frips) <- gsub('^X','',names(frips))

# Save QC table
n <- c(rownames(qc_table), "FRiP", rownames(frips))
qc_table <-  rbind(qc_table, frip, frips[, names(qc_table)])
rownames(qc_table) <- n

plot <- qc_table
qc_table$summary <- sapply(1:nrow(qc_table), function(i) {
  row <- qc_table[i,]
  
  if (any(grepl(")", row))) {
    row <- gsub(".*\\((.*)%)","\\1",row)
    row <- gsub("(.*) \\(.*bp)","\\1",row)
  }
  row <- gsub("bp","",row)
  row <- as.numeric(row)
  
  plot[i,] <<- row
  
  r <- round(range(row, na.rm = TRUE),2) # Range
  m <- round(mean(row, na.rm = TRUE),2) # Mean
  na <- round(mean(is.na(row))*100,0) # Missing 
  
  return(paste0("Mean: ", m, "; Range: ", paste(r, collapse = "-"), "; Missing: ", na, "%"))
} )

qc_table <- qc_table[,names(qc_table)[c(ncol(qc_table), 1:(ncol(qc_table)-1))]]

write.csv(qc_table, paste0(opt$outdir,"/QC_peaks.csv"))

# Plot
print('Plotting QC results')
plot["Treatment",] <- metadata[match(gsub('_.*','',colnames(plot)),metadata$Sample),"Treatment"]
counts$samples$group <- as.factor(as.character(plot["Treatment",paste0(gsub(".*/counts/","",rownames(counts$samples)), '_full')]))

plot <- plot %>% 
  t() %>%
  data.frame(stringsAsFactors = FALSE) %>%
  pivot_longer(1:(nrow(plot)-1), names_to = "perc", values_to ="val")
# Correct value types
plot$perc <- gsub("\\.", " ", plot$perc)
plot$perc <- as.factor(plot$perc)
plot$val <- as.numeric(plot$val)
plot$Treatment <- as.factor(plot$Treatment)
  
plot <- plot[!is.na(plot$val),]
# Set values for each plots
rows <- list(c("BAM records analyzed", "Fragments analyzed", "MAPQ   30",#"Paired alignments",
               "Unpaired alignments"),
             c(#"Paired aln sets duplicates",
               "Singleton aln sets duplicates"),
             c("Peaks identified"),
             c("FRiP"))
labels <- list(c("Read Summary", "", "Reads"),
               c("Duplicate summary", "", "Percentage"),
               c("Peak summary", "Experiments", "Number of Peaks Identified"),
               c("Fraction of reads in peaks", "Experiments", "Fraction"))

# Plot summaries
for (i in 1:length(rows)) {
  print(paste("Plot", labels[[i]][1]))

  idx <- unique(as.numeric(unlist(lapply(rows[[i]], grep, plot$perc))))
  p <- ggplot(plot[idx,] , aes(x=perc, y=val, fill=Treatment)) + 
      geom_boxplot() + 
      labs(title=labels[[i]][1]) + xlab(labels[[i]][2]) + ylab(labels[[i]][3]) +
      theme_classic() + theme( plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  png(paste0(opt$outdir,"/QC_peaks",i,".png"))
  print(p)
  dev.off()
}

# Check samples
print("Samples with most peaks identified:")
peaks <- qc_table["Peaks identified",-1]
peaks <- gsub(" \\(.*","",peaks)
peaks <- as.numeric(peaks)

samples <- metadata[metadata$Sample%in% gsub('_.*','',names(qc_table))[which(peaks > 50000)],c("Strain","Treatment")]
samples <- cbind(samples, "Peaks"=peaks[which(peaks > 50000)])
  
print(samples)

# Check peak sizes
features_original <- read.delim(gsub("gtf","bed", feature), header = FALSE)
features_processed <- read.delim(feature, header = FALSE)

features_original <- features_original$V3 - features_original$V2
features_processed <- features_processed$V5 - features_processed$V4    

png(paste0(opt$outdir,"/QC_peaks_sizes.png"))

p1 <- hist(log(features_original),breaks = seq(0,15,0.1), plot=FALSE)
p2 <- hist(log(features_processed), breaks = seq(0,15,0.1), plot=FALSE)

p1$counts <- p1$counts / sum(p1$counts)
p2$counts <- p2$counts / sum(p2$counts)

plot(0,0,type="n",xlim=c(0,10),ylim = c(0, max(p1$counts, p2$counts)), xlab="log(bp)",ylab="Fraction of peaks",main="Comparison of peaks sizes")

plot(p1,border = FALSE, col=rgb(0,0,1,1/4),add=TRUE)
plot(p2,border = FALSE,col=rgb(0,1,0,1/4),add=TRUE)

abline(v=log(seq(50,250,50)), col = rgb(1,0,0,0.25))

legend('topleft',c('Called peaks','Processed peaks', "50bp intervals"),
       fill = c(rgb(0,0,1,1/4), rgb(0,1,0,1/4), rgb(1,0,0,1/4)), bty = 'n',
       border = NA)
dev.off()

sessionInfo()