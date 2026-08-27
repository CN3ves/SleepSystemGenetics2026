# Script designed to plot and summarise footprint results

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ggplot2")
  library("ggrepel")
  library("ComplexHeatmap")
  library("RColorBrewer")
  library("dendsort")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-f", "--footprints"), type="character", default=NULL, 
              help="Files with the differential footpritns for each BXD line", metavar="character"),
  make_option(c("-a", "--heatmap"), type="character", default=NULL, 
              help="File with differential footprint with all lines", metavar="character"),
  make_option(c("-t", "--test"), type="character", default=NULL, 
              help="Which group of (down)samples are used", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$footprints)){
  print_help(opt_parser)
  stop("Files with the differential footprints for each BXD line (-f) is missing", call.=FALSE)
}
if (is.null(opt$heatmap)){
  print_help(opt_parser)
  stop("File with differential footprint with all lines (-a) is missing", call.=FALSE)
}
if (is.null(opt$test)){
  print_help(opt_parser)
  stop("Which group of (down)samples are used (-t) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/12-BXD_footprints/plots_12.1',opt$test,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Reading files\n")
files <- strsplit(opt$footprints, " ")[[1]]
footprints <- lapply(files, read.delim)

# Get significance
n <- length(footprints)
footprints <- do.call(rbind, footprints)

sig_footprints <- footprints[footprints$P_values < 0.05,]
tf_sig <- table(sig_footprints$Motif)

cat(paste("Footprints are significantly different between SD and CTRL in all strains for these TF:",
            paste(names(which(tf_sig == n)), collapse = "\n"), "\n"))

# Filter TFs for those with a good number of motifs for statistical robustness
footprints <- footprints[footprints$Num >100, c("Motif", "TF_Activity", "P_values")]

# Aggregate statistics meta-analysis
fisher_res <- list()
for (tf in unique(footprints$Motif)) {
  tf_res <- footprints[footprints$Motif == tf,]
  pvals <- tf_res$P_values
  activities <- tf_res$TF_Activity
  
  # If TFs are missing:
  print(paste("TF", tf, "has", nrow(tf_res), "samples after filtering"))
  if(!length(pvals)>0) {
    pvals <- 1
    activities <- 0
  }

  fisher_res[[tf]] <- c("mean_activity" = mean(activities, na.rm= T), 
                        "meta_pval" = pchisq((sum(log(pvals))*-2), df=length(pvals)*2, lower.tail=F))
  
}

meta_footprints <- as.data.frame(do.call(rbind, fisher_res))
meta_footprints$padj <- p.adjust(meta_footprints$meta_pval, method ="BH")

write.csv(meta_footprints, paste0(opt$outdir,"/footprint_analysis", opt$test,".csv"))

print(paste(sum(meta_footprints$padj < 0.05),"out of", nrow(meta_footprints) ,"TFs have statistically differential activity"))

# Make aggregate plot
df <- meta_footprints[order(meta_footprints$mean_activity),]
df$tf <- rownames(df)
df$Sig <- "Not sig."
df$Sig[df$padj <0.05] <- "Sig."
df$Sig[df$tf %in% names(which(tf_sig == n))] <- "All"

cols <- c("Not sig." = "grey", "Sig." = "darkred", "All" = "darkgreen")
df$labs <- df$tf
df$labs[df$Sig == "Not sig."] <- ""
df$labs <- gsub("MA[0-9]*\\.[0-9]\\.", "", df$labs)
df$labs[df$Sig != "Not sig."] <- paste0(df$labs[df$Sig != "Not sig."], " (", tf_sig[df$tf[df$Sig != "Not sig."]], "/",n, ")")
df$labs <- gsub("NA", "0", df$labs)

y_lims <- round(max(abs(df$mean_activity)),2)+0.01

svg(paste0(opt$outdir,"/agreggated_diff",opt$test,".svg"))
ggplot(df, aes(x=tf,y=mean_activity, fill=Sig, col = Sig, label = labs)) +
  geom_point(shape = 21, size=1) +
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) + 
  theme_classic() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title=element_text(size=14,face="bold")) +
  geom_text_repel(size=3,show.legend = FALSE,max.overlaps=20) +
  ylab("Mean activity score\nCTRL <-> SD") + xlab("") + 
  ylim(c(-y_lims,y_lims)) + 
  geom_hline(yintercept=0, color = "black", linewidth=0.5) 
dev.off()

### Heatmap with all samples ran simultaneouslty ###
if(opt$test=="") { # no need to repeat for subsampled data

  df <- read.delim(opt$heatmap)
  samples <- unique(gsub("Protection_Score_","",names(df)[grep("Protection_Score_",names(df))]))

  for(s in samples) {
    df[,s] <- rowSums(df[,grepl(s, names(df))])
  }

  df <- df[,-c(grep("Protection_Score_", names(df)), grep("TC_", names(df)))]

  rownames(df) <- df$Motif
  df$Motif <- NULL
  df_plot <- subset(df, Num > 100)
  df_plot$Num <- NULL

  df_plot_scale <-t(apply(df_plot, 1, scale))
  colnames(df_plot_scale) <- colnames(df_plot)
  rownames(df_plot_scale) <- rownames(df_plot)

  rowlabels <- rownames(df_plot_scale)  
  rowlabels[!rowlabels %in% rownames(meta_footprints)[meta_footprints$padj < 0.05]] <- "" 
    
  rowlabels <- gsub("MA[0-9]*\\.[0-9]\\.","",rowlabels)
  rowlabels <- gsub("\\(var","",rowlabels)

  ra = rowAnnotation(
    labels = anno_text(rowlabels, which = "row", just="left", gp = gpar(fontsize = 5)),
    show_legend = FALSE
   )

  p <- Heatmap(as.matrix(df_plot_scale),
               name = "TF",
               left_annotation = ra,
               cluster_rows = TRUE,
               show_row_names = FALSE,
               row_dend_side = "right",
               row_names_gp = gpar(fontsize = 4),
               column_names_gp = gpar(fontsize = 8, col=ifelse(grepl("_SD", colnames(df_plot_scale)), "purple","darkgreen")),
               show_column_dend = TRUE,
               show_row_dend = TRUE,
               clustering_method_rows = "ward.D2",
               col = rev(brewer.pal(n = 11, name = "RdBu")))


  png(paste0(opt$outdir,"/heatmap_tf.png"), width = 500, height = 500)
  p <- draw(p)
  p
  dev.off()

  # Split by clusters
  r.dend <- cutree(as.hclust(row_dend(p)),10)  

  for (i in unique(r.dend)) {
    
    col_dend = dendsort(hclust(dist(t(df_plot_scale[r.dend!=i,]))))
    
    p <- Heatmap(as.matrix(df_plot_scale[r.dend==i,]),
               name = "TF",
               left_annotation = ra[r.dend==i],
               cluster_columns = col_dend,
               cluster_rows = TRUE,
               show_row_names = FALSE,
               row_names_gp = gpar(fontsize = 4),
               column_names_gp = gpar(fontsize = 8, col=ifelse(grepl("_SD", colnames(df_plot_scale)), "purple","darkgreen")),
               show_column_dend = TRUE,
               show_row_dend = TRUE,
               clustering_method_rows = "ward.D2",
               col = rev(brewer.pal(n = 11, name = "RdBu")))
    
    png(paste0(opt$outdir,"/heatmap_tf", i, ".png"), width = 500, height = 500)
    draw(p)
    dev.off()
  }
}

sessionInfo()