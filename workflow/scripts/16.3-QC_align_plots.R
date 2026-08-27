# Script designed to plot alignment QC summaries

# Redirect all R logs to Snakemake log
log <- file('logs/16-NRF_bam/alignQC_16.3.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("tidyr")
  library("Rsamtools")
  library("ggplot2")
  library("reshape2")
  library("foreach")
  library("doParallel")
  library("dplyr")     
  library("ComplexHeatmap")
  library("optparse")
})
source("workflow/scripts/16.0-QC_align_helper.R")

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="Directory with the RData aligment QC files", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="Processed metadata file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("Directory containing the Fastp .json output (-d) is missing", call.=FALSE)
}
if (is.null(opt$meta)){
  print_help(opt_parser)
  stop("Processed metadata fiile (-m) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
qc_data <- list.files(opt$dir, pattern="RData", full.names= TRUE)
qcs <- lapply(qc_data, function(x) {load(x); return(res)})
names(qcs) <- gsub("Aligned.RData", "", basename(qc_data))

metadata <- read.csv(opt$meta, row.names=1)

# Bam QC plots
mapqc <- lapply(names(qcs), function(sample) {
  x <- qcs[[sample]]$QC$MAPQ 
  x$sample <- sample
  return(x)
  })
mapqc <- do.call(rbind,mapqc)
id <-gsub(".RData","",mapqc$sample)
stopifnot(all(metadata$SampleName[match(id, metadata$SampleName)] == id))
mapqc$group <- metadata$Group[match(id, metadata$SampleName)]

mapqc$alpha <-  0.1
mapqc$Var1 <- as.numeric(mapqc$Var1)

png(paste0(opt$outdir,"/QC_bam_mapq.png"))
ggplot(mapqc , mapping = aes(x = Var1, y = Freq, group=sample, color=group, alpha=alpha)) +
  geom_line() + 
  labs(title="Aligment Quality Distribution") + xlab("Mapping quality") + ylab("Density") +
  theme_classic() + 
  guides(alpha = "none") +
  theme( plot.title = element_text(hjust = 0.5))
dev.off()

# Chromosome distribution
chrs <- lapply(names(qcs), function(sample) {
  x <- qcs[[sample]]$QC$idxstats
  main_chr <- grep("^chr[0-9XYM]*$", x$seqnames) 
  others <- x[-main_chr,]
  x$seqnames <- as.character(x$seqnames)
  x <- rbind(x[main_chr,], c("seqnames" = "others", colSums(others[,-1])))
  x$rate <- as.numeric(x$mapped)/sum(as.numeric(x$mapped))

  x <- data.frame("chr" = x$seqnames, "sample" = sample, "rate"=x$rate)
  return(x)
  })
chrs <- do.call(rbind, chrs)
chrs <- as.data.frame(chrs, stringsAsFactor=FALSE)

id <-gsub(".RData","",chrs$sample)
stopifnot(all(metadata$SampleName[match(id, metadata$SampleName)] == id)) 
chrs$group <- metadata$Group[match(id, metadata$SampleName)]

chrs$chr <- factor(chrs$chr, levels = c(paste0("chr", 1:19), "chrX", "chrY", "chrM", "others"))

png(paste0(opt$outdir,"/QC_bam_chrs.png"))
ggplot(chrs , mapping = aes(x=chr, y = rate, fill = group)) +
  geom_boxplot() + 
  labs(title="Chromosomal alignment rate") + ylab("Proportion of reads") + xlab("") +
  theme_classic() + 
  guides(alpha = 'none') +
  geom_hline(yintercept=0.3) + 
  geom_hline(yintercept=0.5) + 
  theme( plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()

mito <- chrs[chrs$chr=='chrM',]
mito$sample <- gsub('[0-9]*_(.*)_full','\\1',mito$sample)

png(paste0(opt$outdir,"/QC_bam_mitocrondria.png"))
ggplot(mito , mapping = aes(x=sample, y = rate, fill=group)) +
  geom_boxplot() + geom_point() + 
  labs(title="Chromosomal alignment rate") + ylab("Proportion of reads") + xlab("") +
  theme_classic() + 
  guides(alpha = "none") +
  theme( plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()


# Bottleneck
bneck <- lapply(names(qcs), function(sample) {
  x <- qcs[[sample]]$QC$PCRbottleneckCoefficient_1 
  x <- c("coef" = x, "sample"=sample)
  return(x)
  })
bneck <- do.call(rbind,bneck)
bneck <- as.data.frame(bneck, stringsAsFactor=FALSE)

id <-gsub(".RData","",bneck$sample)
stopifnot(all(metadata$SampleName[match(id, metadata$SampleName)] == id)) 
bneck$group <- metadata$Group[match(id, metadata$SampleName)]

bneck$coef <- as.numeric(bneck$coef)

png(paste0(opt$outdir,"/QC_bam_bottleneck.png"))
ggplot(bneck , mapping = aes(y = coef, group=group, color=group)) +
  geom_boxplot() + 
  labs(title="Aligment Bottleneck Distribution") + ylab("Density") +
  theme_classic() + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  guides(alpha = "none") +
  geom_hline(yintercept=0.8) + geom_hline(yintercept=0.9) +  geom_hline(yintercept=0.5) + 
  theme( plot.title = element_text(hjust = 0.5))
dev.off()
 
print("SAMPLES with high bottleneck")
print(bneck[bneck$coef < 0.8, c("coef", "sample")])

# Complexity plots 
comp <- lapply(names(qcs), function(sample) {
  x <- qcs[[sample]]$Complexity 
  x$sample <- sample
  return(x)
  })
comp <- as.data.frame(do.call(rbind, comp))

id <-gsub(".RData","",comp$sample)
stopifnot(all(metadata$SampleName[match(id, metadata$SampleName)] == id)) 
comp$group <- metadata$Group[match(id, metadata$SampleName)]
comp$alpha <-  0.1

png(paste0(opt$outdir,"/QC_bam_libcomplexity.png"))
ggplot(comp , mapping = aes(x = reads, y = values, group=sample, color=group, alpha = alpha)) +
  geom_point(shape = 1, size=2) + 
  geom_line() + 
  labs(title="Estimation of ATAC-seq\nlibrary complexity") + xlab(expression(Putative ~ sequenced ~ 
       fragments ~ x ~ 10^6)) + ylab(expression(Distinct ~  fragments ~ x ~ 10^6)) +
  theme_classic() + scale_alpha_continuous(guide ='none') + 
  theme( plot.title = element_text(hjust = 0.5)) 
dev.off()

# Read size distributions
rsize <- lapply(names(qcs), function(sample) {
  x <- qcs[[sample]]$Length 
  x$sample <- sample
  return(x)
  })
rsize <- as.data.frame(do.call(rbind, rsize))

id <-gsub(".RData","",rsize$sample)
stopifnot(all(metadata$SampleName[match(id, metadata$SampleName)] == id))  
rsize$group <- metadata$Group[match(id, metadata$SampleName)]
rsize$alpha <-  0.1

png(paste0(opt$outdir,"/QC_bam_fragsize.png"))
ggplot(rsize, aes(x=x, y=y, group=sample, color=group, alpha=alpha)) + 
  geom_line() + 
  labs(title="Read sizes") + xlab("Read length (bp)") + ylab(expression(Normalized ~ read ~ density ~ x ~ 10^-3)) +
  theme_classic() +  theme( plot.title = element_text(hjust = 0.5)) +scale_alpha_continuous(guide ='none') 
dev.off()

# PT Scores
pt <- lapply(names(qcs), function(sample) {
  x <- unique(qcs[[sample]]$PT)
  x$sample <- sample
  return(x)
  })
pt <- as.data.frame(do.call(rbind, pt))

id <-gsub(".RData","",pt$sample)
stopifnot(all(metadata$SampleName[match(id, metadata$SampleName)] == id))  
pt$group <- metadata$Group[match(id, metadata$SampleName)]

pt$alpha <-  0.1
pt$alpha[grepl("chrM", pt$chr)] <-  0.5

pt$shape <- sapply(pt$chr, function(x) if(x=="chrM") "Mitocondrial" else "Nuclear")
pt$log <- round(pt$log,1) 
pt$PT <- round(pt$PT,1) 
points <- unique(pt[,-3])

print("PT score:")
print(summary(pt$PT))

png(paste0(opt$outdir,"/QC_bam_pt.png"))
ggplot(points, mapping = aes(x = log, y = PT,  color=group, alpha=alpha, shape = shape)) +
  geom_point(size=1) +  facet_wrap(~group, ncol= 1) + 
  labs(title="PT score") + xlab("log2 mean coverage") + ylab("Promoter vs Transcript")  +
  geom_hline(yintercept=0) + 
  theme_classic() + scale_alpha_continuous(guide ='none') +
  scale_shape_manual(name = "DNA", values = c(4, 1)) +
  theme( plot.title = element_text(hjust = 0.5))
dev.off()

# NFR Scores
nfr <- lapply(names(qcs), function(sample) {
  x <- unique(qcs[[sample]]$NFR)
  x$sample <- sample
  return(x)
  })
nfr <- as.data.frame(do.call(rbind, nfr))

id <-gsub(".RData","",nfr$sample)
stopifnot(all(metadata$SampleName[match(id, metadata$SampleName)] == id))  
nfr$group <- metadata$Group[match(id, metadata$SampleName)]

nfr$alpha <-  0.1
nfr$alpha[grepl("chrM", nfr$chr)] <-  0.5

nfr$shape <- sapply(nfr$chr, function(x) if(x=="chrM") "Mitocondrial" else "Nuclear")
nfr$log <- round(nfr$log,1) 
nfr$NFR <- round(nfr$NFR,1) 
points <- unique(nfr[,-3])

print("NFR score:")
print(summary(nfr$NFR))

png(paste0(opt$outdir,"/QC_bam_nfr.png"))
ggplot(points, mapping = aes(x = log, y = NFR,  color=group, alpha=alpha, shape = shape)) +
  geom_point(size=1) +  facet_wrap(~group, ncol= 1) + 
  labs(title="NFRscore for 200bp flanking TSSs") + xlab("log2 mean coverage") + ylab("Nucleosome Free Regions score") +
  geom_hline(yintercept=0) + 
  theme_classic() + scale_alpha_continuous(guide ='none') +
  scale_shape_manual(name = "DNA", values = c(4, 1)) +
  theme( plot.title = element_text(hjust = 0.5))
dev.off()

# TSSE scores
tsse_line <- lapply(names(qcs), function(sample) {
  x <- qcs[[sample]]$TSS$values
  x <- c(x, "sample" = sample)
  return(x)
  })

tsse_line <- as.data.frame(do.call(rbind, tsse_line))

id <-gsub(".RData","",tsse_line$sample)
stopifnot(all(metadata$SampleName[match(id, metadata$SampleName)] == id))  
tsse_line$group <- metadata$Group[match(id, metadata$SampleName)]
names(tsse_line)[1:(ncol(tsse_line)-2)] <- 100*(-9:10-.5)

plot <- tsse_line %>% 
  data.frame(stringsAsFactors = FALSE) %>%
  pivot_longer(1:(ncol(tsse_line)-2), names_to = "perc", values_to ="val")

plot$perc <- as.numeric(gsub("\\.","-", gsub("X","", plot$perc)))
plot$val <- as.numeric(plot$val)
label <- plot[plot$group=="CTRL" & plot$perc>450 & plot$val>5,]

png(paste0(opt$outdir,"/QC_bam_tsseline.png"))
ggplot(plot, mapping = aes(x = perc, y = val, group=sample, color=group)) +
  geom_line() + geom_point() +  geom_text(data=label, aes(x=perc,y=val+1,label=sample)) + 
  facet_wrap(~group, ncol= 1) + 
  labs(title="Coverage around TSS") + xlab("Distance to TSS") + ylab("Aggregate TSS score") +
  geom_hline(yintercept=0) + 
  theme_classic() + scale_alpha_continuous(guide ='none') +
  scale_shape_manual(name = "DNA", values = c(4, 1)) +
  theme( plot.title = element_text(hjust = 0.5))
dev.off()

tsse_box <- lapply(names(qcs), function(sample) {
  x <- qcs[[sample]]$TSS$TSSEscore
  x <- c(x, "sample" = sample)
  return(x)
  })

tsse_box <- as.data.frame(do.call(rbind, tsse_box))

id <-gsub(".RData","",tsse_box$sample)
stopifnot(all(metadata$SampleName[match(id, metadata$SampleName)] == id))  
tsse_box$group <- metadata$Group[match(id, metadata$SampleName)]
tsse_box$V1 <- as.numeric((tsse_box$V1))

png(paste0(opt$outdir,"/QC_bam_tssebox.png"))
ggplot(tsse_box, mapping = aes(y = V1, group=group, color=group)) +
  geom_boxplot()  + 
  labs(title="Average TSS enrichement") + ylab("TSSE score") +
  geom_hline(yintercept=10, color="darkred") + 
  geom_hline(yintercept=15, color="darkgreen") + 
  theme_classic() +
  theme( plot.title = element_text(hjust = 0.5))
dev.off()

# Sample correlation
mcov <- lapply(names(qcs), function(sample) {
  x <- qcs[[sample]]$VM
  x <- c(x[[1]], "sample" = sample)
  return(x)
  })
mcov <- as.data.frame(do.call(rbind, mcov))
rownames(mcov) <- mcov$sample
mcov <- mcov[,-which(names(mcov)=="sample")]
mcov <- t(mcov)
class(mcov) <- "numeric"

pca <- prcomp(mcov)
plot <- as.data.frame(pca$rotation)

id <- gsub('.RData','',rownames(plot))
stopifnot(all(metadata$SampleName[match(id, metadata$SampleName)] == id))  
plot$group <- metadata$Group[match(id, metadata$SampleName)]

png(paste0(opt$outdir,"/QC_bam_pca.png"))
ggplot(plot, mapping = aes(x = PC1, y = PC2, color=group)) +
  geom_point()  + 
  labs(title="Promoter Coverage PCA") + 
  xlab(paste("PC1 (", round(summary(pca)$importance[2,1]*100,2),"% variace explained)")) + 
  ylab(paste("PC2 (", round(summary(pca)$importance[2,2]*100,2),"% variace explained)")) + 
  theme_classic() + 
  scale_shape_manual(values= rep(1:20,2))+
  stat_ellipse(aes(group = group)) + 
  theme( plot.title = element_text(hjust = 0.5))
dev.off()

id <- gsub('.RData','',colnames(mcov))
stopifnot(all(metadata$SampleName[match(id, metadata$SampleName)] == id))
groups <- metadata$Group[match(id, metadata$SampleName)]

pcor <- cor(mcov, method = "spearman")

gcols <- rainbow(length(unique(groups)))
names(gcols) <- unique(groups)

ha = HeatmapAnnotation(
    Groups = groups,
    col = list(Groups=gcols),
    show_legend = c(TRUE, FALSE)
)

png(paste0(opt$outdir,"/QC_bam_heatmap.png"))
Heatmap(pcor, border_gp = gpar(col = "black", lty = 1),
        top_annotation = ha,
        heatmap_legend_param = list(title = "Correlation"),
        show_row_names = FALSE,
        show_column_names = FALSE)
dev.off()

# Sex phenotype (based on  PureCN)
sex <- lapply(names(qcs), function(sample) {
  x <- qcs[[sample]]$sex[[1]]
  x <- c(x[paste0("chr", c(1:19, "X", "Y"))], "sample" =sample)
  return(x)})
sex <- as.data.frame( t(sapply(sex, unlist)), stringsAsFactors = FALSE)
rownames(sex) <- sex$sample

sex$somatic <- rowMeans(apply(sex[,-which(names(sex) %in% c("chrX", "chrY", "sample"))], 2, as.numeric))
sex$prop <- as.numeric(sex$chrX)/as.numeric(sex$chrY)

png(paste0(opt$outdir,"/QC_bam_sex.png"))
par(mfrow = c(3, 1))
hist(sex$prop, main="Sex genotyping", xlab="Chr X/Y coverage ratio")
plot(sex$chrX ~ sex$chrY, main="", xlab="Average coverage chrY", ylab="Average coverage chrX")
for(i in 1:5) abline(a=0,b=i, lty="dashed", col="grey")
hist(as.numeric(sex$chrX)/sex$somatic, main="",  xlab="Chr X/somatic coverage ratio")
dev.off()

# Write alignment QC table
res <- list()
for (sample in names(qcs)) {
  res[[sample]] <- as.data.frame(qcs[[sample]][['QC']]$PCRbottleneckCoefficient_1)
  res[[sample]]$sample <-  gsub('.RData','',sample)
  res[[sample]]$TSSscore <- qcs[[sample]][['TSSE']]$TSSEscore
}
res <- do.call(rbind,res)

write.csv(res, paste0(opt$outdir,"/QCbam_align.csv"))

sessionInfo()