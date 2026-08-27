# Script designed for differential analyses

# Redirect all R logs to Snakemake log
log <- file('logs/logs/17-NRF_differential/DA.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("edgeR")
  library("EDASeq")
  library("org.Mm.eg.db")
  library("ggplot2")
  library("ggrepel")
  library("tidyverse")
  library('TxDb.Mmusculus.UCSC.mm10.knownGene')
  library("ComplexHeatmap")
  library("circlize")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="Rdata file with the normalized feature counts", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("Rdata file with the normalized feature counts (-c) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
counts <- readRDS(opt$counts)

# Complex interaction contrast
pData(counts)$NRF <- ifelse(gsub("_.*", "", pData(counts)$Group) == "FT", "negative", "positive")

pData(counts)$NRF <- as.factor(pData(counts)$NRF)
pData(counts)$NRF <- relevel(pData(counts)$NRF, ref="positive")
pData(counts)$Treatment <- as.factor(pData(counts)$Treatment)
pData(counts)$Treatment <- relevel(pData(counts)$Treatment, ref="Vehicle")

# design 1 no interactions; design 2 with interactions; design 3 splits group interactions
design1 <- model.matrix(~0+Group, data=pData(counts))
colnames(design1) <- gsub("Group", "", colnames(design1))

pData(counts)$Genotype <- factor(gsub("_.*","", pData(counts)$Group), levels=c("FT", "FV","CT"))
design2 <- model.matrix(~Genotype*Sleep, data=pData(counts))
colnames(design2) <- gsub("Genotype", "", 
                         gsub("Sleep", "", 
                             gsub(":", "_", colnames(design2))))

pData(counts)$Comment <- ifelse(pData(counts)$Genotype == "FT","mutant","wt")
design3 <- model.matrix(~Comment*Sleep, data=pData(counts))
colnames(design3) <- gsub("Comment", "", 
                         gsub("Sleep", "", 
                             gsub(":", "_", colnames(design3))))


y <- DGEList(counts=counts(counts), group=pData(counts)$Group)
y$offset <- -offst(counts)
# Filtering
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

logcpm <- cpm(y, log=TRUE)

colors <- c("CT_NSD"="grey", "FT_NSD"="red", "FV_NSD"="green", "CT_SD"="darkgrey", "FT_SD"="darkred", "FV_SD"="darkgreen")
  
png(paste0(opt$outdir,"/MDS.png")
plotMDS(logcpm, col = colors[pData(counts)$Group], pch = ifelse(pData(counts)$Sleep=="SD", 16,23),main="MDS (PCA)",top=5000, gene.selection="common", dim.plot = c(1,2))
legend("topright", legend=c("SD","CTRL"), pch=c(16,23))
legend("topleft", legend=names(colors), col=colors, pch=19, ncol=1)
dev.off()

#Dispersion
y1 <- estimateDisp(y, design1)
plotBCV(y1)
y2 <- estimateDisp(y, design2)
plotBCV(y2)
y3 <- estimateDisp(y, design3)
plotBCV(y3)

#Differential expression
fit1 <- glmQLFit(y1, design1)
plotQLDisp(fit1)

fit2 <- glmQLFit(y2, design2)
plotQLDisp(fit2)

fit3 <- glmQLFit(y3, design3)
plotQLDisp(fit3)

contrasts1 <- makeContrasts("FT"= FT_SD - FT_NSD,
                           "CT"= CT_SD - CT_NSD,
                           "FV"= FV_SD - FV_NSD, 
                           "T_control"= (FT_SD - FT_NSD) - (CT_SD - CT_NSD),
                           "F_control"= (FT_SD - FT_NSD) - (FV_SD - FV_NSD),
                           "FTcvFV_baseline"= FT_NSD - FV_NSD,
                           "FTcvCT_baseline"= FT_NSD - CT_NSD,
                           levels=design1)
#one-way analysis of deviance (ANODEV) 
# https://f1000research.com/articles/5-1438
# add as "out of range" contrast

contrasts2<- makeContrasts("baseline"= (CT+FV)/2,
                           "baselineCT"= CT,
                           "baselineFV"= FV,
                           "SDFT"= SD,
                           "int"= (CT_SD+FV_SD)/2,
                           "SDCT"= CT_SD,
                           "SDFV"= FV_SD,
                           levels=design2)

contrasts3<- makeContrasts("mutant"= mutant,
                           "mutant*SD"= mutant_SD,
                           levels=design3)

stats <- lapply(1:(ncol(contrasts1)+ncol(contrasts2)+1), function(i) {
 
 if (i > ncol(contrasts1)) {
   if(i > ncol(contrasts1)+ncol(contrasts2)) {
     contrast=contrasts1[,1:3]
     n <- "ANODEV_sleep"
     fit <- fit1
   } else{
     contrast=contrasts2[,i-ncol(contrasts1)]
     n <- colnames(contrasts2)[i-ncol(contrasts1)]
     fit <- fit2
   }
   
 } else {
   contrast=contrasts1[,i]
   n <- colnames(contrasts1)[i]
   fit <- fit1
 }
  
 qlf <- glmQLFTest(fit, contrast=contrast)

 png(paste0(opt$outdir,"/MD_", n,".png"))
 plotMD(qlf)
 abline(h=c(-1, 1), col="blue")
 dev.off()
 
 tab <-  as.data.frame(topTags(qlf,n=30000,p.value=1))
 tab$SYMBOL <- mapIds(org.Mm.eg.db, gsub("\\..*","",rownames(tab)),  "SYMBOL","ENSEMBL")
 tab$NAME <- mapIds(org.Mm.eg.db, gsub("\\..*","",rownames(tab)),  "GENENAME","ENSEMBL")
 
 write.csv(tab, paste0(opt$outdir,"/Diff_", n,".csv"))

 genes <- tab
 
 if (i > ncol(contrasts1)+ncol(contrasts2)) genes$logFC <- rowMeans(genes[,grep("logFC",names(genes))])
 
 genes$color <- "No change"
 genes$color[genes$FDR < 0.05 & genes$logFC > 0] <- "Up"
 genes$color[genes$FDR < 0.05 & genes$logFC < 0] <- "Down"
   
 idx <- unique(c(order(genes$logFC)[1:20],order(genes$logFC, decreasing=TRUE)[1:20],order(genes$FDR)[1:50]))
 genes$label <- NA
 genes$label[idx] <- genes$SYMBOL[idx]
 
 lim <-  round(max(abs(genes$logFC)),1)
    
 g_lines <- c(1.2,1.5,2)
 g_lines <- c(1/g_lines,g_lines)
 g_lines <- log(g_lines)

 png(paste0(opt$outdir,"/Volcano_", n,".png"))

 g <- ggplot(data=genes, aes(x=logFC, y=-log10(FDR), col=color, label = label)) +
   geom_point(size = 0.2, alpha=0.75) + 
   geom_text_repel(show_guide = FALSE, size=3) + 
   theme_classic() + xlim(-lim,lim) + 
   guides(alpha= "none", col=guide_legend(title="DEG", override.aes = list(size=2))) + 
   scale_color_manual(values=c("navy", "grey","#e2725b")) +
   geom_vline(xintercept=g_lines, col="grey", lty="dashed", alpha = 0.5) +
   geom_hline(yintercept=-log10(0.05), col="grey", alpha = 0.75) + 
   ggtitle("Differential expression") +
   theme(plot.title = element_text(hjust = 0.5))
 print(g)
 dev.off()
 
 return(list(result=qlf, summary=summary(decideTests(qlf)), table = tab,contrast = n))
})

names(stats) <- lapply(stats, function(x) unique(x$contrast))

saveRDS(stats, paste0(opt$outdir,"/stats.RData"))


# Correlations
mat <- data.frame("BSL" = rowMeans(stats3[[1]]$table[rownames(fit3),grep("logFC", names(stats3[[1]]$table))]),
                  "SD"=   rowMeans(stats3[[2]]$table[rownames(fit3),grep("logFC", names(stats3[[2]]$table))]),
                  row.names = rownames(fit3))



genes <- AnnotationDbi::select(org.Mm.eg.db, keytype="GOALL", keys=c("GO:0042775","GO:0006390","GO:0140053"), columns=c("SYMBOL","ENSEMBL"))

targets <- c("Nrf1", "Gabrb1", "Tfam",unique(genes$SYMBOL))

labels <- sapply(mapIds(org.Mm.eg.db, gsub("\\..*","",rownames(mat)),  "SYMBOL","ENSEMBL"), function(x) {
  ifelse(x %in% targets, x,NA)})

labels[is.na(labels)] <- ""

ht <- Heatmap(as.matrix(mat), 
             col = colorRamp2(c(-1,  0, 1), c("red", "white", "blue")),
             heatmap_legend_param = list(title = "LogFC SDvsCTRL"),
             row_labels = labels,
             column_title="Condition", row_title="Transcripts") 

png(paste0(opt$outdir,"/heatmap.png"))
draw(ht)
dev.off()


sessionInfo()