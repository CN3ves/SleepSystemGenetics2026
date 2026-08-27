# Script designed for differential analyses

# Redirect all R logs to Snakemake log
log <- file('logs/17-NRF_differential/DA_17.1.log', open = "wt")
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
  library("optparse")

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
pData(counts)$Comment <- factor(pData(counts)$Comment, levels=c("wt","mutant"))
design3 <- model.matrix(~Comment*Sleep, data=pData(counts))
colnames(design3) <- gsub("Comment", "", 
                         gsub("Sleep", "", 
                             gsub(":", "_", colnames(design3))))


y <- DGEList(counts=counts(counts), group=pData(counts)$Group)
y$offset <- -offst(counts)

cat("Filtering counts\n")
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

logcpm <- cpm(y, log=TRUE)

colors <- c("CT_NSD"="grey", "FT_NSD"="red", "FV_NSD"="green", "CT_SD"="darkgrey", "FT_SD"="darkred", "FV_SD"="darkgreen")
  
png(paste0(opt$outdir,"/MDS.png"))
plotMDS(logcpm, col = colors[pData(counts)$Group], pch = ifelse(pData(counts)$Sleep=="SD", 16,23),main="MDS (PCA)",top=5000, gene.selection="common", dim.plot = c(1,2))
legend("topright", legend=c("SD","CTRL"), pch=c(16,23))
legend("topleft", legend=names(colors), col=colors, pch=19, ncol=1)
dev.off()

cat("Estimating disertion\n")
y1 <- estimateDisp(y, design1)
png(paste0(opt$outdir,"/BCV1.png"))
plotBCV(y1)
dev.off()

y2 <- estimateDisp(y, design2)
png(paste0(opt$outdir,"/BCV2.png"))
plotBCV(y2)
dev.off()

y3 <- estimateDisp(y, design3)
png(paste0(opt$outdir,"/BCV3.png"))
plotBCV(y3)
dev.off()

cat("Differential analyses\n")
fit1 <- glmQLFit(y1, design1)
png(paste0(opt$outdir,"/QLDisp1.png"))
plotQLDisp(fit1)
dev.off()

fit2 <- glmQLFit(y2, design2)
png(paste0(opt$outdir,"/QLDisp2.png"))
plotQLDisp(fit2)
dev.off()

fit3 <- glmQLFit(y3, design3)
png(paste0(opt$outdir,"/QLDisp3.png"))
plotQLDisp(fit3)
dev.off()

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

cat("Get statistics\n")
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
   geom_text_repel(show.legend = FALSE, size=3) + 
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

saveRDS(stats, paste0(opt$outdir,"/nrf_stats.RData"))

cat("Heatmap\n")

stats3 <- lapply(1:ncol(contrasts3), function(c) {
  qlf <- glmQLFTest(fit3, contrast=contrasts3[,c])

  tab <-  as.data.frame(topTags(qlf,n=30000,p.value=1))
  tab$SYMBOL <- mapIds(org.Mm.eg.db, gsub("\\..*","",rownames(tab)),  "SYMBOL","ENSEMBL")
  tab$NAME <- mapIds(org.Mm.eg.db, gsub("\\..*","",rownames(tab)),  "GENENAME","ENSEMBL")

  return(tab)
})


mat <- data.frame("BSL" = stats3[[1]][rownames(fit3),grep("logFC", names(stats3[[1]]))],
                  "SD"=   stats3[[2]][rownames(fit3),grep("logFC", names(stats3[[2]]))],
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

cat("Response correlation\n")
mat <- data.frame(FT= stats[['FT']]$table[rownames(fit1),"logFC"],
                  FV= stats[['FV']]$table[rownames(fit1),"logFC"],
                  CT= stats[['CT']]$table[rownames(fit1),"logFC"], 
                  row.names = rownames(fit1))

mat <- data.frame(x=c(mat$FT,mat$FT), 
                  y=c(mat$FV,mat$CT),
                  Genes = rep(gsub("\\..*","",rownames(mat)),2), 
                  Legend=c(rep("FV vs FT", nrow(mat)),
                           rep("CT vs FT", nrow(mat))))

labels <- mat[abs(mat$y-mat$x) > 1,"Genes"]
mat$Labels <- mapIds(org.Mm.eg.db, mat$Genes,  "SYMBOL","ENSEMBL")

mat$Legend[!(mat$Genes %in%labels)] <- "not sig"
mat$Labels[mat$Legend == "not sig"] <- NA

png(paste0(opt$outdir,"/FCcor.png"))
ggplot(mat, mapping = aes(x = x, y = y, color=Legend, label=Labels)) +
  geom_point(alpha=0.2)  + 
  labs(title="Correlation with FT") + 
  xlab("Log2FC FT") + 
  ylab("Log2FC CTRL") + 
  theme_classic() + 
  geom_text_repel(show.legend = FALSE) +
  geom_abline(slope=1, intercept = 0) +
  scale_color_manual(values=c("FV vs FT"="red","CT vs FT"= "blue", "not sig"="grey")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1) ) ) + 
  geom_vline(xintercept = 0, color="grey", linetype="dashed")+
  geom_hline(yintercept = 0, color="grey", linetype="dashed")+
  theme( plot.title = element_text(hjust = 0.5))
dev.off()


sessionInfo()