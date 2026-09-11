# Script to produce figure S1e

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure1Se.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("ggplot2")
  library("ggrepel")
  library("tidyverse")
  library("openxlsx")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--S2"), type="character", default=NULL, 
              help="Table S2", metavar="character"),
  make_option(c("-b", "--S3"), type="character", default=NULL, 
              help="Table S3", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$S2)){
  print_help(opt_parser)
  stop("Table S2 (-a) is missing", call.=FALSE)
}
if (is.null(opt$S3)){
  print_help(opt_parser)
  stop("Table S3 (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load Table S2\n")
# Get differential results
dar <- read.xlsx(opt$S2, sheet="Differential_Accessibility")
deg <- read.xlsx(opt$S2, sheet="Differential_Expression")

sig_rna <- deg$Gene_ID[deg$RNA_FDR < 0.05]

correlations <- dar[dar$ATAC_FDR < 0.05 & dar$Pearson_correlation_all != "-",
                    c("Region_ID","Nearest.Gene","Pearson_correlation_all","Correlation_FDR_all")]

rownames(dar) <- dar$Region_ID
rownames(deg) <- deg$Gene_ID

correlations$RNA <- deg[correlations$Nearest.Gene,"RNA_logFC"]
correlations$ATAC <- dar[correlations$Region_ID,"ATAC_logFC"]

correlations$alpha <- as.factor(ifelse(as.numeric(correlations$Correlation_FDR_all) < 0.05,"Correlated","Not correlated"))

correlations <- correlations %>% select(RNA,ATAC,Nearest.Gene,Pearson_correlation_all) %>%
  group_by(Nearest.Gene) %>% 
  summarize(RNA=unique(RNA), ATAC= log2(mean(2^ATAC)), cor=max(abs(as.numeric(Pearson_correlation_all))))

cat("Determining gene groups\n")
ov <- read.xlsx(opt$S3, "Overlap")
plasticity <- c(grep('ynap',ov$Description),
                grep('memo',ov$Description),
                grep('hor',ov$Description),
                grep('neu',ov$Description))
plasticity <- unique(unlist(strsplit(ov[plasticity,'Genes'],', ')))

IEG <- unique(c("Atf3", "Bhlhe40", "Ccl2", "Ccn1", "Ccn2", "Ccnl1", "Cebpd", "Csrnp1", "Cxcl3", "Dusp1", "Dusp5", "Dusp6", "Egr1", "Egr3", "F3", "Flg", "Fos", "Fosb", "Gadd45b", "Gbp1", "Gem", "Hbegf", "Ier3", "Il6", "Jun", "Junb", "Klf10", "Klf6", "Klhl21", "Ldlr", "Mcl1", "Nfkbia", "Nfkbiz", "Nr2c2", "Nr4a1", "Nr4a2", "Plau", "Pmaip1", "Rcan1", "Sgk1", "Slc2a3", "Srf", "Tnfaip3", "Trib1", "Tsc22d1", "Zfp36","Fos", "Fosb","Fosl1", "Fosl2","Jun", "Junb","Jund", "Egr1","Egr2", "Egr3","Egr4", "Nr4a1","Nr4a2", "Nr4a3","Arc", "Homer1","Rheb", "Rgs2","Plk2", "Ptgs2","Bdnf", "Inhba","Nptx2", "Plat","Nrn1", "Myc","Dusp1", "Dusp5","Dusp6", "Pcdh8","Cyr61", "Gadd45b","Trib1", "Gem","Btg2", "Ier2","Npas4", "Rasd1","Crem", "Mbnl2","Arf4", "Gadd45g","Arih1", "Nup98","Ppp1r15a", "Fbxo33","Per1", "Per2","Maff", "Zfp36","Srf", "Mcl1","Ctgf", "Il6","Atf3", "Rcan1","Ncoa7", "Cxcl2","Bhlhe40", "Slc2a3","Nfkbia", "Ier3","Sgk1", "Klf6","Klf10", "Nfkbiz","Flg", "Gbp2b","Tnfaip3", "Cebpd","Hbegf", "Ldlr","Tsc22d1", "F3","Ccl2", "Csrnp1","Pmaip1", "Zfp36l2","Plau", "Ccl5","Saa3", "Ifnb1","Tnf", "Irf1","Cd83", "Map3k8","Socs3", "Csf2","Il1a", "Cxcl1","Il12b", "Il1b","Sod2", "Pim1","Peli1", "Tlr2","Ccl3", "Noct","Bcl3", "Ifit2","Icam1", "Ifit1","Tnfsf9", "Ccrl2","Cxcl10", "Gbp2","Il10", "Clec4e","Acod1", "Mmp13","Cxcl11", "Il23a","Arhgef3", "Serpine1","Traf1", "Vcam1","Ackr4", "Marcksl1","Nfkbid", "Ikbke","Ccl12", "Ifit3","Cebpb", "Zfp36l1","Txnip", "Nfib","Hes1", "Pias1","Klf2", "Cd69","Dusp2", "Wee1","Thbs1", "Sik1","Gdf15", "Ier5","Rgs1", "Id2","Apold1"))

print("IEG in plasticity")
print(plasticity[plasticity %in% IEG])
#[1] "Arc"    "Homer1" "Rheb"   "Cebpb"  "Fos"    "Npas4"  "Ptgs2"  "Rcan1" 
#[9] "Egr2"   "Nr4a1"  "Nr4a3" 

ov <- unique(unlist(sapply(ov$Genes,strsplit, ', ')))

cat("Set plot parameters\n")
correlations$fill <- "other"
correlations$fill[correlations$Nearest.Gene %in% ov] <- "Enrichment"
correlations$fill[correlations$Nearest.Gene %in% ov[ov %in% IEG]] <- "IEG"
correlations$col <- correlations$fill
correlations$col[correlations$Nearest.Gene %in% plasticity] <- "Neuronal enrichment terms"
correlations$alpha <- 0.25
correlations$alpha[correlations$Nearest.Gene %in% ov] <- 0.7
correlations$pch <- "No differential expression (RNA)"
correlations$pch[correlations$Nearest.Gene %in% sig_rna] <- "Differential expressed (RNA)"
correlations$pch[correlations$Nearest.Gene %in% plasticity] <- "Neuronal enrichment terms"
correlations$size <- abs(correlations$cor)*5

correlations$label <- correlations$Nearest.Gene
correlations$label[correlations$col == "other"] <- ""
correlations$label[correlations$pch =="No differential expression (RNA)"] <-  ""

cat("Save plot\n")
svg(paste0(opt$outdir,"/FigS1e.svg"), width=18, height=7)

# filter non-relevant dots forligher images
params <- correlations$fill == 'other' & correlations$col == 'other' 
rm_random <- sample(which(params),sum(params) *0.5)
                                                                                            
ggplot(data=correlations[-rm_random,], aes(x=ATAC, y=RNA, fill=fill,col=col, label =label))+
  geom_point(aes(shape=pch), alpha=correlations$alpha[-rm_random], size=correlations$size[-rm_random]) + 
  geom_vline(xintercept=0,color="grey") + geom_hline(yintercept=0,color="grey") + 
  geom_text_repel(show.legend = FALSE, size=5,max.overlaps = 50)  + theme_bw() +
  xlab("Average ATAC logFC") + ylab("RNA logFC") +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values=c("other"="black", 'Enrichment'='blue',"IEG"="purple")) +
  scale_color_manual(values=c("other"="black", "Neuronal enrichment terms" ='darkred', 'Enrichment'='blue',"IEG"="purple")) +
  scale_shape_manual(values=c("Differential expressed (RNA)"=19, "Neuronal enrichment terms" = 22,"No differential expression (RNA)"=4)) 

dev.off()

write.csv(correlations, paste0(opt$outdir, "/data/S1e.csv"))

sessionInfo()