# Script to produce table S7

# Redirect all R logs to Snakemake log
log <- file('logs/19-Tables/S7_19.7.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("biomaRt")
  library("edgeR")
  library("GRaNIE")
  library("openxlsx")
  library("org.Mm.eg.db")
  library("GenomicRanges")
  library("readr")
  library("igraph")
  library("optparse")
})

cat("Checking arguments\n")

option_list = list(
  make_option(c("-a", "--granie"), type="character", default=NULL, 
              help="GRaNIE network", metavar="character"),
  make_option(c("-b", "--atac"), type="character", default=NULL, 
              help="ATAC counts", metavar="character"),
  make_option(c("-c", "--rna"), type="character", default=NULL, 
              help="RNA counts", metavar="character"),
  make_option(c("-d", "--footprints"), type="character", default=NULL, 
              help="Detected footprints", metavar="character"),
  make_option(c("-e", "--S2"), type="character", default=NULL, 
              help="Table S2", metavar="character"),
  make_option(c("-f", "--S4"), type="character", default=NULL, 
              help="Table S4", metavar="character"),
  make_option(c("-g", "--S5"), type="character", default=NULL, 
              help="Table S5", metavar="character"),
  make_option(c("-i", "--S6"), type="character", default=NULL, 
              help="Table S6", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$granie)){
  print_help(opt_parser)
  stop("GRaNIE network (-a) is missing", call.=FALSE)
}
if (is.null(opt$atac)){
  print_help(opt_parser)
  stop("ATAC counts (-b) is missing", call.=FALSE)
}
if (is.null(opt$rna)){
  print_help(opt_parser)
  stop("RNA counts (-c) is missing", call.=FALSE)
}
if (is.null(opt$S2)){
  print_help(opt_parser)
  stop("Table S2 (-d) is missing", call.=FALSE)
}
if (is.null(opt$S4)){
  print_help(opt_parser)
  stop("Table S4 (-e) is missing", call.=FALSE)
}
if (is.null(opt$S5)){
  print_help(opt_parser)
  stop("Table S5 (-f) is missing", call.=FALSE)
}
if (is.null(opt$S6)){
  print_help(opt_parser)
  stop("Table S6 (-e) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

# Create a new workbook and add a sheet
wb <- createWorkbook()

cat("Loading GRaNIE network\n")
gr_net <- readRDS(opt$granie)
ATACdiff <- read.xlsx(opt$S2, "Differential_Accessibility")
RNAdiff <- read.xlsx(opt$S2, "Differential_Expression")
footprints <-  read.xlsx(opt$S4, "Footprints")
kea_list <- list(read.xlsx(opt$S5, "Mean_rank"),
  read.xlsx(opt$S5, "Integrated_scaled_rank"))

edges <- rbind(
  data.frame(source=toupper(as.character(gr_net$TF.name)), target=as.character(gr_net$peak.ID), type="regulatory"),
  data.frame(source=as.character(gr_net$peak.ID), 
             target=ifelse(as.character(gr_net$gene.name) != '',
                           as.character(gr_net$gene.name),
                           as.character(gr_net$gene.ENSEMBL)), 
             type="regulatory"))


nodes <- rbind(data.frame(name = toupper(as.character(gr_net$TF.name)), carac= "TF"), 
               data.frame(name = as.character(gr_net$gene.name), carac= "Transcript"), 
               data.frame(name = as.character(gr_net$gene.ENSEMBL), carac= "Transcript"), 
               data.frame(name = as.character(gr_net$peak.ID), carac= "Chromatin"))

nodes <- nodes[nodes$name %in% c(edges[,1],edges[,2]),]

grnet <- graph_from_data_frame(d=unique(edges), vertices=unique(nodes), directed=F) 

saveRDS(grnet ,paste0(opt$outdir,"/data/grnet.RData"))

addWorksheet(wb, "Gene Regulatory (GR) Net")
writeData(wb, "Gene Regulatory (GR) Net", gr_net, rowNames=FALSE)

cat("Loading QTL network\n")
qtls <- list("RNA_baseline" =  read.xlsx(opt$S6,sheet="RNA_QTL"),
             "RNA_response" = read.xlsx(opt$S6,sheet="RNAxSD_QTL"),
             "ATAC_baseline" = read.xlsx(opt$S6,sheet="ATAC_QTL"),
             "ATAC_response" = read.xlsx(opt$S6,sheet="ATACxSD_QTL"),
             "TF_baseline" = read.xlsx(opt$S6,sheet="TF_QTL"),
             "TF_response" = read.xlsx(opt$S6,sheet="TFxSD_QTL"),
             "Pheno_baseline" = read.xlsx(opt$S6,sheet="pheno_QTL"))

edges <- lapply(names(qtls), function(x) {
  temp <- qtls[[x]]
  temp$loci <- paste0('chr',temp$chr,':',temp$bp-10000,'-',temp$bp+10000)
  col <- names(temp)[names(temp) %in% c("Gene", "region", "Motif_TF", "phenotype_ID", "regionID")]
  temp <- temp[,c("loci","snp",col)]
  
  if(ncol(temp) == 4){ #ATAC
    temp <- temp[,c("loci","snp", "region", "regionID")]
  } else {
    temp$ID <- temp[,col]
  }
  temp$type <- ifelse(grepl("baseline",x), "QTL", "QTLxSD")
  names(temp) <- c("source", 'snp',  "target", "ID", "type")
  
  return(temp)
})

edges <- do.call('rbind', edges)

cat("Calculating correlations\n")
mean_cpm <- function(obj) {
  
  mean_counts <- lapply(unique(obj$samples$group), function(group) {
    samples <- obj$samples$group == group
    
    group <- gsub("BXD0*", "BXD", group)
    group <- gsub("96","48a", group)
    group <- gsub("97","65a", group)
    group <- gsub("103","73b", group)
    
    samples_means <- data.frame(rowMeans(cpm(obj)[,samples, drop = FALSE]))
    names(samples_means) <- group
    
    return(samples_means)
  })
  
  as.data.frame(mean_counts)
}

atac_cpm <- mean_cpm(readRDS(opt$atac))
rna_cpm <- mean_cpm(readRDS(opt$rna))
rna_cpm <- as.data.frame(rna_cpm[,colnames(atac_cpm)])
atac_cpm <- as.data.frame(atac_cpm[,colnames(rna_cpm)])

cors <- ATACdiff[,c("Region_ID","Nearest.Gene","Pearson_correlation_all","Correlation_FDR_all")]
cors$region <- paste(GRanges(ATACdiff))
cors <- cors[cors$Correlation_FDR_all != '-',]
cors$Correlation_FDR_all <- as.numeric(cors$Correlation_FDR_all)
cors <- cors[cors$Correlation_FDR_all < 0.1,]

genes <- unique(edges$ID)[unique(edges$ID) %in% rownames(rna_cpm)]

for (gene in genes) {
  # Get regions sharing QTL peak
  loci <- GRanges(unique(edges[edges$ID  %in% gene,"source"]))
  regions <- unique(edges[queryHits(findOverlaps(GRanges(edges$source), loci)),"ID"])
  regions <-  regions[regions %in% rownames(atac_cpm)]
  
  # Get correlated/annotated region
  if(any(cors$Nearest.Gene %in% gene)) regions <- unique(c(regions, cors[cors$Nearest.Gene %in% gene,"Region_ID"]))
  
  for (region in regions) {
    
    message('\r', round(which(genes ==  gene)/length( unique(genes))*100,2), "% done"," - ",round(which(regions ==  region)/length(regions)*100,2),rep(" ",10), appendLF = FALSE)
    
    correlation <- cor.test(as.numeric(rna_cpm[gene,]),as.numeric(atac_cpm[region,]))
    
    if (correlation$p.value < 0.05) {
      id <- ATACdiff[ATACdiff$Region_ID == region,]
      id <- paste0(id$Chromosome,":", id$Start,"-",id$End)
      
      if(correlation$estimate >0 ) {
        edges <- rbind(edges,data.frame("source"= id, snp="", "target"= gene, "ID" = "cor","type" = "+Cor"))
      } 
      if(correlation$estimate < 0 ) {
        edges <- rbind(edges,data.frame("source"= id, snp="", "target"= gene,  "ID" = "cor","type" = "-Cor"))
      } 
    }
  }
  
}

save.image(paste0(opt$outdir,"/data/sgNet.Rimage"))

edges$ID[-grep('^chr[0-9XY]*_', edges$ID)] <- "-"
names(edges)[4] <- "regionID"

addWorksheet(wb, "Systems Genetics (SG) Net")
writeData(wb, "Systems Genetics (SG) Net", unique(edges), rowNames=FALSE)

cat("Make Genetic network\n")
nodes <- lapply(names(qtls), function(n) {
  temp <- qtls[[n]]
  col <- names(temp) %in%  c("Gene", "region", "Motif_TF", "phenotype_ID")
  data.frame(name=temp[,col], carac=gsub("_.*","",n))
})
nodes[['snp']] <- data.frame(name=unique(unlist(sapply(qtls,function(x) x$snp))), carac='Loci')
nodes[['regions']] <- data.frame(name=unique(edges[grep("Cor",edges$type),"source"]), carac='ATAC')

nodes <- do.call(rbind,nodes)

nodes$carac <- gsub("RNA","Transcript",
                    gsub("Pheno","Sleep",
                    gsub("ATAC","Chromatin",nodes$carac)))



nodes$name[nodes$carac=='Loci'] <- edges$source[match(nodes$name[nodes$carac=='Loci'],edges$snp)] 

sgnet <- graph_from_data_frame(d=unique(edges[,c('source', 'target', 'type')]), vertices=unique(nodes), directed=F) 
saveRDS(sgnet,paste0(opt$outdir,"/data/sgnet.RData"))

cat("Merge networks\n")
edges <-  unique(rbind(as_data_frame(sgnet,"edges"),
                as_data_frame(grnet,"edges")))

nodes <-  unique(rbind(as_data_frame(sgnet,"vertices"),
                as_data_frame(grnet,"vertices")))



cat("Make TF names uniform\n") 
rename_nodes <-function(edges,nodes, newnames)  {
  
  edges$from[edges$from %in% rownames(newnames)] <-  newnames[edges$from[edges$from %in% rownames(newnames)],'n']
  
  edges$to[edges$to %in% rownames(newnames)] <-  newnames[edges$to[edges$to %in% rownames(newnames)],'n']
  
  nodes$name[nodes$name %in% rownames(newnames)] <-  newnames[nodes$name[nodes$name %in% rownames(newnames)],'n']
  
  return(list(edges=unique(edges),nodes=unique(nodes)))
  
}

tfnames <- unique(data.frame(n=toupper(gsub("MA.*\\.","",
                                                    gsub("\\(.*","",
                                                         nodes[nodes$carac=='TF','name']))),
                             old= nodes[nodes$carac=='TF','name'],
                             row.names=nodes[nodes$carac=='TF','name']))

net <- rename_nodes(edges,nodes,tfnames)

cat("Connect TF to transcript\n") 
tfs <- unique(unlist(strsplit(net[["nodes"]]$name[net[["nodes"]]$carac == "TF"],"::")))

tfs <- toupper(net[["nodes"]]$name[net[["nodes"]]$carac == "Transcript"]) %in% toupper(tfs)
tfs <- net[["nodes"]]$name[net[["nodes"]]$carac == "Transcript"][tfs]

net[["edges"]] <- unique(rbind(net[["edges"]], 
               unique(data.frame(from= tfs,
                                 to=toupper(tfs), 
                                 type="Translation"))))

net[["nodes"]] <- unique(rbind(net[["nodes"]], data.frame(name = toupper(tfs), carac= "TF")))

cat("Identify TF dimers\n") 
dimers <- lapply(grep("::",net[["nodes"]]$name,value=TRUE), function(dimer) {
  TFs <- strsplit(dimer,"::")[[1]]
  idx <-  TFs %in% net[['nodes']]$name[net[['nodes']]$carac=="TF"]
  
  if (any(idx)) {
    res <- data.frame(from= rep(dimer, 2),to=TFs, type=rep("Heterodimer",sum(idx)))
    
    if (sum(idx)==2) { # add heterodimer binding to regions where both TFs bind
      tf1 <- c(net[['edges']]$from[net[['edges']]$to %in% TFs[1]],
               net[['edges']]$to[net[['edges']]$from %in% TFs[1]])
      tf2 <- c(net[['edges']]$from[net[['edges']]$to %in% TFs[2]],
               net[['edges']]$to[net[['edges']]$from %in% TFs[2]])
      shared <- tf1[tf1 %in% tf2]
      if(length(shared)>0) {
        res <- rbind(res, data.frame(from=dimer,to=tf1[tf1 %in% tf2], type="Binding"))
      }
    }
    return(res)
  }
  
})

dimers <- do.call(rbind,dimers)

net[["edges"]] <- unique(rbind(net[["edges"]],dimers))
net[["nodes"]] <- unique(rbind(net[["nodes"]], data.frame(name = dimers[dimers$type=="Heterodimer","to"], carac= "TF")))

cat("Filter footprints to TFs in network\n") 
footprints <- footprints[toupper(footprints$motif_TF) %in% net[["nodes"]]$name,]
footprints <- GRanges(footprints)

cat("Match footprints to regions in network\n") 
regions <- GRanges(net[["nodes"]]$name[net[["nodes"]]$carac %in% "Chromatin"])

ov <- findOverlaps(footprints,regions)
#TF to motif
net[["edges"]] <- rbind(net[["edges"]], 
               unique(data.frame(from=toupper(gsub(".*\\.","",
                                                   gsub("\\(.*","",footprints$motif_TF[queryHits(ov)]))),
                                 to=footprints$motif_ID[queryHits(ov)],
                                 type="Footprint")))
#Motif to region
net[["edges"]] <- rbind(net[["edges"]], 
               unique(data.frame(from=footprints$motif_ID[queryHits(ov)],
                                 to=paste(regions)[subjectHits(ov)],
                                 type="Footprint")))

net[["nodes"]] <- unique(rbind(net[["nodes"]], data.frame(name = footprints$motif_ID[queryHits(ov)], carac= "Motif")))

cat("Rename overlapping regions\n") 
regions <- GRanges(net[['nodes']]$name[net[['nodes']]$carac=="Chromatin"])
reduced <- reduce(regions)
ov <- findOverlaps(regions,reduced)

regnames <- data.frame(n=paste(reduced[subjectHits(ov)]),
                       old=paste(regions[queryHits(ov)]),
                       row.names=paste(regions[queryHits(ov)]))

net <- rename_nodes(net[["edges"]],net[["nodes"]],regnames)

cat("Rename overlapping loci\n") 
loci <- GRanges(net[['nodes']]$name[net[['nodes']]$carac=="Loci"])
reduced <- reduce(loci)
ov <- findOverlaps(loci,reduced)

regnames <- data.frame(n=paste(reduced[subjectHits(ov)]),
                       old=paste(loci[queryHits(ov)]),
                       row.names=paste(loci[queryHits(ov)]))

net <- rename_nodes(net[["edges"]],net[["nodes"]],regnames)

cat("Match Loci and regions\n") 
regions <- GRanges(net[['nodes']]$name[net[['nodes']]$carac=="Chromatin"])
loci <- GRanges(net[['nodes']]$name[net[['nodes']]$carac=="Loci"])
ov <- findOverlaps(loci,regions)

net[["edges"]] <- rbind(net[["edges"]], 
               unique(data.frame(from=paste(loci[queryHits(ov)]),
                                 to=paste(regions)[subjectHits(ov)],
                                 type="Overlap")))


cat("Add Kinases\n") 
kea <- rbind(kea_list[[1]][1:20,-3],kea_list[[2]][1:20,-3])

kinases <- strsplit(kea$Overlapping.Proteins, ",")
names(kinases) <- kea$Protein

#Get only the shared kinases marked with *
kinase_edges <- lapply(grep("\\*",names(kinases), value=TRUE), function(k) {
  data.frame(from=gsub("\\*","",k),to=kinases[[k]], type="Phosphorylation")
})
kinase_edges <- do.call('rbind',kinase_edges)

kinase_edges <- kinase_edges[kinase_edges$to %in% net[["nodes"]]$name[net[["nodes"]]$carac == 'TF'],]

kinases <- net[["nodes"]]$name[net[["nodes"]]$carac == 'Transcript'][toupper(net[["nodes"]]$name[net[["nodes"]]$carac == 'Transcript']) %in% toupper(kinase_edges$from)]

net[["nodes"]] <- unique(rbind(net[["nodes"]], data.frame(name = unique(kinase_edges$from), carac= "Kinase")))

kinase_edges <- rbind(kinase_edges,
                      data.frame(from=kinases, to=toupper(kinases), type="Translation"))

net[['edges']] <- unique(rbind(net[['edges']],kinase_edges))

fullnet <- graph_from_data_frame(d=unique(net[['edges']]), vertices=unique(net[['nodes']]), directed=F) 
saveRDS(fullnet,paste0(opt$outdir,"/data/fullnet.RData"))

addWorksheet(wb, "Unfitered Net")
writeData(wb, "Unfitered Net", unique(as_data_frame(fullnet,"edges")), rowNames=FALSE)

save.image(paste0(opt$outdir,"/data/fullnet.Rimage"))

cat("Filter network\n") 
network <- delete_vertices(fullnet, grep("FitNormal",V(fullnet)$name, value=TRUE))

nodes <- as_data_frame(network,"vertices")

cat("Filter transcripts affected by SD\n") 
rna <- RNAdiff[RNAdiff$RNA_FDR < 0.05,"Gene_ID"]
keep <- nodes$name[nodes$carac == 'Transcript'] %in% rna

network <- delete_vertices(network,nodes$name[nodes$carac == 'Transcript'][!keep])

cat("Filter regions affected by SD\n") 
atac <- GRanges(ATACdiff[ATACdiff$ATAC_FDR < 0.05,])
regions <- GRanges(nodes$name[nodes$carac == 'Chromatin'])
ov <- findOverlaps(regions,atac)
keep <- nodes$name[nodes$carac == 'Chromatin'] %in% paste(regions[queryHits(ov)])

network <- delete_vertices(network,nodes$name[nodes$carac == 'Chromatin'][!keep])

cat("Filter regions not connected to SNPs or transcripts\n") 
nodes <- as_data_frame(network,"vertices")
edges <- as_data_frame(network, "edges")

regions <- nodes$name[nodes$carac == 'Chromatin']
genes <- nodes$name[nodes$carac == 'Transcript']
loci <- nodes$name[nodes$carac == 'Loci']

r <- edges$to[edges$to %in% regions & edges$from %in% c(genes, loci)]
r <- c(r,edges$from[edges$from %in% regions & edges$to %in% c(genes, loci)])

keep <- regions %in% r

network <- delete_vertices(network,regions[!keep])

cat("Filter motifs with no regions\n") 
nodes <- as_data_frame(network,"vertices")
edges <- as_data_frame(network, "edges")

regions <- nodes$name[nodes$carac == 'Chromatin']
motifs <- nodes$name[nodes$carac == 'Motif']

m <- edges$to[edges$to %in% motifs & edges$from %in% regions]
m <- c(m,edges$from[edges$from %in% motifs & edges$to %in% regions])

keep <- motifs %in% m

network <- delete_vertices(network,motifs[!keep])

cat("Filter footprints with no regulation\n") 
nodes <- as_data_frame(network,"vertices")
edges <- as_data_frame(network, "edges")

motifs <- nodes$name[nodes$carac == 'Motif']

for (m in motifs) {
  n <- names(neighbors(network, m))
  tf <- n[n %in% nodes$name[nodes$carac == 'TF']]

  r <- edges$to[edges$to %in% n & edges$from %in% tf]
  r <- c(r,edges$from[edges$from %in% n & edges$to %in% tf])

  if (length(r) == 0)  network <- delete_vertices(network,m)
}

tfs <- nodes$name[nodes$carac == 'TF']
rm_edges <- list()
for (t in tfs) {
  n <- names(neighbors(network, t))
  
  motif <- n[n %in% nodes$name[nodes$carac == 'Motif']]

  if(length(motif) > 0) {
    r <- n[n %in% nodes$name[nodes$carac == 'Chromatin']]
    if(length(r) > 0) {
      for (m in motif) {
        f <- names(neighbors(network, m))
        rm <- r[r %in% f]

        for (e in rm) {
          rm_edges <- c(rm_edges, list(
            c(which(edges$to == t & edges$from == e),
              which(edges$to == e & edges$from == t))))         
        }
      }
    }
  } else if(length(motif) == 0) { #try remove all tf-chr
    r <- n[n %in% nodes$name[nodes$carac == 'Chromatin']]
    for (e in r) {
      rm_edges <- c(rm_edges, list(
        c(which(edges$to == t & edges$from == e),
         which(edges$to == e & edges$from == t))))    
    }
  }
}

network <- delete_edges(network,E(network)[unlist(rm_edges)])

cat("Filter unconnected nodes\n") 
nodes <- as_data_frame(network,"vertices")
edges <- as_data_frame(network, "edges")
# prune network progressively edge chromatin and motifs
for(o in c('Chromatin', 'Motif')){ 
  d <- names(which(degree(network)==1))
  rm <- d[d %in% nodes$name[nodes$carac %in% o]]
  network <- delete_vertices(network, rm)
}

network <- delete_vertices(network, names(which(degree(network)==0)))

saveRDS(network, paste0(opt$outdir,"/data/SleepNet.RDS"))

cat("Rename edges\n") 
net <- as_data_frame(network, "both")

E(network)$type[E(network)$type == 'QTL'] <- "Genetic"
E(network)$type[E(network)$type == 'QTLxSD'] <- "Genetic interaction"

atac <- net$vertices$name[net$vertices$carac == "Chromatin"]
rna <- net$vertices$name[net$vertices$carac == "Transcript"]
E(network)$type[(net$edges$to %in% atac & net$edges$from %in% rna) | (net$edges$to %in% rna & net$edges$from %in% atac)] <- "Regulatory"

tf <- net$vertices$name[net$vertices$carac == "TF"]
E(network)$type[(net$edges$to %in% tf & net$edges$from %in% rna) | (net$edges$to %in% rna & net$edges$from %in% tf)] <- "Translation"

kinases <- net$vertices$name[net$vertices$carac == "Kinase"]
E(network)$type[(net$edges$to %in% tf & net$edges$from %in% kinases) | (net$edges$to %in% kinases & net$edges$from %in% tf)] <- "Phosphorylation"

E(network)$type[(net$edges$to %in% kinases & net$edges$from %in% rna) | (net$edges$to %in% rna & net$edges$from %in% kinases)] <- "Translation"

E(network)$type[((net$edges$to %in% tf & net$edges$from %in% atac) | (net$edges$to %in% atac & net$edges$from %in% tf)) & net$edges$type != "Footprint"] <- "Regulatory"

motif <- net$vertices$name[net$vertices$carac == "Motif"]
E(network)$type[net$edges$to %in% motif | net$edges$from %in% motif] <- "Binding"

net$edges$simplified <- E(network)$type

addWorksheet(wb, "Final SD-Net edges")
writeData(wb, "Final SD-Net edges", unique(net$edges) , rowNames=FALSE)

addWorksheet(wb, "Final SD-Net nodes")
writeData(wb, "Final SD-Net nodes", unique(net$vertices) , rowNames=FALSE)

cat("Save table\n")
saveWorkbook(wb, paste0(opt$outdir,"/TableS7-GRN.xlsx"), overwrite = TRUE)
   
saveRDS(network, paste0(opt$outdir,"/data/SleepNet_simple.RDS"))

sessionInfo()

