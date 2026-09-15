# Script to produce figure 8c

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure8c.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("edgeR")
  library("GenomicRanges")
  library("RColorBrewer")
  library("tidyverse")
  library("colorspace")
  library("openxlsx")
  library("qgraph")
  library("igraph")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--net"), type="character", default=NULL, 
              help="Sleep net", metavar="character"),
  make_option(c("-b", "--atac"), type="character", default=NULL, 
              help="ATAC-seq counts", metavar="character"),
  make_option(c("-c", "--rna"), type="character", default=NULL, 
              help="RNA-seq counts", metavar="character"),
  make_option(c("-d", "--prints"), type="character", default=NULL, 
              help="BXD footprints", metavar="character"),
  make_option(c("-e", "--phenos"), type="character", default=NULL, 
              help="Sleep phenotypes", metavar="character"),
  make_option(c("-f", "--S2"), type="character", default=NULL, 
              help="Table S2", metavar="character"),
  make_option(c("-g", "--S9"), type="character", default=NULL, 
              help="Table S9", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$net)){
  print_help(opt_parser)
  stop("Sleep net (-a) is missing", call.=FALSE)
}
if (is.null(opt$atac)){
  print_help(opt_parser)
  stop("ATAC-seq counts (-b) is missing", call.=FALSE)
}
if (is.null(opt$rna)){
  print_help(opt_parser)
  stop("RNA-seq counts(-c) is missing", call.=FALSE)
}
if (is.null(opt$prints)){
  print_help(opt_parser)
  stop("BXD footprints (-d) is missing", call.=FALSE)
}
if (is.null(opt$phenos)){
  print_help(opt_parser)
  stop("Sleep phenotypes (-e) is missing", call.=FALSE)
}
if (is.null(opt$S2)){
  print_help(opt_parser)
  stop("Table S2 (-f) is missing", call.=FALSE)
}
if (is.null(opt$S9)){
  print_help(opt_parser)
  stop("Table S9 (-g) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

# Helper functions
net_params <- function(nw, e_cols,v_cols) {
  V(nw)$size <- 1
  V(nw)$color <-  unlist(lapply(V(nw)$carac, function(x) v_cols$color[v_cols$legend == x]))
  
  E(nw)$width <- 1
  E(nw)$color <- unlist(lapply(E(nw)$type, function(x) e_cols$color[e_cols$legend == x]))
  
  V(nw)$shape <- "circle"
  V(nw)$shape[V(nw)$carac == 'Region'] <- "square"
  V(nw)$shape[V(nw)$carac %in% c('SNP','Sleep')] <- "sphere" 
  V(nw)$shape[V(nw)$carac %in% c('Kinase')] <- "square" 
  V(nw)$shape[grep("TF",V(nw)$carac)] <- "circle" 
  
  return(nw)
}

find_paths <- function(group1, group2) {

  paths <- lapply(group1, function(p) {
    paths <- shortest_paths(network,p,group2[group2 != p])$vpath
    idx <- sapply(paths, function(x) any(names(x) %in% V(network)$name[V(network)$carac  == "Motif"]))
    if(!all(idx)) {

      paths <- paths[!idx]
      lengths <- sapply(paths, length)
      return(paths[lengths == min(lengths)])
    }
  })

  paths <- unlist(paths[!sapply(paths, is.null)], recursive=FALSE)

  return(lapply(paths, names))
}
fix_names <- function(group) {
  group <- gsub("BXD0*", "BXD", group)
  group <- gsub("96","48a", group)
  group <- gsub("97","65a", group)
  group <- gsub("103","73b", group)
  
  return(group)
}

mean_fcs <- function(obj) {
  
  means <- lapply(unique(obj$samples$Strain), function(group) {
    ctrl <- obj$samples$group == paste0(group,"_CTRL")
    sd <- obj$samples$group == paste0(group,"_SD")
    
    group <- fix_names(group)
    
    ctrl_means <- rowMeans(cpm(obj)[,ctrl, drop = FALSE])
    sd_means <- rowMeans(cpm(obj)[,sd, drop = FALSE])
    
    fc <- data.frame(sd_means[names(ctrl_means)]/ctrl_means)
    names(fc) <- group
    
    return(as.data.frame(fc[order(rownames(fc)),,drop=FALSE]))
  })

  return(as.data.frame(means))
}

cat("Load Sleep Network \n")
network <- readRDS(opt$net)
components <- igraph::components(network)
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(network)[components$membership == biggest_cluster_id]
network <- induced_subgraph(network, vert_ids)

nodes <- as_data_frame(network, what = "vertices")

cat("Load statistics \n")

ct_int <- read.xlsx(opt$S9, sheet="-Interaction- SD x Genotype")
ct_int <- ct_int[ct_int$PValue < 0.01,]
fv_int <- read.xlsx(opt$S9, sheet="-Interaction- SD x Tamoxifen")
fv_int <- fv_int[fv_int$PValue < 0.01,]

int <- fv_int$SYMBOL[fv_int$ENSEMBLE %in% ct_int$ENSEMBLE]
int <- int[!is.na(int)]

ct_bsl <- read.xlsx(opt$S9, sheet="-Baseline- ciKO vs Genotype")
ct_bsl <- ct_bsl[ct_bsl$FDR < 0.05,]
fv_bsl <- read.xlsx(opt$S9, sheet="-Baseline- ciKO vs Tamoxifen")
fv_bsl <- fv_bsl[fv_bsl$FDR < 0.05,]

bsl <- fv_bsl$SYMBOL[fv_bsl$ENSEMBLE %in% ct_bsl$ENSEMBLE]
bsl <- bsl[!is.na(bsl)]

int <- int[int %in% nodes$name]
bsl <- bsl[bsl %in% nodes$name]
genes <- unique(c(int,bsl))

cat("Load DARs \n")
ATACdiff <- read.xlsx(opt$S2, sheet="Differential_Accessibility" )

cat("Load footprint meta-analysis \n")
prints <- read.delim(opt$prints)
samples <- unique(gsub("Protection_Score_","",names(prints)[grep("Protection_Score_",names(prints))]))
for(s in samples) {
  prints[,s] <- rowSums(prints[,grepl(s, names(prints))])
}
tfs <- prints[,-c(grep("Protection_Score_", names(prints)), grep("TC_", names(prints)))]
tfs$Motif <- toupper(gsub("MA.*\\.","",gsub("\\(.*","",tfs$Motif)))
tfs <- tfs %>% select(!Num) %>% group_by(Motif) %>% summarise_all(mean) %>% as.data.frame
rownames(tfs) <- tfs$Motif
tfs <- tfs[,-1]

cat("Load phenotypes \n")
sleep <- read.delim(opt$phenos, sep=" ")
rownames(sleep) <- fix_names(sleep$ID)
colnames(sleep) <- gsub("quant.out.","",colnames(sleep))
sleep <- t(sleep[,-1])

cat("Prepare Sleep Network \n")
keep_phenos <- c("BXD_TPF.P.freqFD_sws", "BXD_NREMS_EEG_bands.del_fast." )

colors <- readRDS(paste0(opt$outdir,"/graph_cols.RDS"))
v_cols <- colors[['node']]
e_cols <- colors[['edge']]

network <- net_params(network, e_cols,v_cols)

paths <- find_paths(keep_phenos, genes)
paths <-  paths[sapply(paths,length) <5]

sub_nodes <- seeds <- unique(unlist(paths))
edges <- as_edgelist(network)
for (i in 1:3) {
  sub_nodes <- sub_nodes[!sub_nodes %in% V(network)$name[V(network)$carac =="Motif"]]
  sub_edges <- edges[edges[,1] %in% sub_nodes | edges[,2] %in% sub_nodes,]
  sub_nodes <- unique(c(sub_edges[,1], sub_edges[,2]))
}

subnet <- induced_subgraph(network, V(network)[V(network)$name %in% sub_nodes]) 

V(subnet)$label <- V(subnet)$name
V(subnet)$label[V(subnet)$carac %in% "Chromatin"] <- ""

layout <- qgraph.layout.fruchtermanreingold(as_edgelist(subnet, names = FALSE), vcount = vcount(subnet))

lab_genes <- genes[genes %in% unique(unlist(paths))]
lab_phenos <- keep_phenos[keep_phenos %in% unlist(paths)]
lab.cex <- ifelse(V(subnet)$label %in% c(lab_genes,lab_phenos),1.5,1)

e <- as_edgelist(subnet)
for (path in paths) {
  
  for(i in 2:length(path)){
    idx <- which(e[,1] %in% path[(i-1):i] & e[,2] %in% path[(i-1):i])
    E(subnet)$width[idx] <- 3
  }
  
}

cat("Calculate FC responses \n")
atac_fc <- mean_fcs(readRDS(opt$atac))
rna_fc <- mean_fcs(readRDS(opt$rna))

samples <- colnames(sleep)[colnames(sleep) %in% colnames(atac_fc)]
samples <- samples[samples %in% colnames(rna_fc)]

sleep <- as.data.frame(sleep[phenos[phenos %in% V(subnet)$name],samples,drop=FALSE])
rna_fc <- as.data.frame(rna_fc[rownames(rna_fc) %in% V(subnet)$name,samples])

ov <- findOverlaps(GRanges(ATACdiff),GRanges(V(subnet)$name[V(subnet)$carac == 'Chromatin']))

reg_convertion <- data.frame(ID = ATACdiff[queryHits(ov),"Region_ID"], 
                             test = paste(GRanges(ATACdiff[queryHits(ov),])),
                             Region = paste(GRanges(V(subnet)$name[V(subnet)$carac == 'Chromatin'])[subjectHits(ov)]))

atac_fc <- lapply(unique(reg_convertion$Region), function(r) {
  ids <- reg_convertion[reg_convertion$Region == r,"ID"]
  mat <- matrix(colMeans(atac_fc[ids,samples]),nrow=1,  dimnames = list(r,samples))
})

atac_fc <- as.data.frame(do.call('rbind',atac_fc))   

for(group in unique(gsub("_.*","", names(tfs)))) {
  ctrl <- grep(paste0(group,"_CTRL"),names(tfs))
  sd <- grep(paste0(group,"_SD"),names(tfs))
  
  group <- fix_names(group)
  
  fc <- data.frame(tfs[,sd]/tfs[,ctrl])
  names(fc) <- group
  tfs <- cbind(tfs, fc)
  tfs <- tfs[,-c(sd,ctrl)]
    
}

names(tfs)[-grep("BXD", names(tfs))] <- paste0("BXD",names(tfs)[-grep("BXD", names(tfs))])
tfs <- as.data.frame(tfs[rownames(tfs) %in% V(subnet)$name,samples])

nodes <- rbind(rna_fc,atac_fc,tfs)

cors <- sapply(rownames(nodes), function(n) {
  x <- as.numeric(nodes[n,])
  rm <- is.na(x) | x==0 | x ==Inf
  y <- as.numeric(colMeans(sleep[keep_phenos,]))
  
  cor(x[!rm],y[!rm])
})

cuts <- cut(cors, breaks=seq(round(min(cors)-0.04,1),round(max(cors)+0.04,1),0.1))
names(cuts) <- names(cors)
cols <- colorRampPalette(c("darkred", "red", "gold", "cyan","blue"))(length(levels(cuts)))
names(cols) <- levels(cuts)

V(subnet)$label.color <- V(subnet)$color
V(subnet)$color <- as.character(cols[cuts[V(subnet)$name]])
V(subnet)$color[is.na(V(subnet)$color)] <- "darkgray"
V(subnet)$size <- abs(round(cors[V(subnet)$name]*10)+1)
V(subnet)$size[is.na(V(subnet)$size)] <- 2

svg(paste0(opt$outdir, "/Fig8c.svg"),width=10,height=10)
plot(subnet,vertex.label.cex =  lab.cex, vertex.label.font = 2, vertex.frame.width = 0,edge.curved = 0.25, layout=layout)
  
legend("bottomright", legend=v_cols$legend, col = v_cols$color, bty = "n", pch=0 , pt.cex = 1, cex = 0.5, text.col=v_cols$color , horiz = FALSE, inset = c(0.1, 0.1))

legend("topleft", legend=e_cols$legend, col = e_cols$color, bty = "n", lty=1 , pt.cex = 1, cex = 0.5, text.col=e_cols$color, horiz = FALSE, inset = c(0.1, 0.1))

legend("topright", legend=c(paste0("ρ ⊂ ",names(cols)),"No cor"), bty = "n", col =c(cols,"darkgray"), pch=20 , pt.cex = 1, cex = 0.5, text.col=c(cols,"darkgray"), horiz = FALSE, inset = c(0.1, 0.1))

dev.off()

write.csv(as_data_frame(subnet, what = "edges"), paste0(opt$outdir, "/data/8a.csv"))

sessionInfo()
