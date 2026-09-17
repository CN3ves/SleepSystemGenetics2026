# Script to produce figure 8a

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure8a.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("tidyverse")
  library("colorspace")
  library("openxlsx")
  library("qgraph")
  library("RColorBrewer")
  library("igraph")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--net"), type="character", default=NULL, 
              help="Sleep net", metavar="character"),
  make_option(c("-b", "--S9"), type="character", default=NULL, 
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
if (is.null(opt$S9)){
  print_help(opt_parser)
  stop("Table S9 (-b) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load Sleep Network \n")
network <- readRDS(opt$net)
components <- igraph::components(network)
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(network)[components$membership == biggest_cluster_id]
network <- induced_subgraph(network, vert_ids)

nodes <- as_data_frame(network, what = "vertices")
phenos <- nodes[nodes$carac=="Sleep", "name"]

cat("Load statistics \n")

ct_int <- read.xlsx(opt$S9, sheet="-Interaction- SD x Genotype")
ct_int <- ct_int[ct_int$PValue < 0.01,]
fv_int <- read.xlsx(opt$S9, sheet="-Interaction- SD x Tamoxifen")
fv_int <- fv_int[fv_int$PValue < 0.01,]

int <- fv_int$SYMBOL[fv_int$ENSEMBL %in% ct_int$ENSEMBL]
int <- int[!is.na(int)]

ct_bsl <- read.xlsx(opt$S9, sheet="-Baseline- ciKO vs Genotype")
ct_bsl <- ct_bsl[ct_bsl$FDR < 0.05,]
fv_bsl <- read.xlsx(opt$S9, sheet="-Baseline- ciKO vs Tamoxifen")
fv_bsl <- fv_bsl[fv_bsl$FDR < 0.05,]

bsl <- fv_bsl$SYMBOL[fv_bsl$ENSEMBL %in% ct_bsl$ENSEMBL]
bsl <- bsl[!is.na(bsl)]

int <- int[int %in% nodes$name]
bsl <- bsl[bsl %in% nodes$name]
genes <- unique(c(int,bsl))

keep_phenos <- c(grep('the',V(network)$name[V(network)$carac %in% "Sleep"], value=T),
                 grep('TPF',V(network)$name[V(network)$carac %in% "Sleep"], value=T),
                 grep('del',V(network)$name[V(network)$carac %in% "Sleep"], value=T),
                 grep('\\.r',V(network)$name[V(network)$carac %in% "Sleep"], value=T),
                 grep('BXDr',V(network)$name[V(network)$carac %in% "Sleep"], value=T))
keep_phenos <-  grep('ref',keep_phenos, invert = TRUE, value=T)

colors <- readRDS(paste0(opt$outdir,"/graph_cols.RDS"))
v_cols <- colors[['node']]
e_cols <- colors[['edge']]

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

network <- net_params(network, e_cols,v_cols)

paths <- find_paths(keep_phenos, genes)
paths <-  paths[sapply(paths,length) <5]

seeds <- unique(unlist(paths))

sub_nodes <- c(seeds,
  unique(unlist(sapply(genes, function(g) names(unlist(shortest_paths(network, g, phenos)$vpath))))))

subnet <- induced_subgraph(network, V(network)[V(network)$name %in% sub_nodes]) 

#Highlight paths
e <- as_edgelist(subnet)
for (path in paths) {
  
  for(i in 2:length(path)){
    idx <- which(e[,1] %in% path[(i-1):i] & e[,2] %in% path[(i-1):i])
    E(subnet)$width[idx] <- 3
  }
  
}

V(subnet)$label <- V(subnet)$name
V(subnet)$label[!V(subnet)$name %in% seeds] <- ""
V(subnet)$label[!V(subnet)$name %in% genes & !V(subnet)$carac %in% c("TF", "Kinase","Sleep")] <- ""

cat("Save plot \n")
layout <- qgraph.layout.fruchtermanreingold(as_edgelist(subnet, names = FALSE), vcount = vcount(subnet))

svg(paste0(opt$outdir, "/Fig8a.svg"),width=10,height=10)
plot(subnet,vertex.label.cex = 1, vertex.label.font = 2, vertex.frame.width = 0,vertex.label.color=V(subnet)$color,edge.curved = 0.25, layout=layout)
  
legend("bottomright", legend=v_cols$legend, col = v_cols$color, bty = "n", pch=20 , pt.cex = 1, cex = 0.5,  text.col=v_cols$color , horiz = FALSE, inset = c(0.1, 0.1))

legend("topleft", legend=e_cols$legend, col = e_cols$color, bty = "n", lty=1 , pt.cex = 1, cex = 0.5, text.col=e_cols$color, horiz = FALSE, inset = c(0.1, 0.1))

dev.off()

write.csv(as_data_frame(subnet, what = "edges"), paste0(opt$outdir, "/data/8a.csv"))

sessionInfo()
