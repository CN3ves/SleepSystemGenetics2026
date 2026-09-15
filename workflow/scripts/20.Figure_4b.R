# Script to produce figure 4b

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure4b.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("igraph")
  library("qgraph")
  library("RColorBrewer")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--net"), type="character", default=NULL, 
              help="GR network", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$net)){
  print_help(opt_parser)
  stop("GR network (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load SD network\n")
network <- readRDS(opt$net)

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

network <- net_params(network, e_cols,v_cols)
  
#subset component
components <- components(network)
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(network)[components$membership == biggest_cluster_id]
subnet <- induced_subgraph(network, vert_ids)

edges <- as_data_frame(subnet, what = "edges")

seed <-  readRDS(paste0(opt$outdir,"/sg_genes.RDS"))
seed <- seed[seed %in% V(subnet)$name]

#extend net
sub_nodes <- lapply(seed, function(n) names(unlist(shortest_paths(subnet,n,seed)$vpath)))
sub_nodes <- unique(unlist(sub_nodes))

for(i in 1:2) {
  sub_edges <- edges[edges$from %in% sub_nodes | edges$to %in% sub_nodes,]
  sub_nodes <- unique(c(sub_edges$from, sub_edges$to))

}
subnet <- induced_subgraph(subnet,  V(subnet)[V(subnet)$name %in% sub_nodes]) 

V(subnet)$label <- V(subnet)$name
V(subnet)$label[!V(subnet)$carac =='TF'] <- ""

V(subnet)$size <- ifelse(V(subnet)$carac %in% "TF", 3,1)

# prune network progressively edge chromatin and motifs

d <- names(which(degree(subnet)==1))
rm <- d[d %in% V(subnet)$name[V(subnet)$carac %in% 'Chromatin']]
subnet <- delete_vertices(subnet, rm)

layout <- qgraph.layout.fruchtermanreingold(as_edgelist(subnet, names = FALSE), vcount = vcount(subnet))

svg(paste0(opt$outdir,"/Fig4b.svg"),width=10,height=10) 
plot(subnet,vertex.label.cex =1, vertex.label.font = 2, vertex.frame.width = 0,edge.curved = 0.25,vertex.label.color="black", layout=layout)
  
legend("bottomleft", legend=v_cols$legend, col = v_cols$color, bty = "n", pch=20 , pt.cex = 1, cex = 0.5, text.col=v_cols$color , horiz = FALSE, inset = c(0.1, 0.1))

legend("topleft", legend=e_cols$legend, col = e_cols$color, bty = "n", lty=1 , pt.cex = 1, cex = 0.5, text.col=e_cols$color, horiz = FALSE, inset = c(0.1, 0.1))
  
dev.off()

write.csv(as_data_frame(subnet, what = "edges"), paste0(opt$outdir, "/data/4b.csv"))

saveRDS(V(subnet)$name,paste0(opt$outdir,"/gr_genes.RDS"))

sessionInfo()

