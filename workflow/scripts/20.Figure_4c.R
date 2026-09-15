# Script to produce figure 4c

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure4c.log', open = "wt")
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
find_paths<- function(group1, group2) {

  paths <- lapply(group1, function(p) {
    paths <- shortest_paths(subnet,p,group2[group2 != p])$vpath
    lengths <- sapply(paths, length)
    return(paths[lengths == min(lengths)])
  })
  paths <- unlist(paths[!sapply(paths, is.null)], recursive=FALSE)

  return(lapply(paths, names))
}

network <- net_params(network, e_cols,v_cols)
  
#subset component
components <- components(network)
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(network)[components$membership == biggest_cluster_id]
subnet <- induced_subgraph(network, vert_ids)

sg <-  readRDS(paste0(opt$outdir,"/sg_genes.RDS"))
gr <-  readRDS(paste0(opt$outdir,"/gr_genes.RDS"))
sg <- sg[sg %in% V(subnet)$name]
gr <- gr[gr %in% V(subnet)$name]
sg <- sg[!sg %in% V(subnet)$name[V(subnet)$carac == 'Chromatin']]
gr <- gr[!gr %in% V(subnet)$name[V(subnet)$carac == 'Chromatin']]

#find paths between phenotypes
sub_nodes <- c(find_paths(sg, gr),find_paths(gr, sg))
sub_nodes <- unique(unlist(sub_nodes))

subnet <- induced_subgraph(subnet,  V(subnet)[V(subnet)$name %in% sub_nodes]) 

phenos <- V(subnet)$name[V(subnet)$carac=="Sleep"]
keep_phenos <- c(grep('the',phenos, value=T),
                 grep('TPF',phenos, value=T),
                 grep('del',phenos, value=T),
                 grep('\\.r',phenos, value=T),
                 grep('BXDr',phenos, value=T))
keep_phenos <-  grep('ref',keep_phenos, invert = TRUE, value=T)

V(subnet)$label <- V(subnet)$name
V(subnet)$label[!(V(subnet)$carac =='TF' | V(subnet)$name %in% keep_phenos)] <- ""

V(subnet)$size <- ifelse(V(subnet)$carac %in% "TF", 3,1)

d <- names(which(degree(subnet)==1))
rm <- d[d %in% V(subnet)$name[V(subnet)$carac %in% 'Chromatin']]
subnet <- delete_vertices(subnet, rm)

layout <- qgraph.layout.fruchtermanreingold(as_edgelist(subnet, names = FALSE), vcount = vcount(subnet))

svg(paste0(opt$outdir,"/Fig4c.svg"),width=10,height=10) 
plot(subnet,vertex.label.cex =1, vertex.label.font = 2, vertex.frame.width = 0,edge.curved = 0.25,vertex.label.color="black", layout=layout)
  
legend("bottomleft", legend=v_cols$legend, col = v_cols$color, bty = "n", pch=20 , pt.cex = 1, cex = 0.5, text.col=v_cols$color , horiz = FALSE, inset = c(0.1, 0.1))

legend("topleft", legend=e_cols$legend, col = e_cols$color, bty = "n", lty=1 , pt.cex = 1, cex = 0.5, text.col=e_cols$color, horiz = FALSE, inset = c(0.1, 0.1))
  
dev.off()

write.csv(as_data_frame(subnet, what = "edges"), paste0(opt$outdir, "/data/4c.csv"))

sessionInfo()

