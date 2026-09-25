# Script to produce figure S4a

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure4Sa.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("igraph")
  library("qgraph")
  library("colorspace")
  library("RColorBrewer")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--net"), type="character", default=NULL, 
              help="SD network", metavar="character"),
  make_option(c("-b", "--full"), type="character", default=NULL, 
              help="Full SD network", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$net)){
  print_help(opt_parser)
  stop("SD network (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load SD network\n")
network <- readRDS(opt$net)
full <- readRDS(opt$full)


colors <- list(
  edge = data.frame(legend=unique(E(full)$type), 
    color= rev(brewer.pal(n = length(unique(E(full)$type))+1, name = "Set3")[-2])),# bad yellow
  node = data.frame(legend=unique(V(full)$carac), 
    color= rev(brewer.pal(n = length(unique(V(full)$carac)), name = "Dark2")))
  )

v_cols <- colors[['node']]
e_cols <- colors[['edge']] 

saveRDS(colors, paste0(opt$outdir,"/graph_cols.RDS"))

net_params <- function(nw, e_cols,v_cols) {
  V(nw)$size <- 1
  V(nw)$color <-  unlist(lapply(V(nw)$carac, function(x) v_cols$color[v_cols$legend == x]))
  
  E(nw)$width <- 0.5
  E(nw)$color <- unlist(lapply(E(nw)$type, function(x) e_cols$color[e_cols$legend == x]))
  
  V(nw)$shape <- "circle"
  V(nw)$shape[V(nw)$carac == 'Region'] <- "square"
  V(nw)$shape[V(nw)$carac %in% c('SNP','Sleep')] <- "sphere" 
  V(nw)$shape[V(nw)$carac %in% c('Kinase')] <- "square" 
  V(nw)$shape[grep("TF",V(nw)$carac)] <- "circle" 
  
  idx <- V(nw)$carac %in% c('TF','Transcript','Kinase', 'Sleep')
  V(nw)$label[idx] <- V(nw)$name[idx]

  
  return(nw)
}

full <- net_params(full, e_cols,v_cols)
  
#subset component
components <- components(network)
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(network)[components$membership == biggest_cluster_id]
subnet <- induced_subgraph(network, vert_ids)

edges <- as_data_frame(subnet, what = "edges")
sub_nodes <- seeds <- "Nrf1"
d <- 5
for (i in 1:d) {
  sub_edges <- edges[edges$from %in% sub_nodes | edges$to %in% sub_nodes,]
  sub_nodes <- unique(c(sub_edges$from, sub_edges$to))
  
  if (i < d)  sub_nodes <- sub_nodes[!sub_nodes %in% V(network)$name[V(network)$carac %in% c('Motif')]]

  if (i == d) { # extend from motifs on the last loop
    motifs <- sub_nodes[sub_nodes %in% V(network)$name[V(network)$carac %in% c('Motif')]]
    n <- V(network)$name[!V(network)$carac %in% c('Chromatin')]
    extra_edges <- rbind(edges[edges$from %in% motifs & edges$to %in% n,],
                     edges[edges$to %in% motifs & edges$from %in% n,])

    extra_nodes <- c(extra_edges$from, extra_edges$to)
    extra_nodes <- extra_nodes[!extra_nodes %in% sub_nodes]

    extra_edges <- rbind(edges[edges$from %in% extra_nodes & edges$to %in% n,],
                     edges[edges$to %in% extra_nodes & edges$from %in% n,])
    sub_nodes <- unique(c(sub_nodes,extra_edges$from, extra_edges$to))
  }
  
}

net <- induced_subgraph(full,  V(full)[V(full)$name %in% sub_nodes]) 

edges <- as_data_frame(net, what = "edges")

unique(V(network)$carac)[!unique(V(network)$carac) %in% V(net)$carac]

unique(E(network)$type)[!unique(E(network)$type) %in% E(net)$type]

layout <- qgraph.layout.fruchtermanreingold(as_edgelist(net, names = FALSE), vcount = vcount(net))

plot_graph <- function(net, old_edges, new_edges, e_cols,v_cols, file, seeds=NULL) {
  
  plot_net <- delete_edges(net,E(net)[!E(net)$type %in% c(old_edges,new_edges)])

  E(plot_net)$width[E(plot_net)$type %in% new_edges] <- 2

  if(!is.null(seeds)) V(plot_net)$size[V(plot_net)$name %in% seeds] <- 2

  svg(paste0(opt$outdir,"/Fig",file, ".svg"),width=10,height=10)
  plot(plot_net,vertex.label.cex = ifelse(V(plot_net)$name %in% seeds,1, 0.75), vertex.label.font = 2, vertex.frame.width = 0,edge.curved = 0.25,vertex.label.color="black", layout=layout)
  
  legend("bottomleft", legend=v_cols$legend, col = v_cols$color, bty = "n", pch=20 , pt.cex = 1, cex = 0.5, text.col=v_cols$color , horiz = FALSE, inset = c(0.1, 0.1))

  legend("topleft", legend=e_cols$legend, col = e_cols$color, bty = "n", lty=1 , pt.cex = 1, cex = 0.5, text.col=e_cols$color, horiz = FALSE, inset = c(0.1, 0.1))
  
  dev.off()
}

# A) QTL x FC edges
cat("Plot figure S5a\n")
old_edges <- c(""); new_edges <- c("QTLxSD","QTL")
plot_graph(net, old_edges, new_edges, e_cols,v_cols,file="S5a", seeds) 

# B) + Cor edges
cat("Plot figure S5b\n")
old_edges <- c(old_edges,new_edges); new_edges <- c("-Cor","+Cor")
plot_graph(net, old_edges, new_edges, e_cols, v_cols, file="S5b", seeds) 

# C) + GRN edges
cat("Plot figure S5c\n")
old_edges <- c(old_edges,new_edges); new_edges <- c("regulatory")
plot_graph(net, old_edges, new_edges, e_cols, v_cols, file="S5c", seeds) 

# D) + footprint edges
cat("Plot figure S5d\n")
old_edges <- c(old_edges,new_edges); new_edges <- c("Footprint")
plot_graph(net, old_edges, new_edges, e_cols, v_cols,file="S5d", seeds) 

# E) + alls edges
cat("Plot figure S5e\n")
old_edges <- c(old_edges,new_edges); new_edges <- c( "Overlap", "Translation", "Phosphorylation")
plot_graph(net, old_edges, new_edges, e_cols, v_cols, file="S5e", seeds) 

# F) network
cat("Plot figure S5f\n")
subnet <- delete_edges(net,  E(net)[!attr(E(net),'vnames') %in% attr(E(network),'vnames')]) 
old_edges <- c(old_edges,new_edges); new_edges <-""
plot_graph(subnet, old_edges, new_edges, e_cols,v_cols, file="S5f", seeds) 

write.csv(as_data_frame(subnet, what = "edges"), paste0(opt$outdir, "/data/S5.csv"))

sessionInfo()
