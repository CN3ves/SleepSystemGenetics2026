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
  library('visNetwork') 
  library('htmlwidgets')
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

layout <- qgraph.layout.fruchtermanreingold(as_edgelist(subnet, names = FALSE), vcount = vcount(subnet))




links <- as_data_frame(plot_net, what = "edges")
links$value <- 2
links$title <- links$type 

nodes <- as_data_frame(plot_net, what = "vertices")
nodes$id <- nodes$title <- nodes$name 
nodes$label <- substr(nodes$label,1,6)

ln <- unique(nodes[,c('carac','shape','color')])
names(ln)[1] <- 'label'

le <- unique(links[,c('title','color')])
names(le)[1] <- 'label'

data <- toVisNetworkData(igraph_network)
visNetwork(nodes = data$nodes, edges = data$edges, height = "500px")

fig <- visNetwork(nodes, links, width="100%", height="800px", background="#eeefff",
           main="SD-Net", submain="the most beautiful of all networks") %>% 
  visOptions(selectedBy = "carac", 
             highlightNearest = list(enabled=TRUE, hover=TRUE, degree=1),
             nodesIdSelection = TRUE) %>%
  visLegend(useGroups = F, addNodes = ln,  addEdges = le)  %>% 
  visInteraction(navigationButtons = TRUE)

saveWidget(fig, file = "myplot.html")

sessionInfo()

