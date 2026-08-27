# Function to make random colors for each strain

make_colors <- function(group, parents) {
  
  colors <- rainbow(length(levels(group)))
  names(colors) <- levels(group)
  
  if (length(parents) == 2) {
    colors["DBA"]<-parents[["DBA"]]
    colors["C57Bl6"]<-parents[["C57Bl6"]]
  }
  
  colors <- unname(sapply(group, function(x) colors[x]))
  return(colors)
}

plot_helper <- function(col, pch) {
  # Fix DBA and C57 colors
  f0_color <- list("DBA"=rgb(226/255,187/255,144/255), "C57Bl6"=rgb(148/255,151/255,152/255))
  # Get color and shapes
  colors <- make_colors(as.factor(col), parents = f0_color)
  points <- as.numeric(gsub("CTRL", 19,  gsub("SD", 17, pch )))
  
  return(list("colors" = colors, "shapes" = points, "parents"=f0_color))
}