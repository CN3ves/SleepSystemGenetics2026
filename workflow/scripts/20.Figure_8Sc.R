# Script to produce figure S8c

# Redirect all R logs to Snakemake log
log <- file('logs/20-Figures/figure8Sc.log', open = "wt")
sink(log, type = "output")
sink(log, type = "message")

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("VennDiagram")
  library("RColorBrewer")
  library("openxlsx")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')

option_list = list(
  make_option(c("-a", "--S9"), type="character", default=NULL, 
              help="Table S9", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$S9)){
  print_help(opt_parser)
  stop("Table S9 (-a) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Load Table S9 \n")
ct <- read.xlsx(opt$S9, sheet="SD effect -Genotype Control")
fv <- read.xlsx(opt$S9, sheet="SD effect -Tamoxifen Control")
ft <- read.xlsx(opt$S9, sheet="SD effect -ciKO")

stats <- list("Tamoxifen control"=ct, "Genotype control"=fv, "Nrf1 cKO"=ft)
sets <- lapply(stats, function(x) x$ENSEMBL[x$FDR < 0.05])

venn.diagram(
  x = sets,
  category.names = names(sets),
  filename = paste0(opt$outdir,"/FigS8c.svg"),
  output=TRUE,
  disable.logging = TRUE,
  # Output features
  imagetype="svg" ,
  height = 7, 
  width = 7, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = brewer.pal(3, "Pastel2"),
  # Numbers
  cex = .9,
  fontface = "italic",
  fontfamily = "sans",
  # Set names
  cat.cex = 2,
  cat.col = brewer.pal(3, "Pastel2"),
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)
n <- max(sapply(sets, length))
df <- data.frame(c(sets[[1]],rep('',n-length(sets[[1]]))),
    c(sets[[2]],rep('',n-length(sets[[2]]))),
    c(sets[[3]],rep('',n-length(sets[[3]]))))

names(df) <- names(sets)

write.csv(df, paste0(opt$outdir, "/data/S8c.csv"))

sessionInfo()
