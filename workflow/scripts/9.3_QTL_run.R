# Script designed to run r/qtl2 regression model

cat("Loading packages\n")
.libPaths('Rlibs')

suppressMessages({
  library("qtl2")
  library("optparse")
})

cat("Checking arguments\n")
options(bitmapType='cairo')
option_list = list(
  make_option(c("-j", "--json"), type="character", default=NULL, 
              help="json rqtl2 control file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$json)){
  print_help(opt_parser)
  stop("json rqtl2 control file (-j) is missing", call.=FALSE)
}
if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory (-o) is missing", call.=FALSE)
}

cat("Reading inputs\n")
json <- opt$json
exp <- gsub('.json','',basename(json))

# Redirect all R logs to Snakemake log
log <- file(paste0('logs/9-BXD_qtl/qtl_9.3_', exp,'.log'), open = "wt")
sink(log, type = "output")
sink(log, type = "message")

print(paste("Reading", json))
qtl_data <- read_cross2(file =  json)

summary(qtl_data)

#To find QTL at positions between markers (so called “pseudomarkers”), first insert pseudomarkers into the genetic map 
map <- insert_pseudomarkers(map=qtl_data$gmap, step=1) # spaced at 1 cM intervals 

pr <- calc_genoprob(cross=qtl_data, map=map, error_prob=0.0001)

#In calculating the kinship matrix, you can eliminate the effect of varying marker density across the genome
grid <- calc_grid(map = qtl_data$gmap, step=1)
pr_grid <- probs_to_grid(probs = pr, grid = grid)

# Calculating A Kinship Matrix - Linear mixed models (LMMs)
# https://doi.org/10.1093/g3journal/jkab131 suggests Inflation in the RIL population was corrected by an LOCO kinship correction 

kinship <- calc_kinship(probs = pr_grid)

png(paste0(opt$outdir, "/kinship", exp, ".png"))
heatmap(kinship, symm = TRUE)
dev.off()

kinship <- calc_kinship(probs = pr_grid, type="loco")

treatment <- gsub(".*_", "", names(qtl_data$is_female))
treatment <- factor(treatment, levels = c("CTRL", "SD"))

if (all(is.na(treatment)) ) { # There's no treatment for phenotype
  treatment <- rep("1", length(treatment))
}

treatment <- matrix(as.numeric(treatment)-1, ncol=1, dimnames = list(names(qtl_data$is_female), "Treatment"))

print("Running addictive model")
out_addictive <- scan1(genoprobs = pr,  pheno = qtl_data$pheno, addcovar = treatment,  kinship=kinship)

print("Saving results")
file <- paste0(opt$outdir,"/",exp, "_QTL_model.RData")

saveRDS(list(
    qtl_data=qtl_data,
    genoprobs = pr,  
    pheno = qtl_data$pheno,
    addcovar = treatment,
    kinship=kinship, 
    model=out_addictive), 
  file)

sessionInfo()