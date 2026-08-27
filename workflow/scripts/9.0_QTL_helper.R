# Functions to make QTL input files uniform

add_header <- function(original_header, object, shift = 0 ) {
  header <- cbind(original_header, object[1:nrow(original_header), -1])
  header[-1] <- "" 
  header[2+shift,1] <- paste("# nrow", nrow(object)) 
  header[3+shift,1] <- paste("# ncol", ncol(object))
  names(header) <-  names(object)
  return(rbind(header, names(object), object))
}

make_genotype <- function(template, select, sample_names, file) { # make genotype file
  
  geno <-  template[['data']][,c("marker", select)] 
  names(geno) <- c("marker", sample_names)
  geno <- geno[,unique(names(geno))] # remove duplicated names
  
  geno <- add_header(template[['header']], geno)
  
  write.table(geno, file, row.names =FALSE, col.names = FALSE, quote=FALSE, sep=",")
  
  return(geno)
}

make_phenotype <- function(counts, template, file) {

  pheno <- t(counts)
  pheno <- as.data.frame(pheno)
  pheno$id <- rownames(pheno)
  pheno <- pheno[,c("id",grep("^id$", names(pheno), invert = TRUE, value = TRUE))]
  
  pheno <- add_header(template[['header']], pheno)
  
  pheno[is.na(pheno)] <- 0 
  
  write.table(pheno, file , row.names =FALSE, col.names = FALSE, quote=FALSE, sep = ",")
  
  return(pheno)
}

collapse_counts <- function(obj) {
  
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

make_cross <- function(genotype, template, strains, file) {
  
  cross <- data.frame("id" = genotype, "cross_direction"= "BxD")

  cross <- add_header(template[['header']], cross, shift = 1)
  
  write.table(cross, file, row.names =FALSE, col.names = FALSE, quote=FALSE, sep = ",")
  
  return(cross)
}
  
make_covar <- function(phenotype, regions, template, file) {
  covar <- data.frame("id" = phenotype, 
                      "description"= gsub(",","_",regions$description),
                      regions[, -which(names(regions)=="description")])
  covar <- covar[,!sapply(covar, function(col) all(is.na(col)))]
  covar <- add_header(template[['header']], covar)
  
  write.table(covar, file, row.names =FALSE, col.names = FALSE, quote=FALSE, sep = ",")
  
  return(covar)
}

read_template <- function(template_dir) {
  
  template_data <- list("phenotype" = list(),"genotype" = list(),"covariates" = list(),"cross" = list())
  
  # Load phenotype data
  template_data[["phenotype"]][["header"]] <- read.csv(paste0(template_dir,"/bxd_pheno.csv"), nrows=3, header = FALSE)
  template_data[["phenotype"]][["data"]] <- read.csv(paste0(template_dir,"/bxd_pheno.csv"), skip=3)
  
  # Load Genotype data
  template_data[["genotype"]][["header"]] <- read.csv(paste0(template_dir,"/bxd_geno.csv"), nrows=3, header = FALSE)
  geno <- read.csv(paste0(template_dir,"/bxd_geno.csv"), skip=3)
  #geno$C57Bl6 <- geno$C57 <- geno$C571 <- geno$C572 <- "B" #set parent strains
  #geno$DBA <- geno$DBA1 <- geno$DBA2 <- "D"
  template_data[["genotype"]][["data"]] <- geno[, -which(names(geno) == "BXD100")] # Remove problematy strain
  
  # Load covariates data
  template_data[["covariates"]][["header"]] <- read.csv(paste0(template_dir,"/bxd_phenocovar.csv"), nrows=3, header = FALSE)
  template_data[["covariates"]][["data"]] <- read.csv(paste0(template_dir,"/bxd_phenocovar.csv"), skip=3)
  
  # Load cross data 
  template_data[["cross"]][["header"]] <- read.csv(paste0(template_dir,"/bxd_crossinfo.csv"), nrows=4, header = FALSE)
  template_data[["cross"]][["data"]] <- read.csv(paste0(template_dir,"/bxd_crossinfo.csv"),skip=4)
  
  # Load control file
  template_data[["control"]]  <-  fromJSON(file = paste0(template_dir, "/bxd.json"))

  return(template_data)
}

get_regions <- function(exp, genes= NULL) {
  if (grepl("rna", exp)) {    
    gene.coords <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene,single.strand.genes.only=FALSE)
    gene.coords <- as.data.frame(gene.coords)
    
    symb <- as.data.frame(org.Mm.egGENENAME[mappedkeys(org.Mm.egGENENAME)])
    symb$SYMBOL <- select(org.Mm.eg.db, keys = symb$gene_id, columns="SYMBOL", keytype="ENTREZID")$SYMBOL
    
    gene.coords <- merge(gene.coords, symb, all.x=TRUE, by.x="group_name", by.y="gene_id")
    gene.coords <- gene.coords[match(tolower(genes), tolower(gene.coords$SYMBOL) ),]
    print(paste("Proportion of genes matched:", mean(!is.na(gene.coords[,1]))))
    
    regions <- data.frame("description"= gene.coords$gene_name, 
                          "chr"= gene.coords$seqnames,
                          "start"= gene.coords$start,
                          "end"= gene.coords$end)
  } else if(grepl("atac", exp)) {
    features <- read.delim('results/6-BXD_features/list/genomic_features.gtf', header = FALSE)
    features$V9 <- gsub("Peak_ID (.*);", "\\1", features$V9)
    
    print(paste("Proportion of regions matched to feature list:", mean(genes %in% features$V9)))
    feats <- features[match(rownames(mat), features$V9),]
    print(paste("All features selected in data:", all(feats$V9 == rownames(mat))))
    
    regions <- data.frame("description"= gsub("\\..*", "",feats$V9), 
                      "chr"= feats$V1,
                      "start"= feats$V4,
                      "end"= feats$V5)
  } else if(grepl("sleep", exp)) {
    
    metadata <- read.delim("rawdata/sleep_bxd/phenotypes_metadata.txt", header = TRUE)
    
    feature_idx <- sapply(metadata[,1], function(x) {
      grep(paste0("^",gsub("\\.", "\\\\.", gsub("-", ".",gsub("out.", "", x))), "$"), 
           gsub("out.", "",gsub("quant.", "", genes)),
           ignore.case  =TRUE)
    })# names don't match between files
    
    regions <- data.frame("description"= metadata$filename.variable, 
                      "Condition"= metadata$condition,
                      "Type"= paste(metadata$variable.type, metadata$subtype, sep =))
    
    regions$feature <- feature_idx
    regions$feature <- sapply(regions$feature, function(i) if(length(genes[i]) >0) genes[i] else NA)
    print("Phenotypes with no measurements found:")
    print(regions[is.na(regions$feature),"description"])
    regions <- regions[!is.na(regions$feature),]
    
    print("Measurements with no metadata found:")
    print(genes[!genes %in% regions$feature])
    add_mat <- data.frame("description"= NA, 
                      "Condition"= NA,
                      "Type"= NA,
                      "feature" = genes[!genes %in% regions$feature])
    regions <- rbind(regions, add_mat)

  }

  return(regions)
}
