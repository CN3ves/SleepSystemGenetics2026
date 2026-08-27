# Script containing the helper functions to integrate RNA/ATAC results

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

makeQTL_ranges <-function(qtls, cut_off=0.05){

  qtls <- qtls[grep(cut_off,qtls$test),]

  idx <- lapply(qtls, function(x) grep("\\|-\\|", x))

  for(col in names(idx)) {
    if(length(idx[[col]]) == 0) next

    qtls[idx[[col]],col] <- sapply(idx[[col]], function(i) {
      #remove the multiple tests done during QTL
      idx <- grep(cut_off,strsplit(qtls[i,"test"],"\\|-\\|")[[1]])
      value <- strsplit(qtls[i,col],"\\|-\\|")[[1]][idx]

      return(value)
    })
  }

  qtls <- GRanges(seqnames = qtls$chr, 
                  IRanges(start=as.numeric(qtls$ci_lo)*1e6, end=as.numeric(qtls$ci_hi)*1e6), # position in centimorgans, so multiplied to remove decimals and pesky rounding
                  region = qtls$lodcolumn,
                  pos = qtls$pos,
                  snp = qtls$snp)
  
  return(qtls)
}

classify_qtl <- function(data) {
  #filter QTLs for significant opposed effects (inverted) or only significant effect (protected) of SD on the two genotypes. Only significant baseline differences are excluded
  data$type <- "rm"
  data$type[(data$'BB.SDvsCTRL' == "up" & data$'DD.SDvsCTRL' == "down") | (data$'BB.SDvsCTRL' == "down" &   data$'DD.SDvsCTRL' == "up")] <- "inverted"
  data$type[(data$'BB.SDvsCTRL' == "ns" & data$'DD.SDvsCTRL' != "ns") | (data$'BB.SDvsCTRL' != "ns" &   data$'DD.SDvsCTRL' == "ns")]  <- "protected"
  
  data$type[data$'baseline.DDvsBB' != "ns"]  <- "rm"
  
  return(data)
}

match2coords <- function(lst, ranges,filter=FALSE) {
  if(filter) {
    lst <- lst[lst$type != 'rm',]
  }
  ranges <- data.frame(ranges)[,c("seqnames", "start", "end", "region")]

  lst <- merge(lst,ranges, by="region", all.x =TRUE)
  lst <- makeGRangesFromDataFrame(lst, keep.extra.columns = TRUE,seqnames.field="seqnames")
  
  
  return(lst)
}
