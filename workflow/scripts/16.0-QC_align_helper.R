# Script containing the helper functions for alignment QC summary

# Function to read and extract the alignment results
alignQC <- function(logs) {
  log_summary <- data.frame("File"="File", "Reads" = "Reads", "Unaligned" = "Unaligned", "Unaligned.rate" = "Unaligned.rate",  "Unique_map" = "Unique_map", "Unique_map.rate" = "Unique_map.rate", "Multi_map" = "Multi_map","Multi_map.rate"="Multi_map.rate", stringsAsFactors = FALSE)

  # Extract data from logs
  for (i in 1:length(logs)) {
    log <- read.delim(logs[i], header=FALSE, stringsAsFactors = FALSE)
    idx <- c(grep("MULTI-MAPPING READS", log[,1]),
             grep("UNMAPPED READS", log[,1]),
             grep("CHIMERIC READS", log[,1]))
    idx <- list(log[idx[1]:idx[2],2],
             log[idx[2]:idx[3],2])

    sums <- sapply(idx, function(x) sum(as.numeric(grep("%",x, invert=TRUE, value=TRUE)), na.rm=TRUE))
    rates <- sapply(idx, function(x) paste0(sum(as.numeric(gsub("%","",grep("%",x, value=TRUE)))),"%"))
    
    log <- data.frame("File"= gsub(".*/(.*)_Log.final.out","\\1", logs[i]),
                      "Reads" = log[grep("Number of input reads",log[,1]),2],
                      "Unaligned" = sums[2],
                      "Unaligned rate" = rates[2],
                      "Unique_map" = log[grep("Uniquely mapped reads number",log[,1]),2],
                      "Unique_map rate" = log[grep("Uniquely mapped reads %",log[,1]),2],
                      "Multi_map" = sums[1],
                      "Multi_map rate" =rates[1],
                      stringsAsFactors = FALSE)
    
    log_summary <- rbind(log_summary, log)
  }
  names(log_summary) <- log_summary[1,]
  
  # Summarize data
  log_summary[1,] <- c("", apply(log_summary[-1,-1], 2, function(col) {
    millions <- !any(grepl("%", col))
    col <- gsub("%", "", col)
    col <- as.numeric(col)
    col <- summary(col)
    if (millions) col <- round(col/1000000,2)
    col <- paste0("Min: ", col["Min."], "%; Max: ", col["Max."], "%; Mean: ", round(col["Mean"],2),"%")
    if (millions) col <- gsub("%", "M", col)
    return(col)
  }))
  
  return(log_summary)
}


# Function to calculate read duplication (adapted from ATACseqQC)
readsDupFreq <- function (bamFile) {
    lst <- lapply(names(bamFile[[1]]), function(.elt) {
        do.call(c, unname(lapply(bamFile, "[[", .elt)))
    })
    names(lst) <- names(bamFile[[1]])
    df <- do.call(data.frame, lst)
    rm("lst")
    freqByPos <- as.data.frame(table(paste(df$rname, df$strand, 
        df$cigar, df$pos, df$mrnm, df$isize, df$mpos, sep = "_")))
    rm("df")
    freqByDuplication <- as.data.frame(table(freqByPos$Freq))
    freqByDuplication$Var1 <- as.numeric(as.character(freqByDuplication$Var1))
    freqByDuplication <- as.matrix(freqByDuplication)
  
    return(freqByDuplication)
}

# Function to estimate library complexity (adapted from ATACseqQC)
estimateLibComplexity <- function (histFile, times = 100, interpolate.sample.sizes = seq(0.1, 1, by = 0.1), extrapolate.sample.sizes = seq(5, 20, by = 5)) {
    total <- histFile[, 1] %*% histFile[, 2]
    suppressWarnings({
        result = ds.rSAC.bootstrap(histFile, r = 1, times = times)
    })
    sequences <- c(interpolate.sample.sizes, extrapolate.sample.sizes)
    estimates <- data.frame(relative.size = sequences, values = rep(NA, 
        length(sequences)))
    for (i in seq_along(sequences)) {
        estimates$values[i] <- result$f(sequences[i])/10^6
    }
    suppressWarnings(estimates$reads <- estimates$relative.size * 
        total/10^6)
  
    return(invisible(estimates))
}

# Function to estimate mean coverage (adapted from ATACseqQC's plotCorrelation)
mean_cov <- function (obj, txs, seqlev = intersect(seqlevels(obj), seqlevels(txs)), upstream = 2000, downstream = 500)  {
  
  meta <- metadata(obj)
  header <- meta$header
  meta$header <- NULL
  meta$param <- ScanBamParam(flag=meta$param@flag, what=("qname"))
  meta$which <- NULL
  meta$asMates <- NULL
 
  gal <- do.call(readGAlignments, meta)
  metadata(gal) <- list(header=header)

  obj <- as(gal, "GRanges")
  mcols(obj) <- NULL
  obj <- promoters(obj, upstream = 0, downstream = 1)

  cvg <- coverage(obj)
  avg_cvg <- lapply(cvg, function(x) sum(as.numeric(x)))
  cvg <-cvg[names(cvg) %in% seqlev]
  
  txs <- txs[seqnames(txs) %in% seqlev]
  txs <- unique(txs)
    
  sel.gr <- promoters(txs, upstream = upstream, downstream = downstream)
  sel.gr <- split(sel.gr, seqnames(sel.gr))
  
  seqlev <- seqlev[seqlev %in% names(sel.gr)]
  cvg <- cvg[seqlev]
  sel.gr <- sel.gr[seqlev]
  
  vms <- unlist(viewMeans(Views(cvg, sel.gr)))
  sel.gr <- unlist(sel.gr)
  mcols(sel.gr) <- vms
 
  smallNumber <- max(c(1e-06, min(unlist(vms))))
  vms.log2 <-log2(vms + smallNumber)
  return(list("vms"=vms.log2, "sex"= avg_cvg))
}