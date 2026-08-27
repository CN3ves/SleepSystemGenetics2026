# Scriptcontaining the helper functions for annotation plots

# Adapted from ChIPseeker
plot_covfc <- function(peak.gr, weightCol = "logFC", xlab = "Chromosome Size (bp)", 
    ylab = "", title = "Significant peaks location", labels=c("Open","Closed","Chr limit"), legend= "SD vs CTRL") {
  
  # Get signal
  weight <- mcols(peak.gr)[[weightCol]]
  peak.cov <- coverage(peak.gr, weight=weight)
  
  pcov <- lapply(peak.cov, IRanges::slice, lower=0.1)
  ncov <- lapply(peak.cov, IRanges::slice, upper=-0.1)
  cov <- lapply(1:length(peak.cov), function(i) {
    c(pcov[[i]],ncov[[i]])
  })
  names(cov) <- names(pcov)
  
  # Make data frame
  ldf <- lapply(1:length(cov), function(i) {
    
    x <- cov[[i]]
    l <- length(BSgenome.Mmusculus.UCSC.mm10[[names(cov)[i]]])
    
    if (length(x@ranges) == 0) {
      msg <- paste0(names(cov[i]),
                    " doesn't contain significant signal")
      message(msg)
      return(NA)
    }
    
    rval <- runValue(x)
    rval <- sapply(rval@listData, mean)
    
    data.frame(chr   = names(cov[i]),
               start = start(x),
               end   = end(x),
               cnt   = rval,
               len = l
              )
  })
  
  ldf <- ldf[!is.na(ldf)]
  df <- do.call("rbind", ldf)
  
  #Sort chrs
  chr.name <- as.character(unique(df$chr))
  noChr <- suppressWarnings(as.numeric(sub("chr", "", chr.name)))
  ch.idx <- which(is.na(noChr))
  n.idx <- which(!is.na(noChr))
  chr.n <- noChr[n.idx]
  
  chr.sorted <- chr.name[n.idx][order(chr.n)]
  if (length(ch.idx) != 0) {
    chr.sorted <- c(chr.sorted, sort(chr.name[ch.idx]))
  }
  
  df$chr <- factor(df$chr, levels=chr.sorted)

  # Group data for plotting
  tm <- group_by(df, chr, start, end, len) %>% summarise(value=sum(cnt))
  
  
  # Add colors
  colors <- rep("coral", nrow(tm))
  colors[tm$value > 0] <- "cornflowerblue"
  tm$color <- colors
    
  # get chr boundaries
  limits <- unique(tm[,c("chr","len")])
  limits$start <- limits$len # right limit up tick
  limits$end <- limits$len
  limits$color <- "black"
  limits$value <- round(max(abs(tm$value)))
  limits2 <- limits# right limit down tick
  limits2$value <- -round(max(abs(tm$value)))
  limits <- rbind(limits,limits2)
  limits <- limits[,names(tm)] #sort limits to fit plot data
  limits2 <- limits
  limits2$start <- limits2$end <- 0 # length limits
  limits <- rbind(limits,limits2)
  tm <- rbind(tm,limits)
  
  tm$color <- factor(tm$color, levels=c("cornflowerblue",  "coral","black"))
  
  print(paste(table(tm$color)))
  p <- ggplot(tm, aes(start, value, color=color, fill=color))
  p <- p + geom_rect(aes(xmin = start, ymin = 0, xmax = end,  ymax = value))
  p <- p + facet_grid(chr ~ ., scales = "free")
  p <- p + theme_classic()
  p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title)
  p <- p + scale_y_continuous(breaks=seq(-0.5,0.5,0.5))
  p <- p + theme(strip.text.y = element_text(angle = 360),plot.title = element_text(hjust = 0.5))
  p <- p  + scale_color_manual(values=c("cornflowerblue", "coral","black"), labels=labels)
  p <- p  + scale_fill_manual(values=c("cornflowerblue", "coral","black"), labels=labels)
  p <- p + guides(color=guide_legend(title=legend), fill=guide_legend(title=legend))
  p <- p  + geom_hline(yintercept=0, linetype="dashed", color = "grey")
 
  return(p)
}


qtl_plot <- function(peaks, snps) {
  snps <- snps[gsub("\\.[0-9]*","",snps$lodcolumn) %in% peaks$PeakID,]
  snps <- snps[, c("chr","bp", "bp", "lod")]
  names(snps) <- c("seqnames","start", "end", "logFC")
  snps$seqnames <- paste0("chr",snps$seqnames) 
  snps <- GRanges(snps)
   
  merged <- list()
  ov <- findOverlaps(snps,snps)
  for (i in which(!duplicated(snps))) {
 
    subs <- subjectHits(ov)[queryHits(ov)==i] # get replicates
    
    temp <- reduce(snps[subs]) 
    temp$logFC <-  -log2(mean(2^snps[subs]$logFC)) # get log of mean LODs
    
    merged[[i]] <- temp
    
  }
  
  merged <- do.call(c,merged)
  
  return(c(peaks, merged))
}

tag_plots <- function(window, regions, file) {
  promoter <- getPromoters(TxDb=txdb, upstream=window, downstream=window)

  
  SD <- lapply(regions, getTagMatrix, weightCol="SD", windows=promoter)
  names(SD) <- paste0( names(SD), "_SD")

  CTRL <- lapply(regions, getTagMatrix, weightCol="CTRL", windows=promoter)
  names(CTRL) <- paste0( names(CTRL), "_CTRL")
  
  tagMatrixList <- c(SD, CTRL)
  tagMatrixList <- tagMatrixList[sort(names(tagMatrixList))]

  for(test in names(regions)) {
    cat(paste("Tag:",test,"\n"))
 
    png(paste0(file,"/tagheatmap_", test , "_", window ,".png"))
    tagHeatmap(tagMatrixList[grep(paste0("^",test,"_"),names(tagMatrixList))])
    dev.off()
    
  }
 
}
