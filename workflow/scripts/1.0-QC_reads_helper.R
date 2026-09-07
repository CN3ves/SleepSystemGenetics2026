# Script containing the helper functions for read QC summary

# Function to extract QC read summary
qc_summary <- function(report) {
  
  # Get filtering summary 
  summary <- c(report[["summary"]], report[["filtering_result"]])
  res <- unlist(summary)
  
  # Standardize the names to "value (before/after filtering)" 
  names(res) <- gsub("(.*)_filtering.(.*)", "\\2 (\\1 filtering)", names(res))
  
  # Remove GC content (replicated later)
  res <- res[-grep("gc_content", names(res))]
  
  # Format numbers to strings
  res<- format(res, trim = TRUE, scientific = FALSE)
  
  
  return(res)
}

# Function to extract QC duplication summary
qc_duplication <- function(report) {
  
  # Get duplication summary 
  duplications <- report[["duplication" ]]
  
  # Standardize the names
  names(duplications) <- c("Duplication rate")
  
  # Format numbers to strings
  duplications <- format(duplications, trim = TRUE, scientific = FALSE)

  return(duplications)
}

# Function to extract QC adapter summary
qc_adaptor <- function(report) {
  
  # Get adapter summary 
  res <- unlist(report[["adapter_cutting"]])
  
  # Standardize the names if there is information on the report
  if (! is.null(res)) {
    names(res) <- gsub("read1_adapter_counts.(.*)","read1_adapter_counts (\\1)", names(res))
    return(res)
  }
  
  return(NULL)
}

#Function to extract QC quality summary
qc_quality <- function(values) {# called by qc_filter
  # Check bases with lower than 30 phred score
  tryCatch(
    expr = {
      
      data.frame("bases"= which(values<30)) %>% 
        group_by(grp = cumsum(bases-lag(bases,default = 0) > 1)) %>% # Group consecutive bases
        summarise(bases = paste(range(bases),  collapse = "-")) %>% # Make string
        dplyr::select(-grp) %>%
        paste(sep="; ")
      
    }, warning = function(w) { # No bases found then result is -
        "-" 
      }
  )
}

# Function to extract QC read filtering summary
qc_filter <- function(report, stage) {
  
  # Define values name (after/before) 
  col <- paste0("read1_", stage, "_filtering")
  
  # Quality QC
  qual <- sapply(report[[col]][["quality_curves"]], qc_quality)
  names(qual) <- paste0("Cycles with quality < 30 ", time," filtering (", names(qual),")")
  
  # Base content QC
  content <- sapply(report[[col]][["content_curves"]], function(c) mean(c)*100)
  names(content) <- paste0("Average ", names(content), " content % ", time, " filtering")
  
  # Format numbers to strings
  content <- format(content, trim = TRUE, scientific = FALSE)

  res <- c(qual, content)
    
  return(res)
}

# Function to summarize results
summarize <- function(i, tab) {
  
  # Report value (row) to summarize - 1st column is a fake summary placeholder
  row <- trimws(tab[i,-1])
  
  # If the value is numeric
  if (! any(grepl("[^0-9\\.]", row[!is.na(row)])) & ! any(grepl("\\..*\\.", row)) & ! grepl("adapter_counts", rownames(tab)[i])) {
    row <- as.numeric(row)
    
    r <- round(range(row, na.rm = TRUE),2) # Range
    m <- round(mean(row, na.rm = TRUE),2) # Mean
    na <- round(mean(is.na(row))*100,0) # Missing 
    
    value <- paste0("Mean: ", m, "; Range: ", paste(r, collapse = "-"), "; Missing: ", na, "%")
 
  # If the value is range of cycles 
  } else if (grepl("Cycles with quality", rownames(tab)[i])) {
    # Remove string characters  and extract all gaps
    c <- gsub("[c(\" )]", "" , row) 
    c <- sapply(c, function(x) {strsplit(x, ",")[[1]]})
    
    # Get unique list of gaps with low quality
    c <- unlist(c)
    c <- unique(c)
    c <- c[grep("[0-9]", c)]
    
    # If there are any low quality gaps, get "largest range"
    if (length(c)>0) {
      l <- min(as.numeric(sapply(c, function(x) {strsplit(x, "-")[[1]][1]})), na.rm=TRUE)
      r <- max(as.numeric(sapply(c, function(x) {strsplit(x, "-")[[1]][2]})), na.rm=TRUE)
    
      value <- paste("Low-quality gap: ", l,"-", r)
    } else {
      value <- paste("No low-quality bases")
    }
  
  # If the value is categorical/string
  } else {
    row <- as.character(row)
    
    u <- length(unique(row[!is.na(row)])) # Number unique categories 
    na <- round(mean(is.na(row))*100,0) # Missing 
    
    value <- paste0("Number of unique value: ", u,"; Missing: ", na, "%")
  }
  
  return(value)
} 
