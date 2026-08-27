# Script containing the helper functions for enrichment

get_lists_from_genes <- function(genes, universe) {
  genes$ID <- AnnotationDbi::mapIds(org.Mm.eg.db, 
       keys = gsub("\\.[0-9]*","",rownames(genes)),
       column = "ENTREZID",
       keytype = "ENSEMBL")

  genes <- genes[!is.na(genes$ID),] # there's some creepy nesting of a df in a df!
  genes <- genes[!duplicated(genes$ID),]
  
  FCList = genes$logFC
  names(FCList) = genes$ID
  FCList = sort(FCList, decreasing = TRUE)

  up <- genes[genes$logFC > 0 & genes$FDR < 0.05, "ID"]
  down <- genes[genes$logFC < 0 & genes$FDR < 0.05, "ID"]
 
  genes <- genes[genes$FDR < 0.05, "ID"]

  return(list(genes = genes, up = up, down=down, universe =  universe, fc = FCList))  
}

run_enrich <- function(lsts, params = list(ora = c("genes", "up", "down"), gsea = FALSE), type, dir, name=NULL) {
  # ORA ANALYSES
  

  if(type=="GO") {
    ontology=c("all", "BP", "MF","CC")
  } else {
    ontology=''
  }
  
  for (test in params[['ora']]) {
    print(paste("ORA analyses for", type))
    for(ont in ontology) {
      print(paste(test, ont))
      enrich <- enrichment(genes = lsts[[test]], universe = lsts[['universe']], type =type, ontology=ont)
      print(enrich)
      write.csv(enrich, paste0(dir,"/",name,"_ORA_", type, "_", test,"_",ont, "_enrichment.csv"))
    
      # Plots
      if (!is.null(enrich) & (sum(enrich@result$qvalue < 0.05) > 5 & !(type=="GO" & ont=='all'))) {
        print(paste('plotting',name,"ORA", type, test,ont))
        plot_enrich(enrich, lsts[['fc']], paste0(dir,"/",name,"_ORA_", type, "_", test,"_",ont), type)
        
        if(type == 'KEGG') {
          print("Top pathways:")
          print(sapply(1:5, function(i) browseKEGG(enrich, enrich$ID[i])))
        }
      }
    }
  }
  
   # GSEA ANALYSES
  if (params[['gsea']]) {
    print(paste("GSEA analyses for", type))
     for(ont in ontology) {

      enrich <- enrichment(genes = lsts[['fc']], universe = NULL, type =type, ontology=ont)
      print(enrich)
      write.csv(enrich, paste0(dir,"/",type,"_GSEA_", name, "_enrichment.csv"))
    
      # Plots
      if (!is.null(enrich) & (sum(enrich@result$qvalue < 0.05) > 5 & !(type=="GO" & ont=='all'))) {
        print(paste('plotting',name,"GSEA", type))
        plot_enrich(enrich, lsts[['fc']], paste0(dir,"/",name,"_GSEA_", type,"_",ont), type)
        
        if(type == 'KEGG') {
          print("Top pathways:")
          print(sapply(1:5, function(i) browseKEGG(enrich, enrich$ID[i])))
        }
      }
     }
  }
}

enrichment <- function(genes, universe=NULL, type="GO", ontology=NULL) {
  tryCatch({
    if (type == "GO") {
    if (is.null(universe)) {
      res <- gseGO(geneList= genes, 
                   OrgDb = org.Mm.eg.db, ont = ontology, eps=0,
                   pvalueCutoff  = 0.05, verbose      = FALSE)
  
    } else {
      res <- enrichGO(gene = genes, universe = universe,
                      OrgDb = org.Mm.eg.db, ont = ontology,
                      pAdjustMethod = "BH", pvalueCutoff  = 0.05, 
                      qvalueCutoff  = 0.05,readable = TRUE)
  
    }
  }
  
  if (type == "KEGG") {
    if (is.null(universe)) {
      res <- gseKEGG(geneList = genes, organism = 'mmu',
              eps=0, pvalueCutoff = 0.05, verbose = FALSE)

    } else {
      res <- enrichKEGG(gene= genes, universe = universe,
                 organism = 'mmu', pvalueCutoff = 0.05)
  
    }
  }
  
  if (type == "Wikipaths") {
     if (is.null(universe)) {
      res <- gseWP(genes, organism = "Mus musculus")
    } else {
      res <- enrichWP(genes, universe=universe, organism = "Mus musculus") 
    }
  }
  
  if (type == "Reactome") {
    if (is.null(universe)) {
      res <- gsePathway(genes,  pvalueCutoff = 0.2, 
                        pAdjustMethod = "BH", organism = "mouse", verbose = FALSE)

    } else {
      res <- enrichPathway(gene=genes, universe=universe, 
                           organism = "mouse", pvalueCutoff = 0.05, readable=TRUE)
    }
  }},error =function(cond) {
    print("Error skipped")
    print(cond)
    res=NULL
  })
    
  return(res)
}

plot_enrich <- function(enrich, foldchange,test,type="GO") {
  
  analysis <- class(enrich)
  ont=NULL
  
  if (analysis == "enrichResult") {
    ont <- enrich@ontology
    # Ontology plot
    if (type == "GO") {
      print("Ontology plot")
      p <- goplot(enrich, showCategory=25, geom="label") + ggtitle(paste(type, "ontology", ont)) +
        theme(plot.title = element_text(hjust = 0.5))
      png(paste0(test,"_ontology.png"))
      print(p)
      dev.off()
    }
    
    
    #barplot
    print("Barplot")
    p <- barplot(enrich, showCategory=25, font.size = 8, title = paste(type, "enrichment barplot", ont)) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Gene Ratio")
    png(paste0(test,"_barplot.png"))
    print(p)
    dev.off()
  
  }
  if (analysis == "gseaResult") {
    ont <- enrich@setType
    
    # Ridge plot 
    print("Ridgeplot")
    p <- ridgeplot(enrich, showCategory = 25) + theme(axis.text.y = element_text(size = 7), plot.title = element_text(hjust = 0.5)) + ggtitle(paste(type, "ridgeplot", ont)) + xlab("Log counts")
    png(paste0(test,"_ridge.png")) 
    print(p)
    dev.off()
    
    # GSEA
    print("GSEA")
    p <- gseaplot2(enrich,  geneSetID = 1:5, title = paste(type, "GSEA plot", ont)) + theme(plot.title = element_text(hjust = 0.5))
    png(paste0(test,"_geneset.png"))
    print(p)
    dev.off()
    
  }
  
  # Dotplot
  print("Dotplot")
  p <-  dotplot(enrich,color = "p.adjust", showCategory = 25,font.size = 8, title = paste(type, "enchiment dotplot", ont) ,label_format = 3) + guides(size = 'none')+ theme(plot.title = element_text(hjust = 0.5))
  png(paste0(test,"_dotplot.png"))
  print(p)
  dev.off()
  
  #concept net plot
  print("Cnet")
  p <- cnetplot(enrich, foldChange=foldchange, node_label="category", showCategory=10)+ guides(size = 'none') + ggtitle(paste("GO Gene-Concept Network", ont)) + theme(plot.title = element_text(hjust = 0.5))
  png(paste0(test,"_cnetplot.png"))
  print(p)
  dev.off()
  
  # Heatplot
  print("Heatplot")
  p <- heatplot(enrich, foldChange=foldchange, showCategory=25)+ ggtitle(paste(type,"heatplot", ont)) + theme(axis.text.x = element_text(size = 1), plot.title = element_text(hjust = 0.5)) + xlab("Genes")
  tryCatch({
    png(paste0(test,"_heatplot.png"))
    print(p)
    dev.off()
  },error =function(cond) {
    print("Heatplot has problems")
    print(cond)
    res=NULL
  })
  
  # enrichemnt map plot
  print("Emap")
  p <- emapplot(pairwise_termsim(enrich))
  png(paste0(test,"_emap.png"))
  print(p)
  dev.off() 
  
  # upset plot
  print("Upset plot")
  p <- upsetplot(enrich)
  tryCatch({
    png(paste0(test,"_upset.png"))
    print(p)
    dev.off() 
  },error =function(cond) {
    print("Upset has problems")
    print(cond)
    res=NULL
  })
  return(NULL)
}

browseKEGG <- function (x, pathID) {
    url <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?", 
                  pathID, "/", x[pathID, "geneID"])
    return(url)
}
  