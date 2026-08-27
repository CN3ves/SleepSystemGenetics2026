# Functions for plotting

coef_box <- function(coeff, title, file) {
  
  plot <- coeff %>% 
    pivot_longer(1:(ncol(coeff)-1), names_to = "perc", values_to ="val")
  
  png(file)
  
  p <- ggplot(plot, aes(x=perc, y=val, color=sig)) +
    geom_boxplot() + 
    labs(title=title) + ylab("Regression coefficients") +  xlab("") + 
    theme_classic() + 
    guides(color = guide_legend(title="Significant region")) +
    theme( plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  
  print(p)
  dev.off()

}

coef_heatmap <- function(coeff, filter, title, file) {
  coeff <- as.matrix(coeff[filter, -which(names(coeff) %in% c("SD","sig"))])
  
  png(file)
  p <- Heatmap(coeff, border_gp = gpar(col = "black", lty = 1),
          column_title = title,
          heatmap_legend_param = list(title = "Coeff vs CTRL"),
          show_row_names = FALSE)
  
  print(p)
  dev.off()
}
