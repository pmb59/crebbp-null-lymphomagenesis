setwd("/working-dir/")

library(ggplot2)

# Multiple plot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


df0 <- read.table(file="quartiles.txt", header = FALSE, sep = "\t")[,c(1:7)]
# row titles
colnames(df0) <- c('run', 'decile' ,'pathway', 'median_expr_pathway', 'mean_expr_pathway', 'median_expr_crebbp',  'mean_expr_crebbp')

df0$run <- as.character(df0$run)
df0$pathway <- as.character(df0$pathway)
U <- sort(unique(df0$pathway))

pdf("plot_corr_with_Smoothing_mean_line.pdf", width = 11, height = 5.5)

for (k in U){
  
  df <- df0[which(df0$pathway == k),]

  # calculate correlation
  df$corrMedian <- NA
  df$corrMean <- NA
  df$corrMedianp.value <- NA
  df$corrMeanp.value <- NA
  for (i in unique(df$pathway) ){
    ind <- which(df$pathway == i)
    CORREL_mean<- cor.test(df$mean_expr_pathway[ind],df$mean_expr_crebbp[ind], method='pearson') 
    df$corrMean[ind] <- CORREL_mean$estimate
    df$corrMeanp.value[ind] <- CORREL_mean$p.value
    CORREL_median <- cor.test(df$median_expr_pathway[ind],df$median_expr_crebbp[ind], method='pearson') 
    df$corrMedian[ind] <- CORREL_median $estimate
    df$corrMedianp.value[ind] <- CORREL_median $p.value

    data_to_write <- data.frame(
      Pathway = i, 
      CORREL_mean_estimate = as.numeric(CORREL_mean$estimate), 
      CORREL_mean_pvalue = CORREL_mean$p.value, 
      CORREL_median_estimate = as.numeric(CORREL_median$estimate), 
      CORREL_median_pvalue = CORREL_median$p.value
    )
  
    write.table(data_to_write, file = "KEGG_CORRELATION_3.txt", append = TRUE, quote = FALSE, 
              sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
  
  }

  p1 <- ggplot(df, aes(mean_expr_crebbp, mean_expr_pathway, colour = pathway, group=pathway)) +
    geom_point(size = 1, shape = 21, colour = "black", fill = 'red', stroke = 0.5) + 
    geom_smooth(method = "loess", se = TRUE, colour = "#B22222", fill = "#FF9999") +  # Add smooth line and confidence interval
    xlab('CREBBP expression (TMM)') +
    ylab('Mean expression of pathway genes (TMM)') +
    facet_wrap(~ pathway) + 
    ggtitle(paste("Correlation = ", round(df$corrMean,2), ", p = ", round(df$corrMeanp.value,4))) +  
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) +
    theme_cowplot()


  p2 <- ggplot(df, aes(median_expr_crebbp, median_expr_pathway, colour = pathway, group=pathway)) +
    geom_point(size = 1, shape = 21, colour = "black", fill = 'blue', stroke = 0.5) + 
    geom_smooth(method = "loess", se = TRUE, colour = "#000080", fill = "#87CEFA") +  # Add smooth line and confidence interval
    xlab('CREBBP expression (TMM)') +
    ylab('Median expression of pathway genes (TMM)') +
    facet_wrap(~ pathway) + 
    ggtitle(paste("Correlation = ", round(df$corrMedian,2), ", p = ", round(df$corrMedianp.value,4))) +  
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) +
    theme_cowplot()


  multiplot(p1, p2, cols=2)
       
}
dev.off()

