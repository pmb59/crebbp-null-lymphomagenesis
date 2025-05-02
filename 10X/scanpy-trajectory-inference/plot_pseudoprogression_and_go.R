#--------------------------------------------
# Clustering of Pseudotime Trajectories
#--------------------------------------------
setwd("./PAGA")

library(Seurat) 
library(cowplot)
library(dplyr)
set.seed(1)

#  NO contaminants c(0:9,12)
FOLDER ='cluster8'
CLUSTER = 8
ALLOW_CONDS = c("WT" ,  "CRE" , "PRE"  , "MLEU") 
TOP_n_GENES = 1000    # most variable genes
optimal_K = 10
run_select_optimal_K = FALSE
N = 20     # Enriched GOs to plot

# read Seurat Data
Bcells.combined <- readRDS(file ='seurat3_all.rds')   # no TLEU, MT 5%, res=0.4

temp  <- SubsetData(object=Bcells.combined, subset.name = "batch", accept.value = ALLOW_CONDS )  #  NO contaminants c(0:9,12)

WT_Bcells  <- SubsetData(object=temp, subset.name = "seurat_clusters", accept.value = CLUSTER )  #  NO contaminants c(0:9,12)

rm(Bcells.combined, temp)



setwd(paste0( "./PAGA/", FOLDER)  )
######################################################################
# Pseutotime (0-1) obtained with Scanpy
# Norm gene expresion from Seurat 3
######################################################################
ptime <- read.csv( paste0("pseudotime_Cluster_" ,CLUSTER, ".csv" ) , header =TRUE )
colnames(ptime)[1] <- 'cell_name'
x <- t(as.matrix(GetAssayData(object = WT_Bcells, assay = "RNA", slot = "data")) )    # 'data' or 'scaled.data' or 'counts'
x <- as.data.frame( x )


# find most variable genes IN THIS SET OF CELLS
library(transcripTools)
# MostVar: samples in columns, genes in rows
x <- as.data.frame( t( mostVar(data = as.matrix(t(x) ) , n = TOP_n_GENES ) ) )

Ly <- matrix(NA, nrow=nrow(x), ncol= TOP_n_GENES )
ptime$dpt_pseudotime <- ptime$dpt_pseudotime +  rnorm(length(ptime$dpt_pseudotime), mean = 0, sd = 0.0001)  # FClust does not access duplcayed t-values
sortPseudo <- sort( ptime$dpt_pseudotime, decreasing=FALSE, index.return=TRUE)[[2]]

for (j in 1:ncol(x) ) {
  Ly[,j] <- x[sortPseudo,j]
}


if (run_select_optimal_K == TRUE) {

  tw <- c()
  for(d in 1:20){
    print(d)
    tw [d] <- kmeans(t(Ly), centers=d)$tot.withinss
  }
  plot(1:20, tw, type='o')
  library(factoextra)
  library(cluster) 
  library(ggplot2)
  
  # Silhouette method
  fviz_nbclust(t(Ly), kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette method")

  # Elbow method
  fviz_nbclust(t(Ly), kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2) +
    labs(subtitle = "Elbow method")

  # Gap statistic
  gap_stat <- clusGap(t(Ly), FUN = kmeans, nstart = 25, K.max = 10, B = 50)
  print(gap_stat, method = "firstmax")
  fviz_gap_stat(gap_stat)

}

library(fda.usc)
tt <- ptime$dpt_pseudotime [sortPseudo]
tt [ which(tt<0) ] <- 0
fdataobj<- fdata(t(Ly),tt)

# Perform k-means clustering on functional data
km <- kmeans.fd(fdataobj, ncl = optimal_K , draw = FALSE, cluster.size =2 )
cl <- km $cluster

temp1 <- data.frame(colnames(x), cl )
temp2<- temp1[sort(temp1$cl, decreasing=FALSE, index.return=TRUE)$ix, ]
write.csv(temp2, file=paste0("Top_",TOP_n_GENES,"_most_variable_genes_Cluster",CLUSTER,'.csv') )
rm(temp1, temp2)


gene_expression <- c()
gene <- c()
pseudotime <- c()
cluster <- c()

for (i in 1:optimal_K ){
  
  print(i)
  ind <- which( as.numeric(cl) == i  )
  
  for (j in 1:length(ind)){
    gene_expression <- c( gene_expression ,  x[ ,ind[j] ]  )
    gene <- c(gene, rep( colnames(x) [ind[j] ]  , nrow(x)  )  )
    pseudotime <- c( pseudotime ,  ptime$dpt_pseudotime  )
    cluster <- c(cluster,   rep(i,nrow(x)  ) )
  }
  rm(ind)
}

df <- data.frame(gene_expression=gene_expression, gene=gene, cluster=cluster, pseudotime=pseudotime)
df$cluster <- paste0('pattern ' , df$cluster) 

library(ggplot2)

p <- ggplot( df, aes(pseudotime, gene_expression) ) +
  geom_smooth( se = TRUE, span = 1, fill='magenta') + 
  scale_x_continuous(breaks = c(0,0.5,1) ) +
  facet_wrap(~ cluster, nrow = 2 , scales='free')+ ylab('gene expression') +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90) ) + ggtitle( paste0('Cluster',CLUSTER) ) + theme_classic()

p
ggsave(paste0("Patterns_TopGenes_",TOP_n_GENES   ,"_cluster_",CLUSTER,'.pdf')  , height=5, width= 8.5 )

pdf(paste0("Barplot_TopGenes_",TOP_n_GENES   ,"_cluster_",CLUSTER,'.pdf')  , height=6, width= 5.5 )
barplot(table(km$cluster) , xlab='pattern' , ylab='number of genes', fill='blue', col='blue')
dev.off()




library(circlize)
library(ComplexHeatmap)

ha = rowAnnotation(df = as.data.frame( as.character(cl) ), name='cluster', show_annotation_name = FALSE)

Condition = factor( ptime$batch[sortPseudo] )


Pseudotime = ptime$dpt_pseudotime[sortPseudo]
col_fun = colorRamp2(c(0, 0.5, 1), c("white", "gray", "black")) 


hb = HeatmapAnnotation(Condition= Condition , Pseudotime = Pseudotime ,
                       show_annotation_name = TRUE, which = "column", col=list(Pseudotime = col_fun, Condition = c("WT"="lightcoral","CRE"="chartreuse3","PRE"="mediumturquoise","MLEU"="purple") ))


pdf(paste0("Heatmap_TopGenes_",TOP_n_GENES   ,"_cluster_",CLUSTER,'.pdf')  , height=6, width= 7 )
Heatmap(t(Ly), name = "Gene expression", split = cl, 
        cluster_columns=FALSE, cluster_rows=TRUE,
        column_dend_reorder=FALSE , gap = unit(3, "mm"), column_title = "Cell pseudotime ordering" , 
        row_title = "Genes",
        row_title_rot = 0,
        show_row_dend = FALSE , top_annotation = hb, use_raster = TRUE )   #+ ha

dev.off()



############################################################
# GO analysis
############################################################

go <- data.frame(colnames(x), cl )

# Replace by original
go <- read.csv(".../Top_1000_most_variable_genes_Cluster8.csv")[,2:3]


library(limma)
library(org.Mm.eg.db)
library(GO.db)
library(GOstats)
library(biomaRt)
library(ggplot2)

library(biomaRt)
ensembl84 <- useMart(host='http://may2024.archive.ensembl.org',
                     biomart='ENSEMBL_MART_ENSEMBL',
                     dataset='mmusculus_gene_ensembl')

# Find matches of "entry" in the list
A <- sort(listAttributes(ensembl84)$name)
matches <- grep("entrez", A, value = TRUE, ignore.case = TRUE)

for (i in 1:optimal_K){
  print(i)
  ## First, get Entrez Gene identifiers  fot our list of interest
  ens_id_interest <- go[ which(go$cl == i), 1]
  
  my_entrez_gene <- getBM(attributes='entrezgene_id',
                          filters = 'external_gene_name',
                          values = ens_id_interest,
                          mart = ensembl84)
  
  universe_gene <- getBM(attributes='entrezgene_id', mart = ensembl84)
  
  GOanalysis <- goana(de=as.character(my_entrez_gene$entrezgene_id), universe = NULL, species = "Mm", prior.prob = NULL, covariate=NULL, plot=FALSE)
  
  GOanalysis$P.DE <- p.adjust(GOanalysis$P.DE, "BH")
  GOordered <- GOanalysis[order(GOanalysis$P.DE, decreasing = FALSE),]

  ##############################
  #now split between ontologies
  ##############################
  #BP
  GOordered1 <- GOordered[ which(GOordered$Ont=='BP'),  ]
  GOordered1 <- transform(GOordered1,  Term= reorder(Term, -P.DE))
  
  p1 <- ggplot(data=GOordered1[1:N,] , aes(x=Term, y= -log10(P.DE)  ) )  +  # -log10(P.DE )
    geom_bar(stat="identity", fill="royalblue1")+ theme_linedraw() + ggtitle('GO (BP)') + ylab('-log10(adj p-value)') +xlab('') + coord_flip() + theme(axis.text.x = element_text(angle = 90))
  
  #CC
  GOordered2 <- GOordered[ which(GOordered$Ont=='CC'),  ]
  GOordered2 <- transform(GOordered2,  Term= reorder(Term, -P.DE))
  
  p2 <- ggplot(data=GOordered2[1:N,] , aes(x=Term, y= -log10(P.DE)  ) )  +
    geom_bar(stat="identity", fill="palevioletred2")+ theme_linedraw() + ggtitle('GO (CC)') + ylab('-log10(adj p-value)') +xlab('')+ coord_flip()+ theme(axis.text.x = element_text(angle = 90))
  
  #MF
  GOordered3 <- GOordered[ which(GOordered$Ont=='MF'),  ]
  GOordered3 <- transform(GOordered3,  Term= reorder(Term, -P.DE))
  
  p3 <- ggplot(data=GOordered3[1:N,] , aes(x=Term, y= -log10(P.DE)  ) )  +
    geom_bar(stat="identity", fill="seagreen3")+ theme_linedraw() + ggtitle('GO (MF)') + ylab('-log10(adj p-value)') + xlab('') + coord_flip()+ theme(axis.text.x = element_text(angle = 90))
  
  #############################
  #KEGG
  
  KEGG <- kegga(de=as.character(my_entrez_gene$entrezgene_id), universe = NULL, #as.character(universe_gene$entrezgene_id) 
                species = "Mm", species.KEGG = 'mmu', convert = FALSE,
                gene.pathway = NULL, pathway.names = getKEGGPathwayNames(species.KEGG='mmu', remove.qualifier =TRUE) ,
                prior.prob = NULL, covariate=NULL, plot=FALSE)
  KEGG$PathwayID <- sub("path:", "", rownames(KEGG)) 
  
  PATHWAYS <- getKEGGPathwayNames(species.KEGG='mmu', remove.qualifier =TRUE)

  KEGG2 <- merge(KEGG, PATHWAYS, by='PathwayID' )
  
  KEGG2$P.DE <- p.adjust(KEGG2$P.DE, "BH")
  KEGGOrdered <- KEGG2[order(KEGG2$P.DE, decreasing = FALSE),]
  KEGGOrdered$Pathway <- as.character(KEGGOrdered$Description)
  
  write.table(x=KEGGOrdered, file = paste('KEGG', 'Top', N, 'Cluster',CLUSTER, 'Pattern', i, 'csv', sep='.'),
              append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
  KEGGOrdered <- transform(KEGGOrdered,   Pathway = reorder(Pathway, -P.DE))
  

  p4 <- ggplot(data=KEGGOrdered[1:N,] , aes(x=Pathway, y= -log10(P.DE)  ) )  +
    geom_bar(stat="identity", fill="yellow3")+ theme_linedraw() + ggtitle('KEEG pathway') +xlab('') + ylab('-log10(adj p-value)') + coord_flip()+ theme(axis.text.x = element_text(angle = 90))
  
  source("multiplot.R")
  pdf(paste('new_GO', 'Top', N, 'Cluster',CLUSTER, 'Pattern', i, 'pdf', sep='.'), width=11, height=11 )
  multiplot(p1, p2, p3, p4, cols=2)
  dev.off()
  
  rm(ens_id_interest, universe_gene, KEGG2, my_entrez_gene, GOanalysis, GOordered, GOordered1, GOordered2, GOordered3, KEGG , KEGGOrdered, p1, p2, p3, p4 )

}

