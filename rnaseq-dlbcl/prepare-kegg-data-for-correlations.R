
DATA_PATH <- "./working-dir/path-to-rna_seq.star_gene_counts/gdc_download_20200311_174748.005232"
useGenesOfInterest <- TRUE    # if FALSE uses 200 most variable genes
ReadRawFromRobject <- TRUE

# h=23 "Oxidative phosphorylation - Homo sapiens (human)"
# h=196 B cell receptor signaling pathway

library("KEGGREST")
human_pathways <- names(keggList("pathway", "hsa"))
head(human_pathways)
length(human_pathways)

for (h in 1:length(human_pathways)) {
  print(h)
  pathway_id <- human_pathways[h]
  pathway_info <- keggGet(pathway_id)

  if ("GENE" %in% names(pathway_info[[1]])) {
    print("GENE exists in the list")
    genes <- pathway_info[[1]]$GENE
    namePathway <- pathway_info[[1]]$NAME
    genes_list <- genes[seq(2, length(genes), 2)]
    gene_names <- unique(gsub(";.*", "", genes_list))
    cleaned_vector <- gene_names[!grepl(" ", gene_names)]
    gene_names <- cleaned_vector

    print(gene_names)

    GENE_SETS <- paste0('./path-kegg/KEGG/', pathway_id, '.txt')
    
    write.table(x = gene_names, file = GENE_SETS, append = FALSE, quote = FALSE, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = FALSE, qmethod = c("escape", "double"),
                fileEncoding = "")
    
    TITLE <- namePathway
  }
}

FILE_QUARTILES <- 'quartiles.txt'

metada <- read.csv("gdc_sample_sheet.2020-03-11.tsv", sep = "\t")

n <- list.files(path = DATA_PATH,
                pattern = ".tsv.gz", all.files = FALSE,
                full.names = FALSE, recursive = TRUE,
                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

if (ReadRawFromRobject == FALSE) {
  
  for (i in 1:length(n)) {
    print(paste0('reading file..', i, ' of ', length(n)))
    temp <- read.table(paste0(DATA_PATH, "/", n[i]))[, 1:2]
    NEWNAME <- as.character(metada$Case.ID[which(metada$File.Name == strsplit(n[i], split = '/')[[1]][2])])
    colnames(temp) <- c('gene', NEWNAME)
    print(NEWNAME)

    if (i == 1) {
      final <- temp
    }
    if (i > 1) {
      final <- merge(x = final, y = temp, by = 'gene', sort = FALSE)
    }
    rm(temp)
    rm(NEWNAME)
  }
}

if (ReadRawFromRobject == TRUE) {
  final <- readRDS("./working-dir/final.rds")
}

final <- final[5:nrow(final), ]


gene2 <- c()
for (i in 1:length(final$gene)) {
  gene2[i] <- strsplit(x = as.character(final$gene[i]), split = ".", fixed = TRUE)[[1]][1]
}
final$gene <- gene2
rownames(final) <- gene2
final <- final[, 2:ncol(final)]

col_final <- data.frame(sample = as.character(colnames(final)),
                        celltype = 'DLBCL',
                        source = 'NEJM',
                        type = 'tumor',
                        Crebbp = 'NONE',
                        n_Crebbp_mut = 0,
                        ecotyper = 'X',
                        gender = 'X',
                        COO = 'x',
                        CrebbpBin = 'x',
                        Cohort = 'NA',
                        Sample.Type = 'NA',
                        Status.at.last.Follow.up = 'Unknown')

col_final[] <- lapply(col_final, as.character)

coo <- read.csv('/working-dir/mmc2.csv')

for (i in 1:nrow(col_final)) {
  cooi <- which(as.character(coo$Donor.name) == as.character(col_final$sample[i]))
  if (length(cooi) != 1) {
    print(paste0('WARNING - no match in ', i))
  }
  col_final$celltype[i] <- paste0(col_final$celltype[i], "_", as.character(coo$LymphGen.call[cooi]))
  col_final$gender[i] <- as.character(coo$Gender[cooi])
  col_final$COO[i] <- as.character(coo$COO.Class[cooi])
  col_final$Cohort[i] <- as.character(coo$Cohort[cooi])
  col_final$Sample.Type[i] <- as.character(coo$Sample.Type[cooi])
  if (!is.na(coo$Status.at.last.Follow.up[cooi])) {
    col_final$Status.at.last.Follow.up[i] <- as.character(coo$Status.at.last.Follow.up[cooi])
  }
}

setwd("./working-dir")
ecot <- read.csv('ecotyper_cell_state_assignments.tsv', sep = '\t')[, c(1, 3)]

for (i in 1:nrow(col_final)) {
  cooi <- which(as.character(ecot$ID) == as.character(col_final$sample[i]))
  if (length(cooi) != 1) {
    print(paste0('WARNING - no match in ', i))
  }
  col_final$ecotyper[i] <- as.character(ecot$B.cells_state_assignments[cooi])
}
col_final$ecotyper[which(is.na(col_final$ecotyper) == TRUE)] <- 'unclassified'

mut <- read.csv("./working-dir/phs001444/MAF_NCICCR-DLBCL_phs001444.txt", header = TRUE, sep = '\t')

for (i in 1:nrow(col_final)) {
  cooi <- which(as.character(mut$SUBJECT_NAME) == as.character(col_final$sample[i]))
  occurrences <- which(as.character(mut$GENE.SYMBOL[cooi]) == "CREBBP")
  if (length(occurrences) > 0) {
    col_final$n_Crebbp_mut[i] <- length(occurrences)
    for (k in 1:length(occurrences)) {
      if (k == 1) {
        col_final$Crebbp[i] <- as.character(mut$MUTATION_TYPE[cooi[occurrences[k]]])
      } else {
        col_final$Crebbp[i] <- paste(col_final$Crebbp[i], as.character(mut$MUTATION_TYPE[cooi[occurrences[k]]]), sep = "+")
      }
    }
  }
}

col_final$CrebbpBin <- col_final$Crebbp
col_final$CrebbpBin[which(col_final$Crebbp == 'TRUNC+TRUNC')] <- 'TRUNC'
col_final$CrebbpBin[which(col_final$Crebbp == 'MISSENSE+MISSENSE')] <- 'MISSENSE'
col_final$CrebbpBin[which(col_final$Crebbp == 'MISSENSE+MISSENSE+TRUNC')] <- 'mixed'
col_final$CrebbpBin[which(col_final$Crebbp == 'MISSENSE+TRUNC')] <- 'mixed'
col_final$CrebbpBin[which(col_final$Crebbp == 'MISSENSE+TRUNC+MISSENSE')] <- 'mixed'

col_final$gender[which(col_final$gender == "")] <- 'NA'
table(col_final$gender)

write.csv(x = col_final, file = "metadata_dlbcl.csv", quote = FALSE,
          eol = "\n", na = "NA", row.names = TRUE)

table(colnames(final) == col_final$sample)

orden <- sort(col_final$celltype, index.return = TRUE)$ix
col_final <- col_final[orden, ]
final <- final[, orden]

table(colnames(final) == col_final$sample)

M <- final
col_meta <- col_final
rownames(col_meta) <- col_meta$sample

library(ggplot2)
g <- ggplot(col_meta, aes(celltype, fill = ecotyper))
g + geom_bar() + coord_flip()

M <- as.matrix(M)
table(col_meta$sample == colnames(M))

library("NOISeq")
mydata <- readData(data = M, length = NULL, factors = col_meta)
TMM <- tmm(assayData(mydata)$exprs, long = 1000, lc = 0, k = 0)
table(colnames(TMM) == col_meta$sample)

df <- read.table(GENE_SETS, header = TRUE)
colnames(df) <- 'external_gene_name'
length(unique(df$external_gene_name))

library(biomaRt)
ensembl95 <- useMart(host = 'http://jul2023.archive.ensembl.org',
                     biomart = 'ENSEMBL_MART_ENSEMBL',
                     dataset = 'hsapiens_gene_ensembl')
geneinfo <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'external_gene_name', 'gene_biotype'),
                  mart = ensembl95)

df2 <- merge(x = df, y = geneinfo, by = 'external_gene_name')
df3 <- df2[!(grepl(pattern = "CHR_", x = df2$chromosome_name)), ]
df4 <- df3[!(grepl(pattern = "HSCHR", x = df3$chromosome_name)), ]

genesOfInterest <- unique(df4$ensembl_gene_id)

if (useGenesOfInterest == FALSE) {
  row_variability <- apply(TMM, 1, sd)
  n <- 200
  most_variable_indices <- order(row_variability, decreasing = TRUE)[1:n]
  most_variable_rows <- TMM[most_variable_indices, , drop = FALSE]
  genesOfInterest <- rownames(most_variable_rows)
}

genesOfInterest2 <- c()
for (i in genesOfInterest) {
  match_i <- which(rownames(TMM) == i)
  if (length(match_i) != 0) {
    genesOfInterest2 <- c(genesOfInterest2, i)
  }
}
genesOfInterest <- genesOfInterest2

MM0 <- TMM[genesOfInterest, 1:481]
MM0 <- MM0[which(rowSums(MM0) != 0), ]
MM <- t(scale(t(MM0), center = TRUE, scale = TRUE))

library(cluster)
selected_k <- 5
print(paste0("best K could be... ", selected_k))
clusters <- kmeans(x = t(MM), centers = selected_k)$cluster
table(clusters)

####################################################################################
# CREBBP (ENSG00000005339) expression in all clusters
which(rownames(M) == "ENSG00000005339" )  #44111

# Compute quartiles
quartiles <- quantile(TMM[44111,], probs = seq(0, 1, length.out = length(TMM[44111,])+1 ) )

# Assign labels to quartiles
# labels <- c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10")
labels <- paste0("Q", 1:length(TMM[44111,]) )  

# Cut numeric vector into quartiles and assign labels
quartile_labels <- cut(TMM[44111,], breaks = quartiles, labels = labels, include.lowest = TRUE)

col_final$CREBBP_exp <- as.character(quartile_labels)

# Now I use matrix MM0 to get mean values in Q1, Q2,..., Q10
random_string <- paste(sample(c(0:9, letters, LETTERS), 10, replace = TRUE), collapse = "") # run id

library(gtools)
unique_q <- mixedsort(unique(as.character(quartile_labels)))

N <- nrow(MM0)
for (i in unique_q ){
  print(i)
  median_value_inQuartile <- median(  MM0[, which (quartile_labels == i )]  ) 
  mean_value_inQuartile <- mean(  MM0[, which (quartile_labels == i )]  ) 
  
  median_crebbp <-  median(TMM[44111, which(col_final$CREBBP_exp == i )  ])
  mean_crebbp <- mean(TMM[44111, which(col_final$CREBBP_exp == i )  ])
  
  write.table(data.frame(random_string, i,TITLE, median_value_inQuartile, mean_value_inQuartile, median_crebbp,mean_crebbp), file = FILE_QUARTILES, 
              append = TRUE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = FALSE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
}
rm(genes, namePathway,cleaned_vector,  genes_list, gene_names, GENE_SETS, TITLE, metada, n, temp, NEWNAME, final, gene2, col_final, coo, ecot, mut, cooi, occurrences, orden, M, col_meta, g, mydata, TMM, df, ensembl95 , geneinfo, df2, df3, df4, genesOfInterest, MM, MM0,clusters, quartiles, labels, quartile_labels, random_string, unique_q, N, median_value_inQuartile,mean_value_inQuartile, median_crebbp, mean_crebbp  )

  } else {
    print("GENE does not exist in the list")
  }
  
}


