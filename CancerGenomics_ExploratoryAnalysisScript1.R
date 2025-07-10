setwd("C:\\Users\\ADMIN\\MyFirstRproject")
getwd()
#Loading the data
my_data <- read.csv("sample_data (1).csv", header=TRUE )
rna_data <- read.table("GSE231379_fpkm_allsamples.txt", header = TRUE)
my_counts <- rna_data
#PREPARING COUNT MATRIX DATAFRAM

my_counts <- as.data.frame(rna_data) #CHANGING VARIABLE TO DATAFRAME

# Load your count matrix (assumes rows = genes, columns = samples)
dim(my_counts)  # Number of genes (rows), number of samples (columns)

# Summary of total library sizes
summary(colSums(my_counts[ ,-1]))

#Summary with boxplot
png( "Boxplot summary of log2expression.png", width = 400, height = 500)
boxplot(log2(my_counts[ , -1] + 1), 
        las = 2, 
        main = "Log2 Raw Counts Distribution per Sample",
        ylab = "Log2 Expression",
        col = "skyblue")
dev.off()

#SETTING COLUMN ORDER 

column_order <- c("Ctrl", "OTA", "MCPD", "TCP")
#SETTING GENE_ID COLUMN AS ROW NAME
rownames(my_counts) <- my_counts$gene_id

#REMOVING GENE_ID COLUMN SINCE THE CONTENT HAS BEEN SET AS ROWNAMES 
my_counts <- my_counts[ ,-1]

#REORDERING COLUMNS
my_counts <- my_counts[, c("Ctrl", "OTA", "MCPD", "TCP")]
rna_data <- my_counts

#DEFINING EXPERIMENTAL GROUPS 
group <- factor(c("Control", "Treatment", "Treatment", "Treatment"))
levels(group)

#CREATING THE EXPERIMENTAL DESIGN DATA
targets <- data.frame(row.names = colnames(my_counts), group)
targets

#SET ROWNAMES OF THE TARGET DATAFRAME AS THE ROWNAMES 
sample_labels <- rownames(targets)

#Making the DGEList
library("edgeR") 
my_DGElist <- DGEList(my_counts)  #converts raw data into a DGEList object for doownstream analysis
my_DGElist$samples$group <- group


#Filter out lowly expressed genes

keep <- filterByExpr(my_DGElist, group = group)
my_DGElist.filtered <- my_DGElist[keep, , keep.lib.sizes = FALSE]

#calculate cpm and transform 
log2_cpm <- cpm(my_DGElist.filtered, log = TRUE) 

#Subset Variable genes for heatmap
# Calculate variance for each gene
gene_var <- apply(log2_cpm, 1, var)

# Select top 100 most variable genes
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:100]
heatmap_data <- log2_cpm[top_genes, ]
head(top_genes)
length(top_genes)
tail(top_genes)

#plot the heatmap

library(pheatmap)
png("Heatmap of top 100 most variable genes.png", width = 400, height = 500)
pheatmap(heatmap_data,
         scale = "row",  # Normalize each gene across samples
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Top 100 Most Variable Genes")
dev.off()


#Hierarchical clustering
#setting groups 
group <- targets$group
group <- factor(group)
View(distance <- dist(t(log2_cpm), method = "euclidean"))
head(distance)
png("Hierachical Clustering of Samples.png",width = 400, height = 500)
distance <- dist(t(log2_cpm), method = "euclidean")
clusters <- hclust(distance,method = "complete")
plot(clusters, labels = sample_labels)

clusters
dev.off()

trans <- t(log2_cpm)

#Principal Component analysis
pca.res <- prcomp(trans,scale.= F, retx = T)
pca.result <-  pca.res$x

#eigenValues
(pc.var <- pca.res$sdev^2)
(pc.percentage <- round(pc.var/sum(pc.var )*100, 1))
library(tidyverse)

pca.res.df <- as_tibble(pca.res$x)

#Plotting Principal Component Analysis
png("PCA.png", width= 400, height = 500)
ggplot(pca.res.df, aes(x=PC1,y=PC2, label = sample_labels,colour = group)) +
  geom_point(size= 4) + geom_label() +
  xlab(paste0("PC1(",pc.percentage[1],"&",")"))+
  ylab(paste0("PC2(",pc.percentage[2],"&",")"))+
  labs(title = "Principal Component Analysis") +
  theme_bw()
dev.off()

#Histogram of log2expressed genes
png("Histogram_log2_expression.png", width = 500, height = 400)
hist(log2(as.matrix(my_counts) + 1),
     breaks = 100,
     col = "lightgreen",
     main = "Histogram of Log2 Expression Values",
     xlab = "Log2(count + 1)")
dev.off()

#Save filtered data for further analysis

write.csv(log2_cpm,"log2_cpm_filtered.csv")

writeLines(capture.output(sessionInfo()), "session_info.txt")



