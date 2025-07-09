setwd("C:\\Users\\ADMIN\\MyFirstRproject")
getwd()
#Loading the data
my_data <- read.csv("sample_data (1).csv", header=TRUE )
rna_data <- read.table("GSE231379_fpkm_allsamples.txt", header = TRUE)
#PREPARING COUNT MATRIX DATAFRAME

my_counts <- as.data.frame(rna_data) #CHANGING VARIABLE TO DATAFRAME

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

#calculating CPM using the output of DGEList
cpm <- cpm(my_DGElist)
head(cpm)

#confirm that the cpm is correct 
colSums(cpm)

#transform cpm
log2.cpm <- cpm(my_DGElist, log=TRUE)
#Converting data Matrix to dataframe 

library(tibble)
log2.cpm.df <- as_tibble(log2.cpm, rownames= "gene ID")
colnames(log2.cpm.df) 

library("tidyr")
#converting my data to a pivot table so I can make a plot
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, cols = -1, names_to = "sample",
                                  values_to = "expression")
library("ggplot2")
#making a violin plot for all sample 
plot1 <- ggplot(log2.cpm.df.pivot) + aes(x = sample, y= expression, fill = sample)+
  geom_violin(trim = FALSE, show.legend = FALSE) + 
  stat_summary(fun = "median", geom = "point",
  shape = 95, size = 10, colour = "black", show.legend = FALSE)+
  labs(y = "log2_expression", x = "Samples", title = "Log2 Count Per Million of all Samples", 
       subtitle = "Unfiltered and Non-normalized") + theme_bw() + coord_flip()
plot1
png("Log2 Count Per Million of all Samples.png", width= 5000, height= 6000, res= 600)
plot1
dev.off()

#filtering to remove the lowly expressed genes
cpm <- cpm(my_DGElist)
fil_threshold <- rowSums(cpm>1)>=2
my_DGElist.filtered <- my_DGElist[fil_threshold,]
#preparing filtered data for plotting 
log2.cpm.filtered <- cpm(my_DGElist.filtered,log= TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered,rownames = "gene id")

colnames(log2.cpm.filtered.df) <- c("gene id",sample_labels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, cols = -1,
names_to = "Samples", values_to = "Expression" )
#making a violin plot for filtered data 
plot2 <- ggplot(log2.cpm.filtered.df.pivot) + aes(x = Samples, y= Expression,
  fill = Samples)+
  geom_violin(trim = FALSE, show.legend = FALSE) + 
  stat_summary(fun = "median", geom = "point",
               shape = 95, size = 10, colour = "black",
               show.legend = FALSE)+
  labs(y = "log2_expression", x = "Samples", 
       title = "Log2 Count Per Million of all Samples", 
       subtitle = "filtered and Non-normalized") + theme_bw() + coord_flip()
plot2

#save plot2 
png("Log2 Count Per Million of filtered Samples2.png", width= 5000,height= 6000,res= 600)
plot2
dev.off()

#Normalization 
my_DGElist.filtered.Norm <- calcNormFactors(my_DGElist.filtered, method = "TMM")
#preparing normalized data for plotting 
log2.cpm.filtered.Norm <- cpm(my_DGElist.filtered.Norm,log= TRUE)
log2.cpm.filtered.Norm.df <- as_tibble(log2.cpm.filtered.Norm,rownames = "gene id")

colnames(log2.cpm.filtered.df) <- c("gene id",sample_labels)

#Creating pivot table for normalized data 

log2.cpm.filtered.normalized.df.pivot <- pivot_longer(log2.cpm.filtered.Norm.df, cols = -1,
                                           names_to = "Samples", values_to = "Expression" )
plot3 <- ggplot(log2.cpm.filtered.normalized.df.pivot) + aes(x = Samples, y= Expression,
                                                  fill = Samples)+
  geom_violin(trim = FALSE, show.legend = FALSE) + 
  stat_summary(fun = "median", geom = "point",
               shape = 95, size = 10, colour = "black",
               show.legend = FALSE)+
  labs(y = "log2_expression", x = "Samples", 
       title = "Log2 Count Per Million of all Samples", 
       subtitle = "filtered and normalized") + theme_bw() + coord_flip()
plot3
png("Log2 Count Per Million of filtered and Normalized Samples.png", width= 5000,height= 6000,res= 600)
plot3
dev.off()
png("Log2 Count Per Million of filtered Samples2.png", width= 5000,height= 6000,res= 600)
plot2
dev.off
install.packages("cowplot")
library("cowplot")
norm_plot <- plot_grid(plot1, plot2, plot3, labels = c("A", "B", "C"),label_size = 12)
norm_plot
ggsave("normplot.png", norm_plot, width=10, height=9,dpi=300)

#Hierarchical clustering
#setting groups 
group <- targets$group
View(group <- factor(group))
View(distance <- dist(t(log2.cpm.filtered.Norm), method = "euclidean"))
head(distance)
trans <-  t(log2.cpm.filtered.Norm)
png("Hierachical Clustering of Samples.png",width = 400, height = 500)
View(distance <- dist(trans, method = "euclidean"))
clusters <- hclust(distance,method = "complete")
plot(clusters, labels = sample_labels)

clusters
dev.off()

#Principal Component analysis
pca.res <- prcomp(trans,scale.= F, retx = T)
view(pca.res$x)

#eigenValues
(pc.var <- pca.res$sdev^2)
view(pc.percentage <- round(pc.var/sum(pc.var )*100, 1))
pca.res.df <- as.tibble(pca.res$x)

#Plotting Principal Component Analysis
png("PCA.png", width= 400, height = 500)
ggplot(pca.res.df, aes(x=PC1,y=PC2, label = sample_labels,colour = group)) +
  geom_point(size= 4) + geom_label() +
  xlab(paste0("PC1(",pc.percentage[1],"&",")"))+
ylab(paste0("PC2(",pc.percentage[2],"&",")"))+
labs(title = "Principal Component Analysis") +
theme_bw()
dev.off()
png("PCA 2.png", width= 400, height = 500)
ggplot(pca.res.df, aes(x=PC3,y=PC4, label = sample_labels,colour = group)) +
  geom_point(size= 4) + geom_label() +
  xlab(paste0("PC3(",pc.percentage[3],"&",")"))+
  ylab(paste0("PC4(",pc.percentage[4],"&",")"))+
  labs(title = "Principal Component Analysis") +
  theme_bw()
dev.off()
