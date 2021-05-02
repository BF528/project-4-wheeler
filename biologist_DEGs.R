#loading tidyverse so I can use dpylr to filer and select DEGs
library("tidyverse")

#storing the deferentially expressed genes (DEGs) for the clusters identified
gene_list <- read.csv(file="/projectnb/bf528/project_4_scrnaseq/GSM2230760_marker_genes.csv")


#creating a variable to store the DEGs for each cluster
cluster_0_genes <- filter(gene_list, gene_list$cluster == 0)
cluster_1_genes <- filter(gene_list, gene_list$cluster == 1)
cluster_2_genes <- filter(gene_list, gene_list$cluster == 2)
cluster_3_genes <- filter(gene_list, gene_list$cluster == 3)
cluster_4_genes <- filter(gene_list, gene_list$cluster == 4)
cluster_5_genes <- filter(gene_list, gene_list$cluster == 5)
cluster_6_genes <- filter(gene_list, gene_list$cluster == 6)
cluster_7_genes <- filter(gene_list, gene_list$cluster == 7)
cluster_8_genes <- filter(gene_list, gene_list$cluster == 8)
cluster_9_genes <- filter(gene_list, gene_list$cluster == 9)
cluster_10_genes <- filter(gene_list, gene_list$cluster == 10)
cluster_11_genes <- filter(gene_list, gene_list$cluster == 11)
cluster_12_genes <- filter(gene_list, gene_list$cluster == 12)


#filtered for DEGs with a log fold change greater than 0 and an adjusted p-value less than 0.05
c1_DEGs <- filter(cluster_1_genes, avg_logFC > 0 & p_val_adj < 0.05)

#sorting DEGs based upon significance so the smallest p-values are at the top 
#c1_DEGs <- c1_DEGs[order(c1_DEGs$p_val_adj, decreasing=FALSE), ]  

#printed the filtered DEGs for each cluster on a new line that could readily be pasted into Metascape
cat(c1_DEGs$gene,sep="\n")


c0_DEGs <- filter(cluster_0_genes, avg_logFC > 0 & p_val_adj < 0.05)
cat(c0_DEGs$gene, sep="\n")

c2_DEGs <- filter(cluster_2_genes, avg_logFC > 0 & p_val_adj < 0.05)
cat(c2_DEGs$gene, sep="\n")

c3_DEGs <- filter(cluster_3_genes, avg_logFC > 0 & p_val_adj < 0.05)
cat(c3_DEGs$gene, sep="\n")

c4_DEGs <- filter(cluster_4_genes, avg_logFC > 0 & p_val_adj < 0.05)
cat(c4_DEGs$gene, sep="\n")

c5_DEGs <- filter(cluster_5_genes, avg_logFC > 0 & p_val_adj < 0.05)
cat(c5_DEGs$gene, sep="\n")

c6_DEGs <- filter(cluster_6_genes, avg_logFC > 0 & p_val_adj < 0.05)
cat(c6_DEGs$gene, sep="\n")

c7_DEGs <- filter(cluster_7_genes, avg_logFC > 0 & p_val_adj < 0.05)
cat(c7_DEGs$gene, sep="\n")

c8_DEGs <- filter(cluster_8_genes, avg_logFC > 0 & p_val_adj < 0.05)
cat(c8_DEGs$gene, sep="\n")

c9_DEGs <- filter(cluster_9_genes, avg_logFC > 0 & p_val_adj < 0.05)
cat(c9_DEGs$gene, sep="\n")

c10_DEGs <- filter(cluster_10_genes, avg_logFC > 0 & p_val_adj < 0.05)
cat(c10_DEGs$gene, sep="\n")

c11_DEGs <- filter(cluster_11_genes, avg_logFC > 0 & p_val_adj < 0.05)
cat(c11_DEGs$gene, sep="\n")

c12_DEGs <- filter(cluster_12_genes, avg_logFC > 0 & p_val_adj < 0.05)
cat(c12_DEGs$gene, sep="\n")




