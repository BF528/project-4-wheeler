#libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(bseqsc)
library(purrr)
library(ggplot2)

#load the data
cells <- readRDS("/Users/arielxue/Documents/bf528/project4/GSM2230760_seurat.rda")
str(cells)
#find all markers
cells.markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cells.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
#save the markers
saveRDS(cells.markers, file = "/Users/arielxue/Documents/bf528/project4/markers.rds")

#find top10 markers
gen_marker_table <- function(x){
  cells.markers[cells.markers$cluster == x, ] %>%
    head(n=10)
}
top10_markers <- map_dfr(0:12, gen_marker_table)

#trying to sign marker genes to cell types
FeaturePlot(cells, features = c(top10_markers[top10_markers$cluster == 2, "gene"]), label = TRUE)
markers_2vs9 <- FindMarkers(cells, ident.1 = 2, ident.2 = 9, label = TRUE)

VlnPlot(cells, features = c("GCG")) #Alpha
VlnPlot(cells, features = c("SST")) # Delta 
VlnPlot(cells, features = c("INS")) # Beta
VlnPlot(cells, features = c("PPY")) # Gamma 
VlnPlot(cells, features = c("KRT19")) # Ductal
VlnPlot(cells, features = c("CPA1")) # Acinar
VlnPlot(cells, features = c("SDS")) # Macrophage
VlnPlot(cells, features = c('GHRL','RGS5','PDGFRA','VWF','TRAC')) #none consistent clusters

#find marker genes for cluster 7 and 9
cluster7.markers <- FindMarkers(cells, ident.1 = 7, min.pct = 0.25)
cluster9.markers <- FindMarkers(cells, ident.1 = 9, min.pct = 0.25)

#heatmap for top10 marker genes in each cluster
top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(cells, features = top10$gene) + NoLegend()

#log normalize the counts
cells <- NormalizeData(cells, normalization.method = "LogNormalize", scale.factor = 10000)

#heatmap for top2 marker genes in each cluster
top2 <- cells.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
DoHeatmap(cells, features = top2$gene, slot = "counts") + scale_fill_gradientn(colors = c("white","purple"))

#feature plot for each cell type's marker gene
FeaturePlot(cells, features = c("GCG","INS",'SST','PPY','GHRL','KRT19'), label = TRUE)
FeaturePlot(cells, features = c('CPA1','RGS5','PDGFRA','VWF','SDS','TRAC'), label = TRUE)

#signing cell type to each cluster
new.cluster.ids <- c("Delta/Gamma", "Beta", "Alpha", "Delta/Acinar", "Alpha", "Ductal","Beta", "7", "Alpha","9","Ductal","Acinar","Macrophage")
names(new.cluster.ids) <- levels(cells)
cells <- RenameIdents(cells, new.cluster.ids)
DimPlot(cells, reduction = "umap", label = TRUE, pt.size = 0.5)

#umap
cell <- RunUMAP(cells, dims = 1:10)
DimPlot(cell, reduction = "umap", label = TRUE)
