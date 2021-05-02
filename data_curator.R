
# R packages
library(tidyverse)
library(ggplot2)

# identify path to read and save files
proj_data_path = "/projectnb/bf528/project_4/"
my_data_path = "/projectnb/bf528/users/wheeler/project_4/samples"
output_path = "/projectnb/bf528/users/wheeler/project_4/samples"

samples = c("SRR3879604", "SRR3879605", "SRR3879606")

bc_count = seq_along(samples) %>% 
  map_dfr(
    function(i){
      #i=3;
      count = read.table(file.path(my_data_path, paste0(samples[i], "_bc_count.txt")), header=FALSE, col.names=c("count", "bc")) %>%
        mutate(sample = samples[i], mean=mean(count, na.rm=T))

      return(count)
    }
  )

bc_combined_whitelist <- bc_count %>% filter(count > mean) %>% select(bc) %>% filter(!duplicated(bc)) %>% 
  group_by(bc) %>% 
  summarize(count=n()) %>% 
  select(bc)

write.table(bc_combined_whitelist, file.path(output_path, paste0("combined_whitelist.txt")), col.names=F, row.names=F, quote=F)

cdf_plot <- bc_count %>% 
  ggplot(aes(x=count)) + 
  stat_ecdf(size=1, geom="point", color="black") + 
  geom_hline(yintercept=c(0,1), linetype="dashed") +
  facet_wrap(~sample) +
  xlab("barcode counts") + 
  ylab("percent") +
  theme(
    axis.text=element_text(size=12)
  ) +
  theme_classic()
 
cdf_plot

ggsave(filename=paste0("bc_cdf_plot.png"), plot=cdf_plot, path=output_path, width=7, height=4, units="in", dpi = 300)

## Create the tgMap for salmon alevin
transcripts <- read.table(file.path(my_data_path, "gencode.v37.annotation.gtf.gz"), header=FALSE, sep="\t", col.names=c("chromosome_name", "annotation_source", "feature_type", "genomic_start_location", "genomic_end_location", "score", "genomic_strand", "genomic_phase", "key_value_pairs"))

transcripts_tgmap <- transcripts %>% 
  mutate(
    gene_id=gsub("gene_id", "", grep("gene_id", unlist(strsplit(key_value_pairs, "; ", fixed=TRUE)), value=TRUE)) %>% trimws(),
    transcript_id=gsub("transcript_id", "", grep("transcript_id", unlist(strsplit(key_value_pairs, "; ", fixed=TRUE)), value=TRUE)) %>% trimws(),
    gene_name=gsub("gene_name", "", grep("gene_name", unlist(strsplit(key_value_pairs, "; ", fixed=TRUE)), value=TRUE)) %>% trimws(),
    transcript_name=gsub("transcript_name", "", grep("transcript_name", unlist(strsplit(key_value_pairs, "; ", fixed=TRUE)), value=TRUE)) %>% trimws()
  )

## save the transcript out with associate gene symbols and transcript names
write.table(transcripts_tgmap, file.path(output_path, paste0("transcripts.tsv")), sep="\t", col.names=T, row.names=F)

## save the transcript id and gene id only to run salmon alevin
write.table(transcripts_tgmap %>% dplyr::select(transcript_id, gene_id), file.path(output_path, paste0("gencode.v37.transcripts.tgmap.tsv")), sep="\t", col.names=F, row.names=F, quote=F)

## R packages to do QC metrics on the umi counts
library(tidyverse)
library(ggplot2)
library(tximeta)
library(fishpond)
library(SingleCellExperiment)
library(Seurat)
library(patchwork)
library(edgeR)

# identify path to read and save files
proj_data_path = "/projectnb/bf528/project_4/"
my_data_path = "/projectnb/bf528/users/wheeler/project_4/samples"
output_path = "/projectnb/bf528/users/wheeler/project_4/samples"

samples = c("SRR3879604", "SRR3879605", "SRR3879606")

## import salmon alevin counts
files <- file.path(output_path, "salmon_alevin_output", "alevin", "quants_mat.gz")
file.exists(files)

se <- tximeta(files, type="alevin", alevinArgs=list(filterBarcodes=TRUE))

sce <- as(se, "SingleCellExperiment")

# Get the names in the sce assays
assayNames(sce)

# Get the counts in the sce assays
data <- assays(sce)[["counts"]]

nrow(data); ncol(data)

data[1:20,1:4]

colSums(counts(sce)[,1:4])

# get the transcript that has transcript id and gene id
transcripts <- read.table(file.path(my_data_path, "transcripts.tsv"), header=TRUE, sep="\t") 

# get the gene symbol for the assays data
data_gene_symbol <- data.frame(gene_id=rownames(data)) %>% 
  dplyr::left_join(transcripts) %>% 
  group_by(gene_id, gene_name) %>% 
  summarize(count=n())

# map gene id in the count matrix to gene symbols
rownames(data) <- data_gene_symbol$gene_name[match(rownames(data), data_gene_symbol$gene_id)]

## Setup the Seurat Object
umicount <- CreateSeuratObject(counts = data, min.cells = 10, min.features = 200, project = "GSE84133")

nrow(umicount); ncol(umicount)

umicount[["percent.mt"]] <- PercentageFeatureSet(umicount, pattern = "^MT-")

head(umicount@meta.data, 5)

### Standard pre-processing workflow
VlnPlot(umicount, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(umicount, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(umicount, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## filter cells that have unique feature counts over 4000 or less than 200
filter_umicount_feature <- subset(umicount, subset = nFeature_RNA > 200 & nFeature_RNA < 4000)
nrow(filter_umicount_feature); ncol(filter_umicount_feature)

## filter cells that have >10% mitochondrial counts
filter_umicount <- subset(filter_umicount_feature, subset = percent.mt < 10)
nrow(filter_umicount); ncol(filter_umicount)

## Normalize the umicount
normalized_umicount <- NormalizeData(filter_umicount, normalization.method = "LogNormalize", scale.factor = 10000)

## Identification of highly variable features (feature selection)
highly_variable_feature <- FindVariableFeatures(normalized_umicount, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(highly_variable_feature), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(highly_variable_feature)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

## Scaling the data
scale_featurecount <- ScaleData(highly_variable_feature, features = rownames(highly_variable_feature))

##Perform linear dimensional reduction
pca_featurecount <- RunPCA(scale_featurecount)

print(pca_featurecount[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pca_featurecount, dims = 1:2, reduction = "pca")

DimPlot(pca_featurecount, reduction = "pca")

DimHeatmap(pca_featurecount, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pca_featurecount, dims = 1:15, cells = 500, balanced = TRUE)

# percent of variance explained by each PC
percent_variance_explained <- pca_featurecount@reductions[["pca"]]@stdev/sum(pca_featurecount@reductions[["pca"]]@stdev)*100
percent_variance_explained

### Determine the ‘dimensionality’ of the dataset
resample_featurecount <- JackStraw(pca_featurecount, num.replicate = 100)
resample_featurecount <- ScoreJackStraw(resample_featurecount, dims = 1:20)

JackStrawPlot(resample_featurecount, dims = 1:15)

ElbowPlot(resample_featurecount)

## Cluster the cells
cluster_cells <- FindNeighbors(pca_featurecount, dims = 1:15)
cluster_cells <- FindClusters(cluster_cells, resolution = 0.5)

head(Idents(cluster_cells), 5)

cluster_dt <- data.frame(
  cluster=cluster_cells@meta.data[["seurat_clusters"]]
) %>% 
  group_by(cluster) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  mutate(percent=n/sum(n))

cluster_dt %>% 
  ggplot(aes(x=cluster, y=percent, label = scales::percent(percent))) +
  geom_bar(stat="identity", fill="darkgray", color="white") +
  geom_text(position = position_dodge(width = .9), vjust = -0.5, size = 3) +
  theme_classic()




