#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=budankm@ccf.org
#SBATCH --job-name=scrna_ccf
#SBATCH --time=10:00:00
#SBATCH --partition=bigmem
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --output=/home/budankm/isilon/Meghana/slurm/slurm_out_%j.out
#SBATCH --error=/home/budankm/isilon/Meghana/slurm/slurm_err_%j.err
#####################################################################
source /home/budankm/miniforge3/bin/activate
conda activate bioinformatics
#####################################################################
cd /home/budankm/isilon/XieLab_NGSData/PIP_9_12_2024/Demultiplexed

# Create necessary directories
mkdir -p trimmed_reads fastp_qc_output raw_qc_output

# Loop over all R1 FASTQ files in the folder and process their paired R2 files
for r1 in Meghana_*_L008_R1_001.fastq.gz; do
    sample_name=$(basename ${r1} _L008_R1_001.fastq.gz)
    r2="${sample_name}_L008_R2_001.fastq.gz"  

    echo "Processing sample: $sample_name"

    # Run FastQC on raw files
    fastqc -t 4 --outdir=raw_qc_output ${r1} ${r2}

    # Decompress the FASTQ files
    gunzip -c ${r1} > trimmed_reads/${sample_name}_R1.fastq
    gunzip -c ${r2} > trimmed_reads/${sample_name}_R2.fastq

    # Run FastP trimming
    fastp -i trimmed_reads/${sample_name}_R1.fastq -o trimmed_reads/${sample_name}_R1_trimmed.fastq \
          -I trimmed_reads/${sample_name}_R2.fastq -O trimmed_reads/${sample_name}_R2_trimmed.fastq \
          --cut_front --cut_tail --qualified_quality_phred 20 --detect_adapter_for_pe \
          --disable_length_filtering --thread 4 --html fastp_qc_output/fastp_report_${sample_name}.html \
          --json fastp_qc_output/fastp_report_${sample_name}.json

    # Run FastQC on trimmed files
    fastqc -t 4 --outdir=fastp_qc_output trimmed_reads/${sample_name}_R1_trimmed.fastq trimmed_reads/${sample_name}_R2_trimmed.fastq
done

# === Step 3: MultiQC ===
multiqc raw_qc_output fastp_qc_output "$OUTPUT_DIR"

# Check if STAR is available in the Conda environment
echo "Checking if STAR is installed and available in the bioinformatics environment..."
which STAR

# If STAR is available, display its version
if which STAR > /dev/null; then
  echo "STAR is available at: $(which STAR)"
  STAR --version
else
  echo "STAR is not available in this environment. Please check your Conda installation."
  exit 1
fi

# Define directories
STAR_INDEX_DIR="/home/budankm/isilon/Meghana/starindexGRCh38"  
FASTQ_DIR="/home/budankm/isilon/XieLab_NGSData/250318_lh00134_0671_A22MKLCLT4/DemultiplexPIP_ATAC_MB/trimmed_reads" 
OUTPUT_DIR="/home/budankm/isilon/Meghana/star_output"  

# Define fastq files for input
FASTQ_R1="${FASTQ_DIR}/*_R1_trimmed.fastq"
FASTQ_R2="${FASTQ_DIR}/*_R2_trimmed.fastq"

# Loop over each pair of FASTQ files
for R1 in $FASTQ_R1; do
 
  R2="${R1/_R1_trimmed.fastq/_R2_trimmed.fastq}"
  

  SAMPLE_NAME=$(basename "$R1" _R1_trimmed.fastq)
  

  OUTPUT_PREFIX="${OUTPUT_DIR}/${SAMPLE_NAME}_aligned"


# Run STAR genome generation
STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir /home/budankm/isilon/Meghana/starindexGRCh38 \
--genomeFastaFiles /home/budankm/isilon/Meghana/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /home/budankm/isilon/Meghana/Homo_sapiens.GRCh38.113.gtf \
--sjdbOverhang 150 

# Run STAR for alignment
  /home/budankm/miniforge3/envs/bioinformatics/bin/STAR --runThreadN 12 \
    --genomeDir "$STAR_INDEX_DIR" \
    --readFilesIn "$R1" "$R2" \
    --outFileNamePrefix "$OUTPUT_PREFIX" \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outSAMattributes NH HI AS nM MD \
    --sjdbOverhang 150


######################################################################
#Counts were generated in R and downstream analysis was performed. 
#scrna
#####################

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)

# List of sample folders
sample_dirs <- list(
  "Anterior" = "/Users/budankm/Desktop/sensitivity_5",
  "Lateral" = "/Users/budankm/Desktop/sensitivity_5e",
  "Medial" = "/Users/budankm/Desktop/sensitivity_5f",
  "Posterior" = "/Users/budankm/Desktop/sensitivity_5g",
  "Intra_Tumor" = "/Users/budankm/Desktop/sensitivity_5h"
)

output_dir <- "/Users/budankm/Desktop/output_counts"  
dir.create(output_dir, showWarnings = FALSE)

# Initialize list to store Seurat objects and pseudo-bulk counts
seurat_list <- list()
pseudo_bulk_list <- list()

# Loop through each sample
for (sample_name in names(sample_dirs)) {
  data_path <- sample_dirs[[sample_name]]
  
  # Read 10X data
  data <- Read10X(
    data.dir = data_path,
    gene.column = 2,
    cell.column = 1,
    unique.features = TRUE,
    strip.suffix = FALSE
  )
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_name)
  seurat_obj$sample <- sample_name  # Add sample ID
  
  # Save single-cell counts
  single_cell_counts <- as.matrix(seurat_obj@assays$RNA@counts)
  write.csv(pseudo_bulk_df, file = file.path(output_dir, paste0(sample_name, "_single_bulk_counts.csv")), row.names = FALSE)
  # Create pseudo-bulk (sum across all cells)
  pseudo_bulk_counts <- rowSums(single_cell_counts)
  
  # Save pseudo-bulk counts
  pseudo_bulk_df <- data.frame(
    gene = names(pseudo_bulk_counts),
    counts = pseudo_bulk_counts
  )
  write.csv(pseudo_bulk_df, file = file.path(output_dir, paste0(sample_name, "_pseudo_bulk_counts.csv")), row.names = FALSE)
  

  # Store Seurat object and pseudo-bulk vector in lists
  seurat_list[[sample_name]] <- seurat_obj
  pseudo_bulk_list[[sample_name]] <- pseudo_bulk_counts
}

# combine all pseudo-bulk samples into one matrix
pseudo_bulk_matrix <- do.call(cbind, pseudo_bulk_list)
pseudo_bulk_matrix <- as.data.frame(pseudo_bulk_matrix)
pseudo_bulk_matrix$gene <- rownames(pseudo_bulk_matrix)

# Save the full pseudo-bulk matrix
write.csv(pseudo_bulk_matrix, file = file.path(output_dir, "All_Samples_PseudoBulk_Matrix.csv"), row.names = FALSE)
seurat_list <- list()

seurat_list <- list()

for (sample_name in names(sample_dirs)) {
  data_path <- sample_dirs[[sample_name]]
  data <- Read10X(data.dir = data_path, gene.column = 2,
                  cell.column = 1,
                  unique.features = TRUE,
                  strip.suffix = FALSE)
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_name)
  seurat_obj$sample <- sample_name  # Add sample ID as metadata
  seurat_list[[sample_name]] <- seurat_obj
}


combined_seurat <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "AllSamples"
)
#

# Add percent mitochondrial content
combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^MT[-\\.]")

# Create violin plot for pre-QC metrics
preqc_vlnplot <- VlnPlot(
  combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample",  
  ncol = 3,  
  pt.size = 0.1  
) + 
  theme(
    axis.title.x = element_text(size = 12),  
    axis.title.y = element_text(size = 12),  
    axis.text.x = element_text(angle = 45, hjust = 1),  
    strip.text.x = element_text(size = 12)  
  ) +
  scale_y_continuous(labels = scales::comma)  


preqc_vlnplot <- preqc_vlnplot + 
  plot_annotation(title = "Initial Pre-QC Metrics") & 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

# Display the plot
print(preqc_vlnplot)

# Save the plot as a PDF
ggsave("~/Desktop/plots/preqc_vlnplot.pdf", plot = preqc_vlnplot, width = 10, height = 5)


#############################
# Subset to retain high-quality cells:
# - More than 200 detected genes
# - Fewer than 7500 detected genes
# - Less than 10% mitochondrial content
combined_sub <- subset(
  combined,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 7500 &
    percent.mt < 10
)

paste("Before:", ncol(combined), " After:", ncol(combined_sub))

postqc_vlnplot <- VlnPlot(
  combined_sub,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample",  
  ncol = 3,  
  pt.size = 0.1
) + 
  theme(
    axis.title.x = element_text(size = 12),  
    axis.title.y = element_text(size = 12),  
    axis.text.x = element_text(angle = 45, hjust = 1),  
    strip.text.x = element_text(size = 12)  
  ) +
  scale_y_continuous(labels = scales::comma)  

postqc_vlnplot <- postqc_vlnplot + 
  plot_annotation(title = "Post-QC Metrics") & 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))


print(postqc_vlnplot)

ggsave("~/Desktop/plots/postqc_vlnplot.pdf", plot = postqc_vlnplot , width = 10, height = 5)

############################

genecount_vlnplot <- VlnPlot(combined_sub, features = "nFeature_RNA") +
  labs(
    title = "Gene Count per Cell",
    x = "Tumor Samples",
    y = "# of genes per cell"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Title size
    axis.title.x = element_text(size = 14),  # X-axis title size
    axis.title.y = element_text(size = 14),  # Y-axis title size
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # X-axis text
    axis.text.y = element_text(size = 12)   # Y-axis text
  )

print(genecount_vlnplot)

ggsave("~/Desktop/plots/genecount_vlnplot.pdf", plot = genecount_vlnplot , width = 10, height = 5)

###########################
#quantifying cells and genes

gene_counts <- combined_sub@meta.data %>%
  group_by(orig.ident) %>%
  summarise(avg_genes = mean(nFeature_RNA))

p1 <- ggplot(gene_counts, aes(x = orig.ident, y = avg_genes, fill = orig.ident)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(title = "Average Genes per Cell",
       x = "Sample", y = "Avg. Genes") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Total number of cells per sample
cell_counts <- combined_sub@meta.data %>%
  count(orig.ident)

p2 <- ggplot(cell_counts, aes(x = orig.ident, y = n, fill = orig.ident)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(title = "Number of Cells per Sample",
       x = "Sample", y = "Cell Count") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Combine the plots
p1 + p2 + plot_layout(ncol = 2)
combined_qc_plot <- p1 + p2 + plot_layout(ncol = 2)

# Save to Desktop/plots
ggsave("~/Desktop/plots/combined_qc_plot.pdf",
       plot = combined_qc_plot,
       width = 10, height = 5)
###########################
library(Seurat)
library(ggplot2)
library(patchwork)

# First plot: nCount_RNA vs percent.mt
plot1 <- FeatureScatter(
  combined_sub,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
) +
  labs(
    title = "RNA Count vs Mitochondrial Content",
    x = "RNA Count (nCount_RNA)",
    y = "Mitochondrial Percentage (percent.mt)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12)
  )

# Second plot: nCount_RNA vs nFeature_RNA
plot2 <- FeatureScatter(
  combined_sub,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
) +
  labs(
    title = "RNA Count vs Gene Features",
    x = "RNA Count (nCount_RNA)",
    y = "Detected Genes (nFeature_RNA)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12)
  )

# Combine the plots side-by-side with patchwork
combined_scatter <- plot1 + plot2 +
  plot_layout(guides = "collect") & theme(legend.position = "right")

# Print the combined plot
print(combined_scatter)

# Optional: Save to file
ggsave("~/Desktop/plots/qc_scatter_plots.pdf", plot = combined_scatter, width = 12, height = 5)

############################## QC POST scattering
# 1. Subset High-Quality Cells
seurat <- subset(
  combined_sub,
  subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5
)

VlnPlot(seurat, features = c("nFeature_RNA", "percent.mt"))

# 2. Normalize and Identify Variable Features

seurat <- NormalizeData(seurat, verbose = FALSE)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)


# 3. Visualize Top Variable Features

top_features <- head(VariableFeatures(seurat), 10)
var_plot <- VariableFeaturePlot(seurat)


highlighted_plot <- LabelPoints(
  plot = var_plot, 
  points = top_features, 
  repel = TRUE,
  xnudge = 0, 
  ynudge = 0
)
print(var_plot + highlighted_plot + plot_layout(guides = "collect"))


# 4. Scale the Data

all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes, verbose = FALSE)


# 5. PCA (Principal Component Analysis)


seurat <- RunPCA(seurat, features = VariableFeatures(seurat), verbose = FALSE)

# View PCA details
print(seurat[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(seurat, dims = 1:2, reduction = "pca") + 
  ggtitle("PCA Loadings") + 
  theme_minimal()

# Visualize PCA clusters
DimPlot(seurat, reduction = "pca", label = TRUE) +
  ggtitle("PCA Clustering") +
  theme_minimal()

# PCA Heatmaps
DimHeatmap(seurat, dims = 1:5, balanced = TRUE) +
  ggtitle("PCA Heatmap (Top PCs)")
DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat, dims = 1:10, cells = 500, balanced = TRUE)

# Elbow Plot to determine optimal PC count
ElbowPlot(seurat)


# 6. Graph-Based Clustering

seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat, resolution = 0.8)


# 7. UMAP for Dimensionality Reduction

seurat <- RunUMAP(seurat, dims = 1:10)

# Visualize UMAP with Cluster Labels

DimPlot(seurat, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("UMAP: Clusters of Tumor Cells") +
  theme_minimal()
##################

# Create a directory to store plots (optional)
if (!dir.exists("plots")) dir.create("plots")

# ---- 1. Variable Feature Plot ----
var_plot <- VariableFeaturePlot(seurat)
highlighted_plot <- LabelPoints(plot = var_plot, points = top_features, repel = TRUE)
combined_var_plot <- var_plot + highlighted_plot + plot_layout(guides = "collect")
ggsave("plots/variable_features.pdf", plot = combined_var_plot, width = 10, height = 6)

# ---- 2. PCA Loadings Plot ----
pca_loadings <- VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
ggsave("plots/pca_loadings.pdf", plot = pca_loadings, width = 8, height = 6)

# ---- 3. PCA DimPlot ----
pca_plot <- DimPlot(seurat, reduction = "pca", label = TRUE) +
  ggtitle("PCA Clustering") +
  theme_minimal()
ggsave("plots/pca_clusters.pdf", plot = pca_plot, width = 8, height = 6)

# ---- 4. PCA Heatmaps ----
pdf("plots/pca_heatmaps.pdf", width = 10, height = 6)
DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat, dims = 1:10, cells = 500, balanced = TRUE)
dev.off()

# ---- 5. Elbow Plot ----
elbow <- ElbowPlot(seurat)
ggsave("plots/elbow_plot.pdf", plot = elbow, width = 6, height = 5)

# ---- 6. UMAP Plot ----
umap <- DimPlot(seurat, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("UMAP: Clusters of Tumor Cells") +
  theme_minimal()
ggsave("plots/umap_clusters.pdf", plot = umap, width = 8, height = 6)


###################

# 1. Identify markers for each cluster
cl_markers <- FindAllMarkers(
  seurat,
  only.pos = TRUE,               # Only return upregulated markers
  min.pct = 0.25,                # Minimum % of cells expressing the gene in either group
  logfc.threshold = log(1.2)     # Log fold-change threshold
)

head(cl_markers)
table(cl_markers$cluster)
##########################

# 4. View top markers per cluster
library(dplyr)
# Get top 10 markers per cluster
top_markers <- cl_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup() %>%
  dplyr::select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  arrange(cluster, desc(avg_log2FC))
############################

top_genes_all_clusters <- cl_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%  # Select top 10 genes per cluster
  ungroup() %>%
  arrange(desc(avg_log2FC))  # Sort all top genes by log fold-change in descending order

# 2. View the top genes for all clusters
head(top_genes_all_clusters)

write.csv(top_genes_all_clusters, "~/Desktop/plots/top_genes_all_clusters.csv", row.names = FALSE)

# Select the top genes from the 'top_genes_all_clusters' dataframe
top_genes_list <- top_genes_all_clusters$gene

# Create a heatmap of the top genes across all clusters
DoHeatmap(seurat, features = top_genes_list) + 
  theme_minimal() + 
  ggtitle("Top Genes Across All Clusters")


#######################

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
all_marker_genes <- cl_markers$gene
ego_all_markers <- enrichGO(gene = all_marker_genes, 
                            OrgDb = org.Hs.eg.db, 
                            keyType = "SYMBOL", 
                            ont = "BP", 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05)

# Extract the top 10 enriched GO terms based on adjusted p-value
top_go <- ego_all_markers@result[order(ego_all_markers@result$p.adjust), ][1:10, ]

# Bar plot for top 10 enriched GO terms
top10enrichedgenes <- ggplot(top_go, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", fill = "lightpink") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +  
  labs(x = "GO Term", y = "-log10(adjusted p-value)", title = "Top 10 Enriched GO Terms (All Markers)")

print(top10enrichedgenes)

# Dot plot for top 20 enriched GO terms
dotplottop10enrichedgenes <- dotplot(ego_all_markers, showCategory = 20) +
  theme(axis.text.y = element_text(size = 8))  

# Show the dot plot
print(dotplottop10enrichedgenes)
ggsave("~/Desktop/plots/top10_enriched_barplot.pdf", plot = top10enrichedgenes, width = 8, height = 6, dpi = 300)

# Save the dot plot
ggsave("~/Desktop/plots/top20_enriched_dotplot.df", plot = dotplottop10enrichedgenes, width = 8, height = 6, dpi = 300)
############################
# enrich by cluster

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# 1. Make cluster-to-orig.ident mapping
cluster_to_origident <- seurat@meta.data %>%
  dplyr::select(cluster = seurat_clusters, orig.ident) %>%
  distinct()

# 2. Merge back into cl_markers
cl_markers_with_orig <- cl_markers %>%
  left_join(cluster_to_origident, by = "cluster")

# 3. Split by orig.ident
marker_list_by_origident <- split(cl_markers_with_orig$gene, cl_markers_with_orig$orig.ident)

# 4. Loop through each orig.ident and do GO enrichment
dir.create("~/Desktop/plots/GO_by_origident", recursive = TRUE)

for (ident_name in names(marker_list_by_origident)) {
  
  genes <- marker_list_by_origident[[ident_name]]
  
  # Perform GO enrichment
  ego <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  
  # Skip if no enrichment results
  if (is.null(ego) || nrow(ego@result) == 0) {
    message(paste("No enrichment found for", ident_name))
    next
  }
  
  # Extract top 10 GO terms
  top_go <- ego@result %>% arrange(p.adjust) %>% head(10)
  
  # Bar plot
  barplot <- ggplot(top_go, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
    geom_bar(stat = "identity", fill = "lightpink") +
    coord_flip() +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8)) +
    labs(
      x = "GO Term",
      y = "-log10(adjusted p-value)",
      title = paste("Top 10 Enriched GO Terms:", ident_name)
    )
  
  # Dot plot
  dotplot_plot <- dotplot(ego, showCategory = min(20, nrow(ego@result))) +
    theme(axis.text.y = element_text(size = 8)) +
    ggtitle(paste("Top Enriched GO Terms:", ident_name))
  
  # Save the plots
  ggsave(
    filename = paste0("~/Desktop/plots/GO_by_origident/", ident_name, "_barplot.pdf"),
    plot = barplot,
    width = 8, height = 6, dpi = 300
  )
  
  ggsave(
    filename = paste0("~/Desktop/plots/GO_by_origident/", ident_name, "_dotplot.pdf"),
    plot = dotplot_plot,
    width = 8, height = 6, dpi = 300
  )
  
  # Optionally print them
  print(barplot)
  print(dotplot_plot)
}

#########################################################################

###########################
library(Seurat)
library(celldex)
library(ggplot2)
library(SingleR)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
expr_matrix <- seurat@assays$RNA@data
ref_data <- celldex::HumanPrimaryCellAtlasData()
predictions <- SingleR(test = expr_matrix, ref = ref_data, labels = ref_data$label.main)
head(predictions)
str(predictions)
predicted_cell_types <- predictions$pruned.labels
prediction_scores <- predictions$scores
head(predicted_cell_types)

seurat$SingleR_labels <- predicted_cell_types

# Create count table
celltype_counts <- as.data.frame(table(seurat_filtered$orig.ident, seurat_filtered$SingleR_labels))
colnames(celltype_counts) <- c("Sample", "CellType", "Count")

# Set factor levels for nicer legend sorting
celltype_counts$CellType <- factor(celltype_counts$CellType, 
                                   levels = names(sort(tapply(celltype_counts$Count, celltype_counts$CellType, sum), decreasing = TRUE)))

# Color palette
n_types <- length(unique(celltype_counts$CellType))
my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_types)

# Plot
celltype_barplot <- ggplot(celltype_counts, aes(x = Sample, y = Count, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  theme_minimal() +
  labs(title = "Cell Type Composition per Sample",
       x = "Sample",
       y = "Number of Cells",
       fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(celltype_barplot)
ggsave("~/Desktop/plots/Cell_Type_BarPlot_noNA.png", plot = celltype_barplot, width = 10, height = 7, dpi = 300)

# Plot 
Cell_Type_Annotations <- DimPlot(seurat, group.by = "SingleR_labels", 
  label = TRUE, repel = TRUE, label.size = 3, cols = my_colors) +
  ggtitle("Cell Type Annotations") +
  theme_minimal()
print(Cell_Type_Annotations)
ggsave("~/Desktop/plots/Cell_Type_Annotations.png", plot = Cell_Type_Annotations, width = 10, height = 8, dpi = 300)



###############################

# Convert the table to a data frame
annotation_table <- as.data.frame(table(seurat$seurat_clusters, seurat$SingleR_labels))
colnames(annotation_table) <- c("Seurat_Cluster", "SingleR_Label", "Cell_Count")

# Plot the heatmap
library(ggplot2)

celltypes_heatmap <- ggplot(annotation_table, aes(x = Seurat_Cluster, y = SingleR_Label, fill = Cell_Count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of Cell Type Distribution per Seurat Cluster",
       x = "Seurat Cluster", 
       y = "SingleR Cell Type")
ggsave("~/Desktop/plots/celltypes_heatmap.pdf", plot = celltypes_heatmap, width = 10, height = 8, dpi = 300)
###################
Sys.getenv("R_LIBS_USER")
#[1] "/Users/budankm/Library/R/arm64/4.4/library"
dir.exists(Sys.getenv("R_LIBS_USER"))
.libPaths(Sys.getenv("R_LIBS_USER"))
#BiocManager::install("ComplexHeatmap", force = TRUE)
library(ComplexHeatmap)

# Prepare the data
annotation_matrix <- as.matrix(table(seurat$seurat_clusters, seurat$SingleR_labels))
rownames(annotation_matrix) <- paste0("Cluster_", rownames(annotation_matrix))
colnames(annotation_matrix) <- gsub(" ", "_", colnames(annotation_matrix))  # Clean column names
head(annotation_matrix)

##############################

install.packages("pheatmap")
library(pheatmap)

# Normalize the matrix by row to show relative proportions
annotation_matrix_norm <- prop.table(annotation_matrix, margin = 1)
# Plot the heatmap
Heatmap(annotation_matrix, name = "Cell_Count",
        row_title = "Seurat Cluster", column_title = "SingleR Labels",
        cluster_rows = FALSE, cluster_columns = FALSE,
        col = colorRampPalette(c("white", "blue"))(100))


######################################################################################################
pseudo_bulk_counts <- apply(raw_counts, 1, function(x) tapply(x, grouping, sum))

# Transpose pseudo_bulk_counts to match DESeq2 input format (genes as rows, samples as columns)
pseudo_bulk_counts_transposed <- t(pseudo_bulk_counts)

# Create a metadata data frame with sample conditions
sample_info <- data.frame(
  condition = c("anterior", "posterior", "medial", "lateral", "intra_tumor"),
  row.names = colnames(pseudo_bulk_counts_transposed)  # Ensure row names match the column names of pseudo_bulk_counts
)
head(sample_info)
library(DESeq2)
# don't have replicates:
dds <- DESeqDataSetFromMatrix(
  countData = pseudo_bulk_counts_transposed, 
  colData = sample_info, 
  design = ~ 1  # No condition
)

# Run DESeq2 without conditions for simple estimation
dds <- DESeq(dds)
results_dds <- results(dds)

# View the results
head(results_dds)
plotMA(results_dds, ylim = c(-5, 5))
volcano_data_clean <- volcano_data[!is.na(volcano_data$padj) & !is.na(volcano_data$log2FoldChange), ]
# Now plot the volcano plot with cleaned data
ggplot(volcano_data_clean, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
    geom_point() + theme_minimal()

rlog_counts <- rlog(dds, blind = FALSE)

#generate the heatmap for top genes
top_genes <- head(order(results_dds$padj), 20)  
heatmap_data <- assay(rlog_counts)[top_genes, ]
pheatmap(heatmap_data)

pheatmap(heatmap_data, fontsize_row = 6)    

######################################################################################################
# ssgsea
# BiocManager::install("GSVA")
# BiocManager::install("clusterProfiler")

# Load necessary libraries
library(Seurat)
library(GSVA)
library(clusterProfiler)
library(msigdbr)

normalized_counts <- GetAssayData(combined_sub, slot = "data")
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")

unique_gene_sets <- unique(hallmark_sets$gs_name)
# Define the list of selected gene sets
selected_gene_sets <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_HYPOXIA",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_PROTEIN_SECRETION",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_MYC_TARGETS_V2"
)

# Filter the hallmark_sets dataframe for the selected gene sets
filtered_sets <- hallmark_sets[hallmark_sets$gs_name %in% selected_gene_sets, ]

# Create a gene set list by splitting the data based on gs_name
gene_sets_list <- split(filtered_sets$gene_symbol, filtered_sets$gs_name)

# View the gene set list
gene_sets_list



















  # Assuming results_dds is the result from DESeq2 or equivalent
  # Sort by log2FoldChange (you could use other statistics like p-value, but fold change is typical for GSEA)
  gene_list <- results_dds$log2FoldChange
  names(gene_list) <- rownames(results_dds)
  gene_list <- sort(gene_list, decreasing = TRUE)  # GSEA prefers descending order
  library(clusterProfiler)
  # Add small random noise to the gene list to break ties
  set.seed(42)  # Set a seed for reproducibility
  gene_list <- gene_list + rnorm(length(gene_list), mean = 0, sd = 1e-6)
  
  # Now run fgsea again
  fgsea_results <- fgsea(
    pathways = gene_sets,
    stats = gene_list,
    minSize = 15,
    maxSize = 500,
    nPermSimple = 1000
  )
  

fgsea_results_clean <- fgsea_results[!is.na(fgsea_results$NES) & is.finite(fgsea_results$NES), ]
top_fgsea_results <- fgsea_results_clean[order(fgsea_results_clean$NES, decreasing = TRUE), ]
  


library(GSVA)

# Make the parameter object correctly
ssgseaPar <- ssgseaParam(expression_matrix_sparse, geneSets = gene_sets_list)

# Check it
ssgseaPar
library(pheatmap)
> 
> # Scale scores if you want better contrast (optional)
> scaled_scores <- t(scale(t(ssgsea_es)))
> 
> # Draw heatmap
> pheatmap(
+     scaled_scores,
+     cluster_rows = TRUE,
+     cluster_cols = TRUE,
+     show_rownames = TRUE,
+     show_colnames = TRUE,
+     fontsize_row = 8,
+     fontsize_col = 10,
+     color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
+ )
# Pick a sample
> sample_id <- colnames(ssgsea_es)[1]
> 
> # Get top pathways
> top_pathways <- ssgsea_es[, sample_id] %>%
+     sort(decreasing = TRUE) %>%
+     head(20)
> 
> # Plot
> ggplot(data = data.frame(Pathway = names(top_pathways), Score = top_pathways), aes(x = reorder(Pathway, Score), y = Score)) +
+     geom_bar(stat = "identity", fill = "steelblue") +
+     coord_flip() +
+     theme_minimal() +
+     xlab("Pathway") +
+     ylab("ssGSEA Score") +
+     ggtitle(paste("Top 20 Pathways in", sample_id))
> metadata <- data.frame(
+     cell = colnames(combined_sub),
+     original.ident = combined_sub$orig.ident
+ )
> 
> View(metadata)
> # Make sure metadata is a data.frame
> metadata <- as.data.frame(metadata)
> 
> # Now reorder metadata rows to match the order of ssgsea_es columns
> metadata <- metadata[match(colnames(ssgsea_es), metadata$cell), ]
> 
> # Check if they match
> stopifnot(all(metadata$cell == colnames(ssgsea_es)))
>  library(dplyr)
> 
> # Add tumor identity info to your ssGSEA matrix
> ssgsea_df <- as.data.frame(t(ssgsea_es))  # transpose so cells are rows
> ssgsea_df$original.ident <- metadata$original.ident  # add tumor info
> 
> # Average the pathway scores per tumor
> ssgsea_avg <- ssgsea_df %>%
+     group_by(original.ident) %>%
+     summarise(across(everything(), mean)) %>%
+     as.data.frame()
> 
> # Set tumor names as rownames
> rownames(ssgsea_avg) <- ssgsea_avg$original.ident
> ssgsea_avg <- ssgsea_avg[, -1]  # remove the 'original.ident' column
> 
> library(ggplot2)
> # Set how many top pathways you want
> top_n <- 10
> 
for (tumor in rownames(ssgsea_avg)) {
+     
+     # Extract the scores as a numeric named vector
+     tumor_scores <- as.numeric(ssgsea_avg[tumor, ])
+     names(tumor_scores) <- colnames(ssgsea_avg)
+     
+     # Sort to get top pathways
+     top_pathways <- sort(tumor_scores, decreasing = TRUE)[1:top_n]
+     
+     # Create a data frame for plotting
+     plot_df <- data.frame(
+         pathway = names(top_pathways),
+         score = as.numeric(top_pathways)
+     )
+     
+     # Create plot
+     p <- ggplot(plot_df, aes(x = reorder(pathway, score), y = score, fill = score)) +
+         geom_bar(stat = "identity", show.legend = FALSE) +
+         coord_flip() +
+         theme_minimal() +
+         labs(title = paste("Top", top_n, "Pathways in", tumor),
+              x = "Pathways", y = "Mean ssGSEA Score") +
+         scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
+         theme(axis.text = element_text(size = 12),
+               axis.title = element_text(size = 14))
+     
+     print(p)
+     
+     # Save the plot
+     ggsave(filename = paste0("~/Desktop/plots", top_n, "_Pathways_", tumor, ".pdf"), plot = p, width = 8, height = 6)
+     
+ }
> library(pheatmap)
library(dplyr)

# Define tumor sample names
tumor_samples <- c("Anterior", "Lateral", "Medial", "Posterior", "Intra_Tumor")

# Set the number of top genes you want to visualize
top_n <- 20

# Loop through each tumor sample
for (tumor in tumor_samples) {
    
    # Get the cells belonging to this specific tumor sample
    tumor_cells <- metadata$cell[metadata$original.ident == tumor]
    
    # Extract ssGSEA scores for the tumor sample (for those specific cells)
    tumor_ssgsea_scores <- ssgsea_es[, tumor_cells]
    
    # Calculate the mean ssGSEA score for each gene across all cells in this tumor sample
    gene_avg_scores <- rowMeans(tumor_ssgsea_scores, na.rm = TRUE)
    
    # Sort the genes by mean ssGSEA score in descending order
    sorted_genes <- sort(gene_avg_scores, decreasing = TRUE)
    
    # Get the top N enriched genes
    top_enriched_genes <- names(sorted_genes)[1:top_n]
    
    # Check if the top genes are in the ssGSEA scores matrix
    top_enriched_genes <- top_enriched_genes[top_enriched_genes %in% rownames(tumor_ssgsea_scores)]
    
    # Extract the ssGSEA scores for the top genes
    top_genes_ssgsea_scores <- tumor_ssgsea_scores[top_enriched_genes, , drop = FALSE]
    
    # Create and save the heatmap directly using pheatmap
    pheatmap(top_genes_ssgsea_scores, 
              cluster_rows = TRUE, 
              cluster_cols = TRUE, 
              scale = "row", 
              show_rownames = TRUE, 
              show_colnames = TRUE,    # Show cell/sample names on the top
              main = paste("ssGSEA Scores - Top", top_n, "Genes in", tumor),
              color = colorRampPalette(c("blue", "white", "red"))(100),  # Valid color palette
              fontsize = 10,
              fontsize_row = 8,        # Adjust the font size for row labels (gene names)
              fontsize_col = 10,       # Adjust the font size for column labels (cells)
              angle_col = 45,          # Rotate column names to 45 degrees for better readability
              legend = TRUE,           # Show legend for colors
              filename = paste0("~/Desktop/plots/Heatmap_Top_", top_n, "_Genes_", tumor, "_AllGeneset.pdf")) # Save directly as PDF
}
library(pheatmap)
library(dplyr)

# Define tumor sample names
tumor_samples <- c("Anterior", "Lateral", "Medial", "Posterior", "Intra_Tumor")

# Set the number of top genes you want to visualize
top_n <- 20

# Define output directory
output_dir <- "~/Desktop/plots/"
dir.create(output_dir, showWarnings = FALSE)  # Create the directory if it doesn't exist

# Loop through each tumor sample
for (tumor in tumor_samples) {
    
    # Get the cells belonging to this specific tumor sample
    tumor_cells <- metadata$cell[metadata$original.ident == tumor]
    
    # Extract ssGSEA scores for the tumor sample (for those specific cells)
    tumor_ssgsea_scores <- ssgsea_es[, tumor_cells]
    
    # Calculate the mean ssGSEA score for each gene across all cells in this tumor sample
    gene_avg_scores <- rowMeans(tumor_ssgsea_scores, na.rm = TRUE)
    
    # Sort the genes by mean ssGSEA score in descending order
    sorted_genes <- sort(gene_avg_scores, decreasing = TRUE)
    
    # Get the top N enriched genes
    top_enriched_genes <- names(sorted_genes)[1:top_n]
    
    # Check if the top genes are in the ssGSEA scores matrix
    top_enriched_genes <- top_enriched_genes[top_enriched_genes %in% rownames(tumor_ssgsea_scores)]
    
    # Extract the ssGSEA scores for the top genes
    top_genes_ssgsea_scores <- tumor_ssgsea_scores[top_enriched_genes, , drop = FALSE]
    
    # Create a heatmap for this tumor sample
    pheatmap(top_genes_ssgsea_scores, 
             cluster_rows = TRUE, 
             cluster_cols = TRUE, 
             scale = "row", 
             show_rownames = TRUE, 
             show_colnames = FALSE, 
             main = paste("ssGSEA Scores - Top", top_n, "Genes in", tumor),
             color = colorRampPalette(c("blue", "white", "red"))(100),
             fontsize = 10)
    
    # Optionally save the heatmap to a file
    ggsave(filename = paste0(output_dir, "Heatmap_Top_", top_n, "_Genes_", tumor, ".pdf"), 
           width = 8, height = 6)
}

> 
######################################################################################################
#######################
#######################

# Enrichment map

############################

BiocManager::install("scmap")
library(scmap)
# Optionally, normalize the expression matrix (scmap works best with normalized data)
expr_matrix <- NormalizeData(seurat)
expr_matrix <- seurat@assays$RNA@data
expr_matrix <- as.data.frame(expr_matrix)

library(SingleCellExperiment)
sce <- as.SingleCellExperiment(seurat)
rowData(sce)$feature_symbol <- rownames(sce)
# Remove duplicated gene names (required by scmap)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
sce <- selectFea# Make sure you have cluster identities in colData
sce$cluster <- seurat@active.ident  # If not already present
sce <- indexCluster(sce, cluster_col = "cluster")
tures(sce, suppress_plot = TRUE)

sce <- scmap::selectFeatures(sce, suppress_plot = TRUE)
sce <- indexCluster(sce, cluster_col = "seurat_clusters")
sce <- scmapCluster(projection = sce, index_list = list(ref = metadata(sce)$scmap_cluster_index))
head(colData(sce)$scmap_cluster_labs)



######################################################################
#PIPseeker pipeline
# performed pipseeker pipeline on these samples initial for quality meterics.
/home/budankm/pipseeker-v3.3.0-linux/pipseeker full --skip-version-check \
--fastq /home/budankm/isilon/XieLab_NGSData/PIP_9_12_2024/Demultiplexed/D_S4_L008 \
--star-index-path /home/budankm/isilon/Meghana/starindexGRCh38 \
--chemistry V --output-path /home/budankm/isilon/Meghana/Dsample

/home/budankm/pipseeker-v3.3.0-linux/pipseeker full --skip-version-check \
--fastq /home/budankm/isilon/XieLab_NGSData/PIP_9_12_2024/Demultiplexed/E_S5_L008 \
--star-index-path /home/budankm/isilon/Meghana/starindexGRCh38 \
--chemistry V --output-path /home/budankm/isilon/Meghana/Esample

/home/budankm/pipseeker-v3.3.0-linux/pipseeker full --skip-version-check \
--fastq /home/budankm/isilon/XieLab_NGSData/PIP_9_12_2024/Demultiplexed/F_S6_L008 \
--star-index-path /home/budankm/isilon/Meghana/starindexGRCh38 \
--chemistry V --output-path /home/budankm/isilon/Meghana/Fsample

/home/budankm/pipseeker-v3.3.0-linux/pipseeker full --skip-version-check \
--fastq /home/budankm/isilon/XieLab_NGSData/PIP_9_12_2024/Demultiplexed/G_S7_L008 \
--star-index-path /home/budankm/isilon/Meghana/starindexGRCh38 \
--chemistry V --output-path /home/budankm/isilon/Meghana/Gsample

/home/budankm/pipseeker-v3.3.0-linux/pipseeker full --skip-version-check \
--fastq /home/budankm/isilon/XieLab_NGSData/PIP_9_12_2024/Demultiplexed/H_S8_L008 \
--star-index-path /home/budankm/isilon/Meghana/starindexGRCh38 \
--chemistry V --output-path /home/budankm/isilon/Meghana/Hsample


