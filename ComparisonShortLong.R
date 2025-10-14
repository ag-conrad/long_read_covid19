# =========================================================================================================
# Script:     Long-read vs. Short-read Seurat Integration and Label Transfer
# Author:     Kristin Köhler
# Date:       2025-10-14
# Description:
#   This script compares single-cell long-read (ONT) and short-read (Illumina) datasets.
#   It matches cell barcodes, transfers cell-type annotations between datasets using Seurat,
#   and visualizes UMAP projections and annotations and label transfer accuracy.
# Contact:
#   kristin.koehler@bih-charite.de
# =========================================================================================================

# Load required packages -------------------------------------------------------------------
library(Seurat)
library(tidyverse)       
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(svglite)

# ---------------------------------------------------------------------------------------------------------
# 1. Load Seurat objects
# ---------------------------------------------------------------------------------------------------------
lr_NS <- readRDS("Seurat_FinalLongRead.rds")
sr_NS <- readRDS("Seurat_FinalShortRead.rds")

# ---------------------------------------------------------------------------------------------------------
# 2. Match barcodes between short- and long-read data
# ---------------------------------------------------------------------------------------------------------
matched_barcodes_sr <- sub("(_|-).*", "", rownames(sr_NS[[]]))
matched_barcodes_sr <- paste0(matched_barcodes_sr, "_", sr_NS$orig.ident.corrected)
sr_NS$matched_barcodes <- matched_barcodes_sr

matched_barcodes_lr <- sub("(_|-).*", "", rownames(lr_NS[[]]))
matched_barcodes_lr <- paste0(matched_barcodes_lr, "_", lr_NS$orig.ident)
lr_NS$matched_barcodes <- matched_barcodes_lr

# Find matching indices between datasets
match_indices_s2l <- match(sr_NS$matched_barcodes, lr_NS$matched_barcodes)
match_indices_l2s <- match(lr_NS$matched_barcodes, sr_NS$matched_barcodes)

# Transfer annotation labels
sr_NS$lr_anno <- lr_NS$finalAnno[match_indices_s2l] %>% unname()
lr_NS$sr_anno <- sr_NS$finalAnnoRenamed[match_indices_l2s] %>% unname()

# ---------------------------------------------------------------------------------------------------------
# 3. Visualize UMAPs with annotations
# ---------------------------------------------------------------------------------------------------------
options(repr.plot.height = 9, repr.plot.width = 20)
DimPlot(lr_NS, reduction = "umap.harmony", group.by = "finalAnno", pt.size = .6, label = TRUE, cols = paired_palette) +
  ggtitle("Long-read-based UMAP & annotation") |
DimPlot(lr_NS, reduction = "umap.harmony", group.by = "sr_anno", pt.size = .6, label = TRUE, cols = paired_palette) +
  ggtitle("Transferred short-read annotation")

# Comparison of both technologies
options(repr.plot.height = 7, repr.plot.width = 17)
(
  DimPlot(lr_NS, reduction = "umap.harmony", group.by = "finalAnno", pt.size = .3, label = TRUE, repel = TRUE, cols = paired_palette) +
    ggtitle("Long-read-based UMAP & annotation") |
  DimPlot(sr_NS, reduction = "umap.harmony.wrapper2", group.by = "finalAnnoRenamed", pt.size = .3, label = TRUE, repel = TRUE, cols = paired_palette) +
    ggtitle("Short-read-based UMAP & annotation")
)

# Gene detection UMAPs
lr_umap_gene <- FeaturePlot(lr_NS, reduction = "umap.harmony", features = "nFeature_RNA", pt.size = .4, label = TRUE, repel = TRUE) +
  ggtitle("Detected genes (long-read UMAP)")
sr_umap_gene <- FeaturePlot(sr_NS, reduction = "umap.harmony.wrapper2", features = "nFeature_RNA", pt.size = .4, label = TRUE, repel = TRUE) +
  ggtitle("Detected genes (short-read UMAP)")

# Short-read UMAP with projected long-read annotation
sr_umap_proj <- DimPlot(sr_NS, reduction = "umap.harmony.wrapper2", group.by = "lr_anno", pt.size = .3, label = TRUE, repel = TRUE, cols = paired_palette) +
  ggtitle("Transferred long-read annotation") +
  theme(legend.position = "None")

# ---------------------------------------------------------------------------------------------------------
# 4. Label transfer between datasets
# ---------------------------------------------------------------------------------------------------------
lr_NS.anchors <- FindTransferAnchors(
  reference = sr_NS,
  query = lr_NS,
  dims = 1:30,
  reference.reduction = "pca"
)
predictions <- TransferData(
  anchorset = lr_NS.anchors,
  refdata = sr_NS$finalAnnoRenamed,
  dims = 1:30
)

# Add transfer scores to metadata
lr_NS_meta <- lr_NS
lr_NS_meta@meta.data <- lr_NS_meta@meta.data %>% select(!starts_with("prediction.score"))
lr_NS_meta <- AddMetaData(lr_NS_meta, metadata = predictions)

# Prepare transfer score matrix
transfer_score_df <- lr_NS_meta@meta.data %>%
  select(predicted.id, finalAnno, starts_with("prediction.score"))
heat_df <- transfer_score_df %>%
  group_by(finalAnno) %>%
  summarise(across(starts_with("prediction"), mean, na.rm = TRUE)) %>%
  rename(lrAnno = finalAnno)
colnames(heat_df) <- colnames(heat_df) %>% str_replace("prediction.score.", "")
heat_df <- heat_df %>% select(!max)

# Define annotation levels
ct_levels_sr <- c(
  "Neutrophil", "MoMa", "rMa", "cDCs", "Mast", "B",
  "CD4_T", "CD8_T", "T_lowQual",
  "Basal", "Secretory", "IFN_high", "CiliatedDiff_Stress",
  "CiliatedDiff", "Ciliated", "FOXN4", "Ionocyte", "Mitotic", "Squamous"
)
ct_levels_lr <- c(
  "Neutrophil", "Neutrophil_Outlier", "MoMa", "Mast_cell", "B_cell",
  "CD4_T", "CD8_T", "Basal", "Secretory", "CiliatedDiff",
  "Ciliated", "FOXN4", "Ionocyte", "Mitotic", "Squamous"
)

# Melt for plotting
heat_df_long <- melt(heat_df, id.vars = "lrAnno") %>% rename(srAnno = variable)
heat_df_long$lrAnno <- factor(heat_df_long$lrAnno, levels = ct_levels_lr)
heat_df_long$srAnno <- factor(heat_df_long$srAnno, levels = ct_levels_sr)

# ---------------------------------------------------------------------------------------------------------
# 5. Heatmap of label transfer accuracy
# ---------------------------------------------------------------------------------------------------------
palette <- colorRampPalette(brewer.pal(9, "OrRd"))(100)

labeltrans_hm <- ggplot(heat_df_long, aes(x = lrAnno, y = srAnno, fill = value)) +
  geom_tile(color = "#999999", lwd = 0.5) +
  scale_fill_gradientn(colors = palette, name = "Transfer score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 15, margin = margin(t = 10)),
    axis.title.y = element_text(size = 15, margin = margin(r = 10)),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank()
  ) +
  xlab("Long-read annotation") +
  ylab("Short-read annotation")

# ---------------------------------------------------------------------------------------------------------
# 6. Combined UMAP & label transfer visualization
# ---------------------------------------------------------------------------------------------------------
options(repr.plot.height = 7, repr.plot.width = 17)
lr_sr_UMAP <- (
  DimPlot(lr_NS, reduction = "umap.harmony.adapted", group.by = "finalAnno", pt.size = .1, label = TRUE, repel = TRUE, cols = paired_palette) +
    ggtitle("Long-read-based UMAP & annotation") |
  DimPlot(sr_NS, reduction = "umap.harmony.wrapper2", group.by = "finalAnnoRenamed", pt.size = .1, label = TRUE, repel = TRUE, cols = paired_palette) +
    ggtitle("Short-read-based UMAP & annotation")
)

options(repr.plot.height = 8, repr.plot.width = 17)
UMAPproj_labeltrans <- (
  DimPlot(lr_NS, reduction = "umap.harmony.adapted", group.by = "sr_anno", pt.size = .1, label = TRUE, repel = TRUE, cols = paired_palette) +
    ggtitle("Transferred short-read annotation") | labeltrans_hm
) + plot_layout(widths = c(1, 1))

# ---------------------------------------------------------------------------------------------------------
# 7. Save figures
# ---------------------------------------------------------------------------------------------------------
plot_path <- "path"
ggsave(paste0(plot_path, "comp_lr_sr.svg"), plot = UMAPproj_labeltrans, width = 17, height = 8)
ggsave(paste0(plot_path, "lr_sr_UMAP.svg"), plot = lr_sr_UMAP, width = 17, height = 7)

combi_plt <- (lr_umap_gene | sr_umap_gene) / (sr_umap_proj | overlap_bcs)
ggsave(paste0(plot_path, "projections.svg"), plot = combi_plt, width = 17, height = 14)
