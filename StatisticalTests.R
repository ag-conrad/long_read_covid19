# =========================================================================================================
# Script:     Cell Distribution Test and PTPRC Transcript Analysis 
# Author:     Kristin Köhler
# Date:       2025-10-14
# Description:
#   This script evaluates cell-type–specific PTPRC transcript expression within immune cells
#   using long-read single-cell RNA-seq data. It performs Fisher’s exact test for global
#   cell type–condition associations and transcript-level marker testing in immune subsets.
# Contact:
#   kristin.koehler@bih-charite.de
# =========================================================================================================

# ---------------------------------------------------------------------------------------------------------
# 1. Load required libraries
# ---------------------------------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)        # includes dplyr, ggplot2, stringr, etc.
library(patchwork)
library(viridis)
library(scCustomize)
library(ggrepel)
library(RColorBrewer)
library(ggtranscript)
library(ComplexHeatmap)
library(spatstat.explore)
library(gridExtra)
library(grid)

# ---------------------------------------------------------------------------------------------------------
# 2. Load data
# ---------------------------------------------------------------------------------------------------------
dict  <- read.table("transcript2gene.tsv", sep = " ", header = TRUE)
lr_NS <- readRDS("Seurat_FinalLongRead.rds")

# ---------------------------------------------------------------------------------------------------------
# 3. Fisher’s exact test: Association between cell type and condition
# ---------------------------------------------------------------------------------------------------------
# Create contingency table
im_epi_condi_table <- table(lr_NS$celltype, lr_NS$condition)

# Global Fisher’s test (K × 2 table): tests association between cell type and condition
global_fisher <- fisher.test(im_epi_condi_table)
global_fisher$p.value
global_fisher  # Full test output

# ---------------------------------------------------------------------------------------------------------
# 4. Define PTPRC transcript IDs
# ---------------------------------------------------------------------------------------------------------
ptprc_abc <- "ENST00000442510"
ptprc_a   <- "ENST00000413409"
ptprc_o   <- "ENST00000348564"
ptprc_b   <- "ENST00000530727"

# ---------------------------------------------------------------------------------------------------------
# 5. Subset immune cells
# ---------------------------------------------------------------------------------------------------------
lr_NS_immune <- lr_NS %>% subset(celltype == "immune")

# ---------------------------------------------------------------------------------------------------------
# 6. Differential expression (FindMarkers) for individual PTPRC transcripts
# ---------------------------------------------------------------------------------------------------------
# Each test compares a target immune cell type (e.g., B_cell, T_cell subsets, MoMa)
# against the rest within the immune subset. Using assay = "Transcript" allows
# isoform-level comparisons.

# --- B cells vs others ---
FindMarkers(
  object = lr_NS_immune,
  ident.1 = "B_cell",
  assay = "Transcript",
  slot = "data",
  features = ptprc_abc,
  logfc.threshold = 0,
  min.pct = 0
)

# --- CD4+/CD8+ T cells vs others (PTPRC isoforms) ---
FindMarkers(
  object = lr_NS_immune,
  ident.1 = c("CD4_T", "CD8_T"),
  assay = "Transcript",
  slot = "data",
  features = ptprc_a,
  logfc.threshold = 0,
  min.pct = 0
)

FindMarkers(
  object = lr_NS_immune,
  ident.1 = c("CD4_T", "CD8_T"),
  assay = "Transcript",
  slot = "data",
  features = ptprc_o,
  logfc.threshold = 0,
  min.pct = 0
)

FindMarkers(
  object = lr_NS_immune,
  ident.1 = c("CD4_T", "CD8_T"),
  assay = "Transcript",
  slot = "data",
  features = ptprc_b,
  logfc.threshold = 0,
  min.pct = 0
)

# --- Monocytes/Macrophages (MoMa) ---
FindMarkers(
  object = lr_NS_immune,
  ident.1 = "MoMa",
  assay = "Transcript",
  slot = "data",
  features = ptprc_abc,
  logfc.threshold = 0,
  min.pct = 0
)

FindMarkers(
  object = lr_NS_immune,
  ident.1 = "MoMa",
  assay = "Transcript",
  slot = "data",
  features = ptprc_a,
  logfc.threshold = 0,
  min.pct = 0
)

FindMarkers(
  object = lr_NS_immune,
  ident.1 = "MoMa",
  assay = "Transcript",
  slot = "data",
  features = ptprc_b,
  logfc.threshold = 0,
  min.pct = 0
)

FindMarkers(
  object = lr_NS_immune,
  ident.1 = "MoMa",
  assay = "Transcript",
  slot = "data",
  features = ptprc_o,
  logfc.threshold = 0,
  min.pct = 0
)
