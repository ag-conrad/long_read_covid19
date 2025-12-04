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

# ---------------------------------------------------------------------------------------------------------
# 7. Plot transcript statistics
# ---------------------------------------------------------------------------------------------------------



# --- Transcript biotype piechart (Figure 1d) ---

transcript_type = read.table('mart_export.txt', sep = '\t', header = TRUE)
transcript_type$transcript_type_collapsed = ifelse(startsWith(transcript_type$Transcript.type, "IG"), "IG", 
                                                  ifelse(startsWith(transcript_type$Transcript.type, "TR"), "TR", transcript_type$Transcript.type))

DefaultAssay(lr_NS) = "Transcript"
transcript_type =  transcript_type |> 
    subset(Transcript.stable.ID %in% rownames(lr_NS))

transcript_type = transcript_type |> 
    group_by(Transcript.type) |> 
    summarise(count = n()) |> 
    mutate(pct = count/sum(count)) |> arrange(desc(pct))

pie_transtype = ggplot(transcript_type |> head(8) , aes(x="", y = pct, fill=Transcript.type)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_brewer(palette="Set2") +
  ggtitle('Isoform biotypes')

ggsave(file='../plots_manuscript/transcript_biotypes.svg', plot=abundance_line_plt, width=12, height=4)

# --- Barplot- Transcripts per gene (Figure 1e) ---

lr_NS$dummy = 'all'

aggExp = AggregateExpression(lr_NS, group.by = 'dummy')
aggExp = aggExp$Transcript |> data.frame()
aggExp = merge(aggExp, dict[c('gene_name', 'transcript_id')], by.x = 0, by.y = 'transcript_id') |> 
                    rename('transcript_id' = 'Row.names')

nrtranscripts_all = aggExp |> 
    filter(all != 0) |>                 
    group_by(gene_name) |>                    
    summarise(transcript_count = n_distinct(transcript_id))

nrtranscripts_all$transcript_count_collapsed = ifelse(nrtranscripts_all$transcript_count > 5, ">5", as.character(nrtranscripts_all$transcript_count))
nrtranscripts_all$transcript_count_collapsed = factor(nrtranscripts_all$transcript_count_collapsed, levels = c('1', '2', '3', '4', '5', '>5'))

nrtranscripts_all = ggplot(nrtranscripts_all) + 
    geom_bar(aes(x=transcript_count_collapsed, fill = transcript_count)) + 
    scale_fill_manual(values = rev(brewer.pal(9, "Greens"))) + 
    xlab("Isoforms per gene") + ylab('Number of genes') +
    theme_minimal() +
    theme(legend.position = 'none', 
         axis.text = element_text(size=12), 
         axis.title = element_text(size=14)) +
   ggtitle("Number of isoforms per gene")

ggsave(file='../plots_manuscript/transcripts_per_gene.svg', plot=abundance_line_plt, width=12, height=4)


# --- Lineplot- Ranked transcript expression (Figure 1f) ---

transranks_all = aggExp |> 
  group_by(gene_name) |> 
  arrange(gene_name, desc(all)) |> 
  mutate(rank = row_number(), transcript_count = n_distinct(transcript_id)) |> 
  mutate(transcript_pct = all/ sum(all)) |> 
  group_by(gene_name) |> 
  mutate(
    pct_rank1 = transcript_pct[rank == "1"]
  ) |> 
  ungroup()

abundance_line_plt = ggplot(transranks_all %>% subset(rank < 6)) + 
    geom_line(aes(x=rank, y = transcript_pct, group = gene_name, colour = pct_rank1), alpha = 0.05, linewidth = 0.1) + 
    xlab("Isoform rank") + ylab('Expression percentage') +
    theme_minimal() +
    scale_colour_viridis(
        option = "D"
    ) +
    theme(
         axis.text = element_text(size=12), 
         axis.title = element_text(size=14)) + 
   ggtitle("Isoform expression per gene (n = 30,932)")

ggsave(file='../plots_manuscript/iso_abund_lines.svg', plot=abundance_line_plt, width=12, height=4)


