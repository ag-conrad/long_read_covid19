# =========================================================================================================
# Script:     DEG and DTU Summary Statistics and Pathway Enrichment
# Author:     Kristin Köhler
# Date:       2025-10-14
# Description:
#   This script summarizes and visualizes the number of differentially expressed genes (DEGs)
#   and differentially used transcripts (DTUs) across cell types. It additionally performs pathway 
#   enrichment analysis (GSEA) for differentially expressed (DEG)
#   and differentially used transcript (DTU) genes using Reactome, GO, and Hallmark gene sets.
#   Generates per-cell-type enrichment results and publication-ready heatmaps and plots.
#
# Contact:
#   kristin.koehler@bih-charite.de
# =========================================================================================================

# ---------------------------------------------------------------------------------------------------------
# 1. Load required libraries
# ---------------------------------------------------------------------------------------------------------
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(reactome.db)
library(ReactomePA)
library(org.Hs.eg.db)
library(enrichplot)
library(pheatmap)
library(RColorBrewer)


# ---------------------------------------------------------------------------------------------------------
# 2. Load DEG and DTU test statistics
# ---------------------------------------------------------------------------------------------------------

degs_covcon <- read.table("DEG_combined_CovVsCon.tsv", header = TRUE)
dtus_covcon <- read.table("DTU_combined_CovVsCon.tsv", header = TRUE)

degs_modcrit <- read.table("DEG_combined_ModVsCrit.tsv", header = TRUE)
dtus_modcrit <- read.table("DTU_combined_ModVsCrit.tsv", header = TRUE)

# =========================================================================================================
# 3. DEG and DTU numbers — Controls vs. COVID-19
# =========================================================================================================

# Summarize DEGs per cell type
test_stats_covcon <- degs_covcon %>%
  group_by(Celltype) %>%
  summarise(
    total = n(),
    significant = sum(padj < 0.05, na.rm = TRUE)
  ) %>%
  mutate(type = "DEGs")

# Summarize DTUs per cell type
test_stats_covcon <- rbind(
  test_stats_covcon,
  dtus_covcon %>%
    select(Celltype, groupID, padj) %>%
    unique() %>%
    group_by(groupID) %>%
    slice_min(order_by = padj, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(Celltype) %>%
    summarise(
      total = n(),
      significant = sum(padj < 0.05, na.rm = TRUE)
    ) %>%
    mutate(type = "DTUs")
)

# Compute DEG–DTU overlaps per cell type
overlaps <- list()
for (ct in unique(degs_covcon$Celltype)) {
  overlaps[[ct]] <- intersect(
    degs_covcon %>% subset(Celltype == ct, padj < 0.05) %>% pull(gene_name) %>% unique(),
    dtus_covcon %>% subset(Celltype == ct, padj < 0.05) %>% pull(gene_name) %>% unique()
  ) %>% length()
}

# Add overlap statistics
test_stats_covcon <- rbind(
  test_stats_covcon,
  data.frame(
    Celltype = names(overlaps),
    total = 0,
    significant = unname(unlist(overlaps)),
    type = "Overlap"
  )
)

# Define consistent cell type order
test_stats_covcon$Celltype <- factor(
  test_stats_covcon$Celltype,
  levels = c(
    "Secretory", "Ciliated", "BasalDiff", "CiliatedDiff", "Ionocytes",
    "Squamous", "FOXN4", "MoMa", "T_cell", "B_cell",
    "Plasma_B", "Neutrophil", "Mitotic"
  )
)

# =========================================================================================================
# 4. DEG and DTU numbers — Moderate vs. Critical COVID-19
# =========================================================================================================

# Summarize DEGs per cell type
test_stats_modcrit <- degs_modcrit %>%
  group_by(Celltype) %>%
  summarise(
    total = n(),
    significant = sum(padj < 0.05, na.rm = TRUE)
  ) %>%
  mutate(type = "DEGs")

# Summarize DTUs per cell type
test_stats_modcrit <- rbind(
  test_stats_modcrit,
  dtus_modcrit %>%
    select(Celltype, groupID, padj) %>%
    unique() %>%
    group_by(groupID) %>%
    slice_min(order_by = padj, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(Celltype) %>%
    summarise(
      total = n(),
      significant = sum(padj < 0.05, na.rm = TRUE)
    ) %>%
    mutate(type = "DTUs")
)

# Compute DEG–DTU overlaps per cell type
overlaps <- list()
for (ct in unique(degs_modcrit$Celltype)) {
  overlaps[[ct]] <- intersect(
    degs_modcrit %>% subset(Celltype == ct, padj < 0.05) %>% pull(gene_name) %>% unique(),
    dtus_modcrit %>% subset(Celltype == ct, padj < 0.05) %>% pull(gene_name) %>% unique()
  ) %>% length()
}

# Add overlap statistics
test_stats_modcrit <- rbind(
  test_stats_modcrit,
  data.frame(
    Celltype = names(overlaps),
    total = 0,
    significant = unname(unlist(overlaps)),
    type = "Overlap"
  )
)

# Define consistent cell type order
test_stats_modcrit$Celltype <- factor(
  test_stats_modcrit$Celltype,
  levels = c(
    "Secretory", "Ciliated", "BasalDiff", "CiliatedDiff", "Ionocytes",
    "Squamous", "FOXN4", "MoMa", "T_cell", "B_cell",
    "Plasma_B", "Neutrophil", "Neutrophil_Outlier", "Mitotic"
  )
)

# =========================================================================================================
# 5. Visualization — DEG and DTU summary plots (Figure 3a–b)
# =========================================================================================================

# --- Control vs. COVID-19 ---
nr_covcon_plt <- ggplot(test_stats_covcon) +
  geom_bar(aes(x = Celltype, y = significant, fill = type), stat = "identity") +
  theme_minimal() +
  scale_fill_brewer(palette = "Paired", name = "") +
  scale_y_sqrt(breaks = c(10, 100, 250, 500, 750)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "None"
  ) +
  xlab("") + ylab("Number of genes") +
  coord_flip()

# --- Moderate vs. Critical COVID-19 ---
nr_modcrit_plt <- ggplot(test_stats_modcrit) +
  geom_bar(aes(x = Celltype, y = significant, fill = type), stat = "identity") +
  theme_minimal() +
  scale_fill_brewer(palette = "Paired", name = "") +
  scale_y_sqrt(breaks = c(10, 50, 100)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("") + ylab("Number of genes") +
  coord_flip()

# Save figures
ggsave("COVID_manuscript/DEG_DTU_CtrlCov.svg", plot = nr_covcon_plt, width = 9, height = 3)
ggsave("COVID_manuscript/DEG_DTU_ModCrit.svg", plot = nr_modcrit_plt, width = 9, height = 3)



# =========================================================================================================
# 6. Utility functions
# =========================================================================================================

# ---------------------------------------------------------------------------------------------------------
# 6.1 Convert Ensembl IDs to Entrez IDs
# ---------------------------------------------------------------------------------------------------------
translate2entrez <- function(ensembl_id) {
  entrez_id <- mapIds(
    org.Hs.eg.db,
    keys = ensembl_id,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  data.frame(ensembl_gene_id = ensembl_id, entrezgene_id = entrez_id)
}

# ---------------------------------------------------------------------------------------------------------
# 6.2 Reactome GSEA (Gene Set Enrichment Analysis)
# ---------------------------------------------------------------------------------------------------------
reactomeGSEA <- function(de_df, celltype, p_cutoff) {
  id_conversion <- translate2entrez(de_df$Gene)
  de_df_merged <- merge(de_df, id_conversion, by.x = "Gene", by.y = "ensembl_gene_id")

  ranked_genes <- if (celltype == "all") {
    de_df_merged %>%
      select(entrezgene_id, stat) %>%
      na.omit() %>%
      group_by(entrezgene_id) %>%
      slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
      arrange(-stat) %>%
      distinct() %>%
      deframe()
  } else {
    de_df_merged %>%
      subset(Celltype == celltype) %>%
      select(entrezgene_id, stat) %>%
      na.omit() %>%
      group_by(entrezgene_id) %>%
      slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
      arrange(-stat) %>%
      distinct() %>%
      deframe()
  }

  gsePathway(
    gene = ranked_genes,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = p_cutoff,
    nPermSimple = 100000
  )
}

# ---------------------------------------------------------------------------------------------------------
# 6.3 General pathway GSEA (GO / Hallmark)
# ---------------------------------------------------------------------------------------------------------
pathwayGSEA <- function(de_df, pathway_set, celltype, p_cutoff) {
  ranked_genes <- if (celltype == "all") {
    de_df %>%
      select(gene_name, stat) %>%
      arrange(-stat) %>%
      na.omit() %>%
      distinct() %>%
      deframe()
  } else {
    de_df %>%
      subset(Celltype == celltype) %>%
      select(gene_name, stat) %>%
      arrange(-stat) %>%
      na.omit() %>%
      distinct() %>%
      deframe()
  }

  GSEA(
    geneList = ranked_genes,
    TERM2GENE = pathway_set,
    pvalueCutoff = p_cutoff,
    nPermSimple = 10000,
    eps = 1e-20
  )
}

# ---------------------------------------------------------------------------------------------------------
# 6.4 Combined GSEA plot showing DEG ranks and DTU overlaps
# ---------------------------------------------------------------------------------------------------------
combi_gseaplot <- function(gsea_res, dtus, pathway_set, celltype, pathway_id, reactome) {
  plt <- gseaplot2(gsea_res, geneSetID = pathway_id, title = gsea_res$Description[pathway_id])
  plt1 <- plt[[1]] + ylab("Running ES") + theme(axis.title.y = element_text(size = 6))
  plt2 <- plt[[2]] + ylab("Pathway\ngenes") + theme(axis.title.y = element_text(size = 6))
  plt3 <- plt[[3]] + ylab("Ranked\nstatistic") + theme(axis.title.y = element_text(size = 6))

  dtus_sign <- dtus %>%
    subset(padj < 0.1, Celltype == celltype) %>%
    pull(gene_name) %>%
    unique()

  if (reactome) {
    pathway_genes <- pathway_set[[gsea_res$ID[pathway_id]]]
    id_conversion <- bitr(
      pathway_genes, toType = "SYMBOL",
      fromType = "ENTREZID", OrgDb = org.Hs.eg.db
    )
    dtus_sign <- id_conversion %>%
      subset(SYMBOL %in% dtus_sign) %>%
      pull(ENTREZID)
  } else {
    pathway_genes <- pathway_set %>%
      subset(term == gsea_res$Description[pathway_id]) %>%
      pull(gene)
  }

  dtus_sign <- dtus_sign[dtus_sign %in% pathway_genes]
  dtus_sign_ranks <- match(dtus_sign, names(gsea_res@geneList)) %>% na.omit()

  dtu_segm_plt <- ggplot() +
    geom_segment(data = data.frame(x = dtus_sign_ranks),
                 aes(x = x, xend = x, y = -0.1, yend = -1),
                 color = "black", size = 0.5) +
    ylab("DTU\ngenes") +
    scale_x_continuous(limits = c(1, length(gsea_res@geneList)), expand = c(0, 0)) +
    theme(
      axis.title.y = element_text(size = 6),
      axis.title.x = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA)
    )

  (plt1 / plt2 / dtu_segm_plt / plt3) +
    plot_annotation(title = ifelse(reactome, paste0("Reactome GSEA - ", celltype),
                                   paste0("MolSigDB GSEA - ", celltype))) +
    plot_layout(heights = c(2, 0.5, 0.5, 1.5))
}

# ---------------------------------------------------------------------------------------------------------
# 6.5 Extract NES and adjusted p-value matrices from GSEA results
# ---------------------------------------------------------------------------------------------------------
getEnrichmentMatrices <- function(gsea_result_list, celltypes, top_n_pathways, color.by) {
  gsea_stats <- data.frame()
  for (celltype in celltypes) {
    message(celltype)
    if (length(gsea_result_list[[celltype]]@result$ID) > 0) {
      gsea_stats <- rbind(
        gsea_stats,
        data.frame(
          PathwayName = gsea_result_list[[celltype]]@result$Description,
          PathwayID = gsea_result_list[[celltype]]@result$ID,
          NES = gsea_result_list[[celltype]]@result$NES,
          Padj = gsea_result_list[[celltype]]@result$p.adjust,
          LogDirPadj = ifelse(gsea_result_list[[celltype]]@result$NES > 0,
                              -log10(gsea_result_list[[celltype]]@result$p.adjust),
                              log10(gsea_result_list[[celltype]]@result$p.adjust)),
          Celltype = celltype
        )
      )
    }
  }

  pathways <- unlist(lapply(celltypes, function(ct) {
    head(gsea_result_list[[ct]]@result$ID, top_n_pathways)
  }))

  gsea_stats_filtered <- gsea_stats %>%
    subset(PathwayID %in% unique(pathways))

  NES_matrix <- gsea_stats_filtered %>%
    distinct(Celltype, PathwayName, .keep_all = TRUE) %>%
    pivot_wider(id_cols = PathwayName,
                names_from = Celltype,
                values_from = !!sym(color.by),
                values_fill = NA) %>%
    column_to_rownames("PathwayName") %>%
    as.matrix() %>%
    .[rowSums(is.na(.)) == 0, ]

  sign_matrix <- gsea_stats_filtered %>%
    distinct(Celltype, PathwayName, .keep_all = TRUE) %>%
    pivot_wider(id_cols = PathwayName,
                names_from = Celltype,
                values_from = Padj,
                values_fill = NA) %>%
    column_to_rownames("PathwayName") %>%
    as.matrix() %>%
    .[rowSums(is.na(.)) == 0, ]

  sign_matrix <- ifelse(sign_matrix < 0.05, "*", "")
  list(NES_matrix, sign_matrix)
}

# ---------------------------------------------------------------------------------------------------------
# 6.6 Plot enrichment heatmaps
# ---------------------------------------------------------------------------------------------------------
plotEnrichmentHeat <- function(NES_matrix, sign_matrix, cluster_cols, cluster_rows, title) {
  pheatmap(NES_matrix,
           display_numbers = sign_matrix,
           scale = "none",
           cluster_cols = cluster_cols,
           cluster_rows = cluster_rows,
           show_rownames = TRUE,
           show_colnames = TRUE,
           na_col = "darkgrey",
           treeheight_row = 0,
           treeheight_col = 0,
           main = title)
}

# =========================================================================================================
# 7. Pathway enrichment on DEGs and DTUs
# =========================================================================================================

hallmark.pathways <- read.gmt("resources/h.all.v2024.1.Hs.symbols.gmt.txt")
go.pathways <- read.gmt("resources/c5.go.bp.v2024.1.Hs.symbols.gmt.txt")
reactome.pathways <- as.list(reactomePATHID2EXTID)

plots_reactome_ct <- list()
plots_go_ct <- list()
plots_hallmark_ct <- list()

celltypes <- unique(degs_modcrit$Celltype)
results_modcrit <- list()
results_covcon <- list()

# Reactome GSEA for both comparisons
for (celltype in unique(degs_covcon$Celltype)) {
  message(celltype)
  results_covcon[[celltype]] <- reactomeGSEA(degs_covcon, celltype, 1)
}

for (celltype in unique(degs_modcrit$Celltype)) {
  message(celltype)
  results_modcrit[[celltype]] <- reactomeGSEA(degs_modcrit, celltype, 1)
}

# =========================================================================================================
# 8. Summary heatmaps (Figure 3c, 3d)
# =========================================================================================================
ct_covcon <- c("Mitotic", "Ionocytes", "FOXN4", "Ciliated", "CiliatedDiff", "Secretory",
               "BasalDiff", "MoMa", "Neutrophil", "Plasma_B", "B_cell", "T_cell")

matrices_covcon <- getEnrichmentMatrices(results_covcon, ct_covcon, 8, "NES")
nes_mat_covcon <- matrices_covcon[[1]]
sign_mat_covcon <- matrices_covcon[[2]]

pathways_fil <- rownames(nes_mat_covcon)
heatmap_covcon <- plotEnrichmentHeat(nes_mat_covcon[pathways_fil, ct_covcon],
                                     sign_mat_covcon[pathways_fil, ct_covcon],
                                     FALSE, FALSE,
                                     "Reactome GSEA on DEGs - Control vs. COVID-19")
ggsave("heatmap_covcon.svg", plot = heatmap_covcon, width = 8, height = 8)

ct_modcrit <- c("Mitotic", "Ionocytes", "FOXN4", "Ciliated", "CiliatedDiff", "Secretory",
                "BasalDiff", "MoMa", "Neutrophil", "Plasma_B", "B_cell", "T_cell")

matrices_mod <- getEnrichmentMatrices(results_modcrit[names(results_modcrit) %in% ct_modcrit],
                                      ct_modcrit, 8, "NES")
nes_mat_modcrit <- matrices_mod[[1]]
sign_mat_modcrit <- matrices_mod[[2]]

heat_modcrit <- plotEnrichmentHeat(-nes_mat_modcrit, sign_mat_modcrit,
                                   FALSE, FALSE,
                                   "Reactome DEG GSEA - Moderate vs. Severe COVID-19")
ggsave("heatmap_modcrit.svg", plot = heat_modcrit, width = 8.5, height = 8)

# =========================================================================================================
# 9. GSEA plots (Figure 4)
# =========================================================================================================
for (celltype in celltypes) {
  message(celltype)

  reactome_res <- reactomeGSEA(degs_covcon, celltype, 1)
  hallmark_res <- pathwayGSEA(degs_covcon, hallmark.pathways, celltype, 1)
  go_res <- pathwayGSEA(degs_covcon, go.pathways, celltype, 1)

  plots_reactome <- plots_go <- plots_hallmark <- list()

  if (length(hallmark_res@result$ID) != 0) {
    for (i in seq_len(min(9, length(hallmark_res@result$ID)))) {
      plots_hallmark[[i]] <- combi_gseaplot(hallmark_res, dtus_nof, hallmark.pathways, celltype, i, FALSE)
    }
  }

  if (length(reactome_res@result$ID) != 0) {
    for (i in seq_len(min(9, length(reactome_res@result$ID)))) {
      plots_reactome[[i]] <- combi_gseaplot(reactome_res, dtus_nof, reactome.pathways, celltype, i, TRUE)
    }
  }

  if (length(go_res@result$ID) != 0) {
    for (i in seq_len(min(9, length(go_res@result$ID)))) {
      plots_go[[i]] <- combi_gseaplot(go_res, dtus_nof, go.pathways, celltype, i, FALSE)
    }
  }

  plots_reactome_ct[[celltype]] <- wrap_plots(plots_reactome, ncol = 3, nrow = 3) +
    plot_annotation(title = paste0("Reactome GSEA - ", celltype))
  plots_go_ct[[celltype]] <- wrap_plots(plots_go, ncol = 3, nrow = 3) +
    plot_annotation(title = paste0("GO GSEA - ", celltype))
  plots_hallmark_ct[[celltype]] <- wrap_plots(plots_hallmark, ncol = 3, nrow = 3) +
    plot_annotation(title = paste0("Hallmark GSEA - ", celltype))
}

pdf("ReactomePathways_Enrichment_top9_scaled.pdf", width = 12, height = 10)
for (ct in names(plots_reactome_ct)) {
  if (length(plots_reactome_ct[[ct]]) != 0) print(plots_reactome_ct[[ct]])
}
dev.off()

pdf("GOPathways_Enrichment_top9_scaled.pdf", width = 12, height = 10)
for (ct in names(plots_go_ct)) {
  if (length(plots_go_ct[[ct]]) != 0) print(plots_go_ct[[ct]])
}
dev.off()

pdf("HallmarkPathways_Enrichment_top9_scaled.pdf", width = 12, height = 10)
for (ct in names(plots_hallmark_ct)) {
  if (length(plots_hallmark_ct[[ct]]) != 0) print(plots_hallmark_ct[[ct]])
}
dev.off()
