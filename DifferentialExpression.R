# =========================================================================================================
# Script:     Differential Expression (DEG) and Differential Transcript Usage (DTU) Analysis
# Author:     Kristin Köhler
# Date:       2025-10-14
# Description:
#   This script performs DEG and DTU analyses per cell type using DESeq2 and DEXSeq on Seurat-derived
#   count matrices. It aggregates counts, applies filtering, runs gene- and transcript-level tests,
#   and writes per-cell-type results. 
# Contact:
#   kristin.koehler@bih-charite.de
# =========================================================================================================

# ---------------------------------------------------------------------------------------------------------
# 1. Load required libraries
# ---------------------------------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)     
library(DESeq2)
library(DEXSeq)
library(DRIMSeq)

# =========================================================================================================
# 2. Utility functions
# =========================================================================================================

# ---------------------------------------------------------------------------------------------------------
# 2.1 Aggregate expression counts per sample for one cell type
# ---------------------------------------------------------------------------------------------------------
aggregateExpressionMatrix <- function(seurat_obj, cellType, trans2gene, sampleData) {
  # Subset by cell type
  cellTypeSubset <- seurat_obj %>% subset(celltype == cellType)

  # Aggregate counts per sample (gene- and transcript-level)
  counts_geneExp  <- AggregateExpression(cellTypeSubset, group.by = "orig.ident", assay = "RNA")$RNA
  counts_transExp <- AggregateExpression(cellTypeSubset, group.by = "orig.ident", assay = "Transcript")$Transcript

  # Ensure all samples are present in subset
  samples <- colnames(counts_geneExp)[colnames(counts_geneExp) %in% sampleData$sample_id]
  sampleData <- sampleData[samples, ]
  sampleDataCheck <- sampleData %>% group_by(condition) %>% mutate(groupSize = n()) %>% ungroup()
  if (min(sampleDataCheck$groupSize) < 2) stop("Too few samples per condition")

  # Reorder columns and merge with transcript–gene dictionary
  counts_geneExp <- counts_geneExp[, rownames(sampleData)]
  counts_transExp <- counts_transExp[, rownames(sampleData)]
  counts_transExp <- merge(
    trans2gene[c("transcript_id", "gene_id")],
    counts_transExp, by.y = 0, by.x = "transcript_id",
    all.x = FALSE, all.y = TRUE
  ) %>% rename(feature_id = transcript_id)
  rownames(counts_geneExp) <- dict[match(rownames(counts_geneExp), dict$gene_name), "gene_id"]

  return(list(geneExpDf = counts_geneExp, transExpDf = counts_transExp, sampleDat = sampleData))
}

# ---------------------------------------------------------------------------------------------------------
# 2.2 Aggregate across all cells (bulk-like)
# ---------------------------------------------------------------------------------------------------------
aggregateExpressionToBulk <- function(seurat_obj, trans2gene, sampleData) {
  counts_geneExp  <- AggregateExpression(seurat_obj, group.by = "orig.ident", assay = "RNA")$RNA
  counts_transExp <- AggregateExpression(seurat_obj, group.by = "orig.ident", assay = "Transcript")$Transcript

  counts_transExp <- merge(
    trans2gene[c("transcript_id", "gene_id")],
    counts_transExp, by.y = 0, by.x = "transcript_id",
    all.x = FALSE, all.y = TRUE
  ) %>% rename(feature_id = transcript_id)
  rownames(counts_geneExp) <- dict[match(rownames(counts_geneExp), dict$gene_name), "gene_id"]

  return(list(geneExpDf = counts_geneExp, transExpDf = counts_transExp, sampleDat = sampleData))
}

# ---------------------------------------------------------------------------------------------------------
# 2.3 Filter genes based on expression thresholds
# ---------------------------------------------------------------------------------------------------------
filterGenes <- function(counts_geneExp, sampleData, min_samps_gene_expr, min_gene_expr) {
  counts_geneExp %>%
    as.data.frame() %>%
    mutate(minSamplesGeneCount = rowSums(across(everything(), ~ . > min_gene_expr))) %>%
    subset(minSamplesGeneCount >= min_samps_gene_expr) %>%
    select(!minSamplesGeneCount)
}

# ---------------------------------------------------------------------------------------------------------
# 2.4 Filter transcripts based on expression thresholds
# ---------------------------------------------------------------------------------------------------------
filterTranscripts <- function(counts_transExp, sampleData,
                              min_samps_gene_expr, min_gene_expr,
                              min_samps_trans_expr_condition1, min_samps_trans_expr_condition2,
                              min_trans_expr, min_trans_prop) {
  # Get sample names per condition
  condition1_samples <- (sampleData %>% subset(condition == unique(sampleData$condition)[1]))$sample_id
  condition2_samples <- (sampleData %>% subset(condition == unique(sampleData$condition)[2]))$sample_id

  # Keep genes with minimum expression in enough samples
  geneIDs <- counts_transExp %>%
    select(!feature_id) %>%
    aggregate(. ~ gene_id, FUN = sum) %>%
    mutate(minSamplesGeneCount = rowSums(across(-gene_id, ~ . > min_gene_expr))) %>%
    subset(minSamplesGeneCount >= min_samps_gene_expr) %>%
    pull(gene_id)

  # Filter transcripts
  counts_transExp %>%
    subset(gene_id %in% geneIDs) %>%
    mutate(
      minSamplesTransCount_condition1 = rowSums(across(all_of(condition1_samples), ~ . > min_trans_expr)),
      minSamplesTransCount_condition2 = rowSums(across(all_of(condition2_samples), ~ . > min_trans_expr)),
      transcriptCount = rowSums(across(!c("feature_id", "gene_id")))
    ) %>%
    group_by(gene_id) %>%
    mutate(transcriptPerc = transcriptCount / sum(transcriptCount)) %>%
    subset(
      minSamplesTransCount_condition1 >= min_samps_trans_expr_condition1 |
      minSamplesTransCount_condition2 >= min_samps_trans_expr_condition2
    ) %>%
    subset(transcriptPerc > min_trans_prop) %>%
    ungroup() %>%
    select(!c(minSamplesTransCount_condition1, minSamplesTransCount_condition2, transcriptCount, transcriptPerc))
}

# ---------------------------------------------------------------------------------------------------------
# 2.5 Test for DEG and DTU per cell type
# ---------------------------------------------------------------------------------------------------------
testDE <- function(counts_geneExp, counts_transExp, sampleData) {
  # DEG analysis (DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = counts_geneExp, colData = sampleData, design = ~condition)
  dds <- DESeq(dds)

  # DTU analysis (DEXSeq)
  dxd <- DEXSeqDataSet(
    countData = as.matrix(counts_transExp[-c(1, 2)]),
    sampleData = sampleData,
    design = ~sample + exon + condition:exon,
    featureID = counts_transExp$feature_id,
    groupID = counts_transExp$gene_id
  )

  message("Estimating size factors...")
  dxd <- estimateSizeFactors(dxd)
  message("Estimating dispersion...")
  dxd <- estimateDispersions(dxd)
  message("Testing for differential exon usage...")
  dxd <- testForDEU(dxd)
  message("Estimating exon fold changes...")
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")

  return(list(deseq = dds, dexseq = dxd))
}

# ---------------------------------------------------------------------------------------------------------
# 2.6 Combine and export results
# ---------------------------------------------------------------------------------------------------------
combineResults <- function(deseq_obj, dexseq_obj, celltype, out_dir) {
  # Gene-level results
  results_gene <- data.frame(DESeq2::results(deseq_obj))
  results_gene$Gene <- rownames(results_gene)
  write.table(results_gene, paste0(out_dir, "/DESeq2Output_", celltype, ".tsv"),
              row.names = FALSE, quote = FALSE, sep = "\t")

  results_gene <- results_gene %>%
    mutate(Isoform = "Expression") %>%
    select(Gene, Isoform, log2FoldChange, padj) %>%
    rename("ExplogFC/FC" = log2FoldChange, BH = padj)

  # Transcript-level results
  results_transcript <- DEXSeqResults(dexseq_obj) %>%
    as.data.frame() %>%
    select(!genomicData)
  write.table(results_transcript, paste0(out_dir, "/DEXSeqOutput_", celltype, ".tsv"),
              row.names = FALSE, quote = FALSE, sep = "\t")

  logfc <- names(results_transcript)[startsWith(names(results_transcript), "log2fold")]
  results_trans <- results_transcript %>%
    mutate(foldChange = 2^!!sym(logfc)) %>%
    select(groupID, featureID, foldChange, padj) %>%
    rename(Gene = groupID, Isoform = featureID, "ExplogFC/FC" = foldChange, BH = padj)

  # Combine gene and transcript results
  results_combined <- rbind(results_gene, results_trans) %>% arrange(Gene)
  results_combined$`ExplogFC/FC` <- as.numeric(results_combined$`ExplogFC/FC`)
  results_combined$BH <- as.numeric(results_combined$BH)

  write.table(na.omit(results_combined),
              paste0(out_dir, "/isopret_", celltype, "_noNA.tsv"),
              row.names = FALSE, quote = FALSE, sep = "\t")
}

# =========================================================================================================
# 3. DEG and DTU analysis
# =========================================================================================================

# Load Seurat object and transcript–gene mapping
seurat_object <- readRDS("Seurat_FinalLongRead.rds")
dict <- read.table("resources/transcript2gene.tsv", sep = " ", header = TRUE)

# Define sample metadata for designs
severity <- c(
  "COVID0" = "critical", "COVID2" = "critical", "COVID5" = "critical", "COVID6" = "critical",
  "COVID9" = "critical", "COVID12" = "critical",
  "COVID1" = "moderate", "COVID3" = "moderate", "COVID4" = "moderate"
)
sampleData_ModCrit <- data.frame(sample_id = names(severity),
                                 condition = severity[names(severity)],
                                 row.names = names(severity))

sampleData_CovCon <- data.frame(sample_id = unique(seurat_object$orig.ident),
                                condition = ifelse(startsWith(unique(seurat_object$orig.ident), "COVID"), "COVID", "Control"),
                                row.names = unique(seurat_object$orig.ident))

dir <- "DE_results/"

# Filtering parameters
parameters <- list(
  CovVsCon = list(min_samps_gene_expr = 3, min_gene_expr = 10,
                  min_samps_gene_expr_dtu = 9, min_samps_trans_expr_condition1 = 2,
                  min_samps_trans_expr_condition2 = 6, min_trans_expr = 5, min_trans_prop = 0),
  ModVsCrit = list(min_samps_gene_expr = 3, min_gene_expr = 10,
                   min_samps_gene_expr_dtu = 6, min_samps_trans_expr_condition1 = 4,
                   min_samps_trans_expr_condition2 = 2, min_trans_expr = 5, min_trans_prop = 0)
)

# ---------------------------------------------------------------------------------------------------------
# 3.1 Loop over designs and cell types
# ---------------------------------------------------------------------------------------------------------
for (testDesign in c("CovVsCon", "ModVsCrit")) {
  message("Testing ", testDesign)
  if (testDesign == "CovVsCon") {
    sampleDat <- sampleData_CovCon
    outdir <- paste0(dir, "CovVsCon")
  } else {
    sampleDat <- sampleData_ModCrit
    outdir <- paste0(dir, "ModVsCrit")
  }

  for (ct in unique(seurat_object$celltype)) {
    message(ct)

    skipCelltype <- FALSE
    aggData <- tryCatch(
      aggregateExpressionMatrix(seurat_object, ct, dict, sampleDat),
      error = function(e) {
        skipCelltype <<- TRUE
      }
    )
    if (skipCelltype) {
      message("Too few samples per condition.")
      next
    }

    # Gene-level filtering
    counts_geneExp <- filterGenes(aggData$geneExpDf, aggData$sampleDat,
                                  parameters[[testDesign]]$min_samps_gene_expr,
                                  parameters[[testDesign]]$min_gene_expr)

    # Transcript-level filtering
    counts_transExp <- filterTranscripts(aggData$transExpDf, aggData$sampleDat,
                                         parameters[[testDesign]]$min_samps_gene_expr_dtu,
                                         parameters[[testDesign]]$min_gene_expr,
                                         parameters[[testDesign]]$min_samps_trans_expr_condition1,
                                         parameters[[testDesign]]$min_samps_trans_expr_condition2,
                                         parameters[[testDesign]]$min_trans_expr,
                                         parameters[[testDesign]]$min_trans_prop)

    cat("Genes after filtering:", nrow(counts_geneExp), "\n")
    cat("Transcripts after filtering:", nrow(counts_transExp), "\n")

    # Test for DEGs and DTUs
    skipCelltype <- FALSE
    testResults <- tryCatch(testDE(counts_geneExp, counts_transExp, aggData$sampleDat),
                            error = function(e) { skipCelltype <<- TRUE })
    if (skipCelltype) {
      message("Dataset too small for DE(X)Seq.")
      next
    }

    # Write outputs
    combineResults(testResults$deseq, testResults$dexseq, ct, outdir)
  }
}

