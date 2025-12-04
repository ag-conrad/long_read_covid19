# Long-read single-cell transcriptomics in COVID-19

This repository contains the analysis code accompanying the preprint:  
**[Long-read single-cell transcriptomics reveals transcript isoform regulation in COVID-19](https://doi.org/10.1101/2025.09.23.676511)**  

All raw and processed data are available via Zenodo:  
📦 **[Zenodo Record 17135811](https://zenodo.org/records/17135811)**  

---

## Repository Structure

- **TranscriptStatistics.R** – Fisher's exact test for cell type/condition testing and PTPRC differential expression. Aggregation and plotting of transcript statistics. (Figure 1)
- **SeuratWorkflow.ipynb** – Notebook used for preprocessing, QC and annotating long-read data.  (Figure 1)
- **ComparisonShortLong.R** – Barcode matching and label transfer between long- and short-read Seurat datasets. (Figure 2)
- **DifferentialExpression.R** – DEG and DTU testing workflow using *DESeq2* and *DEXSeq* on Seurat-derived counts. (Figure 3)
- **EnrichmentAnalysis.R** – Pathway enrichment (Reactome, GO, Hallmark) and visualization scripts for DEGs and DTUs.  (Figure 3, 4, 5)
- **resources/** – Supporting annotation files (e.g., transcript-to-gene mappings).  

---

## Citation

If you use this code, please cite the preprint:  
> Köhler K., Morris L., Rakszewska A. *et al.* Long-read single-cell transcriptomics reveals transcript isoform regulation in COVID-19. *bioRxiv* (2025).  
> [https://doi.org/10.1101/2025.09.23.676511](https://doi.org/10.1101/2025.09.23.676511)

---
