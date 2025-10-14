library(Seurat)
library(dplyr)
library(DESeq2)
library(DEXSeq)
library(DRIMSeq)

#############################################################################################################################################################
#################################################################### Utility functions ######################################################################
#############################################################################################################################################################

#################################################Aggregate expression counts per sample for one celltype#####################################################

aggregateExpressionMatrix = function(seurat_obj, cellType, trans2gene, sampleData){
    ##subset by cell type
    cellTypeSubset <- seurat_obj %>% subset(celltype==cellType)

    ##aggregate counts per sample on gene- and transcript-level 
    counts_geneExp <- AggregateExpression(cellTypeSubset, group.by = 'orig.ident', assay = 'RNA')$RNA #gene-level
    counts_transExp <- AggregateExpression(cellTypeSubset, group.by = 'orig.ident', assay = 'Transcript')$Transcript #transcript-level

    ##check if all samples are present in cell subset
    samples = colnames(counts_geneExp)[colnames(counts_geneExp) %in% sampleData$sample_id]
    sampleData = sampleData[samples,]
    sampleDataCheck = sampleData |> group_by(condition) |> mutate(groupSize= n()) |> ungroup() 
    if(min(sampleDataCheck$groupSize) < 2) stop()

    ##reorder columns and merge with dict 
    counts_geneExp <- counts_geneExp[,rownames(sampleData)]
    counts_transExp <- counts_transExp[,rownames(sampleData)]
    counts_transExp = merge(trans2gene[c( 'transcript_id','gene_id')], counts_transExp, by.y = 0, by.x = 'transcript_id', all.x = FALSE, all.y = TRUE) |>
                        dplyr::rename(feature_id = transcript_id)
    rownames(counts_geneExp) <- dict[match(rownames(counts_geneExp), dict$gene_name),"gene_id"]

    return(list(geneExpDf = counts_geneExp, transExpDf = counts_transExp, sampleDat = sampleData))
}

##########################Aggregate across all cells to mimic bulk RNAseq##################################### 

aggregateExpressionToBulk = function(seurat_obj, cellType, trans2gene, sampleData){

    counts_geneExp <- AggregateExpression(seurat_obj, group.by = 'orig.ident', assay = 'RNA')$RNA #gene-level
    counts_transExp <- AggregateExpression(seurat_obj, group.by = 'orig.ident', assay = 'Transcript')$Transcript #transcript-level

    counts_transExp = merge(trans2gene[c( 'transcript_id','gene_id')], counts_transExp, by.y = 0, by.x = 'transcript_id', all.x = FALSE, all.y = TRUE) |>
                        dplyr::rename(feature_id = transcript_id)
    rownames(counts_geneExp) <- dict[match(rownames(counts_geneExp), dict$gene_name),"gene_id"]

    return(list(geneExpDf = counts_geneExp, transExpDf = counts_transExp, sampleDat = sampleData))
}

########################################################Filter genes based on expression per condition########################################################

filterGenes = function(counts_geneExp, sampleData, min_samps_gene_expr, min_gene_expr){

    ##keep genes with min_gene_expr counts in at least min_samps_gene_expr
    counts_geneExp = counts_geneExp |> 
            as.data.frame() |>
            mutate(minSamplesGeneCount = rowSums(across(everything(), ~ . > min_gene_expr))) |> 
            subset(minSamplesGeneCount >= min_samps_gene_expr) |> 
            dplyr::select(!c(minSamplesGeneCount)) 

    return(counts_geneExp)
}

###################################################Filter transcripts based on expression per condition########################################################

filterTranscripts = function(counts_transExp, sampleData, min_samps_gene_expr, min_gene_expr, min_samps_trans_expr_condition1, min_samps_trans_expr_condition2, min_trans_expr, min_trans_prop){

    ##get sample names per condition
    condition1_samples = (sampleData |> subset(condition == unique(sampleData$condition)[1]))$sample_id
    condition2_samples = (sampleData |> subset(condition == unique(sampleData$condition)[2]))$sample_id
    #cat("condition 1: ", condition1_samples, "\n")
    #cat("condition 2: ", condition2_samples, "\n")

    ##keep genes with min_gene_expr in min_samps_gene_expr
    geneIDs <- counts_transExp |> dplyr::select(!feature_id) |> 
                aggregate(. ~ gene_id, FUN = sum) |> 
                mutate(minSamplesGeneCount = rowSums(across(-gene_id, ~ . > min_gene_expr))) |> 
                subset(minSamplesGeneCount >= min_samps_gene_expr) |> 
                pull(gene_id)

    ##keep transcripts with min_trans_expr counts in at least min_samps_trans_expr_condition1/min_samps_trans_expr_condition2 
    counts_transExp = counts_transExp |> 
        subset(gene_id %in% geneIDs) |> 
        mutate(minSamplesTransCount_condition1 = rowSums(across(all_of(condition1_samples), ~ . > min_trans_expr))) |> 
        mutate(minSamplesTransCount_condition2 = rowSums(across(all_of(condition2_samples), ~ . > min_trans_expr))) |> 
        mutate(transcriptCount = rowSums(across(!c("feature_id", "gene_id")))) |>
        group_by(gene_id) |>
        mutate(transcriptPerc = transcriptCount/sum(transcriptCount))  |> 
        subset(minSamplesTransCount_condition1 >= min_samps_trans_expr_condition1 | minSamplesTransCount_condition2 >= min_samps_trans_expr_condition2) |> 
        subset(transcriptPerc > min_trans_prop) |>
        ungroup() |>
        dplyr::select(!c(minSamplesTransCount_condition1, minSamplesTransCount_condition2, transcriptCount, transcriptPerc))

    return(counts_transExp)
}

###############################################################Test for DEG and DTU per celltype###############################################################

testDE = function(counts_geneExp, counts_transExp, sampleData){
    ##run DESeq2
    dds <- DESeqDataSetFromMatrix(countData = counts_geneExp,
                                colData = sampleData,
                                design = ~ condition)                                    
    dds <- DESeq(dds)

    ##run DEXSeq
    dxd <- DEXSeqDataSet(countData = as.matrix(counts_transExp[-c(1,2)]), sampleData = sampleData,
                                    design = ~sample + exon + condition:exon, featureID = counts_transExp$feature_id,
                                            groupID = counts_transExp$gene_id)
    ##parallelization possible if subsets are large enough
    BPPARAM = MulticoreParam(32)
    message('estimating size factors..')
    dxd = estimateSizeFactors(dxd)
    message('estimating dispersion..')
    dxd = estimateDispersions(dxd)#, BPPARAM=BPPARAM)
    message('testing for differential exon usage..')
    dxd = testForDEU(dxd)#, BPPARAM=BPPARAM)
    message('estimating exon fold changes..')
    dxd = estimateExonFoldChanges(dxd, fitExpToVar="condition")#,BPPARAM=BPPARAM,  )    

    return(list(deseq = dds, dexseq = dxd))
}

####################################################Combine and write test results into output directory########################################################

combineResults = function(deseq_obj, dexseq_obj, celltype, out_dir){
    results_gene <- data.frame(DESeq2::results(deseq_obj))
    results_gene$Gene <- rownames(results_gene)
    write.table(results_gene, paste0(out_dir, "/DESeq2Output_", celltype ,'.tsv'), row.names = FALSE, quote = FALSE, sep = "\t")
    results_gene <- results_gene |> mutate(Isoform = 'Expression') |>
                                    dplyr::select(Gene, Isoform, log2FoldChange, padj) |>
                                    dplyr::rename('ExplogFC/FC' = log2FoldChange, BH = padj)

    results_transcript <- DEXSeqResults(dexseq_obj) |> 
                            as.data.frame() |>
                            dplyr::select(!genomicData)
    write.table(results_transcript, paste0(out_dir, "/DEXSeqOutput_", celltype ,'.tsv'), row.names = FALSE, quote = FALSE, sep = "\t")

    logfc = names(results_transcript)[startsWith(names(results_transcript), 'log2fold')]
    results_trans <- results_transcript |>
                        mutate(foldChange = 2^!!sym(logfc)) |>
                        dplyr::select(groupID, featureID, foldChange, padj) |>
                        dplyr::rename(Gene = groupID, Isoform = featureID, 'ExplogFC/FC' = foldChange, BH = padj)

    results_combined <- rbind(results_gene, results_trans) |> 
                            arrange(Gene)

    results_combined$`ExplogFC/FC` <- as.numeric(results_combined$`ExplogFC/FC`)
    results_combined$BH <- as.numeric(results_combined$BH)

    write.table(na.omit(results_combined), paste0(out_dir, "/isopret_", celltype ,'_noNA.tsv'), row.names = FALSE, quote = FALSE, sep = "\t")

}




#############################################################################################################################################################
################################################################### DEG and DTU analysis ####################################################################
#############################################################################################################################################################

##read data
seurat_object = readRDS("NS_Integrated.rds")
dict = read.table("resources/transcript2gene.tsv", sep = " ", header = TRUE)

##define test design dataframe
severity = c("COVID0" = "critical", "COVID2" = "critical", "COVID5" = "critical", "COVID6" = "critical", "COVID9" = "critical", "COVID12" = "critical",
            "COVID1" = "moderate", "COVID3" = "moderate", "COVID4" = "moderate")
sampleData_ModCrit <- data.frame(sample_id = names(severity),
                        condition = severity[ names(severity)])
rownames(sampleData_ModCrit) <- sampleData_ModCrit$sample_id
sampleData_CovCon <- data.frame(sample_id = unique(seurat_object$orig.ident),
                        condition = ifelse(startsWith(unique(seurat_object$orig.ident), 'COVID'),'COVID', 'Control'),
                        row.names = unique(seurat_object$orig.ident))
dir = "DE_results/"

##set filtering parameters
parameters = list("CovVsCon" = list(min_samps_gene_expr = 3,
                                min_gene_expr = 10,
                                min_samps_gene_expr_dtu = 9, 
                                min_samps_trans_expr_condition1 = 2,
                                min_samps_trans_expr_condition2 = 6,
                                min_trans_expr = 5,
                                min_trans_prop = 0),
                  "ModVsCrit" = list(min_samps_gene_expr = 3,
                                min_gene_expr = 10,
                                min_samps_gene_expr_dtu = 6, 
                                min_samps_trans_expr_condition1 = 4,
                                min_samps_trans_expr_condition2 = 2,
                                min_trans_expr = 5,
                                min_trans_prop = 0)
                )


##set filtering parameters
parameters_drim = list("CovVsCon" = list(min_samps_gene_expr = 3,
                                min_gene_expr = 10,
                                min_samps_gene_expr_dtu = 11, 
                                min_samps_trans_expr = 3,
                                min_trans_expr = 5, 
                                min_samps_feature_prop = 3, 
                                min_feature_prop = 0.1),
                  "ModVsCrit" = list(min_samps_gene_expr = 3,
                                min_gene_expr = 10,
                                min_samps_gene_expr_dtu = 6, 
                                min_samps_trans_expr = 3,
                                min_trans_expr = 5, 
                                min_samps_feature_prop = 3, 
                                min_feature_prop = 0.1)
                )

for(testDesign in c("CovVsCon", "ModVsCrit")){
    ##specify design and directory
    message("Testing ", testDesign)
    if(testDesign == "CovVsCon"){
        sampleDat = sampleData_CovCon
        outdir = paste0(dir, "CovVsCon")
    } else {
        sampleDat = sampleData_ModCrit
        outdir = paste0(dir, "ModVsCrit")
    }

    ##test per celltype
    for(ct in unique(seurat_object$celltype)){
        message(ct, "\n")

        ##aggregate counts
        skipCelltype <- FALSE
        aggData = tryCatch(aggregateExpressionMatrix(seurat_obj = seurat_object, cellType = ct, trans2gene = dict, sampleData = sampleDat), error = function(e) { skipCelltype <<- TRUE})
        if(skipCelltype) {
            message("Number of samples per condition is too small!") 
            next
        }  
        ##gene-level filtering
        counts_geneExp = filterGenes(counts_geneExp = aggData$geneExpDf, sampleData = aggData$sampleDat, 
                                    min_samps_gene_expr = parameters[[testDesign]]$min_samps_gene_expr, 
                                    min_gene_expr = parameters[[testDesign]]$min_gene_expr)                        
        ##gene- and transcript-level filtering
        counts_transExp = filterTranscripts(counts_transExp = aggData$transExpDf, sampleData = aggData$sampleDat, 
                                    min_samps_gene_expr = parameters[[testDesign]]$min_samps_gene_expr_dtu, min_gene_expr = parameters[[testDesign]]$min_gene_expr, 
                                    min_samps_trans_expr_condition1 = parameters[[testDesign]]$min_samps_trans_expr_condition1, min_samps_trans_expr_condition2 = parameters[[testDesign]]$min_samps_trans_expr_condition2, 
                                    min_trans_expr = parameters[[testDesign]]$min_trans_expr, min_trans_prop = parameters[[testDesign]]$min_trans_prop)
        
        ##gene- and transcript-level filtering  with DrimSeq
        #drim_obj = dmDSdata(counts = aggData$transExpDf, samples = aggData$sampleDat)
        #skipCelltype <- FALSE
        #counts_transExp = tryCatch(counts(dmFilter(drim_obj, min_samps_gene_expr = parameters_drim[[testDesign]]$min_samps_gene_expr_dtu, min_samps_feature_expr = parameters_drim[[testDesign]]$min_samps_trans_expr, 
        #                            min_gene_expr = parameters_drim[[testDesign]]$min_gene_expr, min_feature_expr = parameters_drim[[testDesign]]$min_trans_expr, 
        #                            min_samps_feature_prop = parameters_drim[[testDesign]]$min_samps_feature_prop, min_feature_prop = parameters_drim[[testDesign]]$min_feature_prop)), 
        #                           error = function(e) { skipCelltype <<- TRUE})
        #if(skipCelltype) {
        #    message("No transcript left after DrimSeq filtering!") 
        #    next
        #} 


        cat("Number of genes after filtering: ", nrow(counts_geneExp), "\n")
        cat("Number of transcripts after filtering: ", nrow(counts_transExp), "\n")
        ##test for DEG and DTUs
        skipCelltype <- FALSE
        testResults = tryCatch(testDE(counts_geneExp, counts_transExp, aggData$sampleDat), error = function(e) { skipCelltype <<- TRUE})
        if(skipCelltype) {
            message("Remaining dataset too small for DE(X)Seq!") 
            next
        }     
        ##write output
        combineResults(deseq_obj = testResults$deseq, dexseq_obj = testResults$dexseq, celltype = ct, out_dir = outdir)
    }

}



######################################################### bulk simulation ###############################################################

aggData = aggregateExpressionToBulk(seurat_obj = seurat_object, trans2gene = dict, sampleData = sampleData_CovCon)

##gene-level filtering
counts_geneExp = filterGenes(counts_geneExp = aggData$geneExpDf, sampleData = aggData$sampleDat, 
                            min_samps_gene_expr = parameters[["CovVsCon"]]$min_samps_gene_expr, 
                            min_gene_expr = parameters[["CovVsCon"]]$min_gene_expr)                        
##gene- and transcript-level filtering
counts_transExp = filterTranscripts(counts_transExp = aggData$transExpDf, sampleData = aggData$sampleDat, 
                            min_samps_gene_expr = parameters[["CovVsCon"]]$min_samps_gene_expr_dtu, min_gene_expr = parameters[["CovVsCon"]]$min_gene_expr, 
                            min_samps_trans_expr_condition1 = parameters[["CovVsCon"]]$min_samps_trans_expr_condition1, min_samps_trans_expr_condition2 = parameters[["CovVsCon"]]$min_samps_trans_expr_condition2, 
                            min_trans_expr = parameters[["CovVsCon"]]$min_trans_expr, min_trans_prop = parameters[["CovVsCon"]]$min_trans_prop)

counts_transExp = counts_transExp |> dplyr::select(feature_id, gene_id, aggData$sampleDat$sample_id)
counts_geneExp = counts_geneExp |> dplyr::select(aggData$sampleDat$sample_id)

testResults = testDE(counts_geneExp, counts_transExp, aggData$sampleDat)
combineResults(deseq_obj = testResults$deseq, dexseq_obj = testResults$dexseq, celltype = "allCells", out_dir = paste0(dir, "CovVsCon"))
