
# .libPaths("/home4/2711498i/R_packages")
# Install ggplot2 from CRAN
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyverse")
# Install Bioconductor installer and DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager :: install ("IHW")
BiocManager::install("DESeq2")

library(tidyverse)
library(data.table)
setDTthreads(8)
library(IHW)
library(DEWSeq)
library(ggrepel)
library(BiocParallel)
library(scales)
#library(dplyr)
dir.create("./DEW-seq/Regions_w50s20")

setwd("/home2/2459810a/iCLIP_analysis/NOP2_Innes/")

## Read data
print("read data")
countData <- fread("DEW-seq/Matrix_w50s20/Full_NOP2_iCLIP.matrix.txt.gz", sep = "\t")

print("read annotation")
annotationData <- fread("GENOMEDIR/sliding_window/HS.GRCh38.SINV_attribute-fix.flattened.w50s20.maptoid.txt.gz", sep = "\t")

extract_regions <- function(resultWindows, padjThresh, log2FoldChangeThresh, trackName, fileName) {
  resultWindows <- as.data.frame(resultWindows)
  resultRegions <- extractRegions(windowRes  = resultWindows,
                                  padjCol    = "p_adj_IHW",
                                  padjThresh = padjThresh, 
                                  log2FoldChangeThresh = log2FoldChangeThresh)
  toBED(windowRes = resultWindows,
        regionRes = resultRegions,
        padjCol    = "p_adj_IHW",
        padjThresh = padjThresh,
        log2FoldChangeThresh = log2FoldChangeThresh,
        trackName = trackName,
        fileName  = fileName
  )
  return(resultRegions)
}


print("NOP2_SINV")
##
print("prep data")
NOP2_SINV <- countData %>%
dplyr::select(unique_id, contains("SINV"), -contains("mock")) %>%
column_to_rownames(var = "unique_id")

NOP2_SINV_colData <- tibble(name = colnames(NOP2_SINV)) %>%
  mutate(type = case_when(
    #str_detect(name, "IP") ~ "IP",
    str_detect(name, "SMI") ~ "SMI", # changes to SMI if name contains SMI
    TRUE ~ "IP" # changes to IP if above conditions are not met
  )) %>%
  mutate(type = fct_relevel(type, c("IP", "SMI"))) %>%
  column_to_rownames(var = "name")

##
print("convert to DESeq dataset")
NOP2_SINV_ddw <- DESeqDataSetFromSlidingWindows(countData = NOP2_SINV,
                                           colData = NOP2_SINV_colData,
                                           annotObj = annotationData,
                                           tidy = FALSE,
                                           design = ~type,
                                           start0based = TRUE)

##
print("filter output based on frequency")
keep <- rowSums(counts(NOP2_SINV_ddw)) >= 5
NOP2_SINV_ddw <- NOP2_SINV_ddw[keep,]
rm(list = "keep")

print("estimate size factors")
NOP2_SINV_ddw <- estimateSizeFactors(NOP2_SINV_ddw)

## 
NOP2_SINV_ddw <- estimateDispersions(NOP2_SINV_ddw, fitType = "local", quiet = TRUE)
NOP2_SINV_ddw <- nbinomWaldTest(NOP2_SINV_ddw)

NOP2_SINV_resultWindows <- resultsDEWSeq(NOP2_SINV_ddw,
                                    contrast = c("type", "IP", "SMI"),
                                    parallel = FALSE,
                                    BPPARAM = 12,
                                    tidy = TRUE) %>% as_tibble

NOP2_SINV_resultWindows[,"p_adj_IHW"] <- adj_pvalues(ihw(pSlidingWindows ~ baseMean, 
                                                    data = NOP2_SINV_resultWindows,
                                                    alpha = 0.05,
                                                    nfolds = 10))

saveRDS(NOP2_SINV_ddw, "./DEW-seq/Matrix_w50s20/NOP2_SINV_ddw_ct5.RDS")
saveRDS(NOP2_SINV_resultWindows, "./DEW-seq/Matrix_w50s20/NOP2_SINV_resultWindows.RDS")


pdf("./Plots/DEW-seq_Dispersion_NOP2_SINV_ct5.pdf", width = 7, height = 7)
plotDispEsts(NOP2_SINV_ddw)
dev.off()

## extract regions
extract_regions(resultWindows = NOP2_SINV_resultWindows,
                padjThresh = 0.01,
                log2FoldChangeThresh = 2,
                trackName = "Binding site NOP2 SINV 0.01",
                fileName = "./DEW-seq/Regions_w50s20/NOP2_SINV_enrichedregions.log2FC_2.0_p0.01.bed") %>%
  saveRDS("./DEW-seq/Regions_w50s20/NOP2_SINV_enrichedregions.log2FC_2.0_p0.01.RDS")


extract_regions(resultWindows = NOP2_SINV_resultWindows,
                padjThresh = 0.05,
                log2FoldChangeThresh = 2,
                trackName = "Binding site NOP2 SINV 0.05",
                fileName = "./DEW-seq/Regions_w50s20/NOP2_SINV_enrichedregions.log2FC_2.0_p0.05.bed") %>%
  saveRDS("./DEW-seq/Regions_w50s20/NOP2_SINV_enrichedregions.log2FC_2.0_p0.05.RDS")


###### MOCK
print("NOP2_MOCK")
##
print("prep data")
NOP2_mock <- countData %>%
dplyr::select(unique_id, contains("mock"), -contains("SINV")) %>%
column_to_rownames(var = "unique_id")

NOP2_mock_colData <- tibble(name = colnames(NOP2_mock)) %>%
  mutate(type = case_when(
    #str_detect(name, "IP") ~ "IP",
    str_detect(name, "SMI") ~ "SMI", # changes to SMI if name contains SMI
    TRUE ~ "IP" # changes to IP if above conditions are not met
  )) %>%
  mutate(type = fct_relevel(type, c("IP", "SMI"))) %>%
  column_to_rownames(var = "name")

##
print("convert to DESeq dataset")
NOP2_mock_ddw <- DESeqDataSetFromSlidingWindows(countData = NOP2_mock,
                                           colData = NOP2_mock_colData,
                                           annotObj = annotationData,
                                           tidy = FALSE,
                                           design = ~type,
                                           start0based = TRUE)

##
print("filter output based on frequency")
keep <- rowSums(counts(NOP2_mock_ddw)) >= 5
NOP2_mock_ddw <- NOP2_mock_ddw[keep,]
rm(list = "keep")

print("estimate size factors")
NOP2_mock_ddw <- estimateSizeFactors(NOP2_mock_ddw)

## 
NOP2_mock_ddw <- estimateDispersions(NOP2_mock_ddw, fitType = "local", quiet = TRUE)
NOP2_mock_ddw <- nbinomWaldTest(NOP2_mock_ddw)

NOP2_mock_resultWindows <- resultsDEWSeq(NOP2_mock_ddw,
                                    contrast = c("type", "IP", "SMI"),
                                    parallel = FALSE,
                                    BPPARAM = 12,
                                    tidy = TRUE) %>% as_tibble

NOP2_mock_resultWindows[,"p_adj_IHW"] <- adj_pvalues(ihw(pSlidingWindows ~ baseMean, 
                                                    data = NOP2_mock_resultWindows,
                                                    alpha = 0.05,
                                                    nfolds = 10))

saveRDS(NOP2_mock_ddw, "./DEW-seq/Matrix_w50s20/NOP2_mock_ddw_ct5.RDS")
saveRDS(NOP2_mock_resultWindows, "./DEW-seq/Matrix_w50s20/NOP2_mock_resultWindows.RDS")


pdf("./Plots/DEW-seq_Dispersion_NOP2_mock_ct5.pdf", width = 7, height = 7)
plotDispEsts(NOP2_mock_ddw)
dev.off()

## extract regions
extract_regions(resultWindows = NOP2_mock_resultWindows,
                padjThresh = 0.01,
                log2FoldChangeThresh = 2,
                trackName = "Binding site NOP2 mock 0.01",
                fileName = "./DEW-seq/Regions_w50s20/NOP2_mock_enrichedregions.log2FC_2.0_p0.01.bed") %>%
  saveRDS("./DEW-seq/Regions_w50s20/NOP2_mock_enrichedregions.log2FC_2.0_p0.01.RDS")


extract_regions(resultWindows = NOP2_mock_resultWindows,
                padjThresh = 0.05,
                log2FoldChangeThresh = 2,
                trackName = "Binding site NOP2 mock 0.05",
                fileName = "./DEW-seq/Regions_w50s20/NOP2_mock_enrichedregions.log2FC_2.0_p0.05.bed") %>%
  saveRDS("./DEW-seq/Regions_w50s20/NOP2_mock_enrichedregions.log2FC_2.0_p0.05.RDS")


