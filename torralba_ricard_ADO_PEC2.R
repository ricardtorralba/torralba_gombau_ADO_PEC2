params <-
list(samples = 292L, seed = 1336336L, folder = "./Data", file_targets = "targets.csv", 
    file_counts = "counts.csv")

## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------
library(knitr)
opts_chunk$set(comment = NA,
               echo = F,
               prompt = T,
               message = F,
               warning = F,
               cache = T)


## ----libraries, include = F, message = F, warning = F---------------------------------------------------------------------------------
# Load packages
library(Biobase)
library(ggplot2)
library(DESeq2)
library(reshape2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(limma)


## ---- results = "hide"----------------------------------------------------------------------------------------------------------------
n_nit = 236
n_sfi = 42
n_eli = 14
n = 30
n_group = n/3


## -------------------------------------------------------------------------------------------------------------------------------------
set.seed(params$seed)
targ = read.csv(file = file.path(params$folder, params$file_targets),
                header = T)
targ$Group = factor(targ$Group)
rownames(targ) = gsub("-", ".", targ$Sample_Name)

nit = targ[targ$Group == "NIT",]
sfi = targ[targ$Group=="SFI",]
eli = targ[targ$Group=="ELI",]
obs_nit = sample(nrow(nit), 10, replace = F)
obs_sfi = sample(nrow(sfi), 10, replace = F)
obs_eli = sample(nrow(eli), 10, replace = F)
df_nit = nit[obs_nit,]
df_sfi = sfi[obs_sfi,]
df_eli = eli[obs_eli,]
coldata = rbind(df_nit, df_sfi, df_eli)
val_coldata = coldata$Sample_Name
vars_coldata = gsub("-", ".", val_coldata)

data = read.csv(file = file.path(params$folder, params$file_counts),
                sep = ";", header = T, row.names = 1)
countdata = data[vars_coldata]
rownames(countdata) = gsub("\\..*", "", rownames(countdata), fixed = F)
countdata = countdata[, rownames(coldata)]

ddsMat = DESeqDataSetFromMatrix(countData = countdata,
                                colData = coldata,
                                design = ~ Group)
ddsMat


## -------------------------------------------------------------------------------------------------------------------------------------
summary(ddsMat$Group)


## -------------------------------------------------------------------------------------------------------------------------------------
dds = ddsMat[rowSums(counts(ddsMat)) >1,]
nrow(dds)
removed_rows = nrow(ddsMat) - nrow(dds)


## -------------------------------------------------------------------------------------------------------------------------------------
norm_dds = estimateSizeFactors(dds)


## -------------------------------------------------------------------------------------------------------------------------------------
dds_res = DESeq(norm_dds, parallel = T)


## -------------------------------------------------------------------------------------------------------------------------------------
res_nit.sfi = results(dds_res, contrast = c("Group", "NIT", "SFI"))
resSig_nit.sfi = subset(res_nit.sfi, padj < 0.05)
tab_nit.sfi = table(res_nit.sfi$padj<0.05)
tab_nit.sfi


## -------------------------------------------------------------------------------------------------------------------------------------
kable(head(resSig_nit.sfi[order(resSig_nit.sfi$log2FoldChange),], 10),
      caption = "Most down-regulated transcripts in the NIT-SFI comparison")


## -------------------------------------------------------------------------------------------------------------------------------------
kable(head(resSig_nit.sfi[order(resSig_nit.sfi$log2FoldChange, decreasing = T),], 10),
      caption = "Most up-regulated transcrpits in the NIT-SFI comparison")


## -------------------------------------------------------------------------------------------------------------------------------------
res_nit.eli = results(dds_res, contrast = c("Group", "NIT", "ELI"))
resSig_nit.eli = subset(res_nit.eli, padj < 0.05)
tab_nit.eli = table(res_nit.eli$padj<0.05)
tab_nit.eli


## -------------------------------------------------------------------------------------------------------------------------------------
kable(head(resSig_nit.eli[order(resSig_nit.eli$log2FoldChange),], 10),
      caption = "Most down-regulated transcripts in the NIT-ELI comparison")


## -------------------------------------------------------------------------------------------------------------------------------------
kable(head(resSig_nit.eli[order(resSig_nit.eli$log2FoldChange, decreasing = T),], 10),
      caption = "Most up-regulated transcripts in the NIT-ELI comparison")


## -------------------------------------------------------------------------------------------------------------------------------------
res_sfi.eli = results(dds_res, contrast = c("Group", "SFI", "ELI"))
resSig_sfi.eli = subset(res_sfi.eli, padj < 0.05)
tab_sfi.eli = table(res_sfi.eli$padj<0.05)
tab_sfi.eli


## -------------------------------------------------------------------------------------------------------------------------------------
kable(head(resSig_sfi.eli[order(resSig_sfi.eli$log2FoldChange),], 10),
      caption = "Most down-regulated transcripts in the SFI-ELI comparison")


## -------------------------------------------------------------------------------------------------------------------------------------
kable(head(resSig_sfi.eli[order(resSig_sfi.eli$log2FoldChange, decreasing = T),], 10),
      caption = "Most up-regulated transcripts in the SFI-ELI comparison")


## -------------------------------------------------------------------------------------------------------------------------------------
resSig_nit.sfi$symbol = mapIds(org.Hs.eg.db,
                               keys = rownames(resSig_nit.sfi),
                               column = "SYMBOL",
                               keytype = "ENSEMBL",
                               multiVals = "first")
resSig_nit.sfi$entrez = mapIds(org.Hs.eg.db,
                               keys = row.names(resSig_nit.sfi),
                               column = "ENTREZID",
                               keytype = "ENSEMBL",
                               multiVals = "first")

resSig_nit.eli$symbol = mapIds(org.Hs.eg.db,
                               keys = rownames(resSig_nit.eli),
                               column = "SYMBOL",
                               keytype = "ENSEMBL",
                               multiVals = "first")
resSig_nit.eli$entrez = mapIds(org.Hs.eg.db,
                               keys = row.names(resSig_nit.eli),
                               column = "ENTREZID",
                               keytype = "ENSEMBL",
                               multiVals = "first")

resSig_sfi.eli$symbol = mapIds(org.Hs.eg.db,
                               keys = rownames(resSig_sfi.eli),
                               column = "SYMBOL",
                               keytype = "ENSEMBL",
                               multiVals = "first")
resSig_sfi.eli$entrez = mapIds(org.Hs.eg.db,
                               keys = row.names(resSig_sfi.eli),
                               column = "ENTREZID",
                               keytype = "ENSEMBL",
                               multiVals = "first")


## -------------------------------------------------------------------------------------------------------------------------------------
kable(head(resSig_nit.eli[order(resSig_nit.eli$padj),],10),
      caption = "Most differently expressed genes in the NIT-ELI comparison")


## -------------------------------------------------------------------------------------------------------------------------------------
kable(head(resSig_nit.sfi[order(resSig_nit.sfi$padj),],10),
      caption = "Most differently expressed genes in the NIT-SFI comparison")


## -------------------------------------------------------------------------------------------------------------------------------------
kable(head(resSig_sfi.eli[order(resSig_sfi.eli$padj),],10),
      caption = "Most differently expressed genes in the NIT-SFI comparison")


## -------------------------------------------------------------------------------------------------------------------------------------
common_nit.sfi_nit.eli = intersect(rownames(resSig_nit.sfi), rownames(resSig_nit.eli))
common_nit.sfi_sfi.eli = intersect(rownames(resSig_nit.sfi), rownames(resSig_sfi.eli))
common_nit.eli_sfi.eli = intersect(rownames(resSig_nit.eli), rownames(resSig_sfi.eli))
common_all = intersect(common_nit.sfi_nit.eli, rownames(resSig_sfi.eli))

comp = c(length(common_nit.sfi_nit.eli), length(common_nit.sfi_sfi.eli), length(common_nit.eli_sfi.eli), length(common_all))
multiple_comparison = data.frame("Number of genes" = comp)
rownames(multiple_comparison) = c("NIT-SFI AND NIT-ELI", "NIT-SFI AND SFI-ELI", "NIT-ELI AND SFI-ELI", "ALL")
kable(multiple_comparison, caption = "Number of genes differently expressed in more than one comparison")

