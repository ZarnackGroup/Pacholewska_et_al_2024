### ======================================================================== ###
###                                                                          ###
###                  Run differential binding site detection                 ###
###                                                                          ###
### ======================================================================== ###

### Without binding sites in background
### Total Set

### ============================================================================
### Load packages
### ----------------------------------------------------------------------------
#
library(BindingSiteFinder)
library(DESeq2)
library(AnnotationDbi)
library(rtracklayer)
library(GenomicFeatures)
library(tidyr)
library(dplyr)
library(IHW)
# source("/Users/mirko/Projects/BindingSiteStrength/03_markdowns/02_u2af65_hnrpc/02_prediction/scripts/Functions.R")
source("/Users/mirko/Projects/BindingSiteStrength/03_markdowns/02_u2af65_hnrpc/04_predictionFinal/scripts/Functions.R")

### ============================================================================
### Load gene annotation
### ----------------------------------------------------------------------------
#load("/Users/mirko/Projects/Annotations/human/gencode_29/filtered/gencode_v29_filtered.rda")
#anno.db = loadDb("/Users/mirko/Projects/Annotations/human/gencode_29/filtered/gencode_v29_filtered.sqlite")
load("/Users/mirko/Projects/Annotations/human/gencode_36/filtered/gencode_v36_filtered.rda")
anno.db = loadDb("/Users/mirko/Projects/Annotations/human/gencode_36/filtered/gencode_v36_filtered.sqlite")
gns = genes(anno.db)
idx = match(gns$gene_id, anno$gene_id)
elementMetadata(gns) = cbind(elementMetadata(gns), elementMetadata(anno)[idx,])
names(gns) = sub("\\..*", "", names(gns))
meta = data.frame(gene_id = gns$gene_id, gene_name = gns$gene_name, gene_type = gns$gene_type)
mcols(gns) = meta


### ============================================================================
### Load data
### ----------------------------------------------------------------------------
#
# Load binding sites
load("/Users/mirko/Projects/sf3b1/02_markdowns/01_transcriptome/01_bindingSites/data/bsTranscript.rda")
bindingSites = bsTranscript

# Load clip data
clipFilesWt = "/Users/mirko/Projects/sf3b1/01_data_subsamp/wt/cov/replicate"
clipFilesMut = "/Users/mirko/Projects/sf3b1/01_data_subsamp/mut/cov/replicate"
clipFiles = c(clipFilesWt, clipFilesMut)
clipFiles = list.files(clipFiles, pattern = ".bw$", full.names = TRUE)
clipFilesP = clipFiles[grep(clipFiles, pattern = "Plus")]
clipFilesM = clipFiles[grep(clipFiles, pattern = "Minus")]
# Organize clip data in dataframe
colData = data.frame(
    id = c(1:5),
    condition = factor(c("MUT", "MUT", "WT", "WT", "WT"), levels = c("MUT", "WT")),
    clPlus = clipFilesP,
    clMinus = clipFilesM)
# Make BindingSiteFinder object
bds = BSFDataSetFromBigWig(ranges = bindingSites, meta = colData)

### ============================================================================
### Construct background 
### ----------------------------------------------------------------------------
#
countObj = makeBsBackgroundMatrix(object = bds, geneRanges = gns, offset = 5, minCounts = 100)
m = as.data.frame(mcols(countObj))
m = m %>% select(starts_with("counts"))

### ============================================================================
### Run test with DESeq
### ----------------------------------------------------------------------------
#
# internal modification for col data
colDataMod = rbind.data.frame(colData,colData)
colDataMod$type = c(rep("bs", nrow(colData)),
                    rep("bg", nrow(colData)))
colDataMod$condition = factor(colDataMod$condition, levels = c("MUT", "WT"))
colDataMod$type = factor(colDataMod$type, levels = c("bs", "bg"))

# create SE object
se = SummarizedExperiment(assays = list(counts = as.matrix(m)),
                          rowRanges = granges(countObj), colData = colDataMod)

# set design
dds = DESeqDataSet(se, design = ~ type + condition + type:condition)
dds$condition = relevel(dds$condition, "WT")
dds$type = relevel(dds$type, "bs")
dds = DESeq(dds, fitType = "local", test = "LRT", reduced = ~ type + condition)
res = results(dds, contrast = c("condition", "MUT", "WT"), alpha = 0.05, filterFun = ihw, independentFiltering = TRUE)
res = lfcShrink(dds, res = res, contrast = c("condition", "MUT", "WT"), type = "ashr")


### ============================================================================
### Run GENE test with DESeq
### ----------------------------------------------------------------------------
#
# construct background dataframe for testing of gene expression level
m = as.data.frame(mcols(countObj))
m = m %>% 
    select(starts_with("counts.bg")) %>%
    unique() 
rownames(m) = sapply(strsplit(rownames(m),"\\."), `[`, 1)

# create SE object
seGene = SummarizedExperiment(assays = list(counts = as.matrix(m)), colData = colData)

# DESeq design
ddsGene = DESeqDataSet(seGene, design = ~ condition)
ddsGene$condition = relevel(ddsGene$condition, "WT")
ddsGene = DESeq(ddsGene)
resGene = results(ddsGene, contrast = c("condition", "MUT", "WT"), alpha = 0.05, filterFun = ihw, independentFiltering = TRUE)
resGene = lfcShrink(ddsGene, res = resGene, contrast = c("condition", "MUT", "WT"), type = "ashr")

### ============================================================================
### Export results
### ----------------------------------------------------------------------------
#
# deseq fetch results 
obj = countObj
colnames(res) = paste0("res.", colnames(res))
idx = match(names(obj), rownames(res))
mcols(obj) = cbind(mcols(obj), res[idx,])

# store result
searchRes = list(obj = obj, dds = dds, resGene = resGene, ddsGene = ddsGene)
save(searchRes, file = "/Users/mirko/Projects/sf3b1/02_markdowns/01_transcriptome/11_differentialBinding/data/searchRes.rda")


