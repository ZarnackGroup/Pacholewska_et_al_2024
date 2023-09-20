library(branchpointer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)

# set genome object
g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

# select all exons from the annotation (-> unfiltered)
exns = gtfToExons("/Users/mirko/Projects/Annotations/human/gencode_36/raw/gencode.v36.annotation.gtf")

# generate the search window for all exons (18-40 nt from 3'SS)
geneIds = unique(exns$gene_id)
queryIntrons <- makeBranchpointWindowForExons(geneIds, idType = "exon_id", exons = exns)
queryIntrons = subset(queryIntrons, gene_type == "protein_coding")

# make prediction
nWorker = 3
param = BiocParallel::MulticoreParam(workers = nWorker, exportglobals = TRUE)
register(param)

q = queryIntrons
xSplit = sample(1:nWorker, length(q), replace = TRUE)
qSplit = split(q, xSplit)

pred = bplapply(qSplit, function(x, BPPARAM = param){
    g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    branchpointer::predictBranchpoints(query = x, queryType = "region", BSgenome = g)
})
pred = unlist(GRangesList(pred))

saveRDS(pred, file = "./data/pred.rds")
