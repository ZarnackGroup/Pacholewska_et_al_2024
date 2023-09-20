# ---------------------------------------------------------------------------- #
#                                                                              #
#                         Branchpointer Analysis                               #
#                                                                              #
# ---------------------------------------------------------------------------- #
# The branchpointer analysis takes all putative alternative AGs and calculates
# the strength of potential branchpoints using the branchpointer prediction.
# For the prediction to work a spcific window ranging from -18 to -44nt from the
# splice site has to be calculated. 
#
library(branchpointer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(readr)

### ----------------------------------------------------------------------------
###                                Load data
### ----------------------------------------------------------------------------
ss3PacBioEnd_1 = readRDS("/Users/mirko/Projects/sf3b1/02_markdowns/03_clean/02_splicingRegulation/data/ss3PacBioEnd_1.rds")
g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38


### ----------------------------------------------------------------------------
###                        Caclulate branchpoints 
### ----------------------------------------------------------------------------
### make spcific branchpointer windows
rng.minus = subset(ss3PacBioEnd_1, strand == "-")
rng.minus = granges(rng.minus) 
rng.minus = GenomicRanges::shift(rng.minus, shift = 18)
rng.minus = GenomicRanges::resize(rng.minus, width = 27, fix = "end")

rng.plus = subset(ss3PacBioEnd_1, strand == "+")
rng.plus = granges(rng.plus) 
rng.plus = GenomicRanges::shift(rng.plus, shift = -18)
rng.plus = GenomicRanges::resize(rng.plus, width = 27, fix = "end")

bpQuery = c(rng.plus, rng.minus)
bpQuery = sortSeqlevels(bpQuery)
bpQuery = sort(bpQuery)

export(bpQuery, con = "./bpQuery.bed", format = "BED")

bpQuery$id = ss3PacBioEnd_1$id
bpQuery$to_3prime = 18
bpQuery$to_5prime = width(ss3PacBioEnd_1)

grBpPacBioBp = branchpointer::predictBranchpoints(bpQuery, queryType = "region", BSgenome = g)
saveRDS(grBpPacBioBp, file = "/Users/mirko/Projects/sf3b1/02_markdowns/03_clean/02_splicingRegulation/data/grBpPacBioBp.rds")



