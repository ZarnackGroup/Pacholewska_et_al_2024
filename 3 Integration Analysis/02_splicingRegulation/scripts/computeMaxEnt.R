# ---------------------------------------------------------------------------- #
#                                                                              #
#                              MaxEnt Analysis                                 #
#                                                                              #
# ---------------------------------------------------------------------------- #
# 
# Script that takes canonical and alternative splice sites from the IsoSeq 
# alternative splicing analysis and locates all possible splice sites (AGs) in 
# a range of 23nt (20nt into the intron plus 3nt of the exon) for the MaxEnt
# based splice site strength analysis
#

library(rtracklayer)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(stringi)
source("./02_splicingRegulation/scripts/helperMaxEnt.R")

### ----------------------------------------------------------------------------
###                     Filter PacBio events
### ----------------------------------------------------------------------------
###
ss3PacBioEnd_1 = readRDS("./01_splicingMaps/data/rngEnd.rds") # this is the range based on the longer isoform
ss3PacBioAlt_1 = readRDS("./01_splicingMaps/data/rngAlt.rds") # this is the range based on the shorter isoform 

### ----------------------------------------------------------------------------
###                Reduce PacBio events to 3'splice sites
### ----------------------------------------------------------------------------
###
ss3PacBioEnd_2 = unique(resize(ss3PacBioEnd_1, fix = "end", width = 1))
ss3PacBioAlt_2 = unique(resize(ss3PacBioAlt_1, fix = "end", width = 1))

tmpEnd_plus = GenomicRanges::shift(subset(ss3PacBioEnd_2, strand == "+"), shift = -1) %>% resize(., width = 2, fix = "end")
tmpEnd_minus = GenomicRanges::shift(subset(ss3PacBioEnd_2, strand == "-"), shift = 1) %>% resize(., width = 2, fix = "end")
tmpAlt_plus = GenomicRanges::shift(subset(ss3PacBioAlt_2, strand == "+"), shift = -1) %>% resize(., width = 2, fix = "end")
tmpAlt_minus = GenomicRanges::shift(subset(ss3PacBioAlt_2, strand == "-"), shift = 1) %>% resize(., width = 2, fix = "end")

ss3PacBioEnd_2 = sort(sortSeqlevels(c(tmpEnd_plus, tmpEnd_minus)))
ss3PacBioAlt_2 = sort(sortSeqlevels(c(tmpAlt_plus, tmpAlt_minus)))

### ----------------------------------------------------------------------------
###                Annotate splice sites by deltaPSI
### ----------------------------------------------------------------------------
###
mcols(ss3PacBioEnd_2)$regulation = ifelse(ss3PacBioEnd_2$deltaPSI > 0, "shorter_introns", "longer_introns")
mcols(ss3PacBioAlt_2)$regulation = ifelse(ss3PacBioAlt_2$deltaPSI > 0, "shorter_introns", "longer_introns")
### NOTE:
### -> if deltaPSI goes up (>0), 
### ---> then the longer isoform ends at the canonical 3'SS 
### ---> then the shorter isoform ends in the intron (exon gets longer)

### ----------------------------------------------------------------------------
###               Clean splice sites
### ----------------------------------------------------------------------------
###
# set exon id for event
ss3PacBioEnd_2$exonID = names(ss3PacBioEnd_2)
ss3PacBioAlt_2$exonID = names(ss3PacBioAlt_2)

# keep only 3'splice sites with AG at start - END
ss3PacBioEnd_ag = RNAStringSet(getSeq(Hsapiens, ss3PacBioEnd_2), use.names = TRUE)
countAg = stri_count(str = ss3PacBioEnd_ag, regex = "AG", mode = "all")
ss3PacBioEnd = ss3PacBioEnd_2[countAg == 1]

# keep only 3'splice sites with AG at start - ALT
ss3PacBioAlt_ag = RNAStringSet(getSeq(Hsapiens, ss3PacBioAlt_2), use.names = TRUE)
countAg = stri_count(str = ss3PacBioAlt_ag, regex = "AG", mode = "all")
ss3PacBioAlt = ss3PacBioAlt_2[countAg == 1]

### ----------------------------------------------------------------------------
###               adjust positions to exon ranges
### ----------------------------------------------------------------------------
###
toExonPosition <- function(x) {
    tmp_plus = GenomicRanges::shift(subset(x, strand == "+"), shift = 1) %>% resize(., width = 1, fix = "end")
    tmp_minus = GenomicRanges::shift(subset(x, strand == "-"), shift = -1) %>% resize(., width = 1, fix = "end")
    
    relocated = sort(sortSeqlevels(c(tmp_plus, tmp_minus)))
    return(relocated)
}
ss3PacBioExon = toExonPosition(ss3PacBioEnd)

# exports 
export(granges(ss3PacBioExon), "./02_splicingRegulation/data/ss3PacBioExon.bed", format = "BED")
export(granges(ss3PacBioEnd), "./02_splicingRegulation/data/ss3PacBioEnd.bed", format = "BED")
export(granges(ss3PacBioAlt), "./02_splicingRegulation/data/ss3PacBioAlt.bed", format = "BED")

saveRDS(ss3PacBioEnd_1, file = "./02_splicingRegulation/data/ss3PacBioEnd_1.rds")
saveRDS(ss3PacBioEnd_2, file = "./02_splicingRegulation/data/ss3PacBioEnd_2.rds")
saveRDS(ss3PacBioEnd, file = "./02_splicingRegulation/data/ss3PacBioEnd.rds")


### ----------------------------------------------------------------------------
###                     Find all AGs in PacBio exons
### ----------------------------------------------------------------------------
###
grAgPacBio = makeAgGrangesForExon(x = ss3PacBioExon, wExon = 50, wIntron = 100)
names(grAgPacBio) = grAgPacBio$combID
export(granges(grAgPacBio), "./02_splicingRegulation/data/grAgPacBio.bed", format = "BED")
saveRDS(grAgPacBio, file = "./02_splicingRegulation/data/grAgPacBio.rds")

### ----------------------------------------------------------------------------
###             Locate the matching alternative AG for each 3'SS
### ----------------------------------------------------------------------------
###
grAgPacBio_match = grAgPacBio %>% 
    as.data.frame() %>%
    mutate(matchID = paste0(start, "_", end, "_", exonID))
ss3PacBioAlt_match = ss3PacBioAlt %>%
    as.data.frame() %>%
    mutate(matchID = paste0(start, "_", end, "_", exonID))

isAlternative = match(grAgPacBio_match$matchID, ss3PacBioAlt_match$matchID)
isAlternative[!is.na(isAlternative)] = "Yes"
isAlternative[is.na(isAlternative)] = "No"

grAgPacBio$isAlternative = isAlternative


### ----------------------------------------------------------------------------
###             write fasta file in maxEnt format 
### ----------------------------------------------------------------------------
###
writeMaxEntSeq(grAgPacBio, path = "./02_splicingRegulation/data/seqsAGsPacBio.fasta")





