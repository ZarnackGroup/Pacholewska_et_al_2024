# ---------------------------------------------------------------------------- #
#                                                                              #
#                  MaxEnt specific helper functions                            #
#                                                                              #
# ---------------------------------------------------------------------------- #


### ----------------------------------------------------------------------------
###                   Make GRanges for all AGs near exons
### ----------------------------------------------------------------------------
### For all exons in the given GRanges all putative splices sites are recorded in 
### a GRanges object
### The searched width is defined by wIntron and wExon, describing how the far the
### window reaches into the exon/ intron
###
makeAgGrangesForExon <- function(x, wIntron, wExon){
    
    # get range and find AGs
    flankingRange = flank(x, width = wIntron, start = TRUE)
    flankingRangeSequence = RNAStringSet(getSeq(Hsapiens, flankingRange), use.names = TRUE)
    positionAg = stri_locate(str = flankingRangeSequence, regex = "AG", mode = "all")
    numberOfAGsPerRange = sapply(positionAg, function(x){nrow(as.data.frame(x))})
    
    grFlank = data.frame(
        chr = rep(seqnames(flankingRange), numberOfAGsPerRange),
        strand = rep(strand(flankingRange), numberOfAGsPerRange),
        startRange = rep(start(flankingRange), numberOfAGsPerRange),
        endRange = rep(end(flankingRange), numberOfAGsPerRange),
        exonID = rep(flankingRange$exonID, numberOfAGsPerRange),
        agPosition = unlist(sapply(positionAg, function(x){x[,1]}))
    ) %>%
        mutate(start = ifelse(strand == "+", startRange + agPosition -1, endRange - agPosition )) %>%
        mutate(end = ifelse(strand == "+", startRange + agPosition, endRange - agPosition +1)) %>%
        group_by(exonID) %>%
        mutate(agID = rev(paste0("AG+", 1:n()-1))) %>%
        mutate(agID2 = rev(1:n()-1)) %>%
        mutate(combID = paste0(exonID, "_", agID)) %>%
        makeGRangesFromDataFrame(. , keep.extra.columns = TRUE)
    
    
    exonsRange = resize(x, width = wExon, fix = "start")
    exonSequence = RNAStringSet(getSeq(Hsapiens, exonsRange), use.names = TRUE)
    positionAgExon = stri_locate(str = exonSequence, regex = "AG", mode = "all")
    numberOfAGsPerRangeExon = sapply(positionAgExon, function(x){nrow(as.data.frame(x))})
    
    grExon = data.frame(
        chr = rep(seqnames(exonsRange), numberOfAGsPerRangeExon),
        strand = rep(strand(exonsRange), numberOfAGsPerRangeExon),
        startRange = rep(start(exonsRange), numberOfAGsPerRangeExon),
        endRange = rep(end(exonsRange), numberOfAGsPerRangeExon),
        exonID = rep(exonsRange$exonID, numberOfAGsPerRangeExon),
        agPosition = unlist(sapply(positionAgExon, function(x){x[,1]}))
    ) %>%
        mutate(start = ifelse(strand == "+", startRange + agPosition -1, endRange - agPosition )) %>%
        mutate(end = ifelse(strand == "+", startRange + agPosition, endRange - agPosition +1)) %>%
        group_by(exonID) %>%
        mutate(agID = (paste0("AG-", 1:n()))) %>%
        mutate(agID2 = -(1:n())) %>%
        mutate(combID = paste0(exonID, "_", agID)) %>%
        filter(!is.na(start)) %>%
        makeGRangesFromDataFrame(. , keep.extra.columns = TRUE)
    
    gr = c(grFlank, grExon)
    gr = sort(sortSeqlevels(gr))
    
    return(gr)
}
### ----------------------------------------------------------------------------
###                   Make sequences for MaxEntScan
### ----------------------------------------------------------------------------
### Function that takes the position of AGs as GRanges and make the appropriate
### frame for MaxEntScan and exports it as .fasta
###
writeMaxEntSeq <- function(x, path) {
    # make 20 + 3 window range
    rngFrame = resize(x, width = 23, fix = "end")
    rngFramePlus = subset(rngFrame, strand(rngFrame) == "+")
    rngFrameMinus = subset(rngFrame, strand(rngFrame) == "-")
    rngFramePlus = shift(rngFramePlus, shift = +3)
    rngFrameMinus = shift(rngFrameMinus, shift = -3)
    rngFrame = sort(sortSeqlevels(c(rngFrameMinus, rngFramePlus)))
    
    # get sequence
    rngFrameSequence = DNAStringSet(getSeq(Hsapiens, rngFrame), use.names = TRUE)
    
    seqDf = data.frame(seq = as.character(rngFrameSequence), names = names(rngFrameSequence))
    seqDf$seqSplit = strsplit(seqDf$seq, "")
    
    require(seqinr)
    write.fasta(sequences = seqDf$seqSplit,
                names = seqDf$names,
                file.out = path, as.string = FALSE)
}
