# ---------------------------------------------------------------------------- #
#                                                                              #
#                         General helper functions                             #
#                                                                              #
# ---------------------------------------------------------------------------- #

minMaxNorm <- function(m) {
    nM = apply(m, 1, function(x){
        n = ((x - min(x)) / (max(x) - min(x)))
        return(n)
    })
    nM = t(nM)
    return(nM)
}

smoothing <- function(y, lambda, dim){
    p <- length(y)
    
    which(y != max(y))
    if ( length(which(y != max(y))) != 0 ){
        x = c(1:p)
        y.smooth = smooth.spline(x, y, spar = lambda)
        x_new = seq(1, p, p/dim)
        y_new = predict(y.smooth, x_new)$y
        y_new = y_new / max(y_new)
        return(y_new)
    } 
    
    x_new <- seq(1, p, p / dim)
    y_new <- rep( max(y), length(x_new) )
    return(y_new)
}

selMaxNorm <- function(m, cutoff = 2) {
    nM = apply(m, 1, function(x){
        if (max(x) > cutoff) {
            r = x / (max(x)/cutoff)
            return(r)
        }
        if (max(x) <= cutoff) {
            r = x
            return(r)
        }
    })
    nM = t(nM)
    return(nM)
}



myFormat <- function(x){
    format(x, big.mark = ",", decimal.mark = ".")
}


filterByRegion <- function(region, peaks, keepAbove) {
    # remove peaks not located on selected regions
    ols = findOverlaps(region, peaks)
    queries = unique(queryHits(ols))
    # apply selected cutoff to every region
    filteredPerRegion = sapply(queries, function(x){
        currentHits = subjectHits(ols)[queryHits(ols) == x]
        currentPeaks = peaks[currentHits]
        currentQuantiles = quantile(currentPeaks$score, probs = seq(0,1, by = 0.01))
        names(currentQuantiles) = seq(0,1, by = 0.01)
        filteredPeaks = currentPeaks[currentPeaks$score >= currentQuantiles[names(currentQuantiles) == keepAbove]]
        filteredPeaks
    })
    filteredPerRegion = unlist(GRangesList(filteredPerRegion))
    # remove peaks on multiple regions
    filteredPerRegion = filteredPerRegion[countOverlaps(filteredPerRegion) == 1] 
    return(filteredPerRegion)
}

basicVectorToNiceDf <- function(x){
    df = data.frame(Type = names(table(x$olType)), Freq = as.vector(table(x$olType)))
    df = df[order(df$Freq, decreasing = F),]
    df$Type = factor(df$Type, levels = df$Type)
    df$Frac = df$Freq / sum(df$Freq)
    df$ymax = cumsum(df$Frac)
    df$ymin = c(0, head(df$ymax, n=-1))
    df$labPos = (df$ymax + df$ymin) / 2
    df$NFrac = round(df$Frac * 100)
    df$NFrac2 = round(df$Freq / sum(df$Freq), digits = 4)
    df$NFracNice = df$NFrac2 * 100
    df$labelNice = paste0(format(df$Freq, big.mark = ",", decimal.mark = "."), " (", df$NFracNice, "%)")
    df$labelNice2 = paste0(df$Type, ": ", format(df$Freq, big.mark = ",", decimal.mark = "."), " (", df$NFracNice, "%)")
    return(df)
}



makeHmMatrixList <- function(clipSet, plot.topLim.size = 1){
    ### --------------------------------------------------------------------------
    ### Ranges
    ### --------------------------------------------------------------------------
    # select ranges
    # -> select within 100nt distance
    rngEndCurr = rngEnd
    rngEndCurr = rngEnd[rngEnd$distAltEnd <= 100]
    
    # shuffle ranges around
    # -> to avoid auto sorting by P value
    set.seed(1234)
    reorderIdx = sample(1:length(rngEndCurr))
    rngEndCurr = rngEndCurr[reorderIdx,]
    
    ### --------------------------------------------------------------------------
    ### Coverage
    ### --------------------------------------------------------------------------
    # coverage 3'AS outer position
    bdsSel = setRanges(clipSet, resize(rngEndCurr, fix = "end", width = 1) + 200)
    cov1 = coverageOverRanges(bdsSel, returnOptions = "merge_all_replicates", method = "mean")
    m1 = cov1
    # sorting by distance to 3'AS
    distIdx = order(rngEndCurr$distAltEnd, decreasing = FALSE)
    m1 = m1[distIdx,]
    
    # remove rows with no signal
    m1 = m1[rowSums(m1) > 0,]
    
    # scale values in matrix for plotting
    m1 = minMaxNorm(m1)
    m1 = t(apply(m1, 1, smoothing, lambda=0.05, dim=401))
    
    ### --------------------------------------------------------------------------
    ### Annotations
    ### --------------------------------------------------------------------------
    # set reference range for annotations
    range = rngEndCurr
    currM = m1
    # annotate with distance between alternative and constitutive 3'ss
    dfDist = data.frame(dist = range$distAltEnd[as.numeric(rownames(currM))])
    dfDist$dist[dfDist$dist > 100] = 100
    haDist = rowAnnotation(
        dist = anno_points(dfDist, size = unit(2, "points"), gp = gpar(col = "#595959")),
        width = unit(1, "cm"), annotation_name_gp = gpar(fontsize = 6)
    )
    # annotate with dPSI/ PSI values
    df = data.frame(deltaPSI = range$deltaPSI[as.numeric(rownames(currM))] * -1)
    haDeltaPSI = rowAnnotation(
        deltaPSI = anno_points(df, size = unit(2, "points"), gp = gpar(col = "#4c435a", alpha = 1), ylim = c(-1,1) ), 
        width = unit(2, "cm"), annotation_name_gp = gpar(fontsize = 6)
    )
    df = data.frame(psiWT = range$wtPSI[as.numeric(rownames(currM))] * -1)
    haWtPSI = rowAnnotation(
        psiWT = anno_points(df, size = unit(2, "points"), gp = gpar(col = "#0086b3", alpha = 0.5)), # , ylim = c(-1,0)
        width = unit(2, "cm"), annotation_name_gp = gpar(fontsize = 6)
    )
    df = data.frame(psiMUT = range$mutPSI[as.numeric(rownames(currM))] * -1)
    haMutPSI = rowAnnotation(
        psiMUT = anno_points(df, size = unit(2, "points"), gp = gpar(col = "#990000", alpha = 0.5)), # , ylim = c(-1,0)
        width = unit(2, "cm"), annotation_name_gp = gpar(fontsize = 6)
    )
    # annotate with significant results
    sig = range$sig[as.numeric(rownames(currM))]
    haSig = rowAnnotation(sigRes = sig, col = list(sigRes = c("FALSE" = "white", "TRUE" = "black")),
                          gp = gpar(lwd = 0.1),
                          show_annotation_name=F, annotation_name_gp = gpar(fontsize = 6),
                          annotation_legend_param = list(sigRes = list(nrow = 1, direction = "horizontal"))
    )
    # split by range
    split = factor(ifelse(dfDist$dist >= 12 & dfDist$dist <= 21, "2",
                          ifelse(dfDist$dist <= 12 & dfDist$dist <= 21, "1", "3")),
                   levels = c("1", "2", "3"))
    
    ### --------------------------------------------------------------------------
    ### Plotting
    ### --------------------------------------------------------------------------
    
    # use raw-data for barplot on top
    mm1 = cov1[rownames(cov1) %in% rownames(m1),]
    mm1 = mm1[,c(30:250)]
    
    # formatting of matrix outline
    rownames(m1) = NULL
    colnames(m1) = -200:200
    colnames(m1)[as.numeric(colnames(m1)) %% 100 != 0] = ""
    
    # change heatmap plotting frame
    m1 = m1[,c(30:250)]
    
    # annotate with top col means profile
    df1 = data.frame(sums = colMeans(mm1))
    haMeans1 = HeatmapAnnotation(cov = anno_barplot(df1, gp = gpar(fill = "#595959", col = "#595959"), ylim = c(0, plot.topLim.size)),
                                 height = unit(2, "cm"), show_annotation_name=F, annotation_name_gp = gpar(fontsize = 6))
    # color scale
    custom.col = colorRamp2(c(0,0.05,0.1,0.15,0.2,0.25), viridis(6, option = "mako", direction = -1))
    
    # plot heatmap
    h1 = Heatmap(m1,
                 column_title = "3'AS outer", name = "iCLIP signal",
                 cluster_rows = FALSE, cluster_columns = FALSE,
                 show_column_names = TRUE, show_row_names = FALSE,
                 col = custom.col,
                 row_split = split, row_gap = unit(0.1, "cm"),
                 column_names_gp = gpar(fontsize = 8),
                 top_annotation = haMeans1,
                 right_annotation = c(haSig, haDeltaPSI, haWtPSI, haMutPSI),
                 border = TRUE,
                 use_raster = TRUE, raster_by_magick = TRUE,
                 raster_magick_filter = "Spline",
                 raster_quality = 1, raster_resize_mat = mean,
                 heatmap_legend_param = list(
                     legend_width = unit(6, "cm"),
                     title_position = "topleft", direction = "horizontal"
                 )
    )
    plotList = h1
    matDim = dim(m1)
    l = list(plotList = plotList, matDim = matDim)
    return(l)
}



minMaxNorm <- function(m) {
    nM = apply(m, 1, function(x){
        n = ((x - min(x)) / (max(x) - min(x)))
        return(n)
    })
    nM = t(nM)
    return(nM)
}

smoothing <- function(y, lambda, dim){
    p <- length(y)
    
    which(y != max(y))
    if ( length(which(y != max(y))) != 0 ){
        x = c(1:p)
        y.smooth = smooth.spline(x, y, spar = lambda)
        x_new = seq(1, p, p/dim)
        y_new = predict(y.smooth, x_new)$y
        y_new = y_new / max(y_new)
        return(y_new)
    } 
    
    x_new <- seq(1, p, p / dim)
    y_new <- rep( max(y), length(x_new) )
    return(y_new)
}





makeCoverageAndCompare <- function(x, y, w, Xname, Yname, rngSs3, rngSs5, rngBp) {
    # set frame
    f3ss = rngSs3 + w
    f5ss = rngSs5 + w
    fbp = rngBp + w
    
    # calc 3'ss
    # ----------------------------------------------------------------------------
    # calc coverage
    cObjX = setRanges(x, f3ss)
    cCovX = coverageOverRanges(cObjX, returnOptions = "merge_all_replicates", method = "mean")
    cCovX = cCovX[rowSums(cCovX) > 0,] # This removes rows with all zeros !
    cObjY = setRanges(y, f3ss)
    cCovY = coverageOverRanges(cObjY, returnOptions = "merge_all_replicates", method = "mean")
    cCovY = cCovY[rowSums(cCovY) > 0,] # This removes rows with all zeros !
    dfCov3ssX = data.frame(pos = -w:w, mean = colMeans(cCovX), sd = colSds(cCovX), type = "3'SS", data = Xname) 
    dfCov3ssY = data.frame(pos = -w:w, mean = colMeans(cCovY), sd = colSds(cCovY), type = "3'SS", data = Yname)
    dfCov3ss = rbind(dfCov3ssX, dfCov3ssY)
    # testing
    testRes3 = sapply(1:ncol(cCovX), function(i){
        tr = t.test(x = cCovX[,i], y = cCovY[,i])
        # tr = wilcox.test(x = cCovX[,i], y = cCovY[,i])
        pVal = tr$p.value
        return(pVal)
    })
    dfTest3 = data.frame(pos = -w:w, pval = testRes3, pAdj = p.adjust(testRes3, method = "BH"), type = "3'SS")
    
    # calc 5'ss
    # ----------------------------------------------------------------------------
    # calc coverage
    cObjX = setRanges(x, f5ss)
    cCovX = coverageOverRanges(cObjX, returnOptions = "merge_all_replicates", method = "mean")
    cCovX = cCovX[rowSums(cCovX) > 0,] # This removes rows with all zeros !
    cObjY = setRanges(y, f5ss)
    cCovY = coverageOverRanges(cObjY, returnOptions = "merge_all_replicates", method = "mean")
    cCovY = cCovY[rowSums(cCovY) > 0,] # This removes rows with all zeros !
    dfCov5ssX = data.frame(pos = -w:w, mean = colMeans(cCovX), sd = colSds(cCovX), type = "5'SS", data = Xname) 
    dfCov5ssY = data.frame(pos = -w:w, mean = colMeans(cCovY), sd = colSds(cCovY), type = "5'SS", data = Yname)
    dfCov5ss = rbind(dfCov5ssX, dfCov5ssY)
    # testing
    testRes5 = sapply(1:ncol(cCovX), function(i){
        tr = t.test(x = cCovX[,i], y = cCovY[,i])
        pVal = tr$p.value
        return(pVal)
    })
    dfTest5 = data.frame(pos = -w:w, pval = testRes5, pAdj = p.adjust(testRes5, method = "BH"), type = "5'SS")
    
    # calc BP
    # ----------------------------------------------------------------------------
    # calc coverage
    cObjX = setRanges(x, fbp)
    cCovX = coverageOverRanges(cObjX, returnOptions = "merge_all_replicates", method = "mean")
    cCovX = cCovX[rowSums(cCovX) > 0,] # This removes rows with all zeros !
    cObjY = setRanges(y, fbp)
    cCovY = coverageOverRanges(cObjY, returnOptions = "merge_all_replicates", method = "mean")
    cCovY = cCovY[rowSums(cCovY) > 0,] # This removes rows with all zeros !
    dfCovBPssX = data.frame(pos = -w:w, mean = colMeans(cCovX), sd = colSds(cCovX), type = "BP", data = Xname) 
    dfCovBPssY = data.frame(pos = -w:w, mean = colMeans(cCovY), sd = colSds(cCovY), type = "BP", data = Yname)
    dfCovBP = rbind(dfCovBPssX, dfCovBPssY)
    # testing
    testResBP = sapply(1:ncol(cCovX), function(i){
        tr = t.test(x = cCovX[,i], y = cCovY[,i])
        pVal = tr$p.value
        return(pVal)
    })
    dfTestBP = data.frame(pos = -w:w, pval = testResBP, pAdj = p.adjust(testResBP, method = "BH"), type = "BP")
    
    # output
    # ----------------------------------------------------------------------------
    dfCov = rbind(dfCov3ss, dfCov5ss, dfCovBP)
    dfTest = rbind(dfTest3, dfTest5, dfTestBP)
    d = list(dfCov = dfCov, dfTest = dfTest)
    return(d)
}
