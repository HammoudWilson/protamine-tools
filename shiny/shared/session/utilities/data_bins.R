# handle bin exclusions, including at gc extremes
isExcludedBin <- function(chrom, genome, included, gc){
    !endsWith(chrom, paste0("-", genome)) |
    included == 0 | 
    gc < gcLimits[1] | 
    gc > gcLimits[2]
}
isIncludedAutosomeBin <- function(chrom, genome, included, gc, isAutosome){
    !isExcludedBin(chrom, genome, included, gc) & isAutosome
}
getIncludedAutosomeBins <- function(bins, gc, genome){
    isAutosome <- bins[, !startsWith(chrom, c("chrX", "chrY"))]
    isIncludedAutosomeBin(bins$chrom, genome, bins$included, gc, isAutosome)
}
getGenomeBins <- function(bins, genome){
    bins[, endsWith(chrom, paste0("-", genome))]
}

# similar for score bins from primary genome only
getIncludedAutosomeBins_scores <- function(bins){
    isAutosome <- bins[, !startsWith(chrom, c("chrX", "chrY"))]
    isIncluded <- bins$included == 1
    passesGc   <- bins$gc >= gcLimits[1] & bins$gc <= gcLimits[2]
    isAutosome & isIncluded & passesGc
}
