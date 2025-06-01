# utilities for describing bins

# unpack frac_gc limit values
unpackGcLimits <- function(env) strsplit(env$GC_LIMITS, ",")[[1]]

# handle bin exclusions, including at gc extremes
isExcludedBin <- function(chrom, genome, included, gc, gcLimits){
    !endsWith(chrom, paste0("-", genome)) |
    included == 0 | 
    gc < gcLimits[1] | 
    gc > gcLimits[2]
}
isIncludedAutosomeBin <- function(chrom, genome, included, gc, gcLimits, isAutosome){
    !isExcludedBin(chrom, genome, included, gc, gcLimits) & isAutosome
}
getIncudedAutosomeBins <- function(bins, genome, gcLimits){
    isAutosome <- bins[, !startsWith(chrom, c("chrX-", "chrY-"))]
    isIncludedAutosomeBin(bins$chrom, genome, bins$included, bins$gc, gcLimits, isAutosome)
}
