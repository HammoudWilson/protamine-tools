# utilities for describing bins

# unpack the fractionGC limit values
unpackGcLimits <- function(env) strsplit(env$GC_LIMITS, ",")[[1]]

# handle bin exclusions, including at gc extremes
isExcludedBin <- function(excluded, pct_gc, gcLimits){
    excluded == 1 | pct_gc < gcLimits[1] | pct_gc > gcLimits[2]
}
isIncludedAutosomeBin <- function(excluded, pct_gc, nAlleles, gcLimits){
    !isExcludedBin(excluded, pct_gc, gcLimits) & nAlleles == 2
}
getIncudedAutosomeBins <- function(bins, gcLimits){
    isIncludedAutosomeBin(bins$excluded, bins$pct_gc, bins$nAlleles, gcLimits)
}



