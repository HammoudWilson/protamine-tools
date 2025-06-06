#----------------------------------------------------------------------
# nbinomCountsGC class generic methods, called as method(obj)
#----------------------------------------------------------------------

# prevent errors that arise when the GC fraction is outside the range of the model
suppressGcOutliers.nbinomCountsGC <- function(nb, fractionGC){
    allowedGC <- range(nb$model$fractionGC)
    pmax(allowedGC[1], pmin(allowedGC[2], fractionGC))
}

# returns an object with both x and y data (since x may have been our default) binned fraction GC
predict.nbinomCountsGC <- function(nb, fractionGC=NULL, type=c('mu', 'peak', 'adjustedPeak'), peakThreshold=10){
    rows <- if(is.null(fractionGC)) TRUE else round(fractionGC * nb$nGcSteps, 0) - nb$indexGC_offset
    rows[rows < 1] <- NA_integer_
    if(type[1] == 'adjustedPeak'){
        peak <- nb$model$readsPerAllele_peak[rows]
        mu   <- nb$model$readsPerAllele_mu[rows]
        ifelse(peak < peakThreshold, mu, peak) # revert to mu when peak is unreliable due to collision with zero
    } else {
        col <- paste('readsPerAllele', type[1], sep = '_')
        nb$model[[col]][rows]
    }
}

# calculate the z-score of a set of bin counts given a nbinomCountsGC model
# use mid-p approach similar to edgeR nbinomZScore
zScore.nbinomCountsGC <- function(nb, binCounts, fractionGC, binCN = 2){

    # get the expected bin count based on bin fraction GC and CN
    # expectedReadCount = rpa(reads per allele) * nAlleles
    mu <- predict.nbinomCountsGC(nb, suppressGcOutliers.nbinomCountsGC(nb, fractionGC), type = 'mu') * binCN 

    # get the probability on th higher and lower side of the observed bin counts
    binCounts <- as.integer(ceiling(binCounts + 0.001))
    x <- data.table(
        a = pnbinom(binCounts,      size = nb$theta, mu = mu),
        b = pnbinom(binCounts - 1L, size = nb$theta, mu = mu)
    )

    # project the probability to normal distribution-equivalent z-scores
    qnorm(x[, rowMeans(.SD)])
}
