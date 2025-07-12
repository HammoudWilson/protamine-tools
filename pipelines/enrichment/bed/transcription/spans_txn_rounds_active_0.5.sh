
# a higher RPKM threshold to consider a genome span as actively transcribed
export MIN_RPKM=0.5  # bin CPM is the same as RPKM for the input 1kb bins
export TXED_FILTER=1 # thus, selecting transcribed spans for output BED file
source transcription/txn_spans/txn_spans_active.sh
