
# a lower RPKM threshold to consider a genome span as actively transcribed
export MIN_RPKM=0.1  # bin CPM is the same as RPKM for the input 1kb bins
export TXED_FILTER=0 # thus, selecting non-transcribed spans for output BED file
source transcription/txn_spans/txn_spans_active.sh
