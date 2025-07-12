
# a higher threshold for RPKM to consider a transcript as active
export MIN_RPKM=0.5 # same as CPM on the input 1kb bins
export TXED_FILTER=1
source transcription/txn_spans/txn_spans_active.sh
