#!/bin/bash

if [ -f ${GENOME_METADATA_PREFIX}.primary.tn5_site_freqs_exp.txt.gz ]; then
    echo "file already exists:"
    echo "    "${GENOME_METADATA_PREFIX}.primary.tn5_site_freqs_exp.txt.gz
    echo "skipping Tn5 site counting"
    exit
else
    Rscript ${ACTION_DIR}/sites.R
    checkPipe
fi
