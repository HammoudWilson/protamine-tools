#!/bin/bash

mkdir -p ${PLOTS_DIR}

Rscript ${ACTION_DIR}/tss.R
checkPipe
