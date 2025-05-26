#!/bin/bash

export BIAS_PLOT_DIR=${PLOTS_DIR}/bias
mkdir -p ${BIAS_PLOT_DIR}

Rscript ${ACTION_DIR}/bias.R
checkPipe
