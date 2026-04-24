#!/bin/bash
#SBATCH --account=hammou0
#SBATCH --partition=standard
#SBATCH --job-name=liftover_mm10_to_mm39
#SBATCH --output=output/round2/liftover_%j.out
#SBATCH --error=error/round2/liftover_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# Usage: sbatch liftover_mm10_to_mm39.sbatch <input.bed> <output.bed> <unmapped.bed>

#set -euo pipefail

INPUT=${1:-}
OUTPUT=${2:-}
UNMAPPED=${3:-}
CHAIN_URL="http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz"
CHAIN_FILE="mm10ToMm39.over.chain.gz"
LIFTOVER_BIN="liftOver"

if [[ -z "$INPUT" || -z "$OUTPUT" || -z "$UNMAPPED" ]]; then
    echo "Usage: sbatch liftover_mm10_to_mm39.sbatch <input.bed> <output.bed> <unmapped.bed>"
    exit 1
fi

# Download chain file if not present
if [[ ! -f "$CHAIN_FILE" ]]; then
    echo "Downloading $CHAIN_FILE ..."
    wget "$CHAIN_URL" -O "$CHAIN_FILE"
fi

# Check for liftOver binary
if ! command -v "$LIFTOVER_BIN" &> /dev/null; then
    echo "Error: liftOver binary not found in PATH. Please load the appropriate module or add it to your PATH."
    exit 2
fi

echo "Running liftOver..."
"$LIFTOVER_BIN" "$INPUT" "$CHAIN_FILE" "$OUTPUT" "$UNMAPPED"
echo "Done. Output: $OUTPUT, Unmapped: $UNMAPPED"
