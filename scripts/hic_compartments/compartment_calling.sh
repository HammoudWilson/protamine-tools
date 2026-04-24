#!/bin/bash
#SBATCH --job-name=compartment_calling
#SBATCH --account=hammou0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --output=output/round2/compartment_calling_%j.out
#SBATCH --error=error/round2/compartment_calling_%j.err
#SBATCH --mail-user=zapell@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=10G
#SBATCH --time=0:30:00
#SBATCH --profile=Task

# --- Argument parsing ---
usage() {
  echo "Usage: $0 -o OUTDIR -f HIC_FILE -n NORM -r RES"
  echo "  -o OUTDIR     Output directory (required)"
  echo "  -f HIC_FILE   Input .hic file (required)"
  echo "  -n NORM       Normalization: KR, VC, NONE (default: KR)"
  echo "  -r RES        Resolution in bp (default: 250000)"
  exit 1
}

NORM=KR
RES=250000

while getopts ":o:f:n:r:" opt; do
  case $opt in
    o) OUTDIR="$OPTARG" ;;
    f) HIC_FILE="$OPTARG" ;;
    n) NORM="$OPTARG" ;;
    r) RES="$OPTARG" ;;
    *) usage ;;
  esac
done

if [[ -z "$OUTDIR" || -z "$HIC_FILE" ]]; then
  usage
fi

JUICER_TOOLS=/nfs/turbo/umms-hammou/minjilab/juicer/scripts/juicer_tools.jar

mkdir -p "$OUTDIR"
chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y)

for chr in "${chroms[@]}"; do
  echo "Processing chr$chr..."
  java -Xmx8g -jar "$JUICER_TOOLS" eigenvector "$NORM" "$HIC_FILE" chr$chr BP "$RES" "$OUTDIR/chr${chr}_pc1_raw.txt" -p
  awk -v chr=chr$chr -v res=$RES '{
    start = (NR - 1) * res;
    end = start + res;
    if ($1 != "NaN") print chr"\t"start"\t"end"\t"$1;
  }' "$OUTDIR/chr${chr}_pc1_raw.txt" >> "$OUTDIR/compartments_${RES}_pc1.bedGraph"
  rm "$OUTDIR/chr${chr}_pc1_raw.txt"
done

echo "All Juicebox-compatible compartment annotation files saved in: $OUTDIR"