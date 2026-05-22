#!/usr/bin/env bash
# Run the scGAS benchmark pipeline.
# Usage (from smk-benchmark/):
#   bash run_benchmark.sh           # 8 cores
#   bash run_benchmark.sh 16        # 16 cores
#   bash run_benchmark.sh 8 -n      # dry-run

set -euo pipefail
cd "$(dirname "$0")"

CORES=${1:-8}
shift || true   # remaining args passed through

echo "=== scGAS benchmark pipeline ==="
echo "Working dir : $(pwd)"
echo "Cores       : $CORES"
echo ""

snakemake \
  --use-conda \
  --cores "$CORES" \
  --rerun-incomplete \
  --printshellcmds \
  "$@"
