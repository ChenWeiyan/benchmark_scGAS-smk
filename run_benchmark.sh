#!/usr/bin/env bash
# Run the scGAS benchmark pipeline.
#
# Setup (once):
#   1. Clone both repos side by side:
#        git clone git@github.com:ChenWeiyan/scGAS-R.git
#        git clone git@github.com:ChenWeiyan/benchmark_scGAS-smk.git
#   2. Create and activate the R environment:
#        mamba env create -f benchmark_scGAS-smk/envs/r_benchmark.yaml
#        conda activate r_benchmark
#   3. Install scGAS into that env:
#        Rscript -e 'devtools::install("scGAS-R")'
#      (If scGAS is not installed, run_scgas.R falls back to devtools::load_all()
#       on `scgas.pkg_dir` in config.yaml, default ../scGAS-R.)
#
# Usage (from benchmark_scGAS-smk/):
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
