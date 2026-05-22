"""
cisFormer: Transformer-based prediction of gene expression from scATAC-seq.
Paper: https://doi.org/10.1038/s41467-023-43204-1

Install: pip install cisformer   (or follow paper GitHub)

TODO: Verify the cisFormer package name and API.
      The implementation below is a structural template.
"""

import os
import sys
import numpy as np

outdir   = snakemake.params["outdir"]
out_path = snakemake.output["activity"]
os.makedirs(outdir, exist_ok=True)

log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh

fragment_path = snakemake.input["fragment"]
barcodes_path = snakemake.input["barcodes"]
genome        = snakemake.params["genome"]
batch_size    = snakemake.params["batch_size"]
seq_len       = snakemake.params["seq_len"]
n_cores       = snakemake.params["n_cores"]

with open(barcodes_path) as f:
    barcodes = [line.strip() for line in f]
print(f"[cisformer] Cells: {len(barcodes)}", file=sys.stderr)

# ── Parse fragment file into per-cell accessibility matrix ────────────────────
# TODO: cisFormer may require a sequence-level accessibility representation.
# This step needs to be adapted based on the actual cisFormer input format.

try:
    import cisformer   # TODO: confirm correct package name

    print(f"[cisformer] Running cisFormer (genome={genome}) ...", file=sys.stderr)

    # TODO: Replace with actual cisFormer workflow:
    # 1. Build accessibility profile from fragment file
    # 2. Run model inference to predict gene expression
    # 3. Extract per-cell per-gene prediction scores

    # model = cisformer.load_pretrained(genome=genome)
    # atac_input = cisformer.prepare_input(fragment_path, cells=barcodes, seq_len=seq_len)
    # predictions = model.predict(atac_input, batch_size=batch_size)
    # gene_names, activity_mat = predictions.genes, predictions.matrix   # genes × cells

    raise NotImplementedError("cisFormer API not yet wired — see TODO comments")

except ImportError:
    raise ImportError(
        "[cisformer] Cannot import cisformer. "
        "Install from the cisFormer GitHub repository."
    )

# ── Save as RDS via rpy2 ──────────────────────────────────────────────────────
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
numpy2ri.activate()
base = importr("base")

act_r = ro.r["matrix"](
    ro.FloatVector(activity_mat.flatten()),
    nrow=activity_mat.shape[0],
    dimnames=ro.r["list"](
        ro.StrVector(list(gene_names)),
        ro.StrVector(barcodes)
    )
)
base.saveRDS(act_r, out_path)
print(f"[cisformer] Saved: {out_path}", file=sys.stderr)
log_fh.close()
