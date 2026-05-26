"""
LINGER: Deep learning gene regulatory network inference.
https://github.com/SUwonglab/LINGER

Install: pip install lingerpy   (or follow GitHub README)

Outputs a genes × cells activity matrix (predicted gene expression from ATAC)
saved as an RDS file via rpy2 for compatibility with the R evaluate script.

TODO: Verify LINGER API against the installed version.
"""

import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp

outdir   = snakemake.params["outdir"]
out_path = snakemake.output["activity"]
os.makedirs(outdir, exist_ok=True)

# Redirect stderr to log
log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh

fragment_path  = snakemake.input["fragment"]
barcodes_path  = snakemake.input["barcodes"]
rna_counts_path= snakemake.input["rna_counts"]
genome         = snakemake.params["genome"]
method         = snakemake.params["method"]
n_cores        = snakemake.params["n_cores"]

# ── Load filtered barcodes ────────────────────────────────────────────────────
with open(barcodes_path) as f:
    barcodes = [line.strip() for line in f]
print(f"[linger] Cells: {len(barcodes)}", file=sys.stderr)

# ── Load RNA counts (from RDS via rpy2) ───────────────────────────────────────
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri, pandas2ri
numpy2ri.activate(); pandas2ri.activate()
base = importr("base"); Matrix = importr("Matrix")

rna_r     = base.readRDS(rna_counts_path)
rna_np    = np.array(Matrix.as_matrix(rna_r))
rna_genes = list(ro.r["rownames"](rna_r))
rna_cells = list(ro.r["colnames"](rna_r))
rna_df    = pd.DataFrame(rna_np, index=rna_genes, columns=rna_cells)
rna_df    = rna_df[barcodes] if all(b in rna_df.columns for b in barcodes) else rna_df

# ── Build ATAC peak matrix from fragment file ─────────────────────────────────
# TODO: LINGER may accept fragment file paths directly or a peak matrix.
# The following is a placeholder — adapt to the actual LINGER input format.
# See: https://github.com/SUwonglab/LINGER for the correct input API.

try:
    import linger
    from linger import LINGER

    print(f"[linger] Running LINGER (method={method}) ...", file=sys.stderr)

    # LINGER typically requires:
    #   1. scATAC peak × cell matrix
    #   2. scRNA gene × cell matrix
    #   3. bulk reference (handled internally by LINGER with genome)
    # TODO: replace with actual LINGER call once input format is confirmed
    model = LINGER(
        genome      = genome,
        method      = method,
        n_jobs      = n_cores,
        output_dir  = outdir,
    )
    # model.load_data(atac_matrix, rna_matrix, barcodes)
    # model.train()
    # activity_mat = model.predict_expression()   # genes × cells

    raise NotImplementedError("LINGER input loading not yet implemented — see TODO above")

except ImportError:
    raise ImportError(
        "[linger] Cannot import linger package. "
        "Install: pip install lingerpy  (or follow https://github.com/SUwonglab/LINGER)"
    )

# ── Save as RDS ───────────────────────────────────────────────────────────────
# activity_mat: numpy array (genes × cells) or scipy sparse
# Convert to R sparse matrix and save as RDS
act_r = ro.r["matrix"](
    ro.FloatVector(activity_mat.flatten()),
    nrow=activity_mat.shape[0],
    dimnames=ro.r["list"](
        ro.StrVector(gene_names),
        ro.StrVector(barcodes)
    )
)
base.saveRDS(act_r, out_path)
print(f"[linger] Saved: {out_path}", file=sys.stderr)
log_fh.close()
