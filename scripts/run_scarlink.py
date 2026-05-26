"""
SCARlink: Peak-to-gene link inference from scMultiome co-accessibility.
https://github.com/snehamitra/SCARlink

Install: pip install SCARlink

Gene activity = weighted sum of peak accessibilities using SCARlink link weights.
"""

import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp

outdir   = snakemake.params["outdir"]
out_path = snakemake.output["activity"]
os.makedirs(outdir, exist_ok=True)

log_fh = open(snakemake.log[0], "w")
sys.stderr = log_fh

fragment_path  = snakemake.input["fragment"]
barcodes_path  = snakemake.input["barcodes"]
rna_counts_path= snakemake.input["rna_counts"]
genome         = snakemake.params["genome"]
n_neighbors    = snakemake.params["n_neighbors"]
n_components   = snakemake.params["n_components"]
n_cores        = snakemake.params["n_cores"]

with open(barcodes_path) as f:
    barcodes = [line.strip() for line in f]
print(f"[scarlink] Cells: {len(barcodes)}", file=sys.stderr)

# ── Load RNA via rpy2 ─────────────────────────────────────────────────────────
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri, pandas2ri
numpy2ri.activate(); pandas2ri.activate()
base_r   = importr("base")
Matrix_r = importr("Matrix")

rna_r     = base_r.readRDS(rna_counts_path)
rna_dense = np.array(Matrix_r.as_matrix(rna_r))
rna_genes = list(ro.r["rownames"](rna_r))
rna_cells = list(ro.r["colnames"](rna_r))
shared    = [c for c in barcodes if c in rna_cells]
rna_idx   = [rna_cells.index(c) for c in shared]
rna_mat   = rna_dense[:, rna_idx]   # genes × shared cells

try:
    import scarlink
    from scarlink.src.scarlink import SCARlink

    print(f"[scarlink] Fitting SCARlink model ...", file=sys.stderr)

    # TODO: SCARlink takes an AnnData-like object with ATAC + RNA layers.
    # The block below is a structural template:

    import anndata as ad
    # Build ATAC peak × cell matrix from fragment file (e.g. using pybedtools or pyranges)
    # atac_mat: peaks × cells sparse matrix
    # peak_names: list of "chr:start-end" strings

    # adata = ad.AnnData(X=atac_mat.T)   # cells × peaks
    # adata.obsm["X_lsi"] = ...
    # adata.uns["rna"] = rna_mat.T       # cells × genes
    # adata.uns["gene_names"] = rna_genes

    # model = SCARlink(adata, n_neighbors=n_neighbors, n_components=n_components)
    # model.fit()
    # links = model.get_links()            # DataFrame: peak, gene, weight

    # Gene activity = dot(atac_mat, link_weight_matrix)
    # activity_mat: genes × cells

    raise NotImplementedError("SCARlink ATAC matrix construction not yet wired — see TODO")

except ImportError:
    raise ImportError(
        "[scarlink] Cannot import scarlink. "
        "Install: pip install SCARlink  (https://github.com/snehamitra/SCARlink)"
    )

# ── Save as RDS ───────────────────────────────────────────────────────────────
act_r = ro.r["matrix"](
    ro.FloatVector(activity_mat.flatten()),
    nrow=activity_mat.shape[0],
    dimnames=ro.r["list"](
        ro.StrVector(list(gene_names)),
        ro.StrVector(shared)
    )
)
base_r.saveRDS(act_r, out_path)
print(f"[scarlink] Saved: {out_path}", file=sys.stderr)
log_fh.close()
