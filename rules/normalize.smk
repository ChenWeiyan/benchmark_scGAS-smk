NORM_METHODS = [m for m in METHODS if m in ["signac", "archr", "cicero"]]

if len(NORM_METHODS) > 0:
    rule normalize_methods:
        """
        Add each raw gene activity matrix as a Seurat assay and normalize
        using Seurat NormalizeData (log-CP10k).  Produces activity_lognorm.rds
        for each count-based method (signac, archr, cicero).
        """
        input:
            expand(f"{OUTDIR}/methods/{{method}}/activity_raw.rds", method=NORM_METHODS),
        output:
            expand(f"{OUTDIR}/methods/{{method}}/activity_lognorm.rds", method=NORM_METHODS),
        params:
            methods     = NORM_METHODS,
            methods_dir = f"{OUTDIR}/methods",
        log:   f"{OUTDIR}/logs/normalize_methods.log"
        conda: "../envs/r_benchmark.yaml"
        script: "../scripts/normalize_methods.R"
