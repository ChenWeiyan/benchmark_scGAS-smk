rule evaluate:
    """
    Compute per-gene Pearson correlation (PCC) between each method's activity
    matrix and RNA log-norm expression.  Three gene sets: all expressed genes,
    top-N HVG from RNA, and top-N DEG per cell type (requires cell_type_csv).
    Outputs a long-format CSV and a 4-panel faceted comparison plot.
    """
    input:
        rna_lognorm = f"{OUTDIR}/qc/filtered_rna_lognorm.rds",
        activities  = expand(
            f"{OUTDIR}/methods/{{method}}/activity_lognorm.rds",
            method=METHODS
        ),
    output:
        pcc_csv  = f"{OUTDIR}/evaluation/pcc_results.csv",
        plot_pdf = f"{OUTDIR}/evaluation/pcc_comparison.pdf",
        plot_png = f"{OUTDIR}/evaluation/pcc_comparison.png",
    params:
        methods       = METHODS,
        methods_dir   = f"{OUTDIR}/methods",
        min_expr_frac = config["evaluation"]["min_expr_frac"],
        n_hvg         = config["evaluation"]["n_hvg"],
        n_top_deg     = config["evaluation"]["n_top_deg"],
        cell_type_csv = config["evaluation"]["cell_type_csv"],
        cell_type_col = config["evaluation"]["cell_type_col"],
    log:    f"{OUTDIR}/logs/evaluate.log"
    conda:  "../envs/r_benchmark.yaml"
    script: "../scripts/evaluate.R"
