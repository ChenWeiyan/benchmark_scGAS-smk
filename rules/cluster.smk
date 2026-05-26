if CELL_TYPE_CSV:
    rule cluster_benchmark:
        """
        Per-method clustering benchmark: activity_lognorm → PCA → Leiden
        resolution sweep → ARI / NMI / AMI vs cell-type ground truth.
        Requires evaluation.cell_type_csv to be set in config.yaml.
        """
        input:
            activities    = expand(
                f"{OUTDIR}/methods/{{method}}/activity_lognorm.rds",
                method=METHODS
            ),
            cell_type_csv = CELL_TYPE_CSV,
        output:
            metrics_csv = f"{OUTDIR}/evaluation/cluster_metrics.csv",
            plot_pdf    = f"{OUTDIR}/evaluation/cluster_metrics.pdf",
            plot_png    = f"{OUTDIR}/evaluation/cluster_metrics.png",
        params:
            methods       = METHODS,
            methods_dir   = f"{OUTDIR}/methods",
            cell_type_csv = CELL_TYPE_CSV,
            cell_type_col = config["evaluation"]["cell_type_col"],
            resolutions   = config["clustering"]["resolutions"],
            n_pcs         = config["clustering"]["n_pcs"],
        log:   f"{OUTDIR}/logs/cluster_benchmark.log"
        conda: "../envs/r_benchmark.yaml"
        script: "../scripts/cluster_benchmark.R"
