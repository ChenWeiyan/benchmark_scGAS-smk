rule preprocess:
    """
    Shared QC step.  Reads fragment file + RNA, computes ATAC/RNA QC metrics,
    applies thresholds, and saves filtered barcodes + RNA matrices used by
    every downstream method.
    """
    input:
        fragment = config["input"]["fragment_file"],
        rna      = config["input"]["rna"],
    output:
        barcodes    = f"{OUTDIR}/qc/filtered_barcodes.txt",
        rna_counts  = f"{OUTDIR}/qc/filtered_rna_counts.rds",
        rna_lognorm = f"{OUTDIR}/qc/filtered_rna_lognorm.rds",
        qc_metrics  = f"{OUTDIR}/qc/qc_metrics.rds",
    params:
        genome  = config["input"]["genome"],
        qc      = config["qc"],
        n_cores = config["resources"]["n_cores"],
    threads: config["resources"]["n_cores"]
    log:     f"{OUTDIR}/logs/preprocess.log"
    conda:   "../envs/r_benchmark.yaml"
    script:  "../scripts/preprocess_shared.R"
