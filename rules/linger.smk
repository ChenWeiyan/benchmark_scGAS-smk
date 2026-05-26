if "linger" in METHODS:
    rule run_linger:
        """
        LINGER: deep-learning gene regulatory network inference.
        Uses scMultiome (ATAC + RNA) to predict per-cell gene expression
        from chromatin accessibility.
        Requires Python LINGER package in its own conda env.
        """
        input:
            fragment  = config["input"]["fragment_file"],
            barcodes  = f"{OUTDIR}/qc/filtered_barcodes.txt",
            rna_counts= f"{OUTDIR}/qc/filtered_rna_counts.rds",
        output:
            activity  = f"{OUTDIR}/methods/linger/activity_lognorm.rds",
        params:
            genome    = config["linger"]["genome"],
            method    = config["linger"]["method"],
            n_cores   = config["resources"]["n_cores"],
            outdir    = f"{OUTDIR}/methods/linger",
        threads: config["resources"]["n_cores"]
        log:    f"{OUTDIR}/logs/linger.log"
        conda:  "../envs/linger.yaml"
        script: "../scripts/run_linger.py"
