if "signac" in METHODS:
    rule run_signac:
        """
        Signac GeneActivity: count fragments overlapping gene bodies
        (extended upstream by `upstream` bp).
        """
        input:
            fragment = config["input"]["fragment_file"],
            barcodes = f"{OUTDIR}/qc/filtered_barcodes.txt",
        output:
            activity = f"{OUTDIR}/methods/signac/activity_raw.rds",
        params:
            genome     = config["input"]["genome"],
            upstream   = config["signac"]["upstream"],
            downstream = config["signac"]["downstream"],
        threads: config["resources"]["n_cores"]
        log:    f"{OUTDIR}/logs/signac.log"
        conda:  "../envs/r_benchmark.yaml"
        script: "../scripts/run_signac.R"
