if "scarlink" in METHODS:
    rule run_scarlink:
        """
        SCARlink: infers peak-to-gene links from scMultiome co-accessibility,
        then aggregates peak accessibility into gene-level activity scores.
        Requires Python SCARlink package in its own conda env.
        """
        input:
            fragment  = config["input"]["fragment_file"],
            barcodes  = f"{OUTDIR}/qc/filtered_barcodes.txt",
            rna_counts= f"{OUTDIR}/qc/filtered_rna_counts.rds",
        output:
            activity  = f"{OUTDIR}/methods/scarlink/activity_lognorm.rds",
        params:
            genome      = config["scarlink"]["genome"],
            n_neighbors = config["scarlink"]["n_neighbors"],
            n_components= config["scarlink"]["n_components"],
            n_cores     = config["resources"]["n_cores"],
            outdir      = f"{OUTDIR}/methods/scarlink",
        threads: config["resources"]["n_cores"]
        log:    f"{OUTDIR}/logs/scarlink.log"
        conda:  "../envs/scarlink.yaml"
        script: "../scripts/run_scarlink.py"
