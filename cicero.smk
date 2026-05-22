if "cicero" in METHODS:
    rule run_cicero:
        """
        Cicero co-accessibility gene activity: links peaks to genes via
        co-accessibility scores, then aggregates peak accessibility per gene.
        Requires cicero >= 1.3.
        """
        input:
            fragment = config["input"]["fragment_file"],
            barcodes = f"{OUTDIR}/qc/filtered_barcodes.txt",
        output:
            activity = f"{OUTDIR}/methods/cicero/activity_raw.rds",
        params:
            genome  = config["input"]["genome"],
            cfg     = config["cicero"],
            n_cores = config["resources"]["n_cores"],
            outdir  = f"{OUTDIR}/methods/cicero",
        threads: config["resources"]["n_cores"]
        log:    f"{OUTDIR}/logs/cicero.log"
        conda:  "../envs/r_benchmark.yaml"
        script: "../scripts/run_cicero.R"
