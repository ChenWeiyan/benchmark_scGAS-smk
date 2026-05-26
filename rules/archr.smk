if "archr" in METHODS:
    rule run_archr:
        """
        ArchR GeneScoreMatrix: regression-based gene score from ATAC tiles.
        Requires ArchR >= 1.0.2.
        """
        input:
            fragment = config["input"]["fragment_file"],
            barcodes = f"{OUTDIR}/qc/filtered_barcodes.txt",
        output:
            activity = f"{OUTDIR}/methods/archr/activity_raw.rds",
        params:
            genome    = config["input"]["genome"],
            cfg       = config["archr"],
            n_cores   = config["resources"]["n_cores"],
            archr_dir = f"{OUTDIR}/methods/archr/ArchRProject",
        threads: config["resources"]["n_cores"]
        log:    f"{OUTDIR}/logs/archr.log"
        conda:  "../envs/r_benchmark.yaml"
        script: "../scripts/run_archr.R"
