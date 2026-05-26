if NEED_LIFTOVER:
    rule liftover:
        """
        LiftOver the input fragment file to the scGAS reference genome (hg19).
        Produces a bgzipped, tabix-indexed fragment file.
        Only triggered when input genome != scGAS ref genome.
        """
        input:
            fragment = config["input"]["fragment_file"],
            chain    = config["scgas"]["chain_file"],
        output:
            lifted = f"{OUTDIR}/liftover/fragments_{SCGAS_GENOME}.tsv.gz",
            tbi    = f"{OUTDIR}/liftover/fragments_{SCGAS_GENOME}.tsv.gz.tbi",
        params:
            n_cores = config["resources"]["n_cores"],
            outdir  = f"{OUTDIR}/liftover",
        threads: config["resources"]["n_cores"]
        log:    f"{OUTDIR}/logs/liftover.log"
        conda:  "../envs/r_benchmark.yaml"
        script: "../scripts/liftover.R"
