if "cisformer" in METHODS:
    rule run_cisformer:
        """
        cisFormer: transformer-based prediction of gene expression from
        scATAC-seq chromatin accessibility.
        Requires Python cisFormer package in its own conda env.
        """
        input:
            fragment = config["input"]["fragment_file"],
            barcodes = f"{OUTDIR}/qc/filtered_barcodes.txt",
        output:
            activity = f"{OUTDIR}/methods/cisformer/activity_lognorm.rds",
        params:
            genome     = config["cisformer"]["genome"],
            batch_size = config["cisformer"]["batch_size"],
            seq_len    = config["cisformer"]["seq_len"],
            n_cores    = config["resources"]["n_cores"],
            outdir     = f"{OUTDIR}/methods/cisformer",
        threads: config["resources"]["n_cores"]
        log:    f"{OUTDIR}/logs/cisformer.log"
        conda:  "../envs/cisformer.yaml"
        script: "../scripts/run_cisformer.py"
