def _scgas_fragment(wildcards):
    if NEED_LIFTOVER:
        return f"{OUTDIR}/liftover/fragments_{SCGAS_GENOME}.tsv.gz"
    return config["input"]["fragment_file"]

def _scgas_inputs(wildcards):
    d = {
        "barcodes"   : f"{OUTDIR}/qc/filtered_barcodes.txt",
        "rna_counts" : f"{OUTDIR}/qc/filtered_rna_counts.rds",
        "fragment"   : _scgas_fragment(wildcards),
    }
    return d

if "scgas" in METHODS:
    rule run_scgas:
        """
        Full scGAS pipeline: preprocess → metacell → train Lasso models →
        compute single-cell GAS via network propagation.
        Output is a genes × cells matrix saved as RDS.
        """
        input:
            unpack(_scgas_inputs)
        output:
            activity = f"{OUTDIR}/methods/scgas/activity_lognorm.rds",
        params:
            ref_dir    = config["scgas"]["ref_dir"],
            genome     = config["scgas"]["ref_genome"],
            cfg        = config["scgas"],
            n_cores    = config["resources"]["n_cores"],
            chunk_size = config["scgas"]["chunk_size"],
            pkg_dir    = "../scGAS-R",   # path to local scGAS-R package (relative to smk-benchmark/)
        threads: config["resources"]["n_cores"]
        log:    f"{OUTDIR}/logs/scgas.log"
        conda:  "../envs/r_benchmark.yaml"
        script: "../scripts/run_scgas.R"
