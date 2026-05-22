import os

configfile: "config.yaml"

OUTDIR       = config["outdir"]
METHODS      = config["methods"]
GENOME       = config["input"]["genome"]
SCGAS_GENOME = config["scgas"]["ref_genome"]

# Liftover is needed when scgas is active and input genome differs from scGAS ref genome
NEED_LIFTOVER = "scgas" in METHODS and GENOME != SCGAS_GENOME

# Cluster benchmark is enabled when a cell-type CSV is provided
CELL_TYPE_CSV = config["evaluation"].get("cell_type_csv", "")

include: "rules/preprocess.smk"
include: "rules/liftover.smk"
include: "rules/scgas.smk"
include: "rules/signac.smk"
include: "rules/archr.smk"
include: "rules/cicero.smk"
include: "rules/linger.smk"
include: "rules/cisformer.smk"
include: "rules/scarlink.smk"
include: "rules/normalize.smk"
include: "rules/evaluate.smk"
include: "rules/cluster.smk"


rule all:
    input:
        f"{OUTDIR}/evaluation/pcc_comparison.pdf",
        *([f"{OUTDIR}/evaluation/cluster_metrics.pdf"] if CELL_TYPE_CSV else []),
