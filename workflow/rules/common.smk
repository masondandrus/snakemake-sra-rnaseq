"""Sample-sheet parsing, wildcard constraints, and input helper functions."""

import shlex

import pandas as pd


# ---------------------------------------------------------------------------
# Sample sheet
#
# samples.tsv is one row per SRA run. It can be written by hand or generated
# from the public SRA metadata table in BigQuery:
#     python workflow/scripts/fetch_sra_metadata.py --study SRP123456
# See resources/sra_metadata.sql for the query.
# ---------------------------------------------------------------------------
samples = pd.read_csv(config["samples"], sep="\t", dtype=str).set_index(
    "run", drop=False
)
samples.index.names = ["run_id"]

RUNS = samples["run"].tolist()
PE_RUNS = samples.loc[samples["library_layout"] == "PAIRED", "run"].tolist()
SE_RUNS = samples.loc[samples["library_layout"] == "SINGLE", "run"].tolist()

_unknown = set(RUNS) - set(PE_RUNS) - set(SE_RUNS)
if _unknown:
    raise ValueError(
        f"library_layout must be PAIRED or SINGLE; offending runs: {sorted(_unknown)}"
    )


def constrain(runs):
    """Regex matching exactly these runs (and nothing at all if the list is empty).

    Paired- and single-end rules write some identically-named outputs (the fastp
    JSON/HTML), so pattern uniqueness alone does not disambiguate them. Scoping
    each rule's `run` wildcard to the runs of that layout resolves the DAG
    statically, without a ruleorder that would silently prefer one branch.
    """
    return "|".join(runs) if runs else "^$"


wildcard_constraints:
    run=constrain(RUNS),
    read="1|2",
    contrast="|".join(config["contrasts"]) if config.get("contrasts") else "^$",


# ---------------------------------------------------------------------------
# Input helpers
# ---------------------------------------------------------------------------
def is_paired(run):
    return samples.loc[run, "library_layout"] == "PAIRED"


def trimmed_fastq(wildcards):
    """Trimmed reads for a run, layout-aware."""
    if is_paired(wildcards.run):
        return expand(
            "results/trimmed/{run}_{read}.fastq.gz", run=wildcards.run, read=["1", "2"]
        )
    return ["results/trimmed/{run}.fastq.gz".format(run=wildcards.run)]


def all_fastqc_zips():
    pe = expand("results/qc/fastqc/{run}_{read}_fastqc.zip", run=PE_RUNS, read=["1", "2"])
    se = expand("results/qc/fastqc/{run}_fastqc.zip", run=SE_RUNS)
    return pe + se


def all_fastp_reports():
    return expand("results/qc/fastp/{run}.fastp.json", run=RUNS)


# ---------------------------------------------------------------------------
# Evidence axes for the convergence module (see rules/annotate.smk)
# ---------------------------------------------------------------------------
ANNOT = config.get("annotation") or {}
ENRICH = config.get("enrichment") or {}
if "species" not in ENRICH:
    raise ValueError(
        "enrichment.species is required in the configuration, "
        'for example "Mus musculus" or "Homo sapiens"'
    )
AXES = ANNOT.get("axes", []) if ANNOT.get("enabled") else []
WGCNA_OPTS = ANNOT.get("wgcna") or {}
AXIS_BY_NAME = {a["name"]: a for a in AXES}
INTERNAL_TYPES = {"deseq2", "wgcna"}


def axes_of_type(*types):
    return [a["name"] for a in AXES if a["type"] in types]


def axis_constraint(names):
    return "|".join(names) if names else "^$"


def is_internal(name):
    return AXIS_BY_NAME[name]["type"] in INTERNAL_TYPES


def cli_opt(flag, value):
    """'--flag value' for a set value, an empty string otherwise."""
    return f"{flag} {shlex.quote(str(value))}" if value else ""


def axis_path(name, contrast, suffix="tsv"):
    """Contrast-specific axes live under results/; external resources are
    cached once under resources/ and reused across contrasts."""
    if is_internal(name):
        return f"results/convergence/{contrast}.axis.{name}.{suffix}"
    return f"resources/annotation/axis.{name}.{suffix}"


def axis_inputs(w):
    return [axis_path(a["name"], w.contrast) for a in AXES]


def axis_provenance(w):
    return [axis_path(a["name"], w.contrast, "provenance.json") for a in AXES]


