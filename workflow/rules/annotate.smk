"""Multi-evidence convergence over a config-declared set of evidence axes.

Every axis, whatever its origin, is normalised to the same contract: a TSV with
a gene_id and/or gene_symbol column plus metadata prefixed with the axis name,
and a sibling provenance JSON declaring the join key, the species, and whether
the axis is usable. convergence.R then consumes N of them without knowing what
any of them are, so adding a line of evidence is a config entry.

Dispatch by type uses the same wildcard-scoping pattern as the paired/single
FASTQ rules: each rule constrains {axis} to the axis names of its own type, so
the DAG resolves statically from config with no ruleorder.
"""

wildcard_constraints:
    axis=axis_constraint([a["name"] for a in AXES]),


# ---------------------------------------------------------------------------
# Internal axes, derived from this workflow's own outputs
# ---------------------------------------------------------------------------
rule axis_deseq2:
    input:
        src="results/deseq2/{contrast}.results.tsv",
    output:
        tsv="results/convergence/{contrast}.axis.{axis}.tsv",
        provenance="results/convergence/{contrast}.axis.{axis}.provenance.json",
    log:
        "logs/annotation/{contrast}.axis.{axis}.log",
    wildcard_constraints:
        axis=axis_constraint(axes_of_type("deseq2")),
    params:
        axis=lambda w: w.axis,
        type="deseq2",
        species=lambda w: AXIS_BY_NAME[w.axis].get("species") or ANNOT["species"],
        opts=lambda w: AXIS_BY_NAME[w.axis],
    conda:
        "../envs/annotate.yaml"
    script:
        "../scripts/axis_internal.R"


rule axis_wgcna:
    """Hub genes of the modules associated with this contrast. The association
    is fitted with the DESeq2 design formula and corrected across modules."""
    input:
        src="results/wgcna/modules.tsv",
        eigengenes="results/wgcna/eigengenes.tsv",
        src_provenance="results/wgcna/provenance.json",
        dds="results/deseq2/dds.rds",
    output:
        tsv="results/convergence/{contrast}.axis.{axis}.tsv",
        provenance="results/convergence/{contrast}.axis.{axis}.provenance.json",
        module_trait="results/wgcna/{contrast}.{axis}.module_trait.tsv",
    log:
        "logs/annotation/{contrast}.axis.{axis}.log",
    wildcard_constraints:
        axis=axis_constraint(axes_of_type("wgcna")),
    params:
        axis=lambda w: w.axis,
        type="wgcna",
        species=lambda w: AXIS_BY_NAME[w.axis].get("species") or ANNOT["species"],
        opts=lambda w: AXIS_BY_NAME[w.axis],
        contrast=lambda w: config["contrasts"][w.contrast],
        design=config["deseq2"]["design"],
        kme_threshold=WGCNA_OPTS.get("kme_threshold", 0.8),
        module_alpha=WGCNA_OPTS.get("module_alpha", 0.05),
    conda:
        "../envs/annotate.yaml"
    script:
        "../scripts/axis_internal.R"


rule wgcna:
    """Co-expression modules, eigengenes and kME, built once from all samples.
    Below wgcna.min_samples the rule writes empty tables and records why."""
    input:
        vst="results/deseq2/vst.rds",
    output:
        modules="results/wgcna/modules.tsv",
        eigengenes="results/wgcna/eigengenes.tsv",
        soft_threshold="results/wgcna/soft_threshold.tsv",
        plot=report("results/wgcna/modules.png", category="Co-expression"),
        provenance="results/wgcna/provenance.json",
    log:
        "logs/annotation/wgcna.log",
    params:
        min_samples=WGCNA_OPTS.get("min_samples", 15),
        n_top_genes=WGCNA_OPTS.get("n_top_genes", 5000),
        min_module_size=WGCNA_OPTS.get("min_module_size", 30),
        scale_free_r2=WGCNA_OPTS.get("scale_free_r2", 0.8),
        power=WGCNA_OPTS.get("power", "auto"),
    conda:
        "../envs/annotate.yaml"
    threads: 8
    resources:
        mem_mb=16000,
    script:
        "../scripts/wgcna.R"


# ---------------------------------------------------------------------------
# External axes, cached under resources/annotation/
# ---------------------------------------------------------------------------
rule axis_eqtl:
    """Region-matched cis-eQTL genes. GTEx by default; RatGTEx covers the brain
    regions the addiction literature uses (NAcc, IL, PL, OFC, LHb, pVTA, RMTg,
    BLA) and is the right choice for rat work."""
    output:
        tsv="resources/annotation/axis.{axis}.tsv",
        provenance="resources/annotation/axis.{axis}.provenance.json",
    log:
        "logs/annotation/axis.{axis}.log",
    wildcard_constraints:
        axis=axis_constraint(axes_of_type("eqtl")),
    params:
        script=workflow.source_path("../scripts/fetch_eqtl.py"),
        source=lambda w: AXIS_BY_NAME[w.axis].get("source", "gtex"),
        region=lambda w: AXIS_BY_NAME[w.axis]["region"],
        qval=lambda w: AXIS_BY_NAME[w.axis].get("qval", 0.05),
        url=lambda w: cli_opt("--url", AXIS_BY_NAME[w.axis].get("url")),
        species=lambda w: cli_opt("--species", AXIS_BY_NAME[w.axis].get("species")),
    conda:
        "../envs/annotate.yaml"
    retries: 2
    shell:
        "python {params.script:q} "
        "--axis {wildcards.axis} --source {params.source} "
        "--region {params.region:q} --qval {params.qval} {params.url} {params.species} "
        "--out {output.tsv} --provenance {output.provenance} > {log} 2>&1"


rule axis_gwas:
    """Genes mapped to trait-associated variants in the NHGRI-EBI GWAS Catalog.
    Proximity-based, not fine-mapped — see docs/convergence.md."""
    output:
        tsv="resources/annotation/axis.{axis}.tsv",
        provenance="resources/annotation/axis.{axis}.provenance.json",
    log:
        "logs/annotation/axis.{axis}.log",
    wildcard_constraints:
        axis=axis_constraint(axes_of_type("gwas")),
    params:
        script=workflow.source_path("../scripts/fetch_gwas.py"),
        traits=lambda w: AXIS_BY_NAME[w.axis]["traits"],
        pvalue=lambda w: AXIS_BY_NAME[w.axis].get("pvalue", 5e-8),
        species=lambda w: cli_opt("--species", AXIS_BY_NAME[w.axis].get("species")),
    conda:
        "../envs/annotate.yaml"
    retries: 2
    shell:
        "python {params.script:q} "
        "--axis {wildcards.axis} --traits {params.traits:q} "
        "--pvalue {params.pvalue} {params.species} "
        "--out {output.tsv} --provenance {output.provenance} > {log} 2>&1"


rule axis_file:
    """Any gene list you supply: differential ATAC peaks, CUT&RUN targets, a
    druggable-genome list, cell-type markers, hits from a previous study."""
    output:
        tsv="resources/annotation/axis.{axis}.tsv",
        provenance="resources/annotation/axis.{axis}.provenance.json",
    log:
        "logs/annotation/axis.{axis}.log",
    wildcard_constraints:
        axis=axis_constraint(axes_of_type("file")),
    params:
        script=workflow.source_path("../scripts/axis_from_file.py"),
        path=lambda w: AXIS_BY_NAME[w.axis]["path"],
        species=lambda w: AXIS_BY_NAME[w.axis].get("species") or ANNOT["species"],
    conda:
        "../envs/annotate.yaml"
    shell:
        "python {params.script:q} "
        "--axis {wildcards.axis} --path {params.path:q} "
        "--species {params.species:q} "
        "--out {output.tsv} --provenance {output.provenance} > {log} 2>&1"


# ---------------------------------------------------------------------------
# Convergence
# ---------------------------------------------------------------------------
rule convergence:
    input:
        universe="results/deseq2/{contrast}.results.tsv",
        axes=axis_inputs,
        provenance=axis_provenance,
    output:
        table="results/convergence/{contrast}.convergence.tsv",
        plot=report("results/convergence/{contrast}.upset.png", category="Convergence"),
        provenance="results/convergence/{contrast}.provenance.json",
    log:
        "logs/annotation/convergence.{contrast}.log",
    params:
        species=ANNOT.get("species"),
        orthologs=ANNOT.get("orthologs", True),
        require=ANNOT.get("require", []),
    conda:
        "../envs/annotate.yaml"
    threads: 2
    resources:
        mem_mb=8000,
    script:
        "../scripts/convergence.R"
