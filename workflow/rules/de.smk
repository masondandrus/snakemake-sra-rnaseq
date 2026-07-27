"""Differential expression with DESeq2."""


rule deseq2_init:
    input:
        txi="results/counts/txi.rds",
        samples=config["samples"],
    output:
        dds="results/deseq2/dds.rds",
        vst="results/deseq2/vst.rds",
        pca=report("results/deseq2/pca.png", category="Diagnostics"),
    log:
        "logs/deseq2/init.log",
    params:
        design=config["deseq2"]["design"],
        min_count=config["deseq2"]["min_count"],
        min_samples=config["deseq2"]["min_samples"],
        pca_group=config["deseq2"]["pca_group"],
    conda:
        "../envs/deseq2.yaml"
    threads: 4
    resources:
        mem_mb=16000,
    script:
        "../scripts/deseq2_init.R"


rule deseq2_results:
    input:
        dds="results/deseq2/dds.rds",
    output:
        tsv="results/deseq2/{contrast}.results.tsv",
        ma=report("results/deseq2/{contrast}.ma.png", category="Differential expression"),
    log:
        "logs/deseq2/{contrast}.log",
    params:
        contrast=lambda w: config["contrasts"][w.contrast],
        alpha=config["deseq2"]["alpha"],
        shrink=config["deseq2"]["lfc_shrink"],
    conda:
        "../envs/deseq2.yaml"
    threads: 4
    resources:
        mem_mb=16000,
    script:
        "../scripts/deseq2_results.R"
