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
        tx2gene="resources/reference/tx2gene.tsv",
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


rule fgsea:
    """Rank-based enrichment. Uses the whole gene list, not just the
    significant subset, which is what makes it usable on small designs."""
    input:
        results="results/deseq2/{contrast}.results.tsv",
    output:
        tsv="results/enrichment/{contrast}.gsea.tsv",
        plot=report("results/enrichment/{contrast}.gsea.png", category="Enrichment"),
    log:
        "logs/fgsea/{contrast}.log",
    params:
        species=config["enrichment"]["species"],
        min_size=config["enrichment"]["min_size"],
        max_size=config["enrichment"]["max_size"],
    conda:
        "../envs/report.yaml"
    threads: 4
    resources:
        mem_mb=8000,
    script:
        "../scripts/fgsea.R"


rule report_html:
    input:
        dds="results/deseq2/dds.rds",
        vst="results/deseq2/vst.rds",
        results="results/deseq2/{contrast}.results.tsv",
        gsea="results/enrichment/{contrast}.gsea.tsv",
    output:
        report("results/report/{contrast}.report.html", category="Report"),
    log:
        "logs/report/{contrast}.log",
    params:
        contrast=lambda w: config["contrasts"][w.contrast],
        alpha=config["deseq2"]["alpha"],
    conda:
        "../envs/report.yaml"
    threads: 2
    resources:
        mem_mb=8000,
    script:
        "../scripts/report.Rmd"
