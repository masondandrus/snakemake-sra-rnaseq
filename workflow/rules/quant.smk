"""Transcript quantification with salmon, summarised to gene level with tximport."""


rule salmon_decoys:
    """Build the decoy-aware index input: transcriptome + whole genome as decoy."""
    input:
        txome="resources/reference/transcriptome.fa.gz",
        genome="resources/reference/genome.fa.gz",
    output:
        gentrome=temp("resources/reference/gentrome.fa.gz"),
        decoys=temp("resources/reference/decoys.txt"),
    log:
        "logs/salmon/decoys.log",
    conda:
        "../envs/base.yaml"
    threads: 2
    shell:
        """
        (
        zcat {input.genome} | grep '^>' | cut -d ' ' -f 1 | sed 's/>//g' > {output.decoys}
        cat {input.txome} {input.genome} > {output.gentrome}
        ) > {log} 2>&1
        """


rule salmon_index:
    input:
        gentrome="resources/reference/gentrome.fa.gz",
        decoys="resources/reference/decoys.txt",
    output:
        directory("resources/reference/salmon_index"),
    log:
        "logs/salmon/index.log",
    params:
        k=config["salmon"]["kmer"],
    conda:
        "../envs/salmon.yaml"
    threads: 8
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "salmon index "
        "--transcripts {input.gentrome} "
        "--decoys {input.decoys} "
        "--index {output} "
        "--kmerLen {params.k} "
        "--gencode "
        "--threads {threads} > {log} 2>&1"


rule salmon_quant_pe:
    wildcard_constraints:
        run=constrain(PE_RUNS),
    input:
        index="resources/reference/salmon_index",
        r1="results/trimmed/{run}_1.fastq.gz",
        r2="results/trimmed/{run}_2.fastq.gz",
    output:
        quant="results/salmon/{run}/quant.sf",
        lib="results/salmon/{run}/lib_format_counts.json",
    log:
        "logs/salmon/quant/{run}.log",
    params:
        libtype=config["salmon"]["libtype"],
        extra=config["salmon"]["extra"],
        outdir=lambda w: f"results/salmon/{w.run}",
    conda:
        "../envs/salmon.yaml"
    threads: 8
    resources:
        mem_mb=16000,
    shell:
        "salmon quant "
        "--index {input.index} "
        "--libType {params.libtype} "
        "-1 {input.r1} -2 {input.r2} "
        "--output {params.outdir} "
        "--threads {threads} {params.extra} > {log} 2>&1"


rule salmon_quant_se:
    wildcard_constraints:
        run=constrain(SE_RUNS),
    input:
        index="resources/reference/salmon_index",
        fq="results/trimmed/{run}.fastq.gz",
    output:
        quant="results/salmon/{run}/quant.sf",
        lib="results/salmon/{run}/lib_format_counts.json",
    log:
        "logs/salmon/quant/{run}.log",
    params:
        libtype=config["salmon"]["libtype"],
        extra=config["salmon"]["extra"],
        outdir=lambda w: f"results/salmon/{w.run}",
    conda:
        "../envs/salmon.yaml"
    threads: 8
    resources:
        mem_mb=16000,
    shell:
        "salmon quant "
        "--index {input.index} "
        "--libType {params.libtype} "
        "-r {input.fq} "
        "--output {params.outdir} "
        "--threads {threads} {params.extra} > {log} 2>&1"


rule tx2gene:
    """Transcript-to-gene map parsed from GENCODE-style FASTA headers.

    Header format: >ENST...|ENSG...|OTTHUMG...|...|gene_symbol|length|type|
    If you point the config at an Ensembl (non-GENCODE) FASTA, swap this rule
    for a GTF-derived map instead — the header fields differ.
    """
    input:
        "resources/reference/transcriptome.fa.gz",
    output:
        "resources/reference/tx2gene.tsv",
    log:
        "logs/reference/tx2gene.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/tx2gene.py"


rule tximport:
    input:
        quants=expand("results/salmon/{run}/quant.sf", run=RUNS),
        tx2gene="resources/reference/tx2gene.tsv",
        samples=config["samples"],
    output:
        counts="results/counts/gene_counts.tsv",
        tpm="results/counts/gene_tpm.tsv",
        rds="results/counts/txi.rds",
    log:
        "logs/tximport.log",
    conda:
        "../envs/deseq2.yaml"
    threads: 1
    resources:
        mem_mb=8000,
    script:
        "../scripts/tximport.R"
