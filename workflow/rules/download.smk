"""Retrieval of raw reads from SRA and of the reference transcriptome.

fasterq-dump is run through prefetch first: prefetch is resumable and far more
reliable over long transfers than streaming fasterq-dump directly, which is the
usual cause of silent truncation on large studies.
"""


rule prefetch:
    output:
        temp(directory("resources/sra/{run}")),
    log:
        "logs/prefetch/{run}.log",
    params:
        max_size=config["sra"]["max_size"],
    conda:
        "../envs/sra-tools.yaml"
    retries: 3
    threads: 1
    resources:
        runtime=120,
    shell:
        "prefetch {wildcards.run} "
        "--max-size {params.max_size} "
        "--output-directory resources/sra "
        "> {log} 2>&1"


rule fasterq_dump_pe:
    wildcard_constraints:
        run=constrain(PE_RUNS),
    input:
        "resources/sra/{run}",
    output:
        r1=temp("results/fastq/{run}_1.fastq.gz"),
        r2=temp("results/fastq/{run}_2.fastq.gz"),
    log:
        "logs/fasterq_dump/{run}.log",
    conda:
        "../envs/sra-tools.yaml"
    threads: 6
    resources:
        runtime=120,
    shell:
        """
        (
        mkdir -p results/scratch && tmp=$(mktemp -d -p results/scratch)
        fasterq-dump {wildcards.run} \
            --split-3 \
            --threads {threads} \
            --temp "$tmp" \
            --outdir "$tmp" \
            "{input}/{wildcards.run}.sra"
        pigz -p {threads} -c "$tmp/{wildcards.run}_1.fastq" > {output.r1}
        pigz -p {threads} -c "$tmp/{wildcards.run}_2.fastq" > {output.r2}
        rm -rf "$tmp"
        ) > {log} 2>&1
        """


rule fasterq_dump_se:
    wildcard_constraints:
        run=constrain(SE_RUNS),
    input:
        "resources/sra/{run}",
    output:
        temp("results/fastq/{run}.fastq.gz"),
    log:
        "logs/fasterq_dump/{run}.log",
    conda:
        "../envs/sra-tools.yaml"
    threads: 6
    resources:
        runtime=120,
    shell:
        """
        (
        mkdir -p results/scratch && tmp=$(mktemp -d -p results/scratch)
        fasterq-dump {wildcards.run} \
            --split-3 \
            --threads {threads} \
            --temp "$tmp" \
            --outdir "$tmp" \
            "{input}/{wildcards.run}.sra"
        pigz -p {threads} -c "$tmp/{wildcards.run}.fastq" > {output}
        rm -rf "$tmp"
        ) > {log} 2>&1
        """


rule get_transcriptome:
    output:
        "resources/reference/transcriptome.fa.gz",
    log:
        "logs/reference/transcriptome.log",
    params:
        url=config["reference"]["transcriptome_url"],
    conda:
        "../envs/base.yaml"
    retries: 3
    shell:
        "curl -fsSL {params.url} -o {output} 2> {log}"


rule get_genome:
    """Decoy sequence for selective alignment (Srivastava et al. 2020)."""
    output:
        "resources/reference/genome.fa.gz",
    log:
        "logs/reference/genome.log",
    params:
        url=config["reference"]["genome_url"],
    conda:
        "../envs/base.yaml"
    retries: 3
    shell:
        "curl -fsSL {params.url} -o {output} 2> {log}"
