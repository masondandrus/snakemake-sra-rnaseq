"""Adapter/quality trimming and QC aggregation."""


rule fastp_pe:
    wildcard_constraints:
        run=constrain(PE_RUNS),
    input:
        r1="results/fastq/{run}_1.fastq.gz",
        r2="results/fastq/{run}_2.fastq.gz",
    output:
        r1="results/trimmed/{run}_1.fastq.gz",
        r2="results/trimmed/{run}_2.fastq.gz",
        json="results/qc/fastp/{run}.fastp.json",
        html="results/qc/fastp/{run}.fastp.html",
    log:
        "logs/fastp/{run}.log",
    params:
        extra=config["fastp"]["extra"],
    conda:
        "../envs/fastp.yaml"
    threads: 4
    shell:
        "fastp "
        "--in1 {input.r1} --in2 {input.r2} "
        "--out1 {output.r1} --out2 {output.r2} "
        "--detect_adapter_for_pe "
        "--json {output.json} --html {output.html} "
        "--thread {threads} {params.extra} "
        "> {log} 2>&1"


rule fastp_se:
    wildcard_constraints:
        run=constrain(SE_RUNS),
    input:
        "results/fastq/{run}.fastq.gz",
    output:
        fq="results/trimmed/{run}.fastq.gz",
        json="results/qc/fastp/{run}.fastp.json",
        html="results/qc/fastp/{run}.fastp.html",
    log:
        "logs/fastp/{run}.log",
    params:
        extra=config["fastp"]["extra"],
    conda:
        "../envs/fastp.yaml"
    threads: 4
    shell:
        "fastp "
        "--in1 {input} --out1 {output.fq} "
        "--json {output.json} --html {output.html} "
        "--thread {threads} {params.extra} "
        "> {log} 2>&1"


rule fastqc_pe:
    wildcard_constraints:
        run=constrain(PE_RUNS),
    input:
        "results/trimmed/{run}_{read}.fastq.gz",
    output:
        html="results/qc/fastqc/{run}_{read}_fastqc.html",
        zip="results/qc/fastqc/{run}_{read}_fastqc.zip",
    log:
        "logs/fastqc/{run}_{read}.log",
    conda:
        "../envs/qc.yaml"
    threads: 2
    shell:
        "fastqc --quiet --threads {threads} "
        "--outdir results/qc/fastqc {input} > {log} 2>&1"


rule fastqc_se:
    wildcard_constraints:
        run=constrain(SE_RUNS),
    input:
        "results/trimmed/{run}.fastq.gz",
    output:
        html="results/qc/fastqc/{run}_fastqc.html",
        zip="results/qc/fastqc/{run}_fastqc.zip",
    log:
        "logs/fastqc/{run}.log",
    conda:
        "../envs/qc.yaml"
    threads: 2
    shell:
        "fastqc --quiet --threads {threads} "
        "--outdir results/qc/fastqc {input} > {log} 2>&1"


rule multiqc:
    input:
        fastqc=all_fastqc_zips(),
        fastp=all_fastp_reports(),
        salmon=expand("results/salmon/{run}/quant.sf", run=RUNS),
    output:
        report("results/qc/multiqc_report.html", category="QC"),
    log:
        "logs/multiqc.log",
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc --force "
        "--outdir results/qc --filename multiqc_report.html "
        "results/qc results/salmon > {log} 2>&1"
