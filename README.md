# snakemake-sra-rnaseq

Takes a public SRA study accession and returns a gene-level count matrix, a QC report, and DESeq2 results.

Dependencies are pinned per rule in conda environments. Sample sheets are built by querying the SRA metadata tables in BigQuery instead of being assembled by hand. The DAG is dry-run tested in CI on every push.

<!-- TODO (Mason): two or three sentences here on why you built this and what you use it for.
     Written in your own words, this is the single most valuable paragraph in the file.
     Delete this comment when you've replaced it. -->

```
SRA accession
   ↓  prefetch + fasterq-dump
FASTQ
   ↓  fastp (adapter/quality trimming) → FastQC
trimmed FASTQ
   ↓  salmon (selective alignment, decoy-aware index)
transcript quantifications
   ↓  tximport (lengthScaledTPM)
gene count matrix
   ↓  DESeq2 (+ apeglm LFC shrinkage)
results tables, MA plots, PCA, MultiQC report
```

## Status

CI validates the DAG on every push: lint, plus a full dry run against a mixed
paired/single-end test sample sheet. That confirms the workflow is structurally
sound on a clean machine. It does not confirm the numbers, because nothing
executes in CI.

End-to-end validation against a published study is in progress. This section
will name the accession and report the recovered contrast once that's done.

## Quickstart

```bash
git clone https://github.com/masondandrus/snakemake-sra-rnaseq
cd snakemake-sra-rnaseq

# 1. Build the sample sheet from SRA metadata in BigQuery
python workflow/scripts/fetch_sra_metadata.py --study SRP123456 --audit
python workflow/scripts/fetch_sra_metadata.py --study SRP123456 --out config/samples.tsv

# 2. Point config/config.yaml at the right reference and set your contrasts

# 3. Check the plan, then run
snakemake -n
snakemake --use-conda --cores 16
```

Outputs land in `results/`: `counts/gene_counts.tsv`, `qc/multiqc_report.html`,
`deseq2/*.results.tsv`, and diagnostic PCA and MA plots.

For a small study you can skip step 1 and write `config/samples.tsv` by hand.
It's four columns.

## Building the sample sheet with SQL

SRA metadata lives in BigQuery (`nih-sra-datastore.sra.metadata`), so run
accessions, library layout, and submitter sample attributes can be pulled with a
query rather than clicked out of the Run Selector.
`resources/sra_metadata.sql` does that and normalises the result into the four
columns the workflow needs: `run`, `sample`, `condition`, `library_layout`.

The complication is that sample attributes are submitter-defined key/value
pairs, not fixed columns. What one study calls `treatment_sam` the next calls
`condition_sam` or `genotype_sam`, and some studies bury the design in
`source_name_sam` as prose. Run the audit query before anything else. It lists
every attribute key a study actually used, with example values, so you can edit
the `COALESCE` block to match what's really there.

```bash
python workflow/scripts/fetch_sra_metadata.py --study SRP123456 --audit
```

Get this step wrong and nothing downstream will tell you. A mislabelled group
survives every QC check in the pipeline and comes out the other end as a clean
volcano plot.

## Design decisions

**Selective alignment with a decoy-aware index.** Reads originating from
unannotated genomic sequence get force-assigned to transcripts they didn't come
from, inflating abundance estimates for a subset of them. Including the genome
as decoy sequence gives those reads somewhere else to go (Srivastava et al.,
2020). The cost is a genome-sized index build, roughly 32 GB of memory for
human.

**`countsFromAbundance = "lengthScaledTPM"`.** Raw TPM isn't on a count scale
and violates DESeq2's negative binomial assumptions. Raw counts ignore the fact
that average transcript length shifts between samples when isoform usage
changes. `lengthScaledTPM` corrects for length while keeping values countable.

**prefetch before fasterq-dump.** Streaming `fasterq-dump` straight from NCBI is
the usual cause of silently truncated FASTQ files on large studies. `prefetch`
is resumable, and the rule carries `retries: 3`.

**Layout-aware rule pairs instead of a checkpoint.** Paired- and single-end runs
go through separate rules whose `run` wildcard is scoped to the runs of that
layout, so the DAG resolves statically from the sample sheet. The two fastp
rules write identically named JSON reports, so output patterns alone left the
graph ambiguous. A `ruleorder` would have resolved that too, by silently sending
every single-end run down the paired-end branch.

**Intermediate FASTQs are `temp()`.** Raw and trimmed reads are deleted once
consumed. A 40-sample study stays in the tens of GB rather than the hundreds.

**One conda environment per rule.** The salmon and Bioconductor toolchains never
have to co-resolve, and conda isn't asked to solve one large environment, which
is where reproducibility tends to break a year later. Versions are pinned to the
minor release.

## Repository layout

```
config/
  config.yaml            reference URLs, tool parameters, design, contrasts
  samples.tsv            one row per SRA run
workflow/
  Snakefile              entry point and targets
  rules/                 common, download, qc, quant, de
  envs/                  pinned conda environments, one per toolchain
  scripts/               tx2gene.py, tximport.R, deseq2_init.R, deseq2_results.R,
                         fetch_sra_metadata.py, curate_samples.py
resources/
  sra_metadata.sql       BigQuery sample-sheet query + attribute audit query
docs/
  llm-curation.md        optional LLM-assisted curation helper
.github/workflows/ci.yml lint + dry-run on every push
```

## Configuration

| Key | Purpose |
| --- | --- |
| `reference.transcriptome_url` / `genome_url` | GENCODE FASTAs. Swap both together; the decoy index needs a matching genome build. |
| `salmon.kmer` | 31 by default. Drop to 23–25 for reads under ~50 bp. |
| `salmon.libtype` | `A` infers strandedness. Verify against `lib_format_counts.json`. |
| `deseq2.design` | Any DESeq2 formula over columns present in the sample sheet. |
| `deseq2.min_count` / `min_samples` | Pre-filtering threshold, independent of the design. |
| `contrasts` | `[factor, numerator, denominator]` per contrast. The key becomes the output filename. |

Non-human studies need both reference URLs changed. `tx2gene.py` parses
GENCODE-style pipe-delimited FASTA headers; an Ensembl FASTA needs a GTF-derived
map instead.

## Where to look when it finishes

Three outputs are worth checking before trusting any results table.

`results/salmon/*/lib_format_counts.json` reports the inferred library type per
sample. It should be consistent across the study. Split results usually mean
something is wrong with the deposit, not the pipeline.

`results/qc/multiqc_report.html` catches a collapsed mapping rate or a
duplication spike before it drags a group mean.

`results/deseq2/pca.png` is the important one. If PC1 separates by batch or
sequencing date rather than by condition, stop and add the covariate to
`deseq2.design`. Snakemake will rerun only the DESeq2 steps; quantification
stays cached.

## Running on a cluster

```bash
snakemake --use-conda --profile <your-profile> --jobs 100
```

Rules declare `threads` and `mem_mb`. The salmon index build sets the memory
ceiling.

## Testing

CI runs `snakemake --lint` and a dry run against `.test/config.yaml`, which
catches wildcard, input-function, and DAG errors without downloading anything.
The test sample sheet deliberately mixes paired- and single-end runs so both
branches are exercised. Same checks locally:

```bash
snakemake --lint
snakemake -n --configfile .test/config.yaml
```

## Optional tooling

`workflow/scripts/curate_samples.py` proposes `condition` assignments from
submitter free text using the Claude API, for the cases the SQL can't reach. It
runs outside the workflow and writes a proposal you review and commit. See
[docs/llm-curation.md](docs/llm-curation.md).

## Reference

Srivastava, A., Malik, L., Sarkar, H., Zakeri, M., Almodaresi, F., Soneson, C.,
Love, M. I., Kingsford, C., & Patro, R. (2020). Alignment and mapping
methodology influence transcript abundance estimation. *Genome Biology, 21*(1),
239. https://doi.org/10.1186/s13059-020-02151-8

## License

MIT
