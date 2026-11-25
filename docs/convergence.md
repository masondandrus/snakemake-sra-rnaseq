# Multi-evidence convergence

Optional. The workflow runs without it — set `annotation.enabled: false`.

A differential expression run on a whole transcriptome gives you a ranked list
and no obvious place to start. This module narrows it by asking which genes are
named by more than one kind of evidence.

The axes are declared in config, not hardcoded. Five types ship:

| Type | Source | Question |
| --- | --- | --- |
| `deseq2` | this workflow's DE result | Does expression change with condition? |
| `wgcna` | co-expression on the same counts | Is it a hub of a module associated with this contrast? |
| `eqtl` | GTEx or RatGTEx, region-matched | Is its expression under genetic control in this tissue? |
| `gwas` | NHGRI-EBI GWAS Catalog | Is it mapped to a variant associated with the trait? |
| `file` | anything you supply | Whatever you can express as a gene list |

Each axis is normalised to the same contract — a TSV keyed on `gene_id` and/or
`gene_symbol` with metadata prefixed by the axis name, plus a provenance JSON
declaring the key, the species, and whether the axis is usable. `convergence.R`
consumes N of them without knowing what any of them are.

You can declare the same type more than once. Two brain regions is two `eqtl`
axes; a differential ATAC list and a CUT&RUN target list is two `file` axes.

Output is `results/convergence/{contrast}.convergence.tsv`: one row per gene
with `n_evidence`, an `evidence` string like `DEG|eQTL|GWAS`, the supporting
statistics from every axis, and per-axis flags showing which hits required an
ortholog mapping. Plus an UpSet plot of how the sets intersect.

`annotation.require` is a hard filter — genes must hit those axes to appear at
all. `require: ["DEG"]` restricts the table to differentially expressed genes
and asks what else is true of them.

## What the score is not

`n_evidence` is a triage device. It is not a statistic, and three specific
things follow:

**The axes are not independent.** Expressed genes are likelier to have detected
eQTLs. GWAS mapped genes skew toward gene-dense regions. Counting hits does not
multiply independent probabilities and no p-value is emitted for the count.

**eQTL/GWAS overlap is not colocalisation.** When both axes name a gene, it
means both analyses named that gene — not that the eQTL and the GWAS signal
share a causal variant. Establishing that needs coloc, SMR, or fine-mapping,
none of which this module does. This is the claim the framework is easiest to
misread as making, so it is worth being explicit about in any write-up.

**GWAS "mapped gene" is mostly proximity.** The Catalog assigns genes to
variants largely by position. A gene on this axis has a trait-associated
variant in or near it.

## Sample size and the hub axis

WGCNA needs enough samples for the correlation structure to mean something —
the package authors put the floor at 15 and recommend 20 or more. Below
`annotation.wgcna.min_samples` the rule writes an empty hub table, records
`skipped: true` with the reason, and exits cleanly. Convergence then runs on
the remaining axes and reports how many were available.

This is deliberate. Modules estimated from six samples look exactly as
plausible as real ones and reproduce nothing.

The network is built once from all samples. For each contrast, every module
eigengene is then regressed on the DESeq2 design formula, and the coefficient
for the contrast's numerator against its denominator is tested. P-values are
corrected across modules (Benjamini-Hochberg). A gene enters the hub axis if
its module passes `annotation.wgcna.module_alpha` and its kME reaches
`annotation.wgcna.kme_threshold`. `results/wgcna/{contrast}.{axis}.module_trait.tsv`
holds the estimate, p-value and adjusted p-value of every module.

The soft-threshold power is the lowest power whose scale-free fit reaches
`annotation.wgcna.scale_free_r2`. `results/wgcna/soft_threshold.tsv` holds the
fit at every power tried. When no power passes, the rule takes the value that
the WGCNA FAQ gives for a signed network at that number of samples, and the
provenance file says so.

## Species and orthologs

The DE data and the annotation resources are usually different species: mouse
or rat expression, rat eQTLs, human GWAS. `g:Orth` maps DE gene IDs into each
resource's species, and every hit that depended on a mapping is flagged in the
output (`{axis}_via_ortholog` for each axis, and `any_ortholog`), because
ortholog assignment is lossy in both directions and a reviewer should be able
to see which calls rest on it.

Set `annotation.orthologs: false` to require exact ID matching instead, which
is right when the data and the resource are the same species.

## Adding an axis

Anything expressible as a gene list becomes an axis without touching workflow
code. Write a TSV with a `gene_id` and/or `gene_symbol` column — extra columns
carry through as metadata — and declare it:

```yaml
    - name: ATAC
      type: file
      path: "resources/annotation/atac_differential_peaks.tsv"
      species: "Mus musculus"
```

Validation is strict on purpose: a missing file, a missing key column, or an
empty result stops the run. A silently empty axis makes every gene look like it
has one fewer line of evidence, and that error survives all the way to a
conclusion.

## Regions

The `region` key of each eQTL axis selects the tissue. RatGTEx covers the regions the
addiction literature actually uses:

`NAcc` nucleus accumbens core · `IL` infralimbic cortex · `PL` prelimbic cortex
· `OFC` orbitofrontal cortex · `LHb` lateral habenula · `pVTA` posterior
ventral tegmental area · `RMTg` rostromedial tegmental nucleus · `BLA`
basolateral amygdala · `IC` agranular insular cortex · `Brain` hemisphere

GTEx is the default for human, using tissue IDs like
`Brain_Nucleus_accumbens_basal_ganglia`. The GTEx v2 API base and tissue
vocabulary are documented, but the eGene endpoint path used here has not been
verified against a live response — if the first run fails, the error prints the
actual response body, and you can either set `url` to a downloaded eGene table
or switch that axis to `source: file`.

## Usage

```yaml
annotation:
  enabled: true
  species: "Rattus norvegicus"
  orthologs: true
  require: ["DEG"]
  axes:
    - name: DEG
      type: deseq2
      alpha: 0.05
    - name: HUB
      type: wgcna
    - name: eQTL_NAcc
      type: eqtl
      source: ratgtex
      region: "NAcc"
      species: "Rattus norvegicus"
    - name: eQTL_PL
      type: eqtl
      source: ratgtex
      region: "PL"
      species: "Rattus norvegicus"
    - name: GWAS
      type: gwas
      species: "Homo sapiens"
      traits: ["MONDO_0005186"]
  wgcna:
    min_samples: 15
```

```bash
snakemake --use-conda --cores 8
```

Fetched resources are cached under `resources/annotation/` with a provenance
JSON recording the source URL, retrieval date, thresholds, and row counts.
Commit those alongside a result you intend to defend — the resources are
versioned and they move.

## Reference

Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted
correlation network analysis. *BMC Bioinformatics, 9*, 559.
https://doi.org/10.1186/1471-2105-9-559

Munro, D., Wang, T., Chitre, A. S., Polesskaya, O., Ehsan, N., Gao, J.,
Gusev, A., Woods, L. C. S., Saba, L. M., Chen, H., Palmer, A. A., &
Mohammadi, P. (2022). The regulatory landscape of multiple brain regions in
outbred heterogeneous stock rats. *Nucleic Acids Research, 50*(19),
10882–10895. https://doi.org/10.1093/nar/gkac912

Sollis, E., Mosaku, A., Abid, A., Buniello, A., Cerezo, M., Gil, L.,
Groza, T., Güneş, O., Hall, P., Hayhurst, J., Ibrahim, A., Ji, Y., John, S.,
Lewis, E., MacArthur, J. A. L., McMahon, A., Osumi-Sutherland, D., Panoutsopoulou, K.,
Pendlington, Z., ... Harris, L. W. (2023). The NHGRI-EBI GWAS Catalog:
knowledgebase and deposition resource. *Nucleic Acids Research, 51*(D1),
D977–D985. https://doi.org/10.1093/nar/gkac1010
