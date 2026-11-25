#!/usr/bin/env python3
"""Write a synthetic working directory that starts at salmon output.

Genes 0-79 form a module driven by treatment, 80-159 a module unrelated to any
factor, 160-239 a module driven by sex. Genes 240-269 change with treatment
one by one. Gene 270 has an effect large enough to give a p-value of zero.
Genes 271 onward are background; the last 60 have low counts.

The fixture also holds an ortholog table to human for genes 0-999 (genes 0-9
have two orthologs each, genes 1000 onward have none), a human gene list
KNOWN that names the orthologs of genes 240-269, and a mouse gene list TENTH
that names every tenth gene.
"""

import argparse
import pathlib

import numpy as np
import pandas as pd
import yaml

N_GENES = 1500


def simulate(n_per_group: int, seed: int):
    rng = np.random.default_rng(seed)
    n = 2 * n_per_group
    treated = np.repeat([0, 1], n_per_group)
    sex = np.tile([0, 1], n_per_group)

    log_mu = rng.normal(5.5, 1.0, N_GENES)
    log_mu[-60:] = rng.normal(1.2, 0.4, 60)
    eta = np.tile(log_mu[:, None], (1, n))

    z_a = 1.2 * treated + rng.normal(0, 0.5, n)
    z_b = rng.normal(0, 1.0, n)
    z_c = 1.2 * sex + rng.normal(0, 0.5, n)
    eta[0:80] += rng.uniform(0.5, 0.9, 80)[:, None] * z_a
    eta[80:160] += rng.uniform(0.5, 0.9, 80)[:, None] * z_b
    eta[160:240] += rng.uniform(0.5, 0.9, 80)[:, None] * z_c

    lfc = rng.choice([-1.0, 1.0], 30) * 1.0
    eta[240:270] += lfc[:, None] * treated
    eta[270] = np.log(20000.0) + 5.5 * treated - 2.75

    disp = np.full(N_GENES, 0.05)
    disp[270] = 0.01
    mu = np.exp(eta)
    counts = rng.negative_binomial(1.0 / disp[:, None], 1.0 / (1.0 + mu * disp[:, None]))
    if n_per_group < 7:
        counts[300, 0] = counts[300].max() * 400 + 5000
    return counts, treated, sex


def write(workdir: pathlib.Path, n_per_group: int, seed: int, extra_config: dict):
    counts, treated, sex = simulate(n_per_group, seed)
    n = counts.shape[1]
    runs = [f"SRR90000{i:02d}" for i in range(n)]
    samples = pd.DataFrame(
        {
            "run": runs,
            "sample": [f"{'trt' if t else 'ctl'}_{i:02d}" for i, t in enumerate(treated)],
            "condition": ["treated" if t else "control" for t in treated],
            "sex": ["F" if s else "M" for s in sex],
            "library_layout": "SINGLE",
        }
    )
    (workdir / "config").mkdir(parents=True, exist_ok=True)
    samples.to_csv(workdir / "config" / "samples.tsv", sep="\t", index=False)

    genes = [f"ENSMUSG{g:011d}" for g in range(N_GENES)]
    tx = [f"ENSMUST{g:011d}" for g in range(N_GENES)]
    ref = workdir / "resources" / "reference"
    ref.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {"tx_id": tx, "gene_id": genes, "gene_symbol": [f"Gene{g}" for g in range(N_GENES)]}
    ).to_csv(ref / "tx2gene.tsv", sep="\t", index=False)

    eff = 1300.0
    for j, run in enumerate(runs):
        reads = counts[:, j].astype(float)
        rate = reads / eff
        tpm = rate / rate.sum() * 1e6
        out = workdir / "results" / "salmon" / run
        out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            {
                "Name": [f"{t}.1" for t in tx],
                "Length": 1500,
                "EffectiveLength": eff,
                "TPM": tpm,
                "NumReads": reads,
            }
        ).to_csv(out / "quant.sf", sep="\t", index=False)

    config = {
        "samples": "config/samples.tsv",
        "sra": {"max_size": "1g"},
        "reference": {"transcriptome_url": "unused", "genome_url": "unused"},
        "fastp": {"extra": ""},
        "salmon": {"kmer": 31, "libtype": "A", "extra": ""},
        "deseq2": {
            "design": "~ condition",
            "min_count": 10,
            "min_samples": min(3, n_per_group),
            "alpha": 0.05,
            "lfc_shrink": "none",
            "pca_group": "condition",
        },
        "enrichment": {"species": "Mus musculus", "min_size": 15, "max_size": 500},
        "contrasts": {"treated_vs_control": ["condition", "treated", "control"]},
        "annotation": {
            "enabled": True,
            "species": "Mus musculus",
            "orthologs": False,
            "require": [],
            "axes": [
                {"name": "DEG", "type": "deseq2", "alpha": 0.05},
                {"name": "HUB", "type": "wgcna"},
            ],
            "wgcna": {
                "min_samples": 15,
                "n_top_genes": 5000,
                "min_module_size": 20,
                "kme_threshold": 0.8,
                "power": "auto",
            },
        },
    }
    annot = workdir / "resources" / "annotation"
    annot.mkdir(parents=True, exist_ok=True)
    mapped = list(range(1000))
    pairs = [(genes[g], f"ENSG{g:011d}", f"HGENE{g}") for g in mapped]
    pairs += [(genes[g], f"ENSG9{g:010d}", f"HGENE{g}B") for g in range(10)]
    pd.DataFrame(pairs, columns=["gene_id", "ortholog_gene_id", "ortholog_symbol"]).to_csv(
        annot / "orthologs.mmusculus_to_hsapiens.tsv", sep="\t", index=False)
    (annot / "orthologs.mmusculus_to_hsapiens.provenance.json").write_text(
        '{"service": "fixture", "gprofiler_version": "fixture"}')

    lists = workdir / "resources" / "fixture"
    lists.mkdir(parents=True, exist_ok=True)
    known = [f"ENSG{g:011d}" for g in range(240, 270)] + ["ENSG90000000005"]
    pd.DataFrame({"gene_id": known}).to_csv(lists / "known_human.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_id": [genes[g] for g in range(0, N_GENES, 10)]}).to_csv(
        lists / "every_tenth_mouse.tsv", sep="\t", index=False)
    config["annotation"]["orthologs"] = True
    config["annotation"]["axes"] += [
        {"name": "KNOWN", "type": "file", "path": "resources/fixture/known_human.tsv",
         "species": "Homo sapiens"},
        {"name": "TENTH", "type": "file", "path": "resources/fixture/every_tenth_mouse.tsv",
         "species": "Mus musculus"},
    ]

    config["deseq2"]["design"] = extra_config.get("design", "~ condition")
    if extra_config.get("sex_contrast"):
        config["contrasts"]["F_vs_M"] = ["sex", "F", "M"]
    with open(workdir / "config" / "config.yaml", "w") as fh:
        yaml.safe_dump(config, fh, sort_keys=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, help="directory to create the fixture in")
    ap.add_argument("--n-per-group", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--design", default="~ condition")
    ap.add_argument("--sex-contrast", action="store_true",
                    help="add a second contrast, F against M")
    args = ap.parse_args()
    write(pathlib.Path(args.out), args.n_per_group, args.seed,
          {"design": args.design, "sex_contrast": args.sex_contrast})
    print(f"fixture written to {args.out}")
