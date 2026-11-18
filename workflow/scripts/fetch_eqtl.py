#!/usr/bin/env python3
"""Fetch region-matched cis-eQTL genes and write a normalised gene list.

The point of this step is narrow: for one brain region, which genes have a
credible cis-eQTL? That list becomes one axis of the convergence table. It is
not a colocalisation analysis and makes no claim that the eQTL and any GWAS
signal share a causal variant — see docs/convergence.md.

Sources
-------
ratgtex : RatGTEx (ratgtex.org), outbred HS rats, brain regions NAcc, IL, PL,
          OFC, LHb, pVTA, RMTg, BLA plus non-brain tissues. Region-matched to
          the addiction literature, which is why it is the default here.
gtex    : GTEx v8+ human tissues, via a bulk signif-pairs file.
file    : a local TSV you supply. Use this for any resource not covered above,
          or to pin an exact version for a published analysis.

Every run writes a provenance JSON next to the output recording the source URL,
retrieval date, and row counts, so a convergence table can be traced back to
the exact resource version that produced it.
"""

from __future__ import annotations

import argparse
import datetime
import gzip
import io
import json
import pathlib
import sys

# RatGTEx v3 (rn7). The "conditionally independent cis-eQTLs" table covers all
# tissues in one file with a `tissue` column, which is cheaper than pulling
# per-tissue files. Verified against https://ratgtex.org/download/v3/
RATGTEX_INDEP = "https://ratgtex.org/data/v3/eqtl/eqtls_indep.v3_rn7.txt"

# GTEx Portal API v2 (https://gtexportal.org/api/v2/redoc). The base URL and
# tissue ID vocabulary are documented; the eGene path below is the expected one
# but has NOT been verified against a live response, so it is overridable and
# the code fails loudly with the actual payload if the schema is not what it
# expects. Set annotation eqtl.url to a bulk eGene file to bypass the API.
GTEX_API = "https://gtexportal.org/api/v2"
GTEX_EGENE_PATH = "/association/egene"
RATGTEX_TOP = "https://ratgtex.org/data/v3/eqtl/top_assoc.v3_rn7.txt"
RATGTEX_TISSUES = {
    "NAcc": "nucleus accumbens core",
    "IL": "infralimbic cortex",
    "PL": "prelimbic cortex",
    "OFC": "orbitofrontal cortex",
    "LHb": "lateral habenula",
    "pVTA": "posterior ventral tegmental area",
    "RMTg": "rostromedial tegmental nucleus",
    "BLA": "basolateral amygdala",
    "IC": "agranular insular cortex",
    "Brain": "brain hemisphere",
    "Adipose": "adipose",
    "Eye": "eye",
    "Liver": "liver",
}


def http_get(url: str) -> bytes:
    try:
        import requests
    except ImportError:
        sys.exit("requests is not installed.\n    pip install requests")
    r = requests.get(url, timeout=300)
    r.raise_for_status()
    return r.content


def read_table(raw: bytes, url: str):
    import pandas as pd

    buf = io.BytesIO(raw)
    if url.endswith(".gz"):
        buf = io.BytesIO(gzip.decompress(raw))
    return pd.read_csv(buf, sep="\t", low_memory=False)


def pick(df, *candidates):
    """Return the first column present, or None. Column names drift between
    resource versions; failing loudly beats silently using the wrong one."""
    for c in candidates:
        if c in df.columns:
            return c
    return None


def from_ratgtex(region: str, qval: float, url: str | None):
    if region not in RATGTEX_TISSUES:
        sys.exit(
            f"'{region}' is not a RatGTEx tissue. Available:\n  "
            + "\n  ".join(f"{k:8s} {v}" for k, v in RATGTEX_TISSUES.items())
        )
    src = url or RATGTEX_INDEP
    df = read_table(http_get(src), src)

    tissue_col = pick(df, "tissue", "Tissue")
    gene_col = pick(df, "gene_id", "gene", "phenotype_id")
    q_col = pick(df, "qval", "qvalue", "pval_beta")
    var_col = pick(df, "variant_id", "variant", "top_variant_id")

    if tissue_col is None or gene_col is None:
        sys.exit(
            f"Unexpected RatGTEx schema. Columns found: {list(df.columns)}\n"
            "The file layout has probably changed; set annotation.eqtl.url or "
            "switch to source: file."
        )

    sub = df[df[tissue_col] == region].copy()
    if sub.empty:
        sys.exit(f"No rows for tissue '{region}' in {src}")

    if q_col is not None:
        sub = sub[sub[q_col].astype(float) <= qval]

    out = pd.DataFrame(
        {
            "gene_id": sub[gene_col].astype(str).str.replace(r"\.\d+$", "", regex=True),
            "eqtl_qval": sub[q_col].astype(float) if q_col else float("nan"),
            "eqtl_variant": sub[var_col].astype(str) if var_col else "",
        }
    )
    # One row per gene: keep the strongest association.
    out = out.sort_values("eqtl_qval").drop_duplicates("gene_id", keep="first")
    return out, src, "Rattus norvegicus"


def from_gtex_api(region: str, qval: float):
    """Page the GTEx v2 eGene endpoint for one tissue."""
    import pandas as pd
    import requests

    rows, page, url = [], 0, GTEX_API + GTEX_EGENE_PATH
    while True:
        r = requests.get(
            url,
            params={"tissueSiteDetailId": region, "page": page, "itemsPerPage": 10000},
            headers={"Accept": "application/json"},
            timeout=300,
        )
        if r.status_code != 200:
            sys.exit(
                f"GTEx API returned {r.status_code} for {url}\n"
                f"tissue: {region}\n"
                f"body: {r.text[:400]}\n\n"
                "The eGene endpoint path is unverified. Either correct it via\n"
                "annotation eqtl.url, or set source: file with a downloaded\n"
                "GTEx eGene table."
            )
        payload = r.json()
        chunk = payload.get("data", payload if isinstance(payload, list) else [])
        if not chunk:
            break
        rows.extend(chunk)
        paging = payload.get("paging_info") or {}
        if page + 1 >= int(paging.get("numberOfPages", page + 1)):
            break
        page += 1

    if not rows:
        sys.exit(f"GTEx returned no eGenes for tissue '{region}'.")

    df = pd.DataFrame(rows)
    gene_col = pick(df, "gencodeId", "geneId", "gene_id")
    q_col = pick(df, "qValue", "qval", "pValue")
    if gene_col is None:
        sys.exit(
            f"Unexpected GTEx schema. Keys found: {list(df.columns)[:20]}\n"
            "Set source: file with a downloaded eGene table instead."
        )
    if q_col:
        df = df[df[q_col].astype(float) <= qval]
    out = pd.DataFrame({
        "gene_id": df[gene_col].astype(str).str.replace(r"\.\d+$", "", regex=True),
        "eqtl_qval": df[q_col].astype(float) if q_col else float("nan"),
        "eqtl_variant": df[pick(df, "variantId", "variant_id")].astype(str)
        if pick(df, "variantId", "variant_id") else "",
    })
    out = out.sort_values("eqtl_qval").drop_duplicates("gene_id", keep="first")
    return out, GTEX_API + GTEX_EGENE_PATH, "Homo sapiens"


def from_gtex(region: str, qval: float, url: str | None):
    if not url:
        return from_gtex_api(region, qval)
    df = read_table(http_get(url), url)
    gene_col = pick(df, "gene_id", "phenotype_id", "gene")
    q_col = pick(df, "qval", "qvalue", "pval_beta", "pval_nominal")
    var_col = pick(df, "variant_id", "variant")
    if gene_col is None:
        sys.exit(f"No gene column in {url}. Columns: {list(df.columns)}")
    if q_col is not None:
        df = df[df[q_col].astype(float) <= qval]
    out = pd.DataFrame(
        {
            "gene_id": df[gene_col].astype(str).str.replace(r"\.\d+$", "", regex=True),
            "eqtl_qval": df[q_col].astype(float) if q_col else float("nan"),
            "eqtl_variant": df[var_col].astype(str) if var_col else "",
        }
    )
    out = out.sort_values("eqtl_qval").drop_duplicates("gene_id", keep="first")
    return out, url, "Homo sapiens"


def from_file(path: str, qval: float):
    import pandas as pd

    df = pd.read_csv(path, sep="\t")
    if "gene_id" not in df.columns:
        sys.exit(f"{path} must have a gene_id column. Found: {list(df.columns)}")
    if "eqtl_qval" in df.columns:
        df = df[df["eqtl_qval"].astype(float) <= qval]
    for c in ("eqtl_qval", "eqtl_variant"):
        if c not in df.columns:
            df[c] = "" if c.endswith("variant") else float("nan")
    return df[["gene_id", "eqtl_qval", "eqtl_variant"]], path, "unknown"


def main() -> None:
    import pandas as pd  # noqa: F401  (imported for side effects in helpers)

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--axis", default="eQTL", help="axis name used to prefix columns")
    ap.add_argument("--source", default="gtex", choices=["ratgtex", "gtex", "file"])
    ap.add_argument("--region", required=True, help="tissue/region code, e.g. NAcc")
    ap.add_argument("--qval", type=float, default=0.05)
    ap.add_argument("--url", default=None, help="override or local path for source: file")
    ap.add_argument("--out", required=True)
    ap.add_argument("--provenance", required=True)
    args = ap.parse_args()

    if args.source == "ratgtex":
        table, src, species = from_ratgtex(args.region, args.qval, args.url)
    elif args.source == "gtex":
        table, src, species = from_gtex(args.region, args.qval, args.url)
    else:
        table, src, species = from_file(args.url, args.qval)

    table = table.rename(columns={
        "eqtl_qval": f"{args.axis}_qval",
        "eqtl_variant": f"{args.axis}_variant",
    })

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out, sep="\t", index=False)

    prov = {
        "name": args.axis,
        "type": "eqtl",
        "species": species,
        "key": "gene_id",
        "skipped": False,
        "retrieved_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "source": args.source,
        "url": src,
        "region": args.region,
        "qval_threshold": args.qval,
        "n_genes": int(len(table)),
    }
    pathlib.Path(args.provenance).write_text(json.dumps(prov, indent=2))

    print(f"{len(table)} eQTL genes for region {args.region} -> {out}")
    print(f"resource species: {species}")


if __name__ == "__main__":
    main()
