#!/usr/bin/env python3
"""Fetch GWAS Catalog associations for one or more EFO traits, as a gene list.

Returns the genes the Catalog has mapped to each associated variant. Two things
this deliberately is not:

- It is not a fine-mapped causal gene list. "Mapped gene" in the Catalog is
  largely proximity-based, so a gene appearing here means a trait-associated
  variant sits in or near it, not that the gene mediates the association.
- It is not colocalisation. Overlap with an eQTL gene in the convergence table
  means both signals name the same gene, not that they share a causal variant.

Both caveats matter downstream: convergence is a prioritisation heuristic, and
overstating it is the fastest way to lose an argument about the result.

API: https://www.ebi.ac.uk/gwas/rest/api  (NHGRI-EBI GWAS Catalog)
"""

from __future__ import annotations

import argparse
import datetime
import json
import pathlib
import sys
import time

BASE = "https://www.ebi.ac.uk/gwas/rest/api"


def get_json(url: str, params: dict | None = None, tries: int = 3):
    import requests

    for attempt in range(tries):
        try:
            r = requests.get(
                url, params=params, headers={"Content-Type": "application/json"}, timeout=120
            )
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r.json()
        except Exception:
            if attempt == tries - 1:
                raise
            time.sleep(2 * (attempt + 1))
    return None


def resolve_trait(term: str) -> tuple[str, str] | None:
    """Accept either an EFO id (EFO_0004329) or a free-text trait label."""
    if term.upper().startswith(("EFO_", "MONDO_", "ORPHANET_", "HP_")):
        data = get_json(f"{BASE}/efoTraits/{term}")
        if data:
            return term, data.get("trait", term)
        return None
    data = get_json(f"{BASE}/efoTraits/search/findByEfoTrait", {"trait": term})
    if not data:
        return None
    traits = data.get("_embedded", {}).get("efoTraits", [])
    if not traits:
        return None
    t = traits[0]
    return t.get("shortForm", term), t.get("trait", term)


def associations_for_trait(efo_id: str, pmax: float, page_size: int = 500):
    """Walk the paginated associations endpoint for one EFO id."""
    rows = []
    url = f"{BASE}/efoTraits/{efo_id}/associations"
    params = {"projection": "associationByEfoTrait", "size": page_size}
    while url:
        data = get_json(url, params)
        params = None  # subsequent pages carry their own query string
        if not data:
            break
        assocs = data.get("_embedded", {}).get("associations", [])
        for a in assocs:
            p = a.get("pvalue")
            if p is None or float(p) > pmax:
                continue
            genes = set()
            for locus in a.get("loci", []) or []:
                for ag in locus.get("authorReportedGenes", []) or []:
                    name = (ag.get("geneName") or "").strip()
                    if name and name not in {"NR", "intergenic", "-"}:
                        genes.add(name)
            for snp in a.get("snps", []) or []:
                for ctx in snp.get("genomicContexts", []) or []:
                    g = (ctx.get("gene") or {}).get("geneName")
                    if g:
                        genes.add(g.strip())
            for g in genes:
                rows.append(
                    {
                        "gene_symbol": g,
                        "gwas_pvalue": float(p),
                        "gwas_efo_id": efo_id,
                    }
                )
        url = (data.get("_links", {}).get("next", {}) or {}).get("href")
    return rows


def main() -> None:
    import pandas as pd

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--axis", default="GWAS", help="axis name used to prefix columns")
    ap.add_argument("--traits", nargs="+", required=True,
                    help="EFO ids or trait labels, e.g. EFO_0004329 'cocaine dependence'")
    ap.add_argument("--pvalue", type=float, default=5e-8)
    ap.add_argument("--out", required=True)
    ap.add_argument("--provenance", required=True)
    args = ap.parse_args()

    all_rows, resolved, missing = [], [], []
    for term in args.traits:
        hit = resolve_trait(term)
        if not hit:
            missing.append(term)
            print(f"WARNING: no EFO trait matched '{term}'", file=sys.stderr)
            continue
        efo_id, label = hit
        rows = associations_for_trait(efo_id, args.pvalue)
        print(f"{label} ({efo_id}): {len(rows)} gene-association rows")
        resolved.append({"query": term, "efo_id": efo_id, "label": label,
                         "n_rows": len(rows)})
        all_rows.extend(r | {"gwas_trait": label} for r in rows)

    if not all_rows:
        sys.exit(
            "No associations retrieved. Check the trait terms against\n"
            "https://www.ebi.ac.uk/gwas/efotraits — the Catalog is strict about labels."
        )

    df = pd.DataFrame(all_rows)
    # One row per gene: strongest association, with traits collapsed.
    traits = (
        df.groupby("gene_symbol")["gwas_trait"]
        .apply(lambda s: ";".join(sorted(set(s))))
        .rename("gwas_trait")
    )
    best = df.sort_values("gwas_pvalue").drop_duplicates("gene_symbol", keep="first")
    out_df = best.drop(columns=["gwas_trait"]).merge(traits, on="gene_symbol")
    out_df = out_df[["gene_symbol", "gwas_pvalue", "gwas_trait", "gwas_efo_id"]].rename(
        columns={
            "gwas_pvalue": f"{args.axis}_pvalue",
            "gwas_trait": f"{args.axis}_trait",
            "gwas_efo_id": f"{args.axis}_efo_id",
        }
    )

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out, sep="\t", index=False)

    prov = {
        "name": args.axis,
        "type": "gwas",
        "species": "Homo sapiens",
        "key": "gene_symbol",
        "skipped": False,
        "retrieved_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "api": BASE,
        "pvalue_threshold": args.pvalue,
        "traits_requested": args.traits,
        "traits_resolved": resolved,
        "traits_unmatched": missing,
        "n_genes": int(len(out_df)),
        "note": "mapped/reported genes, largely proximity-based; not fine-mapped",
    }
    pathlib.Path(args.provenance).write_text(json.dumps(prov, indent=2))
    print(f"{len(out_df)} unique GWAS genes -> {out}")


if __name__ == "__main__":
    main()
