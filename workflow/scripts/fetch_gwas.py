#!/usr/bin/env python3
"""Fetch GWAS Catalog associations for one or more EFO traits, as a gene list.

Returns the genes the Catalog has mapped to each associated variant. Two things
this deliberately is not:

- It is not a fine-mapped causal gene list. "Mapped gene" in the Catalog is
  positional (Ensembl genes overlapping or nearest to the variant), so a gene
  appearing here means a trait-associated variant sits in or near it, not that
  the gene mediates the association.
- It is not colocalisation. Overlap with an eQTL gene in the convergence table
  means both signals name the same gene, not that they share a causal variant.

Both caveats matter downstream: convergence is a prioritisation heuristic, and
overstating it is the fastest way to lose an argument about the result.

API: https://www.ebi.ac.uk/gwas/rest/api/v2  (NHGRI-EBI GWAS Catalog, REST v2;
v1 was retired in August 2026)
"""

from __future__ import annotations

import argparse
import datetime
import json
import pathlib
import sys
import time

BASE = "https://www.ebi.ac.uk/gwas/rest/api/v2"
ONTOLOGY_PREFIXES = ("EFO_", "MONDO_", "ORPHANET_", "HP_", "OBA_", "GO_", "CHEBI_")


def get_json(url: str, params: dict | None = None, tries: int = 3):
    import requests

    for attempt in range(tries):
        try:
            r = requests.get(
                url, params=params, headers={"Accept": "application/json"}, timeout=120
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


def is_ontology_id(term: str) -> bool:
    return term.upper().startswith(ONTOLOGY_PREFIXES)


def parse_associations(data: dict, term: str, pmax: float) -> tuple[list[dict], set[str]]:
    """Rows and trait labels from one page of /associations."""
    rows, labels = [], set()
    for a in (data.get("_embedded") or {}).get("associations") or []:
        traits = a.get("efo_traits") or []
        if is_ontology_id(term):
            matched = [t for t in traits if (t.get("efo_id") or "").upper() == term.upper()]
        else:
            matched = [t for t in traits if (t.get("efo_trait") or "").lower() == term.lower()]
        matched = matched or traits
        labels.update(t["efo_trait"] for t in matched if t.get("efo_trait"))

        p = a.get("p_value")
        if p is None or float(p) > pmax:
            continue
        label = ";".join(sorted({t.get("efo_trait", "") for t in matched} - {""})) or term
        efo_id = ";".join(sorted({t.get("efo_id", "") for t in matched} - {""})) or term
        for g in {(g or "").strip() for g in a.get("mapped_genes") or []} - {""}:
            rows.append(
                {
                    "gene_symbol": g,
                    "gwas_pvalue": float(p),
                    "gwas_efo_id": efo_id,
                    "gwas_trait": label,
                }
            )
    return rows, labels


def associations_for_trait(term: str, pmax: float, child_traits: bool, page_size: int):
    """Walk the paginated v2 associations endpoint for one trait."""
    key = "efo_id" if is_ontology_id(term) else "efo_trait"
    params = {
        key: term,
        "show_child_traits": str(child_traits).lower(),
        "size": page_size,
        "page": 0,
    }
    url = f"{BASE}/associations"
    rows, labels, n_seen, pages = [], set(), 0, 0
    while True:
        data = get_json(url, params)
        if not data:
            break
        page_rows, page_labels = parse_associations(data, term, pmax)
        n_page = len((data.get("_embedded") or {}).get("associations") or [])
        n_seen += n_page
        pages += 1
        rows.extend(page_rows)
        labels.update(page_labels)

        nxt = ((data.get("_links") or {}).get("next") or {}).get("href")
        page = data.get("page") or {}
        if nxt:
            url, params = nxt, None
        elif params and page.get("totalPages") and params["page"] + 1 < page["totalPages"]:
            params["page"] += 1
        else:
            break
        if n_page == 0:
            break
    return rows, sorted(labels), n_seen, pages


def main() -> None:
    import pandas as pd

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--axis", default="GWAS", help="axis name used to prefix columns")
    ap.add_argument("--traits", nargs="+", required=True,
                    help="ontology ids or exact trait labels, e.g. EFO_0002610 'cocaine dependence'")
    ap.add_argument("--pvalue", type=float, default=5e-8)
    ap.add_argument("--child-traits", action="store_true",
                    help="also include associations for more specific child traits")
    ap.add_argument("--page-size", type=int, default=100)
    ap.add_argument("--out", required=True)
    ap.add_argument("--provenance", required=True)
    args = ap.parse_args()

    all_rows, resolved, missing = [], [], []
    for term in args.traits:
        rows, labels, n_seen, pages = associations_for_trait(
            term, args.pvalue, args.child_traits, args.page_size
        )
        if n_seen == 0:
            missing.append(term)
            print(f"WARNING: no associations returned for '{term}'", file=sys.stderr)
            continue
        label = "; ".join(labels) or term
        print(f"{label} ({term}): {n_seen} associations over {pages} page(s), "
              f"{len(rows)} gene rows at p <= {args.pvalue}")
        resolved.append({"query": term, "label": label, "n_associations": n_seen,
                         "n_rows": len(rows)})
        all_rows.extend(rows)

    if not all_rows:
        sys.exit(
            "No associations retrieved. Check the trait terms against\n"
            "https://www.ebi.ac.uk/gwas/efotraits — the Catalog is strict about labels."
        )

    df = pd.DataFrame(all_rows)
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
        "show_child_traits": args.child_traits,
        "traits_requested": args.traits,
        "traits_resolved": resolved,
        "traits_unmatched": missing,
        "n_genes": int(len(out_df)),
        "note": "Catalog mapped genes (positional, Ensembl); not fine-mapped",
    }
    pathlib.Path(args.provenance).write_text(json.dumps(prov, indent=2))
    print(f"{len(out_df)} unique GWAS genes -> {out}")


if __name__ == "__main__":
    main()
