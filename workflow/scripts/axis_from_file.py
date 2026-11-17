#!/usr/bin/env python3
"""Validate a user-supplied gene list and emit a normalised evidence axis.

This is the extension point. Anything you can express as "here is a set of
genes, and optionally a number per gene" becomes an axis without touching the
workflow code: differential ATAC peaks, CUT&RUN targets, a druggable-genome
list, cell-type markers, a hit list from a previous study.

Input: a TSV with a `gene_id` and/or `gene_symbol` column. Any other columns
are carried through as axis metadata, prefixed with the axis name so several
axes can be joined without collisions.

The validation is deliberately strict. A silently empty or mis-keyed axis makes
every gene look like it has one fewer line of evidence, which is the kind of
error that survives all the way to a conclusion.
"""

from __future__ import annotations

import argparse
import datetime
import json
import pathlib
import sys

RESERVED = {"gene_id", "gene_symbol"}


def main() -> None:
    import pandas as pd

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--axis", required=True, help="axis name, e.g. ATAC")
    ap.add_argument("--path", required=True, help="input TSV")
    ap.add_argument("--species", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--provenance", required=True)
    args = ap.parse_args()

    src = pathlib.Path(args.path)
    if not src.exists():
        sys.exit(
            f"Axis '{args.axis}' points at {src}, which does not exist.\n"
            "Supply the file or remove the axis from annotation.axes."
        )

    sep = "," if src.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(src, sep=sep)

    keys = [c for c in ("gene_id", "gene_symbol") if c in df.columns]
    if not keys:
        sys.exit(
            f"{src} needs a gene_id and/or gene_symbol column.\n"
            f"Found: {list(df.columns)}"
        )
    key = "gene_id" if "gene_id" in keys else "gene_symbol"

    before = len(df)
    df = df.dropna(subset=[key])
    df = df[df[key].astype(str).str.strip() != ""]
    if "gene_id" in df.columns:
        # Strip Ensembl version suffixes so IDs join across annotation releases.
        df["gene_id"] = df["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    df = df.drop_duplicates(subset=[key], keep="first")
    dropped = before - len(df)

    if df.empty:
        sys.exit(f"{src} has no usable rows after cleaning ({before} rows read).")

    rename = {c: f"{args.axis}_{c}" for c in df.columns if c not in RESERVED}
    df = df.rename(columns=rename)

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)

    prov = {
        "name": args.axis,
        "type": "file",
        "species": args.species,
        "key": key,
        "source_path": str(src),
        "n_genes": int(len(df)),
        "n_dropped": int(dropped),
        "metadata_columns": list(rename.values()),
        "skipped": False,
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    pathlib.Path(args.provenance).write_text(json.dumps(prov, indent=2))

    print(f"axis {args.axis}: {len(df)} genes keyed on {key}"
          + (f" ({dropped} rows dropped)" if dropped else ""))


if __name__ == "__main__":
    main()
