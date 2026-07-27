#!/usr/bin/env python3
"""Generate config/samples.tsv from the public SRA metadata table in BigQuery.

The SRA moved its metadata into BigQuery, so a study's run list, library
layout, and sample attributes can be pulled with SQL instead of scraped from
the SRA Run Selector by hand:
https://www.ncbi.nlm.nih.gov/sra/docs/sra-bigquery/

Usage
-----
    # what attribute keys does this study actually use?
    python workflow/scripts/fetch_sra_metadata.py --study SRP123456 --audit

    # build the sample sheet
    python workflow/scripts/fetch_sra_metadata.py \
        --study SRP123456 --out config/samples.tsv

Requires `pip install google-cloud-bigquery pandas` and an authenticated GCP
project (`gcloud auth application-default login`). Queries against the public
dataset are billed to your own project; the metadata table is large, so keep
filters on sra_study rather than scanning it whole.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

QUERY_FILE = pathlib.Path(__file__).resolve().parents[2] / "resources" / "sra_metadata.sql"

AUDIT_SQL = """
SELECT
  LOWER(a.k) AS attribute_key,
  COUNT(DISTINCT m.acc) AS n_runs,
  ARRAY_AGG(DISTINCT a.v LIMIT 10) AS example_values
FROM `nih-sra-datastore.sra.metadata` AS m, UNNEST(m.attributes) AS a
WHERE m.sra_study = @study
GROUP BY attribute_key
ORDER BY n_runs DESC
"""

ATTR_SQL = """
SELECT m.acc AS run, m.jattr AS attributes
FROM `nih-sra-datastore.sra.metadata` AS m
WHERE m.sra_study = @study AND m.consent = 'public'
ORDER BY m.acc
"""

REQUIRED = ["run", "sample", "condition", "library_layout"]


def load_query() -> str:
    """Read the sample-sheet query, dropping the trailing commented audit block."""
    text = QUERY_FILE.read_text()
    return text.split("-- AUDIT QUERY")[0]


def run_query(sql: str, study: str, project: str | None):
    try:
        from google.cloud import bigquery
    except ImportError:
        sys.exit(
            "google-cloud-bigquery is not installed.\n"
            "    pip install google-cloud-bigquery pandas db-dtypes"
        )

    client = bigquery.Client(project=project)
    job_config = bigquery.QueryJobConfig(
        query_parameters=[bigquery.ScalarQueryParameter("study", "STRING", study)]
    )
    return client.query(sql, job_config=job_config).to_dataframe()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--study", required=True, help="SRA study accession, e.g. SRP123456")
    ap.add_argument("--out", default="config/samples.tsv", help="output TSV path")
    ap.add_argument("--project", default=None, help="GCP billing project ID")
    ap.add_argument(
        "--dump-attributes",
        metavar="PATH",
        help="write raw per-run sample attributes as JSON (input to curate_samples.py)",
    )
    ap.add_argument(
        "--audit",
        action="store_true",
        help="list the sample attribute keys used by this study and exit",
    )
    args = ap.parse_args()

    if args.audit:
        df = run_query(AUDIT_SQL, args.study, args.project)
        if df.empty:
            sys.exit(f"No public runs found for {args.study}.")
        print(df.to_string(index=False))
        return

    if args.dump_attributes:
        import json

        df = run_query(ATTR_SQL, args.study, args.project)
        if df.empty:
            sys.exit(f"No public runs found for {args.study}.")
        out = {}
        for _, row in df.iterrows():
            try:
                out[row["run"]] = json.loads(row["attributes"]) if row["attributes"] else {}
            except (TypeError, ValueError):
                out[row["run"]] = {"_unparsed_jattr": str(row["attributes"])}
        p = pathlib.Path(args.dump_attributes)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(out, indent=2, sort_keys=True))
        print(f"Wrote raw attributes for {len(out)} runs to {p}")
        return

    df = run_query(load_query(), args.study, args.project)

    if df.empty:
        sys.exit(
            f"No public RNA-Seq runs found for {args.study}. "
            "Check the accession, or relax the assay_type filter in the SQL."
        )

    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        sys.exit(f"Query result is missing required columns: {missing}")

    bad_layout = sorted(set(df["library_layout"]) - {"PAIRED", "SINGLE"})
    if bad_layout:
        sys.exit(f"Unexpected library_layout values: {bad_layout}")

    if df["run"].duplicated().any():
        sys.exit("Duplicate run accessions returned — inspect the query result.")

    # Disambiguate repeated sample names (multi-run samples are common).
    if df["sample"].duplicated().any():
        df["sample"] = df["sample"] + "_" + df.groupby("sample").cumcount().add(1).astype(str)

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)

    n_unspecified = int((df["condition"] == "unspecified").sum())
    print(f"Wrote {len(df)} runs to {out}")
    print(df["condition"].value_counts().to_string())
    if n_unspecified:
        print(
            f"\nWARNING: {n_unspecified} runs have no parsed condition. "
            "Re-run with --audit and edit the CASE block in resources/sra_metadata.sql, "
            "or fill the column in by hand before running the workflow."
        )


if __name__ == "__main__":
    main()
