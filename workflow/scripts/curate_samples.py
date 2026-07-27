#!/usr/bin/env python3
"""Propose `condition` assignments for an SRA study from free-text sample attributes.

WHY THIS IS NOT A SNAKEMAKE RULE
--------------------------------
Language model output is not deterministic. Putting a model call inside the DAG
would mean the same accession and the same config could produce different counts
on different days, which destroys the reproducibility guarantee the workflow
exists to provide.

So this script sits outside the workflow. It writes a *proposal* plus a
provenance record. You review the proposal, accept it explicitly, and commit the
resulting `samples.tsv` to version control. The pipeline only ever reads that
committed file. Every downstream number traces back to a reviewed artifact, and
`git log` shows who accepted what and when.

WHAT PROBLEM IT SOLVES
----------------------
SRA sample attributes are submitter free text. One study puts the design in
`treatment_sam` as "5 mg/kg cocaine, 24h withdrawal"; the next buries it in
`source_name_sam`; the next uses `genotype_sam` for what is really a treatment.
The COALESCE block in resources/sra_metadata.sql handles clean key/value pairs
and gives up on prose. Mapping prose onto a factor is a language problem.

USAGE
-----
    export ANTHROPIC_API_KEY=sk-ant-...

    # 1. dump the raw attributes for the study
    python workflow/scripts/fetch_sra_metadata.py \
        --study SRP123456 --dump-attributes resources/SRP123456.attributes.json

    # 2. propose assignments (writes proposal + provenance, changes nothing else)
    python workflow/scripts/curate_samples.py \
        --attributes resources/SRP123456.attributes.json

    # 3. read the proposal, then accept it
    python workflow/scripts/curate_samples.py \
        --attributes resources/SRP123456.attributes.json --accept

Requires `pip install anthropic`.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import pathlib
import sys

DEFAULT_MODEL = "claude-sonnet-5"

SYSTEM_PROMPT = """\
You are helping curate a sample sheet for an RNA-seq differential expression \
analysis. You will be given the free-text sample attributes for each run in a \
public SRA study, exactly as the submitter entered them.

Your job is to propose an experimental `condition` label for each run, and to \
flag anything that could confound the analysis.

Rules:
- Base every assignment ONLY on text present in that run's attributes. Do not \
infer from the study accession, from prior knowledge of the publication, or \
from what a study like this "usually" looks like.
- For each assignment, quote the exact attribute value you based it on. If you \
cannot quote a source, the assignment is not supported and must be UNRESOLVED.
- Use "UNRESOLVED" for condition when the attributes do not clearly state the \
experimental group. An honest UNRESOLVED is far more useful than a plausible \
guess, because a wrong group label silently corrupts every result downstream.
- Condition labels must be lowercase snake_case, and must be identical across \
runs that belong to the same group.
- Confidence is "high" only when the attribute text states the group \
explicitly. Use "medium" when it requires interpretation, "low" when you are \
reading between the lines.
- In candidate_covariates, list attributes that vary across runs and could \
confound the comparison (batch, sex, age, collection date, sequencing run, \
tissue, cell line). Include the attribute key and why it matters.
- In concerns, note anything a human should look at: unbalanced groups, \
apparent technical replicates, runs whose attributes contradict each other, \
groups with fewer than three samples.

Return ONLY a JSON object with this shape, no prose and no markdown fences:

{
  "assignments": [
    {"run": "SRR...", "condition": "...", "evidence": "...",
     "evidence_key": "...", "confidence": "high|medium|low"}
  ],
  "candidate_covariates": [{"key": "...", "reason": "..."}],
  "concerns": ["..."]
}
"""


def load_attributes(path: pathlib.Path) -> dict:
    data = json.loads(path.read_text())
    if not isinstance(data, dict) or not data:
        sys.exit(f"{path} is not a non-empty JSON object of run -> attributes.")
    return data


def call_model(attributes: dict, model: str) -> tuple[dict, dict]:
    """Return (parsed_response, provenance)."""
    try:
        import anthropic
    except ImportError:
        sys.exit("anthropic is not installed.\n    pip install anthropic")

    if not os.environ.get("ANTHROPIC_API_KEY"):
        sys.exit(
            "ANTHROPIC_API_KEY is not set.\n"
            "Create a key at https://platform.claude.com and export it:\n"
            "    export ANTHROPIC_API_KEY=sk-ant-..."
        )

    user_content = json.dumps(attributes, indent=2, sort_keys=True)

    client = anthropic.Anthropic()
    resp = client.messages.create(
        model=model,
        max_tokens=8000,
        temperature=0,  # not a determinism guarantee, but removes gratuitous variance
        system=SYSTEM_PROMPT,
        messages=[{"role": "user", "content": user_content}],
    )

    text = "".join(block.text for block in resp.content if block.type == "text").strip()
    text = text.removeprefix("```json").removeprefix("```").removesuffix("```").strip()

    try:
        parsed = json.loads(text)
    except json.JSONDecodeError as exc:
        sys.exit(f"Model did not return parseable JSON: {exc}\n\n{text[:2000]}")

    provenance = {
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "model": model,
        "temperature": 0,
        "system_prompt_sha256": hashlib.sha256(SYSTEM_PROMPT.encode()).hexdigest(),
        "input_sha256": hashlib.sha256(user_content.encode()).hexdigest(),
        "n_runs_submitted": len(attributes),
        "usage": {
            "input_tokens": resp.usage.input_tokens,
            "output_tokens": resp.usage.output_tokens,
        },
        "raw_response": text,
    }
    return parsed, provenance


def validate(parsed: dict, attributes: dict) -> list[dict]:
    """Structural checks. The model is untrusted input, not an authority."""
    assignments = parsed.get("assignments")
    if not isinstance(assignments, list) or not assignments:
        sys.exit("Response contained no assignments.")

    submitted = set(attributes)
    returned = {a.get("run") for a in assignments}

    if returned - submitted:
        sys.exit(f"Model invented runs not in the input: {sorted(returned - submitted)}")
    if submitted - returned:
        sys.exit(f"Model dropped runs: {sorted(submitted - returned)}")
    if len(returned) != len(assignments):
        sys.exit("Model returned duplicate run entries.")

    for a in assignments:
        if not a.get("condition"):
            sys.exit(f"Missing condition for {a.get('run')}.")
        if a["condition"] != "UNRESOLVED" and not a.get("evidence"):
            sys.exit(f"Assignment for {a['run']} has no supporting evidence quote.")

    return assignments


def report(parsed: dict, assignments: list[dict], attributes: dict) -> int:
    """Print the proposal for human review. Returns count of rows needing attention."""
    from collections import Counter

    print("\nPROPOSED ASSIGNMENTS")
    print("-" * 78)
    width = max(len(a["run"]) for a in assignments)
    needs_attention = 0

    for a in sorted(assignments, key=lambda x: (x["condition"], x["run"])):
        flag = " " if a.get("confidence") == "high" else "!"
        if flag == "!" or a["condition"] == "UNRESOLVED":
            needs_attention += 1
        print(f"{flag} {a['run']:<{width}}  {a['condition']:<24} [{a.get('confidence','?')}]")
        ev = (a.get("evidence") or "").replace("\n", " ")
        if ev:
            key = a.get("evidence_key", "?")
            print(f"  {'':<{width}}  from {key}: {ev[:90]}")

    print("\nGROUP SIZES")
    for cond, n in sorted(Counter(a["condition"] for a in assignments).items()):
        warn = "  <- fewer than 3 replicates" if n < 3 else ""
        print(f"  {cond:<24} {n}{warn}")

    covs = parsed.get("candidate_covariates") or []
    if covs:
        print("\nCANDIDATE COVARIATES (consider adding to deseq2.design)")
        for c in covs:
            print(f"  {c.get('key')}: {c.get('reason')}")

    concerns = parsed.get("concerns") or []
    if concerns:
        print("\nCONCERNS")
        for c in concerns:
            print(f"  - {c}")

    return needs_attention


def write_outputs(
    assignments: list[dict],
    attributes: dict,
    provenance: dict,
    samples_path: pathlib.Path,
    proposal_path: pathlib.Path,
    provenance_path: pathlib.Path,
    accept: bool,
) -> None:
    import pandas as pd

    proposal = pd.DataFrame(assignments)[
        ["run", "condition", "confidence", "evidence_key", "evidence"]
    ]
    proposal.to_csv(proposal_path, sep="\t", index=False)
    provenance_path.write_text(json.dumps(provenance, indent=2))

    print(f"\nProposal written to {proposal_path}")
    print(f"Provenance written to {provenance_path}")

    if not accept:
        print(
            "\nNothing else was modified. Review the proposal above, then re-run "
            "with --accept to write these conditions into the sample sheet."
        )
        return

    if not samples_path.exists():
        sys.exit(f"{samples_path} does not exist. Build it first with fetch_sra_metadata.py.")

    samples = pd.read_csv(samples_path, sep="\t", dtype=str)
    mapping = dict(zip(proposal["run"], proposal["condition"]))

    missing = set(samples["run"]) - set(mapping)
    if missing:
        sys.exit(f"Sample sheet contains runs absent from the proposal: {sorted(missing)}")

    old = samples["condition"].copy() if "condition" in samples else None
    samples["condition"] = samples["run"].map(mapping)

    if old is not None:
        changed = (old != samples["condition"]).sum()
        print(f"\nUpdated {changed} of {len(samples)} condition values.")

    if (samples["condition"] == "UNRESOLVED").any():
        n = int((samples["condition"] == "UNRESOLVED").sum())
        print(
            f"WARNING: {n} runs are UNRESOLVED. Fill these in by hand before "
            "running the workflow — DESeq2 will happily model a nonsense factor level."
        )

    samples.to_csv(samples_path, sep="\t", index=False)
    print(f"Wrote {samples_path}. Commit it alongside {provenance_path.name}.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--attributes", required=True, help="JSON from --dump-attributes")
    ap.add_argument("--samples", default="config/samples.tsv")
    ap.add_argument("--model", default=DEFAULT_MODEL)
    ap.add_argument(
        "--accept",
        action="store_true",
        help="write the proposed conditions into the sample sheet",
    )
    args = ap.parse_args()

    attr_path = pathlib.Path(args.attributes)
    attributes = load_attributes(attr_path)

    stem = attr_path.name.replace(".attributes.json", "")
    proposal_path = attr_path.with_name(f"{stem}.curation_proposal.tsv")
    provenance_path = attr_path.with_name(f"{stem}.curation_provenance.json")

    print(f"Submitting {len(attributes)} runs to {args.model}...")
    parsed, provenance = call_model(attributes, args.model)
    assignments = validate(parsed, attributes)
    needs_attention = report(parsed, assignments, attributes)

    if needs_attention:
        print(
            f"\n{needs_attention} row(s) marked with ! are unresolved or below high "
            "confidence. Check those against the SRA record before accepting."
        )

    write_outputs(
        assignments,
        attributes,
        provenance,
        pathlib.Path(args.samples),
        proposal_path,
        provenance_path,
        args.accept,
    )


if __name__ == "__main__":
    main()
