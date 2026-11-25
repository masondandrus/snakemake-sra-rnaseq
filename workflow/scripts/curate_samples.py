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

    # 2. propose assignments (one model call; writes proposal + provenance)
    python workflow/scripts/curate_samples.py \
        --attributes resources/SRP123456.attributes.json

    # 3. read resources/SRP123456.curation_proposal.tsv. To overrule a label,
    #    edit its condition and set its confidence to "manual".

    # 4. accept the saved proposal (no model call; writes a new sample sheet)
    python workflow/scripts/curate_samples.py \
        --attributes resources/SRP123456.attributes.json \
        --accept --out config/samples.curated.tsv

Step 4 reads the proposal file from step 2. It never contacts the model, so
the labels that enter the sample sheet are the labels that were reviewed. It
stops if the attributes changed after the proposal, and it records whether
the proposal file was edited by hand.

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


def _strings(obj):
    if isinstance(obj, dict):
        for v in obj.values():
            yield from _strings(v)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            yield from _strings(v)
    elif obj is not None:
        yield str(obj)


def _norm(text: str) -> str:
    return " ".join(str(text).split())


def quote_found(evidence: str, value) -> bool:
    """True if the evidence quote occurs verbatim in this attribute value."""
    quote = _norm(evidence).strip("\"'")
    return bool(quote) and any(quote in _norm(v) for v in _strings(value))


def _resolved(a: dict) -> bool:
    return a["condition"] != "UNRESOLVED" and a.get("confidence") != "manual"


def unsupported_quotes(assignments: list[dict], attributes: dict) -> list[str]:
    """Runs whose quote is not in the attribute that evidence_key names."""
    bad = []
    for a in assignments:
        if not _resolved(a):
            continue
        key = a.get("evidence_key")
        run_attrs = attributes[a["run"]]
        if key not in run_attrs or not quote_found(a["evidence"], run_attrs[key]):
            bad.append(a["run"])
    return sorted(bad)


def shared_quotes(assignments: list[dict], attributes: dict) -> list[str]:
    """Quotes that also occur, under the same key, in a run of another condition.

    Such a quote cannot be the reason for the label, because a run with a
    different label carries the same text. Runs that a reviewer labelled by
    hand are not used for the comparison.
    """
    clashes = []
    for a in assignments:
        if not _resolved(a):
            continue
        key = a.get("evidence_key")
        for b in assignments:
            if b["condition"] in (a["condition"], "UNRESOLVED"):
                continue
            if b.get("confidence") == "manual":
                continue
            other = attributes[b["run"]].get(key)
            if other is not None and quote_found(a["evidence"], other):
                clashes.append(f"{a['run']} ({a['condition']}) quotes text also in "
                               f"{b['run']} ({b['condition']})")
                break
    return clashes


def validate(assignments, attributes: dict) -> list[dict]:
    """Structural checks. The model is untrusted input, not an authority."""
    if not isinstance(assignments, list) or not assignments:
        sys.exit("Response contained no assignments.")

    submitted = set(attributes)
    returned = {a.get("run") for a in assignments}

    if returned - submitted:
        sys.exit(f"Assignments name runs not in the input: {sorted(returned - submitted)}")
    if submitted - returned:
        sys.exit(f"Assignments are missing runs: {sorted(submitted - returned)}")
    if len(returned) != len(assignments):
        sys.exit("Assignments contain duplicate run entries.")

    for a in assignments:
        if not a.get("condition"):
            sys.exit(f"Missing condition for {a.get('run')}.")
        if _resolved(a) and not a.get("evidence"):
            sys.exit(f"Assignment for {a['run']} has no supporting evidence quote.")

    bad = unsupported_quotes(assignments, attributes)
    if bad:
        sys.exit(
            "Evidence quote not found in the attribute named by evidence_key for: "
            f"{bad}"
        )

    clashes = shared_quotes(assignments, attributes)
    if clashes:
        sys.exit(
            "Evidence quote does not distinguish the label from another label:\n  "
            + "\n  ".join(clashes)
        )

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


def sha256_file(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def attributes_hash(attributes: dict) -> str:
    return hashlib.sha256(json.dumps(attributes, indent=2, sort_keys=True).encode()).hexdigest()


PROPOSAL_COLUMNS = ["run", "condition", "confidence", "evidence_key", "evidence"]


def write_proposal(assignments, provenance, proposal_path, provenance_path) -> None:
    import pandas as pd

    pd.DataFrame(assignments).reindex(columns=PROPOSAL_COLUMNS).to_csv(
        proposal_path, sep="\t", index=False
    )
    provenance["proposal_sha256"] = sha256_file(proposal_path)
    provenance_path.write_text(json.dumps(provenance, indent=2))
    print(f"\nProposal written to {proposal_path}")
    print(f"Provenance written to {provenance_path}")
    print(
        "\nNothing else was modified. Review the proposal, then run again with "
        "--accept --out <new sample sheet>."
    )


def accept(attributes, proposal_path, provenance_path, acceptance_path,
           samples_path, out_path, allow_unresolved, force) -> None:
    """Write a sample sheet from the saved proposal. Makes no model call."""
    import pandas as pd

    for path in (proposal_path, provenance_path):
        if not path.exists():
            sys.exit(f"{path} does not exist. Run without --accept first.")
    if out_path is None:
        sys.exit("--accept needs --out, the path of the new sample sheet.")
    if out_path.exists() and not force:
        sys.exit(f"{out_path} exists. Choose another --out or pass --force.")
    if not samples_path.exists():
        sys.exit(f"{samples_path} does not exist. Build it first with fetch_sra_metadata.py.")

    provenance = json.loads(provenance_path.read_text())
    if provenance.get("input_sha256") != attributes_hash(attributes):
        sys.exit(
            "The attributes file changed after the proposal was made. "
            "Run without --accept to make a new proposal."
        )

    proposal = pd.read_csv(proposal_path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in PROPOSAL_COLUMNS if c not in proposal.columns]
    if missing:
        sys.exit(f"{proposal_path} is missing columns: {missing}")
    assignments = validate(proposal.to_dict("records"), attributes)

    n_unresolved = sum(a["condition"] == "UNRESOLVED" for a in assignments)
    if n_unresolved and not allow_unresolved:
        sys.exit(
            f"{n_unresolved} run(s) are UNRESOLVED. Set their condition in "
            f"{proposal_path} with confidence 'manual', or pass --allow-unresolved."
        )

    samples = pd.read_csv(samples_path, sep="\t", dtype=str)
    mapping = {a["run"]: a["condition"] for a in assignments}
    absent = set(samples["run"]) - set(mapping)
    if absent:
        sys.exit(f"Sample sheet contains runs absent from the proposal: {sorted(absent)}")

    old = samples["condition"].copy() if "condition" in samples else None
    samples["condition"] = samples["run"].map(mapping)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    samples.to_csv(out_path, sep="\t", index=False)

    proposal_sha = sha256_file(proposal_path)
    record = {
        "accepted_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "proposal": str(proposal_path),
        "proposal_sha256": proposal_sha,
        "proposal_edited_after_model": proposal_sha != provenance.get("proposal_sha256"),
        "attributes_sha256": provenance.get("input_sha256"),
        "model": provenance.get("model"),
        "n_runs": len(assignments),
        "n_manual": sum(a.get("confidence") == "manual" for a in assignments),
        "n_unresolved": n_unresolved,
        "sample_sheet_in": str(samples_path),
        "sample_sheet_out": str(out_path),
        "sample_sheet_out_sha256": sha256_file(out_path),
    }
    acceptance_path.write_text(json.dumps(record, indent=2))

    if old is not None:
        print(f"Changed {(old != samples['condition']).sum()} of {len(samples)} condition values.")
    if record["proposal_edited_after_model"]:
        print("The proposal file was edited after the model wrote it; the record says so.")
    print(f"Wrote {out_path} and {acceptance_path}. {samples_path} was not modified.")
    print(f"Commit {out_path.name}, {proposal_path.name}, {provenance_path.name} "
          f"and {acceptance_path.name} together.")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--attributes", required=True, help="JSON from --dump-attributes")
    ap.add_argument("--samples", default="config/samples.tsv",
                    help="existing sample sheet; read, never modified")
    ap.add_argument("--model", default=DEFAULT_MODEL)
    ap.add_argument("--accept", action="store_true",
                    help="write a sample sheet from the saved proposal; makes no model call")
    ap.add_argument("--out", default=None, help="with --accept: path of the new sample sheet")
    ap.add_argument("--allow-unresolved", action="store_true",
                    help="with --accept: write the sheet even if runs are UNRESOLVED")
    ap.add_argument("--force", action="store_true",
                    help="with --accept: overwrite --out if it exists")
    args = ap.parse_args(argv)

    attr_path = pathlib.Path(args.attributes)
    attributes = load_attributes(attr_path)

    stem = attr_path.name.replace(".attributes.json", "")
    proposal_path = attr_path.with_name(f"{stem}.curation_proposal.tsv")
    provenance_path = attr_path.with_name(f"{stem}.curation_provenance.json")
    acceptance_path = attr_path.with_name(f"{stem}.curation_acceptance.json")

    if args.accept:
        accept(attributes, proposal_path, provenance_path, acceptance_path,
               pathlib.Path(args.samples), pathlib.Path(args.out) if args.out else None,
               args.allow_unresolved, args.force)
        return

    print(f"Submitting {len(attributes)} runs to {args.model}...")
    parsed, provenance = call_model(attributes, args.model)
    assignments = validate(parsed.get("assignments"), attributes)
    needs_attention = report(parsed, assignments, attributes)

    if needs_attention:
        print(
            f"\n{needs_attention} row(s) marked with ! are unresolved or below high "
            "confidence. Check those against the SRA record before accepting."
        )

    write_proposal(assignments, provenance, proposal_path, provenance_path)


if __name__ == "__main__":
    main()
