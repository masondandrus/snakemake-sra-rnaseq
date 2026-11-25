import json

import pandas as pd
import pytest

ATTRS = {
    "SRR1": {"organism": "Mus musculus", "treatment": "cocaine 20 mg/kg", "tissue": "NAc"},
    "SRR2": {"organism": "Mus musculus", "treatment": "cocaine 20 mg/kg", "tissue": "NAc"},
    "SRR3": {"organism": "Mus musculus", "treatment": "saline", "tissue": "NAc"},
    "SRR4": {"organism": "Mus musculus", "treatment": "saline", "tissue": "NAc"},
}


def good():
    return [
        {"run": r, "condition": c, "evidence": e, "evidence_key": "treatment", "confidence": "high"}
        for r, c, e in [("SRR1", "cocaine", "cocaine 20 mg/kg"), ("SRR2", "cocaine", "cocaine"),
                        ("SRR3", "saline", "saline"), ("SRR4", "saline", "saline")]
    ]


def test_a_supported_proposal_passes(curate):
    assert len(curate.validate(good(), ATTRS)) == 4


def test_a_quote_from_another_attribute_is_rejected(curate):
    a = good()
    a[0]["evidence"] = "Mus musculus"
    with pytest.raises(SystemExit) as err:
        curate.validate(a, ATTRS)
    assert "evidence_key" in str(err.value)


def test_a_quote_shared_with_another_condition_is_rejected(curate):
    a = good()
    for row in a:
        row["evidence"], row["evidence_key"] = "Mus musculus", "organism"
    with pytest.raises(SystemExit) as err:
        curate.validate(a, ATTRS)
    assert "does not distinguish" in str(err.value)


def test_invented_dropped_and_duplicated_runs_are_rejected(curate):
    for mutate in (
        lambda a: a.append({**a[0], "run": "SRR9"}),
        lambda a: a.pop(),
        lambda a: a.__setitem__(3, dict(a[2])),
    ):
        a = good()
        mutate(a)
        with pytest.raises(SystemExit):
            curate.validate(a, ATTRS)


@pytest.fixture
def study(tmp_path, curate, monkeypatch):
    attrs = tmp_path / "SRP1.attributes.json"
    attrs.write_text(json.dumps(ATTRS))
    samples = tmp_path / "samples.tsv"
    pd.DataFrame({"run": list(ATTRS), "sample": ["a", "b", "c", "d"],
                  "condition": "unspecified", "library_layout": "SINGLE"}).to_csv(
        samples, sep="\t", index=False)
    calls = []

    def fake_model(attributes, model):
        calls.append(model)
        return {"assignments": good()}, {
            "model": model, "input_sha256": curate.attributes_hash(attributes)}

    monkeypatch.setattr(curate, "call_model", fake_model)
    return {"dir": tmp_path, "attrs": attrs, "samples": samples, "calls": calls}


def test_accept_reads_the_saved_proposal_and_never_calls_the_model(curate, study):
    curate.main(["--attributes", str(study["attrs"]), "--samples", str(study["samples"])])
    assert len(study["calls"]) == 1
    before = study["samples"].read_text()

    out = study["dir"] / "samples.curated.tsv"
    curate.main(["--attributes", str(study["attrs"]), "--samples", str(study["samples"]),
                 "--accept", "--out", str(out)])
    assert len(study["calls"]) == 1
    assert study["samples"].read_text() == before
    sheet = pd.read_csv(out, sep="\t")
    assert list(sheet["condition"]) == ["cocaine", "cocaine", "saline", "saline"]
    record = json.loads((study["dir"] / "SRP1.curation_acceptance.json").read_text())
    assert record["proposal_edited_after_model"] is False


def test_a_hand_edit_is_accepted_and_recorded(curate, study):
    curate.main(["--attributes", str(study["attrs"]), "--samples", str(study["samples"])])
    proposal = study["dir"] / "SRP1.curation_proposal.tsv"
    table = pd.read_csv(proposal, sep="\t", dtype=str)
    table.loc[table["run"] == "SRR2", ["condition", "confidence"]] = ["cocaine_high", "manual"]
    table.to_csv(proposal, sep="\t", index=False)
    out = study["dir"] / "new.tsv"
    curate.main(["--attributes", str(study["attrs"]), "--samples", str(study["samples"]),
                 "--accept", "--out", str(out)])
    record = json.loads((study["dir"] / "SRP1.curation_acceptance.json").read_text())
    assert record["proposal_edited_after_model"] is True
    assert record["n_manual"] == 1
    assert pd.read_csv(out, sep="\t")["condition"].iloc[1] == "cocaine_high"


def test_accept_stops_when_attributes_changed_or_output_exists_or_runs_unresolved(curate, study):
    args = ["--attributes", str(study["attrs"]), "--samples", str(study["samples"])]
    curate.main(args)
    out = study["dir"] / "new.tsv"

    proposal = study["dir"] / "SRP1.curation_proposal.tsv"
    original = proposal.read_text()
    proposal.write_text(original.replace("SRR4\tsaline", "SRR4\tUNRESOLVED"))
    with pytest.raises(SystemExit) as err:
        curate.main(args + ["--accept", "--out", str(out)])
    assert "UNRESOLVED" in str(err.value)
    proposal.write_text(original)

    out.write_text("occupied")
    with pytest.raises(SystemExit) as err:
        curate.main(args + ["--accept", "--out", str(out)])
    assert "exists" in str(err.value)
    out.unlink()

    changed = dict(ATTRS, SRR4={"organism": "Mus musculus", "treatment": "cocaine"})
    study["attrs"].write_text(json.dumps(changed))
    with pytest.raises(SystemExit) as err:
        curate.main(args + ["--accept", "--out", str(out)])
    assert "changed after the proposal" in str(err.value)


def test_accept_without_a_proposal_or_without_out_stops(curate, study):
    args = ["--attributes", str(study["attrs"]), "--samples", str(study["samples"])]
    with pytest.raises(SystemExit):
        curate.main(args + ["--accept", "--out", str(study["dir"] / "x.tsv")])
    curate.main(args)
    with pytest.raises(SystemExit) as err:
        curate.main(args + ["--accept"])
    assert "--out" in str(err.value)
