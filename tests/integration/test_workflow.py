"""Run the workflow from salmon output to the convergence table on synthetic data.

Needs R with DESeq2, tximport and WGCNA on the path. Skips when they are absent.
"""

import json
import pathlib
import re
import subprocess
import sys

import pandas as pd
import pytest

REPO = pathlib.Path(__file__).resolve().parents[2]
RULES = ["tximport", "deseq2_init", "deseq2_results", "wgcna",
         "axis_deseq2", "axis_wgcna", "axis_file", "convergence"]


def r_ready():
    try:
        return subprocess.run(
            ["Rscript", "-e", "library(DESeq2); library(tximport); library(WGCNA)"],
            capture_output=True).returncode == 0
    except FileNotFoundError:
        return False


pytestmark = pytest.mark.skipif(not r_ready(), reason="R with DESeq2, tximport and WGCNA not found")


def run(workdir, fixture_args, contrasts):
    subprocess.run([sys.executable, str(REPO / "tests/fixtures/make_fixture.py"),
                    "--out", str(workdir), *fixture_args], check=True, capture_output=True)
    targets = [f"results/convergence/{c}.convergence.tsv" for c in contrasts]
    done = subprocess.run(
        [sys.executable, "-m", "snakemake", "-s", str(REPO / "workflow/Snakefile"),
         "--directory", str(workdir), "--cores", "2", "--allowed-rules", *RULES, "--", *targets],
        capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-3000:]


def modeled_genes(workdir):
    log = (workdir / "logs/deseq2/init.log").read_text()
    return int(re.search(r"genes after filtering: (\d+)", log).group(1))


def block(gene_ids):
    return gene_ids.str[-5:].astype(int)


@pytest.fixture(scope="module")
def large(tmp_path_factory):
    workdir = tmp_path_factory.mktemp("eight_per_group")
    run(workdir, ["--design", "~ sex + condition", "--sex-contrast"],
        ["treated_vs_control", "F_vs_M"])
    return workdir


@pytest.fixture(scope="module")
def small(tmp_path_factory):
    workdir = tmp_path_factory.mktemp("three_per_group")
    run(workdir, ["--n-per-group", "3"], ["treated_vs_control"])
    return workdir


def results(workdir, contrast="treated_vs_control"):
    return pd.read_csv(workdir / f"results/deseq2/{contrast}.results.tsv", sep="\t")


def test_results_table_holds_every_modeled_gene(large, small):
    for workdir in (large, small):
        table = results(workdir)
        assert len(table) == modeled_genes(workdir)
        assert table["gene_id"].is_unique
        assert {"stat", "padj_status", "log2FoldChange_unshrunk"} <= set(table.columns)


def test_a_count_outlier_keeps_its_row_and_its_reason(small):
    table = results(small)
    row = table[table["gene_id"] == "ENSMUSG00000000300"].iloc[0]
    assert pd.isna(row["padj"]) and pd.isna(row["pvalue"])
    assert row["padj_status"] == "count_outlier"


def test_a_gene_with_a_zero_p_value_has_a_finite_statistic(large):
    row = results(large).set_index("gene_id").loc["ENSMUSG00000000270"]
    assert row["pvalue"] == 0
    assert row["stat"] > 30


def test_hubs_come_from_modules_associated_with_the_contrast(large):
    hub = pd.read_csv(large / "results/convergence/treated_vs_control.axis.HUB.tsv", sep="\t")
    b = block(hub["gene_id"])
    assert ((b < 80).sum()) >= 60
    assert ((b >= 80) & (b < 240)).sum() == 0

    hub_sex = pd.read_csv(large / "results/convergence/F_vs_M.axis.HUB.tsv", sep="\t")
    b = block(hub_sex["gene_id"])
    assert ((b >= 160) & (b < 240)).sum() >= 60
    assert (b < 160).sum() == 0


def test_module_trait_table_is_corrected_across_modules(large):
    trait = pd.read_csv(large / "results/wgcna/treated_vs_control.HUB.module_trait.tsv", sep="\t")
    assert "grey" not in set(trait["module"])
    assert (trait["padj"] >= trait["p"]).all()
    prov = json.loads((large / "results/convergence/treated_vs_control.axis.HUB.provenance.json").read_text())
    assert prov["design"] == "~sex + condition"
    assert prov["n_modules_associated"] == int((trait["padj"] < 0.05).sum())


def test_network_records_its_soft_threshold(large):
    prov = json.loads((large / "results/wgcna/provenance.json").read_text())
    assert prov["skipped"] is False
    assert prov["scale_free_r2_cut"] == 0.8
    fit = pd.read_csv(large / "results/wgcna/soft_threshold.tsv", sep="\t")
    assert prov["soft_power"] in set(fit["Power"])


def test_small_study_skips_the_network_and_says_so(small):
    prov = json.loads((small / "results/convergence/treated_vs_control.axis.HUB.provenance.json").read_text())
    assert prov["skipped"] is True
    assert "below the minimum" in prov["skip_reason"]
    conv = pd.read_csv(small / "results/convergence/treated_vs_control.convergence.tsv", sep="\t")
    assert not conv["HUB"].any()
    assert conv["DEG"].any()
    assert (conv["n_evidence"] >= 1).all()


def convergence(workdir, contrast="treated_vs_control", kind="convergence.tsv"):
    path = workdir / f"results/convergence/{contrast}.{kind}"
    if kind.endswith("json"):
        return json.loads(path.read_text())
    return pd.read_csv(path, sep="\t")


def test_a_known_overlap_is_detected_and_a_chance_overlap_is_not(large):
    overlap = convergence(large, kind="overlap.tsv").set_index("axis")
    known = overlap.loc["KNOWN"]
    assert known["reference"] == "DEG"
    assert known["n_universe"] <= 1000
    assert known["n_both"] >= 25
    assert known["odds_ratio"] > 5
    assert known["padj"] < 1e-6

    tenth = overlap.loc["TENTH"]
    assert tenth["n_universe"] == modeled_genes(large)
    assert 0.3 < tenth["odds_ratio"] < 3
    assert tenth["p_greater"] > 0.01


def test_genes_without_an_ortholog_are_not_evaluable_on_a_human_axis(large):
    prov = convergence(large, kind="provenance.json")
    known = next(a for a in prov["axes"] if a["name"] == "KNOWN")
    assert known["status"] == "ok"
    assert known["n_evaluable"] <= 1000
    assert known["n_one_to_many"] == 10
    assert prov["ortholog_tables"][0]["gprofiler_version"] == "fixture"

    table = convergence(large)
    states = table.set_index("gene_id")["KNOWN_state"]
    unmapped = [g for g in states.index if int(g[-5:]) >= 1000]
    assert unmapped and (states[unmapped] == "not_evaluable").all()


def test_a_gene_is_positive_through_its_second_ortholog(large):
    table = convergence(large).set_index("gene_id")
    row = table.loc["ENSMUSG00000000005"]
    assert row["KNOWN_state"] == "positive"
    assert bool(row["KNOWN_via_ortholog"])


def test_axis_states_and_counts_in_the_table(large, small):
    table = convergence(large)
    for axis in ("DEG", "HUB", "KNOWN", "TENTH"):
        assert set(table[f"{axis}_state"]) <= {"positive", "negative", "not_evaluable"}
        assert (table[axis] == (table[f"{axis}_state"] == "positive")).all()
    assert (table["n_evidence"] == table[["DEG", "HUB", "KNOWN", "TENTH"]].sum(axis=1)).all()

    prov = convergence(small, kind="provenance.json")
    hub = next(a for a in prov["axes"] if a["name"] == "HUB")
    assert hub["status"] == "skipped"
    assert set(convergence(small)["HUB_state"]) == {"axis_skipped"}
    assert "HUB" not in set(convergence(small, kind="overlap.tsv")["axis"])
