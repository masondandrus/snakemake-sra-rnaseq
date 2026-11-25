import json

import pytest

RATGTEX = (
    b"tissue\tgene_id\tqval\tvariant_id\n"
    b"NAcc\tENSRNOG00000000001.5\t0.01\tchr1:100\n"
    b"NAcc\tENSRNOG00000000001.5\t0.001\tchr1:200\n"
    b"NAcc\tENSRNOG00000000002\t0.20\tchr1:300\n"
    b"PL\tENSRNOG00000000003\t0.01\tchr2:100\n"
)
GTEX_BULK = (
    b"gene_id\tqval\tvariant_id\n"
    b"ENSG00000000001.12\t0.001\tchr1_1_A_G_b38\n"
    b"ENSG00000000002.3\t0.30\tchr1_2_A_G_b38\n"
)


def test_ratgtex_source_returns_genes_for_one_tissue(fetch_eqtl, monkeypatch):
    monkeypatch.setattr(fetch_eqtl, "http_get", lambda url: RATGTEX)
    table, src = fetch_eqtl.from_ratgtex("NAcc", 0.05, None)
    assert list(table["gene_id"]) == ["ENSRNOG00000000001"]
    assert table["eqtl_qval"].iloc[0] == 0.001
    assert table["eqtl_variant"].iloc[0] == "chr1:200"


def test_gtex_bulk_file_source_runs(fetch_eqtl, monkeypatch):
    monkeypatch.setattr(fetch_eqtl, "http_get", lambda url: GTEX_BULK)
    table, src = fetch_eqtl.from_gtex("x", 0.05, "https://example.org/egenes.txt")
    assert list(table["gene_id"]) == ["ENSG00000000001"]


def test_species_is_fixed_by_gtex_and_ratgtex(fetch_eqtl):
    assert fetch_eqtl.resolve_species("ratgtex", None) == "Rattus norvegicus"
    assert fetch_eqtl.resolve_species("gtex", "Homo sapiens") == "Homo sapiens"


def test_a_wrong_declared_species_stops_the_script(fetch_eqtl):
    with pytest.raises(SystemExit) as err:
        fetch_eqtl.resolve_species("ratgtex", "Homo sapiens")
    assert "Rattus norvegicus" in str(err.value)


def test_file_source_needs_and_uses_the_declared_species(fetch_eqtl, tmp_path):
    with pytest.raises(SystemExit):
        fetch_eqtl.resolve_species("file", None)
    assert fetch_eqtl.resolve_species("file", "Mus musculus") == "Mus musculus"


def test_main_writes_axis_and_provenance_for_a_file_source(fetch_eqtl, tmp_path, monkeypatch):
    src = tmp_path / "eqtl.tsv"
    src.write_text("gene_id\teqtl_qval\nENSMUSG00000000001\t0.01\nENSMUSG00000000002\t0.5\n")
    out, prov = tmp_path / "axis.tsv", tmp_path / "axis.json"
    monkeypatch.setattr("sys.argv", [
        "fetch_eqtl.py", "--axis", "eQTL_mouse", "--source", "file", "--region", "NAc",
        "--url", str(src), "--species", "Mus musculus",
        "--out", str(out), "--provenance", str(prov),
    ])
    fetch_eqtl.main()
    record = json.loads(prov.read_text())
    assert record["species"] == "Mus musculus"
    assert record["n_genes"] == 1
    assert out.read_text().splitlines()[0].split("\t") == [
        "gene_id", "eQTL_mouse_qval", "eQTL_mouse_variant"]
