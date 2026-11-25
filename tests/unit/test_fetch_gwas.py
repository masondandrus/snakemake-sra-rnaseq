import pytest

PAGE = {
    "_embedded": {
        "associations": [
            {"p_value": 1e-9, "mapped_genes": ["CHRNA5", "CHRNA3"],
             "efo_traits": [{"efo_id": "MONDO_0005186", "efo_trait": "cocaine dependence"}]},
            {"p_value": 1e-5, "mapped_genes": ["DRD2"],
             "efo_traits": [{"efo_id": "MONDO_0005186", "efo_trait": "cocaine dependence"}]},
            {"p_value": 2e-10, "mapped_genes": [],
             "efo_traits": [{"efo_id": "MONDO_0005186", "efo_trait": "cocaine dependence"}]},
        ]
    }
}


def test_only_associations_below_the_threshold_give_genes(fetch_gwas):
    rows, labels = fetch_gwas.parse_associations(PAGE, "MONDO_0005186", 5e-8)
    assert sorted(r["gene_symbol"] for r in rows) == ["CHRNA3", "CHRNA5"]
    assert labels == {"cocaine dependence"}


def test_a_non_human_species_on_the_axis_stops_the_script(fetch_gwas, monkeypatch, tmp_path):
    monkeypatch.setattr("sys.argv", [
        "fetch_gwas.py", "--traits", "MONDO_0005186", "--species", "Mus musculus",
        "--out", str(tmp_path / "a.tsv"), "--provenance", str(tmp_path / "a.json"),
    ])
    with pytest.raises(SystemExit) as err:
        fetch_gwas.main()
    assert "Homo sapiens" in str(err.value)
