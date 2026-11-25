import importlib.util
import pathlib

import pytest

SCRIPTS = pathlib.Path(__file__).resolve().parents[2] / "workflow" / "scripts"


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def fetch_eqtl():
    return load("fetch_eqtl")


@pytest.fixture
def fetch_gwas():
    return load("fetch_gwas")


@pytest.fixture
def curate():
    return load("curate_samples")
