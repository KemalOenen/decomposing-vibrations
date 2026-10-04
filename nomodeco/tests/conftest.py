from pathlib import Path

import pytest

from nomodeco.tests import molecules

REPO = Path(__file__).resolve().parents[2]
# Moves to nomodeco/tests/data in step 0.3 of the restructure plan
ORCA_DATA = REPO / "test_calculations" / "orca_tests"
GAUSSIAN_DATA = REPO / "test_calculations" / "gaussian_tests"


@pytest.fixture(autouse=True)
def _run_in_tmp(tmp_path, monkeypatch):
    """Every test runs in its own tmp dir, so .out/.png/.csv files never land in the repo."""
    monkeypatch.chdir(tmp_path)


@pytest.fixture(params=sorted(molecules.ALL))
def mol_name(request):
    return request.param
