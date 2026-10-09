"""
IC generation must not depend on PYTHONHASHSEED (set iteration order), and every place that
cuts the search short must show up in the search report.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

from nomodeco.libraries import icsel, logfile
from nomodeco.tests import molecules
from nomodeco.tests.helpers import analyse

REPO_ROOT = Path(__file__).resolve().parents[3]

# prints every generated IC set of ethylene (before the dihedral sort, a different order for
# every hash seed)
SCRIPT = """
from nomodeco.tests import molecules
from nomodeco.tests.helpers import generate_sets
ic_dict, _ = generate_sets(molecules.ethylene())
print(repr([ic_dict[k] for k in range(len(ic_dict))]))
"""


def run_with_hash_seed(seed):
    env = dict(os.environ, PYTHONHASHSEED=str(seed))
    result = subprocess.run([sys.executable, "-c", SCRIPT], cwd=REPO_ROOT, env=env,
                            capture_output=True, text=True, check=True)
    return result.stdout.strip().splitlines()[-1]


@pytest.mark.slow
def test_ic_sets_do_not_depend_on_hash_seed():
    assert run_with_hash_seed(1) == run_with_hash_seed(2)


def test_dihedrals_are_in_atom_index_order():
    mol = molecules.acetyl_cyanide()
    index = {a.symbol: i for i, a in enumerate(mol)}
    dihedrals = analyse(mol).dihedrals
    assert dihedrals
    keys = [tuple(index[a] for a in d) for d in dihedrals]
    assert keys == sorted(keys)


def test_subset_cap_is_reported(monkeypatch):
    monkeypatch.setattr(icsel, "_MAX_SUBSETS", 3)
    report = logfile.SearchReport()
    logfile.search_log.addHandler(report)
    try:
        icsel._flat_subsets_of_size([[i] for i in range(6)], 2)
    finally:
        logfile.search_log.removeHandler(report)
    assert len(report.messages) == 1
    assert "stopped at 3 subsets of size 2" in report.messages[0]


def test_untruncated_search_reports_nothing():
    report = logfile.SearchReport()
    logfile.search_log.addHandler(report)
    try:
        icsel._flat_subsets_of_size([[i] for i in range(6)], 2)
    finally:
        logfile.search_log.removeHandler(report)
    assert report.messages == []
