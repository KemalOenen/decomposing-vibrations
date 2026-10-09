"""
topology.LazyIcDict: indexing must agree with iteration across all product blocks.
"""

import itertools
import logging

from nomodeco.libraries import topology
from nomodeco.libraries.topology import LazyIcDict

BLOCK_1 = (["b1"], ["la1"], [["a1"], ["a2"]], [["o1"], ["o2"], ["o3"]], [["d1"], ["d2"]])
BLOCK_2 = (["b2"], [], [["a3"]], [[]], [["d3"], ["d4"]])


def expected_sets(*blocks):
    for bonds, la, a_s, o_s, d_s in blocks:
        for ang, oop, dih in itertools.product(a_s, o_s, d_s):
            yield {"bonds": bonds, "angles": ang, "linear valence angles": la,
                   "out of plane angles": oop, "dihedrals": dih}


def lazy(*blocks):
    d = LazyIcDict()
    for block in blocks:
        d.add_product(*block)
    return d


def test_empty():
    d = LazyIcDict()
    assert len(d) == 0
    assert list(d.values()) == []


def test_length_and_keys():
    d = lazy(BLOCK_1, BLOCK_2)
    assert len(d) == 2 * 3 * 2 + 1 * 1 * 2
    assert list(d.keys()) == list(range(len(d)))


def test_indexing_matches_product_order_across_blocks():
    d = lazy(BLOCK_1, BLOCK_2)
    expected = list(expected_sets(BLOCK_1, BLOCK_2))
    assert [d[k] for k in d.keys()] == expected
    assert list(d.values()) == expected
    assert list(d.items()) == list(enumerate(expected))


def test_empty_block_is_skipped_with_a_warning(caplog):
    d = lazy(BLOCK_1)
    with caplog.at_level(logging.WARNING):
        d.add_product(["b3"], [], [["a4"]], [], [[]])
    assert "Empty IC product block" in caplog.text
    assert len(d) == 12
    assert [d[k] for k in d.keys()] == list(expected_sets(BLOCK_1))


def test_block_with_no_dihedrals_needed():
    """[[]] (choose nothing) is one selection, not zero."""
    d = lazy((["b1"], [], [["a1"]], [[]], [[]]))
    assert len(d) == 1
    assert d[0]["dihedrals"] == []


def test_total_cap_skips_further_blocks(monkeypatch):
    monkeypatch.setattr(topology, "_MAX_IC_SETS_TOTAL", 10)
    d = lazy(BLOCK_1, BLOCK_2)
    assert len(d) == 12
