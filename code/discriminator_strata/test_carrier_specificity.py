from __future__ import annotations

import pytest

import carrier_specificity as cs


def test_every_member_builds_with_declared_counts() -> None:
    expected = {
        "icosahedron": (12, 30),
        "cuboctahedron": (12, 24),
        "truncated_tetrahedron": (12, 18),
        "hexagonal_prism": (12, 18),
        "octahedron": (6, 12),
        "cube": (8, 12),
        "dodecahedron": (20, 30),
        "tetrahedron": (4, 6),
        "petersen": (10, 15),
        "complete_k12": (12, 66),
        "cycle_c12": (12, 12),
        "complete_bipartite_k66": (12, 36),
    }
    for name, (vertices, edges) in expected.items():
        row = cs.probe_member(name)
        assert (row["vertices"], row["edges"]) == (vertices, edges), name


def test_icosahedron_is_the_unique_full_hit() -> None:
    receipt = cs.build_receipt()
    assert receipt["full_hit_members"] == ["icosahedron"]
    assert receipt["unique_full_hit"] is True
    assert receipt["ensemble_size"] == 12


def test_cuboctahedron_shares_multiplicities_and_fails_elsewhere() -> None:
    row = cs.probe_member("cuboctahedron")
    assert row["band_multiplicities"] == [1, 3, 3, 5]
    assert row["probes"]["multiplicities_1_3_3_5"] is True
    assert row["probes"]["degree_five_regular"] is False
    assert row["probes"]["galois_band_pair"] is False
    assert row["probes"]["response_involution"] is False
    assert row["full_hit"] is False


def test_per_probe_counts_are_stable() -> None:
    receipt = cs.build_receipt()
    assert receipt["per_probe_pass_counts"] == {
        "twelve_ports": 7,
        "degree_five_regular": 1,
        "four_bands": 3,
        "galois_band_pair": 3,
        "multiplicities_1_3_3_5": 2,
        "unique_antipode_involution": 7,
        "twenty_triangles": 1,
        "response_involution": 1,
    }


def test_galois_pair_members() -> None:
    for name, expected in (
        ("icosahedron", True),
        ("dodecahedron", True),
        ("cycle_c12", True),
        ("petersen", False),
        ("cube", False),
    ):
        row = cs.probe_member(name)
        assert row["probes"]["galois_band_pair"] is expected, name


def test_response_law_family_counts_and_declared_recovery() -> None:
    laws = cs.response_law_family()
    assert laws["law_count"] == 16
    assert laws["rational_law_count"] == 8
    assert laws["declared_coefficients_recovered"] is True
    declared = [row for row in laws["laws"] if row["declared_law"]]
    assert len(declared) == 1
    assert declared[0]["rational_coefficients"] is True
    assert declared[0]["coefficients_low_to_high"] == [
        "1+0*sqrt5",
        "-1/2+0*sqrt5",
        "-2/5+0*sqrt5",
        "1/10+0*sqrt5",
    ]
    for row in laws["laws"]:
        symmetric = row["sign_pattern"][1] == row["sign_pattern"][2]
        assert row["rational_coefficients"] is symmetric


def test_committed_receipt_is_byte_exact() -> None:
    committed = cs.OUTPUT_PATH.read_bytes()
    assert committed == cs.canonical_json_bytes(cs.build_receipt())


def test_duplicate_edge_fails_closed() -> None:
    with pytest.raises(cs.SpecificityError, match="duplicate edge"):
        cs.adjacency_matrix([("a", "b"), ("b", "a")])


def test_minimal_polynomial_on_known_graph() -> None:
    matrix = cs.adjacency_matrix(cs.complete_k12())
    minimal = cs.minimal_polynomial(matrix)
    assert len(minimal) - 1 == 2
    structure = cs.spectral_structure(matrix)
    assert structure["rational_bands"] == ["-1", "11"]
