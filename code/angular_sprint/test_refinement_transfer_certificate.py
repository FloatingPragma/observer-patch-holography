from __future__ import annotations

import json

import pytest

import refinement_transfer_certificate as rtc


@pytest.fixture(scope="module")
def receipt() -> dict:
    return rtc.build_receipt()


def test_receipt_schema_status_and_self_digest(receipt: dict) -> None:
    assert receipt["schema"] == "oph.refinement_transfer_receipt.v1"
    assert receipt["issue"] == 643
    assert receipt["status"] == (
        "STATIC_BASE_PORT_TRANSFER_NONIDENTIFIABLE__"
        "CANONICAL_BAND_LEVEL_ONE_REFINEMENT_IDENTIFIABLE__"
        "PHYSICAL_TRANSFER_OPEN"
    )
    body = {key: value for key, value in receipt.items() if key != "receipt_sha256"}
    assert receipt["receipt_sha256"] == rtc.tagged_sha256(
        rtc.canonical_json_bytes(body)
    )


def test_committed_receipt_is_byte_exact(receipt: dict) -> None:
    committed = rtc.RECEIPT_PATH.read_bytes()
    assert committed == rtc.canonical_json_bytes(receipt)


def test_harmonic_bases_dimensions_harmonicity_and_rank() -> None:
    for degree, expected in ((0, 1), (1, 3), (2, 5), (3, 7), (6, 13), (10, 21), (12, 25)):
        monos, basis = rtc.harmonic_basis(degree)
        assert len(monos) == (degree + 1) * (degree + 2) // 2
        assert len(basis) == expected
        for vec in basis:
            assert len(vec) == len(monos)
            image = rtc.apply_laplacian(vec, degree)
            assert all(coeff == 0 for coeff in image)
        rows = [[rtc.F5(coeff) for coeff in vec] for vec in basis]
        assert rtc.rank_f5(rows) == expected


def test_vertex_norms_adjacency_counts_and_edge_total() -> None:
    verts = rtc.icosahedron_vertices()
    assert len(verts) == 12
    for vertex in verts:
        assert rtc.dot(vertex, vertex) == rtc.R_VERTEX
    pairs = rtc.adjacent_pairs(verts)
    assert len(pairs) == 30
    for index in range(12):
        assert sum(1 for a, b in pairs if index in (a, b)) == 5


def test_midpoint_radius_identity_and_unit_projection() -> None:
    assert rtc.ONE_PLUS_SQRT5 * rtc.ONE_PLUS_SQRT5 == rtc.R_MIDPOINT
    verts = rtc.icosahedron_vertices()
    units = rtc.unit_midpoints(verts)
    assert len(units) == 30
    assert len(set(units)) == 30
    for point in units:
        assert rtc.dot(point, point) == rtc.ONE


def test_vertex_radius_nonsquare_guard() -> None:
    facts = rtc.verify_vertex_radius_nonsquare(rtc.R_VERTEX)
    assert facts["galois_norm"] == "5"
    assert facts["square_in_q_sqrt5"] is False
    with pytest.raises(rtc.CertificateError):
        rtc.verify_vertex_radius_nonsquare(rtc.F5(4))
    with pytest.raises(rtc.CertificateError):
        rtc.verify_vertex_radius_nonsquare(rtc.F5(6, 2))  # (1+sqrt5)^2


def test_kernel_dimensions_are_pinned(receipt: dict) -> None:
    kernels = receipt["kernels"]
    assert kernels["canonical_band_vertex_only"]["space_dimension"] == 16
    assert kernels["canonical_band_vertex_only"]["condition_rows"] == 24
    assert kernels["canonical_band_vertex_only"]["kernel_dimension"] == 4
    assert kernels["canonical_band_even_part_vertex_only"]["space_dimension"] == 6
    assert kernels["canonical_band_even_part_vertex_only"]["kernel_dimension"] == 0
    assert kernels["canonical_band_odd_part_vertex_only"]["space_dimension"] == 10
    assert kernels["canonical_band_odd_part_vertex_only"]["kernel_dimension"] == 4
    assert kernels["canonical_band_level_one"]["condition_rows"] == 54
    assert kernels["canonical_band_level_one"]["kernel_dimension"] == 0
    assert kernels["comb_band_level_one"]["space_dimension"] == 60
    assert kernels["comb_band_level_one"]["condition_rows"] == 54
    assert kernels["comb_band_level_one"]["kernel_dimension"] == 39
    assert kernels["canonical_band_midpoints_only"]["condition_rows"] == 30
    assert kernels["canonical_band_midpoints_only"]["kernel_dimension"] == 0


def test_parity_split_sanity_is_pinned(receipt: dict) -> None:
    sanity = receipt["parity_split_sanity"]
    assert sanity["canonical_band_zero_vertex_rows"] == 0
    assert sanity["canonical_even_part_zero_vertex_rows"] == 12
    assert sanity["canonical_odd_part_zero_vertex_rows"] == 12
    assert sanity["comb_band_zero_vertex_rows"] == 12


def test_midpoint_corruption_fails_closed() -> None:
    verts = rtc.icosahedron_vertices()
    units = list(rtc.unit_midpoints(verts))
    point = units[0]
    axis = next(k for k in range(3) if not point[k].is_zero())
    corrupted = list(point)
    corrupted[axis] = corrupted[axis] * rtc.F5(2)
    units[0] = (corrupted[0], corrupted[1], corrupted[2])
    with pytest.raises(rtc.CertificateError):
        rtc.midpoint_rows(rtc.CANONICAL_DEGREES, units)


def _assert_no_floats(value: object) -> None:
    if isinstance(value, dict):
        for key, item in value.items():
            assert not isinstance(key, float)
            _assert_no_floats(item)
    elif isinstance(value, list):
        for item in value:
            _assert_no_floats(item)
    else:
        assert not isinstance(value, float)


def test_receipt_contains_no_float_values(receipt: dict) -> None:
    _assert_no_floats(receipt)
    _assert_no_floats(json.loads(rtc.RECEIPT_PATH.read_text(encoding="ascii")))
