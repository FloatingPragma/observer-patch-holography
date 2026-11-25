import json

import numpy as np
import pytest
from scipy import sparse

from . import held, mounts, transfer


def test_independent_mounts_and_muted_material_response():
    board, route, _, patches = mounts.apparatus()
    first = mounts.cuts.native.map_matrix(board, patches[0])
    for p in patches[1:]:
        other = mounts.cuts.native.map_matrix(board, p)
        assert (first - other).count_nonzero() and (first.T - other).count_nonzero()
    # At zero marked-head interaction, both parities have the same mobile
    # generator; the shift and every flux inherit identity on the two heads.
    for charge in (-1, 1):
        bg = mounts.cuts.native.blank_h(board, charge)
        assert (bg - sparse.kron(bg[::4, ::4], sparse.eye(4))).count_nonzero() == 0
    for p in patches:
        s = mounts.cuts.native.map_matrix(board, p)
        edges = mounts.cuts.native.flux_on_edges(board, p)
        blank = mounts.cuts.native.blank_h(board, 1)
        one = sparse.eye(blank.shape[0], format="csr")
        domain = -s @ blank @ (one - s.T @ s)
        image = (one - s @ s.T) @ blank @ s
        assert (
            sum(edges.values()) - domain - domain.T - image - image.T
        ).count_nonzero() == 0
        for op in (s - s.T, domain + domain.T, image + image.T, *edges.values()):
            assert (op - sparse.kron(op[::4, ::4], sparse.eye(4))).count_nonzero() == 0
    # On the two charge-preserving repeated preparations, a scalar block
    # readout independent of parity annihilates every trace-zero tangent.
    assert np.count_nonzero(mounts.source_components([np.eye(4), np.eye(4)])) == 0


def test_source_and_held_chart_certificates():
    r = transfer.certify()
    assert all(
        arm["rank"] == 31 and float(arm["true_modular_error_upper"]) < 3e-5
        for arm in r["instruments"]
    )
    assert float(r["instruments"][-1]["sensor_error_l2_to_source_hs_upper"]) < 2.1e9
    assert all(
        float(row["detector_vector_error_per_unit_HS_upper"]) < 1e-20
        for row in r["transfers"]
    )
    assert all(
        float(row["modular_disagreement_per_unit_HS_upper"]) < 1e-19
        for row in r["transfers"]
    )
    original = transfer.read_native()[0]
    with pytest.raises(ZeroDivisionError):
        rows = original.copy()
        rows[0] = rows[1]
        _ = transfer.center(rows) ** -1


def test_independently_flown_pooled_readers():
    report = json.loads((mounts.HERE / "transfer.json").read_text())
    r = held.run(report)
    for arm in r["instruments"][3:]:
        assert arm["reference"]["source_HS_error"] < 1e-5
        assert arm["reference"]["twelve_directed_range_max_error"] < 1e-6
        for c, src in zip(arm["cases"], r["preparations"]):
            assert c["reconstructed_source_HS_error"] < 1e-5
            assert abs(c["modular"] - src["entropy"]) < 2e-5
            assert abs(c["material_reference_modular"] - src["entropy"]) < 2e-5
            assert abs(c["modular"] - src["ordinary_calibrated_area"]) < 2e-5
            assert all(v["relative_l2_error"] < 1e-5 for v in c["held_readings"])


def test_two_nested_regional_reports_predict_the_unseen_full_mount():
    from . import apertures, regional

    patches = apertures.apparatus()[-1]
    assert [len(p["cells"]) for p in patches] == [1, 3, 1, 2, 1, 3]
    assert all(
        set(patches[i]["cells"]) < set(patches[i + 1]["cells"]) for i in (0, 2, 4)
    )
    result = regional.run()
    for arm in result["instruments"]:
        assert arm["rank"] == 31 and arm["min_singular_lower_in_1e20"] > 10**8
        for case, prep in zip(arm["cases"], result["preparations"]):
            assert (
                abs(case["material_only_modular"] - prep["ordinary_calibrated_area"])
                < 2e-5
            )
            assert abs(case["material_only_modular"] - prep["entropy"]) < 2e-5
            assert all(c["relative_l2_error"] < 1e-6 for c in case["full_held_charts"])
