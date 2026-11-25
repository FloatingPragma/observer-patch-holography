import copy
import json

import numpy as np
import pytest

from .. import station
from . import measure


def record():
    return json.loads((measure.HERE / "results.json").read_text())


def test_exact_rank_and_common_root_control():
    source = record()
    measure.iv.dps = 32
    result = measure.certify(source)
    assert result["minimum_total_product_terms"] == 5
    assert result["chart_subdivisions"] == source["chart_subdivisions"]
    # Duplicated parity effects have a common rank drop at any even root.
    duplicate = copy.deepcopy(source)
    for name in ("area", "volume"):
        duplicate[name][1] = duplicate[name][0]
    with pytest.raises(ValueError, match="uncertified pencil"):
        measure.certify(duplicate)


def test_five_term_reader_against_independent_native_ranging():
    source = record()
    reader = source["signed_reader"]
    stored = json.loads((measure.HERE.parent / "results.json").read_text())
    state = np.array(stored["x_real"]) + 1j * np.array(stored["x_imag"])
    effects = station.radar_menus()[6]
    area, volume = station.range_gradient(station.pair(effects, state))
    gradient = -stored["ordinary_eta"] * (area - float(reader["volume_price"]) * volume)
    panel = np.einsum("p,pbij->bij", gradient, effects)
    sigmas = [
        np.array([[0, 1], [1, 0]]),
        np.array([[0, -1j], [1j, 0]]),
        np.diag([1, -1]),
    ]
    local = [np.kron(p, np.eye(2)) / 2 for p in sigmas] + [
        np.kron(np.eye(2), p) / 2 for p in sigmas
    ]
    assembled = []
    for sector in range(2):
        dyads = reader["dyads"][sector]
        assert len(dyads) == 2 + sector
        block = sum(
            float(t) * axis for t, axis in zip(reader["one_head"][sector], local)
        )
        for term in dyads:
            left, right = [
                sum(float(v) * p for v, p in zip(term[key], sigmas))
                for key in ("left", "right")
            ]
            block += float(term["weight"]) * np.kron(left, right) / 2
        block += (
            ((-1) ** sector) * float(reader["parity_contrast"]) * np.eye(4) / np.sqrt(8)
        )
        assembled.append(block)
    panel -= np.eye(4) * np.trace(panel, axis1=1, axis2=2).sum() / 8
    assert (
        np.linalg.norm(panel - assembled)
        < reader["error_upper_times_scale"] / source["scale"]
    )
    # The state spectrum supplies an independent entropy derivative at the
    # published central preparation. Restrict it to the same volume tangent.
    def traceless(blocks):
        return blocks - np.trace(blocks, axis1=1, axis2=2).sum() * np.eye(4) / 8

    dv = traceless(np.einsum("p,pbij->bij", volume, effects))
    difference = traceless(-station.log_blocks(state) - assembled)
    orthogonal = difference - dv * np.vdot(dv, difference).real / np.vdot(dv, dv).real
    assert np.linalg.norm(orthogonal) < 3.106e-7
