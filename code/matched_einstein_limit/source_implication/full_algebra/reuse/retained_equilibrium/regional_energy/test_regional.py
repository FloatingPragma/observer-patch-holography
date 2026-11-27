import json
from fractions import Fraction

import numpy as np

from . import certificate, native, ward


def test_integer_product_with_both_signs_and_large_coordinates():
    values = np.array(
        [-(2**54) - 1445, -(2**20), -1, 0, 1, 2**29 + 3, 2**55 - 1], dtype=np.int64
    )
    for c in (-(2**46) - 8, -(2**23), 0, 4, 2**46 + 17):
        np.testing.assert_array_equal(
            native.scaled_product(values, c), [int(v) * c // 2**46 for v in values]
        )


def test_independent_energy_and_momentum_certificates():
    data = certificate.verify()
    assert data == json.loads((native.HERE / "results.json").read_text())
    prior = json.loads(
        (certificate.SOURCE / "directional_stress/translation/results.json").read_text()
    )
    ordinary_energy = sum(Fraction(v) for v in (".297562881119", ".297563446805")) / 2
    for q, result in data["tracers"].items():
        assert result["profile_isentropic_upper_both_marks_peak"] < -16_000_000
        assert result["both_mark_profiles_abs_collar_over_ordinary_upper"] < 590_000
        assert result["full_flight_collar"][0] < 0 < result["full_flight_collar"][1]
        for mode in ("ordinary", "pairing"):
            reference = prior["signs"][q]["regional_energy"][mode]
            expected = (
                reference + (float(ordinary_energy) if mode == "ordinary" else 0)
            ) * 1e8
            key = "twelve_cell_ordinary" if mode == "ordinary" else "twelve_cell_collar"
            assert result[key][0] < expected < result[key][1]
    assert data["tracers"]["1"]["profile_isentropic_upper_one_mark_peaks"] < -3700000
    assert data["fixed_vacancy_region"]["excess_energy"][1] < -100000
    assert data["fixed_vacancy_region"]["source_entropy"][0] > 22_900_000
    assert ward.evidence() == json.loads((native.HERE / "ward_counts.json").read_text())
