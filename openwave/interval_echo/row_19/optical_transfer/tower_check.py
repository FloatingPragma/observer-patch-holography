"""Exact Einstein-frame control for the inference from rest energy to null focusing."""

from fractions import Fraction as F
import json
from pathlib import Path

from matched_einstein_limit import frame


def tensor_control():
    # Local tensor data in an orthonormal frame. Set the common density to one.
    # Each diagonal pressure has magnitude at most the density.
    cases = {}
    for sign in (-1, 1):
        diag = (F(1), F(sign, 2), F(-sign, 2), F(0))
        stress = tuple(diag[i] if i == j else F(0) for i, j in frame.PAIRS)
        timelike = frame.apply(frame.design(), stress)
        assert frame.decode(timelike) == stress
        null_x = frame.apply(frame.design(((1, 1, 0, 0),)), stress)[0]
        active_density = sum(diag)
        assert active_density == timelike[0] == 1
        assert max(map(abs, diag[1:])) <= diag[0]
        cases[str(sign)] = dict(rest_density=str(timelike[0]), density_plus_pressure_trace=str(active_density),
                                null_x=str(null_x), boosted_x=str(timelike[1]),
                                ten_tower_readings=list(map(str, timelike)))
    assert [cases[k]["null_x"] for k in ("-1", "1")] == ["1/2", "3/2"]
    assert [cases[k]["boosted_x"] for k in ("-1", "1")] == ["17/9", "11/3"]
    return cases


def native_focusing():
    # Odd forward/reverse differences test flux. Even sums test the proposed
    # longitudinal stress and energy inferred from work on the whole seam.
    path = Path(frame.__file__).parent / "source_implication/jacobi_screen/source_response.json"
    sources = json.loads(path.read_text())["species"]
    odd = []
    even = []
    bands = {}
    for q, probe in sources.items():
        flux = probe["forward_minus_reverse"]["held_focusing_over_current_prediction"]
        work = probe["forward_plus_reverse"]["held_focusing_over_whole_seam_work_prediction"]
        def endpoints(span):
            return [F(t.strip()) for t in span[1:-1].split(",")]
        lo_odd, hi_odd = endpoints(flux)
        lo_even, hi_even = endpoints(work)
        assert 0 < lo_even <= hi_even < 1 < lo_odd <= hi_odd
        odd.append(lo_odd)
        even.append(hi_even)
        bands[q] = dict(focusing_over_flux_prediction=flux, focusing_over_work_prediction=work)
    # For any positive scale a, at least one actual ratio r obeys
    # |r/a - 1| >= (H-L)/(H+L) for H=max(odd lower), L=min(even upper).
    high, low = max(odd), min(even)
    floor_ppm = (10**6 * (high-low)) // (high+low)
    assert floor_ppm >= 142_882
    return dict(probes=bands, shared_scale_relative_error_floor_ppm=floor_ppm)


if __name__ == "__main__":
    print(json.dumps(dict(tensor_control=tensor_control(), native_focusing=native_focusing()), indent=2))
