"""Exact affine span, Hessian sign, and frozen modular work on the native volume level."""

import json
from itertools import combinations
from math import prod

import numpy as np
from mpmath import iv

from ....full_algebra.enclosure import geometry_jet
from ....full_algebra.station import ORDER
from ..geometry import HERE as BALLS, interval_numbers, raw_bounds

# Independent circulations on the four native stations. Reversal negates the exponent.
UPPER = (1, 2, -3, 5, -4, 7)


def certificate():
    ring = np.zeros((4, 4), dtype=int)
    for (i, j), a in zip(combinations(range(4), 2), UPPER, strict=True):
        ring[i, j] = a
        ring[j, i] = -a
    assert np.all(ring.sum(axis=1) == 0)
    exponents = ring[ORDER]
    assert sorted(abs(exponents)) == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 7, 7]
    vandermonde = prod(map(int, exponents)) * prod(
        int(b - a) for i, a in enumerate(exponents) for b in exponents[i + 1 :]
    )
    assert vandermonde != 0

    # The previous state-independent right inverse controls these exact same
    # twelve native effects, already on the 18 marginal-preserving correlations.
    algebra = BALLS.parents[1] / "full_algebra"
    inv = json.loads((algebra / "reuse" / "invisible.json").read_text())
    assert inv["range_rows"] == 12 and inv["exact_range_kernel_dimension"] == 6
    assert inv["correlation_coordinates"] == 18
    prior = json.loads((algebra / "certificate.json").read_text())
    assert prior["multiplier_drift_upper"] < 0.043
    eta = iv.mpf(prior["gain_interval"])
    lam = iv.mpf("51.0886428") + iv.mpf(["-.043", ".043"])

    radar = raw_bounds()[ORDER]
    area, volume = geometry_jet(radar)
    sign = np.array([1, 1, -1, -1])[np.asarray(ORDER[0])]
    velocity = sign * radar
    direct = velocity @ volume.h @ velocity
    target = iv.mpf([str(v) for v in prior["volume_interval"]])
    assert volume.v.a <= target.a and target.b <= volume.v.b
    hess = -3 * target
    assert hess.a >= direct.a and hess.b <= direct.b
    charge = 3 * lam * target - 2 * eta * area.v
    assert charge.b < -230
    return dict(
        native_range_rank=12,
        trace_one_dimension=31,
        range_fiber_dimension=19,
        directed_circulation_exponents=list(map(int, exponents)),
        derivative_vandermonde=str(vandermonde),
        finite_affinely_independent_preparations=32,
        area=interval_numbers(area.v),
        volume=prior["volume_interval"],
        direct_hessian_interval=interval_numbers(direct),
        gauge_hessian_interval=interval_numbers(hess),
        frozen_modular_coefficient=interval_numbers(charge),
    )


if __name__ == "__main__":
    print(json.dumps(certificate(), indent=2))
