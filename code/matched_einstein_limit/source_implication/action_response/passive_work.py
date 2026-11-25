"""Energy attainable by collective passive rearrangement of independent parked heads."""

import argparse
from functools import lru_cache
from math import comb, log
import json

import numpy as np
import sympy as sp
from scipy.optimize import brentq

from .measure import HERE, preparations
from ..common_front import certificate as front


def occupancy_types(copies, levels):
    if levels == 1:
        yield (copies,), 1
        return
    for count in range(copies + 1):
        for tail, multiplicity in occupancy_types(copies - count, levels - 1):
            yield (count, *tail), multiplicity * comb(copies, count)


@lru_cache(maxsize=None)
def type_arrays(copies):
    types, counts = zip(*occupancy_types(copies, 8))
    return np.array(types), np.array(counts)


def sorted_energy(copies, populations, derivatives):
    """Sort eigenvalues against physical energies, retaining multiplicities of both."""
    types, counts = type_arrays(copies)
    logs = types @ np.log(populations)
    derivatives = types @ (derivatives / populations[:, None])
    order = np.argsort(-logs, kind="stable")
    masses = np.exp(logs[order])
    # Cut equal-energy output shells without expanding individual words.
    boundaries = np.cumsum([4**copies * comb(copies, m) for m in range(copies + 1)])
    assert int(boundaries[-1]) == 8**copies and int(counts.sum()) == 8**copies
    starts = np.r_[0, np.cumsum(counts[order])[:-1]]
    stops = starts + counts[order]
    weights = np.zeros(len(order))
    for m in range(1, copies + 1):
        overlap = np.maximum(
            np.minimum(stops, boundaries[m]) - np.maximum(starts, boundaries[m - 1]), 0
        )
        weights += m * overlap / copies
    weighted = weights * masses
    return np.r_[weighted.sum(), weighted @ derivatives[order]]


def run():
    x = preparations()
    chi = x["equilibrium"]
    turn = (x["plus"] - x["minus"]) / 2
    rows = json.loads((HERE.parent / "signed_transfers.json").read_text())["tests"]
    tangent = [
        x[label] - chi + row["angle_per_unit"] * turn
        for label, row in zip(("ordinary", "pairing"), rows)
    ]
    populations = []
    derivatives = []
    for sector, rho in enumerate(chi):
        w, v = np.linalg.eigh(rho)
        populations.extend(w)
        derivatives.extend(
            np.column_stack([np.diag(v.conj().T @ d[sector] @ v).real for d in tangent])
        )
    pop = np.array(populations)
    der = np.array(derivatives)
    slopes = -np.log(pop) @ der
    entropy = -pop @ np.log(pop)

    def h_binary(p):
        return -p * log(p) - (1 - p) * np.log1p(-p)

    q = brentq(lambda p: log(4) + h_binary(p) - entropy, 1e-6, 0.5)
    beta = log((1 - q) / q)
    # Sylvester tests on the exact rational center, widened by the certified
    # radius for the Hermitian source. Failure via a strictly negative minor
    # gives a strictly negative Rayleigh direction.
    mid, _, _, _, radius = front.inputs()
    eye = sp.eye(4)

    def signs(block):
        return [sp.sign(block[:k, :k].det()) for k in range(1, 5)]

    even, odd = mid
    assert signs(even - (sp.Rational(".10133") + radius) * eye) == [1] * 4
    assert -1 in signs(even - (sp.Rational(".10135") - radius) * eye)
    assert signs((sp.Rational(".07730") - radius) * eye - odd) == [1] * 4
    assert -1 in signs((sp.Rational(".07726") + radius) * eye - odd)
    assert -1 in signs((sp.Rational(".29130") + radius) * eye - even)
    guaranteed_one_copy_gap = 0.02402
    guaranteed_two_copy_work = 0.01222
    assert sp.Rational(".10133") - sp.Rational(".07730") > sp.Rational(".02402")
    assert sp.Rational(".29130") * sp.Rational(".07726") - sp.Rational(
        ".10135"
    ) ** 2 > sp.Rational(".01222")
    data = {}
    for n in (1, 2, 4, 8, 12):
        energy, ordinary, held = sorted_energy(n, pop, der)
        predicted_ratio = (held / ordinary) * (slopes[0] / slopes[1])
        data[str(n)] = dict(
            least_excitation_energy_per_copy=energy,
            ordinary_derivative=ordinary,
            pairing_derivative=held,
            ordinary_calibrated_pairing_fraction=predicted_ratio,
        )
    return dict(
        native_excitation_probability=float(chi[1].trace().real),
        equilibrium_entropy=entropy,
        one_copy_population_order_margin_lower=guaranteed_one_copy_gap,
        two_copy_one_swap_work_lower=guaranteed_two_copy_work,
        finite=data,
        infinite=dict(
            least_excitation_energy=q,
            beta=beta,
            ordinary_derivative=slopes[0] / beta,
            pairing_derivative=slopes[1] / beta,
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    record = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "passive_work.json").write_text(record)
    else:
        print(record, end="")
