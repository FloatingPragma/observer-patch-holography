"""Test the flat diamond stress kernel on one family of native spectral transfers."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np

from ..experiment import endpoint_positions
from ..memories import THETA
from ..native_pair import LAW, LOGICAL, NativePair, entropy, sharp_counts

HERE = Path(__file__).resolve().parent
RATIOS = (1, 4, 16)
COEFFICIENTS = (4, -17, 4)
SEPARATION = 4**5
SOURCE_RADII_START = Fraction(1, 512)


def radius(shell):
    return SOURCE_RADII_START / SEPARATION**shell


def kernel_over_pi(ball, distance):
    """Flat rest-frame kernel including its support. Multiply by pi for B."""
    if distance >= ball:
        return 0
    return ball - distance**2/ball


def annular_kernel_over_pi(ball, distance):
    return sum(c*kernel_over_pi(q*ball, distance)
               for c, q in zip(COEFFICIENTS, RATIOS, strict=True))


def normalized_lower_bound():
    """For an entropy slope of one times R**4 on each of the three balls."""
    norm = sum(abs(c)*q**4 for c, q in zip(COEFFICIENTS, RATIOS, strict=True))
    return Fraction(abs(sum(COEFFICIENTS)), norm)


def _tickets(probabilities):
    """Verify and express a spectral or load distribution as integer sixteenths."""
    tickets = np.rint(16*np.asarray(probabilities)).astype(int)
    np.testing.assert_allclose(tickets/16, probabilities, rtol=0, atol=4e-13)
    return tickets.tolist()


def native_transfer(low, side):
    """Replay the two native preparations, their conserved observables and undo read."""
    pair = NativePair(low, side)
    gate = pair.gate(THETA)
    high, low_word = gate[:, 0], gate[:, 3]
    tangent = np.outer(low_word, low_word.conj()) - np.outer(high, high.conj())
    vacuum = pair.state(THETA)
    finite_shift = 1/16
    changed = vacuum + finite_shift*tangent

    energy, axes = np.linalg.eigh(pair.generator)
    levels = np.rint(energy).astype(int)
    np.testing.assert_allclose(energy, levels, rtol=0, atol=4e-13)
    assert set(levels) == {0, 8, 9}
    energy_laws = []
    for state in (vacuum, changed):
        diagonal = np.diag(axes.conj().T @ state @ axes).real
        energy_laws.append(_tickets([sum(diagonal[levels == level]) for level in (0, 8, 9)]))
    assert energy_laws[0] == energy_laws[1] == [4, 6, 6]

    charges = np.array([sum(pair.physical(*word)) for word in LOGICAL])
    distinct_loads = sorted(set(charges))
    load_laws = []
    for state in (vacuum, changed):
        diagonal = np.diag(state).real
        load_laws.append(_tickets([sum(diagonal[charges == load]) for load in distinct_loads]))
    assert load_laws[0] == load_laws[1] == [2, 12, 2]

    eigen, basis = np.linalg.eigh(vacuum)
    log_vacuum = (basis*np.log(eigen)) @ basis.conj().T
    slope = -float(np.trace(tangent @ log_vacuum).real)
    assert math.isclose(slope, math.log(5), abs_tol=4e-13, rel_tol=0)
    finite_gain = entropy(changed) - entropy(vacuum)
    assert 0 < finite_shift*slope - finite_gain <= (48/5)*finite_shift**2

    exposures = [pair.undo(state, THETA) for state in (vacuum, changed)]
    inverse_laws = [_tickets(np.diag(state).real) for state in exposures]
    assert inverse_laws == [[10, 2, 2, 2], [9, 2, 2, 3]]
    seed = 261126 + 4*low + 2*int(side == "right")
    shots = 16384
    counts = [sharp_counts(state, seed+i, shots) for i, state in enumerate(exposures)]
    return dict(low_load=low, reversed_port=side, transfer_load=int(charges[0]),
                eigenwords=[list(pair.physical(*LOGICAL[i])) for i in (0, 3)],
                undo_laws_in_sixteenths=inverse_laws, energy_levels=[0, 8, 9],
                energy_laws_in_sixteenths=energy_laws,
                load_levels=[int(q) for q in distinct_loads], load_laws_in_sixteenths=load_laws,
                entropy_slope=round(slope, 12), finite_entropy_gain=round(finite_gain, 12),
                seed_for_first_bank=seed, shots_per_bank=shots, undo_word_counts=counts)


def finite_bands(depth):
    """Restrict one common product preparation to each nested triple."""
    if depth <= 0:
        raise ValueError("at least one active pair is required")
    shells = [radius(k) for k in range(depth)]
    pairs = [np.linalg.norm(endpoint_positions("rest", float(r/4)), axis=1) for r in shells]
    rows = []
    for k, r in enumerate(shells):
        masks = []
        entropy_weights = []
        for q in RATIOS:
            masks.append([[bool(d < float(q*r)) for d in endpoints] for endpoints in pairs])
            assert all(a == b for a, b in masks[-1])
            inside = [i for i, (a, _) in enumerate(masks[-1]) if a]
            assert inside == list(range(k, depth))
            entropy_weights.append(sum((shells[i]**4 for i in inside), Fraction()))
        assert len(set(entropy_weights)) == 1
        weight = entropy_weights[0] / r**4
        rows.append(dict(radius=str(r), perturbed_pair_scale=str(r/4),
                         contained_active_pairs=[i for i,(a,_) in enumerate(masks[0]) if a],
                         common_slope_in_units_of_R4_log5=float(weight),
                         minimum_normalized_error=float(normalized_lower_bound()*weight)*math.log(5)))
    return rows


def deposit():
    tail = 1/(1-Fraction(1, SEPARATION)**4)
    peak = float(tail)*math.log(5)
    return dict(schema="oph.published_string_nested_stress.v1", native_dimension=LAW.DIM,
                native_transfers=[native_transfer(low, side) for low in (0, 1) for side in ("left", "right")],
                radii_ratios=list(RATIOS), annular_coefficients=list(COEFFICIENTS),
                shell_radius_ratio=SEPARATION, active_weight="R_k^4",
                finite_prefix=finite_bands(4),
                exact_uniform_error_bound_in_units_of_log5=str(normalized_lower_bound()),
                uniform_error_bound=float(normalized_lower_bound())*math.log(5),
                exact_infinite_tail_factor=str(tail),
                entropy_slopes_per_radius_fourth_infinite=[peak/q**4 for q in RATIOS],
                minimum_limsup_error_for_regular_stress=peak*(1-1/16**4)/2,
                minimum_annular_negative_variation_per_R3=9*math.log(5)/(60*math.pi))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = deposit()
    path = HERE/"results.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    elif json.loads(path.read_text()) != data:
        raise ValueError("native transfer or three-ball archive changed")
    else:
        print("Replayed native spectral transfers and four compatible three-ball bands.")


if __name__ == "__main__":
    main()
