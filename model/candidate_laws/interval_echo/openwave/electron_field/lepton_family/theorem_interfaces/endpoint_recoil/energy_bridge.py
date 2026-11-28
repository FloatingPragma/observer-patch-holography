"""Certify shared interacting energy support and a finite Wilson energy read."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import comb, prod
import json
from pathlib import Path

import sympy as sp

from ..charged_endpoints.cube import (
    BETA_SEA,
    CELLS,
    LINKS,
    Word,
    c_action,
    charge,
    gauss,
)
from . import measure as recoil


HERE = Path(__file__).resolve().parent
CONJUGATION = "model/particle_primary_charge_conjugation/PROOF.md"
HOP = F(1, 2)


def common_support(row, dimension):
    sending = row["absolute_source_clock_diagnostics"]
    receiving = row["absolute_receiver_clock_diagnostics"]
    si, sf = sending["variance"], receiving["variance"]
    separation = sending["mean"] - receiving["mean"]
    # Use c=(mean_i+mean_f)/2 and bound the geometric mean by the arithmetic mean.
    moment_product_upper = (si + sf) / 2 + separation**2 / 4
    shared_amplitude = HOP**2 / moment_product_upper
    return {
        "shared_absolute_amplitude_lower": shared_amplitude,
        "sender_energy_in_shared_set_lower": HOP**2 / sf,
        "receiver_energy_in_shared_set_lower": HOP**2 / si,
        "sharp_energy_then_matched_receiver_lower": shared_amplitude**2 / dimension,
        "vacuum_or_particle_energy_assigned": False,
    }


def finite_bit(row, spectral_radius):
    wait = F(1, 8192)
    wait_error = F(1, 2**22)
    kraus_error = F(1, 2**22)
    wait_lo, wait_hi = wait - wait_error, wait + wait_error
    assert wait_hi * spectral_radius < F(3, 2)
    moment = row["absolute_source_clock_diagnostics"]
    transfer = row["fresh_receiver"]
    delta = abs(transfer["sender_to_receiver_diagonal_detuning"])
    remainder = transfer["third_order_amplitude_remainder_coefficient_upper"]
    raw = recoil.probability_band(
        wait_lo / 2,
        wait_hi / 2 + delta * wait_hi**2 / 4,
        remainder * wait_hi**3 + 2 * kraus_error,
    )
    # K=(I±iU)/2 has amplitude g/2 between these orthogonal endpoints.
    no_counterclock = [value / 4 for value in raw]

    low, high = wait_lo / 2, wait_hi / 2
    angle = abs(moment["mean"]) * high
    assert angle < F(1, 16)
    cos_plus_sin = 1 + angle
    cos_minus_sin = 1 - angle**2 / 2 - angle
    # For L=(cos(sH)±sin(sH))/sqrt(2), s=t/2. Bound the complex
    # Taylor remainder after combining ±s before dividing by sqrt(2).
    allowance = (
        delta * high**2 * cos_plus_sin / 4
        + 2 * remainder * high**3
        + F(3, 2) * kraus_error
    )
    lueders = [
        x / 2
        for x in recoil.probability_band(
            low * cos_minus_sin / 2, high * cos_plus_sin / 2, allowance
        )
    ]
    assert no_counterclock[0] > 0 and lueders[0] > 0
    assert moment["variance"] < 14**2
    return {
        "wait": wait,
        "absolute_wait_error_upper": wait_error,
        "each_outcome_Kraus_operator_norm_error_upper": kraus_error,
        "absolute_energy_times_wait_upper": spectral_radius * wait_hi,
        "per_pointer_and_signed_gauge_and_Bell_probability_with_common_wait": no_counterclock,
        "per_pointer_and_signed_gauge_and_Bell_probability_Lueders": lueders,
        "unread_Lueders_half_trace_disturbance_on_source_upper": 14 * high
        + 2 * kraus_error,
    }


def spin_chart(gammas):
    e = sp.eye(4)
    chart = sp.Matrix.hstack(
        e[:, 0] + e[:, 2],
        e[:, 1] + e[:, 3],
        -e[:, 0] + e[:, 2],
        -e[:, 1] + e[:, 3],
    ) / sp.sqrt(2)
    m = chart.H * gammas[1] * chart
    assert m == m.T == m.conjugate() and m**2 == sp.eye(4)
    beta = sp.diag(1, 1, -1, -1)
    assert -m.T * beta.T * m.conjugate() == beta
    for w in recoil.beta_hops(gammas):
        assert -m.T * w.T * m.conjugate() == w.H
    transform = []
    for row in range(4):
        (col,) = (col for col in range(4) if m[row, col])
        coefficient = int(m[row, col])
        assert coefficient in (-1, 1)
        transform.append((col, coefficient))
    return transform


def conjugate_word(word, transform):
    """Use the full CAR charge conjugation and linear link pullback."""
    occupied = (1 << (60 * len(CELLS))) - 1
    phase = 1
    for original in reversed(range(60 * len(CELLS))):
        if not (word.occupied >> original) & 1:
            continue
        site, orbital = divmod(original, 60)
        spin_row, internal = divmod(orbital, 15)
        new_row, coefficient = transform[spin_row]
        orbital = 60 * site + 15 * new_row + internal
        occupied, sign = c_action(occupied, orbital, False)
        phase *= sign * coefficient
    folded = Word(occupied, tuple(-e for e in word.characters))
    assert charge(folded) == {site: -q for site, q in charge(word).items()}
    assert not any(gauss(folded).values())
    return folded, phase


def check_conjugate_path(path, transform, data):
    columns = []
    for spin in (0, 1):
        source, source_phase = recoil.prepare(path, 1, spin)
        sink, raw_hop_phase = recoil.destination(source, 1, spin)
        negative_source, source_conj_phase = conjugate_word(source, transform)
        negative_sink, sink_conj_phase = conjugate_word(sink, transform)
        hole_spin = transform[spin][0] - 2
        proposed_sink, neg_hop_phase = recoil.destination(
            negative_source, -1, hole_spin
        )
        assert proposed_sink == negative_sink
        assert raw_hop_phase * sink_conj_phase == source_conj_phase * neg_hop_phase
        assert conjugate_word(negative_source, transform) == (source, source_conj_phase)
        for left, right in ((source, negative_source), (sink, negative_sink)):
            a, b = [recoil.diagnostic(word, data, F(3)) for word in (left, right)]
            assert all(a[k] == b[k] for k in ("mean", "variance", "electric_mean"))
        columns.append(
            {
                "original_positive_spin": spin,
                "holder_positive_spin_in_negative_source": transform[2][0],
                "negative_mobile_hole_spin_with_same_idle_reference": hole_spin,
                "source_CAR_phase": source_phase * source_conj_phase,
                "receiver_CAR_phase": source_phase * raw_hop_phase * sink_conj_phase,
            }
        )
    return {"path": path, "paired_columns": columns}


def run():
    sources = {
        name: recoil.read_archive(name)
        for name in (
            recoil.CONTRACT,
            recoil.GAMMA,
            recoil.M15,
            recoil.SCALE,
            CONJUGATION,
        )
    }
    member = json.loads(sources[recoil.CONTRACT])["fixed_primary_member"]
    blocks, stays, casimirs = recoil.charge_blocks(sources)
    dims = [len(ids) for ids in blocks.values()]
    volume = len(CELLS)
    link_dimension = 1 + 2 * sum(d**2 for d in dims)
    assert link_dimension == 119 and len(LINKS) == 12
    # H preserves the exact matter number of every native charge block.
    # Physical ghost links have only the Haar-constant word. Each eigenspace
    # in this binomial carrier has one label, regardless of degeneracy.
    dimension = link_dimension ** len(LINKS) * prod(
        comb(4 * d * volume, 2 * d * volume) for d in dims
    )
    gammas = recoil.clifford(sources[recoil.GAMMA])
    source = member["source"]
    data = (
        blocks,
        stays,
        casimirs,
        recoil.beta_hops(gammas),
        F(source["electric_coefficient"]),
        F(source["magnetic_coefficient"]),
    )
    transform = spin_chart(gammas)
    blank = Word(BETA_SEA, (0,) * len(LINKS))
    assert conjugate_word(blank, transform)[0] == blank
    state = member["state"]
    radius = max(abs(F(state[key])) for key in ("spectral_lower", "spectral_upper"))
    input_receipt = recoil.HERE / "results.json"
    evidence = recoil.run()
    assert json.loads(json.dumps(evidence, default=str)) == json.loads(
        input_receipt.read_text()
    )
    rows = [
        {
            "path": row["path"],
            "native_moving_and_holder_hypercharge": row[
                "native_moving_and_holder_hypercharge"
            ],
            **common_support(row, dimension),
            **finite_bit(row, radius),
        }
        for row in evidence["clocked_rows"]
    ]
    for a, b in zip(rows[:2], rows[2:], strict=True):
        assert {
            k: v for k, v in a.items() if k != "native_moving_and_holder_hypercharge"
        } == {k: v for k, v in b.items() if k != "native_moving_and_holder_hypercharge"}
    return {
        "source_revision": recoil.SOURCE_REVISION,
        "original_source": evidence["original_source"],
        "source_sha256": {
            name: sha256(raw).hexdigest() for name, raw in sources.items()
        },
        "endpoint_recoil_receipt_sha256": sha256(
            input_receipt.read_bytes()
        ).hexdigest(),
        "distinct_active_energy_labels_upper": dimension,
        "minimum_required_relation": "full H, sender-to-receiver offdiagonal −1/2",
        "stationary_idle_spin_reference": True,
        "complete_ghost_links_reduced_by_Gauss": True,
        "one_original_stock_acceptance_required": evidence[
            "starting_stock_original_admission_lower"
        ],
        "rows_conditioned_on_that_stock": rows,
        "exact_charge_conjugate_comparison": [
            check_conjugate_path(path, transform, data) for path in recoil.PATHS
        ],
        "ideal_sharp_energy_labels_compiled": False,
        "distinct_exact_energy_gap_lower_computed": False,
        "mass_shell_certified": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = HERE / "energy_bridge.json"
    actual = json.dumps(run(), default=str, sort_keys=True, indent=2) + "\n"
    if args.write:
        receipt.write_text(actual)
    else:
        assert receipt.read_text() == actual
        print(
            "Verified shared Wilson energy support and finite spin/charge receiver bounds."
        )


if __name__ == "__main__":
    main()
