"""Replay the two-path energy and fresh-receiver comparison on the Wilson cube."""

import argparse
from dataclasses import replace
from fractions import Fraction as F
from hashlib import sha256
from itertools import pairwise, product
import json
from pathlib import Path

import sympy as sp

from ..charged_endpoints.cube import (
    A,
    BETA_SEA,
    C,
    CELLS,
    CENTER,
    LINKS,
    MASKS,
    PARTNERS,
    Q6,
    Word,
    addressed_character,
    c_action,
    charge,
    gauss,
    mode,
)
from ..charged_endpoints.measure import (
    CONTRACT,
    GAMMA,
    SOURCE_REVISION,
    clifford,
    read_archive,
)

HERE = Path(__file__).resolve().parent
M15 = "model/observer_particles/results/representation_identification_operators.json"
SCALE = "model/particle_primary_spectral_scale/results.json"
HOLDER = PARTNERS[0]
PATHS = ((HOLDER, A), (HOLDER, C, CENTER, A))
RECEIVER = CENTER
ACTIVE_BOND = tuple(sorted((A, RECEIVER)))


def closed_faces():
    for a, b in ((0, 1), (0, 2), (1, 2)):
        for x in CELLS:
            if x[a] != 0 or x[b] != 0:
                continue
            u = tuple(v + int(j == a) for j, v in enumerate(x))
            v = tuple(w + int(j == b) for j, w in enumerate(x))
            y = tuple(w + int(j == b) for j, w in enumerate(u))
            yield (x, u, y, v, x)


FACES = tuple(closed_faces())


def charge_blocks(sources):
    raw = json.loads(sources[M15])
    packet = raw["even_identification_operators"]
    archived_masks = [
        sum(1 << j for j in raw["basis_subsets"][i]) for i in raw["even_indices"]
    ]
    assert sorted(archived_masks) == list(MASKS)

    def diagonal(label):
        matrix = packet[label]
        assert matrix["shape"] == [15, 15]
        assert all(i == j for i, j, _ in matrix["entries"])
        values = {i: F(k) for i, _, k in matrix["entries"]}
        return [values.get(archived_masks.index(m), F(0)) for m in MASKS]

    charges = diagonal("q")
    weights = list(zip(*(diagonal(k) for k in ("H1", "H2", "T3")), strict=True))
    assert sum(charges) == 0
    assert all(sum(values) == 0 for values in zip(*weights, strict=True))
    blocks = {
        z: tuple(i for i, q in enumerate(charges) if q == z)
        for z in sorted(set(charges))
    }
    assert blocks[6] == (Q6,)
    signatures = {
        z: tuple(sorted(weights[i] for i in ids)) for z, ids in blocks.items()
    }
    signed = (
        signatures
        | {
            -z: tuple(sorted(tuple(-k for k in weight) for weight in sig))
            for z, sig in signatures.items()
        }
        | {0: ((0, 0, 0),)}
    )
    assert len(signed) == 11
    # Haar projection retains all normalized coefficients of an isotypic
    # matrix-coefficient space, including an equivalent conjugate weak doublet.
    survives = {
        (z, e): signed.get(z + 6 * e) == sig
        for z, sig in signed.items()
        if z != 0
        for e in (-1, 0, 1)
    }
    scale = json.loads(sources[SCALE])["actual_local_expansion_channels"]
    casimirs = {
        row["integer_Q6_charge"]: F(row["trace_metric_Casimir"])
        for row in scale["all_native_blocks"]
    }
    assert set(casimirs) == set(blocks)
    return blocks, survives, casimirs


def beta_hops(gammas):
    e = sp.eye(4)
    chart = sp.Matrix.hstack(
        e[:, 0] + e[:, 2], e[:, 1] + e[:, 3], -e[:, 0] + e[:, 2], -e[:, 1] + e[:, 3]
    ) / sp.sqrt(2)
    beta = gammas[3]
    assert chart.H * beta * chart == sp.diag(1, 1, -1, -1)
    hops = [sp.simplify(chart.H * (-beta + beta * g) * chart / 2) for g in gammas[:3]]
    for hop in hops:
        gram = hop.H * hop
        assert gram * gram == gram and sp.trace(gram) == 2
    return hops


def prepare(path, signed, spin):
    """The same two local Wilson legs as the exchange, with a fixed holder spin."""
    x = path[-1]
    j = mode(HOLDER, -1, 0) if signed > 0 else mode(x, -1, spin)
    k = mode(x, 1, spin) if signed > 0 else mode(HOLDER, 1, 0)
    mask, parity1 = c_action(BETA_SEA, j, False)
    mask, parity2 = c_action(mask, k, True)
    word = Word(mask, (0,) * len(LINKS))
    for u, v in pairwise(path):
        word = addressed_character(word, u, v, signed)
    assert not any(gauss(word).values())
    assert charge(word) == {
        c: 6 * signed * (int(c == x) - int(c == HOLDER)) for c in CELLS
    }
    return word, parity1 * parity2


def destination(word, signed, spin):
    """Apply the signed native one-link isometry and keep its CAR phase."""
    tmp, a = c_action(word.occupied, mode(A, signed, spin), signed < 0)
    mask, b = c_action(tmp, mode(RECEIVER, signed, spin), signed > 0)
    final = addressed_character(replace(word, occupied=mask), A, RECEIVER, signed)
    assert not any(gauss(final).values())
    assert charge(final) == {
        c: 6 * signed * (int(c == RECEIVER) - int(c == HOLDER)) for c in CELLS
    }
    return final, a * b


def hopping_norms(word, blocks, stays, hops, cut=False):
    norms = {}
    for (x, y), exponent in zip(LINKS, word.characters, strict=True):
        e = (x, y)
        if cut and e == ACTIVE_BOND:
            continue
        a = next(j for j in range(3) if x[j] != y[j])
        w = hops[a]
        value = F(0)
        for target, origin, m, orientation in ((x, y, w, 1), (y, x, w.H, -1)):
            x0, y0 = (60 * CELLS.index(site) for site in (target, origin))
            for z, idx in blocks.items():
                if not stays[orientation * z, exponent]:
                    continue
                for r, s in product(range(4), repeat=2):
                    if m[r, s] == 0:
                        continue
                    unit = F(sp.simplify(sp.Abs(m[r, s]) ** 2)) / len(idx)
                    for i, j in product(idx, repeat=2):
                        arrival = x0 + 15 * r + i
                        departure = y0 + 15 * s + j
                        if (word.occupied >> departure) & 1 and not (
                            word.occupied >> arrival
                        ) & 1:
                            value += unit
        norms[e] = value
    return norms


def magnetic_norm(word, blocks, stays, coefficient):
    face_live = []
    for path in FACES:
        live = []
        for z in blocks:
            for orientation in (-1, 1):
                if all(
                    stays[
                        orientation * z * (1 if u < v else -1),
                        word.characters[LINKS.index(tuple(sorted((u, v))))],
                    ]
                    for u, v in pairwise(path)
                ):
                    live.append(orientation * z)
        face_live.append(live)
    return sum(map(len, face_live)) * (coefficient / 30) ** 2, face_live


def diagnostic(word, input_data, alpha, cut=False):
    blocks, stays, casimirs, hops, electric, magnetic = input_data
    defects = (BETA_SEA ^ word.occupied).bit_count()
    links = sum(bool(e) for e in word.characters)
    electric_mean = electric / 2 * casimirs[6] * links
    mean = alpha * (-240 + defects) + electric_mean
    link = hopping_norms(word, blocks, stays, hops, cut=cut)
    mag, live = magnetic_norm(word, blocks, stays, magnetic)
    return {
        "mean": mean,
        "electric_mean": electric_mean,
        "variance": sum(link.values()) + mag,
        "hopping_variance_by_bond": list(link.values()),
        "magnetic_variance": mag,
        "surviving_magnetic_character_labels_by_face": live,
    }


def ramsey(source, stock, wait):
    """Rational finite-time bands from the complete unshuttered clock arms."""
    gap = source["mean"] - stock["mean"]
    x = gap * wait
    assert 0 < x <= 1
    sin_lo = x - x**3 / 6
    sin_hi = sin_lo + x**5 / 120
    allowance = wait**2 * (source["variance"] + stock["variance"]) / 2
    # C = <source|U|source><stock|U^*|stock>, P(Y+)=(1+Im C)/2.
    return [(1 - sin_hi - allowance) / 2, (1 - sin_lo + allowance) / 2]


def probability_band(magnitude_lo, magnitude_hi, error):
    return [max(F(0), magnitude_lo - error) ** 2, (magnitude_hi + error) ** 2]


def bell_diagnostics(rows):
    a, b = rows
    for key in (
        "mean",
        "electric_mean",
        "magnetic_variance",
        "surviving_magnetic_character_labels_by_face",
    ):
        assert a[key] == b[key]
    return {
        **a,
        "variance": (a["variance"] + b["variance"]) / 2,
        "unrecombined_sharp_spin_variances": [a["variance"], b["variance"]],
        "hopping_variance_by_bond": [
            (v + w) / 2
            for v, w in zip(
                a["hopping_variance_by_bond"],
                b["hopping_variance_by_bond"],
                strict=True,
            )
        ],
    }


def transfer(source, sink, alpha, string_unit, electric_max, magnetic, cut=False):
    """Finite holding test using a fresh addressed and Bell-checked receiver."""
    delta = sink["mean"] - source["mean"]
    assert abs(delta) == string_unit and source["variance"] < 14**2
    assert sink["variance"] < 14**2
    # One application of H reaches at most four changed links and two further
    # beta defects. The cut-law bound is dominated by the uncut norm bound.
    offdiagonal_bound = 30 * len(LINKS) + magnetic * len(FACES)
    support_bound = abs(delta) + 2 * alpha + 4 * electric_max + offdiagonal_bound
    sink_centered_bound = 22
    assert delta**2 + sink["variance"] < sink_centered_bound**2
    remainder = 14 * (abs(delta) * sink_centered_bound + 14 * support_bound) / 6
    wait = F(1, 1024)
    error = remainder * wait**3
    prob = (
        probability_band(F(0), F(0), error)
        if cut
        else probability_band(wait / 2, wait / 2 + abs(delta) * wait**2 / 4, error)
    )
    phase_wait = F(1, 16384)
    # Use the genuine signed hop reference with one common calibrated phase
    # for both members of a same-centroid pair.
    phase_error = remainder * phase_wait**3
    amplitude_bound = phase_error + phase_wait / 2 + abs(delta) * phase_wait**2 / 4
    # The separately calibrated absolute Ramsey centroid can err by at most
    # one quarter source energy unit without making the held signs overlap.
    center_error = F(1, 4)
    phase_error += center_error * phase_wait * amplitude_bound
    real = F(0) if cut else delta * phase_wait**2 / 4
    return {
        "leading_joint_receiver_probability_coefficient": "0" if cut else "1/4",
        "arrival_wait": wait,
        "same_address_charge_gauge_and_Bell_probability": prob,
        "sender_to_receiver_diagonal_detuning": delta,
        "held_phase_wait": phase_wait,
        "centered_reference_transition_real": [real - phase_error, real + phase_error],
        "allowed_center_calibration_error": center_error,
        "third_order_amplitude_remainder_coefficient_upper": remainder,
    }


def experiment(signed, path, alpha, data, wait):
    _, _, casimirs, hops, electric, magnetic = data
    blank = Word(BETA_SEA, (0,) * len(LINKS))
    stock = diagnostic(blank, data, alpha)
    source_rows = []
    sink_rows = []
    cut_rows = []
    cut_sink_rows = []
    for spin in (0, 1):
        word, initial_phase = prepare(path, signed, spin)
        final, raw_hop_phase = destination(word, signed, spin)
        # The forward and backward diagonal beta blocks are scalar in spin.
        # Their particle and hole CAR transfers need opposite direction/orientation.
        axis = next(j for j in range(3) if A[j] != RECEIVER[j])
        w = hops[axis]
        block = w[:2, :2] if signed > 0 else -w[2:, 2:].T
        assert block == -sp.eye(2) / 2
        assert initial_phase in (-1, 1) and raw_hop_phase in (-1, 1)
        source_rows.append(diagnostic(word, data, alpha))
        sink_rows.append(diagnostic(final, data, alpha))
        cut_rows.append(diagnostic(word, data, alpha, cut=True))
        cut_sink_rows.append(diagnostic(final, data, alpha, cut=True))
    source, sink = bell_diagnostics(source_rows), bell_diagnostics(sink_rows)
    string_unit = electric / 2 * casimirs[6]
    electric_max = electric / 2 * max(casimirs.values())
    # On separate copies, use absolute Y Ramsey at the same short increment.
    # |t*mean| <= 1/2, so the in-domain inverse sine has slope < 8/7.
    dt = F(1, 16384)
    y_error = F(1, 2**18)
    assert abs(dt * source["mean"]) <= F(1, 2)
    sine_error = 2 * y_error + dt**2 * source["variance"] / 2
    half = F(1, 2)
    assert (half - half**3 / 6 + half**5 / 120 + sine_error) ** 2 < F(15, 64)
    center_error = F(8, 7) * sine_error / dt
    assert center_error < F(1, 4)
    return {
        "path": path,
        "native_moving_and_holder_hypercharge": [6 * signed, -6 * signed],
        "onsite_beta_coefficient": alpha,
        "absolute_source_clock_diagnostics": source,
        "absolute_receiver_clock_diagnostics": sink,
        "matched_blank_expectation": stock["mean"],
        "matched_transfer_centroid": source["mean"] - stock["mean"],
        "matched_transfer_variance": source["variance"] + stock["variance"],
        "relative_centroid_with_electric_mean_subtracted": (
            source["mean"] - stock["mean"] - source["electric_mean"]
        ),
        "unshuttered_two_copy_Ramsey_Y_plus_band": ramsey(source, stock, wait),
        "absolute_Ramsey_center_calibration": {
            "wait": dt,
            "required_Y_plus_probability_resolution": y_error,
            "derived_energy_error_upper": center_error,
        },
        "fresh_receiver": transfer(
            source, sink, alpha, string_unit, electric_max, magnetic
        ),
        "cut_bond_sham": transfer(
            bell_diagnostics(cut_rows),
            bell_diagnostics(cut_sink_rows),
            alpha,
            string_unit,
            electric_max,
            magnetic,
            cut=True,
        ),
    }


def run():
    sources = {name: read_archive(name) for name in (CONTRACT, GAMMA, M15, SCALE)}
    member = json.loads(sources[CONTRACT])["fixed_primary_member"]
    source = member["source"]
    for name in (GAMMA, M15):
        assert (
            sha256(sources[name]).hexdigest() == source["implementation_sha256"][name]
        )
    assert source["side"] == 2 and source["per_link_word_degree"] == 1
    assert F(source["spacing"]) == 1 and source["interior_plaquettes"] == len(FACES)
    electric = F(source["electric_coefficient"])
    magnetic = F(source["magnetic_coefficient"])
    blocks, stays, casimirs = charge_blocks(sources)
    assert electric == 105 and magnetic == F(1, 7) and casimirs[6] == F(3, 10)
    data = (
        blocks,
        stays,
        casimirs,
        beta_hops(clifford(sources[GAMMA])),
        electric,
        magnetic,
    )
    blank = diagnostic(Word(BETA_SEA, (0,) * 12), data, F(3))
    assert blank["mean"] == -720 and blank["variance"] == 180 + F(1, 735)
    wait = F(1, 128)
    rows = [
        experiment(sign, path, F(3), data, wait) for sign in (1, -1) for path in PATHS
    ]
    tuned = [
        experiment(sign, PATHS[0], 3 + electric / 2 * casimirs[6], data, wait)
        for sign in (1, -1)
    ]
    for i, control in enumerate(tuned):
        short, detour = rows[2 * i : 2 * i + 2]
        assert short["matched_transfer_centroid"] == F(87, 4)
        assert detour["matched_transfer_centroid"] == F(213, 4)
        assert (
            control["matched_transfer_centroid"] == detour["matched_transfer_centroid"]
        )
        assert (
            detour["unshuttered_two_copy_Ramsey_Y_plus_band"][1]
            < short["unshuttered_two_copy_Ramsey_Y_plus_band"][0]
        )
        for row in (short, detour, control):
            live = row["fresh_receiver"]
            assert (
                live["same_address_charge_gauge_and_Bell_probability"][0]
                > (
                    row["cut_bond_sham"][
                        "same_address_charge_gauge_and_Bell_probability"
                    ][1]
                )
            )
        assert (
            detour["fresh_receiver"]["centered_reference_transition_real"][1]
            < 0
            < control["fresh_receiver"]["centered_reference_transition_real"][0]
        )
    for a, b in zip(rows[:2], rows[2:], strict=True):
        assert all(
            a[k] == b[k]
            for k in (
                "matched_transfer_centroid",
                "relative_centroid_with_electric_mean_subtracted",
                "matched_transfer_variance",
                "fresh_receiver",
                "cut_bond_sham",
            )
        )
    state = member["state"]
    width = F(state["spectral_upper"]) - F(state["spectral_lower"])
    assert F(state["beta"]) * width <= state["dyadic_scale"]
    rho_floor = 480 + 9 * source["owned_links"] + 2 * state["dyadic_scale"]
    return {
        "original_source": {"member_id": member["member_id"], **member["ids"]},
        "source_revision": SOURCE_REVISION,
        "source_sha256": {name: sha256(b).hexdigest() for name, b in sources.items()},
        "starting_stock_original_admission_lower": {
            "base": 2,
            "negative_exponent": rho_floor,
        },
        "clocked_stock_requires_two_independent_acceptances": True,
        "idle_original_references_kept": True,
        "site_charge_units_not_calibrated_to_SI": True,
        "native_terms_removed": [],
        "blank_full_source_clock_diagnostics": blank,
        "path_clock_wait": wait,
        "clocked_rows": rows,
        "counterfactual_onsite_tuned_to_match_detour_centroid": tuned,
        "signed_sites": {
            "sending": A,
            "receiving": RECEIVER,
            "conjugate_holder": HOLDER,
        },
        "oriented_bonds": LINKS,
        "oriented_faces": FACES,
        "physical_refinement_and_endpoint_mass_identified": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = HERE / "results.json"
    actual = json.dumps(run(), default=str, sort_keys=True, indent=2) + "\n"
    if args.write:
        receipt.write_text(actual)
    else:
        assert receipt.read_text() == actual
        print(
            "Verified full Wilson cube moments and finite charge/spin receiver bands."
        )


if __name__ == "__main__":
    main()
