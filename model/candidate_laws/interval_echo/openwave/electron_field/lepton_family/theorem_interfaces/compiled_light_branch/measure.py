"""Exact historical rest bank and matched weak clicks at a held axial mode."""

import argparse
import hashlib
import json
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import isqrt
from pathlib import Path

import sympy as sp

from ..charged_endpoints.measure import SOURCE_REVISION, read_archive
from ..endpoint_weak_partner import representation as rep

JOINT = "model/particle_joint_mass_identity/results.json"
COMPILER = "model/lepton_mass_compiler/results.json"
SOURCES = (
    JOINT,
    COMPILER,
    "model/particle_joint_mass_identity/joint.py",
    "model/lepton_mass_compiler/compiler.py",
    "model/lepton_mass_compiler/PROOF.md",
    "model/coherent_defect_wave/mass_instrument/model.py",
    "model/common_source/chiral_candidate/kinetic.py",
    "Lean/Screen/ChargedCurrentDictionary.lean",
)
X, Z = sp.Matrix([[0, 1], [1, 0]]), sp.diag(1, -1)
I2, I4, I16 = (sp.eye(n) for n in (2, 4, 16))
BETA = sp.kronecker_product(X, I2)
CHIRALITY = sp.kronecker_product(Z, I2)
SPIN = sp.kronecker_product(I2, Z)
READOUTS = {"temporal": I2, "transverse_x": -X}


@lru_cache(maxsize=None)
def annihilator(site, count):
    op = sp.zeros(1 << count)
    bit = 1 << site
    for word in range(1 << count):
        if word & bit:
            op[word ^ bit, word] = (-1) ** ((word & (bit - 1)).bit_count())
    return op


def second_quantize(matrix):
    modes = [annihilator(j, 4) for j in range(4)]
    return sum(
        (
            matrix[i, j] * a.H * b
            for i, a in enumerate(modes)
            for j, b in enumerate(modes)
            if matrix[i, j]
        ),
        sp.zeros(16),
    )


def sea_density(one_particle):
    pairs = list(combinations(range(4), 2))
    sea = sp.zeros(16)
    for a in pairs:
        for b in pairs:
            left, right = (sum(1 << j for j in pair) for pair in (a, b))
            sea[left, right] = one_particle.extract(a, b).det()
    return sea


def reservoir_swap():
    swap = sp.zeros(32)
    order = (4, 1, 2, 3, 0)
    for word in range(32):
        occupied = [order[j] for j in range(5) if word & (1 << j)]
        inversions = sum(a > b for a, b in combinations(occupied, 2))
        swap[sum(1 << j for j in occupied), word] = (-1) ** inversions
    assert swap.H * swap == sp.eye(32)
    for site in range(5):
        assert swap.H * annihilator(site, 5) * swap == annihilator(order[site], 5)
    return swap


@lru_cache(maxsize=None)
def weak_kraus(charge, name, target_neutral_spin, wrong_chirality=False):
    """One calibrated neutral-polarization transfer, before final sea check."""
    row = READOUTS[name][target_neutral_spin, :]
    offset = 2 * wrong_chirality
    neutral = annihilator(4 + target_neutral_spin, 6)
    charged = sum((row[j] * annihilator(j + offset, 6) for j in range(2)), sp.zeros(64))
    swap_leg = neutral.H * charged
    # Rest neutral word is 01 (up occupied). Output is 11 or 00.
    before = 16
    after = 48 if charge == -1 else 0
    if charge == 1:
        swap_leg = swap_leg.H
    return swap_leg.extract(range(after, after + 16), range(before, before + 16))


def rest_bank(archive):
    h = second_quantize(BETA) + 2 * I16
    q = 2 * I16 - second_quantize(I4)
    spin = second_quantize(SPIN)
    sea = sea_density((I4 - BETA) / 2)
    assert sea * sea == sea and sp.trace(sea) == 1 and h * sea == sp.zeros(16)
    assert h * q == q * h and h * spin == spin * h
    effects = []
    for energy in range(5):
        projected = I16
        for other in set(range(5)) - {energy}:
            projected = projected * (h - other * I16) / (energy - other)
        for charge in range(-2, 3):
            for twice_spin in range(-2, 3):
                matches = sp.diag(
                    *[
                        int(q[w, w] == charge and spin[w, w] == twice_spin)
                        for w in range(16)
                    ]
                )
                effect = projected * matches
                if effect != sp.zeros(16):
                    assert effect.H == effect and effect * effect == effect
                    effects.append((energy, charge, twice_spin, effect))
    assert sum((e for *_, e in effects), sp.zeros(16)) == I16
    assert [int(sp.trace(e)) for *_, e in effects] == archive["joint_projector_ranks"]
    swap = reservoir_swap()
    rows = []
    for initial in (1, 0):
        full_input = sp.zeros(32)
        full_input[
            16 * initial : 16 * (initial + 1), 16 * initial : 16 * (initial + 1)
        ] = sea
        output = swap * full_input * swap.H
        for flag in range(2):
            source = output[16 * flag : 16 * (flag + 1), 16 * flag : 16 * (flag + 1)]
            for energy, charge, twice_spin, effect in effects:
                retained = effect * source * effect / 2
                raw = sp.trace(retained)
                row = {
                    "setting": "add" if initial else "remove",
                    "reservoir_out": flag,
                    "native_electric_charge": charge,
                    "spin_z_twice": twice_spin,
                    "source_energy_in_same_m": energy,
                    "probability_exact": str(raw),
                }
                trials = {}
                for name in READOUTS:
                    if energy != 1 or abs(charge) != 1:
                        trials[name] = {"weak_not_admitted": str(raw)}
                        continue
                    outcomes = []
                    for wrong in (False, True):
                        transfer = weak_kraus(charge, name, int(charge == -1), wrong)
                        returned = sea * transfer * effect
                        effect_click = returned.H * returned
                        assert effect * effect_click == effect_click
                        assert effect_click * q == charge * effect_click
                        click = sp.trace(effect_click * retained)
                        assert 0 <= click <= raw
                        outcomes.append(
                            {
                                "transfer_and_charged_sea_restored": str(click),
                                "complement": str(raw - click),
                            }
                        )
                    assert outcomes[0] == outcomes[1]
                    trials[name] = dict(
                        zip(("left", "wrong_right"), outcomes, strict=True)
                    )
                rows.append({"old_record": row, "weak_outcomes_per_bank": trials})
    assert [r["old_record"] for r in rows] == archive["complete_records_per_family"]
    assert sum(Fraction(r["old_record"]["probability_exact"]) for r in rows) == 1
    admitted = [
        r
        for r in rows
        if r["old_record"]["source_energy_in_same_m"] == 1
        and abs(r["old_record"]["native_electric_charge"]) == 1
    ]
    assert sum(
        Fraction(r["old_record"]["probability_exact"]) for r in admitted
    ) == Fraction(1, 2)
    for charge in (-1, 1):
        cohort = [
            r for r in admitted if r["old_record"]["native_electric_charge"] == charge
        ]
        totals = {
            name: sum(
                Fraction(
                    r["weak_outcomes_per_bank"][name]["left"][
                        "transfer_and_charged_sea_restored"
                    ]
                )
                for r in cohort
            )
            for name in READOUTS
        }
        matched = "transverse_x" if charge == -1 else "temporal"
        assert totals == {
            name: Fraction(1, 8) if name == matched else 0 for name in READOUTS
        }
    return rows


def generator_checks():
    raise_w = rep.exterior(rep.unit(3, 4))
    sites = {m: sp.eye(15)[:, rep.MASKS.index(m)] for m in (15, 23, 24)}
    assert raise_w * sites[23] == sites[15]
    assert raise_w.H * sites[15] == sites[23]
    assert raise_w * sites[24] == sp.zeros(15, 1)
    assert {m: rep.PHOTON[m] for m in sites} == {15: 0, 23: -6, 24: 6}
    assert (BETA * CHIRALITY + CHIRALITY * BETA) == sp.zeros(4)
    assert BETA * SPIN == SPIN * BETA
    admitted_control = rep.exterior(rep.unit(3, 4))
    same_family, wrong_family = sp.eye(3), sp.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    for initial in range(3):
        charged = sp.kronecker_product(sites[23], same_family[:, initial])
        for final in range(3):
            neutral = sp.kronecker_product(sites[15], same_family[:, final])
            assert (
                neutral.H
                * sp.kronecker_product(admitted_control, same_family)
                * charged
            )[0] == int(initial == final)
            assert (
                neutral.H
                * sp.kronecker_product(admitted_control, wrong_family)
                * charged
            )[0] == wrong_family[final, initial]
    # A single calibrated exchange pulse has unit transfer between the two
    # singly occupied legs and leaves empty/doubly occupied sectors unchanged.
    charged, neutral = (annihilator(j, 2) for j in range(2))
    transfer = neutral.H * charged
    pulse_generator = transfer + transfer.H
    pulse = I4 - pulse_generator**2 - sp.I * pulse_generator
    assert pulse.H * pulse == I4 and pulse[2, 1] == -sp.I
    assert [pulse[i, i] for i in (0, 3)] == [1, 1]
    pointer = sp.eye(3)
    total_charge = sp.kronecker_product(-charged.H * charged, pointer)
    total_charge += sp.kronecker_product(I4, sp.diag(-1, 0, 1))
    for charge in (-1, 1):
        destination = pointer[:, charge + 1] * pointer[1, :]
        material = transfer if charge == -1 else transfer.H
        gauged = sp.kronecker_product(material, destination)
        gate = gauged + gauged.H
        unitary = sp.eye(12) - gate**2 - sp.I * gate
        assert total_charge * unitary == unitary * total_charge
        assert unitary.H * unitary == sp.eye(12)
        origin_word, final_word = (1, 2) if charge == -1 else (2, 1)
        assert unitary[3 * final_word + charge + 1, 3 * origin_word + 1] == -sp.I


def held_projectors():
    """Use rho=p/E and mu=m/E with rho**2+mu**2=1, before numeric evaluation."""
    rho, mu = sp.symbols("rho mu", positive=True)
    divisor = mu**2 + rho**2 - 1
    rows = []
    for charge in (-1, 1):
        # A positron at +p is the hole of a negative band at -p; flip spin too.
        momentum_sign = -charge
        gamma = -momentum_sign * rho * CHIRALITY * SPIN + mu * BETA
        assert (gamma**2 - I4).applyfunc(lambda x: sp.rem(x, divisor, mu)) == sp.zeros(
            4
        )
        for physical_spin in (-1, 1):
            dirac_spin = momentum_sign * physical_spin
            line = ((I4 + momentum_sign * gamma) / 2) * ((I4 + dirac_spin * SPIN) / 2)
            assert line.H == line and sp.trace(line) == 1
            assert (line * line - line).applyfunc(
                lambda x: sp.rem(x, divisor, mu)
            ) == sp.zeros(4)
            assert (gamma * line - momentum_sign * line).applyfunc(
                lambda x: sp.rem(x, divisor, mu)
            ) == sp.zeros(4)
            # Both receiving excitations have underlying neutral down: empty
            # at +p for a particle, occupied at -p for a hole.
            values = {}
            for name, component in READOUTS.items():
                detector = component[1, :]
                left = detector.row_join(sp.zeros(1, 2))
                right = sp.zeros(1, 2).row_join(detector)
                values[name] = [
                    sp.factor((leg * line * leg.H)[0]) for leg in (left, right)
                ]
            matching = "temporal" if dirac_spin == -1 else "transverse_x"
            other = next(name for name in READOUTS if name != matching)
            preferred = sp.Rational(1, 2) * (1 + charge * physical_spin * rho)
            assert sp.simplify(values[matching][0] - preferred) == 0
            assert sp.simplify(sum(values[matching]) - 1) == 0
            assert values[other] == [0, 0]
            rows.append(
                {
                    "charge": charge,
                    "physical_spin_z_twice": physical_spin,
                    "matched_component": matching,
                    "coefficients_in_rho_equals_p_over_energy": {
                        name: [str(x) for x in v] for name, v in values.items()
                    },
                }
            )
    return rows


def finite_start_records():
    """Propagate the switch and refusal branches on the corresponding spin block."""
    w = sp.symbols("w", positive=True)
    amplitude = sp.sqrt(w * (1 - w))
    filled_band = I2[:, 0]
    outputs = {}
    for charge in (-1, 1):
        # After a reservoir refusal the charged source has the original
        # occupancy for addition, or the original vacancy for removal.
        # Their diagonal weights in the filled-mode reference are both 1-w.
        occupation_amount = 1 - w if charge == -1 else w
        occupancy = sp.Matrix(
            [[occupation_amount, amplitude], [amplitude, 1 - occupation_amount]]
        )
        assert (occupancy**2 - occupancy).applyfunc(sp.simplify) == sp.zeros(2)
        bad = occupancy if charge == -1 else I2 - occupancy
        flip = sp.simplify((filled_band.H * (I2 - bad) * filled_band)[0])
        assert flip == w
        refused = bad * filled_band
        energy_0, energy_2e = (sp.factor(x**2) for x in refused)
        assert sp.simplify(energy_0 + energy_2e + flip - 1) == 0
        assert energy_0 == (w - 1) ** 2 and energy_2e == -w * (w - 1)
        matched_success = {}
        for label, conditioned in (("left", w), ("wrong_right", 1 - w)):
            click = sp.factor(flip * conditioned)
            complement = sp.factor(flip - click)
            assert sp.simplify(energy_0 + energy_2e + click + complement - 1) == 0
            matched_success[label] = {
                "click_per_start": str(click),
                "nonclick_after_charged_admission_per_start": str(complement),
            }
        outputs[str(charge)] = {
            "original_charge_admission_per_start": str(flip),
            "neutral_refusal_energy_0_per_start": str(energy_0),
            "neutral_refusal_energy_2E_per_start": str(energy_2e),
            "matched_component": matched_success,
            "unmatched_component_per_start": {
                "click": "0",
                "nonclick_after_charged_admission": str(w),
            },
        }
    return outputs


def alternating_atan(x, n):
    terms = [(-1) ** j * x ** (2 * j + 1) / (2 * j + 1) for j in range(n + 1)]
    return sorted([sum(terms[:-1]), sum(terms)])


def momentum_bounds(period):
    # Machin's identity plus alternating remainder, followed by sin bounds.
    a0, a1 = alternating_atan(Fraction(1, 5), 14)
    b0, b1 = alternating_atan(Fraction(1, 239), 14)
    p0, p1 = (16 * a0 - 4 * b1), (16 * a1 - 4 * b0)
    x0, x1 = 2 * p0 / period, 2 * p1 / period
    assert x1 < Fraction(1, 64) < 2 * x0
    return (x0 - x0**3 / 6 + x0**5 / 120 - x0**7 / 5040, x1 - x1**3 / 6 + x1**5 / 120)


def sqrt_bounds(x):
    denominator = 10**25
    lower = isqrt(x.numerator * denominator**2 // x.denominator)
    return Fraction(lower, denominator), Fraction(lower + 1, denominator)


def outward(bounds, precision=12):
    base = 10**precision
    low, high = (v * base for v in bounds)
    integers = (
        low.numerator // low.denominator,
        -(-high.numerator // high.denominator),
    )
    return [f"{n // base}.{n % base:0{precision}d}" for n in integers]


def held_families(archive, compiler):
    period = 512
    p0, p1 = momentum_bounds(period)
    assert compiler["fixed_scalar_background_norm_exact"] == "1"
    energies = compiler["same_effect_rest_energies_k0_k1_k2"]
    families = sorted(
        archive["families"],
        key=lambda row: Fraction(row["source_rest_energy_interval"][0]),
    )
    results = []
    preferred_intervals = []
    held_energies = []
    for family in families:
        m0, m1 = map(Fraction, family["source_rest_energy_interval"])
        k = family["character"]
        assert abs(float((m0 + m1) / 2) - energies[k]) < 1e-14
        e0, e1 = sqrt_bounds(p0**2 + m0**2)[0], sqrt_bounds(p1**2 + m1**2)[1]
        ratio0 = p0 / sqrt_bounds(p0**2 + m1**2)[1]
        ratio1 = p1 / sqrt_bounds(p1**2 + m0**2)[0]
        large = ((1 + ratio0) / 2, (1 + ratio1) / 2)
        small = (1 - large[1], 1 - large[0])
        preferred_intervals.append(large)
        held_energies.append((e0, e1))
        results.append(
            {
                "character": k,
                "rest_energy_enclosure": outward((m0, m1)),
                "held_energy_enclosure": outward((e0, e1)),
                "external_resonance_work_abs_enclosure": outward((e0 - p1, e1 - p0)),
                "negative_helicity_negative_charge_or_positive_helicity_positive_charge": outward(
                    large
                ),
                "opposite_helicity_either_charge": outward(small),
                "matched_left_click_per_original_preferred_helicity_start": outward(
                    (large[0] ** 2, large[1] ** 2)
                ),
                "matched_left_click_per_original_opposite_helicity_start": outward(
                    (small[0] ** 2, small[1] ** 2)
                ),
            }
        )
    assert all(
        a[0] > b[1] for a, b in zip(preferred_intervals, preferred_intervals[1:])
    )
    cut = Fraction("0.904651804389")
    margin = Fraction("0.09533")
    assert preferred_intervals[0][0] - cut > margin
    assert cut - max(v[1] for v in preferred_intervals[1:]) > margin
    meter_halfwidth = Fraction(1, 10**7)
    assert held_energies[0][0] - p1 > 2 * meter_halfwidth
    assert all(
        b[0] - a[1] > 2 * meter_halfwidth
        for a, b in zip(held_energies, held_energies[1:])
    )
    return {
        "period": period,
        "axial_mode": "pi/256",
        "sine_momentum_enclosure": outward((p0, p1)),
        "families": results,
        "light_cut_on_preferred_conditional_response": str(cut),
        "conditional_probability_margin_exceeds": str(margin),
        "sufficient_symmetric_energy_meter_halfwidth": str(meter_halfwidth),
    }


def finite_compiler_margin():
    # Apply the archived complete-law bound 108 V² T²/N at T=1 to the
    # declared held periodic board. An identical appended ideal instrument
    # cannot increase trace distance from the limiting outcome law.
    volume = 4 * 4 * 512
    trace_norm_bound = Fraction(1, 1000)
    repetitions = 108 * volume**2 / trace_norm_bound
    assert repetitions.denominator == 1
    event_error = trace_norm_bound / 2
    lower_admission = Fraction("0.523")
    conditioned_error = 2 * event_error / (lower_admission - event_error)
    assert conditioned_error < Fraction(1, 500)
    return {
        "conditional_on_archive_hypotheses_and_identical_appended_instruments": True,
        "illustrative_periods": [4, 4, 512],
        "total_compiler_wait": 1,
        "sufficient_round_count": int(repetitions),
        "complete_trace_norm_bound": str(trace_norm_bound),
        "per_event_absolute_error_bound": str(event_error),
        "preferred_conditional_probability_error_below": "1/500",
        "executed_finite_compilation": False,
    }


def calculate():
    if not __debug__:
        raise RuntimeError(
            "Run without Python optimization so all exact checks execute."
        )
    archive = json.loads(read_archive(JOINT))["compiled_lepton_joint"]
    compiler = json.loads(read_archive(COMPILER))["same_source_kinetic_attachment"]
    generator_checks()
    return {
        "schema": "oph.compiled_light_branch.v1",
        "historical_source_revision": SOURCE_REVISION,
        "historical_blobs_sha256": {
            path: hashlib.sha256(read_archive(path)).hexdigest() for path in SOURCES
        },
        "rest_records": rest_bank(archive),
        "held_line_identities": held_projectors(),
        "held_reservoir_and_weak_denominators_in_w": finite_start_records(),
        "held_weak_residue_intervals": held_families(archive, compiler),
        "conditional_finite_compiler_margin": finite_compiler_margin(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    target = Path(__file__).with_name("results.json")
    body = json.dumps(calculate(), indent=2, sort_keys=True) + "\n"
    if parser.parse_args().record:
        target.write_text(body)
    elif target.read_text() != body:
        raise RuntimeError("Independent replay differs from the archived receipt.")
    print(
        "Exact archived bank preserved; held energy, charge, spin and weak comparisons replayed."
    )


if __name__ == "__main__":
    main()
