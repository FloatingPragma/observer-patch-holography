"""Decode the row 9 scientific gate from its local research receipts."""

from hashlib import sha256
from itertools import product
import json
from math import fsum, isclose, sqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[5]
EVIDENCE = ROOT / "model/candidate_laws/interval_echo/openwave/row_09_neutral"
LEVELS = (0, 8, 9, 32, 33)


def read(name):
    return json.loads((EVIDENCE / name).read_bytes())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def valid_bank(bank, original, size):
    return (
        len(bank) == size
        and all(type(x) is int and x >= 0 for x in bank)
        and sum(bank) == original
    )


def checksum(receipt, path, field):
    require(
        sha256((EVIDENCE / path).read_bytes()).hexdigest() == receipt[field],
        f"Input digest differs for {path}",
    )


def decode_cycle(panel, original):
    references = panel["daughter_minus_cold_contact_records"]
    duty = panel["relative_link_duty"]
    n = panel["dimers"]
    require(len(references) == 2 * n - 1, "Incomplete layer contacts")
    terms = []
    for j, contact in enumerate(references):
        differences = []
        for bank in contact["daughter_and_cold_bins"]:
            require(valid_bank(bank, original, 7), "Incomplete original contact bank")
            accepted = sum(bank[2:])
            require(accepted > 0, "No valid outcomes")
            differences.append(
                sum(x * y for x, y in zip(LEVELS, bank[2:], strict=True)) / accepted
            )
        w = (1.0 if not j % 2 else duty / 2 if j == 1 else duty) / (1 + duty)
        terms.append(w * (differences[0] - differences[1]))
    cycle = fsum(terms)
    stated = panel["observed_energy_and_switch_work"]["cycle_average"]
    require(
        isclose(cycle, stated["sampled"], abs_tol=1e-13, rel_tol=0),
        "Cycle average does not decode",
    )
    require(
        isclose(
            stated["calculated"],
            panel["slope_of_lowest_relative_phase_per_original_period"],
            abs_tol=2e-11,
            rel_tol=0,
        ),
        "Contact value and phase response disagree",
    )
    lo, hi = stated["simultaneous_99pct_interval"]
    require(lo <= cycle <= hi, "Cycle interval excludes its own centre")
    c0, c1 = panel["light_cutoff_from_separate_clock"]
    require(
        panel["source_link_substitution_full_refusal"] == [original, 0, 0, 0, 0, 0, 0],
        "Cut parent has nonzero admissions",
    )
    if panel["schedule"] == "half":
        require(lo > c1, "Half-duty control does not reject the light interpretation")
        return f"{n} half duty: ordinary cycle average {cycle:.8f}"
    require(0 < lo <= hi < c0, "Primary ordinary energy isn't light")
    joining = panel["inherited_same_mode_emission"]
    lower = joining["lower_99pct_from_original_joint_banks"]
    bright = joining["bright_only_probability_upper"]
    require(
        panel["lowest_relative_quasifrequency"]
        < c0
        <= c1
        < joining["first_bright_relative_frequency"],
        "Unresolved low-mode identity",
    )
    require(
        0 < bright < lower,
        "Parent/absorption/vacuum record cannot identify the low mode",
    )
    source_duties = joining["registered_source_duty_and_joint_lower"]
    require(
        source_duties[0] == [0.0, 0.0]
        and [x[0] for x in source_duties] == [0.0, 0.25, 0.5, 1.0]
        and all(x[1] > bright for x in source_duties[1:]),
        "Source coupling controls do not separate",
    )
    mass = (sqrt(lower) - sqrt(bright)) ** 2
    require(
        isclose(
            mass,
            joining["lowest_mode_mass_lower_from_coherent_amplitude_triangle"],
            abs_tol=1e-14,
            rel_tol=0,
        ),
        "Incorrect coherent leakage subtraction",
    )
    raw = panel["guard_only_terminal_charge_control"]["guarded_and_open_guard_counts"]
    require(
        all(valid_bank(bank, original, 5) for bank in raw),
        "Incomplete charge-control bank",
    )
    guard, open_guard = raw
    require(
        not (guard[2] or guard[4]) and open_guard[2] > 0 and open_guard[4] > 0,
        "Same-input charge guard control does not separate",
    )
    return f"{n} native: average {cycle:.8f}; tagged mode ≥{mass:.8g} per original; guard-only charged bins {open_guard[2]}/{open_guard[4]}"


def bare_floor(data):
    states = {
        (p, m): v for p, m, v in data["local_counts_plus_minus_and_twice_shifted_floor"]
    }
    require(
        data["words_checked"] == 3**12 and len(states) == 91,
        "Incomplete exact native cell census",
    )
    require(states.pop((0, 0)) == 0, "Filled cold reference differs")
    require(min(states.values()) == 7, "Local contrast-shift floor differs")
    single = [(p, 0, e) for (p, m), e in states.items() if p == m]
    double = [
        (p + q, 1, states[p, m] + states[q, r])
        for (p, m), (q, r) in product(states, repeat=2)
        if p + q == m + r
    ]
    require(
        min(e for _, _, e in single + double) == 18
        and all(p == 1 for p, _, e in single + double if e == 18),
        "The exact native regional floor or equality case differs",
    )
    require(
        min(22, min(e for p, _, e in single + double if p >= 2)) == 22,
        "The multiple-pair floor differs",
    )
    require(
        data["all_nonvacuum_neutral_formation_bounds"]
        == {
            "boundary_shut_floor_contrast_units": 9,
            "intact_board_floor_contrast_units": 10,
            "two_or_more_pairs_floor_contrast_units": 11,
        },
        "Unexpected native formation findings",
    )
    reads = data["finite_512_unit_period_32768_label_joint_charge_and_energy_read"]
    native, changed = [
        reads[k]["neutral_nontrivial"]["computed_bright_per_uniform_sharp_enrollment"]
        for k in ("8", "0")
    ]
    require(
        0 <= native < 1e-5 and changed > 0.1,
        "The two same-reader light controls do not separate",
    )
    print(
        f"Bare neutral floors: intact 80, owner-shuttered 72, two or more defect pairs 88; bright {native:.8g} vs {changed:.8g}"
    )


def main():
    energy = read("critical_emission/cycle_energy/results.json")
    checksum(energy, "critical_emission/emitter.json", "source_emission_sha256")
    checksum(energy, "critical_emission/cycle_energy/PROTOCOL.md", "protocol_sha256")
    bare = read("bare_neutral_bound/results.json")
    checksum(bare, "bare_neutral_bound/PROTOCOL.md", "protocol_sha256")
    checksum(bare, "passive_decay/records.json", "independent_parent_receipt_sha256")
    require(
        {(x["dimers"], x["schedule"]) for x in energy["panels"]}
        == {(23, "native"), (23, "half"), (45, "native"), (45, "half")},
        "Incomplete pulse controls",
    )
    for panel in energy["panels"]:
        print(decode_cycle(panel, energy["shots_per_contact_and_branch"]))
    bare_floor(bare)
    print(
        "HOLD for full row: native powered release passes the light test; the exact neutral outcome depends on the amount guard, and the candidate differs from the pin."
    )


if __name__ == "__main__":
    main()
