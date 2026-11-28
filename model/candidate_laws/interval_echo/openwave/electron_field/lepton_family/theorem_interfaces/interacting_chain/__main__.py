"""Replay the contact-ring instruments and compare raw reference/receiver records."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial
from pathlib import Path

from ....replay import compare
from . import measure, ring
from .independent import measure as independent


HERE = Path(__file__).resolve().parent
REFUSALS_AND_REPORTS = [
    "ground_refused",
    "pump_refused",
    "energy_window_refused",
    "two_site_charge_and_single_gate_refused",
    "++",
    "+-",
    "-+",
    "--",
]
AXES = ["XX", "YY", "ZZ"]
SOURCE_SIGNS = [1, -1]
RECEIVER = "simultaneous_signed_charge_and_neighbor_spin"


def require(condition, description):
    if not condition:
        raise ValueError(description)


def protocol_digest(folder):
    return sha256((folder / "PROTOCOL.md").read_bytes()).hexdigest()


def contains(interval, value):
    return interval[0] <= value <= interval[1]


def raw_witness(report, starts, sampling_radius):
    """Evaluate X−Y+Z−p on all starts, with the joint event fixed before axes."""
    pauli = report["simulated_local_Pauli_comparison"]["banks"]
    require([bank["axis_pair"] for bank in pauli] == AXES, "Pauli settings differ")
    p = list(map(Fraction, report["joint_gate_probability_per_original_Neel_start"]))
    eta = list(
        map(
            Fraction,
            report["reference_channel_without_global_terminal_drain"][
                "single_qubit_depolarizing_eta"
            ],
        )
    )
    p_center = sum(p) / 2
    eta_center = sum(eta) / 2
    # The simulated bank uses interval midpoints. These allowances also cover
    # the exact-model spectral intervals and the midpoint's floating arithmetic.
    roundoff = Fraction(1, 10**12)
    gate_error = (p[1] - p[0]) / 2 + roundoff
    signal_error = (
        max(abs(q * s - p_center * eta_center) for q in p for s in eta) + roundoff
    )
    require(
        0 < p[0] <= p[1] <= 1 and -1 <= eta[0] <= eta[1] <= 1, "Nonphysical interval"
    )
    signed = []
    admitted_upper = [p[1]]
    for record in pauli:
        counts = record["counts"]
        require(
            len(counts) == len(REFUSALS_AND_REPORTS)
            and all(type(v) is int and v >= 0 for v in counts)
            and sum(counts) == starts,
            "Pauli bank omits or corrupts original-start outcomes",
        )
        pp, pm, mp, mm = counts[4:]
        signed.append(Fraction(pp - pm - mp + mm, starts))
        admitted_upper.append(
            Fraction(pp + pm + mp + mm, starts) + sampling_radius + gate_error
        )
    lower = signed[0] - signed[1] + signed[2] - 3 * (sampling_radius + signal_error)
    return lower - min(admitted_upper)


def joined_observer_records(main_record, independent_record):
    """Join only the original eight-site source and its same-event receiver."""
    require(
        main_record["source_revision"] == measure.SOURCE_REVISION
        and main_record["archived_source_sha256"]
        == independent_record["archived_source_sha256"]
        == ring.ANCESTRY,
        "Contact-chain source ancestry differs",
    )
    require(
        main_record["protocol_sha256"] == protocol_digest(HERE)
        and independent_record["protocol_sha256"]
        == protocol_digest(HERE / "independent"),
        "The corresponding preparation and receiver protocols must be retained",
    )
    first, eight = main_record["ring_measurements"]
    require(
        (first["length"], eight["length"], independent_record["length"]) == (4, 8, 8),
        "The simultaneous charge-and-spin receiver belongs to the eight-site source",
    )
    local = independent_record["signed_branches"]
    require(
        [bank["sign"] for bank in local] == SOURCE_SIGNS,
        "Signed trials are missing or reordered",
    )
    require(
        independent_record["simulated_count_bin_order"] == REFUSALS_AND_REPORTS,
        "Refusal alphabet differs",
    )
    starts = independent_record["simulated_original_Starts_per_Pauli_bank"]
    cohorts = independent_record["simultaneous_local_cohort_count"]
    require(
        starts == independent.TRIALS and cohorts == independent.COMPARISONS,
        "Count resolution differs",
    )
    radius = Fraction(1, 7168)
    exponent = starts * radius**2 / 2
    # Two statistics on three axes per cohort; the length-2 Hoeffding bound
    # is conservative for both the signed report and the admission indicator.
    exp_lower = sum(exponent**k / factorial(k) for k in range(31))
    family_error = Fraction(12 * cohorts, 1) / exp_lower
    require(
        family_error < Fraction(1, 100), "Simultaneous raw count allowance is too small"
    )
    margins = {"half_beat": [], "fixed_wait_149": []}
    checked = 0
    for baseline, bank in zip(eight["sign_banks"], local, strict=True):
        sign = bank["sign"]
        require(baseline["sign"] == sign, "The signed source banks do not match")
        window = baseline["two_energy_window"]
        lines = bank["complete_dark_and_bright_lines"]
        require(
            len(lines) == 4
            and len(window["levels_inside_gate"]) == 4
            and bank["both_bright_translation_parities"] == [1, -1]
            and bank["two_spin_source_and_receiver_intertwiner_checked"]
            and bank[
                "half_beat_T4_on_reachable_support_exact_by_line_count_and_translation"
            ]
            and bank["fixed_wait_direct_propagation_checked"],
            "Same-source two-energy and spin premises differ",
        )
        for row, original in zip(lines, window["levels_inside_gate"], strict=True):
            require(
                contains(row["reservoir_corrected_energy"], original["corrected"])
                and len(row["momenta"]) == original["spin_z_block_rank"]
                and row["exact_spin_selection_excludes_source"]
                == (row["total_spin"] == "3/2"),
                "An energy component has changed, including the dark spectrum",
            )
        require(
            contains(bank["half_beat_native_time"], window["half_beat_wait"])
            and contains(
                bank["joint_window_acceptance_per_Neel_start"],
                window["window_pass_per_original"],
            ),
            "Wait or original preparation denominator differs",
        )
        settings = bank["free_readouts"]
        require(
            set(settings) == {"half_beat", "fixed_wait_149", "parked"},
            "Free-wait controls differ",
        )
        dephased = bank["resolved_energy_then_forgotten_receiver_local"]
        for wait, record in (("resolved", dephased), *settings.items()):
            reports = record[RECEIVER]
            require(
                [(r["charge_site"], r["spin_site"]) for r in reports]
                == [(4, 3), (4, 5)],
                "The charge site and spin neighbor must belong to the same terminal gate",
            )
            for report in reports:
                require(
                    report["joint_native_excess_on_the_two_receiving_sites"] == sign,
                    "Receiver charge differs",
                )
                margin = raw_witness(report, starts, radius)
                channel = report["reference_channel_without_global_terminal_drain"]
                if wait in margins:
                    require(
                        margin > 0
                        and channel["unfiltered_optimal_CHSH_interval"][1] < 2,
                        "Local Bell witness or its CHSH boundary differs",
                    )
                    margins[wait].append(margin)
                else:
                    require(
                        channel["Bell_overlap_interval"][1] < 0.5,
                        "Parked or dephased comparison differs",
                    )
                checked += 1
    require(checked == cohorts, "The recorded comparison family is incomplete")
    return {
        "source": "archived interacting periodic eight-site Weyl/contact ring",
        "separate_four_site_singlet_return_and_CHSH_replayed": True,
        "common_local_gate": "signed charge at 4 and single-occupancy spin at 3 or 5",
        "cohorts": checked,
        "simultaneous_raw_report_error_probability_upper": float(family_error),
        "minimum_raw_X_minus_Y_plus_Z_minus_admission_lower_by_wait": {
            wait: float(min(values)) for wait, values in margins.items()
        },
        "parked_and_resolved_analytical_Bell_overlap_upper_below_half": True,
        "local_unfiltered_CHSH_violation_claimed": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--flux",
        action="store_true",
        help="also replay the separate supplied four-site twist",
    )
    args = parser.parse_args()
    main_record = measure.run()
    compare(
        main_record, json.loads((HERE / "results.json").read_bytes()), "contact ring"
    )
    independent_record = json.loads(json.dumps(independent.run(), allow_nan=False))
    require(
        independent_record
        == json.loads((HERE / "independent/records.json").read_bytes()),
        "Independent eight-site record changed",
    )
    summary = joined_observer_records(main_record, independent_record)
    if args.flux:
        from .flux_continuation import check

        compare(
            check.run(),
            json.loads((HERE / "flux_continuation/results.json").read_bytes()),
            "separate four-site supplied twist",
        )
        summary["separate_four_site_supplied_twist_replayed"] = True
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
