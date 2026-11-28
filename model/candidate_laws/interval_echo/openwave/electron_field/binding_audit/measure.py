"""Compare the released four-face fixture to its exactly factorized parked bridge."""

import argparse
import json
from fractions import Fraction
from hashlib import sha256
from math import floor
from pathlib import Path
from runpy import run_path
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import numpy as np
from scipy.linalg import eigh

HERE = Path(__file__).resolve().parent
FIELD = HERE.parent
GAP = 1e-9


def ranges(operators, columns):
    """Outcome probability range over the specified input amplitudes."""
    result = {}
    for label, operator in operators.items():
        probabilities = np.linalg.eigvalsh(columns.conj().T @ operator @ columns)
        assert min(probabilities) > -GAP and max(probabilities) < 1 + GAP
        result[label] = np.clip(probabilities, 0, 1).tolist()
    return result


def mean_reads(operators, columns):
    """Outcome means for a fair mixture of the same two preparations."""
    return {
        label: float(
            np.trace(columns.conj().T @ operator @ columns).real / columns.shape[1]
        )
        for label, operator in operators.items()
    }


def count_below(h, cut):
    """Count eigenvalues below a rational cut using exact symmetric Bareiss elimination."""
    twice = np.rint(2 * h).astype(int)
    np.testing.assert_array_equal(twice, 2 * h)
    a = (cut.denominator * twice).tolist()
    for i, row in enumerate(a):
        row[i] -= 2 * cut.numerator
    previous = 1
    negative = 0
    for k, row in enumerate(a):
        pivot = row[k]
        assert pivot != 0
        # Each pivot is a leading principal determinant. Consecutive ratios
        # give the diagonal entries of a real LDLᵀ congruence.
        negative += int((pivot > 0) != (previous > 0))
        for i in range(k + 1, len(a)):
            for j in range(i, len(a)):
                next_entry, remainder = divmod(
                    pivot * a[i][j] - row[i] * row[j], previous
                )
                assert remainder == 0
                a[i][j] = next_entry
        previous = pivot
    return negative


def parked_reference(space, field, heads, h):
    """Remove exactly the native gateway term, then check the independent factors."""
    gateway = np.zeros_like(h)
    left, right = (port.ravel() for port in space.portal)
    for i, j in zip(left, right, strict=True):
        assert space.material_step(space.basis[i], (0, 6)) == space.basis[j]
        assert h[i, j] == -4
        gateway[i, i] = gateway[j, j] = 4
        gateway[i, j] = gateway[j, i] = -4
    parked = h - gateway
    cell_zero = np.array([x < 6 for x, _ in space.basis])
    np.testing.assert_array_equal(parked[np.ix_(cell_zero, ~cell_zero)], 0)
    neutral = field.ring_laplacian(field.words[heads])
    np.testing.assert_allclose(
        np.linalg.eigvalsh(neutral), [0, 1, 1, 2], atol=1e-14, rtol=0
    )
    charged_fans = []
    for stay in (cell_zero, ~cell_zero):
        block = parked[np.ix_(stay, stay)]
        size = len(block) // 4
        # The archive orders left words first, within each charge position.
        # Put the opposite, neutral cell's word last on both sides.
        if not stay[0]:
            positions = [x for x, _ in space.basis if x >= 6]
            lengths = np.unique(positions, return_counts=True)[1]
            rows = []
            start = 0
            for length in lengths:
                offset = np.arange(start, start + length).reshape(4, -1)
                rows.extend(offset.T.ravel())
                start += length
            block = block[np.ix_(rows, rows)]
            assert start == len(block)
        charged = block[::4, ::4] - neutral[0, 0] * np.eye(size)
        np.testing.assert_array_equal(
            block, np.kron(charged, np.eye(4)) + np.kron(np.eye(size), neutral)
        )
        charged_fans.append(charged)
    return parked, charged_fans


def run_sector(field, heads, grades, sign, times):
    from binding import other_cell_detectors
    from search import JoinedSpace, vacuum_code

    space = JoinedSpace(field, (heads, heads), sign)
    h = space.hamiltonian([6, *grades[0], 6, *grades[1]])
    parked, charged_fans = parked_reference(space, field, heads, h)
    grounds = [float(eigh(fan, eigvals_only=True)[0]) for fan in charged_fans]
    threshold = min(grounds) + 1
    lower_ground = floor(min(grounds) * 100)
    ground_cuts = [Fraction(lower_ground + offset, 100) for offset in (0, 1)]
    separated_cuts = [cut + 1 for cut in ground_cuts]
    fan_counts = [
        [count_below(fan, cut) for cut in ground_cuts] for fan in charged_fans
    ]
    assert all(low == 0 and high >= 1 for low, high in fan_counts)
    release_counts = [count_below(h, cut) for cut in separated_cuts]
    assert release_counts == [1, 1]
    vacuum, code = vacuum_code(field, heads)
    source = np.zeros((len(h), 2))
    sink = np.zeros_like(source)
    receiver = np.zeros((len(h), 4, 2))
    source[space.portal[0]] = np.einsum("ik,j->ijk", code, vacuum)
    sink[space.portal[1]] = np.einsum("i,jk->ijk", vacuum, code)
    receiver[space.portal[1]] = np.einsum("ia,jk->ijak", np.eye(4), code)
    receiver = receiver.reshape(len(h), 8)
    opposite = other_cell_detectors(space, field, heads)
    arrival = {
        "gateway_charge_and_receiver_doublet": receiver @ receiver.T,
        "gateway_charge_and_receiver_doublet_and_source_vacuum": sink @ sink.T,
    }
    reads = arrival | {key: value["projector"] for key, value in opposite.items()}
    cell_zero = np.array([x < 6 for x, _ in space.basis])
    np.testing.assert_array_equal(source[~cell_zero], 0)
    for read in arrival.values():
        np.testing.assert_array_equal(read[:, cell_zero], 0)
    for description in opposite.values():
        read = description["projector"]
        np.testing.assert_array_equal(read @ source, 0)
        np.testing.assert_array_equal(read @ parked, parked @ read)
    energy, vectors = eigh(h)
    sub = energy < threshold
    assert np.min(abs(energy - threshold)) > GAP
    assert np.count_nonzero(sub) == release_counts[0]
    inside = vectors[:, sub]
    energy_clusters = np.split(
        np.arange(len(energy)), np.flatnonzero(np.diff(energy) > GAP) + 1
    )
    modes = []
    for group in energy_clusters:
        if not sub[group].all():
            continue
        basis = vectors[:, group]
        modes.append(
            {
                "energy": float(energy[group].mean()),
                "multiplicity": len(group),
                "source_probability_range": ranges(
                    {"energy_read": basis @ basis.T}, source
                )["energy_read"],
                "projector_weights_on_eigenspace": ranges(reads, basis),
            }
        )

    def history(steps, columns, time):
        values, bases = steps
        return bases @ (np.exp(-1j * time * values[:, None]) * (bases.T @ columns))

    timed = []
    parked_modes = eigh(parked)
    for time in times:
        full = history((energy, vectors), source, time)
        low = history((energy[sub], inside), source, time)
        high = full - low
        full_read = mean_reads(reads, full)
        low_read = mean_reads(reads, low)
        high_read = mean_reads(reads, high)
        control = mean_reads(reads, history(parked_modes, source, time))
        assert max(map(abs, control.values())) < 1e-12
        timed.append(
            {
                "release_time": time,
                "full_arrival_probability_ranges": ranges(arrival, full),
                "full_mean_probabilities": full_read,
                "below_threshold_mean_probabilities": low_read,
                "above_threshold_mean_probabilities": high_read,
                "interference_contributions_to_full_probabilities": {
                    name: full_read[name] - low_read[name] - high_read[name]
                    for name in full_read
                },
                "parked_control_probabilities_exact": dict.fromkeys(control, 0),
            }
        )
    return {
        "sign": sign,
        "heads_each_side": list(heads),
        "dimension": len(energy),
        "parked_charged_fan_ground_energies": grounds,
        "opposite_neutral_doublet_energy": 1,
        "parked_separated_threshold": threshold,
        "exact_threshold_dimension_certificate": {
            "parked_charged_ground_bracket": [str(cut) for cut in ground_cuts],
            "parked_charged_fan_counts_below_bracket_endpoints": fan_counts,
            "separated_threshold_bracket": [str(cut) for cut in separated_cuts],
            "released_counts_below_bracket_endpoints": release_counts,
        },
        "closest_released_level_to_threshold": float(np.min(abs(energy - threshold))),
        "released_energy_eigenspaces_below_threshold": modes,
        "source_probability_range_below_threshold": ranges(
            {"threshold_filter": inside @ inside.T}, source
        )["threshold_filter"],
        "subspace_detector_probability_ranges": ranges(reads, inside),
        "release_reads": timed,
    }


def run():
    # Verify the inherited experiment inventory and run only its shared fixture builders.
    manifest = dict(
        line.split("  ", maxsplit=1)[::-1]
        for line in (FIELD / "records.sha256").read_text().splitlines()
    )
    model_root = FIELD.parents[4]
    sys.path.insert(0, str(model_root))
    with TemporaryDirectory() as work, ZipFile(FIELD / "records.zip") as archive:
        root = Path(work)
        assert sorted(archive.namelist()) == sorted(manifest)
        assert all(
            sha256(archive.read(name)).hexdigest() == digest
            for name, digest in manifest.items()
        )
        (root / model_root.name).symlink_to(model_root, target_is_directory=True)
        archive.extractall(root)
        package = root / "experiments" / "interval_echo_collective_spin"
        sys.path.insert(0, str(package))
        from search import fan

        record = json.loads((package / "results.json").read_text())
        best = record["headline"]["best_worst_qubit_success"]
        reference = best["release"]
        assert (
            reference["sign"] == 1
            and reference["head_sectors"][0] == reference["head_sectors"][1]
        )
        heads = tuple(reference["head_sectors"][0])
        field = fan(4)
        conjugate = tuple(
            sum(i in e for e in field.edges) - h for i, h in enumerate(heads)
        )
        times = [reference["best_receiver_only_time"], reference["best_worst_time"]]
        assert times == [27.0625, 30.8125]
        positive = run_sector(field, heads, best["degrees"], 1, times)
        negative = run_sector(field, conjugate, best["degrees"], -1, times)
        local, shared = positive["release_reads"]
        np.testing.assert_allclose(
            local["full_mean_probabilities"]["gateway_charge_and_receiver_doublet"],
            reference["receiver_only_mean_joint_success_at_that_time"],
            atol=1e-12,
            rtol=0,
        )
        np.testing.assert_allclose(
            shared["full_arrival_probability_ranges"][
                "gateway_charge_and_receiver_doublet_and_source_vacuum"
            ],
            [
                reference["best_worst_qubit_success"],
                reference["better_input_success_at_same_time"],
            ],
            atol=1e-12,
            rtol=0,
        )
        np.testing.assert_allclose(
            shared["full_mean_probabilities"]["opposite_cell_doublet"],
            sum(
                value
                for name, value in record["selected_witness_controls"].items()
                if name.endswith("at_selected_time") and "doublet" in name
            ),
            atol=1e-12,
            rtol=0,
        )
        return {
            "physical_ports": best["actual_ports"],
            "rim_degrees": best["degrees"],
            "archived_collective_search_sha256": manifest[
                "experiments/interval_echo_collective_spin/search.py"
            ],
            "archived_collective_receipt_sha256": manifest[
                "experiments/interval_echo_collective_spin/results.json"
            ],
            "sectors": [positive, negative],
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = run()
    if arguments.output:
        with arguments.output.open("x") as handle:
            handle.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        compare = run_path(str(FIELD / "replay.py"))["compare"]
        compare(
            result, json.loads((HERE / "results.json").read_text()), "binding_audit"
        )
        print("Verified the matched two-cell binding read.")


if __name__ == "__main__":
    main()
