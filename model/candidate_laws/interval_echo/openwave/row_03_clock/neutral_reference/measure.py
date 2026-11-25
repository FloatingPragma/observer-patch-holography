"""Check which neutral rest phases the native observer grammar can read."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import candidate

HERE = Path(__file__).resolve().parent
SHOTS, SEED = 2**18, 178303


def checked(condition, message):
    if not condition:
        raise ArithmeticError(message)


def close(left, right, message):
    checked(np.max(np.abs(np.asarray(left) - np.asarray(right))) < 8e-13, message)


def primitives():
    echo = law.echo()
    kinetic, contrast = law.local_term_parts()
    h = kinetic + candidate.KAPPA * contrast
    charge = np.tile([a + b - 2 for a, b in law.SRC], 1 + law.N_LAB)
    block = [law.SRC.index((1, 1))] + [
        (1 + law.LABELS.index((coin, 0))) * law.N_SRC + law.SRC.index((1, 1))
        for coin in (0, 1)
    ]
    projector = np.zeros(law.DIM)
    projector[block] = 1
    # Charge conservation of the echo follows over Q(sqrt(2)) directly from
    # the integer write support and the native reversal's integer support.
    wp = law._write_pattern()
    r, _ = law._label_ops()
    for row, col in zip(*np.nonzero(wp), strict=True):
        checked(charge[law.N_SRC + row] == charge[col], "write changes total charge")
    checked(all(law.LABELS[i][1] ** 2 == law.LABELS[j][1] ** 2
                for i, j in zip(*np.nonzero(r), strict=True)), "reversal changes contrast square")
    for name, operator in (("echo", echo), ("exchange", kinetic), ("contrast", contrast)):
        close((charge[:, None] - charge) * operator, 0, name + " changes charge")
        close((projector[:, None] - projector) * operator, 0, name + " mixes C")
    close(h[:, block], 0, "uniform zero-contrast state acquires a nontrivial phase")
    # Every load/memory pointer is diagonal in this basis, hence commutes
    # exactly with both diagonal projectors. No enumeration of a subset is used.
    full = np.eye(law.DIM)
    vacuum = full[:, law.SRC.index((1, 1))]
    pair = (full[:, law.SRC.index((0, 2))] + full[:, law.SRC.index((2, 0))]) / np.sqrt(2)
    close(h @ pair, 4 * candidate.KAPPA * pair, "native symmetric pair energy differs")
    close(h @ vacuum, 0, "filled reference energy differs")
    written_vacuum, written_pair = echo @ vacuum, echo @ pair
    close(echo @ (written_vacuum + written_pair) / np.sqrt(2),
          (vacuum + pair) / np.sqrt(2), "native emission differs")
    return echo, h, charge, block, vacuum, pair


def transcript(echo, h, state):
    """All 11 first memory outputs, feedback, and all 9 terminal load outputs."""
    result = []
    for record in range(1 + law.N_LAB):
        branch = np.zeros(law.DIM, complex)
        region = slice(record * law.N_SRC, (record + 1) * law.N_SRC)
        branch[region] = (echo @ state)[region]
        # Distinct classical first records choose distinct native waits.
        followed = echo @ (expm(-1j * h * ((record + 1) / 128)) @ branch)
        result.extend(np.sum(np.abs(followed.reshape(-1, law.N_SRC)) ** 2, axis=0))
    checked(abs(sum(result) - 1) < 1e-12, "native transcript loses a complement")
    return np.asarray(result)


def sample_counts(rng, probabilities):
    """Retain raw positive predictions; account explicitly for floating norm error."""
    raw = np.asarray(probabilities, dtype=float)
    checked(np.all(np.isfinite(raw)) and np.min(raw) >= 0,
            "negative or nonfinite detector probability")
    norm_error = float(abs(np.sum(raw) - 1))
    checked(norm_error < 2e-12, "detector normalization error exceeds budget")
    # A numerical unitary gives a distribution with a bounded norm residual.
    # Use the largest outcome as the explicit complement to all other raw
    # probabilities. No negative entry is clipped and no vector is rescaled.
    last = int(np.argmax(raw))
    law_for_sampling = raw.copy()
    law_for_sampling[last] = 1 - sum(raw[i] for i in range(len(raw)) if i != last)
    checked(law_for_sampling[last] >= 0, "negative explicit complement")
    bound = float(np.sum(abs(law_for_sampling - raw)))
    checked(bound < 2e-12, "sampling complement error exceeds budget")
    order = [i for i in range(len(raw)) if i != last] + [last]
    counts = np.empty(len(raw), dtype=np.int64)
    counts[order] = rng.multinomial(SHOTS, law_for_sampling[order])
    return counts.tolist(), bound


def branch_boundary(h, charge):
    eye = np.eye(law.DIM)
    state = (eye[:, law.SRC.index((0, 1))] + eye[:, law.SRC.index((1, 0))]) / np.sqrt(2)
    close(h @ state, candidate.KAPPA * state, "control sample is not a native eigenstate")
    close(charge * state, -state, "control sample charge differs")
    t, offset = np.pi / 32, -16
    amplitude = np.vdot(state, expm(-1j * h * t) @ state)
    conditional = np.vdot(state, expm(-1j * (h + np.diag(offset * charge)) * t) @ state)
    # True global reference shift acts on the idle branch too.
    shifted_idle = np.exp(1j * offset * t)
    p = (1 + amplitude.real) / 2
    wrong = (1 + conditional.real) / 2
    global_shift = (1 + (shifted_idle.conjugate() * conditional).real) / 2
    close(p, global_shift, "uniform global convention changed the control record")
    checked(abs(p - wrong) > 0.70, "controlled-vs-idle boundary is insensitive")
    return dict(exposure=float(t), offset=offset, native_plus=float(p),
                shifted_every_branch_plus=float(global_shift), conditional_shift_only_plus=float(wrong),
                distinction="A conditional energy drive is not a uniform global energy convention")


def changed_creation_control(h, charge, filled, pair):
    """A charge- and seam-swap-preserving rival can create the same neutral pair."""
    coupling = 16 * (np.outer(filled, pair) + np.outer(pair, filled))
    close((charge[:, None] - charge) * coupling, 0, "creation rival changes total charge")
    # The added term joins only the same seam's existing load words and has
    # the same amplitude for exchanging its two endpoints.
    swap = np.tile([law.SRC.index((b, a)) for a, b in law.SRC], 1 + law.N_LAB)
    swap += np.repeat(np.arange(1 + law.N_LAB) * law.N_SRC, law.N_SRC)
    close(coupling[np.ix_(swap, swap)], coupling, "creation rival picks an endpoint")
    t = np.pi / (32 * np.sqrt(2))
    pairs = np.tile([int(a + b == 2 and a != b) for a, b in law.SRC], 1 + law.N_LAB)
    native = expm(-1j * h * t) @ filled
    changed = expm(-1j * (h + coupling) * t) @ filled
    probabilities = [float(pairs @ abs(v)**2) for v in (native, changed)]
    close(probabilities, [0, 0.5], "same-reader neutral-creation sensitivity fails")
    return dict(exposure=float(t), added_blank_seam_pair_mixing=16,
                native_pair_census_probability=probabilities[0],
                changed_law_pair_census_probability=probabilities[1],
                scope="Changed local law on the same99 alphabet; not Interval Echo")


def run():
    echo, h, charge, block, vacuum, pair = primitives()
    rng = np.random.default_rng(SEED)
    rows, largest = [], 0.0
    for t in np.linspace(0, 0.5, 17):
        # Two checked native eigenstates give the exact reduced free evolution.
        evolved_pair = np.exp(-4j * candidate.KAPPA * t) * pair
        pure = (vacuum + evolved_pair) / np.sqrt(2)
        coherent = transcript(echo, h, pure)
        dephased = (transcript(echo, h, vacuum) + transcript(echo, h, evolved_pair)) / 2
        largest = max(largest, float(np.max(np.abs(coherent - dephased))))
        close(coherent, dephased, "primitive native records resolve cross-C phase")
        # A supplied analyzer outside the primitive grammar can interfere C
        # with its complement. This comparison tests the missing operation.
        plus, minus = (vacuum + pair) / np.sqrt(2), (vacuum - pair) / np.sqrt(2)
        ramsey = [float(abs(np.vdot(v, pure)) ** 2) for v in (plus, minus)]
        close(ramsey, [(1 + np.cos(32 * t)) / 2, (1 - np.cos(32 * t)) / 2],
              "supplied analyzer misses the prepared neutral phase")
        native_counts, native_bound = sample_counts(rng, coherent)
        ramsey_counts, ramsey_bound = sample_counts(rng, ramsey)
        rows.append(dict(exposure=float(t), complete_native_probabilities=coherent.tolist(),
                         complete_native_counts=native_counts,
                         native_sampling_l1_correction_bound=native_bound,
                         supplied_ramsey_probabilities=ramsey,
                         supplied_ramsey_counts=ramsey_counts,
                         ramsey_sampling_l1_correction_bound=ramsey_bound,
                         dephased_ramsey_probabilities=[0.5, 0.5]))
    # The interferometer in branch_boundary retains a coherent control, which
    # is why its conditional-sample shift has to be separated from a convention.
    boundary = branch_boundary(h, charge)
    return dict(
        source="Interval Echo native 99-state seam at kappa=8; isolated addressed schedule",
        scope="Primitive-grammar reachability and global charge-reference identifiability",
        shots_per_appointment=SHOTS, seed=SEED,
        primitive_uniform_subspace_indices=block,
        neutral_pair_energy_above_filled=32,
        maximum_native_coherent_minus_dephased_probability=largest,
        uniform_to_pair_probability_all_primitive_histories=0,
        native_outcomes="memory record 0..10, then load word index 0..8, in row-major order",
        records=rows, coherent_control_boundary=boundary,
        changed_law_creation_control=changed_creation_control(h, charge, vacuum, pair),
        charged_reference_conventions=[
            dict(offset=0, vacancy=35, excess=67, neutral_pair=102),
            dict(offset=-16, vacancy=51, excess=51, neutral_pair=102),
        ],
        exact_two_convention_minimax_absolute_error_lower_bound=8,
        official_status="OW-03 untested at OpenWave 0df73dc6; local full criterion partial",
        obligation="Physical rest-energy selection and a self-starting internal rest clock are unresolved",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = run()
    if args.output:
        with args.output.open("x") as stream:
            json.dump(report, stream, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        from ...electron_field.replay import compare
        compare(report, json.loads((HERE / "results.json").read_text()), "neutral_reference")
        print("Replayed native neutral-record indistinguishability and the conditional-drive boundary.")
