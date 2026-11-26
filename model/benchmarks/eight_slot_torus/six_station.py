"""Six local ports on the free one-constituent sector of the pinned torus law.

The labels 0 through 5 refer to +x, -x, +y, -y, +z, -z on a periodic
5 by 5 by 5 carrier with one private position. A station can launch its
normalized equal-slot mode or swap it with an empty isolated latch. Reading
latch occupation distinguishes a hit from the orthogonal no-click state.
A hit leaves the carrier vacant. A reverse swap releases the same slot mode.
These apparatus operations and synchronized address counters are supplied.
They do not use the two-constituent collision vertex.

All times are integer updates. Each window lasts 12 updates. The marked
address (0,0,0,0) applies a sign to its eight slots immediately after
updates 3 and 15, before queries. The sham applies the identity on the
same schedule. B queries after each update in the first window. On a hit
it parks the excitation and records j. In the range menu B releases at
tick 12 after its last query. A queries at ticks 13 through 24 and records
the relative tick k. A's elapsed time minus the logged dwell at B is
(12+k)-(12-j)=j+k. Two successful local records define the range (j+k)/2.
No per-trial result travels ahead of the wave to enable a remote service.

In the challenge menu A launches at tick 0 or 6. B decodes its first hit
at 1..6 or 7..12 and releases at tick 12 or 18. A decodes its own first
click by the same two slots relative to tick 12. A silent B keeps a hit
parked. A disabled B makes no query. In every menu the actual no-click
wave continues into A's window. A zero on either time axis is its timeout.
All clocks, controls and challenges are explicit inputs. Comparing a pair
of local records after the run requires their classical collection.

The review's historical outcomes were known when this finite design was
chosen. Coordinates, horizons, phases, trial counts, controls and seed are
fixed before this acquisition. The range certificate concerns a separate
fixed 5-percent question on these measurements. The station addresses do
not enter the distance reader.
"""
from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction
from functools import lru_cache
from io import BytesIO
from itertools import permutations
from pathlib import Path

import numpy as np

from .law import EQUAL, LAWS, step

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "six_station_records.npz"
DIGEST = HERE / "six_station_records.sha256"
SOURCE = "3db6e6d74416c26ca0f2e1be4bb3da2f6ad9d5fa"
SOURCE_DIGESTS = {
    "law.py": "d0f1b4221671a4867cdcde097f484f6085bf334221b2a068be53f9cb5fd0dbd9",
    # The relocated benchmark specification is pinned separately from the executable law.
    "six_station_source/SPEC.md": "466b1dc20f278d83ec79796cacf75e2e3441c21e0680bcbc1a15b455745519e6",
}
CARRIER = (5, 5, 5, 1)
MARK = (0, 0, 0, 0)
STATIONS = ((1, 0, 0, 0), (4, 0, 0, 0), (0, 1, 0, 0),
            (0, 4, 0, 0), (0, 0, 1, 0), (0, 0, 4, 0))
HORIZON = 12
SLOT = 6
PULSE_TICK = 3
ACTIVE, SILENT, DISABLED = range(3)
GRID = 2**48
RANGE_SHOTS = 2**34
TAG_SHOTS = 2**20
SEED = 122303
SOURCE_TOLERANCE = 2e-13
# Two-sided Bernstein bound over 420 menus with 169 cells each:
# 2 * 420 * 169 * exp(-27) < 2.67e-7 under the independent-trial model.
COUNT_LOG_TAIL = 27


def check(condition, message):
    if not condition:
        raise ValueError(message)


def check_source():
    for name, digest in SOURCE_DIGESTS.items():
        data = (HERE / name).read_bytes().replace(b"\r\n", b"\n")
        check(hashlib.sha256(data).hexdigest() == digest, f"pinned source changed: {name}")


def literal_history(sender, receiver, marked, control=ACTIVE, bit=0, tagged=False):
    """Run the literal update with B's time record and its local latch.

    In the first window a hit transfers a complex port amplitude into an
    orthogonal time record and an isolated latch. In the second window
    the latch swap puts this amplitude back into B's launch mode. The
    source, queries and phase gates below use the absolute lab tick.
    """
    source, target = STATIONS[sender], STATIONS[receiver]
    history = np.zeros((HORIZON + 1, HORIZON + 1))
    vacant = np.zeros(CARRIER + (8,), complex)
    wave = vacant.copy()
    held = {}

    def update(state, clock):
        advanced = step(state, LAWS["torus"])
        if marked and clock in (PULSE_TICK, HORIZON+PULSE_TICK):
            advanced[MARK] *= -1
        return advanced

    def swap_out(state, station):
        parked = np.vdot(EQUAL, state[station])
        state[station] -= EQUAL * parked
        return parked

    for tick in range(1, HORIZON + 1):
        if tick == 1 + bit*SLOT:
            wave[source] = EQUAL
        wave = update(wave, tick)
        if control != DISABLED:
            amplitude = swap_out(wave, target)
            if amplitude != 0:
                held[tick] = amplitude
    # B's no-click outcome carries its unnormalized outgoing wave into A's window.
    branches = {0: wave}
    for j, parked in held.items():
        if control == ACTIVE:
            branches[j] = vacant.copy()
        else:
            history[j, 0] = abs(parked)**2
    for k in range(1, HORIZON+1):
        for j, state in branches.items():
            if j:
                release = SLOT if tagged and j > SLOT else 0
                if k <= release:
                    continue
                if k == release+1:
                    state[target] = EQUAL * held[j]
            state = update(state, HORIZON+k)
            history[j, k] = abs(swap_out(state, source))**2
            branches[j] = state
    for j, residual in branches.items():
        history[j, 0] = np.vdot(residual, residual).real
    return history


class Dyadic:
    """Independent integer coin, reflected edges, and equal-port instrument.

    The state is x/(sqrt(8)*d). After the coin and routing the denominator
    is 4d. The equal-port hit amplitude is the local slot sum divided by
    32d. Subtract its vector from the wave at denominator 32d. All weights
    including the final survivor norm are rational.
    """

    @staticmethod
    def launch(site):
        x = np.zeros(CARRIER + (8,), object)
        x[site] = 1
        return x, 1

    @staticmethod
    def norm(wave):
        x, denominator = wave
        return Fraction(np.sum(x*x), 8 * denominator**2)

    @staticmethod
    def monitor(wave, site, start, marked):
        x, denominator = wave
        hits = np.full(HORIZON + 1, Fraction(0), object)
        for tick in range(start + 1, HORIZON + 1):
            totals = x.sum(axis=-1)
            routed = np.empty_like(x)
            for axis in range(4):
                for slot, direction in ((2*axis, 1), (2*axis+1, -1)):
                    routed[..., slot ^ 1] = np.roll(totals - 4*x[..., slot], direction, axis)
            if tick == PULSE_TICK and marked:
                routed[MARK] = -routed[MARK]
            denominator *= 32
            x = 8*routed
            if site is not None:
                slot_sum = sum(routed[site])
                hits[tick] = Fraction(slot_sum, denominator)**2
                x[site] -= slot_sum
        return hits, (x, denominator)


def exact_menus(sender, receiver, marked):
    """Full time histories per range or (switch, challenge) menu."""
    source, target = STATIONS[sender], STATIONS[receiver]
    check(source != target, "two distinct stations required")

    @lru_cache(maxsize=None)
    def return_law(delay):
        hits, residual = Dyadic.monitor(Dyadic.launch(target), source, delay, marked)
        hits[0] = Dyadic.norm(residual)
        return hits

    @lru_cache(maxsize=None)
    def outbound(bit, enabled):
        hits, missed = Dyadic.monitor(Dyadic.launch(source), target if enabled else None,
                                      bit*SLOT, marked)
        false_hits, survivor = Dyadic.monitor(missed, source, 0, marked)
        false_hits[0] = Dyadic.norm(survivor)
        return hits, false_hits

    def histories(control, bit, tagged):
        hits, false_hits = outbound(bit, control != DISABLED)
        result = np.zeros((HORIZON + 1, HORIZON + 1), dtype=hits.dtype)
        result[0] = false_hits
        for j in range(1, HORIZON + 1):
            if control == ACTIVE:
                delay = SLOT if tagged and j > SLOT else 0
                result[j] = hits[j] * return_law(delay)
            else:
                result[j, 0] = hits[j]
        return result

    radar = histories(ACTIVE, 0, False)
    tag = np.stack([[histories(control, bit, True) for bit in range(2)]
                    for control in (ACTIVE, SILENT, DISABLED)])
    return radar, tag


def integer_weights(law):
    """Largest remainders on the exact rational law; row-major ties."""
    flat = law.ravel()
    check(all(q >= 0 for q in flat) and sum(flat) == 1, "exact history law incomplete")
    scaled = [q*GRID for q in flat]
    weights = np.asarray([int(value) for value in scaled], dtype=np.int64)
    left = GRID - int(weights.sum())
    remainders = [value - int(weight) for value, weight in zip(scaled, weights)]
    ranked = sorted(range(len(flat)), key=remainders.__getitem__, reverse=True)
    weights[ranked[:left]] += 1
    return weights.reshape(law.shape)


def source_weights():
    """Compare every inclusive history to the literal torus update before quantizing it.

    The maximum absolute rounding error per categorical cell is 1/GRID,
    so TV from the exact law is at most (HORIZON+1)^2/(2*GRID) per menu.
    """
    check_source()
    axes = (2, len(STATIONS), len(STATIONS))
    times = (HORIZON+1, HORIZON+1)
    weights = {"range_time_weights": np.zeros(axes + times, np.int64),
               "tag_time_weights": np.zeros(axes + (3, 2) + times, np.int64)}
    worst = 0.
    for sender, receiver in permutations(range(len(STATIONS)), 2):
        for marked in (0, 1):
            exact = exact_menus(sender, receiver, marked)
            literal = (literal_history(sender, receiver, marked),
                       np.array([[literal_history(sender, receiver, marked, control, bit, True)
                                  for bit in (0, 1)] for control in (ACTIVE, SILENT, DISABLED)]))
            for key, rational, actual in zip(weights, exact, literal):
                error = float(np.max(abs(rational.astype(float) - actual)))
                check(error <= SOURCE_TOLERANCE, "literal and integer history laws disagree")
                worst = max(worst, error)
                for index in np.ndindex(rational.shape[:-2]):
                    weights[key][(marked, sender, receiver) + index] = integer_weights(rational[index])
    return weights, worst


def sample(weights):
    """One seeded multinomial of the entire time histogram per menu.

    Integer weights divided by 2^48 are exact binary64 probabilities.
    Finite-count coverage uses the independent categorical sampling model.
    The TV bound concerns the input law; this experiment supplies no
    probability bound on deviations of a deterministic pseudorandom source.
    """
    counts = {}
    for number, (key, shots) in enumerate((("range_time_weights", RANGE_SHOTS),
                                           ("tag_time_weights", TAG_SHOTS))):
        laws = weights[key]
        values = np.zeros_like(laws)
        for index in np.ndindex(laws.shape[:-2]):
            if index[1] == index[2]:
                continue
            rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED, number, *index])))
            drawn = rng.multinomial(shots, laws[index].ravel() / GRID)
            values[index] = drawn.reshape(laws.shape[-2:])
        counts[key.replace("weights", "counts")] = values
    return counts


def validate_counts(weights, counts):
    """Check conservation, support, and simultaneous binomial deviations.

    After acquisition, this diagnostic tests the archived counts against the
    reconstructed sampling law. Each multinomial cell has marginal variance
    N*p*(1-p). Bernstein's two-sided tail at deviation d is bounded by
    2 exp[-d²/(2 N*p*(1-p)+2d/3)]. The integer inequality below compares
    this exponent to 27 without rounding a rare probability or a residual.
    """
    for key, shots in (("range_time_weights", RANGE_SHOTS), ("tag_time_weights", TAG_SHOTS)):
        laws = weights[key]
        values = counts[key.replace("weights", "counts")]
        check(values.shape == laws.shape and values.dtype == np.dtype("int64"),
              "invalid count histogram")
        check(np.all((values >= 0) & (values <= shots)) and not np.any((laws == 0) & (values != 0)),
              "count outside source support")
        check(np.array_equal(values.sum(axis=(-2, -1)),
                             (laws.sum(axis=(-2, -1)) != 0).astype(np.int64)*shots),
              "selection or loss of launches")
        w = laws.astype(object)
        error = abs(values.astype(object)*GRID - shots*w)
        variance = shots*w*(GRID-w)
        check(np.all(3*error**2 <= 2*COUNT_LOG_TAIL*(3*variance + GRID*error)),
              "histogram cell outside the simultaneous source-count bound")


def acquire(path):
    weights, error = source_weights()
    counts = sample(weights)
    validate_counts(weights, counts)
    with Path(path).open("xb") as dest:
        np.savez_compressed(dest, **weights, **counts)
    return error


def verify_archive(path=ARCHIVE, seed_replay=False):
    """Regenerate both source distributions and score the archived counts."""
    path = Path(path)
    expected_digest = DIGEST.read_text(encoding="ascii").split()[0]
    payload = path.read_bytes()
    check(hashlib.sha256(payload).hexdigest() == expected_digest,
          "acquisition checksum mismatch")
    with np.load(BytesIO(payload), allow_pickle=False) as archive:
        expected = {f"{kind}_time_{stat}" for kind in ("range", "tag") for stat in ("weights", "counts")}
        check(len(archive.files) == len(expected) and set(archive.files) == expected,
              "unexpected or duplicate record menu")
        weights, _ = source_weights()
        for key, value in weights.items():
            check(archive[key].dtype == value.dtype and np.array_equal(archive[key], value),
                  f"archived integer inputs differ from source: {key}")
        counts = {key: archive[key] for key in expected if key.endswith("counts")}
        validate_counts(weights, counts)
        if seed_replay:
            for key, values in sample(weights).items():
                check(np.array_equal(values, counts[key]), "multinomial draw differs on this backend")
    return {**weights, **counts}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--acquire", type=Path, help="create an acquisition at an unused path")
    parser.add_argument("--records", type=Path, default=ARCHIVE, help="archive to verify")
    parser.add_argument("--seed-replay", action="store_true", help="also repeat the backend's original draws")
    args = parser.parse_args()
    if args.acquire is None:
        verify_archive(args.records, args.seed_replay)
        print("Verified six-station counts and integer source weights")
    else:
        error = acquire(args.acquire)
        print(f"Maximum absolute literal/integer probability difference {error:.3g}")


if __name__ == "__main__":
    main()
