"""A one-port intensity clock on a heralded native pair in two finite rings.

The two ordinary inputs occupy disjoint one- and two-link launch sites.
Each is in the inverse-coin image of its inward route. At the second
vertex appointment the local wire head heralds one birthplace. On rail
zero the vertex writes a scalar times that site's equal-slot vector, so
the complete conditional material law is independent of incoming width.
The wire census distinguishes every other birthplace and the unmarked
sector. A refusal ends that trial. Parallel preparations at each age avoid
intermediate measurement of any retained material.

After the 128 opening updates, a fixed receiver one link away tests the
matched axial difference (|4,4>-|5,5>)/sqrt(2) at its own address. It
records this click, an orthogonal local click or a spatial miss. The
launch amplitudes, wire queries, phase-sensitive local slot mixer, unit
service clock, and finite torus are supplied. The mature clock's intensity
comes from the native evolution; the apparatus makes no spectral selection.

At N=5 or 6 the local odd port kills the zero and Nyquist modes. Reflection
combines the two signed orders within each of the remaining absolute orders
1 and 2. Both energy signs enter for both groups. Four nonzero angular
frequencies can therefore enter the intensity. Same-sign pairs at the
difference frequency are indistinguishable in intensity. A uniform random
translation of both launches and the chosen herald, keeping the detector
fixed, erases interference between different momenta at fixed modal power.

The native waveforms and pole magnitudes on these rings were explored
before choosing the holdout. The count acquisition and tests below are
fixed before drawing its independent counts. The width theorem concerns
the registered birth and finite ring. It supplies no width-independent
absolute reception, infinite-volume spectral line or particle assignment.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from functools import lru_cache
from io import BytesIO
from pathlib import Path

import numpy as np

from ..formation import ANGLE, capture
from ..law import EQUAL, LAWS, coin, step_contact, step_maturing, walk_matrix

SIZES = (5, 6)
FAR_POWER = (0., .125, .25, .5)
MENUS = ("ordinary", "reflected", "translated", "vertex_off", "reversed_collision", "competing_birth")
OUTCOMES = ("earlier_birth", "other_birthplace", "unmarked", "pair_away",
            "pair_local_complement", "pair_local_odd")
OPENING = 128
AGES = np.arange(OPENING, OPENING + 128)
SHOTS = 2**32
GRID = 2**44
SEED = 225031
NATIVE = LAWS["torus"]
ACQUISITION_SHA256 = "aa3845fd025d27a4cb7a6654bfeb93a8ae2ea167c34fd0a1bc4b36699686a8b2"


def check(value, message):
    if not value:
        raise ValueError(message)


def ordinary_birth(size, far, vertex=True, reverse=False, competing=False):
    """Apply the two scheduled literal vertices to separated one-rail trains."""
    check(size in SIZES and 0 <= far <= 1, "invalid ordinary train")
    carrier = (1, 1, size, 1, 8)
    one = np.eye(8)
    trains = []
    for direction, slot in ((1, 4), (-1, 5)):
        wave = np.zeros(carrier, complex)
        for distance, power in ((1, 1-far), (2, far)):
            z = -direction*distance
            if competing and direction == 1 and distance == 2:
                z = 0
            wave[0, 0, z % size, 0] = np.sqrt(power)*coin(one[slot])
        trains.append(wave)
    if reverse:
        trains.reverse()
    free = np.tensordot(trains[0], trains[1], axes=0)
    angle = ANGLE if vertex else 0.
    free, early = capture(free, angle)
    free = step_contact(free, NATIVE.rail_contact_phases[0])
    free, born = capture(free, angle)
    selected = np.zeros(carrier, complex)
    elsewhere = 0.
    for site, branch in born:
        if site == 0:
            selected = branch
        else:
            elsewhere += float(np.vdot(branch, branch).real)
    refused = [sum(float(np.vdot(b, b).real) for _, b in early),
               elsewhere, float(np.vdot(free, free).real)]
    check(abs(sum(refused) + np.vdot(selected, selected).real - 1) < 2e-13, "wire census incomplete")
    return np.array(refused), selected


def slot_read(wave):
    """Return absolute selected-branch local incidence and odd-port probabilities."""
    slots = wave[0, 0, :, 0]
    return np.sum(abs(slots)**2, axis=1), abs((slots[:, 4]-slots[:, 5])/np.sqrt(2))**2


def native_words(size, far, menu):
    """At each age return six unconditional terminal POVM probabilities."""
    check(menu in MENUS, "unknown observation menu")
    refused, wave = ordinary_birth(size, far, menu != "vertex_off", menu == "reversed_collision",
                                  menu == "competing_birth")
    weight = float(np.vdot(wave, wave).real)
    words = []
    for n in range(1, int(AGES[-1])+1):
        wave = step_maturing(wave, NATIVE, n-1)
        if n >= OPENING:
            total, port = slot_read(wave)
            site = -1 if menu == "reflected" else 1
            if menu == "translated":
                incidence, click = total.mean(), port.mean()
            else:
                incidence, click = total[site], port[site]
            words.append([*refused, weight-incidence, incidence-click, click])
    result = np.array(words)
    check(result.min() >= -2e-14 and np.max(abs(result.sum(axis=1)-1)) < 2e-12,
          "incomplete point detector")
    return result


@lru_cache(maxsize=2)
def poles(size):
    """Resolve the signed mature poles of the actual opening into one local port."""
    check(size in SIZES, "two enrolled axial rings required")
    _, born = ordinary_birth(size, 0.)
    wave = born / np.linalg.norm(born)
    for age in range(OPENING):
        wave = step_maturing(wave, NATIVE, age)
    transform = np.fft.fft(wave[0, 0, :, 0], axis=0)
    phi = NATIVE.rail_contact_phases[0]
    energy, amplitudes = [], np.zeros((2, 2), complex)
    dark_orders = [0] + ([size//2] if size % 2 == 0 else [])
    dark_leak = 0.
    for order in range(size):
        k = 2*np.pi*order/size
        w = np.exp(.5j*phi)*walk_matrix([0, 0, k, 0], np.exp(-1j*phi))
        state = transform[order]
        # At self-inverse orders inversion symmetry fixes both axial slots
        # to the same amplitude, before and after each future update.
        if order in dark_orders:
            for vector in (state, w@state):
                dark_leak = max(dark_leak, float(abs(vector[4]-vector[5])))
            continue
        h = min(order, size-order)
        e = np.arccos((3+np.cos(k))*np.cos(phi/2)/4)
        if order == h:
            energy.append(e)
        for index, s in enumerate((1, -1)):
            root = np.exp(-1j*s*e)
            v = (w@state-root.conjugate()*state)/(root-root.conjugate())
            check(np.linalg.norm(w@v-root*v) < 1e-12, "opening outside the two bright poles")
            amplitudes[h-1, index] += np.exp(1j*k)*(v[4]-v[5])/(size*np.sqrt(2))
    check(dark_leak < 1e-12, "self-inverse order entered odd port")
    energy = np.asarray(energy)
    # Rows are complex Fourier coefficients at positive angular frequency.
    a, b = amplitudes
    lines = np.array([b[0]*a[0].conj()+a[1]*b[1].conj(),
                      b[0]*a[1].conj()+a[0]*b[1].conj(),
                      a[0]*a[1].conj(), b[0]*b[1].conj()])
    return energy, amplitudes, lines


def predicted_intensity(size):
    e, amplitudes, _ = poles(size)
    phase = e[:, None]*np.array([1, -1])[None, :]
    return abs(np.exp(-1j*np.outer(AGES-OPENING, phase.ravel())) @ amplitudes.ravel())**2


def categorical_weights(words):
    """Round each six-way distribution on a shared dyadic grid."""
    scaled = np.maximum(words, 0.) * GRID
    weights = np.floor(scaled).astype(np.int64)
    for n, row in enumerate(weights):
        left = GRID - int(row.sum())
        check(0 <= left <= len(row), "grid apportionment failed")
        order = np.argsort(-(scaled[n]-row), kind="stable")
        row[order[:left]] += 1
    check(not np.any((words == 0) & (weights != 0)), "rounded a forbidden outcome into support")
    return weights


def acquire():
    """Return all input weights and independently sampled counts in published axis order."""
    shape = (len(SIZES), len(FAR_POWER), len(MENUS), len(AGES), len(OUTCOMES))
    weights = np.zeros(shape, np.int64)
    for i, size in enumerate(SIZES):
        for j, far in enumerate(FAR_POWER):
            for m, menu in enumerate(MENUS):
                weights[i, j, m] = categorical_weights(native_words(size, far, menu))
    return weights, draw_counts(weights)


def draw_counts(weights):
    """Draw from the retained integer input laws with a separate seed per appointment."""
    draws = np.zeros_like(weights)
    for index in np.ndindex(weights.shape[:-1]):
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED, *index])))
        draws[index] = rng.multinomial(SHOTS, weights[index]/GRID)
    return draws


def prefix_fit(counts):
    """One blind frequency from the first 64 appointments; forecast all 128."""
    y = counts[:64, -1] / SHOTS
    t = np.arange(64)
    candidates = np.arange(1, 3142)/1000
    def scan(frequencies):
        phase = frequencies[:, None]*t
        design = np.stack([np.ones_like(phase), np.cos(phase), np.sin(phase)], axis=2)
        coefficients = np.linalg.pinv(design) @ y
        residual = (design @ coefficients[:, :, None])[:, :, 0] - y
        chosen = np.argmin(np.sum(residual**2, axis=1))
        return frequencies[chosen], coefficients[chosen]
    coarse, _ = scan(candidates)
    frequency, fit = scan(coarse + np.arange(-50, 51)/50000)
    phase = np.arange(len(AGES))*frequency
    return float(frequency), fit[0]+fit[1]*np.cos(phase)+fit[2]*np.sin(phase)


def frequency_intervals(counts, remainder, q):
    """Conservative all-alias bands for a constant plus one line with bounded remainder.

    For each lag h, differences of a single sinusoid satisfy d2+d0 =
    2 cos(h*w) d1. Bound all four probabilities by the same absolute
    error, including the subleading poles. Every frequency in [0,pi]
    lies within half a cell of its mesh point. Enlarging each half-space
    by the Lipschitz variation h*half_cell preserves that entire interval.
    """
    radius = np.sqrt(14/SHOTS) + 2e-12
    # A union bound on 6144 odd-port frequencies gives 12288 exp(-28).
    budget = radius + q*remainder
    y = counts[:, -1]/SHOTS
    edges = np.linspace(0, np.pi, 65537)
    grid = (edges[:-1]+edges[1:])/2
    keep = np.ones(len(grid), bool)
    for h in (1, 3, 7, 11, 17, 23, 31, 39, 41):
        d = y[h:]-y[:-h]
        x = 2*d[h:-h]
        target = d[2*h:]+d[:-2*h]
        usable = x != 0
        if not np.any(usable):
            continue
        center = target[usable]/x[usable]
        width = 8*budget/abs(x[usable])
        mesh = h*(edges[1]-edges[0])/2 + 1e-12
        cosine = np.cos(h*grid)
        keep &= (cosine >= np.max(center-width)-mesh) & (cosine <= np.min(center+width)+mesh)
    boundaries = np.diff(np.r_[False, keep, False].astype(int))
    return [[float(edges[a]), float(edges[b])] for a, b in
            zip(np.flatnonzero(boundaries == 1), np.flatnonzero(boundaries == -1))]


def empirical_wash(ordinary, translated, q):
    """Post-acquisition contrast of temporal spreads at the fixed local port.

    Every odd-click frequency lies within sqrt(14/SHOTS) of its sampling
    probability except with family probability 2*6144*exp(-28). The source
    comparison adds the declared 2e-12 allowance per probability. This
    simultaneous event also covers selection of the largest and smallest
    probabilities at different ages in each cohort. Divide by the known
    native herald mass q only after comparing original-start probabilities.
    """
    check(ordinary.shape == translated.shape == (len(AGES), len(OUTCOMES)) and 0 < q <= 1,
          "both complete point-receiver cohorts required")
    epsilon = np.sqrt(14/SHOTS) + 2e-12
    ordinary_range = float(np.ptp(ordinary[:, -1]/SHOTS))
    washed_range = float(np.ptp(translated[:, -1]/SHOTS))
    lower = ordinary_range-washed_range-4*epsilon
    return {"ordinary_original_start_spread": ordinary_range, "translated_original_start_spread": washed_range,
            "conditional_difference_lower": float(lower/q),
            "conditional_translated_upper": float((washed_range+2*epsilon)/q),
            "positive_lower_bound": bool(lower > 0)}


def score(weights, counts, *, exact_weights=True):
    """Read held waveforms and all-alias intervals with gates fixed before sampling."""
    shape = (2, 4, 6, 128, 6)
    check(weights.shape == counts.shape == shape and weights.dtype == counts.dtype == np.int64,
          "unexpected receipt axes")
    check(np.all((weights >= 0) & (weights <= GRID)) and np.all(weights.sum(axis=-1) == GRID),
          "incomplete input laws")
    check(np.all((counts >= 0) & (counts <= SHOTS)) and np.all(counts.sum(axis=-1) == SHOTS)
          and not np.any((weights == 0) & (counts != 0)), "incomplete or unsupported counts")
    rows = []
    for i, size in enumerate(SIZES):
        e, _, lines = poles(size)
        prediction = float(e[1]-e[0])
        remainder = 2*abs(lines[1:]).sum() + 5e-12
        for j, far in enumerate(FAR_POWER):
            q = (1-far)**2/2
            for m, menu in enumerate(MENUS):
                c = counts[i,j,m]
                law = native_words(size, far, menu)
                if exact_weights:
                    check(np.array_equal(categorical_weights(law), weights[i,j,m]),
                          "integer acquisition weights differ on this backend")
                check(np.max(abs(law-weights[i,j,m]/GRID)) <= 2e-12,
                      "recorded sampling inputs differ from the native source")
                row = dict(size=size, far_power=far, menu=menu,
                           totals=[int(v) for v in c.sum(axis=0)])
                if menu in ("ordinary", "reflected", "competing_birth"):
                    check(np.max(abs(law[:,-1]-q*predicted_intensity(size))) < 3e-12,
                          "point spectrum differs from complete native probabilities")
                    frequency, forecast = prefix_fit(c)
                    hold_error = float(np.max(abs(c[64:,-1]/SHOTS-forecast[64:])))
                    bands = frequency_intervals(c, remainder, q)
                    unique = (len(bands) == 1 and bands[0][0] <= prediction <= bands[0][1]
                              and bands[0][1]-bands[0][0] <= .0048)
                    row.update(predicted=prediction, prefix=frequency, held_max_error=hold_error,
                               alias_intervals=bands, conditional_off_line_bound=float(remainder),
                               resolved=bool(unique and abs(frequency-prediction) <= .0005
                                             and hold_error <= .00045))
                elif menu == "translated":
                    row.update(conditional_odd_spread=float(np.ptp(law[:,-1])/q),
                               wash_gate=bool(np.ptp(law[:,-1]) <= q*.00025))
                    row["empirical_spread"] = empirical_wash(counts[i,j,MENUS.index("ordinary")], c, q)
                else:
                    row["no_birth"] = bool(np.all(c[:,2] == SHOTS))
                rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--acquire", type=Path, help="create a count archive at an unused path")
    parser.add_argument("--records", type=Path, default=Path(__file__).with_name("formed_clock_counts.npz"))
    parser.add_argument("--details", action="store_true")
    parser.add_argument("--seed-replay", action="store_true",
                        help="also require exact generated integer laws and seeded counts on this backend")
    args = parser.parse_args()
    if args.acquire:
        weights, counts = acquire()
        with args.acquire.open("xb") as target:
            np.savez_compressed(target, weights=weights, counts=counts)
        path = args.acquire
    else:
        path = args.records
        payload = path.read_bytes()
        check(hashlib.sha256(payload).hexdigest() == ACQUISITION_SHA256, "acquisition checksum mismatch")
        with np.load(BytesIO(payload), allow_pickle=False) as archive:
            check(len(archive.files) == 2 and set(archive.files) == {"counts", "weights"},
                  "unexpected acquisition members")
            weights, counts = archive["weights"], archive["counts"]
    rows = score(weights, counts, exact_weights=bool(args.acquire or args.seed_replay))
    if args.seed_replay:
        check(np.array_equal(draw_counts(weights), counts), "original seeded counts differ on this backend")
    result = {"sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
              "trials": int(counts.sum()),
              "passes_by_gate": {name: [sum(r[name] for r in rows if name in r),
                                        sum(name in r for r in rows)]
                                 for name in ("resolved", "wash_gate", "no_birth")}}
    empirical = [r["empirical_spread"] for r in rows if "empirical_spread" in r]
    result["post_acquisition_wash"] = {
        "positive_comparisons": sum(v["positive_lower_bound"] for v in empirical),
        "comparisons": len(empirical),
        "minimum_conditional_difference_lower": min(v["conditional_difference_lower"] for v in empirical),
        "simultaneous_failure_probability_at_most": float(12288*np.exp(-28))}
    if args.details:
        result["menus"] = rows
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
