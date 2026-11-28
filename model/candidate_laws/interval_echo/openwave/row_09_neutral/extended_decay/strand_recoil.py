"""Read residual strand and spectator recoil together after the real heavy-source flight."""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
from scipy.special import jv
from . import strand_exterior as exterior, weak_strand as ws
from .measure import HERE
from ..passive_decay.measure import compare_records, hamiltonians

WAITS = (0, 4, 8, 16, 32)
DESIGN = "STRAND_RECOIL.md"


def global_recoil(source_gaps):
    target = source_gaps[2] - source_gaps[1]
    result = []
    for level in range(3):
        b = exterior.federation(level)
        a, z = np.asarray(b.seams).T
        ports = np.arange(b.n)
        h = sparse.coo_matrix(
            (
                np.r_[0.5 * b.degree, np.full(2 * len(a), -0.5)],
                (np.r_[ports, a, z], np.r_[ports, z, a]),
            ),
            shape=(b.n, b.n),
        ).tocsc()
        # Fix the initial solver vector independent of the solution. The
        # returned residual also bounds distance to an exact physical eigenline.
        trial = np.sin(ports + 1.0)
        values, v = eigsh(h, k=4, sigma=target, which="LM", v0=trial, tol=2e-13)
        i = np.argmin(abs(values - target))
        state = v[:, i]
        ev = values[i]
        error = np.linalg.norm(h @ state - ev * state)
        assert error < 2e-11
        if level == 0:
            all_values = eigh(h.toarray(), eigvals_only=True)
            assert abs(np.min(abs(all_values - target)) - abs(ev - target)) < 1e-12
        result.append(
            dict(
                owners=b.cells,
                ports=b.n,
                target=target,
                nearest=ev,
                mismatch=ev - target,
                residual=error,
                per_port_density_transition_norm=float(
                    np.linalg.norm(state / np.sqrt(b.n))
                ),
            )
        )
    return result


def embedded_rows(words, terminal, cell_words):
    lookup = {frozenset(w): i for i, w in enumerate(words)}
    return np.array(
        [
            [lookup[frozenset((*pair, terminal[1 - e]))] for pair in cell_words]
            for e in (0, 1)
        ]
    )


def inject(medium, cold, rows, m):
    side = medium.shape[0] // 2
    wave = np.zeros((side, m, medium.shape[1]), dtype=np.result_type(medium, cold))
    for e in (0, 1):
        wave[:, rows[e], :] = (
            medium.reshape(side, 2, -1)[:, e, None, :] * cold[None, :, None]
        )
    return wave.reshape(-1, medium.shape[1])


def severed_flight(
    board, dimers, cell, path, interior, cell_words, cold, words, known, gaps
):
    # With the terminal link severed the strand is an independent system.
    # Only the original two exterior units can reach the still-open approach.
    sites = tuple(sorted(set(path[1:]) | set(cell)))
    outside, small = exterior.outlet.binary_generator(
        board, sites, [*interior, *zip(path[1:], path[2:])], 2
    )
    initial = np.zeros(len(outside))
    index = {w: i for i, w in enumerate(outside)}
    initial[[index[w] for w in cell_words]] = cold
    ev, v = eigh(small)
    amp = v.T @ initial
    rows = embedded_rows(words, dimers[-1], outside)

    def at(t):
        strand = known * np.exp(-1j * t * gaps)
        stock = v @ (np.exp(-1j * t * ev) * amp)
        return inject(strand, stock, rows, len(words))

    return at


def final_joint(wave, rows, strand_frame, cell_frame, ngroup, cgroup):
    # Only words with exactly one terminal unit and exactly the two cell
    # residents enter this double-energy measurement. All other words refuse.
    n = len(strand_frame)
    prepared = wave.reshape(n // 2, -1, wave.shape[-1])[:, rows, :].reshape(
        n, len(cell_frame), -1
    )
    after = np.einsum(
        "ai,awc,wb->ibc", strand_frame, prepared, cell_frame, optimize=True
    )
    results = np.empty(
        (len(ngroup), len(cgroup), wave.shape[-1], wave.shape[-1]), complex
    )
    for i, sg in enumerate(ngroup):
        for j, cg in enumerate(cgroup):
            block = after[np.ix_(sg, cg)].reshape(-1, wave.shape[-1])
            results[i, j] = block.conj().T @ block
    complete = np.einsum("nbc,nbd->cd", prepared.conj(), prepared, optimize=True)
    np.testing.assert_allclose(results.sum(axis=(0, 1)), complete, atol=2e-12)
    return results


def open_flights(h, initial, times):
    """Full-state Chebyshev expansion at all requested clocks on one recurrence."""
    diagonal = h.diagonal()
    radius_by_row = np.asarray(abs(h).sum(axis=1)).ravel() - abs(diagonal)
    lo = min(diagonal - radius_by_row)
    hi = max(diagonal + radius_by_row)
    centre = (hi + lo) / 2
    radius = (hi - lo) / 2
    if radius == 0:
        return [np.exp(-1j * t * centre) * initial for t in times]
    normalized = (h - centre * sparse.eye(h.shape[0], format="csr")) / radius
    coefficients = []
    for t in times:
        x = t * radius
        if x == 0:
            coefficients.append(np.ones(1))
            continue
        k = int(np.ceil(x)) + 1
        while True:
            # Bessel's generating function on |z|=exp(a) bounds the whole
            # omitted tail by 2 exp(x sinh(a)-k a)/(1-exp(-a)).
            a = np.arccosh(k / x)
            if np.log(2) + x * np.sinh(a) - k * a - np.log(-np.expm1(-a)) < np.log(
                1e-13
            ):
                break
            k += 1
        weights = 2 * jv(np.arange(k), x)
        weights[0] /= 2
        coefficients.append(weights)
    results = [weights[0] * initial.astype(complex) for weights in coefficients]
    previous = initial
    current = normalized @ initial
    twice = normalized * 2
    for k in range(1, max(map(len, coefficients))):
        sign = (-1.0) ** ((k + 1) // 2)
        for weights, output in zip(coefficients, results, strict=True):
            if k < len(weights):
                channel = output.imag if k % 2 else output.real
                channel += sign * weights[k] * current
        previous, current = current, twice @ current - previous
    return [
        np.exp(-1j * t * centre) * out for t, out in zip(times, results, strict=True)
    ]


def groups(gaps, budget):
    indices = np.flatnonzero(gaps <= budget)
    starts = np.r_[
        0, np.flatnonzero(abs(np.diff(gaps[indices])) > 1e-9) + 1, len(indices)
    ]
    resolved = [indices[a:z] for a, z in zip(starts[:-1], starts[1:], strict=True)]
    return [*resolved, np.flatnonzero(gaps > budget)]


def neutral_orbit(cell, internal, cell_words):
    # Every binary site's private phase is read directly or by a reversible
    # in-cell swap to a degree-five standard. These signatures separate words.
    words = {frozenset(w) for w in cell_words}
    assert len({tuple(int(p in w) for p in cell) for w in words}) == len(words)
    first = next(iter(words))
    seen = {first}
    active = [first]
    while active:
        w = active.pop()
        for p, q in internal:
            if (p in w) != (q in w):
                nxt = w ^ {p, q}
                if nxt not in seen:
                    assert nxt in words
                    seen.add(nxt)
                    active.append(nxt)
    assert seen == words
    return len(seen)


def summarize(joint, rho):
    # j=0 is the independently commissioned cold-strand control at exactly
    # the same admitted launch rate. The remaining columns carry actual modes.
    retained = np.einsum("ijab,ba->ij", joint[:, :, 1:, 1:], rho).real
    reset = joint[:, :, 0, 0].real * np.trace(rho).real
    return dict(
        keep=retained.tolist(),
        reset=reset.tolist(),
        inventory_refused_keep=float(np.trace(rho).real - retained.sum()),
        inventory_refused_reset=float(np.trace(rho).real - reset.sum()),
    )


def one(n):
    b, d, heads, _ = ws.ring_layout(n, True, False)
    path, cell, interior = exterior.geometry(b, d, heads)
    sectors = exterior.outlet.receiving_cell(b, cell, interior)
    cell_words, energy, cv = sectors[2]
    _, ordinary = exterior.outlet.binary_generator(b, cell, interior, 2)
    assert np.linalg.norm(ordinary @ cv - cv * energy, ord=2) < 1e-10
    e, sv = eigh(ws.strand(n).toarray())
    core = eigh(hamiltonians(4)[0], eigvals_only=True)
    budget = core[-1] - core[2]
    sg, cg = e - e[0], energy - energy[0]
    affordable = (sg > 1e-10) & (sg <= budget)
    sr = groups(sg, budget)
    cr = groups(cg, budget)
    # The published binary pointer routes each rail to an actual private
    # standard. This independent charge comparison resolves both signs relative
    # to two resident units without exposing individual port addresses.
    _, pointer = exterior.outlet.relative_inventory_reader(
        b,
        cell,
        {count: (list(combinations(cell, count)), None, None) for count in (1, 2, 3)},
        8,
        2,
    )
    assert pointer == [7, 0, 1]
    assert neutral_orbit(cell, interior, cell_words) == 66
    standard = dict(
        strand_low=[sg[a[0]] for a in sr[:-1]],
        cell_ground=energy[0],
        cell_low=[[cg[a[0]], len(a)] for a in cr[:-1]],
        first_higher_cell=cg[cr[-1][0]],
        destination_pointers_for_one_two_three=pointer,
    )
    if len(cr) == 2:
        assert min(cg[1:]) > budget
        return dict(
            standards=standard,
            reason="no positive two-resident level inside the complete incoming budget",
        )
    known = np.column_stack((sv[:, 0], sv[:, affordable]))
    stocks = exterior.source_stocks(b, d, heads, sv, affordable)
    source_probs = {k: np.diag(rho).real.tolist() for k, rho in stocks.items()}
    all_matrices = {}
    for setting in ("open", "link_cut", "parked"):
        words, h = exterior.mobile_generator(
            b, d, heads, path, cell, interior, setting, number=3
        )
        all_matrices[setting] = h
    rows = embedded_rows(words, d[-1], cell_words)
    initial = inject(known, cv[:, 0], rows, len(words))
    isolated_flight = severed_flight(
        b,
        d,
        cell,
        path,
        interior,
        cell_words,
        cv[:, 0],
        words,
        known,
        sg[: len(sr) - 1],
    )
    parked = all_matrices["parked"]
    initial_h = initial.T @ (parked @ initial)
    np.testing.assert_allclose(
        initial_h - initial_h[0, 0] * np.eye(len(initial_h)),
        np.diag(sg[: len(initial_h)]),
        atol=1e-10,
    )
    np.testing.assert_allclose(parked @ initial, initial @ initial_h, atol=1e-11)
    trials = {}
    for key in ("open", "link_cut"):
        h = all_matrices[key]
        exchange = h - parked
        start = initial.T @ (exchange @ initial)
        pages = []
        flights = (
            (isolated_flight(t) for t in WAITS)
            if key == "link_cut"
            else open_flights(h, initial, WAITS)
        )
        for t, wave in zip(WAITS, flights, strict=True):
            np.testing.assert_allclose(
                wave.conj().T @ wave, np.eye(wave.shape[1]), atol=2e-10
            )
            joint = final_joint(wave, rows, sv, cv, sr, cr)
            stop = -(wave.conj().T @ (exchange @ wave))
            net = wave.conj().T @ (parked @ wave) - initial_h
            np.testing.assert_allclose(start + stop, net, atol=3e-9)
            resolved = np.diagonal(joint, axis1=2, axis2=3).real
            # The per-level input comparison is independent of any heavy-source
            # normalization; multiply its columns by the recorded original weights.
            original = {}
            for label, rho in stocks.items():
                rec = summarize(joint, rho)
                target = joint[1, 1]
                rec["incoming_second_residual_first_cell_doublet"] = float(
                    rho[1, 1].real * target[2, 2].real
                )
                rec["incoming_second_reset_same_joint"] = float(
                    rho[1, 1].real * target[0, 0].real
                )
                rec["mean_net_work_per_original"] = float(
                    np.trace(net[1:, 1:] @ rho).real
                )
                original[label] = rec
            pages.append(
                dict(
                    wait=t,
                    conditional_by_incoming_level=resolved.tolist(),
                    open_work=np.diag(start).real.tolist(),
                    close_work=np.diag(stop).real.tolist(),
                    net_work=np.diag(net).real.tolist(),
                    source=original,
                )
            )
        trials[key] = pages
    return dict(
        standards=standard,
        path=path,
        dimension=initial.shape[0],
        original_mode_probabilities=source_probs,
        far_port_transition=float(
            abs((1 - 2 * ws.fields(n)[0][:, -1]) * sv[:, 1] @ sv[:, 2])
        ),
        bare_board_separate_standard=global_recoil(sg),
        trials=trials,
    )


def acquire():
    return dict(
        protocol_sha256=sha256((HERE / DESIGN).read_bytes()).hexdigest(),
        tested={str(n): one(n) for n in (6, 8)},
    )


def calibrate_receiver():
    """Read the known doublet from a cold two-load cell using a local phase."""
    b, d, h, _ = ws.ring_layout(8, True, False)
    path, cell, edges = exterior.geometry(b, d, h)
    words, energy, v = exterior.outlet.receiving_cell(b, cell, edges)[2]
    gaps = energy - energy[0]
    lines = np.unique(np.round(gaps, 10))
    low, high = lines[1] / 2, (lines[1] + lines[2]) / 2
    clock = exterior.outlet.clock(energy, energy[0], high, labels=32768, low_edge=low)[
        :, 1
    ]
    levels = np.isclose(gaps, lines[1], atol=1e-9)
    # Uniformly sample all 66 sharp two-resident words. Three physical,
    # non-demolition zero-bin comparisons prepare a cold receiver before use.
    # No input word or positive-energy refusal is omitted in this rate.
    ratio = gaps / exterior.outlet.PERIOD
    zero = (np.sinc(32768 * ratio) / np.sinc(ratio)) ** 2
    dirty = float((zero[1:] ** 3).sum())
    result = {}
    for port in cell:
        phase = 1 - 2 * np.array([int(port in w) for w in words])
        amplitudes = v.T @ (phase * v[:, 0])
        result[str(port)] = dict(
            true=float(np.dot(amplitudes[levels], amplitudes[levels])),
            read=float(amplitudes**2 @ clock),
        )
    return dict(
        entrance=int(path[-1]),
        cold_preparation_admission=(1 + dirty) / len(words),
        cold_preparation_non_ground_upper=dirty / (1 + dirty),
        lines=[low, high],
        true_clock_min=float(clock[levels].min()),
        false_clock_max=float(clock[~levels].max()),
        three_read_false_max=float(clock[~levels].max() ** 3),
        three_read_true_min=float(clock[levels].min() ** 3),
        cold_false=float(clock[0]),
        site_phase=result,
    )


def decode(data):
    for n, rec in data["tested"].items():
        print(n, rec["standards"])
        for key, pages in rec.get("trials", {}).items():
            for row in pages:
                cells = np.array(row["conditional_by_incoming_level"])
                pulse = row["source"]["native_pulse"]
                print(
                    key,
                    row["wait"],
                    "conditional cold/1/2/3, tagged original, reset, coherent, work2",
                    np.round(
                        [
                            *cells[1, 1],
                            pulse["incoming_second_residual_first_cell_doublet"],
                            pulse["incoming_second_reset_same_joint"],
                            pulse["keep"][1][1],
                            row["net_work"][2],
                        ],
                        9,
                    ),
                )
        if "bare_board_separate_standard" in rec:
            print("independent bare board", rec["bare_board_separate_standard"])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquire", type=Path)
    p.add_argument("--receiver-only", action="store_true")
    args = p.parse_args()
    if args.receiver_only:
        print(json.dumps(calibrate_receiver(), indent=2))
        raise SystemExit(0)
    data = acquire()
    decode(data)
    if args.acquire:
        with args.acquire.open("x") as f:
            json.dump(data, f, allow_nan=False, separators=(",", ":"))
            f.write("\n")
    else:
        compare_records(json.loads((HERE / "strand_recoil.json").read_text()), data)
