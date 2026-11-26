"""Two material sign probes in 132- and 660-word native interacting sectors."""

import argparse
import json
from itertools import combinations, permutations
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import cell

HERE = Path(__file__).resolve().parent
STATIONS = (0, 1, 3)
DONOR, MARK = 5, 2
STOCK = tuple(i for i in range(12) if i not in (*STATIONS, DONOR))
MATERIAL = tuple(permutations(STOCK, 2))
ROUTES = tuple(permutations(STATIONS, 2))
TICKS, STEP, SHOTS = 12, 0.5, 2**20
GROUPS = ("near", "quiet", "vacancy_at_mark")


def stock_source():
    size = len(MATERIAL)
    at = {pair: i for i, pair in enumerate(MATERIAL)}
    matrix = np.zeros((size, size), dtype=np.int64)
    for a, b in combinations(STOCK, 2):
        change = {a: b, b: a}
        for i, pair in enumerate(MATERIAL):
            image = tuple(change.get(j, j) for j in pair)
            matrix[at[image], i] += 1
    eye = np.eye(size, dtype=np.int64)
    if not np.all((matrix - 20 * eye) @ (matrix - 14 * eye) @ (matrix - 12 * eye) == 32):
        raise ArithmeticError("Stock mixing polynomial failed")
    if (matrix != matrix.T).any() or (matrix.sum(axis=0) != 28).any():
        raise ArithmeticError("Stock transition is not a reversible probability schedule")
    launch = at[STOCK[0], STOCK[-1]]
    return np.linalg.matrix_power(matrix / 28, 64)[:, launch]


def words(vacancies, excesses):
    out = []
    for vs in combinations(range(12), vacancies):
        rest = sorted(set(range(12)) - set(vs))
        for es in combinations(rest, excesses):
            out.append(tuple(0 if i in vs else 2 if i in es else 1 for i in range(12)))
    return out


def native_generator(board, inventory, k):
    at = {w: j for j, w in enumerate(inventory)}
    h = np.zeros((len(inventory), len(inventory)), dtype=np.int64)
    for j, w in enumerate(inventory):
        diag = 2 * k * sum((6 - int(deg)) * v*v for deg, v in zip(board.degree, w, strict=True))
        for a, b in board.seams:
            if w[a] != w[b]:
                v = list(w)
                v[a], v[b] = v[b], v[a]
                h[at[tuple(v)], j] -= 1
                diag += 1 + 2 * k * (w[a] - w[b])**2
        h[j, j] = diag
    if not np.array_equal(h, h.T):
        raise ArithmeticError("The native word contact generator is not Hermitian")
    return h


def launch(item, sender, bit):
    v, e = item
    w = [1] * 12
    w[v], w[e], w[sender] = 0, 2, bit
    return tuple(w)


def first_hits(raw, target_returns):
    """Renewal decomposition at the first recorded target hit."""
    found = [None]
    for tick in range(1, TICKS + 1):
        value = raw[tick].copy()
        for prior in range(1, tick):
            value -= target_returns[tick-prior] @ found[prior]
        found.append(value)
    return found


def flight(board, k, inserted):
    basis = words(1 + inserted, 1)
    paired = words(1, 1 + inserted)
    by_word = {word: j for j, word in enumerate(basis)}
    conjugate_at = {word: j for j, word in enumerate(paired)}
    conjugate = [conjugate_at[tuple(2-n for n in w)] for w in basis]
    h = native_generator(board, basis, k)
    reverse = native_generator(board, paired, k)[np.ix_(conjugate, conjugate)]
    if not np.array_equal(reverse - h, 8 * k * inserted * np.eye(len(h), dtype=np.int64)):
        raise ArithmeticError("Two exact physical generators fail charge reversal up to a scalar")
    e, vec = eigh(h / 2)
    e -= e[0]
    raw_basis = np.array(basis)
    sites = {s: np.flatnonzero(raw_basis[:, s] == 0) for s in STATIONS}
    labels = {s: np.array([by_word[launch(m, s, 1 - inserted)] for m in MATERIAL]) for s in STATIONS}
    propagators = {(i, j): [None] for i in STATIONS for j in STATIONS}
    sources = {route: [None] for route in ROUTES}
    activity = {}
    for tick in range(1, TICKS + 1):
        z = np.exp(-1j * tick * STEP * e)
        for dest in STATIONS:
            left = vec[sites[dest]] * z
            for orig in STATIONS:
                propagators[dest, orig].append(left @ vec[sites[orig]].T)
                if dest != orig:
                    sources[orig, dest].append(left @ vec[labels[orig]].T)
        if not inserted and tick in (2, 4):
            inside = (vec[raw_basis[:, MARK] != 1] * z) @ vec[labels[STATIONS[0]]].T
            activity[str(tick * STEP)] = np.sum(abs(inside)**2, axis=0)
    routes = {}
    for s, r in ROUTES:
        outward = first_hits(sources[s, r], propagators[r, r])
        back = first_hits(propagators[s, r], propagators[s, s])
        result = np.zeros((len(MATERIAL), TICKS + 2))
        for a in range(1, TICKS):
            for b in range(1, TICKS + 1 - a):
                result[:, a+b] += np.sum(abs(back[b] @ outward[a])**2, axis=0)
        result[:, -1] = 1 - result.sum(axis=1)
        if result.min() < -2e-12:
            raise ArithmeticError("A completion or timeout has negative probability")
        routes[s, r] = result
    return routes, activity


def masks():
    v = np.array([a == MARK for a, _ in MATERIAL])
    e = np.array([b == MARK for _, b in MATERIAL])
    return (v, e, ~(v | e))


def summarize(raw):
    # The three original categories resolve the herald; charge-blind near pools two.
    rows = (raw[0] + raw[1], raw[2], raw[0])
    clock = np.minimum(np.arange(TICKS + 2), TICKS) * STEP
    return {label: {"mean": float(row @ clock / row.sum()), "return": float(1 - row[-1] / row.sum())}
            for label, row in zip(GROUPS, rows, strict=True)}


def count_intervals(bank):
    # Conditional on the included herald count, these copies have its conditional law.
    rows = (bank[0] + bank[1], bank[2], bank[0])
    clock = np.minimum(np.arange(TICKS + 2), TICKS) * STEP
    span = (TICKS - 2) * STEP
    log_budget = np.log(4 * 72 / 0.01)  # 6 routes, 2 signs, 2 arms, 3 heralds
    enclosures = {}
    for name, row in zip(GROUPS, rows, strict=True):
        n = int(row.sum())
        mean = float(row @ clock / n)
        delay = clock.max() - mean
        coarse = span * np.sqrt(log_budget / (2 * n))
        variance = span * (delay + coarse)
        radius = np.sqrt(2 * variance * log_budget / n) + 2 * span * log_budget / (3 * n)
        enclosures[name] = [mean - radius, mean + radius]
    return enclosures


def energy_meter(board, k, source, archive=None, save_counts=False):
    bonds = [edge for edge in board.seams if MARK in edge]
    levels = sorted({0, k, k + 1, 4*k, 4*k + 1})
    conditional = np.zeros((3, len(levels)))
    heralds = masks()
    for j, m in enumerate(MATERIAL):
        category = next(i for i, group in enumerate(heralds) if group[j])
        w = launch(m, STATIONS[0], 1)
        for a, b in bonds:
            d = abs(w[a] - w[b])
            for energy in (k*d*d, k*d*d + (d != 0)):
                conditional[category, levels.index(energy)] += source[j] / (2 * len(bonds))
    bank = (np.random.default_rng(812064).multinomial(SHOTS, conditional.ravel()).reshape(conditional.shape)
            if save_counts else archive["seam_energy"] if archive is not None else None)
    def means(table):
        rows = (table[0] + table[1], table[2], table[0])
        return {g: float(len(bonds) * (r @ levels) / r.sum()) for g, r in zip(GROUPS, rows, strict=True)}
    if bank is not None and (bank.sum() != SHOTS or (bank < 0).any() or bank.shape != conditional.shape):
        raise ArithmeticError("Incomplete local seam-energy cohort")
    if archive is not None:
        np.testing.assert_allclose(conditional, archive["seam_energy_rate"], rtol=0, atol=2e-12)
    return means(conditional), means(bank) if bank is not None else {}, bank, conditional


def cycle(data, herald):
    def mean(s, r):
        return data[f"{s}-{r}"]["vacancy"]["probe"][herald]["mean"]
    return (mean(0, 1) * mean(1, 3) * mean(3, 0)
            / (mean(1, 0) * mean(3, 1) * mean(0, 3)))


def cycle_interval(intervals, herald):
    forward = [intervals[k]["vacancy"]["probe"][herald] for k in ("0-1", "1-3", "3-0")]
    backward = [intervals[k]["vacancy"]["probe"][herald] for k in ("1-0", "3-1", "0-3")]
    return [float(np.prod([r[i] for r in forward]) / np.prod([r[1-i] for r in backward])) for i in (0, 1)]


def run(k, archive=None, save_counts=False):
    board = cell()
    source = stock_source()
    signed_masks = masks()
    permutation = [MATERIAL.index((e, v)) for v, e in MATERIAL]
    panels = {}
    for arm, bit in (("sham", False), ("probe", True)):
        panels[arm], active = flight(board, k, bit)
        if not bit:
            activity = active
    rng = np.random.default_rng(812063)
    full, read, counts, intervals = {}, {}, {}, {}
    largest_even_tv = 0
    min_prob = [sum(source[m]) for m in signed_masks]
    for s, r in ROUTES:
        key = f"{s}-{r}"
        full[key], read[key], intervals[key] = {}, {}, {}
        signed = {}
        for q in ("vacancy", "excess"):
            full[key][q], read[key][q], intervals[key][q] = {}, {}, {}
            for arm, panel in panels.items():
                words_panel = panel[s, r] if q == "vacancy" else panel[s, r][permutation]
                table = np.stack([source[m] @ words_panel[m] for m in signed_masks])
                signed[q, arm] = table
                full[key][q][arm] = summarize(table)
                if save_counts or archive is not None:
                    name = f"{s}_{r}_{q}_{arm}"
                    bank = (rng.multinomial(SHOTS, table.ravel()).reshape(table.shape)
                            if archive is None else archive[name])
                    if bank.shape != table.shape or bank.dtype.kind not in "iu" or bank.sum() != SHOTS or bank.min() < 0:
                        raise ArithmeticError("Bad original cohort for the material herald")
                    if np.max(abs(bank / SHOTS - table)) > np.sqrt(np.log(2 * 1008 / .01) / (2 * SHOTS)):
                        raise ArithmeticError("A joint count differs from its free-propagation probability")
                    if archive is not None:
                        np.testing.assert_allclose(table, archive[name + "_rate"], rtol=0, atol=2e-12)
                    counts[name] = bank
                    counts[name + "_rate"] = table
                    read[key][q][arm] = summarize(bank)
                    intervals[key][q][arm] = count_intervals(bank)
        for arm in panels:
            v, e = signed["vacancy", arm], signed["excess", arm]
            for idx in ([0, 1], [2]):
                vout, eout = v[idx].sum(axis=0), e[idx].sum(axis=0)
                largest_even_tv = max(largest_even_tv, .5 * abs(vout/vout.sum() - eout/eout.sum()).sum())
    energy, read_energy, energy_counts, energy_rates = energy_meter(board, k, source, archive, save_counts)
    if energy_counts is not None:
        counts["seam_energy"] = energy_counts
        counts["seam_energy_rate"] = energy_rates
    selected = [signed_masks[0] | signed_masks[1], signed_masks[2], signed_masks[0]]
    def conditional(values):
        return {g: float(source[m] @ np.array(values)[m] / source[m].sum()) for g, m in zip(GROUPS, selected, strict=True)}
    return {
        "strength": k, "source_TV_from_uniform": .5 * float(abs(source - 1/len(source)).sum()),
        "source_TV_bound": np.sqrt(55) / 2 * (5/7)**64,
        "probabilities_vacancy_excess_unit_at_mark": min_prob,
        "local_seam_energy": energy, "measured_local_seam_energy": read_energy,
        "unmonitored_nonunit_at_mark": {t: conditional(v) for t, v in activity.items()},
        "maximum_two_charge_history_TV_on_sign_blind_heralds": largest_even_tv,
        "clock": full, "clock_products_forward_over_reverse": {g: cycle(full, g) for g in GROUPS},
        "clock_product_intervals_99_percent": {g: cycle_interval(intervals, g) for g in GROUPS} if counts else {},
        "count_clock": read if counts else {},
        "clock_mean_intervals_simultaneous_99_percent": intervals if counts else {},
    }, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    path = HERE / "counts.npz"
    if args.write:
        native, counts = run(int(candidate.KAPPA), save_counts=True)
        np.savez_compressed(path, **counts)
    else:
        with np.load(path, allow_pickle=False) as bank:
            native, _ = run(int(candidate.KAPPA), archive=bank)
    control, _ = run(0)
    report = {"native": native, "exchange_only": control, "shots_per_original_cohort": SHOTS}
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    else:
        saved = json.loads(path.read_text())
        def compare(a, b):
            if isinstance(a, dict):
                if a.keys() != b.keys():
                    raise ArithmeticError("Cohort or summary fields changed")
                for name, value in a.items():
                    compare(value, b[name])
            else:
                np.testing.assert_allclose(a, b, rtol=2e-10, atol=2e-11)
        compare(report, saved)
    print("Verified native paired radar histories on the equilibrated complete cell.")


if __name__ == "__main__":
    main()
