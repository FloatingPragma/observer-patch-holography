"""Native marked-face sampling and its port-weighted compact-family read."""

import argparse
from collections import deque
from functools import cache
from itertools import combinations
import json
from math import asin, log, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from model.candidate_laws.interval_echo.candidate import KAPPA, native_law as law
from model.executable_model.carrier import CELL_SEAMS, federation, isomorphism

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SHOTS = 2**22
N = 12
EDGES = set(CELL_SEAMS)
WORDS = [
    (d, a, b) for d in range(N) for a, b in combinations(range(N), 2) if d not in (a, b)
]
INDEX = {w: i for i, w in enumerate(WORDS)}
FACES = [f for f in combinations(range(N), 3) if set(combinations(f, 2)) <= EDGES]
SEED = 130910
CLASSIFIER = ROOT / "code/b14_jacobi"


def path():
    def extend(w):
        if len(w) == N:
            return w
        for v in range(N):
            if v not in w and tuple(sorted((w[-1], v))) in EDGES:
                found = extend(w + (v,))
                if found:
                    return found

    return extend((0,))


def routes(level):
    car = federation(level)
    gates = [(a, b) if a < 12 else (b, a) for a, b in car.external if min(a, b) < 12][
        :2
    ]
    if len(gates) != 2:
        raise ArithmeticError("Two distinct external references required")
    rows = []
    for gate, reference in gates:
        paths = {gate: (gate,)}
        left = deque([gate])
        while left:
            i = left.popleft()
            for v in sorted(j for j in car.neighbours[i] if j < 12):
                if v not in paths:
                    paths[v] = paths[i] + (v,)
                    left.append(v)
        entries = []
        for v in range(N):
            route = paths[v][::-1]
            edges = list(zip(route, route[1:]))
            swaps = edges + list(reversed(edges[:-1]))
            labels = list(range(N))
            for i, j in swaps:
                labels[i], labels[j] = labels[j], labels[i]
            expected = list(range(N))
            expected[gate], expected[v] = expected[v], expected[gate]
            if labels != expected:
                raise ArithmeticError("A transit register changed")
            entries.append(route)
        rows.append(dict(reference=reference, owner=reference // 12, paths=entries))
    return rows


def meter():
    echo = law.echo()
    bright = np.tile([v == (1, 1) for v in law.SRC], law.DIM // 9)
    tags = []
    for value in range(3):
        reads = []
        for ref in (0, 2):
            i = law.SRC.index((value, ref))
            amplitude = echo[:, i]
            expected = value + ref == 2
            selected = amplitude * (bright == expected)
            np.testing.assert_allclose(
                echo.conj().T @ selected, np.eye(law.DIM)[:, i], atol=1e-13
            )
            reads.append(expected)
        tags.append(2 if reads[0] else 0 if reads[1] else 1)
    return tags


def read_word(word, reference_routes, tags):
    """Read each source address through two distinct guarded external seams."""
    source = [0] * N
    d, a, b = word
    source[d], source[a], source[b] = 2, 1, 1
    values = []
    for port in range(N):
        for ref, routed in zip((0, 2), reference_routes, strict=True):
            path_to_gate = routed["paths"][port]
            forward = list(zip(path_to_gate, path_to_gate[1:]))
            exchanges = forward + forward[:-1][::-1]
            for u, v in exchanges:
                source[u], source[v] = source[v], source[u]
            level = source[path_to_gate[-1]]
            bright = tags[level] + ref == 2
            for u, v in reversed(exchanges):
                source[u], source[v] = source[v], source[u]
            if bright:
                break
        values.append(2 - ref if bright else 1)
    return (values.index(2), *(i for i, x in enumerate(values) if x == 1))


@cache
def swap_indices(u, v):
    def exchange(i):
        return v if i == u else u if i == v else i

    return np.array(
        [INDEX[(exchange(d), *sorted([exchange(a), exchange(b)]))] for d, a, b in WORDS]
    )


@cache
def rates(num, den, speed):
    e, c = law.local_term_parts()
    full = speed * e + KAPPA * c
    h = full[:9, :9]
    np.testing.assert_allclose(full[9:, :9], 0, atol=2e-13)
    tau = 2 * asin(sqrt(num / den))
    pulse = expm(-1j * tau * h)
    swap = np.eye(9)[:, [law.SRC.index((b, a)) for a, b in law.SRC]]
    np.testing.assert_allclose(expm(-1j * pi * h), swap, atol=2e-13)
    expected = np.sin(speed * tau / 2) ** 2
    values = []
    for i, (a, b) in enumerate(law.SRC):
        j = law.SRC.index((b, a))
        sq = abs(pulse[:, i]) ** 2
        pred = np.zeros(9)
        pred[i] += 1 - expected
        pred[j] += expected
        np.testing.assert_allclose(sq, pred, atol=2e-13)
        if a != b:
            values.append(sq[j])
    return float(np.mean(values))


def evolve(seed, offers, speed, integer=False):
    f = np.zeros(len(WORDS), dtype=object if integer else float)
    f[INDEX[seed]] = 1
    for u, v, num, den in offers:
        p = rates(num, den, speed)
        if integer:
            if speed != 1:
                raise ValueError("Integer insertion proof belongs to native exchange")
            f = (den - num) * f + num * f[swap_indices(u, v)]
        else:
            f = (1 - p) * f + p * f[swap_indices(u, v)]
    if integer:
        if f[0] <= 0 or any(value != f[0] for value in f):
            raise ArithmeticError("The complete exact terminal measure is not uniform")
        return str(f[0])
    return f / f.sum()


def calibration_inputs():
    # These two sparse coefficient tables and their coordinate order are the
    # same family parametrizations used by the exact classifier.
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "oph_oriented_selector", CLASSIFIER / "oriented_face_bracket_selector.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    manifest, raw, _, system, _ = mod.load_inputs()
    tables = np.zeros((mod.M, len(mod.COORDS)))
    coord = {key: i for i, key in enumerate(mod.COORDS)}
    for i, r in enumerate(raw["basis"]):
        for o, a, b, n, d in r["entries"]:
            tables[i, coord[o, a, b]] = n / d
    inverse = np.zeros((mod.M, mod.M))
    for i, pieces in enumerate(
        system["fixed_line_reduction"]["channel_decomposition"][
            "inverse_transform_rows"
        ]
    ):
        for j, an, ad, bn, bd in pieces:
            inverse[i, j] = an / ad + sqrt(5) * bn / bd
    channel = inverse.T @ tables
    class_specs = {
        "P": [("t_pp_to_p",), ("t_mm_to_m",)],
        "F": [("t_pp_to_p", "t_pf_to_f", -sqrt(5)), ("t_ff_to_p",), ("t_mm_to_m",)],
        "G": [("t_mm_to_m", "t_mf_to_f", sqrt(5)), ("t_ff_to_m",), ("t_pp_to_p",)],
    }
    bases = {
        name: np.array(
            [
                channel[mod.CH[items[0]]]
                + (items[2] * channel[mod.CH[items[1]]] if len(items) == 3 else 0)
                for items in cols
            ]
        ).T
        for name, cols in class_specs.items()
    }
    pin = manifest["carrier"]
    source_edges = [tuple(pin["ports"].index(p) for p in e) for e in pin["edges"]]
    to_pin = isomorphism(list(EDGES), source_edges)
    return mod, manifest, bases, to_pin


def classify(weights, coords, target, bases):
    o, i, j = np.array(coords).T
    factors = np.sqrt(weights[o] / (weights[i] * weights[j]))
    fits = {}
    for name, matrix in bases.items():
        coeff, _, rank, _ = np.linalg.lstsq(
            matrix * factors[:, None], target * factors, rcond=None
        )
        if rank != matrix.shape[1]:
            raise ArithmeticError("Unidentified family coordinates")
        fits[name] = dict(
            coefficients=coeff.tolist(),
            squared_distance=float(np.sum((factors * (target - matrix @ coeff)) ** 2)),
        )
    return fits


def measure():
    rng = np.random.default_rng(SEED)
    p = path()
    banks = []
    meter_tags = meter()
    rows = [routes(level) for level in range(3)]
    readbacks = [[INDEX[read_word(w, r, meter_tags)] for w in WORDS] for r in rows]
    if any(table != list(range(len(WORDS))) for table in readbacks):
        raise ArithmeticError("Terminal instrument changed the word")
    for speed in (1, 3):
        for seed in [
            (FACES[0][0], *FACES[0][1:]),
            next(w for w in WORDS if tuple(sorted(w)) not in FACES),
        ]:
            for reverse in (False, True):
                chain = p[::-1] if reverse else p
                ops = [
                    (chain[i], chain[i + 1], i + 1, i + 2)
                    for j in range(1, N)
                    for i in reversed(range(j))
                ]
                prob = evolve(seed, ops, speed)
                exact = evolve(seed, ops, 1, integer=True) if speed == 1 else None
                counts = rng.multinomial(SHOTS, prob)
                banks.append(
                    dict(
                        speed=speed,
                        seed=seed,
                        reverse=reverse,
                        counts=counts.tolist(),
                        max_uniform_error=float(max(abs(prob - 1 / len(prob)))),
                        numerator_per_word=exact,
                    )
                )
    # These independently seeded launches run the three deployed read maps.
    native_prob = evolve(
        (0, 1, 2),
        [
            (p[i], p[i + 1], i + 1, i + 2)
            for j in range(1, N)
            for i in reversed(range(j))
        ],
        1,
    )
    host_counts = []
    for readback in readbacks:
        emitted = rng.multinomial(SHOTS, native_prob)
        recorded = np.zeros(len(WORDS), dtype=int)
        np.add.at(recorded, readback, emitted)
        host_counts.append(recorded.tolist())
    uniform = np.ones(len(WORDS)) / len(WORDS)
    passive = []
    half = [(u, v, 1, 2) for u, v in sorted(EDGES)]
    for sweeps in (1, 4, 16):
        law_row = evolve(WORDS[0], half * sweeps, 1)
        counts = rng.multinomial(SHOTS, law_row)
        passive.append(
            dict(
                sweeps=sweeps,
                TV=float(abs(law_row - uniform).sum() / 2),
                counts=counts.tolist(),
            )
        )
    mod, manifest, bases, to_pin = calibration_inputs()
    data = np.array(banks[0]["counts"])
    flags = np.array([tuple(sorted(w)) in FACES for w in WORDS])
    pcounts = np.bincount(
        [to_pin[w[0]] for w, b in zip(WORDS, flags) if b],
        weights=data[flags],
        minlength=N,
    )
    weights = pcounts / pcounts.sum()
    # Handedness is a convention. The measured source supplies the support
    # for the alternating face tensor used in the registered classifier.
    seen = {
        tuple(sorted(to_pin[v] for v in w))
        for w, c, b in zip(WORDS, data, flags)
        if c and b
    }
    oriented = [f for f in mod.canonical_faces(manifest) if tuple(sorted(f)) in seen]
    T = mod.face_tensor(oriented)
    target = np.array([T[o][i][j] for o, i, j in mod.COORDS])
    fitted = classify(weights, mod.COORDS, target, bases)
    reference = classify(np.full(N, 1 / N), mod.COORDS, target, bases)
    f_control = bases["F"] @ fitted["F"]["coefficients"]
    foreign = classify(weights, mod.COORDS, f_control, bases)
    allowance = sqrt(log(2 * N / 0.01) / (2 * pcounts.sum()))
    lo, hi = min(weights) - allowance, max(weights) + allowance
    intervals = {
        k: [
            val["squared_distance"] * lo / (N * hi**2),
            val["squared_distance"] * hi / (N * lo**2),
        ]
        for k, val in reference.items()
    }
    return dict(
        shots=SHOTS,
        path=p,
        seed=SEED,
        words=WORDS,
        routes=rows,
        ternary_tags=meter_tags,
        native_and_controls=banks,
        host_counts=host_counts,
        half_contacts=passive,
        material_to_pinned=to_pin,
        classification=dict(
            admissions=int(pcounts.sum()),
            port_counts=list(map(int, pcounts)),
            fitted=fitted,
            reference=reference,
            score_envelopes=intervals,
            f_control=foreign,
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    write = parser.parse_args().write
    result = measure()
    content = json.dumps(result, separators=(",", ":")) + "\n"
    file = HERE / "records.json"
    if write:
        file.write_text(content)
    else:
        archived = json.loads(file.read_text())
        for key in (
            "routes",
            "words",
            "ternary_tags",
            "path",
            "seed",
            "shots",
            "material_to_pinned",
            "host_counts",
        ):
            if json.dumps(archived[key], sort_keys=True) != json.dumps(
                result[key], sort_keys=True
            ):
                raise SystemExit(f"Changed source or receipt: {key}")
        for key in ("native_and_controls", "half_contacts"):
            for before, after in zip(archived[key], result[key], strict=True):
                for field in before:
                    if field in ("TV", "max_uniform_error"):
                        np.testing.assert_allclose(
                            before[field], after[field], atol=1e-13
                        )
                    elif json.dumps(before[field]) != json.dumps(after[field]):
                        raise SystemExit(f"Changed record: {key}.{field}")
        for key in ("admissions", "port_counts"):
            if archived["classification"][key] != result["classification"][key]:
                raise SystemExit("Observer bins differ")
        for group in ("fitted", "reference", "f_control"):
            for name, fit in result["classification"][group].items():
                old = archived["classification"][group][name]
                for field in fit:
                    np.testing.assert_allclose(
                        fit[field], old[field], atol=1e-10, rtol=1e-11
                    )
        for name, bounds in result["classification"]["score_envelopes"].items():
            np.testing.assert_allclose(
                bounds, archived["classification"]["score_envelopes"][name], atol=1e-10
            )
    c = result["classification"]
    print(
        c["admissions"],
        {k: v["squared_distance"] for k, v in c["fitted"].items()},
        c["score_envelopes"],
    )
