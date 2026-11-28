"""Compare the two-clock overlap records with all static mobile-load profiles."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import linprog

from . import joint_records as joint

HERE = joint.HERE
TIMES = (1, 2, 3, 4, 5, 6, 8, 12)
SHOTS = 2**30
SEED = 193241


def clock_laws(stiffness):
    board, seam, _, _ = joint.source.matched.fp.make_board()
    nodes = list(range(board.n))
    states, bright, dark, other = joint.source.preparation()
    w = (4 * 0.68 - 1) / 3
    v = np.sqrt(w) * dark + 1j * np.sqrt(1 - w) * other
    sources = {
        "paired": states["P"],
        "matched": states["FG0"],
        "unmatched": states["FG4"],
        "interior": np.outer(bright, bright) / 4 + 0.75 * np.outer(v, v.conj()),
    }
    read = {k: np.zeros((len(TIMES), 4)) for k in ("empty", *nodes, *sources)}
    for n in (0, 1):
        terms, starts, _ = joint.flight.parts(n, stiffness)
        inputs = starts[None]
        systems = {}
        for par in range(8):
            h4 = terms[0].copy()
            for k in range(3):
                if par & (1 << k):
                    h4 = h4 + terms[k + 1]
            h4 = h4.toarray()
            if stiffness == 8:
                np.testing.assert_array_equal(h4, joint.source.timing.model(n, par)[0])
            systems[par] = eigh(h4 / 4)
        for ti, t in enumerate(TIMES):
            U = {p: (V * np.exp(-1j * t * E)) @ V.T for p, (E, V) in systems.items()}
            for src_par in (0, 1):
                for near, far in np.ndindex(2, 2):
                    # Coherently read both sign bits; the source sign is traced.
                    block = (
                        sum(
                            (-1) ** (near * i + far * j) * U[src_par + 2 * i + 4 * j]
                            for i, j in np.ndindex(2, 2)
                        )
                        / 4
                    )
                    out = 2 * near + far
                    diagonal = np.sum(abs(block) ** 2, axis=0)
                    if n == 0:
                        read["empty"][ti, out] += diagonal.mean() / 2
                        continue
                    port_means = diagonal.reshape(len(nodes), 8).mean(axis=1) / 2
                    for p in nodes:
                        read[p][ti, out] += port_means[p]
                    b = block[:, inputs]
                    gram = (b.conj().T @ b).reshape(4, 4, 4, 4)
                    effect = np.trace(gram, axis1=1, axis2=3) / 8
                    for k, rho in sources.items():
                        read[k][ti, out] += np.trace(rho @ effect).real

    def mix(x, y):
        return (x + y) / 2

    laws = {"baseline": mix(read["empty"], read["paired"])}
    laws.update({k: mix(read["empty"], read[k]) for k in sources if k != "paired"})
    ports = [f"port_{p}" for p in nodes]
    laws.update(
        {k: mix(read["paired"], read[p]) for k, p in zip(ports, nodes, strict=True)}
    )
    laws["planted"] = mix(laws["baseline"], mix(laws[ports[0]], laws[ports[-1]]))
    np.testing.assert_allclose(
        laws["interior"],
        w * laws["matched"] + (1 - w) * laws["unmatched"],
        rtol=0,
        atol=1e-12,
    )
    for a in laws.values():
        assert a.min() >= -1e-13
        np.testing.assert_allclose(a.sum(axis=1), 1, rtol=0, atol=1e-12)
    return laws, ports, seam


def separator(load, endpoints, cone):
    # Scores have l1 norm at most one. A physical lottery includes zero load.
    # With free amplitude require opposite signs and maximize the smaller margin.
    size = load.shape[1]
    rows = []
    for v in endpoints:
        rows.append(np.r_[-v, v, 1.0, int(not cone)])
    for v in load:
        rows.append(np.r_[v, -v, int(cone), -int(not cone)])
    rows.append(np.r_[np.ones(2 * size), 0, 0])
    result = linprog(
        np.r_[np.zeros(2 * size), -1, 0],
        A_ub=rows,
        b_ub=[0] * (len(rows) - 1) + [1],
        bounds=[(0, None)] * (2 * size + 2),
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    return np.rint(100000 * (result.x[:size] - result.x[size : 2 * size])).astype(int)


def run():
    rng = np.random.default_rng(SEED)
    results = {}
    for stiffness in (8, 0):
        law, ports, seam = clock_laws(stiffness)
        raw = {
            k: np.array([rng.multinomial(SHOTS, row / row.sum()) for row in a])
            for k, a in law.items()
        }
        count = 2 * len(law) * len(TIMES) * 4
        eps = np.sqrt(np.log(2 * count / 0.01) / (2 * SHOTS))

        def changes(banks):
            return {k: (a - banks["baseline"]).ravel() for k, a in banks.items()}

        predicted = changes(law)
        observed = changes({k: a / SHOTS for k, a in raw.items()})
        archive = json.loads((HERE / "joint_records.json").read_bytes())
        shutter = [TIMES.index(t) for t in (2, 4)]
        pi = 2 / 33
        unit = (1 - pi) * law["baseline"] + pi * (
            law[ports[seam[0]]] + law[ports[seam[1]]]
        ) / 2
        deviations = []
        for k, a in (
            ("baseline", law["baseline"]),
            ("unit", unit),
            ("matched", law["matched"]),
            ("unmatched", law["unmatched"]),
        ):
            bank = np.asarray(archive["laws"][str(stiffness)]["counts"][k][0])
            deviations.extend(
                abs(a[shutter].mean(axis=0) - bank.sum(axis=2).ravel() / bank.sum())
            )
        archived_eps = np.sqrt(np.log(2 * 32 / 0.01) / (2 * archive["trials_per_bank"]))
        assert max(deviations) < archived_eps
        certificates = {}
        load = np.array([predicted[k] for k in ports])
        for band, (low, high) in {"all": (0.25, 1), "inward": (0.67, 0.69)}.items():
            weights = (4 * np.array([low, high]) - 1) / 3
            ends = np.array(
                [
                    w * predicted["matched"] + (1 - w) * predicted["unmatched"]
                    for w in weights
                ]
            )
            for cone in (True, False):
                coeff = separator(load, ends, cone)
                l1 = abs(coeff).sum()
                score = {k: float(coeff @ v / max(l1, 1)) for k, v in observed.items()}
                lower = (
                    min(
                        w * score["matched"] + (1 - w) * score["unmatched"]
                        for w in weights
                    )
                    - 2 * eps
                )
                upper = max(score[k] for k in ports) + 2 * eps
                passed = (lower > 0 and upper <= 0) if cone else lower > max(0, upper)
                certificates[f"{band}_{'cone' if cone else 'lottery'}"] = dict(
                    coefficients=coeff.reshape(len(TIMES), 4).tolist(),
                    l1=int(l1),
                    source_lower=float(lower),
                    max_port_upper=float(upper),
                    scores=score,
                    separated=bool(passed),
                )
        control = (observed[ports[0]] + observed[ports[-1]]) / 4
        planted_error = float(max(abs(observed["planted"] - control)))
        assert planted_error < 2 * eps
        results[str(stiffness)] = dict(
            counts={k: a.tolist() for k, a in raw.items()},
            certificates=certificates,
            planted_error=planted_error,
            planted_bound=2 * eps,
            shutter_max_error=max(deviations),
            archived_bound=archived_eps,
        )
    return dict(
        protocol_sha256=sha256(
            (HERE / "PORT_CLOCK_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        seed=SEED,
        shots=SHOTS,
        settings=TIMES,
        outcomes=["++", "+-", "-+", "--"],
        simultaneous_single_frequency_half_width=eps,
        laws=results,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    out = json.loads(json.dumps(run(), allow_nan=False))
    path = HERE / "port_clocks.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2) + "\n")
    else:
        joint.source.region.check(json.loads(path.read_bytes()), out)
    for k, result in out["laws"].items():
        print(
            k,
            {
                band: (c["source_lower"], c["max_port_upper"], c["separated"])
                for band, c in result["certificates"].items()
            },
        )
