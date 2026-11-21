"""Read an off-diagonal load-and-coin matrix element after unassisted flight."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .experiment import ARMS, HERE, canonical, local_check, sector


SHOTS = 1_048_576
TIMES = (0.25, 0.5, 1.0)


def terminal_gates():
    law = candidate.native_law
    h0, h1 = law.local_term_parts()
    kappa = candidate.KAPPA
    h = h0 + kappa * h1
    v = np.zeros((law.DIM, 2))
    for col, (coin, loads) in enumerate([(1, (1, 0)), (0, (0, 1))]):
        for d in (-1, 1):
            idx = (1 + law.LABELS.index((coin, d))) * law.N_SRC + law.SRC.index(loads)
            v[idx, col] = 1 / np.sqrt(2)
    gates = {}
    for name, angle in (("plus", np.pi / 2), ("minus", 3 * np.pi / 2)):
        # Remove the common phase of a spectator configuration.
        operated = np.exp(1j * angle) * expm(-1j * angle * h / kappa) @ v
        u = v.T @ operated
        prediction = np.eye(2) + (np.exp(1j * angle) - 1) * np.ones((2, 2)) / 2
        if np.linalg.norm(operated - v @ prediction) > 1e-11:
            raise ArithmeticError(
                "native full-column mixer does not implement the joint quadratures"
            )
        gates[name] = u
    ref = np.eye(law.DIM)[:, [law.SRC.index((0, 0)), law.SRC.index((1, 0))]]
    if np.linalg.norm(expm(-1j * np.pi * h1) @ ref - ref @ np.diag((1, -1))) > 1e-11:
        raise ArithmeticError(
            "target dephaser differs from the native blank-reference phase"
        )
    return gates


def run():
    local_check()
    gates = terminal_gates()
    car = federation(0)
    gateway = min(e for e in car.external if e[0] < 12)[1]
    targets = [
        q
        for q in car.neighbours[gateway]
        if q // 12 == gateway // 12 and car.degree[q] == 5
    ]
    ids = list(range(car.n))
    orders = [ids, sorted(ids, key=lambda p: (p // 12, (-p) % 12)), ids[::-1]]
    comparisons = len(orders) * len(TIMES) * len(targets) * 2
    margin = np.sqrt(np.log(2 * comparisons / 1e-6) / SHOTS)
    trials = []
    for arm, (label, order) in enumerate(zip(ARMS, orders, strict=True)):
        h, sites, inventory = sector(car, order)
        if inventory["omitted_states"]:
            raise ValueError(
                "the phase reader needs closure under every selected terminal contact"
            )
        offdiag = h.tocoo()
        source = np.zeros(len(sites))
        source[0] = 1.0
        for i, t in enumerate(TIMES):
            state = expm_multiply(-1j * t * h, source)
            for j, target in enumerate(targets):
                match = (sites[offdiag.row] == target) & (sites[offdiag.col] == gateway)
                x, y = offdiag.row[match], offdiag.col[match]
                if len(set(x)) != len(x) or len(set(y)) != len(y):
                    raise ArithmeticError(
                        "a joint coin state has two partners on the same seam"
                    )
                ordinary = float(np.sum(abs(state[sites == target]) ** 2))
                # States populate only their own sector. Both terminal contacts conserve it.
                spectator = float(
                    np.sum(
                        abs(
                            state[
                                (sites == target) & ~np.isin(np.arange(len(sites)), x)
                            ]
                        )
                        ** 2
                    )
                )
                result = {}
                for control, dephased in enumerate((False, True)):
                    readings = {}
                    for k, name in enumerate(("sham", "plus", "minus")):
                        if name == "sham":
                            prob = ordinary
                        else:
                            a, b = gates[name][0]
                            p0 = float(np.sum(abs(a * state[x] + b * state[y]) ** 2))
                            p1 = float(np.sum(abs(-a * state[x] + b * state[y]) ** 2))
                            prob = spectator + ((p0 + p1) / 2 if dephased else p0)
                        rng = np.random.default_rng([231310, arm, i, j, control, k])
                        clicks = int(rng.binomial(SHOTS, prob))
                        readings[name] = dict(
                            probability=prob, clicks=clicks, quiet=SHOTS - clicks
                        )
                    contrast = (
                        readings["plus"]["probability"]
                        - readings["minus"]["probability"]
                    )
                    predicted = (
                        0.0 if dephased else -2 * np.vdot(state[x], state[y]).imag
                    )
                    if abs(contrast - predicted) > 1e-10:
                        raise ArithmeticError(
                            "terminal frequencies do not select the expected joint quadrature"
                        )
                    observed = (
                        readings["plus"]["clicks"] - readings["minus"]["clicks"]
                    ) / SHOTS
                    result["dephased" if dephased else "coherent"] = dict(
                        settings=readings,
                        contrast=contrast,
                        observed_contrast=observed,
                        resolved_nonzero=bool(abs(observed) > margin),
                    )
                trials.append(
                    dict(
                        orientation=label,
                        t=t,
                        receiver=target,
                        gateway=gateway,
                        compatible_coin_links=len(x),
                        joint_overlap_imaginary=float(np.vdot(state[x], state[y]).imag),
                        arms=result,
                    )
                )
    return canonical(
        dict(
            candidate=candidate.NAME,
            scope="one load, 630 sharp even-sign retained collars",
            shots_per_setting=SHOTS,
            comparisons=comparisons,
            family_type_i_bound=1e-6,
            margin=margin,
            receiver_cell=gateway // 12,
            receiving_ports=targets,
            trials=trials,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="save a fresh receipt")
    output = parser.parse_args().output
    data = run()
    if output:
        with output.open("x") as f:
            json.dump(data, f, indent=2, allow_nan=False)
            f.write("\n")
    elif data != json.loads((HERE / "phase.json").read_text()):
        raise SystemExit("Joint phase receipt differs from archive")
    for label in ARMS:
        rows = [r for r in data["trials"] if r["orientation"] == label]
        print(
            label,
            sum(r["arms"]["coherent"]["resolved_nonzero"] for r in rows),
            "coherent contrasts resolved;",
            sum(r["arms"]["dephased"]["resolved_nonzero"] for r in rows),
            "dephased contrasts resolved",
        )
