"""Observed first rejections and terminal addresses after two ordered projected currents."""

import argparse
import json
from functools import cache
from itertools import product
from math import log, pi, sqrt

import numpy as np

from model.experiments.interval_echo_consensus.signed_effects.sequential import (
    grid_counts,
)
from model.experiments.interval_echo_consensus.signed_effects.full_address import (
    measure as address,
)
from model.executable_model import federation
from .apparatus import LAW
from .code_frame import (
    RAILS,
    ACCEPT,
    CHARGES,
    FLUXES,
    N,
    evaluate,
    inverse,
    membership,
    phase_pair,
    target,
    logical_frame,
)
from .measure import HERE

SHOTS, SEED, PIECES = 262144, 962611, 16
SOURCES = (0, 3, 5, 8)


def keys():
    return product(range(2), SOURCES, range(2), (1, 3))


def instruments():
    targets = [target(k) for k in FLUXES]
    word = logical_frame()
    forward, back = evaluate(word), evaluate(inverse(word))
    mask, _ = membership()
    native = [
        ([back[s] @ mask[s] @ forward[s]], back[s] @ (np.eye(N) - mask[s]) @ forward[s])
        for s in range(2)
    ]
    supplied = []
    fine = []
    for s, (_, r, _, P) in enumerate(targets):
        p = np.zeros((N, N), complex)
        p[:12, :12] = P
        p[-1, -1] = int(r < 3)  # vacancy spare goes to the unused third pass slot
        passed = [np.outer(forward[s, a].conj(), forward[s, a]) for a in ACCEPT]
        supplied.append(([p], np.eye(N) - p))
        fine.append((passed, np.eye(N) - p))
    return native, supplied, fine


@cache
def terminal(sign):
    """Full address-pointer tape plus a signed head distinguishes spare from port zero."""
    calls = address.schedule()
    spare = RAILS[-1]
    witness = next(
        q
        for q in federation(0).neighbours[spare]
        if q not in {*RAILS, *address.RESERVED}
    )
    reading = np.zeros((N + 1, N))
    for col, site in enumerate(RAILS):
        loads, _ = address.execute_basis(calls, (-1, 1)[sign], site)
        raw = address.read_word(loads)
        assert loads[witness] == 1
        source = LAW.SRC.index((loads[spare], loads[witness]))
        for i in np.flatnonzero(LAW.write()[:, source]):
            _, difference = LAW.LABELS[i // 9]
            label = (
                (N - 1 if raw == 0 else N) if difference else (raw if raw < 12 else N)
            )
            reading[label, col] += abs(LAW.write()[i, source]) ** 2
    np.testing.assert_allclose(reading, np.eye(N + 1, N), rtol=0, atol=1e-14)
    return reading


def probabilities(key, instrument, gains, actual=False, reverse=False):
    s, source, order, angle = key
    # No measurement can distinguish a slot's common phase from a frame convention.
    response = tuple(-q * angle * pi / (6 * PIECES) for q in gains)
    z = phase_pair(*response)[s] if actual else np.exp(1j * response[s])
    yes, no = instrument[s]
    rho = np.diag((np.arange(N) == source).astype(complex))
    ports = (0, 1) if bool(order) == reverse else (1, 0)
    out = []
    for p in (None,) + tuple(port for port in ports for _ in range(PIECES)):
        if p is not None:
            rho[p, :] *= z
            rho[:, p] *= z.conjugate()
        out.append(np.trace(no @ rho @ no.conj().T).real)
        rho = sum(K @ rho @ K.conj().T for K in yes)
    effects = terminal(s) if actual else np.eye(N + 1, N)
    out.extend(effects @ rho.diagonal().real)  # includes spare and malformed records
    prob = np.maximum(0, np.array(out))
    np.testing.assert_allclose(prob.sum(), 1, atol=3e-10, rtol=0)
    return prob / prob.sum()


def read(raw, targets, fine):
    menus = {
        "declared": (targets, CHARGES, False),
        "unit_load": (targets, (-1, 1), False),
        "exchanged_gains": (targets, (2, -3), False),
        "zero": (targets, (0, 0), False),
        "reversed_order": (targets, CHARGES, True),
        "fine": (fine, CHARGES, False),
    }
    rows = {}
    for *k, c in raw:
        k = tuple(k)
        if (
            k in rows
            or len(c) != 47
            or any(type(v) is not int or v < 0 for v in c)
            or sum(c) != SHOTS
        ):
            raise ValueError("Invalid two-port enrollment")
        rows[k] = np.array(c) / SHOTS
    if set(rows) != set(keys()):
        raise ValueError("Incomplete two-port enrollment")
    errors = {
        name: max(
            float(max(abs(counts - probabilities(key, inst, g, reverse=rev))))
            for key, counts in rows.items()
        )
        for name, (inst, g, rev) in menus.items()
    }
    eps = sqrt(log(2 * len(rows) * 47 / 0.01) / (2 * SHOTS))
    assert (
        errors["declared"]
        < eps
        < min(value for name, value in errors.items() if name != "declared")
    )
    return dict(
        banks=len(rows),
        starts=SHOTS * len(rows),
        maximum_raw_frequency_errors=errors,
        simultaneous_frequency_radius=eps,
        malformed_or_spare_counts=sum(
            int((v[-1] + v[-2]) * SHOTS) for v in rows.values()
        ),
    )


def run(write=False):
    native, targets, fine = instruments()
    rng = np.random.default_rng(SEED)
    rows = []
    for key in keys():
        prob = probabilities(key, native, CHARGES, actual=True)
        rows.append([*key, grid_counts(rng, SHOTS, prob, tolerance=2.3e-11).tolist()])
    summary = read(rows, targets, fine)
    path = HERE / "code_records.json"
    if write:
        path.write_text(json.dumps(rows, separators=(",", ":")) + "\n")
    elif json.loads(path.read_text()) != rows:
        raise ValueError("Held two-port replay differs")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    print(json.dumps(run(parser.parse_args().write), indent=2))
