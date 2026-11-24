"""Neutral-pair labels, conserved memory records and a coherent reversal test."""

import argparse
import json
from math import comb
from pathlib import Path

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import candidate, currents
from .. import registers
from ..written_junction import polarization


def local():
    a, b = registers.split(registers.native_generators()[0])
    head = np.zeros((11, 11))
    head[1:, 1:] = b / 2
    mat = np.diag([int(x == 0) + int(y == 0) for x, y in law.SRC])
    charge = np.kron(np.eye(11), mat) + np.kron(head, np.eye(9))
    np.testing.assert_array_equal(a, np.diag([2, 0, 0]))
    np.testing.assert_allclose(head @ head, head, atol=1e-15)
    hist = np.kron(np.eye(11), mat)
    start = np.zeros(9)
    start[[law.SRC.index((0, 2)), law.SRC.index((2, 0))]] = 1 / np.sqrt(2)
    balanced = np.eye(9)[law.SRC.index((1, 1))]
    difference = start - balanced
    clean = np.kron(np.eye(11), np.eye(9) - np.outer(difference, difference))
    echo = law.echo()
    h0, h1 = law.local_term_parts()

    def comm(x, y):
        return float(np.linalg.norm(x @ y - y @ x))

    diagrams = {}
    for name, gate in [("native", echo), ("no_mark", clean)]:
        half = (np.eye(99) - 1j * gate) / np.sqrt(2)
        initial = np.kron(np.eye(11)[0], start)
        v = half @ initial
        rho = v.reshape(11, 9).T @ v.reshape(11, 9).conj()
        material = np.stack((start, balanced))
        reduced = material @ rho @ material.T
        survived = float(reduced[0, 0].real)
        np.testing.assert_allclose(np.diag(reduced), [0.5, 0.5], atol=1e-12)
        back = half.conj().T
        erased = 0
        for proj in (head, np.eye(11) - head):
            state = np.kron(proj, np.eye(9)) @ v
            erased += abs(np.vdot(initial, back @ state)) ** 2
        target = (start - 1j * balanced) / np.sqrt(2)
        diagrams[name] = dict(
            commutator=comm(charge, gate),
            squared_change_from_initial_count=np.linalg.norm((charge - np.eye(99)) @ v)
            ** 2,
            surviving_pair=survived,
            offdiagonal_material=abs(reduced[0, 1]),
            coherent_material_read=float(np.vdot(target, rho @ target).real),
            full_reversal=abs(np.vdot(initial, back @ v)) ** 2,
            reversal_after_mark_read=erased,
        )
    return diagrams, [comm(charge, h) for h in (h0, h1, echo)], comm(hist, echo)


def run():
    checks, exact, bare = local()
    path, _, _ = polarization.layout()
    x = polarization.X
    aperture = abs(x) <= 18

    def score(s):
        return int((np.roll(s, 1) + s + np.roll(s, -1) - 4)[aperture].sum())

    out = []
    for shift in (1, -2):
        bulk, wall, _, _ = polarization.translated(path, shift)
        ref, source = 1 + bulk, 1 + wall
        sites = []
        for pairs in range(5):
            rows = []
            for survivors in range(pairs + 1):
                rows.append(
                    dict(
                        vacancies=survivors,
                        compensating_marks=pairs - survivors,
                        probability=comb(pairs, survivors) / 2**pairs,
                    )
                )
            out.append(
                dict(
                    shift=shift,
                    pairs=pairs,
                    addresses=list(sites),
                    aperture_units=(score(source) - score(ref)) / 3,
                    joint=rows,
                    reverse_after_read=checks["native"]["reversal_after_mark_read"]
                    ** pairs,
                )
            )
            if pairs == 4:
                break
            i = next(
                j
                for j in np.argsort(abs(x))
                if abs(x[j]) < 15 and source[j] == source[j + 1] == 1
            )
            sites.append([int(path[i]), int(path[i + 1])])
            for pattern in [(0, 2), (2, 0)]:
                check = source.copy()
                check[i : i + 2] = pattern
                assert score(check) == score(source)
            source[i : i + 2] = [0, 2]
    return currents._portable(
        dict(
            gates=checks,
            native_charge_commutators=exact,
            bare_vacancy_echo_commutator=bare,
            coupling=candidate.KAPPA,
            core_reads=out,
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    result = run()
    with args.output.open("x") as f:
        f.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
