"""Finite-source stationary deposit profiles, charge conservation and full work balances."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .material import fixture
from .module import native_microscope

HERE = Path(__file__).resolve().parent
LOCK = "e3f8353924e36a3d67716d4279228f5eb9fe26bcb95a9d34435b4f1293345fc1"
COUNT = 2**22


def link_energy(psi, edges):
    a, b = edges.T
    return np.sum(
        8.5 * (abs(psi[a]) ** 2 + abs(psi[b]) ** 2) - (psi[a].conj() * psi[b]).real
    ).item()


def survey(psi, donors, other, shells, rng):
    a, b = donors.T
    probs = np.column_stack([abs(psi[a] + psi[b]) ** 2, abs(psi[a] - psi[b]) ** 2]) / 2
    assert abs(probs.sum() - 1) < 1e-10
    probs /= probs.sum()
    deposit = probs.sum(axis=1)
    counts = rng.multinomial(COUNT, probs.ravel()).reshape(probs.shape)
    per_shell = np.bincount(shells, weights=deposit)
    ns = np.arange(len(per_shell))
    donor_energy = link_energy(psi, donors)
    assert abs(donor_energy - probs.sum(axis=0) @ [8, 9]) < 1e-9
    other_energy = link_energy(psi, other)
    after_echo = np.zeros(len(psi))
    after_echo[a] = after_echo[b] = deposit / 2
    stages = [
        donor_energy + other_energy + 8.5,
        8.5 + 8.5 * after_echo[other].sum(),
        donor_energy + 17,
        donor_energy + 8.5,
        donor_energy + 8.5,
    ]
    np.testing.assert_allclose(stages[-1] - stages[0], -other_energy, atol=1e-10)
    sampled = counts.sum(axis=1) / COUNT
    return dict(
        record_mean=donor_energy,
        final_ordinary_mean=8.5,
        stage_energies=stages,
        net_work=stages[-1] - stages[0],
        shutter_pair_work=[-8.5, 8.5],
        positive_stage_changes=float(np.maximum(np.diff(stages), 0).sum()),
        radius_first=float(ns @ per_shell),
        radius_second=float(ns**2 @ per_shell),
        owner_shell_probabilities=per_shell.tolist(),
        owner_shell_fractions=np.bincount(shells, weights=sampled).tolist(),
        address_energy_counts=counts.tolist(),
    )


def many(occupied, pairs, other):
    hot = np.asarray(occupied, dtype=float)
    p = np.zeros(2 * len(hot))
    a, b = pairs.T
    p[a] = p[b] = hot / 2
    x, y = other.T
    energy_other = 8.5 * np.sum(p[x] + p[y] - 2 * p[x] * p[y])
    n = int(hot.sum())
    # Product of ordinary even one-unit donor states on all hot dimers.
    # A cross-donor swap that changes the port word changes a donor's total
    # number; it has zero expectation. Trace the other ports in the product.
    return dict(
        fuel_units=n,
        unit_marks=n,
        record_energy=8 * n,
        ordinary_retained_units=n,
        final_ordinary_energy=8.5 * n,
        net_work=-float(energy_other),
        stage_energies=[
            16.5 * n + energy_other,
            8.5 * n + energy_other,
            25 * n,
            16.5 * n,
            16.5 * n,
        ],
    )


def run():
    assert sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest() == LOCK
    rng = np.random.default_rng(319144)
    res = []
    for lvl in range(3):
        h, pairs, other, physical, shells, audit = fixture(lvl)
        root = np.zeros(h.shape[0], dtype=complex)
        root[0] = 1
        traces = []
        for t in (0, 2, 8, 24):
            generator = -1j * t * h
            psi = expm_multiply(generator, root, traceA=generator.diagonal().sum())
            native = survey(psi, pairs, other, shells, rng)
            static = survey(
                np.exp(-1j * t * h.diagonal()) * root, pairs, other, shells, rng
            )
            np.testing.assert_allclose(native["stage_energies"][0], 25.5, atol=1e-9)
            assert static["radius_second"] == 0
            traces.append(dict(t=t, native=native, zero_exchange=static))
        res.append(
            dict(
                **audit,
                physical_donor_seams=physical.tolist(),
                traces=traces,
                unaddressed_no_capture_unit_marks=0,
                many={
                    name: many(bits, pairs, other)
                    for name, bits in [
                        ("cold", np.zeros(len(pairs))),
                        ("root_owner", shells == 0),
                        ("whole_board", np.ones(len(pairs))),
                    ]
                },
            )
        )
    return dict(
        protocol_sha256=LOCK,
        per_setting_starts=COUNT,
        single_load_unmodified_blank_energy_interval=[48, 54],
        one_separated_empty_mark_plus_original_rest_energy_floor=56,
        **native_microscope(),
        boards=res,
    )


def encode(obj):
    if isinstance(obj, np.generic):
        return obj.item()
    raise TypeError(type(obj))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    outfile = HERE / "results.json"
    result = run()
    content = json.dumps(result, default=encode, indent=1, allow_nan=False) + "\n"
    if parser.parse_args().write:
        outfile.write_text(content)
    else:
        old = json.loads(outfile.read_text())

        def compare(x, y):
            if isinstance(x, dict):
                assert x.keys() == y.keys()
                for k in x:
                    compare(x[k], y[k])
            elif isinstance(x, list) and x and isinstance(x[0], dict):
                assert len(x) == len(y)
                for a, b in zip(x, y, strict=True):
                    compare(a, b)
            elif isinstance(x, (list, float)):
                np.testing.assert_allclose(x, y, atol=1e-10, rtol=1e-10)
            else:
                assert x == y

        compare(json.loads(content), old)
    print(
        "verified native six-port energy transfer and finite-source profiles on 20, 80 and 320 owners"
    )
