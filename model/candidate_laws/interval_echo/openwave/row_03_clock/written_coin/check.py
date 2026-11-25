"""Exact native columns and a separate closed-form check of the count banks."""

import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.stats import binom

from model.candidate_laws.contrast_memory_spherical_v1 import law

HERE = Path(__file__).parent


def native_columns():
    # Use the integer native write and reversal tables, without the producer's
    # floating echo, diagonalization, propagation or Fourier transform.
    w = sp.Matrix(law._write_pattern().tolist()) / sp.sqrt(2)
    z = sp.zeros(law.N_SRC)
    u = z.row_join(w.T).col_join(w.row_join(sp.eye(w.rows) - w * w.T))
    r, d = law._label_ops()
    r, d = sp.Matrix(r.tolist()), sp.Matrix(d.tolist())
    basis = [law.SRC.index((1, 1))] + [
        (1 + law.LABELS.index((c, delta))) * law.N_SRC + law.SRC.index((1, 1))
        for c in (0, 1)
        for delta in (2, -2)
    ]
    j = sp.eye(law.DIM)[:, basis]
    proj = sp.Matrix([[1, -1], [-1, 1]]) / 2
    for mem, small in [((sp.eye(r.rows) - r) / 2, proj), (d, 4 * sp.eye(2))]:
        g = sp.diag(z, sp.kronecker_product(mem, sp.eye(law.N_SRC)))
        expected = sp.diag(0, sp.kronecker_product(proj, small))
        assert u * g * u * j == j * expected
    for sign in (1, -1):
        blank = (
            sp.eye(law.DIM)[:, law.SRC.index((2, 0))]
            + sign * sp.eye(law.DIM)[:, law.SRC.index((0, 2))]
        )
        assert u * blank / sp.sqrt(2) == j * sp.Matrix([0, 1, sign, 1, sign]) / 2


def check_sample(count, n, p):
    count = np.asarray(count)
    p = np.asarray(p)
    # An exact-zero theoretical bin cannot have a count. Protect tiny positive
    # sinc tails against finite precision on the analytic verification only.
    alpha = 0.01 / (2 * 18000)
    assert np.all(
        (count >= binom.ppf(alpha, n, p)) & (count <= binom.ppf(1 - alpha, n, p))
    )


def check():
    native_columns()
    data = json.loads((HERE / "records.json").read_text())
    n, m = data["starts"], data["timer_starts"]
    times = np.linspace(0.02, 12.02, 193)
    cols = np.arange(512)
    for tag, panel in data["panels"].items():
        strength, kind = tag.split("_")
        k = int(strength)
        even = {"even": 1.0, "odd": 0.0, "sharp": 0.5}[kind]
        rows = np.asarray(panel["counts"])
        assert np.all(rows[:, 1:].sum(axis=1) == n) and np.all(rows[:, -1] == 0)
        flip = (
            even * np.sin(2 * k * times) ** 2
            + (1 - even) * np.sin((4 * k + 1) * times / 2) ** 2
        )
        check_sample(rows[:, 2], n, flip)
        check_sample(rows[:, 0], m, np.sin(times / 16) ** 2)
        e = panel["energy"]
        assert sum(e["fourier_bins"]) == sum(e["local_blank_even_odd"]) == n
        check_sample(e["energy_increment_calibration"], m, np.sin(0.5) ** 2)
        check_sample(e["local_blank_even_odd"], n, [0.5, 0.5 * even, 0.5 * (1 - even)])
        qs = 0.5 * (cols == 0)
        for frequency, weight in ((4 * k, even / 2), (4 * k + 1, (1 - even) / 2)):
            offset = frequency / 16 / (2 * np.pi) - cols / 512
            qs += weight * (np.sinc(512 * offset) / np.sinc(offset)) ** 2
        check_sample(e["fourier_bins"], n, qs)
        if k and kind != "sharp":
            assert abs(panel["fit"]["frequency"] - e["peak"]) < 0.15
            assert panel["fit"]["held_rms"] < 0.035
        if k and kind == "sharp":
            assert panel["fit"]["held_rms"] > 0.035
    dephase = np.asarray(data["dephased_same_flip"])
    assert np.all(dephase.sum(axis=1) == n)
    check_sample(dephase[:, 1], n, 0.5)
    assert data["uniform_blank_sham_counts"] == [0, 0, n]
    print(
        "Exact invariant columns, both writes, 1,737 clock appointments, nine 512-bin energy banks and controls passed."
    )


if __name__ == "__main__":
    check()
