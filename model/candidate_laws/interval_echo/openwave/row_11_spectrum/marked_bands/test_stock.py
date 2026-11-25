from math import pi

import numpy as np
from scipy.linalg import expm

from model.executable_model import federation
from model.experiments.material_inertia.dark_collar import label
from . import material, paired, stock


def test_each_real_pair_contributes_one_distinct_side_head():
    for bench, levels, number in (("ring", 2, 48), ("paired", 0, 12)):
        board = federation(levels)
        expected = material.geometry()[1] if bench == "ring" else paired.side_contacts()
        for sign in (-1, 1):
            row = stock.verify(bench, sign)
            deposits = row["stored_even_side_marks"]
            inputs = [p for x in deposits for p in x["stock"]]
            assert len(set(inputs)) == 2 * number
            assert row["initial_vacancies"] == row["initial_excesses"] == number + 1
            assert row["unspent_vacancies"] == row["unspent_excesses"] == 1
            assert [x["seam"] for x in deposits] == [sorted(e) for e in expected]
            assert {x["dwell_pi_over_64"] for x in deposits} == {1, 3}
            assert all(board.degree[x["reference"]] == 5 for x in deposits)
            assert row["private_zero_writes"] == (
                48
                if bench == "ring"
                else len(paired.rail.preparation()["silent_private"])
            )
            assert row["total_full_swaps"] > 2 * number


def test_reversing_the_phase_reference_requires_the_other_positive_dwell():
    law = material.LAW
    exchange, contrast = law.local_term_parts()
    half = expm(-1j * pi / 2 * (exchange + 8 * contrast)[:9, :9])[
        :, law.SRC.index((0, 2))
    ]
    target = label(2, 1)
    odd = label(2, -1)
    for endpoint, dwell in [(0, pi / 64), (1, 3 * pi / 64)]:
        for theta, head in [(dwell, target), (pi / 16 - dwell, odd)]:
            private = np.array([w[endpoint] ** 2 for w in law.SRC])
            out = (law.echo()[:, :9] @ (np.exp(-8j * theta * private) * half)).reshape(
                11, 9
            )
            overlap = abs(np.vdot(head, out[:, law.SRC.index((1, 1))]))
            assert abs(1 - overlap) < 3e-13
