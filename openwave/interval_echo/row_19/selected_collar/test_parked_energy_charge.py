"""Verify the sharp seam spectra and the degree bound independently."""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations_with_replacement

import numpy as np
import sympy as sp

from openwave.interval_echo.row_19.selected_collar import (
    parked_energy_charge as experiment,
)


def test_sharp_energies_with_vacuum_or_one_travelling_load():
    """Apply the native echo and generator over Q(sqrt(2)) on each target."""
    volume = experiment.reference.local.survey.routes.volume
    law = volume.candidate.native_law
    eye = sp.eye(law.N_SRC)
    W = sp.Matrix(law._write_pattern().tolist()) / sp.sqrt(2)
    reversal, square = (sp.Matrix(a.tolist()) for a in law._label_ops())
    memory = (sp.eye(law.N_LAB) - reversal) / 2 + sp.Rational(
        str(volume.candidate.KAPPA)
    ) * square
    generator = sp.kronecker_product(memory, eye)

    def echo(blank, written):
        read = W.T * written
        return read, W * blank + written - W * read

    for magnitude, sign, energy in ((None, 1, 0), (0, 1, 0), (2, 1, 32), (2, -1, 33)):
        letters = sp.zeros(law.N_LAB, 1)
        if magnitude is not None:
            for i, (_, contrast) in enumerate(law.LABELS):
                if abs(contrast) == magnitude:
                    letters[i] = sign if contrast < 0 else 1
            native = volume.label(magnitude, sign)
            normalized = np.asarray(letters / letters.norm(), dtype=float).ravel()
            np.testing.assert_allclose(native, [0, *normalized])
        for loads in ((0, 0), (1, 0), (0, 1), (2, 0), (0, 2), (1, 1)):
            endpoints = eye[:, law.SRC.index(loads)]
            blank = endpoints if magnitude is None else sp.zeros(law.N_SRC, 1)
            written = sp.kronecker_product(letters, endpoints)
            _, inner = echo(blank, written)
            output = echo(sp.zeros(law.N_SRC, 1), generator * inner)
            if magnitude is None:
                # A blank seam stays blank and only the mobile load moves.
                swapped = eye[:, law.SRC.index(loads[::-1])]
                kappa = sp.Rational(str(volume.candidate.KAPPA))
                mobile = (endpoints - swapped) / 2 + kappa * (
                    loads[0] - loads[1]
                ) ** 2 * endpoints
                assert output == (mobile, written)
            elif loads == (1, 1):
                assert output == (
                    sp.zeros(law.N_SRC, 1),
                    sp.zeros(law.N_SRC * law.N_LAB, 1),
                )
            else:
                assert output == (energy * blank, energy * written)


def test_no_cubic_and_exact_quadratic_class_read():
    result = experiment.reference.local.inputs("parked_energy_charge")
    for file, digest in result["inputs"].items():
        assert (
            sha256((experiment.reference.local.HERE / file).read_bytes()).hexdigest()
            == digest
        )
    symbols = np.array(result["sharp_spectra_in_calibration_order"], dtype=int)
    values = np.zeros((9, 2), dtype=object)
    for record in result["quadratic_class_read"]:
        monomial = symbols[:, record["factors"]].prod(axis=1)
        values += np.outer(monomial, list(map(Fraction, record["weights"])))
    np.testing.assert_array_equal(
        values, np.eye(3, dtype=int)[result["classes_in_calibration_order"]][:, 1:]
    )

    joined = np.hstack([np.repeat(symbols, 9, axis=0), np.tile(symbols, (9, 1))])
    delta = np.outer(result["w"], result["w"]).ravel()
    for degree in range(4):
        for factors in combinations_with_replacement(range(joined.shape[1]), degree):
            assert delta @ joined[:, factors].prod(axis=1) == 0
    quartic = joined[:, [0, 2, 4, 6]].prod(axis=1)
    assert delta @ quartic != 0
    assert result["charge_change_nats"]["fourth_order"]["exact"][0] > 0.00231337
    previous = experiment.reference.local.inputs("dust_source_read")["arms"]["shared"]
    current = result["charge_change_nats"]["shared"]["exact"]
    assert (
        previous["full_charge"][0]
        <= current[0]
        < current[1]
        <= previous["full_charge"][1]
    )
