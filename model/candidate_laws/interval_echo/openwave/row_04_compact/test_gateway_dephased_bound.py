from fractions import Fraction
import json

import numpy as np
import pytest
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

from model.executable_model import federation
from . import gateway_dephased_bound as bound
from .experiment import gateway_table, geometry
from .gateway_bound import SCORE, arbitrary_sweep_floor
from .gateway_confirmation import PROTOCOL
from .gateway_transport_bound import swap_home_ceiling


def test_continuous_dephased_window():
    midpoint, error, ceiling = bound.dephased_home_ceiling()
    assert midpoint == Fraction(2999545820205203, 4503599627370496)
    assert error < Fraction(208, 10**6)
    assert ceiling < Fraction(916448, 10**6)
    assert -(-ceiling.numerator * 10**12 // ceiling.denominator) == 916447948118
    assert list(bound.CENTERS_TWICE) == list(range(31, 76, 2))
    gates = json.loads(PROTOCOL.read_text(encoding="utf-8"), parse_float=Fraction)["gates"]
    charts = {tuple(row) for row in gateway_table(federation(0))}
    coherent = Fraction(min(arbitrary_sweep_floor(tips) for tips in charts), SCORE)
    assert coherent > gates["exact_joint_home_strictly_above"]
    assert coherent > gates["coherent_terminal_home_strictly_above"]
    assert coherent-ceiling >= gates["minimum_coherent_minus_dephased"]
    assert swap_home_ceiling()[2] < gates["stochastic_swap_home_strictly_below"]


def test_column_error_and_coherent_exception():
    board = federation(0)
    step = bound._exchange(board)
    # Compare an initial sharp column at the first half-integer center with
    # an unrounded complex Krylov action of the same time generator.
    re = np.eye(board.n, 1, dtype=np.int64)*bound.SCALE
    im = np.zeros_like(re)
    exact = expm_multiply(31j*step/4, re/bound.SCALE)
    for _ in range(31):
        re, im = bound._half_step(step, re, im)
    distance = np.linalg.norm((re+1j*im)/bound.SCALE-exact)
    assert distance < 31*float(bound._half_step_error(board.n))

    # Restoring phase between sharp ports can give home probability exactly
    # one at every exposure. The same interval cannot bound coherent inputs.
    _, antipodes = geometry()
    tips = set(gateway_table(board)[0])
    pairs = [pair for pair in antipodes if tips.isdisjoint(pair)]
    signed = np.zeros(board.n, dtype=np.int64)
    signed[list(pairs[0])] = 1
    signed[list(pairs[1])] = -1
    product = step @ signed
    eigenvalue = product[pairs[0][0]]
    assert np.array_equal(product, eigenvalue*signed)
    assert np.dot(signed[:12], signed[:12]) == np.dot(signed, signed) == 4


def test_nonuniform_phase_cannot_be_removed(monkeypatch):
    board = federation(0)
    native = bound.candidate.generator(board)
    change = diags([1, *[0]*(board.n-1)], format="csr", dtype=np.int64)
    monkeypatch.setattr(bound.candidate, "generator", lambda b: native+change)
    with pytest.raises(ArithmeticError, match="native one-load evolution"):
        bound._exchange(board)
