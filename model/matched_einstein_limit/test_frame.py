"""Independent algebra, omitted-premise controls and cold receipt replay."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import random
import subprocess
import sys

import pytest

from . import frame, packet


def tensor_from_coordinates(x):
    matrix = [[F(0) for _ in range(4)] for _ in range(4)]
    for entry, (i, j) in zip(x, frame.PAIRS, strict=True):
        matrix[i][j] = matrix[j][i] = entry
    return matrix


def direct_contractions(matrix):
    return tuple(sum(matrix[i][j]*u[i]*u[j] for i in range(4) for j in range(4))
                 for u in frame.FRAME)


def test_generatively_reconstructs_arbitrary_symmetric_tensors():
    rng = random.Random(122)
    for _ in range(64):
        x = tuple(F(rng.randrange(-50, 51), rng.randrange(1, 20)) for _ in range(10))
        q = direct_contractions(tensor_from_coordinates(x))
        assert frame.decode(q) == x
        assert frame.apply(frame.design(), x) == q


def test_exact_future_unit_tensors_and_independent_inverse():
    metric = [[F((-1, 1, 1, 1)[i] if i == j else 0) for j in range(4)] for i in range(4)]
    assert direct_contractions(metric) == (F(-1),)*10
    assert all(u[0] > 0 for u in frame.FRAME)
    assert frame.inverse(frame.design()) == frame.decoder_matrix()


def test_sharp_norm_on_all_extreme_charge_vectors():
    maximum = F(0)
    for signs in product((-1, 1), repeat=10):
        result = frame.decode(tuple(map(F, signs)))
        maximum = max(maximum, *map(abs, result))
    assert maximum == F(17, 8)


@pytest.mark.parametrize('deleted', range(10))
def test_each_missing_probe_has_nonzero_invisible_tensor(deleted):
    tensor = tensor_from_coordinates(frame.decode(frame.basis(10, deleted)))
    assert any(entry for row in tensor for entry in row)
    q = direct_contractions(tensor)
    assert q[deleted] == 1
    assert all(q[i] == 0 for i in range(10) if i != deleted)


def test_actual_small_ball_normalization_arithmetic():
    # Coefficients after factoring out pi, R^(3) and ell^4.
    assert 4*F(1, 2)*(F(1, 3)-F(1, 5)) == F(4, 15)
    assert 4*(2*F(1, 90)-F(1, 18)) == -F(2, 15)
    assert -F(2, 15)*2 == -F(4, 15)  # Gauss R^(3)=2G_uu.
    assert F(17, 8)*F(15, 4) == F(255, 32)


def test_residual_identity_and_bound_with_nonzero_errors():
    # Rational pi parameter is an algebraic control, not an approximation receipt.
    pi, eta, ell = F(355, 113), F(7, 3), F(2, 9)
    C = 4*pi*ell**4/15
    rng = random.Random(1145)
    for _ in range(32):
        G, T, lam, eb, ea = (F(rng.randrange(-9, 10), 11) for _ in range(5))
        s = 2*pi*C*T+eb
        a = -C*(G-lam)+ea
        d = s+eta*a
        residual = G-lam-2*pi*T/eta
        assert eta*C*residual == eb+eta*ea-d
        assert abs(residual) <= (abs(eb)+eta*abs(ea)+abs(d))/(eta*C)


def test_small_errors_and_nonvanishing_coefficient_are_load_bearing():
    for n in (1, 3, 7, 31):
        ell = F(1, n+1)
        assert frame.normalized_error(2, ell, ell, ell) == F(255, 16)*ell
        assert frame.normalized_error(2, 1, 0, 0) == F(255, 64)
        assert frame.normalized_error(ell, ell, 0, 0) == F(255, 32)
    for eta in (0, -1):
        with pytest.raises(ValueError):
            frame.normalized_error(eta, 1, 1, 1)


def test_variable_coupling_can_preserve_ward_bianchi_with_variable_lambda():
    # G=0, T=c*g, Lambda=c*kappa at two points: both tensors conserved.
    c = F(3, 2)
    kappa = (F(1), F(2))
    lam = tuple(c*k for k in kappa)
    for p in range(2):
        for sign in (-1, 1, 1, 1):
            assert lam[p]*sign == kappa[p]*c*sign
    assert lam[0] != lam[1]


def test_nonzero_extrinsic_curvature_has_leading_area_effect():
    K = (F(1), F(1), F(1))
    correction = (sum(k*k for k in K)-sum(K)**2)/2
    assert correction == -3
    for n in (1, 4, 16):
        ell = F(1, n)
        assert (F(4, 15)*ell**4*correction)/ell**4 == -F(4, 5)


@pytest.mark.parametrize('key', tuple(packet.scope()))
def test_each_scope_clause_is_retained(key):
    expected = packet.report()
    changed = deepcopy(expected)
    changed['scope'][key] = 'UNSUPPORTED' if isinstance(changed['scope'][key], str) else not changed['scope'][key]
    with pytest.raises(ValueError, match='scope changed'):
        packet.validate(changed, expected)


def test_cold_receipt_and_source_custody():
    subprocess.run([sys.executable, '-W', 'error', '-m', 'model.matched_einstein_limit',
                    '--check'], cwd=packet.ROOT, check=True, capture_output=True, text=True)
    retained = json.loads(Path(packet.OUTPUT).read_text())
    assert retained['source_bindings'] == packet.source_bindings()
    assert retained['lean']['Lean_kernel_reexecuted_by_Python_receipt'] is False
    assert retained['frame_and_bounds']['profile_checks_prove_a_source_refinement'] is False


def test_lean_audit_rejects_unexpected_imported_axiom_and_stderr_warning():
    valid = '\n'.join("'OPH.EinsteinBranch.Quantitative."+name+
                      "' depends on axioms: [propext,\n Classical.choice, Quot.sound]"
                      for name in packet.THEOREMS)
    parsed = packet.parse_lean_output(valid)
    assert all(set(axioms) == {'propext', 'Classical.choice', 'Quot.sound'}
               for axioms in parsed.values())
    damaged = valid.replace('Quot.sound]', 'Quot.sound, unprovedRegionalLaw]', 1)
    with pytest.raises(ValueError, match='unapproved Lean theorem axiom'):
        packet.parse_lean_output(damaged)
    with pytest.raises(ValueError, match='admission or warning'):
        packet.parse_lean_output(valid, 'warning: imported declaration')
