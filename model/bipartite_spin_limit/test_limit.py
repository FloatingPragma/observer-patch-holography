"""Independent finite controls for the analytic graph theorem."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
from scipy.linalg import eigh, expm

from tools import strict_json
from . import algebra, bounds, packet


def jordan_wigner(sites, edges, c, U):
    """Build the full tensor Fock action, without the producer's hop routine."""
    z = np.diag([1., -1.])
    annihilator = np.array([[0., 1.], [0., 0.]])
    modes = []
    for j in range(2*sites):
        operator = np.array([[1.]])
        for k in range(2*sites):
            operator = np.kron(operator, z if k < j else annihilator if k == j else np.eye(2))
        modes.append(operator)
    number = [a.T@a for a in modes]
    H = sum(U*number[2*x]@number[2*x+1] for x in range(sites))
    for x, y, sign in edges:
        for spin in (0, 1):
            hop = modes[2*x+spin].T@modes[2*y+spin]
            H += c*sign/2*(hop+hop.T)
    configurations = list(product((0, 1), repeat=2*sites))
    half = [i for i, bits in enumerate(configurations) if sum(bits) == sites]
    single = [i for i in half if all(sum(configurations[i][2*x:2*x+2]) == 1 for x in range(sites))]
    P = np.eye(len(half))[:, [half.index(i) for i in single]]
    H = H[np.ix_(half, half)]
    # Compute the independent leading spin action from actual virtual hopping.
    # The zero diagonal of H on P makes HP equal TP.
    HP = H@P
    effective = -HP.T@HP/U
    return H, P, effective


@pytest.mark.parametrize('sides', [(2, 1, 1), (3, 1, 1), (2, 2, 1), (2, 2, 2)])
def test_exact_parity_on_complete_half_filled_source(sides):
    result = algebra.exact_audit(sides)
    assert all(result['identities'].values())
    assert result['virtual_charge_states_retained']


@pytest.mark.parametrize('sites,edges,U', [
    (2, ((0, 1, 1),), 4),
    (3, ((0, 1, -1), (1, 2, 1)), 8),
    (4, ((0, 1, -1), (0, 2, 1), (1, 3, 1), (2, 3, 1)), 16),
    (4, ((0, 1, -1), (0, 2, 1), (1, 3, 1), (2, 3, 1)), 64),
])
def test_independent_full_fock_spectrum_canonical_generator_and_output(sites, edges, U):
    H, P, effective = jordan_wigner(sites, edges, 1, U)
    values, vectors = eigh(H)
    m = 2**sites
    projection = vectors[:, :m]@vectors[:, :m].T
    values_gram, basis = eigh(P.T@projection@P)
    W = projection@P@(basis@np.diag(values_gram**(-.5))@basis.T)
    certificate = bounds.finite_graph(1, U, len(edges), F(3, 2))
    assert np.linalg.norm(W.T@W-np.eye(m), 2) < 1e-12
    assert values[m] >= U-len(edges)-1e-11
    assert values[m-1] < 1e-11
    assert max(abs(values[:m]-eigh(effective, eigvals_only=True))) < float(certificate['ordered_energy_error_upper'])+1e-11
    assert np.linalg.norm(W.T@H@W-effective, 2) < float(certificate['canonical_generator_error_upper'])+1e-11
    actual = np.linalg.norm(expm(-1.5j*H)@P-P@expm(-1.5j*effective), 2)
    assert actual < float(certificate['complete_dynamics_error_upper'])+1e-11
    # Reversing the sign loses the physical second-order energy shift.
    if U >= 16*len(edges):
        assert np.linalg.norm(W.T@H@W+effective, 2) > float(certificate['canonical_generator_error_upper'])


def test_dimer_exact_quartic_correction_is_nonzero():
    # Exact low singlet: (U-sqrt(U^2+4c^2))/2; its difference from -c^2/U
    # is 4c^4/[U(U+sqrt(U^2+4c^2))^2], avoiding cancellation.
    for U in (8, 16, 64, 256):
        correction = 4/(U*(U+np.sqrt(U*U+4))**2)
        assert .96 < correction*U**3 < 1
        assert correction < float(bounds.finite_graph(1, U, 1)['ordered_energy_error_upper'])


def test_odd_cycle_is_outside_parity_domain():
    with pytest.raises(ValueError, match='bipartite'):
        algebra.bipartition(3, ((0, 1, 1), (1, 2, 1), (0, 2, 1)))


@pytest.mark.parametrize('c,U,edges,time', [(True, 8, 1, 0), (1., 8, 1, 0), (1, 3, 1, 0),
                                         (1, 8, True, 0), (1, 8, 0, 0), (1, 8, 1, -1)])
def test_invalid_or_uncertified_scales_rejected(c, U, edges, time):
    with pytest.raises(ValueError):
        bounds.finite_graph(c, U, edges, time)


def test_refinement_threshold_and_frozen_exchange_are_distinct():
    assert not bounds.refinement_exponents(F(16, 3))['global_error_proved_to_vanish']
    assert bounds.refinement_exponents(F(11, 2))['global_error_proved_to_vanish']
    assert bounds.refinement_exponents(F(11, 2))['exchange_power'] < 0
    assert not bounds.refinement_exponents(2)['global_error_proved_to_vanish']


def test_receipt_rejects_scope_bound_and_custody_mutations():
    value = strict_json.load(packet.OUTPUT)
    packet.validate(value)
    for key, changed in [('scope', {'integrated_gate_accepted': True}),
                         ('bounds', {}), ('source_bindings', {}), ('all_pass', 1)]:
        mutation = deepcopy(value)
        mutation[key] = changed
        with pytest.raises(ValueError):
            packet.validate(mutation)


def test_warm_import_cannot_certify_changed_source(monkeypatch):
    changed = dict(packet._IMPORTED_BINDINGS)
    changed['model/bipartite_spin_limit/algebra.py'] = '0'*64
    monkeypatch.setattr(packet, '_IMPORTED_BINDINGS', changed)
    with pytest.raises(ValueError, match='changed after import'):
        packet.report()
